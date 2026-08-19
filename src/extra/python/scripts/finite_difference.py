"""
Routines associated with finite differencing on the sphere.
"""
import numpy as np

# Static parameters
RAD = np.pi / 180.0
EARTH_R = 6.371e6


class NDSlicer(object):
    """N-Dimensional slice class for numpy arrays."""

    def __init__(self, axis, ndim, start=None, stop=None, step=None):
        """
        Create an n-dimensional slice list.

        Parameters
        ----------
        axis : integer
            Axis on which to apply the slice
        ndim : integer
            Total number of dimensions of array to be sliced
        start, stop, step : integer, optional
            Index of beginning, stop and step width of the slice [start:stop:step]
            default for each is None.

        """
        self.axis = axis
        self.ndim = ndim
        self.start = start
        self.stop = stop
        self.step = step
        self.slicer = None
        self.__getitem__(slice(start, stop, step))

    def __getitem__(self, key):
        """
        Create an n-dimensional slice list.

        Parameters
        ----------
        axis : integer
            Axis on which to apply the slice
        ndim : integer
            Total number of dimensions of array to be sliced
        start, stop, step : integer, optional
            Index of beginning, stop and step width of the slice [start:stop:step]
            default for each is None.

        Returns
        -------
        slicer : list
            list of slices such that all data at other axes are kept, one axis is sliced

        Examples
        --------
        Create random array, slice it::

            x = np.random.randn(5, 3)

            # Create slicer equivalent to [1:-1, :]
            slc = NDSlicer(0, x.ndim)
            print(x)
            [[ 0.68470539  0.87880216 -0.45086367]
             [ 1.06804045  0.63094676 -0.76633033]
             [-1.69841915  0.35207064 -0.4582049 ]
             [-0.56431067  0.62833728 -0.04101542]
             [-0.02760744  2.02814338  0.13195714]]
            print(x[slc[1:-1]])
            [[ 1.06804045  0.63094676 -0.76633033]
             [-1.69841915  0.35207064 -0.4582049 ]
             [-0.56431067  0.62833728 -0.04101542]]

        """
        if isinstance(key, slice):
            self.start = key.start
            self.stop = key.stop
            self.step = key.step
        elif isinstance(key, int):
            self.start = key
            self.stop = key + 1
            self.step = None

        self.slicer = [slice(None)] * self.ndim
        self.slicer[self.axis] = slice(self.start, self.stop, self.step)
        return self.slicer

    def slice(self, start=None, stop=None, step=None):
        """Legacy compatibility method, calls `__getitem__`."""
        self.__getitem__(slice(start, stop, step))


def cfd(data, lons, lats, cyclic=False):
    """
    Vectorized central finite difference for N-D array data by lon / lat.

    Parameters
    ----------
    data : array_like
        N-dimensional array to be differentiated
    lons, lats : array_like
        1D coordinate arrays for longitude and latitude dimensions
    cyclic : Boolean, optional
        Data is cyclic on <axis> if true

    Returns
    -------
    diff_lon, diff_lat : array_like
        ND arrays of longitudinal, latitudinal centered finite differences of `data`
        respectively, with same dimensionality as `data`

    """
    # Find the axis where `data` shape matches longitudes
    axis_x = np.where(np.array(data.shape) == lons.shape[0])[0][0]

    # Find the axis where `data` shape matches latitudes
    axis_y = np.where(np.array(data.shape) == lats.shape[0])[0][0]

    dlong, dlatg = dlon_dlat(lons, lats, cyclic=cyclic)
    diff_x = diff_cfd(data, axis_x, cyclic=cyclic) / dlong
    diff_y = diff_cfd(data, axis_y, cyclic=False) / dlatg

    return diff_x, diff_y


def convert_radians_latlon(lat, lon):
    """
    Convert input lat/lon array to radians if input is degrees, do nothing if radians.

    Parameters
    ----------
    lat : array_like
        ND array of latitude
    lon : array_like
        ND array of longitude

    Returns
    ----------
    lat : array_like
        ND array of latitude in radians
    lon : array_like
        ND array of longitude in radians

    """
    if (np.max(np.abs(lat)) - np.pi / 2.0) > 1.0:
        lat_out = lat * RAD
    else:
        lat_out = lat

    if(np.min(lon) < 0 and np.max(lon) > 0 and
       np.abs(np.max(np.abs(lon)) - np.pi) > np.pi):
        lon_out = lon * RAD
    elif np.abs(np.max(np.abs(lon)) - np.pi * 2) > np.pi:
        lon_out = lon * RAD
    else:
        lon_out = lon

    return lat_out, lon_out


def dlon_dlat(lon, lat, cyclic=True):
    """
    Compute horizontal center finite differences of latitude/longitude on Earth [m].

    Parameters
    ----------
    lon, lat : array_like
        1D array of longitudes and latitudes respectively
    cyclic : boolean
        If True, longitudes are cyclic

    Returns
    -------
    dlong, dlatg : array_like
        Longitudinal, latitudinal centered finite differences in m
        Shapes are (lat.shape[0], lon.shape[0])

    """
    # Check that lat/lon are in radians
    lat, lon = convert_radians_latlon(lat, lon)

    # Calculate centre finite difference of lon / lat
    dlon = lon[2:] - lon[:-2]
    dlat = lat[2:] - lat[:-2]

    # If we want cyclic data, repeat dlon[0] and dlon[-1] at edges
    if cyclic:
        dlon = np.append(dlon[0], dlon)      # cyclic boundary in East
        dlon = np.append(dlon, dlon[-1])     # cyclic boundary in West
        _, lat2d = np.meshgrid(lon, lat)
    else:
        _, lat2d = np.meshgrid(lon[1:-1], lat)

    dlat = np.append(lat[1] - lat[0], dlat)    # boundary in South
    dlat = np.append(dlat, lat[-1] - lat[-2])  # boundary in North
    dlong, dlatg = np.meshgrid(dlon, dlat)

    # Lon/Lat differences in spherical coords
    dlong *= EARTH_R * np.cos(lat2d)
    dlatg *= EARTH_R

    return dlong, dlatg


def diff_cfd(data, axis=-1, cyclic=False, ndiff=1):
    """
    Calculate centered finite difference of a field along an axis with *even* spacing.

    Parameters
    ----------
    data : array_like
        ND array of data of which to calculate the differences
    axis : integer
        Axis of `data` on which differences are calculated
    cyclic : bool
        Flag to indicate whether `data` is cyclic on `axis`

    Returns
    -------
    diff : array_like
        ND array of central finite differences of `data` along `axis`

    Notes
    -----
    The output of this must be divided by 2x the grid spacing. For example, say `x` is
    the coordinate, and y is the `data`

        >>> x = np.array([-3, -2, -1, 0, 1, 2])
        >>> y = np.array([0.0, -0.9, -0.9, 0.0, 0.9, 0.9])
        >>> diff_y = diff_cfd(y, cyclic=True)
        >>> dydx = diff_y / (x[2] - x[0])

    Since the coordinate itself isn't cyclic, and the spacing is constant using
    the difference across a point works (shown is the difference across index == 1)

    Also, note that if data is *not* cyclic, then differences on the boundaries are
    forward / backward differences, thus the coordinate must be differenced in the
    same fashion (e.g. using this function on the coordinate variable)

    """
    # Calculate centred differences along longitude direction
    # Equivalent to: diff = data[..., 2:] - data[..., :-2] for axis == -1
    slc = NDSlicer(axis, data.ndim)
    if ndiff == 1:
        diff = data[slc[2:]] - data[slc[:-2]]
    elif ndiff == 2:
        diff = data[slc[2:]] - 2 * data[slc[1:-1]] + data[slc[:-2]]

    if cyclic:
        # Cyclic boundary in "East"
        if ndiff == 1:
            # Equiv to diff[..., 0] = data[..., 1:2] - data[..., -1:]
            d_1 = data[slc[1:2]] - data[slc[-1:]]
            # Cyclic boundary in "West"
            # Equiv to diff[..., -1] = data[..., 0:1] - data[..., -2:-1]
            d_2 = data[slc[0:1]] - data[slc[-2:-1]]

        elif ndiff == 2:
            d_1 = data[slc[1:2]] - 2 * data[slc[0:1]] + data[slc[-1:]]
            d_2 = data[slc[0:1]] - 2 * data[slc[-1:]] + data[slc[-2:-1]]

    else:
        if ndiff == 1:
            # Otherwise edges are forward/backward differences
            # Boundary in "South", (data[..., 1:2] - data[..., 0:1])
            d_1 = data[slc[1:2]] - data[slc[0:1]]

            # Boundary in "North" (data[..., -1:] - data[..., -2:-1])
            d_2 = data[slc[-1:]] - data[slc[-2:-1]]

        elif ndiff == 2:
            d_1 = data[slc[2:3]] - 2 * data[slc[1:2]] + data[slc[0:1]]
            d_2 = data[slc[-3:-2]] - 2 * data[slc[-2:-1]] + data[slc[-1:None]]

    diff = np.append(d_1, diff, axis=axis)
    diff = np.append(diff, d_2, axis=axis)

    return diff


def diffz(data, vcoord, axis=None):
    """
    Calculate vertical derivative for data on uneven vertical levels.

    Parameters
    ----------
    data : array_like
        N-D array of input data to be differentiated, where
        data.shape[axis] == vcoord.shape[0]
    vcoord : array_like
        Vertical coordinate, 1D
    axis : integer
        Axis where data.shape[axis] == vcoord.shape[0]

    Returns
    -------
    dxdz : array_like
        N-D array of d(data)/d(vcoord), same shape as input `data`

    """
    if axis is None:
        # Find matching axis between data and vcoord
        axis = np.where(np.array(data.shape) == vcoord.shape[0])[0][0]

    # Create array to hold vertical derivative
    dxdz = np.ones(data.shape)

    slc = NDSlicer(axis, data.ndim)
    # Create an n-dimensional broadcast along matching axis, same as [None, :, None, None]
    # for axis=1, ndim=4
    bcast = [np.newaxis] * data.ndim
    bcast[axis] = slice(None)

    dz1 = vcoord[1:-1] - vcoord[:-2]  # z[i] - z[i - 1]
    dz2 = vcoord[2:] - vcoord[1:-1]   # z[i + 1] - z[i]
    dz1 = dz1[bcast]
    dz2 = dz2[bcast]

    dxdz[slc[1:-1]] = ((dz1**2 * data[slc[2:]] + (dz2**2 - dz1**2) * data[slc[1:-1]] -
                        dz2**2 * data[slc[:-2]]) / (dz1 * dz2 * (dz2 + dz1)))

    # Do forward difference at 0th level [:, 1, :, :] - [:, 0, :, :]
    i = 0
    dz1 = vcoord[i + 1] - vcoord[i]
    dz2 = vcoord[i + 2] - vcoord[i + 1]
    dxdz[slc[i]] = (-dz1**2 * data[slc[i + 2]] + (dz1 + dz2)**2 * data[slc[i + 1]]
                    - (dz2**2 + 2 * dz1 * dz2) * data[slc[i]]) / (dz1 * dz2 * (dz1 + dz2))

    # Do backward difference at Nth level [:, -1, :, :] - [:, -2, :, :]
    i = data.shape[axis] - 1
    dz1 = vcoord[i - 1] - vcoord[i - 2]
    dz2 = vcoord[i] - vcoord[i - 1]
    dxdz[slc[i]] = (((dz1**2 + 2 * dz1 * dz2) * data[slc[i]] -
                     (dz1 + dz2)**2 * data[slc[i - 1]] + dz2**2 * data[slc[i - 2]]) /
                    (dz1 * dz2 * (dz1 + dz2)))

    return dxdz


def diff2z(data, vcoord, axis=None):
    """
    Calculate 2nd order vertical derivative for data on uneven vertical levels.

    Parameters
    ----------
    data : array_like
        N-D array of input data to be differentiated, where
        data.shape[axis] == vcoord.shape[0]
    vcoord : array_like
        Vertical coordinate, 1D
    axis : integer
        Axis where data.shape[axis] == vcoord.shape[0]

    Returns
    -------
    d2xdz2 : array_like
        N-D array of d**2(data)/d(vcoord)**2, same shape as input `data`

    """
    if axis is None:
        # Find matching axis between data and vcoord
        axis = np.where(np.array(data.shape) == vcoord.shape[0])[0][0]

    # Create array to hold vertical derivative
    d2xdz2 = np.ones(data.shape)

    slc = NDSlicer(axis, data.ndim)
    # Create an n-dimensional broadcast along matching axis,
    # this is the same as [None, :, None, None] for axis=1, ndim=4
    bcast = [np.newaxis] * data.ndim
    bcast[axis] = slice(None)

    dz1 = vcoord[1:-1] - vcoord[:-2]  # z[i] - z[i - 1]
    dz2 = vcoord[2:] - vcoord[1:-1]   # z[i + 1] - z[i]
    dz1 = dz1[bcast]
    dz2 = dz2[bcast]

    d2xdz2[slc[1:-1]] = (2 * (dz1 * data[slc[2:]] - (dz1 + dz2) * data[slc[1:-1]]
                              + dz2 * data[slc[:-2]]) / (dz1 * dz2 * (dz1 + dz2)))

    i = 0
    dz1 = vcoord[i + 1] - vcoord[i]
    dz2 = vcoord[i + 2] - vcoord[i + 1]
    d2xdz2[slc[i]] = 2 * (dz1 * data[slc[i + 2]] - (dz1 + dz2) * data[slc[i + 1]]
                          + dz2 * data[slc[i]]) / (dz1 * dz2 * (dz1 + dz2))

    i = data.shape[axis] - 1
    dz1 = vcoord[i - 1] - vcoord[i - 2]
    dz2 = vcoord[i] - vcoord[i - 1]
    d2xdz2[slc[i]] = 2 * (dz1 * data[slc[i]] - (dz1 + dz2) * data[slc[i - 1]]
                          + dz2 * data[slc[i - 2]]) / (dz1 * dz2 * (dz1 + dz2))

    return d2xdz2


def dp(pres):
    """
    Calculate centered finite difference on pressure field with evenly spaced levels.

    Note
    ----
    Difference will be negative if input pressure is monotonic decreasing
    (surface is 0th element)

    Parameters
    ----------
    pres : array_like
        1D array of pressure level data

    Returns
    -------
    dp : array_like
        1D array of same shape as `pres`, that is the centered finite difference for
        elements [1:-1] and forward/backward difference for elements [0] and [-1] resp.

    """
    d_p = pres[2:] - pres[:-2]
    d_p = np.append(pres[1] - pres[0], dp)    # boundary at surface
    d_p = np.append(dp, pres[-1] - pres[-2])  # boundary aloft

    return dp