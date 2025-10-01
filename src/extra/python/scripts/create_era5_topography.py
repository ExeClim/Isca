import numpy as np
import xarray as xr
import pyshtools as pysh
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
from scipy.special import j1
import copy


# Physical constants
g = 9.80 # m/s^2, Earth gravity
a = 6376.0e3 # m, Earth radius

# You will need to download invariant ERA-5 fields (such as topography) either
# directly from the CDS archive or from CEDA (login required):
#     https://catalogue.ceda.ac.uk/uuid/2c8f38fac04945b89cf12d6e9c928c6f/
#
# The following files will be required (if downloading from CEDA):
#     Land-sea mask: ecmwf-era5_oper_an_sfc_200001010000.lsm.inv
#     Topography:    ecmwf-era5_oper_an_sfc_200001010000.z.inv

file_dir = 'era5_source_files/'
file_lsm = 'ecmwf-era5_oper_an_sfc_200001010000.lsm.inv.nc'
file_z = 'ecmwf-era5_oper_an_sfc_200001010000.z.inv.nc'

# Open files and drop time dimension as not required
era5_lsm = xr.open_dataset(file_dir + file_lsm).squeeze()
era5_z = xr.open_dataset(file_dir + file_z).squeeze()

if era5_lsm.lsm.shape != era5_z.z.shape:
    raise ValueError(
        f"Shape mismatch between geopotential (z: {era5_z.z.shape}) "
        f"and land-sea mask (lsm: {era5_lsm.lsm.shape})"
    )

invar = xr.merge([era5_lsm, era5_z], compat='no_conflicts')
invar = invar.assign(
    zsurf=invar.z / g,
)
invar = invar.drop_vars(['z'])

# Print statistics
print("Dataset Statistics:")
nlat = len(invar['latitude'])
nlon = len(invar['longitude'])
print(f"Shape (lat, lon): {nlat} x {nlon}")
print()
print(f"{'Variable':<6} {'Min':>8} {'Max':>8}")
print("-" * 26)
print(f"{'z':<10} {invar.zsurf.min():6.1f}   {invar.zsurf.max():6.1f}")
print(f"{'lsm':<10} {invar.lsm.min():6.1f}   {invar.lsm.max():6.1f}")
print()


def gaussian_grid(n_lat, n_lon):
    """
    Compute a Gaussian latitude-longitude grid.

    Gaussian latitudes are the roots of the Legendre polynomial of degree n_lat,
    converted to degrees. Longitudes are equally spaced.

    Parameters
    ----------
    n_lat : int
        Number of Gaussian latitude points.
    n_lon : int
        Number of longitude points.

    Returns
    -------
    lats : np.ndarray
        Gaussian latitudes in degrees, ordered from north to south.
    lons : np.ndarray
        Regular longitudes in degrees, ranging from 0 to 360 (exclusive).
    """
    gauss_points, _ = leggauss(n_lat)
    lats = (np.arcsin(gauss_points) * (180.0 / np.pi))[::-1]
    lons = np.linspace(0,360,n_lon, endpoint=False)
    return lats, lons

def dh_sh_filter(ds, tnum=None, n_target_lat=None, n_target_lon=None, gaussian=False):
"""
    Perform spherical harmonic filtering on an input dataset and return fields
    on a specified latitude-longitude grid.

    The input dataset is first interpolated to a Driscoll-Healy (DH) grid,
    expanded into spherical harmonics, coefficients are filtered, and then the
    field is reconstructed and interpolated onto the user-specified target grid.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset with dimensions 'latitude' (ascending or descending) and
        'longitude'. Can contain multiple variables.
    tnum : int, optional
        Maximum spherical harmonic degree/order for the expansion.
        Determines truncation and effective filtering strength.
    n_target_lat : int, optional
        Number of latitude points for the output grid. If None, must be supplied
        when `gaussian=True`. For regular grids, defaults to the largest even
        number ≤ input latitude count.
    n_target_lon : int, optional
        Number of longitude points for the output grid. If None, defaults to
        input longitude count.
    gaussian : bool, optional
        If True, output grid will use Gaussian quadrature latitudes with
        `n_target_lat` points. If False (default), output grid will be
        regularly spaced in latitude and longitude.

    Returns
    -------
    filtered_ds : xarray.Dataset
        Dataset of filtered fields on the specified target grid. Output has
        dimensions ('lat', 'lon'). If the input contains a land-sea mask variable,
        it is rounded to a binarised mask'.
    """
    nlat_orig = len(ds.latitude)
    nlon_orig = len(ds.longitude)
    # pyshtools expects input grids to match certain formats, one of which is
    # a DH grid, which has shape (x, x/2) rather than the supplied (x, x/2 + 1)
    lats_dh = np.linspace(90, -90, (nlat_orig // 2) * 2)
    lons_dh = np.linspace(0, 360, nlon_orig, endpoint=False)

    # Define target lat-lon arrays (can be uniform or Gaussian)
    if gaussian:
        target_lats, target_lons = gaussian_grid(n_target_lat, n_target_lon)
    else:
        target_lats = np.linspace(90, -90, n_target_lat)
        target_lons = np.linspace(0, 360, n_target_lon, endpoint=False)

    # Ensure latitude ascending for xarray interp
    ds_sorted = ds.sortby('latitude')

    filtered_dict = {}
    for var in ds_sorted.data_vars:
        da = ds_sorted[var]
        # Interpolate input to DH grid (roughly matching input resolution)
        da_interp = da.interp(
            latitude=xr.DataArray(lats_dh, dims='latitude'),
            longitude=xr.DataArray(lons_dh, dims='longitude'),
            method='linear'
        )

        # Create SHGrid object and expand coefficients upto tnum
        data_dh = da_interp.values
        grid = pysh.SHGrid.from_array(data_dh)
        clm = grid.expand(lmax_calc=tnum)

        # Apply filter
        clm_filtered = copy.deepcopy(clm)
        Theta_opt = 3.8317 / (tnum + 0.5)
        for l in range(tnum + 1):
            for m in range(l + 1):
                factor = 2 * j1(l * Theta_opt) / (l * Theta_opt) if l > 0 else 1.0
                clm_filtered.coeffs[0, l, m] *= factor
                clm_filtered.coeffs[1, l, m] *= factor

        # Expand back to DH grid (PySh decides shape based on tnum)
        filtered_grid = clm_filtered.expand()
        filtered_array = filtered_grid.to_array()
        lat_filtered = np.linspace(90, -90, filtered_array.shape[0])
        lon_filtered = np.linspace(0, 360, filtered_array.shape[1], endpoint=False)

        # Wrap as xarray DataArray and interpolate to target grid
        da_filtered = xr.DataArray(filtered_array,
                              dims=('latitude', 'longitude'),
                              coords={'latitude': lat_filtered,
                                      'longitude': lon_filtered})
        da_filtered_interp = da_filtered.interp(latitude=target_lats,
                                                longitude=target_lons,
                                                method='linear')

        # Round values in land mask
        if var == 'lsm':
            da_filtered_interp = np.rint(da_filtered_interp)

        filtered_dict[var] = (('latitude', 'longitude'), da_filtered_interp.values)

    coords_dict = {'latitude': target_lats, 'longitude': target_lons}
    filtered_ds = xr.Dataset(filtered_dict, coords=coords_dict)

    # Sort latitude ascending and rename
    filtered_ds = filtered_ds.sortby('latitude', ascending=True)
    filtered_ds = filtered_ds.rename({'latitude':'lat', 'longitude':'lon',
                                      'lsm':'land_mask'})

    return filtered_ds



#ds_t42 = dh_sh_filter(invar, tnum=42, n_target_lat=64, n_target_lon=128, gaussian=True)
#ds_t42.to_netcdf('era-spectral_T42_64x128.nc')
