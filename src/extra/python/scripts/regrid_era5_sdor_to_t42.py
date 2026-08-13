"""
Regrids ERA5 "standard deviation of sub-gridscale orography" (sdor) to Isca's
T42 Gaussian grid (64 lat x 128 lon, same grid as era-spectral7_T42_64x128.out.nc
elsewhere in this repo), producing a direct mg_drag "ghprime" input field.
Writes exp/test_cases/mg_drag/input/mg_drag.res.nc, the input file used by
exp/test_cases/mg_drag/mg_drag_socrates_topo_test_case.py.

Source data: ERA5 single-level invariant parameter "Standard deviation of
orography" (sdor, ECMWF param 160), downloaded from the Copernicus Climate
Data Store (CDS) -- see https://cds.climate.copernicus.eu, reanalysis-era5-single-levels,
variable "standard_deviation_of_orography". Not committed to git (large binary,
regeneratable from CDS) -- edit ERA5_SDOR_PATH below to point at your own copy.

Why 'input' mode instead of 'computed': mg_drag_nml's source_of_sgsmtn='computed'
path (src/atmos_param/mg_drag/mg_drag.f90 -> src/shared/topography/topography.F90)
expects *raw high-resolution elevation* and computes subgrid standard deviation
itself by binning many fine-grid points into each destination grid cell. ERA5's
sdor is already a precomputed standard-deviation field (ECMWF derives it the same
way, from a much finer DEM), so feeding it through 'computed' would compute the
"stdev of a stdev field" -- meaningless. Instead this script regrids sdor
directly onto the model grid and mg_drag reads it verbatim via
source_of_sgsmtn='input' (mg_drag_init's INPUT/mg_drag.res.nc restart-style
read, variable name 'ghprime').

Regridding: bilinear interpolation onto Isca's Gaussian T42 latitudes (via
numpy.polynomial.legendre.leggauss). Interpolation is done by hand with two
separable 1D numpy.interp passes (longitude then latitude) rather than
xarray's/scipy's interp -- this works cleanly here because both the ERA5
source grid and the T42 target grid are regular in longitude and monotonic in
latitude, so a separable linear interpolation is exact for bilinear
interpolation on a rectilinear grid. Spherical-harmonic filtering is
deliberately not used: sdor is a standard-deviation field and must stay >= 0
everywhere, and spectral truncation can ring negative near sharp gradients
(e.g. the Himalayas); linear interpolation of a non-negative field stays
non-negative.

Latitude convention: mg_drag_init reads this file via FMS's read_data(), which
loads the raw array by index position into the model's internal field with no
coordinate-based matching -- it just trusts the array is already ordered to
match the model's own domain-decomposed grid, exactly like every other Isca
input file (era-spectral7_T42_64x128.out.nc is ascending, -90 -> 90). Writing
descending here would silently mirror the whole field north-south when
mg_drag reads it at runtime, so gaussian_grid() below returns ascending
latitudes to match.
"""
import os

import numpy as np
import xarray as xr
from numpy.polynomial.legendre import leggauss

from isca import GFDL_BASE

ERA5_SDOR_PATH = '/home/links/sit204/gwd_tests/69eb30085e940a6d80b1ce305b746b0b.nc'
OUT_PATH = os.path.join(GFDL_BASE, 'exp/test_cases/mg_drag/input/mg_drag.res.nc')

NLAT_T42 = 64
NLON_T42 = 128


def gaussian_grid(n_lat: int, n_lon: int) -> tuple[np.ndarray, np.ndarray]:
    """Ascending (-90 -> 90) Gaussian latitudes, matching Isca's own grid
    convention -- see the 'Latitude convention' note above for why this matters."""
    gauss_points, _ = leggauss(n_lat)
    lats = np.arcsin(gauss_points) * (180.0 / np.pi)
    lons = np.linspace(0, 360, n_lon, endpoint=False)
    return lats, lons


def bilinear_regrid(data, lat_src, lon_src, lat_tgt, lon_tgt):
    """Separable linear regrid (no scipy). lat_src must be ascending;
    lon_src must be ascending and span >= [lon_tgt.min(), lon_tgt.max()]."""
    # interpolate along longitude for every source latitude row
    by_lon = np.empty((data.shape[0], len(lon_tgt)))
    for i in range(data.shape[0]):
        by_lon[i, :] = np.interp(lon_tgt, lon_src, data[i, :])
    # then interpolate along latitude for every target longitude column
    out = np.empty((len(lat_tgt), len(lon_tgt)))
    for j in range(len(lon_tgt)):
        out[:, j] = np.interp(lat_tgt, lat_src, by_lon[:, j])
    return out


def main():
    era5 = xr.open_dataset(ERA5_SDOR_PATH).squeeze()
    sdor = era5['sdor'].sortby('latitude')  # ERA5 ships lat descending; need ascending for np.interp

    lat_src = sdor['latitude'].values
    lon_src = sdor['longitude'].values
    data_src = sdor.values

    target_lat, target_lon = gaussian_grid(NLAT_T42, NLON_T42)
    regridded = bilinear_regrid(data_src, lat_src, lon_src, target_lat, target_lon)

    assert regridded.shape == (NLAT_T42, NLON_T42)
    assert regridded.min() >= 0.0, 'regridded sdor went negative -- unexpected for linear interp of a >=0 field'

    ghprime = xr.DataArray(
        regridded.astype('float64'),
        dims=('lat', 'lon'),
        coords={'lat': target_lat, 'lon': target_lon},
        name='ghprime',
        attrs={'long_name': 'subgrid-scale orography standard deviation (from ERA5 sdor)', 'units': 'm'},
    )
    ghprime.to_netcdf(OUT_PATH)
    print(f'Wrote {OUT_PATH}: ghprime min={float(ghprime.min()):.1f} max={float(ghprime.max()):.1f} '
          f'mean={float(ghprime.mean()):.2f} m, shape={ghprime.shape}')


if __name__ == '__main__':
    main()
