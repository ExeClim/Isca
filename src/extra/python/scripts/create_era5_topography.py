import numpy as np
import xarray as xr
import pyshtools as pysh
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
import copy, os
from scipy.special import j1
from resolutions import get_grid_for_truncation
from rich.console import Console
from rich.table import Table

# Physical constants
g = 9.80  # m/s^2, Earth gravity
a = 6376.0e3  # m, Earth radius

# Global console object for Rich output
console = Console()

def load_era5_invariant_fields(file_dir: str, file_lsm: str, file_z: str) -> xr.Dataset:
    """Load and merge ERA5 invariant fields (land-sea mask and topography)."""
    lsm_path = os.path.join(file_dir, file_lsm)
    z_path = os.path.join(file_dir, file_z)

    if not os.path.exists(lsm_path) or not os.path.exists(z_path):
        console.print("[bold red]Error:[/bold red] ERA5 invariant files not found.")
        console.print(f"Expected paths:\n- {lsm_path}\n- {z_path}")
        console.print("\n[bold]Note:[/bold] You must download the ERA5 invariant fields manually.")
        console.print("Please refer to the README.md for download instructions.\n")
        raise FileNotFoundError("ERA5 invariant files not found.")

    era5_lsm = xr.open_dataset(file_dir + file_lsm).squeeze()
    era5_z = xr.open_dataset(file_dir + file_z).squeeze()

    if era5_lsm.lsm.shape != era5_z.z.shape:
        raise ValueError(
            f"Shape mismatch between geopotential (z: {era5_z.z.shape}) "
            f"and land-sea mask (lsm: {era5_lsm.lsm.shape})"
        )

    invar = xr.merge([era5_lsm, era5_z], compat='no_conflicts')
    invar = invar.assign(zsurf=invar.z / g)
    invar = invar.drop_vars(['z'])

    return invar

def print_era5_statistics(ds: xr.Dataset) -> None:
    """Print statistics for ERA5 invariant fields."""
    nlat = len(ds['latitude'])
    nlon = len(ds['longitude'])

    console.print("[bold]ERA5 Invariant Fields Statistics[/bold]")
    console.print(f"[bold]Grid Shape (lat, lon):[/bold] {nlat} x {nlon}\n")

    table = Table(title="Variable Statistics", show_header=True, header_style="bold")
    table.add_column("Variable", style="cyan", no_wrap=True)
    table.add_column("Minimum")
    table.add_column("Maximum")
    table.add_row("zsurf", f"{ds.zsurf.min():.1f}", f"{ds.zsurf.max():.1f}")
    table.add_row("lsm", f"{ds.lsm.min():.1f}", f"{ds.lsm.max():.1f}")

    console.print(table)

def gaussian_grid(n_lat: int, n_lon: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a Gaussian latitude-longitude grid.
    Gaussian latitudes are the roots of the Legendre polynomial of degree n_lat,
    converted to degrees. Longitudes are equally spaced.
    """
    gauss_points, _ = leggauss(n_lat)
    lats = (np.arcsin(gauss_points) * (180.0 / np.pi))[::-1]
    lons = np.linspace(0, 360, n_lon, endpoint=False)
    return lats, lons

def dh_sh_filter(
    ds: xr.Dataset,
    tnum: int,
    n_target_lat: int = None,
    n_target_lon: int = None,
    gaussian: bool = False,
) -> xr.Dataset:
    """
    Perform spherical harmonic filtering on an input dataset and return fields
    on a specified latitude-longitude grid.
    """
    if tnum is None:
        raise ValueError("`tnum` must be provided as a positive, non-zero integer.")
    if not isinstance(tnum, int) or tnum <= 0:
        raise ValueError("`tnum` must be a positive, non-zero integer.")
    if (n_target_lat is None) != (n_target_lon is None):
        raise ValueError("`n_target_lat` and `n_target_lon` must be provided as a pair or not at all.")

    # If tnum is provided but n_target_lat and n_target_lon are not, use get_grid_for_truncation
    if tnum is not None and n_target_lat is None and n_target_lon is None:
        grid_params = get_grid_for_truncation(tnum)
        n_target_lat = grid_params['nlat']
        n_target_lon = grid_params['nlon']
        console.print(f"[bold]Note:[/bold] Using grid parameters for T{tnum}: nlat={n_target_lat}, nlon={n_target_lon}")

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
    filtered_ds = filtered_ds.astype({v: "float32" for v in filtered_ds.data_vars})
    filtered_ds = filtered_ds.rename({'latitude':'lat', 'longitude':'lon', 'lsm':'land_mask'})

    return filtered_ds

def generate_spectral_grids(ds: xr.Dataset, truncations: list[int]) -> None:
    """Generate and save spectral grids for specified truncations."""
    for t in truncations:
        gstat = get_grid_for_truncation(t)
        grid_name = f'era-spectral_T{gstat["nfou"]}_{gstat["nlat"]}x{gstat["nlon"]}'
        ds_tgrid = dh_sh_filter(ds, tnum=gstat["nfou"], gaussian=True)
        ds_tgrid.to_netcdf(f"{grid_name}.nc")
        console.print(f"Saved new T{gstat['nfou']} topography file [bold]{grid_name}.nc[/bold]\n")

def main():
    # File paths
    file_dir = 'era5_invariant_fields/'
    file_lsm = 'ecmwf-era5_oper_an_sfc_200001010000.lsm.inv.nc'
    file_z = 'ecmwf-era5_oper_an_sfc_200001010000.z.inv.nc'

    # Load and print statistics
    invar = load_era5_invariant_fields(file_dir, file_lsm, file_z)
    print_era5_statistics(invar)

    # Generate spectral grids
    truncations = [10, 21, 42, 85, 170, 341]
    generate_spectral_grids(invar, truncations)

if __name__ == "__main__":
    main()
