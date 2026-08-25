# -*- coding: utf-8 -*-
"""Creates an input file for hs_forcing_nml's 'from_file' local_heating_option,
combining a polar-cap heating perturbation with a zonal-wavenumber
midlatitude heating perturbation, for use with equilibrium_t_option='Polvani_Kushner'.

Ported from Regan Mudhar's polar_heating()/heat_perturb()/combo_heat1()
functions in rgnmudhar/Isca_rmudhar polvani_kushner_rm:input/input_files.py,
adapted to build the horizontal grid from the repo's standard T42 grid file
and the vertical grid from the repo's ozone climatology file (rather than
hardcoded arrays), so no new binary files need to be committed.
"""
import os

import numpy as np
from netCDF4 import Dataset

__author__ = 'Regan Mudhar'


def _get_grid(gfdl_base):
    grid_file = Dataset(os.path.join(gfdl_base, 'src/extra/python/scripts/gfdl_grid_files/t42.nc'), 'r', format='NETCDF3_CLASSIC')
    lons = grid_file.variables['lon'][:]
    lats = grid_file.variables['lat'][:]
    lonbs = grid_file.variables['lonb'][:]
    latbs = grid_file.variables['latb'][:]

    ozone_file = Dataset(os.path.join(gfdl_base, 'input/rrtm_input_files/ozone_1990_notime.nc'), 'r', format='NETCDF3_CLASSIC')
    p_full = ozone_file.variables['pfull'][:]
    p_half = ozone_file.variables['phalf'][:]

    return lats, lons, latbs, lonbs, p_full, p_half


def polar_heating_field(lats, lons, p_full, y_wid=15., th_mag=4., p_th=50., p_top=600., p_ref=800.):
    """Polar-cap heating perturbation (K/s), following Orlanski and Solman (2010).

    y_wid: latitudinal decay of heating away from the pole (degrees)
    th_mag: magnitude of forcing (K/day) when centred on the p_ref level
    p_th: sets the vertical gradient of the forcing at its cap (hPa)
    p_top: depth of forcing - top pressure level below which heating decays to zero (hPa)
    p_ref: reference pressure level defining the forcing magnitude (hPa)
    """
    heat_lat = np.exp(-((lats - 90.) / y_wid) ** 2.)  # (lat,)
    heat_vert = 0.5 * th_mag * (1000. - p_ref) / (1000. - p_top) * (1. + np.tanh((p_full - p_top) / p_th)) / 86400.  # (pfull,)

    field = heat_vert[:, np.newaxis, np.newaxis] * heat_lat[np.newaxis, :, np.newaxis] * np.ones((1, 1, len(lons)))
    return field  # (pfull, lat, lon)


def midlat_wave_heating_field(lats, lons, p_full, q_0=6., m=2, y_cen=45., p_0=800., p_t=200.):
    """Zonal-wavenumber-m tropospheric diabatic heating perturbation (K/s) to
    induce NH winter-like wave activity, following Lindgren et al. (2018).

    q_0: magnitude of forcing (K/day)
    m: zonal wavenumber
    y_cen: latitudinal centre of the heating (degrees)
    p_0, p_t: lower and upper pressure bounds of the forcing (hPa)
    """
    y_wid = 0.175 * 360. / (2. * np.pi)

    heat_lat = np.exp(-0.5 * ((lats - y_cen) / y_wid) ** 2.)  # (lat,)
    heat_lon = np.cos(m * np.deg2rad(lons))  # (lon,)
    heat_p = np.sin(np.pi * np.log(p_full / p_0) / np.log(p_t / p_0))  # (pfull,)

    field = (q_0 / 86400.) * heat_p[:, np.newaxis, np.newaxis] * heat_lat[np.newaxis, :, np.newaxis] * heat_lon[np.newaxis, np.newaxis, :]
    field[(p_full < p_t) | (p_full > p_0), :, :] = 0.

    return field  # (pfull, lat, lon)


def create_polar_and_midlat_heating_file(output_dir, file_name='w15a4p600f800g50_q6m2y45l800u200.nc',
                                          polar_kwargs=None, midlat_kwargs=None):
    """Writes the sum of the polar-cap and midlatitude wavenumber-m heating
    fields (K/s) to output_dir/file_name, for use as hs_forcing_nml's
    local_heating_file with local_heating_option='from_file'.

    Reproduces Regan Mudhar's combo_heat1() with its default parameters
    (a polar heating perturbation combined with a zonal wavenumber-2 midlatitude
    heating perturbation centred at 45N) unless overridden via polar_kwargs/midlat_kwargs."""

    gfdl_base = os.environ['GFDL_BASE']
    lats, lons, latbs, lonbs, p_full, p_half = _get_grid(gfdl_base)

    polar = polar_heating_field(lats, lons, p_full, **(polar_kwargs or {}))
    midlat = midlat_wave_heating_field(lats, lons, p_full, **(midlat_kwargs or {}))
    heating_field = polar + midlat

    file_path = os.path.join(output_dir, file_name)
    output_file = Dataset(file_path, 'w', format='NETCDF3_CLASSIC')

    output_file.createDimension('lat', len(lats))
    output_file.createDimension('lon', len(lons))
    output_file.createDimension('latb', len(latbs))
    output_file.createDimension('lonb', len(lonbs))
    output_file.createDimension('pfull', len(p_full))
    output_file.createDimension('phalf', len(p_half))

    latitudes = output_file.createVariable('lat', 'd', ('lat',))
    longitudes = output_file.createVariable('lon', 'd', ('lon',))
    latitudebs = output_file.createVariable('latb', 'd', ('latb',))
    longitudebs = output_file.createVariable('lonb', 'd', ('lonb',))
    pfulls = output_file.createVariable('pfull', 'd', ('pfull',))
    phalfs = output_file.createVariable('phalf', 'd', ('phalf',))

    latitudes.units = 'degrees_N'
    latitudes.cartesian_axis = 'Y'
    latitudes.long_name = 'latitude'
    latitudes.edges = 'latb'

    longitudes.units = 'degrees_E'
    longitudes.cartesian_axis = 'X'
    longitudes.long_name = 'longitude'
    longitudes.edges = 'lonb'

    latitudebs.units = 'degrees_N'
    latitudebs.cartesian_axis = 'Y'
    latitudebs.long_name = 'latitude edges'

    longitudebs.units = 'degrees_E'
    longitudebs.cartesian_axis = 'X'
    longitudebs.long_name = 'longitude edges'

    pfulls.units = 'hPa'
    pfulls.cartesian_axis = 'Z'
    pfulls.positive = 'down'
    pfulls.long_name = 'full pressure level'

    phalfs.units = 'hPa'
    phalfs.cartesian_axis = 'Z'
    phalfs.positive = 'down'
    phalfs.long_name = 'half pressure level'

    var_name = os.path.splitext(file_name)[0]
    heating_out = output_file.createVariable(var_name, 'f4', ('pfull', 'lat', 'lon'))

    latitudes[:] = lats
    longitudes[:] = lons
    latitudebs[:] = latbs
    longitudebs[:] = lonbs
    pfulls[:] = p_full
    phalfs[:] = p_half
    heating_out[:] = heating_field

    output_file.close()

    return file_path, var_name


if __name__ == "__main__":
    create_polar_and_midlat_heating_file(os.getcwd())
