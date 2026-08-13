# -*- coding: utf-8 -*-
"""Creates an input file for hs_forcing_nml's 'from_file' local_heating_option.

By default this generates a file where the heating rate is zero everywhere - a
starting point that a user can edit (see the 'DEFINE A HEATING RATE HERE' comment
below) to prescribe their own local heating field.
"""
import os

import numpy as np
from netCDF4 import Dataset

import create_timeseries as cts

__author__ = 'Stephen Thomson'


def create_zero_heating_rate_file(output_dir, file_name='heating_rate.nc', num_levels=25, surface_pressure_hpa=1000., manual_grid_option=False):
    """Writes a 12-month climatology of a local heating rate (K/s) to
    output_dir/file_name, on the horizontal grid selected by manual_grid_option
    (see create_timeseries.create_grid) and on num_levels evenly-spaced pressure
    levels up to surface_pressure_hpa.

    The heating rate itself is all zeros - edit the array below this docstring
    to prescribe a non-zero heating field."""

    lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(manual_grid_option)

    p_full = np.linspace(0., surface_pressure_hpa, num_levels + 1)[1:]

    ntime = 12
    time_arr = np.arange(ntime) * 30.  # 12 months, 30 days apart, matching Isca's thirty_day calendar

    # DEFINE A HEATING RATE HERE. Units are K/s. Defaults to zero everywhere.
    heating_rate = np.zeros((ntime, num_levels, nlat, nlon))

    file_path = os.path.join(output_dir, file_name)
    output_file = Dataset(file_path, 'w', format='NETCDF3_CLASSIC')

    output_file.createDimension('lat', nlat)
    output_file.createDimension('lon', nlon)
    output_file.createDimension('latb', nlatb)
    output_file.createDimension('lonb', nlonb)
    output_file.createDimension('pfull', num_levels)
    output_file.createDimension('time', 0)  # unlimited time axis

    latitudes = output_file.createVariable('lat', 'd', ('lat',))
    longitudes = output_file.createVariable('lon', 'd', ('lon',))
    latitudebs = output_file.createVariable('latb', 'd', ('latb',))
    longitudebs = output_file.createVariable('lonb', 'd', ('lonb',))
    pfulls = output_file.createVariable('pfull', 'd', ('pfull',))
    times = output_file.createVariable('time', 'd', ('time',))

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

    times.units = 'days since 0000-01-01 00:00:00.0'
    times.calendar = 'THIRTY_DAY_MONTHS'
    times.calendar_type = 'THIRTY_DAY_MONTHS'
    times.cartesian_axis = 'T'

    heating_rate_out = output_file.createVariable('heating_rate', 'f4', ('time', 'pfull', 'lat', 'lon'))

    latitudes[:] = lats
    longitudes[:] = lons
    latitudebs[:] = latbs
    longitudebs[:] = lonbs
    pfulls[:] = p_full
    times[:] = time_arr
    heating_rate_out[:] = heating_rate

    output_file.close()

    return file_path


if __name__ == "__main__":
    create_zero_heating_rate_file(os.getcwd())
