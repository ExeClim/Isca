# -*- coding: utf-8 -*-
"""
Module for handling of CMIP5 date/time, wraps complicated stuff with fake and/or
real datetime/DatetimeIndex objects
"""
import netCDF4 as nc
import datetime as dt
import pandas as pd
import numpy as np


class FakeDT(object):
    """
    An object created to mimic the behavior of a pandas DatetimeIndex object, but
    one that allows for dates from non-standard calendars (e.g. 360 day or no leap)
    Parameters
    ----------
    dates : array_like
            Array or list of dates to be converted into a FakeDT object
    units : string
            Time units (e.g. hours since...)
    calendar : string
            Calendar to which <dates> belong
    """

    def __init__(self, dates, units='hours since 1800-01-01 00:00:00',
                 calendar='standard'):

        self.dates = np.array(dates)
        self.units = units
        self.calendar = calendar
        self.ndates = len(self.dates)
        self.dtype = type(dates[0])
        if self.ndates == 1:
            self.year = dates.year
            self.month = dates.month
            self.day = dates.day
            self.hour = dates.hour
            self.minute = dates.minute
            try:
                    self.dayofyear = dates.timetuple().tm_yday
            except AttributeError:
                    self.dayofyear = dates.timetuple()[7]
        else:
            self.year = np.array([dk.year for dk in dates])
            self.month = np.array([dk.month for dk in dates])
            self.day = np.array([dk.day for dk in dates])
            self.hour = np.array([dk.hour for dk in dates])
            self.minute = np.array([dk.minute for dk in dates])
            try:
                    self.dayofyear = np.array([dk.timetuple().tm_yday for dk in dates])
            except AttributeError:
                    self.dayofyear = np.array([dk.timetuple()[7] for dk in dates])

    def __getitem__(self, idx):
        # If <idx> is array_like, return a new FakeDT object restricted to those
        # indicies, if not, just return the member at a particular location
        if isinstance(idx, (list, np.ma.MaskedArray, np.ndarray)):
            return FakeDT(self.dates[idx], self.units, self.calendar)
        else:
            return self.dates[idx]

    def __str__(self):
        if self.ndates == 1:
            return "[ {}, dtype={} ]".format(self.dates, type(self.dates))
        else:
            out_s = "[ "
            for k in range(self.ndates - 1):
                if k % 5 == 0:
                    out_s += "{},\n".format(self.dates[k])
                else:
                    out_s += "{}, ".format(self.dates[k])
            out_s += "{}, dtype={} ]".format(self.dates[-1], type(self.dates))
        return out_s

    def __reduce__(self):
        """ Special method for pickle to output in binary format """
        return (self.__class__, (self.dates, self.units, self.calendar))

    def __len__(self):
        return self.ndates

    def get_loc(self, date):
        """
        FakeDT class method for returning the index of a particular date
        raises KeyError if the date is not found. Uses bisection method
        Parameters
        ----------
        date : scalar_like
                netcdftime.datetime or datetime.datetime date for which to search
        Returns
        -------
        c : scalar_like
                Index of <date> in <self.dates>
        """
        a, b = 0, len(self.dates)-1
        niter = 0
        # Compare dates using the .timetuple() method, since this works if <date> is
        # a datetime or netcdftime .datetime, otherwise only == works, not > or <
        while True and niter < len(self.dates):
            if self.dates[a] == date:
                return a
            elif self.dates[b] == date:
                return b
            elif self.dates[a].timetuple() < date.timetuple()\
                    and self.dates[b].timetuple() > date.timetuple():
                c = a + (b-a)/2
                if self.dates[c] == date:
                    return c
                elif self.dates[c].timetuple() > date.timetuple():
                    b = c
                elif self.dates[c].timetuple() < date.timetuple():
                    a = c
            else:
                # First error string only raised if 'c' has been assigned
                if 'c' in locals():
                    raise KeyError('Date not found {}, a({}): {}, b({}): {}, '
                                   'c({}):{}'.format(date, a, self.dates[a],
                                                     b, self.dates[b],
                                                     c, self.dates[c]))
                else:
                    raise KeyError('Date not found {}, a({}): {},'
                                   ' b({}): {}'.format(date, a, self.dates[a],
                                                       b, self.dates[b]))
            niter += 1
        return c


def num2date_wrap(intimes):
    """
    Wrapper method for num2date, since num2date doesn't return a
    real datetime object if the calendar is non-standard
    (e.g. 360 day or noleap)
    Parameters
    ----------
        intimes : array_like
                netCDF4 time array, has .units and .calendar attributes
    Returns
    -------
        dates_out : array_like
                Pandas DatetimeIndex array, if calendar is standard or a FakeDT
                array if the calendar is non-standard
    """

    from netCDF4 import num2date
    n2d = num2date(intimes[:], intimes.units, calendar=intimes.calendar)
    try:
        return pd.DatetimeIndex(n2d)
    except (AttributeError, TypeError):
        return FakeDT(n2d, intimes.units, intimes.calendar)


def add_timedelta(date, delta, units='hours since 1800-01-01 00:00',
                  calendar='standard'):
    """
    Wrapper method to add a timedelta to a netcdftime.datetime or datetime.datetime
    or pd.Timestap object

    Parameters
    ----------
    date : scalar_like
            Date to add <delta> time to
    delta : scalar_like
            datetime.timedelta object to be added to <date>
    units : string
            Time units, only relavent if <date> is netcdf.datetime object
    calendar : string
            Calendar of <date>, only relavent if <date> is netcdf.datetime object

    Returns
    -------
    date + delta : scalar_like
            <date> + <delta> as either a datetime.datetime or netcdf.datetime object
    """
    if isinstance(date, pd.Timestamp) or isinstance(date, dt.datetime):
        return date + delta
    elif isinstance(date, nc.netcdftime.datetime):
        tot_hrs = delta.total_seconds()/(60.0**2)  # Calc number of hours in delta
        time = nc.date2num(date, units, calendar)  # Convert date to hours since...
        # Returns a real datetime object if calendar is std, otherwise it's netcdftime
        return nc.num2date(time + tot_hrs, units, calendar)
    else:
        raise TypeError("No method for date {} type {}".format(date, type(date)))


def sub_timedelta(date, delta, units='hours since 1800-01-01 00:00',
                  calendar='standard'):
    """
    Wrapper method to subtract a timedelta to a netcdftime.datetime
    or datetime.datetime or pd.Timestap object

    Parameters
    ----------
    date : scalar_like
            Date to subtract <delta> time from
    delta : scalar_like
            datetime.timedelta object to be subracted from <date>
    units : string
            Time units, only relavent if <date> is netcdf.datetime object
    calendar : string
            Calendar of <date>, only relavent if <date> is netcdf.datetime object

    Returns
    -------
    date - delta : scalar_like
            <date> - <delta> as either a datetime.datetime or netcdf.datetime object
    """

    if isinstance(date, pd.Timestamp) or isinstance(date, dt.datetime):
        return date - delta
    elif isinstance(date, nc.netcdftime.datetime):
        tot_hrs = delta.total_seconds()/(60.0**2)  # See add_timedelta above
        time = nc.date2num(date, units, calendar)
        return nc.num2date(time - tot_hrs, units, calendar)
    else:
        raise TypeError("No method for date {} type {}".format(date, type(date)))


def sub_ncdate(date1, date0, units='hours since 1900-01-01 00:00',
               calendar='standard'):
    """
    Wrapper method to add two dates together for two netcdftime.datetime
    or datetime.datetime or pd.Timestap objects

    Parameters
    ----------
    date1 : scalar_like
            Date to subtract <date0> from
    date0 : scalar_like
            Date to subtract from <date1>
    units : string
            Time units, only relavent if <date0> and <date1> are
            netcdf.datetime objects
    calendar : string
            Calendar of <date>, only relavent if <date0> and <date1> are
            netcdf.datetime objects

    Returns
    -------
    date1 - date0 : scalar_like
            <date0> + <date1> as a pandas.Timedelta object
    """

    if isinstance(date0, (pd.Timestamp, pd.DatetimeIndex)) and \
       isinstance(date1, (pd.Timestamp, pd.DatetimeIndex)):
        return date1 - date0
    else:
        # Try each time conversion seperately, that way if it fails, we know which one
        try:
            time0 = nc.date2num(date0, units, calendar)
        except (ValueError, AttributeError) as e:
            print('ERROR {}'.format(e))
            raise KeyError(date0)

        try:
            time1 = nc.date2num(date1, units, calendar)
        except (ValueError, AttributeError) as e:
            print('ERROR {}'.format(e))
            raise KeyError(date1)

        # Determine units of dates, since TimedeltaIndex
        # needs this to output correctly
        if 'hours' in units:
            out_units = 'h'
        elif 'minutes' in units:
            out_units = 'm'
        elif 'days' in units:
            out_units = 'd'

        # Try to output as a 'list' (maybe input dates are FakeDT objects?) if not
        # fall back on outputing a single pandas Timedelta, if it fails because
        # the date is after 2262-04-11 or April 11th, 2262, since pandas stores dates
        # in nanosecond resolution, so after this it overflows, return a
        # datetime.timedelta object, after re-converting so that units are in days
        try:
            return pd.TimedeltaIndex((time1 - time0), unit=out_units)

        except ValueError:
            return pd.Timedelta((time1 - time0), unit=out_units)

        except OverflowError:
            units_new = 'days since 2000-01-01 00:00'
            time0 = nc.date2num(date0, units_new, calendar)
            time1 = nc.date2num(date1, units_new, calendar)
            return dt.timedelta(time1-time0)
