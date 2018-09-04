from cftime import utime
from datetime import  datetime
from cmip_time import FakeDT
import numpy as np
import pdb

__author__='Stephen Thomson'

def day_number_to_datetime_array(time_in, calendar_type, units_in):

    cdftime = utime(units_in, calendar = calendar_type)

    date_out = cdftime.num2date(time_in)

    return date_out

def day_number_to_date(time_in, calendar_type = '360_day', units_in = 'days since 0001-01-01 00:00:00'):
    """
    Aim is to make the time array have attributes like .month, or .year etc. This doesn't work with
    normal datetime objects, so Mike's FakeDT does this for you. First step is to turn input times
    into an array of datetime objects, and then FakeDT makes the array have the attributes of the
    elements themselves.
    """

    time_in = day_number_to_datetime_array(time_in, calendar_type, units_in)

    cdftime = FakeDT( time_in, units=units_in,
                 calendar=calendar_type)

    return cdftime

def month_to_season(months_in, avg_or_daily):

    seasons=np.zeros(len(months_in))
    idx_djf=(months_in==1)|(months_in==2)|(months_in==12)    
    idx_mam=(months_in==3)|(months_in==4)|(months_in==5)
    idx_jja=(months_in==6)|(months_in==7)|(months_in==8)
    idx_son=(months_in==9)|(months_in==10)|(months_in==11)

    seasons[idx_djf]=0
    seasons[idx_mam]=1
    seasons[idx_jja]=2
    seasons[idx_son]=3

    return seasons


def month_to_two_months(months_in, avg_or_daily):

    two_months=np.zeros(len(months_in))
    idx_jf=(months_in==1)|(months_in==2)
    idx_ma=(months_in==3)|(months_in==4)
    idx_mj=(months_in==5)|(months_in==6)
    idx_ja=(months_in==7)|(months_in==8)
    idx_so=(months_in==9)|(months_in==10)
    idx_nd=(months_in==11)|(months_in==12)
        
    two_months[idx_jf]=0
    two_months[idx_ma]=1
    two_months[idx_mj]=2
    two_months[idx_ja]=3
    two_months[idx_so]=4
    two_months[idx_nd]=5    

    return two_months

def recurring_to_sequential(time_in):

    seq_time=np.zeros_like(time_in)

    for time in np.arange(len(seq_time)-1)+1:        
        if time_in[time]!=time_in[time-1]:
            seq_time[time]=seq_time[time-1]+1
        else:
            seq_time[time]=seq_time[time-1]

    return seq_time
        


if __name__ == "__main__":
    import numpy as np
    from datetime import  datetime

    cdftime = utime('hours since 0001-01-01 00:00:00')
    date = datetime.now()
    print(date)
    t = cdftime.date2num(date)
    print(t)
    date = cdftime.num2date(t)
    print(date)

    time = 1772.5
    date_equivalent = day_number_to_date(time)
    print(date_equivalent)





