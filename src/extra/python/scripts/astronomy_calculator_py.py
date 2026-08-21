import numpy as np
import copy as cp
import matplotlib.pyplot as plt
import pyshtools as pysh
import pdb

def daily_mean_solar(time_since_ae, lat):
    """ Function to calculate the daily-mean incoming solar radiation cos(zen) as a function of the time-of year.
        Taken from updated Isca astronomy module on 9th Jan 2017 with commit-id 46988a7.
    """
    
    PI = np.pi

    ang = angle (time_since_ae)
    dec = declination(ang)
    h   = half_day    (lat, dec)

    cosz = np.sin(lat)*np.sin(dec)*(h)/(PI) + np.cos(lat)*np.cos(dec)*np.sin(h)/(PI)
    cosz[np.where(cosz < 0.)] = 0.
    
    return cosz
    
def half_day (latitude, dec):

    PI = np.pi
    eps = 1.0E-05
    tan_dec = np.tan(dec)

    latitude[np.where(latitude ==  0.5*PI)]  = latitude[np.where(latitude ==  0.5*PI)] - eps
    latitude[np.where(latitude ==  -0.5*PI)] = latitude[np.where(latitude ==  -0.5*PI)] + eps

    cos_half_day = -np.tan(latitude)*tan_dec
    h = np.arccos(cos_half_day)
    h[np.where(cos_half_day <= -1.0)] = PI
    h[np.where(cos_half_day >= +1.0)] = 0.0
    
    return h


def orbit():

    num_angles=3600.
    ecc = 0.0
    twopi = 2.* np.pi

    orb_angle = np.zeros(int(num_angles+1))

    orb_angle[0] = 0.0
    dt = twopi/float(num_angles)
    norm = np.sqrt(1.0 - ecc**2)
    dt = dt*norm

    for n in range(1,int(num_angles+1)):
        d1 = dt*r_inv_squared(orb_angle[n-1])
        d2 = dt*r_inv_squared(orb_angle[n-1]+0.5*d1)
        d3 = dt*r_inv_squared(orb_angle[n-1]+0.5*d2)
        d4 = dt*r_inv_squared(orb_angle[n-1]+d3)
        d5 = d1/6.0 + d2/3.0 +d3/3.0 +d4/6.0
        orb_angle[n] = orb_angle[n-1] + d5

    return orb_angle

def angle(t):

    num_angles=3600.
    twopi = np.pi * 2.
    
    norm_time = t*num_angles/twopi
    
    norm_time_keep = cp.copy(norm_time)
    
    norm_time_floor = np.floor(norm_time)
    norm_time_floor_mod = np.mod(norm_time_floor,num_angles)
    
    x = norm_time_keep - np.floor(norm_time_keep)
    
    orb_angle = orbit()
    
    angle = (1.0 -x)*orb_angle[int(norm_time_floor_mod)] + x*orb_angle[int(norm_time_floor_mod+1)]
    angle = np.mod(angle, twopi)

    return angle
    
def r_inv_squared (ang):
    """
    !--------------------------------------------------------------------
    !    define the earth-sun distance (r) and then return the inverse of
    !    its square (r_inv_squared) to the calling routine.
    !--------------------------------------------------------------------
    """
    
    deg_to_rad = np.pi / 180.
    
    per   = 102.932
    ecc   = 0.

    rad_per       = per*deg_to_rad
    r             = (1. - ecc**2)/(1. + ecc*np.cos(ang - rad_per))
    r_inv_squared = r**(-2)

    return r_inv_squared    



def declination (ang):

    obliq = 23.439
    deg_to_rad = np.pi / 180.

    rad_obliq   =   obliq*deg_to_rad
    sin_dec     = - np.sin(rad_obliq)*np.sin(ang)
    declination =   np.arcsin(sin_dec)

    return declination

def area_average(in_arr, lat_arr):

    area_average = 0.5 * np.mean(in_arr * np.cos(lat_arr)* np.pi) #Has factor of pi as doing a mean means summing and dividing by number of elements. Whereas doing it properly would be done by summing over the latitude axis with delta_lat = np.pi / nlat. So to convert the first way to be like the second way, we multiply by np.pi, as dividing by nlat in the mean will be like multiplying by delta_lat.
    
    return area_average

def frierson_sw(lat_arr, del_sol = 1.4, del_sw = 0.):
    """Include Frierson incoming-solar radiation formulation for comparison"""
    
    
    p2     = (1. - 3.*(np.sin(lat_arr)**2.))/4.        
    s = 0.25 * (1.0 + del_sol*p2 + del_sw * np.sin(lat_arr))

    return s

def mod_frierson_sw(lat_arr, p2_coeff, p4_coeff, p6_coeff):
    
    p2     = (1. - 3.*(np.sin(lat_arr)**2.))/2.        
    p4     = (1./8.)*(35.*(np.sin(lat_arr)**4.) - 30.*(np.sin(lat_arr)**2.) + 3.)
    p6     = (1./16.)*(231.*(np.sin(lat_arr)**6.) - 315.*(np.sin(lat_arr)**4.) + 105.*(np.sin(lat_arr)**2.) - 5.)
    
    s = 0.25 + p2_coeff * p2 + p4_coeff * p4 + p6_coeff*p6

    return s

def spherical_decompostion(input_arr, lon_arr, lat_arr):
    """Routine that takes evenly-spaced lat-lon data and transforms to spherical harmonics"""


    coeffs = pysh.expand.SHExpandDH(input_arr, norm=3) #norm=3 returns un-normalised coefficients.
    
#     check_grid = pysh.expand.MakeGridDH(coeffs) #Inverse transforms to check that input_arr and check_grid match closely.
    
    return coeffs

if __name__=="__main__":

    num_points_in_year = 360. #How many points to split the year into
    
    #Create lat-lon grid
    lons = np.arange(0.,360., 1.)* np.pi/ 180.
    lats = np.arange(-90.,90., 1.)* np.pi/ 180.
    lon_arr, lat_arr = np.meshgrid(lons,lats)
    
    #Create storage arrays
    cosz_time_arr = np.zeros((int(num_points_in_year), lats.shape[0], lons.shape[0]))
    area_av_time_arr = np.zeros(int(num_points_in_year))

    #Run calculation separately for each point in the year
    
    for tick in range(int(num_points_in_year)):

        if tick % 10 == 0.:
            print(tick)
        time_since_ae = 2.*np.pi * tick / num_points_in_year
        cosz_time_arr[tick,...] = daily_mean_solar(time_since_ae, lat_arr)
        area_av_time_arr[tick] = area_average(cosz_time_arr[tick,...], lat_arr)
    
    #Perform annual average and spherical decomposition
    annual_average = np.mean(cosz_time_arr, axis=0)
    annual_average_coeffs = spherical_decompostion(annual_average, lon_arr, lat_arr)
    
    plt.plot(annual_average[:,0], label='ann_av')

    #Compare with Frierson sw scheme with various parameters
    for del_sw_value in [0.8, 0.9, 1.0, 1.1]:
        frierson_sw_arr = frierson_sw(lat_arr, del_sol = del_sw_value)
        plt.plot(frierson_sw_arr[:,0], label='frierson del_sol='+str(del_sw_value))
    
    #Compare with 3 coefficient expansion in Legendre polynomials
    mod_frierson_sw_arr = mod_frierson_sw(lat_arr, -annual_average_coeffs[0,2,0], -annual_average_coeffs[0,4,0], -annual_average_coeffs[0,6,0] )
    
    plt.plot(mod_frierson_sw_arr[:,0], label='p2 p4 p6 expansion')
        
    plt.legend()
    plt.show()
    
    
    pdb.set_trace()
