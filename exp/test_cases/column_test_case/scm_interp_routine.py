import xarray as xr 
import numpy as np 
from isca import GFDL_BASE
import os 



def vinterp(data, vcoord, vlevels):
    
    """ vertical linear interpolation, credit ExeClim/ShareCode"""

    assert (vcoord.ndim == data.ndim or vcoord.ndim == 1 and data.ndim == 4 or
            vcoord.ndim == 4 and data.ndim == 1)
    if vcoord.ndim == 1 and data.ndim > 1:
        # This handles the case where vcoord is 1D and data is N-D
        v_dim = int(np.where(np.array(data.shape) == vcoord.shape[0])[0])
        # numpy.broadcast_to only works for the last axis of an array,
        # swap our shape around so that vertical dimension is last, broadcast
        # vcoord to it, then swap the axes back so vcoord.shape == data.shape
        data_shape = list(data.shape)
        data_shape[-1], data_shape[v_dim] = data_shape[v_dim], data_shape[-1]

        vcoord = np.broadcast_to(vcoord, data_shape)
        vcoord = np.swapaxes(vcoord, -1, v_dim)

    vcoord_shape = list(vcoord.shape)
    vcoord_shape.pop(1)
    vcoord_shape = tuple(vcoord_shape)

    valid = np.min([np.prod(vcoord_shape) -
                    np.sum(np.isnan(vcoord[:, 0, ...])),
                    np.prod(vcoord_shape) -
                    np.sum(np.isnan(vcoord[:, -1, ...]))])

    if np.sum(vcoord[:, 0, ...] > vcoord[:, -1, ...]) / valid > 0.80:
        # Vcoord data is decreasing on interpolation axis, (at least 80% is)
        idx_gt = 1
        idx_lt = 0
    else:
        # Data is increasing on interpolation axis
        idx_gt = 0
        idx_lt = 1

    if data.ndim >= vcoord.ndim:
        # Handle case where data has the same dimensions or data has more
        # dimensions compared to vcoord (e.g. vcoord is 4D, data is 4D,
        # or vcoord is 1D, data is 4D)
        out_shape = list(data.shape)
    else:
        # Handle case where data has fewer dimensions than vcoord
        # (e.g. data is 1-D vcoord is N-D)
        out_shape = list(vcoord.shape)
    out_shape[1] = vlevels.shape[0]

    out_shape = tuple(out_shape)
    out_data = np.zeros(out_shape) + np.nan

    for lev_idx, lev in enumerate(vlevels):
        if idx_gt == 0:
            # Case where vcoord data is increasing, find index where
            # vcoord below [:-1] is equal or less than desired lev, and
            # vcoord above [1:] is greater than lev, this means <data> for lev
            # is between these points, use weight to determine exactly where
            idx = np.where(np.logical_and(vcoord[:, :-1, ...] <= lev,
                                          vcoord[:, 1:, ...] > lev))

        else:
            # This does the same, but where vcoord is decreasing with index,
            # so find where vcoord below [:-1] is greater, and vcoord above
            # [1:] is less or equal
            idx = np.where(np.logical_and(vcoord[:, :-1, ...] > lev,
                                          vcoord[:, 1:, ...] <= lev))
        # Reduce diminsions of `idx`
        idx = np.squeeze(idx)
        # Create copies of this index, so they can be modified for
        # weighting functions and output array
        idx_abve = idx.copy()
        idx_belw = idx.copy()
        out_idx = idx.copy()

        # The interpolation axis index (1) for output
        # is the level index (lev_idx)
        out_idx[1, :] = lev_idx

        # Weighting function 'above' is index +1 for decreasing,
        # or index +0 for decr.
        idx_abve[1, :] += idx_gt
        # Weighting function 'below' is index +0 for decreasing,
        # or index +1 for decr.
        idx_belw[1, :] += idx_lt

        # Change indicies back into tuples so
        # numpy.array.__getitem__ understands them
        idx_abve = tuple(idx_abve)
        idx_belw = tuple(idx_belw)
        out_idx = tuple(out_idx)

        # Weighting function for distance above lev
        wgt1 = ((lev - vcoord[idx_belw]) /
                (vcoord[idx_abve] - vcoord[idx_belw]))

        # Weighting function for distance below lev
        wgt0 = 1.0 - wgt1

        if data.ndim >= vcoord.ndim:
            # Handle case where data has same or more dimensions than vcoord
            out_data[out_idx] = (wgt0 * data[idx_belw] + wgt1 * data[idx_abve])
        else:
            # Handle case where data has fewer dimensions than vcoord
            out_data[out_idx] = (wgt0 * data[idx_belw[1]] +
                                 wgt1 * data[idx_abve[1]])

    return np.squeeze(out_data)

def global_average_lat_lon(ds_in, var_name, radius=6371.e3):

    try:
        ds_in['area_array']
    except KeyError:
        cell_area(ds_in, radius)

    
    weighted_data = ds_in[var_name]*ds_in['area_array']

    area_average = weighted_data.mean(('lat', 'lon')) / ds_in['area_array'].mean(('lat','lon'))

    var_in_dims = ds_in[var_name].dims

    var_out_dims = tuple(x for x in var_in_dims if x!='lat' and x!='lon')

    ds_in[var_name+'_area_av'] = (var_out_dims, area_average)

def cell_area(dataset_in, radius = 6371.e3):

    lonb = dataset_in['lonb']
    latb = dataset_in['latb']

    lonb_1 = lonb[1::].values
    lonb_2 = lonb[0:-1].values

    delta_lon = lonb_1 - lonb_2

    latb_1 = latb[1::].values
    latb_2 = latb[0:-1].values

    delta_lat = latb_1 - latb_2

    dataset_in['delta_lon'] = (('lon'), delta_lon)
    dataset_in['delta_lat'] = (('lat'), delta_lat)

    dataset_in['latb_1'] = (('lat'), latb_1)
    dataset_in['latb_2'] = (('lat'), latb_2)

    xsize = radius*np.absolute(np.deg2rad(dataset_in['delta_lon']))*(np.sin(np.deg2rad(dataset_in['latb_1']))-np.sin(np.deg2rad(dataset_in['latb_2'])))
    ysize = radius

    area_array = xsize*ysize

    dataset_in['area_array'] = (('lat','lon'), area_array.transpose('lat','lon'))

def pkbk(coord_option, nlevels, surf_res=.1, exponent=2.5, scale_heights=4.):

    if coord_option == 'even_sigma':
        pk = np.zeros(nlevels+1)
        bk = np.zeros(nlevels+1)
        for lvl in range(0, nlevels):
            bk[lvl] = float(lvl) / float(nlevels)
        bk[-1] = 1.0
    elif coord_option == 'uneven_sigma':
        pk = np.zeros(nlevels+1)
        bk = np.zeros(nlevels+1)
        for lvl in range(0, nlevels):
            zeta = (1. - (float(lvl)/float(nlevels)))
            z = surf_res*zeta + (1. - surf_res)*(zeta**exponent)
            bk[lvl] = np.exp(-z * scale_heights)
        bk[-1] = 1.0
        bk[0] = 0.0
    else: 
        print('pkbk: '+coord_option+' is NOT a coordinate option supported by this script')
        
    return pk, bk 

def calc_pfull(pk, bk, psurf, diff_option):

    nhalflevels = len(pk)
    phalf = np.zeros(nhalflevels)
    ln_top_level_factor = -1.0

    for lvl in range(0, nhalflevels):
        phalf[lvl] = pk[lvl] + bk[lvl] * psurf 
    
    lnphalf = np.zeros(nhalflevels)
    lnpfull = np.zeros(nhalflevels-1)
    if diff_option == 'simmons_and_burridge':
        for lvl in range(1, nhalflevels):
            lnphalf[lvl] = np.log(phalf[lvl])
        
        for lvl in range(1, nhalflevels-1):
            alpha = 1 - phalf[lvl]*(lnphalf[lvl+1] - lnphalf[lvl])/(phalf[lvl+1]-phalf[lvl])
            lnpfull[lvl] = lnphalf[lvl+1] - alpha
        lnpfull[0] = lnphalf[1] + ln_top_level_factor 
        lnphalf[0] = 0.0
        pfull = np.exp(lnpfull)
    else:
        print('calc_pfull: '+diff_option+' is NOT a vertical differencing option supported by this script')


    return pfull 

def scm_interp(filename, varname='ozone_1990', vcoord_option='even_sigma', nlevels=None, 
                 vert_difference_option = 'simmons_and_burridge', psurf=1.e3, pk_input=None, bk_input=None):

    # psurf in hPa 


    ds = xr.open_dataset(filename, decode_times=False)

    if vcoord_option != 'input':
        pk, bk = pkbk(vcoord_option, nlevels)
    else:
        pk = pk_input
        bk = bk_input

    pfull = calc_pfull(pk, bk, psurf, vert_difference_option)
    

    data = vinterp(ds[varname].values, ds.pfull.values, pfull)


    ds.coords['pfull_new'] = pfull 
    ds[varname+'_interp'] = (('time', 'pfull_new', 'lat'), data)


    return ds 

if __name__ == "__main__":

    scm_interp(filename=os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'), 
               varname='ozone_1990', 
               nlevels=31)

    