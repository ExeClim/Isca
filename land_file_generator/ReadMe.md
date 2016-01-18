# ReadMe #

The planet's surface in MiMA is a mixed-layer ocean with no dynamics, meaning that it has a simple heat capacity, and exchanges sensible and latent heat fluxes with the atmosphere. In MiMA it is possible to add areas of `land` to the surface, where `land` is differentiated from the surrounding mixed-layer ocean by differences in surface properties. The surface properties that can be modified in areas defined as `land` are:

* Mixed-layer heat capacity
* The surface albedo
* The surface roughness length
* Latent heat fluxes

It is to be noted that MiMA makes a distinction between `land` and `topography`, with `land` being areas where the surface properties are different, and `topography` being areas where the surface height is different. In principle, areas of `land` and `topography` do not have to be the same.

In what follows `land` will be added first, with the extension to adding `topography` discussed later.

### Adding land ###


The first step is to generate a file that will tell the model where the `land` is. This is done using a NetCDF input file, which can be generated using `land_file_generator.py`. See section on **_Using land_file_generator.py_** for instructions on how to do this. The result will be a file named `land.nc`, whose only varibales are `land_mask`, `latitude` and `longitude`. The `land_mask` array contains `0's` and `1's` defining `not land` and `land`, respectively. 

The second step is to tell the model you want to add land. This is done **_first_** by placing the relevant `land.nc` file in the `input` folder of `./exp/YOUR_EXPERIMENTS_NAME_HERE/input/`. **_Second_**, the variable `land_nc_file = land.nc` must be set in `./exp/YOUR_EXPERIMENTS_NAME_HERE/runscript`. Then **_thirdly_** `land_option = 'input'` in `idealized_moist_phys_nml` which is part of the `phys.nml` file. This will read in the `land_mask` array, and turn it into a logical array called `land`. 

The final step is then to decide which of the surface properties you wish to modify. These are done as follows:

* Mixed-layer heat capacity:
    * Set `land_option = 'input'` in `mixed_layer_nml`
    * Change `land_heat_capacity_prefactor` to a value other than `1.0`, which is its default. Setting a value of `0.8`, for example, will make the heat capacity over land 80% of its value elsewhere. 
* The surface albedo
    * Make sure `land_option = 'input'` in `mixed_layer_nml`
    * Change `land_albedo_prefactor` to a value other than `1.0`, which is its default. Setting a value of `0.8`, for example, will make the surface albedo over land 80% of its value elsewhere. 
* The surface roughness length
    * Make sure `land_option = 'input'` in `idealized_moist_phys_nml`, which it should be anyway if you're using `land`. 
    * Change `land_roughness_prefactor` to a value other than `1.0`, which is its default. Note that all 3 of the roughness lengths (i.e. `roughness_mom`, `roughness_heat` and `roughness_moist`) will be multiplied by this prefactor.
* Latent heat fluxes
    * Make sure `land_option = 'input'` in `idealized_moist_phys_nml`, which it should be anyway if you're using `land`. 
    * Change `land_humidity_prefactor` to a value other than `1.0`, which is its default. This prefactor is put into the equation for the evaporative flux E \propto (q_a - land_humidity_prefactor * q_s) where `q_a` is the specific humidity at the lowest model level, and `q_s` is the equivalent specific humidity of the surface, as in equation (11) of Frierson et al. 2006. This prefacor, when it is < 1.0, has the effect of making the surface drier, meaning latent heat fluxes will be lower than over areas of ocean. 
    * **_Note_** that there may well be some problems with this approach. Including the possibility of latent heat fluxes into the atmosphere being **_negative_**, i.e. the surface is withdrawing heat from the atmosphere and putting it into the surface. The realism of these negative fluxes needs to be discussed. In the future it would be best to go to a `bucket` surface model. 

### Adding topography ###

Adding `topography` is done by essentially the same process as adding `land`.

The first step is to generate a file that will tell the model where the `topography` is, and what shape / height the `topography` takes. This is done using a NetCDF input file, which can be generated using `topo_land_file_generator.py`. See section on **_Using topo_land_file_generator.py_** for instructions on how to do this. The result will be a file named `land.nc`, whose varibales are the same as before: `land_mask`, `latitude` and `longitude`, but now also include `zsurf`, which is the height of the topography (in metres).  

The second step is to tell the model you want to add topography. This is done by setting `topography_option = 'input'` and `topog_file_name = 'land.nc'` in `spectral_init_cond_nml` which is part of the `phys.nml` file. This will read in the `zsurf` array from `land.nc' and use this as the topography. **_It is to be noted that the model currently outputs data on sigma levels, and so using topography will require the standard interpolation onto pressure levels._**

### Using land_file_generator.py ###

In `land_file_generator.py` the program first requires the definition of the resolution of the grid to be used. This is set using the variable `t_res`, and is set to `42` in the current example. The program then looks in `./gfdl_grid_files/` for a file named `t42.nc`, from which it will read in the latitude and longitude grid. When you want to use a resolution other than `t42` a new resolution file must be made, and it must then be placed in `./gfdl_grid_files/`. I have written a program to do this called `grid_file_generator.py`, instructions for which are below.

Having defined the grid, `land_file_generator.py` then creates `land_array` of size `(nlat, nlon)`. To define where the land is, a list of array indices called `idx` is created based on some criteria for land, and then `land_array[idx] = 1.0`. In the current example

```
idx = (20. < lat_array) & (lat_array < 60) & (20. < lon_array) & (lon_array < 60)
```

which makes `idx` a list of indices where `20. < latitude < 60.` and `20. < longitude < 60`. Latitude goes from `-90 -> 90`. and longitude goes from `0 -> 360`. In principle any shape could be defined this way. 

Having defined `land_array` the results are outputted to `land.nc` and this file can then be used as an input for the model.

### Using grid_file_generator.py ###

The `grid_file_generator.py` program takes an input file with a name of the form `t42_atmos_daily.nc`, then reads in its latitude and longitude values, and outputs a file with name `t42.nc`. Such a file can then be used as input to the `*_file_generator.py` scripts.


### Using topo_land_file_generator.py ###

Defining `land` in `topo_land_file_generator.py` is exactly the same as in `land_file_generator.py`, and in principle can be over different areas than the `topography`. In the example below `topography` and `land` are specified together. 

To generate a gaussian mountain the following code is used:

```
central_lat = 40. 
central_lon = 40.
radius_degrees = 20.
std_dev = 10.
height = 3500.


rsqd_array = np.sqrt((lon_array - central_lon)**2.+(lat_array - central_lat)**2.)

idx = (rsqd_array < radius_degrees) 

topo_array[idx] = height* np.exp(-(rsqd_array[idx]**2.)/(2.*std_dev**2.))
land_array[idx] = 1.0
```

This firsts generates `rsqd_array`, which measures the distance from a central point. The `idx` list is then generated using `rsqd_array`, and `topo_array` and `land_array` are defined using `idx`. 

To specify `land` differently, a different `idx` could be generated and used in `land_array`. 


