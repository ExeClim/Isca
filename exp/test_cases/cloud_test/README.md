# Using Isca with the Socrates radiation scheme

The Socrates radiation scheme is a highly-flexible scheme written and maintained by the UK Met Office. It has many significant advantages over RRTM, notably the ease with the properties of the radiation scheme can be changed using input-files, e.g. different numbers of spectral bands, different atmospheric compositions, etc. 

To read more about Socrates itself, see here:
[Socrates project](https://code.metoffice.gov.uk/trac/socrates)


## Getting started

### 1. Download Socrates source code

All the code needed to interface Isca and the Socrates radiation scheme is provided in the folder 

`$GFDL_BASE/src/atmos_param/socrates/interface`,

and this is packed with Isca itself. However, the **Socrates source code is not packaged with Isca**, as the Met Office maintains its own version control system for Socrates that they wish to keep as the difinitive version. The latest version of Socrates is freely available to download from their website:

`https://code.metoffice.gov.uk/trac/socrates`.

Access to this website requires a user account on the Met Office's Science Repository Service. Instructions for how to get an account are given at 

`https://code.metoffice.gov.uk/trac/home`

To download a packaged release of Socrates, click on the `tar.xz` file of the latest version. Once the folder has downloaded, unzip it using

`tar -xf INSERT_FILE_NAME_HERE.xz`

Once you have unzipped the file, the resulting folder should contain folders like `data`, `docs` and `src`. 

### 2. Move the Socrates source code to the correct folder within the Isca directory

Navigate to the following Isca directory:

`$GFDL_BASE/src/atmos_param/socrates/src`

make a folder called `trunk` and then put the contents of your downloaded Socrates code into the `trunk` folder.

You should then have the following directory structure:

`$GFDL_BASE/src/atmos_param/socrates/src/trunk/src/radiance_core/`

### 3. Edit the number of angles in the phase function

From this point, the Isca test case for Socrates should run without issue. However, it is strongly advised to first make the following changes to the Socrates source code.

Open the file

`$GFDL_BASE/src/atmos_param/socrates/src/trunk/src/modules_core/dimensions_spec_ucf.F90`

and make the following changes:

```
npd_k_term=14
npd_scale_variable = 4
npd_continuum = 2
npd_drop_type = 5
npd_cloud_parameter = 30
npd_phase_term = 1
```

The final one of these is the most important, as a large value for this term significantly increases Socrates' memory usage, and will make Isca much slower. 

### 4. Run the Socrates test case.

Navigate to 

`$GFDL_BASE/exp/test_cases/socrates_test/` 

and run the test-case `socrates_aquaplanet.py`. This will compile and run Isca with a default Earth-like aquaplanet.


### 5. If Socrates seems slow, here are some hints

* Isca is set up to pass socrates a certain number of vertical profiles for each time the Socrates is called. This number is set as `chunk_size` in the `socrates_rad_nml`. A value of 16 was found to be optimal on large linux-server-type machines in Exeter, but it is worth playing with this number to find the optimal number on your machine.

* Socrates reads external input files that tell it the number of spectral bands to use, with one file setting the short-wave options, and another file setting the long-wave options. Some spectral files have lots of bands, which will make the model run slowly. The default files used in the Met Office's Unified Model-GA7, and also in Isca, can be found here:
	* `$GFDL_BASE/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7` for the long-wave
	* `$GFDL_BASE/src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7` for the short-wave
* Other options are available within this folder, and a useful set of other spectral files are provided via this [webpage](https://simplex.giss.nasa.gov/gcm/ROCKE-3D/).


	