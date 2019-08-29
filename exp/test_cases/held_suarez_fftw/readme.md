# Held-Suarez using FFTW

This configuration is idential to the Held-Suarez configuration found in `../held_suarez`, however it uses the FFTW library to perform the FFT instead of the default Temperton's FFT. 

## Running this configuration
To run this configuration, you need to ensure that you have FFTW installed on your system, and that it can be linked to the Isca code. This usually involves loading the relevant module file. For BCP4 at Bristol, this is done by adding the following command to the environment file:

```
module load libs/fftw/3.3.6
```

Once this has been done, you can simply run from within this folder: 

```
python held_suarez_test_case_fftw.py
```

## Other configurations
In order to use FFTW with other model configurations there are two steps you must take:

### 1. Add the correct preprocessor directive 

You must add the `-DFFTW3` preprocessor directive to the list of compiler flags that will be used to compile Isca. This can be done by adding the following line to your python script:

```
cb.compile_flags.append('-DFFTW3')
``` 
This must be done **before** calling `cb.compile()`. An example of this is given in the `held_suarez_test_case_fftw.py` file inside this folder. 

### 2. Add the fftw.F90 module

To use FFTW, a wrapper module has been created (`Isca/src/shared/fft/fftw.F90`). This file must be added to the list of files that are to be compiled using your configuration, as done in `src/extra/model/isca_fftw/path_names`.



Having followed both steps, the FFTW3 library should be configured, however feel free to send me an email at qv18258@bristol.ac.uk if you get stuck. 

