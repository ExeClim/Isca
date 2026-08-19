Experiment configuration scripts
==============

Summary
-------
Once you have managed to run a test case, you may want to run your own experiment. The model configuration, diagnostics, and simulation parameters
of your experiment can be configured with a python-based run script.
Find the test case in ``Isca/exp/test_cases`` that most resembles your experiment, copy and modify it. The important lines are described below.

Setup
-------
::

    from isca import DryCodeBase
imports the Isca python module that should have been installed. 
``DryCodeBase`` refers to the part of the code needed to run the dry model. The codebase object is defined in ``Isca/src/extra/python/isca/codebase.py``. The different default codebases are designed to avoid the model compiling un-needed parts of the code, which speeds up compilation. ``DryCodeBase`` can be replaced by ``IscaCodeBase``, ``GreyCodeBase``, or ``SocratesCodeBase`` (provided you setup Socrates, see Socrates documentation).

::

    NCORES = 16
specifies how many cores should be used when running the model.

::

    RESOLUTION = T42, 25
specifies the horizontal (in this case, spectral T42) and vertical (25 pressure levels) 
resolution. Runs at T42 by default. If the horizontal resolution is not T42, make sure the 
input files (such as land masks) are changed and add ``exp.set_resolution(*RESOLUTION)`` after exp is defined.
Other common spectral resolutions are T21 and T85.

::

    cb = DryCodeBase.from_directory(GFDL_BASE)
tells Isca to find the model code at the location ``GFDL_BASE``, which is where the source code may be changed. You can replace ``from_directory(GFDL_BASE)`` with ``from_repo(repo='https://github.com/isca/isca', commit='isca1.1')`` to point to a specific git repository and commit id. This should ensure future, independent, reproducibility of results. But, the compilation depends on computer specific settings.  The ``$GFDL_ENV`` environment variable is used to determine which ``$GFDL_BASE/src/extra/env`` file is used to load the correct compilers.  The env file is always loaded from ``$GFDL_BASE`` and not the checked out git repo.

::

    cb.compile()
compiles the source code.

::

    exp = Experiment('EXPERIMENT_NAME', codebase=cb)
creates an experiment object which will handle the configuration of model parameters and output diagnostics. The output files will be found at ``$GFDL_DATA/EXPERIMENT_NAME``.


Diagnostics
-------

::

    diag = DiagTable()
creates a DiagTable object which we can configure to tell Isca which variables to output.

::

    diag.add_file('atmos_monthly', 30, 'days', time_units='days')
creates an ``atmos_monthly.nc`` file every 30 days, which is an average of the output over the previous 30 days. The output files can be found at ``$GFDL_DATA/EXPERIMENT_NAME/run####/*``.

::

    diag.add_field(MODULE_NAME, VARIABLE_NAME, time_avg=True)
determines which fields will be written in ``atmos_monthly.nc``. Find the available VARIABLE_NAMEs by going to the MODULE_NAME documentation or by finding the relevant source code (``cd Isca/src/ & find . -name "MODULE_NAME*"``).


Namelist
-------

::

    namelist = Namelist({...})
defines a namelist object, which lets us configure the science options. 
It is only necessary to set values that are different from the default parameters, which are defined 
in the relevant module documentation (for example, ``atmosphere_nml`` parameters can be found in the ``atmosphere`` 
module documentation or at the beginning of the ``atmosphere.F90`` source file).

Running the experiment
-------

::

    exp.run(...)
will make the model run for the amount of time specified in ``main_nml`` (usually 30 days). 

The ``use_restart`` option can be set to ``False`` to start from scratch (isothermal atmosphere) or can point to a restart file (``use_restart = $GFDL_DATA/exp_name/run####/restarts/*``) to initialize the run from the output of a previous run. If unspecified, it will start from where the previous run left off or from an isothermal atmosphere in the absence of a previous run.

Output
-------

Output from the experiment can be found at ``$GFDL_DATA/EXP_NAME``. The atmospheric output is provided on 
sigma levels where sigma is the pressure normalized by the surface pressure. For a planet with no topography, sigma and pressure levels are quite similar. 
If there is topography present (such as in the ``realistic_continents`` test case), you need to interpolate the 
data onto pressure levels before analyzing it. Top of atmosphere and surface values are not affected, but in-atmosphere values are.

The details and code for interpolation to pressure levels can be found at https://github.com/ExeClim/Isca/tree/master/postprocessing/plevel_interpolation

In the python code, there is a convenient function which can be used to call the interpolation code: https://github.com/ExeClim/Isca/blob/master/src/extra/python/isca/util.py  (line 134).

For example::

    from isca.util import interpolate_output
    for run in ["EXPERIMENT_NAME"]: 
        print(run)    
        for i in range(121, 241):
            try:
                infile = '/data_directory/' + run + '/run%04d/atmos_monthly.nc' % i   
                outfile = '/data_directory/' + run + '/run%04d/plev_monthly.nc' % i
                interpolate_output(infile, outfile, p_levs='EVEN', var_names=['slp', 'height'])
            except:
                print(i)

Authors
-------

This documentation was written by Matthew Henry (heavily inspired from document written by Neil Lewis), peer reviewed by Will Seviour, and quality controlled by Brett McKim.
