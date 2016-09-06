## Using the `gfdl` library for running experiments

`gfdl` is a Python module that provides an easy way to configure, compile and run an experiment based on the `GFDLmoistModel` source code.

The primary advantages of using it are:

1. Each experiment is completely isolated - there is no chance another run will overwrite existing data, even if both are run in parallel.
2. Experiments can be tied to a specific git commit, from which the code will be compiled and run.  This ensures repeatability of results when the same script is run again in the future.
2. Namelists and diagnostic tables can be created and manipulated in Python.
3. Experiment scripts can be killed at any point. In the future when you want to resume, simply rerunning the script will continue from the last completed month.
4. If you use the UNIX screen command, the script will update the screen terminal name to show the current status i.e. the current day that is being run.  The screen continues to update when switching to another tab (see example screenshot below)

![](http://g.recordit.co/ax5h0Hw9IE.gif)

### Installation
1. Activate your preferred Python environment. I use Anaconda and have an env called `jpgfdl` on all servers.
```
jp492@emps-gv4:~$ module load python/anaconda
jp492@emps-gv4:~$ source activate jpgfdl
discarding /usr/local/anaconda-2.1.0/bin from PATH
prepending /scratch/jp492/envs/jpgfdl/bin to PATH
```
2. Set the environment variables to point at the location of your source, where you want temporary working directories to be created, and where you want run data to be stored.
```
(jpgfdl)jp492@emps-gv4:~$  export GFDL_WORK=/scratch/jp492/gfdl_work
(jpgfdl)jp492@emps-gv4:~$  export GFDL_BASE=/scratch/jp492/GFDLmoistModel
(jpgfdl)jp492@emps-gv4:~$  export GFDL_DATA=/scratch/jp492/gfdl_data
```
3. Navigate to `src/extra/python` and install the `gfdl` library in development mode:
```
(jpgfdl)jp492@emps-gv4:~$ cd $GFDL_BASE/src/extra/python
(jpgfdl)jp492@emps-gv4:.../python$ pip install -e .
Obtaining file:///scratch/jp492/GFDLmoistModel/src/extra/python
...
Installing collected packages: GFDL
  Running setup.py develop for GFDL
Successfully installed GFDL-0.1
(jpgfdl)jp492@emps-gv4:.../python$
```

The module will now be available in your Python environment.  The advantage of installing in this way is that any changes you make to the module source code in `src/extra/python/gfdl` will be made available immediately in your experiment scripts.

### Usage

The most basic experiment would be:
```
from gfdl.experiment import Experiment, DiagTable

exp = Experiment('playground', overwrite_data=True)

diag = DiagTable()

# create one or more output files
diag.add_file('6hourly', 6*60*60, 'seconds', time_units='days')

# add diag fields to the output files
diag.add_field('dynamics', 'ucomp')
diag.add_field('dynamics', 'vcomp')
diag.add_field('dynamics', 'temp')
diag.add_field('dynamics', 'vor')
diag.add_field('dynamics', 'div')

diag.add_field('two_stream', 'olr')
diag.add_field('two_stream', 'flux_sw')
diag.add_field('two_stream', 'flux_lw')

exp.use_diag_table(diag)

# compile the source code to $work_dir/exec
exp.disable_rrtm()	# when using two-stream gray rad we don't need rrtm
exp.compile()

exp.clear_rundir()

# set some values in the namelist
# overwrite the whole main_nml section
exp.namelist['main_nml'] = {
    'dt_atmos': 900,
    'seconds': 86400.0*30,
    'calendar': 'no_calendar'
}

# update specific values in some namelist sections
exp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True
exp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
exp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True
exp.namelist['spectral_dynamics_nml']['num_levels'] = 25

exp.runmonth(1, use_restart=False)
for i in range(2, 13):
  exp.runmonth(2)  # use the restart i-1 by default
```

For examples of creating experiments with custom `diag_tables`, changing namelists and looping over several namelist configurations see `exp/python_gfdl/example....py`.

For further documentation, dig into the source! It should be documented and fairly obvious how it works.  The entire `Experiment` and `DiagTables` objects are found in `src/extra/python/gfdl/experiment.py`.
