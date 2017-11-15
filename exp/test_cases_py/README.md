## Using the `isca` library for running experiments

`isca` is a Python module that provides an easy way to configure, compile and run an Isca experiment.

The primary advantages of using it are:

1. Each experiment is completely isolated - there is no chance another run will overwrite existing data, even if both are run in parallel.
2. Experiments can be tied to a specific git commit, from which the code will be compiled and run.  This ensures repeatability of results when the same script is run again in the future.
2. Namelists and diagnostic tables can be created and manipulated in Python.
3. Experiment scripts can be killed at any point. In the future when you want to resume, simply rerunning the script will continue from the last completed run.

### Installation
1. Activate your preferred Python environment. I use Anaconda and have an env called `jpisca` on all servers.
```
jp492@emps-gv4:~$ module load python/anaconda
jp492@emps-gv4:~$ source activate jpisca
discarding /usr/local/anaconda-2.1.0/bin from PATH
prepending /scratch/jp492/envs/jpisca/bin to PATH
```
2. Set the environment variables to point at the location of your source, where you want temporary working directories to be created, and where you want run data to be stored.
```
(jpisca)jp492@emps-gv4:~$  export GFDL_WORK=/scratch/jp492/gfdl_work
(jpisca)jp492@emps-gv4:~$  export GFDL_BASE=/scratch/jp492/isca
(jpisca)jp492@emps-gv4:~$  export GFDL_DATA=/scratch/jp492/gfdl_data
(jpisca)jp492@emps-gv4:~$  export GFDL_ENV=emps-gv
```
3. Navigate to `src/extra/python` and install the `isca` library in development mode:
```
(jpisca)jp492@emps-gv4:~$ cd $GFDL_BASE/src/extra/python
(jpisca)jp492@emps-gv4:.../python$ pip install -e .
Obtaining file:///scratch/jp492/isca/src/extra/python
...
Installing collected packages: Isca
  Running setup.py develop for Isca
Successfully installed Isca-0.1
(jpisca)jp492@emps-gv4:.../python$
```

The module will now be available in your Python environment.  The advantage of installing in this way is that any changes you make to the module source code in `src/extra/python/isca` will be made available immediately in your experiment scripts.

### Usage

There are two main components to running an experiment:

1. A CodeBase
2. An Experiment

The `CodeBase` object configures which version of the model you want to run - moist, dry, shallow etc, and handles compilation.

The `Experiment` object configures the parameters (Namelist) and diagnostic output (DiagTable) of a specific set of runs of the model.

A basic example experiment script setup can be seen in `exp/test_case_py/held_suarez.py`, which runs the Held Suarez Newtonian cooling experiment at T42 resolution.

For further documentation, dig into the source! It should be documented and fairly obvious how it works.  The entire `Experiment` and `DiagTables` objects are found in `src/extra/python/isca/<experiment, codebase, ...>.py`.
