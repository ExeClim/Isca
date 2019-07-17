# Instructions for running Isca on GW4's Isambard Supercomputer (GNU compiler)

These instructions are intended to get you up-and-running with a simple Held-Suarez test case. They assume you are starting with a default user environment on BlueCrystal 3, so some changes might be needed if you have already modified your environment.

Isamabard is a cray supercomputer comprised of Arm processor nodes. This means that configuration is slightly different to that of the usual Intel-based machines. 

First, you need to clone the Isca repository into your home directory

```{bash}
$ git clone git@github.com:ExeClim/Isca.git
$ cd Isca
```

Before you can run Isca, you'll need to load the CrayPython module:

```{bash}
$ module load cray-python/3.6.5.6
```

Next we'll make a Python environment for Isca (this means it will have all the right versions of the various packages on which it depends):

```{bash}
$ python3 -m venv ~/isca_env 
$ source ~/isca_env/bin/activate
(isca_env) $ cd Isca/src/extra/python
(isca_env) $ pip install -r requirements.txt

Successfully installed MarkupSafe-1.0 f90nml jinja2-2.9.6 numpy-1.13.3 pandas-0.21.0 python-dateutil-2.6.1 pytz-2017.3 sh-1.12.14 six-1.11.0 xarray-0.9.6

(isca_env) $ pip install -e .
...
Successfully installed Isca
```

Finally, we'll need to update the `~/.bashrc` file. Add the following lines:

```{bash}
# directory of the Isca source code
export GFDL_BASE=$HOME/Isca
# "environment" configuration for bc4
export GFDL_ENV=bristol-bc4
# temporary working directory used in running the model
export GFDL_WORK=$HOME/Isca_work
# directory for storing model output
export GFDL_DATA=$HOME/Isca_data

#Cray-specific environment variables 
export EXECUTION_TYPE=APRUN 
#export GFDL_MKMF_TEMPLATE=isambard_gfortran                                                                        #export GFDL_ENV=isambard_gfortran
module unload PrgEnv-cray/6.0.5
```

Then make the `Isca_work` and `Isca_data` directories:

```{bash}
(isca_env) $ mkdir -p $HOME/isca_home/isca_work
(isca_env} $ mkdir -p $HOME/isca_home/isca_data
(isca_env) $ bash
```

Now everything should be set up and we can try a test run. The following should compile and run 12 months of a Held-Suarez test case, at T85 resolution spread over 64 cores. 

```{bash}
(isca_env) $ cd Isca/exp/site-specific/gw4_isambard
(isca_env) $ qsub isca_pbspro.job
```

This should produce a two outputfiles: isca_pbspro.job.e<jobid> and isca_pbspro.job.e<jobid> reporting the standard output and error output of the job.

You can view details about the job using 

```{bash}
(isca_env) $ qstat -f <jobid>
```

All being well, after about 60 minutes the job should complete, and you'll find some output files in `$GFDL_DATA`.

```{bash}
(isca_env) $ cd $GFDL_DATA
```
