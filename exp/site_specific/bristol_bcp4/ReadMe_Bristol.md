
# Instructions for running Isca on BlueCrystal 4 at Bristol

These instructions are intended to get you up-and-running with a simple Held-Suarez test case. They assume you are starting with a default user environment on BlueCrystal 4, so some changes might be needed if you have already modified your environment. Please don't hesistate to get in touch (w.seviour@bristol.ac.uk) if you have any questions.

First, to enable `git`, run:

```{bash}
$ module load git
```

Now you can clone the Isca repository:

```{bash}
$ git clone git@github.com:ExeClim/Isca.git
$ cd Isca
```

Before you can run Isca, you'll need to load the Anaconda python distribution:

```{bash}
$ module load languages/anaconda3
```

Next we'll make a conda environment for Isca (this means it will have all the right versions of the various packages on which it depends):

```{bash}
$ conda create -n isca_env python ipython
$ source activate isca_env
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
export GFDL_WORK=/mnt/storage/home/$USER/scratch/Isca_work
# directory for storing model output
export GFDL_DATA=/mnt/storage/home/$USER/scratch/Isca_data
```

Then make the `Isca_work` and `Isca_data` directories:

```{bash}
(isca_env) $ mkdir -p /mnt/storage/home/$USER/scratch/Isca_work
(isca_env} $ mkdir -p /mnt/storage/home/$USER/scratch/Isca_data
(isca_env) $ bash
```
Now everything should be set up and we can try a test run. The following should compile and run 12 months of a Held-Suarez test case, at T42 resolution spread over 16 cores. 

```{bash}
(isca_env) $ cd Isca/exp/site-specific
(isca_env) $ sbatch isca_slurm_job.sh
```

This should produce a slurm file showing the progress as it compiles and runs. To track the progress:

```{bash}
(isca_env) $ tail -f slurm_*.o
```

All being well, after about 20 minutes the job should complete, and you'll find some output files in `/mnt/storage/home/$USER/scratch/Isca_data`.