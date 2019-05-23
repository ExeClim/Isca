#!/bin/bash -l

#SBATCH --job-name=held_suarez_test_case
#SBATCH --partition=veryshort
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#number of tasks ~ processes per node
#SBATCH --ntasks-per-node=16
#number of cpus (cores) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_%j.o

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`


module purge
source $HOME/.bashrc
source $GFDL_BASE/src/extra/env/bristol-bc4
source activate isca_env



$HOME/.conda/envs/isca_env/bin/python $GFDL_BASE/exp/test_cases/held_suarez/held_suarez_test_case.py
