#!/bin/bash -l

#SBATCH --job-name=held_suarez_test_case

#SBATCH --qos=debug
#SBATCH --constraint=knl
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1

#SBATCH --time=00:30:00
#SBATCH --license=project,SCRATCH
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end



echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`



source /etc/profile.d/modules.sh

module load $HOME/.bashrc
module load cray-netcdf
module load git
module load python/3.6-anaconda-5.2
source activate isca_env

python $GFDL_BASE/exp/test_cases/held_suarez/held_suarez_test_case.py
