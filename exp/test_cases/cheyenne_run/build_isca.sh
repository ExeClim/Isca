#!/bin/bash -f
#PBS -A UHAR0008
#PBS -N fixqflux
#PBS -q regular
#PBS -l select=1:ncpus=16:mpiprocs=16:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -S /bin/csh -V
#PBS -M kangwanying1992@gmail.com
#PBS -m abe
#############
# environment:
source /glade/u/home/${USER}/.bashrc
module purge
module load ncarenv/1.2
module load intel/17.0.1
module load ncarcompilers
module load mpt/2.15f
module load netcdf/4.4.1.1
module load git
source /glade/p/work/${USER}/${USER}_python_package/bin/activate
export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

# passing info
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

export casename=${1:-"test"}
export ncores=${2:-"16"}
export todo=${3:-"compile"}
echo casename=$casename, ncores=$ncores, todo=$todo

# running model
python build_isca.py $casename $ncores $todo

