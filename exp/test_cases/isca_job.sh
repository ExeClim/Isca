#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -d . # set working directory to .
#PBS -q ptq # submit to the parallel test queue
#PBS -l nodes=1:ppn=8 # nodes=number of nodes required. ppn=number of processors per node
#PBS -l walltime=0:30:00 # Maximum wall time for the job.
#PBS -A Research_Project-PROJECTNUMBER # research project to submit under.
#PBS -m e -M USER@exeter.ac.uk # email me at job completion

# This is an example of how to submit an experiment to the queue on the Isca supercomputer
# in Exeter.  USER and PROJECTNUMBER in the metadata above need to be
# updated to appropriate values for your experiment.
# For more information about running on the Exeter cluster, see the wiki
# here: https://wiki.exeter.ac.uk/display/ISCA/Isca+Wiki
module load Anaconda3
source activate gfdl
python held_suarez/parameter_sweep.py
