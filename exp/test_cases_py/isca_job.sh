
#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -d . # set working directory to .
#PBS -q ptq # submit to the parallel test queue
#PBS -l nodes=1:ppn=8 # nodes=number of nodes required. ppn=number of processors per node
#PBS -l walltime=0:30:00 # Maximum wall time for the job.
#PBS -A Research_Project-161613 # research project to submit under. 
#PBS -m e -M jp492@exeter.ac.uk # email me at job completion
python example_basic.py
