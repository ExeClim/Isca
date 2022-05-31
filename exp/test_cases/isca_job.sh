# This is an example of how to submit an experiment to the queue on the Isca supercomputer
# in Exeter.  USER and PROJECTNUMBER in the metadata above need to be
# updated to appropriate values for your experiment.
# For more information about running on the Exeter cluster, see the wiki
# here: https://universityofexeteruk.sharepoint.com/sites/ExeterARC

#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=03:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-PROJECTNUMBER # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=USER@exeter.ac.uk # email address

module load Anaconda3
source activate isca_env
python held_suarez/held_suarez_test_case.py
