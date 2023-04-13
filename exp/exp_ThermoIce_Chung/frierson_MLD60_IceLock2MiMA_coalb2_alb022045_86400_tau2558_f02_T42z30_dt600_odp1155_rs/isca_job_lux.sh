#!/bin/bash
#SBATCH --job-name=frierson_MLD60_IceLock2MiMA_coalb2_alb022045_86400_tau2558_f02_T42z30_dt600_odp1155_rs # Job name
#SBATCH --partition=cpuq             # queue for job submission
#SBATCH --account=cpuq               #
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=pochung@ucsc.edu   # Where to send mail
#SBATCH --ntasks=160                  # Number of MPI ranks
#SBATCH --nodes=4                    # Number of nodes
#SBATCH --ntasks-per-node=40         # How many tasks on each node
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=mpi_test_%j.log     # Standard output and error log

pwd; hostname; date

echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."

source $HOME/.bashrc
source $GFDL_BASE/src/extra/env/ucsc-lux-gfortran
source $HOME/Isca/venvs/Isca/bin/activate

python $GFDL_BASE/exp/exp_ThermoIce_Chung/frierson_MLD60_IceLock2MiMA_coalb2_alb022045_86400_tau2558_f02_T42z30_dt600_odp1155_rs/frierson_MLD60_IceLock2MiMA_coalb2_alb022045_86400_tau2558_f02_T42z30_dt600_odp1155_rs.py

date
