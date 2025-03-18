#!/bin/bash
# This script is an example used to build and run the held suarez test case
# on the NSF NCAR Derecho Supercomputer. Following the previous steps in the
# README is required for this script to work correctly. You must substitute
# the [ACCOUNT CODE] with your own account code before submitting the job.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=128:mpiprocs=128
#PBS -q main
#PBS -A [ACCOUNT CODE]
#PBS -N Isca-test
#PBS -k eod
#PBS -j oe
#PBS -o Isca-test.out

module --force purge
module load ncarenv/23.09 intel/2023.2.1 craype/2.7.31 cray-mpich/8.1.27 ncarcompilers/1.0.0 hdf5/1.12.2 netcdf/4.9.2 conda/latest
conda activate isca_env

cd $GFDL_BASE/exp/test_cases/held_suarez
python held_suarez_test_case.py
