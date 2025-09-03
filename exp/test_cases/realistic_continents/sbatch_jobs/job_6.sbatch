#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-00:00
#SBATCH --job-name=RT42_sst_E_6
#SBATCH --output=/scratch/philbou/outerr/%x-%j.out
#SBATCH --error=/scratch/philbou/outerr/%x-%j.err
#SBATCH --account=def-rfajber

export GFDL_BASE=/home/philbou/Isca
export GFDL_ENV=narval.ifort
export GFDL_WORK=/scratch/philbou/isca_work
export GFDL_DATA=/scratch/philbou/isca_data
source /home/philbou/.bashrc
conda activate isca_env
cd $GFDL_BASE/exp/test_cases/realistic_continents

python RT42_sst_m1.py
