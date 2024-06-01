#!/bin/bash
#SBATCH --job-name=soc_aquaplanet_co2_test1  # 作业名称
#SBATCH --output=output_%j.out  # 输出文件名，%j表示当前作业ID
#SBATCH --error=error_%j.err    # 错误输出文件名
#SBATCH -N 1 # how many nodes are using 1 node = 32 cores
#SBATCH -c 16 # how many cores are using
#SBATCH --partition=wzhcnormal
#SBATCH --exclusive
module purge
ulimit -s unlimited
ulimit -l unlimited
CaseName=soc_aquaplanet_co2_test1
source /work/home/ac9b0k6rio/miniconda3/bin/activate isca_soc2207_env
source /work/home/ac9b0k6rio/isca_17/Isca/src/extra/env/sugon

python ${GFDL_BASE}/exp/run_isca/${CaseName}.py

