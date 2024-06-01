#!/bin/bash
#SBATCH --job-name=socrates_aquaplanet   # 作业名称
#SBATCH --output=output_%j.out  # 输出文件名，%j表示当前作业ID
#SBATCH --error=error_%j.err    # 错误输出文件名
#SBATCH -N 1 # how many nodes are using 1 node = 64 cores
#SBATCH -c 16 # how many cores are using
#SBATCH --partition=wzhcnormal
module purge
CaseName=socrates_aquaplanet
/work/home/ac9b0k6rio/miniconda3/bin/conda activate isca_soc2207_env
source ~/isca_17/env.sh

python ${GFDL_BASE}/exp/test_cases/socrates_test/${CaseName}.py

