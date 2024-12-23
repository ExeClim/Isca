#!/bin/bash
#SBATCH --job-name=55cnce_tb_test1 # 作业名称
#SBATCH --output=output_%j.out  # 输出文件名，%j表示当前作业ID
#SBATCH --error=error_%j.err    # 错误输出文件名
#SBATCH -N 1 # how many nodes are using 1 node = 32 cores
#SBATCH -n 32 # how many cores are using
#SBATCH --partition=wzhctdnormal
module purge
ulimit -c unlimited
ulimit -s unlimited
ulimit -l unlimited
CaseName=55cnce_tb_test1
source /work/home/ac9b0k6rio/miniconda3/bin/activate isca_soc2207_env
source /work/home/ac9b0k6rio/isca_17/Isca/src/extra/env/sugon

python ${GFDL_BASE}/exp/run_isca/55cnce_var/${CaseName}.py

