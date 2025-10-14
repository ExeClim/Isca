#!/bin/bash

# === CONFIGURATION ===
SCRIPT="RT_sst.py"
SBATCH_TEMPLATE_DIR="./sbatch_jobs"
mkdir -p "$SBATCH_TEMPLATE_DIR"

# === ENVIRONMENT SETUP BLOCK ===
read -r -d '' ENV_BLOCK << 'EOF'
export GFDL_BASE=/home/philbou/Isca
export GFDL_ENV=narval.ifort
export GFDL_WORK=/scratch/philbou/isca_work
export GFDL_DATA=/scratch/philbou/isca_data
source /home/philbou/.bashrc
conda activate isca_env
cd $GFDL_BASE/exp/test_cases/realistic_continents
EOF

# === GENERATE AND SUBMIT JOBS FOR RANGE -4 TO 4 (INCLUDE 0) ===
i_values=(-1)
newvals=(295)

for idx in "${!i_values[@]}"; do
    i="${i_values[$idx]}"
    newval="${newvals[$idx]}"
    if (( i < 0 )); then
        abs=${i#-}            # absolute value of negative number
        JOB_SUFFIX="m${abs}"  # e.g. m4, m3, m2, m1
    else
        JOB_SUFFIX="${i}"     # e.g. 0, 1, 2, 3, 4
    fi

    JOB_NAME="RT42_sst_${JOB_SUFFIX}"
    JOB_FILE="$SBATCH_TEMPLATE_DIR/job_${JOB_SUFFIX}.sbatch"

    cat <<EOF > "$JOB_FILE"
#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=3G
#SBATCH --time=2-00:00
#SBATCH --job-name=$JOB_NAME
#SBATCH --output=/scratch/philbou/outerr/%x-%j.out
#SBATCH --error=/scratch/philbou/outerr/%x-%j.err
#SBATCH --account=def-rfajber
#SBATCH --mail-user=philippe.boulanger@mail.mcgill.ca
#SBATCH --mail-type=ALL

$ENV_BLOCK

python $SCRIPT $i $newval
EOF

    sbatch "$JOB_FILE"
    echo "Submitted $JOB_FILE with argument $i (job name: $JOB_NAME)"
done
