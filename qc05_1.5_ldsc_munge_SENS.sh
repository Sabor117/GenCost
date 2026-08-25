#!/bin/bash

#SBATCH --job-name=qc05_1.5_LDSC_munge_SENS_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --array=1-9
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_01.5_ldsc_munge_SENS_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_01.5_ldsc_munge_SENS_%a.log

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

### Load LDSC Tykky environment
### Unsure if needed to load Tykky

module purge
module load tykky

### Path to environment

export PATH="/projappl/project_2007428/software/tykky_envs/ldsc_env/bin/:$PATH"

### Main directories

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir="${mainDir}processing/ldsc_intermediate_files/SENSITIVITY_TEST/"
intermediateDir="${mainDir}processing/ldsc_intermediate_files/SENSITIVITY_TEST/"

### File selection

index_no=$((${SLURM_ARRAY_TASK_ID} - 1))

allFiles=("${sumstatDir}"*_ldsc_input.txt.gz)

currFile="${allFiles[$index_no]}"

printf "\n\nFile selected.\n\n"
printf "==========================\n\n"

echo "Now working on:"
echo "${currFile}"
echo

### Extract filename information

file_only=$(basename "${currFile}")

if [[ "$file_only" =~ ^([^_]+)_([^_]+)_(ALL|NO_ZERO|RAW_COST)_ldsc_input\.txt\.gz$ ]]; then
    cohort="${BASH_REMATCH[1]}"
    phenotype="${BASH_REMATCH[2]}"
    stratif="${BASH_REMATCH[3]}"
else
    echo "ERROR: Filename does not match expected format:"
    echo "${file_only}"
    exit 1
fi

echo "File: ${file_only}"
echo "Cohort: ${cohort}"
echo "Phenotype: ${phenotype}"
echo "Stratification: ${stratif}"

printf "\n\n==========================\n\n"

output="${intermediateDir}${cohort}_${phenotype}_${stratif}_ldsc_munged"

echo "Output:"
echo "${output}.sumstats.gz"
echo

/projappl/project_2007428/software/ldsc/munge_sumstats.py \
    --out "${output}" \
    --merge-alleles /projappl/project_2007428/software/ldsc/w_hm3.snplist \
    --sumstats "${currFile}" \
    --snp rsid \
    --a1 a1 \
    --a2 a0 \
    --p p \
    --N-col n \
    --signed-sumstats zscore1,0 \
    --a1-inc

status=$?

if [[ $status -ne 0 ]]; then
    echo
    echo "ERROR: LDSC munge_sumstats.py failed with exit code ${status}"
    exit $status
fi

echo
echo "Sumstats munged."
echo "==============="
echo
echo "Script complete."
date



