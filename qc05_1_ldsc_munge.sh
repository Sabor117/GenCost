#!/bin/bash

#SBATCH --job-name=qc05_1_LDSC_munge_INOUT_ALL_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --array=1-8
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_01_ldsc_munge_INOUT_ALL_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_01_ldsc_munge_INOUT_ALL_%a.log

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

phenotype="INOUT_ALL"

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir="${mainDir}processing/ldsc_intermediate_files/${phenotype}/"
intermediateDir="${mainDir}processing/ldsc_intermediate_files/${phenotype}/"

### File selection

index_no=$((${SLURM_ARRAY_TASK_ID} - 1))

allFiles=(${sumstatDir}*_ldsc_input.txt.gz)

currFile=${allFiles[index_no]}

printf "\n\nFile selected.\n\n"
printf "\n\n==========================\n\n"

printf "\nNow working on:\n\n"
echo "${currFile}"
printf "Phenotype is: ${phenotype}\n"
printf "\n\n"

### Extract basename
file_only=$(basename "${currFile}")

### Extract the first element
#cohort=$(echo "$file_only" | sed "s/_${phenotype}.*//g")
cohort=$(echo "$file_only" | sed "s/_ldsc.*//g")

### Account for snipar FINNGEN variants
if [[ "$cohort" == "FINNGEN" ]]
then
    gwas_format=$(echo "$file_only" | awk -F "${phenotype}_|_ldsc" '{print $2}')

    cohort="${cohort}_${gwas_format}"
fi

### Print the results
echo "File: ${file_only}"
echo "Cohort: $cohort"
printf "\n\n==========================\n\n"

/projappl/project_2007428/software/ldsc/munge_sumstats.py \
    --out "${intermediateDir}${cohort}_${phenotype}_ldsc_munged" \
    --merge-alleles /projappl/project_2007428/software/ldsc/w_hm3.snplist \
    --sumstats ${currFile} \
    --snp rsid \
    --a1 a1 \
    --a2 a0 \
    --p p \
    --N-col n \
    --signed-sumstats zscore1,0 \
    --a1-inc


printf "\nSumstats munged.\n==============\n\n"
printf "Script complete.\n\n"

date





