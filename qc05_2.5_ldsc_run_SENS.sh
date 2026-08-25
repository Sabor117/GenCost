#!/bin/bash

#SBATCH --job-name=qc05_02.5_LDSC_run_SENS
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_02.5_ldsc_run_SENS.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_02.5_ldsc_run_SENS.log
printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

### Load LDSC Tykky environment
module purge
module load tykky

### Path to LDSC environment
export PATH="/projappl/project_2007428/software/tykky_envs/ldsc_env/bin/:$PATH"

### Main directories
mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir="${mainDir}processing/ldsc_intermediate_files/SENSITIVITY_TEST/"

### Get all munged summary statistics
allFiles=("${sumstatDir}"*munged.sumstats.gz)

### Iterate through each file as the starting file
for startingFile in "${allFiles[@]}"; do

    ### Extract basename
    file_only=$(basename "${startingFile}")

    ### Only UKB EUR
    cohort=$(echo "$file_only" | cut -d'_' -f2,3)

    ### Build list of all other files
    restFiles=()
    for file in "${allFiles[@]}"; do
        [[ "$file" != "$startingFile" ]] && restFiles+=("$file")
    done

    ### Join files into comma-separated string
    filesString=$(IFS=,; echo "${startingFile},${restFiles[*]}")

    echo "Cohort = ${cohort}"
    echo "File = ${startingFile}"
    printf "\n\n"

    ### Run LDSC
    /projappl/project_2007428/software/ldsc/ldsc.py \
        --rg ${filesString} \
        --ref-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --w-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --out "${sumstatDir}${cohort}_SENSITIVITY_TEST_ldsc_corr"

    printf "\n==============\n\n"

done

printf "SCRIPT COMPLETE.\n\n"
date
printf "\n\n==========================\n\n"




