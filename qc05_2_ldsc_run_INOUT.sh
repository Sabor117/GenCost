#!/bin/bash

#SBATCH --job-name=qc05_02_LDSC_run_INOUT_ALL
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_02_ldsc_run_INOUT_ALL.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc05_02_ldsc_run_INOUT_ALL.log

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

### Load LDSC Tykky environment
### Unsure if needed to load Tykky

module purge
module load tykky

### Path to environment

export PATH="/projappl/project_2007428/software/tykky_envs/ldsc_env/bin/:$PATH"

phenotype="INOUT_ALL"

### Main directories

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir="${mainDir}processing/ldsc_intermediate_files/${phenotype}/"

allFiles=(${sumstatDir}*.sumstats.gz)

### Iterate through each file as the starting file

for startingFile in "${allFiles[@]}"; do

    ### Extract basename
    file_only=$(basename "${startingFile}")

    ### Extract the first element
    cohort=$(echo "$file_only" | cut -d'_' -f1)

    ### Add statement to include UKB multiple ancestries

    if [[ $cohort == *UKB* ]]; then

        echo "Cohort is UKB."

        population=$(echo "$file_only" | cut -d'_' -f2)

        cohort="${cohort}_${population}"

    else

        printf "Cohort is not UKB. These people have no ancestry.\n\n"

    fi

    ### Account for snipar FINNGEN variants

    if [[ "$cohort" == "FINNGEN" ]]

    then

        gwas_format=$(echo "$file_only" | awk -F "_${phenotype}" '{print $1}')

        cohort="${gwas_format}"
    
    else

        printf "Cohort is not FINNGEN. These people have no parents.\n\n"
    
    fi

    ### Remove the starting file from the array

    restFiles=(${allFiles[@]/$startingFile})

    ### Join the files into a comma-separated string

    filesString=$(IFS=,; echo "${startingFile},${restFiles[*]}")

    echo "Cohort = ${cohort}"
    echo "File = ${startingFile}"
    printf "\n\n"

    ### Run the specific command with the current filesString

    /projappl/project_2007428/software/ldsc/ldsc.py \
        --rg ${filesString} \
        --ref-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --w-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --out /scratch/project_2007428/projects/prj_001_cost_gwas/processing/ldsc_intermediate_files/INOUT_ALL/${cohort}_INOUT_ALL_ldsc_corr

    printf "\n==============\n\n"

done

printf "SCRIPT COMPLETE.\n\n"
date
printf "\n\n==========================\n\n"





