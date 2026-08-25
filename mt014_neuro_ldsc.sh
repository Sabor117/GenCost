#!/bin/bash

#SBATCH --job-name=mt014_neuro_ldsc
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt14_neuro_ldsc.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt14_neuro_ldsc.log

##### =========================== #####

### Setup environment

##### =========================== #####

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

phenotype="neuroticism"

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
intermediateDir="${mainDir}processing/ldsc_intermediate_files/misc/"
metaDir="${mainDir}processing/ldsc_intermediate_files/META_v4/"

currFile="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/neuroticism_ldsc_input.txt.gz"


##### =========================== #####

### Munge run

##### =========================== #####

printf "\nNow working on:\n\n"
echo "${currFile}"
printf "Phenotype is: ${phenotype}\n"
printf "\n\n"

### Extract basename

file_only=$(basename "${currFile}")

### Print for debugging
echo "File: ${file_only}"
echo "Phenotype: ${phenotype}"
printf "\n\n==========================\n\n"

/projappl/project_2007428/software/ldsc/munge_sumstats.py \
    --out "${intermediateDir}/${phenotype}_ldsc_munged" \
    --merge-alleles /projappl/project_2007428/software/ldsc/w_hm3.snplist \
    --sumstats "${currFile}" \
    --snp snpid \
    --a1 a1 \
    --a2 a0 \
    --p p \
    --N-col n \
    --signed-sumstats zscore1,0 \
    --a1-inc


printf "\nSumstats munged.\n==============\n\n"

##### =========================== #####

### LDSC run

##### =========================== #####

munged_file="${intermediateDir}/${phenotype}_ldsc_munged.sumstats.gz"

for meta_file in "${metaDir}"/*_META_v4_META_v4_ldsc_munged.sumstats.gz; do

    meta_base=$(basename "${meta_file}" _META_v4_META_v4_ldsc_munged.sumstats.gz)

    printf "\nRunning LDSC rg:\n"
    printf "Trait 1: %s\n" "${munged_file}"
    printf "Trait 2: %s\n\n" "${meta_file}"

    /projappl/project_2007428/software/ldsc/ldsc.py \
        --rg "${munged_file},${meta_file}" \
        --ref-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --w-ld-chr /projappl/project_2007428/software/ldsc/Ref/eur_w_ld_chr/ \
        --out "${intermediateDir}/${phenotype}_${meta_base}"

done

printf "Nobody expects the Seb Inquistion!\n\n"

printf "SCRIPT COMPLETE.\n\n"
date
printf "\n\n==========================\n\n"




