#!/bin/bash

#SBATCH --job-name=mt013_find_ld_snps_%a
#SBATCH --time=0-24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --array=1-31
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt013_find_ld_snps_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt013_find_ld_snps_%a.log

##### =========================== #####

### Setup environment

##### =========================== #####

module load r-env
module load plink/1.90b7.7

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron


##### =========================== #####

### Start script

##### =========================== #####

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
outDir="${mainDir}outputs/table/"
refDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/white_british_30k_reference_panel/"

snp_list="${mainDir}processing/misc_data/novel_snps_to_check_ld.txt"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"
printf "\n\n"
printf "Specifying SNP to search for others in LD.\n\n"
printf "\n\n"
printf "\n\n==========================\n\n"

run_no=${SLURM_ARRAY_TASK_ID}

rsid=$(tail -n+2 ${snp_list} | sed -n ${run_no}p | cut -f1)
chrom=$(tail -n+2 ${snp_list} | sed -n ${run_no}p | cut -f2)

printf "\nSNP selected for run #${run_no} is: ${rsid} on chromosome ${chrom}\n\n"
printf "Command to run is:\n"
printf "...\n\n"
echo "plink --bfile ${refDir}ukb22828_chr${chrom}_30k_random_unrelated_white_british --ld-snp ${rsid} --ld-window-r2 0.1 --r2 inter-chr --out ${outDir}${rsid}_ld_snp"

plink --bfile ${refDir}ukb22828_chr${chrom}_30k_random_unrelated_white_british \
        --ld-snp ${rsid} \
        --ld-window-r2 0.1 \
        --r2 inter-chr \
        --out ${outDir}${rsid}_ld_snp

printf "\n\n"
printf "SNPS IN LD DEFINED.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
