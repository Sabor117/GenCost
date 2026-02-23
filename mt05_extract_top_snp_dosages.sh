#!/bin/bash

#SBATCH --job-name=mt05_get_costs_per_snp_%a
#SBATCH --time=1-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt05_get_costs_per_snp_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt05_get_costs_per_snp_%a.log


##### =========================== #####

### Setup environment

### Step 5.1 runs Rscript inputting GCTA-annotated top SNPs and meta-analysis
### Outputs the top SNPs and Z-scores + P-vals for each SNP

### Step 5.2 inputs top SNP list into Plink command to extract dosages

##### =========================== #####

module load r-env
module load plink/2.00a5

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron


##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "Extracting top SNPs:\n\n"
printf "First command is:\n"
echo "Rscript ${mainDir}scripts/meta_step5_1_get_top_snps.R"
printf "\n\n"
printf "==========================================\n\n"

##### ====== #####
### Step 5.1
##### ====== #####

#Rscript ${mainDir}scripts/meta_step5_1_get_top_snps.R

printf "\nTop SNPs extracted. Running Plink.\n==============\n\n"


##### ====== #####
### Step 5.2
##### ====== #####

phenotype="DRUG_ALL"
maf_filter="MAF_0_001"

snpListDir="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/cost_per_snp_processing/"
snpListPattern="${phenotype}_top_snps_${maf_filter}_list_chr"

pgenDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/plink2/"
pgenPrefix="ukb22828_chr_"

outDosages="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/cost_per_snp_processing/${phenotype}_top_snps_${maf_filter}_dosages_"

### Loop through each SNP list file

for snpListFile in ${snpListDir}/${snpListPattern}*.txt; do

    ### Extract the chromosome number from the file name

    chrom=$(basename ${snpListFile} | sed -E "s/${snpListPattern}([0-9]+)\.txt/\1/")
    
    ### Construct the full path to the .pgen files

    pfilePath="${pgenDir}${pgenPrefix}${chrom}"
    
    ### Construct the output file prefix

    outFile="${outDosages}chr${chrom}"

    ### Run the plink2 command
    
    printf "Plink2 command is:\n\n"
    echo "plink2 --pfile ${pfilePath} --export A --extract ${snpListFile} --out ${outFile}"
    printf "____\n\n"
    plink2 --pfile ${pfilePath} --export A  --extract ${snpListFile} --out ${outFile}

done

printf "==========================================\n\nScript complete.\n\n===========\n\n"
date





