#!/bin/bash

#SBATCH --job-name=mt11_1_coloc_preprocess
#SBATCH --time=3:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt11_1_coloc_preprocess.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt11_1_coloc_preprocess.log

##### =========================== #####

### Setup environment

##### =========================== #####

module load r-env
module load samtools

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron

### Directories

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir="${mainDir}outputs/METAL_v4/"
miscDataDir="${mainDir}processing/misc_data/"


##### =========================== #####

### Starting pre-process

##### =========================== #####

### Download SNP information files from NCBI (build 37 and 38) 

# wget -P /scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz 
# wget -P /scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz

builg37_vcf="${miscDataDir}All_20180423.vcf.gz"
builg38_vcf="${miscDataDir}GCF_000001405.40.gz"

### Run script to combine the main GWAS SNPs into one file for getting position data

printf "First script starts:\n\n"
printf "NOBODY EXPECTS THE SEB INQUISTION:\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/meta_step11_1_get_gwas_snps.R"
printf "\n\n"
printf "==========================================\n\n"

#Rscript ${mainDir}scripts/meta_step11_1_get_gwas_snps.R

printf "\n\n"
printf "==========================================\n\n"
echo "Pre-process script 11-1 run."
date
printf "\n\n==========================\n\n"

### Filter VCF using the NCBI files to get build 37 positions for rsIDs
### First use tabix to make an index file for the VCF

echo -e "\n============\nSTARTING VCF TRANSFER.\n\n================================\n\n"

echo "tabix -p vcf ${builg37_vcf}"
#tabix -p vcf ${builg37_vcf}

echo -e "\n\nTabix file created.\n......\n\n"

echo "bcftools query -i 'ID=@'"${metaDir}ALL_METAL_output_RSIDs_only.txt" -f '%ID\t%CHROM\t%POS\n' ${builg37_vcf} > ${metaDir}ALL_METAL_output_RSIDs_to_hg37.tsv"
#bcftools query -i 'ID=@'"${metaDir}ALL_METAL_output_RSIDs_only.txt" -f '%ID\t%CHROM\t%POS\n' ${builg37_vcf} > ${metaDir}ALL_METAL_output_RSIDs_to_hg37.tsv

echo -e "\n\nBCFTOOLs run.\n......\n\n"

### Created BED file for SNPs with the Rscript
### Testing with pos-1/pos, pos/pos and pos/pos+1
### Hopefully limited difference between files: pos-1/pos is strictly correct

echo -e "\n============\nSTARTING LIFTOVER.\n\n================================\n\n"

### Testing shows that liftOver functions with any of: pos-1/pos, pos/pos, pos/pos+1

#pos0_pos_bed="${metaDir}ALL_METAL_output_SNPIDs.bed"
#pos_pos_bed="${metaDir}ALL_METAL_output_SNPIDs_pos_pos.bed"
pos_pos1_bed="${metaDir}ALL_METAL_output_SNPIDs_pos_pos1.bed"

chain_file="/projappl/project_2007428/software/liftOver/chain_files/hg38ToHg19.over.chain.gz"

#echo "/projappl/project_2007428/software/liftOver/liftOver ${pos0_pos_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos0_pos_liftOver_output.out ${metaDir}liftOver_pos0_pos.err"
#/projappl/project_2007428/software/liftOver/liftOver ${pos0_pos_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos0_pos_liftOver_output.out ${metaDir}liftOver_pos0_pos.err
#echo -e "\n\nPos0-pos complete.\n\n......\n\n"

#echo "/projappl/project_2007428/software/liftOver/liftOver ${pos_pos_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos_pos_liftOver_output.out ${metaDir}liftOver_pos_pos.err"
#/projappl/project_2007428/software/liftOver/liftOver ${pos_pos_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos_pos_liftOver_output.out ${metaDir}liftOver_pos_pos.err
#echo -e "\n\nPos-pos complete.\n\n......\n"

echo "/projappl/project_2007428/software/liftOver/liftOver ${pos_pos1_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos_pos1_liftOver_output.out ${metaDir}liftOver_pos_pos1.err"
#/projappl/project_2007428/software/liftOver/liftOver ${pos_pos1_bed} ${chain_file} ${metaDir}ALL_METAL_output_SNPIDs_pos_pos1_liftOver_output.out ${metaDir}liftOver_pos_pos1.err
echo -e "\n\nPos-pos1 complete.\n\n......\n"

printf "\n\n"
printf "==========================================\n\n"
echo "LiftOver and conversion phase 1 run."
date
printf "\n\n==========================\n\n"

