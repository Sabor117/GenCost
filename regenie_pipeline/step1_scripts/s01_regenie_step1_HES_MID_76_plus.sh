#!/bin/bash

#SBATCH --job-name=REGENIE_step_1_GWAS_HES_MID_76_plus
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=32G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/a01_regenie_step1_HES_MID_76_plus.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/a01_regenie_step1_HES_MID_76_plus.log

### Genotype Plink file

genotypeDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/full_plink_array/"
genotype_file="ukb22418_allChrs"

### Phenotype file

phenotypeDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/"
phenotype_file="FINAL_hes_based_phenotype_frame_REGENIE.tsv"

### Location of covariate file

covarDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/"

### Location of QC'd SNPs based on REGENIE tutorial

extract_file="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/plink_superpop_qc/MID_qc_snps_for_regenie.snplist"

### Location of file of IIDs to keep

keep_file="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_76_plus_iids_only.tsv"

### tmpdir location if using low-mem argument

tmpDir_for_lowmem="/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/regenie_tmp_preds/HES_MID_76_plus/"
mkdir -p ${tmpDir_for_lowmem}

### Output prefix

out_prefix="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/regenie/step1_output/HES_MID_76_plus/"
mkdir -p ${out_prefix}

printf "Script settup complete. Starting run.\n\n"
date
printf "\n\n==========================\n\n"
printf "Nobody expects the Seb Inqusition...\n\n"
printf "\n\n==========================\n\n"

analysis_name="HES_MID_76_plus"

/projappl/project_2007428/software/regenie/regenie_*_Centos7_mkl \
	--step 1 \
	--bed ${genotypeDir}${genotype_file} \
	--phenoFile ${phenotypeDir}${phenotype_file} \
	--covarFile ${covarDir}gwas_cost_covariates_population_MID.tsv \
	--extract ${extract_file} \
	--keep ${keep_file} \
	--qt \
	--bsize 1000 \
	--lowmem \
	--lowmem-prefix ${tmpDir_for_lowmem} \
	--out ${out_prefix}${analysis_name}
