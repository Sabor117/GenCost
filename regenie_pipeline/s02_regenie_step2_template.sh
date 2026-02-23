#!/bin/bash

#SBATCH --job-name=a02_REGENIE_step_2_GWAS_<OUT_PREFIX>_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=32G
#SBATCH --account=project_2007428
#SBATCH --array=1-23
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/regenie_logs/a02_regenie_step2_<OUT_PREFIX>_%A_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/regenie_logs/a02_regenie_step2_<OUT_PREFIX>_%A_%a.log

analysis_name=<OUT_PREFIX>

### Genotype bgen file

genotypeDir="/scratch/project_2007428/data/base_data/ukbb_78537/genotypes/bgen/"

if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]

    then

        genotype_file="ukb22828_chr_X.bgen"
        index_file="ukb22828_chr_X.bgen.bgi"

        ### Other genotype files

        sample_file="ukb22828_chr_X.sample"

    else

        genotype_file="ukb22828_chr_"${SLURM_ARRAY_TASK_ID}".bgen"
        index_file="ukb22828_chr_"${SLURM_ARRAY_TASK_ID}".bgen.bgi"

        ### Other genotype files

        sample_file="ukb22828.sample"

fi

### Phenotype file

phenotypeDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/"
phenotype_file="<PHENOTYPE_FILE>"

### Location of covariate file

covarDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/"
covariate_file=<COVARIATE_FILE>

### Predcitors file

predictorsDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/regenie/step1_output/"${analysis_name}"/"
predictors_file=${analysis_name}"_pred.list"

### Output prefix

out_prefix="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/regenie/step2_output/"${analysis_name}"/"
mkdir -p ${out_prefix}

printf "Script settup complete. Starting run.\n\n"
date
printf "\n\n==========================\n\n"
printf "Nobody expects the Seb Inqusition!\n\n"
printf "\n\n==========================\n\n"

/projappl/project_2007428/software/regenie/regenie_*_Centos7_mkl \
	--step 2 \
	--bgen ${genotypeDir}${genotype_file} \
    --bgi ${genotypeDir}${index_file} \
    --ref-first \
    --sample ${genotypeDir}${sample_file} \
	--phenoFile ${phenotypeDir}${phenotype_file} \
	--covarFile ${covarDir}${covariate_file} \
	--qt \
    --firth --approx --pThresh 0.01 \
    --pred ${predictorsDir}${predictors_file} \
    --bsize 400 \
    --out ${out_prefix}${analysis_name}"_regenie_step2_chr"${SLURM_ARRAY_TASK_ID}
