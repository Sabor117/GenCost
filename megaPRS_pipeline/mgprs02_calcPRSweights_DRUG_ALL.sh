#!/bin/sh
#SBATCH --job-name=mgPRS_02_calcPRSweights_DRUG_ALL
#SBATCH --account=project_2007428
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_02_calcPRSweights_DRUG_ALL.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_02_calcPRSweights_DRUG_ALL.log

printf "\n\n"
printf "==========================================\n\n"
date
printf "\n\n"
printf "==========================================\n\n"

### load ldak

LDAK=/projappl/project_2007428/software/LDAK/ldak6.seb

### Files and names

outDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/"
outName="DRUG_ALL_meta_V4_hapmap_plus_BAYESR"

sumstats="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/hapmap_plus_sumstats/DRUG_ALL_meta_v4_hapmap_plus_mega_sumstats.txt"
snplist="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/hapmap_plus_sumstats/DRUG_ALL_meta_v4_hapmap_plus_mega_snplist.txt"

### Bed files used for version 5.1

ref_genotypes="/scratch/project_2007428/users/FAHagenbeek/data/1000Greferencepanel/ref"

### High LD gene predictors - used for version 5.2

highld_pred="/scratch/project_2007428/users/FAHagenbeek/data/1000Greferencepanel/highld/genes.predictors.used"

### Construct prediction model
### For prior to version 6 include the commented arguments (high-LD argument, bfile argument and window-cm)

printf "\nRunning -mega-prs\n\n===============\n\n"

$LDAK --mega-prs ${outDir}${outName} \
    --model bayesr \
    --summary ${sumstats} \
    --power -0.25 \
    --max-threads 8 \
    --cors /scratch/project_2007428/data/processing/1000G/megaPRS_cors_v6/cors \
    --allow-ambiguous YES \
    --extract ${snplist} #\
#    --bfile ${ref_genotypes} # \ --high-LD ${highld_pred} \ # --window-cm 1 \




printf "\n==========================================\n\n"
printf "Script complete.\n\n"
date
printf "==========================================\n\n"
