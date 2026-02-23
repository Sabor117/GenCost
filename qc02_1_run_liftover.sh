#!/bin/bash

#SBATCH --job-name=Liftover_run
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc02_1_liftover_run.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc02_1_liftover_run.log

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cohort_sumstats/"

chain_file="/projappl/project_2007428/software/liftOver/chain_files/hg19ToHg38.over.chain.gz"

for curr_file in ${mainDir}*/*_liftOver_input.tsv;
    
    do
    echo -e "\n\n========================\n\nRunning LiftOver for cohort:\n"
    currCohort=$(basename ${curr_file%_liftOver_input.tsv} | tr '[a-z]' '[A-Z]')
    echo ${currCohort}
    outDir="$(dirname "${curr_file}")""/"
    echo -e "\n\nCurrent input file:\n"
    echo ${curr_file}
    echo -e "\n\nLiftOver command is:\n"
    echo "/projappl/project_2007428/software/liftOver/liftOver ${curr_file} ${chain_file} ${outDir}${currCohort}_liftOver_output.out ${outDir}liftOver.err"

    /projappl/project_2007428/software/liftOver/liftOver ${curr_file} ${chain_file} ${outDir}${currCohort}_liftOver_output.out ${outDir}liftOver.err

    echo -e "\n\nRun complete.\n\n================================\n\n"

done






