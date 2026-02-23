#!/bin/bash

#SBATCH --job-name=METAL_GWAS_<OUT_PREFIX>
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/METAL/a01_METAL_<OUT_PREFIX>.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/METAL/a01_METAL_<OUT_PREFIX>.log

### METAL run file for <CURRENT_TRAIT>

curr_analysis_file=<CURRENT_FILE>

curr_analysis_file_heritability=<CURRENT_FILE_HERIT>

outfile=<CURRENT_OUTFILE>

outfile_heritability=<CURRENT_OUTFILE_HERIT>

printf "Script settup complete. Starting run.\n\n"
date
printf "\n\n==========================\n\n"
printf "Nobody expects the Seb Inqusition!\n\n"
printf "\n\n==========================\n\n"

/projappl/project_2007428/software/development-metal/METAL/build/bin/metal ${curr_analysis_file}

/projappl/project_2007428/software/development-metal/METAL/build/bin/metal ${curr_analysis_file_heritability}

printf "METAL run complete.\n\n"
date

gzip ${outfile}
gzip ${outfile_heritability}

printf "Files Gzipped.\n\n"
printf "Script complete.\n\n"
date
printf "\n\n==========================\n\n"



