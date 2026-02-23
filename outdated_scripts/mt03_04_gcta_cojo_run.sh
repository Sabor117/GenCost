#!/bin/bash

#SBATCH --job-name=mt03_gctacojo_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --array=1-40
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_gctacojo_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_gctacojo_%a.log

##### =========================== #####

### Setup environment

##### =========================== #####

module load r-env

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron

### Directories and files

metaDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL/"
fileList=(${metaDir}*.gz)

gctaDir="/scratch/project_2007428/users/Zhiyu/Tool/GCTA/"

### Array job

run_no=${SLURM_ARRAY_TASK_ID}


##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

### Reduce run_no by 1 for selection

fileselect=$((run_no - 1))

### Check if the selected number is within the range of available files and select it

if [ $fileselect -ge 0 ] && [ $fileselect -lt ${#fileList[@]} ]; then
    
    currMeta=${fileList[$fileselect]}
    echo "Selected file: $currMeta"
    printf "=======\n\n"
    
else

    printf "Invalid selection. Please choose a number within the range.\n\n"

fi

### Run GCTA on a per-chromosome basis

${gctaDir}gcta-1.94.1 


printf "\n\n==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

