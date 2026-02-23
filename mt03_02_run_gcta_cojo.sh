#!/bin/bash

#SBATCH --job-name=mt03_02_gcta_run_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --array=1-32
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_02_gcta_run_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_02_gcta_run_%a.log

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

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir="${mainDir}outputs/METAL_v4/"
intermediateDir="${mainDir}processing/gcta_intermediate_files_MAF_0_001/"
outDir="${mainDir}outputs/gcta_cojo_v4_MAF_0_001/"

### Reference genotypes

refDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/white_british_30k_reference_panel/"

##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

index_no=$((${SLURM_ARRAY_TASK_ID} - 1))

allFiles=(${metaDir}*.TBL.gz)
currFile=${allFiles[index_no]}

#gzip ${currFile}

fileName=$(basename "$currFile")
filePrefix=${fileName%1.TBL.gz}
fileoutPrefix=${fileName%_metal_output_1.TBL.gz}

printf "Command is:\n"
echo "Rscript ${mainDir}scripts/meta_step3_2_make_gcta_input.R $fileName"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step3_2_make_gcta_input.R $fileName



##### =========================== #####

### Running GCTA

##### =========================== #####

workingfile="${intermediateDir}${filePrefix}gcta_input.txt"

for curr_chr in {1..23}; do

    if [ "$curr_chr" -eq 23 ];
    then
        chr_name="X"
    else
        chr_name=$curr_chr
    fi

    echo -e "\n\nLiftOver command is:\n"
    echo "/scratch/project_2007428/users/Zhiyu/Tool/GCTA/gcta-1.94.1 --bfile ${refDir}/ukb22828_chr${curr_chr}_30k_random_unrelated_white_british  --chr ${curr_chr} --maf 0.001 --cojo-file ${workingfile} --cojo-slct --out ${outDir}${fileoutPrefix}_chr${chr_name}_out"

    /scratch/project_2007428/users/Zhiyu/Tool/GCTA/gcta-1.94.1 \
        --bfile ${refDir}/ukb22828_chr${curr_chr}_30k_random_unrelated_white_british  \
        --chr ${curr_chr} \
        --maf 0.001 \
        --cojo-file ${workingfile} \
        --cojo-slct \
        --out ${outDir}${fileoutPrefix}_chr${chr_name}_out

    echo -e "\n\nChromosome ${chr_name} complete.\n\n=============\n\n"

done

printf "\n\n==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"




