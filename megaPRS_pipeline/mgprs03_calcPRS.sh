#!/bin/sh
#SBATCH --job-name=mgPRS_03_calcPRS
#SBATCH --account=project_2007428
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_03_calcPRS_noUKB.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_03_calcPRS_noUKB.log

module load plink

printf "\n\n"
printf "==========================================\n\n"
date
printf "\n\n"
printf "==========================================\n\n"

genoDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/plink2/"
weightsFile="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/DS_bayes.IN_ALL_noUKB.effects"
outDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/scores/"

for chrom in {1..22}

    do

        echo "Running on chromosome ${chrom}. Command is:"
        echo "plink2 --pfile $genoDir/ukb22828_chr_${chrom} --score $weightsFile 6 2 5 ignore-dup-ids --out $outDir"

        plink2 --pfile ${genoDir}/ukb22828_chr_${chrom} --score ${weightsFile} 6 2 5 ignore-dup-ids --out ${outDir}IN_ALL_noUKB_megaprs_chr${chrom}

        echo "MegaPRS scores made."
        printf "\n\n=========\n\n"

        echo "Running LDpred2 scores."
        echo "Command is:"
        echo "plink2 --pfile ${genoDir}/ukb22828_chr_${chrom} --score ${weightsFile} 1 4 6 ignore-dup-ids --out ${outDir}IN_ALL_noUKB_ldpred2_chr${chrom}"

        plink2 --pfile ${genoDir}/ukb22828_chr_${chrom} --score ${weightsFile} 1 4 6 ignore-dup-ids --out ${outDir}IN_ALL_noUKB_ldpred2_chr${chrom}


    done


printf "\n==========================================\n\n"
printf "Script complete.\n\n"
date
printf "==========================================\n\n"
