#!/bin/sh
#SBATCH --job-name=mt09_3_calcPRS_%a
#SBATCH --account=project_2007428
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=8-11
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt09_3_calcPRS_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt09_3_calcPRS_%a.log

module load plink

printf "\n\n"
printf "==========================================\n\n"
date
printf "\n\n"
printf "==========================================\n\n"

genoDir="/scratch/project_2007428/data/processing/ukbb_78537/genotypes/plink2/"

prsFiles=("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/DRUG_ALL_noUKB_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/IN_ALL_noUKB_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/PRIM_ALL_noUKB_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/INOUT_ALL_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/DRUG_ALL_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/IN_ALL_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/PRIM_ALL_meta_V4_hapmap_plus_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights_outdated/DS_bayes.IN_ALL.effects"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights_outdated/DS_bayes.IN_ALL_noUKB.effects"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights_outdated/DS_bayes.DRUG_ALL.effects"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights_outdated/DS_bayes.INOUT_ALL.effects"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/DRUG_ALL_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/IN_ALL_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/INOUT_ALL_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/PRIM_ALL_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/IN_ALL_noUKB_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/DRUG_ALL_noUKB_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/megaprs_weights/PRIM_ALL_noUKB_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt"
            )

prs_index_no=$((${SLURM_ARRAY_TASK_ID} - 1))

prs_select="${prsFiles[prs_index_no]}"

outDir="/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/prs/scores/"

for chrom in {1..22}

    do

        echo "Running on chromosome ${chrom}. Command is:"
        echo "plink2 --pfile $genoDir/ukb22828_chr_${chrom} --score $prs_select 9 2 5 ignore-dup-ids --out ${outDir}"

        analysis_name=$(basename "$prs_select" | sed -E 's/.*FUSION_(.+)\.txt/\1/')
        analysis_name=$(basename "$prs_select" .effects)

        if [[ "$analysis_name" == *DS_* ]]; then
            snp=6
            allele=2
            weight=5
        else
            snp=9
            allele=2
            weight=5
        fi

        echo "plink2 --pfile $genoDir/ukb22828_chr_${chrom} \
                --score $prs_select $snp $allele $weight ignore-dup-ids \
                --out ${outDir}${analysis_name}_chr${chrom}"

        plink2 --pfile "$genoDir/ukb22828_chr_${chrom}" \
                --score "$prs_select" $snp $allele $weight ignore-dup-ids \
                --out "${outDir}${analysis_name}_chr${chrom}"

        echo "${analysis_name} scores made for chromosome ${chrom}."
        printf "\n\n=========\n\n"

    done


printf "\n==========================================\n\n"
printf "Script complete.\n\n"
date
printf "==========================================\n\n"
