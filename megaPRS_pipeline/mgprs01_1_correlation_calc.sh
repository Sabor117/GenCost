#!/bin/bash

#SBATCH --job-name=mgPRS_01_1_correlation_calc
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_01_1_correlation_calc.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_01_1_correlation_calc.log

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


##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts. NOBODY EXPECTS THE SEB INQUISITION:\n\n"
date
printf "\n\n==========================\n\n"

### load ldak

LDAK5=/projappl/project_2007428/software/LDAK/ldak5.1.linux
LDAK6=/projappl/project_2007428/software/LDAK/ldak6.seb

### Consistent files and directories

ref_genotype=/scratch/project_2007428/users/FAHagenbeek/data/1000Greferencepanel/ref
ldak_sections=/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/sections/
ldak_annotations=/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/megaPRS_annotations/
summary_Dir=/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/

snplist=/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/IN_ALL_meta_v4_hapmap_plus_mega_snplist.txt

### Construct correlations model for LDAK VERSION 6.0

$LDAK6 --calc-cors /scratch/project_2007428/data/processing/1000G/megaPRS_cors_v6/cors --bfile ${ref_genotype} --extract ${snplist}

### Cut LDAK sections for MegaPRS - VERSION 5.1

#$LDAK5 --cut-weights ${ldak_sections} --bfile ${ref_genotype}
#$LDAK5 --calc-weights-all ${ldak_sections} --bfile ${ref_genotype}
#mv ${ldak_sections}weights.short bld65

### Used BLD annotations and new bld65 (from cut sections) to make summary matrix - VERSION 5.1

#$LDAK5 --calc-tagging ${ldak_annotations}ldak_bld_matrix \
#          --bfile ${ref_genotype} \
#          --power -.25 \
#          --annotation-number 65 \
#          --ignore-weights YES \
#          --window-kb 1 \
#          --annotation-prefix ${ldak_annotations}bld \
#          --save-matrix YES

#$LDAK5 --sum-hers ${ldak_annotations}ldak_ind_hers_DRUG_ALL --tagfile ${ldak_annotations}ldak_bld_matrix.tagging --summary ${summary_Dir}DRUG_ALL_meta_sumstats_for_megaPRS.txt --matrix ${ldak_annotations}ldak_bld_matrix.matrix --check-sums NO

#$LDAK5 --sum-hers ${ldak_annotations}ldak_ind_hers_IN_ALL --tagfile ${ldak_annotations}ldak_bld_matrix.tagging --summary ${summary_Dir}IN_ALL_meta_sumstats_for_megaPRS.txt --matrix ${ldak_annotations}ldak_bld_matrix.matrix --check-sums NO

#$LDAK5 --sum-hers ${ldak_annotations}ldak_ind_hers_IN_ALL_noFINNGEN --tagfile ${ldak_annotations}ldak_bld_matrix.tagging --summary ${summary_Dir}IN_ALL_noFINNGEN_meta_sumstats_for_megaPRS.txt --matrix ${ldak_annotations}ldak_bld_matrix.matrix --check-sums NO

#$LDAK5 --sum-hers ${ldak_annotations}ldak_ind_hers_IN_ALL_noUKB --tagfile ${ldak_annotations}ldak_bld_matrix.tagging --summary ${summary_Dir}IN_ALL_noUKB_meta_sumstats_for_megaPRS.txt --matrix ${ldak_annotations}ldak_bld_matrix.matrix --check-sums NO

#$LDAK5 --sum-hers ${ldak_annotations}ldak_ind_hers_INOUT_ALL --tagfile ${ldak_annotations}ldak_bld_matrix.tagging --summary ${summary_Dir}INOUT_ALL_meta_sumstats_for_megaPRS.txt --matrix ${ldak_annotations}ldak_bld_matrix.matrix --check-sums NO


printf "\n\n"
printf "Correlations created.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
