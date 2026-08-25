#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library(dplyr)

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Gonna make REGENIE scripts like magic snort snort...")

input_step1_template_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/regenie_pipeline/s01_regenie_step1_template.sh"
input_step2_template_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/regenie_pipeline/s02_regenie_step2_template.sh"

step1Out = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/regenie_pipeline/step1_scripts/"
step2Out = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/regenie_pipeline/step2_scripts/"


##############################

### Making output scripts

##############################

project_name = "GenCOST"

### Based on all of the REGENIE input files

analysis_list = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/*REGENIE*")

### Exclude specific files based on requirements

exclude_string = c("CPI|frame")
analysis_list = analysis_list[!(grepl(exclude_string, analysis_list, ignore.case = TRUE))]

### Define the PanUKB populations used

ukb_pops = c("EUR", "AFR", "CSA", "EAS")

input_step1_template = readLines(input_step1_template_file)
input_step2_template = readLines(input_step2_template_file)

for (i in 1:length(analysis_list)){

    curr_input_step1_template = input_step1_template
    curr_input_step2_template = input_step2_template

    curr_analysis = analysis_list[i]

    analysis_name = basename(curr_analysis)

    analysis_name_basic = case_when(grepl("INPAT", analysis_name) ~ "IN",
                                    grepl("DRUG", analysis_name) ~ "DRUG",
                                    grepl("PRIM", analysis_name) ~ "PRIM",
                                    .default = analysis_name)

    stratif_name = case_when(grepl("non_zero", analysis_name) ~ "NON_ZERO",
                                    grepl("raw_cost", analysis_name) ~ "RAW_COST",
                                    .default = NA)

    if (is.na(stratif_name)){

        analysis_name = paste0(analysis_name_basic, "_ALL")

    } else {

        analysis_name = paste0(analysis_name_basic, "_", stratif_name)

    }

    ### Phenotype file = input file

    curr_input_step1_template = gsub("<PHENOTYPE_FILE>", basename(curr_analysis), curr_input_step1_template)
    curr_input_step2_template = gsub("<PHENOTYPE_FILE>", basename(curr_analysis), curr_input_step2_template)

    for (j in 1:length(ukb_pops)){

        currPop = ukb_pops[j]

        analysis_prefix = paste0(analysis_name, "_", currPop)

        ### Out prefix = newly defined analysis name

        curr_input_step1_template = gsub("<OUT_PREFIX>", analysis_prefix, curr_input_step1_template)
        curr_input_step2_template = gsub("<OUT_PREFIX>", analysis_prefix, curr_input_step2_template)

        ### Keep IIDs based on population
        ### This would need to be adjusted if re-running other stratifications

        keep_file = paste0("/scratch/project_2007428/data/processing/ukbb_78537/genotypes/plink_superpop_qc/", currPop, "_sample.tsv")

        curr_input_step1_template = gsub("<IIDS_TO_KEEP>", keep_file, curr_input_step1_template)

        ### Extract SNPs also based on population

        extract_snps = paste0(currPop, "_qc_snps_for_regenie.snplist")

        curr_input_step1_template = gsub("<SNPLIST_FILE>", extract_snps, curr_input_step1_template)

        ### Covariate file also based on population
        ### This would need to be adjusted if re-running other stratifications

        covar_file = paste0("gwas_cost_covariates_population_", currPop, ".tsv")

        curr_input_step1_template = gsub("<COVARIATE_FILE>", covar_file, curr_input_step1_template)
        curr_input_step2_template = gsub("<COVARIATE_FILE>", covar_file, curr_input_step2_template)

        ### Writing

        writeLines(curr_input_step1_template, paste0(step1Out, "s01_regenie_step1_", analysis_prefix, ".sh"))
        writeLines(curr_input_step2_template, paste0(step2Out, "s02_regenie_step2_", analysis_prefix, ".sh"))


    }
}


