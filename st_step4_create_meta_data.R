#####################

### Script set-up

#####################

library(data.table)
library(stringr)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
phenoDir = paste0(mainDir, "outputs/cost_phenotypes/")
outDir = paste0(mainDir, "processing/ukb_gwas_meta_data/")

### Primary files

PRIM_frame_file = paste0(phenoDir, "FINAL_gp_based_phenotype_frame.tsv")
IN_frame_file = paste0(phenoDir, "FINAL_hes_based_phenotype_frame.tsv")
DRUG_frame_file = paste0(phenoDir, "FINAL_prescription_based_phenotype_frame.tsv")

### IID files

iids_36_55_file = paste0(phenoDir, "gwas_cost_covariates_36_55.tsv")
iids_56_75_file = paste0(phenoDir, "gwas_cost_covariates_56_75.tsv")
iids_76_plus_file = paste0(phenoDir, "gwas_cost_covariates_76_plus.tsv")
iids_M_file = paste0(phenoDir, "gwas_cost_covariates_male.tsv")
iids_F_file = paste0(phenoDir, "gwas_cost_covariates_female.tsv")

### Population IIDs

populations = c("EUR", "AFR", "CSA", "EAS")

covars_file = paste0(phenoDir, "gwas_cost_covariates_population_%%%.tsv")

### Reading IIDs

iids_36_55 = fread(iids_36_55_file, data.table = FALSE)
iids_56_75 = fread(iids_56_75_file, data.table = FALSE)
iids_76_plus = fread(iids_76_plus_file, data.table = FALSE)
iids_M = fread(iids_M_file, data.table = FALSE)
iids_F = fread(iids_F_file, data.table = FALSE)

### Stratum set up

strata = c("36_55", "56_75", "76_plus", "males", "females")
covars = c("sex", "age_at_end_of_followup")

### Values function

meta_values = function(phenotype_data, colselect){

    ### Gets mean, median, SD and percentage female based on input
    ### Outputs a row

    out_mean = mean(phenotype_data[,colselect])
    out_med = median(phenotype_data[,colselect])
    out_sd = sd(phenotype_data[,colselect])
    out_percent_F = nrow(phenotype_data[phenotype_data$sex == 0,]) / nrow(phenotype_data) * 100
    out_mad = mad(phenotype_data[,colselect])
    out_iqr = IQR(phenotype_data[,colselect])

    outrow = data.frame(mean = out_mean,
                        median = out_med,
                        sd = out_sd,
                        percent_female = out_percent_F,
                        mad = out_mad,
                        iqr = out_iqr)

    return(c(outrow))

}

### Overarching function

get_meta_data = function(phenotype_data, covars){

    ### Set up output

    pheno_outFrame = data.frame(matrix(ncol = 9, nrow = 0))
    colnames(pheno_outFrame) = c("strata", "phenotype", "mean", "median", "sd", "percent_female", "mad", "iqr", "n")

    age_outFrame = data.frame(matrix(ncol = 9, nrow = 0))
    colnames(age_outFrame) = c("strata", "phenotype", "mean", "median", "sd", "percent_female", "mad", "iqr", "n")

    for (currCol in 4:ncol(phenotype_data)){

        currPheno_name = colnames(phenotype_data)[currCol]

        currPheno = phenotype_data[,c("eid", c(covars), currPheno_name)]

        over_meta_values = meta_values(currPheno, currPheno_name)
        over_meta_age_values = meta_values(currPheno, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "All",
                                            phenotype = currPheno_name,
                                            mean = over_meta_values$mean,
                                            median = over_meta_values$median,
                                            sd = over_meta_values$sd,
                                            percent_female = over_meta_values$percent_female,
                                            mad = over_meta_values$mad,
                                            iqr = over_meta_values$iqr,
                                            n = nrow(currPheno))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "All",
                                            phenotype = currPheno_name,
                                            mean = over_meta_age_values$mean,
                                            median = over_meta_age_values$median,
                                            sd = over_meta_age_values$sd,
                                            percent_female = over_meta_age_values$percent_female,
                                            mad = over_meta_age_values$mad,
                                            iqr = over_meta_age_values$iqr,
                                            n = nrow(currPheno))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)

        currPheno_36_55 = currPheno[currPheno$eid %in% iids_36_55$IID,]
        currPheno_56_75 = currPheno[currPheno$eid %in% iids_56_75$IID,]
        currPheno_76_plus = currPheno[currPheno$eid %in% iids_76_plus$IID,]
        currPheno_M = currPheno[currPheno$eid %in% iids_M$IID,]
        currPheno_F = currPheno[currPheno$eid %in% iids_F$IID,]

        ###

        meta_values_36_55 = meta_values(currPheno_36_55, currPheno_name)
        meta_age_values_36_55 = meta_values(currPheno_36_55, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "36_55",
                                            phenotype = currPheno_name,
                                            mean = meta_values_36_55$mean,
                                            median = meta_values_36_55$median,
                                            sd = meta_values_36_55$sd,
                                            percent_female = meta_values_36_55$percent_female,
                                            mad = meta_values_36_55$mad,
                                            iqr = meta_values_36_55$iqr,
                                            n = nrow(currPheno_36_55))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "36_55",
                                            phenotype = currPheno_name,
                                            mean = meta_age_values_36_55$mean,
                                            median = meta_age_values_36_55$median,
                                            sd = meta_age_values_36_55$sd,
                                            percent_female = meta_age_values_36_55$percent_female,
                                            mad = meta_age_values_36_55$mad,
                                            iqr = meta_age_values_36_55$iqr,
                                            n = nrow(currPheno_36_55))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)

        ###

        meta_values_56_75 = meta_values(currPheno_56_75, currPheno_name)
        meta_age_values_56_75 = meta_values(currPheno_56_75, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "56_75",
                                            phenotype = currPheno_name,
                                            mean = meta_values_56_75$mean,
                                            median = meta_values_56_75$median,
                                            sd = meta_values_56_75$sd,
                                            percent_female = meta_values_56_75$percent_female,
                                            mad = meta_values_56_75$mad,
                                            iqr = meta_values_56_75$iqr,
                                            n = nrow(currPheno_56_75))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "56_75",
                                            phenotype = currPheno_name,
                                            mean = meta_age_values_56_75$mean,
                                            median = meta_age_values_56_75$median,
                                            sd = meta_age_values_56_75$sd,
                                            percent_female = meta_age_values_56_75$percent_female,
                                            mad = meta_age_values_56_75$mad,
                                            iqr = meta_age_values_56_75$iqr,
                                            n = nrow(currPheno_56_75))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)

        ###

        meta_values_76_plus = meta_values(currPheno_76_plus, currPheno_name)
        meta_age_values_76_plus = meta_values(currPheno_76_plus, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "76_plus",
                                            phenotype = currPheno_name,
                                            mean = meta_values_76_plus$mean,
                                            median = meta_values_76_plus$median,
                                            sd = meta_values_76_plus$sd,
                                            percent_female = meta_values_76_plus$percent_female,
                                            mad = meta_values_76_plus$mad,
                                            iqr = meta_values_76_plus$iqr,
                                            n = nrow(currPheno_76_plus))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "76_plus",
                                            phenotype = currPheno_name,
                                            mean = meta_age_values_76_plus$mean,
                                            median = meta_age_values_76_plus$median,
                                            sd = meta_age_values_76_plus$sd,
                                            percent_female = meta_age_values_76_plus$percent_female,
                                            mad = meta_age_values_76_plus$mad,
                                            iqr = meta_age_values_76_plus$iqr,
                                            n = nrow(currPheno_76_plus))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)

        ###

        meta_values_M = meta_values(currPheno_M, currPheno_name)
        meta_age_values_M = meta_values(currPheno_M, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "males",
                                            phenotype = currPheno_name,
                                            mean = meta_values_M$mean,
                                            median = meta_values_M$median,
                                            sd = meta_values_M$sd,
                                            percent_female = meta_values_M$percent_female,
                                            mad = meta_values_M$mad,
                                            iqr = meta_values_M$iqr,
                                            n = nrow(currPheno_M))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "males",
                                            phenotype = currPheno_name,
                                            mean = meta_age_values_M$mean,
                                            median = meta_age_values_M$median,
                                            sd = meta_age_values_M$sd,
                                            percent_female = meta_age_values_M$percent_female,
                                            mad = meta_age_values_M$mad,
                                            iqr = meta_age_values_M$iqr,
                                            n = nrow(currPheno_M))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)

        ###

        meta_values_F = meta_values(currPheno_F, currPheno_name)
        meta_age_values_F = meta_values(currPheno_F, "age_at_end_of_followup")

        pheno_outFrame_row = data.frame(strata = "females",
                                            phenotype = currPheno_name,
                                            mean = meta_values_F$mean,
                                            median = meta_values_F$median,
                                            sd = meta_values_F$sd,
                                            percent_female = meta_values_F$percent_female,
                                            mad = meta_values_F$mad,
                                            iqr = meta_values_F$iqr,
                                            n = nrow(currPheno_F))

        pheno_outFrame = rbind(pheno_outFrame, pheno_outFrame_row)

        age_outFrame_row = data.frame(strata = "females",
                                            phenotype = currPheno_name,
                                            mean = meta_age_values_F$mean,
                                            median = meta_age_values_F$median,
                                            sd = meta_age_values_F$sd,
                                            percent_female = meta_age_values_F$percent_female,
                                            mad = meta_age_values_F$mad,
                                            iqr = meta_age_values_F$iqr,
                                            n = nrow(currPheno_F))

        age_outFrame = rbind(age_outFrame, age_outFrame_row)
    }

    return(c(pheno_outFrame, age_outFrame))

}


#####################

### Script run

#####################

### Read main three files

PRIM_cols = c("eid", "gp_costs_year", "total_costs_year")
PRIM_frame = fread(PRIM_frame_file, data.table = FALSE, select = PRIM_cols)
colnames(PRIM_frame)[2] = "Prim"

IN_cols = c("eid", "hes_costs_year")
IN_frame = fread(IN_frame_file, data.table = FALSE, select = IN_cols)
colnames(IN_frame)[2] = "Inpatient"

DRUG_cols = c("eid", "prescription_costs_year")
DRUG_frame = fread(DRUG_frame_file, data.table = FALSE, select = DRUG_cols)
colnames(DRUG_frame)[2] = "Drug"

### Set up output tables

IN_meta_data = data.frame(matrix(ncol = 10, nrow = 0))
IN_meta_age_data = data.frame(matrix(ncol = 10, nrow = 0))
DRUG_meta_data = data.frame(matrix(ncol = 10, nrow = 0))
DRUG_meta_age_data = data.frame(matrix(ncol = 10, nrow = 0))
PRIM_meta_data = data.frame(matrix(ncol = 10, nrow = 0))
PRIM_meta_age_data = data.frame(matrix(ncol = 10, nrow = 0))

output_cols = c("cohort", "phenotype", "strata", "percent_female", "mean", "sd", "median", "mad", "iqr", "n")

colnames(IN_meta_data) = output_cols
colnames(IN_meta_age_data) = output_cols
colnames(DRUG_meta_data) = output_cols
colnames(DRUG_meta_age_data) = output_cols
colnames(PRIM_meta_data) = output_cols
colnames(PRIM_meta_age_data) = output_cols

### Run loop for each population

for (i in 1:length(populations)){

    ### Select population

    curr_pop = populations[i]

    curr_covar_file = gsub("%%%", curr_pop, covars_file)

    ### Read covariates for population

    total_covars = fread(curr_covar_file, data.table = FALSE)

    ### Merge with cost frames

    curr_PRIM_frame = merge(PRIM_frame, total_covars[,c("IID", "sex", "age_at_end_of_followup")], by.x = "eid", by.y = "IID")
    curr_PRIM_frame = curr_PRIM_frame[,c("eid", covars, "Prim", "total_costs_year")]

    curr_IN_frame = merge(IN_frame, total_covars[,c("IID", "sex", "age_at_end_of_followup")], by.x = "eid", by.y = "IID")
    curr_IN_frame = curr_IN_frame[,c("eid", covars, "Inpatient")]

    curr_DRUG_frame = merge(DRUG_frame, total_covars[,c("IID", "sex", "age_at_end_of_followup")], by.x = "eid", by.y = "IID")
    curr_DRUG_frame = curr_DRUG_frame[,c("eid", covars, "Drug")]

    ### Creating medians, means etc
    ### For inpatient

    IN_meta_stats = get_meta_data(curr_IN_frame, covars)

    IN_meta_cost_frame = data.frame(cohort = paste0("UKB_", curr_pop),
                                    phenotype = IN_meta_stats[2],
                                    strata = IN_meta_stats[1],
                                    percent_female = IN_meta_stats[6],
                                    mean = IN_meta_stats[3],
                                    sd = IN_meta_stats[5],
                                    median = IN_meta_stats[4],
                                    mad = IN_meta_stats[7],
                                    iqr = IN_meta_stats[8],
                                    n = IN_meta_stats[9])

    IN_meta_age_frame = data.frame(cohort = paste0("UKB_", curr_pop),
                                    phenotype = IN_meta_stats[11], 
                                    strata = IN_meta_stats[10],
                                    percent_female = IN_meta_stats[15],
                                    mean = IN_meta_stats[12],
                                    sd = IN_meta_stats[14],
                                    median = IN_meta_stats[13],
                                    mad = IN_meta_stats[16],
                                    iqr = IN_meta_stats[17],
                                    n = IN_meta_stats[18])

    ### Merge into output

    IN_meta_data = rbind(IN_meta_data, IN_meta_cost_frame)
    IN_meta_age_data = rbind(IN_meta_age_data, IN_meta_age_frame)

    ### For DRUG

    DRUG_meta_stats = get_meta_data(curr_DRUG_frame, covars)

    DRUG_meta_cost_frame = data.frame(cohort = paste0("UKB_", curr_pop),
                                        phenotype = DRUG_meta_stats[2],
                                        strata = DRUG_meta_stats[1],
                                        percent_female = DRUG_meta_stats[6],
                                        mean = DRUG_meta_stats[3],
                                        sd = DRUG_meta_stats[5],
                                        median = DRUG_meta_stats[4],
                                        mad = DRUG_meta_stats[7],
                                        iqr = DRUG_meta_stats[8],
                                        n = DRUG_meta_stats[9])


    DRUG_meta_age_frame = data.frame(cohort = paste0("UKB_", curr_pop),
                                        phenotype = DRUG_meta_stats[11], 
                                        strata = DRUG_meta_stats[10],
                                        percent_female = DRUG_meta_stats[15],
                                        mean = DRUG_meta_stats[12],
                                        sd = DRUG_meta_stats[14],
                                        median = DRUG_meta_stats[13],
                                        mad = DRUG_meta_stats[16],
                                        iqr = DRUG_meta_stats[17],
                                        n = DRUG_meta_stats[18])

    ### Merge into output

    DRUG_meta_data = rbind(DRUG_meta_data, DRUG_meta_cost_frame)
    DRUG_meta_age_data = rbind(DRUG_meta_age_data, DRUG_meta_age_frame)

    ### For PRIM

    PRIM_meta_stats = get_meta_data(curr_PRIM_frame, covars)

    PRIM_meta_cost_frame = data.frame(cohort = paste0("UKB_", curr_pop),
                                        phenotype = PRIM_meta_stats[2],
                                        strata = PRIM_meta_stats[1],
                                        percent_female = PRIM_meta_stats[6],
                                        mean = PRIM_meta_stats[3],
                                        sd = PRIM_meta_stats[5],
                                        median = PRIM_meta_stats[4],
                                        mad = PRIM_meta_stats[7],
                                        iqr = PRIM_meta_stats[8],
                                        n = PRIM_meta_stats[9])


    PRIM_meta_age_frame = data.frame(ccohort = paste0("UKB_", curr_pop),
                                        phenotype = PRIM_meta_stats[11], 
                                        strata = PRIM_meta_stats[10],
                                        percent_female = PRIM_meta_stats[15],
                                        mean = PRIM_meta_stats[12],
                                        sd = PRIM_meta_stats[14],
                                        median = PRIM_meta_stats[13],
                                        mad = PRIM_meta_stats[16],
                                        iqr = PRIM_meta_stats[17],
                                        n = PRIM_meta_stats[18])

    ### Merge into output

    PRIM_meta_data = rbind(PRIM_meta_data, PRIM_meta_cost_frame)
    PRIM_meta_age_data = rbind(PRIM_meta_age_data, PRIM_meta_age_frame)

}

### Writing files

fwrite(IN_meta_data, paste0(outDir, "IN_meta_cost_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(IN_meta_age_data, paste0(outDir, "IN_meta_age_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(DRUG_meta_data, paste0(outDir, "DRUG_meta_cost_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(DRUG_meta_age_data, paste0(outDir, "DRUG_meta_age_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(PRIM_meta_data, paste0(outDir, "PRIM_meta_cost_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(PRIM_meta_age_data, paste0(outDir, "PRIM_meta_age_frame.txt"), row.names = FALSE, quote = FALSE, sep = "\t")