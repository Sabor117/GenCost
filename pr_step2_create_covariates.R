#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library("lubridate")

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("She turned me into a newt!")


#####################

### File locations and script prep

#####################

ukb_base_data_file = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.tsv"
ukb_additional_data_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/ukb_data/ukb671404_cost_phenos.csv"

ukb_translation_file = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_translations.tsv"

death_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/death/death.txt"

ancestry_file = "/scratch/project_2007428/projects/prj_002_breast_cancer_prs/PCSsAncestry/Output/RF_ntree100_nPCS10_nodesize_5_prob_AFR_0.45,AMR_0.42,CSA_0.3,EAS_0.3,EUR_0.3,MID_0.25,OCE_0.3/ukb22828_ancestry.txt"

ukb_translations = fread(ukb_translation_file, data.table = FALSE)
death_reads = fread(death_file, data.table = FALSE)
ukb_ancestry = fread(ancestry_file, data.table = FALSE)

#####################

### Get age at end of follow-up

#####################

### Start date defined by Padraig

start_date = as.Date("01/06/2006", format = "%d/%m/%Y")

### End date defined as last date in death register

end_date = as.Date("19/12/2022", format = "%d/%m/%Y")

### All individuals

age_calc = fread(ukb_additional_data_file, data.table = FALSE, select = c("eid", "34-0.0", "52-0.0"))

colnames(age_calc)[c(1:3)] = c("FID", "year_of_birth", "month_of_birth")

age_calc$dob = as.Date(paste0(1, "/", age_calc$month_of_birth, "/", age_calc$year_of_birth), format = "%d/%m/%Y")

death_reads$date_of_death = as.Date(death_reads$date_of_death, format = "%d/%m/%Y")

age_calc$final_date = death_reads$date_of_death[match(age_calc[["FID"]], death_reads$FID)]

age_calc$final_date = ifelse(is.na(age_calc$final_date) == TRUE, as.Date(end_date), age_calc$final_date)

### This next line is needed for some reason

age_calc$final_date = as.Date(age_calc$final_date, origin = "1970-01-01")

age_calc$age_at_end_of_followup = floor(as.numeric(age_calc$final_date - age_calc$dob) / 365.25)

age_calc = age_calc[complete.cases(age_calc),]


#####################

### Other phenotypes

#####################

pc_fields = paste0("22009-0.", c(1:20))

covariates = fread(ukb_base_data_file, data.table = FALSE, select = c("eid", "31-0.0", "22000-0.0", pc_fields))

colnames(covariates)[c(1:ncol(covariates))] = c("FID", "sex", "batch", paste0("pc", c(1:20)))

covariates$array = ifelse(covariates$batch < 0, 1, ifelse(covariates$batch > 0, 2, NA))

covariates = covariates[complete.cases(covariates),]

covariates = merge(covariates, age_calc[,c("FID", "year_of_birth", "age_at_end_of_followup")])

covariates$age_end_squared = covariates$age_at_end_of_followup ^ 2

covariates$age_end_times_sex = covariates$age_at_end_of_followup * covariates$sex

covariates$IID = covariates$FID

covariates = covariates[,c(1, ncol(covariates), 2:(ncol(covariates)-1))]


#####################

### Getting super-population

#####################

covariates = merge(covariates, ukb_ancestry[,c("IID", "SuperPop")], by = "IID", all = TRUE)

covariates$SuperPop[is.na(covariates$SuperPop)] = "NA"

covariates = covariates[!(is.na(covariates$sex)),]

covariates = covariates[,c(2, 1, 3:ncol(covariates))]

### This line converts the SuperPop types into a number
### As the original file for SuperPops was purged, the values here were reverse-engineered based on superpop counts
### 1 = EUR
### 2 = CSA
### 3 = OTHER
### 4 = AFR
### 5 = EAS
### 6 = NA
### 7 = AMR
### 8 = MID

unique_pops = data.frame(population = unique(covariates$SuperPop),
                            replacement = c(1:length(unique(covariates$SuperPop)))
                        )

for (i in 1:nrow(unique_pops)){

    currPop = unique_pops$population[i]

    covariates$SuperPop[covariates$SuperPop == currPop] = unique_pops$replacement[unique_pops$population == currPop]

}

covariates$SuperPop = as.numeric(covariates$SuperPop)


#####################

### Saving

#####################

fwrite(covariates, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

#covariates$SuperPop = NULL

covariates_male = covariates[covariates$sex == 1,]
covariates_female = covariates[covariates$sex == 0,]

covariates_male$sex = NULL
covariates_male$age_end_times_sex = NULL
covariates_female$sex = NULL
covariates_female$age_end_times_sex = NULL

fwrite(covariates_male, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_male.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

fwrite(covariates_male[,c("IID", "FID")], "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_male_iids_only.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", col.names = FALSE)

fwrite(covariates_female, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_female.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

fwrite(covariates_female[,c("IID", "FID")], "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_female_iids_only.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", col.names = FALSE)

covariates_18 = covariates[covariates$age_at_end_of_followup <= 18,] # NONE in UKB
covariates_35 = covariates[covariates$age_at_end_of_followup <= 35 & covariates$age_at_end_of_followup > 18,] # NONE in UKB
covariates_55 = covariates[covariates$age_at_end_of_followup <= 55 & covariates$age_at_end_of_followup > 35,]
covariates_75 = covariates[covariates$age_at_end_of_followup <= 75 & covariates$age_at_end_of_followup > 55,]
covariates_100 = covariates[covariates$age_at_end_of_followup > 75,]

fwrite(covariates_55, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_36_55.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

fwrite(covariates_75, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_56_75.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

fwrite(covariates_100, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_76_plus.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

fwrite(covariates_55[,c("IID", "FID")], "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_36_to_55_iids_only.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", col.names = FALSE)

fwrite(covariates_75[,c("IID", "FID")], "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_56_to_75_iids_only.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", col.names = FALSE)

fwrite(covariates_100[,c("IID", "FID")], "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/ukb_76_plus_iids_only.tsv",
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA", col.names = FALSE)

for (i in 1:nrow(unique_pops)){

    currPop = as.numeric(unique_pops$replacement[i])

    currPop_name = unique_pops$population[i]

    curr_covariates = covariates[covariates$SuperPop == currPop,]

    curr_covariates$SuperPop = NULL

    fwrite(curr_covariates, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_population_", currPop_name, ".tsv"),
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

}

for (i in 1:nrow(unique_pops)){

    currPop = as.numeric(unique_pops$replacement[i])

    currPop_name = unique_pops$population[i]

    curr_covariates = covariates[covariates$SuperPop == currPop,]

    curr_covariates$SuperPop = NULL
    curr_covariates$sex = NULL
    curr_covariates$age_end_times_sex = NULL

    fwrite(curr_covariates, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cost_phenotypes/gwas_cost_covariates_population_", currPop_name, "_sex.tsv"),
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

}


