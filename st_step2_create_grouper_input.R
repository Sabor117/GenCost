#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library(Hmisc)
library(dplyr)
library(stringr)
library(tidyr)

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("She turned me into a newt!")


#####################

### File locations and script prep

#####################

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"

ukb_base_cost_phenotype_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/ukb_data/ukb671404_cost_phenos.csv" # Extracted from UKB encoded file
ukb_hesin_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin.txt" # Direct UKB download
ukb_hesin_diagnoses_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin_diag.txt" # Direct UKB download
ukb_translation_file = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_translations.tsv" #Extracted from UKB encoded file
ukb_hesin_operations_file =  "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin_oper.txt" # Direct UKB download
ukb_crit_days_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin_critical.txt" # Direct UKB download

ukb_base_cost_phenotypes = fread(ukb_base_cost_phenotype_file, data.table = FALSE)
ukb_hesin = fread(ukb_hesin_file, data.table = FALSE)
ukb_hesin_diag = fread(ukb_hesin_diagnoses_file, data.table = FALSE)
ukb_translations = fread(ukb_translation_file, data.table = FALSE)
ukb_hesin_oper = fread(ukb_hesin_operations_file, data.table = FALSE)
ukb_crit_days = fread(ukb_crit_days_file, data.table = FALSE)


#####################

### Start data parsing

#####################

#####################

### Data formatting for base data

#####################

### Removing multiples of attending assessment centre, death and location of centre from columns

keep_cols = str_split_fixed(colnames(ukb_base_cost_phenotypes), "\\.", n = 2)
keep_cols = str_split_fixed(keep_cols[,1], "-", n = 2)
keep_cols = colnames(ukb_base_cost_phenotypes)[keep_cols[,2] == "" | keep_cols[,2] == 0]

ukb_base_cost_phenotypes = ukb_base_cost_phenotypes[,keep_cols]

field_names = c("eid",
                "sex",
                "birth_year",
                "birth_month",
                "date_assessed",
                "assessment_centre",
                "age_recruited",
                "genetic_sex",
                "death_age",
                "inpatient_record_0",
                "inpatient_record_1",
                "inpatient_record_2",
                "inpatient_record_3",
                "inpatient_record_4",
                "inpatient_record_5"
                )

colnames(ukb_base_cost_phenotypes) = field_names

### Add DOB column
### No birth day, so assume 01 to estimate

ukb_base_cost_phenotypes$dob = as.Date(paste0("01/", ukb_base_cost_phenotypes$birth_month, "/", ukb_base_cost_phenotypes$birth_year), format = "%d/%m/%Y")


#####################

### Data formatting for HESIN data

#####################

### Remove data from prior to 2006
### Remove data with NA in ID column

ukb_hesin$epistart = as.Date(ukb_hesin$epistart, format = "%d/%m/%Y")
ukb_hesin$admidate = as.Date(ukb_hesin$admidate, format = "%d/%m/%Y")

ukb_hesin = ukb_hesin[ukb_hesin$epistart > "2006-06-01",]
ukb_hesin = ukb_hesin[ukb_hesin$admidate > "2006-06-01",]
ukb_hesin = ukb_hesin[!(is.na(ukb_hesin$eid)),]

### Add sex and DOB to HESIN data

ukb_hesin = merge(ukb_hesin, ukb_base_cost_phenotypes[,c("eid", "sex", "dob")], by = "eid")

### Calculate start age based on estimated DOB and admittance date

ukb_hesin$startage = floor(as.numeric(difftime(ukb_hesin$admidate, ukb_hesin$dob, unit = "weeks")) / 52.25)

### Convert sex from UKB format (1/0) to NHS format (1/2)
### Men number 1

ukb_hesin$sex = ifelse(ukb_hesin$sex == 0, 2, ifelse(ukb_hesin$sex == 1, 1, NA))


#####################

### Data formatting for ICD10 data

#####################

### Remove HES entries based on UZ05 errors from first Grouper run
### May not be necessary

#diag_tidy_list_0 = fread(paste0(mainDir, "processing/nhs_grouper/diag_tidy_list_2_rems.txt"), data.table = FALSE)

#ukb_hesin_diag = ukb_hesin_diag[!(ukb_hesin_diag$diag_icd10 %in% diag_tidy_list_0$code_to_remove),]

### Order data based first on EID, then instance index and then array index
### Organises diagnoses in order

ukb_hesin_diag = ukb_hesin_diag[with(ukb_hesin_diag, order(ukb_hesin_diag$eid, ukb_hesin_diag$ins_index, ukb_hesin_diag$arr_index)),]

### Remove duplicated rows based on EID + instance

ukb_hesin_diag_duplicates = data.frame(iid = paste0(ukb_hesin_diag$eid, "_", ukb_hesin_diag$ins_index),
                                        diag = ukb_hesin_diag$diag_icd10)

ukb_hesin_diag = ukb_hesin_diag[!(duplicated(ukb_hesin_diag_duplicates)),]                                        

### Convert long data format into wide

ukb_hesin_diag_aggr = ukb_hesin_diag %>%
  group_by(eid, ins_index) %>%
  dplyr::summarize(diag_icd10 = paste(diag_icd10, collapse = ",")) %>%
  ungroup()

ukb_hesin_diag_out = ukb_hesin_diag_aggr %>%
  separate(col = diag_icd10, sep = ",", into = paste0("diag_", seq_len(max(str_count(ukb_hesin_diag_aggr$diag_icd10, ","))+1)), fill = "right")

ukb_hesin_diag_out = as.data.frame(ukb_hesin_diag_out)

ukb_hesin_diag_out = ukb_hesin_diag_out %>%
  rename_with(~ sprintf("diag_%02d", as.numeric(str_extract(., "\\d+"))), starts_with("diag_"))


#####################

### Data formatting for operations data

#####################

### Remove OPCS entries from X998, X999 and X64

oper_tidy_list_0 = c("X998", "X999", "X64")

ukb_hesin_oper = ukb_hesin_oper[!(ukb_hesin_oper$oper4 %in% oper_tidy_list_0),]

### Order data based first on EID, then instance index and then array index
### Organises operations in order

ukb_hesin_oper = ukb_hesin_oper[with(ukb_hesin_oper, order(ukb_hesin_oper$eid, ukb_hesin_oper$ins_index, ukb_hesin_oper$arr_index)),]

### Convert long data format into wide
### Only use oper4 column as this is OPCS4

ukb_hesin_oper_aggr = ukb_hesin_oper %>%
  group_by(eid, ins_index) %>%
  dplyr::summarize(oper = paste(oper4, collapse = ",")) %>%
  ungroup()

ukb_hesin_oper_out = ukb_hesin_oper_aggr %>%
  separate(col = oper, sep = ",", into = paste0("oper_", seq_len(max(str_count(ukb_hesin_oper_aggr$oper, ","))+1)), fill = "right")

ukb_hesin_oper_out = as.data.frame(ukb_hesin_oper_out)

ukb_hesin_oper_out = ukb_hesin_oper_out %>%
  rename_with(~ sprintf("oper_%02d", as.numeric(str_extract(., "\\d+"))), starts_with("oper_"))


#####################

### Data formatting for critical care days

#####################

ukb_crit_days$ccstartdate = as.Date(ukb_crit_days$ccstartdate, format = "%d/%m/%Y")
ukb_crit_days$ccdisdate = as.Date(ukb_crit_days$ccdisdate, format = "%d/%m/%Y")

ukb_crit_days$ccdays = as.numeric(difftime(ukb_crit_days$ccdisdate, ukb_crit_days$ccstartdate, unit = "days"))


#####################

### Main data output

#####################

### Convert EID + INS Index into "IID"

ukb_hesin$iid = paste0(ukb_hesin$eid, "_", ukb_hesin$ins_index)
ukb_hesin_diag_out$iid = paste0(ukb_hesin_diag_out$eid, "_", ukb_hesin_diag_out$ins_index)
ukb_hesin_oper_out$iid = paste0(ukb_hesin_oper_out$eid, "_", ukb_hesin_oper_out$ins_index)
ukb_crit_days$iid = paste0(ukb_crit_days$eid, "_", ukb_crit_days$ins_index)

### Convert critical care days to just days

ukb_crit_days = ukb_crit_days[,c("iid", "ccdays")]

### Some IIDs are duplicated - for those, sum the ccdays (usually instances of individuals being discharged and readmitted on the same day)

ukb_crit_days = ukb_crit_days %>% 
                        group_by(iid) %>% 
                        dplyr::summarize(total_ccdays = sum(ccdays)) %>% 
                        distinct()

ukb_crit_days = as.data.frame(ukb_crit_days)

### For those IIDs without any critical care days, add 0 value                        

ukb_crit_days_miss = data.frame(iid = ukb_hesin$iid[!(ukb_hesin$iid %in% ukb_crit_days$iid)],
                                  total_ccdays = 0)

ukb_crit_days = rbind(ukb_crit_days, ukb_crit_days_miss)

### Remove EID + ins_index columns from diag and oper data for merging

ukb_hesin_diag_out = ukb_hesin_diag_out[,-c(1,2)]
ukb_hesin_oper_out = ukb_hesin_oper_out[,-c(1,2)]

### Add additional NAs for IIDs in HESIN not present in either operations or diagnosis data set

ukb_hesin_diag_out_miss = data.frame(iid = ukb_hesin$iid[!(ukb_hesin$iid %in% ukb_hesin_diag_out$iid)])
ukb_hesin_diag_out_miss[,c(2:ncol(ukb_hesin_diag_out))] = NA
colnames(ukb_hesin_diag_out_miss)[2:ncol(ukb_hesin_diag_out_miss)] = colnames(ukb_hesin_diag_out)[1:(ncol(ukb_hesin_diag_out)-1)]

ukb_hesin_diag_out = rbind(ukb_hesin_diag_out, ukb_hesin_diag_out_miss)

ukb_hesin_oper_out_miss = data.frame(iid = ukb_hesin$iid[!(ukb_hesin$iid %in% ukb_hesin_oper_out$iid)])
ukb_hesin_oper_out_miss[,c(2:ncol(ukb_hesin_oper_out))] = NA
colnames(ukb_hesin_oper_out_miss)[2:ncol(ukb_hesin_oper_out_miss)] = colnames(ukb_hesin_oper_out)[1:(ncol(ukb_hesin_oper_out)-1)]

ukb_hesin_oper_out = rbind(ukb_hesin_oper_out, ukb_hesin_oper_out_miss)

### Merge operations and diagnoses datasets with HESIN

ukb_hesin = merge(ukb_hesin, ukb_hesin_diag_out, by = "iid")
print(nrow(ukb_hesin))
print(length(unique(ukb_hesin$eid)))
ukb_hesin = merge(ukb_hesin, ukb_hesin_oper_out, by = "iid")
print(nrow(ukb_hesin))
print(length(unique(ukb_hesin$eid)))
ukb_hesin = merge(ukb_hesin, ukb_crit_days, by = "iid")
print(nrow(ukb_hesin))
print(length(unique(ukb_hesin$eid)))

### Create grouper output

grouper_output = data.frame(procodet = "ZZZ",
                              provspno = ukb_hesin$iid,
                              epiorder = ukb_hesin$epiorder,
                              startage = ukb_hesin$startage,
                              sex = ukb_hesin$sex,
                              classpat = ukb_hesin$classpat,
                              admisorc = ukb_hesin$admisorc,
                              admimeth = ukb_hesin$admimeth,
                              disdest = ukb_hesin$disdest,
                              dismeth = ukb_hesin$dismeth,
                              epidur = ukb_hesin$epidur,
                              mainspef = ukb_hesin$mainspef,
                              neocare = 8,
                              tretspef = ukb_hesin$tretspef,
                              diag_01 = ukb_hesin$diag_01,
                              diag_02 = ukb_hesin$diag_02,
                              diag_03 = ukb_hesin$diag_03,
                              diag_04 = ukb_hesin$diag_04,
                              diag_05 = ukb_hesin$diag_05,
                              diag_06 = ukb_hesin$diag_06,
                              diag_07 = ukb_hesin$diag_07,
                              diag_08 = ukb_hesin$diag_08,
                              diag_09 = ukb_hesin$diag_09,
                              diag_10 = ukb_hesin$diag_10,
                              diag_11 = ukb_hesin$diag_11,
                              diag_12 = ukb_hesin$diag_12,
                              diag_13 = ukb_hesin$diag_13,
                              diag_14 = ukb_hesin$diag_14,
                              diag_15 = ukb_hesin$diag_15,
                              diag_16 = ukb_hesin$diag_16,
                              diag_17 = ukb_hesin$diag_17,
                              diag_18 = ukb_hesin$diag_18,
                              diag_19 = ukb_hesin$diag_19,
                              diag_20 = ukb_hesin$diag_20,
                              diag_21 = ukb_hesin$diag_21,
                              oper_01 = ukb_hesin$oper_01,
                              oper_02 = ukb_hesin$oper_02,
                              oper_03 = ukb_hesin$oper_03,
                              oper_04 = ukb_hesin$oper_04,
                              oper_05 = ukb_hesin$oper_05,
                              oper_06 = ukb_hesin$oper_06,
                              oper_07 = ukb_hesin$oper_07,
                              oper_08 = ukb_hesin$oper_08,
                              oper_09 = ukb_hesin$oper_09,
                              oper_10 = ukb_hesin$oper_10,
                              oper_11 = ukb_hesin$oper_11,
                              oper_12 = ukb_hesin$oper_12,
                              oper_13 = ukb_hesin$oper_13,
                              oper_14 = ukb_hesin$oper_14,
                              oper_15 = ukb_hesin$oper_15,
                              oper_16 = ukb_hesin$oper_16,
                              oper_17 = ukb_hesin$oper_17,
                              oper_18 = ukb_hesin$oper_18,
                              oper_19 = ukb_hesin$oper_19,
                              oper_20 = ukb_hesin$oper_20,
                              oper_21 = ukb_hesin$oper_21,
                              oper_22 = ukb_hesin$oper_22,
                              oper_23 = ukb_hesin$oper_23,
                              oper_24 = ukb_hesin$oper_24,
                              criticalcaredays = ukb_hesin$total_ccdays,
                              rehabilitationdays = 0,
                              spcdays = 0
                              )

grouper_output = unique(grouper_output)

fwrite(grouper_output, paste0(mainDir, "processing/nhs_grouper/grouper_unfiltered_data.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
grouper_output = fread(paste0(mainDir, "processing/nhs_grouper/grouper_unfiltered_data.tsv"), data.table = FALSE)


#####################

### Begin correcting Grouper names

#####################

grouper_output_filtered = grouper_output

### Read in code fixes

diag_tidy_list_1 = fread(paste0(mainDir, "processing/nhs_grouper/input/diag_tidy_list_1_PD.txt"), data.table = FALSE)
diag_tidy_list_2 = fread(paste0(mainDir, "processing/nhs_grouper/input/diag_tidy_list_addX_PD.txt"), data.table = FALSE)
diag_tidy_list_3 = fread(paste0(mainDir, "processing/nhs_grouper/input/diag_tidy_list_addX_SMW.txt"), data.table = FALSE)
oper_tidy_list = fread(paste0(mainDir, "processing/nhs_grouper/input/oper_tidy_list_1_PD.txt"), data.table = FALSE)

oper_cols = colnames(grouper_output)[(15+21):(length(colnames(grouper_output))-3)]
diag_cols = colnames(grouper_output)[15:(length(colnames(grouper_output))-3-24)]

### Fix HES codes based on Padraig Dixon curated list

for (i in 1:nrow(diag_tidy_list_1)){

  curr_old_code = diag_tidy_list_1$old_code[i]

  curr_new_code = diag_tidy_list_1$new_code[i]

  grouper_output_filtered[, diag_cols][grouper_output_filtered[, diag_cols] == curr_old_code] = curr_new_code

}

### Add "X" to end of HES codes based on Padraig Dixon + SMW curated list

diag_tidy_list_2 = rbind(diag_tidy_list_2, diag_tidy_list_3)

for (i in 1:nrow(diag_tidy_list_2)){

  curr_old_code = diag_tidy_list_2$old_code[i]

  curr_new_code = paste0(curr_old_code, "X")

  grouper_output_filtered[, diag_cols][grouper_output_filtered[, diag_cols] == curr_old_code] = curr_new_code

}

### Fix OPCS4 codes based on Padraig Dixon curated list

for (i in 1:nrow(oper_tidy_list)){

  curr_old_code = oper_tidy_list$old_code[i]

  curr_new_code = oper_tidy_list$new_code[i]

  grouper_output_filtered[, oper_cols][grouper_output_filtered[, oper_cols] == curr_old_code] = curr_new_code

}

### Remove any IIDs which have all NAs for operations + HESIN

colchecks = c(diag_cols, oper_cols)

grouper_output_filtered = grouper_output_filtered[!(apply(is.na(grouper_output_filtered[,colchecks]), 1, all)),]

### Remove any IIDs which have all NAs for only HESIN (operations but no diagnoses)

grouper_output_filtered = grouper_output_filtered[!(apply(is.na(grouper_output_filtered[,diag_cols]), 1, all)),]

### Remove any IIDs which have all "" for only HESIN (operations but no diagnoses)

grouper_output_filtered = grouper_output_filtered[!(apply(grouper_output_filtered[,diag_cols] == "", 1, all)),]

### Fix column names to capitals

colnames(grouper_output_filtered) = toupper(colnames(grouper_output_filtered))

#####################

### Final output

#####################

fwrite(grouper_output_filtered, paste0(mainDir, "processing/nhs_grouper/grouper_filtered_data.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

fwrite(grouper_output_filtered, paste0(mainDir, "processing/nhs_grouper/grouper_filtered_data.csv"), row.names = FALSE, sep = ",")

end_time = Sys.time()

time_taken = end_time - start_time

heading(paste0("Total time taken is: ", time_taken))


heading("I got better...")

