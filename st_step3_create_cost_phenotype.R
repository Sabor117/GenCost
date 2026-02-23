#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Come and see the payments inherent in the system. Help! Help! We're being repressed!!")


#####################

### File locations and script prep

#####################

### Overarching directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
dataDir = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/"

### Specific files

### UKB phenptypes

ukb_phenos_file = paste0(mainDir, "processing/ukb_data/ukb671404_cost_phenos.csv")

### Grouper files

grouper_output_file = paste0(mainDir, "processing/nhs_grouper/output/grouper_filtered_output_FCE.csv")
grouper_classpat_recode_list_file = paste0(mainDir, "processing/misc_grouper_processing/grouper_classpat_recoding.txt")

### Care cost files (NHS national schedule) and prescription costs
### 2017 costs only used for excess days - converted to 2022 in this script
### Prescription costs were manually extracted from the drug tariff for 2022

nhs_care_costs_file_2022 = paste0(mainDir, "processing/cost_coding_values/2_National_schedule_of_NHS_costs_FY21-22_v2.xlsx")
nhs_care_costs_file_2017 = paste0(mainDir, "processing/cost_coding_values/2_-_National_schedule_of_reference_costs_-_the_main_schedule.xlsx")
prescriptions_costs_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_prescriptions_sheet.csv")

### Prescription translations files

prescriptions_AMPVMP_translations_file_1 = paste0(mainDir, "processing/cost_coding_values/prescription_costing/AMP_to_VMP_translation_1.xlsx")
prescriptions_AMPVMP_translations_file_2 = paste0(mainDir, "processing/cost_coding_values/prescription_costing/AMP_to_VMP_translation_2.csv")
prescriptions_AMPVMP_translations_file_3 = paste0(mainDir, "processing/cost_coding_values/prescription_costing/AMP_to_VMP_translation_3.csv")
prescriptions_incremental_costs_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_price_per_cost_prescriptions.csv")
prescriptions_curated_translations_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_prescriptions_translations_curated.csv")
prescriptions_appliances_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_appliances_costs.csv")
outdated_prescriptions_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_prescriptions_outdated_costs.csv")

### Records files

hes_data_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin.txt"
death_file = paste0(dataDir, "healthcare_phenotypes/death/death.txt")
prescriptions_file = paste0(dataDir, "healthcare_phenotypes/gp_data/gp_scripts.txt.gz")
gp_visits_file = paste0(dataDir, "healthcare_phenotypes/gp_data/gp_clinical.txt.gz")

### File reading

grouper_output = readLines(grouper_output_file)
grouper_classpat_recode_list = fread(grouper_classpat_recode_list_file, data.table = FALSE)
care_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "APC")
consultant_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "OP")
unbun_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "Total HRGs")
elective_excess_costs_2017 = read_excel(nhs_care_costs_file_2017, sheet = "EL_XS")
non_elective_excess_costs_2017 = read_excel(nhs_care_costs_file_2017, sheet = "NEL_XS")
death_reads = fread(death_file, data.table = FALSE)
gp_prescriptions = fread(prescriptions_file, data.table = FALSE)
prescriptions_costs = fread(prescriptions_costs_file, data.table = FALSE)
hes_data = fread(hes_data_file, data.table = FALSE, select = c("eid", "ins_index", "epistat", "epistart", "epiend"))
gp_visits = fread(gp_visits_file, data.table = FALSE)
prescriptions_AMPVMP_translations_1 = read_excel(prescriptions_AMPVMP_translations_file_1, sheet = "AMP to VMP")
prescriptions_exclude_1 = read_excel(prescriptions_AMPVMP_translations_file_1, sheet = "Never prescribable")
prescriptions_exclude_2 = read_excel(prescriptions_AMPVMP_translations_file_1, sheet = "No idea")
prescriptions_appliances_1 = as.data.frame(read_excel(prescriptions_AMPVMP_translations_file_1, sheet = "Medical devices"))
prescriptions_appliances_2 = fread(prescriptions_appliances_file, data.table = FALSE)
prescriptions_AMPVMP_translations_2 = fread(prescriptions_AMPVMP_translations_file_2, data.table = FALSE)
prescriptions_AMPVMP_translations_3 = fread(prescriptions_AMPVMP_translations_file_3, data.table = FALSE)
prescriptions_incremental_costs = fread(prescriptions_incremental_costs_file, data.table = FALSE)
prescriptions_curated_translations = fread(prescriptions_curated_translations_file, data.table = FALSE)
outdated_prescriptions = fread(outdated_prescriptions_file, data.table = FALSE)
ukb_phenos = fread(ukb_phenos_file, data.table = FALSE, select = c("eid", "53-0.0"))

### Global variables

cpi_index_2017_to_2022 = 1.18 # Value extracted from: https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator

### Start date defined by Padraig

start_date = as.Date("01/06/2006", format = "%d/%m/%Y")

### End date defined as last date in death register

end_date = as.Date("19/12/2022", format = "%d/%m/%Y")

### Get "initial date" for HES follow-up

colnames(ukb_phenos)[2] = "date_of_entry"

#####################
#####################

### COST CODE parsing

#####################
#####################

### Care costs read from National Schecule Excel sheet (one sheet of several)
### Sheet is converted into data frame

care_costs = as.data.frame(care_costs_2022)

### Superfluous rows removed

care_costs = care_costs[-c(1:4),]

### Columns renamed correctly

colnames(care_costs) = c("department_code",
                            "department_desc",
                            "currency_code",
                            "currency_desc",
                            "n_fce",
                            "national_avg_unit_cost",
                            "total_cost",
                            "n_submissions",
                            "apc_sector")

### Ensuring columns are numeric correctly                            

care_costs$n_fce = as.numeric(care_costs$n_fce)
care_costs$n_fce[is.na(care_costs$n_fce)] = 0

care_costs$national_avg_unit_cost = as.numeric(care_costs$national_avg_unit_cost)
care_costs$national_avg_unit_cost[is.na(care_costs$national_avg_unit_cost)] = 0

care_costs$n_submissions = as.numeric(care_costs$n_submissions)

### Some negative values (for some reason)

care_costs$national_avg_unit_cost[care_costs$national_avg_unit_cost < 0] = care_costs$national_avg_unit_cost[care_costs$national_avg_unit_cost < 0] * -1

heading("Care costs parsed and ready.")

print(head(care_costs, n = 10))


#####################

### Calculating costs for NEL/EL/NES units

#####################

### Extracting EL, NES and NEL coding "Elective Inpatients", "Non-Elective Inpatient - Short Stay" and "Non-Elective Inpatient - Long Stay"

care_costs_el = care_costs[care_costs$department_code == "EL", c(1,3,5:9)]
care_costs_nes = care_costs[care_costs$department_code == "NES", c(1,3,5:9)]
care_costs_nel = care_costs[care_costs$department_code == "NEL", c(1,3,5:9)]

### Re-naming columns for merging

colnames(care_costs_el) = paste0(colnames(care_costs_el), "_el")
colnames(care_costs_nes) = paste0(colnames(care_costs_nes), "_nes")
colnames(care_costs_nel) = paste0(colnames(care_costs_nel), "_nel")

colnames(care_costs_el)[2] = "currency_code"
colnames(care_costs_nes)[2] = "currency_code"
colnames(care_costs_nel)[2] = "currency_code"

care_costs_el[,1] = NULL
care_costs_nes[,1] = NULL
care_costs_nel[,1] = NULL

### Creating additional weighted column by adding community costs
### DONE FOR EACH SET

weighted_avg_maker = function(cost_table,
                                currency_code_field,
                                unit_cost_field,
                                num_visits_field,
                                healthcare_area){

    ### Place table in order so Acute visit is always first

    cost_table = cost_table[order(cost_table[,healthcare_area]),]

    ### Creating new column for weighted average of Acute + Community care costs
    ### Also create summed N

    cost_table$weighted_avg_all_unit_cost = cost_table[,unit_cost_field]
    cost_table$total_fce = cost_table[,num_visits_field]

    ### Identify duplicate currency_code values

    duplicate_codes = cost_table[duplicated(cost_table[,currency_code_field]), currency_code_field]

    ### Calculate weighted average for duplicates

    for (code in duplicate_codes) {
        
        subset_cost_data = cost_table[cost_table[,currency_code_field] == code,]

        weighted_avg_cost = sum(subset_cost_data[,num_visits_field] * subset_cost_data[,unit_cost_field]) / sum(subset_cost_data[,num_visits_field])

        cost_table[cost_table[,currency_code_field] == code, "weighted_avg_all_unit_cost"] = weighted_avg_cost

        cost_table[cost_table[,currency_code_field] == code, "total_fce"] = sum(subset_cost_data[,num_visits_field])

    }

    cost_table = cost_table[!duplicated(cost_table$currency_code),]

    return(cost_table)

}

### Running function for EL

care_costs_el = weighted_avg_maker(cost_table = care_costs_el,
                                currency_code_field = "currency_code",
                                unit_cost_field = "national_avg_unit_cost_el",
                                num_visits_field = "n_fce_el",
                                healthcare_area = "apc_sector_el")

care_costs_el[,"apc_sector_el"] = NULL
colnames(care_costs_el)[6:7] = paste0(colnames(care_costs_el)[6:7], "_el")

### Running function for NEL

care_costs_nel = weighted_avg_maker(cost_table = care_costs_nel,
                                currency_code_field = "currency_code",
                                unit_cost_field = "national_avg_unit_cost_nel",
                                num_visits_field = "n_fce_nel",
                                healthcare_area = "apc_sector_nel")

care_costs_nel[,"apc_sector_nel"] = NULL
colnames(care_costs_nel)[6:7] = paste0(colnames(care_costs_nel)[6:7], "_nel")

### Running function for NES

care_costs_nes = weighted_avg_maker(cost_table = care_costs_nes,
                                currency_code_field = "currency_code",
                                unit_cost_field = "national_avg_unit_cost_nes",
                                num_visits_field = "n_fce_nes",
                                healthcare_area = "apc_sector_nes")

care_costs_nes[,"apc_sector_nes"] = NULL
colnames(care_costs_nes)[6:7] = paste0(colnames(care_costs_nes)[6:7], "_nes")

### Merging data and keeping all codes

combined_cost_frame = merge(care_costs_el, care_costs_nel, by = "currency_code", all = TRUE)
combined_cost_frame = merge(combined_cost_frame, care_costs_nes, by = "currency_code", all = TRUE)

### Creating unit costs based on NEL/NES ratio

combined_cost_frame$nel_nes_sum = rowSums(combined_cost_frame[, c("n_fce_nel", "n_fce_nes")], na.rm = TRUE)
combined_cost_frame$nel_nes_share_nel = combined_cost_frame$n_fce_nel / combined_cost_frame$nel_nes_sum
combined_cost_frame$nel_nes_share_nes = combined_cost_frame$n_fce_nes / combined_cost_frame$nel_nes_sum

combined_cost_frame$nel_nes_unit_cost = (combined_cost_frame$national_avg_unit_cost_nel * combined_cost_frame$nel_nes_share_nel) + (combined_cost_frame$national_avg_unit_cost_nes * combined_cost_frame$nel_nes_share_nes)

### Creating unit costs based on EL/NES ratio

combined_cost_frame$el_nes_sum = rowSums(combined_cost_frame[, c("n_fce_el", "n_fce_nes")], na.rm = TRUE)
combined_cost_frame$el_nes_share_el = combined_cost_frame$n_fce_el / combined_cost_frame$el_nes_sum
combined_cost_frame$el_nes_share_nes = combined_cost_frame$n_fce_nes / combined_cost_frame$el_nes_sum

combined_cost_frame$el_nes_unit_cost = (combined_cost_frame$national_avg_unit_cost_el * combined_cost_frame$el_nes_share_el) + (combined_cost_frame$national_avg_unit_cost_nes * combined_cost_frame$el_nes_share_nes)

### REPEATING COST CREATION WITH THE WEIGHTED AVERAGE
### Creating unit costs based on NEL/NES ratio

combined_cost_frame$weighted_nel_nes_sum = rowSums(combined_cost_frame[, c("total_fce_nel", "total_fce_nes")], na.rm = TRUE)
combined_cost_frame$weighted_nel_nes_share_nel = combined_cost_frame$total_fce_nel / combined_cost_frame$total_fce_nes
combined_cost_frame$weighted_nel_nes_share_nes = combined_cost_frame$total_fce_nes / combined_cost_frame$weighted_nel_nes_sum

combined_cost_frame$weighted_nel_nes_unit_cost = (combined_cost_frame$weighted_avg_all_unit_cost_nel * combined_cost_frame$weighted_nel_nes_share_nes) + (combined_cost_frame$weighted_avg_all_unit_cost_nes * combined_cost_frame$weighted_nel_nes_share_nes)

### Creating unit costs based on EL/NES ratio

combined_cost_frame$weighted_el_nes_sum = rowSums(combined_cost_frame[, c("total_fce_el", "total_fce_nes")], na.rm = TRUE)
combined_cost_frame$weighted_el_nes_share_nel = combined_cost_frame$total_fce_el / combined_cost_frame$total_fce_nes
combined_cost_frame$weighted_el_nes_share_nes = combined_cost_frame$total_fce_nes / combined_cost_frame$weighted_el_nes_sum

combined_cost_frame$weighted_el_nes_unit_cost = (combined_cost_frame$weighted_avg_all_unit_cost_el * combined_cost_frame$weighted_el_nes_share_nes) + (combined_cost_frame$weighted_avg_all_unit_cost_nes * combined_cost_frame$weighted_el_nes_share_nes)

### Creating merge-status column

combined_cost_frame$el_nel_nes_merge = 0

combined_cost_frame = combined_cost_frame %>%
    mutate(el_nel_nes_merge = case_when(
        !is.na(n_fce_nel) & is.na(n_fce_el) & is.na(n_fce_nes) ~ 1, # Only in NEL
        is.na(n_fce_nel) & !is.na(n_fce_el) & is.na(n_fce_nes) ~ 2, # Only in EL
        is.na(n_fce_nel) & is.na(n_fce_el) & !is.na(n_fce_nes) ~ 3, # Only in NES
        !is.na(n_fce_nel) & !is.na(n_fce_el) & is.na(n_fce_nes) ~ 4, # In EL + NEL not NES
        !is.na(n_fce_nel) & is.na(n_fce_el) & !is.na(n_fce_nes) ~ 5, # In NEL + NES not NES
        is.na(n_fce_nel) & !is.na(n_fce_el) & !is.na(n_fce_nes) ~ 6, # In EL + NES
        !is.na(n_fce_nel) & !is.na(n_fce_el) & !is.na(n_fce_nes) ~ 7, # In all three
        TRUE ~ combined_cost_frame$el_nel_nes_merge
        ))

heading("Combined cost frame for EL/NEL/NES created and merged.")

print(head(combined_cost_frame, n = 10))


#####################

### Calculating costs for consultant and non-consultant led

#####################

### Consultant costs read from National Schecule Excel sheet (OP sheet)
### Sheet is converted into data frame

consultant_costs = as.data.frame(consultant_costs_2022)

### Superfluous rows removed

consultant_costs = consultant_costs[-c(1:4),]

### Columns renamed correctly

colnames(consultant_costs) = c("department_code",
                            "department_desc",
                            "service_code",
                            "service_desc",
                            "currency_code",
                            "currency_desc",
                            "n_attendees",
                            "national_avg_unit_cost",
                            "total_cost",
                            "n_submissions"
                            )

### Removing cost rows relating to Paediatrics

consultant_costs = consultant_costs[!(grepl("Paediatric", consultant_costs)),]

#consultant_costs = consultant_costs[, c(5, 7:10)]

### Making sure columns are numeric

consultant_costs$n_attendees = as.numeric(consultant_costs$n_attendees)

consultant_costs$national_avg_unit_cost = as.numeric(consultant_costs$national_avg_unit_cost)
consultant_costs$national_avg_unit_cost[is.na(consultant_costs$national_avg_unit_cost)] = 0 # Replace NAs in unit costs with 0

consultant_costs$total_cost = as.numeric(consultant_costs$total_cost)
consultant_costs$n_submissions = as.numeric(consultant_costs$n_submissions)

### If any values are negative, make them positive

consultant_costs$national_avg_unit_cost[consultant_costs$national_avg_unit_cost < 0] = consultant_costs$national_avg_unit_cost[consultant_costs$national_avg_unit_cost < 0] * -1

### Restrict cost rows to CL (consultant led)

consultant_costs_cl = consultant_costs[consultant_costs$department_code == "CL",]

### Grouping currency codes by the total number who have attended each code
### I.e. multiple instances of each code, with numerous attendees, so combine n_submissions

consultant_costs_cl = consultant_costs_cl %>% 
    group_by(currency_code) %>% 
    mutate(total_attend = sum(n_submissions))

### Divide individual attendences by total overall to get proportions

consultant_costs_cl$attend_share = consultant_costs_cl$n_submissions / consultant_costs_cl$total_attend

### Weighted weighted cost by multiplying unit cost by proportion

consultant_costs_cl$weighted_unit_cost = consultant_costs_cl$attend_share * consultant_costs_cl$national_avg_unit_cost

### Grouped total of weighted unit costs

consultant_costs_cl = consultant_costs_cl %>% 
    group_by(currency_code) %>% 
    mutate(cl_weighted_unit_cost = sum(weighted_unit_cost))

### Keep only data of unique codes + unique costs

consultant_costs_cl = as.data.frame(distinct(consultant_costs_cl, currency_code, .keep_all = TRUE))

### Restrict cost rows to NCL (consultant led)
### Repeat process for CL on NCL

consultant_costs_ncl = consultant_costs[consultant_costs$department_code == "NCL",]

### Grouping currency codes by the total number who have attended each code
### I.e. multiple instances of each code, with numerous attendees, so combine n_submissions

consultant_costs_ncl = consultant_costs_ncl %>% 
    group_by(currency_code) %>% 
    mutate(total_attend = sum(n_submissions))

### Divide individual attendences by total overall to get proportions

consultant_costs_ncl$attend_share = consultant_costs_ncl$n_submissions / consultant_costs_ncl$total_attend

### Weighted weighted cost by multiplying unit cost by proportion

consultant_costs_ncl$weighted_unit_cost = consultant_costs_ncl$attend_share * consultant_costs_ncl$national_avg_unit_cost

### Grouped total of weighted unit costs

consultant_costs_ncl = consultant_costs_ncl %>% 
    group_by(currency_code) %>% 
    mutate(ncl_weighted_unit_cost = sum(weighted_unit_cost))

### Keep only data of unique codes + unique costs

consultant_costs_ncl = as.data.frame(distinct(consultant_costs_ncl, currency_code, .keep_all = TRUE))

### Merge NCL + CL datasets

colnames(consultant_costs_ncl)[-5] = paste0(colnames(consultant_costs_ncl)[-5], "_ncl")
colnames(consultant_costs_cl)[-5] = paste0(colnames(consultant_costs_cl)[-5], "_cl")

merged_consultant_costs = merge(consultant_costs_cl[,c(5, 7:14)],
                                consultant_costs_ncl[,c(5, 7:14)],
                                by = "currency_code")

### Creating total attendances column

merged_consultant_costs$total_attendances = merged_consultant_costs$total_attend_ncl + merged_consultant_costs$total_attend_cl

### Creating proportional shares of visits

merged_consultant_costs$ncl_total_share = merged_consultant_costs$total_attend_ncl / merged_consultant_costs$total_attendances
merged_consultant_costs$cl_total_share = merged_consultant_costs$total_attend_cl / merged_consultant_costs$total_attendances

### Creating final weighted CL/NCL unit cost

merged_consultant_costs$weighted_cl_ncl_cost = merged_consultant_costs$ncl_total_share * merged_consultant_costs$ncl_weighted_unit_cost_ncl + 
                                                merged_consultant_costs$cl_total_share * merged_consultant_costs$cl_weighted_unit_cost_cl

### Data frame with only code and weighted cost

consultant_cost_fees = merged_consultant_costs[, c("currency_code", "weighted_cl_ncl_cost")]

heading("Consultant cost frame for CL/NCL created and merged.")

print(head(consultant_cost_fees, n = 10))


#####################

### Creating Zero-cost sheet

#####################

zero_cost_fees = data.frame(currency_code = c("DZ13A", "DZ13B", "LA97A", "LA97B", "SB97Z", "SC97Z", "RD97Z", "RN97Z"),
                                            unit_cost = 0)


#####################

### Extracting day-care costs for DC/RP units

#####################

### Extracting DC coding "Day care" and "Regular Day or Night Admissions"

care_costs_dc = care_costs[care_costs$department_code == "DC", c(1,3,5:9)]
care_costs_rp = care_costs[care_costs$department_code == "RP", c(1,3,5:9)]

### Re-naming columns for merging

colnames(care_costs_dc) = paste0(colnames(care_costs_dc), "_dc")
colnames(care_costs_rp) = paste0(colnames(care_costs_rp), "_rp")

colnames(care_costs_dc)[2] = "currency_code"
colnames(care_costs_rp)[2] = "currency_code"

### Creating community-weighted costs for both

care_costs_dc$department_code_dc = NULL

care_costs_dc = weighted_avg_maker(cost_table = care_costs_dc,
                                currency_code_field = "currency_code",
                                unit_cost_field = "national_avg_unit_cost_dc",
                                num_visits_field = "n_fce_dc",
                                healthcare_area = "apc_sector_dc")

care_costs_dc[,"apc_sector_dc"] = NULL
colnames(care_costs_dc)[6:7] = paste0(colnames(care_costs_dc)[6:7], "_dc")

care_costs_rp$department_code_rp = NULL

care_costs_rp = weighted_avg_maker(cost_table = care_costs_rp,
                                currency_code_field = "currency_code",
                                unit_cost_field = "national_avg_unit_cost_rp",
                                num_visits_field = "n_fce_rp",
                                healthcare_area = "apc_sector_rp")

care_costs_rp[,"apc_sector_rp"] = NULL
colnames(care_costs_rp)[6:7] = paste0(colnames(care_costs_rp)[6:7], "_rp")

heading("Day care, regular visits and non-cost frames created and merged.")

print(head(care_costs_rp, n = 10))


#####################

### Extracting excess bed day costs for EL/NEL patients

#####################

### Care costs read from National Schecule Excel sheet (one sheet of several)
### Sheet is converted into data frame

excess_el_care_costs = as.data.frame(elective_excess_costs_2017)
excess_nel_care_costs = as.data.frame(non_elective_excess_costs_2017)

### Superfluous rows removed

excess_el_care_costs = excess_el_care_costs[-c(1:4),]
excess_nel_care_costs = excess_nel_care_costs[-c(1:4),]

### Columns renamed correctly

colnames(excess_el_care_costs) = c("currency_code",
                            "currency_desc",
                            "excess_bed_days",
                            "national_avg_unit_cost",
                            "lower_quartile_unit_cost",
                            "upper_quartile_unit_cost",
                            "n_submissions",
                            "total_cost")

colnames(excess_nel_care_costs) = c("currency_code",
                            "currency_desc",
                            "excess_bed_days",
                            "national_avg_unit_cost",
                            "lower_quartile_unit_cost",
                            "upper_quartile_unit_cost",
                            "n_submissions",
                            "total_cost")

### Ensuring columns are numeric correctly
### Remove NAs
### Multiply by 2017 -> 2022 index in order to bring up to 2022 costs                        

excess_el_care_costs$excess_bed_days = as.numeric(excess_el_care_costs$excess_bed_days)
excess_el_care_costs$excess_bed_days[is.na(excess_el_care_costs$excess_bed_days)] = 0
excess_el_care_costs$excess_bed_days = excess_el_care_costs$excess_bed_days * cpi_index_2017_to_2022

excess_el_care_costs$national_avg_unit_cost = as.numeric(excess_el_care_costs$national_avg_unit_cost)
excess_el_care_costs$national_avg_unit_cost[is.na(excess_el_care_costs$national_avg_unit_cost)] = 0
excess_el_care_costs$national_avg_unit_cost = excess_el_care_costs$national_avg_unit_cost * cpi_index_2017_to_2022

excess_nel_care_costs$excess_bed_days = as.numeric(excess_nel_care_costs$excess_bed_days)
excess_nel_care_costs$excess_bed_days[is.na(excess_nel_care_costs$excess_bed_days)] = 0
excess_nel_care_costs$excess_bed_days = excess_nel_care_costs$excess_bed_days * cpi_index_2017_to_2022

excess_nel_care_costs$national_avg_unit_cost = as.numeric(excess_nel_care_costs$national_avg_unit_cost)
excess_nel_care_costs$national_avg_unit_cost[is.na(excess_nel_care_costs$national_avg_unit_cost)] = 0
excess_nel_care_costs$national_avg_unit_cost = excess_nel_care_costs$national_avg_unit_cost * cpi_index_2017_to_2022

heading("Excess costs for EL-XS/NEL-XS created..")

print(head(excess_nel_care_costs, n = 10))

#####################
#####################

### END OF COST CODE parsing

#####################
#####################


#####################

### Grouper data modulation

#####################

### Grouper data read as readLines - convert into data frame

grouper_columns = strsplit(grouper_output[1], ",")[[1]]

### Split the lines by comma

split_grouper_output = strsplit(grouper_output[-1], ",")

### Calculate the maximum number of columns

max_grouper_cols = max(lengths(split_grouper_output))

### Re-split the grouper file with the max number of columns:

split_grouper_output = str_split_fixed(grouper_output[-1], ",", max_grouper_cols)

### Re-name columns of new table

unbundled_columns = grouper_columns[length(grouper_columns)]
unbundled_columns = paste(unbundled_columns, c((length(grouper_columns)-(length(grouper_columns)-1)):(max_grouper_cols-(length(grouper_columns)-1))), sep = "_")

colnames(split_grouper_output) = c(grouper_columns[-length(grouper_columns)], unbundled_columns)
grouper_output = as.data.frame(split_grouper_output)

### Recode classpats of some Grouper outputs based on current ClassPat and currency code

for (i in 1:nrow(grouper_classpat_recode_list)){

    currCode = grouper_classpat_recode_list$code[i]

    currClasspat_ins = grouper_classpat_recode_list$classpat_in[i]
    currClasspat_recode = grouper_classpat_recode_list$classpat_out[i]

    currClasspat_ins = str_split_fixed(currClasspat_ins, ",", 2)

    curr_grouper_index = which(grouper_output$FCE_HRG == currCode)

    lines_to_change = curr_grouper_index[which(grouper_output$CLASSPAT[curr_grouper_index] == currClasspat_ins[,1] | grouper_output$CLASSPAT[curr_grouper_index] == currClasspat_ins[,2])]

    if (length(lines_to_change) > 0){

        heading("The following lines need to be fixed in the Grouper output:")

        print(grouper_output[lines_to_change,])

        grouper_output$CLASSPAT[lines_to_change] = currClasspat_recode

        heading(paste0("This has been changed from ", currClasspat_ins[,1], " OR ", currClasspat_ins[,2], " to ", currClasspat_recode))

        print(grouper_output[lines_to_change,])

    }
}

### Create elective column
### Non-elective as standard, change to 1 based on ADMIMETH

grouper_output$Elective = 0

grouper_output$Elective[which(grouper_output$ADMIMETH == "11" | grouper_output$ADMIMETH == "12" | grouper_output$ADMIMETH == "13" | grouper_output$ADMIMETH == "81")] = 1

### Removing instances where the care was not completed based on EPISTAT (from original HESIN data)

hes_data$iid = paste0(hes_data$eid, "_", hes_data$ins_index)

hes_data = hes_data[hes_data$iid %in% grouper_output$PROVSPNO,]

iids_to_remove = hes_data$iid[which(hes_data$epistat == 1)]

grouper_output = grouper_output[!(grouper_output$PROVSPNO %in% iids_to_remove),]

### Creating IN/OUT patient column
### THIS IS NOT IN USE

hes_data$epistart = as.Date(hes_data$epistart, format = "%d/%m/%Y")
hes_data$epiend = as.Date(hes_data$epiend, format = "%d/%m/%Y")

hes_data$healthcare_format = ifelse(hes_data$epistart == hes_data$epiend, "OUT", "IN")

colnames(hes_data)[6] = "PROVSPNO"

### Adding epistart and date_of_entry columns to Grouper data

hes_data = merge(hes_data, ukb_phenos, by = "eid")

grouper_output = merge(grouper_output, hes_data[,c("PROVSPNO", "healthcare_format", "epistart", "epiend", "date_of_entry")], by = "PROVSPNO")

### Grouper outputs where the start of episode is prior to the entry to UKB are removed from the dataset
### Grouper outputs where the end of episode occurs after the defined end date are also removed
### END DATE DEFINED AS - max date of HES data minus one year (to account for time for registry data to catch up)

grouper_output = grouper_output[grouper_output$date_of_entry <= grouper_output$epistart,]

maxHesDate = max(hes_data$epiend[!(is.na(hes_data$epiend))]) - 365

heading(paste0("Last date for HES data is counted as: ", maxHesDate))

grouper_output = grouper_output[grouper_output$epiend < maxHesDate,]

heading("Grouper output data formatting complete.")

print(head(grouper_output, n = 10))


#####################

### Obtaining final APC cost

#####################

### Remove invalid Grouper code from output

error_individuals = grouper_output[grouper_output$FCE_HRG == "UZ01Z",]
error_iids = unique(str_split_fixed(error_individuals$PROVSPNO, "_", 2)[,1])

grouper_output = grouper_output[!(grouper_output$FCE_HRG == "UZ01Z"),]

### Creating data frame of IIDs + Classpat + all the costing values

costing_frame = data.frame(iid = grouper_output$PROVSPNO,
                            currency_code = grouper_output$FCE_HRG,
                            elective = grouper_output$Elective,
                            classpat = grouper_output$CLASSPAT)

costing_frame = costing_frame[costing_frame$classpat == 1,]

costing_frame = merge(costing_frame, combined_cost_frame[,c("currency_code",
                                                            "national_avg_unit_cost_el",
                                                            "national_avg_unit_cost_nel",
                                                            "national_avg_unit_cost_nes",
                                                            "el_nes_unit_cost",
                                                            "nel_nes_unit_cost",
                                                            "weighted_avg_all_unit_cost_el",
                                                            "weighted_avg_all_unit_cost_nel",
                                                            "weighted_avg_all_unit_cost_nes",
                                                            "weighted_el_nes_unit_cost",
                                                            "weighted_nel_nes_unit_cost",
                                                            "el_nel_nes_merge"
                                                            )],
                        by = "currency_code", all = TRUE)

### Merge introduces NA IIDs where no codes were present, remove those

costing_frame = costing_frame[which(!(is.na(costing_frame$iid))),]

### Create final APC cost columns - normal + weighted-average to include the Community costs

costing_frame$unit_cost = 0
costing_frame$weighted_unit_cost = 0

### Create cost values based on merge status column created earlier
### Some currency codes are not in EL/NEL/NES, so a merge value is created for them (0)

costing_frame$el_nel_nes_merge[which(is.na(costing_frame$el_nel_nes_merge))] = 0

costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 1] = costing_frame$national_avg_unit_cost_nel[costing_frame$el_nel_nes_merge == 1]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 2] = costing_frame$national_avg_unit_cost_el[costing_frame$el_nel_nes_merge == 2]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 3] = costing_frame$national_avg_unit_cost_nes[costing_frame$el_nel_nes_merge == 3]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 5] = costing_frame$nel_nes_unit_cost[costing_frame$el_nel_nes_merge == 5]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 6] = costing_frame$el_nes_unit_cost[costing_frame$el_nel_nes_merge == 6]

costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 1] = costing_frame$weighted_avg_all_unit_cost_nel[costing_frame$el_nel_nes_merge == 1]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 2] = costing_frame$weighted_avg_all_unit_cost_el[costing_frame$el_nel_nes_merge == 2]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 3] = costing_frame$weighted_avg_all_unit_cost_nes[costing_frame$el_nel_nes_merge == 3]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 5] = costing_frame$weighted_nel_nes_unit_cost[costing_frame$el_nel_nes_merge == 5]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 6] = costing_frame$weighted_el_nes_unit_cost[costing_frame$el_nel_nes_merge == 6]

### Where BOTH EL and NEL values are present (I.e. merge status 4 or 7) the cost value is based on the elective column

costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 0] = costing_frame$nel_nes_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 0]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 1] = costing_frame$el_nes_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 1]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 0] = costing_frame$national_avg_unit_cost_nel[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 0]
costing_frame$unit_cost[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 1] = costing_frame$national_avg_unit_cost_el[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 1]

costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 0] = costing_frame$weighted_nel_nes_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 0]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 1] = costing_frame$weighted_el_nes_unit_cost[costing_frame$el_nel_nes_merge == 7 & costing_frame$elective == 1]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 0] = costing_frame$weighted_avg_all_unit_cost_nel[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 0]
costing_frame$weighted_unit_cost[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 1] = costing_frame$weighted_avg_all_unit_cost_el[costing_frame$el_nel_nes_merge == 4 & costing_frame$elective == 1]

### Adding costs from CL/NCL data
### Only WF01A and WF02A present in UKB data (plus two 0 costs)

costing_frame$unit_cost[costing_frame$currency_code == "WF01A"] = consultant_cost_fees$weighted_cl_ncl_cost[consultant_cost_fees$currency_code == "WF01A"]
costing_frame$unit_cost[costing_frame$currency_code == "WF02A"] = consultant_cost_fees$weighted_cl_ncl_cost[consultant_cost_fees$currency_code == "WF02A"]

costing_frame$weighted_unit_cost[costing_frame$currency_code == "WF01A"] = consultant_cost_fees$weighted_cl_ncl_cost[consultant_cost_fees$currency_code == "WF01A"]
costing_frame$weighted_unit_cost[costing_frame$currency_code == "WF02A"] = consultant_cost_fees$weighted_cl_ncl_cost[consultant_cost_fees$currency_code == "WF02A"]

costing_frame = costing_frame[,c("iid", "currency_code", "unit_cost", "weighted_unit_cost")]

heading("APC-costs per IID created.")

print(head(costing_frame, n = 10))


#####################

### Obtaining DC + RP cost

#####################

### Restrict IIDs to classpat == 2

day_care_costing_frame = data.frame(iid = grouper_output$PROVSPNO,
                            currency_code = grouper_output$FCE_HRG,
                            classpat = grouper_output$CLASSPAT)

day_care_costing_frame = day_care_costing_frame[day_care_costing_frame$classpat == 2,]

### Cost value

day_care_costing_frame = merge(day_care_costing_frame, care_costs_dc[,c("currency_code",
                                                                        "national_avg_unit_cost_dc",
                                                                        "weighted_avg_all_unit_cost_dc")],
                                by = "currency_code", all = TRUE)

day_care_costing_frame = day_care_costing_frame[which(!(is.na(day_care_costing_frame$iid))),]


### Restrict IIDs to classpat == 3 | 4

regular_care_costing_frame = data.frame(iid = grouper_output$PROVSPNO,
                            currency_code = grouper_output$FCE_HRG,
                            classpat = grouper_output$CLASSPAT)

regular_care_costing_frame = regular_care_costing_frame[regular_care_costing_frame$classpat == 3 | regular_care_costing_frame$classpat == 4,]

### Cost value

regular_care_costing_frame = merge(regular_care_costing_frame, care_costs_rp[,c("currency_code",
                                                                        "national_avg_unit_cost_rp",
                                                                        "weighted_avg_all_unit_cost_rp")],
                                by = "currency_code", all = TRUE)

regular_care_costing_frame = regular_care_costing_frame[which(!(is.na(regular_care_costing_frame$iid))),]


#####################

### Obtaining Maternity cost

#####################


### Restrict IIDs to classpat == 5
### Original code also used classpat == 8 but there are no such codes in the data

maternity_care_costing_frame = data.frame(iid = grouper_output$PROVSPNO,
                            currency_code = grouper_output$FCE_HRG,
                            classpat = grouper_output$CLASSPAT)

maternity_care_costing_frame = maternity_care_costing_frame[maternity_care_costing_frame$classpat == 5,]

### Costs are taken from the HRG care costs (as opposed to APC)
### Care costs read from National Schecule Excel sheet (one sheet of several)
### Sheet is converted into data frame

hrg_care_costs = as.data.frame(unbun_costs_2022)

### Superfluous rows and columns removed

maternal_care_costs = hrg_care_costs[-c(1:5),c(1:5)]

### Columns renamed correctly

trio_col_names = c("_activity", "_unit_cost", "_total_cost")

colnames(maternal_care_costs) = c("currency_code",
                            "currency_desc",
                            paste0("total", trio_col_names))

### Making sure total cost is numeric correctly

maternal_care_costs$total_unit_cost = as.numeric(maternal_care_costs$total_unit_cost)
maternal_care_costs$total_unit_cost[is.na(maternal_care_costs$total_unit_cost)] = 0
maternal_care_costs$total_unit_cost[maternal_care_costs$total_unit_cost < 0] = maternal_care_costs$total_unit_cost[maternal_care_costs$total_unit_cost < 0] * -1

### Extracting cost

maternity_care_costing_frame$unit_cost_mat = maternal_care_costs$total_unit_cost[match(maternity_care_costing_frame[["currency_code"]], maternal_care_costs$currency_code)]


#####################

### Obtaining Unbundled HRGs cost

#####################

### Define unbun frame

unbun_costing_frame = data.frame(iid = grouper_output$PROVSPNO,
                            epidur = grouper_output$EPIDUR,
                            classpat = grouper_output$CLASSPAT, # Needed to pick between day cases/RP etc
                            elective = grouper_output$Elective, # To differentiate between EL/NEL
                            total_cost = 0,
                            temp_cost = 0,
                            grouper_output[,93:114]) # All unbundledHRG columns

### Sortinmg the unbun care costs
### Care costs read from National Schecule Excel sheet (one sheet of several)
### Sheet is converted into data frame

hrg_care_costs = as.data.frame(unbun_costs_2022)

### Superfluous rows and columns removed

unbun_care_costs = hrg_care_costs[-c(1:5),c(1:20)]

### Columns renamed correctly

trio_col_names = c("_activity", "_unit_cost", "_total_cost")

colnames(unbun_care_costs) = c("currency_code",
                            "currency_desc",
                            paste0("total", trio_col_names),
                            paste0("elective", trio_col_names),
                            paste0("non_elective_ls", trio_col_names),
                            paste0("non_elective_ss", trio_col_names),
                            paste0("day_case", trio_col_names),
                            paste0("day_night_admi", trio_col_names))

### Ensuring columns are numeric correctly                            

unbun_care_costs$total_unit_cost = as.numeric(unbun_care_costs$total_unit_cost)
unbun_care_costs$elective_unit_cost = as.numeric(unbun_care_costs$elective_unit_cost)
unbun_care_costs$non_elective_ls_unit_cost = as.numeric(unbun_care_costs$non_elective_ls_unit_cost)
unbun_care_costs$non_elective_ss_unit_cost = as.numeric(unbun_care_costs$non_elective_ss_unit_cost)
unbun_care_costs$day_case_unit_cost = as.numeric(unbun_care_costs$day_case_unit_cost)
unbun_care_costs$day_night_admi_unit_cost = as.numeric(unbun_care_costs$day_night_admi_unit_cost)

unbun_care_costs$total_unit_cost[is.na(unbun_care_costs$total_unit_cost)] = 0
unbun_care_costs$elective_unit_cost[is.na(unbun_care_costs$elective_unit_cost)] = 0
unbun_care_costs$non_elective_ls_unit_cost[is.na(unbun_care_costs$non_elective_ls_unit_cost)] = 0
unbun_care_costs$non_elective_ss_unit_cost[is.na(unbun_care_costs$non_elective_ss_unit_cost)] = 0
unbun_care_costs$day_case_unit_cost[is.na(unbun_care_costs$day_case_unit_cost)] = 0
unbun_care_costs$day_night_admi_unit_cost[is.na(unbun_care_costs$day_night_admi_unit_cost)] = 0

### Making negative values positive

unbun_care_costs$total_unit_cost[unbun_care_costs$total_unit_cost < 0] = unbun_care_costs$total_unit_cost[unbun_care_costs$total_unit_cost < 0] * -1
unbun_care_costs$elective_unit_cost[unbun_care_costs$elective_unit_cost < 0] = unbun_care_costs$elective_unit_cost[unbun_care_costs$elective_unit_cost < 0] * -1
unbun_care_costs$non_elective_ls_unit_cost[unbun_care_costs$non_elective_ls_unit_cost < 0] = unbun_care_costs$non_elective_ls_unit_cost[unbun_care_costs$non_elective_ls_unit_cost < 0] * -1
unbun_care_costs$non_elective_ss_unit_cost[unbun_care_costs$non_elective_ss_unit_cost < 0] = unbun_care_costs$non_elective_ss_unit_cost[unbun_care_costs$non_elective_ss_unit_cost < 0] * -1
unbun_care_costs$day_case_unit_cost[unbun_care_costs$day_case_unit_cost < 0] = unbun_care_costs$day_case_unit_cost[unbun_care_costs$day_case_unit_cost < 0] * -1
unbun_care_costs$day_night_admi_unit_cost[unbun_care_costs$day_night_admi_unit_cost < 0] = unbun_care_costs$day_night_admi_unit_cost[unbun_care_costs$day_night_admi_unit_cost < 0] * -1

### Separating unbundled frame by various types of cost
### MAY NOT BE NECESSARY

#unbun_costing_frame_el = unbun_costing_frame[unbun_costing_frame$elective == 1,]
#unbun_costing_frame_nel_ss = unbun_costing_frame[unbun_costing_frame$elective == 0 & !(unbun_costing_frame$epidur > 0),]
#unbun_costing_frame_nel_ls = unbun_costing_frame[unbun_costing_frame$elective == 0 & unbun_costing_frame$epidur > 0,]
#unbun_costing_frame_dc = unbun_costing_frame[unbun_costing_frame$classpat == 2,]
#unbun_costing_frame_rp = unbun_costing_frame[unbun_costing_frame$classpat == 3 | unbun_costing_frame$classpat == 4,]

### Creating summed column of costs
### Iterate over each HRG column

for (i in 7:ncol(unbun_costing_frame)) {

  ### Get the column name

  curr_code_column = names(unbun_costing_frame)[i]
  
  ### Replace empty strings with 0 in the code column

  unbun_costing_frame[,curr_code_column][unbun_costing_frame[,curr_code_column] == ""] = 0
  
  ### Match the codes with costs and calculate the total cost

  unbun_costing_frame$temp_cost = unbun_care_costs$total_unit_cost[match(unbun_costing_frame[[curr_code_column]], unbun_care_costs$currency_code)]

  unbun_costing_frame$total_cost = rowSums(unbun_costing_frame[,c("total_cost", "temp_cost")], na.rm = TRUE)

}

heading("Unbundled costs and maternity costs frames crated.")

print(head(unbun_costing_frame, n = 10))


#####################

### Creating final cost frame

#####################

### Merge each cost frame by IID (individual ID + installment at hospital)

final_cost_frame = data.frame(iid = grouper_output$PROVSPNO)

final_cost_frame = merge(final_cost_frame, costing_frame[,c(1,3,4)], by = "iid", all = TRUE)
colnames(final_cost_frame)[2:3] = c("apc_cost", "weighted_apc_cost")

final_cost_frame = merge(final_cost_frame, regular_care_costing_frame[,c(2,4,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[4:5] = c("rp_cost", "weighted_rp_cost")

final_cost_frame = merge(final_cost_frame, day_care_costing_frame[,c(2,4,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[6:7] = c("dc_cost", "weighted_dc_cost")

final_cost_frame = merge(final_cost_frame, maternity_care_costing_frame[,c(1,4)], by = "iid", all = TRUE)
colnames(final_cost_frame)[8] = c("mat_cost")

final_cost_frame = merge(final_cost_frame, unbun_costing_frame[,c(1,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[9] = c("unbun_cost")

### Sum of all costs and sum of weighted costs

final_cost_frame$hes_costs = rowSums(final_cost_frame[,c(2,4,6,8,9)], na.rm = TRUE)
final_cost_frame$weighted_hes_costs = rowSums(final_cost_frame[,c(3,5,7:9)], na.rm = TRUE)

### Finally create sum of each hospital instance for each individual

final_cost_frame$eid = str_split_fixed(final_cost_frame$iid, "_", 2)[,1]

### Making the final phenotype data frame: final_costs_individual

final_costs_individual = aggregate(cbind(hes_costs, weighted_hes_costs) ~ eid, data = final_cost_frame, FUN = sum)

### Add in individuals who are not present in HES data at all as 0 cost

hes_zero_cost_individuals = data.frame(eid = ukb_phenos$eid[!(ukb_phenos$eid %in% unique(hes_data$eid))])
hes_zero_cost_individuals$hes_costs = 0
hes_zero_cost_individuals$weighted_hes_costs = 0

hes_zero_cost_individuals = hes_zero_cost_individuals[!(hes_zero_cost_individuals$eid %in% error_iids),]

final_costs_individual = rbind(final_costs_individual, hes_zero_cost_individuals)

### Adding in end of follow-up date

final_costs_individual$hes_end_of_followup = maxHesDate

### Adding last date for HES data recorded
### SUPERCEDED

#hes_dates = grouper_output[,c("PROVSPNO", "epiend")]
#
#hes_dates$eid = str_split_fixed(grouper_output$PROVSPNO, "_", 2)[,1]
#hes_dates$PROVSPNO = NULL
#
#most_recent_dates = hes_dates %>%
#                      group_by(eid) %>%
#                        summarize(epiend = max(epiend))
#
#colnames(most_recent_dates)[2] = "hes_end_of_followup"
#
#final_costs_individual = merge(final_costs_individual, most_recent_dates, by = "eid")

fwrite(final_costs_individual, paste0(mainDir, "processing/intermediate_files/hes_cost_phenotypes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

heading("Final cost fame for each individual being parsed.")

print(head(final_costs_individual, n = 10))

#####################

### Calculate GP visit costs

#####################

final_costs_individual = fread(paste0(mainDir, "processing/intermediate_files/hes_cost_phenotypes.tsv"), data.table = FALSE)

gp_visits = gp_visits[,c("eid", "event_dt")]

### GP visits kept between start and end date

gp_visits$event_dt = as.Date(gp_visits$event_dt, format = "%d/%m/%Y")
gp_visits = gp_visits[complete.cases(gp_visits$event_dt),]
gp_visits = gp_visits[gp_visits$event_dt < end_date,]
gp_visits = unique(gp_visits) # Remove multiple visits occuring on the same day - seems unlikely

### Get GP start of follow-up
### Group by EID and calculate the earliest date for each individual

gp_start = gp_visits %>%
	group_by(eid) %>%
	summarise(gp_start_of_followup = min(event_dt))

### Merge with date_of_entry to determine the start of follow-up

gp_start = merge(gp_start, ukb_phenos[, c("eid", "date_of_entry")], by = "eid", all.x = TRUE)

### Replace gp_start_of_followup with date_of_entry if it's earlier

gp_start$gp_start_of_followup = ifelse(gp_start$gp_start_of_followup < gp_start$date_of_entry, 
                                                gp_start$date_of_entry, 
                                                gp_start$gp_start_of_followup)

gp_start$gp_start_of_followup = as.Date(gp_start$gp_start_of_followup)

### Select required columns

gp_start = gp_start[, c("eid", "gp_start_of_followup")]
gp_start = unique(gp_start)

### Merge gp_visits with UKB data
### Remove all GP visits which occur prior to entry into UKB

gp_visits = merge(gp_visits, gp_start, by = "eid")

gp_visits = gp_visits[gp_visits$event_dt >= gp_visits$gp_start_of_followup,]

maxGPDate = max(gp_visits$event_dt[!(is.na(gp_visits$event_dt))])

heading(paste0("Last date for GP data is counted as: ", maxGPDate))

### Tally total number of visits per individual

gp_visit_counts = table(gp_visits$eid)

gp_visit_counts = data.frame(gp_visit_counts)
colnames(gp_visit_counts) = c("eid", "count")

gp_visit_counts$eid = as.character(gp_visit_counts$eid)
gp_visit_counts$count = as.numeric(gp_visit_counts$count)

### Add extra IIDs

gp_missing_eids = gp_start$eid[!(gp_start$eid %in% gp_visit_counts$eid)]

gp_missing_counts = data.frame(eid = gp_missing_eids,
                                count = 0)

gp_visit_counts = rbind(gp_visit_counts, gp_missing_counts)

### Cost calculated as £42 per visit
### Based on: https://www.kingsfund.org.uk/insight-and-analysis/data-and-charts/key-facts-figures-nhs

gp_visit_counts$gp_costs = gp_visit_counts$count * 42

### Combined with outpatient healthcare cost

individual_all_costs = merge(gp_visit_counts[,c(1,3)], final_costs_individual, by = "eid", all = TRUE)

individual_all_costs$gp_end_of_followup = maxGPDate
individual_all_costs = merge(individual_all_costs, gp_start, by = "eid", all = TRUE)

### Adding last date for GP data recorded
### SUPERCEDED

#gp_dates = gp_visits[,c("eid", "event_dt")]
#
#most_recent_gp_dates = gp_dates %>%
#                      group_by(eid) %>%
#                        summarize(event_dt = max(event_dt))
#
#colnames(most_recent_gp_dates)[2] = "gp_end_of_followup"
#
#individual_all_costs = merge(individual_all_costs, most_recent_gp_dates, by = "eid", all = TRUE)

### Adding first date for GP data recorded
### SUPERCEDED

#gp_dates = gp_visits[,c("eid", "event_dt")]
#
#oldest_gp_dates = gp_dates %>%
#                      group_by(eid) %>%
#                        summarize(event_dt = min(event_dt))
#
#colnames(oldest_gp_dates)[2] = "gp_start_of_followup"

individual_all_costs = merge(individual_all_costs, ukb_phenos, by = "eid", all = TRUE)

individual_all_costs$eid = as.character(individual_all_costs$eid)

### Some individuals were not added in as zeroes in the previous round
### These are individuals who are not in the costed grouper output but WERE in the HES data
### They are corrected here

individual_all_costs$hes_costs[is.na(individual_all_costs$hes_costs)] = 0
individual_all_costs$weighted_hes_costs[is.na(individual_all_costs$weighted_hes_costs)] = 0

heading("Added GP visits to final cost frame.")

print(head(individual_all_costs, n = 10))

fwrite(individual_all_costs, paste0(mainDir, "processing/intermediate_files/gp_hes_cost_phenotypes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


#####################

### Creating prescription costs per person

#####################

#####################

### Initial data parsing

#####################

### Remove rows which are all NAs or have NA dates

gp_prescriptions = gp_prescriptions[complete.cases(gp_prescriptions[,c("eid", "issue_date", "drug_name", "quantity")]),]

### Remove rows which have blanks

gp_prescriptions = gp_prescriptions[!(gp_prescriptions$issue_date == ""),]
gp_prescriptions = gp_prescriptions[!(gp_prescriptions$drug_name == ""),]
gp_prescriptions = gp_prescriptions[!(gp_prescriptions$quantity == ""),]

### Convert to date format

gp_prescriptions$issue_date = as.Date(gp_prescriptions$issue_date, format = "%d/%m/%Y")

### Get prescription start of follow-up
### Group by EID and calculate the earliest date for each individual

prescript_start = gp_prescriptions %>%
	group_by(eid) %>%
	summarise(prescript_start_of_followup = min(issue_date))

### Remove all entries prior to 1st of April, 2006

gp_prescriptions = gp_prescriptions[gp_prescriptions$issue_date > "2006-06-01",] # This is done just to trim down the dataset to make later date parsing easier
gp_prescriptions = gp_prescriptions[gp_prescriptions$issue_date < end_date,] # Done due to some erroneous dates occurring in the future

maxPrescriptDate = max(gp_prescriptions$issue_date)

heading(paste0("Last date for Prescription data is counted as: ", maxPrescriptDate))

### Merge prescription start with date_of_entry to determine the start of follow-up

prescript_start = merge(prescript_start, ukb_phenos[, c("eid", "date_of_entry")], by = "eid", all.x = TRUE)

### Replace gp_start_of_followup with date_of_entry if it's earlier

prescript_start$prescript_start_of_followup = ifelse(prescript_start$prescript_start_of_followup < prescript_start$date_of_entry, 
                                         prescript_start$date_of_entry, 
                                         prescript_start$prescript_start_of_followup)

### Select required columns

prescript_start = prescript_start[, c("eid", "prescript_start_of_followup")]
prescript_start = unique(prescript_start)

prescript_start$prescript_start_of_followup = as.Date(prescript_start$prescript_start_of_followup)

### Remove all prescriptions for individuals prior to entry into UKB

gp_prescriptions_dates = gp_prescriptions[,c("eid", "issue_date")]

gp_prescriptions_dates = merge(gp_prescriptions_dates, ukb_phenos, by = "eid")

gp_prescriptions = gp_prescriptions[gp_prescriptions_dates$issue_date > gp_prescriptions_dates$date_of_entry,]

heading("Begin parsing GP prescrptions.")

print(head(gp_prescriptions, n = 10))


#####################

### Translating prescription names in gp_prescriptions into more consistent names for costing

#####################

### Some instances of prescriptions containining unusual space characters

gp_prescriptions$drug_name = gsub(" ", " ", gp_prescriptions$drug_name)

### Also instances of double spaces, removing those may fix other names

gp_prescriptions$drug_name = gsub("  ", " ", gp_prescriptions$drug_name)

### Doing the same with the comparison list

remove_naughty_spaces = function(df_vector){

  df_vector = gsub(" ", " ", df_vector)
  df_vector = gsub("  ", " ", df_vector)

  return(df_vector)

}

prescriptions_appliances_1[,1] = remove_naughty_spaces(prescriptions_appliances_1[,1])
prescriptions_appliances_2$product_part1 = remove_naughty_spaces(prescriptions_appliances_2$product_part1)
prescriptions_costs$drug = remove_naughty_spaces(prescriptions_costs$drug)
prescriptions_incremental_costs$product = remove_naughty_spaces(prescriptions_incremental_costs$drug)
outdated_prescriptions$outdated_name = remove_naughty_spaces(outdated_prescriptions$outdated_name)
prescriptions_curated_translations$alt_name = remove_naughty_spaces(prescriptions_curated_translations$alt_name)
prescriptions_curated_translations$translation = remove_naughty_spaces(prescriptions_curated_translations$translation)

### Testing the number of prescription name matches before and after name translation

costed_prescription_list = unique(c(prescriptions_appliances_1[,1],
                                    prescriptions_appliances_2$product_part1,
                                    prescriptions_costs$drug,
                                    prescriptions_incremental_costs$product,
                                    outdated_prescriptions$outdated_name))

matched_names_first = nrow(gp_prescriptions[tolower(gp_prescriptions$drug_name) %in% tolower(costed_prescription_list),])

### Begin translation of prescriptions
### Start with hand-curated names

gp_prescriptions = gp_prescriptions %>%
  mutate(drug_name = ifelse(
    tolower(drug_name) %in% tolower(prescriptions_curated_translations$alt_name),
    prescriptions_curated_translations$translation[match(tolower(drug_name), tolower(prescriptions_curated_translations$alt_name))],
    drug_name
  ))

### Then Sean Harrison names

prescriptions_AMPVMP_translations_1 = as.data.frame(prescriptions_AMPVMP_translations_1)

### Some prescription translation names are the same in both columns - odd, remove them just in case

prescriptions_AMPVMP_translations_1 = prescriptions_AMPVMP_translations_1[-c(which(prescriptions_AMPVMP_translations_1$old == prescriptions_AMPVMP_translations_1$new)),]

gp_prescriptions = gp_prescriptions %>%
  mutate(drug_name = ifelse(
    tolower(drug_name) %in% tolower(prescriptions_AMPVMP_translations_1$old),
    prescriptions_AMPVMP_translations_1$new[match(tolower(drug_name), tolower(prescriptions_AMPVMP_translations_1$old))],
    drug_name
  ))

### Remove prescriptions which cannot be costed (Exclude_1) or cannot be priced (Exclude_2)

removal_list = c(as.data.frame(prescriptions_exclude_1)[,1], as.data.frame(prescriptions_exclude_2)[,1])

gp_prescriptions = gp_prescriptions[!(tolower(gp_prescriptions$drug_name) %in% tolower(removal_list)),]

### Further fixes for AMPs and VMPs

gp_prescriptions = gp_prescriptions %>%
  mutate(drug_name = ifelse(
    tolower(drug_name) %in% tolower(prescriptions_AMPVMP_translations_2$AMP),
    prescriptions_AMPVMP_translations_2$VMP[match(tolower(drug_name), tolower(prescriptions_AMPVMP_translations_2$AMP))],
    drug_name
  ))

gp_prescriptions = gp_prescriptions %>%
  mutate(drug_name = ifelse(
    tolower(drug_name) %in% tolower(prescriptions_AMPVMP_translations_3$AMP),
    prescriptions_AMPVMP_translations_3$VMP[match(tolower(drug_name), tolower(prescriptions_AMPVMP_translations_3$AMP))],
    drug_name
  ))

### Checking the difference with how many matching prescriptions there are after all the name translations

matched_names_translated = nrow(gp_prescriptions[tolower(gp_prescriptions$drug_name) %in% tolower(costed_prescription_list),])

heading(paste0("There are now ", matched_names_translated - matched_names_first, " more matching rows than were originally in the data. Original number: ", matched_names_first, " New number: ", matched_names_translated, "\n\nThis corresponds to ", matched_names_translated/nrow(gp_prescriptions) * 100, "% of coverage."))

missing_scripts = table(gp_prescriptions$drug_name[!(tolower(gp_prescriptions$drug_name) %in% tolower(costed_prescription_list))])

heading("Name matching complete.")

print(head(sort(missing_scripts, decreasing = TRUE), n = 40))


#####################

### Obtaining true quantity of prescriptions from messy quantity column

#####################

### Find prescriptions which are multiples (e.g. "2 packs of 10 pills" or "2*10 pills")

gp_prescriptions$multiples_status = grepl(" of ", gp_prescriptions$quantity) | grepl("\\*", gp_prescriptions$quantity)

### Split this into 2 data frames

gp_prescriptions_nomulti = gp_prescriptions[gp_prescriptions$multiples_status == FALSE,]
gp_prescriptions_multi = gp_prescriptions[gp_prescriptions$multiples_status == TRUE,]

### Create simple column for extracting numbers
### First we try to turn strings into numbers

gp_prescriptions_nomulti$basic_column = as.numeric(gp_prescriptions_nomulti$quantity)

### Then we attempt to remove all letters and punctuation and then turn the resultant strings into numbers

gp_prescriptions_nomulti$simple_column = gsub("[^0-9.]", " ", gp_prescriptions_nomulti$quantity)

gp_prescriptions_nomulti$basic_column_2 = as.numeric(gp_prescriptions_nomulti$simple_column)

### If this still does not produce a final number then we replace all introduced spaces with a "_" and then use the first number

gp_prescriptions_nomulti$simple_column = gsub(" ", "_", gp_prescriptions_nomulti$simple_column)

### Outdated versions of the regex - kept for posterity

#gp_prescriptions_nomulti$simple_column = gsub("[[:alpha:]]+|[[:punct:]]", "", gp_prescriptions_nomulti$quantity)
#gp_prescriptions_nomulti$simple_column = gsub("(?<=\\d)\\s+(?=\\d)", "_", gp_prescriptions_nomulti$simple_column, perl = TRUE)

### For columns that do not have a "multiples" status, extract the first value from any numbers with a "_" as being the correct quantity
### This may not be perfect

nomulti_split = as.numeric(str_split_fixed(gp_prescriptions_nomulti$simple_column, "_", 2)[,1])

gp_prescriptions_nomulti$simple_column_nomulti = nomulti_split

### For columns that DO have a "multiples" status, split into a data frame of multiple values
### Most prescriptions follow the naming scheme of "3 packs of 10, 0.1%"
### I.e. the true value following splitting will be the first * second value and ignore the third (which will be the prescription strength)
### ~27k instances of prescriptions have more than 3 "_" in them - in these instances it is often: "56 tablet(s) - 1.25 grams + 10 micrograms (400 units)"
### I.e. not a multiple, as classified by the multiples column
### For ~3.5k multiples status is true, however, the same formula can be applied because the name is stuff like "1 pack of 120 tablets (28 oestrogen + 12 norgestrel) (3 x 40)"

### THIS SECTION WAS RUN ONCE

cat(paste0("Number of rows of previous multiple gp_scripts: \n\n"))
print(nrow(gp_prescriptions_multi))
cat(paste0("\n\n"))

extract_numbers = function(string, verbose = FALSE){

  if (grepl("\\*", string)){

    ### Extract numbers on either side of "*"

    split_string = strsplit(string, "\\*")[[1]]

    ### Extract the two multipliers

    multiplier_1 = sub(".*?(\\d+)$", "\\1", split_string[1])
    
    multiplier_2 = sub("^(\\d+).*", "\\1", split_string[2])

  } else if (grepl(" of ", string)){

    ### Extract numbers on either side of " of "
  
    split_string = strsplit(string, " of ")[[1]]

    ### Extract the two multipliers

    multiplier_1 = gsub(".*?(\\d+)(?!.*\\d).*", "\\1", split_string[1], perl = TRUE)

    multiplier_2 = sub("\\D*(\\d+).*", "\\1", split_string[2])

  }

  if (verbose == TRUE){

    cat(paste0("Input string was: ", string, "\nOutput multipliers are: ", multiplier_1, " and ", multiplier_2, "\n________\n\n"))

  }

  return(c(multiplier_1, multiplier_2))

}

if (file.exists(paste0(mainDir, "processing/intermediate_files/prescription_costs_multiple_quantity.tsv"))){


  ### OUTPUT OF CURRENT SECTION READ HERE
  ### ~16 hours to run + ~1.5k NA results in 3.3M output

  gp_prescriptions_multi = fread(paste0(mainDir, "processing/intermediate_files/prescription_costs_multiple_quantity.tsv"), data.table = FALSE)


} else {

  gp_prescriptions_multi$multi_1 = 0
  gp_prescriptions_multi$multi_2 = 0

  function_start_time = Sys.time()

  for (i in 1:nrow(gp_prescriptions_multi)){

    if ((i %% 10000) == 0){

      curr_time = Sys.time()

      function_time_taken = curr_time - function_start_time

      cat(paste0("Currently working on row: ", i, " of ", nrow(gp_prescriptions_multi), ". I.e.: ", i/nrow(gp_prescriptions_multi)*100, "% complete.\n__________\n\nTime is now: ", curr_time, "\n\nAnd total function time so far is: "))
      
      print(function_time_taken)

      cat("\n\n=======\n\n")

    }

    curr_string = gp_prescriptions_multi$quantity[i]

    if ((i %% 10000) == 0){

      out_strings = extract_numbers(curr_string, verbose = TRUE)

    } else {

      out_strings = extract_numbers(curr_string)

    }

    gp_prescriptions_multi$multi_1[i] = out_strings[1]
    gp_prescriptions_multi$multi_2[i] = out_strings[2]

  }

  fwrite(gp_prescriptions_multi, paste0(mainDir, "processing/intermediate_files/prescription_costs_multiple_quantity.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

  heading("Prescription multi-split function complete and written.")

}

### Making quantity columns numerics

gp_prescriptions_multi$multi_1 = as.numeric(gp_prescriptions_multi$multi_1)
gp_prescriptions_multi$multi_2 = as.numeric(gp_prescriptions_multi$multi_2)

### Creating true quantity column by multiplying the two multipliers

gp_prescriptions_multi = gp_prescriptions_multi[complete.cases(gp_prescriptions_multi[,c("multi_1", "multi_2")]),]
gp_prescriptions_multi$true_quantity = gp_prescriptions_multi$multi_1 * gp_prescriptions_multi$multi_2

### Creating true quantity column for no multiples
### If the number is in basic then that is used, otherwise then use basic_2 otherwise then use simplified number

gp_prescriptions_nomulti$simple_column_nomulti = as.numeric(gp_prescriptions_nomulti$simple_column_nomulti)

gp_prescriptions_nomulti$true_quantity = ifelse(!(is.na(gp_prescriptions_nomulti$basic_column)),
                                                    gp_prescriptions_nomulti$basic_column,
                                                    ifelse(!(is.na(gp_prescriptions_nomulti$basic_column_2)), gp_prescriptions_nomulti$basic_column_2, gp_prescriptions_nomulti$simple_column_nomulti))

gp_prescriptions_nomulti = gp_prescriptions_nomulti[!(is.na(gp_prescriptions_nomulti$true_quantity)),]

gp_prescription_costing = rbind(gp_prescriptions_multi[,c("eid", "issue_date", "drug_name", "true_quantity")],
                                  gp_prescriptions_nomulti[,c("eid", "issue_date", "drug_name", "true_quantity")])

heading("GP prescription data frame created.")

fwrite(gp_prescription_costing, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_data.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

gp_prescription_costing = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_data.tsv"), data.table = FALSE)

print(head(gp_prescription_costing, n = 10))

#####################

### Obtaining the ACTUAL cost of the bloody prescriptions

#####################

### Finalises to-lower status of drug names (remove all capital letters)

gp_prescription_costing$drug_lower = tolower(gp_prescription_costing$drug)

### Also for the prescription costs themselves

prescriptions_costs$drug_lower = tolower(prescriptions_costs$drug)

gp_prescription_costing_1 = gp_prescription_costing[gp_prescription_costing$drug_lower %in% prescriptions_costs$drug_lower,]

### Merge data frames based on drug column
### This takes ~10 minutes

if (file.exists(paste0(mainDir, "processing/intermediate_files/merged_prescription_costing_1.csv"))){

  merged_prescription_costing_1 = fread(paste0(mainDir, "processing/intermediate_files/merged_prescription_costing_1.csv"), data.table = FALSE)

  print(head(merged_prescription_costing_1, n = 10))

} else {

  heading("Merging prescriptions with costs data frame.")

  merged_prescription_costing_1 = merge(gp_prescription_costing_1, prescriptions_costs, by.x = "drug_lower", by.y = "drug_lower", all.x = TRUE)

  merged_prescription_costing_1$drug_name = NULL
  merged_prescription_costing_1$drug = NULL

  ### Duplicates added somewhere, removing them
  ### Some may be genuinely prescribed twice, but seems unlikely

  merged_prescription_costing_1 = merged_prescription_costing_1[-(which(duplicated(merged_prescription_costing_1))),]

  fwrite(merged_prescription_costing_1, paste0(mainDir, "processing/intermediate_files/merged_prescription_costing_1.csv"), quote = TRUE, row.names = FALSE, sep = ",")

  heading("Merge and write complete.")

  print(head(merged_prescription_costing_1, n = 10))

}

### Merged data frame used to create long data frame of costs to select out fitting cost per prescription
### Takes A HUGE AMOUNT of memory

if (file.exists(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_1.csv"))){

  gp_prescription_costing_1 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_1.csv"), data.table = FALSE)

} else {

  ### Picking out which columns of quantity/cost contain values
  ### Creates index of which quantity value is closest to true quantity value (1-4)

  quantity_cols = paste0("quantity_", 1:4)
  price_cols = paste0("basic_price_", 1:4)

  price_column_index = apply(abs(merged_prescription_costing_1[, quantity_cols] - merged_prescription_costing_1$true_quantity), 1, which.min)

  ### Calculate true_price
  ### Initialise column and then calculate for each row

  merged_prescription_costing_1 = as.data.table(merged_prescription_costing_1)

  ### Convert to long format

  long_merged_prescription_costing_1 = melt(merged_prescription_costing_1, measure.vars = grep("quantity_|basic_price_", colnames(merged_prescription_costing_1), value = T), 
                                              value.factor = FALSE, variable.factor = FALSE)


  ### Split index from variable

  long_merged_prescription_costing_1[,index := substr(variable, nchar(variable), nchar(variable))]

  ### Remove suffixes from variable

  long_merged_prescription_costing_1[,variable := ifelse(grepl("quantity_", variable), "quantity", "basic_price")]

  ### Remove missing values

  long_merged_prescription_costing_1 = long_merged_prescription_costing_1[!is.na(value)]

  print(head(long_merged_prescription_costing_1))

  ### Cast to wide format

  wide_merged_prescription_costing_1 = dcast(long_merged_prescription_costing_1, drug_lower + eid + issue_date + true_quantity + index ~ variable, value.var = "value")

  ### Find closest match

  closest_prescription_costing_1 = wide_merged_prescription_costing_1[,.SD[which.min(abs(true_quantity-quantity))],
                                                                                by = c("drug_lower", "eid", "issue_date", "true_quantity")]

  ### Create true cost column                                                                              

  closest_prescription_costing_1$true_cost = closest_prescription_costing_1$basic_price * (closest_prescription_costing_1$true_quantity / closest_prescription_costing_1$quantity)

  ### Write file and remove temporary big objects

  fwrite(closest_prescription_costing_1, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_1.csv"), quote = TRUE, row.names = FALSE, sep = ",")

  wide_merged_prescription_costing_1 = 0
  long_merged_prescription_costing_1 = 0

}

####################
####################

### Repeat above process for INCREMENTAL costs

prescriptions_incremental_costs$drug_lower = tolower(prescriptions_incremental_costs$product)

gp_prescription_costing_2 = gp_prescription_costing[gp_prescription_costing$drug_lower %in% prescriptions_incremental_costs$drug_lower,]

heading("Merging prescriptions with incremental costs data frame.")

merged_prescription_costing_2 = merge(gp_prescription_costing_2, prescriptions_incremental_costs, by.x = "drug_lower", by.y = "drug_lower", all.x = TRUE)

merged_prescription_costing_2$drug_name = NULL
merged_prescription_costing_2$product = NULL

### Duplicates added somewhere, removing them
### Some may be genuinely prescribed twice, but seems unlikely

if(any(duplicated(merged_prescription_costing_2))){

  merged_prescription_costing_2 = merged_prescription_costing_2[-(which(duplicated(merged_prescription_costing_2))),]

}

heading("Merge complete.")

print(head(merged_prescription_costing_2, n = 10))

merged_prescription_costing_2$true_cost = 0

for (i in 1:nrow(merged_prescription_costing_2)){

  currRow = merged_prescription_costing_2[i,]

  if (currRow$true_quantity == currRow$basic_quantity){

    merged_prescription_costing_2$true_cost[i] = currRow$basic_price

  } else if (currRow$true_quantity < currRow$basic_quantity){

    merged_prescription_costing_2$true_cost[i] = currRow$basic_price * (currRow$true_quantity / currRow$basic_quantity)

  } else if (currRow$true_quantity > currRow$basic_quantity){

    merged_prescription_costing_2$true_cost[i] = currRow$basic_price + ((currRow$true_quantity - currRow$basic_quantity) * currRow$price_per_unit_increase)

  }
}

fwrite(merged_prescription_costing_2, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_2.csv"), quote = TRUE, row.names = FALSE, sep = ",")

heading("True cost calculation and writing for merged_prescription_costing_2 complete.")

####################
####################

### Repeat above process for OUTDATED costs

outdated_prescriptions$drug_lower = tolower(outdated_prescriptions$outdated_name)

gp_prescription_costing_3 = gp_prescription_costing[gp_prescription_costing$drug_lower %in% outdated_prescriptions$drug_lower,]

heading("Merging prescriptions with incremental costs data frame.")

merged_prescription_costing_3 = merge(gp_prescription_costing_3, outdated_prescriptions, by.x = "drug_lower", by.y = "drug_lower", all.x = TRUE)

merged_prescription_costing_3$drug_name = NULL
merged_prescription_costing_3$outdated_name = NULL
merged_prescription_costing_3$year = NULL

### Duplicates added somewhere, removing them
### Some may be genuinely prescribed twice, but seems unlikely

if(any(duplicated(merged_prescription_costing_3))){

  merged_prescription_costing_3 = merged_prescription_costing_3[-(which(duplicated(merged_prescription_costing_3))),]

}

heading("Merge complete.")

print(head(merged_prescription_costing_3, n = 10))

merged_prescription_costing_3 = as.data.table(merged_prescription_costing_3)

### Convert to long format

long_merged_prescription_costing_3 = melt(merged_prescription_costing_3, measure.vars = grep("quantity_|basic_price_", colnames(merged_prescription_costing_3), value = T), 
                                            value.factor = FALSE, variable.factor = FALSE)


### Split index from variable

long_merged_prescription_costing_3[,index := substr(variable, nchar(variable), nchar(variable))]

### Remove suffixes from variable

long_merged_prescription_costing_3[,variable := ifelse(grepl("quantity_", variable), "quantity", "basic_price")]

### Remove missing values

long_merged_prescription_costing_3 = long_merged_prescription_costing_3[!is.na(value)]

print(long_merged_prescription_costing_3)

### Cast to wide format

wide_merged_prescription_costing_3 = dcast(long_merged_prescription_costing_3, drug_lower + eid + issue_date + true_quantity + index + cpi_index_to_2022 ~ variable, value.var = "value")

### Find closest match

closest_prescription_costing_3 = wide_merged_prescription_costing_3[,.SD[which.min(abs(true_quantity-quantity))],
                                                                              by = c("drug_lower", "eid", "issue_date", "true_quantity")]

### Create true cost column                                                                              

closest_prescription_costing_3$true_cost = (closest_prescription_costing_3$basic_price * (closest_prescription_costing_3$true_quantity / closest_prescription_costing_3$quantity)) * closest_prescription_costing_3$cpi_index_to_2022

### Write file and remove temporary big objects

fwrite(closest_prescription_costing_3, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_3.csv"), quote = TRUE, row.names = FALSE, sep = ",")

wide_merged_prescription_costing_3 = 0
long_merged_prescription_costing_3 = 0

heading("True cost calculation and writing for closest_prescription_costing_3 complete.")

####################
####################

### Repeat above process for APPLIANCES ONE costs

prescriptions_appliances_1$drug_lower = tolower(prescriptions_appliances_1[,1])

gp_prescription_costing_4 = gp_prescription_costing[gp_prescription_costing$drug_lower %in% prescriptions_appliances_1$drug_lower,]

heading("Merging prescriptions with incremental costs data frame.")

merged_prescription_costing_4 = merge(gp_prescription_costing_4, prescriptions_appliances_1[,c(2,3,5)], by.x = "drug_lower", by.y = "drug_lower", all.x = TRUE)

merged_prescription_costing_4$drug_name = NULL

### Duplicates added somewhere, removing them
### Some may be genuinely prescribed twice, but seems unlikely

if(any(duplicated(merged_prescription_costing_4))){

  merged_prescription_costing_4 = merged_prescription_costing_4[-(which(duplicated(merged_prescription_costing_4))),]

}

heading("Merge complete.")

print(head(merged_prescription_costing_4, n = 10))

merged_prescription_costing_4$quantity = as.numeric(gsub("g", "", merged_prescription_costing_4$quantity))

merged_prescription_costing_4$true_cost = merged_prescription_costing_4$cost * (merged_prescription_costing_4$true_quantity / merged_prescription_costing_4$quantity)

### Write file

fwrite(merged_prescription_costing_4, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_4.csv"), quote = TRUE, row.names = FALSE, sep = ",")

heading("True cost calculation and writing for merged_prescription_costing_4 complete.")

####################
####################

### Repeat above process for APPLIANCES TWO costs

prescriptions_appliances_2$part1_lower = tolower(prescriptions_appliances_2$product_part1)
prescriptions_appliances_2$part2_lower = tolower(prescriptions_appliances_2$product_part2)

### Reduce gp_prescription_costing to gp_prescription_costing_refined to reduce computational time

appliance_match_term = unique(prescriptions_appliances_2$part1_lower)[1]

appliance_match_vector = grepl(appliance_match_term, gp_prescription_costing$drug_lower)

for (i in 2:length(unique(prescriptions_appliances_2$part1_lower))){

  currterm = unique(prescriptions_appliances_2$part1_lower)[i]

  cat(paste0("Currently working on row: ", i, " of ", length(unique(prescriptions_appliances_2$part1_lower)), ". I.e.: ", currterm, "\n\n"))

  curr_match_vector = grepl(currterm, gp_prescription_costing$drug_lower)

  appliance_match_vector = appliance_match_vector | curr_match_vector

}

### Filter the data frame based on the appliance_match_vector

gp_prescription_costing_refined = gp_prescription_costing[appliance_match_vector, ]
gp_prescription_costing_refined = gp_prescription_costing_refined[-(which(duplicated(gp_prescription_costing_refined))),]

merged_prescription_costing_5 = as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(merged_prescription_costing_5) = c("drug_lower", "eid", "issue_date", "true_quantity", "index", "basic_price", "quantity")

for (i in 1:nrow(prescriptions_appliances_2)){

  cat(paste0("Currently working on row: ", i, " of ", nrow(prescriptions_appliances_2), ". I.e.: \n\n"))
  print(prescriptions_appliances_2[i,])
  cat("\n\n")

  if (prescriptions_appliances_2$part2_lower[i] == ""){

    matched_rows = gp_prescription_costing_refined[gp_prescription_costing_refined$drug_lower %in% prescriptions_appliances_2$part1_lower[i],]

  } else {

    matched_rows = gp_prescription_costing_refined[grepl(prescriptions_appliances_2$part1_lower[i], gp_prescription_costing_refined$drug_lower) & grepl(prescriptions_appliances_2$part2_lower[i], gp_prescription_costing_refined$drug_lower),]

  }

  if (nrow(matched_rows) == 0){

    next

  }

  ### Calculate true_price
  ### Initialise column and then calculate for each row

  matched_rows$quantity_1 = prescriptions_appliances_2[i,"quantity_1"]
  matched_rows$quantity_2 = prescriptions_appliances_2[i,"quantity_2"]
  matched_rows$quantity_3 = prescriptions_appliances_2[i,"quantity_3"]
  matched_rows$quantity_4 = prescriptions_appliances_2[i,"quantity_4"]

  matched_rows$basic_price_1 = prescriptions_appliances_2[i,"basic_price_1"]
  matched_rows$basic_price_2 = prescriptions_appliances_2[i,"basic_price_2"]
  matched_rows$basic_price_3 = prescriptions_appliances_2[i,"basic_price_3"]
  matched_rows$basic_price_4 = prescriptions_appliances_2[i,"basic_price_4"]

  merged_prescription_costing_appli = as.data.table(matched_rows)

  ### Convert to long format

  long_merged_prescription_costing_appli = melt(merged_prescription_costing_appli, measure.vars = grep("quantity_|basic_price_", colnames(merged_prescription_costing_appli), value = T), 
                                              value.factor = FALSE, variable.factor = FALSE)


  ### Split index from variable

  long_merged_prescription_costing_appli[,index := substr(variable, nchar(variable), nchar(variable))]

  ### Remove suffixes from variable

  long_merged_prescription_costing_appli[,variable := ifelse(grepl("quantity_", variable), "quantity", "basic_price")]

  ### Remove missing values

  long_merged_prescription_costing_appli = long_merged_prescription_costing_appli[!is.na(value)]

  print(long_merged_prescription_costing_appli)

  ### Cast to wide format

  wide_merged_prescription_costing_appli = dcast(long_merged_prescription_costing_appli, drug_lower + eid + issue_date + true_quantity + index ~ variable, value.var = "value")

  ### Find closest match

  closest_prescription_costing_appli = wide_merged_prescription_costing_appli[,.SD[which.min(abs(true_quantity-quantity))],
                                                                                by = c("drug_lower", "eid", "issue_date", "true_quantity")]

  merged_prescription_costing_5 = rbind(merged_prescription_costing_5, closest_prescription_costing_appli)

}

### For some reason "hypromellose 0.5% eye drops" was not aggregating correctly - not apparent why

merged_prescription_costing_5$basic_price[merged_prescription_costing_5$drug_lower == "hypromellose 0.5% eye drops"] = 99
merged_prescription_costing_5$quantity[merged_prescription_costing_5$drug_lower == "hypromellose 0.5% eye drops"] = 10

### Create true cost column                                                                              

merged_prescription_costing_5$true_cost = merged_prescription_costing_5$basic_price * (merged_prescription_costing_5$true_quantity / merged_prescription_costing_5$quantity)

fwrite(merged_prescription_costing_5, paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_5.csv"), quote = TRUE, row.names = FALSE, sep = ",")

heading("True cost calculation and writing for merged_prescription_costing_5 complete.")


#####################

### Five prescription cost files made, read them back in and create final prescription costs

#####################

gp_prescription_costing_1 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_1.csv"), data.table = FALSE, select = c("eid", "true_cost"))
gp_prescription_costing_2 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_2.csv"), data.table = FALSE, select = c("eid", "true_cost"))
gp_prescription_costing_3 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_3.csv"), data.table = FALSE, select = c("eid", "true_cost"))
gp_prescription_costing_4 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_4.csv"), data.table = FALSE, select = c("eid", "true_cost"))
gp_prescription_costing_5 = fread(paste0(mainDir, "processing/intermediate_files/gp_prescription_costing_5.csv"), data.table = FALSE, select = c("eid", "true_cost"))

gp_prescription_costing_total = rbind(gp_prescription_costing_1, gp_prescription_costing_2, gp_prescription_costing_3, gp_prescription_costing_4, gp_prescription_costing_5)

### Drug costs are in penny format - convert to pounds

gp_prescription_costing_total$true_cost = gp_prescription_costing_total$true_cost/100

### Cost of providing a prescription is £1.27, add this to each instance

gp_prescription_costing_total$true_cost = gp_prescription_costing_total$true_cost + 1.27

final_individual_prescription_costs = aggregate(true_cost ~ eid, data = gp_prescription_costing_total, FUN = sum)

colnames(final_individual_prescription_costs)[2] = "prescription_costs"

### Get final date for prescriptions
### SUPERCEDED

#gp_prescriptions_dates = gp_prescriptions[c("eid", "issue_date")]
#
#most_recent_gp_prescriptions_dates = gp_prescriptions_dates %>%
#                                      group_by(eid) %>%
#                                        summarize(issue_date = max(issue_date))
#
#colnames(most_recent_gp_prescriptions_dates)[2] = "prescription_end_of_followup"
#
#most_recent_gp_prescriptions_dates = most_recent_gp_prescriptions_dates[most_recent_gp_prescriptions_dates$eid %in% final_individual_prescription_costs$eid,]
#
#final_individual_prescription_costs = merge(final_individual_prescription_costs, most_recent_gp_prescriptions_dates, by = "eid", all = TRUE)

### Get first date for prescriptions
### SUPERCEDED

#oldest_gp_prescriptions_dates = gp_prescriptions_dates %>%
#                                      group_by(eid) %>%
#                                        summarize(issue_date = min(issue_date))
#
#colnames(oldest_gp_prescriptions_dates)[2] = "prescription_start_of_followup"
#
#oldest_gp_prescriptions_dates = oldest_gp_prescriptions_dates[oldest_gp_prescriptions_dates$eid %in% final_individual_prescription_costs$eid,]
#
#final_individual_prescription_costs = merge(final_individual_prescription_costs, oldest_gp_prescriptions_dates, by = "eid", all = TRUE)

### Merging everything

final_individual_prescription_costs = final_individual_prescription_costs[final_individual_prescription_costs$eid %in% individual_all_costs$eid,]

heading("Data frames pre-rbind.")

print(head(individual_all_costs))

individual_all_costs = merge(individual_all_costs, final_individual_prescription_costs, by = "eid", all = TRUE)

individual_all_costs$prescription_end_of_followup = maxPrescriptDate

individual_all_costs = merge(individual_all_costs, prescript_start, by = "eid", all.x = TRUE)

heading("Added prescription costs to final cost frame.")

print(head(individual_all_costs, n = 10))


#####################

### Get days of follow-up

#####################

### MOST OF THIS WAS SUPERCEDED WITH BY JUST INCORPORATING DATE OF JOINING UKB AS START OF FOLLOW-UP
##################
#years_of_followup = data.frame(eid = individual_all_costs$eid)
#
#ukb_phenos = ukb_phenos[ukb_phenos$eid %in% years_of_followup$eid,]
#
#years_of_followup = merge(years_of_followup, ukb_phenos, by = "eid")
#
### Start of follow-up is EITHER when the individual was admitted to the UKB or the start date (1st April 2006 - 2006-06-01)
#
#years_of_followup$initial_date[years_of_followup$initial_date < start_date] = start_date
#
### Adding death dates for individuals who have died (I.e. end of follow-up)
##################

death_reads$date_of_death = as.Date(death_reads$date_of_death, format = "%d/%m/%Y")

death_dates = death_reads[,c("eid", "date_of_death")]
death_dates = unique(death_dates[death_dates$eid %in% individual_all_costs$eid,])

### Adding date of death to the follow-up

individual_all_costs = merge(individual_all_costs, death_dates, by = "eid", all = TRUE)

#years_of_followup$final_date = ifelse(is.na(years_of_followup$final_date) == TRUE, as.Date(end_date), years_of_followup$final_date)

### This next line is needed for some reason

#years_of_followup$final_date = as.Date(years_of_followup$final_date, origin = "1970-01-01")

#years_of_followup$days_of_followup = as.numeric(years_of_followup$final_date - years_of_followup$initial_date)

#individual_all_costs = merge(individual_all_costs, years_of_followup, by = "eid")

heading("Added death dates to final cost frame.")

fwrite(individual_all_costs, paste0(mainDir, "outputs/cost_phenotypes/ALL_cost_phenotypes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

print(head(individual_all_costs, n = 10))


#####################

### SCRIPT COMPLETE

#####################


heading("HE'S ONLY GONE AND BLOWN THE BLOODY DOORS OFF!")

heading("Script complete.")
