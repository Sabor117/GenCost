
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

nhs_care_costs_file_2022 = paste0(mainDir, "processing/cost_coding_values/2_National_schedule_of_NHS_costs_FY21-22_v2.xlsx")
nhs_care_costs_file_2017 = paste0(mainDir, "processing/cost_coding_values/2_-_National_schedule_of_reference_costs_-_the_main_schedule.xlsx")
prescriptions_costs_file = paste0(mainDir, "processing/cost_coding_values/prescription_costing/GenCOST_prescriptions_sheet.csv")

### Records files

hes_data_file = "/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin.txt"

### Sense check file

sense_check_file = paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_grouper_BB.csv")

### File reading

grouper_output = readLines(grouper_output_file)
grouper_classpat_recode_list = fread(grouper_classpat_recode_list_file, data.table = FALSE)
care_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "APC")
consultant_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "OP")
unbun_costs_2022 = read_excel(nhs_care_costs_file_2022, sheet = "Total HRGs")
elective_excess_costs_2017 = read_excel(nhs_care_costs_file_2017, sheet = "EL_XS")
non_elective_excess_costs_2017 = read_excel(nhs_care_costs_file_2017, sheet = "NEL_XS")
ukb_phenos = fread(ukb_phenos_file, data.table = FALSE, select = c("eid", "53-0.0", "34-0.0"))
sense_check = fread(sense_check_file, data.table = FALSE)
hes_data = fread(hes_data_file, data.table = FALSE, select = c("eid", "ins_index", "epistat", "epistart", "epiend"))

### Global variables

cpi_index_2017_to_2022 = 1.18 # Value extracted from: https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator

### Start date defined by Padraig

start_date = as.Date("01/06/2006", format = "%d/%m/%Y")

### End date defined as last date in death register

end_date = as.Date("19/12/2022", format = "%d/%m/%Y")


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


#####################

### Data comparison

#####################

grouper_output_test = grouper_output

grouper_output_test$iid = str_split_fixed(grouper_output_test$PROVSPNO, "_", 2)[,1]

colnames(ukb_phenos) = c("iid", "date_of_admission", "year_of_birth")

grouper_output_test = merge(grouper_output_test, ukb_phenos[,c("iid", "year_of_birth")], by = "iid")

colnames(grouper_output_test) = tolower(colnames(grouper_output_test))

colnames(grouper_output_test) = gsub("diag_0", "diag_", colnames(grouper_output_test))
colnames(grouper_output_test) = gsub("oper_0", "oper_", colnames(grouper_output_test))

grouper_output_test = grouper_output_test[grouper_output_test$year_of_birth %in% unique(sense_check$year_of_birth),]
grouper_output_test = grouper_output_test[grouper_output_test$diag_1 %in% unique(sense_check$diag_1),]

sense_check_refined = sense_check[,c("epiorder",
                                        "startage",
                                        "sex",
                                        "classpat",
                                        "admisorc",
                                        "admimeth",
                                        "disdest",
                                        "dismeth",
                                        "epidur",
                                        "mainspef",
                                        "neocare",
                                        "tretspef",
                                        paste0("diag_", 1:20),
                                        paste0("oper_", 1:24),
                                        "year_of_birth",
                                        "spellhrg",
                                        "fce_pbc")]

check_frame = data.frame(matrix(ncol = ncol(grouper_output_test), nrow = 0))

colnames(check_frame) = colnames(grouper_output_test)

for (i in 37:nrow(sense_check_refined)){

    currRow = sense_check_refined[i,]

    cat(paste0("\nCurrent row being examined is: \n\n"))
    print(currRow)
    cat(paste0("\n================\n\n"))

    check_one = grouper_output_test[grouper_output_test$diag_1 == currRow$diag_1 &
                                        grouper_output_test$oper_1 == currRow$oper_1 &
                                        grouper_output_test$year_of_birth == currRow$year_of_birth,]

    cat(paste0("\nFirst check findings: \n\n"))
    print(check_one)
    cat(paste0("\n================\n\n"))

    if (nrow(check_one) > 1){

        check_two = check_one[check_one$diag_2 == currRow$diag_2 &
                                check_one$oper_2 == currRow$oper_2,]

        if (nrow(check_two) > 1){

            check_three = check_two[check_two$diag_3 == currRow$diag_3 &
                                check_two$oper_3 == currRow$oper_3,]

            if (nrow(check_three) > 1){

                check_four = check_three[check_three$diag_4 == currRow$diag_4 &
                                check_three$diag_5 == currRow$diag_5 &
                                check_three$diag_6 == currRow$diag_6 &
                                check_three$oper_4 == currRow$oper_4 &
                                check_three$oper_5 == currRow$oper_5 &
                                check_three$oper_6 == currRow$oper_6 &
                                check_three$classpat == currRow$classpat &
                                check_three$tretspef == currRow$tretspef &
                                check_three$admisorc == currRow$admisorc &
                                check_three$epidur == currRow$epidur,]


                outRow = check_four

                if (nrow(check_four) > 1){

                    break

                }
            } else {

                outRow = check_three

            }
        } else {

            outRow = check_two

        }
    } else {

        outRow = check_one

    }

    cat(paste0("\nFinal findings: \n\n"))
    print(outRow)
    cat(paste0("\n================\n\n"))

    check_frame = rbind(check_frame, outRow)

}

fwrite(check_frame, paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_grouper_output_smw.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

check_frame = fread(paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_grouper_output_smw.tsv"), data.table = FALSE)

output_comparison = data.frame(iid = check_frame$provspno,
                                padraig_spellhrg = sense_check_refined$spellhrg,
                                padraig_fce_pbc = sense_check_refined$fce_pbc,
                                seb_spellhrg = check_frame$spellhrg,
                                seb_fce_pbc = check_frame$fce_pbc)

fwrite(output_comparison, paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_output_code_comparison.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

grouper_output = grouper_output[grouper_output$PROVSPNO %in% check_frame$provspno,]


#####################

### Data parsing

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

excess_el_care_costs$excess_bed_days = as.numeric(excess_el_care_costs$excess_bed_days)
excess_el_care_costs$excess_bed_days[is.na(excess_el_care_costs$excess_bed_days)] = 0

excess_el_care_costs$national_avg_unit_cost = as.numeric(excess_el_care_costs$national_avg_unit_cost)
excess_el_care_costs$national_avg_unit_cost[is.na(excess_el_care_costs$national_avg_unit_cost)] = 0

excess_nel_care_costs$excess_bed_days = as.numeric(excess_nel_care_costs$excess_bed_days)
excess_nel_care_costs$excess_bed_days[is.na(excess_nel_care_costs$excess_bed_days)] = 0

excess_nel_care_costs$national_avg_unit_cost = as.numeric(excess_nel_care_costs$national_avg_unit_cost)
excess_nel_care_costs$national_avg_unit_cost[is.na(excess_nel_care_costs$national_avg_unit_cost)] = 0

heading("Excess costs for EL-XS/NEL-XS created..")

print(head(excess_nel_care_costs, n = 10))

#####################

### Obtaining final APC cost

#####################

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

final_cost_frame = merge(final_cost_frame, costing_frame, by = "iid", all = TRUE)
colnames(final_cost_frame)[3:4] = c("apc_cost", "weighted_apc_cost")

final_cost_frame = merge(final_cost_frame, regular_care_costing_frame[,c(2,4,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[5:6] = c("rp_cost", "weighted_rp_cost")

final_cost_frame = merge(final_cost_frame, day_care_costing_frame[,c(2,4,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[7:8] = c("dc_cost", "weighted_dc_cost")

final_cost_frame = merge(final_cost_frame, maternity_care_costing_frame[,c(1,4)], by = "iid", all = TRUE)
colnames(final_cost_frame)[9] = c("mat_cost")

final_cost_frame = merge(final_cost_frame, unbun_costing_frame[,c(1,5)], by = "iid", all = TRUE)
colnames(final_cost_frame)[10] = c("unbun_cost")

### Remove invalid Grouper code from output

final_cost_frame[final_cost_frame$currency_code == "UZ01Z", c(3,4)] = 0

### Sum of all costs and sum of weighted costs

final_cost_frame$hes_costs = rowSums(final_cost_frame[,c(3,5,7,9,10)], na.rm = TRUE)
final_cost_frame$weighted_hes_costs = rowSums(final_cost_frame[,c(4,6,8:10)], na.rm = TRUE)

### Adding the codes for Padraig

grouper_codes = grouper_output[,c("PROVSPNO", "FCE_HRG")]
colnames(grouper_codes)[1] = "iid"

final_cost_frame = merge(final_cost_frame, grouper_codes, by = "iid")

final_cost_frame$currency_code = final_cost_frame$FCE_HRG

final_cost_frame$FCE_HRG = NULL

fwrite(final_cost_frame, paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_output_instance_values.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

### Finally create sum of each hospital instance for each individual

final_cost_frame$eid = str_split_fixed(final_cost_frame$iid, "_", 2)[,1]

### Making the final phenotype data frame: final_costs_individual

final_costs_individual = aggregate(cbind(hes_costs, weighted_hes_costs) ~ eid, data = final_cost_frame, FUN = sum)

heading("Final cost fame for each individual being parsed.")

print(head(final_costs_individual, n = 10))

fwrite(final_costs_individual, paste0(mainDir, "processing/nhs_grouper/sense_check/sense_check_output_individual_values.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
