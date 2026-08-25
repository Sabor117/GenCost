### ---------------------------
###
### Script name: meta_step15_1_gwas_cat_info.R
###
### Purpose of script: Create GWAS Catalog upload metadata table for meta-analysis summary statistics
###
### Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2026-06-09
###
### ---------------------------
###

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(stringr)
library(dplyr)
library(purrr)

options(scipen = 999)

### Start timer

script_start = Sys.time()

### Classic memes

heading = function(text){
  
  cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))
  
}

### Helper function for compact metadata strings

collapse_unique = function(x){
  
  x |>
    unique() |>
    sort() |>
    paste(collapse = "|")
  
}

### Helper function to extract maximum sample size from a cohort file

get_max_n = function(file){
  
  curr_data = fread(file, select = "n", showProgress = FALSE)
  
  if (nrow(curr_data) == 0 || all(is.na(curr_data$n))){
    return(NA_real_)
  }
  
  max(curr_data$n, na.rm = TRUE)
  
}


##### =========================== #####

### Directories and input files

##### =========================== #####

heading("Setup directories and input files")

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "processing/meta_sumstats/")
gwasDir = paste0(mainDir, "outputs/METAL_v4_upload/")

array_table_file = paste0(mainDir, "processing/misc_data/cohort_array_table.txt")

allMetas = Sys.glob(paste0(metaDir, "*"))

### Remove ALL-level meta folders/files if present

allMetas = allMetas[!(grepl("/ALL_", allMetas))]

cat("Number of meta-analysis folders/files found:\n")
print(length(allMetas))

cat("\nFirst few meta-analysis names:\n")
print(head(basename(allMetas)))


##### =========================== #####

### Read supporting metadata

##### =========================== #####

heading("Read cohort array table")

array_table = fread(array_table_file, data.table = FALSE)

array_table$Study = toupper(array_table$Study)

cat("Dimensions of array table:\n")
print(dim(array_table))

cat("\nFirst few studies in array table:\n")
print(head(array_table$Study))


##### =========================== #####

### Define fixed metadata

##### =========================== #####

heading("Define fixed metadata")

country_map = c(
  FINNGEN = "Finland",
  UKB = "UK",
  GNH = "England",
  GE = "England",
  GS20K = "Scotland",
  ESTBB = "Estonia",
  CHB = "Denmark",
  QGP = "Qatar",
  MGBB = "USA",
  GBP = "Australia",
  AGDS = "Australia"
)

ancestry_description = c(
  EUR = "EUR - multiple cohorts containing individuals of white European descent (UK, Finland, Estonia, Australia, etc, based on cohorts present). Sample size: ",
  CSA = "CSA - Central/South Asian ancestry from two possible cohorts (UK Biobank individuals, defined by PanUKBB methodology, and Genes and Health). Sample size: ",
  AFR = "AFR - Subset of UK Biobank individuals of African descent as defined by PanUKBB methodology. Sample size: ",
  EAS = "EAS - Subset of UK Biobank individuals of East Asian ancestry descent as defined by PanUKBB methodology. Sample size: ",
  MID = "MID - Individuals of Middle Eastern descent from Qatar Biobank. Sample size: "
)


##### =========================== #####

### Prepare output table

##### =========================== #####

heading("Prepare output table")

outtable = data.frame(matrix(ncol = 16, nrow = 0))

colnames(outtable) = c("study",
                       "genotyping",
                       "manufacturer",
                       "software",
                       "imputation",
                       "variant_n",
                       "covars",
                       "n",
                       "trait",
                       "sumstats",
                       "md5",
                       "cohorts",
                       "sex_strat",
                       "description",
                       "ancestry",
                       "country")


##### =========================== #####

### Run metadata extraction

##### =========================== #####

heading("Run metadata extraction")

for (i in seq_along(allMetas)){
  
  loop_start = Sys.time()
  
  currMeta = allMetas[i]
  meta_name = basename(currMeta)
  
  heading(paste0("Working on meta ", i, " / ", length(allMetas), ": ", meta_name))
  
  ##### --------------------------- #####
  ### Step 1. Identify cohort files
  ##### --------------------------- #####
  
  all_cohorts = Sys.glob(paste0(currMeta, "/*txt.gz"))
  
  cat("Number of cohort summary statistic files found:\n")
  print(length(all_cohorts))
  
  if (length(all_cohorts) == 0){
    
    warning(paste0("No cohort files found for ", meta_name, ". Skipping this meta-analysis."))
    
    loop_end = Sys.time()
    
    cat("\nTime spent on skipped meta-analysis:\n")
    print(round(difftime(loop_end, loop_start, units = "mins"), 2))
    
    next
    
  }
  
  cohorts = basename(all_cohorts) |>
    str_extract("^[^.]+") |>
    toupper()
  
  cat("\nCohorts present:\n")
  print(sort(unique(cohorts)))
  
  cohortString = cohorts |>
    unique() |>
    sort() |>
    paste(collapse = "|")
  
  
  ##### --------------------------- #####
  ### Step 2. Match cohorts to metadata table
  ##### --------------------------- #####
  
  presentStudies = array_table[array_table$Study %in% cohorts, ]
  
  missingStudies = setdiff(unique(cohorts), presentStudies$Study)
  
  cat("\nNumber of cohorts matched in array table:\n")
  print(paste0(nrow(presentStudies), " / ", length(unique(cohorts))))
  
  if (length(missingStudies) > 0){
    
    cat("\nWARNING: These cohorts were not found in the array table:\n")
    print(missingStudies)
    
  }
  
  assocString = presentStudies |>
    pull(SoftwareCollapsed) |>
    collapse_unique()
  
  arrayString = presentStudies |>
    pull(ArrayCollapsed) |>
    collapse_unique()
  
  imputationString = presentStudies |>
    pull(ImputationCollapsed) |>
    collapse_unique()
  
  
  ##### --------------------------- #####
  ### Step 3. Extract sample sizes
  ##### --------------------------- #####
  
  cat("\nReading files for sample sizes.\n")
  
  sample_size_start = Sys.time()
  
  cohort_info = tibble(file = all_cohorts,
                       cohort = cohorts,
                       ancestry = str_split_fixed(basename(all_cohorts), "\\.", 8)[, 6]) |>
    mutate(max_n = map_dbl(file, get_max_n))
  
  total_n = sum(cohort_info$max_n, na.rm = TRUE)
  
  cat("\nSample size summary by ancestry:\n")
  print(cohort_info |>
          group_by(ancestry) |>
          summarise(n_files = n(),
                    total_max_n = sum(max_n, na.rm = TRUE),
                    .groups = "drop"))
  
  cat("\nTotal sample size across cohort files:\n")
  print(total_n)
  
  sample_size_end = Sys.time()
  
  cat("\nTime spent extracting sample sizes:\n")
  print(round(difftime(sample_size_end, sample_size_start, units = "mins"), 2))
  
  
  ##### --------------------------- #####
  ### Step 4. Count variants in GWAS Catalog formatted file
  ##### --------------------------- #####
  
  variant_file_path = paste0(gwasDir, meta_name, "_gwas_cat_format.tsv.gz")
  
  cat("\nReading GWAS Catalog formatted file to count variants:\n")
  print(variant_file_path)
  
  if (!file.exists(variant_file_path)){
    
    warning(paste0("GWAS Catalog formatted file not found for ", meta_name))
    variant_n = NA_integer_
    
  } else {
    
    variant_file = fread(variant_file_path,
                         select = "variant_id",
                         showProgress = FALSE)
    
    variant_n = nrow(variant_file)
    
  }
  
  cat("\nNumber of variants:\n")
  print(variant_n)
  
  
  ##### --------------------------- #####
  ### Step 5. Build ancestry description
  ##### --------------------------- #####
  
  ancestry_present = unique(cohort_info$ancestry)
  
  cat("\nAncestries present:\n")
  print(ancestry_present)
  
  ancestry_string = ""
  
  for (anc in names(ancestry_description)){
    
    if (anc %in% ancestry_present){
      
      anc_n = sum(cohort_info$max_n[cohort_info$ancestry == anc], na.rm = TRUE)
      
      ancestry_string = paste0(ancestry_string,
                               ancestry_description[[anc]],
                               anc_n,
                               "|")
      
    }
    
  }
  
  unknown_ancestries = setdiff(ancestry_present, names(ancestry_description))
  
  if (length(unknown_ancestries) > 0){
    
    cat("\nWARNING: These ancestries do not have a predefined description:\n")
    print(unknown_ancestries)
    
  }
  
  
  ##### --------------------------- #####
  ### Step 6. Build country string
  ##### --------------------------- #####
  
  countries = sort(unique(unname(country_map[names(country_map) %in% cohorts])))
  
  country_string = paste(countries, collapse = "|")
  
  cat("\nCountries represented:\n")
  print(country_string)
  
  
  ##### --------------------------- #####
  ### Step 7. Define covariates and stratification labels
  ##### --------------------------- #####
  
  sex_strat = ifelse(grepl("_F", meta_name) | grepl("_M", meta_name), TRUE, FALSE)
  
  if (isTRUE(sex_strat)){
    
    covars = "Age-at-end-of-follow-up|Age^2|PC1-20|Study-specific-covariates"
    
    sex_type = ifelse(grepl("_F", meta_name), "F", "M")
    
  } else {
    
    covars = "Age-at-end-of-follow-up|Sex|Age*Sex|Age^2|PC1-20|Study-specific-covariates"
    
    sex_type = "Combined"
    
  }
  
  overall_strat = case_when(
    grepl("_ALL", meta_name) ~ "",
    grepl("_F", meta_name) ~ "Female",
    grepl("_M", meta_name) ~ "Male",
    grepl("_5_18", meta_name) ~ "5-18 years old",
    grepl("_19_35", meta_name) ~ "19-35 years old",
    grepl("_36_55", meta_name) ~ "36-55 years old",
    grepl("_56_75", meta_name) ~ "56-75 years old",
    grepl("_76_95", meta_name) ~ "76+ years old",
    .default = ""
  )
  
  overall_trait = case_when(
    grepl("IN_", meta_name) ~ "Inpatient cost",
    grepl("DRUG_", meta_name) ~ "Prescription drug cost",
    grepl("INOUT_", meta_name) ~ "Inpatient + outpatient cost",
    grepl("PRIM_", meta_name) ~ "Primary care cost",
    .default = "Cost"
  )
  
  cat("\nTrait and stratification summary:\n")
  print(data.frame(trait = overall_trait,
                   sex_strat = sex_type,
                   description = overall_strat))
  
  
  ##### --------------------------- #####
  ### Step 8. Add row to output table
  ##### --------------------------- #####
  
  sumstatfile = paste0(meta_name, "_gwas_cat_format.tsv.gz")
  md5file = paste0(meta_name, "_gwas_cat_format.md5")
  
  outrow = data.frame(study = meta_name,
                      genotyping = arrayString,
                      manufacturer = "Thermo Fisher Scientific (Applied Biosystems)|Illumina|deCODE genetics",
                      software = assocString,
                      imputation = imputationString,
                      variant_n = variant_n,
                      covars = covars,
                      n = total_n,
                      trait = overall_trait,
                      sumstats = sumstatfile,
                      md5 = md5file,
                      cohorts = cohortString,
                      sex_strat = sex_type,
                      description = overall_strat,
                      ancestry = ancestry_string,
                      country = country_string)
  
  outtable = rbind(outtable, outrow)
  
  cat("\nCurrent number of rows in output table:\n")
  print(nrow(outtable))
  
  
  ##### --------------------------- #####
  ### Step 9. Print loop timing
  ##### --------------------------- #####
  
  loop_end = Sys.time()
  
  cat("\nFinished meta-analysis:\n")
  print(meta_name)
  
  cat("\nTime spent on this meta-analysis:\n")
  print(round(difftime(loop_end, loop_start, units = "mins"), 2))
  
}


##### =========================== #####

### Write output

##### =========================== #####

heading("Write output table")

outfile = paste0(mainDir, "processing/misc_data/gwas_cat_upload_info.tsv")

fwrite(outtable,
       outfile,
       quote = FALSE,
       row.names = FALSE,
       sep = "\t",
       na = "NA")

cat("Output written to:\n")
print(outfile)

cat("\nFinal dimensions of output table:\n")
print(dim(outtable))


##### =========================== #####

### Finish script

##### =========================== #####

heading("Finished script")

script_end = Sys.time()

cat("Total script runtime:\n")
print(round(difftime(script_end, script_start, units = "mins"), 2))
