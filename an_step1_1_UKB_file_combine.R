##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)
library(R.utils)

options(scipen = 999) # prevent scientific notation in output

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gwasDir = paste0(mainDir, "outputs/regenie/step2_output/")

all_gwas = Sys.glob(paste0(gwasDir, "*"))

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. PathWAS incomplete. Nae args.")
  
}

### Test example:
### args = c("1")

run_no = as.numeric(args[1])

currGWAS = all_gwas[run_no]

cat(paste0("\nCurrent run# = ", run_no, "\n\n"))
cat(paste0("\nCorresponds to ", currGWAS, "\n\n==================\n\n"))

all_gwas_files = Sys.glob(paste0(currGWAS, "/*.regenie"))

if (length(all_gwas_files) == 0){

    cat(paste0("\nGWAS has no valid files for combining. CHECK 3 COMPLETE.\n\n==================\n\n"))

    stop()

}

phenotypes = basename(all_gwas_files)
phenotypes = gsub(".regenie", "", phenotypes)
phenotypes_split = str_split_fixed(phenotypes, "_log_", 2)
phenotypes = unique(phenotypes_split[,2])

population = strsplit(basename(currGWAS), "_")[[1]][2]

for (i in 1:length(phenotypes)){

    currPheno = phenotypes[i]

    cat(paste0("\nStarting run on ", currPheno, ". CHECK 1 COMPLETE.\n\n==================\n\n"))

    curr_file_set = all_gwas_files[grepl(currPheno, all_gwas_files)]

    curr_output = fread(curr_file_set[1], data.table = FALSE)

    total_snps = nrow(curr_output)

    cat(paste0("File 1 read.\n\n"))
    print(head(curr_output))

    for (j in 2:length(curr_file_set)){

        cat(paste0("Reading file ", j,".\n\n"))

        curr_file = fread(curr_file_set[j], data.table = FALSE)

        curr_output = rbind(curr_output, curr_file)

        total_snps = nrow(curr_file) + total_snps

    }

    cat(paste0("\nOutput has been made. CHECK 2 COMPLETE.\n\n"))

    cat(paste0("File has: ", nrow(curr_output), " corresponding with a total nrow of individual files = ", total_snps, "\n\n"))

    if (length(strsplit(basename(currGWAS), population)[[1]]) == 1){

        age_strat = "ALL"
        gender_strat = "ALL"

    } else {

        if (grepl("_to_", strsplit(basename(currGWAS), population)[[1]][2]) | grepl("plus", strsplit(basename(currGWAS), population)[[1]][2])){

            gender_strat = "ALL"

            age_strat = strsplit(strsplit(basename(currGWAS), population)[[1]][2], "_")

            if (length(age_strat[[1]]) == 3){

                age_strat = paste0(age_strat[[1]][2], "_99")

            } else {

                age_strat = paste0(age_strat[[1]][2], "_", age_strat[[1]][4])

            }
        } else {

            age_strat = "ALL"

            gender_strat = ifelse(grepl("female", strsplit(basename(currGWAS), population)[[1]][2]), "F", "M")

        }
    }

    sample_size = max(unique(curr_output$N))

    nice_pheno = case_when(currPheno == "gp_costs_year" ~ "PRIM",
                            currPheno == "total_costs_year" ~ "TOTAL",
                            currPheno == "hes_costs_year" ~ "IN",
                            currPheno == "prescription_costs_year" ~ "DRUG",
                            .default = as.character(currPheno)
                            )

    file_output_name = paste0("UKB.SMW.", nice_pheno, ".", age_strat, ".", gender_strat, ".", population, ".", sample_size, ".REGENIE.20240305.txt.gz")

    if (grepl("HES", currGWAS)){

        file_output_name = gsub("TOTAL", "total_costs_year", file_output_name)

    } else {

        file_output_name = gsub("TOTAL", "ALL", file_output_name)

    }

    cat(paste0("\nFile name is: ", file_output_name, "\n\n"))

    fwrite(curr_output, paste0(currGWAS, "/", file_output_name),
            row.names = FALSE, quote = FALSE, sep = "\t", na = "NA", compress = "gzip")

    cat(paste0("\nFile written and zipped. CHECK 3 COMPLETE.\n\n==================\n\n"))

}

