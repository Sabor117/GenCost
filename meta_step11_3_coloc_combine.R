## ---------------------------
##
## Script name: meta_step11_3_coloc_combine.R
##
## Purpose of script: Combine the outputs of coloc_run
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-06-02
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(dplyr)
library(stringr)

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Might go somewhere sunny. Sit on beach, look at ocean, collect sea shells. Might run tests on the sea shells.")

sessionInfo()
start_time = Sys.time()

### Files and directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
colocDir = paste0(mainDir, "outputs/coloc_MAF_0_001/")

all_coloc_files = Sys.glob(paste0(colocDir, "*.tsv"))
all_coloc_files = all_coloc_files[!(grepl("duplicated_var", all_coloc_files))]
all_coloc_files = all_coloc_files[!(grepl("all_phenos", all_coloc_files))]


##### =========================== #####

### Combine

##### =========================== #####

output = fread(all_coloc_files[1], data.table = FALSE)

current_analysis = gsub(".tsv", "", basename(all_coloc_files[1]))
curr_pheno = strsplit(current_analysis, "_")[[1]][3]
curr_query = strsplit(current_analysis, "_")[[1]][4]

output$phenotype = curr_pheno
output$query = curr_query

for (i in 2:length(all_coloc_files)){

    curr_file = fread(all_coloc_files[i], data.table = FALSE)

    current_analysis = gsub(".tsv", "", basename(all_coloc_files[i]))
    curr_pheno = strsplit(current_analysis, "_")[[1]][3]

    if (length(strsplit(current_analysis, "_")[[1]]) == 4){

        curr_query = strsplit(current_analysis, "_")[[1]][4]

    } else if (length(strsplit(current_analysis, "_")[[1]]) == 5){

        curr_query = paste0(
                            strsplit(current_analysis, "_")[[1]][4],
                            "_",
                            strsplit(current_analysis, "_")[[1]][5]
        )

    }

    curr_file$phenotype = curr_pheno
    curr_file$query = curr_query

    output = rbind(output, curr_file)

}

fwrite(output, paste0(colocDir, "full_coloc_res_all_phenos_v4_1.tsv"),
            quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

output_pruned = na.omit(output)
output_pruned = output[output$PP4 >= 0.9,]

fwrite(output_pruned, paste0(colocDir, "strong_coloc_res_all_phenos_v4_1.tsv"),
            quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)
