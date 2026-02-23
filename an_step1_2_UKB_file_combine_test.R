##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)
library(R.utils)

options(scipen = 999) # prevent scientific notation in output

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gwasDir = paste0(mainDir, "outputs/regenie/step2_output/HES_CPI_test/")

all_gwas_files = Sys.glob(paste0(gwasDir, "*.regenie"))

if (length(all_gwas_files) == 0){

    cat(paste0("\nGWAS has no valid files for combining. CHECK 3 COMPLETE.\n\n==================\n\n"))

    stop()

}

gwas_set = c("INPAT_0_96", "INPAT_1_02", "INPAT")

for (ngwas in 1:length(gwas_set)){

    currpheno = gwas_set[ngwas]

    if (currpheno == "INPAT"){

        gwas_files = all_gwas_files[!(grepl("1_02", all_gwas_files))]
        gwas_files = all_gwas_files[!(grepl("0_96", all_gwas_files))]

    } else {

        gwas_files = all_gwas_files[grepl(currpheno, all_gwas_files)]

    }

    gwas = fread(gwas_files[1], data.table = FALSE)

    for (nfile in 2:length(gwas_files)){

        currfile = fread(inpat_96[nfile], data.table = FALSE)

        gwas = rbind(gwas, currfile)

        cat(paste0("\nCurrent file #", nfile, " added.\n...\n\n"))

    }

    fwrite(gwas, paste0(gwasDir, "CPI_test_", currpheno, "_regenie.txt.gz"),
            quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

    heading(paste0(currpheno, " GWAS written."))

    rm("gwas")

}

