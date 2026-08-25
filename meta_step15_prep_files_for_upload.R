## ---------------------------
##
## Script name: meta_step15_prep_files_for_upload.R
##
## Purpose of script: Convert GWAS files and PGS weights to format for upload to GWAS/PGS catalogue
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-03-19
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gwasDir = paste0(mainDir, "outputs/METAL_v4/")
gwasCatDir = paste0(mainDir, "outputs/METAL_v4_upload/")
pgsDir = paste0(mainDir, "outputs/prs/megaprs_weights/")
pgsCatDir = paste0(mainDir, "outputs/prs/pgs_catalog/")

all_meta_files = Sys.glob(paste0(gwasDir, "*.TBL.gz"))
all_meta_files = all_meta_files[!(grepl("_no", all_meta_files))]
### Reorder: TRUE (contains _ALL_) comes first
all_meta_files = all_meta_files[order(!grepl("_ALL_", all_meta_files))]

all_pgs_files = Sys.glob(paste0(pgsDir, "*BAYESR_effects_for_pgs.txt"))
all_pgs_files = all_pgs_files[!(grepl("_no", all_pgs_files))]

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("I am the very model of a scientist salarian!")


###### =========================== #####

### Running

##### =========================== #####

for (i in 1:length(all_meta_files)){

    curr_name = gsub("_metal_output_1.TBL.gz", "", basename(all_meta_files[i]))

    if (file.exists(paste0(gwasCatDir, curr_name, "_gwas_cat_format.tsv.gz"))){

        cat(paste0("\nGWAS file ", curr_name, " already exists. Skipping.\n...\n\n"))

        next

    }

    cat(paste0("\nInitalising GWAS run on ", curr_name, "\n...\n\n"))

    currFile = fread(all_meta_files[i], data.table = FALSE)

    ### Add betas + SEs

    currFile$beta1 = currFile$Zscore / sqrt((2 * currFile$Freq1) * (1 - currFile$Freq1) * (currFile$Weight + (currFile$Zscore^2)))
    currFile$se = 1 / sqrt((2 * currFile$Freq1) * (1 - currFile$Freq1) * (currFile$Weight + (currFile$Zscore^2)))

    ### Remove very low MAF SNPs

    currFile = currFile[currFile$Freq1 >= 0.001,]
    currFile = currFile[currFile$Freq1 <= 0.999,]

    ### Remove SNPs which are present in fewer than 1/3 cohorts

    third_count = floor(str_count(currFile$Direction)[1] / 3)

    if (third_count == 1){

        third_count = 2

    }

    currFile = currFile %>%
        filter((str_count(Direction, "\\-") + str_count(Direction, "\\+")) >= third_count)

    ### Replace chrX with 23

    currFile$MarkerName = gsub("^X:", "23:", currFile$MarkerName)
    currFile$Chromosome = gsub("X", 23, currFile$Chromosome)

    ### Converting into GWAS catalogue format

    currFile = data.frame(chromosome = currFile$Chromosome,
                            base_pair_location = as.numeric(currFile$Position),
                            effect_allele = toupper(currFile$Allele1),
                            other_allele = toupper(currFile$Allele2),
                            beta = as.numeric(currFile$beta1),
                            standard_error = as.numeric(currFile$se),
                            effect_allele_frequency = as.numeric(currFile$Freq1),
                            p_value = as.numeric(currFile[,"P-value"]),
                            variant_id = paste(currFile$Chromosome,
                                                currFile$Position,
                                                toupper(currFile$Allele1),
                                                toupper(currFile$Allele2),
                                                sep = "_"),
                            rs_id = ifelse(grepl("rs", currFile$MarkerName), currFile$MarkerName, "."),
                            n = as.numeric(currFile$Weight),
                            i_squared = as.numeric(currFile$HetISq),
                            het_pval = as.numeric(currFile$HetPVal),
                            zscore = as.numeric(currFile$Zscore)
                            )

    currFile$base_pair_location = format(currFile$base_pair_location, scientific = FALSE, trim = TRUE)

    fwrite(currFile, paste0(gwasCatDir, curr_name, "_gwas_cat_format.tsv.gz"),
            sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, compress = "gzip")

    rm(currFile)

}

heading("GWAS files processed. Now doing PGS files.")

for (i in 1:length(all_pgs_files)){

    curr_name = gsub("_meta_V4_hapmap_plus_BAYESR_effects_for_pgs.txt", "", basename(all_pgs_files[i]))

    if (file.exists(paste0(pgsCatDir, curr_name, "_pgs_cat_format.tsv.gz"))){

        cat(paste0("\nPGS file ", curr_name, " already exists. Skipping.\n...\n\n"))

        next

    }

    cat(paste0("\nInitalising PGS run on ", curr_name, "\n...\n\n"))

    currFile = fread(all_pgs_files[i], data.table = FALSE)

    a1_col = which(grepl("_a1", colnames(currFile)))
    a0_col = which(grepl("_a0", colnames(currFile)))
    effect_col = which(grepl("_hapmap_plus_BAYESR_megaprs", colnames(currFile)))

    currFile = data.frame(rsID = currFile$rsid,
                            chr_name = currFile$chr,
                            chr_position = currFile$hg38_pos,
                            effect_allele = currFile[,a1_col],
                            other_allele = currFile[,a0_col],
                            effect_weight = currFile[,effect_col]
                            )
   
    fwrite(currFile, paste0(pgsCatDir, curr_name, "_pgs_cat_format.tsv.gz"),
            sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, compress = "gzip")

}

heading("ALL FILES PROCESSED. Rscript finished.")


