### ---------------------------
###
### Script name: qc_step4_3_add_info
###
### Purpose of script: Adds info scores to GNH GWAS
###
### Author: Dr. Sebastian May-Wilson
### Contact: sebastian.may-wilson@helsinki.fi
###
### Date Created: 2025-07-04
###
### ---------------------------
###
##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(dplyr)
library(stringr)
options(scipen = 999)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gwasDir = paste0(mainDir, "outputs/cohort_sumstats/gnh/")
tempDir = paste0(mainDir, "tmpdir/")

### Files

allFiles = Sys.glob(paste0(gwasDir, "*20240107.txt.gz"))
allInfoFiles = Sys.glob(paste0(gwasDir, "gnh_snp_infos/gnh_snp_infos_chr*"))

### Reading files

allInfo = fread(allInfoFiles[1], data.table = FALSE, tmpdir = tempDir)

for (i in 2:length(allInfoFiles)){

    currfile = fread(allInfoFiles[i], data.table = FALSE, tmpdir = tempDir)

    allInfo = rbind(allInfo, currfile)

}

fwrite(allInfo, paste0(gwasDir, "gnh_snp_infos_ALLchrs.txt.gz"), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}


##############################

### Adding INFO and rsIDs to files

##############################

allInfo$snpid = paste0(gsub("chr", "", allInfo$V1), "_", allInfo$V2, "_", allInfo$V3, "_", allInfo$V4)
allInfo = allInfo[,c("V5", "V6", "snpid")]
allInfo$snpid = as.character(allInfo$snpid)

heading("Starting run of adding INFO data.")

for (i in 1:length(allFiles)){

    currfile = fread(allFiles[i], data.table = FALSE, tmpdir = tempDir)

    currfile$ID = as.character(currfile$ID)

    currfile = left_join(currfile, allInfo, by = c("ID" = "snpid"))

    currfile$ID[grepl("rs", currfile$V5)] = currfile$V5[grepl("rs", currfile$V5)]

    colnames(currfile)[which(colnames(currfile) == "V6")] = "INFO"
    currfile$V5 = NULL

    print(summary(currfile))
    print(currfile[sample(1:nrow(currfile), 20),])

    fwrite(currfile, allFiles[i], quote = FALSE, row.names = FALSE, sep = "\t", compress = "gzip", na = "NA")

    heading(paste0("Run ", i, " of ", length(allFiles), " (", allFiles[i], ") complete."))

}

