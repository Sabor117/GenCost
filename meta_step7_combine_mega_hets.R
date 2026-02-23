## ---------------------------
##
## Script name: meta_step7_combine_mega_hets.R
##
## Purpose of script: Extract heterogeneities from MR-MEGA outputs for lead SNPs
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-09-25
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(dplyr)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
processDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")
mrmegaDir = paste0(mainDir, "outputs/mr_mega/")
gwasDir = paste0(mainDir, "outputs/METAL_v4/")
outDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/top_snp_matrix/")

### File list

main_analysis_files = Sys.glob(paste0(processDir, "*ALL*annotated*"))
main_analysis_files = main_analysis_files[!(grepl("_no", main_analysis_files))]

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("For those who come after.")


##### =========================== #####

### Start run

##### =========================== #####

for (i in 1:length(main_analysis_files)){

    currfile = fread(main_analysis_files[i], data.table = FALSE)

    phenotype = strsplit(basename(main_analysis_files[i]), "_gcta_")[[1]][1]

    heading(paste0("Now running on ", phenotype))

    all_mr_mega_files = Sys.glob(paste0(mrmegaDir, phenotype, "*result"))
    gwas_file = paste0(gwasDir, phenotype, "_metal_output_1.TBL.gz")

    output = fread(gwas_file, tmpdir = paste0(mainDir, "tmpdir/"))

    output = output[output$MarkerName %in% currfile$SNP,]
    output = output[,c("MarkerName", "HetChiSq", "HetPVal", "HetDf")]
    colnames(output) = c("MarkerName", "gwas_het", "gwas_het_pval", "gwas_het_ndf")

    for (nfile in 1:length(all_mr_mega_files)){

        curr_mr_mega_file = all_mr_mega_files[nfile]
        curr_covar = strsplit(basename(curr_mr_mega_file), "_")[[1]][3]
        curr_covar = gsub("\\.result", "", curr_covar)

        curr_mr_mega = fread(curr_mr_mega_file, tmpdir = paste0(mainDir, "tmpdir/"))

        heading(paste0(curr_covar, " file to merge (", all_mr_mega_files[nfile], ")"))
        print(head(curr_mr_mega))

        colskeep = c("MarkerName", colnames(curr_mr_mega)[grepl("ancestry_het", colnames(curr_mr_mega))])
        print(colskeep)
        curr_mr_mega = curr_mr_mega[, ..colskeep]
        print(curr_mr_mega)

        colnames(curr_mr_mega) = c("MarkerName", paste0(curr_covar, "_het"), paste0(curr_covar, "_het_ndf"), paste0(curr_covar, "_het_pval"))

        output = left_join(output, curr_mr_mega, by = "MarkerName")

        cat(paste0("\n...\n\n"))

    }

    fwrite(output, paste0(outDir, phenotype, "_mr_mega_het_res.txt"), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")

    cat(paste0("\n", phenotype, " complete. Results:\n\n"))
    print(head(output))
    cat(paste0("\n..................\n\n"))

}

heading("0.001 MAF run complete.")

processDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_01/jma_combined/")
outDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_01/top_snp_matrix/")

for (i in 1:length(main_analysis_files)){

    currfile = fread(main_analysis_files[i], data.table = FALSE)

    phenotype = strsplit(basename(main_analysis_files[i]), "_gcta_")[[1]][1]

    heading(paste0("Now running on ", phenotype))

    all_mr_mega_files = Sys.glob(paste0(mrmegaDir, phenotype, "*result"))
    gwas_file = paste0(gwasDir, phenotype, "_metal_output_1.TBL.gz")

    output = fread(gwas_file, tmpdir = paste0(mainDir, "tmpdir/"))

    output = output[output$MarkerName %in% currfile$SNP,]
    output = output[,c("MarkerName", "HetChiSq", "HetPVal", "HetDf")]
    colnames(output) = c("MarkerName", "gwas_het", "gwas_het_pval", "gwas_het_ndf")

    for (nfile in 1:length(all_mr_mega_files)){

        curr_mr_mega_file = all_mr_mega_files[nfile]
        curr_covar = strsplit(basename(curr_mr_mega_file), "_")[[1]][3]
        curr_covar = gsub("\\.result", "", curr_covar)

        curr_mr_mega = fread(curr_mr_mega_file, tmpdir = paste0(mainDir, "tmpdir/"))

        heading(paste0(curr_covar, " file to merge (", all_mr_mega_files[nfile], ")"))
        print(head(curr_mr_mega))

        colskeep = c("MarkerName", colnames(curr_mr_mega)[grepl("ancestry_het", colnames(curr_mr_mega))])
        print(colskeep)
        curr_mr_mega = curr_mr_mega[, ..colskeep]
        print(curr_mr_mega)

        colnames(curr_mr_mega) = c("MarkerName", paste0(curr_covar, "_het"), paste0(curr_covar, "_het_ndf"), paste0(curr_covar, "_het_pval"))

        output = left_join(output, curr_mr_mega, by = "MarkerName")
        cat(paste0("\n...\n\n"))

    }

    fwrite(output, paste0(outDir, phenotype, "_mr_mega_het_res.txt"), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")

    cat(paste0("\n", phenotype, " complete. Results:\n\n"))
    print(head(output))
    cat(paste0("\n..................\n\n"))

}

heading("0.01 MAF run complete.")

heading("Both runs complete.")


