### ---------------------------
###
### Script name: meta_step11_2_postmunge_combine.R
###
### Purpose of script: Gets list of SNPs with build37 positions and combines them
###
### Author: Dr. Sebastian May-Wilson
### Contact: sebastian.may-wilson@helsinki.fi
###
### Date Created: 2025-05-15
###
### ---------------------------
###
##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(dplyr)
options(scipen = 999)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gwasDir = paste0(mainDir, "outputs/METAL_v4/")
tempDir = paste0(mainDir, "tmpdir/")

### Files

rsid_file = paste0(gwasDir, "ALL_METAL_output_RSIDs_to_hg37.tsv")
snpid_file = paste0(gwasDir, "ALL_METAL_output_SNPIDs_pos_pos1_liftOver_output.out")
prev_file = paste0(gwasDir, "ALL_METAL_output_SNPs.txt.gz")

### Reading

rsids = fread(rsid_file, data.table = FALSE, tmpdir = tempDir)
snpids = fread(snpid_file, data.table = FALSE, tmpdir = tempDir)
hg38_build = fread(prev_file, data.table = FALSE, tmpdir = tempDir)

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

### Combine files

colnames(rsids) = c("MarkerName", "Chr", "hg37_pos")
colnames(snpids) = c("chr", "hg37_pos", "hg37_pos_1", "MarkerName")

rsids = rsids[,c("MarkerName", "hg37_pos")]
snpids = snpids[,c("MarkerName", "hg37_pos")]

hg37_build = rbind(rsids, snpids)

hg38_build = hg38_build %>% left_join(
                                hg37_build,
                                by = "MarkerName"
                                )

fwrite(hg38_build, paste0(gwasDir, "ALL_METAL_output_SNPs_mapped.txt.gz"), row.names = FALSE, sep = "\t", quote = FALSE, compress = "gzip")


