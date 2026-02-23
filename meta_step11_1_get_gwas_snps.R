### ---------------------------
###
### Script name: meta_step11_1_premunge_combine.R
###
### Purpose of script: Combine list of SNPs from meta-analyses and get build 37 positions for them
###
## Author: Dr. Sebastian May-Wilson
### Contact: sebastian.may-wilson@helsinki.fi
###
### Date Created: 2025-05-06
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

### Files

mainFiles = Sys.glob(paste0(gwasDir, "*_ALL_metal_output_1.TBL.gz"))

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}


##### =========================== #####

### Run combine

##### =========================== #####

heading("'Tis but a scratch!")

sessionInfo()
start_time = Sys.time()

combined_file = fread(mainFiles[1], data.table = FALSE, select = c("MarkerName", "Chromosome", "Position", "Allele1", "Allele2"))

combined_file_rs = combined_file[grepl("rs", combined_file$MarkerName),]
combined_file_snpid = combined_file[!(grepl("rs", combined_file$MarkerName)),]

print(head(combined_file))
cat(paste0("\nFile ", 1, " of ", length(mainFiles), " read and combined.\n",
                "Current file: ", mainFiles[1], "\n",
                "Currently there are: ", nrow(combined_file_rs), " rsIDs.\n",
                "Currently there are: ", nrow(combined_file_snpid), " rsIDs.\n",
                ".....\n\n"))

for (i in 2:length(mainFiles)){

    loop_time = Sys.time()

    curr_file = fread(mainFiles[i], data.table = FALSE, select = c("MarkerName", "Chromosome", "Position", "Allele1", "Allele2"))

    curr_file_rs = curr_file[grepl("rs", curr_file$MarkerName),]
    curr_file_snpid = curr_file[!(grepl("rs", curr_file$MarkerName)),]

    combined_file_rs = rbind(combined_file_rs, curr_file_rs)
    combined_file_rs = distinct(combined_file_rs)

    combined_file_snpid = rbind(combined_file_snpid, curr_file_snpid)
    combined_file_snpid = distinct(combined_file_snpid)

    run_time = Sys.time()

    cat(paste0("\nFile ", i, " of ", length(mainFiles), " read and combined.\n",
                "Current file: ", mainFiles[i], "\n",
                "Currently there are: ", nrow(combined_file_rs), " rsIDs.\n",
                "Currently there are: ", nrow(combined_file_snpid), " rsIDs.\n",
                "This file took: ", run_time - loop_time, " to run.\n",
                "Total time taken: ", run_time - start_time, " to run.\n",
                ".....\n\n"))

    gc()

}

heading("Combine done.")

colnames(combined_file_snpid)[3] = "hg38_pos"
colnames(combined_file_rs)[3] = "hg38_pos"

### For file with only rsIDs - for getting pos 37

rsids_frame = data.frame(rsid = combined_file_rs$MarkerName)

### For bed file with SNP IDs - for LiftOver

snpids_bed = data.frame(chr = paste0("chr", combined_file_snpid$Chromosome),
                        pos = as.numeric(combined_file_snpid$hg38_pos - 1),
                        pos1 = as.numeric(combined_file_snpid$hg38_pos),
                        snpid = combined_file_snpid$MarkerName)

### For scientific purposes

snpids_bed_1 = data.frame(chr = paste0("chr", combined_file_snpid$Chromosome),
                        pos = as.numeric(combined_file_snpid$hg38_pos),
                        pos1 = as.numeric(combined_file_snpid$hg38_pos),
                        snpid = combined_file_snpid$MarkerName)

snpids_bed_2 = data.frame(chr = paste0("chr", combined_file_snpid$Chromosome),
                        pos = as.numeric(combined_file_snpid$hg38_pos),
                        pos1 = as.numeric(combined_file_snpid$hg38_pos + 1),
                        snpid = combined_file_snpid$MarkerName)

combined_file = rbind(combined_file_rs, combined_file_snpid)
combined_file$Allele1 = toupper(combined_file$Allele1)
combined_file$Allele2 = toupper(combined_file$Allele2)

fwrite(combined_file_snpid, paste0(gwasDir, "ALL_METAL_output_SNPIDs.txt.gz"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(combined_file_rs, paste0(gwasDir, "ALL_METAL_output_RSIDs.txt.gz"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(rsids_frame, paste0(gwasDir, "ALL_METAL_output_RSIDs_only.txt"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)
#fwrite(snpids_bed, paste0(gwasDir, "ALL_METAL_output_SNPIDs.bed"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)

#fwrite(snpids_bed_1, paste0(gwasDir, "ALL_METAL_output_SNPIDs_pos_pos.bed"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)
fwrite(snpids_bed_2, paste0(gwasDir, "ALL_METAL_output_SNPIDs_pos_pos1.bed"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)

fwrite(combined_file_snpid, paste0(gwasDir, "ALL_METAL_output_SNPs.txt.gz"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

heading("Script complete.")

