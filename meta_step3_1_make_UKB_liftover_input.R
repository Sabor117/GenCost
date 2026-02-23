##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)

options(scipen = 999) # prevent scientific notation in output

dataDir = "/scratch/project_2007428/data/processing/ukbb_78537/genotypes/white_british_30k_reference_panel/"
outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/"

all_bims = Sys.glob(paste0(dataDir, "ukb22828_chr*_30k_random_unrelated_white_british.bim"))

all_bims = all_bims[!(grepl("chr23", all_bims))]


##### =========================== #####

### Make full BIM and input for LiftOver

##### =========================== #####

output_bim = fread(all_bims[1], data.table = FALSE)

for (i in 2:length(all_bims)){

    cat(paste0("\nWorking on file ", i, " out of ", length(all_bims), ". Current file is: ", basename(all_bims[i]), "\n\n"))

    currBim = fread(all_bims[i], data.table = FALSE)

    output_bim = rbind(output_bim, currBim)

}

fwrite(output_bim, paste0(dataDir, "ukb22828_ALLSNPs_30k_random_unrelated_white_british.bim"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

cat(paste0("\nFull BIM file created and written.\n\n=====================\n\n"))

liftover_input = data.frame(chr = paste0("chr", output_bim[,1]),
                                pos = as.numeric(output_bim[,4]),
                                pos1 = as.numeric(output_bim[,4]) + 1,
                                rsid = output_bim[,2],
                                snpid = paste0(output_bim[,1], ":", output_bim[,4], ":", output_bim[,5], ":", output_bim[,6]))

liftover_input$chr = gsub("chr23", "chrX", liftover_input$chr)

fwrite(liftover_input, paste0(outDir, "ukb_ALLSNPs_hg19_liftover_input_bim.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

cat(paste0("\nLiftOver input created and written.\n\n=====================\n\n"))

gwas_files = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cohort_sumstats/ukb/UKB.SMW.*.txt.gz")

liftover_input_two = data.frame(matrix(ncol = 5, nrow = 0))
colnames(liftover_input_two) = c("chr", "pos", "pos1", "rsid", "snpid")

for (i in 1:length(gwas_files)){

    curr_file = fread(gwas_files[i], select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1"))

    cat(paste0("\nWorking on file ", i, " out of ", length(gwas_files), ". Current file is: ", basename(gwas_files[i]), "\n\n"))

    curr_file$GENPOS = as.numeric(curr_file$GENPOS)
    curr_file$CHROM = gsub("23", "X", curr_file$CHROM)
    curr_file$CHROM = gsub(23, "X", curr_file$CHROM)

    curr_snp_set = data.frame(chr = paste0("chr", curr_file$CHROM),
                                pos = curr_file$GENPOS,
                                pos1 = curr_file$GENPOS + 1,
                                rsid = curr_file$ID,
                                snpid = paste0(curr_file$CHROM, ":", curr_file$GENPOS, ":", curr_file$ALLELE1, ":", curr_file$ALLELE0)
                                )

    print(head(curr_snp_set))

    curr_snp_set$chr = gsub("chr23", "chrX", curr_snp_set$chr)

    liftover_input_two = rbind(liftover_input_two, curr_snp_set)

    liftover_input_two = unique(liftover_input_two)

    cat(paste0("\nFile complete. Next.\n\n"))

}

cat(paste0("\nFull BIM file created. Unique rows selected. Ordering.\n\n=====================\n\n"))

liftover_input_two$chr_order = NA

liftover_input_two$chr_order[grepl("X", liftover_input_two$chr)] = 23
liftover_input_two$chr_order[grepl("Y", liftover_input_two$chr)] = 24
liftover_input_two$chr_order[grepl("MT|M", liftover_input_two$chr)] = 25

idx = is.na(liftover_input_two$chr_order)
liftover_input_two$chr_order[idx] = as.numeric(gsub("^chr", "", liftover_input_two$chr[idx]))

# Sort by chr and pos
liftover_input_two = liftover_input_two[order(liftover_input_two$chr_order, liftover_input_two$pos), ]

liftover_input_two$chr_order = NULL

fwrite(liftover_input_two, paste0(outDir, "ukb_ALLSNPs_hg19_liftover_input.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

cat(paste0("\nLiftOver input created and written.\n\n=====================\n\n"))

