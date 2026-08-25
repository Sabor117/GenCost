##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
options(scipen = 999)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
outDir = paste0(mainDir, "processing/ldsc_intermediate_files/SENSITIVITY_TEST/")

file_list = Sys.glob(paste0(mainDir, "outputs/cohort_sumstats/ukb/UKB.SMW.*.ALL.ALL.*EUR.*"))
file_list = file_list[!(grepl("UKB.SMW.ALL.ALL.ALL.EUR.", file_list))]

for (i in 1:length(file_list)){

    phenotype = strsplit(basename(file_list[i]), "\\.")[[1]][3]

    hm3_snps = fread(paste0(mainDir, "processing/ldsc_intermediate_files/hm3_snp_list.txt"), data.table = FALSE)

    hm3_snps = hm3_snps[,c(2, 3, 4, 5)]
    colnames(hm3_snps) = c("chr", "pos", "pos1", "rsid")

    hm3_snps$snpid = paste0(gsub("chr", "", hm3_snps$chr), "_", hm3_snps$pos)
    hm3_snps$snpid2 = paste0(gsub("chr", "", hm3_snps$chr), "_", hm3_snps$pos1)

    hm3_snps = unique(hm3_snps)

    currFile = fread(file_list[i], data.table = FALSE)

    outFile = currFile[currFile$ID %in% hm3_snps$rsid,]

    outFile$zscore1 = outFile$BETA / outFile$SE

    outFile = outFile[,c("ID", "ALLELE1", "ALLELE0", "N", "LOG10P", "zscore1")]

    outFile$p = 10^(-outFile$LOG10P)

    colnames(outFile)[1:4] = c("rsid", "a1", "a0", "n")
    outFile = outFile[,c("rsid", "a1", "a0", "n", "p", "zscore1")]

    stratif = case_when(grepl("RAW", basename(file_list[i])) ~ "RAW_COST",
                        grepl("ZERO", basename(file_list[i])) ~ "NO_ZERO",
                        .default = "ALL")

    outname = paste0(outDir, "UKB_", phenotype, "_", stratif, "_ldsc_input.txt.gz")

    fwrite(outFile, outname, row.names = FALSE, na = "NA", quote = FALSE, sep = "\t", compress = "gzip")

}
