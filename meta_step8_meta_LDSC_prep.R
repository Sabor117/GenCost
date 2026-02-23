##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
options(scipen = 999)

phenotype = "META_v4"

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/METAL_v4/")
outDir = paste0(mainDir, "processing/ldsc_intermediate_files/", phenotype, "/")

hm3_snps = fread(paste0(mainDir, "processing/ldsc_intermediate_files/hm3_snp_list.txt"), data.table = FALSE)

hm3_snps = hm3_snps[,c(2, 3, 4, 5)]
colnames(hm3_snps) = c("chr", "pos", "pos1", "rsid")

hm3_snps = unique(hm3_snps)

sumstat_files = Sys.glob(paste0(sumstatDir, "*.TBL.gz"))
sumstat_files = sumstat_files[grepl("_ALL_", sumstat_files)]
sumstat_files = sumstat_files[!(grepl("_no", sumstat_files))]

for (i in 1:length(sumstat_files)){

    currFile = fread(sumstat_files[i], data.table = FALSE)

    currPheno = strsplit(basename(sumstat_files[i]), "_metal")[[1]][1]

    outFile = currFile[currFile$MarkerName %in% hm3_snps$rsid,]

    outFile = outFile[,c("MarkerName", "Allele1", "Allele2", "Weight", "P-value", "Zscore")]

    colnames(outFile) = c("snpid", "a1", "a0", "n", "p", "zscore1")
    outFile$a1 = toupper(outFile$a1)
    outFile$a0 = toupper(outFile$a0)

    outname = paste0(outDir, currPheno, "_META_v4_ldsc_input.txt.gz")

    fwrite(outFile, outname, row.names = FALSE, na = "NA", quote = FALSE, sep = "\t", compress = "gzip")

}