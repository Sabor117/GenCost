##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
options(scipen = 999)

pheno_list = c("IN_ALL", "INOUT_ALL", "PRIM_ALL", "DRUG_ALL")

for (i in 1:length(pheno_list)){

    phenotype = pheno_list[i]

    mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
    sumstatDir = paste0(mainDir, "processing/meta_sumstats/", phenotype, "/")
    outDir = paste0(mainDir, "processing/ldsc_intermediate_files/", phenotype, "/")

    hm3_snps = fread(paste0(mainDir, "processing/ldsc_intermediate_files/hm3_snp_list.txt"), data.table = FALSE)

    hm3_snps = hm3_snps[,c(2, 3, 4, 5)]
    colnames(hm3_snps) = c("chr", "pos", "pos1", "rsid")

    hm3_snps$snpid = paste0(gsub("chr", "", hm3_snps$chr), "_", hm3_snps$pos)
    hm3_snps$snpid2 = paste0(gsub("chr", "", hm3_snps$chr), "_", hm3_snps$pos1)

    hm3_snps = unique(hm3_snps)

    sumstat_files = Sys.glob(paste0(sumstatDir, "*.txt.gz"))

    for (i in 1:length(sumstat_files)){

        currFile = fread(sumstat_files[i], data.table = FALSE)

        currCohort = strsplit(basename(sumstat_files[i]), "\\.")[[1]][1]

        if (currCohort == "UKB"){

            currCohort = paste0(currCohort, "_", strsplit(basename(sumstat_files[i]), "\\.")[[1]][6])

        }

        #currFile_1 = merge(currFile, hm3_snps[,c("snpid", "rsid")], by = "snpid")
        #currFile_2 = merge(currFile, hm3_snps[,c("snpid2", "rsid")], by.x = "snpid", by.y = "snpid2")

        #outFile = rbind(currFile_1, currFile_2)

        outFile = currFile[currFile$snpid %in% hm3_snps$rsid,]

        outFile$zscore1 = outFile$beta1 / outFile$se

        outFile = outFile[,c("snpid", "a1", "a0", "n", "p", "zscore1")]

        outname = paste0(outDir, currCohort, "_", phenotype, "_ldsc_input.txt.gz")

        fwrite(outFile, outname, row.names = FALSE, na = "NA", quote = FALSE, sep = "\t", compress = "gzip")

    }

}
