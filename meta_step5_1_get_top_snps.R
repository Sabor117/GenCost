##### =========================== #####

### Setup environment

### Inputs: GCTA-COJO annotated top SNPs + meta-analysis sumstats
### Outputs: top SNPs from GCTA-COJO, annotated with betas, P-vals, best gene and presence in each meta-analysis

### Output of SNP list for extraction from UKB for cost "QTLs"

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)
set.seed(117)

phenotype = "IN_ALL"
maf_filter = "MAF_0_01"

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gcta_cojo_Dir = paste0(mainDir, "outputs/gcta_cojo_v4_", maf_filter, "/jma_combined/")
metaDir = paste0(mainDir, "outputs/METAL_v4/")
outDir = paste0(mainDir, "processing/cost_per_snp_processing/")

top_snps_file_list = Sys.glob(paste0(gcta_cojo_Dir, phenotype, "*jma_out_annotated.txt"))
top_snps_file_list = top_snps_file_list[!(grepl("_no", top_snps_file_list))]


##### =========================== #####

### Run script

##### =========================== #####

topsnps = data.frame(matrix(ncol = 3, nrow = 0))
colnames(topsnps)[1:3] = c("snpid", "chr", "gene")

for (i in 1:length(top_snps_file_list)){

    curr_gcta = fread(top_snps_file_list[i], data.table = FALSE)

    curr_set = data.frame(snpid = curr_gcta$SNP,
                            chr = curr_gcta$Chr,
                            gene = curr_gcta$nearest_ncbi_gene)

    topsnps = rbind(topsnps, curr_set)

}

topsnps = unique(topsnps)

for (i in 1:length(top_snps_file_list)){

    meta_name = strsplit(basename(top_snps_file_list[i]), "_gcta")[[1]][1]

    meta_file_name = paste0(metaDir, meta_name, "_metal_output_1.TBL.gz")

    cat(paste0("\nReading file: ", meta_file_name, "\n\n"))

    meta_file = fread(meta_file_name, data.table = FALSE)

    meta_file = meta_file[meta_file$MarkerName %in% topsnps$snpid,]

    meta_file = meta_file[,c("MarkerName", "Zscore", "P-value")]
    colnames(meta_file) = c("snpid", paste0(meta_name, "_Z"), paste0(meta_name, "_P"))

    topsnps = merge(topsnps, meta_file, by = "snpid", all.x = TRUE)

}

fwrite(topsnps, paste0(outDir, phenotype, "_top_snps_", maf_filter, "_comparison.txt"), row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")

for (i in 1:length(unique(topsnps$chr))){

    currChr = unique(topsnps$chr)[i]

    topsnps_currChr = topsnps[topsnps$chr == currChr,]

    topsnps_list = data.frame(rsid = topsnps_currChr$snpid)

    fwrite(topsnps_list, paste0(outDir, phenotype, "_top_snps_", maf_filter, "_list_chr", currChr, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t", na = "NA", col.names = FALSE)

}














all_metas = Sys.glob(paste0(metaDir, "DRUG_*TBL.gz"))
all_metas = all_metas[!(grepl("noUKB", all_metas))]
all_metas = all_metas[!(grepl("19_", all_metas))]
all_metas = all_metas[!(grepl("_18_", all_metas))]
all_metas = all_metas[!(grepl("_36_", all_metas))]

meta_snps = fread(all_metas[1], data.table = FALSE)

meta_name = strsplit(basename(all_metas[1]), "_metal")[[1]][1]

meta_snps = meta_snps[meta_snps$MarkerName %in% topsnps,]
meta_snps = meta_snps[,c("MarkerName", "Zscore", "Weight", "Freq1")]
colnames(meta_snps) = c("snpid", paste0(meta_name, "_Z"), paste0(meta_name, "_N"), paste0(meta_name, "_AF"))

for (i in 2:length(all_metas)){

    cat(paste0("Now running on ", i, "\n"))

    curr_meta = fread(all_metas[i], data.table = FALSE)

    meta_name = strsplit(basename(all_metas[i]), "_metal")[[1]][1]

    curr_meta = curr_meta[curr_meta$MarkerName %in% topsnps,]
    curr_meta = curr_meta[,c("MarkerName", "Zscore", "Weight", "Freq1")]
    colnames(curr_meta) = c("snpid", paste0(meta_name, "_Z"), paste0(meta_name, "_N"), paste0(meta_name, "_AF"))

    meta_snps = merge(meta_snps, curr_meta, by = "snpid")

}

fwrite(meta_snps, "/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/beta_comparison_frame_DRUG.txt", row.names = FALSE, quote = FALSE, sep = "\t")


calculate_sig_diffs = function(data, zscore1, zscore2, n1, n2){

    se_1 = 1 / sqrt(data[,zscore1] ^ 2 + data[,n1])
    beta_1 = data[,zscore1] * se_1

    se_2 = 1 / sqrt(data[,zscore2] ^ 2 + data[,n2])
    beta_2 = data[,zscore2] * se_2

    beta_diff = beta_1 - beta_2
    se_diff = sqrt(se_1^2 + se_2^2 - 2 * 0.9 * se_1 * se_2)

    z_diff = beta_diff / se_diff

    p_diff = 2 * pnorm(-abs(z_diff))

    return(p_diff)

}

c(which(calculate_sig_diffs(meta_snps_drug, "DRUG_M_Z", "DRUG_F_Z", "DRUG_M_N", "DRUG_F_N") < 0.05/105))

c(which(calculate_sig_diffs(meta_snps_drug, "DRUG_56_75_Z", "DRUG_76_95_Z", "DRUG_56_75_N", "DRUG_76_95_N") < 0.05/105),
    which(calculate_sig_diffs(meta_snps_drug, "DRUG_56_75_Z", "DRUG_36_55_Z", "DRUG_56_75_N", "DRUG_36_55_N") < 0.05/105),
    which(calculate_sig_diffs(meta_snps_drug, "DRUG_56_75_Z", "DRUG_19_35_Z", "DRUG_56_75_N", "DRUG_19_35_N") < 0.05/105))
