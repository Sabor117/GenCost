################################################################################
#
# Set up script
#
################################################################################

### Packages

library(data.table)
library(dplyr)

### Files

meta_sumstats_file_list = c("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/PRIM_ALL_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/DRUG_ALL_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/IN_ALL_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/PRIM_ALL_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/IN_ALL_noFINNGEN_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/IN_ALL_noUKB_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/DRUG_ALL_noFINNGEN_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/DRUG_ALL_noUKB_metal_output_1.TBL.gz",
                            "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/INOUT_ALL_noFINNGEN_metal_output_1.TBL.gz"
                            )

gencost_allSNPs_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/gencost_allSNPs_markerid_rsids.txt.gz"
hapmap3_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/hm3_plus_snp_list.txt"
hg37_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt"
hg38_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out"

### Locations

outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/"

### Load reference genome - Fiona Hagenbeek - 1000G

ref_gen = fread("/scratch/project_2007428/users/FAHagenbeek/data/1000Greferencepanel/ref.bim", data.table = FALSE, header = FALSE)
ref_gen = ref_gen[,c(2,5,6)] # retain relevant columns only
colnames(ref_gen) = c("Predictor","A1","A2") # rename columns

### Add IUPAC allele coding - based on Fiona's code

ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='A'] = 'W'
ref_gen$IUPAC[ref_gen$A1 == 'C' & ref_gen$A2 =='G' | ref_gen$A1 == 'G' & ref_gen$A2 =='C'] = 'S'
ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='G' | ref_gen$A1 == 'G' & ref_gen$A2 =='A'] = 'R'
ref_gen$IUPAC[ref_gen$A1 == 'C' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='C'] = 'Y'
ref_gen$IUPAC[ref_gen$A1 == 'G' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='G'] = 'K'
ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='C' | ref_gen$A1 == 'C' & ref_gen$A2 =='A'] = 'M'

### Read hg37 build positions from UKB

hg37_snps = fread(hg37_snps_file, data.table = FALSE)
hg37_snps = hg37_snps[,c("V4", "V2")]
colnames(hg37_snps) = c("MarkerName", "hg37_pos")

### Read HapMap SNPs

hapmap3_snps = readRDS(hapmap3_snps_file)

cat(paste0("\nNrow of all HapMap SNPs: ", nrow(hapmap3_snps), "\n\n"))
print(head(hapmap3_snps))

for (nfile in 1:length(meta_sumstats_file_list)){

    ################################################################################
    #
    # Load Data
    #
    ################################################################################

    meta_sumstats_file = meta_sumstats_file_list[nfile]

    ### Analysis name

    analysis_name = strsplit(basename(meta_sumstats_file), "_metal_output_")[[1]][1]

    ### Read sumstats and reduce to those from HapMap

    meta_sumstats = fread(meta_sumstats_file, data.table = FALSE)

    cat(paste0("\nNrow of all meta SNPs: ", nrow(meta_sumstats), "\n\n"))
    cat(paste0("\nMeta analysed: ", meta_sumstats_file, "\n\n"))
    print(head(meta_sumstats, n = 10))

    start_snp_no = nrow(meta_sumstats)

    ### Reduce to HapMap SNPs

    meta_sumstats = meta_sumstats[meta_sumstats$MarkerName %in% hapmap3_snps$rsid,]

    cat(paste0("\nNrow of meta SNPs in HapMap: ", nrow(meta_sumstats), "\n\n"))
    print(head(meta_sumstats, n = 10))
    cat(paste0("\nCorresponding to ", start_snp_no - nrow(meta_sumstats), " SNPs lost when refining to HapMap.\n\n"))

    step1_snp_no = nrow(meta_sumstats)

    ### Remove any potential duplicates here

    if (any(duplicated(meta_sumstats$MarkerName))){

        duplicated_snp_ids = meta_sumstats$MarkerName[duplicated(meta_sumstats$MarkerName)]

        duplicated_snps = meta_sumstats[meta_sumstats$MarkerName %in% duplicated_snp_ids,]

        cat(paste0("\n", length(duplicated_snp_ids), " SNPs removed for having duplicated IDs post-merging.\n\nCorresponds to ", nrow(duplicated_snps), " rows of SNPs.\n\n"))
        print(head(duplicated_snps, n = 20))

        meta_sumstats = meta_sumstats[!(meta_sumstats$MarkerName %in% duplicated_snp_ids),]

    }

    ### Convert to build 37 positions

    meta_sumstats = merge(meta_sumstats, hg37_snps, by = "MarkerName")

    cat("\n\n======================\n\n")

    ### Make GenCOST PGS data frame

    megaprs_sumstats = data.frame(Predictor = paste0(meta_sumstats$Chromosome, ":", meta_sumstats$hg37_pos),
                                A1 = toupper(meta_sumstats$Allele1),
                                A2 = toupper(meta_sumstats$Allele2),
                                n = as.numeric(meta_sumstats$Weight),
                                Z = as.numeric(meta_sumstats$Zscore))

    megaprs_sumstats$Predictor = gsub("_", ":", megaprs_sumstats$Predictor)

    ### Keep only predictor names and IUPAC allele coding

    ref_gen = ref_gen[,c("Predictor","IUPAC")]

    ### Add IUPAC to sumstats

    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'A' & megaprs_sumstats$A2 =='T' | megaprs_sumstats$A1 == 'T' & megaprs_sumstats$A2 =='A'] = 'W'
    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'C' & megaprs_sumstats$A2 =='G' | megaprs_sumstats$A1 == 'G' & megaprs_sumstats$A2 =='C'] = 'S'
    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'A' & megaprs_sumstats$A2 =='G' | megaprs_sumstats$A1 == 'G' & megaprs_sumstats$A2 =='A'] = 'R'
    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'C' & megaprs_sumstats$A2 =='T' | megaprs_sumstats$A1 == 'T' & megaprs_sumstats$A2 =='C'] = 'Y'
    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'G' & megaprs_sumstats$A2 =='T' | megaprs_sumstats$A1 == 'T' & megaprs_sumstats$A2 =='G'] = 'K'
    megaprs_sumstats$IUPAC[megaprs_sumstats$A1 == 'A' & megaprs_sumstats$A2 =='C' | megaprs_sumstats$A1 == 'C' & megaprs_sumstats$A2 =='A'] = 'M'

    ### retain only R, Y, K, and M IUPAC coded alleles - to remove ambiguous alleles

    megaprs_sumstats = megaprs_sumstats[(megaprs_sumstats$IUPAC %in% c('R', 'Y', 'K', 'M')),]

    ### Only retain SNPs where IUPAC alleles in summary statistics and reference file match

    megaprs_sumstats = inner_join(megaprs_sumstats, ref_gen, by = c("Predictor", "IUPAC"))

    predictor_set = data.frame(predictor = megaprs_sumstats$Predictor)

    ### Write summary statistics to file

    fwrite(megaprs_sumstats, paste0(outDir, analysis_name, "_meta_sumstats_for_megaPRS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    fwrite(predictor_set, paste0(outDir, analysis_name, "_meta_snplist_for_megaPRS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

}

