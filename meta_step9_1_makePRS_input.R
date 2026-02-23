##### =========================== #####

### Setup environment

### Script is designed to make pruned sumstats for PRS creation
### Format is designed for LDpred2 + PGSFusion

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)
library(PathWAS, lib.loc = "/projappl/project_2007428/RPackages_421/")

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
tempDir = paste0(mainDir, "tmpdir/")
outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/megaPRS_processing/"

hg37_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt"

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Had to be me. Someone else would have gotten it wrong.")


##### =========================== #####

### Data setup

##### =========================== #####

### Selecting and reading meta-analysis based on list of metas

analysis_list = Sys.glob(paste0(metaDir, "*ALL*noFINNGEN*TBL.gz"))
analysis_list = analysis_list[!(grepl("pruned", analysis_list))]

### Reading hapmap3+ SNPs

hapmap_set = "HMPLUS"

if (hapmap_set == "HM3"){

    hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_snp_list.txt")
    hapmap3_snps = fread(hapmap3_snps_file)

    outDir = paste0(outDir, "hapmap_3_sumstats/")

    hapmap_name = "hapmap3"

} else if (hapmap_set == "HMPLUS"){

    hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_plus_snp_list.rds")
    hapmap3_snps = readRDS(hapmap3_snps_file)

    outDir = paste0(outDir, "hapmap_plus_sumstats/")

    hapmap_name = "hapmap_plus"

}


### Read hg37 file

hg37_snps = fread(hg37_snps_file, data.table = FALSE)

hg37_snps = hg37_snps[,c(2,4)]
colnames(hg37_snps) = c("hg37_pos", "MarkerName")

### Load reference genome - Fiona Hagenbeek - 1000G

ref_gen = fread("/scratch/project_2007428/users/FAHagenbeek/data/data_UKB_GxE_SESDisease/1000Greferencepanel/ref.bim", data.table = FALSE, header = FALSE)
ref_gen = ref_gen[,c(2,5,6)] # retain relevant columns only
colnames(ref_gen) = c("Predictor","A1","A2") # rename columns

### Add IUPAC allele coding - based on Fiona's code

ref_gen$IUPAC = 0

ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='A'] = 'W'
ref_gen$IUPAC[ref_gen$A1 == 'C' & ref_gen$A2 =='G' | ref_gen$A1 == 'G' & ref_gen$A2 =='C'] = 'S'
ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='G' | ref_gen$A1 == 'G' & ref_gen$A2 =='A'] = 'R'
ref_gen$IUPAC[ref_gen$A1 == 'C' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='C'] = 'Y'
ref_gen$IUPAC[ref_gen$A1 == 'G' & ref_gen$A2 =='T' | ref_gen$A1 == 'T' & ref_gen$A2 =='G'] = 'K'
ref_gen$IUPAC[ref_gen$A1 == 'A' & ref_gen$A2 =='C' | ref_gen$A1 == 'C' & ref_gen$A2 =='A'] = 'M'

### Keep only predictor names and IUPAC allele coding

ref_gen_comparison = ref_gen

ref_gen = ref_gen[,c("Predictor", "IUPAC")]


##### =========================== #####

### Starting prune

##### =========================== #####

for (nfile in 1:length(analysis_list)){

    ### Select file based on loop

    curr_analysis = analysis_list[nfile]

    heading(paste0("Creating PRS sumstat file for: ", basename(curr_analysis)))

    phenotype = strsplit(basename(curr_analysis), "_metal")[[1]][1]

    curr_sumstats = fread(curr_analysis, data.table = FALSE, tmpdir = tempDir)

    ### Convert sumstats to PRS format and reduce to HapMap3+
    ### Also add hg37 positions

    curr_sumstats = curr_sumstats[curr_sumstats$MarkerName %in% hapmap3_snps$rsid,]

    hg38_positions = curr_sumstats[,c("MarkerName", "Position")]
    colnames(hg38_positions) = c("rsid", "hg38_pos")

    curr_sumstats = merge(curr_sumstats, hg37_snps, by = "MarkerName")

    curr_sumstats$se = sqrt( 1 / ((curr_sumstats$Zscore^2) + curr_sumstats$Weight)) / sqrt(2 * curr_sumstats$Freq1 * (1 - curr_sumstats$Freq1))
    curr_sumstats$beta1 = curr_sumstats$Zscore * curr_sumstats$se

    curr_sumstats = data.frame(chr = curr_sumstats$Chromosome,
                                hg37_pos = curr_sumstats$hg37_pos,
                                hg38_pos = curr_sumstats$Position,
                                rsid = curr_sumstats$MarkerName,
                                a1 = curr_sumstats$Allele1,
                                a0 = curr_sumstats$Allele2,
                                n_eff = curr_sumstats$Weight,
                                MAF = curr_sumstats$Freq1,
                                p = curr_sumstats[,"P-value"],
                                beta = curr_sumstats$beta1,
                                se = curr_sumstats$se
                                )

    curr_sumstats$a1 = toupper(curr_sumstats$a1)
    curr_sumstats$a0 = toupper(curr_sumstats$a0)

    cat(paste0("\nCurrent sumstats (", phenotype, ") are these:\n\n"))
    print(head(curr_sumstats))
    print(summary(curr_sumstats))

    if (any(is.infinite(curr_sumstats$se))){

        cat(paste0("\nInf values detected. Removing them.\n\n"))

        numINF = nrow(curr_sumstats[is.infinite(curr_sumstats$se),])

        print(curr_sumstats[is.infinite(curr_sumstats$se),])

        curr_sumstats = curr_sumstats[!(is.infinite(curr_sumstats$se)),]

        cat(paste0("\n", numINF, " rows removed.\n\n"))
        print(summary(curr_sumstats))

    }

    fwrite(curr_sumstats, paste0(outDir, phenotype, "_meta_v4_", hapmap_name, ".txt.gz"), row.names = FALSE, quote = FALSE, sep = "\t", compress = "gzip")

    ### Add Predictor

    curr_sumstats$Predictor = paste0(curr_sumstats$chr, ":", curr_sumstats$hg37_pos)

    ### Add IUPAC to sumstats

    curr_sumstats$IUPAC = 0

    curr_sumstats$IUPAC[curr_sumstats$a1 == 'A' & curr_sumstats$a0 =='T' | curr_sumstats$a1 == 'T' & curr_sumstats$a0 =='A'] = 'W'
    curr_sumstats$IUPAC[curr_sumstats$a1 == 'C' & curr_sumstats$a0 =='G' | curr_sumstats$a1 == 'G' & curr_sumstats$a0 =='C'] = 'S'
    curr_sumstats$IUPAC[curr_sumstats$a1 == 'A' & curr_sumstats$a0 =='G' | curr_sumstats$a1 == 'G' & curr_sumstats$a0 =='A'] = 'R'
    curr_sumstats$IUPAC[curr_sumstats$a1 == 'C' & curr_sumstats$a0 =='T' | curr_sumstats$a1 == 'T' & curr_sumstats$a0 =='C'] = 'Y'
    curr_sumstats$IUPAC[curr_sumstats$a1 == 'G' & curr_sumstats$a0 =='T' | curr_sumstats$a1 == 'T' & curr_sumstats$a0 =='G'] = 'K'
    curr_sumstats$IUPAC[curr_sumstats$a1 == 'A' & curr_sumstats$a0 =='C' | curr_sumstats$a1 == 'C' & curr_sumstats$a0 =='A'] = 'M'

    ### retain only R, Y, K, and M IUPAC coded alleles - to remove ambiguous alleles

    curr_sumstats = curr_sumstats[(curr_sumstats$IUPAC %in% c('R', 'Y', 'K', 'M')),]

    ### Only retain SNPs where IUPAC alleles in summary statistics and reference file match

    megaprs_sumstats = inner_join(curr_sumstats, ref_gen, by = c("Predictor", "IUPAC"))

    megaprs_sumstats = data.frame(Predictor = paste0(megaprs_sumstats$chr, ":", megaprs_sumstats$hg37_pos),
                                A1 = toupper(megaprs_sumstats$a1),
                                A2 = toupper(megaprs_sumstats$a0),
                                n = as.numeric(megaprs_sumstats$n_eff),
                                Z = as.numeric(megaprs_sumstats$beta) / as.numeric(megaprs_sumstats$se))

    megaprs_sumstats = unique(megaprs_sumstats)
    
    predictor_set = data.frame(predictor = megaprs_sumstats$Predictor)

    ### Write summary statistics to file

    fwrite(megaprs_sumstats, paste0(outDir, phenotype, "_meta_v4_", hapmap_name, "_mega_sumstats.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

    fwrite(predictor_set, paste0(outDir, phenotype, "_meta_v4_", hapmap_name, "_mega_snplist.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    cat(paste0("\nWritten PGS sumstats for (", phenotype, ")\n.....\n\n"))

}
