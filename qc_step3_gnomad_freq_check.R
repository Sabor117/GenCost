##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
library(PathWAS, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(ggplot2)
options(scipen = 999)
set.seed(117)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")

all_cohorts_files = c(paste0(sumstatDir, "agds/AGDS.LIND.ALL.All.ALL.EUR.11161.SAIGE.20230504.txt.gz"),
						paste0(sumstatDir, "chb/CHB.LOULOUDIS.INPATIENTS.ALL.ALL.EUR.191803.REGENIE.20231211.tsv.gz"),
						paste0(sumstatDir, "finngen/FINNGEN.LEE.IN.5_95.FM.EUR.483907.REGENIE.20231014.txt.gz"),
						paste0(sumstatDir, "gbp/GBP.LIND.ALL.All.ALL.EUR.2833.SAIGE.20230504.txt.gz"),
						paste0(sumstatDir, "gs20k/GS20K.RICHMOND.ALL.ALL.ALL.EUR.14750.BOLT.20230804.txt.gz"),
						paste0(sumstatDir, "mgbb/MGBB.LEE.IN.ALL.ALL.EUR.43862.REGENIE.20231210.txt.gz"),
						paste0(sumstatDir, "ntr/NTR.LAAN.IN.15_102.ALL.EUR.14572.fastGWA.20210702.txt.gz"),
						paste0(sumstatDir, "ukb/UKB.SMW.IN.ALL.ALL.EUR.402248.REGENIE.20240305.txt.gz"),
						paste0(sumstatDir, "ukb/UKB.SMW.IN.ALL.ALL.AFR.7204.REGENIE.20240305.txt.gz"),
						paste0(sumstatDir, "ukb/UKB.SMW.IN.ALL.ALL.CSA.9798.REGENIE.20240305.txt.gz"),
						paste0(sumstatDir, "ukb/UKB.SMW.IN.ALL.ALL.EAS.2348.REGENIE.20240305.txt.gz"),
						paste0(sumstatDir, "qgp/QGP.Mbarek.IN.ALL.ALL.MID.1.SAIGE.20231214.txt.gz"),
						paste0(sumstatDir, "aou/AOU.LEE.IN.ALL.ALL.EUR.22493.REGENIE.20231223.txt.gz"),
						paste0(sumstatDir, "estbb/EstBB.Abner.INPATIENT.ALL.ALL.EUR.133544.REGENIE.20231128.txt.gz"),
						paste0(sumstatDir, "gnh/GNH.SK.IN.ALL.ALL.CSA.51167.REGENIE.20240107.txt.gz"),
						paste0(sumstatDir, "ge/GE.KZ.IN.ALL.ALL.EUR.16135.SAIGE.20240822.txt.gz")
						)


args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
	print("Script obtained arguments from command line incorrectly. Please enter run number.")
	stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("1")

run_no = as.numeric(args[1])

gnomad_file = "/scratch/project_2007428/data/processing/gnomad/v3/gnomad_v3_b38_ref_nfe.gz"

gnomad = fread(gnomad_file, data.table = FALSE)

gnomad = gnomad[,c(1:5)]
colnames(gnomad) = c("chr", "pos", "gno_ref", "gno_alt", "af_alt")

#gnomad = gnomad[gnomad$af_alt >= 0.0001,]
gnomad = na.omit(gnomad)
gnomad$snpid = paste0(gnomad$chr, "_", gnomad$pos)

##############################

### Running frequncy correlations

##############################

freq_corrs = data.frame(matrix(ncol = 2, nrow = 0))
colnames(freq_corrs) = c("cohort", "gnomad_freq_correlation")

cat(paste0("\n\n==================\n\nRunning on ", basename(all_cohorts_files[run_no]), ". CHECK", run_no, ".1 COMPLETE.\n\n"))

cohort = strsplit(basename(all_cohorts_files[run_no]), "\\.")[[1]][1]

chrom_col = case_when(cohort == "UKB" ~ "CHROM",
						cohort == "GNH" ~ "CHROM",
						.default = "CHR")
                        
pos_col = case_when(cohort == "UKB" ~ "GENPOS",
						cohort == "QGP" ~ "POS",
						cohort == "EstBB" ~ "POS",
						cohort == "GNH" ~ "GENPOS",
						.default = "BP")

a1_col = case_when(cohort == "UKB" ~ "ALLELE1",
					cohort == "CHB" ~ "ALLELE1",
					cohort == "QGP" ~ "Allele2",
					cohort == "EstBB" ~ "ALLELE1",
					cohort == "GNH" ~ "ALLELE1",
					.default = "A1")

a2_col = case_when(cohort == "UKB" ~ "ALLELE0",
					cohort == "CHB" ~ "ALLELE0",
					cohort == "QGP" ~ "Allele1",
					cohort == "EstBB" ~ "ALLELE0",
					cohort == "GNH" ~ "ALLELE0",
					.default = "A2")
                    
freq_col = case_when(cohort == "UKB" ~ "A1FREQ",
						cohort == "CHB" ~ "A1FREQ",
						cohort == "QGP" ~ "AF_Allele2",
						cohort == "EstBB" ~ "A1FREQ",
						cohort == "GNH" ~ "A1FREQ",
						.default = "EAF")

if (cohort != "QGP"){

	curr_file = fread(all_cohorts_files[run_no], data.table = FALSE, select = c(chrom_col, pos_col, a1_col, a2_col, freq_col))

} else {

	curr_file = fread(all_cohorts_files[run_no], data.table = FALSE)

	realColumns = colnames(curr_file)[-(which(colnames(curr_file) == "SNPID"))]

	curr_file[,ncol(curr_file)] = NULL

	colnames(curr_file) = realColumns

	print(head(curr_file))

	curr_file = curr_file[,c(chrom_col, pos_col, a1_col, a2_col, freq_col)]

}


cat(paste0("\nBase file has ", nrow(curr_file), " rows.\n\n"))

### Copenhagen data contains ! alleles and Qatar contains * - these should be removed

curr_file = curr_file[curr_file[,a1_col] != "!",]
curr_file = curr_file[curr_file[,a2_col] != "!",]

curr_file = curr_file[curr_file[,a1_col] != "*",]
curr_file = curr_file[curr_file[,a2_col] != "*",]

cat(paste0("\n'!' and '*' alleles removed. File has ", nrow(curr_file), " rows.\n\n"))

### NTR EAF column is not numeric - fix it with this, however this shouldn't affect the others

curr_file[,freq_col] = as.numeric(curr_file[,freq_col])

if (any(is.na(curr_file[,freq_col]))){

	curr_file = curr_file[!(is.na(curr_file[,freq_col])),]

	cat(paste0("\nNA allele freq SNPs removed. File has ", nrow(curr_file), " rows.\n\n"))

}

### Liftover is not necessary for Copenhagen/FINNGEN/MGBB
### These sumstats are in Build 38

hg38_list = c("CHB", "FINNGEN", "MGBB", "QGP", "AOU", "EstBB", "GNH", "GE")

if (cohort %in% hg38_list){

	cat(paste0("\nCohort is ", cohort, ". Do not need LiftOver steps. CHECK", run_no, ".2+3 COMPLETE.\n\n"))

	curr_file$snpid = paste0(curr_file[,chrom_col], "_", curr_file[,pos_col])

	curr_file = merge(curr_file[,c("snpid", a1_col, a2_col, freq_col)], gnomad[,c("snpid", "gno_ref", "gno_alt", "af_alt")], by = "snpid")

	cat(paste0("\nFile merged with gnoMAD. File has ", nrow(curr_file), " rows.\n\n"))
	cat(paste0("\nSumstats and gnoMAD merged successfully. CHECK", run_no, ".4 COMPLETE.\n\n"))
	print(head(curr_file))

	curr_file$FLIP = flip.data.frame(curr_file[,c("snpid", a1_col, a2_col, "gno_alt", "gno_ref")])

	### If FLIP == 0 it is due to difference in alleles
	### E.g.: G/A vs G/GTA
	### Actually not possible to compare this?

	curr_file = curr_file[curr_file$FLIP != 0,]

	cat(paste0("\nSumstats and gnoMAD flipped. CHECK", run_no, ".4 COMPLETE.\n\n"))
	cat(paste0("\n'0' flip status alleles removed. File has ", nrow(curr_file), " rows.\n\n"))
	print(head(curr_file))

	multi_allelic_snps = unique(curr_file$snpid[which(duplicated(curr_file$snpid))])

	### Finding duplicated snpids
	### These have been matched as having the same allele in both directions in gnoMAD
	### Separate curr_file into duplicated and unique IDs

	curr_file_duplicated = curr_file[curr_file$snpid %in% multi_allelic_snps,]
	curr_file = curr_file[!(curr_file$snpid %in% multi_allelic_snps),]

	### Now parsing the duplicated ones
	### First get rid of completely duplicated lines

	curr_file_duplicated = unique(curr_file_duplicated)

	### Compare AFs and keep closest match
	### I.e. when there is multiple matches of gnoMAD + SNP, with different allele flips, the closest allele flip frequency is used

	colnames(curr_file_duplicated)[4] = "file_freq" 

	output_curr_file_duplicated = curr_file_duplicated %>%
										group_by(snpid) %>%
										mutate(abs_diff = abs(af_alt - file_freq)) %>%
										filter(abs_diff == min(abs_diff)) %>%
										select(snpid, file_freq, af_alt, FLIP) %>%
										ungroup()

	### Remove those with double appearances of same AF

	multi_allelic_snps_2 = output_curr_file_duplicated$snpid[which(duplicated(output_curr_file_duplicated$snpid))]

	output_curr_file_duplicated = output_curr_file_duplicated[!(output_curr_file_duplicated$snpid %in% multi_allelic_snps_2),]

	### This new output table is used to select the SNPs from allele_flips_duplicated

	output_curr_file_duplicated$select = paste0(output_curr_file_duplicated$snpid, ":", output_curr_file_duplicated$FLIP)
	curr_file_duplicated$select = paste0(curr_file_duplicated$snpid, ":", curr_file_duplicated$FLIP)

	curr_file_duplicated = curr_file_duplicated[curr_file_duplicated$select %in% output_curr_file_duplicated$select,]

	colnames(curr_file_duplicated)[4] = freq_col
	curr_file_duplicated$select = NULL

	### Add it back in

	curr_file = rbind(curr_file, curr_file_duplicated)

	### Check by full_snpid
	### Remove rows which have multiple different frequencies for the same SNP
	### In either file

	curr_file$full_snpid = paste(curr_file$snpid, curr_file[,a1_col], curr_file[,a2_col], sep = ":")

	multi_allelic_snps_3 = unique(curr_file$full_snpid[which(duplicated(curr_file$full_snpid))])

	curr_file = curr_file[!(curr_file$full_snpid %in% multi_allelic_snps_3),]

	### Perform flipping comparison

	curr_file$freq1 = ifelse(curr_file$FLIP == -1, 1 - curr_file[,freq_col], curr_file[,freq_col])

	cat(paste0("\nTrue frequency calculated. CHECK", run_no, ".5 COMPLETE.\n\n"))
	print(head(curr_file))

	curr_cor = cor(curr_file$freq1, curr_file$af_alt)

	cat(paste0("\nCorrelation calculated for ", cohort, " = ", curr_cor, ". CHECK", run_no, ".6 COMPLETE.\n\n==================\n\n"))

	outrow = data.frame(cohort = cohort, gnomad_freq_correlation = curr_cor)

} else {

	liftOver_file = fread(paste0(dirname(all_cohorts_files[run_no]), "/", cohort, "_liftOver_output.out"), data.table = FALSE)
	colnames(liftOver_file) = c("chr", "pos", "pos1", "snpid")

	cat(paste0("\nSumstats file and LiftOver file read. CHECK", run_no, ".2 COMPLETE.\n\n"))
	cat(paste0("\nLiftover file has ", nrow(liftOver_file), " rows.\n\n"))

	curr_file$snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col])

	curr_file = merge(curr_file, liftOver_file[,c("pos", "snpid")], by = "snpid")

	cat(paste0("\nFile merged with LiftOver. File has ", nrow(curr_file), " rows.\n\n"))
	cat(paste0("\nSumstats file and LiftOver file merged successfully. CHECK", run_no, ".3 COMPLETE.\n\n"))
	print(head(curr_file))

	curr_file$snpid = paste0(curr_file[,chrom_col], "_", curr_file$pos)

	curr_file = merge(curr_file[,c("snpid", a1_col, a2_col, freq_col)], gnomad[,c("snpid", "gno_ref", "gno_alt", "af_alt")], by = "snpid")

	cat(paste0("\nFile merged with gnoMAD. File has ", nrow(curr_file), " rows.\n\n"))
	cat(paste0("\nSumstats and gnoMAD merged successfully. CHECK", run_no, ".4 COMPLETE.\n\n"))
	print(head(curr_file))

	curr_file$FLIP = flip.data.frame(curr_file[,c("snpid", a1_col, a2_col, "gno_alt", "gno_ref")])

	### If FLIP == 0 it is due to difference in alleles
	### E.g.: G/A vs G/GTA
	### Actually not possible to compare this?

	curr_file = curr_file[curr_file$FLIP != 0,]

	cat(paste0("\nSumstats and gnoMAD flipped. CHECK", run_no, ".4 COMPLETE.\n\n"))
	cat(paste0("\n'0' flip status alleles removed. File has ", nrow(curr_file), " rows.\n\n"))
	print(head(curr_file))

	multi_allelic_snps = unique(curr_file$snpid[which(duplicated(curr_file$snpid))])

	### Finding duplicated snpids
	### These have been matched as having the same allele in both directions in gnoMAD
	### Separate curr_file into duplicated and unique IDs

	curr_file_duplicated = curr_file[curr_file$snpid %in% multi_allelic_snps,]
	curr_file = curr_file[!(curr_file$snpid %in% multi_allelic_snps),]

	### Now parsing the duplicated ones
	### First get rid of completely duplicated lines

	curr_file_duplicated = unique(curr_file_duplicated)

	### Compare AFs and keep closest match
	### I.e. when there is multiple matches of gnoMAD + SNP, with different allele flips, the closest allele flip frequency is used

	colnames(curr_file_duplicated)[4] = "file_freq" 

	output_curr_file_duplicated = curr_file_duplicated %>%
										group_by(snpid) %>%
										mutate(abs_diff = abs(af_alt - file_freq)) %>%
										filter(abs_diff == min(abs_diff)) %>%
										select(snpid, file_freq, af_alt, FLIP) %>%
										ungroup()

	### Remove those with double appearances of same AF

	multi_allelic_snps_2 = output_curr_file_duplicated$snpid[which(duplicated(output_curr_file_duplicated$snpid))]

	output_curr_file_duplicated = output_curr_file_duplicated[!(output_curr_file_duplicated$snpid %in% multi_allelic_snps_2),]

	### This new output table is used to select the SNPs from allele_flips_duplicated

	output_curr_file_duplicated$select = paste0(output_curr_file_duplicated$snpid, ":", output_curr_file_duplicated$FLIP)
	curr_file_duplicated$select = paste0(curr_file_duplicated$snpid, ":", curr_file_duplicated$FLIP)

	curr_file_duplicated = curr_file_duplicated[curr_file_duplicated$select %in% output_curr_file_duplicated$select,]

	colnames(curr_file_duplicated)[4] = freq_col
	curr_file_duplicated$select = NULL

	### Add it back in

	curr_file = rbind(curr_file, curr_file_duplicated)

	### Check by full_snpid
	### Remove rows which have multiple different frequencies for the same SNP
	### In either file

	curr_file$full_snpid = paste(curr_file$snpid, curr_file[,a1_col], curr_file[,a2_col], sep = ":")

	multi_allelic_snps_3 = unique(curr_file$full_snpid[which(duplicated(curr_file$full_snpid))])

	curr_file = curr_file[!(curr_file$full_snpid %in% multi_allelic_snps_3),]

	### Perform flipping comparison

	curr_file$freq1 = ifelse(curr_file$FLIP == -1, 1 - curr_file[,freq_col], curr_file[,freq_col])

	cat(paste0("\nTrue frequency calculated. CHECK", run_no, ".5 COMPLETE.\n\n"))
	print(head(curr_file))

	curr_cor = cor(curr_file$freq1, curr_file$af_alt)

	cat(paste0("\nCorrelation calculated for ", cohort, " = ", curr_cor, ". CHECK", run_no, ".6 COMPLETE.\n\n==================\n\n"))

	if (cohort == "UKB"){

		cohort = paste0(strsplit(basename(all_cohorts_files[run_no]), "\\.")[[1]][1], "_", strsplit(basename(all_cohorts_files[run_no]), "\\.")[[1]][6])

	}

	outrow = data.frame(cohort = cohort, gnomad_freq_correlation = curr_cor)

}

freq_corrs = rbind(freq_corrs, outrow)

fwrite(freq_corrs, paste0(mainDir, paste0("outputs/cohort_sumstats/gnomad_frequency_correlations_", cohort, ".txt")), row.names = FALSE, quote = FALSE, sep = "\t")

curr_file = curr_file[sample(nrow(curr_file), nrow(curr_file)/100),]

cat(paste0("\n1 in 100 SNPs selected. File has ", nrow(curr_file), " rows.\n\n"))
cat(paste0("\nData sampled for. CHECK", run_no, ".7 COMPLETE.\n\n==================\n\n"))

af_plotDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/af_plots/"

png(paste0(af_plotDir, cohort, "_gnomad_af_comparison.png"), width = 1000, height = 1000)

ggplot(curr_file, aes(x = af_alt, y = freq1)) +
    geom_point() +
    theme_linedraw() + 
    labs(title = "AF comparison", x = "gnoMAD frequency", y = paste0(cohort, " A1 frequency")) +
    theme(axis.text.x = element_text(size = 12),   # Change x-axis label size
            axis.text.y = element_text(size = 12),   # Change y-axis label size
            axis.title.x = element_text(size = 14),  # Change x-axis title size
            axis.title.y = element_text(size = 14),  # Change y-axis title size
            plot.title = element_text(size = 16)) +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_abline(intercept = 0, slope = 1, color = "red")

dev.off()

cat(paste0("\nAF plot made. CHECK", run_no, ".8 COMPLETE.\n\n==================\n\n"))





