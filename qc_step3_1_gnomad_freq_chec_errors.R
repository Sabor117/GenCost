##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
library(PathWAS, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(ggplot2)
set.seed(117)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")

all_cohorts_files = c(paste0(sumstatDir, "copenhagen/CHB.LOULOUDIS.ALL.ALL.ALL.ALL.EUR.222466.REGENIE.20230704.txt.gz"),
                        paste0(sumstatDir, "finngen/FINNGEN.LEE.ALL.5_95.F.EUR.275659.REGENIE.20231014.txt.gz"),
                        paste0(sumstatDir, "mgbb/MGBB.LEE.IN.19_95.M.EUR.4008.REGENIE.20231005.txt.gz")
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

gnomad = fread(gnomad_file, data.table = FALSE, nrow = 5000000)

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

chrom_col = ifelse(cohort == "UKB",  "CHROM", "CHR")
pos_col = ifelse(cohort == "UKB",  "GENPOS", "BP")
a1_col = ifelse(cohort == "CHB" | cohort == "UKB",  "ALLELE1", "A1")
a2_col = ifelse(cohort == "CHB" | cohort == "UKB",  "ALLELE0", "A2")
freq_col = ifelse(cohort == "CHB" | cohort == "UKB",  "A1FREQ", "EAF")

curr_file = fread(all_cohorts_files[run_no], data.table = FALSE, select = c(chrom_col, pos_col, a1_col, a2_col, freq_col))

cat(paste0("\nBase file has ", nrow(curr_file), " rows.\n\n"))

### Copenhagen data contains ! alleles - these should be removed

curr_file = curr_file[curr_file[,a1_col] != "!",]
curr_file = curr_file[curr_file[,a2_col] != "!",]

cat(paste0("\n'!' alleles removed. File has ", nrow(curr_file), " rows.\n\n"))

### NTR EAF column is not numeric - fix it with this, however this shouldn't affect the others

curr_file[,freq_col] = as.numeric(curr_file[,freq_col])

if (any(is.na(curr_file[,freq_col]))){

  curr_file = curr_file[!(is.na(curr_file[,freq_col])),]

  cat(paste0("\nNA allele freq SNPs removed. File has ", nrow(curr_file), " rows.\n\n"))

}

### Liftover is not necessary for Copenhagen/FINNGEN/MGBB
### These sumstats are in Build 38

cat(paste0("\nCohort is ", cohort, ". Do not need LiftOver steps. CHECK", run_no, ".2+3 COMPLETE.\n\n"))

curr_file$snpid = paste0(curr_file[,chrom_col], "_", curr_file[,pos_col])

curr_file = merge(curr_file[,c("snpid", a1_col, a2_col, freq_col)], gnomad[,c("snpid", "gno_ref", "gno_alt", "af_alt")], by = "snpid")

cat(paste0("\nFile merged with gnoMAD. File has ", nrow(curr_file), " rows.\n\n"))
cat(paste0("\nSumstats and gnoMAD merged successfully. CHECK", run_no, ".4 COMPLETE.\n\n"))
print(head(curr_file))

curr_file$FLIP = flip.data.frame(curr_file[,c("snpid", a1_col, a2_col, "gno_alt", "gno_ref")])

cat(paste0("\nSumstats and gnoMAD flipped. CHECK", run_no, ".5 COMPLETE.\n\n"))
cat(paste0("\nFile still has ", nrow(curr_file), " rows.\n\n"))
print(head(curr_file))

### If FLIP == 0 it is due to difference in alleles
### E.g.: G/A vs G/GTA
### Actually not possible to compare this?

curr_file = curr_file[curr_file$FLIP != 0,]

cat(paste0("\n'0' flip status alleles removed. File has ", nrow(curr_file), " rows.\n\n"))

curr_file$freq1 = ifelse(curr_file$FLIP == -1, 1 - curr_file[,freq_col], curr_file[,freq_col])

cat(paste0("\nTrue frequency calculated. CHECK", run_no, ".6 COMPLETE.\n\n"))
print(head(curr_file))

curr_file_output = curr_file[curr_file$freq1 >= 0.77 & curr_file$freq1 <= 0.90,]
curr_file_output = curr_file_output[curr_file_output$af_alt >= 0.1 & curr_file_output$af_alt <= 0.22,]

fwrite(curr_file_output, paste0(mainDir, paste0("tmpdir/af_comparison_frame_", cohort, ".txt")), row.names = FALSE, quote = FALSE, sep = "\t")

cat(paste0("\nAF plot made. CHECK", run_no, ".8 COMPLETE.\n\n==================\n\n"))





