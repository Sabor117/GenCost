##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
library(PathWAS, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(R.utils)
options(scipen = 999)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")
tempDir = paste0(mainDir, "tmpdir/")

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("1")

run_no = as.numeric(args[1])

all_cohorts_dirs = Sys.glob(paste0(sumstatDir, "*/"))

currDir = all_cohorts_dirs[run_no]
cohort_name = toupper(basename(all_cohorts_dirs[run_no]))

### Several cohorts in build 38

hg38_list = c("CHB", "FINNGEN", "MGBB", "QGP", "AOU", "ESTBB", "GNH", "GE", "NTR")

### For build 38 cohorts, read Liftover input file
### For build 37, read the Liftover output file

if (cohort_name %in% hg38_list){

	file_select = Sys.glob(paste0(currDir, "*_liftOver_input*"))

	cat(paste0("\nRunning on ", cohort_name, ", so file is: ", file_select, "\n\n"))

	build38 = TRUE

	cohort = fread(file_select, data.table = FALSE, tmpdir = tempDir)

} else {

	file_select = Sys.glob(paste0(currDir, "*_liftOver_output*"))

	cat(paste0("\nRunning on ", cohort_name, ", so file is: ", file_select, "\n\n"))

	build38 = FALSE

	cohort = fread(file_select, data.table = FALSE, tmpdir = tempDir)

}

### Set up chromosome list to remove weird chromosomes

chr_list = c(paste0("chr", 1:23), "chrX")

cohort$V1 = gsub("chr23", "chrX", cohort$V1)
cohort = cohort[cohort$V1 %in% chr_list,]

### Cohorts with 6 columns have both frequency + rsID data
### 5 columns is just frequency

if (ncol(cohort) == 5){

	colnames(cohort) = c("chr", "pos", "pos1", "full_snpid", "freq")

	rsid_presence = FALSE

} else {

	colnames(cohort) = c("chr", "pos", "pos1", "full_snpid", "freq", "rsid")

	rsid_presence = TRUE

}

### Fix chromosomes AGAIN

cohort$full_snpid = gsub("^23:", "X:", cohort$full_snpid)

### Specifically UKB is the only cohort with rsID + build 37
### Liftover only kept the first digit of the rsID
### Needs to be merged from Input

if (isTRUE(rsid_presence) & !(cohort_name %in% hg38_list)){

	rsid_match_file = Sys.glob(paste0(currDir, "*_liftOver_input*"))

	cat(paste0("\nrsIDs need to be added back in.\n\n"))
	print(head(cohort))

	rsid_match = fread(rsid_match_file, data.table = FALSE, tmpdir = tempDir, select = c(4, 6))
	colnames(rsid_match) = c("full_snpid", "rsid")

	rsid_match$full_snpid = gsub("^23:", "X:", rsid_match$full_snpid)

	cohort$rsid = rsid_match$rsid[match(cohort$full_snpid, rsid_match$full_snpid)]

}

cat(paste0("\nRunning on current file: ", file_select,
			".\nFor cohort: ", cohort_name,
			"\nCHECK", run_no / 2, ".2 COMPLETE.\n\n"))
print(head(cohort))

### Now read in gnoMAD comparison

gnomad_file = "/scratch/project_2007428/data/processing/gnomad/v3/gnomad_v3_b38_ref_nfe.gz"

gnomad = fread(gnomad_file, data.table = FALSE, tmpdir = tempDir, select = c(1:5))

colnames(gnomad) = c("chr", "pos", "gno_ref", "gno_alt", "af_alt")

print(head(gnomad[gnomad$chr == 23,], n = 10))

gnomad$chr = gsub(23, "X", gnomad$chr)
gnomad$chr = gsub("23", "X", gnomad$chr)

print(head(gnomad[gnomad$chr == "X",], n = 10))

gnomad = na.omit(gnomad)
gnomad$snpid = paste0(gnomad$chr, "_", gnomad$pos)

### Need a check for the presence of chrX (trying to be consistent with this)

cat(paste0("\nGnoMAD set of chromosomes:\n"))
print(unique(gnomad$chr))
cat(paste0("\n\n"))

for (i in 1:length(unique(gnomad$chr))){

	test_chr = unique(gnomad$chr)[i]

	print(gnomad[gnomad$chr == test_chr,][1,])

}



##############################

### Running beta correlations

##############################

### Get SNPID string split

cohort_id_split = str_split_fixed(cohort$full_snpid, ":", 4)

### Chr_pos but now with hg38 positions for hg37 files

cohort$chr = gsub("chr23", "chrX", cohort$chr)

cohort$snpid = paste0(gsub("chr", "", cohort$chr), "_", cohort$pos)
cohort$a1 = cohort_id_split[,3]
cohort$a0 = cohort_id_split[,4]

cat(paste0("\nCohort snpid created. CHECK", run_no, ".3 COMPLETE.\n\n"))
print(nrow(cohort))
print(summary(cohort))
print(head(cohort))

### Merge with gnomad

curr_file = inner_join(cohort, gnomad[,c("snpid", "gno_ref", "gno_alt", "af_alt")], by = "snpid")

cat(paste0("\nSumstats and gnoMAD merged successfully. CHECK", run_no, ".4 COMPLETE.\n\n"))
print(nrow(curr_file))
print(summary(curr_file))
print(head(curr_file))

### First check flip status
### 1 : a1 == gno_alt and a0 == gno_ref
### 2 : a1 == gno_ref and a0 == gno_alt
### 0 : a1 and a0 do not match gnoMAD alleles

curr_file$FLIP = flip.data.frame(curr_file[,c("snpid", "a1", "a0", "gno_alt", "gno_ref")])

cat(paste0("\nSumstats and gnoMAD flipped. CHECK", run_no, ".5 COMPLETE.\n\n"))
print(nrow(curr_file))
print(summary(curr_file))
print(head(curr_file))

### For non-matching alleles, remove them

curr_file = curr_file[curr_file$FLIP != 0,]

### There will be duplicates

if (any(duplicated(curr_file$full_snpid))){

	pre_dup_nrow = nrow(curr_file)

	### Finding duplicated full_snpids in allele_flips
	### These have been matched as having the same allele in both directions in gnoMAD
	### Separate allele_flips into duplicated and unique IDs

	allele_flips_dupli_snps = curr_file$full_snpid[which(duplicated(curr_file$full_snpid))]

	allele_flips_duplicated = curr_file[curr_file$full_snpid %in% allele_flips_dupli_snps,]
	allele_flips = curr_file[!(curr_file$full_snpid %in% allele_flips_dupli_snps),]

	cat(paste0(cohort_name, " test for allele_flips.\n\n"))
	print(head(allele_flips))

	### Now parsing the duplicated ones
	### First get rid of completely duplicated lines

	allele_flips_duplicated = unique(allele_flips_duplicated)

	colnames(allele_flips_duplicated)[6] = "file_freq"
	allele_flips_duplicated$file_freq = as.numeric(allele_flips_duplicated$file_freq)

	### Compare AFs and keep closest match
	### I.e. when there is multiple matches of gnoMAD + SNP, with different allele flips, the closest allele flip frequency is used

	output_duplicated_snps_frame = allele_flips_duplicated %>%
	                                  group_by(full_snpid) %>%
	                                  mutate(abs_diff = abs(af_alt - file_freq)) %>%
	                                  filter(abs_diff == min(abs_diff)) %>%
	                                  select(full_snpid, file_freq, af_alt, FLIP) %>%
	                                  ungroup()

	### Some SNPs REMAIN duplicated now - somehow have the same AF for both gnoMAD variations
	### Pick FLIP = 1
	### Repeat process of extraction and adding back in to create a separate table

	allele_flips_dupli_snps_2 = output_duplicated_snps_frame$full_snpid[which(duplicated(output_duplicated_snps_frame$full_snpid))] # List of SNPs which are duplicated

	output_flips_duplicated = output_duplicated_snps_frame[output_duplicated_snps_frame$full_snpid %in% allele_flips_dupli_snps_2,] # Data frame of duplicated SNPs
	output_flips_duplicated = output_flips_duplicated[output_flips_duplicated$FLIP == 1,] # Restrict to FLIP == 1

	output_duplicated_snps_frame = output_duplicated_snps_frame[!(output_duplicated_snps_frame$full_snpid %in% allele_flips_dupli_snps_2),] # Remove duplicated SNPs from output
	output_duplicated_snps_frame = rbind(output_duplicated_snps_frame, output_flips_duplicated) # Add restricted flips back in

	### This new output table is used to select the SNPs from allele_flips_duplicated

	output_duplicated_snps_frame$select = paste0(output_duplicated_snps_frame$full_snpid, ":", output_duplicated_snps_frame$FLIP) # Use SELECT status to confirm flips in original data frame
	allele_flips_duplicated$select = paste0(allele_flips_duplicated$full_snpid, ":", allele_flips_duplicated$FLIP)

	### Showing some of the removed SNPs

	allele_flips_removed = allele_flips_duplicated[!(allele_flips_duplicated$select %in% output_duplicated_snps_frame$select),]

	allele_flips_duplicated = allele_flips_duplicated[allele_flips_duplicated$select %in% output_duplicated_snps_frame$select,]

	if (nrow(allele_flips_removed) != 0){

		cat(paste0("\nExamples of removed SNPs:\n\n"))
		print(head(allele_flips_removed[sample(nrow(allele_flips_removed), 100),], n = 20))
		cat(paste0("\n=============================\n\nThese were removed.\n=============================\n\n"))

	}

	### Bind all of this back with allele_flips
	### So that flip status can be merged with the sumstats

	colnames(allele_flips_duplicated)[6] = "freq"

	cat(paste0(cohort_name, " pre-bind test for allele_flips.\n\n"))
	print(head(allele_flips))
	cat(paste0(cohort_name, " pre-bind test for allele_flips_duplicated.\n\n"))
	print(head(allele_flips_duplicated))

	if (nrow(allele_flips_duplicated) > 0) {

		allele_flips = rbind(allele_flips, allele_flips_duplicated[,colnames(allele_flips)])

	}

	post_dup_nrow = nrow(allele_flips)

	### Doing checks of duplications

	#cat(paste0("\nDuplicated full_snpids removed. Nrow = ", post_dup_nrow, ", ", pre_dup_nrow - post_dup_nrow, " SNPs removed.\n\n"))
	print(head(allele_flips))
	print(summary(allele_flips))
	#cat(paste0("\n=====\n"))
	#cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
	#          any(duplicated(allele_flips$full_snpid)),
	#          "\n Are there any duplicated snpids? ",
	#          any(duplicated(allele_flips$snpid)), "\n\n"))
	#cat(paste0("\n=====\n"))

} else {

	allele_flips = curr_file

	cat(paste0(cohort_name, " no duplicates test for allele_flips.\n\n"))
	print(head(allele_flips))

	pre_dup_nrow = nrow(allele_flips)
	post_dup_nrow = nrow(allele_flips)

	#cat(paste0("\nDuplicated full_snpids removed. Nrow = ", nrow(allele_flips), ", No SNPs removed.\n\n"))
	print(head(allele_flips))
	print(summary(allele_flips))
	#cat(paste0("\n=====\n"))
	#cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
	#          any(duplicated(allele_flips$full_snpid)),
	#          any(duplicated(allele_flips$snpid)), "\n\n"))
	#          "\n Are there any duplicated snpids? ",
	#cat(paste0("\n=====\n"))

}

if (isFALSE(rsid_presence)){

	out_file = allele_flips[,c("snpid", "full_snpid", "chr", "pos", "a1", "a0", "gno_ref", "gno_alt", "af_alt", "freq", "FLIP")]

} else {

	out_file = allele_flips[,c("snpid", "full_snpid", "rsid", "chr", "pos", "a1", "a0", "gno_ref", "gno_alt", "af_alt", "freq", "FLIP")]

}

out_file$chr = gsub("chr", "", out_file$chr)

cat(paste0("\nFlip column added. CHECK", run_no, ".6 COMPLETE.\n\n"))
print(nrow(out_file))
print(summary(out_file))
print(head(out_file))

out_file_name = paste0(dirname(file_select), "/", cohort_name, "_alleles_aligned.txt.gz")

fwrite(out_file, out_file_name,
        row.names = FALSE, quote = FALSE, sep = "\t", na = "NA", compress = "gzip")





