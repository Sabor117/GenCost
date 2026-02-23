## ---------------------------
##
## Script name: meta_step11_coloc_run.R
##
## Purpose of script: Run coloc for gencost data and curated list of GWASs
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-04-30
##
## Adapted from coloc_run.R by Dr. Tomoko Nakanishi
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Run once
# install.packages("coloc", lib="/projappl/project_2007428/RPackages_421/")

### Packages

library(optparse)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(coloc, lib = "/projappl/project_2007428/RPackages_421/")
library(arrow)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
processDir = paste0(mainDir, "processing/phewas_gwas_sumstats/")

### Read from command line

option_list = list(
	make_option(c('--pheno', '-p' ), help = "phenotype", default = ""),
	make_option(c('--run_no', '-n'), help = "the number of phenotype in the GWAS list", default = 1),
	make_option(c('--out_dir', '-d'), help = "*output prefix", default = "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/coloc_MAF_0_001/"),
	make_option(c('--out_prefix'), help = "*output dir", default = "coloc_res")
)

option.parser = OptionParser(option_list=option_list)
opt = parse_args(option.parser)

pheno = opt$pheno
run_no = opt$run_no
outDir = opt$out_dir
out_prefix = opt$out_prefix

### FOR TESTING

manual_test = FALSE

if (isTRUE(manual_test)){

	pheno = "IN"
	run_no = 52 # 20 = BMI, 35 = MVP Cardiac, 27 = ALT
	outDir = paste0(mainDir, "outputs/coloc_MAF_0_001/")
	out_prefix = "coloc_res"

	gwas_phenotype_file = paste0(mainDir, "outputs/METAL_v4/pruned_meta/", pheno, "_ALL_metal_output_1_pruned.TBL.gz")

}


### Creating out directory for GWAS phenotype

if (!(file.exists(outDir))){

	dir.create(file.path(outDir))

}

### Files

phewas_summary_list = paste0(processDir, "phewas_summary_table.txt")
variant_list_file = paste0(mainDir, "outputs/table/messy_top_variants_ALL_analysis_MAF_0_001_250kb.tsv")
gwas_phenotype_file = paste0(mainDir, "outputs/METAL_v4/", pheno, "_ALL_metal_output_1.TBL.gz")
gwas_position_map_file = paste0(mainDir, "outputs/METAL_v4/ALL_METAL_output_SNPs_mapped.txt.gz")

### Read files

summary_list = fread(phewas_summary_list, data.table = FALSE)
variant_list = fread(variant_list_file, data.table = FALSE)
gwas_pheno = fread(gwas_phenotype_file, data.table = FALSE, tmpdir = paste0(mainDir, "tmpdir/"))

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Might go somewhere sunny. Sit on beach, look at ocean, collect sea shells. Might run tests on the sea shells.")

sessionInfo()
start_time = Sys.time()


##### =========================== #####

### Pre-processing

##### =========================== #####

### List of GWAS variants - reduce to selected phenotype

sig_variants = variant_list %>% filter(phenotype == pheno)

cat(paste0("\n\nPrinting: head(sig_variants)\n\n.........\n\n"))
print(head(sig_variants))
summary(sig_variants)

cat(paste0("\nPrinting: all(sig_variants$group == 'ALL')\n\n"))
print(all(sig_variants$group == 'ALL'))

### Define analysis type

analysis_type = summary_list$analysis[run_no]

### Modulating variant columns to clarify hg37/38 positions

colnames(sig_variants)[which(colnames(sig_variants) == "bp")] = "hg37_pos"

### SUPERCEDED/OUTDATED CODE
### Table now includes alleles and position data from start

#pos_bp = str_split_fixed(sig_variants$hg38_snpid, "_", 4)
#
#sig_variants$a1 = pos_bp[,3]
#sig_variants$a0 = pos_bp[,4]
#sig_variants$hg38_pos = pos_bp[,2]

### Selection of query trait - based on command line
### Read and then remove "chr" from chromosomes where necessary

query_trait_file = Sys.glob(paste0(processDir, "munged_input/", summary_list$file[run_no], "*"))

query_trait = fread(query_trait_file, data.table = FALSE, tmpdir = paste0(mainDir, "tmpdir/"))

### Fix MAF of query trait files

query_trait$MAF = as.numeric(query_trait$MAF)
query_trait = query_trait[query_trait$MAF > 0 & query_trait$MAF < 1,]

### Some query files may contains "chr" from munging
### Otherwise ensure it is numeric

if (any(grepl("chr", query_trait$chr))){

	query_trait = query_trait %>% mutate(chr = gsub("chr", "", chr))

	query_trait$chr = gsub("X", 23, query_trait$chr)
	query_trait$chr = as.numeric(query_trait$chr)

} else {

	query_trait$chr = gsub("X", 23, query_trait$chr)

	query_trait$chr = as.numeric(query_trait$chr)

}

### Initialise output data frmae

result = data.frame(matrix(nrow = nrow(sig_variants), ncol = 11))
colnames(result) = c("locus_snp", "MarkerName", "NSNP", "NPHENO", "NQUERY", "PP0", "PP1", "PP2", "PP3", "PP4", "minP_outcome")

### Make sure alleles are correct and add beta + ses

gwas_pheno$Allele1 = toupper(gwas_pheno$Allele1)
gwas_pheno$Allele2 = toupper(gwas_pheno$Allele2)

gwas_pheno = gwas_pheno %>% mutate(beta = Zscore / sqrt((2 * Freq1) * (1 - Freq1) * (Weight + (Zscore^2))),
									se = 1 / sqrt((2 * Freq1) * (1 - Freq1) * (Weight + (Zscore^2))))

### Remove silly allele frequencies

gwas_pheno = gwas_pheno[gwas_pheno$Freq1 > 0,]
gwas_pheno = gwas_pheno[gwas_pheno$Freq1 < 1,]

### Replace X with 23

gwas_pheno$Chromosome = gsub("X", 23, gwas_pheno$Chromosome)
gwas_pheno$Chromosome = as.numeric(gwas_pheno$Chromosome)

### Setting up table for duplicated SNPs in GWAS

dupli_gwas_snps = as.data.frame(matrix(ncol = ncol(gwas_pheno), nrow = 0))
colnames(dupli_gwas_snps) = colnames(gwas_pheno)

### For build 37 GWAS queries extra steps are required

if (summary_list$build_hg[run_no] == 19){

	### Read map and merge it with GWAS - potentially memory intensive, but probably simplest way

	gwas_position_map = fread(gwas_position_map_file, data.table = FALSE, tmpdir = paste0(mainDir, "tmpdir/"), select = c("MarkerName", "hg37_pos"))

	gwas_pheno = left_join(gwas_pheno, gwas_position_map, by = "MarkerName")

	### RENAME COLUMNS - Later stages rely on "Position" column and re-doing the script would be a pain

	colnames(gwas_pheno)[which(colnames(gwas_pheno) == "Position")] = "hg38_pos"
	colnames(gwas_pheno)[which(colnames(gwas_pheno) == "hg37_pos")] = "Position"

	### Remove the big file to save memory

	rm(list = c("gwas_position_map"))

}

### Test print

cat(paste0("\nStarting COLOC run. GWAS phenotype = ", pheno,
			"\nRun# = ", run_no, 
			"\nQuery trait = ", summary_list$name[run_no],
			"\nInput file = ", processDir, "munged_input/", summary_list$file[run_no], ".tsv.gz",
			"\nOutput file = ", outDir, "coloc_", pheno, "_", summary_list$name[run_no], ".tsv",
			"\n-----\n"))

print(summary_list[run_no,])

cat("\nHead of GWAS file:\n\n")
print(head(gwas_pheno))
summary(gwas_pheno)

cat("\nHead of query trait:\n\n")
print(head(query_trait))
summary(query_trait)

heading("Starting COLOC loop.")


##### =========================== #####

### Begin COLOC loop

##### =========================== #####

for (i in 1:nrow(result)){

	### Analysis locus defined by sentinel SNP

	curr_analysis_locus = sig_variants[i,]

	### Sentinel SNP position defined based on genome build (37 or 38)
	### "19" is INCORRECT nomenlcature and refers to hg37

	pos_col = case_when(summary_list$build_hg[run_no] == 19 ~ "position_37",
						summary_list$build_hg[run_no] == 38 ~ "position_38",
						.default = "pos")

	cat(paste0("\n", i, "/", nrow(result), ":\n\nStart COLOC on locus ",
					sig_variants$SNP[i], " and ", summary_list$name[run_no], "\n===\n\n"))

	cat("Memory used before iteration", i, ":", pryr::mem_used() / (1024^3), "GB\n\n")

	### Reduce GWAS and query sumstats to 1 MB window around sentinel SNP

	temp_gwas = gwas_pheno %>%
				filter(Chromosome == curr_analysis_locus$Chr &
						Position >= as.numeric(curr_analysis_locus[,pos_col]) - 500000 &
						Position <= as.numeric(curr_analysis_locus[,pos_col]) + 500000)

	temp_query = query_trait %>% filter(chr == curr_analysis_locus$Chr &
											pos >= as.numeric(curr_analysis_locus[,pos_col]) - 500000 &
											pos <= as.numeric(curr_analysis_locus[,pos_col]) + 500000)

	### Merge into analysis SNP set

	analysis_set = inner_join(temp_gwas, temp_query, by = c("Chromosome" = "chr", "Position" = "pos"))

	### Positions/SNPs can be duplictaed between datasets
	### Set up a method of tracking these

	if (any(duplicated(analysis_set$Position))){

		### Define duplicates based on position

		dupli_positions = unique(analysis_set$Position[which(duplicated(analysis_set$Position))])

		### Duplicated query SNPs usually due to multi-allelic SNPs
		### E.g. 1_36316763_C_CT vs 1_36316763_C_CTT
		### Just remove these

		duplicated_query_positions = temp_query[temp_query$pos %in% temp_query$pos[temp_query$pos %in% dupli_positions],]
		duplicated_query_positions = unique(duplicated_query_positions$pos[which(duplicated(duplicated_query_positions$pos))])
		temp_query = temp_query[!(temp_query$pos %in% duplicated_query_positions),]

		### Duplicated SNPs in GWAS can be due to multiple reasons:
		### Multi-allelic SNPs as well
		### Also multiple names for same SNP (e.g. rs112706077 and rs61069669)
		### Keep list of SNPs which are a problem?

		duplicated_gwas_positions = temp_gwas[temp_gwas$Position %in% temp_gwas$Position[temp_gwas$Position %in% dupli_positions],]
		duplicated_gwas_positions = unique(duplicated_gwas_positions$Position[which(duplicated(duplicated_gwas_positions$Position))])

		intermediate_gwas_positions = temp_gwas[temp_gwas$Position %in% duplicated_gwas_positions,]

		temp_gwas = temp_gwas[!(temp_gwas$Position %in% duplicated_gwas_positions),]

		dupli_gwas_snps = rbind(dupli_gwas_snps, intermediate_gwas_positions)

		analysis_set = inner_join(temp_gwas, temp_query, by = c("Chromosome" = "chr", "Position" = "pos"))

	}

	### Remove any NA values

	analysis_set = analysis_set[complete.cases(analysis_set[, c("MarkerName", "Freq1", "beta.x", "beta.y", "MAF", "se", "varbeta")]),]

	### Running COLOC if analysis set actually has SNPs

	cat("\nNumber of GWAS SNPs: ", nrow(temp_gwas),
		"\nNumber of query SNPs: ", nrow(temp_query),
		"\nNumber of overlapping SNPs: ", nrow(analysis_set),
		"\nHead of current ", i, "/", nrow(result), "analysis file:\n\n")
	print(head(analysis_set))
	summary(analysis_set)
	
	if(nrow(analysis_set) > 1) {

		if (analysis_type == "quant"){

			coloc.res = coloc.abf(dataset1 = list(varbeta = (analysis_set$se ^ 2),
													snp = analysis_set$MarkerName, 
													MAF = analysis_set$Freq1, 
													beta = analysis_set$beta.x, 
													N = analysis_set$Weight,
													type = "quant"),
									dataset2 = list(varbeta = analysis_set$varbeta,
														snp = analysis_set$MarkerName, 
														MAF = analysis_set$MAF, 
														beta = analysis_set$beta.y,
														N = summary_list$n[run_no],
														type = "quant")
									)

		} else if (analysis_type == "cc"){

			coloc.res = coloc.abf(dataset1 = list(varbeta = (analysis_set$se ^ 2),
													snp = analysis_set$MarkerName, 
													MAF = analysis_set$Freq1, 
													beta = analysis_set$beta.x, 
													N = analysis_set$Weight,
													type = "quant"),
									dataset2 = list(varbeta = analysis_set$varbeta,
														snp = analysis_set$MarkerName, 
														MAF = analysis_set$MAF, 
														beta = analysis_set$beta.y,
														N = summary_list$n[run_no],
														type = "cc")
									)
			
		}

		### Output

		result[i, 1] = sig_variants$snpid[i]
		result[i, 2] = sig_variants$SNP[i]
		result[i, 3] = nrow(analysis_set)
		result[i, 4] = nrow(temp_gwas)
		result[i, 5] = nrow(temp_query)
		result[i, 6:10] = c(coloc.res$summary)[2:6]
		result[i, 11] = min(analysis_set[,"P-value"])
		
	} else {

		cat(paste0("\nRun ", i, "/", nrow(result), " (query pheno = ", summary_list$name[run_no], ") has no overlapping SNPs",
					"\n===\n\n"))

		next

	}

	### Clearing memory

	rm(list = c("analysis_set", "temp_query", "temp_gwas", "coloc.res"))
	gc(reset=TRUE)

	cat(paste0("\nRun ", i, "/", nrow(result), " complete.",
					"\n...\n\n"))

}

heading("Loop complete.")

cat(paste0("\nWriting ", paste0(outDir, out_prefix, "_", pheno, "_", summary_list$name[run_no], ".tsv"), "\n\n...\n"))
print(head(result))

fwrite(result, paste0(outDir, out_prefix, "_", pheno, "_", summary_list$name[run_no], ".tsv"),
		quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

if (nrow(dupli_gwas_snps) > 0){

	fwrite(dupli_gwas_snps, paste0(outDir, out_prefix, "_", pheno, "_", summary_list$name[run_no], "_duplicated_variants.tsv"),
		quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

}

heading("Write complete. Tests on sea shells indicated they were by the sea shore.")



