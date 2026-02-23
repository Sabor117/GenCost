##### =========================== #####

### Setup environment

### Script is for running CoPheScan
### Input = list of phenotypes + input GWAS
### CoPheScan phenos are munged first
### Also need a table of sample sizes

##### =========================== #####

library(data.table)
library(stringr)
library(cophescan, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(stats)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
mungedphenoDir = paste0(mainDir, "processing/phewas_gwas_sumstats/munged_input/")
gctaDir = paste0(mainDir, "outputs/gcta_cojo_v4/jma_combined/")
metaDir = paste0(mainDir, "outputs/METAL_v4/meta_v4_pruned/")
outDir = paste0(mainDir, "outputs/cophescan/")

### Files

phewaspheno_files = Sys.glob(paste0(mungedphenoDir, "*")) 

### Selecting and reading meta-analysis based on list of metas

analysis_list = Sys.glob(paste0(gctaDir, "*ALL*annotated*"))
analysis_list = analysis_list[!(grepl("_no", analysis_list))]

### N and analysis table

analysis_table_file = paste0(mainDir, "processing/phewas_gwas_sumstats/phewas_summary_table.txt")
analysis_table = fread(analysis_table_file, data.table = FALSE)

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

### Read run from command line

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("51")

run_no = as.numeric(args[1])

phenotype_file = analysis_list[run_no]
phenotype = strsplit(basename(phenotype_file), "_gcta")[[1]][1]

### Read GCTA file for top SNPs
### And pruned METAL file for sumstats - MAY NEED TO READ FULL SUMSTATS?

topsnps = fread(phenotype_file, data.table = FALSE)
sumstats = fread(paste0(metaDir, phenotype, "_metal_output_1_pruned.TBL.gz"), data.table = FALSE)


##### =========================== #####

### Starting analysis

##### =========================== #####

heading("'Tis but a scratch!")

sessionInfo()
start_time = Sys.time()

### For GenCOST ~40 phenotype files with >50 SNPs per input GWAS
### Makes more sense to read PheWAS files once and then do the PheWAS for each SNP
### So: run over phenotypes and THEN SNPs rather than the other way around

outfile = data.frame(matrix(ncol = 19, nrow = 0))
colnames(outfile) = c("rsid",
						"chr",
						"hg19_pos",
						"hg38_snpid",
						"beta",
						"se",
						"p",
						"querysnp",
						"phenotype",
						"querytrait",
						"nsnps",
						"PP.Hn",
						"PP.Ha",
						"PP.Hc",
						"lBF.Ha",
						"lBF.Hc",
						"prior_pn",
						"prior_pa",
						"prior_pc")

for (nfile in 1:length(phewaspheno_files)){

	### First read in analysis GWAS sumstats
	### And get file name

	curr_pheno_file = phewaspheno_files[nfile]
	basefile_name = basename(curr_pheno_file)
	basefile_name = gsub(".txt.gz", "", basefile_name)
	basefile_name = gsub(".tsv.gz", "", basefile_name)

	currPheno = fread(curr_pheno_file, data.table = FALSE)

	### Use file name to read from table:
	### N, analysis name, analysis type

	phenotype_meta_data = analysis_table[analysis_table$file == basefile_name,]
    
	phenotype_n = as.numeric(phenotype_meta_data$n)
	phenotype_name = phenotype_meta_data$name
	phenotype_type = phenotype_meta_data$analysis
	phenotype_build = as.numeric(phenotype_meta_data$build_hg)

	if (phenotype_name == "Melanoma"){

		cat(paste0("\nMelanoma sumstats are scuffed. Skipping.\n\n"))

		print(summary(currPheno))

		next

	}

	if (phenotype_name == "RheumArthritis" | phenotype_name == "Albumin"){

		cat(paste0("\nRheumArthritis and Albumin keeps failing during run. RUN MANUALLY.\n\n"))

		print(summary(currPheno))

		next

	}

	cat(paste0("\nBeginning CoPheScan for ", phenotype, " with current phenotype: ", phenotype_name,
			"\nFile ", nfile, " of ", length(phewaspheno_files), " files.",
			"\nAnalysis type: ", phenotype_type,
			"\nN: ", phenotype_n,
			"\nExample of data:\n\n"))

	print(head(currPheno))
	print(summary(currPheno))

	### Many files contain large number of MAF == NA
	### Remove them ahead of time
	### Some files also having missing chr:pos data

	currPheno = currPheno[!(is.na(currPheno$MAF)),]
	
	cat(paste0("\nMAF values of NA have been removed\n\n"))

	print(head(currPheno))
    print(summary(currPheno))

	currPheno = currPheno[complete.cases(currPheno[ , c("chr", "pos")]), ]
	
	cat(paste0("\nMissing position data has been removed\n\n"))

	print(head(currPheno))
    print(summary(currPheno))

	### Start loop for SNPs in own GWAS

	for (nsnp in 1:nrow(topsnps)){

		### SNPID required in EITHER rsid OR chr_pos (depending on SNPIDs of sumstats)

		currsnp = topsnps$SNP[nsnp]

		currsnp_hg19 = paste0(topsnps$Chr[nsnp], "_", topsnps$bp[nsnp])
		currsnp_hg38 = paste0(strsplit(topsnps$hg38_snpid[nsnp], "_")[[1]][1], "_", strsplit(topsnps$hg38_snpid[nsnp], "_")[[1]][2])

		cat(paste0("\nCoPheScan for ", nsnp, " of ", nrow(topsnps), " (", phenotype, ")",
						"\nSNP: ", currsnp,
						"\nhg19: ", currsnp_hg19,
						"\nhg38: ", currsnp_hg38,
						"\n======\n\n"))

		### SNP positions also required for file pruning
		### Requires both hg19 and hg38 checking

		currChr = topsnps$Chr[nsnp]

		if (phenotype_build == 19){

			currPos = topsnps$bp[nsnp]

		} else if (phenotype_build == 38){

			currPos = as.numeric(strsplit(topsnps$hg38_snpid[nsnp], "_")[[1]][2])

		}

		### Upper and lower limits defined as +/- 500k bp

		pos_upper = currPos + 500000
		pos_lower = currPos - 500000

		### Defining test set of SNPs as the ones which surrounded the current SNP
		### Defined by upper and lower limits

		test_set = currPheno[currPheno$chr == currChr,]
		test_set = test_set[test_set$pos <= pos_upper & test_set$pos >= pos_lower,]

		### Input files have multiple SNP formats: rsids, chr_pos_a1_a2, chr:pos:a1:a2, NA
		### Ideally use rsid - if it is another SNP format then MAKE a new "snp" column of "chr_pos" - no alleles
		### Also define test_snp ID based on this

		if (any(grepl("rs", currPheno$snp))){

			test_snp = currsnp

			test_set$snp[is.na(test_set$snp)] = paste0(test_set$chr[is.na(test_set$snp)], "_",
														test_set$pos[is.na(test_set$snp)])

		} else {

			test_snp = paste0(currChr, "_", currPos)

			test_set$snp = paste0(test_set$chr, "_",
									test_set$pos)

		}

		### Check whether test SNP is actually present in data-set

		if (!(test_snp %in% test_set$snp)){

			cat(paste0("\nSNP ", nsnp, " of ", nrow(topsnps), " (", currsnp, ") not present in data.\n...\n\n"))

			outrow = data.frame(rsid = currsnp,
							chr = currChr,
							hg19_pos = topsnps$bp[nsnp],
							hg38_snpid = topsnps$hg38_snpid[nsnp],
							beta = topsnps$b[nsnp],
							se = topsnps$se[nsnp],
							p = topsnps$p[nsnp],
							querysnp = NA,
							phenotype = phenotype,
							querytrait = phenotype_name,
							nsnps = NA,
							PP.Hn = NA,
							PP.Ha = NA,
							PP.Hc = NA,
							lBF.Ha = NA,
							lBF.Hc = NA,
							prior_pn = NA,
							prior_pa = NA,
							prior_pc = NA
							)

			outfile = rbind(outfile, outrow)

			next

		}

		### A little bit of data QC required
		### Remove duplicated SNPs, MAF == 0/1, varbeta = NA
		### Track number of SNPs removed

		number_rows1 = nrow(test_set)

		cat(paste0("\nCoPheScan test set created. Has ", nrow(test_set), " SNPs in data.\n...\n\n"))
		print(head(test_set))
		print(summary(test_set))

		test_set = test_set[!(duplicated(test_set$snp)),]
		number_rows2 = nrow(test_set)

		test_set = test_set[test_set$MAF != 0 & test_set$MAF != 1,]
		number_rows3 = nrow(test_set)

		test_set = test_set[complete.cases(test_set[ , c("beta", "varbeta")]), ]
		number_rows4 = nrow(test_set)

		test_set = test_set[test_set$beta != 0,]
		number_rows5 = nrow(test_set)

		cat(paste0("\nTest set cleaned.\nRemoved duplicates: ", number_rows1 - number_rows2,
					"\nRemoved unacceptable MAFs (1 or 0): ", number_rows2 - number_rows3,
					"\nRemoved missing varbetas (NAs): ", number_rows3 - number_rows4,
					"\nRemoved betas of exactly 0: ", number_rows4 - number_rows5,
					"\nFinal number of SNPs for CoPheScan: ", number_rows5, "\n...\n\n"))

		
		### Re-check whether this filtering somehow removed your SNP

		if (!(test_snp %in% test_set$snp)){

			cat(paste0("\nSNP ", nsnp, " of ", nrow(topsnps), " (", currsnp, ") was removed by filters.\n...\n\n"))

			outrow = data.frame(rsid = currsnp,
							chr = currChr,
							hg19_pos = topsnps$bp[nsnp],
							hg38_snpid = topsnps$hg38_snpid[nsnp],
							beta = topsnps$b[nsnp],
							se = topsnps$se[nsnp],
							p = topsnps$p[nsnp],
							querysnp = NA,
							phenotype = phenotype,
							querytrait = phenotype_name,
							nsnps = NA,
							PP.Hn = NA,
							PP.Ha = NA,
							PP.Hc = NA,
							lBF.Ha = NA,
							lBF.Hc = NA,
							prior_pn = NA,
							prior_pa = NA,
							prior_pc = NA
							)

			outfile = rbind(outfile, outrow)

			next

		}

		### Convert input set into a list for cophescan (whyyyyyy?)

		test_set = list(snp = test_set$snp,
						beta = test_set$beta,
						varbeta = test_set$varbeta,
						MAF = test_set$MAF,
						N = phenotype_n,
						type = phenotype_type)

		cat(paste0("\nCoPheScan input for ", phenotype_name, ":\n\n"))
		print(str(test_set))

		### Define default fixed priors
		### May need to be adjusted if there are many SNPs

		pa_default = 3.82e-5
		pc_default = 1.82e-3
		test_set_nsnps = length(test_set$snp)

		if (((pa_default*(test_set_nsnps-1)) - pc_default) >= 1) {

			cat(paste0("\nYou're the priors? Well, I didn't vote for you...\n"))

			priors = adjust_priors(nsnps = test_set_nsnps, pa = pa_default, pc = pc_default)
			pa_adjust = priors[["pa"]]
			pc_adjust = priors[["pc"]]

			print(priors)

			cat(paste0("\n...\n"))

			res.single = cophe.single(test_set, querysnpid = test_snp, querytrait = phenotype_name, pa = pa_adjust, pc = pc_adjust)

		} else {

			cat(paste0("\nI guess we take those priors then...\n"))

			res.single = cophe.single(test_set, querysnpid = test_snp, querytrait = phenotype_name, pa = pa_default, pc = pc_default)

		}

		### Run CoPheScan
		### Convert output into a table for merging with full output table
		### Priors input should be default unless corrected

		outrow = data.frame(rsid = currsnp,
							chr = currChr,
							hg19_pos = topsnps$bp[nsnp],
							hg38_snpid = topsnps$hg38_snpid[nsnp],
							beta = topsnps$b[nsnp],
							se = topsnps$se[nsnp],
							p = topsnps$p[nsnp],
							querysnp = res.single$querysnp,
							phenotype = phenotype,
							querytrait = phenotype_name,
							nsnps = res.single$summary$nsnps,
							PP.Hn = res.single$summary$PP.Hn,
							PP.Ha = res.single$summary$PP.Ha,
							PP.Hc = res.single$summary$PP.Hc,
							lBF.Ha = res.single$summary$lBF.Ha,
							lBF.Hc = res.single$summary$lBF.Hc,
							prior_pn = res.single$priors[1],
							prior_pa = res.single$priors[2],
							prior_pc = res.single$priors[3]
							)

		outfile = rbind(outfile, outrow)

		cat(paste0("\nSNP ", nsnp, " of ", nrow(topsnps), " complete.\n...\n\n"))

	}

	heading(paste0("CoPheScan run ", nfile, " of ", length(phewaspheno_files), " complete."))

	cat("\n......\n......\n......\n......\n......\n......\n......\n")

}

fwrite(outfile, paste0(outDir, phenotype, "_cophescan_table_res.txt"), quote = FALSE, row.names = FALSE, na = "NA")

heading(paste0("ALL CoPheScan runs for phenotype", phenotype, " complete. Launching missiles."))



