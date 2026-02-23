##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)

options(scipen = 999) # prevent scientific notation in output

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
outDir = paste0(mainDir, "processing/gcta_intermediate_files_MAF_0_001/")

#SNP_translation_file = paste0(mainDir, "processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out")

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

##### =========================== #####

### Make GCTA input

##### =========================== #####

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

currMeta_file = args[1]

### Script takes base file name from command line
### Read files

heading(paste0("Working on file: ", currMeta_file, " (taken from command line)"))

currMeta = fread(paste0(metaDir, currMeta_file), data.table = FALSE)
#snp_translate = fread(SNP_translation_file, data.table = FALSE)

cat(paste0("\nCurrent nrow of file: ", nrow(currMeta), "\n\n"))
print(head(currMeta))

### Remove low AF and P-vals

pval_col = which(colnames(currMeta) == "P-value")

currMeta = currMeta[currMeta$Freq1 >= 0.001,]
currMeta = currMeta[currMeta$Freq1 <= 0.999,]
currMeta = currMeta[currMeta[,pval_col] <= 0.05,]

### Remove SNPs which are present in fewer than 1/3 cohorts

third_count = floor(str_count(currMeta$Direction)[1] / 3)

if (third_count == 1){

  third_count = 2

}

currMeta = currMeta %>%
	filter((str_count(Direction, "\\-") + str_count(Direction, "\\+")) >= third_count)

### Replace chr23 with "X"

currMeta$MarkerName = gsub("^23:", "X:", currMeta$MarkerName)


##### =========================== #####

### Translation to rsIDs using UKB SNPs is no longer necessary as it happens prior to meta-analysis

##### =========================== #####

### Add marker name column to snp_translate and remove non-rsIDs

#snp_translate = snp_translate[grepl("rs", snp_translate$V4),]
#snp_translate$MarkerName = paste0(gsub("chr", "", snp_translate$V1), "_", snp_translate$V2)
#colnames(snp_translate)[4] = "rsid"

#cat(paste0("\nMerging starts.\n\n"))

### Merge based on MarkerName

#currMeta = merge(currMeta, snp_translate[,c("MarkerName", "rsid")], by = "MarkerName", all.x = TRUE)
#cat(paste0("\nNrow of file post merge: ", nrow(currMeta), "\n\n"))
#print(head(currMeta))

### Remove duplicated SNPs

#num_dups = length(which(duplicated(currMeta$MarkerName)))

#duplicated_snpids = currMeta$MarkerName[which(duplicated(currMeta$MarkerName))]

#cat(paste0("\nNumber of duplicated SNP IDs: ", num_dups, "\n\n"))
#print(currMeta[currMeta$MarkerName %in% duplicated_snpids,])

#currMeta = currMeta[!(currMeta$MarkerName %in% duplicated_snpids),]

### Add additional SNPIDs (in case)

#currMeta$rsid[is.na(currMeta$rsid)] = paste0(currMeta$MarkerName[is.na(currMeta$rsid)], "_", toupper(currMeta$Allele1[is.na(currMeta$rsid)]), "_", toupper(currMeta$Allele2[is.na(currMeta$rsid)]))




##### =========================== #####

### REST OF SCRIPT

##### =========================== #####

currMeta_snp_frame = str_split_fixed(currMeta$MarkerName, ":", 4)

currMeta$rsid = ifelse(currMeta_snp_frame[,2] == "",
						currMeta_snp_frame[,1],
						paste0(currMeta_snp_frame[,1], ":", currMeta_snp_frame[,2], "_", currMeta_snp_frame[,3], "_", currMeta_snp_frame[,4])
						)

### Get betas and SEs

currMeta$beta1 = currMeta$Zscore / sqrt((2 * currMeta$Freq1) * (1 - currMeta$Freq1) * (currMeta$Weight + (currMeta$Zscore^2)))

currMeta$se = 1 / sqrt((2 * currMeta$Freq1) * (1 - currMeta$Freq1) * (currMeta$Weight + (currMeta$Zscore^2)))

cat(paste0("\nCalculations complete, creating output.\n\n"))

### Simplify and write

output_file = data.frame(SNP = currMeta$rsid,
                            A1 = toupper(currMeta$Allele1),
                            A2 = toupper(currMeta$Allele2),
                            freq = currMeta$Freq1,
                            b = currMeta$beta1,
                            se = currMeta$se,
                            p = currMeta[,10],
                            N = currMeta$Weight
                            )

print(head(output_file, n = 10))

fwrite(output_file, paste0(outDir, gsub("1.TBL.gz", "gcta_input.txt", currMeta_file)), row.names = FALSE, sep = "\t", quote = FALSE)


