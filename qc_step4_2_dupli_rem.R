##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)
library(R.utils)
options(scipen = 999)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")
tempDir = paste0(mainDir, "tmpdir/")

all_files = Sys.glob(paste0(sumstatDir, "/*/*_alleles_aligned.txt.gz"))

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("1")

run_no = as.numeric(args[1])

curr_file_name = all_files[run_no]


##### =========================== #####

### Run script

##### =========================== #####

cat(paste0("\nStarting run.\n========\n\n"))

curr_file = fread(curr_file_name, data.table = FALSE, tmpdir = tempDir)

cat(paste0("\nNow running on file: ", curr_file_name, "\n\nRemoving duplicates.\n========================\n\n"))

pre_dup_nrow = nrow(curr_file)

curr_dupli_snps = unique(curr_file$full_snpid[which(duplicated(curr_file$full_snpid))])

rem_file = curr_file[curr_file$full_snpid %in% curr_dupli_snps,]
curr_file = curr_file[!(curr_file$full_snpid %in% curr_dupli_snps),]

rem_file_fixed = rem_file %>%
                    group_by(full_snpid) %>%
                    mutate(abs_diff = abs(af_alt - freq)) %>%
                    filter(abs_diff == min(abs_diff)) %>%
                    ungroup()

rem_file_fixed = rem_file_fixed[rem_file_fixed$FLIP == 1,]

rem_file_fixed = as.data.frame(rem_file_fixed)

curr_file = rbind(curr_file, rem_file_fixed[,(colnames(curr_file))])

post_dup_nrow = nrow(curr_file)

fwrite(curr_file, curr_file_name, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip")

cat(paste0("\nFile written. ", pre_dup_nrow - post_dup_nrow, " lines removed.\n========================\n\n"))
