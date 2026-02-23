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
outDir = paste0(mainDir, "processing/meta_sumstats/")
tempDir = paste0(mainDir, "tmpdir/")

sumstat_list = Sys.glob(paste0(sumstatDir, "*/*.gz"))
sumstat_list = sumstat_list[!(grepl("gnh_snp_infos", sumstat_list))]
sumstat_list = sumstat_list[!(grepl("aligned", sumstat_list))]
sumstat_list = sumstat_list[!(grepl("NTR", sumstat_list))]
sumstat_list = sumstat_list[!(grepl("AOU", sumstat_list))]

ukb_hg38_file = paste0(mainDir, "processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out")

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Had to be me. Someone else would have gotten it wrong.")

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("51")

run_no = as.numeric(args[1])

run_start = ((run_no - 1) * 5) + 1
run_end = run_no * 5


##### =========================== #####

### Run script

##### =========================== #####

### All phenotype combinations

phenotypes = c("ALL", "PRIM", "DRUG", "IN", "INOUT")

stratifs = c("ALL", "M", "F", "5_18", "19_35", "36_55", "56_75", "76_95")


### Create combined versions for each combination

combined_phenos = sapply(phenotypes, function(p) {
                            sapply(stratifs, function(s) {
                                paste(p, s, sep = "_")
                            })
                        }, simplify = "matrix")

### Flatten the matrix to get a single vector

combined_phenos = unique(c(combined_phenos))

### Make directories

for (i in 1:length(combined_phenos)){

    cat(paste0("\n\nCommand is:\nmkdir -p ", outDir, "/", combined_phenos[i], "/"))

    system(paste0("mkdir -p ", outDir, "/", combined_phenos[i], "/"))

}

ukb_hg38 = fread(ukb_hg38_file, data.table = FALSE, tmpdir = tempDir)
colnames(ukb_hg38) = c("chr", "pos", "pos1", "rsid", "hg19_snpid")

ukb_alleles = str_split_fixed(ukb_hg38$hg19_snpid, ":", 4)

ukb_snp_comparison = data.frame(rsid = ukb_hg38$rsid,
                                hg19_snpid1 = paste(ukb_alleles[,1], ukb_alleles[,2], ukb_alleles[,3], ukb_alleles[,4], sep = ":"),
                                hg19_snpid2 = paste(ukb_alleles[,1], ukb_alleles[,2], ukb_alleles[,4], ukb_alleles[,3], sep = ":"),
                                hg38_snpid1 = paste(gsub("chr", "", ukb_hg38$chr), ukb_hg38$pos, ukb_alleles[,3], ukb_alleles[,4], sep = ":"),
                                hg38_snpid2 = paste(gsub("chr", "", ukb_hg38$chr), ukb_hg38$pos, ukb_alleles[,4], ukb_alleles[,3], sep = ":")
                                )

ukb_snp_comparison$hg19_snpid1 = gsub("^23:", "X:", ukb_snp_comparison$hg19_snpid1)
ukb_snp_comparison$hg19_snpid2 = sub("^23:", "X:", ukb_snp_comparison$hg19_snpid2)
ukb_snp_comparison$hg38_snpid1 = sub("^23:", "X:", ukb_snp_comparison$hg38_snpid1)
ukb_snp_comparison$hg38_snpid2 = sub("^23:", "X:", ukb_snp_comparison$hg38_snpid2)

hg38_list = c("CHB", "FINNGEN", "MGBB", "QGP", "AOU", "EstBB", "GE", "GNH", "NTR")

### Begin munging

for (i in run_start:run_end){

    ### If the run-number exceeds the number of files, we end the script.

    if (i > length(sumstat_list)){

        cat(paste0("\nFinal run completed.\n"))

        end_time = Sys.time()

        total_time = end_time - start_time

        cat(paste0("\nSession complete! Run-time was: ", total_time, ".\n\n=========================\n\n"))

        break

    }

    ### Read files

    curr_sumstats = fread(sumstat_list[i], data.table = FALSE, tmpdir = tempDir)

    ### Split filename to extract phenotype + stratifications

    filename_split = gsub("FM", "ALL", sumstat_list[i])
    filename_split = gsub("5_95", "ALL", filename_split)
    filename_split = gsub("All", "ALL", filename_split)
    filename_split = gsub("19_95", "ALL", filename_split)
    filename_split = gsub("15_102", "ALL", filename_split)
    filename_split = gsub("18_35", "19_35", filename_split)
    filename_split = gsub("76_99", "76_95", filename_split)
    filename_split = gsub("FEMALES", "F", filename_split)
    filename_split = gsub("MALES", "M", filename_split)
    filename_split = gsub("FEMALE", "F", filename_split)
    filename_split = gsub("MALE", "M", filename_split)
    filename_split = gsub("DRUGS", "DRUG", filename_split)
    filename_split = gsub("PRIMARY", "PRIM", filename_split)
    filename_split = gsub("INPATIENTS", "IN", filename_split)
    filename_split = gsub("75_96", "76_95", filename_split)
    filename_split = gsub("TOTAL", "ALL", filename_split)
    filename_split = gsub("INPATIENT", "IN", filename_split)
    filename_split = gsub("75plus", "76_95", filename_split)

    filename = gsub(".txt.gz", "", basename(filename_split))
    filename = gsub(".tsv.gz", "", basename(filename))

    filename_split = strsplit(basename(filename_split), "\\.")

    cohort = filename_split[[1]][1]
    pheno = filename_split[[1]][3]
    stratif_1 = filename_split[[1]][4]
    stratif_2 = filename_split[[1]][5]
    population = filename_split[[1]][6]
    ssize = filename_split[[1]][7]

    stratif_name = case_when(stratif_1 == "ALL" & stratif_2 == "ALL" ~ "ALL",
                                stratif_1 != "ALL" & stratif_2 == "ALL" ~ stratif_1,
                                stratif_1 == "ALL" & stratif_2 != "ALL" ~ stratif_2,
                                .default = as.character(stratif_1))

    ### Check splitting has worked correctly

    cat(paste0("\nCurrent run number (", i, ") will be on: ", cohort,
                "\nPhenotype: ", pheno,
                "\nStratified by: ", stratif_name,
                "\nOutput goes to: ", outDir, pheno, "_", stratif_name,
                "\nN = ", ssize,
                "\nFile = ", basename(sumstat_list[i]),
                "\n\n"))

    if (cohort == "QGP"){

        realColumns = colnames(curr_sumstats)[-(which(colnames(curr_sumstats) == "SNPID"))]

        curr_sumstats[,ncol(curr_sumstats)] = NULL

        colnames(curr_sumstats) = realColumns

    }

    print(head(curr_sumstats))
    print(summary(curr_sumstats))
    print(nrow(curr_sumstats))
    cat(paste0("\n=====\n"))

    ### Number of rows at start of scripts

    init_nrow = nrow(curr_sumstats)

    ### Define where to output and column names

    curr_outDir = paste0(outDir, pheno, "_", stratif_name, "/")

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

    beta_col = case_when(cohort == "UKB" ~ "BETA",
                            cohort == "CHB" ~ "BETA",
                            cohort == "QGP" ~ "BETA",
                            cohort == "EstBB" ~ "BETA",
                            cohort == "GNH" ~ "BETA",
                            .default = "beta")

    pval_col = case_when(cohort == "UKB" ~ "LOG10P",
                            cohort == "CHB" ~ "P",
                            cohort == "QGP" ~ "p.value",
                            cohort == "EstBB" ~ "p.value",
                            cohort == "GNH" ~ "LOG10P",
                            cohort == "GE" ~ "Pvalue",
                            .default = "Pval")

    se_col = case_when(cohort == "UKB" ~ "SE",
                            cohort == "CHB" ~ "SE",
                            cohort == "QGP" ~ "SE",
                            cohort == "EstBB" ~ "SE",
                            cohort == "GNH" ~ "SE",
                            .default = "se")

    info_col = case_when(cohort == "UKB" ~ "INFO",
                            cohort == "QGP" ~ "imputationInfo",
                            .default = "INFO")

    if ("N" %in% colnames(curr_sumstats)){ 

        col_list = c(chrom_col, pos_col, a1_col, a2_col, freq_col, beta_col, se_col, pval_col, "N")

    } else {

        col_list = c(chrom_col, pos_col, a1_col, a2_col, freq_col, beta_col, se_col, pval_col)

    }

    if (cohort == "EstBB"){

        colnames(curr_sumstats)[which(colnames(curr_sumstats) == "SNPID")] = "SNPID1"

        col_list = c(col_list, "SNPID1")


    }

    if (cohort == "QGP"){

        colnames(curr_sumstats)[which(colnames(curr_sumstats) == "rsid")] = "SNPID1"

        col_list = c(col_list, "SNPID1")

    }

    if (cohort == "UKB"){

        colnames(curr_sumstats)[which(colnames(curr_sumstats) == "ID")] = "SNPID1"

        col_list = c(col_list, "SNPID1")

    }

    if (cohort == "GNH"){

        colnames(curr_sumstats)[which(colnames(curr_sumstats) == "ID")] = "SNPID1"

        col_list = c(col_list, "SNPID1")

    }

    curr_sumstats[,pos_col] = as.numeric(curr_sumstats[,pos_col])
    curr_sumstats[,chrom_col] = as.numeric(curr_sumstats[,chrom_col])    
    curr_sumstats[,beta_col] = as.numeric(curr_sumstats[,beta_col])
    curr_sumstats[,freq_col] = as.numeric(curr_sumstats[,freq_col])
    curr_sumstats[,se_col] = as.numeric(curr_sumstats[,se_col])
    curr_sumstats[,pval_col] = as.numeric(curr_sumstats[,pval_col])

    if (cohort != "NTR" & cohort != "AOU"){

        curr_sumstats[,info_col] = as.numeric(curr_sumstats[,info_col])

    }

    curr_sumstats = curr_sumstats[complete.cases(curr_sumstats[,c(col_list)]),]

    ### Check current state of file

    complete_nrow = nrow(curr_sumstats)

    cat(paste0("\nColumns selected. Nrow = ", complete_nrow, " (", init_nrow - complete_nrow, " rows removed for missing data):\n\n"))
    print(head(curr_sumstats[,c(col_list)]))
    print(summary(curr_sumstats))
    cat(paste0("\n=====\n"))

    complete_del = init_nrow - complete_nrow

    ### Remove low INFO SNPs
    ### NTR does not have INFO

    if (cohort != "NTR" & cohort != "AOU"){

        curr_sumstats = curr_sumstats[!(is.na(curr_sumstats[,info_col])),]

        curr_sumstats = curr_sumstats[curr_sumstats[,info_col] >= 0.6,]

        cat(paste0("\nData refined. Nrow = ", nrow(curr_sumstats), " (", complete_nrow - nrow(curr_sumstats), " rows removed for having low or missing INFO):\n\n"))
        print(head(curr_sumstats[,c(col_list)]))
        print(summary(curr_sumstats))
        cat(paste0("\n=====\n"))

        info_del = complete_nrow - nrow(curr_sumstats)

    } else {

        info_del = 0

    }

    ### Create full SNPID to megre with aligned alleles

    curr_sumstats[,chrom_col] = gsub(23, "X", curr_sumstats[,chrom_col])
    curr_sumstats[,chrom_col] = gsub("23", "X", curr_sumstats[,chrom_col])

    curr_sumstats$full_snpid = paste(curr_sumstats[,chrom_col],
                                curr_sumstats[,pos_col],
                                curr_sumstats[,a1_col],
                                curr_sumstats[,a2_col],
                                sep = ":")

    ### Some files may contain duplicated SNPIDs
    ### Same SNP (position/location etc) but non-matching betas/SEs/etc
    ### These must simply be removed

    duplicated_snps = curr_sumstats$full_snpid[which(duplicated(curr_sumstats$full_snpid))]

    dupli_snps_rem = nrow(curr_sumstats[curr_sumstats$full_snpid %in% duplicated_snps,])

    the_naughty_list_of_snps = curr_sumstats[curr_sumstats$full_snpid %in% duplicated_snps,]

    curr_sumstats = curr_sumstats[!(curr_sumstats$full_snpid %in% duplicated_snps),]

    cat(paste0("\nDuplicated SNPs removed. Nrow = ", nrow(curr_sumstats), ", ", dupli_snps_rem, " SNPs removed.\n\n"))
    print(head(curr_sumstats))
    cat(paste0("\n=====\n"))

    ### Read the alleles aligned file - contains correct positions and FLIP statuses based on gnoMAD

    if (cohort == "EstBB") {

        allele_flips = fread(paste0(dirname(sumstat_list[i]), "/", toupper(cohort), "_alleles_aligned.txt.gz"), data.table = FALSE, tmpdir = tempDir)

    } else {

        allele_flips = fread(paste0(dirname(sumstat_list[i]), "/", cohort, "_alleles_aligned.txt.gz"), data.table = FALSE, tmpdir = tempDir)

    }

    cat(paste0("\nAllele flips read. Nrow = ", nrow(allele_flips), ":\n\n"))
    print(head(allele_flips))
    cat(paste0("\n=====\n"))

    ### Remove instances of non-matching alleles

    zero_flips = nrow(allele_flips[allele_flips$FLIP == 0,])
    allele_flips = allele_flips[allele_flips$FLIP != 0,]

    cat(paste0("\n'0' flips removed. Nrow = ", nrow(allele_flips), ", ", zero_flips, " SNPs removed.\n\n"))
    print(head(allele_flips))
    cat(paste0("\n=====\n"))

    if (any(duplicated(allele_flips$full_snpid))){

        pre_dup_nrow = nrow(allele_flips)

        ### Finding duplicated full_snpids in allele_flips
        ### These have been matched as having the same allele in both directions in gnoMAD
        ### Separate allele_flips into duplicated and unique IDs

        allele_flips_dupli_snps = allele_flips$full_snpid[which(duplicated(allele_flips$full_snpid))]

        allele_flips_duplicated = allele_flips[allele_flips$full_snpid %in% allele_flips_dupli_snps,]

        allele_flips = allele_flips[!(allele_flips$full_snpid %in% allele_flips_dupli_snps),]

        ### Now parsing the duplicated ones
        ### First get rid of completely duplicated lines

        allele_flips_duplicated = unique(allele_flips_duplicated)

        ### Merge with the sumstats frequency

        allele_flips_duplicated = merge(allele_flips_duplicated, curr_sumstats[,c("full_snpid", freq_col)], by = "full_snpid")

        ### Do a second round of unique() as the merge could have re-introduced duplicates

        allele_flips_duplicated = unique(allele_flips_duplicated)

        colnames(allele_flips_duplicated)[ncol(allele_flips_duplicated)] = "file_freq"

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

        allele_flips_dupli_snps_2 = output_duplicated_snps_frame$full_snpid[which(duplicated(output_duplicated_snps_frame$full_snpid))]

        output_flips_duplicated = output_duplicated_snps_frame[output_duplicated_snps_frame$full_snpid %in% allele_flips_dupli_snps_2,]
        output_flips_duplicated = output_flips_duplicated[output_flips_duplicated$FLIP == 1,]

        output_duplicated_snps_frame = output_duplicated_snps_frame[!(output_duplicated_snps_frame$full_snpid %in% allele_flips_dupli_snps_2),]

        output_duplicated_snps_frame = rbind(output_duplicated_snps_frame, output_flips_duplicated)

        ### This new output table is used to select the SNPs from allele_flips_duplicated

        output_duplicated_snps_frame$select = paste0(output_duplicated_snps_frame$full_snpid, ":", output_duplicated_snps_frame$FLIP)
        allele_flips_duplicated$select = paste0(allele_flips_duplicated$full_snpid, ":", allele_flips_duplicated$FLIP)

        allele_flips_duplicated = allele_flips_duplicated[allele_flips_duplicated$select %in% output_duplicated_snps_frame$select,]

        ### Bind all of this back with allele_flips
        ### So that flip status can be merged with the sumstats

        allele_flips = rbind(allele_flips, allele_flips_duplicated[,colnames(allele_flips)])

        post_dup_nrow = nrow(allele_flips)

        ### Doing checks of duplications

        cat(paste0("\nDuplicated full_snpids removed. Nrow = ", post_dup_nrow, ", ", pre_dup_nrow - post_dup_nrow, " SNPs removed.\n\n"))
        print(head(allele_flips))
        print(summary(allele_flips))
        cat(paste0("\n=====\n"))
        cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
                    any(duplicated(allele_flips$full_snpid)),
                    "\n Are there any duplicated snpids? ",
                    any(duplicated(allele_flips$snpid)), "\n\n"))
        cat(paste0("\n=====\n"))

    } else {

        pre_dup_nrow = nrow(allele_flips)
        post_dup_nrow = nrow(allele_flips)

        cat(paste0("\nDuplicated full_snpids removed. Nrow = ", nrow(allele_flips), ", No SNPs removed.\n\n"))
        print(head(allele_flips))
        print(summary(allele_flips))
        cat(paste0("\n=====\n"))
        cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
                    any(duplicated(allele_flips$full_snpid)),
                    "\n Are there any duplicated snpids? ",
                    any(duplicated(allele_flips$snpid)), "\n\n"))
        cat(paste0("\n=====\n"))

    }

    full_snpid_del = pre_dup_nrow - post_dup_nrow

    ### Need to do additional round of duplication removal as some basic SNPIDs may still be present
    ### Likely multi-allelic SNPs (I.e. the full_snpid cannot catch them)
    ### These are also just removed

    if (any(duplicated(allele_flips$snpid))){

        pre_dup_nrow = nrow(allele_flips)

        duplicated_snps = allele_flips$snpid[which(duplicated(allele_flips$snpid))]

        print(head(allele_flips[allele_flips$snpid %in% duplicated_snps[1:5],], n = 10))

        allele_flips = allele_flips[!(allele_flips$snpid %in% duplicated_snps),]

        post_dup_nrow = nrow(allele_flips)

        cat(paste0("\nDuplicated basic snpids removed. Nrow = ", post_dup_nrow, ", ", pre_dup_nrow - post_dup_nrow, " SNPs removed.\n\n"))
        print(head(allele_flips))
        print(summary(allele_flips))
        cat(paste0("\n=====\n"))
        cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
                    any(duplicated(allele_flips$full_snpid)),
                    "\n Are there any duplicated snpids? ",
                    any(duplicated(allele_flips$snpid)), "\n\n"))
        cat(paste0("\n=====\n"))

    } else {

        pre_dup_nrow = nrow(allele_flips)
        post_dup_nrow = nrow(allele_flips)

        cat(paste0("\nDuplicated basic snpids removed. Nrow = ", nrow(allele_flips), ", No SNPs removed.\n\n"))
        print(head(allele_flips))
        print(summary(allele_flips))
        cat(paste0("\n=====\n"))
        cat(paste0("\nAre there currently any duplicated full_snpids in allele_flips? ",
                    any(duplicated(allele_flips$full_snpid)),
                    "\n Are there any duplicated snpids? ",
                    any(duplicated(allele_flips$snpid)), "\n\n"))
        cat(paste0("\n=====\n"))

    }

    half_snpid_del = pre_dup_nrow - post_dup_nrow

    pre_merge_nrow = nrow(curr_sumstats)

    curr_sumstats = inner_join(curr_sumstats[,c("full_snpid", col_list)], allele_flips, by = "full_snpid")

    cat(paste0("\nSumstats merged with allele flips. Nrow = ", nrow(curr_sumstats), ".\n\n"))
    print(head(curr_sumstats))
    cat(paste0("\n=====\n"))

    snpid_merge_del = pre_merge_nrow - nrow(curr_sumstats)

    ### Now do checks for duplicated SNPs once more
    ### And removed once again in case

    cat(paste0("\nAre there currently any duplicated full_snpids in the merged sumstats? ",
                any(duplicated(curr_sumstats$full_snpid)),
                "\n Are there any duplicated basic snpids? ",
                any(duplicated(curr_sumstats$snpid)), "\n\n"))

    pre_dup_nrow = nrow(curr_sumstats)

    duplicated_snps = curr_sumstats$snpid[which(duplicated(curr_sumstats$snpid))]

    print(head(curr_sumstats[curr_sumstats$snpid %in% duplicated_snps[1:5],], n = 10))

    curr_sumstats = curr_sumstats[!(curr_sumstats$snpid %in% duplicated_snps),]

    post_dup_nrow = nrow(curr_sumstats)

    cat(paste0("\nAny merged duplicates are removed. Nrow = ", post_dup_nrow, ", ", pre_dup_nrow - post_dup_nrow, " SNPs removed.\n\n"))
    print(summary(curr_sumstats))
    print(head(curr_sumstats))

    dupli_snps_rem_2 = pre_dup_nrow - post_dup_nrow

    if (cohort != "EstBB" & cohort != "UKB" & cohort != "QGP" & cohort != "GNH"){

        ### Adding in UKB rsIDs

        if (cohort %in% hg38_list){

            if (any(curr_sumstats$full_snpid %in% ukb_snp_comparison$hg38_snpid1)){

                id_comparison_frame_1 = ukb_snp_comparison[ukb_snp_comparison$hg38_snpid1 %in% curr_sumstats$full_snpid,][,c("rsid", "hg38_snpid1")]

            } else {

                id_comparison_frame_1 = data.frame(matrix(ncol = 2, nrow = 0))

            }

            if (any(curr_sumstats$full_snpid %in% ukb_snp_comparison$hg38_snpid2)){

                id_comparison_frame_2 = ukb_snp_comparison[ukb_snp_comparison$hg38_snpid2 %in% curr_sumstats$full_snpid,][,c("rsid", "hg38_snpid2")]

            } else {

                id_comparison_frame_2 = data.frame(matrix(ncol = 2, nrow = 0))

            }

            colnames(id_comparison_frame_1) = c("rsid", "full_snpid")
            colnames(id_comparison_frame_2) = c("rsid", "full_snpid")

            id_comparison_frame = rbind(id_comparison_frame_1, id_comparison_frame_2)

        } else {

            if (any(curr_sumstats$full_snpid %in% ukb_snp_comparison$hg19_snpid1)){

                id_comparison_frame_1 = ukb_snp_comparison[ukb_snp_comparison$hg19_snpid1 %in% curr_sumstats$full_snpid,][,c("rsid", "hg19_snpid1")]

            } else {

                id_comparison_frame_1 = data.frame(matrix(ncol = 2, nrow = 0))

            }

            if (any(curr_sumstats$full_snpid %in% ukb_snp_comparison$hg19_snpid2)){

                id_comparison_frame_2 = ukb_snp_comparison[ukb_snp_comparison$hg19_snpid2 %in% curr_sumstats$full_snpid,][,c("rsid", "hg19_snpid2")]

            } else {

                id_comparison_frame_2 = data.frame(matrix(ncol = 2, nrow = 0))

            }

            colnames(id_comparison_frame_1) = c("rsid", "full_snpid")
            colnames(id_comparison_frame_2) = c("rsid", "full_snpid")

            id_comparison_frame = rbind(id_comparison_frame_1, id_comparison_frame_2)

        }

        pre_rsid_n = nrow(curr_sumstats)

        curr_sumstats = merge(curr_sumstats, id_comparison_frame, by = "full_snpid", all.x = TRUE)

        cat(paste0("\nrsIDs added into sumstats. Nrow = ", nrow(curr_sumstats), ", ", pre_rsid_n - nrow(curr_sumstats), " SNPs removed.\n\n"))
        print(summary(curr_sumstats))
        print(head(curr_sumstats))

        output_sumstats = data.frame(snpid = ifelse(is.na(curr_sumstats$rsid),
                                paste(curr_sumstats[,chrom_col], curr_sumstats$pos, sep = ":"),
                                curr_sumstats$rsid)
                            )

    } else {

        pre_rsid_n = nrow(curr_sumstats)

        output_sumstats = data.frame(snpid = curr_sumstats$SNPID1)

    }

    rsid_merge_del = pre_rsid_n - nrow(curr_sumstats)

    cat(paste0("\nDeletions at stages of script for file: ", basename(sumstat_list[i]),
                "\n========\n\nInitial number of SNPs: ", init_nrow,
                "\nIncomplete SNPs removed (NAs in data): ", complete_del,
                "\nSNPs removed for low INFO (< 0.6): ", info_del,
                "\nSNPs removed with duplicated full chrosomoe-position-alleles: ", dupli_snps_rem,
                "\nSNPs removed following merge with allele flips (allele flips removed duplicated full/half SNPIDs): ", snpid_merge_del,
                "\nSNPs removed for duplicated full SNPIDs following merge: ", dupli_snps_rem_2,
                "\nSNPs removed following merge of UKB rsIDs: ", rsid_merge_del),
                "\nShe turned me into a newt! I got better...\n\n======\n\n")

    output_sumstats$chr = curr_sumstats[,chrom_col]
    output_sumstats$pos = curr_sumstats$pos
    output_sumstats$a1 = ifelse(curr_sumstats$FLIP == 1, curr_sumstats[,a1_col], curr_sumstats[,a2_col])
    output_sumstats$a0 = ifelse(curr_sumstats$FLIP == 1, curr_sumstats[,a2_col], curr_sumstats[,a1_col])
    output_sumstats$beta1 = curr_sumstats[,beta_col] * curr_sumstats$FLIP
    output_sumstats$freq1 = ifelse(curr_sumstats$FLIP == 1, curr_sumstats[,freq_col], (1 - curr_sumstats[,freq_col]))
    output_sumstats$se = curr_sumstats[,se_col]

    missing_rsid_index = which(!(grepl("rs", output_sumstats$snpid)))

    output_sumstats$snpid[missing_rsid_index] = paste(output_sumstats$chr[missing_rsid_index],
                                                        output_sumstats$pos[missing_rsid_index],
                                                        output_sumstats$a1[missing_rsid_index],
                                                        output_sumstats$a0[missing_rsid_index],
                                                        sep = ":")

    if (cohort == "UKB" | cohort == "GNH"){

        curr_sumstats$p_val = 10 ^ -(curr_sumstats[,pval_col])

        output_sumstats$p = curr_sumstats$p_val

    } else {

        output_sumstats$p = curr_sumstats[,pval_col]

    }

    output_sumstats$n = ifelse("N" %in% colnames(curr_sumstats), as.numeric(curr_sumstats$N), ssize)

    output_sumstats = unique(output_sumstats)
                                        
    final_nrow = nrow(output_sumstats)

    cat(paste0("\nOutput sumstats created. Nrow = ", nrow(output_sumstats), ".\n\n"))
    cat(paste0("\nDifference between start and end sumstats. ", init_nrow - final_nrow, " SNPs removed.\n\n"))
    print(summary(output_sumstats))
    print(head(output_sumstats))

    fwrite(output_sumstats, paste0(curr_outDir, filename, "_metal_input_sumstats.txt.gz"),
			row.names = FALSE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

    curr_run_time = Sys.time()

    curr_time = curr_run_time - start_time

    cat(paste0("\nOutput written and zipped. Run-time is: ", curr_time, ".\n\n=========================\n\n"))

}

end_time = Sys.time()

total_time = end_time - start_time

cat(paste0("\nSession complete! Run-time was: ", total_time, ".\n\n=========================\n\n"))
heading("You're the King? Well, I didn't vote for you...")



