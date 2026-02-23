##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(stringr)

options(scipen = 999) # prevent scientific notation in output

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")

all_cohorts = Sys.glob(paste0(sumstatDir, "*/"))
all_files = Sys.glob(paste0(sumstatDir, "*/*.gz"))
all_files = all_files[!(grepl("aligned", all_files))]


##### =========================== #####

### Making LiftOver Input

##### =========================== #####

for (i in 1:length(all_cohorts)){
#for (i in 4:4){

    ### Get cohort name (in caps)

    curr_cohort = toupper(basename(all_cohorts[i]))

    cat(paste0("\nStarting run on ", curr_cohort, ". CHECK 1.", i, " COMPLETE.\n\n"))

    ### Cohort names are pre-defined

    chrom_col = case_when(curr_cohort == "UKB" ~ "CHROM",
                            curr_cohort == "GNH" ~ "CHROM",
                            .default = "CHR")
                            
    pos_col = case_when(curr_cohort == "UKB" ~ "GENPOS",
                            curr_cohort == "QGP" ~ "POS",
                            curr_cohort == "ESTBB" ~ "POS",
                            curr_cohort == "GNH" ~ "GENPOS",
                            .default = "BP")

    a1_col = case_when(curr_cohort == "UKB" ~ "ALLELE1",
                        curr_cohort == "CHB" ~ "ALLELE1",
                        curr_cohort == "QGP" ~ "Allele2",
                        curr_cohort == "ESTBB" ~ "ALLELE1",
                        curr_cohort == "GNH" ~ "ALLELE1",
                        .default = "A1")

    a2_col = case_when(curr_cohort == "UKB" ~ "ALLELE0",
                        curr_cohort == "CHB" ~ "ALLELE0",
                        curr_cohort == "ESTBB" ~ "ALLELE0",
                        curr_cohort == "QGP" ~ "Allele1",
                        curr_cohort == "GNH" ~ "ALLELE0",
                        .default = "A2")

    freq_col = case_when(curr_cohort == "UKB" ~ "A1FREQ",
                            curr_cohort == "CHB" ~ "A1FREQ",
                            curr_cohort == "QGP" ~ "AF_Allele2",
                            curr_cohort == "ESTBB" ~ "A1FREQ",
                            curr_cohort == "GNH" ~ "A1FREQ",
                            .default = "EAF")

    ### Need specific rules for selecting the correct files
    ### EstBB == ESTBB
    ### Only use EUR cohort for the LiftOver file

    if (curr_cohort == "ESTBB"){

        curr_file_set = all_files[grepl("EstBB", all_files)]

    } else if (curr_cohort == "UKB"){

        curr_file_set = all_files[grepl(curr_cohort, all_files)]
        curr_file_set = curr_file_set[grepl("EUR", curr_file_set)]

    } else if (curr_cohort == "GE"){

        curr_file_set = all_files[grepl("GE.KZ", all_files)]

    } else {

        curr_file_set = all_files[grepl(curr_cohort, all_files)]

    }

    ### Get max sample size from file name

    ssize_set = as.numeric(str_split_fixed(basename(curr_file_set), "\\.", 11)[,7])
    ssize_indices = order(ssize_set, decreasing = TRUE)

    curr_file_set = curr_file_set[ssize_indices]

    ### Start LiftOver input .bed

    cohort_output = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(cohort_output) = c("chr", "pos", "pos1", "snpid", "freq1")

    if (curr_cohort == "ESTBB" | curr_cohort == "UKB"){

        cohort_output = data.frame(matrix(ncol = 6, nrow = 0))
        colnames(cohort_output) = c("chr", "pos", "pos1", "snpid", "freq1", "rsid")

    }

    ### Starting main loop
    ### Merges all files selected into on .bed file for input into LiftOver

    for (j in 1:length(curr_file_set)){

        ### Reads current file

        cat(paste0("\nStarting run on ", basename(curr_file_set[j]), ". CHECK ", i, ".", j, ".1 COMPLETE.\n\n"))

        curr_file = fread(curr_file_set[j], data.table = FALSE, tmpdir = paste0(mainDir, "tmpdir/"))

        ### QGP has an extra column which needs to be dealt with

        if (curr_cohort == "QGP"){

            realColumns = colnames(curr_file)[-(which(colnames(curr_file) == "SNPID"))]

            curr_file[,ncol(curr_file)] = NULL

            colnames(curr_file) = realColumns

            print(head(curr_file))

        }

        ### Example of the first file in the set

        if (j == 1){

            cat(paste0("\n=========\nHere is the example header of the first file:\n\n\n"))

            print(head(curr_file))

        }

        cat(paste0("\nFile read. CHECK ", i, ".", j, ".2 COMPLETE.\n\n"))

        ### For each file fix chromosome handling here

        curr_file[,chrom_col] = gsub("23", "X", curr_file[,chrom_col])
        curr_file[,chrom_col] = gsub(23, "X", curr_file[,chrom_col])

        ### For each file it combines:
        ### chrom (chrXXX), pos, pos+1, SNPID, allele frequency (for gnoMAD) and rsID (if available)

        if (curr_cohort == "ESTBB"){

            curr_snp_set = data.frame(chr = paste0("chr", curr_file[,chrom_col]),
                        pos = as.numeric(curr_file[,pos_col]),
                        pos1 = as.numeric(curr_file[,pos_col]) + 1,
                        snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col]),
                        freq1 = curr_file[,freq_col],
                        rsid = curr_file$SNPID)

        } else if (curr_cohort == "UKB"){

            curr_snp_set = data.frame(chr = paste0("chr", curr_file[,chrom_col]),
                        pos = as.numeric(curr_file[,pos_col]),
                        pos1 = as.numeric(curr_file[,pos_col]) + 1,
                        snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col]),
                        freq1 = curr_file[,freq_col],
                        rsid = curr_file$ID)

        } else if (curr_cohort == "QGP"){

            curr_snp_set = data.frame(chr = paste0("chr", curr_file[,chrom_col]),
                        pos = as.numeric(curr_file[,pos_col]),
                        pos1 = as.numeric(curr_file[,pos_col]) + 1,
                        snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col]),
                        freq1 = curr_file[,freq_col],
                        rsid = curr_file$rsid)

        } else if (curr_cohort == "GE"){

            curr_snp_set = data.frame(chr = paste0("chr", curr_file[,chrom_col]),
                        pos = as.numeric(curr_file[,pos_col]),
                        pos1 = as.numeric(curr_file[,pos_col]) + 1,
                        snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col]),
                        freq1 = curr_file[,freq_col])

        } else {

            curr_snp_set = data.frame(chr = paste0("chr", curr_file[,chrom_col]),
                            pos = as.numeric(curr_file[,pos_col]),
                            pos1 = as.numeric(curr_file[,pos_col]) + 1,
                            snpid = paste0(curr_file[,chrom_col], ":", curr_file[,pos_col], ":", curr_file[,a1_col], ":", curr_file[,a2_col]),
                            freq1 = curr_file[,freq_col])

        }

        ### Catching chr23 if there is still any

        if ("chr23" %in% unique(curr_snp_set$chr)){

            curr_snp_set$chr = gsub("chr23", "chrX", curr_snp_set$chr)

        }

        ### rbinding here

        cat(paste0("\nSNP set created. CHECK ", i, ".", j, ".3 COMPLETE.\n\n"))

        cat(paste0("\nCurrent read nrow = ", nrow(curr_snp_set), "\n\n"))

        cohort_output = rbind(cohort_output, curr_snp_set)

        cat(paste0("\nCurrent output nrow = ", nrow(cohort_output), "\n\n"))

        print(summary(cohort_output))

        ### Trying to remove any examples of duplicated SNPs
        ### More done later

        if (any(duplicated(cohort_output$snpid))){

            cohort_output = cohort_output[!(duplicated(cohort_output$snpid)),]

        }

        cat(paste0("\nCurrent unique output nrow = ", nrow(cohort_output), "\n\n"))

        print(summary(cohort_output))

        cat(paste0("\nSNP set bound and uniqued. CHECK ", i, ".", j, ".4 COMPLETE.\n\n"))

    }

    ### Create chr_order column

    cohort_output$chr_order = NA

    ### Assigning values to string chromosomes for ordering

    cohort_output$chr_order[grepl("X", cohort_output$chr)] = 23
    cohort_output$chr_order[grepl("Y", cohort_output$chr)] = 24
    cohort_output$chr_order[grepl("MT|M", cohort_output$chr)] = 25

    ### For numeric chromosomes just just value

    idx = is.na(cohort_output$chr_order)
    cohort_output$chr_order[idx] = as.numeric(gsub("^chr", "", cohort_output$chr[idx]))

    ### Sort by chr and pos in that order

    cohort_output = cohort_output[order(cohort_output$chr_order, cohort_output$pos), ]

    cat(paste0("\nThis is the ordered file output. CHECK ", i, ".5 COMPLETE.\n\n"))
    print(cohort_output[sample(1:nrow(cohort_output), 20)],)

    ### Remove chr_order column

    cohort_output$chr_order = NULL

    fwrite(cohort_output, paste0(all_cohorts[i], "/", curr_cohort, "_liftOver_input.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

    cat(paste0("\nLiftOver file written. LAST CHECK ", i, " COMPLETE.\n\n==================\n\n"))

}
