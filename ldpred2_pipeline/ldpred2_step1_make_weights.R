##### =========================== #####

### Setup environment

### Script is designed to make PGS weights using LDpred2
### This is the primary version of step1 which calculates LD blocks from a reference genotype
### Script has been designed to be manually altered to run either genome-wide or chromosome-by-chromosome

##### =========================== #####

library(data.table)
library(magrittr)
library(bigsnpr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(bigreadr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(rmio, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(bigutilsr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(stringr)
library(dplyr)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
weightDir = paste0(mainDir, "outputs/prs/ldpred2_weights/")
tempDir = paste0(mainDir, "tmpdir/ldpred_tmps/")

refGenome = "/scratch/project_2007428/data/processing/ukbb_78537/genotypes/white_british_30k_reference_panel/ukb22828_chr%%%_30k_random_unrelated_white_british"
refGenomeALLCHR = "/scratch/project_2007428/data/processing/ukbb_78537/genotypes/white_british_30k_reference_panel/pruned_binary/ukb22828_ALLCHR_pruned"

hg37_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt"

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


##### =========================== #####

### Data setup

##### =========================== #####

### Selecting and reading meta-analysis based on list of metas

analysis_list = Sys.glob(paste0(metaDir, "*ALL*.gz"))
analysis_list = analysis_list[!(grepl("pruned", analysis_list))]

curr_analysis = analysis_list[run_no]

heading(paste0("Running LDpred2 on current meta: ", basename(curr_analysis)))

phenotype = strsplit(basename(curr_analysis), "_metal")[[1]][1]

curr_sumstats = fread(curr_analysis, data.table = FALSE, tmpdir = tempDir)

#colexample = fread(curr_analysis, data.table = FALSE, nrow = 5)
#names(curr_sumstats) = names(colexample)

### Reading hapmap3+ SNPs

hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_plus_snp_list.rds")
hapmap3_snps = readRDS(hapmap3_snps_file)

### Read hg37 file

hg37_snps = fread(hg37_snps_file, data.table = FALSE)

hg37_snps = hg37_snps[,c(2,4)]
colnames(hg37_snps) = c("hg37_pos", "MarkerName")


##### =========================== #####

### Starting analysis

##### =========================== #####

### Convert sumstats to LDpred2 format and reduce to HapMap3+
### Also add hg37 positions

curr_sumstats = curr_sumstats[curr_sumstats$MarkerName %in% hapmap3_snps$rsid,]

hg38_positions = curr_sumstats[,c("MarkerName", "Position")]
colnames(hg38_positions) = c("rsid", "hg38_pos")

curr_sumstats = merge(curr_sumstats, hg37_snps, by = "MarkerName")

curr_sumstats$se = sqrt( 1 / ((curr_sumstats$Zscore^2) + curr_sumstats$Weight)) / sqrt(2 * curr_sumstats$Freq1 * (1 - curr_sumstats$Freq1))
curr_sumstats$beta1 = curr_sumstats$Zscore * curr_sumstats$se

curr_sumstats = data.frame(chr = curr_sumstats$Chromosome,
                            pos = curr_sumstats$hg37_pos,
                            rsid = curr_sumstats$MarkerName,
                            a1 = curr_sumstats$Allele1,
                            a0 = curr_sumstats$Allele2,
                            n_eff = curr_sumstats$Weight,
                            MAF = curr_sumstats$Freq1,
                            p = curr_sumstats[,"P-value"],
                            beta_se = curr_sumstats$se,
                            beta = curr_sumstats$beta1
                            )

curr_sumstats$a1 = toupper(curr_sumstats$a1)
curr_sumstats$a0 = toupper(curr_sumstats$a0)

cat(paste0("\nCurrent sumstats (", phenotype, ") are these:\n\n"))
print(head(curr_sumstats))
print(summary(curr_sumstats))

if (any(is.infinite(curr_sumstats$beta_se))){

    cat(paste0("\nInf values detected. Removing them.\n\n"))

    numINF = nrow(curr_sumstats[is.infinite(curr_sumstats$beta_se),])

    print(curr_sumstats[is.infinite(curr_sumstats$beta_se),])

    curr_sumstats = curr_sumstats[!(is.infinite(curr_sumstats$beta_se)),]

    cat(paste0("\n", numINF, " rows removed.\n\n"))
    print(summary(curr_sumstats))

}


##### =========================== #####

### LDpred2 analysis

##### =========================== #####

### Chromosome-by-chromosome or total

chr_by_chr = FALSE

# Get maximum amount of cores

NCORES = nb_cores()

cat(paste0("\nWe have ", NCORES, " to work with.\n\n\n"))

### Set up output

ldpred2_output = as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(ldpred2_output) = c("rsid", "chr", "hg37_pos", "a1", "a0", "ldpred2_beta_auto")

if (isTRUE(chr_by_chr)){

    cat(paste0("\nWe are running chromosome by chromosome.\n\n\n"))

    for (chrom in 1:22){

        # Open a temporary file

        temp_file_ref = tempfile(tmpdir = tempDir) # Required for bigsnpr

        ### Reading genotypes into bigsnp file

        if (chrom != 23){

            curr_bed_chr = gsub("%%%", chrom, refGenome)

        } else if (chrom == 23){

            curr_bed_chr = gsub("%%%", "X", refGenome)

        }

        ### Reduce sumstats per chromosome

        curr_sumstats_chrom = curr_sumstats[curr_sumstats$chr == chrom,]

        cat(paste0("\nRunning on chromosome: ", chrom))
        head(curr_sumstats_chrom)
        summary(curr_sumstats_chrom)

        ### LDpred2 uses beta/se for LDSD - this cannot be negative

        cat(paste0("\nBeta / beta_se = \n\n"))
        print(head(((curr_sumstats_chrom$beta / curr_sumstats_chrom$beta_se)^2)))

        if (any(((curr_sumstats_chrom$beta / curr_sumstats_chrom$beta_se)^2) < 0)){

            chi_2 = (curr_sumstats_chrom$beta / curr_sumstats_chrom$beta_se)^2

            neg_chi_2 = which(chi_2 < 0)

            cat(paste0("\nThe following chi2 were negative. These will be removed. There were negative ", length(neg_chi_2), " SNPs.\n\n"))
            print(head(curr_sumstats_chrom[neg_chi_2,]))

            curr_sumstats_chrom = curr_sumstats_chrom[-c(neg_chi_2),]

        }

        ### Read bed file for chromosome

        snp_readBed2(paste0(curr_bed_chr, ".bed"),
                    backingfile = temp_file_ref)

        bigsnp_ref = snp_attach(paste0(temp_file_ref, ".rds"))

        ### From bigsnp file extract
        ### Genotypes, chromosomes, positions

        ref_G = bigsnp_ref$genotypes
        ref_CHR = as.numeric(bigsnp_ref$map$chromosome)
        ref_POS = bigsnp_ref$map$physical.pos
        ref_y = bigsnp_ref$fam$affection - 1

        cat(paste0("\nLDpred2 setup complete.\n======\n\n"))

        ### Extract required information from bigsnp file into a map

        ref_map = bigsnp_ref$map[-(2:3)]
        names(ref_map) = c("chr", "pos", "a1", "a0")
        ref_map$chr = as.numeric(ref_map$chr)
        ref_info_snp = snp_match(curr_sumstats_chrom, ref_map)
        ref_POS2 = snp_asGeneticPos(ref_CHR, ref_POS)

        cat(paste0("\nMap created.\n======\n\n"))

        ### Indices in info_snp

        ref_df_beta = ref_info_snp[, c("beta", "beta_se", "n_eff")]

        ### Indices in G

        ref_ind.chr2 = ref_info_snp$`_NUM_ID_`

        ### Correlation matrix

        ref_corr0 = snp_cor(ref_G,
                            ind.col = ref_ind.chr2,
                            infos.pos = ref_POS2[ref_ind.chr2],
                            size = 3 / 1000,
                            ncores = NCORES)

        ### Testing for NAs in correlation matrix

        row.names(ref_corr0) = 1:nrow(ref_corr0)

        if (any(is.na(ref_corr0))){

            heading("NAs found in correlation matrix, removing SNPs.")

            na_rem = TRUE

            first_na_test = apply(is.na(ref_corr0), 2, which)

            na_index = as.numeric(names(which(lengths(first_na_test) == (nrow(ref_corr0) - 1))))

            if (any(is.na(ref_corr0[-c(na_index), -c(na_index)]))){

                testing_corr0 = ref_corr0[-c(na_index), -c(na_index)]

                repeat{

                    second_na_test = apply(is.na(testing_corr0), 2, which)

                    row_to_remove = as.numeric(names(sort(table(unlist(second_na_test)), decreasing = TRUE)[1]))
                    
                    testing_corr0 = testing_corr0[-row_to_remove, -row_to_remove]
                    
                    if (!(any(is.na(testing_corr0)))){

                        break

                    }
                }

                na_index = as.numeric(row.names(ref_corr0)[!(row.names(ref_corr0) %in% row.names(testing_corr0))])
            }

            ref_ind.chr2 = ref_ind.chr2[-na_index]

            ref_corr0 = ref_corr0[-na_index, -na_index]
                    
            ref_df_beta = ref_df_beta[-na_index,]

            heading(paste0("NAs removed from correlation matrix. Removed: ", length(na_index), " SNPs."))

        } else {

            na_rem = FALSE

        }

        tmp_file_corr = tempfile(tmpdir = tempDir)

        ref_corr = bigsparser::as_SFBM(ref_corr0,
                                        tmp_file_corr,
                                        compact = TRUE)
        ref_ld = Matrix::colSums(ref_corr0^2)

        cat(paste0("\nCorrelations made.\n======\n\n"))


        ###============

        ### Preparations complete
        ### Running LDpred2

        ###============

        ### LDSC

        ldsc = with(ref_df_beta,
                    snp_ldsc(ref_ld,
                                length(ref_ld),
                                chi2 = ((beta / beta_se)^2),
                    sample_size = n_eff, blocks = NULL))

        h2_est = ldsc[["h2"]]

        cat(paste0("\n\nh2 estimate for trait, ", phenotype, ", on chromosome ", chrom, ": ", h2_est, "\n\n"))

        if (h2_est < 0){

            warning(paste0("Chromosome ", chrom, " produces -ve h2. This chromosome will be skipped.\n======\n\n"))

            next

        }

        ### LDPred2 Auto model

        multi_auto = snp_ldpred2_auto(corr = ref_corr,
                                        df_beta = ref_df_beta,
                                        h2_init = h2_est,
                                        vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                                        allow_jump_sign = FALSE,
                                        shrink_corr = 0.95,
                                        ncores = 1)

        beta_auto = rowMeans(sapply(multi_auto, function(auto) auto$beta_est))

        beta_auto = data.frame(rsid = ref_info_snp$rsid,
                                chr = ref_info_snp$chr,
                                hg37_pos = ref_info_snp$pos,
                                a1 = ref_info_snp$a1,
                                a0 = ref_info_snp$a0,
                                ldpred2_beta_auto = beta_auto)
        
        ### LDpred2 Auto complete
        ### Adding values to output

        heading(paste0("LDpred2_auto weights made for chromosome: ", chrom))
        print(head(beta_auto))

        ldpred2_output = rbind(ldpred2_output, beta_auto)

        system(paste0("rm ", temp_file_ref, "*"))
        system(paste0("rm ", tmp_file_corr, "*"))

        heading("All chromosomes complete.")
    }
} else if (isFALSE(chr_by_chr)){

    cat(paste0("\nWe are NOT running chromosome by chromosome.\n\n\n"))

    ### Open a temporary file

    temp_file_ref = tempfile(tmpdir = tempDir) # Required for bigsnpr

    cat(paste0("\nRunning on ALL CHROMOSOMES."))
    head(curr_sumstats)
    summary(curr_sumstats)

    ### LDpred2 uses beta/se for LDSD - this cannot be negative

    cat(paste0("\nBeta / beta_se = \n\n"))
    print(head(((curr_sumstats$beta / curr_sumstats$beta_se)^2)))

    if (any(((curr_sumstats$beta / curr_sumstats$beta_se)^2) < 0)){

        chi_2 = (curr_sumstats$beta / curr_sumstats$beta_se)^2

        neg_chi_2 = which(chi_2 < 0)

        cat(paste0("\nThe following chi2 were negative. These will be removed. There were negative ", length(neg_chi_2), " SNPs.\n\n"))
        print(head(curr_sumstats[neg_chi_2,]))

        curr_sumstats = curr_sumstats[-c(neg_chi_2),]

    }

    curr_sumstats = curr_sumstats[order(curr_sumstats$chr, curr_sumstats$pos),]

    ### Read bed file for allchrs

    snp_readBed2(paste0(refGenomeALLCHR, ".bed"),
                backingfile = temp_file_ref)

    bigsnp_ref = snp_attach(paste0(temp_file_ref, ".rds"))

    ### From bigsnp file extract
    ### Genotypes, chromosomes, positions

    ref_G = bigsnp_ref$genotypes
    ref_CHR = as.numeric(bigsnp_ref$map$chromosome)
    ref_POS = bigsnp_ref$map$physical.pos
    ref_y = bigsnp_ref$fam$affection - 1

    cat(paste0("\nLDpred2 setup complete.\n======\n\n"))

    ### Extract required information from bigsnp file into a map

    ref_map = bigsnp_ref$map[-(2:3)]
    names(ref_map) = c("chr", "pos", "a1", "a0")
    ref_map$chr = as.numeric(ref_map$chr)
    ref_info_snp = snp_match(curr_sumstats, ref_map)
    ref_POS2 = bigsnpr::snp_asGeneticPos(ref_CHR, ref_POS)

    cat(paste0("\nMap created.\n======\n\n"))

    ### Indices in info_snp

    ref_df_beta = ref_info_snp[, c("beta", "beta_se", "n_eff")]

    ### Indices in G

    ref_ind.chr2 = ref_info_snp$`_NUM_ID_`

    ord_ind.chr2 = order(ref_POS2[ref_ind.chr2])

    cat(paste0("\nORD_IND created.\n======\n\n"))
    print(head(ord_ind.chr2))

    ### Correlation matrix

    ref_corr0 = snp_cor(ref_G,
                        ind.col = ref_ind.chr2[ord_ind.chr2],
                        infos.pos = ref_POS2[ref_ind.chr2[ord_ind.chr2]],
                        size = 3 / 1000,
                        ncores = NCORES)

    cat(paste0("\nref_corr0 created.\n======\n\n"))
    print(names(ref_corr0))
    print(head(ref_corr0$sbk))

    ### Testing for NAs in correlation matrix

    row.names(ref_corr0) = 1:nrow(ref_corr0)

    if (any(is.na(ref_corr0))){

        heading("NAs found in correlation matrix, removing SNPs.")

        na_rem = TRUE

        first_na_test = apply(is.na(ref_corr0), 2, which)

        na_index = as.numeric(names(which(lengths(first_na_test) == (nrow(ref_corr0) - 1))))

        if (any(is.na(ref_corr0[-c(na_index), -c(na_index)]))){

            testing_corr0 = ref_corr0[-c(na_index), -c(na_index)]

            repeat{

                second_na_test = apply(is.na(testing_corr0), 2, which)

                row_to_remove = as.numeric(names(sort(table(unlist(second_na_test)), decreasing = TRUE)[1]))
                
                testing_corr0 = testing_corr0[-row_to_remove, -row_to_remove]
                
                if (!(any(is.na(testing_corr0)))){

                    break

                }
            }

            na_index = as.numeric(row.names(ref_corr0)[!(row.names(ref_corr0) %in% row.names(testing_corr0))])
        }

        ref_ind.chr2 = ref_ind.chr2[-na_index]

        ref_corr0 = ref_corr0[-na_index, -na_index]
                
        ref_df_beta = ref_df_beta[-na_index,]

        heading(paste0("NAs removed from correlation matrix. Removed: ", length(na_index), " SNPs."))

    } else {

        na_rem = FALSE

    }

    tmp_file_corr = tempfile(tmpdir = tempDir)

    ref_corr = bigsparser::as_SFBM(ref_corr0,
                                    tmp_file_corr,
                                    compact = TRUE)
    ref_ld = Matrix::colSums(ref_corr0^2)

    cat(paste0("\nCorrelations made.\n======\n\n"))


    ###============

    ### Preparations complete
    ### Running LDpred2

    ###============

    ### LDSC

    ldsc = with(ref_df_beta,
                snp_ldsc(ref_ld,
                            length(ref_ld),
                            chi2 = ((beta / beta_se)^2),
                sample_size = n_eff, blocks = NULL))

    h2_est = ldsc[["h2"]]

    cat(paste0("\n\nh2 estimate for trait, ", phenotype, ", on chromosome ", chrom, ": ", h2_est, "\n\n"))

    if (h2_est < 0){

        warning(paste0("Chromosome ", chrom, " produces -ve h2. This chromosome will be skipped.\n======\n\n"))

        next

    }

    ### LDPred2 Auto model

    multi_auto = snp_ldpred2_auto(corr = ref_corr,
                                    df_beta = ref_df_beta,
                                    h2_init = h2_est,
                                    vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                                    allow_jump_sign = FALSE,
                                    shrink_corr = 0.95,
                                    ncores = 1)

    beta_auto = rowMeans(sapply(multi_auto, function(auto) auto$beta_est))

    beta_auto = data.frame(rsid = ref_info_snp$rsid,
                            chr = ref_info_snp$chr,
                            hg37_pos = ref_info_snp$pos,
                            a1 = ref_info_snp$a1,
                            a0 = ref_info_snp$a0,
                            ldpred2_beta_auto = beta_auto)
    
    ### LDpred2 Auto complete
    ### Adding values to output

    heading(paste0("LDpred2_auto weights made for ALLCHRS."))
    print(head(beta_auto))

    system(paste0("rm ", temp_file_ref, "*"))
    system(paste0("rm ", tmp_file_corr, "*"))

    heading("All chromosomes complete.")

}

#fwrite(ldpred2_output, paste0(weightDir, phenotype, "_ldpred2_weights.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

print(head(hg38_positions))
print(summary(hg38_positions))

print(head(ldpred2_output))
print(summary(ldpred2_output))

ldpred2_output = merge(ldpred2_output, hg38_positions, by = "rsid")

fwrite(ldpred2_output, paste0(weightDir, phenotype, "_ldpred2_weights_pos38_ALLCHRs.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

heading("R SCRIPT COMPLETE.")

