##### =========================== #####

### Setup environment

### Script is designed to make PGS weights using LDpred2
### This is an alternative version of step1 which uses pre-computed LD blocks

##### =========================== #####

library(data.table)
library(bigsnpr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(bigreadr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(rmio, lib.loc = "/projappl/project_2007428/RPackages_421/")

Sys.info()

start_time = Sys.time()

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
weightDir = paste0(mainDir, "outputs/prs/ldpred2_weights/")
ldblockDir = paste0(mainDir, "processing/misc_data/ldpred_ld_ref/")
tempDir = paste0(mainDir, "tmpdir/ldpred_tmps/")

hg37_snps_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt"

NCORES = nb_cores()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Had to be me. Someone else would have gotten it wrong.")

### Setting up meta-analysis list of metas

analysis_list = Sys.glob(paste0(metaDir, "*ALL*.gz"))
analysis_list = analysis_list[!(grepl("pruned", analysis_list))]

### Information for the variants provided in the LD reference

hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_plus_snp_list.rds")
hapmap3_snps = readRDS(hapmap3_snps_file)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  
  print("Script obtained arguments from command line incorrectly. Please enter run number.")
  stop("SCRIPT ERROR 0. Nae args.")
  
}

### Test example:
### args = c("1")

run_no = as.numeric(args[1])

### Selection and reading meta-analysis

curr_analysis = analysis_list[run_no]
phenotype = strsplit(basename(curr_analysis), "_metal")[[1]][1]

curr_sumstats = fread(curr_analysis, data.table = FALSE, tmpdir = tempDir)

heading(paste0("Running LDpred2 on current meta: ", basename(curr_analysis)))

### Convert sumstats to LDpred2 format and reduce to HapMap3+
### Also add hg37 positions

curr_sumstats = curr_sumstats[curr_sumstats$MarkerName %in% hapmap3_snps$rsid,]

hg38_positions = curr_sumstats[,c("MarkerName", "Position")]
colnames(hg38_positions) = c("rsid", "hg38_pos")

curr_sumstats = merge(curr_sumstats, hapmap3_snps[,c("rsid", "pos")], by.x = "MarkerName", by.y = "rsid") # hapmap3+ list has hg37 positions
curr_sumstats = unique(curr_sumstats) # may introduce duplicates, remove them

### Get betas and SEs

curr_sumstats$beta1 = curr_sumstats$Zscore / sqrt((2 * curr_sumstats$Freq1) * (1 - curr_sumstats$Freq1) * (curr_sumstats$Weight + (curr_sumstats$Zscore^2)))
curr_sumstats$se = 1 / sqrt((2 * curr_sumstats$Freq1) * (1 - curr_sumstats$Freq1) * (curr_sumstats$Weight + (curr_sumstats$Zscore^2)))

### Rename to input columns

curr_sumstats = data.frame(chr = curr_sumstats$Chromosome,
                            pos = curr_sumstats$pos,
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

cat(paste0("\nCurrent sumstats (", phenotype, ") are these (", nrow(curr_sumstats), " rows):\n\n"))
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

setup_time = Sys.time()

heading("Script set-up complete")
cat(paste0("\nTime taken:\n", setup_time - start_time, "\n\n"))


##### =========================== #####

### Starting analysis

##### =========================== #####

### Match SNPs between sumstats and hapmap+
### Columns are named the same

info_snp = snp_match(curr_sumstats, hapmap3_snps)

### LDpred2 example script suggests filtering based on AF
### Unsure how important this is
### For GenCost would result in removal of >50% of SNPs

filter_mode = FALSE

if (isTRUE(filter_mode)){

    ### Filters based on standard deviation of AF and betas

    ### better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)

    sd_ldref = with(info_snp, sqrt(2 * MAF * (1 - MAF)))
    sd_ss = with(info_snp, 2 / sqrt(n_eff * beta_se^2))

    is_bad = sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

    #library(ggplot2)
    #qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
    #theme_bigstatsr() +
    #coord_equal() +
    #scale_color_viridis_d(direction = -1) +
    #geom_abline(linetype = 2, color = "red") +
    #labs(x = "Standard deviations derived from allele frequencies of the LD reference",
    #    y = "Standard deviations derived from the summary statistics",
    #    color = "Removed?")

    df_beta = info_snp[!is_bad, ]

} else {

    df_beta = info_snp
    
}

cat(paste0("\nLDpred2 input sumstats (", phenotype, ") are these (", nrow(df_beta), " rows):\n\n"))
print(head(df_beta))
print(summary(df_beta))

### Open a temporary file

temp_file_ref = tempfile(tmpdir = tempDir) # Required for bigsnpr

for (chr in 1:22) {

    cat(paste0("\nWorking on chromosome ", chr, ".\n\n"))

    ## indices in df_beta
    
    ind.chr = which(df_beta$chr == chr)

    ## indices in 'hapmap3_snps'

    ind.chr2 = df_beta$`_NUM_ID_`[ind.chr]

    ## indices in 'corr_chr'

    ind.chr3 = match(ind.chr2, which(hapmap3_snps$chr == chr))

    cat(paste0("\nIndices created.\n\n"))

    corr_chr = readRDS(paste0(ldblockDir, "LD_with_blocks_chr", chr, ".rds"))
    corr_chr = corr_chr[ind.chr3, ind.chr3]

    cat(paste0("\nCorrelations read.\n\n"))

    if (chr == 1) {

        corr = as_SFBM(corr_chr, temp_file_ref, compact = TRUE)

    } else {

        corr$add_columns(corr_chr, nrow(corr))

    }

    cat(paste0("\nCorrelations written.\n..........\n\n"))
    
}

corr_time = Sys.time()
cat(paste0("\nCorrelation time taken:\n", corr_time - start_time, "\n\n"))

### Heritability estimation of LD score regression
### to be used as a starting value in LDpred2-auto

ldsc = with(df_beta, bigsnpr::snp_ldsc(ld, ld_size = nrow(hapmap3_snps),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                blocks = NULL, # skip parallelisation BUT does not calculate LDSC SEs - default is 200
                                ncores = NCORES))

h2_est = ldsc[["h2"]]

cat(paste0("\n\nh2 estimate for trait, ", phenotype, ": ", h2_est, "\n\n"))

if (h2_est < 0){

    warning(paste0("This phenotype produces a -ve h2. Ending here.\n======\n\n"))

    stop()

}

### LDpred2-auto

multi_auto = snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES) # 5 min

cat("\nMulti_auto summary:\n")
print(summary(multi_auto))

### Filter for best chains and average remaining ones -> the effects sizes of your polygenic score

range = sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep = which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

cat("\n-------\n")
print(summary(range))
cat("\n-------\n")
print(summary(keep))
cat("\n-------\n")

beta_auto = rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

heading("LDpred2 auto weights made")

weight_time = Sys.time()
cat(paste0("\nWeight time taken:\n", weight_time - start_time, "\n\n"))

beta_auto = data.frame(rsid = df_beta$rsid,
                        chr = df_beta$chr,
                        hg37_pos = df_beta$pos,
                        a1 = df_beta$a1,
                        a0 = df_beta$a0,
                        ldpred2_beta_auto = beta_auto)

### LDpred2 Auto complete
### Adding values to output

heading(paste0("LDpred2_auto weights made for phenotype: ", phenotype))
print(head(beta_auto))

ldpred2_output = rbind(ldpred2_output, beta_auto)

system(paste0("rm ", temp_file_ref, "*"))

print(head(hg38_positions))
print(summary(hg38_positions))

print(head(ldpred2_output))
print(summary(ldpred2_output))

ldpred2_output = merge(ldpred2_output, hg38_positions, by = "rsid")

fwrite(ldpred2_output, paste0(weightDir, phenotype, "_ldpred2_precomputed_weights_pos38_ALLCHRs.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

end_time = Sys.time()
cat(paste0("\nTotsl time taken:\n", end_time - start_time, "\n\n"))

heading("R SCRIPT COMPLETE.")







             