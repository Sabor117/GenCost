##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(wesanderson, lib.loc = "/projappl/project_2007428/RPackages_421/")
options(scipen = 5)
set.seed(117)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gctaDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")
all_sumstatDir = paste0(mainDir, "processing/meta_sumstats/")

outDir = paste0(mainDir, "outputs/snp_forest_inputs_v4/")
figureDir = paste0(mainDir, "outputs/figures/top_snp_forest_plots/")

sumstat_list = Sys.glob(paste0(gctaDir, "*gcta_jma_out.txt"))
sumstat_list = sumstat_list[!(grepl("_no", sumstat_list))]
sumstat_list = sumstat_list[grepl("ALL", sumstat_list)]

snp_translation_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out"

for (i in 1:length(sumstat_list)){

    analysis_name = gsub("_gcta_jma_out.txt", "", basename(sumstat_list[i]))

    if (file.exists(paste0(outDir, analysis_name, "_top_snp_metas.txt"))){

        cat(paste0("\n====================\n\nPhenotype: ", gsub("_gcta_jma_out.txt", "", basename(sumstat_list[i])), " already run\n\n==========\n\n"))

        next

    }

    cat(paste0("\n====================\n\nPhenotype: ", gsub("_gcta_jma_out.txt", "", basename(sumstat_list[i])), "\n\n==========\n\n"))

    curr_sumstats = fread(sumstat_list[i], data.table = FALSE)

#    if (analysis_name == "IN_F" | analysis_name == "IN_M" | analysis_name == "ALL_F" | analysis_name == "ALL_M"){
#
#        cat(paste0("\nThis will take too long.\n"))
#
#        next
#
#    }

    sumstatDir = paste0(all_sumstatDir, analysis_name, "/")

    all_sumstat_files = Sys.glob(paste0(sumstatDir, "*_metal_input_sumstats.txt.gz"))

    colnames_set = c("snpid", "chr", "pos", "a1", "a0", "beta1", "freq1", "se", "p", "n", "cohort")

    snp_meta_frame = data.frame(matrix(ncol = 11, nrow = 0))
    colnames(snp_meta_frame) = colnames_set

    cat(paste0("\nSumstats available for this phenotype: ", basename(all_sumstat_files)))

    for (nsnp in 1:nrow(curr_sumstats)){

        currsnp = curr_sumstats$SNP[nsnp]

        cat(paste0("\n========\n\nNow working on SNP ", currsnp, " which is ", nsnp, " of ", nrow(curr_sumstats), "\n\n========\n\n"))

        for (nfile in 1:length(all_sumstat_files)){

            curr_cohort = strsplit(basename(all_sumstat_files[nfile]), "\\.")[[1]][1]

            if (curr_cohort == "UKB"){

                curr_cohort = paste0(curr_cohort, "_", strsplit(basename(all_sumstat_files[nfile]), "\\.")[[1]][6])

            }

            cat(paste0("\nChecking ", curr_cohort, " for ", currsnp, "\n\n"))

            curr_row = system(paste0("zcat ", all_sumstat_files[nfile], " | grep ", currsnp), intern = TRUE)

            if (length(curr_row) == 0){

                cat(paste0(currsnp, " not in cohort.\n"))

                next

            }

            curr_row = str_split_fixed(curr_row, "\t", 10)

            curr_row = as.data.frame(curr_row)

            curr_row$V11 = curr_cohort
            colnames(curr_row) = colnames_set

            snp_meta_frame = rbind(snp_meta_frame, curr_row)

        }

        snp_meta_row = data.frame(snpid = currsnp,
                                    chr = snp_meta_frame$chr[1],
                                    pos = snp_meta_frame$chr[2],
                                    a1 = curr_sumstats$refA[nsnp],
                                    a0 = NA,
                                    beta1 = curr_sumstats$b[nsnp],
                                    freq1 = curr_sumstats$freq[nsnp],
                                    se = curr_sumstats$se[nsnp],
                                    p = curr_sumstats$p[nsnp],
                                    n = curr_sumstats$n[nsnp],
                                    cohort = "Meta")

        snp_meta_frame = rbind(snp_meta_frame, snp_meta_row)

    }

    snp_meta_frame = snp_meta_frame[snp_meta_frame$snpid %in% curr_sumstats$SNP,]

    fwrite(snp_meta_frame, paste0(outDir, analysis_name, "_top_snp_MAF_0_001_meta_forests.txt"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

}