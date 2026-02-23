################################################################################
#
# Set up script
#
################################################################################

### Packages

library(data.table)
library(dplyr)

### Files

meta_sumstats_file_list = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/*_metal_output_1.TBL.gz")

heritability_file_list = meta_sumstats_file_list[grepl("herit", meta_sumstats_file_list)]
meta_sumstats_file_list = meta_sumstats_file_list[!(grepl("herit", meta_sumstats_file_list))]
meta_sumstats_file_list = meta_sumstats_file_list[grepl("_noUKB", meta_sumstats_file_list)]


################################################################################
#
# Running
#
################################################################################

for (i in 1:length(meta_sumstats_file_list)){

    currfile_1 = heritability_file_list[i]

    curr_pheno = strsplit(basename(currfile_1), "_herit")[[1]][1]

    currfile_2 = meta_sumstats_file_list[grepl(curr_pheno, meta_sumstats_file_list)]
    outfile_name = currfile_2

    cat(paste0("\nWorking on phenotype: ", curr_pheno,
                "\nHeritability file selected = ", heritability_file_list[i],
                "\nPosition file selected = ", meta_sumstats_file_list[grepl(curr_pheno, meta_sumstats_file_list)],
                "\nReading files.\n\n"
                ))

    currfile_1 = fread(currfile_1, data.table = FALSE)
    currfile_2 = fread(currfile_2, data.table = FALSE)

    cat(paste0("\nHead of selected files:\n\n"))
    print(head(currfile_1))
    print(head(currfile_2))

    merge_cols = c("MarkerName", "Chromosome", "Position")
    currfile_2 = currfile_2[,merge_cols]

    missing_snps = length(which(currfile_2$MarkerName %in% currfile_1$MarkerName))

    cat(paste0("\nMerging files.",
                "\nNumber of rows in heritability file which are in position file = ", missing_snps,
                "\nNumber of rows in heritability file which are missing = ", nrow(currfile_2) - missing_snps,
                "\n\n"
                ))

    currfile_1 = merge(currfile_1, currfile_2, by = "MarkerName", all.x = TRUE)

    cat(paste0("\nFiles merged.\n\n"))
    print(head(currfile_1))

    fwrite(currfile_1, outfile_name, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip")

    pval = "P-value"

    currfile_pruned = currfile_1[currfile_1[,pval] < 0.05,]

    fwrite(currfile_pruned, gsub("metal_output_1", "metal_output_1_pruned", outfile_name),
            quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip")

    cat(paste0("\nFile pruned and written.",
                "\nRun ", i, " of ", length(meta_sumstats_file_list), " complete.\n\n========\n\n"))

}

system("rm /scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/*herit*")

cat("\n\nINTERNAL CHECK. ALL R RUNS COMPLETE. YOU'RE THE KING? WELL I DIDN'T VOTE FOR YOU.\n\n==============================\n========\n")

