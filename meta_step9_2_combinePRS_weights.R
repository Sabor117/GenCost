##### =========================== #####

### Setup environment

### Script is combine the polygenic weights into one file
### Format is designed for MegaPRS

##### =========================== #####

### Packages

library(data.table)
library(dplyr)
library(stringr)
library(PathWAS, lib.loc = "/projappl/project_2007428/RPackages_421/")

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
weightDir = paste0(mainDir, "outputs/prs/megaprs_weights/")
tempDir = paste0(mainDir, "tmpdir/")

### PRS weight files

all_megaprs = Sys.glob(paste0(weightDir, "*effects"))
all_megaprs = all_megaprs[grepl("noFINN", all_megaprs) & grepl("BAYESR", all_megaprs)]

### Reading hapmap SNPs

hapmap_set = "HMPLUS"

if (hapmap_set == "HM3"){

    hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_snp_list.txt")
    hapmap3_snps = fread(hapmap3_snps_file)

    hapmap_name = "hapmap3"

    sumstat_dir = paste0(mainDir, "processing/megaPRS_processing/hapmap_3_sumstats/")

    prs = as.data.frame(hapmap3_snps[,c("rsid")])

} else if (hapmap_set == "HMPLUS"){

    hapmap3_snps_file = paste0(mainDir, "processing/misc_data/hm3_plus_snp_list.rds")
    hapmap3_snps = readRDS(hapmap3_snps_file)

    hapmap_name = "hapmap_plus"

    sumstat_dir = paste0(mainDir, "processing/megaPRS_processing/hapmap_plus_sumstats/")

    prs = as.data.frame(hapmap3_snps[,c("rsid")])

}



##### =========================== #####

### Starting combine

##### =========================== #####

cat(paste0("\nWake up Mister Freeman. Wake up and smell the ashes...\n\n==============\n\n"))

all_sumstat_files = Sys.glob(paste0(sumstat_dir, "*.txt.gz"))

individ_pgs = TRUE

for (n_prs in 1:length(all_megaprs)){

    curr_PRS = fread(all_megaprs[n_prs], data.table = FALSE)
    curr_set = strsplit(basename(all_megaprs[n_prs]), "\\.")[[1]][1]

    cat(paste0("\nThe right file... In the wrong place... ", curr_set, " can make all the difference in the world.\n====\n\n"))
    print(head(curr_PRS))
    print(summary(curr_PRS))

    sumstats_file = all_sumstat_files[grepl(tolower(gsub("_BAYESR", "", curr_set)), tolower(all_sumstat_files))]
    sumstats = fread(sumstats_file, data.table = FALSE)

    sumstats$Predictor = paste0(sumstats$chr, ":", sumstats$hg37_pos)

    curr_PRS = left_join(curr_PRS, sumstats[,c(1:4,12,7,9)], by = "Predictor")

    colnames(curr_PRS)[c(2,3,5)] = c(paste0(curr_set, "_a1"), paste0(curr_set, "_a0"), paste0(curr_set, "_megaprs"))

    if (isTRUE(individ_pgs)){

        weights_idx = which(grepl("_megaprs", colnames(curr_PRS)))
        curr_PRS = curr_PRS[!(is.na(curr_PRS[,weights_idx])),]

        fwrite(curr_PRS, paste0(weightDir, curr_set, "_effects_for_pgs.txt"),
                quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

        next

    }

    #curr_PRS$pos_hg37 = str_split_fixed(curr_PRS$Predictor, ":", 2)[,2]

    prs = left_join(prs, curr_PRS[,c("rsid", paste0(curr_set, "_a1"), paste0(curr_set, "_a0"), paste0(curr_set, "_megaprs"), "chr", "hg37_pos", "hg38_pos")], by = "rsid")

    if (n_prs > 1){

        ### Now fill in missing positions in prs using coalesce

        prs = prs %>% mutate(
                        chr = coalesce(chr.x, chr.y),
                        hg37_pos = coalesce(hg37_pos.x, hg37_pos.y),
                        hg38_pos = coalesce(hg38_pos.x, hg38_pos.y)
                    ) %>%
                    select(-matches("\\.x$"), -matches("\\.y$"))

    }

    cat(paste0("\nSuck it G-man!\n\n"))
    print(head(prs))
    print(summary(prs))

}

megaprs_cols = grepl("_megaprs$", names(prs))
prs = prs[!apply(is.na(prs[, megaprs_cols, drop = FALSE]), 1, all), ]

### Change column order

first_cols = c("rsid", "chr", "hg37_pos", "hg38_pos")

prs = prs[,c(first_cols, colnames(prs)[!(colnames(prs) %in% first_cols)])]

### Manipulation complete
### Writing

fwrite(prs, paste0(mainDir, "outputs/prs/combined_weights/megaprs_gencost_META_v4_BAYESR_weights.txt"),
        quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")




