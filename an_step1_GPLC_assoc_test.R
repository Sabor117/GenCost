##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(speedglm)
library(dplyr)
library(stringr)
library(ggplot2)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
phenoDir = paste0(mainDir, "outputs/cost_phenotypes/")
prsDir = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/"


##### =========================== #####

### Data setup

##### =========================== #####

### Read GPLC PRS file

prsFile = paste0(prsDir, "ukbb_78537_gplc_prs.csv")
prs_translationsFile = paste0(prsDir, "ukbb_78537_gplc_translations.tsv")

prs_read = fread(prsFile, data.table = FALSE)

### Extract PRS names from translations file

prs_names = fread(prs_translationsFile, data.table = FALSE)

### Split based on ()
### Some PRS have two sets of () so split twice

prs_names_split = str_split_fixed(prs_names$Description, "\\(", 3)
prs_names_split[,2][prs_names_split[,3] != ""] = prs_names_split[,3][prs_names_split[,3] != ""]
prs_names_split[,2] = gsub("\\)", "", prs_names_split[,2])

prs_names$shortname = prs_names_split[,2]

### Rename columns based on new abbreviation names

colnames(prs_read)[-1] = prs_names$shortname[match(colnames(prs_read)[-1], prs_names$UDI)]

### Reading phenotype files

gp_phenosFile = paste0(phenoDir, "FINAL_gp_based_phenotype_frame.tsv")
hes_phenosFile = paste0(phenoDir, "FINAL_hes_based_phenotype_frame.tsv")
prescript_phenosFile = paste0(phenoDir, "FINAL_prescription_based_phenotype_frame.tsv")

gp_phenos = fread(gp_phenosFile, data.table = FALSE)
hes_phenos = fread(hes_phenosFile, data.table = FALSE)
prescript_phenos = fread(prescript_phenosFile, data.table = FALSE)

### Get all IIDs across phenotype files

all_phenos = data.frame(eid = unique(c(hes_phenos$eid, gp_phenos$eid, prescript_phenos$eid)))

### Get log phenotypes + actual cost phenotypes (for later)

all_phenos = base::merge(all_phenos, hes_phenos[,c("eid", "hes_costs_year", "log_hes_costs_year")], by = "eid", all = TRUE)
all_phenos = base::merge(all_phenos, gp_phenos[,c("eid", "log_gp_costs_year", "log_total_costs_year", "gp_costs_year", "total_costs_year")], by = "eid", all = TRUE)
all_phenos = base::merge(all_phenos, prescript_phenos[,c("eid", "log_prescription_costs_year", "prescription_costs_year")], by = "eid", all = TRUE)

### Keep only those IIDs which overlap between phenotypes + PRS

all_phenos_overlap = all_phenos[all_phenos$eid %in% prs_read$eid,]
prs_overlap = prs_read[prs_read$eid %in% all_phenos_overlap$eid,]

### Read covariates file

covariatesFile = paste0(phenoDir, "gwas_cost_covariates.tsv")

### Remove the unncessary columns

covariates_read = fread(covariatesFile, data.table = FALSE)
covariates_read$age_end_squared = NULL
covariates_read$age_end_times_sex = NULL
covariates_read$SuperPop = NULL
covariates_read$year_of_birth = NULL
covariates_read$FID = NULL

### Rename IID for merging purposes
### Limit IIDs to those in PRS + phenotypes

colnames(covariates_read)[1] = "eid"

covariates_overlap = covariates_read[covariates_read$eid %in% all_phenos_overlap$eid,]


##### =========================== #####

### Running analysis

##### =========================== #####

for (nprs in 2:ncol(prs_overlap)){

    currPRS = colnames(prs_overlap)[nprs]

    prs_overlap[,nprs] = scale(prs_overlap[,nprs])

    colnames(prs_overlap)[nprs] = currPRS

}


##### =========================== #####

### Running analysis

##### =========================== #####

### Get list of covariates for formula

covariates = colnames(covariates_overlap)[-1]
covariates = paste(covariates, collapse = " + ")

### Set up base frame for output

outFrame = data.frame(matrix(ncol = 8, nrow = 0))
colnames(outFrame) = c("trait", "score", "median_cost", "top_percentile_cost", "n", "beta", "se", "p")

for (ntrait in which(grepl("log", colnames(all_phenos_overlap)))){

    currtrait = colnames(all_phenos_overlap)[ntrait]

    cat(paste0("Running on trait ", ntrait, " of ", ncol(all_phenos_overlap), ": ", currtrait, "\n\n\n================\n\n\n"))

    nolog_trait = gsub("log_", "", currtrait)
    nolog_trait_num = which(colnames(all_phenos_overlap) == nolog_trait)

    curr_pheno_frame = all_phenos_overlap[,c(1, ntrait, nolog_trait_num)]
    curr_pheno_frame = na.omit(curr_pheno_frame)

    median_cost = median(curr_pheno_frame[,nolog_trait])

    phenotype_frame = merge(curr_pheno_frame, covariates_overlap, by = "eid")

    for (nprs in 2:ncol(prs_overlap)){

        currPRS = colnames(prs_overlap)[nprs]

        cat(paste0("Running on trait ", nprs-1, " of ", ncol(prs_overlap)-1, ": ", currPRS, "\n======\n\n"))

        curr_prs_frame = prs_overlap[,c(1, nprs)]

        assoc_frame = merge(phenotype_frame, curr_prs_frame, by = "eid")

        assoc_frame$prs_twnt = cut(assoc_frame[,currPRS],
                                    breaks = quantile(assoc_frame[,currPRS], prob = 0:20/20, names = FALSE), 
                                    labels = 1:20, include = TRUE)

        assoc_frame$prs_ten = cut(assoc_frame[,currPRS],
                                    breaks = quantile(assoc_frame[,currPRS], prob = 0:10/10, names = FALSE), 
                                    labels = 1:10, include = TRUE)

        top_percentile = assoc_frame[assoc_frame$prs_ten == 10,]
        middle_percentile = assoc_frame[assoc_frame$prs_ten == 4 | assoc_frame$prs_ten == 5 | assoc_frame$prs_ten == 6,]

        top_percentile_cost = median(top_percentile[,gsub("log_", "", nolog_trait)])
        middle_percentile_cost = median(middle_percentile[,gsub("log_", "", nolog_trait)])
        median_cost = median(assoc_frame[,gsub("log_", "", nolog_trait)])

        top_percentile_mean_cost = mean(top_percentile[,gsub("log_", "", nolog_trait)])
        middle_percentile_mean_cost = mean(middle_percentile[,gsub("log_", "", nolog_trait)])
        mean_cost = mean(assoc_frame[,gsub("log_", "", nolog_trait)])

        assoc_formula = as.formula(paste0(nolog_trait, " ~ ", covariates, " + ", currPRS))

        reg_model = glm(formula = assoc_formula, data = assoc_frame)

        reg_summary = summary(reg_model)

        reg_coeff = reg_summary$coefficients

        sample_size = nrow(assoc_frame)

        outRow = data.frame(trait = nolog_trait,
                    score = currPRS,
                    median_cost = median_cost,
                    median_middle_decile = middle_percentile_cost,
                    median_top_decile = top_percentile_cost,
                    mean_cost = mean_cost,
                    mean_middle_decile = middle_percentile_mean_cost,
                    mean_top_decile = top_percentile_mean_cost,
                    n = sample_size,
                    beta = reg_coeff[26,1],
                    se = reg_coeff[26,2],
                    p = reg_coeff[26,4])

        outFrame = rbind(outFrame, outRow)

    }
}

outFrame$trait = case_when(outFrame$trait == "hes_costs_year" ~ "Inpatient costs",
                              outFrame$trait == "gp_costs_year" ~ "Primary care costs",
                              outFrame$trait == "total_costs_year" ~ "Total costs",
                              outFrame$trait == "prescription_costs_year" ~ "Prescription costs",
                              .default = as.character(outFrame$trait))

outFrame$signif = ifelse(outFrame$p < 0.05/length(unique(outFrame$score)), "*", "")

fwrite(outFrame, paste0(mainDir, "outputs/prs/gplc_cost_assoc/cost_gplc_prs_assoc_final_v2.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

legend_colours = c("#006d2c", "#2166ac", "#ff7b00")
legend_labels = c("Primary costs", "Inpatient costs", "Prescription costs")
legend_labels = paste0(legend_labels, ", n = ", c(43712, 102989, 38400))
legend_labels = paste0(legend_labels, ", median cost = €", c(281.90, 233.40, 34.46))

png(paste0(mainDir, "outputs/figures/misc/gplc_assoc_barplot_final.png"), height = 2000, width = 4000, res = 300)

ggplot(outFrame, aes(x = beta, y = reorder(factor(score), beta), fill = factor(trait))) +
    geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" to group by PRS
    labs(x = "Beta", y = "PRS") +
    theme_linedraw() +
    scale_fill_manual(values = legend_colours, name = "Phenotype", labels = legend_labels) +  # Set custom colors for phenotypes
    geom_errorbar(aes(xmin = beta - se, xmax = beta + se, color = signif, group = factor(trait)), position = position_dodge(width = 0.9), size = 0.2) +
    scale_color_manual(values = c("black", "red"), name = "Significance") +
    theme(legend.position = "right",
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()

dev.off()
