### ---------------------------
###
### Script name: meta_step13_1_ukb_pgs_diseases.R
###
### Purpose of script: Test GenCost PGS against number of UKB diagnoses
###
### Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2026-06-09
###
### ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

options(scipen = 999)


### Start timer

script_start = Sys.time()

### Classic memes

heading = function(text){
  
  cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))
  
}

heading("Setup directories and input files")

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
pgsDir = paste0(mainDir, "outputs/prs/scores/")
figureDir = paste0(mainDir, "outputs/figures/misc/")

allPGS = Sys.glob(paste0(pgsDir, "*noUKB*BAYESR*allchrs.tsv"))
allPGS = c(allPGS, paste0(pgsDir, "INOUT_ALL_meta_V4_hapmap_plus_BAYESR_megaprs_allchrs.tsv"))

ukb_base_data_file = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv"
ukb_base_data = fread(ukb_base_data_file, data.table = FALSE)

covariatesFile = paste0(mainDir, "outputs/cost_phenotypes/gwas_cost_covariates.tsv")
covariates = fread(covariatesFile, data.table = FALSE)


#####################

### Obtain diagnosis counts

#####################

### Get long data frame

long_icd = ukb_base_data %>%
                select(eid, starts_with("41270")) %>%
                pivot_longer(
                    cols = starts_with("41270"),
                    names_to = "field",
                    values_to = "ICD10"
                ) %>%
                filter(!is.na(ICD10), ICD10 != "")

### Count of unique ICD10 codes per individual

full_counts = long_icd %>%
                    group_by(eid) %>%
                    summarise(
                        n_unique_icd10 = n_distinct(ICD10),
                        .groups = "drop"
                    )

### Some individuals have no diagnoses at all
### These were removed by the filter - add them back in

healthy_individuals = data.frame(eid = ukb_base_data$eid[!(ukb_base_data$eid %in% full_counts$eid)],
                                    n_unique_icd10 = 0)

full_counts = rbind(full_counts, healthy_individuals)

### Then get version of ICD10 counts for broader diagnoses

broad_icd10 = long_icd %>%
                    mutate(ICD10_3 = str_sub(ICD10, 1, 3)) %>%
                    group_by(eid) %>%
                    summarise(
                        n_unique_icd10_3 = n_distinct(ICD10_3),
                        .groups = "drop"
                    )

colnames(healthy_individuals)[2] = "n_unique_icd10_3"
broad_icd10 = rbind(broad_icd10, healthy_individuals)

### Combine the diagnoses

all_icd10_counts = merge(full_counts, broad_icd10, by = "eid")

### Add in covariates

all_icd10_counts = left_join(all_icd10_counts, covariates, by = c("eid" = "IID"))
all_icd10_counts = all_icd10_counts[!(is.na(all_icd10_counts$FID)),]
all_icd10_counts$FID = NULL


#####################

### Run PRS analysis

#####################

### Covariates to include in the regression model

covars = setdiff(colnames(covariates), c("IID", "FID"))

### --- Main results table (one row per PRS) ---

outFrame = data.frame(matrix(ncol = 11, nrow = 0))
colnames(outFrame) = c("score", "median_diag", "median_diag_middle_decile", "median_diag_top_decile",
                       "median_diag3", "median_diag_middle_decile3", "median_diag_top_decile3",
                       "n", "beta", "se", "p")

### Decile summary table (10 rows per PRS)

decile_summary = data.frame(score = character(),
                                decile = integer(),
                                n = integer(),
                                mean_diag = numeric(),
                                median_diag = numeric(),
                                mean_diag3 = numeric(),
                                median_diag3 = numeric(),
                                stringsAsFactors = FALSE)

for (i in seq_along(allPGS)) {

    ### --- Read and standardise PRS ---

    currPgs = fread(allPGS[i], data.table = FALSE)

    currPgs[,2] = scale(currPgs[,2])

    currPRS = colnames(currPgs)[2]

    ### --- Merge PRS with phenotype data and remove missing diagnosis counts ---

    assoc_frame = left_join(currPgs, all_icd10_counts, by = c("IID" = "eid"))
    assoc_frame = assoc_frame[!is.na(assoc_frame$n_unique_icd10), ]

    ### Assign individuals to PRS deciles

    assoc_frame$prs_ten = cut(
        assoc_frame[, currPRS],
        breaks = quantile(
            assoc_frame[, currPRS],
            prob = 0:10/10,
            names = FALSE
        ),
        labels = 1:10,
        include.lowest = TRUE
    )

    assoc_frame$prs_ten = as.integer(as.character(assoc_frame$prs_ten))

    ### --- Calculate diagnosis statistics for each PRS decile ---

    decile_stats = aggregate(
        cbind(n_unique_icd10, n_unique_icd10_3) ~ prs_ten,
        data = assoc_frame,
        FUN = function(x) c(mean = mean(x), median = median(x), n = length(x))
    )

    decile_out = data.frame(
        score = currPRS,
        decile = decile_stats$prs_ten,
        n = decile_stats$n_unique_icd10[, "n"],
        mean_diag = decile_stats$n_unique_icd10[, "mean"],
        median_diag = decile_stats$n_unique_icd10[, "median"],
        mean_diag3 = decile_stats$n_unique_icd10_3[, "mean"],
        median_diag3 = decile_stats$n_unique_icd10_3[, "median"]
    )

    decile_summary = rbind(decile_summary, decile_out)

    ### --- Calculate summary statistics for the whole cohort and selected deciles ---

    top_percentile = assoc_frame[assoc_frame$prs_ten == 10, ]
    middle_percentile = assoc_frame[assoc_frame$prs_ten %in% c(4, 5, 6), ]

    top_percentile_diagn = median(top_percentile$n_unique_icd10)
    middle_percentile_diagn = median(middle_percentile$n_unique_icd10)
    median_diagn = median(assoc_frame$n_unique_icd10)

    top_percentile_diag3 = median(top_percentile$n_unique_icd10_3)
    middle_percentile_diag3 = median(middle_percentile$n_unique_icd10_3)
    median_diag3 = median(assoc_frame$n_unique_icd10_3)

    ### --- Linear regression of diagnosis count against PRS ---

    assoc_formula = as.formula(
        paste0(
            "n_unique_icd10 ~ ",
            currPRS,
            " + ",
            paste(covars, collapse = " + ")
        )
    )

    reg_model = glm(
        formula = assoc_formula,
        data = assoc_frame
    )

    reg_coeff = summary(reg_model)$coefficients

    ### --- Store regression results ---

    outRow = data.frame(
        score = currPRS,
        median_diag = median_diagn,
        median_diag_middle_decile = middle_percentile_diagn,
        median_diag_top_decile = top_percentile_diagn,
        median_diag3 = median_diag3,
        median_diag3_middle_decile = middle_percentile_diag3,
        median_diag3_top_decile = top_percentile_diag3,
        n = nrow(assoc_frame),
        beta = reg_coeff[2,1],
        se = reg_coeff[2,2],
        p = reg_coeff[2,4]
    )

    outFrame = rbind(outFrame, outRow)

}

### --- Colours for each phenotype ---

pheno_colours = c("IN" = "#0077BB",
                    "DRUG" = "#EE7733",
                    "PRIM" = "#009988",
                    "INOUT" = "#EE3377")

### --- Prepare plotting dataframe ---

plot_df = decile_summary %>% mutate(
                                phenotype = case_when(
                                        str_detect(score, "^IN_")   ~ "IN",
                                        str_detect(score, "^DRUG_") ~ "DRUG",
                                        str_detect(score, "^PRIM_") ~ "PRIM",
                                        str_detect(score, "^INOUT_") ~ "INOUT",
                                        TRUE ~ NA_character_
                                        )
                                ) %>%
                                filter(!is.na(phenotype))

# -------------------------------------------------------------------------
# Plot mean number of diagnoses by PRS decile
# -------------------------------------------------------------------------

deciles_p = ggplot(plot_df,
        aes(x = decile,
            y = mean_diag,
            colour = phenotype,
            group = phenotype)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    geom_point(size = 3) +
    scale_colour_manual(values = pheno_colours) +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "PRS (excluding UK Biobank) decile",
        y = "Mean number of ICD-10 diagnoses\nin UK Biobank",
        colour = "Phenotype"
    ) +
    theme_linedraw() +
    theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold")
    )

ggsave(filename = paste0(figureDir, "icd10_count_pgs_plot.png"),
       deciles_p, width = 10, height = 10, dpi = 300)



