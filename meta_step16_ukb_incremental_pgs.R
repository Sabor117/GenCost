### ---------------------------
###
### Script name: meta_step16_ukb_incremental_pgs.R
###
### Purpose of script: Test GenCost PGS against incremental addition of UKB traits
###
### Author: Dr. Sebastian May-Wilson
### Contact: sebastian.may-wilson@helsinki.fi
###
### Date Created: 2026-08-10
###
### ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(data.table)
library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

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
ukbDir = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/"
figureDir = paste0(mainDir, "outputs/figures/misc/")

allPGS = Sys.glob(paste0(pgsDir, "*noUKB*BAYESR*allchrs.tsv"))
allPGS = c(allPGS, paste0(pgsDir, "INOUT_ALL_meta_V4_hapmap_plus_BAYESR_megaprs_allchrs.tsv"))

ukb_base_data_file = paste0(ukbDir, "ukbb_78537_base_phenotypes.csv")
ukb_blood_biochem_file = paste0(ukbDir, "ukbb_78537_blood_biochem_phenotypes.tsv")

### Cost files

cost_covariatesFile = paste0(mainDir, "outputs/cost_phenotypes/gwas_cost_covariates.tsv")
cost_covariates = fread(cost_covariatesFile, data.table = FALSE)

in_costsFile = paste0(mainDir, "outputs/cost_phenotypes/FINAL_hes_based_phenotype_frame.tsv")
drug_costsFile = paste0(mainDir, "outputs/cost_phenotypes/FINAL_prescription_based_phenotype_frame.tsv")
prim_costsFile = paste0(mainDir, "outputs/cost_phenotypes/FINAL_gp_based_phenotype_frame.tsv")

in_cost = fread(in_costsFile, data.table = FALSE)
drug_cost = fread(drug_costsFile, data.table = FALSE)
prim_cost = fread(prim_costsFile, data.table = FALSE)

### R2 functions

get_r2_glm = function(model) { 1 - (model$deviance / model$null.deviance) }

get_adj_r2_glm = function(model) {
                    n = nobs(model)
                    p = length(coef(model)) - 1  # subtract intercept
                    r2 = 1 - (model$deviance / model$null.deviance)
                    1 - (1 - r2) * (n - 1) / (n - p - 1)
                    }


#####################

### Obtain covariates

#####################

### Main covariates = sex + age

ukb_primary_covariates = c("sex", # sex
                            "age_at_end_of_followup" # age
                            )

### "Phase 1" covariates are those from Jiwoo's paper

ukb_phase1_covariates = c("30720-0.0" = "cystatin_c",
                            "30760-0.0" = "hdl_cholesterol",
                            "48-0.0" = "waist_circumference",
                            "4080-0.0" = "systolic_blood_pressure",
                            "21001-0.0" = "bmi",
                            "30870-0.0" = "triglycerides"
                            )

### "Phase 2" covariates is smoking status

ukb_phase2_covariates = c("20116-0.0" = "smoking_status"
                            )

ukb_covariates = c(ukb_phase1_covariates, ukb_phase2_covariates)

### Read the files

ukb_base_data = fread(ukb_base_data_file, data.table = FALSE, select = c("eid", names(ukb_covariates)))
ukb_blood_biochem = fread(ukb_blood_biochem_file, data.table = FALSE, select = c("eid", names(ukb_covariates)))

### Change the names

ukb_base_data = setNames(ukb_base_data, ifelse(names(ukb_base_data) == "eid",
                                                "IID",
                                                ukb_covariates[names(ukb_base_data)])
                                                )

ukb_blood_biochem = setNames(ukb_blood_biochem, ifelse(names(ukb_blood_biochem) == "eid",
                                                        "IID",
                                                        ukb_covariates[names(ukb_blood_biochem)])
                                                        )

### Create a data frame with the various final costs

cost_frame = merge(in_cost[,c("eid", "log_hes_costs_year")],
                    drug_cost[,c("eid", "log_prescription_costs_year")],
                    by = "eid",
                    all = TRUE)
cost_frame = merge(cost_frame,
                    prim_cost[,c("eid", "log_gp_costs_year")],
                    by = "eid",
                    all = TRUE)

colnames(cost_frame) = c("IID", "IN", "DRUG", "PRIM")


#####################

### Run incremental analysis

#####################

### Covariates to include in the regression model

m0_covars = setdiff(colnames(cost_covariates), c("IID", "FID"))
m1_covars = c(m0_covars, unname(ukb_phase1_covariates))
m2_covars = c(m1_covars, unname(ukb_phase2_covariates))


### --- Main results table ---

colset = c("score", "model", "r2", "adj_r2",
                       "incremental_r2_over_no_pgs", "beta", "se",
                       "p")

outFrame = data.frame(matrix(ncol = length(colset), nrow = 0))
colnames(outFrame) = colset

for (i in seq_along(allPGS)) {

    ### --- Read and standardise PRS ---

    currPgs = fread(allPGS[i], data.table = FALSE)

    currPgs[,2] = scale(currPgs[,2])

    currPRS_col = colnames(currPgs)[2]

    currPRS = case_when(grepl("^DRUG_", currPRS_col) ~ "DRUG",
                            grepl("^IN_",   currPRS_col) ~ "IN",
                            grepl("^INOUT_", currPRS_col) ~ "INOUT",
                            grepl("^PRIM_", currPRS_col) ~ "PRIM",
                            TRUE ~ NA_character_)

    if (currPRS == "INOUT") {

        cat(("\nInpatient + outpatient not present. Skipping.\n...\n\n"))

        next

    }
    
    ### --- Merge PRS with phenotype data and covariates ---

    assoc_frame = left_join(currPgs, cost_covariates[,-1], by = "IID")
    assoc_frame = left_join(assoc_frame, ukb_base_data, by = "IID")
    assoc_frame = left_join(assoc_frame, ukb_blood_biochem, by = "IID")
    assoc_frame = merge(assoc_frame, cost_frame[,c("IID", currPRS)], by = "IID", all = TRUE)

    ### --- Linear regression of models of incremental covariates ---

    m0_assoc_formula = as.formula(
        paste0(currPRS, " ~ ",
            paste(m0_covars, collapse = " + ")
        )
    )

    m1_assoc_formula = as.formula(
        paste0(currPRS, " ~ ",
            paste(m1_covars, collapse = " + ")
        )
    )

    m2_assoc_formula = as.formula(
        paste0(currPRS, " ~ ",
            paste(m2_covars, collapse = " + ")
        )
    )

    ### --- Linear regression of models of incremental covariates PLUS PRS ---

    m0_assoc_formula_prs = as.formula(
        paste0(currPRS, " ~ ",
            currPRS_col, " + ",
            paste(m0_covars, collapse = " + ")
        )
    )

    m1_assoc_formula_prs = as.formula(
        paste0(currPRS, " ~ ",
            currPRS_col, " + ",
            paste(m1_covars, collapse = " + ")
        )
    )

    m2_assoc_formula_prs = as.formula(
        paste0(currPRS, " ~ ",
            currPRS_col, " + ",
            paste(m2_covars, collapse = " + ")
        )
    )

    ### --- Run all models ---
    
    m0_model = glm(formula = m0_assoc_formula, data = assoc_frame)
    m1_model = glm(formula = m1_assoc_formula, data = assoc_frame)
    m2_model = glm(formula = m2_assoc_formula, data = assoc_frame)
    m0_model_prs = glm(formula = m0_assoc_formula_prs, data = assoc_frame)
    m1_model_prs = glm(formula = m1_assoc_formula_prs, data = assoc_frame)
    m2_model_prs = glm(formula = m2_assoc_formula_prs, data = assoc_frame)

    ### --- Get coefficients --- 

    m0_coeff = summary(m0_model)$coefficients
    m1_coeff = summary(m1_model)$coefficients
    m2_coeff = summary(m2_model)$coefficients
    m0_prs_coeff = summary(m0_model_prs)$coefficients
    m1_prs_coeff = summary(m1_model_prs)$coefficients
    m2_prs_coeff = summary(m2_model_prs)$coefficients

    ### --- R2 extraction ---

    r2_table = tibble(model = c("M0", "M1", "M2", "M0_PGS", "M1_PGS", "M2_PGS"),
                        r2 = c(get_r2_glm(m0_model),
                                get_r2_glm(m1_model),
                                get_r2_glm(m2_model),
                                get_r2_glm(m0_model_prs),
                                get_r2_glm(m1_model_prs),
                                get_r2_glm(m2_model_prs)),
                        adj_r2 = c(get_adj_r2_glm(m0_model),
                                    get_adj_r2_glm(m1_model),
                                    get_adj_r2_glm(m2_model),
                                    get_adj_r2_glm(m0_model_prs),
                                    get_adj_r2_glm(m1_model_prs),
                                    get_adj_r2_glm(m2_model_prs))
                        )

    ### --- Incremental R² (PGS contribution at each phase) ---

    r2_table = r2_table %>% mutate(incremental_r2_over_no_pgs = case_when(
                                                                    model == "M0_PGS" ~ r2 - get_r2_glm(m0_model),
                                                                    model == "M1_PGS" ~ r2 - get_r2_glm(m1_model),
                                                                    model == "M2_PGS" ~ r2 - get_r2_glm(m2_model),
                                                                    TRUE ~ NA_real_
                                                                    )
                                    )

    print(r2_table)

    ### --- PGS coefficient extraction ---
    ### Adjust "PGS" to match the exact variable name in your models

    pgs_varname = currPRS_col

    get_pgs_row = function(coeff_matrix, model_label) {

    if (!(pgs_varname %in% rownames(coeff_matrix))) {

        return(tibble(model = model_label, estimate = NA, se = NA, t = NA, p = NA))

    }

    row = coeff_matrix[pgs_varname, ]
    tibble(model = model_label,
            beta = row["Estimate"],
            se = row["Std. Error"],
            p = row["Pr(>|t|)"])

    }

    pgs_coeff_table = bind_rows(get_pgs_row(m0_prs_coeff, "M0_PGS"),
                                get_pgs_row(m1_prs_coeff, "M1_PGS"),
                                get_pgs_row(m2_prs_coeff, "M2_PGS"))

    print(pgs_coeff_table)

    ### --- Combined summary table ---
    
    summary_table = left_join(r2_table, pgs_coeff_table, by = "model")
    print(summary_table)

    summary_table$score = currPRS

    outFrame = rbind(outFrame, summary_table)

}

outFrame = outFrame[, colset]

fwrite(outFrame, paste0(mainDir, "outputs/misc_results/incremental_pgs_ukb_res.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")


#####################

### Plotting

#####################

### --- Colours for each phenotype ---

pheno_colours = c(
    "IN"   = "#0077BB",
    "DRUG" = "#EE7733",
    "PRIM" = "#009988",
    "INOUT" = "#EE3377"
)

### --- Data for total R2 ---
r2_plot = outFrame %>%
    filter(model %in% c("M0", "M1", "M2", "M0_PGS", "M1_PGS", "M2_PGS")) %>%
    mutate(
        model_base = factor(
            case_when(
                model %in% c("M0", "M0_PGS") ~ "M0",
                model %in% c("M1", "M1_PGS") ~ "M1",
                model %in% c("M2", "M2_PGS") ~ "M2"
            ),
            levels = c("M0", "M1", "M2")
        ),
        PGS = grepl("_PGS$", model)
    )

### --- Total R2 ---

p_r2 = ggplot(r2_plot,
        aes(x = model_base,
            y = r2,
            colour = score,
            group = score)) +

    geom_line(
        data = filter(r2_plot, !PGS),
        linewidth = 1
    ) +

    geom_point(
        data = filter(r2_plot, !PGS),
        size = 3
    ) +

    geom_point(
        data = filter(r2_plot, PGS),
        size = 3,
        shape = 17
    ) +

    geom_line(
        data = filter(r2_plot, PGS),
        linewidth = 1,
        linetype = "dashed"
    ) +

    scale_colour_manual(values = pheno_colours) +

    labs(x = "Covariate model",
        y = expression(R^2),
        colour = "Phenotype",
        title = "Variance explained by covariate and PGS models"
    ) +

    theme_linedraw() +
    theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold")
    )


pgs_plot = outFrame %>%
    filter(grepl("_PGS$", model)) %>%
    mutate(
        model_base = factor(
            sub("_PGS$", "", model),
            levels = c("M0", "M1", "M2")
        )
    )

p_incremental = ggplot(pgs_plot,
        aes(x = model_base,
            y = incremental_r2_over_no_pgs,
            colour = score)) +

    geom_line(
        aes(group = score),
        linewidth = 1
    ) +

    geom_point(size = 3) +

    scale_colour_manual(values = pheno_colours) +

    labs(x = "Covariate model",
        y = expression("Incremental " * R^2),
        colour = "Phenotype",
        title = "Additional variance explained by PGS"
    ) +

    theme_linedraw() +

    theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold")
    )

out_plot = p_r2 / p_incremental

ggsave(filename = paste0(figureDir, "incremental_r2_pgs_plot.png"),
       out_plot, width = 10, height = 10, dpi = 300)

