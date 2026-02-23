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
prsDir = paste0(mainDir, "outputs/prs/scores/")
assocDir = paste0(mainDir, "outputs/prs/association/")
baseDir = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/"


##### =========================== #####

### Data setup

##### =========================== #####

### First create full PRS from all chorosome files
### If it already exists, read it instead

all_chromosome_files = Sys.glob(paste0(prsDir, "*.sscore"))

unique_sets = unique(gsub("_chr[0-9]+\\.sscore", "", basename(all_chromosome_files)))

for (i in 1:length(unique_sets)){

    currPRS = unique_sets[i]

    if (file.exists(paste0(prsDir, currPRS, "_allchrs.tsv"))){

        cat(paste0("\n", prsDir, currPRS, "_allchrs.tsv exists. Can be read.\n\n...\n\n"))

        #prs_out = fread(paste0(prsDir, currPRS, "_allchrs.tsv"), data.table = FALSE)

    } else {

        cat(paste0("\n", prsDir, currPRS, "_allchrs.tsv doesn't exist. Making.\n\n...\n\n"))

        allPRS_files = Sys.glob(paste0(prsDir, currPRS, "*.sscore"))

        prs_out = fread(allPRS_files[1], data.table = FALSE)

        ### Use AVG score

        prs_out = prs_out[,c("IID", "SCORE1_AVG")]
        colnames(prs_out) = c("IID", "avg_score_1")

        for (i in 2:length(allPRS_files)){

            currPrs_chr = fread(allPRS_files[i], data.table = FALSE)

            print(basename(allPRS_files[i]))
            print(head(currPrs_chr))

            currPrs_chr = currPrs_chr[,c("IID", "SCORE1_AVG")]
            colnames(currPrs_chr)[2] = paste0("avg_score_", i)

            prs_out = merge(prs_out, currPrs_chr, by = "IID")

        }

        ### Sum of each choromosome score

        prs_out$PRS_out = rowSums(prs_out[,2:ncol(prs_out)])

        prs_out = prs_out[,c("IID", "PRS_out")]
        colnames(prs_out)[2] = paste0(currPRS)

        ### Write full PRS for phenotype

        fwrite(prs_out, paste0(prsDir, currPRS, "_allchrs.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

    }
}

### Select the phenotype for analysis

cost_phenos = c("log_hes_costs_year", "log_gp_costs_year", "log_prescription_costs_year")

### Read cost phenotypes and covariates for association testing

assoc_phenos = fread(paste0(phenoDir, "FINAL_hes_based_phenotype_frame.tsv"), data.table = FALSE)
assoc_phenos = assoc_phenos[,c("IID", cost_phenos)]

covars = fread(paste0(phenoDir, "gwas_cost_covariates.tsv"), data.table = FALSE)

### Restrict to Europeans

ancest_restrict = FALSE

if (isTRUE(ancest_restrict)){

    covars = covars[covars$SuperPop == 1,]

}

covariate_names = c("age_at_end_of_followup", "batch", "array", "sex", paste0("pc", 1:20))

covars = covars[,c("IID", covariate_names)]

### Other PCs

pc_fields = paste0("22009-0.", c(21:40))
additional_pcs = fread(paste0(baseDir, "ukbb_78537_base_phenotypes.csv"), data.table = FALSE, select = c("eid", pc_fields, "21001-0.0"))

colnames(additional_pcs) = c("IID", paste0("pc", 21:40), "BMI")

covars = merge(covars, additional_pcs, by = "IID")

covariate_names = c("age_at_end_of_followup", "batch", "array", "sex", paste0("pc", 1:40))


##### =========================== #####

### Run analysis

##### =========================== #####

prs_files = Sys.glob(paste0(prsDir, "*allchrs.tsv"))

prs_comparison = data.frame(matrix(ncol = 16, nrow = 0))
colnames(prs_comparison) = c("trait",
                "score",
                "median_cost",
                "median_top_percent",
                "median_middle_percent",
                "median_top_percent_bmi",
                "median_middle_percent_bmi",
                "n",
                "beta_cost",
                "se_cost",
                "beta_bmi",
                "se_bmi",
                "p_cost",
                "p_bmi",
                "r2_cost",
                "r2_bmi")

prs_percentiles = data.frame(matrix(ncol = 4, nrow = 0))
colnames(prs_percentiles) = c("percentile", "median_cost", "PGS", "cost")

prs_deciles = data.frame(matrix(ncol = 4, nrow = 0))
colnames(prs_deciles) = c("decile", "median_cost", "PGS", "cost")

for (ncost in 1:length(cost_phenos)){

    ### Select the current cost phenotype

    cost_pheno = cost_phenos[ncost]

    ### Select the correct phenotype in data frame

    curr_phenos = assoc_phenos[,c("IID", cost_pheno)]

    for (nprs in 1:length(prs_files)){

        prs_out = fread(prs_files[nprs], data.table = FALSE)
        currPRS = gsub("_allchrs.tsv", "", basename(prs_files[nprs]))

        colnames(prs_out)[2] = currPRS

        ### Making the analysis data frame

        cat(paste0("\nNumber of PRS individuals: ", nrow(prs_out), "\n_____\n"))

        assoc_frame = merge(prs_out, curr_phenos, by = "IID")

        cat(paste0("\nNumber of individuals with PRS + costs: ", nrow(na.omit(assoc_frame)), "\n_____\n"))

        assoc_frame = merge(assoc_frame, covars, by = "IID")

        cat(paste0("\nNumber of individuals with PRS + costs + covariates: ", nrow(na.omit(assoc_frame)), "\n_____\n"))

        ### Making the analysis formulae

        covariates = paste(covariate_names, collapse = " + ")

        ### Removing NAs from analysis

        assoc_frame = na.omit(assoc_frame)

        ### Starting to run analysis

        median_cost = median(assoc_frame[,cost_pheno])

        ### Get the 10 quantiles of PRS

        assoc_frame$prs_ten = ntile(assoc_frame[[currPRS]], 10)

        ### Get the 100 percentiles of PRS

        assoc_frame$prs_percentile = ntile(assoc_frame[[currPRS]], 100)

        ### Get the top percentile of cost
        ### And the middle percentiles

        top_percentile = assoc_frame[assoc_frame$prs_ten == 10,]
        top_percentile_median = median(top_percentile[,cost_pheno])
        top_percentile_median_bmi = median(na.omit(top_percentile[,"BMI"]))

        middle_percentile = assoc_frame[assoc_frame$prs_ten %in% c(4, 5, 6),]
        middle_percentile_median = median(middle_percentile[,cost_pheno])
        middle_percentile_median_bmi = median(na.omit(middle_percentile[,"BMI"]))

        ### Create the association formula

        assoc_formula = as.formula(paste0(cost_pheno, " ~ ", currPRS, "_scaled + ", covariates))
        bmi_formula = as.formula(paste0("BMI ~ ", currPRS, "_scaled + ", covariates))

        ### For getting the R2 of the PGS - additional formulae

        covar_only_formula = as.formula(paste0(cost_pheno, " ~ ", covariates))
        prs_only_formula = as.formula(paste0(cost_pheno, " ~ ", currPRS, "_scaled"))

        ### Scale PRS prior to association

        assoc_frame$scaled_PRS = scale(assoc_frame[,currPRS])
        rename_integer = which(colnames(assoc_frame) == "scaled_PRS")

        colnames(assoc_frame)[rename_integer] = paste0(currPRS, "_scaled")

        ### glm of formula
        ### Get summary of this and coefficients

        bmi_frame = assoc_frame[!(is.na(assoc_frame$BMI)),]

        ### Make models

        reg_model = glm(formula = assoc_formula, data = assoc_frame)
        bmi_model = glm(formula = bmi_formula, data = bmi_frame)
        covar_model = glm(formula = covar_only_formula, data = assoc_frame)
        prs_model   = glm(formula = prs_only_formula, data = assoc_frame)

        reg_summary = summary(reg_model)
        bmi_summary = summary(bmi_model)

        reg_coeff = reg_summary$coefficients
        bmi_coeff = bmi_summary$coefficients

        ### Multiple R2s

        reg_model_r2 = 1 - sum(residuals(reg_model)^2) / sum((assoc_frame[,cost_pheno] - mean(assoc_frame[,cost_pheno]))^2)
        bmi_model_r2 = 1 - sum(residuals(bmi_model)^2) / sum((bmi_frame[,"BMI"] - mean(bmi_frame[,"BMI"]))^2)
        covar_model_r2 = 1 - sum(residuals(covar_model)^2) / sum((assoc_frame[,cost_pheno] - mean(assoc_frame[,cost_pheno]))^2)
        prs_model_r2 = 1 - sum(residuals(prs_model)^2) / sum((assoc_frame[,cost_pheno] - mean(assoc_frame[,cost_pheno]))^2)

        sample_size = nrow(assoc_frame)

        ### Make the output

        outRow = data.frame(trait = cost_pheno,
                    score = currPRS,
                    median_cost = median_cost,
                    median_top_percent = top_percentile_median,
                    median_middle_percent = middle_percentile_median,
                    median_top_percent_bmi = top_percentile_median_bmi,
                    median_middle_percent_bmi = middle_percentile_median_bmi,
                    n = sample_size,
                    beta_cost = reg_coeff[paste0(currPRS, "_scaled"),1],
                    se_cost = reg_coeff[paste0(currPRS, "_scaled"),2],
                    beta_bmi = bmi_coeff[paste0(currPRS, "_scaled"),1],
                    se_bmi = bmi_coeff[paste0(currPRS, "_scaled"),2],
                    p_cost = reg_coeff[paste0(currPRS, "_scaled"),4],
                    p_bmi = bmi_coeff[paste0(currPRS, "_scaled"),4],
                    r2_model = reg_model_r2,
                    r2_covars = covar_model_r2,
                    r2_cost = prs_model_r2,
                    r2_cost_incremental = reg_model_r2 - covar_model_r2,
                    r2_bmi = bmi_model_r2)

        cat(paste0("Median log(yearly cost) for the top 10% of PRS is: ", top_percentile_median, ". This is equivalent to: €", exp(top_percentile_median), "\n\n"))
        cat(paste0("Median log(yearly cost) for the middle 40-60% of PRS is: ", middle_percentile_median, ". This is equivalent to: €", exp(middle_percentile_median), ".\n\n"))

        prs_comparison = rbind(prs_comparison, outRow)

        decile_table = assoc_frame %>%
                            group_by(decile = prs_ten) %>%
                                summarise(
                                median_cost = median(.data[[cost_pheno]], na.rm = TRUE),
                                PGS  = currPRS,
                                cost = cost_pheno,
                                .groups = "drop"
                                ) %>%
                            arrange(decile)

        prs_deciles = rbind(prs_deciles, decile_table)

        percentile_table = assoc_frame %>%
                            group_by(percentile = prs_percentile) %>%
                                summarise(
                                median_cost = median(.data[[cost_pheno]], na.rm = TRUE),
                                PGS  = currPRS,
                                cost = cost_pheno,
                                .groups = "drop"
                                ) %>%
                            arrange(percentile)

        prs_percentiles = rbind(prs_percentiles, percentile_table)

    }
}

fwrite(prs_comparison, paste0(assocDir, "gencost_UKB_prs_association_frame.txt"), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
fwrite(prs_percentiles, paste0(assocDir, "gencost_UKB_prs_percentiles_frame.txt"), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
fwrite(prs_deciles, paste0(assocDir, "gencost_UKB_prs_deciles_frame.txt"), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")





fwrite(ldpred_frame, "../../processing/misc_data/ldpred_prs_quintiles.txt", quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

assoc_frame$prs_ten = factor(assoc_frame$prs_ten, levels = unique(assoc_frame$prs_ten))

assoc_frame$cost_year = exp(assoc_frame$log_hes_costs_year)

# Create a combined plot showing distributions of `log_hes_costs_year` and `bmi`
ggplot(your_dataframe, aes(x = prs_ten, y = cost_year)) +  
    # Box plot for `log_hes_costs_year`
    geom_boxplot(aes(fill = "cost_year"), position = position_dodge(width = 0.8)) +
    # Optional: If you want to plot `bmi` in the same boxplot
    geom_boxplot(aes(y = bmi, fill = "BMI"), position = position_dodge(width = 0.8)) +

    # Labels and title
    labs(x = "PRS Quintiles", y = "Value", fill = "Variable", 
        title = "Distribution of Log HES Costs and BMI Across PRS Quintiles") +

    # Custom theme and adjustments
    theme_minimal() +
    theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14)) +

    # Set custom colors for the boxplots
    scale_fill_manual(values = c("log_hes_costs_year" = "#1b9e77", "BMI" = "#d95f02"))



