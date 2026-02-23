##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
#library(ggplot2)
library(tidyverse)
library(gt, lib.loc = "/projappl/project_2007428/RPackages_421/")
#library(patchwork)
library(forestplot, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(dplyr)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
snpDir = paste0(mainDir, "outputs/snp_forest_inputs_v4/")
figureDir = paste0(mainDir, "outputs/figures/snp_forests/")

all_files = Sys.glob(paste0(snpDir, "*_top_snp_metas.txt"))

file_select = 12


##### =========================== #####

### Data processing

##### =========================== #####

snp_meta_frame = fread(all_files[file_select], data.table = FALSE)

for (nsnp in 1:length(unique(snp_meta_frame$snpid))){

    #nsnp = nsnp + 1 # for running interactively

    print(nsnp)

    ### Select SNP and get table

    curr_snp = unique(snp_meta_frame$snpid)[nsnp]

    curr_snp_frame = snp_meta_frame[snp_meta_frame$snpid == curr_snp,]

    if (length(unique(curr_snp_frame$snpid)) > 1){

        curr_snp_frame = curr_snp_frame[curr_snp_frame$snpid == curr_snp_frame$snpid[curr_snp_frame$cohort == "Meta"],]

    }

    ### Create upper and lower "CI" bounds

    curr_snp_frame$conf_low = curr_snp_frame$beta1 - curr_snp_frame$se
    curr_snp_frame$conf_high = curr_snp_frame$beta1 + curr_snp_frame$se

    plot_min = min(curr_snp_frame$conf_low)
    plot_max = max(curr_snp_frame$conf_high)

    ### Get plot edges

    rounded_min = floor(plot_min / 0.05) * 0.05
    rounded_max = ceiling(plot_max / 0.05) * 0.05

    ### make plotting beta

    curr_snp_frame$beta_plot = curr_snp_frame$beta1

    ### Create labels in table

    curr_snp_frame = curr_snp_frame |>
    # round estimates and 95% CIs to 2 decimal places for journal specifications
    mutate(across(
        c(beta_plot, conf_low, conf_high),
        ~ str_pad(
        round(.x, 3),
        width = 4,
        pad = "0",
        side = "right"
        )
    ),
    # add an "-" between HR estimate confidence intervals
    estimate_lab = paste0(beta_plot, " (", conf_low, "-", conf_high, ")")) |>
    # round p-values to two decimal places, except in cases where p < .001
    mutate(p_lab = case_when(
        p < .001 ~ "<0.001",
        round(p, 2) == .05 ~ as.character(round(p, 3)),
        p < .01 ~ str_pad( # if less than .01, go one more decimal place
        as.character(round(p, 3)),
        width = 4,
        pad = "0",
        side = "right"
        ),
        TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
        as.character(round(p, 2)),
        width = 4,
        pad = "0",
        side = "right"
        )
    ))

    ### Remove the "Meta" row for the plotting table 

    forest_plot_frame = curr_snp_frame[curr_snp_frame$cohort != "Meta",]

    forest_plot_frame = tibble::tibble(mean  = as.numeric(forest_plot_frame$beta_plot),
                                    lower = as.numeric(forest_plot_frame$conf_low),
                                    upper = as.numeric(forest_plot_frame$conf_high),
                                    cohort = forest_plot_frame$cohort,
                                    estimate_lab = forest_plot_frame$estimate_lab,
                                    p_val = forest_plot_frame$p_lab)

    ### Extract snpid from table and then extract this from Meta

    analysis_name = gsub("_top_snp_metas.txt", "", basename(all_files[file_select]))
    curr_snpid = unique(curr_snp_frame$snpid)

    select_meta_file = paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/", analysis_name, "_metal_output_1.TBL.gz")

    meta_snp_row = system(paste0("zcat ", select_meta_file, " | grep ", curr_snpid), intern = TRUE)

    ### Split METAL row into 15 columns

    meta_snp_row = as.data.frame(str_split_fixed(meta_snp_row, "\t", 17))

    if (nrow(meta_snp_row) > 1){

        meta_snp_row = meta_snp_row[meta_snp_row[,1] == curr_snpid,]

    }

    ### I^2 stat == col 12
    ### I^2 P-val == col 15

    i_sq_stat = as.numeric(meta_snp_row[,12])
    i_sq_pval = as.numeric(meta_snp_row[,15])

    i_sq_pval = round(i_sq_pval, digits = 2)

    ### Create label from values

    meta_label = paste0("Meta\nI^2 = ", i_sq_stat, " (P = ", i_sq_pval,")")

    ### Make outdir

    outDir = paste0(figureDir, analysis_name, "/")
    system(paste0("mkdir -p ", outDir))

    if (length(unique(curr_snp_frame$a1[curr_snp_frame$cohort != "Meta"])) > 1){

        cat(paste0("PROBLEM SNP: ", curr_snp, ". Has multiple alleles: ", 
                paste(unique(curr_snp_frame$a1[curr_snp_frame$cohort != "Meta"]), collapse = ", "), ". Position: ", unique(curr_snp_frame$snpid), "\n\n\n"))

    }


    ### Get meta-values for plot

    if (all(curr_snp_frame$a1[curr_snp_frame$cohort == "Meta"] != unique(curr_snp_frame$a1[curr_snp_frame$cohort != "Meta"]))){

        curr_snp_frame$beta1[curr_snp_frame$cohort == "Meta"] = as.numeric(curr_snp_frame$beta1[curr_snp_frame$cohort == "Meta"]) * -1

    }

    meta_beta = as.numeric(curr_snp_frame$beta1[curr_snp_frame$cohort == "Meta"])
    meta_se = as.numeric(curr_snp_frame$se[curr_snp_frame$cohort == "Meta"])
    meta_lower = meta_beta - meta_se
    meta_higher = meta_beta + meta_se
    meta_pval =  as.numeric(curr_snp_frame$p[curr_snp_frame$cohort == "Meta"])

    meta_pval_lab = sprintf("%.1e", meta_pval)
    meta_est_lab = paste0(round(meta_beta, 3), "  (", round(meta_lower, 3), "-", round(meta_higher, 3), ")")

    plot_title = paste0(curr_snp, ", N = ", curr_snp_frame$n[curr_snp_frame$cohort == "Meta"])

    if (length(unique(curr_snp_frame$a1[curr_snp_frame$cohort != "Meta"])) > 1){

        plot_title = paste0(plot_title, " PROBLEM SNP")

    }


    ##### =========================== #####

    ### Make forest plot

    ##### =========================== #####

    fig_name = paste0(analysis_name, "_", curr_snpid)

    png(paste0(outDir, fig_name, ".png"), width = 3000, height = 2000, res = 300)

    p = forest_plot_frame |> 
        forestplot(labeltext = c(cohort, estimate_lab, p_val), 
                    clip = c(rounded_min, rounded_max),
                    vertices = TRUE,
                    title = plot_title) |> 
        fp_add_lines() |> 
        fp_set_style(box = "royalblue",
                    line = "darkblue",
                    summary = "darkred",
                    align = "lrrr",
                    hrz_lines = "#999999") |> 
        fp_add_header(cohort = c("", "Cohort"),
                        estimate_lab = c("", fp_align_center("Beta (SE)")),
                        p_val = c("", "P-value")) |>
        fp_append_row(mean  = meta_beta,
                        lower = meta_lower,
                        upper = meta_higher,
                        cohort = meta_label,
                        estimate_lab = meta_est_lab,
                        p_val = meta_pval_lab,
                        is.summary = TRUE)

    print(p)

    dev.off()

}
