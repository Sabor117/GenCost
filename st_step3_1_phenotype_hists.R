#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("It's just a flesh wound...")


#####################

### File locations and script prep

#####################

final_costs_individual_plotting = fread("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/cost_phenotypes/cost_phenotypes.tsv", data.table = FALSE)

final_costs_individual_plotting$basic_cost_phenotype = log((final_costs_individual_plotting$total_cost / final_costs_individual_plotting$days_of_followup) * 365.25)
final_costs_individual_plotting$weighted_cost_phenotype = log((final_costs_individual_plotting$weighted_total_cost / final_costs_individual_plotting$days_of_followup) * 365.25)

final_costs_individual_plotting$basic_total_cost_phenotype = log((rowSums(final_costs_individual_plotting[,c("total_cost", "gp_costs")], na.rm = FALSE) / final_costs_individual_plotting$days_of_followup) * 365.25)
final_costs_individual_plotting$weighted_total_cost_phenotype = log((rowSums(final_costs_individual_plotting[,c("weighted_total_cost", "gp_costs")], na.rm = FALSE) / final_costs_individual_plotting$days_of_followup) * 365.25)

final_costs_individual_plotting$basic_cost_phenotype_year = (final_costs_individual_plotting$total_cost / final_costs_individual_plotting$days_of_followup) * 365.25
final_costs_individual_plotting$weighted_cost_phenotype_year = (final_costs_individual_plotting$weighted_total_cost / final_costs_individual_plotting$days_of_followup) * 365.25





#####################

### Making histograms

#####################

final_costs_individual_plotting = final_costs_individual_plotting[final_costs_individual_plotting$total_cost > 0 & final_costs_individual_plotting$weighted_total_cost > 0,]

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/basic_costs_histogram.pdf"),
    height = 8, width = 8)

ggplot(final_costs_individual_plotting, aes(x = total_cost)) +
  geom_histogram(binwidth = 1000, fill = "blue", color = "black") +
  labs(title = "Basic cost", x = "Cost", y = "Frequency") +
  theme_classic()

dev.off()

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/basic_cost_phenotype_histogram.pdf"),
    height = 8, width = 8)

ggplot(final_costs_individual_plotting, aes(x = basic_cost_phenotype)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +
  labs(title = "Basic cost phenotype", x = "Log cost", y = "Frequency") +
  theme_classic()

dev.off()

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/weighted_costs_histogram.pdf"),
    height = 8, width = 8)

ggplot(final_costs_individual_plotting, aes(x = weighted_total_cost)) +
  geom_histogram(binwidth = 1000, fill = "red", color = "black") +
  labs(title = "Weighted cost", x = "Weighted cost", y = "Frequency") +
  theme_classic()

dev.off()

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/weighted_cost_phenotype_histogram.pdf"),
    height = 8, width = 8)

ggplot(final_costs_individual_plotting, aes(x = weighted_cost_phenotype)) +
  geom_histogram(binwidth = 0.5, fill = "red", color = "black") +
  labs(title = "Weighted cost phenotype", x = "Log weighted cost", y = "Frequency") +
  theme_classic()

dev.off()

final_costs_individual_plotting$gp_costs[is.na(final_costs_individual_plotting$gp_costs)] = 0

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/gp_costs_histogram.pdf"),
    height = 8, width = 8)

ggplot(final_costs_individual_plotting, aes(x = gp_costs)) +
  geom_histogram(binwidth = 500, fill = "green", color = "black") +
  labs(title = "GP cost", x = "Cost", y = "Frequency") +
  theme_classic()

dev.off()

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/outputs/figures/overlaid_total_cost_phenotype_histogram.pdf"),
    height = 8, width = 8)

ggplot() +
  geom_histogram(data = final_costs_individual_plotting, aes(x = basic_total_cost_phenotype), fill = "blue", binwidth = 0.5, alpha = 0.5) +
  geom_histogram(data = final_costs_individual_plotting, aes(x = weighted_total_cost_phenotype), fill = "red", binwidth = 0.5, alpha = 0.5) +
  labs(title = "Overlaid total costs", x = "Log cost", y = "Frequency") +
  theme_classic()

dev.off()


