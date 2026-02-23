library(corrplot)
library(data.table)
library(reshape2)
library(reshape)
library(ggplot2)
library(tidyverse)

phenotype = "META_v4"
n_cohorts = length(Sys.glob(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/ldsc_intermediate_files/", phenotype, "/*_ldsc_input.txt.gz"))) - 1

all_corrs = Sys.glob(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/ldsc_intermediate_files/", phenotype, "/*ldsc_corr.log"))

corr_file = readLines(all_corrs[1])

start_line = grep("p1", corr_file)

outFile = fread(all_corrs[1], data.table = FALSE, skip = (start_line - 1), nrow = n_cohorts)

for (i in 2:length(all_corrs)){

	corr_file = readLines(all_corrs[i])

	start_line = grep("p1", corr_file)

	curr_corr = fread(all_corrs[i], data.table = FALSE, skip = (start_line - 1), nrow = n_cohorts)

	outFile = rbind(outFile, curr_corr)

}

outFile$cohort1 = basename(outFile$p1)
outFile$cohort2 = basename(outFile$p2)

outFile$cohort1 = gsub(paste0("_", phenotype, "_ldsc_munged.sumstats.gz"), "", outFile$cohort1)
outFile$cohort2 = gsub(paste0("_", phenotype, "_ldsc_munged.sumstats.gz"), "", outFile$cohort2)

outFile$cohort1 = gsub("_ldsc_input.txt.gz", "", outFile$cohort1)
outFile$cohort2 = gsub("_ldsc_input.txt.gz", "", outFile$cohort2)

### Create custom null hypothesis p-value (for null hypothesis = 1)

outFile$Pval = pchisq((abs(outFile$rg) - 1 / outFile$se) ^ 2,
							df = 1,
							lower = F)

colnames(outFile)[which(colnames(outFile) == "Pval")] = "p_h1"
colnames(outFile)[which(colnames(outFile) == "p")] = "p_h0"

# Function to remove duplicates with opposite cohorts
remove_opposite_duplicates = function(data) {

	### Create a new column with sorted cohort names

	data$sorted_cohorts = apply(data[, c("cohort1", "cohort2")], 1, function(x) paste(sort(x), collapse = " "))

	### Remove duplicates with opposite cohorts

	unique_data = data[!duplicated(data$sorted_cohorts), ]

	### Remove the temporary column

	unique_data$sorted_cohorts = NULL

	return(unique_data)
}

filtered_outFile = remove_opposite_duplicates(outFile)

fwrite(filtered_outFile, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/ldsc/cohort_ldsc_correlations_", phenotype, "_v04.txt"),
		quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
fwrite(outFile, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/ldsc/cohort_ldsc_correlations_nonfiltered_", phenotype, "_v04.txt"),
		quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")

outFile_rgs = outFile

### Removes cohorts manually is they are all NA

#outFile_rgs = outFile[!(outFile$cohort1 == "NTR" | outFile$cohort1 == "CHB" | outFile$cohort1 == "QGP"),]
#outFile_rgs = outFile_rgs[!(outFile_rgs$cohort2 == "NTR" | outFile_rgs$cohort2 == "CHB" | outFile_rgs$cohort2 == "QGP"),]

### Create label column for rg

outFile_rgs$rg_lab = ifelse(is.na(outFile_rgs$rg), "NA", ifelse(outFile_rgs$rg > 1, ">1", round(outFile_rgs$rg, digits = 2)))
outFile_rgs$label = paste0("rg: ", outFile_rgs$rg_lab, "\np = ", outFile_rgs$p_h0)

### Creates long data-frame

filtered_rgMat = melt(outFile_rgs[,c("cohort1", "cohort2", "rg")])
filtered_rgLabel = melt(outFile_rgs[,c("cohort1", "cohort2", "label")])
filtered_pMat = melt(outFile_rgs[,c("cohort1", "cohort2", "p_h0")])

### Colour gradient from white to a deep blue in gradient_int steps (gradient_int == 10)

gradient_int = 10
gradient_colors = colorRampPalette(c("white", "#19579E"))(gradient_int)

filtered_rgMat$value[filtered_rgMat$value > 1] = 1

filtered_rgMat$colour_code = gradient_colors[findInterval(filtered_rgMat$value, seq(0, 1, length.out = gradient_int))]

filtered_rgMat$colour_code[outFile_rgs$p < 0.05] = "#820E31"

filtered_pMat$fdr = p.adjust(filtered_pMat$value, method = "fdr")

### Convert filtered_rgMat into matrix

#filtered_rgMat = cast(filtered_rgMat, cohort1 ~ cohort2 )
#
#rownames(filtered_rgMat) = filtered_rgMat[,1]
#filtered_rgMat = filtered_rgMat[,-1]
#
#for (i in 1:ncol(filtered_rgMat)) {
#
#  filtered_rgMat[i, i] = 1
#
#}
#
#for (i in 1:nrow(filtered_rgMat)) {
#
#  is.na(filtered_rgMat[i,]) = 1
#
#}
#

like_by_like_rg = data.frame(cohort1 = unique(filtered_rgMat$cohort1),
							cohort2 = unique(filtered_rgMat$cohort1),
							variable = "rg",
							value = 1,
							colour_code = "#ADADAD")

like_by_like_label = data.frame(cohort1 = unique(filtered_rgMat$cohort1),
							cohort2 = unique(filtered_rgMat$cohort1),
							label = "rg: 1")

like_by_like_p = data.frame(cohort1 = unique(filtered_rgMat$cohort1),
							cohort2 = unique(filtered_rgMat$cohort1),
							variable = "p_h1",
							value = 1,
							fdr = 1)

filtered_rgMat = rbind(filtered_rgMat, like_by_like_rg)
filtered_rgLabel = rbind(filtered_rgLabel, like_by_like_label)
filtered_pMat = rbind(filtered_pMat, like_by_like_p)

nice_phenotype = case_when(phenotype == "META_v4" ~ "Meta-analysis",
							phenotype == "PRIM" ~ "Primary care",
							phenotype == "DRUG" ~ "Prescription drugs",
							phenotype == "IN" ~ "Inpatient",
							phenotype == "INOUT" ~ "Inpatient + outpatient",
							.default = phenotype)

png(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/ldsc_corrplot_", phenotype, "_v04.png"), width = 1500, height = 1500)

ggplot(filtered_rgMat, aes(x = cohort1, y = cohort2, fill = colour_code)) +
    geom_tile() +
    theme_void() +
	scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
			axis.text.y = element_text(size = 16),
        	plot.title = element_text(size = 20)) +
    geom_text(aes(label = filtered_rgLabel$label), size = 7, color = "black") +
    labs(title = paste0(nice_phenotype, " LDSC genetic correlations")) +
    theme(axis.title = element_blank())

dev.off()

corplot_filtered_rgMat = filtered_rgMat

corplot_filtered_rgMat = %>%
	select(-variable) %>%
	pivot_wider(names_from = cohort1, values_from = value) %>%
	arrange(cohort2) %>%
	select(-cohort2) %>%
	mutate(across(everything(), ~ifelse(is.na(.x), 0, .x))) %>%
	as.data.frame() %>%
	`row.names<-`(names(.)) %>%
	as.matrix() %>%
	corrplot::corrplot(is.corr = FALSE)

geom_point(data = subset(filtered_pMat, value < 0.05), 
			aes(x = cohort1, y = cohort2), 
				size = 0.08) +
geom_point(data = subset(filtered_pMat, fdr < 0.05),
			aes(x = cohort1, y = cohort2),
				size = 1.5,
				shape = 1) +