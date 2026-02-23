##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(ggplot2)
library(dplyr)
library(wesanderson, lib.loc = "/projappl/project_2007428/RPackages_421/")
options(scipen = 5)
set.seed(117)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL_v4/")
sumstatDir = paste0(mainDir, "processing/meta_sumstats/")

figureDir = paste0(mainDir, "outputs/figures/se_n_plots/")

sumstat_list = Sys.glob(paste0(metaDir, "*TBL.gz"))

wes_selected_pallete = c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"), wes_palette("BottleRocket2"))


##### =========================== #####

### Create plots

##### =========================== #####

cohort_colours = c("UKB_EUR" = wes_selected_pallete[1],
				"UKB_AFR" = wes_selected_pallete[1],
				"UKB_CSA" = wes_selected_pallete[1],
				"UKB_EAS" = wes_selected_pallete[1],
				"FINNGEN" = wes_selected_pallete[2],
				"CHB" = wes_selected_pallete[3],
				"GS20K" = wes_selected_pallete[4],
				"MGBB" = wes_selected_pallete[5],
				"AGDS" = wes_selected_pallete[6],
				"AOU" = wes_selected_pallete[7],
				"QGP" = wes_selected_pallete[8],
				"GBP" = wes_selected_pallete[9],
				"NTR" = wes_selected_pallete[10],
				"EstBB" = wes_selected_pallete[11],
				"GNH" = wes_selected_pallete[12],
				"GE" = wes_selected_pallete[13]
)

for (i in 1:length(sumstat_list)){

	##### =========================== #####

	### Data refinement

	##### =========================== #####

	curr_meta_analysis = fread(sumstat_list[i], data.table = FALSE)
	colnames(curr_meta_analysis)[10] = "Pval"

	curr_meta_analysis = curr_meta_analysis[curr_meta_analysis$Pval < 5e-8,]

	analysis_name = strsplit(basename(sumstat_list[i]), "_metal")[[1]][1]

	curr_sumstatDir = paste0(sumstatDir, analysis_name, "/")
	curr_fileSet = Sys.glob(paste0(curr_sumstatDir, "*.txt.gz"))

	plot_snp_set = data.frame(matrix(ncol = 6, nrow = 0))
	colnames(plot_snp_set) = c("snpid", "beta1", "se", "n", "p", "cohort")

	for (j in 1:length(curr_fileSet)){

		curr_cohort_file = curr_fileSet[j]

		curr_cohort = fread(curr_cohort_file, data.table = FALSE)

		cohort_name = strsplit(basename(curr_cohort_file), "\\.")[[1]][1]

		if (cohort_name == "UKB"){

			cohort_name = paste0(cohort_name, "_", strsplit(basename(curr_cohort_file), "\\.")[[1]][6])

		}

		curr_cohort = curr_cohort[curr_cohort$snpid %in% curr_meta_analysis$MarkerName,]

		currsnp_set = data.frame(snpid = curr_cohort$snpid,
									beta1 = curr_cohort$beta1,
									se = curr_cohort$se,
									n = curr_cohort$n,
									p = curr_cohort$p,
									cohort = cohort_name)

		plot_snp_set = rbind(plot_snp_set, currsnp_set)

	}

	plot_snp_set$sqrtn = sqrt(plot_snp_set$n)

	plot_snp_set$se_dev_one = 1 / plot_snp_set$se

	shared_snps = plot_snp_set %>%
					group_by(snpid) %>%
					filter(n_distinct(cohort) == n_distinct(plot_snp_set$cohort))

	select_snp = sample(unique(shared_snps$snpid), 1)
	select_10snp = sample(unique(shared_snps$snpid), 10)

	randomly_selected_1snp = shared_snps[shared_snps$snpid == select_snp,]
	randomly_selected_10snp = shared_snps[shared_snps$snpid %in% select_10snp,]

	top_snps_meta = curr_meta_analysis[curr_meta_analysis$MarkerName %in% unique(shared_snps$snpid),]

	# Sorting the data frame by P-value
	top_snps_meta = top_snps_meta[order(top_snps_meta$Pval),]

	# Create an empty data frame to store selected SNPs
	selected_snps = data.frame()

	# Initialize the loop to select 1-5 SNPs
	for (j in 1:nrow(top_snps_meta)) {
	
	# Check if we've already selected 5 SNPs
		if (nrow(selected_snps) >= 5) {
			break
		}
		
		# Get the SNP to be tested
		current_snp = top_snps_meta[j, ]
		
		# If no SNPs have been selected yet, add the first one
		if (nrow(selected_snps) == 0) {

			selected_snps <- rbind(selected_snps, current_snp)

		} else {

			# Check if the current SNP is on a different chromosome or is at least 500,000 base pairs away
			valid = TRUE
			
			for (k in 1:nrow(selected_snps)) {

				selected_snp = selected_snps[k, ]
			
				# If the SNP is on the same chromosome and within 500,000 base pairs, it's not valid
				if (selected_snp$Chromosome == current_snp$Chromosome && abs(selected_snp$Position - current_snp$Position) < 500000) {

					valid = FALSE
					break

				}
			}
			
			# If the SNP is valid, add it to the selected SNPs list
			if (valid) {

				selected_snps = rbind(selected_snps, current_snp)

			}
		}
	}

	snp_list = selected_snps$MarkerName
	top_select_snps = shared_snps[shared_snps$snpid %in% snp_list,]


	##### =========================== #####

	### Plot names and files

	##### =========================== #####

	box_plotname = paste0("SE x N plot of ", analysis_name, ". Total sig SNPs = ", length(unique(plot_snp_set$snpid)))
	beta_plotname = paste0("Beta x N plot of ", analysis_name, ". SNP = ", select_snp)
	se_plotname = paste0("SE x N plot of ", analysis_name, ". SNP = ", select_snp)
	beta_10plotname = paste0("Beta x N plot of ", analysis_name)
	se_10plotname = paste0("SE x N plot of ", analysis_name)

	boxplot_savename = paste0(figureDir, analysis_name, "_se_n_boxplot_ALLSNP.png")
	beta_plot_savename = paste0(figureDir, analysis_name, "_beta_n_plot_onesnp.png")
	se_plot_savename = paste0(figureDir, analysis_name, "_se_n_plot_onesnp.png")
	beta_10plot_savename = paste0(figureDir, analysis_name, "_beta_n_plot_10snp.png")
	se_10plot_savename = paste0(figureDir, analysis_name, "_se_n_plot_10snp.png")

	beta_topSNPplotname = paste0("Beta x N plot of ", analysis_name)
	se_topSNPplotname = paste0("SE x N plot of ", analysis_name)

	beta_topSNPplot_savename = paste0(figureDir, analysis_name, "_beta_n_plot_topsnp.png")
	se_topSNPplot_savename = paste0(figureDir, analysis_name, "_se_n_plot_topsnp.png")


	##### =========================== #####

	### SE boxplot

	##### =========================== #####

	png(boxplot_savename, width = 1000, height = 1000)

	print(
		ggplot(plot_snp_set, aes(x = sqrtn, y = se, color = cohort)) + 
			geom_boxplot() +
			theme_minimal() +
			labs(title = box_plotname, x = "sqrt(N)", y = "SE") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours) +
			coord_trans(y = "log")
	)

	dev.off()


	##### =========================== #####

	### Beta 1 SNP dot plot

	##### =========================== #####

	png(beta_plot_savename, width = 1000, height = 1000)

	print(
		ggplot(randomly_selected_1snp, aes(x = sqrtn, y = beta1, color = cohort)) + 
			geom_point(size = 4) +
			geom_errorbar(aes(ymin = beta1 - se, ymax = beta1 + se), width = 0.1) +
			theme_minimal() +
			labs(title = beta_plotname, x = "sqrt(N)", y = "Beta") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours)
	)

	dev.off()


	##### =========================== #####

	### SE 1 SNP dot plot

	##### =========================== #####

	png(se_plot_savename, width = 1000, height = 1000)

	print(
		ggplot(randomly_selected_1snp, aes(x = sqrtn, y = se)) + 
			geom_point(size = 4, aes(color = cohort)) +
			geom_smooth(method = "lm", se = FALSE) +
			theme_minimal() +
			labs(title = se_plotname, x = "sqrt(N)", y = "SE") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours)
	)

	dev.off()


	##### =========================== #####

	### Beta 10 SNP dot plot

	##### =========================== #####

	png(beta_10plot_savename, width = 1000, height = 1000)

	print(
		ggplot(randomly_selected_10snp, aes(x = sqrtn, y = beta1, color = cohort)) + 
			geom_point(size = 4) +
			theme_minimal() +
			labs(title = beta_10plotname, x = "sqrt(N)", y = "Beta") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours)
	)

	dev.off()


	##### =========================== #####

	### SE 10 SNP dot plot

	##### =========================== #####

	png(se_10plot_savename, width = 1000, height = 1000)

	print(
		ggplot(randomly_selected_10snp, aes(x = sqrtn, y = se)) + 
			geom_point(size = 4, aes(color = cohort)) +
			geom_smooth(method = "lm", se = FALSE) +
			theme_minimal() +
			labs(title = se_10plotname, x = "sqrt(N)", y = "SE") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours)
	)

	dev.off()

	##### =========================== #####

	### Beta top SNP dot plot

	##### =========================== #####

	png(beta_topSNPplot_savename, width = 1000, height = 1000)

	print(
		ggplot(top_select_snps, aes(x = sqrtn, y = beta1, color = cohort, shape = snpid)) + 
			geom_point(size = 4) +
			theme_classic() +
			labs(title = beta_topSNPplotname, x = "sqrt(N)", y = "Beta") + 
			theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
					axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
					plot.title = element_text(size = 30, face = "bold"),
					legend.text = element_text(size = 20),  # Adjust the size of legend text
					legend.title = element_text(size = 24, face = "bold"),
					axis.text.x = element_text(size = 20),
					axis.text.y = element_text(size = 20)) +
			scale_color_manual(values = cohort_colours)
	)

	dev.off()


	##### =========================== #####

	### SE top SNP dot plot

	##### =========================== #####

	png(se_topSNPplot_savename, width = 1000, height = 1000)

	print(
		ggplot(top_select_snps, aes(x = sqrtn, y = se_dev_one, shape = snpid)) + 
		geom_point(size = 4, aes(color = cohort)) +
		geom_smooth(method = "lm", se = FALSE, aes(group = 1)) +
		theme_classic() +
		labs(title = se_topSNPplotname, x = "sqrt(N)", y = "1/SE") + 
		theme(axis.title.x = element_text(size = 24, face = "bold"),  # Change x-axis title size
				axis.title.y = element_text(size = 24, face = "bold"),  # Change y-axis title size
				plot.title = element_text(size = 30, face = "bold"),
				legend.text = element_text(size = 20),  # Adjust the size of legend text
				legend.title = element_text(size = 24, face = "bold"),
				axis.text.x = element_text(size = 20),
				axis.text.y = element_text(size = 20)) +
		scale_color_manual(values = cohort_colours)
	)

	dev.off()

}

