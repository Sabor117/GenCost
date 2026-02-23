#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)
library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)


sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Come and see the payments inherent in the system. Help! Help! We're being repressed!!")


#####################

### File locations and script prep

#####################

### Overarching directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"

### Unprocessed phenotype file

unprocessed_phenos_file = paste0(mainDir, "outputs/cost_phenotypes/ALL_cost_phenotypes.tsv")

unprocessed_phenos = fread(unprocessed_phenos_file, data.table = FALSE)

unprocessed_phenos = unprocessed_phenos[!(is.na(unprocessed_phenos$eid)),]

### Various IID lists

hes_iids = fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/hesin/hesin.txt", data.table = FALSE)

gp_iids = fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/gp_data/gp_clinical.txt.gz", data.table = FALSE)
gp_iids = unique(gp_iids$eid)

prescriptions_iids = fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/gp_data/gp_scripts.txt.gz", data.table = FALSE)
prescriptions_iids = unique(prescriptions_iids$eid)

### Getting UKB phenotypes for plotting

ukb_phenos = fread(paste0(mainDir, "processing/ukb_data/ukb671404_cost_phenos.csv"), data.table = FALSE)
colnames(ukb_phenos)[which(colnames(ukb_phenos) == "31-0.0")] = "sex"
colnames(ukb_phenos)[which(colnames(ukb_phenos) == "34-0.0")] = "year_of_birth"
colnames(ukb_phenos)[which(colnames(ukb_phenos) == "52-0.0")] = "month_of_birth"
colnames(ukb_phenos)[which(colnames(ukb_phenos) == "53-0.0")] = "date_of_entry"
ukb_phenos$dob = as.Date(paste0("01/", ukb_phenos$month_of_birth, "/", ukb_phenos$year_of_birth), format = "%d/%m/%Y")

death_reads = fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/healthcare_phenotypes/death/death.txt", data.table = FALSE)
death_reads$date_of_death = as.Date(death_reads$date_of_death, format = "%d/%m/%Y")

### Getting list of withdrawals to be removed after final processing

withdrawn = fread("/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/ukbb_78537_withdrawals.csv", data.table = FALSE)


#####################

### Setting up calculations

#####################

cost_years_phenotype = function(phenotype_frame, phenos_to_include, years_col){

    if (length(phenos_to_include) > 1){

        cost_vector = rowSums(phenotype_frame[,phenos_to_include], na.rm = TRUE)

    } else {

        cost_vector = phenotype_frame[,phenos_to_include]

    }

    cost_calc = cost_vector/phenotype_frame[,years_col] * 365.25

}


#####################

### Creating days of follow-up

#####################

### Fixing date of entry mishaps

unprocessed_phenos$date_of_entry = NULL

unprocessed_phenos = merge(unprocessed_phenos, ukb_phenos[,c("eid", "date_of_entry")], by = "eid", all.x = TRUE)

### End of follow-up for GP/prescription and HES has been defined as a set date
### HES date defined as 1 year prior to most up-to-date HES record: "2021-10-31"
### GP date defined as same date as most recent GP record: "2017-09-30"
### Prescription date defined as same date as most recent prescription record: "2017-07-31"
### Have to make sure that if an individual has died prior to that, the date is replaced by death date

hes_end_date = as.Date("2021-10-31")
gp_end_date = as.Date("2017-09-30")
prescript_end_date = as.Date("2017-07-31")

### Adding lines to fix NAs in end_of_followup for HES and GP

unprocessed_phenos$hes_end_of_followup = hes_end_date
unprocessed_phenos$gp_end_of_followup = gp_end_date

### Remove dates of death after most recent end-date (I.e. HES)

unprocessed_phenos$date_of_death = ifelse(unprocessed_phenos$date_of_death <= hes_end_date, unprocessed_phenos$date_of_death, NA)
unprocessed_phenos$date_of_death = as.Date(unprocessed_phenos$date_of_death) # For some reason above line converts the date into numeric - change it back

### Replace end-of-follow-up with death dates for those remaining

unprocessed_phenos$hes_end_of_followup[!(is.na(unprocessed_phenos$date_of_death))] = unprocessed_phenos$date_of_death[!(is.na(unprocessed_phenos$date_of_death))]

### Remove dates of death after next most recent end-date (GP)

unprocessed_phenos$date_of_death = ifelse(unprocessed_phenos$date_of_death <= gp_end_date, unprocessed_phenos$date_of_death, NA)
unprocessed_phenos$date_of_death = as.Date(unprocessed_phenos$date_of_death) # For some reason above line converts the date into numeric - change it back

### Replace end-of-follow-up with death dates for those remaining

unprocessed_phenos$gp_end_of_followup[!(is.na(unprocessed_phenos$date_of_death))] = unprocessed_phenos$date_of_death[!(is.na(unprocessed_phenos$date_of_death))]

### Remove dates of death after next most recent end-date (prescript)

unprocessed_phenos$date_of_death = ifelse(unprocessed_phenos$date_of_death <= prescript_end_date, unprocessed_phenos$date_of_death, NA)
unprocessed_phenos$date_of_death = as.Date(unprocessed_phenos$date_of_death) # For some reason above line converts the date into numeric - change it back

### Replace end-of-follow-up with death dates for those remaining

unprocessed_phenos$prescription_end_of_followup[!(is.na(unprocessed_phenos$date_of_death))] = unprocessed_phenos$date_of_death[!(is.na(unprocessed_phenos$date_of_death))]

### Make sure everything is date format

unprocessed_phenos$gp_end_of_followup = as.Date(unprocessed_phenos$gp_end_of_followup)
unprocessed_phenos$gp_start_of_followup = as.Date(unprocessed_phenos$gp_start_of_followup)
unprocessed_phenos$prescription_end_of_followup = as.Date(unprocessed_phenos$prescription_end_of_followup)
unprocessed_phenos$prescript_start_of_followup = as.Date(unprocessed_phenos$prescript_start_of_followup)
unprocessed_phenos$hes_end_of_followup = as.Date(unprocessed_phenos$hes_end_of_followup)
unprocessed_phenos$date_of_entry = as.Date(unprocessed_phenos$date_of_entry)

### Start calculating days of follow-up

unprocessed_phenos$gp_days_of_followup = as.numeric(unprocessed_phenos$gp_end_of_followup - unprocessed_phenos$gp_start_of_followup)

unprocessed_phenos$prescript_days_of_followup = as.numeric(unprocessed_phenos$prescription_end_of_followup - unprocessed_phenos$prescript_start_of_followup)

unprocessed_phenos$hes_days_of_followup = as.numeric(unprocessed_phenos$hes_end_of_followup - unprocessed_phenos$date_of_entry)

### Some people have costs for after they died
### Fix this

negative_gp_days = which(unprocessed_phenos$gp_days_of_followup < 0)
negative_script_days = which(unprocessed_phenos$prescript_days_of_followup < 0)

unprocessed_phenos$gp_days_of_followup[negative_gp_days] = NA
unprocessed_phenos$prescript_days_of_followup[negative_script_days] = NA

final_phenotypes = unprocessed_phenos[,c("eid",
                                        "gp_costs",
                                        "gp_days_of_followup",
                                        "prescription_costs",
                                        "prescript_days_of_followup",
                                        "hes_costs",
                                        "weighted_hes_costs",
                                        "hes_days_of_followup")]

fwrite(final_phenotypes, paste0(mainDir, "outputs/cost_phenotypes/ALL_cost_phenotypes_and_followup_2022_pounds.tsv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

### Correction to Euros and to 2019 CPI

cpi_index_2022_to_2019 = 0.89 # Value from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator
pound_to_euro_2019 = 1.1405 # Value based on average exchange rate for GBR -> EUR in 2019 https://www.exchangerates.org.uk/GBP-EUR-spot-exchange-rates-history-2019.html

final_phenotypes$gp_costs = final_phenotypes$gp_costs * cpi_index_2022_to_2019 * pound_to_euro_2019
final_phenotypes$prescription_costs = final_phenotypes$prescription_costs * cpi_index_2022_to_2019 * pound_to_euro_2019
final_phenotypes$hes_costs = final_phenotypes$hes_costs * cpi_index_2022_to_2019 * pound_to_euro_2019


#####################

### Running calculations

#####################

### Calculating individual cost per year for each of three phenotypes

final_phenotypes$gp_costs_year = cost_years_phenotype(final_phenotypes, c("gp_costs"), "gp_days_of_followup")
final_phenotypes$hes_costs_year = cost_years_phenotype(final_phenotypes, c("hes_costs"), "hes_days_of_followup")
final_phenotypes$prescription_costs_year = cost_years_phenotype(final_phenotypes, c("prescription_costs"), "prescript_days_of_followup")

### Get list of IIDs with < 5 years follow-up

followup_minimum = 365.25 * 5

hes_less_4_years_iids = final_phenotypes$eid[final_phenotypes$hes_days_of_followup < followup_minimum & !(is.na(final_phenotypes$hes_days_of_followup))]
gp_less_4_years_iids = final_phenotypes$eid[final_phenotypes$gp_days_of_followup < followup_minimum & !(is.na(final_phenotypes$gp_days_of_followup))]
presript_less_4_years_iids = final_phenotypes$eid[final_phenotypes$prescript_days_of_followup < followup_minimum & !(is.na(final_phenotypes$prescript_days_of_followup))]

### HES frame

hes_based_frame = final_phenotypes[,c("eid", "gp_costs_year", "hes_costs_year", "prescription_costs_year")]

### Then remove IIDs with < 4 years follow-up

hes_based_frame = hes_based_frame[!(hes_based_frame$eid %in% hes_less_4_years_iids),]

### Then calculate total costs and summaries

hes_based_frame$total_costs_year = rowSums(hes_based_frame[,c("gp_costs_year", "hes_costs_year", "prescription_costs_year")], na.rm = TRUE)

### Remove NAs and -ve values
### Replace with 1 for creatng natural log

hes_based_frame$gp_costs_year[hes_based_frame$gp_costs_year < 1] = 1
hes_based_frame$prescription_costs_year[hes_based_frame$prescription_costs_year < 1] = 1
hes_based_frame$hes_costs_year[is.na(hes_based_frame$hes_costs_year)] = 1

hes_based_frame$gp_costs_year[hes_based_frame$gp_costs_year < 1] = 1
hes_based_frame$hes_costs_year[hes_based_frame$hes_costs_year < 1] = 1
hes_based_frame$prescription_costs_year[hes_based_frame$prescription_costs_year < 1] = 1
hes_based_frame$total_costs_year[hes_based_frame$total_costs_year < 1] = 1

hes_based_frame$log_hes_costs_year = log(hes_based_frame$hes_costs_year)
hes_based_frame$log_gp_costs_year = log(hes_based_frame$gp_costs_year)
hes_based_frame$log_prescription_costs_year = log(hes_based_frame$prescription_costs_year)
hes_based_frame$log_total_costs_year = log(hes_based_frame$total_costs_year)


### GP frame

gp_based_frame = final_phenotypes[,c("eid", "gp_costs_year", "hes_costs_year", "prescription_costs_year")]

gp_based_frame = gp_based_frame[gp_based_frame$eid %in% gp_iids,]

### Adding in missing IIDs based on GP data
### Has been done at previous stage

if (any(!(gp_iids %in% gp_based_frame$eid))){

    missing_gp_frame = data.frame(eid = unique(gp_iids[which(!(gp_iids %in% gp_based_frame$eid))]),
                            gp_costs_year = -1,
                            hes_costs_year = -1,
                            prescription_costs_year = -1
                            )

    gp_based_frame = rbind(gp_based_frame, missing_gp_frame)

}

### Remove those with insufficient follow-up

gp_based_frame = gp_based_frame[!(gp_based_frame$eid %in% gp_less_4_years_iids),]

### Then calculate total costs and summaries

gp_based_frame$total_costs_year = rowSums(gp_based_frame[,c("gp_costs_year", "hes_costs_year", "prescription_costs_year")], na.rm = TRUE)
gp_based_frame$gp_and_prescript_year = rowSums(gp_based_frame[,c("gp_costs_year", "prescription_costs_year")], na.rm = TRUE)

### Remove NAs and -ve values
### Replace with 1 for creatng natural log

gp_based_frame$gp_costs_year[is.na(gp_based_frame$gp_costs_year)] = 1
gp_based_frame$hes_costs_year[is.na(gp_based_frame$hes_costs_year)] = 1
gp_based_frame$prescription_costs_year[is.na(gp_based_frame$prescription_costs_year)] = 1

gp_based_frame$gp_costs_year[gp_based_frame$gp_costs_year < 1] = 1
gp_based_frame$hes_costs_year[gp_based_frame$hes_costs_year < 1] = 1
gp_based_frame$prescription_costs_year[gp_based_frame$prescription_costs_year < 1] = 1
gp_based_frame$total_costs_year[gp_based_frame$total_costs_year < 1] = 1
gp_based_frame$gp_and_prescript_year[gp_based_frame$gp_and_prescript_year < 1] = 1

gp_based_frame$log_hes_costs_year = log(gp_based_frame$hes_costs_year)
gp_based_frame$log_gp_costs_year = log(gp_based_frame$gp_costs_year)
gp_based_frame$log_prescription_costs_year = log(gp_based_frame$prescription_costs_year)
gp_based_frame$log_total_costs_year = log(gp_based_frame$total_costs_year)
gp_based_frame$log_gp_and_prescript_year = log(gp_based_frame$gp_and_prescript_year)


### Prescription frame

prescript_based_frame = final_phenotypes[,c("eid", "prescription_costs_year")]

prescript_based_frame = prescript_based_frame[prescript_based_frame$eid %in% prescriptions_iids,]

### Remove those with insufficient follow-up

prescript_based_frame = prescript_based_frame[!(prescript_based_frame$eid %in% presript_less_4_years_iids),]

### Remove NAs and -ve values
### Replace with 1 for creatng natural log

prescript_based_frame$prescription_costs_year[is.na(prescript_based_frame$prescription_costs_year)] = 1

prescript_based_frame$prescription_costs_year[prescript_based_frame$prescription_costs_year < 1] = 1

prescript_based_frame$log_prescription_costs_year = log(prescript_based_frame$prescription_costs_year)


### WRITING

hes_based_frame$FID = hes_based_frame$eid
hes_based_frame$IID = hes_based_frame$eid

hes_based_frame = hes_based_frame[!(hes_based_frame$IID %in% withdrawn$withdrawn),]

gp_based_frame$FID = gp_based_frame$eid
gp_based_frame$IID = gp_based_frame$eid

gp_based_frame = gp_based_frame[!(gp_based_frame$IID %in% withdrawn$withdrawn),]

prescript_based_frame$FID = prescript_based_frame$eid
prescript_based_frame$IID = prescript_based_frame$eid

prescript_based_frame = prescript_based_frame[!(prescript_based_frame$IID %in% withdrawn$withdrawn),]

fwrite(hes_based_frame,
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_hes_based_phenotype_frame.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
fwrite(gp_based_frame,
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_gp_based_phenotype_frame.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
fwrite(prescript_based_frame,
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_prescription_based_phenotype_frame.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

fwrite(hes_based_frame[,c("FID", "IID", "log_hes_costs_year", "log_total_costs_year")],
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_hes_based_phenotype_frame_REGENIE.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
fwrite(gp_based_frame[,c("FID", "IID", "log_gp_costs_year", "log_total_costs_year")],
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_gp_based_phenotype_frame_REGENIE.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
fwrite(prescript_based_frame[,c("FID", "IID", "log_prescription_costs_year")],
        paste0(mainDir, "outputs/cost_phenotypes/FINAL_prescription_based_phenotype_frame_REGENIE.tsv"),
        row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")



#####################

### Making densities

#####################

gp_based_plot_frame = melt(gp_based_frame[,c("log_gp_costs_year", "log_total_costs_year", "log_gp_and_prescript_year")])
gp_based_plot_frame$Phenotype = ifelse(gp_based_plot_frame$variable == "log_gp_costs_year", "GP costs", ifelse(gp_based_plot_frame$variable == "log_total_costs_year", "Total costs", "GP + prescription costs"))

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/gp_costs_densities_EUR_2019.pdf"),
    height = 12, width = 12)

ggplot(gp_based_plot_frame, aes(x = value, fill = Phenotype)) +
    geom_density(color = "black", alpha = .6) + 
    scale_fill_manual(breaks = unique(gp_based_plot_frame$Phenotype), values = c("#FAD510", "#CB2314", "#3B9AB2")) +
    labs(title = "GP-based costs density", x = "Log cost (EUR, 2019)", y = "Frequency") +
    theme_classic()

dev.off()


hes_based_plot_frame = melt(hes_based_frame[,c("log_hes_costs_year", "log_total_costs_year")])
hes_based_plot_frame$Phenotype = ifelse(hes_based_plot_frame$variable == "log_total_costs_year", "Total costs", "HES costs")

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/hes_costs_densities_EUR_2019.pdf"),
    height = 12, width = 12)

ggplot(hes_based_plot_frame, aes(x = value, fill = Phenotype)) +
    geom_density(color = "black", alpha = .6) + 
    scale_fill_manual(breaks = unique(hes_based_plot_frame$Phenotype), values = c("#354823", "#C085E6", "grey")) +
    labs(title = "HES-based costs density", x = "Log cost (EUR, 2019)", y = "Frequency") +
    theme_classic()

dev.off()


prescript_based_plot_frame = prescript_based_frame
prescript_based_plot_frame$Phenotype = "Prescription costs"

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/prescript_costs_densities_EUR_2019.pdf"),
    height = 12, width = 12)

ggplot(prescript_based_plot_frame, aes(x = log_prescription_costs_year, fill = Phenotype)) +
    geom_density(color = "black", alpha = .6) + 
    scale_fill_manual(breaks = unique(prescript_based_plot_frame$Phenotype), values = c("grey")) +
    labs(title = "Prescription costs density", x = "Log cost (EUR, 2019)", y = "Frequency") +
    theme_classic()

dev.off()


#gp_based_plot_frame$Phenotype = ifelse(gp_based_plot_frame$Phenotype == "Total costs", "Total costs (GP)", gp_based_plot_frame$Phenotype)
#hes_based_plot_frame$Phenotype = ifelse(hes_based_plot_frame$Phenotype == "Total costs", "Total costs (HES)", hes_based_plot_frame$Phenotype)

gp_based_plot_frame_all = gp_based_plot_frame[!(gp_based_plot_frame$Phenotype == "GP + prescription costs"),]
hes_based_plot_frame_all = hes_based_plot_frame[hes_based_plot_frame$Phenotype == "HES costs",]
prescript_based_plot_frame_all = data.frame(variable = "prescription_costs_year_log",
                                            value = prescript_based_plot_frame$log_prescription_costs_year,
                                            Phenotype = prescript_based_plot_frame$Phenotype)

total_plot_frame = rbind(hes_based_plot_frame_all, gp_based_plot_frame_all, prescript_based_plot_frame_all)

total_plot_frame_tails = total_plot_frame[total_plot_frame$value > 10,]

for (i in 1:length(unique(total_plot_frame$Phenotype))){

    curr_pheno = unique(total_plot_frame$Phenotype)[i]

    n_pheno = length(which(total_plot_frame$Phenotype == curr_pheno))

    total_plot_frame$Phenotype = gsub(curr_pheno, paste0(curr_pheno, ", n = ", n_pheno), total_plot_frame$Phenotype)

}


for (i in 1:length(unique(total_plot_frame_tails$Phenotype))){

    curr_pheno = unique(total_plot_frame_tails$Phenotype)[i]

    n_pheno = length(which(total_plot_frame_tails$Phenotype == curr_pheno))

    total_plot_frame_tails$Phenotype = gsub(curr_pheno, paste0(curr_pheno, ", n = ", n_pheno), total_plot_frame_tails$Phenotype)

}


pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/all_costs_densities_EUR_2019.pdf"),
    height = 12, width = 12)

ggplot(total_plot_frame, aes(x = value, fill = Phenotype)) +
    geom_density(color = "black", alpha = .6) + 
    scale_fill_manual(breaks = unique(total_plot_frame$Phenotype), values = c("#C085E6", "grey", "#FAD510", "#CB2314")) +
    labs(title = "Costs density", x = "Log cost", y = "Frequency") +
    theme_classic()

dev.off()

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/all_costs_densities_tails_EUR_2019.pdf"),
    height = 12, width = 12)

ggplot(total_plot_frame_tails, aes(x = value, fill = Phenotype)) +
    geom_density(color = "black", alpha = .6) + 
    scale_fill_manual(breaks = unique(total_plot_frame_tails$Phenotype), values = c("#C085E6", "grey", "#FAD510", "#CB2314")) +
    labs(title = "Costs density", x = "Log cost", y = "Frequency") +
    theme_classic()

dev.off()


#####################

### Making Jiwoo-densities

#####################

jiwoo_style_plot_frame = hes_based_frame

jiwoo_style_phenotype = "log_hes_costs_year"

selected_end_of_followup = hes_end_date

options(scipen=999)

a = ggplot() +
	geom_histogram(data = jiwoo_style_plot_frame, mapping = aes(x = log_hes_costs_year), fill = "black", alpha = 0.25) +
	geom_density(data = jiwoo_style_plot_frame, mapping = aes(x = log_hes_costs_year, y = ..density..*nrow(jiwoo_style_plot_frame)*0.4), fill = "black", alpha = 0.5, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 140000), breaks = seq(0, 140000, 10000)) +
	labs(x = "", y = "", fill = "", tag = "A") + 
	theme(plot.tag = element_text(size = 25),
            legend.position = "bottom",
            axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"),
            axis.title.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title.y = element_text(size = 15, color = "black"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 25),
            panel.grid.major = element_line(size = 0.1, color = "gray"),
            panel.grid.minor = element_line(size = 0.1, color = "gray"),
            panel.background = element_blank(),
            strip.background = element_rect(fill = "gray95"),
            axis.line = element_line(colour = "black"))

jiwoo_style_plot_frame = merge(jiwoo_style_plot_frame, ukb_phenos[,c("eid", "sex")], by = "eid")

b = ggplot() +
	geom_histogram(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$sex == 0),], mapping = aes(x = log_hes_costs_year, fill = "Female"), alpha = 0.25) +
	geom_density(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$sex == 0),], mapping = aes(x = log_hes_costs_year, y = ..density..*length(unique(jiwoo_style_plot_frame$eid[which(jiwoo_style_plot_frame$sex == 0)]))*0.4, fill = "Female"), alpha = 0.5, adjust = 5) +
	geom_histogram(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$sex == 1),], mapping = aes(x = log_hes_costs_year, fill = "Male"), alpha = 0.25) +
	geom_density(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$sex == 1),], mapping = aes(x = log_hes_costs_year, y = ..density..*length(unique(jiwoo_style_plot_frame$eid[which(jiwoo_style_plot_frame$sex == 1)]))*0.4, fill = "Male"), alpha = 0.5, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 140000), breaks = seq(0, 140000, 10000)) +
	scale_fill_manual(limits = c("Female", "Male"), values = c("dark green", "yellowgreen")) +
	labs(x = "", y = "", fill = "", tag = "B") + 
	theme(plot.tag = element_text(size = 25),
            legend.position = "bottom",
            axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"),
            axis.title.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title.y = element_text(size = 15, color = "black"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 25),
            panel.grid.major = element_line(size = 0.1, color = "gray"),
            panel.grid.minor = element_line(size = 0.1, color = "gray"),
            panel.background = element_blank(),
            strip.background = element_rect(fill = "gray95"),
            axis.line = element_line(colour = "black"))

jiwoo_style_plot_frame$frame_end_of_followup = selected_end_of_followup

death_reads_jiwoo = death_reads[death_reads$date_of_death <= selected_end_of_followup,]
death_reads_jiwoo = death_reads_jiwoo[death_reads_jiwoo$eid %in% jiwoo_style_plot_frame$eid, c("eid", "date_of_death")]
death_reads_jiwoo = unique(death_reads_jiwoo)

jiwoo_style_plot_frame = merge(jiwoo_style_plot_frame, death_reads_jiwoo, by = "eid", all = TRUE)

jiwoo_style_plot_frame$frame_end_of_followup[!(is.na(jiwoo_style_plot_frame$date_of_death))] = jiwoo_style_plot_frame$date_of_death[!(is.na(jiwoo_style_plot_frame$date_of_death))]

jiwoo_style_plot_frame = merge(jiwoo_style_plot_frame, ukb_phenos[,c("eid", "dob")], by = "eid")

jiwoo_style_plot_frame$age_end_of_followup = as.numeric((jiwoo_style_plot_frame$frame_end_of_followup - jiwoo_style_plot_frame$dob)/365.25)

c = ggplot() +
	geom_histogram(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followup > 60),], mapping = aes(x = log_hes_costs_year, fill = ">60 years old"), alpha = 0.5) +
	geom_density(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followup > 60),], mapping = aes(x = log_hes_costs_year, y = ..density..*length(unique(jiwoo_style_plot_frame$eid[which(jiwoo_style_plot_frame$age_end_of_followup > 60)]))*0.5, fill = ">60 years old"), alpha = 0.25, adjust = 5) +
	geom_histogram(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followup < 60 & jiwoo_style_plot_frame$age >= 30),], mapping = aes(x = log_hes_costs_year, fill = "30-60 years old"), alpha = 0.5) +
	geom_density(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followup < 60 & jiwoo_style_plot_frame$age >= 30),], mapping = aes(x = log_hes_costs_year, y = ..density..*length(unique(jiwoo_style_plot_frame$eid[which(jiwoo_style_plot_frame$age_end_of_followup < 60 & jiwoo_style_plot_frame$age_end_of_followup >= 30)]))*0.5, fill = "30-60 years old"), alpha = 0.25, adjust = 5) +
	geom_histogram(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followupage < 30),], mapping = aes(x = log_hes_costs_year, fill = "<30 years old"), alpha = 0.5) +
	geom_density(data = jiwoo_style_plot_frame[which(jiwoo_style_plot_frame$age_end_of_followup < 30),], mapping = aes(x = log_hes_costs_year, y = ..density..*length(unique(jiwoo_style_plot_frame$eid[which(jiwoo_style_plot_frame$age_end_of_followup < 30)]))*0.5, fill = "<30 years old"), alpha = 0.25, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 140000), breaks = seq(0, 140000, 10000)) +
	scale_fill_manual(limits = c("<30 years old", "30-60 years old", ">60 years old"), values = c("coral", "red", "red4")) +
	labs(x = "", y = "", fill = "", tag = "C") + 
	theme(plot.tag = element_text(size = 25),
            legend.position = "bottom",
            axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"),
            axis.title.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title.y = element_text(size = 15, color = "black"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 25),
            panel.grid.major = element_line(size = 0.1, color = "gray"),
            panel.grid.minor = element_line(size = 0.1, color = "gray"),
            panel.background = element_blank(),
            strip.background = element_rect(fill = "gray95"),
            axis.line = element_line(colour = "black"))

pdf(paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/ukb_phenotype_plots/inpatient_cost_densities_jiwoo_style_EUR_2019.pdf"),
    height = 24, width = 32)

grid.arrange(a, b, c, nrow = 1, left = textGrob("Number of Individuals",
                rot = 90, gp = gpar(fontsize = 25)),
                bottom = textGrob("Log Annual Healthcare Cost in Pounds",
                gp = gpar(fontsize = 25)))

dev.off()



