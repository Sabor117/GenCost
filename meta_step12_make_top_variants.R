## ---------------------------
##
## Script name: meta_step12_make_top_variants.R
##
## Purpose of script: Extract lead SNPs and make overlapping "loci" for manuscript
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-09-01
##
## Adapted from coloc_run.R by Dr. Tomoko Nakanishi
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(httr)
library(purrr)
library(fuzzyjoin)

set.seed(117)

### GCTA directory and phenotype lists

outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/gcta_cojo_v4_MAF_0_001/jma_combined/"

pheno = c("IN", "INOUT", "DRUG", "PRIM")
group = c("ALL")

### Define kb-window

kbwindow = 500

kb_name = paste0(kbwindow, "kb")
kb_size = kbwindow * 1000

### Create data frame of all top SNPs from GCTA-COJO

snp_frame = c()

for(i in seq(1,length(pheno))){

	for(j in seq(1, length(group))){

		if(file.exists(paste0(outDir, pheno[i],"_",group[j],"_gcta_jma_out_annotated.txt"))){

			curr_set = fread(paste0(outDir, pheno[i],"_",group[j],"_gcta_jma_out_annotated.txt")) %>% mutate(phenotype = pheno[i]) %>% mutate(group = group[j])
			snp_frame = bind_rows(snp_frame, curr_set) 

		}
	}
}

### Convert to data table

setDT(snp_frame)

### Create data frame of ALL SNPs from GCTA-COJO

all_files = Sys.glob(paste0(outDir, "*_gcta_jma_out_annotated.txt"))
all_files = all_files[!(grepl("_ALL_", all_files))]

all_snp_frame = c()

for(i in 1:length(all_files)){

	curr_set = fread(paste0(all_files[i])) %>% mutate(phenotype = gsub("_gcta_jma_out_annotated.txt", "", basename(all_files[i])))
	all_snp_frame = bind_rows(all_snp_frame, curr_set) 

}

### Convert to data table

setDT(all_snp_frame)


##################### =========================== #####################

### Part 1 = Main phenotypes

##################### =========================== #####################

##### =========================== #####

### Loci sorting

##### =========================== #####

### Sort

snp_frame = snp_frame[order(Chr, bp)]

### Assign raw locus IDs (±250kb grouping within chromosome)

snp_frame[, raw_locus_id := cumsum(c(1, diff(bp) > kb_size)), by = Chr]

### Reassign clean locus IDs (unique across whole dataset, not just per chromosome)

snp_frame[, locus_id := .GRP, by = .(Chr, raw_locus_id)]


##### =========================== #####

### Assign locus names

##### =========================== #####

### First, take the nearest gene as the base locus name

snp_frame[, locus_gene := nearest_ncbi_gene]

### Run two tests to define different loci
### Count how many distinct locus_id each gene appears in

test_gene_loci = snp_frame %>%
					distinct(locus_id, locus_gene) %>%
					group_by(locus_gene) %>%
					summarise(n_loci = n_distinct(locus_id), .groups = "drop") %>%
					arrange(desc(n_loci))

### Show only problematic cases (genes in >1 locus)

problem_genes_loci = test_gene_loci %>% filter(n_loci > 1)

print(problem_genes_loci)

MAF_filter = 0.001

if (MAF_filter == 0.01){

	maf_name = "MAF_0_01"

} else if (MAF_filter == 0.001){

	maf_name = "MAF_0_001"

}

### One instance of a gene in two different loci

if (MAF_filter == 0.001 & kbwindow == 250){

	idx = which(snp_frame$locus_gene == "TTC28")
	snp_frame$locus_gene[idx] = c("TTC28_1", "TTC28_1", "TTC28_2")

}

### Repeat test with locus_ids - find out other loci with multiple genes

test_gene_loci = snp_frame %>%
					distinct(locus_id, locus_gene) %>%
					group_by(locus_id) %>%
					summarise(n_genes = n_distinct(locus_gene), 
							genes = paste(unique(locus_gene), collapse = ", "),
							.groups = "drop") %>%
					arrange(desc(n_genes))

### Show only problematic cases (genes in >1 locus)

problem_genes_names = test_gene_loci %>% filter(n_genes > 1)

print(problem_genes_names)

### Want uniqueness based on locus ID and not just gene name
### I.e. if multiple loci ids share same gene name - distinct locus
### OR multiple gene names for one locus ID
### Manually address loci based on the previous two tests

### MAF FILTER == 0.01

if (MAF_filter == 0.01 & kbwindow == 250){

	snp_frame$locus_gene[snp_frame$locus_id == 111] =  "HLA-DRB5,HLA-DRB6"
	snp_frame$locus_gene[snp_frame$locus_id == 109] =  "MICA-AS1"
	snp_frame$locus_gene[snp_frame$locus_id == 14] =  "MAGI3,PTPN22,HIPK1"
	snp_frame$locus_gene[snp_frame$locus_id == 53] =  "RHOA,BSN,IP6K1"
	snp_frame$locus_gene[snp_frame$locus_id == 100] =  "BTN3A2"
	snp_frame$locus_gene[snp_frame$locus_id == 2] =  "RSRP1"
	snp_frame$locus_gene[snp_frame$locus_id == 10] =  "LOC124904227,LOC124900404"
	snp_frame$locus_gene[snp_frame$locus_id == 40] =  "FIGN,GRB14"
	snp_frame$locus_gene[snp_frame$locus_id == 102] =  "ZSCAN31"
	snp_frame$locus_gene[snp_frame$locus_id == 115] =  "CNPY3,PTK7"
	snp_frame$locus_gene[snp_frame$locus_id == 155] =  "NEK6"
	snp_frame$locus_gene[snp_frame$locus_id == 159] =  "MLLT10"
	snp_frame$locus_gene[snp_frame$locus_id == 169] =  "METTL15"
	snp_frame$locus_gene[snp_frame$locus_id == 178] =  "KCNC2"
	snp_frame$locus_gene[snp_frame$locus_id == 181] =  "MMAB,FAM222A"
	snp_frame$locus_gene[snp_frame$locus_id == 213] =  "FLJ40194,ZNF652"
	snp_frame$locus_gene[snp_frame$locus_id == 216] =  "BPTF"
	snp_frame$locus_gene[snp_frame$locus_id == 226] =  "GDF15"
	snp_frame$locus_gene[snp_frame$locus_id == 234] =  "ZNF831,EDN3"
	snp_frame$locus_gene[snp_frame$locus_id == 237] =  "ADARB1"

}

### MAF FILTER == 0.001

if (MAF_filter == 0.001 & kbwindow == 250){

	snp_frame$locus_gene[snp_frame$locus_id == 110] =  "HLA-DRB5,HLA-DRB6"
	snp_frame$locus_gene[snp_frame$locus_id == 108] =  "MICA-AS1"
	snp_frame$locus_gene[snp_frame$locus_id == 14] =  "MAGI3,PTPN22,HIPK1"
	snp_frame$locus_gene[snp_frame$locus_id == 53] =  "RHOA,BSN,IP6K1"
	snp_frame$locus_gene[snp_frame$locus_id == 2] =  "RSRP1"
	snp_frame$locus_gene[snp_frame$locus_id == 10] =  "LOC124904227,LOC124900404"
	snp_frame$locus_gene[snp_frame$locus_id == 40] =  "FIGN,GRB14"
	snp_frame$locus_gene[snp_frame$locus_id == 100] =  "LOC101928743,LOC124901288"
	snp_frame$locus_gene[snp_frame$locus_id == 102] =  "TOB2P1,ZSCAN31"
	snp_frame$locus_gene[snp_frame$locus_id == 113] =  "CNPY3,PTK7"
	snp_frame$locus_gene[snp_frame$locus_id == 154] =  "NEK6"
	snp_frame$locus_gene[snp_frame$locus_id == 158] =  "MLLT10"
	snp_frame$locus_gene[snp_frame$locus_id == 168] =  "METTL15"
	snp_frame$locus_gene[snp_frame$locus_id == 177] =  "KCNC2"
	snp_frame$locus_gene[snp_frame$locus_id == 180] =  "MMAB,FAM222A"
	snp_frame$locus_gene[snp_frame$locus_id == 212] =  "FLJ40194,ZNF652"
	snp_frame$locus_gene[snp_frame$locus_id == 215] =  "BPTF"
	snp_frame$locus_gene[snp_frame$locus_id == 225] =  "GDF15"
	snp_frame$locus_gene[snp_frame$locus_id == 233] =  "ZNF831,EDN3"
	snp_frame$locus_gene[snp_frame$locus_id == 236] =  "ADARB1"

}

### MAF FILTER == 0.001 + kbwindow == 500

if (MAF_filter == 0.001 & kbwindow == 500){

	snp_frame$locus_gene[snp_frame$locus_id == 110] =  "HLA-DRB5,HLA-DRB6"
	snp_frame$locus_gene[snp_frame$locus_id == 108] =  "MICA-AS1"
	snp_frame$locus_gene[snp_frame$locus_id == 14] =  "MAGI3,PTPN22,HIPK1"
	snp_frame$locus_gene[snp_frame$locus_id == 53] =  "RHOA,BSN,IP6K1"
	snp_frame$locus_gene[snp_frame$locus_id == 2] =  "RSRP1"
	snp_frame$locus_gene[snp_frame$locus_id == 10] =  "LOC124904227,LOC124900404"
	snp_frame$locus_gene[snp_frame$locus_id == 40] =  "FIGN,GRB14"
	snp_frame$locus_gene[snp_frame$locus_id == 100] =  "LOC101928743,LOC124901288"
	snp_frame$locus_gene[snp_frame$locus_id == 102] =  "TOB2P1,ZSCAN31"
	snp_frame$locus_gene[snp_frame$locus_id == 113] =  "CNPY3,PTK7"
	snp_frame$locus_gene[snp_frame$locus_id == 154] =  "NEK6"
	snp_frame$locus_gene[snp_frame$locus_id == 158] =  "MLLT10"
	snp_frame$locus_gene[snp_frame$locus_id == 168] =  "METTL15"
	snp_frame$locus_gene[snp_frame$locus_id == 177] =  "KCNC2"
	snp_frame$locus_gene[snp_frame$locus_id == 180] =  "MMAB,FAM222A"
	snp_frame$locus_gene[snp_frame$locus_id == 212] =  "FLJ40194,ZNF652"
	snp_frame$locus_gene[snp_frame$locus_id == 215] =  "BPTF"
	snp_frame$locus_gene[snp_frame$locus_id == 225] =  "GDF15"
	snp_frame$locus_gene[snp_frame$locus_id == 233] =  "ZNF831,EDN3"
	snp_frame$locus_gene[snp_frame$locus_id == 236] =  "ADARB1"

}


### Locus definition table

locus_table = snp_frame %>%
				group_by(locus_id, locus_gene, Chr) %>%
					summarise(
						start_bp = min(bp),
						end_bp   = max(bp),
						.groups = "drop"
					) %>%
				arrange(Chr, start_bp)

print(locus_table)

### Now for full set of SNPs

### --- Step 1. Expand locus_table with multiple-defined windows ---

locus_windows = locus_table %>%
					mutate(window_start = start_bp - kb_size,
							window_end = end_bp + kb_size)

### --- Step 2. Assign new SNPs to existing loci if within 250kb ---

all_snp_frame$unique_id = paste0(all_snp_frame$SNP, "_", all_snp_frame$phenotype) # Do this because some SNPs appear in multiple analyses

assigned = fuzzy_inner_join(all_snp_frame,
							locus_windows,
							by = c("Chr" = "Chr", "bp" = "window_start", "bp" = "window_end"),
							match_fun = list(`==`, `>=`, `<=`)
							) %>% mutate(locus_id = locus_id, locus_gene = locus_gene)

### Some SNPs are close to two distinct loci - just select closest one
### Select either nearest window, or the window the SNP is inside of

closest_locus_table = assigned %>%
						filter(unique_id %in% unique_id[duplicated(unique_id)]) %>%
						mutate(dist_to_window = case_when(hg38_pos < window_start ~ window_start - hg38_pos,
															hg38_pos > window_end   ~ hg38_pos - window_end,
															TRUE ~ 0 # If inside the window, distance = 0
															)
								) %>%
						arrange(unique_id, dist_to_window) %>%
						group_by(unique_id) %>%
						slice(1) %>%   # keep the closest locus per SNP
						ungroup()

### Removed the duplicated SNPs and then rbind them back in

assigned = assigned[!(assigned$unique_id %in% closest_locus_table$unique_id),]
assigned = rbind(assigned, closest_locus_table[,colnames(assigned)])

### --- Step 3. Identify SNPs that did NOT overlap any locus ---

unassigned = all_snp_frame[!(all_snp_frame$unique_id %in% assigned$unique_id),]

### --- Step 4. For unassigned SNPs, define new loci ---

### Sort

unassigned = unassigned[order(Chr, bp)]

### Assign raw locus IDs (±250kb grouping within chromosome)

unassigned[, raw_locus_id := cumsum(c(1, diff(bp) > kb_size)), by = Chr]

### Reassign clean locus IDs (unique across whole dataset, not just per chromosome)
### But also add to the max value of exisitng locus IDs

unassigned[, locus_id := .GRP, by = .(Chr, raw_locus_id)]
unassigned$locus_id = unassigned$locus_id + max(locus_table$locus_id)

### Locus name is the nearest gene

unassigned[, locus_gene := nearest_ncbi_gene]

### --- Step 5. Combine everything ---

colnames(assigned)[1] = "Chr"
assigned$group = NA
assigned$raw_locus_id = NA
unassigned$group = NA

combined_snp_frame = bind_rows(snp_frame,
								assigned %>% select(names(snp_frame)),
								unassigned %>% select(names(snp_frame))
								)

### Once again perform tests of loci
### Run two tests to define different loci
### Count how many distinct locus_id each gene appears in

test_gene_loci = combined_snp_frame %>%
					distinct(locus_id, locus_gene) %>%
					group_by(locus_gene) %>%
					summarise(n_loci = n_distinct(locus_id), .groups = "drop") %>%
					arrange(desc(n_loci))

### Show only problematic cases (genes in >1 locus)

problem_genes_loci = test_gene_loci %>% filter(n_loci > 1)

print(problem_genes_loci)

### Repeat test with locus_ids - find out other loci with multiple genes

test_gene_loci = combined_snp_frame %>%
					distinct(locus_id, locus_gene) %>%
					group_by(locus_id) %>%
					summarise(n_genes = n_distinct(locus_gene), 
							genes = paste(unique(locus_gene), collapse = ", "),
							.groups = "drop") %>%
					arrange(desc(n_genes))

### Show only problematic cases (genes in >1 locus)

problem_genes_names = test_gene_loci %>% filter(n_genes > 1)

print(problem_genes_names)

### One new locus with multiple genes

if (MAF_filter == 0.01){

	combined_snp_frame$locus_gene[combined_snp_frame$locus_id == 273] =  "HOTTIP"

}

if (MAF_filter == 0.001){

	combined_snp_frame$locus_gene[combined_snp_frame$locus_id == 293] =  "HOTTIP"

}



##### =========================== #####

### Formatting

##### =========================== #####

### Creating tables in format for supplementary data
### ST3

output_frame = snp_frame %>%
	arrange(Chr, bp, locus_id, phenotype) %>%
	mutate(
	phenotype = case_when(
		phenotype == "DRUG"  ~ "Prescription drugs",
		phenotype == "INOUT" ~ "Inpatient + Outpatient",
		phenotype == "PRIM"  ~ "Primary Care",
		phenotype == "IN"    ~ "Inpatient"
	),
	rsid = SNP,
	snpid = paste0(gsub(23, "X", Chr), "_", hg38_pos, "_", toupper(Allele1), "_", toupper(Allele2)),
	chromosome = gsub(23, "X", Chr),
	position_37 = bp,
	position_38 = hg38_pos,
	freq = freq,
	a1 = toupper(Allele1),
	a0 = toupper(Allele2),
	beta = b,
	se = se,
	pval = p,
	n = n
	)

fwrite(output_frame, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/messy_top_variants_main_analysis_", maf_name, "_", kb_name, ".tsv"),
	quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

output_frame = output_frame[,c("rsid", "locus_gene", "snpid", "chromosome", "position_37", "position_38", "a1", "a0", "beta", "se", "pval", "phenotype", "n")]

fwrite(output_frame, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/top_variants_main_analysis_", maf_name, "_", kb_name, ".tsv"),
	quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

### ST4

output_frame_full = combined_snp_frame %>%
	arrange(Chr, bp, locus_id, phenotype) %>%
	mutate(
	nice_phenotype = case_when(
		grepl("DRUG", phenotype) ~ "Prescription drugs",
		grepl("INOUT", phenotype) ~ "Inpatient + Outpatient",
		grepl("PRIM", phenotype) ~ "Primary Care",
		grepl("IN", phenotype) ~ "Inpatient"
	),
	rsid = SNP,
	snpid = paste0(gsub(23, "X", Chr), "_", hg38_pos, "_", toupper(Allele1), "_", toupper(Allele2)),
	chromosome = gsub(23, "X", Chr),
	position_37 = bp,
	position_38 = hg38_pos,
	freq = freq,
	a1 = toupper(Allele1),
	a0 = toupper(Allele2),
	beta = b,
	se = se,
	pval = p,
	n = n
	)

output_frame_full$stratification = case_when(grepl("_M", output_frame_full$phenotype) ~ "(males)",
												grepl("_F", output_frame_full$phenotype) ~ "(females)",
												grepl("_5_18", output_frame_full$phenotype) ~ "(5-18 yo)",
												grepl("_19_35", output_frame_full$phenotype) ~ "(19-35 yo)",
												grepl("_36_55", output_frame_full$phenotype) ~ "(36-55 yo)",
												grepl("_56_75", output_frame_full$phenotype) ~ "(56-75 yo)",
												grepl("_76_95", output_frame_full$phenotype) ~ "(76+ yo)",
												.default = "")

output_frame_full$nice_phenotype = paste0(output_frame_full$nice_phenotype, " ", output_frame_full$stratification)

fwrite(output_frame_full, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/messy_top_variants_ALL_analysis_", maf_name, "_", kb_name, ".tsv"),
	quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

output_frame_full = output_frame_full[,c("nice_phenotype",
											"locus_gene",
											"chromosome",
											"rsid",
											"position_37",
											"a1",
											"a0",
											"freq",
											"beta",
											"se",
											"pval",
											"n",
											"freq_geno",
											"bJ",
											"bJ_se",
											"pJ",
											"LD_r",
											"snpid",
											"l2g_gene",
											"l2g_score",
											"nearest_ncbi_gene",
											"severest_consequence")]

fwrite(output_frame_full, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/top_variants_ALL_analysis_", maf_name, "_", kb_name, ".tsv"),
	quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")

### Also save locus table

fwrite(locus_table, paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/locus_table_main_variants_", maf_name, "_", kb_name, ".tsv"),
	quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")


