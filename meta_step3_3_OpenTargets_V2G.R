##### =========================== #####

### Setup environment

##### =========================== #####

### Run once

#library(devtools)
#install_github("amirfeizi/otargen", lib="/projappl/project_2007428/RPackages_421/")
library(data.table)
library(stringr)
library(httr)
library(jsonlite)
library(purrr)
library(dplyr)
library(otargen, lib.loc = "/projappl/project_2007428/RPackages_421/")

options(scipen = 999) # prevent scientific notation in output

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gctaDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/")
outDir = paste0(gctaDir, "jma_combined/")

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

sessionInfo()


##### =========================== #####

### Combine GCTA outputs

##### =========================== #####

RUN_INTERACTIVE = FALSE

if (isTRUE(RUN_INTERACTIVE)){

	system(paste0("mkdir -p ", gctaDir, "jma_combined/"))
	system(paste0("mkdir -p ", gctaDir, "jma_files/"))
	system(paste0("mkdir -p ", gctaDir, "ldr_cma/"))
	system(paste0("mkdir -p ", gctaDir, "logs/"))

	all_analysis_files = Sys.glob(paste0(gctaDir, "*jma.cojo"))

	all_analyses = unique(str_split_fixed(basename(all_analysis_files), "_chr", 2)[,1])

	system(paste0("mv ", gctaDir, "*.log ", gctaDir, "logs/"))
	system(paste0("mv ", gctaDir, "*.badsnps ", gctaDir, "logs/"))

	for (i in 1:length(all_analyses)){

		curr_analysis = all_analyses[i]

		cat(paste0("\nWorking on ", curr_analysis, "\n\n"))

		analysis_files = all_analysis_files[grepl(curr_analysis, all_analysis_files)]

		if (curr_analysis == "IN_ALL" | curr_analysis == "DRUG_ALL" | curr_analysis == "INOUT_ALL"){

			analysis_files = analysis_files[!(grepl("noUKB", analysis_files))]
			analysis_files = analysis_files[!(grepl("noFINNGEN", analysis_files))]

		}

		cat(paste0("\n", curr_analysis, " has ", length(analysis_files), " files.\n\n"))
		print(analysis_files)

		if (length(analysis_files) == 0){

			cat("\nNo .jma files. Next.")

			next

		}

		curr_analysis_full_jma = fread(analysis_files[1], data.table = FALSE)

		if (length(analysis_files) > 1){

			for (j in 2:length(analysis_files)){

			currFile = fread(analysis_files[j], data.table = FALSE)

			curr_analysis_full_jma = rbind(curr_analysis_full_jma, currFile)

			}

		}

		fwrite(curr_analysis_full_jma, paste0(outDir, curr_analysis, "_gcta_jma_out.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

	}

	system(paste0("mv ", gctaDir, "*cma.cojo ", gctaDir, "ldr_cma/"))
	system(paste0("mv ", gctaDir, "*ldr.cojo ", gctaDir, "ldr_cma/"))
	system(paste0("mv ", gctaDir, "*.jma.cojo ", gctaDir, "jma_files/"))

}

all_full_jma = Sys.glob(paste0(outDir, "*_gcta_jma_out.txt"))


##### =========================== #####

### Get V2G and closest genes for each SNP via OTG

### V2G was deprecated in late 2024 - this was replaced by L2G

##### =========================== #####

### Set base URL of GraphQL API endpoint

base_url = "https://api.platform.opentargets.org/api/v4/graphql"

### Define whether to run interactively - for testing purposes

RUN_INTERACTIVE = TRUE

if (isTRUE(RUN_INTERACTIVE)){

	### Loop over every GCTA file

	for (i in 1:length(all_full_jma)){

		### Read file and define output to extract

		current_analysis = fread(all_full_jma[i], data.table = FALSE)

		heading(paste0("Currently running V2G on ", basename(all_full_jma[i])))

		query_frame = data.frame(matrix(ncol = 5, nrow = 0))
		colnames(query_frame) = c("hg38_snpid", "l2g_gene", "l2g_score", "nearest_gene", "severest_consequence")

		### Loop over every SNP in GCTA file

		for (j in 1:nrow(current_analysis)){

			### Define current SNP

			query_rsID = current_analysis$SNP[j]

			### Get OpenTargets SNP ID

			id_query_string = "query SearchQuery($queryString: String!, $index: Int!, $entityNames: [String!]!) {
									search(
										queryString: $queryString
										entityNames: $entityNames
										page: {index: $index, size: 10}
									) {
										total
										hits {
										id
										highlights
										object {
											... on Variant {
											id
											variantDescription
											referenceAllele
											alternateAllele
											rsIds
											__typename
											}
										}
										}
									}
									}"

			### Set variables object of arguments to be passed to endpoint

			id_variables = list(
								queryString = query_rsID,
								index = 0,
								entityNames = c("Variant")
								)

			### Construct POST request body object with query string and variables
			id_search_body = list(query = id_query_string, variables = id_variables)

			### Perform OpenTargets search request
			id_search_out = POST(url = base_url, body = toJSON(id_search_body, auto_unbox = TRUE),
									add_headers("Content-Type" = "application/json"),
									encode = "raw")

			if (is.null(content(id_search_out)$data$search$total)) {

				cat(paste0("Results are NULL for ", basename(all_full_jma[i])))

				next

			}

			### If content of id_search_out has a length of 0 - the rsID could not be found on OpenTargets
			### Skip to next SNP

			if (content(id_search_out)$data$search$total == 0){

				query_result = data.frame(hg38_snpid = NA,
											l2g_gene = NA,
											l2g_score = NA,
											nearest_gene = NA,
											severest_consequence = NA)

				query_frame = rbind(query_frame, query_result)

				variantID = 0

				next

			}

			### Extract the hits list

			query_hits = content(id_search_out)$data$search$hits

			snp_table = as.data.frame(matrix(ncol = 2, nrow = 0))
			colnames(snp_table) = c("rsid", "snpid")

			### Extraction from OTG also extracts similarly named SNPs
			### Convert the GraphQL data from OTG into a data frame and then extract the CORRECT SNP ID

			for (nsnp in 1:length(query_hits)){

				if (is.null(query_hits[[nsnp]]$object$rsIds)){

					next

				}

				newrow = data.frame(rsid = query_hits[[nsnp]]$object$rsIds[[1]],
									snpid = query_hits[[nsnp]]$object$id)

				snp_table = rbind(snp_table, newrow)

			}

			### If there are SNPs in SNP table - extract the SNP ID
			### Otherwise, skip to NEXT SNP
			
			if (query_rsID %in% snp_table$rsid){

				### Variant ID output

				variantID = snp_table$snpid[snp_table$rsid == query_rsID]

			} else {

				query_result = data.frame(hg38_snpid = NA,
											l2g_gene = NA,
											l2g_score = NA,
											nearest_gene = NA,
											severest_consequence = NA)

				query_frame = rbind(query_frame, query_result)

				variantID = 0

				next

			}

			### With variant ID defined, now run L2G

			v2g_query_string = "query($variantId: String!) {
									variant(variantId: $variantId) {
										id
									    transcriptConsequences {
									    	distanceFromFootprint
									    	target {
									        	id
									        	approvedSymbol       
									    	}
									    }
									    mostSevereConsequence {
									    	label
									    }
									    credibleSets(studyTypes: [gwas]) {
									    	rows {
									    		l2GPredictions {
									    			rows {
									    				score
									    				target {
									    					id
									    					approvedSymbol
									    				}
													}
												}
											}
										}
									}
								}"
					
			### Set variables object of arguments to be passed to V2G

			v2g_variables = list("variantId" = variantID)

			### Construct POST request body object with query string and variables

			v2g_post_body = list(query = v2g_query_string, variables = v2g_variables)

			# Perform POST request

			currSNP_v2g = POST(url = base_url, body = v2g_post_body, encode = "json")

			### As before, if the content of the data is empty it can be skipped
			### Otherwise we now extract the required data
			### Done piece by piece

			query_result = data.frame(hg38_snpid = variantID,
											l2g_gene = NA,
											l2g_score = NA,
											nearest_gene = NA,
											severest_consequence = NA
										)

			### Extract closest gene based on distance from footprint

			if (length(content(currSNP_v2g)$data$variant$transcriptConsequences) != 0){

				closest_gene_table = as.data.frame(matrix(ncol = 2, nrow = 0))
				colnames(closest_gene_table) = c("gene", "distance")

				for (nconseq in 1:length(content(currSNP_v2g)$data$variant$transcriptConsequences)){

					newrow = data.frame(gene = content(currSNP_v2g)$data$variant$transcriptConsequences[[nconseq]]$target$approvedSymbol,
										distance = content(currSNP_v2g)$data$variant$transcriptConsequences[[nconseq]]$distanceFromFootprint)

					closest_gene_table = rbind(closest_gene_table, newrow)

				}

				closest_gene = closest_gene_table$gene[which.min(closest_gene_table$distance)]

				query_result$nearest_gene = closest_gene

			}

			### Extract most severe consequence

			if (length(content(currSNP_v2g)$data$variant$mostSevereConsequence$label) != 0){

				most_severe_consequence = content(currSNP_v2g)$data$variant$mostSevereConsequence$label

				query_result$severest_consequence = most_severe_consequence

			}

			### Extract locus-2-gene scores and genes
			### Previously used V2G

			if (length(content(currSNP_v2g)$data$variant$credibleSets$rows) != 0){

				l2g_gene_table = as.data.frame(matrix(ncol = 2, nrow = 0))
				colnames(l2g_gene_table) = c("gene", "l2g_score")

				for (ngene in 1:length(content(currSNP_v2g)$data$variant$credibleSets$rows)){

					if (length(content(currSNP_v2g)$data$variant$credibleSets$rows[[ngene]]$l2GPredictions$rows) != 0){

						for (kn_gene in 1:length(content(currSNP_v2g)$data$variant$credibleSets$rows[[ngene]]$l2GPredictions$rows)){

							newrow = data.frame(gene = content(currSNP_v2g)$data$variant$credibleSets$rows[[ngene]]$l2GPredictions$rows[[kn_gene]]$target$approvedSymbol,
												l2g_score = content(currSNP_v2g)$data$variant$credibleSets$rows[[ngene]]$l2GPredictions$rows[[kn_gene]]$score)

							l2g_gene_table = rbind(l2g_gene_table, newrow)

						}
					}
				}

				if (nrow(l2g_gene_table) > 0){

					highest_l2g = which.max(l2g_gene_table$l2g_score)

					query_result$l2g_gene = l2g_gene_table$gene[highest_l2g]
					query_result$l2g_score = l2g_gene_table$l2g_score[highest_l2g]

				}
			}

			### Combine into overall data frame

			query_frame = rbind(query_frame, query_result)
		}

		current_analysis_v2g = cbind(current_analysis, query_frame)

		fwrite(current_analysis_v2g, gsub(".txt", "_annotated.txt", all_full_jma[i]), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

	}
}



##### =========================== #####

### PheWAS and associations for each SNP via OTG

##### =========================== #####

for (i in 1:length(all_full_jma)){

	current_analysis = fread(all_full_jma[i], data.table = FALSE)

	heading(paste0("Currently running PheWAS on ", basename(all_full_jma[i])))

	current_analysis_phewas = data.frame(matrix(ncol = 10, nrow = 0))
	colnames(current_analysis_phewas) = c("studyID", "beta", "pval", "nTotal", "nCases", "trait", "traitCategory", "or", "rsid", "snpid")

	for (j in 1:nrow(current_analysis)){

		cat(paste0("Currently running on SNP ", j, " of ", nrow(current_analysis), "\n\n\n"))

		query_rsID = current_analysis$SNP[j]

		### Get OpenTargets SNP ID

		id_query_string = "query useSearchToConvertRSIDIntoIDFormat($query_rsID: String!) {
							search(queryString: $query_rsID) {
								variants {
									id
									rsId
									nearestGene {
										id
										start
										symbol
										tss
										description
										chromosome
										exons
									}
									nearestGeneDistance
								}
							}
						}"

		### Set variables object of arguments to be passed to endpoint
		id_variables = list("query_rsID" = query_rsID)

		### Construct POST request body object with query string and variables
		id_search_body = list(query = id_query_string, variables = id_variables)

		### Perform OpenTargets search request
		id_search_out = POST(url = base_url, body = id_search_body, encode = "json")

		if (length(content(id_search_out)$data$search$variants) == 0){

			query_result = data.frame(studyID = NA,
											beta = NA,
											pval = NA,
											nTotal = NA,
											nCases = NA,
											trait = NA,
											traitCategory = NA,
											or = NA,
											rsid = query_rsID,
											snpid = NA)

			current_analysis_phewas = rbind(current_analysis_phewas, query_result)

			variantID = 0

			next

		}

		### Variant ID output
		variantID = content(id_search_out)$data$search$variants[[1]]$id

		v2d_query_string = "query pheWASsearch($variantId: String!) {
								pheWAS(variantId: $variantId) {
									associations{
										studyId
										eaf
										beta
										pval
										nTotal
										nCases
										study {
											traitReported
											traitCategory
											pmid
										}
										oddsRatio
									}
								}
							}"
				
		### Set variables object of arguments to be passed to V2G
		v2d_variables = list("variantId" = variantID)

		### Construct POST request body object with query string and variables
		v2d_post_body = list(query = v2d_query_string, variables = v2d_variables)

		# Perform POST request
		currSNP_v2d = POST(url = base_url, body = v2d_post_body, encode = 'json')

		if (length(content(currSNP_v2d)$data$pheWAS$associations) > 0){

			studyID = content(currSNP_v2d)$data$pheWAS$associations[[1]]$studyId
			beta = content(currSNP_v2d)$data$pheWAS$associations[[1]]$beta
			pval = content(currSNP_v2d)$data$pheWAS$associations[[1]]$pval
			nTotal = content(currSNP_v2d)$data$pheWAS$associations[[1]]$nTotal
			nCases = ifelse(is.null(content(currSNP_v2d)$data$pheWAS$associations[[1]]$nCases), NA, content(currSNP_v2d)$data$pheWAS$associations[[1]]$nCases)
			trait = content(currSNP_v2d)$data$pheWAS$associations[[1]]$study$traitReported
			or = ifelse(is.null(content(currSNP_v2d)$data$pheWAS$associations[[1]]$oddsRatio), NA, content(currSNP_v2d)$data$pheWAS$associations[[1]]$oddsRatio)
			traitCategory = content(currSNP_v2d)$data$pheWAS$associations[[1]]$study$traitCategory

			currSNP_v2d_frame = data.frame(studyID = studyID,
											beta = beta,
											pval = pval,
											nTotal = nTotal,
											nCases = nCases,
											trait = trait,
											traitCategory = traitCategory,
											or = or,
											rsid = query_rsID,
											snpid = variantID
											)

			if (length(content(currSNP_v2d)$data$pheWAS$associations) > 1){

				for (ntrait in 2:length(content(currSNP_v2d)$data$pheWAS$associations)){

					studyID = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$studyId
					beta = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$beta
					pval = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$pval
					nTotal = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$nTotal
					nCases = ifelse(is.null(content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$nCases), NA, content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$nCases)
					trait = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$study$traitReported
					traitCategory = content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$study$traitCategory
					or = ifelse(is.null(content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$oddsRatio), NA, content(currSNP_v2d)$data$pheWAS$associations[[ntrait]]$oddsRatio)

					currRow = data.frame(studyID = studyID,
											beta = beta,
											pval = pval,
											nTotal = nTotal,
											nCases = nCases,
											trait = trait,
											traitCategory = traitCategory,
											or = or,
											rsid = query_rsID,
											snpid = variantID
											)

					currSNP_v2d_frame = rbind(currSNP_v2d_frame, currRow)

				}
			}
		} else {

			query_result = data.frame(studyID = NA,
											beta = NA,
											pval = NA,
											nTotal = NA,
											nCases = NA,
											trait = NA,
											traitCategory = NA,
											or = NA,
											rsid = query_rsID,
											snpid = variantID)

			current_analysis_phewas = rbind(current_analysis_phewas, query_result)

			variantID = 0

			next

		}

		current_analysis_phewas = rbind(current_analysis_phewas, currSNP_v2d_frame)
	
	}

	fwrite(current_analysis_phewas, gsub(".txt", "_otg_phewas.txt", all_full_jma[i]), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")

}
