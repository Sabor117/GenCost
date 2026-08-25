library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(httr2)
library(jsonlite)
library(purrr)
library(patchwork)
library(forcats)

### PGS files and locations

univar_prs_files = Sys.glob("Documents/Finland PostDoc/GenCOST/Data/FinnGen_PGS_analysis_data/finngen_R13_outdated/*_univar_pgs_analysis.txt")
multivar_prs_files = Sys.glob("Documents/Finland PostDoc/GenCOST/Data/FinnGen_PGS_analysis_data/finngen_R13/*multivar*")

### PGS meta-data

meta_data_file = "Documents/Finland PostDoc/GenCOST/Data/FinnGen_PGS_analysis_data/scores_metadata_filtered.tsv"
meta_data = fread(meta_data_file, data.table = FALSE)

### Pre-existing categories

category_file = "Documents/Finland PostDoc/GenCOST/Data/PGSC_categories_table.tsv"
original_categories = fread(category_file, data.table = FALSE)

##### =========================== #####

### Set up OpenAI API

##### =========================== #####

Sys.setenv(OPENAI_API_KEY = # THIS IS NOT FOR SHARING
)
### OpenAI API caller function
### By ChatGPT

openai_request = function(messages,
                           model = "gpt-4.1-mini",
                           max_tokens = 4096) {
  api_key <- Sys.getenv("OPENAI_API_KEY")
  if (api_key == "") stop("Set OPENAI_API_KEY in your environment first.")
  
  req <- request("https://api.openai.com/v1/chat/completions") |>
    req_auth_bearer_token(api_key) |>
    req_body_json(list(
      model = model,
      messages = messages,
      temperature = 0,
      max_tokens = max_tokens,
      response_format = list(type = "json_object")
    ))
  
  resp <- req_perform(req)
  content <- resp_body_json(resp, simplifyVector = FALSE)
  content
}


##### =========================== #####

### Functions

##### =========================== #####

### Function to try and parse the json more safely
### This did not work and is currently NOT IN USE

parse_json_safe = function(text) {
  # Strip leading/trailing whitespace just in case
  text <- trimws(text)
  
  # If model ever ignored response_format and added fences, strip them
  text <- sub("^```json\\s*", "", text)
  text <- sub("^```\\s*", "", text)
  text <- sub("```\\s*$", "", text)
  
  # Optional: keep only from the first '{' if something else is before JSON
  first_brace <- regexpr("\\{", text)
  if (first_brace > 1) {
    text <- substr(text, first_brace, nchar(text))
  }
  
  # Check if there are any characters that shouldn't be in valid JSON
  invalid_chars <- grep("[^[:print:][:space:],:{}[\\]]", text, value = TRUE)
  if (length(invalid_chars) > 0) {
    cat("Warning: Found potential invalid characters in JSON output.\n")
    cat(paste(invalid_chars, collapse = " "), "\n")
  }
  
  # Quick diagnostic if invalid
  if (!jsonlite::validate(text)) {
    cat("JSON appears invalid. First 300 chars:\n",
        substr(text, 1, 300), "\n...\n")
    stop("JSON from model is invalid / truncated.")
  }
  
  jsonlite::fromJSON(text, simplifyVector = TRUE)
}


##### =========================== #####

### Create category and canonical names table

##### =========================== #####

currfile = fread(univar_prs_files[3], data.table = FALSE)

### First remove PGS which include FinnGen samples

nonFinngen_pgs = meta_data$PGSID[!(meta_data$`FinnGen samples`)]
currfile = currfile[currfile$pgs %in% nonFinngen_pgs,]

### Add official phenotypes

currfile$phenotype = original_categories$`Reported trait`[match(currfile$pgs, original_categories$PGSID)]

### Remove duplicates by phenotype, keeping lowest P-value

currfile_filtered = currfile %>%
  group_by(phenotype) %>%
  filter(pval == min(pval)) %>%
  ungroup()

### Convert to a table

final_table = data.frame(pgs = currfile_filtered$pgs,
                          phenotype = currfile_filtered$phenotype)

### Now use API for the duplicates with the same names

phenos = unique(final_table$phenotype)
phenos = phenos[order(phenos, decreasing = FALSE)]


##### =========================== #####

### First API call - duplicates

##### =========================== #####

### API prompt

system_prompt = "You are a meticulous data curator for genetic epidemiology results. You will receive a list of phenotype names. You must group names that refer to the same underlying phenotype concept (e.g. 'BMI', 'Body mass index', 'Body mass index (BMI)'). Treat spelling variants, abbreviations, and synonyms as the same phenotype if a human genetic epidemiologist would reasonably consider them equivalent.

Return ONLY valid JSON - no explanations or additional text - ENSURE THIS IS VALID. Format:
{
  \"clusters\": [
    {
      \"canonical_name\": \"string\",
      \"aliases\": [\"string\", ...]
    },
    ...
  ]
}

Rules:
- 'canonical_name' should be the most clear and informative label (e.g. 'Body mass index').
- 'canonical_name' should NOT be an abbreviation
- 'canonical_name' should NOT contain multiple names (e.g. 'Body mass index (BMI)') and keep just one per entry (e.g. 'Body mass index')
- Every input name must appear in exactly one 'aliases' list.
- Do NOT invent phenotypes that were not in the input.
- Keep phenotypes distinct if they clearly represent different concepts (e.g. 'Systolic blood pressure' vs 'Diastolic blood pressure')."

### Split phenotypes into chunks of up to 200

chunk_size = 60
idx = split(seq_along(phenos), ceiling(seq_along(phenos) / chunk_size))

### For each chunk:
all_cluster_dfs = lapply(idx, function(idxs) {
  
  # 1) Phenotypes in this chunk
  phs = phenos[idxs]
  
  # 2) Build the list string for this chunk
  pheno_list_str = paste0(seq_along(phs), ". ", phs, collapse = "\n")
  
  # 3) Chunk-specific user prompt
  user_prompt = paste(
    "Here is the list of phenotype names, one per line with an index.",
    "Please group them as described.\n\n",
    pheno_list_str
  )
  
  # 4) Call API
  resp = openai_request(
    messages = list(
      list(role = "system", content = system_prompt),
      list(role = "user",   content = user_prompt)
    )
  )
  
  # 5) Extract JSON string and parse safely
  json_text = resp$choices[[1]]$message$content
  clusters = jsonlite::fromJSON(json_text)
  
  #clusters = parse_json_safe(json_text)
  
  # 6) 'clusters$clusters' is typically a data.frame with columns:
  #    - canonical_name (character)
  #    - aliases (list-column of character vectors)
  #
  # We turn this into a long data frame: one row per alias, with its canonical name.
  chunk_clusters = as.data.frame(clusters$clusters, stringsAsFactors = FALSE)
  
  # Make sure 'aliases' is a list-column (it should be)
  # Then unnest:
  chunk_df = chunk_clusters |>
    tidyr::unnest(aliases) |>
    dplyr::rename(
      phenotype = aliases,
      canonical_phenotype = canonical_name
    )
  
  chunk_df
})

### Combine all chunks into one mapping table

cluster_df = bind_rows(all_cluster_dfs)
cluster_df = unique(cluster_df)

### Merge the canonical_phenotype into the final_table by matching 'phenotype' with 'phenotype' in the cluster_df
final_table_merged = final_table |>
  left_join(cluster_df, by = "phenotype")

final_table_merged = unique(final_table_merged)
final_table_merged = final_table_merged[order(final_table_merged$canonical_phenotype, decreasing = FALSE),]


##### =========================== #####

### Second API call - categories

##### =========================== #####

### List of categories you want to group phenotypes into
### Adjust the number of categories as needed
category_list = c(
  "Anthropometric/body measurement",
  "Cardiovascular disease",
  "Metabolic disorder",
  "Digestive system disorder",
  "Immune system disorder",
  "Neurological disorder",
  "Lipid or lipoprotein measurement",
  "Cardiovascular measurement",
  "Hematological measurement",
  "Inflammatory measurement",
  "Liver enzyme measurement",
  "Cancer",
  "Other lab measurements",
  "Other disease",
  "Other trait"
)
### Make sure the system message is aligned with what the AI needs to do
system_prompt_cat = paste0("You are an expert in genetic epidemiology and PGS phenotypes.",
                           " You will receive a list of phenotype names.",
                           " For each phenotype provided, assign exactly ONE of the following categories:\n",
                           paste0("- ", category_list, collapse = "\n"), "\n\n",
                           "Return ONLY VALID JSON, no explanations. Format:\n",
                           "{\n  \"assignments\": [\n",
                           "    {\"phenotype\": \"...\", \"category\": \"...\"},\n",
                           "    ...\n  ]\n}\n\n",
                           "Rules:\n",
                           "- 'phenotype' must match exactly one of the input names.\n",
                           "- 'category' must be one of the allowed categories above.\n",
                           "- If uncertain, choose the most reasonable category.\n",
                           "- Only use 'Other trait' as last resort if no reasonable cartegory is present.\n")

### The phenotypes list
canonical_phenos = unique(final_table_merged$canonical_phenotype)

### Define chunk size (e.g., 50 phenotypes per chunk)
chunk_size = 60
idx_cat = split(seq_along(canonical_phenos), ceiling(seq_along(canonical_phenos) / chunk_size))

### Process each chunk separately
all_category_dfs = lapply(idx_cat, function(idxs) {
  # Select phenotypes for this chunk
  phs = canonical_phenos[idxs]
  
  # Build the string for the phenotype list in this chunk
  pheno_list_str_cat = paste0(seq_along(phs), ". ", phs, collapse = "\n")
  
  # Generate the user prompt for the chunk
  user_prompt_cat = paste(
    "Here is the list of phenotype names, one per line with an index.",
    "Please assign a category to each.\n\n",
    pheno_list_str_cat
  )
  
  # Send the API request for this chunk
  resp_cat <- openai_request(
    messages = list(
      list(role = "system", content = system_prompt_cat),
      list(role = "user",   content = user_prompt_cat)
    )
  )
  
  # Extract JSON string and parse safely
  json_text_cat <- resp_cat$choices[[1]]$message$content
  clusters_cat <- jsonlite::fromJSON(json_text_cat)
  
  # Convert this chunk's categories to a data frame
  chunk_cat_df <- as.data.frame(clusters_cat$assignments, stringsAsFactors = FALSE)
  
  chunk_cat_df
})

### Combine the results of all chunks into one final data frame
category_df = bind_rows(all_category_dfs)

### Final merging
final_table_merged = left_join(final_table_merged, category_df, by = c("canonical_phenotype" = "phenotype"))
final_table_merged = unique(final_table_merged)

fwrite(final_table_merged, "Documents/Finland PostDoc/GenCOST/Data/PGSC_categories_OpenAI_noFinnGen.tsv",
       row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")


##### =========================== #####

### Getting FINAL categories file

##### =========================== #####

final_table_merged = fread("Documents/Finland PostDoc/GenCOST/Data/PGSC_categories_OpenAI_noFinnGen.tsv", data.table = FALSE)
starter_inpat_file = fread(univar_prs_files[grepl("inpat", univar_prs_files)], data.table = FALSE)

### Read current file
### Select based on phenotype
  
inpatfile = fread(univar_prs_files[grepl("inpat", univar_prs_files)], data.table = FALSE)

### Restrict to non-FinnGen PGS
### Then add names and categories

inpatfile = inpatfile[inpatfile$pgs %in% nonFinngen_pgs,]
inpatfile$phenotype = original_categories$`Reported trait`[match(inpatfile$pgs, original_categories$PGSID)]

inpatfile = inpatfile %>%
  group_by(phenotype) %>%
  filter(pval == min(pval)) %>%
  ungroup()
  
inpatfile = left_join(inpatfile, final_table_merged[,c("pgs", "canonical_phenotype", "category")], by = "pgs")

### For duplicated phenotypes, keep row with smallest pval
  
inpatfile = inpatfile |>
    group_by(canonical_phenotype) |>
    slice_min(pval, with_ties = FALSE) |>
    ungroup()

final_pgs_file = inpatfile[,c("pgs", "canonical_phenotype", "category")]

fwrite(final_pgs_file, "Documents/Finland PostDoc/GenCOST/Data/PGSC_categories_OpenAI_noFinnGen_FINAL.tsv",
       row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")

final_pgs_file = fread("Documents/Finland PostDoc/GenCOST/Data/PGSC_categories_OpenAI_noFinnGen_FINAL.tsv", data.table = FALSE)


##### =========================== #####

### Being plotting sections

##### =========================== #####
  
### Add a column for significance
inpatfile$neglog10p = -log10(inpatfile$pval)
  
### If any ngelog10p are INF convert them to 310
inpatfile$neglog10p = ifelse(is.infinite(inpatfile$neglog10p), 310, inpatfile$neglog10p)
  
### Color-blind friendly palette
cb_palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#000000", "#E6AB02",
               "#7570B3", "#66A61E", "#E7298A", "#A6761D", "#666666",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")

pal_cat = setNames(cb_palette[seq_along(category_list)], category_list) # names = categories
pal_all = c(pal_cat, "Non-significant" = "grey80") # add non-sig grey

# Significance threshold
sig_threshold = -log10(0.05/nrow(inpatfile))


##### =========================== #####

### Manhattan-style plot

##### =========================== #####

p_manhat = ggplot(inpatfile, aes(x = category, y = neglog10p, color = category)) +
  geom_jitter(aes(alpha = ifelse(neglog10p > sig_threshold, 1, 0.3)), 
              width = 0.3, height = 0.3, size = 2) +  # Adjust jitter width/height
  geom_text_repel(data = top_hits, aes(label = canonical_phenotype), size = 3, max.overlaps = 20) +
  scale_color_manual(values = cb_palette) +  # Use the color palette
  geom_hline(yintercept = sig_threshold, color = "red", linetype = "dashed") + 
  theme_linedraw() +
  labs(x = "PGS Category", y = "-log10(P-value)", title = "Phenotype associations with inpatient cost") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = margin(l = 30),
        legend.position = "none")  # Remove legend


##### =========================== #####

### Creating INPAT + DRUG merged file for Z-score comparison plot

##### =========================== #####

drugfile = fread(univar_prs_files[grepl("kela", univar_prs_files)], data.table = FALSE)

drugfile = drugfile[drugfile$pgs %in% nonFinngen_pgs,]
drugfile$phenotype = original_categories$`Reported trait`[match(drugfile$pgs, original_categories$PGSID)]

drugfile = left_join(drugfile, final_table_merged[,c("pgs", "canonical_phenotype", "category")], by = "pgs")

### Add a column for significance
drugfile$neglog10p = -log10(drugfile$pval)

### If any ngelog10p are INF convert them to 310
drugfile$neglog10p = ifelse(is.infinite(drugfile$neglog10p), 310, drugfile$neglog10p)

inpatfile$zscore = inpatfile$est / inpatfile$ste
drugfile$zscore = drugfile$est / drugfile$ste

drugfile = drugfile[drugfile$pgs %in% inpatfile$pgs,]

merged_pgs_file = left_join(inpatfile, drugfile,
                            by = "pgs", suffix = c("_inpat", "_drug"))

merged_pgs_file$sig_col = (merged_pgs_file$neglog10p_drug > sig_threshold | merged_pgs_file$neglog10p_inpat > sig_threshold)

### Make a point-colour column: category for sig, "non-significant" for non-sig

merged_pgs_file$point_col = ifelse(merged_pgs_file$sig_col,
                                    as.character(merged_pgs_file$category_inpat),
                                    "Non-significant")

merged_pgs_file$z_diff = abs(merged_pgs_file$zscore_drug) - abs(merged_pgs_file$zscore_inpat)

top_hits_inpat = merged_pgs_file %>% filter(neglog10p_inpat > 250)  # adjust threshold
top_hits_drug = merged_pgs_file %>% filter(neglog10p_drug > 250)  # adjust threshold

top_hits_combined = merged_pgs_file[merged_pgs_file$pgs %in% unique(c(top_hits_inpat$pgs, top_hits_drug$pgs)),]

### Selected manually from top_hits_combined

selected_pgs = c("Abdominal pain",
                 "Age at first live birth",
                 "Age at first sexual intercourse",
                 "Asthma",
                 "Autoimmune disease",
                 "Body mass index",
                 "Cardiovascular disease",
                 "College education",
                 "Diabetes",
                 "Hypertension",
                 "Lifetime Major Depressive Disorder",
                 "Number of medications taken",
                 "Respiratory and ear-nose-throat diseases",
                 "Average weekly alcohol consumption",
                 
                 ### Adding some with drug prevalance
                 
                 "Average weekly alcohol consumption",
                 "Time spent outdoors in summer",
                 "AST to ALT ratio",
                 "HDL cholesterol",
                 "Lung function (FEV1)",
                 "Hand grip strength"
                 )

top_hits_plotting = merged_pgs_file[merged_pgs_file$canonical_phenotype_drug %in% selected_pgs,]


##### =========================== #####

### MANUAL LABELS fo Z-score comparison plot

##### =========================== #####

manual_labels = top_hits_plotting %>%
  dplyr::transmute(
    canonical_phenotype_drug,
    # the TRUE point coordinates (must exist in top_hits_combined)
    x_point = zscore_drug,
    y_point = zscore_inpat,
    category_inpat,
    # MANUAL label coordinates (EDIT THESE)
    label_x = x_point + 10,
    label_y = y_point + 10
  )

manual_labels$label_x[1] = -46
manual_labels$label_y[1] = -4

manual_labels$label_x[2] = 19
manual_labels$label_y[2] = 39

manual_labels$label_x[3] = -45
manual_labels$label_y[3] = -50

manual_labels$label_x[4] = -45
manual_labels$label_y[4] = -28

manual_labels$label_x[5] = 55
manual_labels$label_y[5] = 10

manual_labels$label_x[6] = 43
manual_labels$label_y[6] = 5.5

manual_labels$label_x[7] = -25
manual_labels$label_y[7] = 20

manual_labels$label_x[8] = 31
manual_labels$label_y[8] = 49

manual_labels$label_x[9] = 71
manual_labels$label_y[9] = 15

manual_labels$label_x[10] = -2
manual_labels$label_y[10] = -43

manual_labels$label_x[11] = 64
manual_labels$label_y[11] = 3

manual_labels$label_x[12] = -39
manual_labels$label_y[12] = -23

manual_labels$label_x[13] = -29
manual_labels$label_y[13] = 13

manual_labels$label_x[14] = 75
manual_labels$label_y[14] = 60

manual_labels$label_x[15] = 45
manual_labels$label_y[15] = 55

manual_labels$label_x[16] = -50
manual_labels$label_y[16] = 7.5

manual_labels$label_x[17] = 67
manual_labels$label_y[17] = 43

manual_labels$label_x[18] = 50
manual_labels$label_y[18] = -5

manual_labels$label_x[19] = -5
manual_labels$label_y[19] = 35

manual_labels$label_x_line = ifelse(manual_labels$label_x < 0, manual_labels$label_x + 0.7, manual_labels$label_x - 0.9)
manual_labels$label_y_line = ifelse(manual_labels$label_y < 0, manual_labels$label_y + 0.7, manual_labels$label_y - 0.9)

manual_labels$label_x_line[1] = -36.7
manual_labels$label_y_line[1] = -4

manual_labels$label_y_line[4] = -28.7

manual_labels$label_y_line[9] = 15.3

manual_labels$label_y_line[11] = 3.3


# Plot the Z-scores
p_scatter = ggplot(merged_pgs_file, aes(x = zscore_drug, y = zscore_inpat)) +
  # points: use point_col for colour + your existing alpha mapping
  geom_point(aes(color = point_col,
                 alpha = ifelse(sig_col, 1, 0.3)),
             size = 3) +
  
  # segments from label -> true point
  geom_segment(
    data = manual_labels,
    aes(x = label_x_line,
        y = label_y_line,
        xend = x_point, yend = y_point),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  
  # label text placed exactly where you want
  geom_text(
    data = manual_labels,
    aes(x = label_x, y = label_y, label = canonical_phenotype_drug, color = category_inpat),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # palette: your cb_palette + add NS grey
  scale_color_manual(
    values = pal_all,
    breaks = names(pal_all), # keep legend to categories only (hide NS)
    name = "PGS category"
  ) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  labs(x = "Z-score in FinnGen (prescription drug costs)",
       y = "Z-score in FinnGen (inpatient costs)") +
  guides(alpha = "none") +
  theme_linedraw() +
  theme(
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.margin  = margin(l = 20)
  ) +
  xlim(-60, 80) +
  ylim(-60, 80) +
  coord_fixed(ratio = 1)

p_combo = p_manhat + p_scatter

ggsave("Documents/Finland PostDoc/GenCOST/Figures/misc/pgs_example_plot.png",
       p_combo, width = 24, height = 12, dpi = 300)


##### =========================== #####

### Create supplementary table

##### =========================== #####

inoutfile = fread(univar_prs_files[grepl("^hilmo", basename(univar_prs_files))], data.table = FALSE)
primfile = fread(univar_prs_files[grepl("avohilmo", basename(univar_prs_files))], data.table = FALSE)

inoutfile = inoutfile[inoutfile$pgs %in% nonFinngen_pgs,]
primfile = primfile[primfile$pgs %in% nonFinngen_pgs,]

inoutfile$phenotype = original_categories$`Reported trait`[match(inoutfile$pgs, original_categories$PGSID)]
primfile$phenotype = original_categories$`Reported trait`[match(primfile$pgs, original_categories$PGSID)]

inoutfile = left_join(inoutfile, final_table_merged, by = "pgs")
primfile = left_join(primfile, final_table_merged, by = "pgs")

inoutfile$zscore = inoutfile$est / inoutfile$ste
primfile$zscore = primfile$est / primfile$ste

inoutfile = inoutfile[inoutfile$pgs %in% inpatfile$pgs,]
primfile = primfile[primfile$pgs %in% inpatfile$pgs,]

inpatfile$query = "Inpatient"
drugfile$query = "Prescription drugs"
inoutfile$query = "Inpatient + outpatient"
primfile$query = "Primary care"

output_colset = c("pgs", "query", "est", "ste", "zscore", "pval", "rsq", "canonical_phenotype", "category")

output_supp = rbind(inpatfile[,output_colset], drugfile[,output_colset])
output_supp = rbind(output_supp, inoutfile[,output_colset])
output_supp = rbind(output_supp, primfile[,output_colset])

fwrite(output_supp, "Documents/Finland PostDoc/GenCOST/Data/supps_finngen_pgs_res_all_phenos.tsv",
       quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)


##### =========================== #####

### Create bar plots table

##### =========================== #####

inpat_set = inpatfile %>% slice_max(abs(zscore), n = 20) %>% pull(canonical_phenotype)
drug_set = drugfile %>% slice_max(abs(zscore), n = 20) %>% pull(canonical_phenotype)

inpat_idx = c(1,4,5,6,9,10,11,13,15,18)
drug_idx = c(1,2,3,4,5,9,11,14,15,18)

drug_set = drug_set[drug_idx]
inpat_set = inpat_set[inpat_idx]

bar_plot_inpat = inpatfile[inpatfile$canonical_phenotype %in% inpat_set,]
bar_plot_drug = drugfile[drugfile$canonical_phenotype %in% drug_set,]

bar_plot_inpat$est = ifelse(bar_plot_inpat$est < 0, bar_plot_inpat$est * -1, bar_plot_inpat$est)

bar_plot_inpat$canonical_phenotype[bar_plot_inpat$canonical_phenotype == "College education"] = "Less college education"
bar_plot_inpat$canonical_phenotype[bar_plot_inpat$canonical_phenotype == "Smoking status"] = "Smoking"
bar_plot_inpat$canonical_phenotype[bar_plot_inpat$canonical_phenotype == "Age at first live birth"] = "Earlier age at first\nlive birth"
bar_plot_inpat$canonical_phenotype[bar_plot_inpat$canonical_phenotype == "Age at first sexual intercourse"] = "Earlier age first had\nsexual intercourse"

bar_plot_drug$canonical_phenotype[bar_plot_drug$canonical_phenotype == "Vascular/heart problems diagnosed by doctor - High blood pressure"] = "High blood pressure\ndiagnosed by doctor"

bar_plot_inpat = bar_plot_inpat %>% mutate(canonical_phenotype = fct_reorder(canonical_phenotype, est, .desc = TRUE))
bar_plot_drug = bar_plot_drug %>% mutate(canonical_phenotype = fct_reorder(canonical_phenotype, est, .desc = TRUE))

bar_plot_inpat$est = (exp(bar_plot_inpat$est) - 1) * 100
bar_plot_inpat$ste = (exp(bar_plot_inpat$ste) - 1) * 100

bar_plot_drug$est = (exp(bar_plot_drug$est) - 1) * 100
bar_plot_drug$ste = (exp(bar_plot_drug$ste) - 1) * 100

## 3) First plot: inpatient
inpat_bar = ggplot(bar_plot_inpat, aes(x = canonical_phenotype, y = est, fill = category)) +
  geom_col() +
  geom_errorbar(aes(ymin = est - ste, ymax = est + ste), width = 0.2) +
  scale_fill_manual(values = pal_all, guide = "none") +  # no legend
  labs(
    x = NULL,
    y = "% change in FinnGen healthcare cost per 1 SD of PGS",
    title = "Inpatient"
  ) +
  ylim(0, 20) +
  theme_linedraw() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 60, vjust = 1, hjust = 1),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none",
    plot.background = element_rect(fill = "#C8E0D7", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(linewidth = 1),
    axis.ticks = element_line(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.8)
  )

## 4) Second plot: drugs
drug_bar = ggplot(bar_plot_drug, aes(x = canonical_phenotype, y = est, fill = category)) +
  geom_col() +
  geom_errorbar(aes(ymin = est - ste, ymax = est + ste), width = 0.2) +
  scale_fill_manual(values = pal_all, guide = "none") +
  labs(
    x = "PGS phenotype",
    y = NULL,
    title = "Prescription drugs"
  ) +
  ylim(0, 20) +
  theme_linedraw() +
  theme(axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 60, vjust = 1, hjust = 1), # size 12 for paper, 20 + bold for poster
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none",
    plot.background = element_rect(fill = "#C8E0D7", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(linewidth = 1),
    axis.ticks = element_line(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.8)
  )

pgs_plot_set = (inpat_bar / drug_bar) | p_scatter

outDir = "Documents/Finland PostDoc/GenCOST/Manuscript/Figures/"

ggsave(filename = paste0(outDir, "GenCOST_Fig_7_updated.png"), plot = pgs_plot_set,
       width = 20, height = 20, units = "in", dpi = 300)

pgs_plot_poster = inpat_bar / drug_bar

ggsave(filename = paste0(outDir, "GenCOST_Fig_7_POSTER.png"), plot = pgs_plot_poster,
       width = 20, height = 20, units = "in", dpi = 300)


##### =========================== #####

### Create INOUT + PRIM supplementary figure

##### =========================== #####

inoutfile = fread(univar_prs_files[grepl("^hilmo", basename(univar_prs_files))], data.table = FALSE)
primfile = fread(univar_prs_files[grepl("avohilmo", basename(univar_prs_files))], data.table = FALSE)

inoutfile$phenotype = original_categories$`Reported trait`[match(inoutfile$pgs, original_categories$PGSID)]
primfile$phenotype = original_categories$`Reported trait`[match(primfile$pgs, original_categories$PGSID)]

### Add a column for significance
inoutfile$neglog10p = -log10(inoutfile$pval)
primfile$neglog10p = -log10(primfile$pval)

### If any ngelog10p are INF convert them to 310
inoutfile$neglog10p = ifelse(is.infinite(inoutfile$neglog10p), 310, inoutfile$neglog10p)
primfile$neglog10p = ifelse(is.infinite(primfile$neglog10p), 310, primfile$neglog10p)

### Color-blind friendly palette
cb_palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#000000", "#E6AB02",
               "#7570B3", "#66A61E", "#E7298A", "#A6761D", "#666666",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")

pal_cat = setNames(cb_palette[seq_along(category_list)], category_list) # names = categories
pal_all = c(pal_cat, "Non-significant" = "grey80") # add non-sig grey

# Significance threshold
sig_threshold = -log10(0.05/nrow(inpatfile))

# Adding Z-scores

inoutfile$zscore = inoutfile$est / inoutfile$ste
primfile$zscore = primfile$est / primfile$ste

# Restrict to PGS in INPAT

inoutfile = inoutfile[inoutfile$pgs %in% inpatfile$pgs,]
primfile = primfile[primfile$pgs %in% inpatfile$pgs,]

merged_pgs_file = left_join(inoutfile, primfile,
                            by = "pgs", suffix = c("_inout", "_prim"))

merged_pgs_file$sig_col = (merged_pgs_file$neglog10p_prim > sig_threshold | merged_pgs_file$neglog10p_inout > sig_threshold)

### Get categories

merged_pgs_file = left_join(merged_pgs_file, final_table_merged[,c("pgs", "category")], by = "pgs")

### Make a point-colour column: category for sig, "non-significant" for non-sig

merged_pgs_file$point_col = ifelse(merged_pgs_file$sig_col,
                                   as.character(merged_pgs_file$category),
                                   "Non-significant")

merged_pgs_file$z_diff = abs(merged_pgs_file$zscore_prim) - abs(merged_pgs_file$zscore_inout)

# Plot the Z-scores
p_scatter_inout = ggplot(merged_pgs_file, aes(x = zscore_inout, y = zscore_prim)) +
  # points: use point_col for colour + your existing alpha mapping
  geom_point(aes(color = point_col,
                 alpha = ifelse(sig_col, 1, 0.3)),
             size = 3) +
  
  # palette: your cb_palette + add NS grey
  scale_color_manual(
    values = pal_all,
    breaks = names(pal_all), # keep legend to categories only (hide NS)
    name = "PGS category"
  ) +
  
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  labs(x = "Z-score in FinnGen\n(inpatient + outpatient costs)",
       y = "Z-score in FinnGen\n(primary care costs)") +
  guides(alpha = "none") +
  theme_linedraw() +
  theme(
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.margin  = margin(l = 60)
  ) +
  xlim(-60, 80) +
  ylim(-60, 80) +
  coord_fixed(ratio = 1)

ggsave(filename = paste0(outDir, "GenCOST_Supp_Fig7.png"), plot = p_scatter_inout,
       width = 12, height = 12, units = "in", dpi = 300)

