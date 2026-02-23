library(data.table)
library(stringr)

outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/processing/"
dataDir = "/scratch/project_2007428/data/base_data/ukbb_78537/phenotypes/"

translationFile = fread(paste0(dataDir, "ukbb_78537_translations.tsv"), data.table = FALSE)

### Subset of UKB phenotypes including basic phenotypes, covariates and costing data

ukb_costing_fields_list = c(33, # date of birth
                        31, # sex
                        21022, # age at recruitment
                        41253, # inpatient record format
                        54, # UKB assessment centre
                        40007, # age at death
                        34, # year of birth
                        52, # month of birth
                        53, # date of attending assessment centre
                        22001 # genetic sex
                        )

### Write field list

ukb_costing_fields_frame = data.frame(field = ukb_costing_fields_list)

fwrite(ukb_costing_fields_frame, paste0(outDir, "costing_fields_list.txt"),
        row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)