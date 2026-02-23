#####################

### Script set-up

#####################

cat("\n=============\n\nRscript starts. NOBODY expects the Seb Inquisition.\n\n===============\n")

library(data.table)

sessionInfo()

start_time = Sys.time()

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Wow. We're gonna meta-analyse this stuff? That's so fucking METAL...")

input_analysis_template_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/metal_pipeline/metal_analysis_template.txt"
input_script_template_file = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/metal_pipeline/metal_run_template.sh"

outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/scripts/metal_pipeline/scripts/"


##############################

### Making output scripts

##############################

meta_list = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/meta_sumstats/*")

input_analysis_template = readLines(input_analysis_template_file)
input_script_template = readLines(input_script_template_file)

analysis_name = "GenCOST"

for (i in 1:length(meta_list)){

    curr_input_analysis_template = input_analysis_template
    curr_input_script_template = input_script_template

    curr_meta = meta_list[i]

    meta_name = basename(curr_meta)

    curr_meta_files = Sys.glob(paste0(curr_meta, "/*txt.gz"))

    curr_meta_files = curr_meta_files[!(grepl("AOU", curr_meta_files))] # Remove AOU for meta v0.3+
    curr_meta_files = curr_meta_files[!(grepl("NTR", curr_meta_files))] # Remove NTR for meta v0.3+

    ### Base meta off UKB EUR and if UKB not present use FINNGEN

    if (any(grepl("UKB", curr_meta_files) & grepl("EUR", curr_meta_files))){

        matching_files = curr_meta_files[grepl("UKB", curr_meta_files) & grepl("EUR", curr_meta_files)]
        non_matching_files = curr_meta_files[!(grepl("UKB", curr_meta_files) & grepl("EUR", curr_meta_files))]

        curr_meta_files = c(matching_files, non_matching_files)

    } else if (any(grepl("FINNGEN", curr_meta_files))){

        curr_meta_files_FINN = curr_meta_files[grepl("FINNGEN", curr_meta_files)]

        curr_meta_files = curr_meta_files[-which(grepl("FINNGEN", curr_meta_files))]

        curr_meta_files = c(curr_meta_files_FINN, curr_meta_files)

    }

    curr_input_analysis_template = gsub("<CURRENT_TRAIT>", paste0(analysis_name, "_", meta_name), curr_input_analysis_template)
    curr_input_script_template = gsub("<CURRENT_TRAIT>", paste0(analysis_name, "_", meta_name), curr_input_script_template)

    curr_input_analysis_template = gsub("<OUT_PREFIX>", meta_name, curr_input_analysis_template)
    curr_input_script_template = gsub("<OUT_PREFIX>", meta_name, curr_input_script_template)

    process_lines = paste("PROCESS", curr_meta_files, sep = " ")

    for (nfile in 1:length(process_lines)){

        process_index = grep("<PROCESS_FILES>", curr_input_analysis_template)

        curr_input_analysis_template = gsub("<PROCESS_FILES>", process_lines[nfile], curr_input_analysis_template)

        n_lines_template = length(curr_input_analysis_template)

        if (nfile < length(process_lines)){

            curr_input_analysis_template = c(curr_input_analysis_template[1:process_index], "", "<PROCESS_FILES>", 
                                    curr_input_analysis_template[(process_index + 1) : n_lines_template])

        }
    }

    output_file = paste0(meta_name, "_metal_output_")
    metal_outDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/"

    curr_input_analysis_template = gsub("<OUTFILE_NAME>", paste0(metal_outDir, output_file), curr_input_analysis_template)

    analysis_file_location = paste0(meta_list[i], "/", meta_name, "_METAL_ANALYSIS.txt")

    writeLines(curr_input_analysis_template, analysis_file_location)

    curr_input_analysis_template = gsub("TRACKPOSITIONS ON", "#TRACKPOSITIONS OFF", curr_input_analysis_template)
    curr_input_analysis_template = gsub("CHROMOSOMELABEL", "#CHROMOSOMELABEL", curr_input_analysis_template)
    curr_input_analysis_template = gsub("POSITIONLABEL", "#POSITIONLABEL", curr_input_analysis_template)
    curr_input_analysis_template = gsub(paste0(meta_name, "_metal_output_"), paste0(meta_name, "_heritability_metal_output_"), curr_input_analysis_template)

    writeLines(curr_input_analysis_template, gsub("_METAL_ANALYSIS", "_METAL_ANALYSIS_heritability", analysis_file_location))

    curr_input_script_template = gsub("<CURRENT_FILE>", analysis_file_location, curr_input_script_template)
    curr_input_script_template = gsub("<CURRENT_FILE_HERIT>", gsub("_METAL_ANALYSIS", "_METAL_ANALYSIS_heritability", analysis_file_location), curr_input_script_template)

    curr_input_script_template = gsub("<CURRENT_OUTFILE>", paste0(metal_outDir, output_file, "1.TBL"), curr_input_script_template)
    curr_input_script_template = gsub("<CURRENT_OUTFILE_HERIT>", paste0(metal_outDir, gsub("_metal_output", "_heritability_metal_output", paste0(output_file, "1.TBL"))), curr_input_script_template)
    
    writeLines(curr_input_script_template, paste0(meta_list[i] , "/", meta_name, "_METAL_ANALYSIS_RUN.sh"))

}


