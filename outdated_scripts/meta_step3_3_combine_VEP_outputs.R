library(data.table)
library(stringr)

workingDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/VEP_outputs/"

all_forwards = Sys.glob(paste0(workingDir, "gencost_allSNPs_VEP_out_forward_*.txt"))
all_forwards = all_forwards[!(grepl("warn", all_forwards))]

all_backwards = Sys.glob(paste0(workingDir, "gencost_allSNPs_VEP_out_backward_*.txt"))
all_backwards = all_backwards[!(grepl("warn", all_backwards))]

limit_to_rsid = function(vep_frame){

    vep_frame = unique(vep_frame)

    multi_names = str_split_fixed(vep_frame[,2], ",", 2)

    vep_frame[,2] = multi_names[,1]

    vep_frame = vep_frame[grepl("rs", vep_frame[,2]),]

    vep_frame = unique(vep_frame)

    return(vep_frame)

}

forward_file_header = readLines(all_forwards[1], n = 300)
forward_file_start = grep("Uploaded", forward_file_header)[2] - 1

forward = fread(all_forwards[1], data.table = FALSE,
                skip = forward_file_start, select = c(1, 13))

forward = limit_to_rsid(forward)

overall_start_time = Sys.time()

for (i in 2:length(all_forwards)){

    function_start_time = Sys.time()

    cat(paste0("\nNow working on run ", i, " of ", length(all_forwards)))

    curr_file_header = readLines(all_forwards[i], n = 300)
    curr_file_header = grep("Uploaded", curr_file_header)[2] - 1

    curr_forwards = fread(all_forwards[i], data.table = FALSE,
                            skip = curr_file_header, select = c(1, 13))

    curr_forwards = limit_to_rsid(curr_forwards)

    forward = rbind(forward, curr_forwards)

    function_end_time = Sys.time()

    time_taken = function_end_time - function_start_time

    total_time = function_end_time - overall_start_time

    cat(paste0("\nRun ", i, " of ", length(all_forwards), " complete. Time taken: ", time_taken, " and total time taken: ", total_time, "\n\n\n"))

}

forward = unique(forward)

fwrite(forward, paste0(workingDir, "gencost_VEP_markerid_rsids_output_forward.txt"), row.names = FALSE, quote = FALSE, sep = "\t")



backward_file_header = readLines(all_backwards[1], n = 300)
backward_file_header = grep("Uploaded", file_header)[2] - 1

backward = fread(all_forwards[1], data.table = FALSE,
                skip = backward_file_header, select = c(1, 13))

backward = limit_to_rsid(backward)

overall_start_time = Sys.time()

for (i in 2:length(all_backwards)){

    function_start_time = Sys.time()

    cat(paste0("\nNow working on run ", i, " of ", length(all_backwards)))

    curr_file_header = readLines(all_backwards[i], n = 300)
    curr_file_header = grep("Uploaded", curr_file_header)[2] - 1


    curr_backwards = fread(all_backwards[i], data.table = FALSE,
                            skip = forward_file_start, select = c(1, 13))

    curr_backwards = limit_to_rsid(curr_backwards)

    backward = rbind(backward, curr_backwards)

    function_end_time = Sys.time()

    time_taken = function_end_time - function_start_time

    total_time = function_end_time - overall_start_time

    cat(paste0("\nRun ", i, " of ", length(all_backwards), " complete. Time taken: ", time_taken, " and total time taken: ", total_time, "\n\n\n"))

}

backward = unique(backward)

fwrite(backward, paste0(workingDir, "gencost_VEP_markerid_rsids_output_backward.txt"), row.names = FALSE, quote = FALSE, sep = "\t")




