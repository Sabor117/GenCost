aou_files = Sys.glob("*/*AOU*")

for (i in 1:length(aou_files)){

    currfile_in = aou_files[i]

    currfile_out = gsub("IN", "ALL", currfile_in)

    cat(paste0("cp ", currfile_in, " ", currfile_out, "\n\n\n"))

    system(paste0("cp ", currfile_in, " ", currfile_out))

}

gs20k_files = Sys.glob("*/*GS20K*")

for (i in 1:length(gs20k_files)){

    currfile_in = gs20k_files[i]

    currfile_out = gsub("ALL_", "IN_", currfile_in)
    currfile_out = gsub("RICHMOND.ALL", "RICHMOND.IN", currfile_out)

    cat(paste0("cp ", currfile_in, " ", currfile_out, "\n\n\n"))

    system(paste0("cp ", currfile_in, " ", currfile_out))

}

mgbb_files = Sys.glob("*/*MGBB*")

for (i in 1:length(mgbb_files)){

    currfile_in = mgbb_files[i]

    currfile_out = gsub("IN", "ALL", currfile_in)

    cat(paste0("cp ", currfile_in, " ", currfile_out, "\n\n\n"))

    system(paste0("cp ", currfile_in, " ", currfile_out))

}

qgp_files = Sys.glob("*/*QGP*INOUT*")

for (i in 1:length(qgp_files)){

    currfile_in = qgp_files[i]

    currfile_out = gsub("INOUT", "ALL", currfile_in)

    cat(paste0("cp ", currfile_in, " ", currfile_out, "\n\n\n"))

    system(paste0("cp ", currfile_in, " ", currfile_out))

}

