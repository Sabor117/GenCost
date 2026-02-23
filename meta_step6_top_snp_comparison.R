##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
library(dplyr)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
gctaDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")
metaDir = paste0(mainDir, "outputs/METAL_v4/")
tempDir = paste0(mainDir, "tmpdir/")
outDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/top_snp_matrix/")

all_gcta_files = Sys.glob(paste0(gctaDir, "*_out_annotated.txt"))
all_gcta_files = all_gcta_files[!(grepl("_no", all_gcta_files))]

all_meta_files = Sys.glob(paste0(metaDir, "*TBL.gz"))
all_meta_files = all_meta_files[!(grepl("_no", all_meta_files))]
all_meta_files = all_meta_files[!(grepl("tar", all_meta_files))]

phenotype_list = unique(str_split_fixed(basename(all_gcta_files), "_", 5)[,1])

main_analyses_files = all_gcta_files[grepl("_ALL_", all_gcta_files)]

start_time = Sys.time()


##### =========================== #####

### Files and 

##### =========================== #####

for (nanalysis in 1:length(phenotype_list)){

    analysis_time = Sys.time()

    curr_pheno = phenotype_list[nanalysis]

    curr_gcta_file = main_analyses_files[grepl(paste0(curr_pheno, "_"), main_analyses_files)]
    curr_top_snps = fread(curr_gcta_file, data.table = FALSE)

    all_snps = unique(curr_top_snps$SNP)

    all_other_files = all_meta_files[!(grepl(paste0(curr_pheno, "_ALL"), all_meta_files))]

    curr_output = data.frame(rsid = curr_top_snps$SNP,
                                chr = curr_top_snps$Chr,
                                pos = curr_top_snps$hg38_pos,
                                beta_ALL = curr_top_snps$b,
                                se_ALL = curr_top_snps$se,
                                p_ALL = curr_top_snps$p
                                )

    ### From previous iteration of code

    manual_replace = FALSE

    if (isTRUE(manual_replace)){

        ### Certain positions (hg38) seem to be missing (what the fuck?)

        curr_output$pos[curr_output$rsid == "rs1064173"] = 32659703
        curr_output$pos[curr_output$rsid == "rs776564339"] = 29045377
        curr_output$pos[curr_output$rsid == "rs202075303"] = 32528549
        curr_output$pos[curr_output$rsid == "rs28366334"] = 32594118
        curr_output$pos[curr_output$rsid == "rs28366358"] = 32597387
        curr_output$pos[curr_output$rsid == "rs2760980"] = 32597424
        curr_output$pos[curr_output$rsid == "rs9274447"] = 32665505
        curr_output$pos[curr_output$rsid == "rs28746841"] = 32666317
        curr_output$pos[curr_output$rsid == "rs199826652"] = 32597424
        curr_output$pos[curr_output$rsid == "rs34017414"] = 32586308
        curr_output$pos[curr_output$rsid == "rs9272327"] = 32636411
        curr_output$pos[curr_output$rsid == "rs9272353"] = 32636679
        curr_output$pos[curr_output$rsid == "rs28746823"] = 32665596
        curr_output$pos[curr_output$rsid == "rs17211937"] = 32669531
        curr_output$pos[curr_output$rsid == "rs4380799"] = 32604087
        curr_output$pos[curr_output$rsid == "rs144225452"] = 32671243
        curr_output$pos[curr_output$rsid == "rs529970560"] = 32468979
        curr_output$pos[curr_output$rsid == "rs28366215"] = 32589533
        curr_output$pos[curr_output$rsid == "rs9270949"] = 32605198
        curr_output$pos[curr_output$rsid == "rs149490268"] = 32608942

    }
    
    cat(paste0("\nCheck for missing positions:\n\n"))

    print(curr_output[curr_output$pos == "",])

    cat(paste0("\n==========================\n\n"))

    ### Loop reading other results

    for (nfile in 1:length(all_other_files)){

        file_time = Sys.time()

        cat(paste0("\nChecking FILE #", nfile, " of ", length(all_other_files), " which is: ", basename(all_other_files[nfile]), "\n\n===============\n\n"))

        curr_file = fread(all_other_files[nfile], data.table = FALSE)

        #curr_file_cols = fread(all_other_files[nfile], data.table = FALSE, nrow = 5)
        #colnames(curr_file) = colnames(curr_file_cols)

        curr_comparison = strsplit(basename(all_other_files[nfile]), "_metal")[[1]][1]

        snp_table = curr_file[curr_file$MarkerName %in% all_snps, c("MarkerName", "Chromosome", "Position", "Zscore", "P-value", "Weight", "Freq1")]

        snp_table$beta1 = snp_table$Zscore / sqrt((2 * snp_table$Freq1) * (1 - snp_table$Freq1) * (snp_table$Weight + (snp_table$Zscore^2)))
        snp_table$se = 1 / sqrt((2 * snp_table$Freq1) * (1 - snp_table$Freq1) * (snp_table$Weight + (snp_table$Zscore^2))) 

        snp_table = snp_table[,c("MarkerName", "Chromosome", "Position", "beta1", "se", "P-value")]

        snp_table$LOCUS = NA
        snp_table$LOCUS_SNP = NA

        curr_colnames = c("rsid", "chr", "pos",
                            paste0("beta_", curr_comparison),
                            paste0("se_", curr_comparison),
                            paste0("p_", curr_comparison),
                            paste0("locus_", curr_comparison),
                            paste0("locus_snp_", curr_comparison)
                            )

        colnames(snp_table) = curr_colnames

        if (any(!(curr_top_snps$SNP %in% snp_table$rsid))){

            missing_snps = data.frame(rsid = curr_top_snps$SNP[!(curr_top_snps$SNP %in% snp_table$rsid)],
                                    chr = curr_top_snps$Chr[!(curr_top_snps$SNP %in% snp_table$rsid)],
                                    pos = curr_top_snps$bp[!(curr_top_snps$SNP %in% snp_table$rsid)],
                                    V1 = NA,
                                    V2 = NA,
                                    V3 = NA,
                                    V4 = NA,
                                    V5 = NA)

            cat(paste0("\nSome SNPs not in overlapping analysis: ", nrow(missing_snps), "\n\n===============\n\n"))

            colnames(missing_snps) = curr_colnames

            snp_table = rbind(snp_table, missing_snps)

        }

        snp_start_time = Sys.time()

        for (nsnp in 1:nrow(snp_table)){

            currSnp = snp_table$rsid[nsnp]
            currPos = as.numeric(curr_output$pos[curr_output$rsid == currSnp])
            currChr = curr_output$chr[curr_output$rsid == currSnp]

            cat(paste0("\nChecking SNP #", nsnp, " of ", nrow(snp_table), ". Which is:\nrsid = ",
                        currSnp, "\nChrom = ",
                        currChr, "\nPos = ",
                        currPos, "\n---\n"))

            topPos = currPos + 1000000
            lowPos = currPos - 1000000

            curr_file_pruned = curr_file[curr_file$Chromosome == currChr,]
            curr_file_pruned = curr_file_pruned[curr_file_pruned$Position <= topPos & curr_file_pruned$Position >= lowPos,]

            if (any(curr_file_pruned[,"P-value"] < 5e-8)){

                snp_table[snp_table$rsid == currSnp, paste0("locus_", curr_comparison)] = TRUE
                snp_table[snp_table$rsid == currSnp, paste0("locus_snp_", curr_comparison)] = curr_file_pruned$MarkerName[which(curr_file_pruned[,"P-value"] == min(curr_file_pruned[,"P-value"]))][1]

                cat(paste0("\nLocus is significant in file ", basename(all_other_files[nfile]), "\nrsid = ",
                        curr_file_pruned$MarkerName[which(curr_file_pruned[,"P-value"] == min(curr_file_pruned[,"P-value"]))][1], "\nPval = ",
                        min(curr_file_pruned[,"P-value"]), "\n---\n"))

            } else {

                snp_table[snp_table$rsid == currSnp, paste0("locus_", curr_comparison)] = FALSE
                snp_table[snp_table$rsid == currSnp, paste0("locus_snp_", curr_comparison)] = NA

            }
        }

        snp_end_time = Sys.time()

        snp_time_run = snp_end_time - snp_start_time
        file_time_run = snp_end_time - file_time
        total_analysis_time = snp_end_time - analysis_time

        snp_table = snp_table[,-which(colnames(snp_table) %in% c("chr", "pos"))]
        curr_output = merge(curr_output, snp_table, by = "rsid")

        cat(paste0("\nOutput has been updated to include analysis: ", basename(all_other_files[nfile]),
                    "\n\nTime taken to search for loci = ", snp_time_run,
                    "\nTime taken to run this file = ", file_time_run,
                    "\nTotal time taken for current primary analysis = ", total_analysis_time,
                    "\n\n===============\n\n"))

        print(head(curr_output))
        print(summary(curr_output))

    }

    fwrite(curr_output, paste0(outDir, curr_pheno, "_top_snp_matrix_v4.txt"), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")

}


