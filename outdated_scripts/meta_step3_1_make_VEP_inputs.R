##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(stringr)
options(scipen = 999)

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/METAL/")
outDir = paste0(mainDir, "processing/misc_data/VEP_inputs/")

sumstat_list = Sys.glob(paste0(metaDir, "*TBL.gz"))


##### =========================== #####

### Prune for GCTA

##### =========================== #####

### The following is run once in order to create a master list of SNPs from the meta
### Used in conjunction with Goldenpath to create chr_pos:rsid match

output_snps = fread(sumstat_list[1], data.table = FALSE, select = c("MarkerName", "Allele1", "Allele2", "Freq1"))

output_snps = output_snps[output_snps$Freq1 > 0.001 & output_snps$Freq1 < 0.999,]
output_snps$fullid = paste(output_snps$MarkerName, output_snps$Allele1, output_snps$Allele2, sep = "_")
output_snps = data.frame(snpid = output_snps$fullid)

for (i in 2:length(sumstat_list)){

    cat(paste0("\nWorking on file ", i, " which is ", strsplit(basename(sumstat_list[i]), "_metal")[[1]][1], "\n"))

    curr_sumstats = fread(sumstat_list[i], data.table = FALSE, select = c("MarkerName", "Allele1", "Allele2", "Freq1"))

    curr_sumstats = curr_sumstats[curr_sumstats$Freq1 > 0.001 & curr_sumstats$Freq1 < 0.999,]
    curr_sumstats$fullid = paste(curr_sumstats$MarkerName, curr_sumstats$Allele1, curr_sumstats$Allele2, sep = "_")
    curr_sumstats = data.frame(snpid = curr_sumstats$fullid)

    cat(paste0("\nNrow of file ", i, " = ", nrow(curr_sumstats), "\n"))
    print(head(curr_sumstats))

    output_snps = rbind(output_snps, curr_sumstats)

    cat(paste0("\nNrow of output = ", nrow(output_snps), "\n"))

    output_snps = unique(output_snps)

    cat(paste0("\nNrow of output after unique = ", nrow(output_snps), "\n"))
    print(tail(output_snps))

}

ident_frame = str_split_fixed(output_snps$snpid, "_", 4)

ident_frame = ident_frame[order(ident_frame[,1], ident_frame[,2]),]
ident_frame[,1] = gsub("23", "X", ident_frame[,1])

division_index = ceiling(nrow(ident_frame)/1500000)

for (j in 1:division_index){

    if (j != division_index){

        start_row = ((j - 1) * 1500000) + 1
        end_row = j * 1500000

    } else {

        start_row = ((j - 1) * 1500000) + 1
        end_row = nrow(ident_frame)

    }

    curr_set_out = ident_frame[c(start_row:end_row),]

    output_snp_frame_forward = data.frame(chr = curr_set_out[,1],
                                    pos = curr_set_out[,2],
                                    pos1 = curr_set_out[,2],
                                    allele = toupper(paste0(curr_set_out[,3], "/", curr_set_out[,4])),
                                    dir = "+")

    fwrite(output_snp_frame_forward, paste0(outDir, "gencost_allSNPs_for_VEP_forward_", j, ".txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

    output_snp_frame_backward = data.frame(chr = curr_set_out[,1],
                                    pos = curr_set_out[,2],
                                    pos1 = curr_set_out[,2],
                                    allele = toupper(paste0(curr_set_out[,3], "/", curr_set_out[,4])),
                                    dir = "-")

    fwrite(output_snp_frame_backward, paste0(outDir, "gencost_allSNPs_for_VEP_backward_", j, ".txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

}

