##### =========================== #####

### Setup environment

##### =========================== #####

library(data.table)
library(locuszoomr, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(ensembldb, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(EnsDb.Hsapiens.v75, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(stringr)
library(dplyr)
library(AnnotationHub)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
snpDir = paste0(mainDir, "outputs/snp_forest_inputs_v4/")
figureDir = paste0(mainDir, "outputs/figures/snp_locuszoom_plots/")

### Main file lists

all_snp_files = Sys.glob(paste0(snpDir, "*_top_snp_metas.txt"))
all_meta_files = Sys.glob(paste0(mainDir, "outputs/METAL_v4/meta_v4_pruned/*.TBL.gz"))
all_meta_files = all_meta_files[!(grepl("_no", all_meta_files))]

### Analysis selection

file_select = 2 # 19 = IN_ALL, 12 = DRUG_ALL
analysis_name = strsplit(basename(all_snp_files[file_select]), "_top")[[1]][1]

### Getting AH database for LD linkage

ah = AnnotationHub(cache = "/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/annotation_hub_cache/")
ensDb_v111 = ah[["AH116291"]]

### Reading files

curr_snp_file = fread(all_snp_files[file_select], data.table = FALSE)
curr_meta_file = fread(all_meta_files[grepl(analysis_name, all_meta_files)], data.table = FALSE)

#rsid_translations = fread(paste0(mainDir, "processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out"), data.table = FALSE, header = FALSE)


##### =========================== #####

### Begin analysis

##### =========================== #####

### Remove low frequency alleles

curr_meta_file = curr_meta_file[curr_meta_file$Freq1 >= 0.001,]
curr_meta_file = curr_meta_file[curr_meta_file$Freq1 <= 0.999,]
colnames(curr_meta_file)[10] = "p"

curr_meta_file = curr_meta_file %>%
	filter((str_count(Direction, "\\-") + str_count(Direction, "\\+")) > 3)

### rsID translations has chrX - convert in meta

#curr_meta_file$Chromosome = gsub(23, "X", curr_meta_file$Chromosome)

### Get chromosome and position

#chr_pos = str_split_fixed(curr_meta_file$MarkerName, "_", 2)

#curr_meta_file$chr = chr_pos[,1]
#curr_meta_file$pos = as.numeric(chr_pos[,2])

### Add MarkerName (chr_pos) to SNP translations file

#rsid_translations$MarkerName = paste0(gsub("chr", "", rsid_translations$V1), "_", rsid_translations$V2)

### Keep only matches

#curr_meta_file = curr_meta_file[curr_meta_file$MarkerName %in% rsid_translations$MarkerName,]

### Merge rsID into meta file

#curr_meta_file = merge(curr_meta_file, rsid_translations[,c("V4", "MarkerName")], by = "MarkerName", all.x = TRUE)

#colnames(curr_meta_file)[18] = "rsid"

### Make LocusZoom plots

for (nsnp in 18:length(unique(curr_snp_file$snpid))){

    currrsid = unique(curr_snp_file$snpid)[nsnp]

    snpPos = as.numeric(unique(curr_snp_file$pos[curr_snp_file$snpid == currrsid]))
    snpChr = unique(curr_snp_file$chr[curr_snp_file$snpid == currrsid])

    plot_frame = curr_meta_file[curr_meta_file$Chromosome == snpChr,]
    plot_frame = plot_frame[plot_frame$Position <= snpPos + 6e5 & plot_frame$Position >= snpPos - 6e5,]

    plot_frame = plot_frame[grepl("rs", plot_frame$MarkerName),]

    analysis_locus = locus(index_snp = currrsid, data = plot_frame, flank = 5e5,
                            labs = "MarkerName",
                            chrom = "Chromosome",
                            pos = "Position",
                            p = "p",
                            ens_db = ensDb_v111,
                            LD = "r2")

    analysis_locus_ld = link_LD(analysis_locus, token = "2030cb076fd3")

    save_name = paste0(analysis_name, "_", currrsid, ".png")
    saveDir = paste0(figureDir, analysis_name, "/")

    system(paste0("mkdir -p ", saveDir))

    #analysis_locus_ld = link_recomb(analysis_locus_ld, genome = "hg38")

    png(paste0(saveDir, save_name), width = 3000, height = 2000, res = 300)

        locus_plot(analysis_locus_ld, labels = "index")

    dev.off()
  
}
