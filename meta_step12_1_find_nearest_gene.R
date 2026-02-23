## ---------------------------
##
## Script name: meta_step12_1_find_nearest_gene.R
##
## Purpose of script: More specifically find the nearest gene for given variant
##
## Author: Dr. Sebastian May-Wilson
## Contact: sebastian.may-wilson@helsinki.fi
##
## Date Created: 2025-09-02
##
## ---------------------------

##### =========================== #####

### Setup environment

##### =========================== #####

### Packages

library(SNPlocs.Hsapiens.dbSNP155.GRCh38, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(TxDb.Hsapiens.UCSC.hg38.knownGene, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(org.Hs.eg.db, lib.loc = "/projappl/project_2007428/RPackages_421/")
library(GenomicRanges)
library(data.table)
library(optparse)
library(dplyr)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
processDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")
gwasDir = paste0(mainDir, "outputs/METAL_v4/")

### Read from command line

option_list = list(

    make_option(c('--run_no', '-n'), help = "the file number from the list of GCTA files", default = 1),
    make_option(c('--working_dir', '-d'), help = "*output prefix", default = "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")

)

option.parser = OptionParser(option_list=option_list)
opt = parse_args(option.parser)

run_no = opt$run_no
outDir = opt$out_dir

### Classic memes

heading = function(text){

    cat(paste0("\n=============\n\n", text, "\n=========================\n\n"))

}

heading("Might go somewhere sunny. Sit on beach, look at ocean, collect sea shells. Might run tests on the sea shells.")

sessionInfo()
start_time = Sys.time()

### File list

all_files = Sys.glob(paste0(processDir, "*_out_annotated.txt"))

### Define file and phenotype

file_select = all_files[run_no]

pheno = gsub("_gcta_jma_out_annotated.txt", "", basename(file_select))

gcta_read = fread(file_select, data.table = FALSE)

### Define and read GWAS file

gwas_file = paste0(gwasDir, pheno, "_metal_output_1.TBL.gz")

gwas = fread(gwas_file, select = c("MarkerName", "Allele1", "Allele2", "Chromosome", "Position"))

### Reading complete

heading(paste0("File reading complete for: ", pheno, " (file read: ", file_select, ")"))


##### =========================== #####

### Start run

##### =========================== #####

### rsIDs from GCTA-COJO output

#rsids = c("rs756711137", "rs9273078") # Example SNP from DRUG
rsids = gcta_read$SNP

cat(paste0("\nSNPs selected for ", pheno, "\n\n"))
print(rsids)
cat(paste0("\n...\n\n"))


### Retrieve SNP positions for these rsIDs from the GWAS file

#snp_table = data.frame(rsid = rsids,
#                        chr = c(12, 6), # dbSNP chromosome for rs756711137
#                        pos = c(15363658, 32644532) # dbSNP position for rs756711137
#                        )

snp_table = data.frame(rsid = rsids,
                        chr = gwas$Chromosome[match(rsids, gwas$MarkerName)],
                        pos = as.numeric(gwas$Position[match(rsids, gwas$MarkerName)])
                        )

### Convert the table to Genomic ranges

snp_gr = GRanges(seqnames = paste0("chr", snp_table$chr),   # ensure "chr" prefix matches genes_gr
                    ranges = IRanges(start = snp_table$pos, end = snp_table$pos + 1),
                    rsid = snp_table$rsid
                    )

cat(paste0("\nSNPs positions retrieved\n\n"))
print(snp_gr)
cat(paste0("\n...\n\n"))

### Get gene table

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr = genes(txdb)

### Get nearest gene

nearest_hits = nearest(snp_gr, genes_gr)

##Build a results table

gene_results = data.frame(
                    rsid = snp_gr$rsid,
                    chr  = as.character(seqnames(snp_gr)),
                    hg38_pos = start(snp_gr),
                    gene_id = mcols(genes_gr)$gene_id[nearest_hits],
                    gene_chr = as.character(seqnames(genes_gr[nearest_hits])),
                    gene_start = start(genes_gr[nearest_hits]),
                    gene_end = end(genes_gr[nearest_hits]),
                    strand = as.character(strand(genes_gr[nearest_hits]))
                )

### Get entrz IDs and then convert the gene IDs in results table

entrez_ids = as.character(gene_results$gene_id)
gene_symbols = mapIds(org.Hs.eg.db, keys = entrez_ids, 
                       column = "SYMBOL", keytype = "ENTREZID", 
                       multiVals = "first")

gene_results$nearest_ncbi_gene = gene_symbols[as.character(gene_results$gene_id)]

### Merge with GCTA file

gcta_read = left_join(gcta_read, gene_results[,c("rsid", "hg38_pos", "nearest_ncbi_gene")], by = c("SNP" = "rsid"))
gcta_read = left_join(gcta_read, gwas[,c("MarkerName", "Allele1", "Allele2")], by = c("SNP" = "MarkerName"))

fwrite(gcta_read, file_select, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")




