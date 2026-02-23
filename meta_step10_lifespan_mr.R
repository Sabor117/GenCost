##### =========================== #####

### Setup environment

### Script to run MR
### Format is designed for Meta V.04 + LifeGen sumstats

##### =========================== #####

library(data.table)
library(MendelianRandomization, lib.loc = "/projappl/project_2007428/RPackages_421/")

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
metaDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_01/jma_combined/")
tempDir = paste0(mainDir, "tmpdir/")

phenotypes = c("IN_ALL", "DRUG_ALL", "INOUT_ALL", "PRIM_ALL")

lifegen_gwas_file = paste0(mainDir, "processing/misc_data/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz")
lifegen_gwas = fread(lifegen_gwas_file, data.table = FALSE)


##### =========================== #####

### Run script

##### =========================== #####

for (i in 1:length(phenotypes)){

    phenotype = phenotypes[i]

    meta_file = paste0(metaDir, phenotype, "_gcta_jma_out_annotated.txt")
    meta_snps = fread(meta_file, data.table = FALSE)

    lifegen_snps = lifegen_gwas[lifegen_gwas$rsid %in% meta_snps$SNP,]
    meta_snps = meta_snps[meta_snps$SNP %in% lifegen_snps$rsid,]

    lifegen_snps = lifegen_snps[order(lifegen_snps$rsid),]
    meta_snps = meta_snps[order(meta_snps$SNP),]

    mr_input_obj = mr_input(bx = meta_snps$b,
                                bxse = meta_snps$se,
                                by = lifegen_snps$beta1,
                                byse = lifegen_snps$se,
                                exposure = phenotype,
                                outcome = "Lifespan",
                                snps = meta_snps$SNP)

    mr_ivw_stats = mr_ivw(mr_input_obj)
    #mr_lasso_stats = mr_lasso(mr_input_obj)

    cat(paste0("\n", phenotype, " lifespan association is ", round(mr_ivw_stats$Estimate, 3), " [", round(mr_ivw_stats$CILower, 2), "-", round(mr_ivw_stats$CIUpper, 2), "], p = ", signif(mr_ivw_stats$Pvalue, 3), "\n"))

}


mr_plot(
        mr_input_obj,
        error = TRUE,
        line = "ivw",
        orientate = FALSE,
        interactive = TRUE
        )



