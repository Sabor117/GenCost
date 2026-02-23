##### =========================== #####

### Setup environment

### Script is for munging sumstats prior to use with CoPheScan
### GWAS sumstats have 3 sources: GWAS catalog, INTERVENE and MVP+UKB
### Munged sumstats must have: SNPID (rsid ideally), beta, varbeta (square of SE), MAF, chr, pos
### Must also supply cc or quant and N - separate table

##### =========================== #####

### Packages

library(data.table)
library(dplyr)
library(stringr)

### Directories

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
mungedDir = paste0(mainDir, "processing/phewas_gwas_sumstats/munged_input/")

### INTERVENE GWAS

intervene_gwas_dir = Sys.glob("/scratch/project_2007428/users/FAHagenbeek/data/data_UKB_GxE_SESDisease/GWASsummstats/")

intervene_gwas_files = c(Epilepsy = paste0(intervene_gwas_dir, "Abou-Khalil_2018_epilepsy_BuildHG19"),
                            AllCancers = paste0(intervene_gwas_dir, "AllCancers.tsv"),
                            RheumArth = paste0(intervene_gwas_dir, "Ha_2020_RA_BuildHG19.tsv"),
                            Asthma = paste0(intervene_gwas_dir, "PreMunge/HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC_buildHG19.txt"),
                            BreastCancer = paste0(intervene_gwas_dir, "Michailidou_2017_BreastCancer_buildHG19.tsv"),
                            Gout = paste0(intervene_gwas_dir, "Tin_2019_urate_buildHG19_csv.txt"),
                            ProstateCancer = paste0(intervene_gwas_dir, "Schumacher_2018_ProstateCancer_buildHG19.tsv"),
                            T1D = paste0(intervene_gwas_dir, "PreMunge/Robertson_2021_T1D_buildHG19.txt"),
                            T2D = paste0(intervene_gwas_dir, "Mahajan.NatGenet2018b.T2D.European_buildHG19.txt"),
                            AtrialFib = paste0(intervene_gwas_dir, "Roselli_2018_AF_buildHG19.tsv"),
                            CAD = paste0(intervene_gwas_dir, "Nelson_2017_CAD_BuildHG19.txt"),
                            Melanoma = paste0(intervene_gwas_dir, "Rashkin_2020_Melanoma_buildHG19.tsv"),
                            MAjDepDis = paste0(intervene_gwas_dir, "Wray_2018_MDD_buildHG19.csv"),
                            LungCancer = paste0(intervene_gwas_dir, "Sakaue_2021_LungCancer_BuildHG19.tsv"),
                            HOsteoarth = paste0(intervene_gwas_dir, "Tachamazidou_2019_HipOA_buildHG19.tsv"),
                            KOsteoarth = paste0(intervene_gwas_dir, "Tachamazidou_2019_KneeOA_buildHG19.tsv"),
                            Appendicitis = paste0(intervene_gwas_dir, "Jiang_2021_Appendicitis_BuildHG19.tsv"),
                            CRC = paste0(intervene_gwas_dir, "Law_2019_colorectal_cancer.txt")
                            )

### MVP GWAS

mvp_gwas_files = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/phewas_gwas_sumstats/mvp_ukb_fin/*")

mvp_gwas_names = str_split_fixed(basename(mvp_gwas_files), "_meta_out", 2)[,1]

### GWAS catalog

gwas_catalog_files = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/phewas_gwas_sumstats/gwas_catalog/*")

gwas_catalog_names = ifelse(str_detect(gwas_catalog_files[-c(1:3)], "_buildGRCh\\d{2}_"),
                            str_extract(gwas_catalog_files[-c(1:3)], "(?<=_buildGRCh\\d{2}_)[^\\.]+"), # Extract after "_buildGRCh##_"
                            str_extract(gwas_catalog_files[-c(1:3)], "(?<=GCST\\d{8}_)[^\\.]+")        # Extract after "GCST########_"
                            )
gwas_catalog_names = c("ChronicPain", "BMI", "WHRadjBMI", gwas_catalog_names)

names(gwas_catalog_files) = gwas_catalog_names

### iPSYCH files

ipsych_gwas_files = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/phewas_gwas_sumstats/ipsych/*")

ipsych_names = c("ADHD2022_iPSYCH", "ASD2017_iPSYCH", "SCZ2022_iPSYCH", "MDD2025_iPSYCH")

names(ipsych_gwas_files) = ipsych_names


##### =========================== #####

### Beginning munge prune

##### =========================== #####

### Start with MVP as all files are in the same format

cat(paste0("\nSet-up complete. Starting munge.\n===============\n\n"))

cat(paste0("\nMunging MVP files\n======\n\n"))

for (nfile in 1:length(mvp_gwas_files)){

    currfile = mvp_gwas_files[nfile]

    cat(paste0("\nRunning on file ", nfile, " of ", length(mvp_gwas_files), " files (", basename(currfile), ").\n...\n\n"))

    if (file.exists(paste0(mungedDir, basename(currfile)))){

        cat(paste0("\nFile exists. Next.\n...\n\n"))

        next

    }

    curr_gwas = fread(currfile, data.table = FALSE)

    outfile = data.frame(snp = curr_gwas$rsid,
                            chr = curr_gwas[,1],
                            pos = curr_gwas$POS,
                            beta = curr_gwas$all_inv_var_meta_beta,
                            varbeta = (curr_gwas$all_inv_var_meta_sebeta ^ 2),
                            MAF = curr_gwas$ukbb_af_alt)

    outfile$snp[is.na(outfile$snp)] = paste0(curr_gwas[is.na(outfile$snp),1], "_", curr_gwas$POS[is.na(outfile$snp)], "_", curr_gwas$ALT[is.na(outfile$snp)], "_", curr_gwas$REF[is.na(outfile$snp)])

    fwrite(outfile, paste0(mungedDir, basename(currfile)), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

    cat(paste0("\nFile ", nfile, " munged and written.\n...\n"))
    print(head(outfile))
    print(summary(outfile))

    curr_gwas = NULL
    outfile = NULL

}

cat(paste0("\nMVP complete. Munging other files\n======\n\n"))

other_gwas_files = c(intervene_gwas_files, gwas_catalog_files, ipsych_gwas_files)

snp_id_cols = c("snp", "id", "rsid", "markername", "variant_id", "snptestid", "rs_id")
chr_cols = c("chromosome", "chr", "chrom", "#chr", "#chrom")
pos_cols = c("pos", "position", "base_pair_location", "bp", "pos_b37", "bp_hg19")
beta_cols = c("beta", "effect")
or_columns = c("odds_ratio", "or")
log_or_columns = c("logor")
se_cols = c("se", "standard_error", "stderr", "se_gc", "sebeta")
zscore_cols = c("z", "zscore")
maf_cols = c("eaf", "a1freq", "maf", "effect_allele_frequency", "frq_a_59851", "coded_af", "freq1", "af_alt", "effect_allele_freq")
sample_cols = c("n", "neff", "weight", "n_total_sum")
a1_cols = c("a1", "allele1", "allele_a", "effect_allele", "alt")
a0_cols = c("a0", "a2", "allele0", "allele2", "allele_b", "other_allele", "ref")

exclude_positions = "HanY_prePMID_asthma|T1D|fat-distn|GCST90000618_buildGRCh37_vitd|PLACEHOLDER"

silly_snp_files = "BreastCancer|fat-distn|PLACHOLDER"

for (nfile in 38:length(other_gwas_files)){

    currfile = other_gwas_files[nfile]

    cat(paste0("\nRunning on file ", nfile, " of ", length(other_gwas_files), " files (", basename(currfile), ").\n...\n\n"))

    outname = tools::file_path_sans_ext(basename(currfile))

    if (basename(currfile) == "Wray_2018_MDD_buildHG19.csv"){

        cat(paste0("\nCurrent file is Wray GWAS. Filling strange rows.\n\n"))

        curr_gwas = fread(currfile, data.table = FALSE, fill = TRUE)

        print(head(curr_gwas, n = 1))
        print(curr_gwas[9533405:9533415,])

    } else if (basename(currfile) == "pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz") {

        cat(paste0("\nCurrent file is MDD GWAS. Skipping top rows.\n\n"))

        curr_gwas = fread(currfile, data.table = FALSE, skip = "#CHROM")

    } else if (basename(currfile) == "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz") {

        cat(paste0("\nCurrent file is MDD GWAS. Skipping top rows.\n\n"))

        curr_gwas = fread(currfile, data.table = FALSE, skip = "CHROM")
        curr_gwas = na.omit(curr_gwas)

    } else {

        curr_gwas = fread(currfile, data.table = FALSE)

    }

    colnames(curr_gwas) = tolower(colnames(curr_gwas))
    gwas_cols = tolower(colnames(curr_gwas))

    ### Select CHR column
    ### CHR + POS necessary only for files without rsID
    ### Can be excluded IF rsID is present (checked manually based on exclude_positions)

    if (any(gwas_cols %in% chr_cols)){

        chr_col = which(gwas_cols %in% chr_cols)

    } else if (grepl(exclude_positions, currfile)) {

        cat(paste0("\nChromosome column excluded.\n......\n\n"))

    } else {

        cat(paste0("\nNO CHR COL FOR: ", currfile, "\n......\n\n"))

        break

    }

    ### Select POS column
    ### CHR + POS necessary only for files without rsID
    ### Can be excluded IF rsID is present (checked manually based on exclude_positions)

    if (any(gwas_cols %in% pos_cols)){

        pos_col = which(gwas_cols %in% pos_cols)

    } else if (grepl(exclude_positions, currfile)) {

        cat(paste0("\nPosition column excluded.\n......\n\n"))

    } else {

        cat(paste0("\nNO POS COL FOR: ", currfile, "\n......\n\n"))

        break

    }

    ### Select a1 column

    if (any(gwas_cols %in% a1_cols)){

        a1_col = which(gwas_cols %in% a1_cols)

    }

    ### Select a0 column

    if (any(gwas_cols %in% a0_cols)){

        a0_col = which(gwas_cols %in% a0_cols)

    }

    ### Select MAF column

    ### MDD and SCZ GWAS contains a frequency for cases + frequency for controls
    ### Very silly - create a total MAF based on this

    if (basename(currfile) == "pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz" | basename(currfile) == "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz"){

        curr_gwas$maf = (curr_gwas$fcas * curr_gwas$ncas + curr_gwas$fcon * curr_gwas$ncon) / (curr_gwas$ncas + curr_gwas$ncon)

        gwas_cols = tolower(colnames(curr_gwas))

    }

    if (basename(currfile) == "ADHD2022_iPSYCH_deCODE_PGC.meta.gz"){

        curr_gwas$maf = (curr_gwas$frq_a_38691 * curr_gwas$nca + curr_gwas$frq_u_186843 * curr_gwas$nco) / (curr_gwas$nca + curr_gwas$nco)

        gwas_cols = tolower(colnames(curr_gwas))

    }

    if (any(gwas_cols %in% maf_cols)){

        maf_col = which(gwas_cols %in% maf_cols)

    } else {

        cat(paste0("\nNO MAF COL FOR: ", currfile, "\n......\n\n"))

        break

    }

    ### Select N column
    ### Only necessary for Z-score. Proper sample sizes will be collected manually

    ### iPSYCH ADHD GWAS contains n-controls + n-cases and not total N

    if (basename(currfile) == "ADHD2022_iPSYCH_deCODE_PGC.meta.gz"){

        curr_gwas$n = curr_gwas$nca + curr_gwas$nco

        gwas_cols = tolower(colnames(curr_gwas))

    }

    ### iPSYCH SCZ GWAS contains NEFFDIV2 (n-effective divided by 2) for some reason

    if (basename(currfile) == "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz"){

        curr_gwas$neff = curr_gwas$neffdiv2 * 2

        gwas_cols = tolower(colnames(curr_gwas))

    }

    if (any(gwas_cols %in% sample_cols)){

        sample_col = which(gwas_cols %in% sample_cols)

    } else {

        cat(paste0("\nNO N COL FOR: ", currfile, "\n......\n\n"))

    }

    ### Effect and SE column
    ### Must be bundled together as effect is sometimes Z-score
    ### Effect can be beta/Z/OR/logOR - have to account for each one
    ### OR has SE as normal

    if (any(gwas_cols %in% beta_cols)){ ### First look for BETA

        effect_col = which(gwas_cols %in% beta_cols)

        ### Internal statement to get SE column

        if (any(gwas_cols %in% se_cols)){

            se_col = which(gwas_cols %in% se_cols)

        } else {

            cat(paste0("\nNO SE COL FOR: ", currfile, "\n......\n\n"))

            break

        }

    } else if (any(gwas_cols %in% or_columns)){ ### THEN CHECK FOR OR

        effect_col = which(gwas_cols %in% or_columns)

        ### Internal statement to get SE column

        if (any(gwas_cols %in% se_cols)){

            se_col = which(gwas_cols %in% se_cols)

        } else {

            cat(paste0("\nNO SE COL FOR: ", currfile, "\n......\n\n"))

            break

        }        
    
    } else if (any(gwas_cols %in% log_or_columns)){ ### THEN CHECK FOR LOG OR

        effect_col = which(gwas_cols %in% log_or_columns)

        cat(paste0("\nLOG OR detected for: ", currfile, " Converting to OR\n......\n\n"))

        curr_gwas[,effect_col] = exp(curr_gwas[,effect_col])

        ### Internal statement to get SE column
        ### Remember to exp(se) as you are correcting the logOR as well

        if (any(gwas_cols %in% se_cols)){

            se_col = which(gwas_cols %in% se_cols)

            curr_gwas[,se_col] = exp(curr_gwas[,se_col])

        } else {

            cat(paste0("\nNO SE COL FOR: ", currfile, "\n......\n\n"))

            break

        }        
    
    } else if (any(gwas_cols %in% zscore_cols)){ ### LASTLY CHECK FOR ZSCORE

        effect_col = which(gwas_cols %in% zscore_cols)

        cat(paste0("\nZ-score detected for: ", currfile, " Converting to beta + se\n......\n\n"))

        curr_gwas$beta1 = curr_gwas[,effect_col] / sqrt((2 * curr_gwas[,maf_col]) * (1 - curr_gwas[,maf_col]) * (curr_gwas[,sample_col] + (curr_gwas[,effect_col]^2)))

        curr_gwas$se = 1 / sqrt((2 * curr_gwas[,maf_col]) * (1 - curr_gwas[,maf_col]) * (curr_gwas[,sample_col] + (curr_gwas[,effect_col]^2)))

        effect_col = which(colnames(curr_gwas) %in% "beta1")
        se_col = which(colnames(curr_gwas) %in% "se")
    
    } else {

        cat(paste0("\nNO EFFECT COL FOR: ", currfile, "\n......\n\n"))

        break

    }

    ### Select SNP column
    ### Occasionally file contains both snpid + rsid
    ### Have to account for that and select rsid
    ### Some GWAS have no SNP column at all - must be created

    if (any(gwas_cols %in% snp_id_cols)){

        snp_col = which(gwas_cols %in% snp_id_cols)

        if (length(snp_col) > 1){

            rsid_checker = curr_gwas[,snp_col]

            columns_with_rs = which(apply(rsid_checker, 2, function(column) any(grepl("rs", column))))

            snp_col = snp_col[columns_with_rs]

            cat(paste0("\nMULTIPLE SNP COLUMNS FOR: ", currfile, " CHOOSING THIS COLUMN:\n\n"))
            print(head(curr_gwas[,snp_col]))

        }

    } else {

        cat(paste0("\nNO SNP COL FOR: ", currfile, " Making one.\n......\n\n"))

        curr_gwas$snpid = paste0(curr_gwas[,chr_col], "_", curr_gwas[,pos_col], "_", curr_gwas[,a1_col], "_", curr_gwas[,a0_col])

        snp_col = which(colnames(curr_gwas) %in% "snpid")

    }

    if (grepl(silly_snp_files, currfile)){

        curr_gwas$snpid_fixed = str_split_fixed(curr_gwas[,snp_col], ":", 4)[,1]

        id_fix2_vector = which(!(grepl("rs", curr_gwas$snpid_fixed)))

        curr_gwas$snpid_fixed[id_fix2_vector] = curr_gwas[id_fix2_vector, snp_col]

        snp_col = which(colnames(curr_gwas) %in% "snpid_fixed")

    }

    ### Create outfile based on selected columns
    ### Also based 

    if (grepl(exclude_positions, currfile)) {

        outfile = data.frame(snp = curr_gwas[,snp_col],
                            a1 = curr_gwas[,a1_col],
                            a0 = curr_gwas[,a0_col],
                            beta = curr_gwas[,effect_col],
                            varbeta = (curr_gwas[,se_col] ^ 2),
                            MAF = curr_gwas[,maf_col])

    } else {

        outfile = data.frame(snp = curr_gwas[,snp_col],
                            a1 = curr_gwas[,a1_col],
                            a0 = curr_gwas[,a0_col],
                            chr = curr_gwas[,chr_col],
                            pos = curr_gwas[,pos_col],
                            beta = curr_gwas[,effect_col],
                            varbeta = (curr_gwas[,se_col] ^ 2),
                            MAF = curr_gwas[,maf_col])

    }

    fwrite(outfile, paste0(mungedDir, outname, ".txt.gz"), quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

    cat(paste0("\nMunge completed for file ", nfile, " (", basename(currfile), ").\n\n"))
    print(head(outfile))
    print(summary(outfile))

    curr_gwas = NULL
    outfile = NULL

}


cat(paste0("\nScript complete. Munge over. Destroy all humans\n===============\n\n"))

### This secondary step is required for some of the files - selected manually
### No chr/pos columns which are required for CoPheScan input SNP pruning
### Have to be added manually

no_position_files = c(paste0(mungedDir, "HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC_buildHG19.txt.gz"),
                        paste0(mungedDir, "Robertson_2021_T1D_buildHG19.txt.gz"), 
                        Sys.glob(paste0(mungedDir, "fat-distn.giant.ukbb.meta*"))
                        paste0(mungedDir, "GCST90000618_buildGRCh37_vitd.txt.gz")
                        )
                        
### Uses UKB SNPs to add hg19 chr + pos columns
### Based on rsID
### Without rsID this is much more complicated

ukb_snps = fread("/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt",
                    data.table = FALSE, header = FALSE)

colnames(ukb_snps) = c("chr", "pos", "pos1", "snp", "snpid")

ukb_snps = ukb_snps[,c("snp", "chr", "pos")]

### For missing chr/pos files instead of using this LiftOver file - will use bcftools based on Tomoko's stuff

missing_snp_positions_file = fread(no_position_files[1], data.table = FALSE)

missing_snp_positions = data.frame(rsid = missing_snp_positions_file$snp)

for (nfile in 2:length(no_position_files)){

    missing_snp_positions_file = fread(no_position_files[nfile], data.table = FALSE)
    missing_snp_positions_file = data.frame(rsid = missing_snp_positions_file$snp)

    missing_snp_positions = rbind(missing_snp_positions, missing_snp_positions_file)

    missing_snp_positions = unique(missing_snp_positions)

}

fwrite(missing_snp_positions, paste0(mainDir, "processing/phewas_gwas_sumstats/missing_chr_pos_rsids.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", col.names = FALSE)


#######################
###
### SUPERCEDED CODE ###
###
#######################

#for (nfile in 1:length(no_position_files)){
#
#    curr_file = no_position_files[nfile]
#
#    file_read = fread(curr_file, data.table = FALSE)
#
#    file_read = left_join(file_read, ukb_snps, by = "snp")
#
#    fwrite(file_read, curr_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")
#
#}

#cat(paste0("\nNOW IT IS TRULY OVER.\n===============\n\n"))

#######################
###
### SUPERCEDED CODE ###
###
#######################

### Another tertiary munging issue
### One file with chr + pos columns but completely void of data

no_position_files = c(paste0(mungedDir, "GCST90000618_buildGRCh37_vitd.txt.gz"))

for (nfile in 1:length(no_position_files)){

    curr_file = no_position_files[nfile]

    cat(paste0("\nWorking on file: ", curr_file, "\n\n"))

    file_read = fread(curr_file, data.table = FALSE)

    if (all(is.na(file_read$pos))){

        file_read = file_read[,c("snp", "beta", "varbeta", "MAF")]

    }

    file_read = left_join(file_read, ukb_snps, by = "snp")

    file_read$chr = as.numeric(gsub("chr", "", file_read$chr))

    fwrite(file_read, curr_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

}

cat(paste0("\nNOW IT IS WELL AND TRULY OVER.\n===============\n\n"))



