##### =========================== #####
###
### build_canonical_rsid_lookup.R
###
### PURPOSE: builds ONE canonical, conflict-checked chr:pos:allele -> rsID
### lookup table combining the UKB hg38 reference with any native rsIDs
### that survived qc_step4_1 for any cohort, so meta_step1 can assign the
### SAME rsID to the same physical variant regardless of which cohort it
### came from. See chat discussion for the original vulnerability this
### fixes (split rsID identity across cohorts -> METAL treats one locus
### as multiple SNPs).
###
### MEMORY NOTES (read this before running - previous version crashed
### at 200GB):
### The original version materialised the full UKB reference TWICE (once
### per allele order) and deduplicated with base R's unique.data.frame(),
### which is not memory-efficient at genome scale. This version:
###   1. Uses an ORDER-INDEPENDENT allele key (alleles sorted before
###      pasting) so each variant needs only ONE row, not two - halves
###      the UKB reference immediately.
###   2. Stays in data.table the whole way through and uses data.table's
###      radix-sort-based unique()/duplicated(), which is dramatically
###      lower-memory than the base R data.frame equivalents.
###   3. Reads only the 5 columns it needs from each cohort's
###      alleles_aligned file, not the full sumstats width.
###   4. Processes chromosome-by-chromosome, writing a small per-chr
###      lookup to disk and rm()/gc()-ing between chromosomes, so peak
###      memory is bounded by ONE chromosome's worth of data rather than
###      the whole genome at once. Chromosome 1 (~8% of the genome) is
###      the effective peak, not all 22+X.
###
### If this still runs out of memory on chr1/chr2 specifically, the
### remaining lever is to also chunk WITHIN a chromosome (fread's `skip`/
### line-count args or process by position bins), but that's unlikely to
### be needed after the above - a single chromosome's worth of UKB
### reference data should be in the low single-digit GB range, not
### hundreds.
###
##### =========================== #####

library(data.table)
setDTthreads(4)   # tune to your allocated cores - fread/dedup parallelise, they don't need to be single-threaded
options(scipen = 999)

##### =========================== #####
### Progress helpers
###
### These print a FULL LINE per call rather than using carriage-return
### (\r) redraws - \r doesn't render sensibly once stdout is redirected
### to a SLURM/nohup log file (you just get a wall of \r characters), so
### a new timestamped line per update is far more useful for tailing a
### log than a "real" terminal progress bar.
##### =========================== #####

script_start_time = Sys.time()

format_elapsed = function(start_time) {
    secs = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    sprintf("%dm%02ds", floor(secs / 60), round(secs %% 60))
}

print_progress = function(current, total, label = "", bar_width = 30, start_time = NULL) {

    pct = current / total
    filled = round(bar_width * pct)
    bar = paste0("[", strrep("=", filled), strrep(" ", bar_width - filled), "]")

    eta_str = ""
    if (!is.null(start_time) && current > 0) {
        secs_elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        secs_per_unit = secs_elapsed / current
        secs_remaining = secs_per_unit * (total - current)
        eta_str = sprintf(" | elapsed %s | ETA %s", format_elapsed(start_time),
                            sprintf("%dm%02ds", floor(secs_remaining / 60), round(secs_remaining %% 60)))
    }

    cat(sprintf("[%s] %s %s %d/%d (%.0f%%)%s\n", format(Sys.time(), "%H:%M:%S"),
                 label, bar, current, total, pct * 100, eta_str))

}

cat(paste0("\n===========================================================\n",
            "build_canonical_rsid_lookup.R starting at ", format(script_start_time), "\n",
            "===========================================================\n\n"))

mainDir    = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
sumstatDir = paste0(mainDir, "outputs/cohort_sumstats/")
tempDir    = paste0(mainDir, "tmpdir/")
logDir     = paste0(mainDir, "outputs/qc_logs/")
chunkDir   = paste0(mainDir, "processing/misc_data/rsid_lookup_by_chr/")

system(paste0("mkdir -p ", logDir, " ", chunkDir))

chrom_list = c(1:22, "X")

##### =========================== #####
### 1. UKB hg38 liftover-derived reference - read once, split by
###    chromosome in memory (cheap - it's just an index split), then
###    process one chromosome at a time.
##### =========================== #####

ukb_hg38_file = paste0(mainDir, "processing/misc_data/ukb_ALLSNPs_hg38_liftOver_output.out")

ukb_hg38 = fread(ukb_hg38_file, tmpdir = tempDir,
                  col.names = c("chr", "pos", "pos1", "rsid", "hg19_snpid"))

ukb_hg38[, chr := gsub("chr", "", chr)]
ukb_hg38[, chr := gsub("^23$", "X", chr)]

## Pre-locate every cohort's alleles_aligned file with a native rsid
## column ONCE (cheap header-only peek), so we don't re-glob per chromosome

aligned_files = Sys.glob(paste0(sumstatDir, "*/*_alleles_aligned.txt.gz"))

cat(paste0("\nScanning ", length(aligned_files), " alleles_aligned files for a surviving native rsid column...\n\n"))

rsid_carrying_files = c()
cohorts_seen = c()

prescan_start = Sys.time()

for (i in seq_along(aligned_files)) {

    f = aligned_files[i]
    cohort_name = strsplit(basename(f), "_")[[1]][1]

    if (i %% 10 == 0 || i == length(aligned_files)) print_progress(i, length(aligned_files), label = "  header scan", start_time = prescan_start)

    if (cohort_name %in% cohorts_seen) next

    hdr = fread(f, tmpdir = tempDir, nrows = 1)
    if ("rsid" %in% colnames(hdr)) {
        rsid_carrying_files = c(rsid_carrying_files, f)
        cohorts_seen = c(cohorts_seen, cohort_name)
        cat(paste0("    -> ", cohort_name, " carries a native rsid column\n"))
    }

}

cat(paste0("\n", length(rsid_carrying_files), " cohorts have a surviving native rsid column: ",
            paste(cohorts_seen, collapse = ", "), "\n\n"))

##### =========================== #####
### 2. Process one chromosome at a time
##### =========================== #####

chr_loop_start = Sys.time()

for (chr_i in seq_along(chrom_list)) {

    this_chr = chrom_list[chr_i]
    this_chr_start = Sys.time()

    print_progress(chr_i - 1, length(chrom_list), label = "chromosomes", start_time = chr_loop_start)
    cat(paste0("\n--- [", format(Sys.time(), "%H:%M:%S"), "] Processing chr", this_chr,
                " (", chr_i, "/", length(chrom_list), ") ---\n\n"))

    ukb_chr = ukb_hg38[chr == this_chr]

    if (nrow(ukb_chr) == 0) { cat(paste0("No UKB rows for chr", this_chr, ", skipping.\n")); next }

    cat(paste0("  UKB reference rows for chr", this_chr, ": ", format(nrow(ukb_chr), big.mark = ","), "\n"))

    alleles = tstrsplit(ukb_chr$hg19_snpid, ":", fixed = TRUE)

    ## Order-independent allele key: sort the two alleles alphabetically
    ## before pasting, so a1/a2 vs a2/a1 collapse to the SAME key. This is
    ## the single biggest memory saving - no more "both orders" doubling.

    a_sorted = matrix(c(alleles[[3]], alleles[[4]]), ncol = 2)
    a_sorted = t(apply(a_sorted, 1, sort))

    ukb_chr_long = data.table(
        full_snpid = paste(this_chr, ukb_chr$pos, a_sorted[,1], a_sorted[,2], sep = ":"),
        rsid = ukb_chr$rsid,
        source = "ukb_hg38_ref"
    )

    rm(ukb_chr, alleles, a_sorted); gc(verbose = FALSE)

    ## Native rsIDs for this chromosome, read with `select` so we never
    ## pull in the columns we don't need (gno_ref, gno_alt, af_alt, freq,
    ## FLIP, snpid)

    native_list = list()

    if (length(rsid_carrying_files) > 0) cat(paste0("  Scanning ", length(rsid_carrying_files), " native-rsid cohort file(s) for chr", this_chr, "...\n"))

    for (j in seq_along(rsid_carrying_files)) {

        f = rsid_carrying_files[j]
        cohort_label = strsplit(basename(f), "_")[[1]][1]

        dt = fread(f, tmpdir = tempDir, select = c("chr", "pos", "a1", "a0", "rsid"))
        dt = dt[chr == this_chr]

        cat(paste0("    [", j, "/", length(rsid_carrying_files), "] ", cohort_label, ": ",
                    format(nrow(dt), big.mark = ","), " rows on chr", this_chr, "\n"))

        if (nrow(dt) == 0) next

        a_sorted = t(apply(matrix(c(dt$a1, dt$a0), ncol = 2), 1, sort))

        native_list[[length(native_list) + 1]] = data.table(
            full_snpid = paste(this_chr, dt$pos, a_sorted[,1], a_sorted[,2], sep = ":"),
            rsid = dt$rsid,
            source = paste0(cohort_label, "_native_survived_liftover")
        )

        rm(dt, a_sorted)

    }

    native_chr = if (length(native_list) > 0) rbindlist(native_list) else data.table(full_snpid = character(0), rsid = character(0), source = character(0))
    rm(native_list); gc(verbose = FALSE)

    native_chr = native_chr[grepl("^rs[0-9]+$", rsid)]

    combined_chr = unique(rbindlist(list(ukb_chr_long, native_chr)))
    rm(ukb_chr_long, native_chr); gc(verbose = FALSE)

    ## Conflict detection & resolution - same logic as before, just scoped
    ## to one chromosome at a time now

    conflict_ids = unique(combined_chr[duplicated(combined_chr, by = "full_snpid")]$full_snpid)

    if (length(conflict_ids) > 0) {

        conflicts = combined_chr[full_snpid %in% conflict_ids]
        fwrite(conflicts, paste0(logDir, "rsid_lookup_conflicts.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE, append = (this_chr != chrom_list[1]))

        combined_chr[, rsid_num := suppressWarnings(as.numeric(gsub("rs", "", rsid)))]
        combined_chr = combined_chr[!is.na(rsid_num)]
        setorder(combined_chr, full_snpid, rsid_num)
        combined_chr = combined_chr[, .SD[1], by = full_snpid]
        combined_chr[, rsid_num := NULL]

        cat(paste0(length(conflict_ids), " conflicts on chr", this_chr,
                    " - logged to rsid_lookup_conflicts.txt, resolved by lowest rsID number.\n"))

    }

    chr_lookup = unique(combined_chr[, .(full_snpid, rsid)])
    stopifnot(!any(duplicated(chr_lookup$full_snpid)))

    fwrite(chr_lookup, paste0(chunkDir, "canonical_rsid_lookup_chr", this_chr, ".txt.gz"),
           sep = "\t", quote = FALSE, row.names = FALSE, compress = "gzip")

    cat(paste0("  chr", this_chr, " done: ", format(nrow(chr_lookup), big.mark = ","),
                " unique rows written, ", length(conflict_ids), " conflicts resolved, took ",
                format_elapsed(this_chr_start), " (total elapsed ", format_elapsed(script_start_time), ")\n"))

    rm(combined_chr, chr_lookup); gc(verbose = FALSE)

}

print_progress(length(chrom_list), length(chrom_list), label = "chromosomes", start_time = chr_loop_start)
cat(paste0("\nAll chromosomes processed in ", format_elapsed(chr_loop_start), ".\n\n"))

##### =========================== #####
### 3. Concatenate the per-chromosome lookups into the final file
###    (each chunk is already deduplicated, so this is a cheap append)
##### =========================== #####

chr_files = paste0(chunkDir, "canonical_rsid_lookup_chr", chrom_list, ".txt.gz")
chr_files = chr_files[file.exists(chr_files)]

cat(paste0("\nConcatenating ", length(chr_files), " per-chromosome lookup files...\n\n"))

concat_start = Sys.time()
chr_pieces = list()

for (k in seq_along(chr_files)) {

    chr_pieces[[k]] = fread(chr_files[k], tmpdir = tempDir)
    print_progress(k, length(chr_files), label = "  concatenating", start_time = concat_start)

}

final_lookup = rbindlist(chr_pieces)
rm(chr_pieces); gc(verbose = FALSE)

stopifnot(!any(duplicated(final_lookup$full_snpid)))

fwrite(final_lookup, paste0(mainDir, "processing/misc_data/canonical_rsid_lookup.txt.gz"),
       sep = "\t", quote = FALSE, row.names = FALSE, compress = "gzip")

cat(paste0("\n===========================================================\n",
            "Done. Final canonical rsID lookup: ", format(nrow(final_lookup), big.mark = ","),
            " unique chr:pos:allele(sorted) -> rsID mappings.\n",
            "Saved to: ", mainDir, "processing/misc_data/canonical_rsid_lookup.txt.gz\n",
            "Total run time: ", format_elapsed(script_start_time), " (finished ", format(Sys.time()), ")\n",
            "===========================================================\n\n"))