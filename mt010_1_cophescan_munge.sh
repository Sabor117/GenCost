#!/bin/bash

#SBATCH --job-name=mt10_1_cophescan_munge
#SBATCH --time=24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt10_1_cophescan_munge.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt10_1_cophescan_munge.log

##### =========================== #####

### Setup environment

##### =========================== #####

module load r-env
module load samtools

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron


##### =========================== #####

### Start script

##### =========================== #####

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
workingDir="${mainDir}processing/phewas_gwas_sumstats/"
mungedDir="${mainDir}processing/phewas_gwas_sumstats/munged_input/"

printf "Script starts:\n\n"
printf "NOBODY EXPECTS THE SEB INQUISTION:\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/meta_step10_1_cophescan_munge.R"
printf "\n\n"
printf "==========================================\n\n"

#Rscript ${mainDir}scripts/meta_step10_1_cophescan_munge.R

printf "\n\n"
printf "MUNGEING MADE.\n\n"
printf "\n\n"
printf "Adding missing chr/pos.\n\n"

missing_snp_positions="${workingDir}missing_chr_pos_rsids.txt"
vcf_file="${mainDir}processing/misc_data/GCF_000001405.40.gz"

tabix -p vcf ${vcf_file}

### Extract chromosome and position from VCF only once (if not already done)

if [ ! -f "$workingDir/missing_rsid_to_chrpos_b38.tsv" ]; then

    echo "Generating missing_rsid_to_chrpos_b38.tsv"
    echo 'bcftools query -i "ID=@$missing_snp_positions" -f '%ID\t%CHROM\t%POS\n' "$vcf_file" > "$workingDir/missing_rsid_to_chrpos_b38.tsv"'
    printf "\n....\n\n"
    bcftools query -i "ID=@$missing_snp_positions" -f '%ID\t%CHROM\t%POS\n' "$vcf_file" > "$workingDir/missing_rsid_to_chrpos_b38.tsv"
    printf "\n....\n\n"

fi

### Join with sequence report if needed (adjust if this varies per file)

if [ ! -f "$workingDir/missing_rsid_to_chrpos_b38_with_chr.tsv" ]; then

    echo "Joining with sequence report."
    printf "\n....\n\n"
    awk 'NR==FNR {seq[$9] = $12; next} $2 in seq {print $1, seq[$2], $3}' \
    "$workingDir/sequence_report.tsv" \
    "$workingDir/missing_rsid_to_chrpos_b38.tsv" > "$workingDir/missing_rsid_to_chrpos_b38_with_chr.tsv"
    printf "\n....\n\n"

fi

### List of files which have no chromosome or position data

file_array=(
#			"fat-distn.giant.ukbb.meta-analysis.bmi.combined.tbl.txt.gz" 
#			"HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC_buildHG19.txt.gz" 
#			"Robertson_2021_T1D_buildHG19.txt.gz" 
#			"fat-distn.giant.ukbb.meta-analysis.whradjbmi.combined.tbl.txt.gz"
            "GCST90000618_buildGRCh37_vitd.txt.gz"
			)


### Process each file

for currfile in "${file_array[@]}"; do

    echo "Processing $currfile."
    printf "\n....\n\n"

    input_path="$mungedDir/$currfile"
    output_path="$input_path.with_chrposb38.tsv"

	### Add chr and pos columns
	zcat "$input_path" | awk '
	NR==FNR {chr[$1] = $2; pos[$1] = $3; next}
	FNR==1 {OFS="\t"; print "chr", "pos", $0; next}
	$1 in pos {OFS="\t"; print chr[$1], pos[$1], $0}
	' "$workingDir/missing_rsid_to_chrpos_b38_with_chr.tsv" - > "$output_path"

    echo "Current file: $currfile processed."
    printf "\n....\n\n"

done

date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

