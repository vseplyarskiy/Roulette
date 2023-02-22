#!/bin/bash

vcf_header=$(cat << 'END' 
##fileformat=VCFv4.3
##contig=<ID=1,length=248956422,assembly=GRCh38>
##contig=<ID=2,length=242193529,assembly=GRCh38>
##contig=<ID=3,length=198295559,assembly=GRCh38>
##contig=<ID=4,length=190214555,assembly=GRCh38>
##contig=<ID=5,length=181538259,assembly=GRCh38>
##contig=<ID=6,length=170805979,assembly=GRCh38>
##contig=<ID=7,length=159345973,assembly=GRCh38>
##contig=<ID=8,length=145138636,assembly=GRCh38>
##contig=<ID=9,length=138394717,assembly=GRCh38>
##contig=<ID=10,length=133797422,assembly=GRCh38>
##contig=<ID=11,length=135086622,assembly=GRCh38>
##contig=<ID=12,length=133275309,assembly=GRCh38>
##contig=<ID=13,length=114364328,assembly=GRCh38>
##contig=<ID=14,length=107043718,assembly=GRCh38>
##contig=<ID=15,length=101991189,assembly=GRCh38>
##contig=<ID=16,length=90338345,assembly=GRCh38>
##contig=<ID=17,length=83257441,assembly=GRCh38>
##contig=<ID=18,length=80373285,assembly=GRCh38>
##contig=<ID=19,length=58617616,assembly=GRCh38>
##contig=<ID=20,length=64444167,assembly=GRCh38>
##contig=<ID=21,length=46709983,assembly=GRCh38>
##contig=<ID=22,length=50818468,assembly=GRCh38>
##FILTER=<ID=low,Description="Low quality regions as determined by gnomAD sequencing metrics. Mappability<0.5;overlap with>50nt simple repeat;ReadPosRankSum>1;0 SNVs in 100bp window.">
##FILTER=<ID=SFS_bump,Description="Pentamer context with abnormal SFS. The fraction of high-frequency SNVS [0.0005<MAF<=0.2] is greater than 1.5x mutation rate controlled average. Tends to be repetitive contexts.">
##FILTER=<ID=TFBS,Description="Transcription factor binding site as determined by overlap with ChIP-seq peaks.">
##FILTER=<ID=high,Description="Nothing suspicious yet identified.">
##INFO=<ID=PN,Number=1,Type=String,Description="Pentanucleotide context">
##INFO=<ID=MR,Number=1,Type=Float,Description="Roulette mutation rate estimate">
##INFO=<ID=AR,Number=1,Type=Float,Description="Adjusted Roulette mutation rate estimate">
##INFO=<ID=MG,Number=1,Type=Float,Description="gnomAD mutation rate estimate (Karczewski et al. 2020)">
##INFO=<ID=MC,Number=1,Type=Float,Description="Carlson mutation rate estimate (Carlson et al. 2018)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
END
)

f_in=$1
roulette_dir=$2
tmp_dir=".tmp/"
mkdir -p "$tmp_dir"
## Create region file
f_region="$tmp_dir/"regions.tsv
tail -n +2 "$f_in" | cut -f1,2 > "$f_region"
## Create the vcf file
tmp_vcf="$tmp_dir/"`basename $f_in .tsv`".vcf"
echo -e "$vcf_header" > "$tmp_vcf"
tail -n +2 "$f_in" | awk -F'\t' -v OFS='\t' '{print $1, $2, ".", $3, $4, ".", ".", "."}' >> "$tmp_vcf"
bgzip "$tmp_vcf"
bcftools index "$tmp_vcf".gz
## Annotate using all chromosomes
out_vcf="$tmp_dir/"`basename $f_in .tsv`".1.vcf"
bcftools annotate --no-version -a "$roulette_dir/"1_rate_v5.2_TFBS_correction_all.vcf.bgz \
    -c FILTER,INFO -R "$f_region" "$tmp_vcf".gz > "$out_vcf"
bgzip "$out_vcf"
bcftools index "$out_vcf".gz
rm "$tmp_vcf".gz "$tmp_vcf".gz.csi 
for ii in {2..22}
do
    out_vcf="$tmp_dir/"`basename $f_in .tsv`".$ii.vcf"
    out_vcf_prev="$tmp_dir/"`basename $f_in .tsv`".$((ii-1)).vcf.gz"
    bcftools annotate --no-version -a "$roulette_dir/$ii"_rate_v5.2_TFBS_correction_all.vcf.bgz \
	-c FILTER,INFO -R "$f_region" "$out_vcf_prev" > "$out_vcf"
    bgzip "$out_vcf"
    bcftools index "$out_vcf".gz
    rm "$out_vcf_prev" "$out_vcf_prev".csi
done
## Extract info back to tsv format
f_out=`basename $f_in .tsv`.roulette.tsv
echo -e 'CHROM\tPOS\tREF\tALT\tFILTER\tPN\tMR\tAR\tMG\tMC' > "$f_out"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/PN\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
    "$out_vcf".gz >> "$f_out"
rm "$out_vcf".gz "$out_vcf".gz.csi
rm "$f_region"
