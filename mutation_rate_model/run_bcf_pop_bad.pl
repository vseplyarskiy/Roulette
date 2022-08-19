use strict;

my $chr=$ARGV[0];
`bcftools query -R /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/SFS/bad_fbcf_$chr  -f '%CHROM %POS %REF %ALT %AC_eas %AC_sas %AC_afr %AC_nfe %AC_amr\n' /net/data/gnomAD/gnomad.genomes.r3.0.sites_hg38.vcf.bgz -o /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/SFS/bad_result_$chr`
