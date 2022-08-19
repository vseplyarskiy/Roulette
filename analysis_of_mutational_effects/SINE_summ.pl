use strict;

for my $chr (1..22){
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/hash/$chr\_ALUtestis_ind_denovo");
    while (<IN>){print "$_"}}
