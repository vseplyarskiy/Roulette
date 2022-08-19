use strict;

my %pair;
open (IN,"rate_evan_penta_outliers_v3");
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $p="$a[2] $a[0]";
    $pair{$p}=1;
}


#my $chr=$ARGV[0];
for my $chr (1..22){
    print "$chr\n";
open (IN,"../mut_model/SFS/input_$chr");
open (OUT,">../mut_model/SFS/bad_fbcf_$chr");
    open (OUT1,">../mut_model/SFS/good_fbcf_$chr");
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $p=int($a[2]*50)/50;
    $p="$a[1] $p";
    if ($pair{$p}){
	print OUT "chr$chr\t$a[0]\n";}

    else { print OUT1 "chr$chr\t$a[0]\n";}
}}

