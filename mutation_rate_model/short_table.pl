use strict;

my $chr=$ARGV[0];

my %qc; my %c; my %dp;
open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/data/all_variants_chr$chr.csv.gz|") ||die;
open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol") ||die;

my $win=100;
while (<IN>){
    chomp;
    my @a=split/\,/;
#    print "$a[8]\n";
#    last;
#    next unless $a[8] eq 'snv' || $a[8] eq 'multi-snv';
    next unless $a[4]<=30;
    my $m="$a[1]\_$a[2]";
    my $l=length($m); 
    next unless $l ==3;
   print OUT "$a[0] $m\n";
}

`gzip  /net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol > /net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol.gz`;
