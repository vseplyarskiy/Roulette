use strict;

my $chr=$ARGV[0];

my %qc; my %c; my %dp;
open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/data/all_variants_chr$chr.csv.gz|") ||die;
my $win=100;
while (<IN>){
    chomp;
    my @a=split/\,/;
    print "$a[8] $a[9] $a[16]   $_\n";
  #  last;
    next unless $a[8] eq 'snv';
    my $pos=int($a[0]/$win)*$win;
    $qc{$pos}+=$a[16];
    $dp{$pos}+=$a[9];
    print "$a[9]\n";
    $c{$pos}++;
#    my $c=($a[0]/$win);
}


#open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_qc_dp_v2") ||die;


foreach my $c (sort {$a<=>$b} keys %c){
    my $pos=$c;
    print OUT "$c $qc{$pos} $dp{$pos} $c{$pos}\n";
}
