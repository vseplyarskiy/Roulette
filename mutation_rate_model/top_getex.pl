use strict;



open (IN,"/net/home/vseplyarsky/TOPMEDs/data/introns_plus_v4");
my %int;
while (<IN>){
    chomp;
    my @a=split/\s/;
    $int{$a[4]}=1;
}
open (IN,"/net/home/vseplyarsky/TOPMEDs/tracks_to_correlate/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")||die;
my %exp;
<IN>;<IN>;<IN>;
while (<IN>){
    chomp; my @a=split/\t/;
    next unless $int{$a[1]};
    $exp{$a[1]}=$a[51];
  #  for my $i (0..$#a){	print "$i $a[$i]\n";}	last;
}


my @values = sort {$a <=> $b} values %exp; 
print "$#values\n";
# Print 95% percentile
print $values[sprintf("%.0f",(0.2*($#values)))];

