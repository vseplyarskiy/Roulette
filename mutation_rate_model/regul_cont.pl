use strict;

my %ra;
my %cor;

open (IN,"../mut_model/validation/cont_distr");
open (OUT,">../mut_model/validation/QQ_norm");

while (<IN>){
    chomp;
    my @a=split/\s/;
    my @site=split/\//,$a[3];
    next unless $site[1]>=1000;
    my $v=int($a[2]*100 +0.5)/100;
    ${$ra{$a[0]}}{$v}=1;
    if ($a[4] >1.2 || $a[4] <1/1.2){${$cor{$a[0]}}{$a[1]}=$v}
}

foreach my $con (sort keys %ra){
    my @mm=sort {$a<=>$b} keys %{$ra{$con}};
    my $mm="$mm[0] $mm[-1]";
    print OUT "$con 100 100 $mm\n";
    foreach my $cor (sort {$a<=>$b} keys %{$cor{$con}}){
	print OUT "$con $cor ${$cor{$con}}{$cor} $mm\n";

    }}
