use strict;

my @ar=</net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/*>;

open(OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/catalogue");

for my $i (0..$#ar){
    my @name=split/\//,$ar[$i];
    print OUT "$name[-1]\n";
}
