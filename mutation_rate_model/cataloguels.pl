use strict;

my @ar=</net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/Logit/*>;

open(OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/script/cataloguels");
print OUT "$#ar\n";
for my $i (0..$#ar){
    my @name=split/\//,$ar[$i];
    print OUT "$name[-1]\n";
}
