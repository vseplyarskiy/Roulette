use strict;


my $chr=$ARGV[0];

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_input");

while (<IN>){
    chomp;
    my @a=split/\t/;
    my $out=$a[2];
    $out =~ s/\s/\_/g;
    open (OUT,">>/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$out");
    print  OUT "$chr $_\n";

}


