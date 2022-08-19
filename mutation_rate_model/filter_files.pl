use strict;

my $fi=$ARGV[0]; 

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$fi")||die;
open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$fi\_cl");
my $out;
while (<IN>){
    chomp;
    my @a=split/\t/;
    next unless $#a==12;
     $out.="$_\n";

}
print OUT "$out";
