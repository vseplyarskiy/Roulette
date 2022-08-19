use strict;

open (IN,"../mut_model/validation/all_ra");
open (OUT,">../mut_model/validation/all_ra_cl");

while (<IN>){
    chomp;
    my @a=split/\s/;
    next unless $#a==5;
    print OUT "$_\n";
}
