use strict;
for my $di (0..2){
    foreach my $fi ("GCCGA_T_pro_4_$di\_met1", "ACCGA_T_pro_4_$di\_met1"){

open (IN,"../mut_model/for_LASSO/context/$fi")||die "../mut_model/for_LASSO/context/$fi";
print "$fi\n";
open (OUT,">test_$fi");
while (<IN>){
    chomp;
    my @a=split/\t/;
    next unless $#a==12;
    print OUT "$_\n";

}
		}
}
