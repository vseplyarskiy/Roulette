use strict;

my %cat;

open (IN,"cataloguels")||die;
while (<IN>){
    chomp;
    my @a=split/\./;
    substr($a[0],-3)='';
    $cat{$a[0]}=1;
#    print "$a[0]\n";
}


open (IN,"../mut_model/for_LASSO/catalogue")||die;
while (<IN>){
    chomp;
   if (!($cat{$_})){if (/CG/){next}print "$_\n"};
#print "XXXXX $_\n"
}


