use strict;

#my $fi=$ARGV[0]; 

my @fi=</net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/TCACT*>;
my @fi1=</net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/TCTGG*>;

    push (@fi,@fi1);
foreach my $fi (@fi){
open (IN,"$fi");#"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$fi")||die;
open  (OUT,">$fi\_cl");#(OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$fi\_cl");
print "$fi\n";
my $out;
while (<IN>){
    chomp;
    my @a=split/\t/;
    next unless $#a==12;
     $out.="$_\n";

}
print OUT "$out";
}
