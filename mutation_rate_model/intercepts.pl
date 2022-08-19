use strict;
use List::Util qw( min max );

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/intercepts");

for my $ex (1..4){
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/gen\ $ex\_penta");
    while (<IN>){
	chomp;
    my @a=split/\s/;
	next unless /X/;
    my $lt="$a[0] gen $ex met\t$a[2]";
	print OUT "$lt\n";    
    }

    for my $di (0..2){
	open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/pro\ $ex\ $di\_penta");
	while (<IN>){
	    chomp;
	    next unless /X/;

	    my @a=split/\s/;
	    my $lt="$a[0] pro $ex $di met\t$a[2]";
	    print OUT "$lt\n";
	}
    }
}
open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/inter_penta");
while (<IN>){
    chomp;
    next unless /X/;
    my @a=split/\s/;
    my $lt="$a[0] inter met\t$a[2]";
    print OUT "$lt\n";
}
 

