use strict;

my $run=$ARGV[0];

my $st=($run-1)*1013;
my $end=($run)*1013+1;
my %cat;
open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/catalogue");
my $i; 
while (<IN>){
    $i++;
    chomp;
    if ($i>$st && $i<$end){$cat{$_}=$i}
}



foreach my $name (sort keys %cat){
#    print "$name\n";
    `perl /net/home/vseplyarsky/GNOMAD_Tufts/script/filter_files.pl $name`;
    `Rscript /net/home/vseplyarsky/GNOMAD_Tufts/script/rscript.R $name\_cl`;
    `rm /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$name\_cl`;

}
