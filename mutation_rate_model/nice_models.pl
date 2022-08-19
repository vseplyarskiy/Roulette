use strict;
open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/script/nice_models");

my %fi;
open (cat,"/net/home/vseplyarsky/GNOMAD_Tufts/script/cataloguels");
<cat>;<cat>;
while (<cat>){
    chomp;
    my $in =$_;
    substr($in,-9);
    my $fi="/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/Logit/$in";
    open (h,"$fi");
    my $h=<h>;chomp($h);

    my @a=split/\t/,$h;
    my $mod;
    if ($#a==1){$mod=$a[1]}
    else{
	if ($a[0]=~/\_cl/){$mod=$a[1];
			   for my $i (2..$#a){if ($a[$i] eq 'NA'){$a[$i]=1};$mod.="\t$a[$i]"}
	}}
    $fi{$fi}=$mod;
    print OUT "$fi,$mod\n";
    <cat>;
}


