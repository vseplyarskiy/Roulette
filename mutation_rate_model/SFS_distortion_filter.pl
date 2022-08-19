use strict;


my $chr=$ARGV[0];
open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/script/rate_evan_penta_outliers_v2")||die;
my %filter;
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $fi="$a[2] $a[0]";
    $filter{$fi}=1;
}


open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/script/rate_evan_penta_outliers_v4")||die;
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $fi="$a[2] $a[0]";
    print "$fi\n";
    $filter{$fi}=1;
}
#my $fi='TTTAA_A 0.02';
#$filter{$fi}=1;

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v4")||die;
open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5")||die;

while (<IN>){
    chomp;
#    next unless /TTTAA_A/;
    if (/low/){print OUT "$_\n";next};
    my @a=split/\s/;
    my $val=$a[3];
    if ($a[5]){$val=$a[6]};
    if ($val<0.21){$val=int($val*50)/50}
    if ($val>0.21){$val=int($val*10)/10}
    if ($val>1.2){$val=1.4}
    my $fi="$a[1] $val";
#    print "$fi $filter{$fi} rate_v4\n";
    if (!($filter{$fi})){print OUT "$_\n";next}
    else {print OUT "$a[0] $a[1] $a[2] $a[3] SFS_bump $a[5] $a[6]\n";
#	  print "$a[0] $a[1] $a[2] $a[3] SFS_bump $a[5] $a[6]\n";
next};
}




