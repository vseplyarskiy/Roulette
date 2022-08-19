use strict;

my $chr=$ARGV[0];

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v3");
open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v4");

my %cor;
open (CTCF,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/CTCF_corection/hash/all_predict");
while (<CTCF>){
    chomp;
    my @a=split/\t/;
    my @b=split/\s/,$a[0];
    next unless $chr==$b[0];
    $a[2]=int($a[2]*100+0.5)/100;
    my $v=$a[2];
    if ($v==0){$v=0.33;$a[2]=0.33;}

    my $lam=-log($v);
    $lam=int($lam*100 +0.5)/100;

    if ($a[2]==0.99){$a[2]=0.9874; $lam=-log($a[2]); $lam=int($lam*1000 +0.5)/1000}
    if ($a[2]==1){$a[2]=0.996; $lam=-log($a[2]); $lam=int($lam*1000 +0.5)/1000}

    ${$cor{$b[1]}}{$a[1]}="$a[2] $lam";
}


open (CTCF,"/net/home/vseplyarsky/GNOMAD_Tufts/script/hypermutable/datah/all_predict_v2");
while (<CTCF>){
    chomp;
    my @a=split/\t/;
    my @b=split/\s/,$a[0];
    next unless $chr==$b[0];

    $a[2]=int($a[2]*100+0.5)/100;
    my $v=$a[2];
    if ($v==0){$v=0.33;$a[2]=0.33;}
    my $lam=-log($v);
    $lam=int($lam*100 +0.5)/100;

    if ($a[2]==0.99){$a[2]=0.9874; $lam=-log($a[2]); $lam=int($lam*1000 +0.5)/1000}
    if ($a[2]==1){$a[2]=0.996; $lam=-log($a[2]); $lam=int($lam*1000 +0.5)/1000}

    ${$cor{$b[1]}}{$a[1]}="$a[2] $lam";
}


my %low;
my %lowh;

open (low,"/net/home/vseplyarsky/GNOMAD_Tufts/qc/data/filters/$chr\_low");
while (<low>){
chomp;
my @a=split/\s/;
if ($a[1] eq 'hn_qc'){$lowh{$a[0]}=1}
else {$low{$a[0]}=1}
}


while (<IN>){
    chomp;
    my @a=split/\s/;
    my $h='high';
    my $pos=int($a[0]/100)*100; 
    if ($lowh{$pos}){ $h='low';}
    if ($low{$a[0]}){ $h='low';}
    my $add;
    if ($cor{$a[0]}){
        my $mu=$a[1];
        substr($mu,0,2)='';
        substr($mu,1,2)='';
	$add=${$cor{$a[0]}}{$mu};
    }
   my $lam=-log($a[2]); $lam=int($lam*1000 +0.5)/1000;
    print OUT "$a[0] $a[1] $a[2] $lam $h $add\n";
}

