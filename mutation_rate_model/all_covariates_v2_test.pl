use strict;
use List::Util qw( min max );
my $chr=$ARGV[0];
my $sc=100;
my $sc1=10000;


open (IN,"/net/home/vseplyarsky/TOPMEDs/data/introns_plus_v4");
my %int;
#my %int1;

my %cord_g;
while (<IN>){
    chomp;
    my @a=split/\s/;
    $int{$a[4]}=1;
}
open (IN,"/net/home/vseplyarsky/TOPMEDs/tracks_to_correlate/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")||die;
my %exp;
<IN>;<IN>;<IN>;
while (<IN>){
    chomp; my @a=split/\t/;
    next unless $int{$a[1]};
    $exp{$a[1]}=$a[51];
  #  for my $i (0..$#a){        print "$i $a[$i]\n";}   last;
}

my @values = sort {$a <=> $b} values %exp;

my %bin;
for my $i (1..4){
# Print 95% percentile
    $bin{$i}=$values[sprintf("%.0f",(0.2*$i*($#values)))];
#    print "$bin{$i}\n";
}

my %dir;
open (IN4,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/OK-seq/rpe_short_hg38.txt")||die;
while (<IN4>){
    chomp; my @a=split/\t/;
    my $c=substr($a[0],3);
    next unless $c eq $chr;
    my $pos=(int($a[1]/$sc1))*$sc1;
    $dir{$pos}=$a[3];}

my %met;
open (met, "/net/home/vseplyarsky/TOPMEDs/tracks_to_correlate/Testis_methylation/data/$chr\_met_short")||die;
while (<met>){
    chomp;
    my @a=split/\s/;
    next unless $a[2]>=5;
    my  $met=$a[3]/$a[2];
    my $met=(int(($met*5)));
    if ($met==5){$met=4}
    $met+=1;
    for my $i ($a[0]..$a[0]+$a[1]-1){
        $met{$i}=$met;
    }
}

my %coef;
for my $ex (1..4){
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/gen\ $ex\_penta");
    while (<IN>){
	chomp;
    my @a=split/\s/;
    my $lt="gen $ex";
	${${$coef{$lt}}{$a[0]}}{$a[1]}=int($a[2]*10000+0.5)/10000;
	print "${${$coef{$lt}}{$a[0]}}{$a[1]} $lt $a[0] $a[1]\n";
    }

    for my $di (0..2){
	open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/pro\ $ex\ $di\_penta");
	while (<IN>){
	    chomp;
	    my @a=split/\s/;
	    my $lt="pro $ex $di";
	    ${${$coef{$lt}}{$a[0]}}{$a[1]}=int($a[2]*10000+0.5)/10000}
    }
}

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/inter_penta");
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $lt="inter";
    ${${$coef{$lt}}{$a[0]}}{$a[1]}=int($a[2]*10000+0.5)/10000}

my %spat; my %spats;
open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/$chr\_mut_count_10kb")||die;
my $sp=<IN>; chomp($sp); my @sp=split/\s/,$sp;
while (<IN>){
    chomp;
    my @a=split/\s/;
    for my $i (1..$#a){${$spat{$a[0]}}{$sp[$i]}=$a[$i]}
}

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/$chr\_site_count_10kb")||die;
my $sp=<IN>; chomp($sp); my @sp=split/\s/,$sp;
while (<IN>){
    chomp;
    my @a=split/\s/;
    for my $i (1..$#a){${$spats{$a[0]}}{$sp[$i]}=$a[$i]}
}

my %aves;
my %ave;
for my $ch (1..22){
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/$ch\_mut_count_10kb")||die;
    my $sp=<IN>; chomp($sp); my @sp=split/\s/,$sp;
    while (<IN>){
	chomp;
	my @a=split/\s/;
	for my $i (1..$#a){$ave{$sp[$i]}+=$a[$i]}
    }

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/data/$ch\_site_count_10kb")||die;
    my $sp=<IN>; chomp($sp); my @sp=split/\s/,$sp;
    while (<IN>){
	chomp;
	my @a=split/\s/;
	for my $i (1..$#a){$aves{$sp[$i]}+=$a[$i]}
    }
}

my %rav;
foreach my $tri (sort keys %aves){
#    print "$tri  ($aves{$tri}-$ave{$tri})/$aves{$tri}\n";
    my $r=($aves{$tri}-$ave{$tri})/$aves{$tri};
    $rav{$tri}=$r;
}



my %spat5; my %spats5;
foreach my $pos (sort keys %spats){

    for my $ii (-2..2){
	my $pp=$pos+$ii*$sc1;
	next unless $spats{$pp};
	for my $i (1..$#sp){${$spats5{$pos}}{$sp[$i]}+=${$spats{$pp}}{$sp[$i]};
			   ${$spats5{$pos}}{'X'}+=${$spats{$pp}}{$sp[$i]};
			   ${$spat5{$pos}}{$sp[$i]}+=${$spat{$pp}}{$sp[$i]};
                           ${$spat5{$pos}}{'X'}+=${$spat{$pp}}{$sp[$i]};
	}
   }
#    print "$pos ${$spats5{$pos}}{'X'}\n";

}




open (IN,"/net/home/vseplyarsky/TOPMEDs/data/introns_plus_v4");
#my %int1;
my %cord_qg;
while (<IN>){
    chomp;
    my @a=split/\s/;
    next unless $a[0] == $chr;
#    print "$exp{$a[4]} $a[4]\n";
    next unless $exp{$a[4]} || $exp{$a[4]}==0;
    my $st=int(($a[1]+$sc/2)/$sc)*$sc;
    my $end=int(($a[2]+$sc/2)/$sc)*$sc;
    my $b=0;
 #   print "$_\n";    
    for my $bin (1..4){
	if ($bin{$bin}<=$exp{$a[4]}){$b=$bin}
    }
    my $len=int(($end-$st)/100);
    next unless $len;
    for my $i (0..$len) {
	my $pos=$st+$i*$sc;
	my $c=int(($i-1)/10);
	my $c1=3;
	if ($c<10){$c1=1;}
        if ($c<30){$c1=2;}

	$cord_qg{$pos}="gen $b $a[3]";
    }
    if ($a[3] eq '-'){
	for my $i (1..50) {
	    my $pos=$end+$i*$sc;
	    my $c=int(($i-1)/17);
	    $cord_qg{$pos}="pro $b $a[3] $c";
#	    print "$pos prom $b $a[3] $c $exp{$a[4]}\n";
	}}

    if ($a[3] eq '+'){
        for my $i (1..50) {
            my $pos=$st-$i*$sc;
            my $c=int(($i-1)/17);
            $cord_qg{$pos}="pro $b $a[3] $c";
        }}
}


my $F;
open (IN1,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_qc_dp_v2");
my %filter; 
    while (<IN1>){
	chomp;
	my @a=split/\s/;
	my $a=$a[1]/$a[3];
	my $lq='low';
	if($a>1 && $a[3]>=3){$lq='high';}
	$filter{$a[0]}=$lq;
	$F+=100;
}
#print "$F\n";


my %excl;
my $F1;
open (IN,"/net/home/vseplyarsky/TOPMEDs/data/knownGene.txt")||die;
my %pass; my $all_p;
while (<IN>){
    chomp;
    my @a=split/\t/;
    my $c=substr($a[1],3);
    next unless $c ==$chr;
    my @b=split/\./,$a[0]; 
   my $stc=$a[5]+1;
    my $endc=$a[6];

  #  print "${$cor{$va}}{'7'} $a[8] $a[9]\n";
#    my @st= split/\s/,${$cor{$va}}{'7'};
my @stc=split/\,/,$a[8];
my @endc=split/\,/,$a[9];
my $ll;
for my $s (0..$#endc){
    $stc[$s]+=1;
    next if $stc>=$endc[$s];
    next if  $endc<=$stc[$s];
    $stc[$s]=max($stc[$s],$stc);
    $endc[$s]=min($endc[$s],$endc);
  #  my $d=$stc[$s]-int($stc[$s]/$sc)*$sc;
  #  if ($d>20){my $p=int($stc[$s]/$sc)*$sc-100}
 #   $excl{$p}=1;
    my $f=int($stc[$s]/$sc);
    my $l=int($endc[$s]/$sc);
#    print "$f $l $s $F1    DFDF $_\n";
    for my $p ($f..$l){my $pos=$p*$sc;  if (!($excl{$pos})){$F1+=100;};$excl{$pos}=1;}

}}
#print "$F1\n";


sub rc {
    my $str = $_[0];
    # reverse complement
    $str =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $str;
}


sub rvc{
    my $str = $_[0];
    my @a=split/\_/,$str;
    my $rv=&rc($a[0]);
    my $rd=&rc($a[1]);
    my $out="$rv".'_'."$rd";
}



my %hg19;
open (hg,"/net/home/vseplyarsky/hg38/chr$chr.line")||die;
<hg>;
my $hg19=<hg>; chomp ($hg19);$hg19{$chr}=uc($hg19);


my %excp;
open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_common.pol")||die;
my %pol;
while (<pol>){
    my @a=split/\s/;
    ${$excp{$a[0]}}{$a[1]}=1;
}


open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol")||die;
my %pol;
while (<pol>){
    my @a=split/\s/;
${$pol{$a[0]}}{$a[1]}=1;
}



my (%regm, %regs, %out, %site);
my @xr=sort {$a<=>$b} keys %filter;
my $P;
#open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_input");

foreach my $pos (sort {$a<=>$b} keys %filter){
    my $lq=$filter{$pos};
 #   next unless $pos >47000000;  
    #next if $excl{$pos};
    my $PP=int($pos/$sc1)*$sc1;
    #next unless ${$spats5{$PP}}{'X'}>75000; 

   for my $p ($pos..$pos+99){

my $shar=$p-7;
my $sseq=substr($hg19{$chr},$shar,13);
next if $sseq=~/N/;
my $si='+';
#print "$type\n";
my $type='inter';

if ($cord_qg{$pos}){$type=$cord_qg{$pos}; $si=substr($type,6,1);substr($type,5,2)='';};
#print "$type AAA$si SI\n";

my $tri=substr($sseq,4,5);
#    print "$_ $a[3] $sseq\n";
my $trit=substr($sseq,5,3);

substr($sseq,4,5)='';
#my $tri=$a[3];
my $anc=substr($tri,2,1);
my $tritr= $trit;

if ($anc ne 'T' && $anc ne 'C'){$tritr=&rvc($trit);}
my $pair=substr($tritr,1,2);
my $CG;
if ($pair eq 'CG'){$CG=$met{$p}; if (!($met{$p})){$CG=4}}
my $pos1=int($pos/$sc1)*$sc1;
#${$regs{$pos1}}{$tritr}++;
my $dir=$dir{$pos1};
if (!($dir)){$dir=4}

for my $nuc ('A', 'C','G','T'){
next if $nuc eq $anc; 

my $mut="$anc\_$nuc";
my $in="$tri\_$nuc";
my $int="$trit\_$nuc";
#next if ${$excp{$p}}{$mut};

my $pol;
if (${$pol{$p}}{$mut}){$pol=1;}

if ($anc ne 'T' && $anc ne 'C'){$int=&rvc($int)}
#${$regm{$pos1}}{$int}+=$pol;
#${$regs{$pos1}}{$int}++;
#print "$PP $int $pos (((${$spats5{$PP}}{$int}-1)-(${$spat5{$PP}}{$int}-$pol)-0.5)/(${$spats5{$PP}}{$int}-1))/$rav{$tri}\n";
my $sp_par=1;
if (${$spats5{$PP}}{$int}>=10){(((${$spats5{$PP}}{$int}-1)-(${$spat5{$PP}}{$int}-$pol))/(${$spats5{$PP}}{$int}-1))/$rav{$int};}
if ($sp_par<=0){$sp_par=1}

$sp_par=int(($sp_par)*10000+0.5)/10000;
#$sp_par=log($sp_par);
if ($type eq 'inter' && ($anc ne 'T' && $anc ne 'C')){$in=&rvc($in);$sseq=&rc($sseq); if ($dir){$dir=8-$dir;}}
if ($si eq '-'){$in=&rvc($in);$sseq=&rc($sseq);if ($dir){$dir=8-$dir;}}

my $type1="$type met$CG";
#if ($type=~/pro/){    $type1="$type $CG";}
$P+=$pol; 
#my @f=keys %{$pol{$p}};
#my $f;
#for my $i (0..$#f){$f.="$f[$i] ";}
#print "$p $in $mut $trit $type1 $int FF $f P $pol $xr[-1] pop $P\n";
my $coeff;
for my $l (1..8){
    my $ls=$l-1;
    my $ll=substr($sseq,$ls,1);
    my $lo="$ll$l";
#    if ($type=~/pro/){print "$type\XXX\n${${$coef{$type}}{$in}}{$lo}XXX XYX$type AAAA$in $lo\n";}
#    print "$type\XXX\n${${$coef{$type}}{$in}}{$lo}XXX XYX$type AAAA$in $lo\n";
    my $x=1;
    if (${${$coef{$type}}{$in}}{$lo}){$x=${${$coef{$type}}{$in}}{$lo}};
    $coeff.="$x\t";
#    ${${$out{$type1}}{$in}}{$lo}+=$pol;
#    ${${$site{$type1}}{$in}}{$lo}++;
   # print "${${$out{$type1}}{$in}}{$lo} ${${$site{$type1}}{$in}}{$lo} $p $in $mut $trit $type1 $int\n";
}

#   $out{$in}++;

my $o=1-$pol;
#open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_input");
print OUT "$p\t$o\t$in $type1\t$coeff\RT_dir$dir\t$sp_par\t$lq\n";
}
   }}
