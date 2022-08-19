use strict;

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


my %sol;

my %bins;

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/bins")||die;
while (<IN>){

    chomp;
    my @a=split/\s/;
    $bins{$a[0]}=1;
}

my @bin = sort {$a<=>$b} keys %bins;


open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/for_reg_all_TFBS_prom_solution");
while (<IN>){
    chomp;
    my @a=split/\t/;
    my $in="$a[0] $a[1] $a[2] $a[3] $a[4]";
    my $r=$a[-1]/$a[-2];
    my $lam=-log($r);

   my $m=5; my $ci; for my $i1 (0..$#bin){my $lm=abs($bin[$i1]-$lam);
				  if ($m>$lm){$ci=$i1; $m=$lm}}


    $r=int($r*1000+0.5)/1000;
    $lam=int($lam*1000+0.5)/1000;
#    print "$_ IN $in V $bin[$ci] Lam $lam\n";
    $sol{$in}="$bin[$ci]";
#    print "$_ VVVV $in XXXX $r $lam\n";
}

my $CHR=$ARGV[0];
my %out;
my %outs;
my %oute;

my %f;
my %fac;
my $wi=20;

for my $chr (1..22){
    next unless $chr == $CHR;
    my $chrf="chr$chr";
    my $sum;
    my %pol;
    my %exep;
    my %tss;
    my %RTder;
    my $pl; my $min;
    my %tssan;
    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/ncbiRefSeq.txt.gz|");
    while (<IN>){
        chomp;
        my @a=split/\t/;
        next unless $a[2] eq $chrf;
        next unless $a[-3]='cmpl';
        next unless $a[3] eq '-' || $a[3] eq '+';
        if ($a[3] eq '-' ){$a[4]=$a[5]}
        my $c=int($a[4]/10)*10;
        my $c1=-1;
        if ($a[3] eq '-'){$c1=1;};
        for my $i (1..200){my $cor=$c+$c1*$i*10; if ($tss{$cor}){$tss{$cor}='Mult';next}
                           $tss{$cor}='Sing'
        }
    }

    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Cistrome/by_chr_v2/$chr\_peaks_bed")||die;
    while (<fac>){
        chomp;
        my $test='other';
        if( /Testis/){$test='Testis';}
        my @a=split/\s/;
#       next unless $id{$a[0]};
     #   substr($a[0],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#       print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
   #     next unless $a[1]==$chr;
        #my $na="$a[1] $a[2]";
        for my $j (-50..50){
            my $i=int(($a[0]+$wi/2)/$wi)*$wi+$j*$wi;
            my $j1=abs($j*20);
#           print "$test\n";
	    if ($fac{$i}=~/other/){$fac{$i}="$j1 $test";
#                                 print "XXX\n";
				   next;}
            $fac{$i}="$j1 $test";
        }
    }


    my %hg19;
    open (hg,"/net/home/vseplyarsky/hg38/chr$chr.line")||die;
    <hg>;
    my $hg19=<hg>; chomp ($hg19);
    $hg19=uc($hg19);





#    print "P2\n";

    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5_filter_ns.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
    open OUT, '>:gzip',  "/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5.2_TFBS_correction.gz"; 
 #  open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/$chr\_for_reg_all_TFBS");

   while (<IN>){
	chomp;
	my @a=split/\s/;
	my $m=5;

	if (!($bins{$a[3]})){my $ci; for my $i1 (0..$#bin){my $lm=abs($bin[$i1]-$a[3]);
							  if ($m>$lm){$ci=$i1; $m=$lm}}

			    $a[3]=$bin[$ci];
	}
	
	

	if ($#a>4 && !($bins{$a[6]})){my $ci; for my $i1 (0..$#bin){my $lm=abs($bin[$i1]-$a[3]);
								   if ($m>$lm){$ci=$i1; $m=$lm}}
				      $a[5]=$bin[$ci];}
				     

	my $mu=$a[1];
	my $lc=$a[0]-1;
        my $n=substr($hg19,$lc,1);
        if ($n ne 'T' && $n ne 'C'){$mu=&rvc($mu);}
        my $in="$a[0] $mu"; for my $i (3..$#a){next if $i==2 || $i==6; $in.=" $a[$i]"}

        my $pos=int(($a[0]+$wi/2)/$wi)*$wi;
#	print "$_ $in\n";
	if (!(/high/) || !($fac{$pos})){print OUT "$in\n";next}
        my $mu2=$a[1];
        substr($mu2,0,1)='';
        substr($mu2,3,1)='';
        my $tssc=int($a[0]/10)*10;
        my $tssv='non-promoter';
	if ($tss{$tssc}){$tssv='promoter';}
	my $in1="$mu2 $a[3] $fac{$pos} $tss{$tssc}";
	print OUT "$a[0] $mu $a[3] TFBS $sol{$in1}\n";
#	print "$_ $a[0] $mu $a[3] TFBS $sol{$in1}\n";
#        print "$a[0] $a[1] $a[2] TFBS $sol{$in}\n";

#	if ($mu2=~/CG/){substr($mu2,0,1)='';} 
#	else {substr($mu2,0,1)=''; substr($mu2,1,1)='';}
     
#	print OUT "$mu2 $a[3] $fac{$pos} $pol\n";
#        print OUT "$a[0] $a[1] $a[2] TFBS $r $lam\n"


   }


}
