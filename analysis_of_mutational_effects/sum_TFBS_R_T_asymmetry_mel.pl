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
                           $tss{$cor}=$a[3]
        }
    }
#	$tss{$c}=$a[5];$tss{$c1}=$a[5]}
#    print "$pl $min\n";
    open (RT,"/net/home/vseplyarsky/TOPMEDs/tracks_to_correlate/signatures_positive2/intensities_rec_corrected.txt.mod.epi");
    while (<RT>){
    chomp;
    my @a=split/\s/;
#print "$a[-5]\n"; last; 
   next unless $a[0] eq $chr;
    next unless ($a[-5]>0.5 || $a[-5]< -0.5 );
    my $or='le';
    if ($a[-5]< -0.5){$or='la';$pl++;}
        $min++;
 
   $RTder{$a[1]}=$or;
    }
#    print "$pl $min\n";




    open (pol,"/net/home/vseplyarsky/all_icgc/icgc_all_data/results/SKCM/$chr\_non-clusters_cont.hg38")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
	my $mut=substr($a[-2],1,1);
	$mut="$mut\_$a[-1]";
#	print "$a[1] $mut\n";
        ${$pol{$a[1]}}{$mut}=1;
    }






    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Cistrome/by_chr_v2/$chr\_motif_mel")||die;
    while (<fac>){
        chomp;
	next unless /mel/;
        my @a=split/\s/;
#	next unless $id{$a[0]};
     #   substr($a[0],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#	print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
   #     next unless $a[1]==$chr;
#	my $na="$a[1] $a[2]";
	my $mean=($a[1]+$a[2])/2;
	for my $j (-50..50){
            my $i=int(($mean+$wi/2)/$wi)*$wi+$j*$wi;
	    my $j1=abs($j*20);
	   # if ($fac{$i}){$fac{$i}="$j1";next;}
            $fac{$i}="$j1";
        }
    }

#    print "P2\n";



    my %hg19;
    open (hg,"/net/home/vseplyarsky/hg38/chr$chr.line")||die;
    <hg>;
    my $hg19=<hg>; chomp ($hg19);
    $hg19=uc($hg19);







    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5_filter_ns.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
    open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/$chr\_TFBS_as_mel_tri");


    my (%site,%obs,%exp);

   while (<IN>){
	chomp;
	next unless /high/;
	my @a=split/\s/;
	my $pos=int(($a[0]+$wi/2)/$wi)*$wi;
        next unless $fac{$pos};
        my $lc=$a[0]-1;
        my $mu=$a[1];
        my $n=substr($hg19,$lc,1);
        if ($n ne 'T' && $n ne 'C'){$mu=&rvc($mu);}

	my $mu2=$mu;
      #  my $mu=$a[1];
        substr($mu,0,2)='';
        substr($mu,1,2)='';
        my $mu1=&rvc($mu);
#	next if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1});
	substr($mu2,0,1)='';
        substr($mu2,3,1)='';
	my $rtc=int($a[0]/10000)*10000;
	my $tssc=int($a[0]/10)*10;
	my $tssv='NA';
        my $RTv='NA';
	if ($tss{$tssc}){$tssv=$tss{$tssc}};
        if ($RTder{$rtc}){$RTv=$RTder{$rtc}};
	my $out="$mu2 $fac{$pos} $tssv $RTv";

 
#	if ($mu2=~/CG/){$out="$mu2 $fac{$pos} $tssv $RTv";} 
#	else {substr($mu2,0,1)=''; substr($mu2,1,1)='';}
        my $pol=0;

 #       if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1}){$pol='NA'};# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
#        my $pol=1;
	if ( ${$pol{$a[0]}}{$mu1} || ${$pol{$a[0]}}{$mu}){$pol=1;}
	$site{$out}++;
        $exp{$out}+=1-$a[2];
        $obs{$out}+=$pol;
}


    foreach my $out (sort keys %site){
	print OUT "$out $obs{$out} $exp{$out} $site{$out}\n";

    }
}
