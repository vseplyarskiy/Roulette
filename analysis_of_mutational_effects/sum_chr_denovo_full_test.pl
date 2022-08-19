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
    my $sum;
    my %pol;
    my %exep;


    open (pol,"/net/home/vseplyarsky/GCPR/mutational_model/ageeffect/data/all_mut")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
        next unless $a[0] eq $chr;
        my $mu=$a[2];
        substr($mu,0,1)='';
        substr($mu,1,1)='';
        ${$pol{$a[1]}}{$mu}++;
    }



    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Cistrome/by_chr_v2/$chr\_peaks_bed")||die;
    while (<fac>){
        chomp;

        my @a=split/\s/;
#	next unless $id{$a[0]};
     #   substr($a[0],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#	print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
   #     next unless $a[1]==$chr;
	my $na="$a[1] $a[2]";
	for my $j (-10..10){
            my $i=int(($a[0]+$wi/2)/$wi)*$wi+$j*$wi;
	    my $j1=abs($j*20);
#	    print "$fac{$i} $a[2] eq 'Testis'\n";
	    my $o=$fac{$i};	 
	    if ($fac{$i}){if  ($a[2] eq 'Testis'){$fac{$i}="$a[2] $j1";}}
            else{$fac{$i}="$a[2] $j1";}
#	    print "{$fac{$i}\n";
	    if ($o){ print "$o $fac{$i}\n";}

        }
    }
 #   foreach my $k (keys %fac){print "$fac{$k}\n";    }

#    print "P2\n";

 #   open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5_filter_ns.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
 #   open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/$chr\_denovo_TFBS");

    my %out; my %site; my %exp;
   while (<IN>){
	chomp;
	next unless /high/;
	my @a=split/\s/;
	my $pos=int(($a[0]+$wi/2)/$wi)*$wi;
        next unless $fac{$pos};
	my $fac=$fac{$pos};
	my $mu2=$a[1];
        my $mu=$a[1];
        substr($mu,0,2)='';
        substr($mu,1,2)='';
        my $mu1=&rvc($mu);
	substr($mu2,0,1)='';
        substr($mu2,3,1)='';
	if ($mu2=~/CG/){substr($mu2,0,1)='';} 
	else {substr($mu2,0,1)=''; substr($mu2,1,1)='';}
        my $pol=0;

#        if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1}){$pol='NA'};# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
#        my $pol=1;
	if ( ${$pol{$a[0]}}{$mu1} || ${$pol{$a[0]}}{$mu}){$pol=1;}
	${$out{$fac}}{$mu2}+=$pol;
	${$site{$fac}}{$mu2}++;
        ${$exp{$fac}}{$mu2}+=$a[-2];

#	print OUT "$mu2 $a[3] $fac{$pos} $pol\n"}
   }

    foreach my $fac (sort keys %site){
	foreach my $mu (sort keys %{$site{$fac}}){
	    print OUT "$fac $mu ${$exp{$fac}}{$mu} ${$out{$fac}}{$mu} ${$site{$fac}}{$mu}\n";
}

    }}
