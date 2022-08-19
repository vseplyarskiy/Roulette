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
my $wi=100;
for my $chr (1..22){
    next unless $chr == $CHR;
    my $sum;
    my %pol;
    my %exep;
    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_common.pol")||die;
    while (<pol>){
        chomp;
	my @a=split/\s/; 
	${$exep{$a[0]}}{$a[1]}=1;
    }

    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
	${$pol{$a[0]}}{$a[1]}++;
    }
    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Cistrome/by_chr_v2/$chr\_peaks_bed")||die;
    while (<fac>){
        chomp;

        my @a=split/\s/;
     #   substr($a[0],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#	print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
      #  next unless $a[0]==$chr;
	my $na="$a[1]\_$a[2]";
	for my $j (-10..10){
            my $i=int(($a[0]+$wi/2)/$wi)*$wi+$j*$wi;
            $fac{$i}.="$na\_$j ";
        }
    }

#    print "P2\n";

    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5_filter_ns.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next unless /high/;
	my @a=split/\s/;
	my $pos=int(($a[0]+$wi/2)/$wi)*$wi;
        next unless $fac{$pos};
	my $mu2=$a[1];
        my $mu=$a[1];
        substr($mu,0,2)='';
        substr($mu,1,2)='';
        my $mu1=&rvc($mu);
	substr($mu2,0,1)='';
        substr($mu2,3,1)='';
#	if ($mu2=~/CG/){substr($mu2,0,1)='';} 
#	else {substr($mu2,0,1)=''; substr($mu2,1,1)='';}
        next if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1});# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
        my $pol=0;
	if ( ${$pol{$a[0]}}{$mu1} || ${$pol{$a[0]}}{$mu}){$pol=1;}
	my @fac=split/\s/,$fac{$pos};
	foreach my $f (@fac){
	    ${$oute{$f}}{$mu2}+=(1-$a[2]);
            ${$out{$f}}{$mu2}+=$pol;
	    ${$outs{$f}}{$mu2}++;
#	    my $val=$f;   print "$val $mu ${$oute{$f}}{$mu}  ${$out{$f}}{$mu}  ${$outs{$f}}{$mu} $a[2] $pol $a[0]\n";
	}
   }

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/$chr\_hash_all_TFBS_C13_14");
    foreach my $val (sort keys %outs){
	foreach my $mu (sort keys %{$outs{$val}}){
	    
    if (!(${$out{$val}}{$mu})){${$out{$val}}{$mu}=0}
    if (!(${$oute{$val}}{$mu})){${$oute{$val}}{$mu}=0}
    
    print OUT "$val $mu ${$outs{$val}}{$mu} ${$oute{$val}}{$mu} ${$out{$val}}{$mu}\n"}
	
    }
    
}
