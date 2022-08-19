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

    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
	${$pol{$a[0]}}{$a[1]}++;
    }
    my %id;
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Histones/topID")||die;
    while (<IN>){
        chomp;
        my @a=split/\s/;
        $id{$a[0]}++;
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
	    if ($fac{$i}){$fac{$i}="Mult $a[2] $j1";next;}
            $fac{$i}="$a[1] $a[2] $j1";
        }
    }

#    print "P2\n";

    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5_filter_ns.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
    open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/DHS/results/$chr\_for_reg_all_TFBS");

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
        my $pol=1;

        if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1}){$pol='NA'};# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
#        my $pol=1;
	if ( ${$pol{$a[0]}}{$mu1} || ${$pol{$a[0]}}{$mu}){$pol=0;}
	print OUT "$mu2 $a[3] $fac{$pos} $pol\n"}
}
