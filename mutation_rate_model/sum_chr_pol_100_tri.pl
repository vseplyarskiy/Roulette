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
my %f;

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
	${$pol{$a[0]}}{$a[1]}=1;
    }

    my %hg19;
    open (hg,"/net/home/vseplyarsky/hg38/chr$chr.line")||die;
    <hg>;
    my $hg19=<hg>; chomp ($hg19);
    $hg19=uc($hg19);

    my %tri;
    my %tris;

    open (tri,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/cont_distr")||die;
    while (<tri>){
        chomp;
        my @a=split/\s/;
	my @b=split/\//,$a[3];
	my $mu1=$a[0];
        substr($mu1,0,1)='';
        substr($mu1,3,1)='';

	$tri{$mu1}+=$b[1]-$b[0];
        $tris{$mu1}+=$b[1];
    }



    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v4")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next if /low/;
	my @a=split/\s/;
	my $lc=$a[0]-1;
#	my $val=$a[2];
	my $val=int($a[0]/100)*100;;
	my $mu=$a[1];
	substr($mu,0,2)='';
        substr($mu,1,2)='';
        my $mu1=$a[1];
        substr($mu1,0,1)='';
        substr($mu1,3,1)='';

#	my $anc=substr($mu,0,1);
        my $n=substr($hg19,$lc,1);
	if ($n ne 'T' && $n ne 'C'){$mu=&rvc($mu);}
	next if (${$exep{$a[0]}}{$mu});# eq $mu[0]||$exep{$a[0]} eq $mu[1]);



	my $pol=0;
        if ( ${$pol{$a[0]}}{$mu}){$pol=1;}
#	print "$_ $n $pol\n";

        $out{$val}+=$pol;
        $outs{$val}+=int(($tri{$mu1}/$tris{$mu1})*10000+0.5)/10000;



#print "$_ $pol $val $exep{$a[0]} "
#	my $r=rand();
#	if ($r<1/2000){print OUT1 "$_ $pol\n"}
   }}

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/hash\_100_$CHR\_tri");

print OUT "pos all_e all\n";


foreach my $val (keys %outs){
print OUT "$val $outs{$val} $out{$val}\n";
}  
 



