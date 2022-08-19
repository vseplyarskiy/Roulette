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
    my %den;
    my %tri; my %tris;
    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_freq.pol")||die;
    while (<pol>){
        chomp;
	my @a=split/\s/; 
	if ($a[2]<=0.0001){${$pol{$a[0]}}{$a[1]}=1}
        if ($a[2]>0.0001){${$exep{$a[0]}}{$a[1]}=1}
    }

    open (pol,"/net/home/vseplyarsky/GCPR/mutational_model/ageeffect/data/all_mut")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
        next unless $a[0] eq $chr;
        my $mu=$a[2];
        substr($mu,0,1)='';
        substr($mu,1,1)='';
        ${$den{$a[1]}}{$mu}=1;
}

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
	my $val=$a[3];
	if ($a[5]){$val=$a[6]}
	my $mu=$a[1];
	substr($mu,0,2)='';
        substr($mu,1,2)='';
	my $mu1=&rvc($mu);

        my $tr=$a[1];
        substr($tr,0,1)='';
        substr($tr,3,1)='';


	next if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1});# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
	$val=int(($tri{$tr}/$tris{$tr})*10000+0.5)/10000;
	my $pol=0;

#	$out{$val}++;
	if ( ${$pol{$a[0]}}{$mu} || ${$pol{$a[0]}}{$mu1}){$pol=1;}

	my $den;
        if ( ${$den{$a[0]}}{$mu} || ${$den{$a[0]}}{$mu1}){$den=1;}
	${$outs{$val}}{$pol}++;
        ${$out{$val}}{$pol}+=$den;



#print "$_ $pol $val $exep{$a[0]} "
#	my $r=rand();
#	if ($r<1/2000){print OUT1 "$_ $pol\n"}
   }}

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/variance_tri_$CHR");
foreach my $val (sort keys %outs){
    if (!(${$outs{$val}}{'0'})){${$outs{$val}}{'0'}=0}
    if (!(${$outs{$val}}{'1'})){${$outs{$val}}{'1'}=0}

    if (!(${$out{$val}}{'0'})){${$out{$val}}{'0'}=0}
    if (!(${$out{$val}}{'1'})){${$out{$val}}{'1'}=0}

    print OUT "$val ${$out{$val}}{'0'} ${$outs{$val}}{'0'} ${$out{$val}}{'1'} ${$outs{$val}}{'1'}\n"}

