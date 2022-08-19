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


my %bin;
$bin{'1'}=0.00005;
$bin{'2'}=0.0005;
$bin{'3'}=0.005;
$bin{'4'}=0.05;
$bin{'5'}=0.2;


my $CHR=$ARGV[0];
my %out;
my %outs;
my %f;

for my $chr (1..22){
    next unless $chr == $CHR;
    my $sum;
    my %pol;
    my %exep;
  
    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_freq.pol")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
	my $bin=6;
	for my $b (1..5){
	    my $inv=6-$b;
	    if ($a[2]<$bin{$inv}){$bin=$inv}
	}
	${$pol{$a[0]}}{$a[1]}=$bin;
#	print "$_ $bin\n";
    }

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v4")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
    print "GO\n";
   while (<IN>){
	chomp;
	next if /low/;
	my @a=split/\s/;
	my $val=int(-log($a[2])*1000+0.5)/1000;

	if ($a[5]){$val=int(-log($a[5])*1000+0.5)/1000}
	my $mu=$a[1];
	substr($mu,0,2)='';
        substr($mu,1,2)='';
	my $mu1=&rvc($mu);
	my @mu;
	$mu[0]=$mu;
        $mu[1]=$mu1;

#	next if (${$pol{$a[0]}}{$mu[0]}==0 || ${$pol{$a[0]}}{$mu[1]}==0);# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
	#$outs{$val}++;
	my $pol=0;


	if ( ${$pol{$a[0]}}{$mu[0]}){$pol=${$pol{$a[0]}}{$mu[0]};}
        if ( ${$pol{$a[0]}}{$mu[1]}){$pol=${$pol{$a[0]}}{$mu[1]};}
#	print "$pol $val ${$out{$val}}{$pol} X $_\n";
	${${$out{$a[1]}}{$val}}{$pol}++;

#print "$_ $pol $val $exep{$a[0]} "
#	my $r=rand();
#	if ($r<1/2000){print OUT1 "$_ $pol\n"}
   }}

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/hash\_evan_penta_$CHR");
foreach my $cont (sort keys %out){
    foreach my $val (sort keys %{$out{$cont}}){
    my $out="$cont $val";
    for my $b (0..6){
	if (!(${${$out{$cont}}{$val}}{$b})){${${$out{$cont}}{$val}}{$b}=0}
	$out.=" ${${$out{$cont}}{$val}}{$b}";
}


    print OUT "$out\n"}

}
