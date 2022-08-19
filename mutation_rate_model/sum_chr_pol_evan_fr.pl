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

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;
    print "GO\n";
   while (<IN>){
	chomp;
	next if /low/;next if /SFS_bump/;
	my @a=split/\s/;
	my $val=$a[3];
	if ($a[5]){$val=$a[6]};

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
	${$out{$val}}{$pol}++;

#print "$_ $pol $val $exep{$a[0]} "
#	my $r=rand();
#	if ($r<1/2000){print OUT1 "$_ $pol\n"}
   }}

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/hash\_evan_$CHR");
foreach my $val (sort keys %out){
    my $out="$val";
    for my $b (0..6){
	if (!(${$out{$val}}{$b})){${$out{$val}}{$b}=0}
	$out.=" ${$out{$val}}{$b}";
}


    print OUT "$out\n"}

