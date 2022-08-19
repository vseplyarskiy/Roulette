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
	next unless $a[2] >=0.00005 && $a[2] <=0.0001;
	
	${$pol{$a[0]}}{$a[1]}=$a[2];
#	print "$_ $bin\n";
    }

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v4")||die;
    open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/SFS/input_$CHR"); # 

    print "GO\n";
   while (<IN>){
	chomp;
	next if /low/;
	my @a=split/\s/;
	next unless $pol{$a[0]};
	my $val=int(-log($a[2])*1000+0.5)/1000;
	if ($a[5]){$val=int(-log($a[5])*1000+0.5)/1000}
	next unless $val > 0.05 && $val <0.15;
	my $mu=$a[1];
	substr($mu,0,2)='';
        substr($mu,1,2)='';
	my $mu1=&rvc($mu);
	my $fr=${$pol{$a[0]}}{$mu};my $M=$mu;
	if(${$pol{$a[0]}}{$mu1}){$fr=${$pol{$a[0]}}{$mu1};$M=$mu1;};
	next unless $fr;
	print OUT "$a[0] $a[1] $val $fr $M\n";
# 	${$pol{$a[0]}}{$mu[1]}==0);# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
   }}


