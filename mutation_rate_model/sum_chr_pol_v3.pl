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

    open (IN,"gunzip -c /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate.gz|")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next if /low/;
	my @a=split/\s/;
	my $val=(int($a[2]*1000 +0.5))/1000;;
	my @mu=split/\//,$a[1];
	next if (${$exep{$a[0]}}{$mu[0]} || ${$exep{$a[0]}}{$mu[1]});# eq $mu[0]||$exep{$a[0]} eq $mu[1]);
	$outs{$val}++;
	my $pol=1;
	$out{$val}++;
	if ( ${$pol{$a[0]}}{$mu[0]} || ${$pol{$a[0]}}{$mu[1]}){$out{$val}--;$pol=0;}
#print "$_ $pol $val $exep{$a[0]} "
	my $r=rand();
	if ($r<1/2000){print OUT1 "$_ $pol\n"}
   }}

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/hash\_$CHR");
foreach my $val (sort keys %outs){
    if (!($out{$val})){$out{$val}=0}
    print OUT "$val $outs{$val} $out{$val}\n"}

