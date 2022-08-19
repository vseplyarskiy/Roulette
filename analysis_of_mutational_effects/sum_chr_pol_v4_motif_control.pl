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
for my $chr (1..22){
    next unless $chr == $CHR;
    my $sum;
    my %pol;
    my %exep;
    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_common.pol")||die;
    while (<pol>){
        chomp;
	my @a=split/\s/; 
	$exep{$a[0]}=1;
    }

    open (pol,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr.pol")||die;
    while (<pol>){
        chomp;
        my @a=split/\s/;
	$pol{$a[0]}++;
    }
    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/Cistrome/by_chr/$chr\_motif");
    while (<fac>){
        chomp;
        my @a=split/\s/;
	$a[2]+=100;
        $a[1]+=100;

	my $l=$a[2]-$a[1]+1;
	for my $i ($a[1]..$a[2]){
	    my $p=$i-$a[1]+1;
	    if ($a[-1] eq '-'){$p=$l-$p+1}
        $fac{$i}.="$a[0]\_$p ";
	}
    }

#    print "P2\n";

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next unless /high/;
	my @a=split/\s/;
	next if ($exep{$a[0]});
	next unless $fac{$a[0]};
	my @fac=split/\s/,$fac{$a[0]};
	foreach my $f (@fac){
	    $oute{$f}+=(1-$a[2]);
            $out{$f}+=($pol{$a[0]}/3);
            $outs{$f}++;
#	    my $val=$f;   print "$val $outs{$val} $oute{$val} $out{$val} $a[2] $pol{$a[0]} $a[0]\n";
	}
   }

open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/hash/$chr\_Cistrome_control_hash");
foreach my $val (sort keys %outs){
    if (!($out{$val})){$out{$val}=0}
    if (!($oute{$val})){$oute{$val}=0}
  
  print OUT "$val $outs{$val} $oute{$val} $out{$val}\n"}

}
