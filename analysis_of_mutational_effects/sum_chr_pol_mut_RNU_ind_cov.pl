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
    my %rnu; 

    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/RNU_corrrect_list.txt")||die;
    while (<fac>){ 
	chomp;
	my @b=split/\,/;
	my @a=split/\s+/,$b[0];
#        print "$a[1]\n";
	next unless $a[1]=~/RNU/;
#	print "$a[1]\n"; 
	$rnu{$a[1]}="$a[1] gene";
    }
    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/RNU_pseudo.txt")||die;
    while (<fac>){
        chomp;
        my @b=split/\,/;
        my @a=split/\s+/,$b[0];
#        print "$a[1]\n";
        next unless $a[1]=~/RNU/;
#       print "$a[1]\n";
        $rnu{$a[1]}="$a[1] pseudogene";
    }





=pod
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
=cut


    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/gencodeV36.bed")||die;
    while (<fac>){
        chomp;

        my @a=split/\t/;
#	print "$a[17]\n";
	next unless  $rnu{$a[17]}; 
	substr($a[0],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#       print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
#       next unless $a[9] eq 'PLS,CTCF-bound';
        my $mean=int(($a[1]+$a[2])/2);
        next unless $a[0]==$chr;
	my $cor=int($mean/100)*100;
	for my $i (-10..10){
	    my $j=$cor+$i*100;
	    $fac{$j}="$rnu{$a[17]} $i";
	}    
    }




    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/results/$chr\_qc_dp_v3")||die;
    while (<IN>){
        chomp;
	my @a=split/\s/;
	next unless $fac{$a[0]};
	my $cov=$a[2]/$a[3];
	my $an=$a[4]/$a[3];
	my $qc=$a[1]/$a[3];
	next if ($qc<1 || $an<140000||$cov<30);
        my $f=$fac{$a[0]};
        $outs{$f}++;
        $out{$f}+=$cov;
        $oute{$f}+=$qc;

    }




open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/hash/$chr\_RNU_ind_denovo_cov");

    foreach my $val (sort keys %outs){	    
    if (!($out{$val})){$out{$val}=0}
    if (!($oute{$val})){$oute{$val}=0}
    
    print OUT "$val $outs{$val} $out{$val} $oute{$val}\n"}
}
