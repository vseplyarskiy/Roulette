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
        my $c=1; if ($a[5] eq '-'){$c=-1;}
#	print "$mean $c $a[5] $a[17]\n";
#	print "$a[6] $mean $a[1] $c\n";
        for my $j ($a[1]..$a[2]){
	    my $i=$mean+$j;
	    $fac{$j}=$rnu{$a[17]};
        }
    }







    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v5")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next unless /high/;
	my @a=split/\s/;
        next unless $fac{$a[0]};
	my $pol;
	my $mu=$a[1];
	substr($mu,0,2)='';
        substr($mu,1,2)='';
	my $mu1=&rvc($mu);
	next if (${$exep{$a[0]}}{$mu} || ${$exep{$a[0]}}{$mu1});
        if (${$pol{$a[0]}}{$mu} || ${$pol{$a[0]}}{$mu1}){$pol=1;}
	my $f=$fac{$a[0]};
	$oute{$f}+=$a[3];
	$outs{$f}++;
	$out{$f}+=$pol;
   }


open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/hash/$chr\_RNU_ind_denovo");

    foreach my $val (sort keys %outs){	    
    if (!($out{$val})){$out{$val}=0}
    if (!($oute{$val})){$oute{$val}=0}
    
    print OUT "$val $outs{$val} $oute{$val} $out{$val}\n"}
	
    
    
}
