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
  



    open (fac,"/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/tRNAs.txt")||die;
    while (<fac>){
        chomp;

        my @a=split/\t/;
        substr($a[1],0,3)='';
#        print "$a[0] $a[1] $a[2] $a[9]  GGGG $_\n\n";
#       print "$a[9] $a[0]==$chr ($a[2]+$a[1])/2\n";
#       next unless $a[9] eq 'PLS,CTCF-bound';
        my $mean=int(($a[2]+$a[3])/2);
        next unless $a[1]==$chr;
        my $c=1; if ($a[6] eq '-'){$c=-1;}

#	print "$a[6] $mean $a[1] $c\n";
        for my $j (-1000..1000){
            my $i=$j*$c+$mean;
            $fac{$i}.="$j ";
        }
    }







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


    my %hg19;
    open (hg,"/net/home/vseplyarsky/hg38/chr$chr.line")||die;
    <hg>;
    my $hg19=<hg>; chomp ($hg19);
    $hg19=uc($hg19);



#    print "P2\n";

    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate_v3")||die;
 #   open (OUT1,">/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/$chr\_rand")||die;

   while (<IN>){
	chomp;
	next if /low/;
	my @a=split/\s/;
        next unless $fac{$a[0]};

	my $mu=$a[1];
        substr($mu,0,2)='';
        substr($mu,1,2)='';
#       my $anc=substr($mu,0,1);
	my $lc=$a[0]-1;
        my $n=substr($hg19,$lc,1);
        if ($n ne 'T' && $n ne 'C'){$mu=&rvc($mu);}
        next if (${$exep{$a[0]}}{$mu});
	my $pol;
        if (${$pol{$a[0]}}{$mu}){$pol=1;}
	my @fac=split/\s/,$fac{$a[0]};
        foreach my $f (@fac){
            ${$oute{$f}}{$mu}+=(1-$a[2]);
            ${$out{$f}}{$mu}+=$pol;
            ${$outs{$f}}{$mu}++;
	}
   }


open (OUT,">/net/home/vseplyarsky/GNOMAD_Tufts/data/CHIPseq/data/hash/$chr\_hash_tRNA");
    foreach my $val (sort keys %outs){
	foreach my $mu (sort keys %{$outs{$val}}){
	    
    if (!(${$out{$val}}{$mu})){${$out{$val}}{$mu}=0}
    if (!(${$oute{$val}}{$mu})){${$oute{$val}}{$mu}=0}
    
    print OUT "$val $mu ${$outs{$val}}{$mu} ${$oute{$val}}{$mu} ${$out{$val}}{$mu}\n"}
	
    }
    
}
