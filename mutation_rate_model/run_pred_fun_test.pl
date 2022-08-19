use strict;




sub pred {
    my @a = @_;
    my $lv;
    if ($#a==25){
	$lv=$a[0];
	my %fac;$fac{'1'}=0;
	my $d=substr($a[25],-1);
	for my $i (2..7){my $c=9+15; my $j=7+$i;$fac{$i}=$a[$j];}
	for my $i(1..8){my $c=$i+15; $lv+=$a[$c]*$a[$i]}
	$lv+=$fac{$d};
	$lv+=$a[15]*$a[25];
#	print "$lv $a[15]*$a[25]\n";
    }

    if ($#a==115){
        $lv=$a[0];
        my %fac;$fac{'1'}=0;
        my $d=substr($a[114],-1);
        for my $i (2..7){my $c=9+105; my $j=7+$i;$fac{$i}=$a[$j];}
        for my $i(1..8){my $c=$i+105; $lv+=$a[$c]*$a[$i]}
        $lv+=$fac{$d};
        $lv+=$a[15]*$a[115];
#       print "$lv $a[15]*$a[25]\n";
	my $p=15;
for my $f (1..8){
    for my $s (1..10){
	next if $s<=$f;
	my $c1=105+$f;
	if ($s==9){my $p1=$p+$d-1;if (!($d==1)){$lv+=$a[$p1]*$a[$c1];}$p+=6;
#		   print "$a[114] $d $p $f $s $p1 $a[$p1]*$a[$c1] $a[$p]\n";
	}
	else{$p++;  my $c2=105+$s;$lv+=$a[$p]*$a[$c1]*$a[$c2]}
    }
}
	if (!($d==1)){my $p1=98+$d;$lv+=$a[$p1]*$a[115];}


}    


    my $pr=exp($lv)/(exp($lv)+1);
    if (!($lv)){$pr='NA'} 
    if ($a[1] eq 'independ'){$pr=log($a[0])+$a[11]; for my $i (2..9){$pr+=$a[$i]}; $pr=exp($pr); if ($pr>=1){ $pr=0.995}} 
    if ($#a==10){$pr=$a[0]}   
    return ($pr);

}


sub tr_name {
    my @a = @_;
#    my $na=$a[0];
    my $out;
    if ($a[1] eq 'pro'){my $n=&rvc($a[0]); $out="$n gen $a[2] $a[4]"}
    if ($a[1] eq 'gen'){my $an=substr($a[0],2,1);  my $n=$a[0];if ($an eq 'G' || $an eq 'A'){$n=&rvc($a[0]);}
			$out="$n inter $a[3]"}
    my @out=split/\s/,$out;
    return($out);
}


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
    return($out);
}


my $chr=$ARGV[0];

#open (OUT,"| gzip -c > /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/$chr\_rate.gz|");

my %coef;
for my $i (1..64){
open (cat,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/cat/nm.$i")||die;
#<cat>;<cat>;
while (<cat>){
    chomp;
    my @a=split/\t/;
    substr($a[0],-3)='';
    my $out="$a[1]";
    if ($#a>=2){
	for my $i (2..$#a){
	    $out.="\t$a[$i]";
	}}
   
    $coef{$a[0]}=$out;
#    print "$a[0] $out\n";
}
}




open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/AACGT_T_inter_met5_cl")||die;
my $h;
while (<IN>){
    chomp;
    $h++;
    my @a=split/\t/;
    my $out=$a[2];
    $out =~ s/\s/\_/g;
    my @out=split/\s/,$a[2]; 
    my $in=$coef{$out};
   # next unless $out=~/pro/;
#    print "UP $out $in\n";
    if ($in eq 'too_short' || !($in)){$out=&tr_name(@out);@out=split/\s/,$out;$out =~ s/\s/\_/g;
				     $in=$coef{$out};}
    
    if ($in eq 'too_short' || !($in)){$out=&tr_name(@out);@out=split/\s/,$out;$out =~ s/\s/\_/g;
				      $in=$coef{$out};}
#    print "DOWN $out $in\n";


    my @pr=split/\t/,$in; 
 # print "$in $#pr Coeff\n";
    for my $i (3..12){if ($i==11){$in.="\t$a[$i]";next;};if ($a[$i]<0.01){$a[$i]=1};my $x=log($a[$i]);$in.="\t$x"}
    my @pr=split/\t/,$in;
 #  print "$#pr len\n";
    my $v=&pred(@pr);
 #   print "$v $h\n\n\n";

    my $lam=-log($v);
    $v=int($v*10000 +0.5)/10000;
    $lam=int($lam*10000 +0.5)/10000;
    my $mut=substr($out,2,1);
    my $mutd=substr($out,6,1);
    $mut="$mut\_$mutd";
    my $mut1=&rvc($mut);
    print  "$a[0] $mut/$mut1 $v $lam $a[1]\n";
#    if (/met5/ || $v<0.3){print "$a[0] $mut/$mut1 $v $lam $a[13]\n";}


}
