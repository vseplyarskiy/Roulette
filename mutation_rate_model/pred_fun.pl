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
    if ($a[1] eq 'independ'){$pr=log($a[0])+$a[11]; for my $i (2..9){$pr+=$a[$i]}; $pr=exp($pr)} 


    return ($pr);



}


open (coef,"test_AATTT_A_gen_4_met.LL")||die;
my %coef;
while (<coef>){
    chomp;
    my @a=split/\t/;
    $coef{$a[0]}="$a[1]";
    for my $i (2..$#a){ $coef{$a[0]}.="\t$a[$i]";}
}



open (IN,"test_AATTT_A_gen_4_met")||die;
my $h;
while (<IN>){
    chomp;
    $h++;
    my @a=split/\t/;
    my $out=$a[2];
    $out =~ s/\s/\_/g;
 
    my $in=$coef{$out};
    my @pr=split/\t/,$in;

 #print "$in $#pr Coeff\n";
    for my $i (3..12){if ($i==11){$in.="\t$a[$i]";next;};if ($a[$i]<0.01){$a[$i]=1};my $x=log($a[$i]);$in.="\t$x"}
    my @pr=split/\t/,$in;
#   print "$#pr len\n";
    my $v=&pred(@pr);
    print "$v $h\n";
}
