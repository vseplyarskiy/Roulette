use strict;

my (%s1,%s2,%s3,%s4);
for my $chr (1..22){
#    next unless $chr ==2;
    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/variance_v4_$chr");
#    open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/validation/variance_tri_$chr");

    while (<IN>){
	my @a=split/\s/;
	if ($a[0]>=0.03){$a[0]=int($a[0]*100+0.5)/100;}
	$s1{$a[0]}+=$a[1];
        $s2{$a[0]}+=$a[2];
        $s3{$a[0]}+=$a[3];
        $s4{$a[0]}+=$a[4];
print ""


    }}


foreach my $r (sort {$a<=>$b} keys %s2){
    next unless $s4{$r} && $s2{$r};
    my $r1= $s1{$r}/$s2{$r};
    my $r2= $s3{$r}/$s4{$r};
    my $r3='NA'; if ($r2){$r3=$r1/$r2;}
  print "$r $r1 $r2 $r3 $s1{$r}/$s2{$r} $s3{$r}/$s4{$r}\n";
}
