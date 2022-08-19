use strict;
#my $chr=$ARGV[0];
my (%sum1,%sum2,%sum3);
for my $chr (1..22){
    open (IN,"gunzip -c ../mut_model/for_LASSO/$chr\_hash.gz|")||next;
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $v=$a[1];
#    $a[0]='xxx';
    $sum1{$a[0]}=$v;
    ${$sum2{$a[0]}}{$v}+=$a[2];
    ${$sum3{$a[0]}}{$v}+=$a[3];

}}
#print "$sum3{TTTTT_A}\n";
foreach my $mm (sort keys %sum2){
	        foreach my $v (sort {$a<=>$b} keys %{$sum2{$mm}}){

  next unless (${$sum2{$mm}}{$v}+${$sum3{$mm}}{$v})>1000;
    my $r=${$sum3{$mm}}{$v}/${$sum2{$mm}}{$v};
#  print "$mm  $r=$sum3{$mm}/$sum2{$mm}\n";
  my $c1=$v;
  my $r1v=-log($v);
  my $rr=-log($r);

  #print "$mm  $r=${$sum3{$mm}}{$v}/${$sum2{$mm}}{$v} $c1 =$v\n";
   
 my $r1=($r1v+0.01)/($rr+0.01);

  my $co=1;
#    if (($r1>=1.5 || $r1<=1/1.5)){print"$mm $v  $r ${$sum3{$mm}}{$v}/${$sum2{$mm}}{$v} $r1 ($r1v)/$rr;\n"}
  print"$mm $v $r ${$sum3{$mm}}{$v}/${$sum2{$mm}}{$v} $r1 ($r1v)/$rr;\n";
 # if ($v<=0.999){print"$mm $v  $r ${$sum3{$mm}}{$v}/${$sum2{$mm}}{$v} $r1 ($r1v)/$rr;\n"}

	     }
}
