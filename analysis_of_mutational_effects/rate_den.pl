use strict;

my $coef=0.000653909515334231;
my (%s1,%s2,%c1);
for my $chr (1..22){
#    next unless $chr ==2;
    open (IN,"hash/$chr\_hash_mut_denovo");
    while (<IN>){
	my @a=split/\s/;
	$s1{$a[0]}+=$a[2];
        $s2{$a[0]}+=$a[3];
	$c1{$a[0]}+=$a[1];
    }}


foreach my $r (sort {$s2{$b}/$s1{$b}<=>$s2{$a}/$s1{$a}} keys %s1){

    next unless $s1{$r}>10000;
    my $r1= $s2{$r}/$s1{$r};
    my $r3=int($s1{$r}*$coef +0.5);
    my @r=split/\./,$r;
    print "$r[1].$r[2] $r1 $s2{$r} $s1{$r} $c1{$r} observed\n";
    print "$r[1].$r[2] $coef $r3 $s1{$r} $c1{$r} expected\n";

}
