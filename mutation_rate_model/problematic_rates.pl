use strict;

my %out;
for my $r1 (2..100){
    for my $r2 (2..100){
	next if $r2>$r1; 
	my $rn1=$r1*0.01;
        my $rn2=$r2*0.01;
	if ($rn1==0.99){$rn1=0.9874;}
	if ($rn1==1){$rn1=0.996}

        if ($rn2==0.99){$rn2=0.9874;}
        if ($rn2==1){$rn2=0.996}
	my $rr=$rn2*$rn1;
	$rr=int($rr*10000+0.5)/10000;
	$out{$rr}.="$rn2,$rn1 ";
    }}

my %c;
open (IN,"coeff_bad");
while (<IN>){
    chomp;
    my @a=split/\s/;
    for my $i (0..$#a){
	$c{$a[$i]}=1;
    }

}
my %h;
foreach my $rr (sort {$a<=>$b} keys %out){
    next unless $c{$rr}; 
    my @a=split/\s/,$out{$rr};

    for my $i (0..$#a){
	my @b=split/\,/,$a[$i];
	for my $i (0..$#b){
	    $h{$b[$i]}++}}
   print "$rr $out{$rr}\n"; 
}

foreach my $h (sort {$h{$b}<=>$h{$a}} keys %h){
    print "$h $h{$h}\n";
}




