use strict;



my %site;
my %siten;
my %site1;
my %siten1;


for my $chr (1..22){
    open (IN,"../mut_model/validation/hash_$chr");

while (<IN>){
    chomp;
    my @a=split/\s/;
    my $v=int($a[0]*100+0.5)/100;
    $site{$v}+=$a[1];
    $siten{$v}+=$a[2];
}


open (IN,"../mut_model/validation/hash_v2_$chr");
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $v=$a[0];
    $site1{$v}+=$a[1];
    $siten1{$v}+=$a[2];
}
}

foreach my $r (sort {$a<=>$b} keys %siten1){

    my $v=$r;
   
#    print "$v $a[2]/$a[1] $siten{$v}/$site{$v}\n";
    my $r1=$siten{$v}/$site{$v};
    my $r=$siten1{$v}/$site1{$v};

print "$v $r $r1 $siten1{$v}/$site1{$v} $siten{$v}/$site{$v}\n";
} 

