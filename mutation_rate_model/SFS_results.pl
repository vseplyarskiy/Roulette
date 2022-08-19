use strict;

open (IN,"/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/SFS/bad_result_11")||die;
my $cross; my $priv;
while (<IN>){
    chomp;
    my @a=split/\s/;
    my $l="$a[2] $a[3]";
    $l=length($l);
    next unless $l ==3;
    next unless ($a[4]+$a[5]+$a[6]+$a[7]+$a[8])>5 && ($a[4]+$a[5]+$a[6]+$a[7]+$a[8])<20;
    if ($a[6] && $a[7]){$cross++}
    if (($a[6] && (!($a[7])))||($a[7] && (!($a[6])))){$priv++}
}

my $r=$priv/$cross;
print "$r=$priv/$cross\n";
