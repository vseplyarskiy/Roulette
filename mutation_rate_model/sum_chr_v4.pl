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
my %f;

for my $chr (1..22){
    next unless $chr == $CHR;
    my $sum;
    my %pol;
    open (pol,"/net/home/vseplyarsky/GCPR/mutational_model/ageeffect/data/all_mut")||die;
    while (<pol>){
        chomp;
	my @a=split/\s/;	
	next unless $a[0] == $CHR;
	my $tri=&rvc($a[2]);
	${$pol{$a[1]}}{$a[2]}=1;
        ${$pol{$a[1]}}{$tri}=1;
	$f{'1'}++;
#	print "$_    ${$pol{$a[1]}}{$tri} $a[1] $tri\n";
    }

    open (IN,"/net/home/vseplyarsky/GCPR/mutational_model/script/results/extended_simple_sv3\_$CHR")||die;
 
   while (<IN>){
	chomp;
	my @a=split/\s/;
	my $val=int($a[2]/100);
	${$outs{$a[1]}}{$val}++;
	if ( ${$pol{$a[0]}}{$a[1]}){my $f=${$pol{$a[0]}}{$a[1]};
				    ${${$out{$a[1]}}{$val}}{$f}++;}
   }
}

open (OUT,">/net/home/vseplyarsky/GCPR/mutational_model/script/results/extended_hash_simple_v3\_$CHR");
my @fa=sort {$a<=>$b} keys %f;
my $o='freq';
for my $i (0..$#fa){$o.=" $fa[$i]"}
print OUT "$o\n";
foreach my $cont (sort keys %outs){
    foreach my $val (sort {$a<=>$b} keys %{$outs{$cont}}){
	my $o="$cont $val ${$outs{$cont}}{$val}";
	for my $i (0..$#fa){ if (!(${${$out{$cont}}{$val}}{$fa[$i]})){${${$out{$cont}}{$val}}{$fa[$i]}=0;} $o.=" ${${$out{$cont}}{$val}}{$fa[$i]}"}
	print OUT "$o\n";
    }}




