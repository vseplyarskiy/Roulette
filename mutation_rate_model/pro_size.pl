use strict;

my @fi=<../mut_model/for_LASSO/context/*pro*met*>;
for my $i (0..$#fi){
    my $h;
    `wc -l $fi[$i]>>ll`;
}
