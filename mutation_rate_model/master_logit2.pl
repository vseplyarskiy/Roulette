use strict;

my $run=$ARGV[0];


my $name='TTTTT_A_inter_met';
if ($run==2){$name='TTTTT_C_inter_met'}


#    print "$name\n";
    `perl /net/home/vseplyarsky/GNOMAD_Tufts/script/filter_files1.pl $name`;
    `Rscript /net/home/vseplyarsky/GNOMAD_Tufts/script/rscript.R $name\_cl`;
    `rm /net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/$name\_cl`;


