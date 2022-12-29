#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

# You can change hostname to any command you would like to run
python add_scaled_rates_new.py --vcf_dir /home/djl34/scratch/genetics.bwh.harvard.edu/downloads/Vova/Roulette --input /home/djl34/lab_pd/roulette/code/pseudo_r2_intergenic/noncoding_polymorphic_chr*.tsv --output_dir /home/djl34/scratch/trial_runs/intergenic --background_type 0