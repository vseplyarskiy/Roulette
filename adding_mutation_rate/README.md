<h1 align="center"> Adding Roulette Rates </h1>

With Roulette, we have a basepair-resolution mutation rate estimates for the human genome. To facilitate analyses using Roulette, we provide a python script to convert unscaled Roulette esimates into the probablility of observing a particular mutation in a site. Scaling rates is important because the probability of actually observing mutations can differ based on the type and composition of the sample. Samples of de novo mutations may differ depending on the age distribution of parents, and the probability of observing a particular mutation in a population sample strongly depends on the sample size.

First, download the raw mutation rate files from this [link](http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/). 

If you simply wish to convert the provided estimates to approximate per-generation rates here are some scaling factors you may use. For Roulette, you can multiply the raw rates by $1.015*10^-7$ to get the per generation mutation rate. For Carlson rates, you can multiply by $2.086 * 10^-9$. The gnomAD rates are already scaled to be approximately per-generation.

In population samples, we cannot use the raw mutation rates directly because we need to linearly scale the mutation rate to match the effective population size of the population sequencing data (more samples increases the likelihood of observing a mutation). Therefore, we need to scale the mutation rate so that the number of expected mutations equals the number of observed mutations for a set of neutral sites, which will refer to as background sites. After scaling the mutation rate we can use the formula $p = 1 - e^{-\mu}$, where $\mu$ is the scaled rate, to get the probability of the site being polymorphic.

We provide three ways to choose a set of background sites. The users can use synonymous sites, the whole genome, or manually provide a set of neutral sites. While the first two options are easier, we think that manually choosing neutral sites will be the most accurate way to scale the mutation rate properly.

Below are instructions on how to run the script.

## Instruction for packages

We use python3 for our script. Please install python packages pandas, scipy, and numpy.

## Instruction for input files

From a list of sites where the user wants a probability of observing a mutation, create a tsv file with with CHROM, POS, REF, and ALT as the first four columns. You may have additional columns after ALT. Please sort the file so that CHROM and POS are in ascending order. For the first row please provide a header. Note that we also do not have estimates for X and Y chromosomes.

Following is the first 10 lines of an example input file:
```sh
CHROM	POS	REF	ALT
1	924437	G	A	
1	924437	G	C	
1	924440	G	A	
1	924440	G	C	
1	924440	G	T	
1	924443	C	A	
```

## Synonymous variants as Background Sites

```sh
  python add_scaled_rates.py --vcf_dir Roulette_vcf_dir --input_filename input_filename --output_header output_header --syn synonymous_variants_filename
```
Roulette_vcf_dir is the directory where the vcf of Roulette rates is located. input_filename is the filename for the input file described above. synonymous_variants_filename is the name for the tsv file that contains a list of observed synonymous variants in the sample. The file should follow the same format as the input file described above.

## Using manual set of Background Sites

To use a user-chosen manual set of sites, as background sites, you must first run add_raw_rates.py.

```sh
  python add_raw_rates.py --vcf_dir Roulette_vcf_dir --input_filename background_sites_filename --output_header output_header
```
Here, note that background_sites_filename follows the same format as the input file described above. Also, background_sites_filename must include both observed variants and potential variants.

The script will output two files output_header.tsv file, which is the same tsv as background_sites_filename, but with the raw Roulette rates added. The second file is output_header_binned.tsv file, which is the summary of the number of sites for each mutation rate bin. output_header_binned.tsv will be used as an input for the add_scaled_rates.py.

```sh
  python add_scaled_rates.py --vcf_dir Roulette_vcf_dir --input_filename input_filename --output_header output_header --background_sites output_header_binned.tsv --polymorphic_count number_of_mutations_in_background_sites
```
Please count up the number of observed mutations in the background sites, and pass as an argument to --polymorphic_count.
