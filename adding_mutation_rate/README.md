<h1 align="center"> Adding Roulette Rates </h1>

With Roulette, we have a basepair-resolution mutation rate estimates for the human genome. We wrote some code in python to make it easier for people to use the mutation rate estimates for their own analyses.

First, download the raw mutation rate files from this [link](http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/).

Unfortunately, we cannot use the raw mutation rates directly because we need the linearly scale the mutation rate for each population sequencing data, since the mutation rate is linearly dependent on the effective population size of the sample. Therefore, we need to scale the mutation rate so that the number of expected mutations equals the number of observed mutations for a set of neutral (background) sites. We provide three ways to choose a set of background sites. The users can use synonymous sites, the whole genome, or manually provide a set of neutral sites. While the first two options are easier, we think that manually choosing a neutral sites will be the most accurate way to scale the mutation rate properly.

Below is instructions on how to get a probability of a mutation mutation under neutrality.

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

## Whole Genome as Background Sites

This script is currently in progress.
