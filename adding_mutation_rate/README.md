<h1 align="center"> Adding Roulette Rates </h1>

With Roulette, we have a basepair-resolution mutation rate estimates for the human genome. We wrote some code in python to make it easier for people to use the mutation rate estimates for their own analyses.

First, download the raw mutation rate files from this [link](http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/).

Unfortunately, we cannot use the raw mutation rates directly because we need the linearly scale the mutation rate for each population sequencing data, since the mutation rate is linearly dependent on the effective population size of the sample. Therefore, we need to scale the mutation rate so that the number of expected mutations equals the number of observed mutations for a set of neutral (background) sites. We provide three ways to choose a set of background sites. The users can use synonymous sites, the whole genome, or manually provide a set of neutral sites. While the first two options are easier, we think that manually choosing a neutral sites will be the most accurate way to scale the mutation rate properly.

Below is instructions on how to get a probability of a mutation mutation under neutrality.

## Instruction for input files

From a list of sites where the user wants a probability of observing a mutation, create a tsv file with with CHROM, POS, REF, and ALT as the first four columns. You may have additional columns after ALT. Please sort the file so that CHROM and POS are in ascending order. For the first row please provide a header. Note that we also do not have estimates for X and Y chromosomes.

Following is the first 10 lines of an example input file:
```sh
CHROM	POS	REF	ALT	AC	polymorphic	mu	mu_quality
1	924437	G	A	0.0	0	0.4	TFBS
1	924437	G	C	0.0	0	0.18600000000000003	TFBS
1	924437	G	T	0.0	0	0.139	TFBS
1	924440	G	A	0.0	0	0.342	TFBS
1	924440	G	C	0.0	0	0.18600000000000003	TFBS
1	924440	G	T	0.0	0	0.073	TFBS
1	924443	C	A	0.0	0	0.03	TFBS
```

# Synonymous variants as Background Sites

```sh
python add_scaled_rates.py --vcf_dir Roulette_vcf_dir --input_filename input_filename --output_header output_header --syn synonymous_variants_filename

```

# Manual Background Sites

Example:
  python add_raw_rates.py --vcf_dir directory_for_Roulette_rates_vcf --input_filename input_filename_to_add_raw_rates --output_header output_header
  
The input filename must be a tsv file with the first column CHROM, second column POS, third column REF, and fourth column ALT. You can have other additional columns. The input filename must have a header column.

This code will output two files. First a file called output_header + ".tsv" file with an added mutation rate. The second file is called output_header + "_binned.tsv", which is a summary of the number of sites for each mutation rate bins.


# To scale the Roulette rates, please use add_scaled_rates.py.

Example:
  python add_scaled_rates.py --vcf_dir directory_for_Roulette_rates_vcf --input_filename input_filename_to_add_scaled_rates --output_header output_header --background_sites binned_background_sites_filename --polymorphic_count number_of_mutations_in_background_sites
  
The input filename must be a tsv file with the first column CHROM, second column POS, third column REF, and fourth column ALT. You can have other additional columns. The input filename must have a header column.

In addition to input filename and the directory for the Roulette rates, you must provide binned_background_sites_filename, which is the "_binned.tsv" output from add_raw_rates.py. In order to scale the mutation rate properly, we need information on the neutral set of sites. Please select a set of sites that you consider neutral and run this on add_raw_rates.py.

For --polymorphic_count please provide the number of sites with mutations from the set of background sites.

This code will output two files. First a file called output_header + ".tsv" file with two columns: scaled mutation rate and a probability of observing a mutation for each site.
