<h1 align="center"> Adding Roulette Rates </h1>

With Roulette, we have a basepair-resolution mutation rate estimates for the human genome. To facilitate analyses using Roulette, we provide a python script to convert the unscaled Roulette estimates into the probablility of observing a particular mutation in a site. Scaling rates is important because the probability of actually observing mutations can differ based on the type and composition of the sample. Samples of de novo mutations may differ depending on the age distribution of parents, and the probability of observing a particular mutation in a population sample strongly depends on the sample size.

First, download the raw mutation rate files from this [link](http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/). Then please unzip the all_hq_synonymous_variants.tsv.gz file to all_hq_synonymous_variants.tsv using the command ``` gunzip all_hq_synonymous_variants.tsv.gz``` 

If you simply wish to convert the provided estimates to approximate per-generation rates here are some scaling factors you may use. For Roulette, you can multiply the raw rates by $1.015*10^-7$ to get the per generation mutation rate. For Carlson rates, you can multiply by $2.086 * 10^-9$. The gnomAD rates are already scaled to be approximately per-generation.

In population samples, we cannot use the raw mutation rates directly because we need to linearly scale the mutation rate to match the effective population size of the population sequencing data (more samples increases the likelihood of observing a mutation). Therefore, we need to scale the mutation rate so that the number of expected mutations equals the number of observed mutations for a set of neutral sites, which will refer to as background sites. After scaling the mutation rate we can use the formula $p = 1 - e^{-\mu}$, where $\mu$ is the scaled rate, to get the probability of the site being polymorphic.

We provide two ways to choose a set of background sites. The users can use synonymous sites or manually provide a set of neutral sites. While the first options is easier, we think that manually choosing neutral sites will be the most accurate way to scale the mutation rate properly.

Below are instructions on how to run the script.

## Instruction for packages

We use python3 for our script. Please install python packages pandas, scipy, and numpy.

## Instruction for input files

Form a list of observed mutations, create a tsv file with with CHROM, POS, REF, and ALT as the first four columns. You may have additional columns after ALT. Please sort the file so that CHROM and POS are in ascending order. For the first row please provide a header. Note that we also do not have estimates for X and Y chromosomes.

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
Please do not zip the tsv file.

The script will output a new tsv file with the columns "mu_roulette_original" for the original roulette rate and "mut_prob" for the probability of observing a mutation (scaled to the user's dataset). 

We give two options for scaling the mutation rate. First, is the option to choose whether the dataset is a population sequencing dataset or denovo sequencing dataset. Second, the user may choose to scale to match the number of observed synonymous variants (for whole exome sequencing) or to number of observed non-coding variants (for whole genome sequencing).


## Synonymous variants as Background Sites

```sh
  python add_scaled_rates.py --vcf_dir Roulette_vcf_dir --input input_filename --output_dir output_directory
```
Roulette_vcf_dir is the directory where the vcf of Roulette rates is located. Input_filename is the filename for the input file described above. If the dataset is denovo sequencing, please include the option ``` --denovo 1```. For population sequencing dataset, you do not need to include any extra options.

## Nonc-coding variants as Background Sites

```sh
  python add_scaled_rates.py --vcf_dir Roulette_vcf_dir --input input_filename --output_dir output_header --background_type 0
```

Similarly to synonymous variants, for denovo sequencing, please include the option ```--denovo 1```.

