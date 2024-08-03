<h1 align="center"> Adding Roulette Rates </h1>

With Roulette, we have a basepair-resolution mutation rate estimates for the human genome. To facilitate analyses using Roulette, we provide a python script to convert the unscaled Roulette estimates into the probablility of observing a particular mutation in a site. Scaling rates is important because the probability of actually observing mutations can differ based on the type and composition of the sample. Samples of de novo mutations may differ depending on the age distribution of parents, and the probability of observing a particular mutation in a population sample strongly depends on the sample size.

First, download the raw mutation rate files from this [link](http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/). 

To directly annotate a set of variants using the provided set of VCF files, we recommend using `bcftools` and the `annotate` function. We have also provided a simple bash script to efficiently annotate relatively small sets of variants: `roulette_annotate.sh`. You can run this as `./roulette_annotate.sh roulette_input.tsv /wherever/roulette/lives`. Using an input file in the following format:

CHROM	POS	REF	ALT<br>
1	905251	A	C<br>
1	905251	A	G<br>
1	905251	A	T<br>
1	905692	T	C<br>
1	954695	C	T<br>

If you simply wish to convert the provided estimates to approximate per-generation rates you can multiply by different scaling factors. For Roulette, multiply the raw rates by $1.015*10^{-7}$ to get the per generation mutation rate. For Carlson estimates, use $2.086 * 10^{-9}$. The gnomAD rate estimates are already scaled to be approximately per-generation. Please note that this is a mutation rate per **diploid** genome.

In population samples, we cannot use the raw mutation rates directly because we need to linearly scale the mutation rate to match the size and overall genetic diversity of the sample (more individuals increases the likelihood of observing a mutation). Therefore, we need to use a set of putative neutral potential mutations to scale the mutation rate so that the number of expected mutations equals the number of observed mutations. We will refer to the chosen set of neutral mutations as background sites even though each corresponds to one of three possible single nucleotide mutations at a site. After scaling the mutation rate we  use the formula $p = 1 - e^{-\mu}$, where $\mu$ is the scaled rate, to get the probability a given site is polymorphic.

We provide two ways to choose a set of background sites: synonymous sites or a user-provided set of neutral sites. While the first options is easier, we think that manually choosing neutral sites will be the most accurate way to scale the mutation rate properly. Download the all_hq_synonymous_variants.tsv.gz containing synonymous variants from gnomAD v2.1.1 and decompress using the command ``` gunzip all_hq_synonymous_variants.tsv.gz```.

Below are instructions on how to run the script.

## Instruction for packages

We use python3 for our script. The environment.yml file is provided in this directory.

## Instruction for input files

Form a list of observed mutations, create a tsv file with with CHROM, POS, REF, and ALT as the first four columns. You may have additional columns after ALT. Sort the file so that CHROM and POS are in ascending order. For the first row please provide a header. Note that we also do not have estimates for X and Y chromosomes.

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

