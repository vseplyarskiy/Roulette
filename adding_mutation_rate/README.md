# To add the raw Roulette mutation rates please use add_raw_rates.py.

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
