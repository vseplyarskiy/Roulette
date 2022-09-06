To add the raw Roulette mutation rates please use add_raw_rates.py.

Example:
  python add_raw_rates.py --vcf_dir directory_for_Roulette_rates_vcf --input_filename input_filename_to_add_raw_rates --output_header output_header
  
The input filename must be a tsv file with the first column CHROM, second column POS, third column REF, and fourth column ALT. You can have other additional columns. The input filename must have a header column.

This code will output two files. First a file called output_header + ".tsv" file with an added mutation rate. The second file is called output_header + "_binned.tsv", which is a summary of the number of sites for each mutation rate bins.


To scale the Roulette rates, please use scale_mutation_rate.py.

Example:
  python scale_mutation_rate.py --vcf_dir directory_for_Roulette_rates_vcf --input_filename input_filename_to_add_raw_rates --output_header output_header
