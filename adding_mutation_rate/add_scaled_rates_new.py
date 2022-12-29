import pandas as pd
import gzip
import os
import csv
import optparse
import sys
import numpy as np
from scipy.optimize import minimize

from dask.distributed import Client
import dask.dataframe as dd

def scale_by_background(vcf_dir, input_filename, background, denovo):
    
    if background == 1:
        df = dd.read_csv(vcf_dir + "/all_hq_synonymous_variants.tsv", sep = "\t")
    else:
        df = dd.read_csv(vcf_dir + "/noncoding/*.tsv.gz", sep = "\t")

    df_observed = dd.read_csv(input_filename, sep = "\t")
    df_observed["polymorphic"] = 1

    df_merged = df[["CHROM", "POS", "REF", "ALT", "mu"]].merge(df_observed[["CHROM", "POS", "REF", "ALT", "polymorphic"]], how = "left", on = ["CHROM", "POS", "REF", "ALT"])
    
    if denovo == 0:
        df_group = pd.DataFrame(df_merged.groupby("mu").agg({"polymorphic": ['size', 'sum']}).compute())
        df_group.columns = df_group.columns.droplevel()

        df_group = df_group.rename({'size': "total", "sum": "polymorphic"}, axis = 1)
        df_group = df_group.reset_index()
        
        df_group["poisson_lambda"] = 1 - np.exp(-1*df_group["mu"])
        df_group["expected_polymorphic"] = df_group["poisson_lambda"] * df_group["total"]

        def poisson_function(x):    
            return abs(sum((1 - np.exp(-1 * x * df_group["mu"]))*df_group["total"]) - sum(df_group["polymorphic"]))

        x0 = np.array([1])
        res = minimize(poisson_function, x0, method='Nelder-Mead', options={'xatol': 1e-3, 'disp': True})

        # this is the proper scaling factor for the mutation rate
        k = res.x[0]
    else:
        mu_sum = df["mu"].sum().compute()
        polymorphic_sum = df_merged["polymorphic"].sum().compute()
        k = polymorphic_sum/mu_sum
        
    return k


def main(vcf_dir, input_filename, output_dir, quality_filter, background, denovo):
    # step 1: find the proper scaling
    client = Client()
    k = scale_by_background(vcf_dir, input_filename, background, denovo)
        
    # step 2: output new file with mutation rate estimates
    chrom_set = [str(x) for x in range(1, 23)]
    
    for chrom in chrom_set:    
        output_filename = chrom + "_scaled_rate_v5.2_TFBS_correction_all.tsv.gz"
        f_out = gzip.open(os.path.join(output_dir, output_filename), "wt", newline="")
        out_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")

        mut_fname_v4 = chrom + "_rate_v5.2_TFBS_correction_all.vcf.gz"
        mut_file_v4 = gzip.open(os.path.join(vcf_dir, mut_fname_v4), "rt")
        mut_reader_v4 = csv.reader(mut_file_v4, delimiter="\t")
        
        #define header for output file
        header = ["CHROM", "POS", "REF", "ALT", "FILTER", "mu_roulette_original", "mut_prob"]
        out_writer.writerow(header)
        
        #iterate over the mut rate
        for mut_current_v4 in mut_reader_v4:
            
            #skip over the comments in the vcf file
            if mut_current_v4[0][0] == "#":
                continue
                
            if "MR" in mut_current_v4[7]:            
                mu = float(mut_current_v4[7].split(";")[1][3:])
                
                if denovo == 0:
                    out_writer.writerow([mut_current_v4[0], mut_current_v4[1], mut_current_v4[3], mut_current_v4[4], 
                                        mut_current_v4[6], mu, 1 - np.exp(-1 * k * mu)])
                else:
                    out_writer.writerow([mut_current_v4[0], mut_current_v4[1], mut_current_v4[3], mut_current_v4[4], 
                                        mut_current_v4[6], mu, k*mu])
            else:
                out_writer.writerow([mut_current_v4[0], mut_current_v4[1], mut_current_v4[3], mut_current_v4[4], 
                                    mut_current_v4[6], None, None])                
        
        mut_file_v4.close()
        f_out.close()


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("--vcf_dir", type="str", dest="vcf_dir", default="", help="location of directory of mutation rate VCF file")
    parser.add_option("--input", type="str", dest="input_filename", help="specify input filename, this is a tsv file with each row being a site with mutation; the tsv file must have a CHROM and POS column for the genomic coordinate, in GRCh38, and also REF and ALT column. The input file should be sorted by CHROM and POS, in ascending order. The CHROM column should have strings such as 1, 2, and 3 instead of chr1, chr2, and chr3.")
    parser.add_option("--output_dir", type="str", dest="output_dir", default="", help="specify directory for output files")
    parser.add_option("--quality", type="int", default=1, dest="quality_filter", help="specify quality filter you want. 0 for no filter, 1 for filtering low-quality regions.")
    parser.add_option("--background_type", type="int", default=1, dest="background", help="1 for synonymous variants (whole exome) as background. 0 for intergenic region as background (whole genome).")
    parser.add_option("--denovo", type="int", default=0, dest="denovo", help="1 for denovo, 0 for population sequencing (0 is default)")
    
    opts, args = parser.parse_args()

    main(opts.vcf_dir, opts.input_filename, opts.output_dir, opts.quality_filter, opts.background, opts.denovo)
    
    