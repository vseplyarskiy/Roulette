import pandas as pd
import gzip
import os
import csv
import optparse

def main(vcf_dir, input_filename, output_header, input_zip, output_zip, chromosome, quality_filter):
    
    mutation_rate_list = []
        
    if chromosome == 0:
        output_file = output_header
    else:
        output_file = output_header +"_chr" + str(chromosome)
        
    output_filename = output_file + ".tsv"
    
    if input_zip == True:
        f_in = gzip.open(input_filename, "rt")
    else:
        f_in = open(input_filename, "rt")
        
    if output_zip == True:
        f_out = gzip.open(output_filename, "wt", newline="")
    else:
        f_out = open(output_filename, "wt", newline="")

    mut_fname_v4 = "1_rate_v5.2_TFBS_correction_all.vcf.gz"
    mut_file_v4 = gzip.open(os.path.join(vcf_dir, mut_fname_v4), "rt")
    mut_reader_v4 = csv.reader(mut_file_v4, delimiter="\t")

    mut_current_v4 = next(mut_reader_v4)
    while mut_current_v4[0][0] == "#":
        mut_current_v4 = next(mut_reader_v4)

    pos_current_v4 = int(mut_current_v4[1])

    in_reader = csv.reader(f_in, delimiter="\t")
    out_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")

    header = next(in_reader)
    out_writer.writerow(header + ["mu", "mu_quality"])

    quality = None

    chrom_col = 0
    pos_col = 1
    ref_col = 2
    alt_col = 3

    for row in in_reader:

        chrom = row[chrom_col]
        
        if chromosome != 0:
            if int(chrom) != chromosome:
                continue
        try:
            pos = int(row[pos_col])
        except ValueError:
            continue

        ref = row[ref_col]
        alt = row[alt_col]
        # print(chrom, pos, ref, alt)

        ## Update mutation rate reader if not on the right chromosome
        if chrom != mut_fname_v4.split("_")[0]:
            mut_file_v4.close()
            mut_fname_v4 = "{}_rate_v5.2_TFBS_correction_all.vcf.gz".format(chrom)

            print(mut_fname_v4)

            mut_file_v4 = gzip.open(os.path.join(vcf_dir, mut_fname_v4), "rt")
            mut_reader_v4 = csv.reader(mut_file_v4, delimiter="\t")
            mut_current_v4 = next(mut_reader_v4)

            while mut_current_v4[0][0] == "#":
                mut_current_v4 = next(mut_reader_v4)
                            
            pos_current_v4 = int(mut_current_v4[1])

        # Locate current mutation in v4
        if pos_current_v4 != pos:
            pos_current_v4 = int(mut_current_v4[1])
            while pos_current_v4 < pos:
                mut_current_v4 = next(mut_reader_v4)
                pos_current_v4 = int(mut_current_v4[1])
            ref_current_v4 = mut_current_v4[3]
            pos_ind_v4 = pos_current_v4

            if ref_current_v4 != ref and (pos_ind_v4 == pos):
                print("reference doesn't match, {}:{}".format(chrom, pos))
                    
            alt_dict_v4 = {}

            while (pos_ind_v4 == pos) and (ref_current_v4 == ref):
                quality = mut_current_v4[6]
                alt_dict_v4[mut_current_v4[4]] = float(mut_current_v4[7].split(";")[1][3:])

                mut_current_v4 = next(mut_reader_v4)
                pos_ind_v4 = int(mut_current_v4[1])
                
        if quality_filter == 1:
            if quality == "TFBS" or quality == "high":
                try:
                    out_writer.writerow(row + [alt_dict_v4[alt], quality])
                    mutation_rate_list.append(alt_dict_v4[alt])
                except KeyError:
                    out_writer.writerow(row + ["NA", "NA"])            
        else:        
            try:
                out_writer.writerow(row + [alt_dict_v4[alt], quality])
                mutation_rate_list.append(alt_dict_v4[alt])
            except KeyError:
                out_writer.writerow(row + ["NA", "NA"])

    f_in.close()
    f_out.close()
    
    df_mu_list = pd.DataFrame(mutation_rate_list, columns = ["mu"])
    df_mu_list = pd.DataFrame(df_mu_list.groupby("mu").size())
    df_mu_list = df_mu_list.reset_index()
    df_mu_list.to_csv(output_file + "_binned.tsv", sep = "\t", index = None)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("--vcf_dir", type="str", dest="vcf_dir", default="", help="location of directory of mutation rate VCF file")
    parser.add_option("--input_filename", type="str", dest="input_filename", help="specify input filename, this is a tsv file with each row being a site with mutation; the tsv file must have a CHROM and POS column for the genomic coordinate, in GRCh38, and also REF and ALT column. The input file should be sorted by CHROM and POS, in ascending order. The CHROM column should have strings such as 1, 2, and 3 instead of chr1, chr2, and chr3.")
    parser.add_option("--output_header", type="str", dest="output_header", default="roulette_raw", help="specify header for output filenames")
    parser.add_option("--input_zip", action="store_true", default=False, help="specify whether input file is zipped")
    parser.add_option("--output_zip", action="store_true", default=False, help="specify whether output file is zipped")
    parser.add_option("--chr", type="int", default=0, dest="chromosome", help="specify if you want to get output for only one chromosome")
    parser.add_option("--quality", type="int", default=1, dest="quality_filter", help="specify quality filter you want. 0 for no filter, 1 for filtering low-quality regions.")
    opts, args = parser.parse_args()

    main(opts.vcf_dir, opts.input_filename, opts.output_header, opts.input_zip, opts.output_zip, opts.chromosome, opts.quality_filter)
    
    