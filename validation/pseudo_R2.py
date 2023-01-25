import argparse
from argparse import RawTextHelpFormatter

import pandas as pd
import numpy as np
import scipy.optimize
import scipy.stats as stats

def exp_likelihood_model(mut_p, bin_p):
    bin_ll = mut_p * np.log(mut_p) + ((1 - mut_p) * np.log(1 - mut_p))
    return np.sum(bin_ll * bin_p)

def exp_likelihood_null(mut_p, bin_p):
    mean_p = np.sum(mut_p * bin_p)
    return mean_p * np.log(mean_p) + (1 - mean_p) * np.log(1 - mean_p)

def avg_likelihood_sample(sites, mutations, mut_p):
    return np.sum(mutations * np.log(mut_p) + (sites - mutations) * np.log(1 - mut_p)) / np.sum(sites)

def null_likelihood_sample(sites, mutations, mut_p, bin_p):
    mean_p = np.sum(mut_p * bin_p)
    return np.sum(mutations * np.log(mean_p) + (sites - mutations) * np.log(1 - mean_p)) / np.sum(sites)

def calculate_R2(ll, ll_null):
    return 1 - np.exp(-2* (ll - ll_null))

def sample_sites(sites, mutations):
    sites_sample = stats.multinomial.rvs(np.sum(sites), sites / np.sum(sites))
    mutations_sample = stats.binom.rvs(sites_sample, mutations / sites)
    return sites_sample, mutations_sample

def calc_pR2(sites, mutations, rates):
    bin_p = sites / np.sum(sites)

    # Use a Poisson link to scale rates to the probability of observing a mutation
    # This is relevant if analyzing population data where recurrent mutation is likely
    def poisson_link(x):
        return -np.sum(stats.binom.logpmf(k=mutations,
                                          n=sites,
                                          p=(1 - np.exp(-x * rates))))
    # Find a non-awful starting place for optimization
    # Replace with ML estimate at some point
    x0 = np.array([-np.log(0.99) / np.max(rates)])
    scaling_fit = scipy.optimize.minimize(poisson_link, x0, method='Nelder-Mead',
                                          options={'xatol': 1e-3, 'disp': False})
    scaling_factor = scaling_fit.x[0]

    mut_p = 1 - np.exp(-scaling_factor * rates)

    # Calculate pseudo-R2 for the observed data
    e_ll = exp_likelihood_model(mut_p, bin_p)
    e_ll_null = exp_likelihood_null(mut_p, bin_p)
    a_ll = avg_likelihood_sample(sites, mutations, mut_p)
    a_ll_null = null_likelihood_sample(sites, mutations, mut_p, bin_p)

    max_R2 = calculate_R2(e_ll, e_ll_null)
    R2 = calculate_R2(a_ll, a_ll_null)
    pR2 = R2 / max_R2
    return max_R2, R2, pR2

def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     epilog="Calculates an adjusted pseudo-R^2 for estimated mutation rates.\nReturns the estimate and 95% CI calculated by bootstrap resampling of sites:\nestimated pseudo-R^2\t0.025Q\t0.975Q")
    parser.add_argument("--mutation_table",
                        help="tsv file with columns for \n- mutation rate (\"rate\")\n- total number of sites with that rate (\"sites\")\n- number of observed mutations at sites with that rate (\"mutations\")")
    parser.add_argument("--n_bootstrap",
                        help="Number of bootstrap resamples", type=int)
    parser.add_argument("--out_bootstrap",
                        help="Location to save bootstrap pseudo-R^2 values")
    args = parser.parse_args()

    mut_fname = args.mutation_table
    mut_table = pd.read_csv(mut_fname, sep="\t")

    # Group together rows with exactly the same mutation rate
    mut_table = mut_table[["rate", "sites", "mutations"]].groupby("rate").sum().reset_index()

    rates = mut_table["rate"].to_numpy(dtype=float)
    sites = mut_table["sites"].to_numpy(dtype=int)
    mutations = mut_table["mutations"].to_numpy(dtype=int)

    max_R2, R2, pR2 = calc_pR2(sites, mutations, rates)

    max_R2_s = []
    R2_s = []
    pR2_s = []
    # Calculate pseudo-R2 for a set of boostrap samples
    for _ in range(args.n_bootstrap):
        sites_sample, mutations_sample = sample_sites(sites, mutations)

        max_R2_sample, R2_sample, pR2_sample = calc_pR2(sites_sample, mutations_sample, rates)
        max_R2_s.append(max_R2_sample)
        R2_s.append(R2_sample)
        pR2_s.append(pR2_sample)

    print(pR2, np.quantile(pR2_s, 0.025), np.quantile(pR2_s, 0.975), sep="\t")

    if args.out_bootstrap:
        pd.DataFrame({"pR2":pR2_s}).to_csv(args.out_bootstrap, index=False)

if __name__ == "__main__":
    main()
