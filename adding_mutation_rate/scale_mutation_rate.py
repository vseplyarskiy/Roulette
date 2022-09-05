import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import dask.dataframe as dd

def main():
    args = sys.argv[1:]
    
    input_filename = args[0]
    output_filename = args[1]

    # each row is a mutation rate bin
    # the "polymorphic" column indicating whether a site is polymorphic (1) or not (0)
    # "mu" column is the raw mutation rate
    
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster()

    with Client() as client:

        df = dd.read_csv(input_filename, sep = "\t")

        #make new dataframe, where you group by mutation rate
        df_mu = pd.DataFrame(df.groupby("mu").size().compute(), columns = ["sites"])
        df_mu["polymorphic"] = df.groupby("mu")["polymorphic"].sum()
        df_mu = df_mu.reset_index()

        df_mu["poisson_lambda"] = 1 - np.exp(-1*df_mu["mu"])
        df_mu["expected_polymorphic"] = df_mu["poisson_lambda"] * df_mu["sites"]

        def poisson_function(x):    
            return abs(sum((1 - np.exp(-1 * x * df_mu["mu"]))*df_mu["sites"]) - sum(df_mu["polymorphic"]))

        #find best fit k that scales mutation rate to the set of sites
        x0 = np.array([1])
        res = minimize(poisson_function, x0, method='Nelder-Mead', options={'xatol': 1e-3, 'disp': True})

        # this is the proper scaling factor for the mutation rate
        k = res.x[0]

        df["scaled_mu"] = k * df["mu"] 
        df["polymorphic_prob"] = 1 - np.exp(-1 * k * df["mu"])

        df.to_csv(output_filename, sep = "\t", index = None, single_file = True)
    
if __name__ == "__main__":
    main()