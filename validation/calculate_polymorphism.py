import pandas as pd
import numpy as np


def simulate(df):
    
    df_sample = df[["rate", "polymorphic"]].sample(n = int(df["polymorphic"].sum()), weights = 'polymorphic', replace = True)

    df_sample = pd.DataFrame(df_sample.groupby("rate").size()).reset_index()
    df_sample.rename({0: "polymorphic_sample"}, axis = 1, inplace = True)

    df_merged = df.merge(df_sample, on = "rate", how = "left")
    df_merged["polymorphic_sample"] = df_merged["polymorphic_sample"].fillna(0)
    
    return df_merged

def bernoulli_likelihood(mu):
    return (mu * np.log(mu)) + ((1 - mu) * np.log(1 - mu))

def likelihood_max(df):
    df["bernoulli_likelihood"] = df["p"].apply(lambda x: bernoulli_likelihood(x))
    l_beta_max = sum(df["bernoulli_likelihood"] * df["prob_bin"])
    return l_beta_max

def likelihood_null(df, column = "polymorphic_sample"):
    
    mean_mu = sum(df["p"] * df["sites"])/sum(df["sites"])

    log_likelihood = sum(df[column] * np.log(mean_mu)) + sum(df["monomorphic"] * np.log(1 - mean_mu))

    l_null = log_likelihood/df["sites"].sum()
    
    return l_null

def calculate_r2(ll, l_null):
    r2 = 1 - np.exp(-2* (ll - l_null))
    return r2

def calculate_loglikelihood(df, column = "polymorphic_sample"):
    log_likelihood = sum(df[column] * np.log(df["p"])) + sum(df["monomorphic"] * np.log(1 - df["p"]))
    
    ll = log_likelihood/df["sites"].sum()
        
    return ll