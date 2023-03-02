import math
import numpy as np
import scipy.special as sp
import scipy.stats as stats


from numba import jit

# Algorithm from Jerome Kelleher
# https://jeromekelleher.net/generating-integer-partitions.html
def accel_asc(n):
    """
    Fast function to generate all partitions of an integer n
    """
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def mult_part(partition):
    """
    Convert a partition from a list of integers to multiplicity representation

    Parameters
    ----------
    partition : list
        [1, 1, 2, 2, 3]

    Returns
    ----------
    list
        [(1, 2), (2, 2), (3, 1)]
    """
    return [(ii, partition.count(ii)) for ii in set(partition)]


def mult_index(partition, max_count, max_recur):
    """
    Convert a partition from a list of integers to a boolean vector

    Parameters
    ----------
    partition : list
        [1, 1, 2, 2, 3]
    max_count : integer
        The maximum integer considered in the analysis,
        determines overall length of the vector
    max_recur : integer
        The maximum partition size to be considered

    Returns
    ----------
    Numpy array (boolean)
        Start index with 0:(max_count-1) for each unique integer in the partition.
        Add to each index max_recur times the nubmer of that integer occurs in the partition.
    """
    result = np.zeros(max_count*max_recur, dtype=np.bool)
    for ii in set(partition):
        result[(partition.count(ii)-1)*max_count + ii - 1] = 1
    return result

def bool_index_to_counts(part, val_vec=None, recur_vec=None):
    vals = part * val_vec
    recur = part * recur_vec
    return vals[vals!=0], recur[recur!=0]

def bool_index_counts_slow(part, max_count, max_recur, val_vec=None, recur_vec=None):
    val_vec = np.tile(np.arange(1, max_count+1), max_recur) if not val_vec else val_vec
    recur_vec = np.repeat(np.arange(1, max_recur+1), max_count) if not recur_vec else recur_vec
    return bool_index_to_counts(part, val_vec, recur_vec)

def max_secondary(part, val_vec, recur_vec):
    if len(part.shape) == 1:
        vals, _ = bool_index_to_counts(part, val_vec, recur_vec)
        return np.sort(vals)[...,-2]
    else:
        result = np.zeros(part.shape[0])
        for ii in range(part.shape[0]):
            vals, _ = bool_index_to_counts(part[ii,:], val_vec, recur_vec)
            result[ii] = np.sort(vals)[-2] if len(vals)>1 else 0
        return result
    
def remove_largest_part_pm(part_matrix, max_count, max_recur):
    """
    Remove the largest integer from each partition in the matrix
    
    Parameters
    ----------
    part_matrix: numpy array (boolean)
        [num. partitions, num. possible integer recurrencs (max_recur*max_count)]
    max_count : integer
        The maximum integer considered in the analysis,
        determines overall length of the vector
    max_recur : integer
        The maximum partition size to be considered
        
    Returns: 
        Numpy array (boolean)
        Integer combinations in boolean form
        NOTE: these are not partitions of a single integer
    """
    val_vec = np.tile(np.arange(1, max_count+1), max_recur)
    max_inds = np.argmax(part_matrix * val_vec, axis=-1)
    result = part_matrix.copy()
    result[np.arange(part_matrix.shape[0]), max_inds] = False
    down_inds = max_inds - max_count
    result[np.arange(part_matrix.shape[0])[down_inds>=0], down_inds[down_inds>=0]] = True
    return result
    
def partition_matrices(max_count, max_recur, fixed_recur=None):
    """
    Generate boolean index partition matrices for each integer up to max_count

    Parameters
    ----------
    max_count : integer
        The maximum integer considered in the analysis,
        determines overall length of the vector
    max_recur : integer
        The maximum partition size to be considered
    fixed_recur : integer
        Use only this recurrence count, still write vector as max_recur

    Returns
    ----------
    list
        A list of numpy arrays. Each array contains all partitions
        of the integer indexing it.
        [num. partitions, num. possible integer recurrencs (max_recur*max_count)]
    """
    if fixed_recur is not None:
        return [np.array([mult_index(part, max_count, max_recur)
                          for part in accel_asc(ii) if len(part) == fixed_recur])
                for ii in range(1, max_count+1)]
    else:
        return [np.array([mult_index(part, max_count, max_recur)
                          for part in accel_asc(ii) if len(part) <= max_recur])
                for ii in range(1, max_count+1)]

def sfs_exp(max_count, Theta, Upsilon):
    """
    The expected SFS of rare neutral alleles under exponential growth
    (NO RECURRENCE)

    Parameters
    ----------
    max_count : integer
        The maximum count to consider in the SFS
    Theta : float
        2N0 * mu, where N0 is the current population size and mu is the mutation rate per generation
    Upsilon : float
        2N0 * gamma / m, where gamma is the exponential growth rate and m is the sample size

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (without recurrence)
    """
    jj = np.arange(1, max_count+1, dtype=np.float64)
    sfs = (1 + jj - jj*Upsilon*sp.hyp2f1(1, 1+jj, 2+jj, 1-Upsilon))/(jj**2 + jj)
    return Theta*sfs

def sfs_exp_v(ac, Theta, Upsilon):
    """
    The expected SFS of rare neutral alleles under exponential growth
    (NO RECURRENCE)

    Parameters
    ----------
    ac      : Numpy array (integer)
        Allele counts to calculate the SFS for
    Theta   : float
        2N0 * mu, where N0 is the current population size and mu is the mutation rate per generation
    Upsilon : float
        2N0 * gamma / m, where gamma is the exponential growth rate and m is the sample size

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (without recurrence)
    """
    sfs = (1 + ac - ac*Upsilon*sp.hyp2f1(1, 1+ac, 2+ac, 1-Upsilon))/(ac**2 + ac)
    return Theta*sfs

def sfs_const(max_count, Theta):
    """
    The expected SFS of rare neutral alleles in a constant-size population
    (NO RECURRENCE)

    Parameters
    ----------
    max_count : integer
        The maximum count to consider in the SFS
    Theta : float
        2N * mu, where N is the population size and mu is the mutation rate per generation

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (without recurrence)
    """
    jj = np.arange(1, max_count+1, dtype=np.float64)
    return Theta/jj

def recur_expand(sfs, max_recur):
    """
    Compute (not yet normalized) Poisson probabilities for each entry in the expected SFS
    up to a maximum value

    Parameters
    ----------
    sfs : numpy array (float)
        The expected SFS
    max_recur : integer
        The maximum number of independent mutations we want to consider occuring at a single site
    """
    return np.concatenate([sfs**kk/sp.factorial(kk) for kk in range(1, max_recur+1)])

def recur_sfs_exp(partition_mats, max_count, max_recur, Theta, Upsilon, count_recurrence=False):
    """
    The expected SFS of rare neutral alleles under exponential growth
    (WITH RECURRENCE)

    Parameters
    ----------
    partition_mats : list of numpy arrays (boolean)
        Partition matrices give the combinations of SFS entries one must consider to compute
        an observed allele count with recurrence
    max_count : integer
        The maximum count to consider in the SFS
    max_recur : integer
        The maximum number of independent mutations we want to consider occuring at a single site
    Theta : float
        2N0 * mu, where N0 is current population size and mu is the mutation rate per generation
    Upsilon : float
        2N0 * gamma / m, where gamma is the exponential growth rate and m is the sample size

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (with recurrence)
    """
    result = np.zeros((max_count+1, max_recur)) if count_recurrence else np.ones(max_count+1)
    if count_recurrence:
        recur_vec = np.repeat(np.arange(1, max_recur+1), max_count)
        recur_vals = [np.matmul(pm, recur_vec) for pm in partition_mats]
        result[0,0] = 1

    sfs = sfs_exp(max_count, Theta, Upsilon)
    recur_probs = recur_expand(sfs, max_recur)
    for ii in range(1, max_count+1):
        if len(partition_mats[ii-1]) == 0:
            print("no partitions given")
            result[ii] = 0
        else:
            if count_recurrence:
                for jj in range(max_recur):
                    result[ii, jj] = np.sum(np.exp(np.matmul(partition_mats[ii-1][recur_vals[ii-1]==(jj+1),:], 
                                                             np.log(recur_probs))))
            else:
                result[ii] = np.sum(np.exp(np.matmul(partition_mats[ii-1], np.log(recur_probs))))
    return result

def recur_sfs_exp_approx(ac, down_mat, max_count, max_recur, Theta, Upsilon):
    n_down_partitions = down_mat.shape[0]
    val_vec = np.tile(np.arange(1, max_count+1), max_recur)
    recur_vec = np.repeat(np.arange(1, max_recur+1), max_count)
    down_sizes = np.sum((down_mat*val_vec) * (down_mat*recur_vec), axis=-1)
    # (len(ac), n_down_partitions)
    big_parts = ac[:,np.newaxis] - down_sizes
    sfs_down = sfs_exp_v(np.arange(1, max_count+1), Theta, Upsilon)
    recur_probs_down = recur_expand(sfs_down, max_recur)
    # (n_down_partitions,)
    part_probs_down = np.matmul(down_mat, np.log(recur_probs_down))
    # Calculate the expected nonrecurrent sfs for the "big" parts 
    sfs = sfs_exp_v(big_parts, Theta, Upsilon)
    part_probs_full = np.exp(np.log(sfs) + part_probs_down)
    # Add additional 2-recurrence partitions not captured by down_mat
    extra_two = [np.arange(np.ceil(max_count/2), np.floor(ac_i/2)) for ac_i in ac]
    extra_probs = np.array([np.sum(sfs_exp_v(extra_two[ii], Theta, Upsilon) * 
                                   sfs_exp_v(ac[ii] - extra_two[ii], Theta, Upsilon)) 
                            for ii in range(len(ac))])
    return np.sum(part_probs_full, axis=1) + extra_probs

def recur_sfs_Ti_approx(ac,  down_mat, max_count, max_recur, ac_set, Ti, Theta=1):
    interp_sfs = lambda x: np.exp(np.interp(np.log(x+1), np.log(ac_set+1), np.log(Ti)))
    n_down_partitions = down_mat.shape[0]
    val_vec = np.tile(np.arange(1, max_count+1), max_recur)
    recur_vec = np.repeat(np.arange(1, max_recur+1), max_count)
    down_sizes = np.sum((down_mat*val_vec) * (down_mat*recur_vec), axis=-1)
    # (len(ac), n_down_partitions)
    big_parts = ac[:,np.newaxis] - down_sizes
    sfs_down = interp_sfs(np.arange(1, max_count+1)) * Theta
    recur_probs_down = recur_expand(sfs_down, max_recur)
    # (n_down_partitions,)
    part_probs_down = np.matmul(down_mat, np.log(recur_probs_down))
    # Calculate the expected nonrecurrent sfs for the "big" parts 
    sfs = interp_sfs(big_parts) * Theta
    part_probs_full = np.exp(np.log(sfs) + part_probs_down)
    # Add additional 2-recurrence partitions not captured by down_mat
    extra_two = [np.arange(np.ceil(max_count/2), np.floor(ac_i/2)) for ac_i in ac]
    extra_probs = np.array([np.sum(interp_sfs(extra_two[ii]) * interp_sfs(ac[ii]))*Theta**2
                            for ii in range(len(ac))])
    return np.sum(part_probs_full, axis=1) + extra_probs
    
def recur_sfs_Ti_combo(ac, partition_mats, down_mat, max_count, max_recur, ac_set, Ti, exact_cutoff, Theta=1):
    ac_low = ac[ac<=exact_cutoff]
    pred_sfs_low = rr.recur_sfs_Ti(part_mats_compare, max_count_compare, 
                                    max_recur_compare, Ti=Ti, Theta=theta_high)

def recur_sfs_const(partition_mats, max_count, max_recur, Theta):
    """
    The expected SFS of rare neutral alleles in a constant-size population
    (WITH RECURRENCE)

    Parameters
    ----------
    partition_mats : list of numpy arrays (boolean)
        Partition matrices give the combinations of SFS entries one must consider to compute
        an observed allele count with recurrence
    max_count : integer
        The maximum count to consider in the SFS
    Theta : float
        2N * mu, where N is the population size and mu is the mutation rate per generation

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (with recurrence)
    """
    result = np.ones(max_count+1)
    sfs = sfs_const(max_count, Theta)
    recur_probs = recur_expand(sfs, max_recur)
    for ii in range(1, max_count+1):
        if len(partition_mats[ii-1]) == 0:
            result[ii] = 0
        else:
            result[ii] = np.sum(np.exp(np.matmul(partition_mats[ii-1], np.log(recur_probs))))
    return result

def recur_sfs_Ti(partition_mats, max_count, max_recur, Ti, Theta, count_recurrence=False):
    """
    The expected SFS of rare neutral alleles with an arbitrary SFS
    (WITH RECURRENCE)

    Parameters
    ----------
    partition_mats : list of numpy arrays (boolean)
        Partition matrices give the combinations of SFS entries one must consider to compute
        an observed allele count with recurrence
    max_count : integer
        The maximum count to consider in the SFS
    Ti : numpy array (float)
        Values proportional to the expected SFS 
    Theta : float
        2N * mu, where N is the population size and mu is the mutation rate per generation

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (with recurrence)
    """
    result = np.zeros((max_count+1, max_recur)) if count_recurrence else np.ones(max_count+1)
    if count_recurrence:
        recur_vec = np.repeat(np.arange(1, max_recur+1), max_count)
        recur_vals = [np.matmul(pm, recur_vec) for pm in partition_mats]
        result[0,0] = 1
    sfs = Ti * Theta
    recur_probs = recur_expand(sfs, max_recur)
    for ii in range(1, max_count+1):
        if len(partition_mats[ii-1]) == 0:
            print("no partitions given")
            result[ii] = 0
        else:
            if count_recurrence:
                for jj in range(max_recur):
                    result[ii, jj] = np.sum(np.exp(np.matmul(partition_mats[ii-1][recur_vals[ii-1]==(jj+1),:], 
                                                             np.log(recur_probs))))
            else:
                result[ii] = np.sum(np.exp(np.matmul(partition_mats[ii-1], np.log(recur_probs))))
    return result

def sfs_llhood_Ti(data_count, data_n, Ti, Theta, max_count, max_recur, partition_mats):
    sfs = recur_sfs_Ti(partition_mats, max_count, max_recur, Ti, Theta)
    sfs_log_p = np.log(sfs) - np.log(np.sum(sfs))
    return np.sum(sfs_log_p[data_count] * data_n)

def sfs_llhood_exp(data_count, data_n, Theta, Upsilon, max_count, max_recur, partition_mats, nz=False):
    if nz:
        sfs = recur_sfs_exp(partition_mats, max_count, max_recur, Theta, Upsilon)[1:]
    else:
        sfs = recur_sfs_exp(partition_mats, max_count, max_recur, Theta, Upsilon)
    sfs_log_p = np.log(sfs) - np.log(np.sum(sfs))
    return np.sum(sfs_log_p[data_count] * data_n)

def make_bins(jmin, bb, nn, incl_zero=False):
    bins = np.arange(1, jmin + 1, dtype=int)
    bb_next = bb
    while bins[-1] >= bb_next: # catch up past starting bins
        bb_next *= bb
    bb_next = math.ceil(bb_next) if math.floor(bb_next)<=bins[-1] else math.floor(bb_next)
    bins = np.append(bins, bb_next)
    while bins[-1] < nn/2:
        bb_next *= bb
        bb_next = math.ceil(bb_next) if math.floor(bb_next)==bins[-1] else math.floor(bb_next)
        bins = np.append(bins, bb_next)
    bins = np.concatenate((bins[:-1], [math.ceil(nn/2), nn]))
    if incl_zero:
        bins = np.concatenate(([0], bins))
    return bins

def bin_means(bins):
    result = np.zeros(len(bins)-1)
    for ii in range(1, len(bins)):
        result[ii-1] = np.mean(np.arange(bins[ii-1], bins[ii]))
    return result

def bin_sizes(bins):
    result = np.zeros(len(bins)-1)
    for ii in range(1, len(bins)):
        result[ii-1] = bins[ii] - bins[ii-1]
    return result

def bin_data(ac, nn, bins):
    result = np.zeros(len(bins)-1)
    for ii in range(1, len(bins)):
        sfs_entries = ((ac >= bins[ii-1]) & (ac < bins[ii]))
        result[ii-1] = np.sum(nn[sfs_entries])
    return result

def Ti_down(m1, m2, Ti_set, ii_vals=None):
    if ii_vals is not None:
        result = np.zeros_like(ii_vals, dtype=float)
        for jj, ii in enumerate(ii_vals):
            result[jj] = np.sum(stats.binom.pmf(k=ii+1, 
                                                n=np.arange(ii+1, len(Ti_set)+1), p=m2/m1) * Ti_set[ii:])
    else:
        result = np.zeros_like(Ti_set, dtype=float)
        for ii, Ti in enumerate(Ti_set):
            result[ii] = np.sum(stats.binom.pmf(k=ii+1, n=np.arange(ii+1, len(Ti_set)+1), p=m2/m1) * Ti_set[ii:])
    return result

def combine_mu(df, mu_comb, ac="ac.nfe"):
    result_count = np.sort(np.array(list(set(np.concatenate([df[ac].to_numpy()[np.where(df.mu==mu)] for mu in mu_comb])))))
    result_n = np.zeros_like(result_count, dtype=np.float64)
    for mu in mu_comb:
        test_inds = np.where((df.mu==mu))
        test_count = df[ac].to_numpy()[test_inds]
        test_n = df["n"].to_numpy()[test_inds]
        for ii, count in enumerate(test_count):
            result_n[result_count == count] += test_n[ii]
    return result_count, result_n
