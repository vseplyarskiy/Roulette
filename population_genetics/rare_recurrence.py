import numpy as np
import scipy.special as sp

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
    Compute (no yet normalized) Poisson probabilities for each entry in the expected SFS
    up to a maximum value

    Parameters
    ----------
    sfs : numpy array (float)
        The expected SFS
    max_recur : integer
        The maximum number of independent mutations we want to consider occuring at a single site
    """
    return np.concatenate([sfs**kk/sp.factorial(kk) for kk in range(1, max_recur+1)])

def recur_sfs_exp(partition_mats, max_count, max_recur, Theta, Upsilon):
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
    Theta : float
        2N0 * mu, where N0 is current population size and mu is the mutation rate per generation
    Upsilon : float
        2N0 * gamma / m, where gamma is the exponential growth rate and m is the sample size

    Returns
    ----------
    Numpy array (float64)
        The expected SFS (with recurrence)
    """
    result = np.ones(max_count+1)
    sfs = sfs_exp(max_count, Theta, Upsilon)
    recur_probs = recur_expand(sfs, max_recur)
    for ii in range(1, max_count+1):
        if len(partition_mats[ii-1]) == 0:
            result[ii] = 0
        else:
            result[ii] = np.sum(np.exp(np.matmul(partition_mats[ii-1], np.log(recur_probs))))
    return result

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

def recur_sfs_Ti(partition_mats, max_count, max_recur, Ti, Theta):
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
    result = np.ones(max_count+1)
    sfs = Ti * Theta
    recur_probs = recur_expand(sfs, max_recur)
    for ii in range(1, max_count+1):
        if len(partition_mats[ii-1]) == 0:
            result[ii] = 0
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