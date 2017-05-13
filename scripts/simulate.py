# TODO
# Power analysis for number of samples required for certain power and FDR
#   Use LR test to determine the theoretical optimum power
# Compare runtime of R version of script with Python version

"""
Pick a subset of data and push data up: increase by 0.25 to 5 fold
Bin peptide intensity: randomly choose 10-20 percent, spike upreg, downreg. Abs mean value comparison between the two peptides. Use control samples to estimate
    Randomly choose 10-20 percent, run simulation with many parameters
    Generate many datasets using this program with known up and down regulation.
    Need to define variance, variance models using common variance vs inverse gamma variance
"""
import math
import numpy as np
import pandas as pd
import rpy2
import scipy as sp

from rpy2.robjects import r
from rpy2.robjects import pandas2ri

from roc import *


# Setup R environment
pandas2ri.activate()
r['source']('modT.r')
r['source']('bayesreg.R')
# Call functions r['modT_test'] and r['bayesT']
# modT_test(data, "placeholder", dframe=TRUE, data.col=[DATA COLS])
# bayesT(data, numC, numE)


def read_data(filename='../data/Neurogranin_KD_PSM_021517.ssv'):
    # Reads in the PSM report and returns raw intensities
    cols_to_use = ['sequence', 'TMT_126', 'TMT_127', 'TMT_128', 'TMT_129', 'TMT_130', 'TMT_131']
    df = pd.read_table(filename,
            sep=';',
            header=0,
            usecols=cols_to_use)
    df = df.dropna(axis=0, how='any')
    sequence = df['sequence']
    data = np.log2(df[cols_to_use[1:]])
    # data = df[cols_to_use[1:]]
    data.columns = ['E1', 'E2', 'E3', 'C1', 'C2', 'C3']
    data['sequence'] = sequence
    return data


def generate_samples(ctrl_data, num_to_change, fold_change, num_samples='same', use_var=None, is_log=True):
    """
    Creates a new dataframe of simulated intensities using empirical intensities
    Some subset of peptides have their mean fold change adjusted accordingly by MULTIPLYING the control mean
        TODO check if variance in controls is correlated with variance in treatment groups

    - ctrl_data: dataframe of raw intensities
        is_log flags whether the data is already log transformed or not
    - num_to_change: number of peptides which should be perturbed (integer > 0)
    - fold_change: fold change of peptides which are perturbed
    - num_samples: 'same' (same number as ctrl_data) or integer > 0
    - use_var: None (default: use average variance for all peptides)
               integer > 0 (rolling window variance for peptides w/ similar intensity)
               'gamma' (sample variances from gamma distribution)

    Returns: len(ctrl_data) x num_samples df containing simulated exp log2 intensities

    Compare modT (one ctrl, many exp) vs cyberT (half and half ctrl vs exp)
    """
    if num_samples == 'same':
        num_samples = ctrl_data.shape[1]
    if num_samples <= 0:
        raise ValueError

    if not is_log:
        ctrl_data = np.log2(ctrl_data)

    n = ctrl_data.shape[0]

    # Create noisy data from distribution
    if use_var is None:
        # Default: uses same variance for all samples, normal dist
        print "Using average variance across all samples"

        variances = np.var(ctrl_data, axis=1)
        avg_var = np.mean(variances)
        print "Average variance: %.2f" % avg_var
        noise = np.random.normal(0.0, avg_var**0.5, (n, num_samples))

    elif isinstance(use_var, int):
        # If use_var is an integer, use a rolling window of average variances
        # around it in terms of intensities
        # Then sample each sample according to normal with given dist
        if use_var <= 0 or use_var > len(variances):
            raise ValueError
        print "Using rolling window of variance of size %d" % use_var

        # Sort in order of mean intensity
        mean_int = np.mean(ctrl_data, axis=1)
        reorder = np.argsort(mean_int)
        ctrl_data = ctrl_data.iloc[reorder]
        # Rolling window average of variance
        variances = pd.Series(np.var(ctrl_data, axis=1))
        roll_var = pd.rolling_mean(variances, use_var)
        roll_var[:use_var-1] = roll_var[use_var-1]
        # Create noise 
        noise = np.array([
            np.random.normal(0.0, std, num_samples)
            for std in roll_var**0.5])

    elif use_var == 'gamma':
        # TODO use inverse gamma distribution for noise?
        # Use average variance and variance of variances
        # To interpolate inverse gamma distribution
        pass
        
    # Sample subset of peptides to change randomly
    sample_idxs = set(np.random.choice(n, num_to_change, replace=False))
    offset = math.log(fold_change, 2)
    avg_ctrl = np.mean(ctrl_data, axis=1)
    tmp = avg_ctrl + np.array([offset if i in sample_idxs else 0 for i in xrange(n)])
    out = pd.DataFrame(noise + tmp[:,np.newaxis])

    # Truncate measurements with negative log values to arbitrary small values
    out[out <= 1] = 1.
    
    out.columns = [('E%d' % (i+1)) for i in xrange(num_samples)]
    is_changed = np.array([1 if i in sample_idxs else 0 for i in xrange(n)])

    return out, is_changed


def generate_samples_no_ctrl_uniform(n, num_to_change, fold_change, var):
    noise = None
    ctrl_data = None
    sample_idxs = set(np.random.choice(n, num_to_change, replace=False))
    offset = math.log(fold_change, 2)
    avg_ctrl = np.mean(ctrl_data, axis=1)
    tmp = avg_ctrl + np.array([offset if i in sample_idxs else 0 for i in xrange(n)])
    out = pd.DataFrame(noise + tmp[:,np.newaxis])


def modT(ctrl, exp):
    # NOTE: only defined if ncol(ctrl) == 1 or ncol(ctrl) = ncol(exp)
    # Be sure ctrl and exp are presented in the same order
    # both ctrl and exp should be log2 ratios
    if ctrl.shape[0] != exp.shape[0]:
        print ctrl.shape
        print exp.shape
        raise ValueError('Length of exp and ctrl data frames not identical')

    if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
        data = pd.DataFrame(exp.values - ctrl.values)
    else:
        raise ValueError('Not valid number of ctrl columns, see documentation')
    # Also create id col for convenience
    data.columns = ['Log2_Ratio_%d' % i for i in xrange(data.shape[1])]
    data_cols = data.columns
    data['id'] = np.arange(len(data))

    res = r['modT_test'](data, "placeholder", id_col='id', data_col=data_cols, dframe=True)
    res = pandas2ri.ri2py(res)

    # TODO filter this data by the relevant columns
    return res


def cyberT(ctrl, exp):
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')
    
    df = pd.concat([ctrl, exp], axis=1)
    df.columns = (['C%d' % i for i in xrange(ctrl.shape[1])] +
                  ['E%d' % i for i in xrange(exp.shape[1])])
    res = r['bayesT'](df, numC = ctrl.shape[1], numE = exp.shape[1], doMulttest=True)
    res = pandas2ri.ri2py(res)

    # TODO filter this data by the relevant columns
    return res


def setup():
    # Convenience functions which encapsulate the ctrl data for quicker access
    ctrl_data = read_data()[['C1', 'C2', 'C3']]

    def quick_modT(n, fold_change, num_to_change=None, **kwargs):
        # If n <= 0, use entire ctrl data field
        # If num_to_change = None: 1/10 of input data
        # Kwargs are passed through to simulation code
        print "Running ModT on cached data: recall setup() to refresh dataset"

        if num_to_change is None:
            num_to_change = n // 10

        subset = np.random.choice(len(ctrl_data), n, replace=False)
        ctrl_data_sub = ctrl_data.iloc[subset]
        print "isnan: %d" % np.sum(np.isnan(ctrl_data_sub).values)
        exp_data, is_changed = generate_samples(ctrl_data_sub, num_to_change, fold_change, **kwargs)
        modT_res = modT(ctrl_data_sub, exp_data)

        p_vals = modT_res['P.Value']

        return (is_changed, p_vals, modT_res)

    def quick_cyberT(n, fold_change, num_to_change=None, **kwargs):
        # If n <= 0, use entire ctrl data field
        # If num_to_change = None: 1/10 of input data
        # Kwargs are passed through to simulation code
        print "Running ModT on cached data: recall setup() to refresh dataset"

        if num_to_change is None:
            num_to_change = n // 10

        subset = np.random.choice(len(ctrl_data), n, replace=False)
        ctrl_data_sub = ctrl_data.iloc[subset]
        print "isnan: %d" % np.sum(np.isnan(ctrl_data_sub).values)
        exp_data, is_changed = generate_samples(ctrl_data_sub, num_to_change, fold_change, **kwargs)
        cyberT_res = cyberT(ctrl_data_sub, exp_data)

        p_vals = cyberT_res['pVal']

        return (is_changed, p_vals, cyberT_res)


    print "Finished setting up with default data"
    return quick_modT, quick_cyberT

quick_modT, quick_cyberT = setup()

if __name__ == "__main__":
    # main()
    pass
