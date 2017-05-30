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
import numbers

import numpy as np
import pandas as pd
import rpy2
import scipy as sp

from rpy2.robjects import r
from rpy2.robjects import pandas2ri


#################
#     SETUP     #
#################

# Constants and globals
MEAN_INTENSITY = 16.0
VAR_INTENSITY = 4.0
DEFAULT_LABELS = [
    'ModT',
    'CyberT P-value',
    'Two sample t-test',
    'One sample t-test',
    'Absolute fold change'
]
PROTEIN_LABELS = [
    'ModT median p-value',
    'CyberT median p-value',
    'Median fold change'
]


# Setup R environment
pandas2ri.activate()
r['source']('modT.r')
r['source']('bayesreg.R')
r['source']('protein.R')
# INSTRUCTIONS FOR USE
# Call functions r['modT_test'], r['bayesT'], r['find_protein_medians']
# r['modT_test'](data, "placeholder", dframe=TRUE, data.col=[DATA COLS])
# r['bayesT'](data, numC, numE)
# r['find_protein_medians'](pepdf)


#################
#  READ  DATA   #
#################

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


#################
#  SIMULATION   #
#################

def generate_samples(ctrl_data, num_to_change, fold_change, num_samples='same', use_var=None, is_log=True):
    """Simulate data based on real ctrl_data

    Creates a new dataframe of simulated intensities using empirical intensities
    Some subset of peptides have their mean fold change adjusted accordingly by MULTIPLYING the control mean
        TODO check if variance in controls is correlated with variance in treatment groups

    Args:
    ctrl_data: dataframe of raw intensities
      is_log flags whether the data is already log transformed or not
    num_to_change: number of peptides which should be perturbed (integer > 0)
    fold_change: fold change of peptides which are perturbed
    num_samples: 'same' (same number as ctrl_data) or integer > 0
    use_var: None (default: use average variance for all peptides)
               integer > 0 (rolling window variance for peptides w/ similar intensity)
               'gamma' (sample variances from gamma distribution)

    Returns: 
        len(ctrl_data) x num_samples df containing simulated exp log2 intensities

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


def _perturb_exp(exp_noise, num_to_change, fold_changes):
    """ Applies artificial perturbation to simulated experimental intensities
    """
    n = exp_noise.shape[0]
    if isinstance(fold_changes, numbers.Number):
        fold_changes = [fold_changes]
    if num_to_change * len(fold_changes) >= n:
        raise ValueError('Too many fold changes for number of peptides!')

    # For convenience, always perturb the first num_to_change peptides
    # This should not affect the analysis: can always randomize order
    perturbs = np.zeros(n, dtype=float)
    perturbs[:len(fold_changes)*num_to_change] = np.repeat(
            np.log2(fold_changes), num_to_change)
    is_changed = np.zeros(n, dtype=int)
    is_changed[:len(fold_changes)*num_to_change] = 1

    exp = exp_noise + perturbs[:,np.newaxis]

    return exp, is_changed



def sample_no_ctrl_uniform(n, num_to_change, fold_changes, var=0.06, nctrl=3, nexp=3, use_var=np.random.normal):
    """ Simulate data from uniform variance background

    Samples log intensity data from normal distribution with uniform variance
    The mean intensity for each peptide sampled from N(16, 4)
    The individual channel intensities are sampled around the mean from N(0, var)

    Args:
        n: integer number of rows to sample
        num_to_change: integer number of rows to which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of positive fold changes
                        Each fold_change will be applied num_to_change times
        var: number, fixed variance of noise
        nctrl: optional, number of control channels
        nexp: optional number of experimental channels
        use_var: sampling function to use. Should take arguments loc, scale, size

    Returns:
        ctrl, exp, is_changed
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: 0-1 Numpy vector, 1 if peptide was perturbed
    """
    # 'Spike-in' data for controls? Other paper did, not sure if useful here
    # Why not do both?

    # Draw noise and means seperately
    noise = use_var(0, var**0.5, (n, nctrl+nexp))
    # Note: always draw means from the same distribution
    avg_ctrl = np.random.normal(16, 2, n)

    background = noise + avg_ctrl[:,np.newaxis]
    ctrl, exp = background[:,:nctrl], background[:,nctrl:]
    exp, is_changed = _perturb_exp(exp, num_to_change, fold_changes)

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed


def sample_no_ctrl_gamma(n, num_to_change, fold_changes, alpha=3, beta=0.1, nctrl=3, nexp=3, use_var=np.random.normal):
    """ Simulate data with variances drawn from inverse gamma distribution

    Samples log intensity data from normal distribution with variance
    sampled from inverse gamma distribution
    The mean intensity for each peptide sampled from N(16, 4)
    The individual channel intensities are sampled around the mean from N(0, var)

    Args:
        n: integer number of rows to sample
        num_to_change: integer number of rows to which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of positive fold changes
                        Each fold_change will be applied num_to_change times
        var: number, [[beta parameter of inverse gamme]]
        nctrl: optional, number of control channels
        nexp: optional number of experimental channels
        use_var: sampling function to use. Should take arguments loc, scale, size

    Returns:
        ctrl, exp, is_changed
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: 0-1 Numpy vector, 1 if peptide was perturbed
    """

    # Sample len(n) variance vector according to inv gamma
    # InvGamma(alpha, beta) ~ 1 / Gamma(alpha, 1/beta)
    variances = 1 / np.random.gamma(alpha, 1./beta, n)
    noise = np.array([
        use_var(loc=0, scale=v**0.5, size=nctrl + nexp)
        for v in variances
    ])

    avg_ctrl = np.random.normal(16, 2, n)

    background = noise + avg_ctrl[:,np.newaxis]
    ctrl, exp = background[:,:nctrl], background[:,nctrl:]
    exp, is_changed = _perturb_exp(exp, num_to_change, fold_changes)

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed


def sample_proteins(m, num_to_change, fold_changes, peps_per_prot, var=0.06, nctrl=3, nexp=3, use_var = np.random.normal, prot_var=0.5, pep_var=3.5):
    """Simulate data from protein-level generation

    Samples log intensity data from protein-level generation
    Protein mean variance fixed at 0.5 - N(16, 0.5)
    Peptide mean variance fixed at 3.5 within protein - N(0, 3.5)
    Peptide variance sampled from NORMAL distribution - N(0, var)
    with fixed variance

    Args:
        M: integer number of proteins to sample
        num_to_change: integer number of proteins for which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of positive fold changes
                        Each fold_change will be applied num_to_change times
        peps_per_prot: int Each protein will have peps_per_prot
                OR len-n vector of ints
        var: number, square of scale parameter of variance distribution
                FOR INDIVIDUAL PEPTIDE INTENSITIES
        nctrl: optional, number of control channels
        nexp: optional, number of experimental channels
        use_var: optional, sampling function to use. 
            Should take arguments loc, scale, size
            Default is normal distribution
        prot_var
        pep_var

    Returns:
        ctrl, exp, is_changed
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: 0-1 Numpy vector, 1 if peptide comes from perturbed protein
        protein_id: 0...(m-1) Numpy vector of ints
    """
    # TODO have inverse gamma noise background for each peptide

    # Set up protein sampling by determining number of peptides for each protein
    if isinstance(peps_per_prot, int):
        peps_per_prot = [peps_per_prot] * m

    # n is total number of peptides to generate
    n = sum(peps_per_prot)
    # Create array of protein indices
    protein = np.repeat(np.arange(m), peps_per_prot)
    # indices[prot] is first index of protein
    # indices[prot+1] is last index of protein
    indices = np.cumsum(peps_per_prot)
    indices = np.insert(indices, 0, 0)
    assert len(protein) == n

    # Constant variance noise from specified background
    noise = use_var(loc=0, scale=var**0.5, size=(n, nctrl + nexp))
    
    # Draw protein and peptide averages (length-n vectors)
    # Protein means are N(16, 0.5)
    protein_means = np.repeat(np.random.normal(16, 0.5**0.5, size=m), peps_per_prot)
    # Peptide means are N(prot_mean, 3.5)
    peptide_means = np.random.normal(0, 3.5**0.5, size=n) + protein_means
    avg_ctrl = peptide_means

    background = noise + avg_ctrl[:,np.newaxis]
    ctrl, exp = background[:,:nctrl], background[:,nctrl:]

    # Validate fold changes
    if isinstance(fold_changes, numbers.Number):
        fold_changes = [fold_changes]
    if num_to_change * len(fold_changes) >= m:
        raise ValueError('Too many fold changes for number of proteins!')

    # Randomly sample set of proteins which to perturb
    prot_to_change = set(np.random.permutation(m)[:len(fold_changes)*num_to_change])
    offsets = np.repeat(np.log2(fold_changes), num_to_change)
    assert len(prot_to_change) == len(offsets)

    is_changed = np.zeros(n)
    perturb = np.zeros(n)
    # Iterate through proteins to change
    # For each protein, set the corresponding section of is_changed to 1
    # and the offset to the fold_change
    for i, prot in enumerate(prot_to_change):
        start, stop = indices[prot], indices[prot+1]
        is_changed[start:stop] = 1
        perturb[start:stop] = offsets[i]

    exp = exp + perturb[:,None]

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed, protein
     


########################
#  STATISTICAL TESTS   #
########################

def modT(ctrl, exp):
    """
    NOTE: only defined if ncol(ctrl) == 1 or ncol(ctrl) = ncol(exp)
    Be sure ctrl and exp are presented in the same order
    both ctrl and exp should be log2 ratios
    """
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

    return res


def cyberT(ctrl, exp):
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')
    ctrl.reset_index(drop=True, inplace=True)
    exp.reset_index(drop=True, inplace=True)

    df = pd.concat([ctrl, exp], axis=1, ignore_index=True)
    df.columns = (['C%d' % (i+1) for i in xrange(ctrl.shape[1])] +
                  ['E%d' % (i+1) for i in xrange(exp.shape[1])])
    res = r['bayesT'](df, numC = ctrl.shape[1], numE = exp.shape[1], doMulttest=True)
    res = pandas2ri.ri2py(res)

    return res


def t_test(ctrl, exp, ratio_test = False):
    """
    Ratio test flag indicates if t test should be performed as 1-sample t test
    """

    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    if ratio_test:
        if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
            data = np.array(exp.values - ctrl.values)
            _, pvals = sp.stats.ttest_1samp(data, 0, axis=1)
        else:
            raise ValueError('Not valid number of ctrl columns, see documentation')
    else:
        _, pvals = sp.stats.ttest_ind(ctrl, exp, axis=1)

    return pvals


def protein_rollup(protein_df):
    """
    Runs cyberT on peptides and rolls up to protein rollup using R script

    Args:
        ctrl, exp should be dataframes of PEPTIDE level measurements
        protein is vector which designates the protein for each row of ctrl, exp

    Returns (WIP): 
        
    """
    protein_df = r['find_protein_medians'](protein_df, use_isoform=True)
    out = pandas2ri.ri2py(protein_df)

    return out


###############
#    MISC     #
###############

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


def do_stat_tests(ctrl, exp):
    """
    Runs modT, cyberT, t-test, fold_change analysis

    Returns modT_pvals,
            cyberT_pvals, 
            cyberT_ppde (+), 
            ttest_pvals,
            ttest_ratio_pvals,
            fold_change (+)
    Any of these may be None if number of channels is not suitable
    (+) = Proper measure is inverted: i.e. x* = max(x) - x
    """

    do_ratio = (ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1])
    do_t = (ctrl.shape[1] > 1)
   
    if do_ratio:
        modT_res = modT(ctrl, exp)
        modT_pvals = modT_res['P.Value']
        print "Ran moderated T test"

        ttest_ratio_pvals = t_test(ctrl, exp, ratio_test=True)
        print "Ran one sample t test"
    else:
        modT_pvals = None
        print "Skipped moderated T test, dimensions not suitable"

        ttest_ratio_pvals = None
        print "Skipped one sample t test, dimensions not suitable"

    if do_t:
        cyberT_res = cyberT(ctrl, exp)
        cyberT_pvals = cyberT_res['pVal']
        print "Ran cyberT test"

        ttest_pvals = t_test(ctrl, exp)
        print "Ran two sample t test"
    else:
        cyberT_pvals = None
        cyberT_ppde = None
        print "Skipped cyberT, too few channels"

        ttest_pvals = None
        print "Skipped two sample t test, too few channels"

    fold_change = np.abs(np.mean(ctrl.values - exp.values, axis=1))
   
    return (modT_pvals,
            cyberT_pvals,
            ttest_pvals,
            ttest_ratio_pvals,
            np.max(fold_change) + 0.01 - fold_change)


def do_stat_tests_protein(ctrl, exp, protein):
    """
    Runs modT, cyberT, t-test, fold_change on protein level

    Returns:
        (modT_pvals (median)
         cyberT_pvals (median)
         fold_change (median)
         
    (+) = Proper measure is inverted: i.e. x* = max(x) - x
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')
    
    ctrl.columns = ['C%d' % (i+1) for i in xrange(ctrl.shape[1])]
    exp.columns = ['E%d' % (i+1) for i in xrange(exp.shape[1])]

    (modT_pvals,
     cyberT_pvals,
     ttest_pvals,
     ttest_ratio_pvals,
     fold_change) = do_stat_tests(ctrl, exp)

    # Revert fold change inversion
    fold_change = np.max(fold_change) + 0.01 - fold_change

    pval_df = pd.DataFrame({
        'accession_number': protein,
        'fold_change': fold_change,
        'modT_PVal': modT_pvals,
        'cyberT_PVal': cyberT_pvals,
        'ttest_PVal': ttest_pvals
    })

    ctrl.reset_index(drop=True, inplace=True)
    exp.reset_index(drop=True, inplace=True)
    pval_df.reset_index(drop=True, inplace=True)

    meanC = np.mean(ctrl.values, axis=1)
    meanE = np.mean(exp.values, axis=1)

    df = pd.concat([ctrl,exp, pval_df], axis=1)
    df['meanC'] = meanC
    df['meanE'] = meanE

    out = protein_rollup(df)
    fc = out['fold_change.med']
    fc = np.max(fc) + 0.01 - fc

    return (
        out['modT_PVal.med'],
        out['cyberT_PVal.med'],
        fc
    ), out


if __name__ == "__main__":
    # quick_modT, quick_cyberT = setup()
    pass
