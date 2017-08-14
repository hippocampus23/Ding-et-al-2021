import math
import numbers

import numpy as np
import pandas as pd
import rpy2
import scipy as sp
import statsmodels.api as sm

from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from scipy import stats
from statsmodels.tools import add_constant


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
DEFAULT_LABELS_MODT_2SAMP = DEFAULT_LABELS + ['ModT (2-sample)']
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
# r['proteinBayesT'](data, numC, numE, pool_intensity=True)


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


def _perturb_exp(exp_noise, num_to_change, fold_changes, binary_labels):
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
            fold_changes, num_to_change)
    is_changed = np.zeros(n, dtype=(int if binary_labels else float))
    is_changed[:len(fold_changes)*num_to_change] = (1 if binary_labels else
            np.repeat(fold_changes, num_to_change))

    exp = exp_noise + perturbs[:,np.newaxis]

    return exp, is_changed



def sample_no_ctrl_uniform(n, num_to_change, fold_changes, var=0.06, nctrl=3, nexp=3, use_var=np.random.normal, binary_labels=True):
    """ Simulate data from uniform variance background

    Samples log intensity data from normal distribution with uniform variance
    The mean intensity for each peptide sampled from N(16, 4)
    The individual channel intensities are sampled around the mean from N(0, var)

    Args:
        n: integer number of rows to sample
        num_to_change: integer number of rows to which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of log2 fold changes
                        Each fold_change will be applied num_to_change times
        var: number, fixed variance of noise
        nctrl: optional, number of control channels
        nexp: optional number of experimental channels
        use_var: sampling function to use. Should take arguments loc, scale, size
        binary_labels: bool, whether to return binary indicators (0-1) or scalar
            (log2 true fc) for each peptide

    Returns:
        ctrl, exp, is_changed
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: length-n Numpy array indicating perturbation of peptide
            (see binary_labels argument)
    """

    # Draw noise and means seperately
    noise = use_var(0, var**0.5, (n, nctrl+nexp))
    # Note: always draw means from the same distribution
    avg_ctrl = np.random.normal(16, 2, n)

    background = noise + avg_ctrl[:,np.newaxis]
    ctrl, exp = background[:,:nctrl], background[:,nctrl:]
    exp, is_changed = _perturb_exp(exp, num_to_change, fold_changes, binary_labels)

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed


def sample_no_ctrl_gamma(n, num_to_change, fold_changes, alpha=3, beta=0.1, nctrl=3, nexp=3, use_var=np.random.normal, binary_labels=True):
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
        fold_changes: 1D array of log2 fold changes
                        Each fold_change will be applied num_to_change times
        alpha, beta: parameters controlling shape of inv gamma. Default 3, 0.1
        nctrl: optional, number of control channels
        nexp: optional number of experimental channels
        use_var: sampling function to use. Should take arguments loc, scale, size
        binary_labels: bool, whether to return binary indicators (0-1) or scalar
            (log2 true fc) for each peptide

    Returns:
        ctrl, exp, is_changed
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: length-n Numpy array indicating perturbation of peptide
            (see binary_labels argument)
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
    exp, is_changed = _perturb_exp(exp, num_to_change, fold_changes, binary_labels)

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed


def sample_proteins(m, num_to_change, fold_changes, peps_per_prot, var=0.06, nctrl=3, nexp=3, use_var = np.random.normal, prot_var=0.5, pep_var=3.5, background = "G", alpha=3, beta=0.1, binary_labels=True):
    """Simulate data from protein-level generation

    Samples log intensity data from protein-level generation
    Protein mean variance defaults to 0.5 - N(16, 0.5)
    Peptide mean variance defaults to 3.5 within protein - N(0, 3.5)
    Peptide variance sampled from NORMAL distribution - N(0, var)
    with fixed variance

    Args:
        m: integer number of proteins to sample
        num_to_change: integer number of proteins for which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of log2 fold changes
                        Each fold_change will be applied num_to_change times
        peps_per_prot: int Each protein will have peps_per_prot
                OR len-n vector of ints
        var: number, square of scale parameter of variance distribution
                FOR INDIVIDUAL PEPTIDE INTENSITIES (only if uniform background)
        nctrl: optional, number of control channels
        nexp: optional, number of experimental channels
        use_var: optional, sampling function to use.
            Should take arguments loc, scale, size
            Default is normal distribution
        prot_var: optional, variance of protein means.
        pep_var: optional, variance of peptide means given protein mean
        background: optional, how variance is sampled (must be "G" or "U")
        binary_labels: bool, whether to return binary indicators (0-1) or scalar
            (log2 true fc) for each peptide

    Returns:
        ctrl, exp, is_changed, protein_id
        ctrl: ctrl intensities (n x nctrl Pandas df)
        exp: experimental intensities (n x nexp Pandas df)
        is_changed: 0-1 Numpy vector, 1 if peptide comes from perturbed protein
        protein_id: 0...(m-1) Numpy vector of ints
    """

    # Set up protein sampling by determining number of peptides for each protein
    if isinstance(peps_per_prot, int):
        peps_per_prot = [peps_per_prot] * m
    elif m % len(peps_per_prot) == 0:
        peps_per_prot = np.repeat(peps_per_prot, m // len(peps_per_prot))

    # n is total number of peptides to generate
    n = sum(peps_per_prot)
    # Create array of protein indices
    protein = np.repeat(np.arange(m, dtype=int), peps_per_prot)
    # indices[prot] is first index of protein
    # indices[prot+1] is last index of protein
    indices = np.cumsum(peps_per_prot)
    indices = np.insert(indices, 0, 0)
    assert len(protein) == n

    if background == "U":
        # Constant variance noise
        noise = use_var(loc=0, scale=var**0.5, size=(n, nctrl + nexp))
    elif background == "G":
        # Inverse gamma background
        variances = 1 / np.random.gamma(alpha, 1./beta, n)
        noise = np.array([
            use_var(loc=0, scale=v**0.5, size=nctrl + nexp)
            for v in variances
            ])
    else:
        raise ValueError("Invalid setting of background")

    # Draw protein and peptide averages (length-n vectors)
    # Protein means are N(16, 0.5)
    protein_means = np.repeat(
            np.random.normal(16, prot_var**0.5, size=m), peps_per_prot)
    # Peptide means are N(prot_mean, 3.5)
    peptide_means = np.random.normal(0, pep_var**0.5, size=n) + protein_means
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
    offsets = np.repeat(fold_changes, num_to_change)
    assert len(prot_to_change) == len(offsets)

    is_changed = np.zeros(n)
    perturb = np.zeros(n)
    # Iterate through proteins to change
    # For each protein, set the corresponding section of is_changed to 1
    # and the offset to the fold_change
    for i, prot in enumerate(prot_to_change):
        start, stop = indices[prot], indices[prot+1]
        is_changed[start:stop] = 1 if binary_labels else 2**offsets[i]
        perturb[start:stop] = offsets[i]

    exp = exp + perturb[:,None]

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed, protein


def sample_phospho(num_peps, num_to_change, fold_changes, protein_df, use_var=np.random.normal, var=0.06, alpha=3, beta=0.1, nctrl=3, nexp=3, background="U", binary_labels=True):
    """ Sample phosphoproteins based on protein fold changes
    
    Samples num_peps with nctrl and nexp channels

    Args:
        num_peps: Number of phosphopeptides to generate
        num_to_change: integer number of peptides for which each fold change
                        perturbation should be applied
                        NOTE: num_to_change * len(fold_changes) must be less than n
        fold_changes: 1D array of log2 fold changes
                        Each fold_change will be applied num_to_change times
        protein_df: DataFrame summarizing protein level information, with cols
                'protein_id': unique key for proteins
                'fold_change': True fold change between ctrl and exp protein
        var: number, square of scale parameter of variance distribution
                FOR INDIVIDUAL PEPTIDE INTENSITIES (only if uniform background)
        nctrl: optional, number of control channels
        nexp: optional, number of experimental channels
        use_var: optional, sampling function to use.
            Should take arguments loc, scale, size
            Default is normal distribution
        background: optional, how variance is sampled (must be "G" or "U")
        binary_labels: bool, whether to return binary indicators (0-1) or scalar
            (log2 true fc) for each peptide

    return ctrl, exp, is_changed, peptide_prot_indices
    """
    # m is number of proteins
    m = protein_df.shape[0]
    # if isinstance(num_peps, int):
    #     num_peps = [num_peps] * m
    # elif m % len(num_peps) == 0:
    #     num_peps = np.repeat(num_peps, m // len(num_peps))
    # n = sum(num_peps)

    # For now, generate peptides and randomly assign to proteins
    # TODO: should mean of peptide depend on mean of protein?
    # Really shouldn't make a difference
    # Would have to pass in variance/mean information from protein samples
    n = num_peps
    peptide_prot_indices = np.random.choice(m, size=n, replace=True)
    avg_ctrl = np.random.normal(16, 2, n)    

    if background == "U":
        # Constant variance noise
        noise = use_var(loc=0, scale=var**0.5, size=(n, nctrl + nexp))
    elif background == "G":
        # Inverse gamma background
        variances = 1 / np.random.gamma(alpha, 1./beta, n)
        noise = np.array([
            use_var(loc=0, scale=v**0.5, size=nctrl + nexp)
            for v in variances
            ])

    background = noise + avg_ctrl[:,np.newaxis]
    ctrl, exp = background[:,:nctrl], background[:,nctrl:]
    # Adjust exp for protein level fold changes
    exp += protein_df['fold_change'].values[peptide_prot_indices,np.newaxis]
    exp, is_changed = _perturb_exp(exp, num_to_change, fold_changes, binary_labels)

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed, peptide_prot_indices


########################
#  STATISTICAL TESTS   #
########################

def modT(ctrl, exp, robust=True):
    """
    NOTE: only defined if ncol(ctrl) == 1 or ncol(ctrl) = ncol(exp)
    Be sure ctrl and exp are presented in the same order
    both ctrl and exp should be log2 ratios
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
        data = pd.DataFrame(exp.values - ctrl.values)
    else:
        raise ValueError('Not valid number of ctrl columns, see documentation')
    # Also create id col for convenience
    data.columns = ['Log2_Ratio_%d' % i for i in xrange(data.shape[1])]
    data_cols = data.columns
    data['id'] = np.arange(len(data))

    res = r['modT_test'](data, "placeholder", id_col='id', data_col=data_cols, dframe=True, robust=robust)
    res = pandas2ri.ri2py(res)
    res['std_err'] = (res['CI.R'] - res['CI.L'])/3.92

    return res


def modT_2sample(ctrl, exp, robust=True):
    """
    Runs moderated T with 2 sample t test
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    data = pd.DataFrame(np.concatenate((ctrl.values, exp.values), axis=1))
    data.columns = (['C%d' % (i+1) for i in xrange(ctrl.shape[1])] +
            ['E%d' % (i+1) for i in xrange(exp.shape[1])])
    data_cols= data.columns
    data['id'] = np.arange(len(data))

    design = np.array(([-1] * ctrl.shape[1]) + ([1] * exp.shape[1]), dtype=int)

    res = r['modT_test'](data, "placeholder", id_col='id', data_col=data_cols, dframe=True, design=design, robust=robust)
    res = pandas2ri.ri2py(res)
    res['std_err'] = (res['CI.R'] - res['CI.L'])/3.92

    return res


def cyberT(ctrl, exp, **kwargs):
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')
    ctrl.reset_index(drop=True, inplace=True)
    exp.reset_index(drop=True, inplace=True)

    df = pd.concat([ctrl, exp], axis=1, ignore_index=True)
    df.columns = (['C%d' % (i+1) for i in xrange(ctrl.shape[1])] +
            ['E%d' % (i+1) for i in xrange(exp.shape[1])])
    res = r['bayesT'](df, numC = ctrl.shape[1], numE = exp.shape[1], ppde=False, doMulttest=True, **kwargs)
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
            _, pvals = stats.ttest_1samp(data, 0, axis=1)
        else:
            raise ValueError('Not valid number of ctrl columns, see documentation')
    else:
        _, pvals = stats.ttest_ind(ctrl, exp, axis=1)

    return pvals


def protein_wls_test(ctrl, exp, protein, onesample=True, use_bayes=False, reg=10, adj_var=False):
    # Do cyberT
    cyberT_peps = cyberT(ctrl, exp)
    cyberT_peps['protein'] = protein
    return protein_wls_test_cyberT(
            cyberT_peps, onesample=onesample, use_bayes=use_bayes, reg=reg, adj_var=adj_var)


def protein_wls_test_cyberT(cyberT_peps, onesample=True, use_bayes=False, reg=10, adj_var=False):
    """ TODO documentation
    Especially what columns are required from cyberT_peps
    and what columns are in output

    reg = number >= 1: Default 10, prior weight on regularization std
                Set to None for no regularization
    """
    # Group by protein
    grouped = cyberT_peps.groupby('protein', sort=False)
    res = pd.DataFrame(
            columns=['protein_id','fold_change','raw_std_err', 'pval_no_reg', 'n'],
            index=np.arange(len(grouped)))
    res['protein_id'] = res['protein_id'].astype(str)
    res['fold_change'] = res['fold_change'].astype(float)
    res['raw_std_err'] = res['raw_std_err'].astype(float)
    res['pval_no_reg'] = res['pval_no_reg'].astype(float)
    res['n'] = 0
    for j, (name, group) in enumerate(grouped):
        # Do WLS
        if onesample:
            fold_change, std_err, pval = _wls_onesample(
                    group, use_bayes=use_bayes, adj_var=adj_var)
        else:
            fold_change, std_err, pval = _wls(
                    group, use_bayes=use_bayes)

        nC, nE = cyberT_peps['nC'].astype(float), cyberT_peps['nE'].astype(float)
        if use_bayes:
            stdC, stdE = cyberT_peps['bayesSDC'], cyberT_peps['bayesSDE']
        else:
            stdC, stdE = cyberT_peps['stdC'], cyberT_peps['stdE']
        res.ix[j] = [name, fold_change, std_err, pval, group.shape[0]]
        # If we regularize, adjust the std
    
    # Add regularization
    # Smooth the protein level variances by moving toward mean
    if reg is not None:
        # Rolling average of the protein level standard deviations
        WINSIZE = 101  # Window size: must be odd positive integer
        m = res.shape[0]
        if m <= WINSIZE:
            p_sd = np.ones(m, dtype=float)*np.mean(res['raw_std_err'].values)
        else:
            # Rolling average sorted in order of fold change
            sorted_std_err = res['raw_std_err'].reindex(
                    np.argsort(res['fold_change'].values))
            p_sd = sorted_std_err.rolling(window=WINSIZE, center=True).mean()
            # Pad ends to avoid NaNs
            pad_size = (WINSIZE-1) / 2
            p_sd.iloc[:pad_size] = p_sd.iloc[pad_size]
            p_sd.iloc[-pad_size:] = p_sd.iloc[-pad_size-1]
            # Now make sure we are in protein id order again
            p_sd.sort_index(inplace=True)

        std_err = ((reg*(p_sd**2) + (res['n']-1)*(res['raw_std_err']**2)) / 
                   (reg+res['n']-1))**0.5
        # Two-tailed test
        p_val = 2*stats.t.cdf(-abs(res['fold_change']) / std_err, df=res['n']-1)
        res['std_err'] = std_err
        res['pval'] = p_val

        # Copy over degenerate pvalues if n=1
        res.loc[res['n'] == 1,'pval'] = res.loc[res['n'] == 1,'pval_no_reg']
        res.loc[res['n'] == 1,'std_err'] = res.loc[res['n'] == 1,'raw_std_err']
    else:
        res['std_err'] = res['raw_std_err']
        res['pval'] = res['pval_no_reg']

    return res


def _wls_degenerate(data, use_bayes):
    """ Weighted least squares with just one peptide

    Peptide mean directly estimates protein mean
    Valid for both one and two sample approach
    data should be df with one row
    """
    data = data.ix[0]
    nC, nE = data['nC'], data['nE']
    if use_bayes:
        stdC, stdE = data['bayesSDC'], data['bayesSDE']
    else:
        stdC, stdE = data['stdC'], data['stdE']
    return (data['meanE'] - data['meanC'],
            (((nC-1)*(stdC**2) + (nE-1)*(stdE**2)) /
                (nC+nE-2)) * (float(nC+nE) / (nC*nE)),
            data['pVal'])

def _wls(data, use_bayes=False):
    """ Weighted least squares for peptides in protein.
    DEPRECATED, do not use
    Args:
        data - df with columns including
            meanC, meanE, stdC, stdE, bayesSDC, bayesSDE
        use_bayes - bool, def False. If True, use bayesSDC/bayesSDE
        reg - (mean, n_obs) mean, pseudocount of of prior observation
            n_obs > 0
    Returns:
        (beta, std_err, pval, n)
    """
    # Degenerate case, only one peptide
    # Regularization doesn't affect answer here
    if data.shape[0] == 1:
        return _wls_degenerate(data, use_bayes)

    y = data['meanC'].append(data['meanE']).values
    if use_bayes:
        w = data['bayesSDC'].append(data['bayesSDE']).values**2
    else:
        w = data['stdC'].append(data['stdE']).values**2
    x = np.ones(data.shape[0]*2)
    x[:data.shape[0]] = 0

    mod_wls = sm.WLS(y, add_constant(x, prepend=False), weights=1./w)
    res_wls = mod_wls.fit()
    return (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0], data.shape[0])

def _wls_onesample(data, use_bayes=False, adj_var=True):
    """ Weighted least squares for one-sample peptides in protein
    Args:
        data - df with columns including
            meanC, meanE, stdC, stdE, bayesSDC, bayesSDE
        use_bayes - bool, def False. If True, use bayesSDC/bayesSDE
        reg - (std, n_obs) standard dev, pseudocount of of prior on std dev
            std, n_obs > 0
    Returns:
        (beta, std_err, pval)
    """
    n = data.shape[0]
    # Degenerate case, only one peptide
    # Regularization doesn't affect answer here
    if n == 1:
        return _wls_degenerate(data, use_bayes)

    y = data['meanE'].values - data['meanC'].values
    x = np.ones(data.shape[0])
    nC, nE = data['nC'].values.astype(float), data['nE'].values.astype(float)
    if use_bayes:
        stdC, stdE = data['bayesSDC'].values, data['bayesSDE'].values
    else:
        stdC, stdE = data['stdC'].values, data['stdE'].values
    # Variance is additive for control and experimental conditions
    w = 1. / (stdC**2 + stdE**2)
    # Weighted least squares with only constant term
    mod_wls = sm.WLS(y, x, weights=w)
    res_wls = mod_wls.fit()
    beta_hat, std_err, p_val = (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0])
    if adj_var:
        # Adjust the standard error for peptide level uncertainty
        # = reciprocal of sum of reciprocals of peptide variances
        std_err += 1. / sum(w)
        p_val = 2*stats.t.cdf(-abs(beta_hat) / std_err, df=(n-1))

    return (beta_hat, std_err, p_val) 


###############
#    MISC     #
###############

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


def do_stat_tests(ctrl, exp, run_modT_2sample=False):
    """
    Runs modT, cyberT, t-test, fold_change analysis

    Returns modT_pvals,
            cyberT_pvals,
            ttest_pvals,
            ttest_ratio_pvals,
            fold_change (+),
            [modT_2sample]
    Any of these may be None if number of channels is not suitable
    (+) = Proper measure is inverted: i.e. x* = max(x) - x
    """

    do_ratio = (ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]) and exp.shape[1] > 1
    do_t = (ctrl.shape[1] > 1)

    if do_ratio:
        modT_res = modT(ctrl, exp)
        modT_pvals = modT_res['P.Value'].values
        print "Ran moderated T test"

        ttest_ratio_pvals = t_test(ctrl, exp, ratio_test=True)
        print "Ran one sample t test"
    else:
        modT_pvals = np.nan
        print "Skipped moderated T test, dimensions not suitable"

        ttest_ratio_pvals = np.nan
        print "Skipped one sample t test, dimensions not suitable"

    if do_t:
        ttest_pvals = t_test(ctrl, exp)
        print "Ran two sample t test"
    else:
        ttest_pvals = np.nan
        print "Skipped two sample t test, too few channels"

    cyberT_res = cyberT(ctrl, exp)
    cyberT_pvals = cyberT_res['pVal'].values
    print "Ran cyberT test"

    fold_change = np.abs(np.mean(ctrl.values, axis=1) - np.mean(exp.values, axis=1))

    res = pd.DataFrame({'modT': modT_pvals,
        'cyberT':cyberT_pvals,
        't-test': ttest_pvals,
        't test (1-sample)': ttest_ratio_pvals,
        'fold change': np.max(fold_change) + 0.01 - fold_change})
    if run_modT_2sample:
        modT_2sample_pvals = modT_2sample(ctrl, exp)['P.Value'].values
        res['modT (2-sample)'] = modT_2sample_pvals

    return res


def do_stat_tests_protein(ctrl, exp, protein):
    """
    Runs modT, cyberT, t-test, fold_change on protein level

    Returns:
        (modT_pvals (median)
         cyberT_pvals (median)
         fold_change (median)
         ), df

    pvals are all series, sorted in increasing order of protein id
    df has columns for pvals. Can pass directly into plot visualization

    TODO implement protein level rollup using variance background

    (+) = Proper measure is inverted: i.e. x* = max(x) - x
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    ctrl.columns = ['C%d' % (i+1) for i in xrange(ctrl.shape[1])]
    exp.columns = ['E%d' % (i+1) for i in xrange(exp.shape[1])]

    ## Do stat tests individually
    # (modT_pvals,
    #  cyberT_pvals,
    #  ttest_pvals,
    #  ttest_ratio_pvals,
    #  fold_change) = do_stat_tests(ctrl, exp)
    # Revert fold change inversion
    # fold_change = np.max(fold_change) + 0.01 - fold_change
    ctrl_mean = ctrl.groupby(protein).median()
    exp_mean = exp.groupby(protein).median()

    modT_pvals = modT(ctrl_mean, exp_mean)['P.Value']
    cyberT_pvals = cyberT(ctrl_mean, exp_mean)['pVal']
    ttest_pvals = t_test(ctrl_mean, exp_mean)
    fold_change = np.abs(
            np.mean(ctrl_mean.values, axis=1)
            - np.mean(exp_mean.values, axis=1)
    )

    pval_df = pd.DataFrame({
        'protein_id': ctrl_mean.index,
        'fold_change_med': fold_change,
        'modT_PVal_med': modT_pvals,
        'cyberT_PVal_med': cyberT_pvals,
        'ttest_PVal_med': ttest_pvals,
        })

    # Do naive aggregates of peptide level tests
    # Sorts in increasing order of id
    out = pval_df.groupby('protein_id').median()

    # ctrl.reset_index(drop=True, inplace=True)
    # exp.reset_index(drop=True, inplace=True)
    # meanC, meanE = cyberT_res['meanC'].values, cyberT_res['meanE'].values
    meanC = np.mean(ctrl.values, axis=1)
    meanE = np.mean(exp.values, axis=1)

    # Note that all these are returned in ascending order of id
    # CyberT by peptide
    protein_cyberT_bypep = r['bayesT.pair'](
            pd.DataFrame.from_items([
                ('ratio', meanE - meanC),
                ('index', meanE + meanC),
            ]),
            1,
            aggregate_by=protein)
    protein_cyberT_bypep = pandas2ri.ri2py(protein_cyberT_bypep)
    protein_ttest_bypep = r['bayesT.pair'](
            pd.DataFrame.from_items([
                ('ratio', meanE - meanC),
                ('index', meanE + meanC),
            ]),
            1,
            bayes=False,
            aggregate_by=protein)
    protein_ttest_bypep = pandas2ri.ri2py(protein_ttest_bypep)

    cyberT_res = cyberT(ctrl, exp)
    cyberT_res['protein'] = protein
    # Testing different regularization coefficients
    # for reg in [2, 4, 8, 12, 16]:
    #     out['wls_pval_reg_%02d' % reg] = protein_wls_test_cyberT(
    #             cyberT_res, onesample=True, reg=reg)['pval'].values
    out['wls'] = protein_wls_test_cyberT(
            cyberT_res, use_bayes=True, onesample=True)['pval'].values

    # out['wls_pval'] = wls_pval_reg['pval']
    out['cyberT_bypep'] = protein_cyberT_bypep['pVal'].values
    out['ttest_bypep'] = protein_ttest_bypep['pVal'].values
    out.reset_index(inplace=True)
    return out


def do_stat_tests_phospho(ctrl, exp, protein_labels, protein_df):
    """ Do statistical tests for phosphoproteins

    protein_df should have columns 'protein_id', 'fold_change', 'std_err'
        If 'std_err' not provided, will default to 0
        Optional: 'fold_change_nwls'
    """
    if 'std_err'in protein_df:
        std_err = protein_df['std_err'][protein_labels].values
    else:
        std_err = 0
        print "Warning! No standard errors were included in protein_df for phosphoprotein statistical tests, defaulting to zero"

    no_norm = cyberT(ctrl, exp)
    norm_exp = exp.subtract(
            protein_df['fold_change'][protein_labels].reset_index(drop=True),
            axis=0)
    mean_norm = cyberT(ctrl, norm_exp)
    norm_with_err = cyberT(ctrl, norm_exp, base_vars=std_err)

    res = pd.DataFrame({
        'no_adj': no_norm['pVal'],
        'mean_adj': mean_norm['pVal'],
        'var_adj': norm_with_err['pVal'],
        })

    if 'fold_change_nwls' in protein_df:
        # TODO what about non-wls mean???
        norm_non_wls_exp = exp.subtract(
                protein_df['fold_change_nwls'][protein_labels].reset_index(drop=True),
                axis=0)
        mean_norm_nwls = cyberT(ctrl, norm_non_wls_exp)
        res['norm_nwls_adj'] = mean_norm_nwls['pVal']

    return res


## TODO this is VERY JANKY
## Convenience function for protein pval labels
def protein_pval_labels():
    ctrl, exp, is_changed, protein = sample_proteins(500, 0, 2, 2)

    tmp = do_stat_tests_protein(ctrl, exp, protein)
    return tmp.columns

## Convenience function for peptide pval labels
def peptide_pval_labels(run_modT_2sample=True):
    ctrl, exp, is_changed = sample_no_ctrl_uniform(500, 0, 2)

    tmp = do_stat_tests(ctrl, exp, run_modT_2sample)
    return tmp.columns

def phospho_pval_labels(nowls=True):
    ctrl, exp, is_changed, protein = sample_proteins(500, 0, 2, 2)
    protein_df = pd.DataFrame({'protein_id': np.arange(500),
                               'fold_change': np.ones(500)})
    ctrl_p, exp_p, is_changed, mapping = sample_phospho(500, 0, 0, protein_df)

    tmp = do_stat_tests_phospho(ctrl_p, exp_p, mapping, protein_df)
    if nowls:
        return list(tmp.columns) + [u'mean_adj_nowls']
    else:
        return tmp.columns
