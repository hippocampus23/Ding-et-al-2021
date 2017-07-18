"""
Scratchpad for random functions and fragments which are useful

Most of these functions deal with setting up and running long running tasks
such as multiple rounds of simulations
"""
from contextlib import contextmanager
from math import pi
import sys, os, time
sys.dont_write_bytecode = True  # Avoid caching problems

from format_results import *
from heavy_tailed import *
from roc import *
from simulate import *


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout


DEF_FOLD_CHANGES = np.arange(0.1, 1.1, 0.1)
STD = 0.3039  ## Std such that 10% of peptides have |log2(FC)| >= 0.5
THRESH = 0.5  ## Set all FC st |log2(FC)| < THRESH to 0

### Functions for running multiple rounds of trials on peptide level ###

def err_bars_peptide_saveall(fold_changes, num_to_change, background="G", n_runs=500, **kwargs):
    """Test function for determining size of dataframe if all results are saved"""
    start = time.time()

    # Hardcoded params
    N_PEPS = 10000

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid background specification")

    res = np.zeros((n_runs, len(DEFAULT_LABELS), N_PEPS), dtype=np.float32)
    is_changed_all = np.zeros((n_runs, N_PEPS), dtype=int)

    for i in xrange(n_runs):
        if i % 25 == 0:
            print "At iteration %d" % i
        try:
            with suppress_stdout():
                ctrl, exp, is_changed = sampler(
                    N_PEPS,
                    num_to_change,
                    fold_changes,
                    **kwargs)
                p_vals = np.array(do_stat_tests(ctrl, exp))
                res[i,:,:] = p_vals
                is_changed_all[i,:] = is_changed

        except rpy2.rinterface.RRuntimeError:
            print "R error!"
            res[i,:,:] = np.nan
            is_changed_all[i,:] = is_changed

    end = time.time()
    print end - start

    return res, is_changed_all



def err_bars_peptide(fold_changes, num_to_change, background = "U", n_runs=500, labels=DEFAULT_LABELS_MODT_2SAMP, run_modT_2samp = True, **kwargs):
    """ Runs multiple rounds of simulations using given sampler 
        Summarizes overall ROC scores
   
    Args:
        fold_changes: non-negative float or list of floats
        num_to_change: int or list of ints
        background: "U" or "G" for uniform or inverse gamma variance
        filename: optional, if present will save results in binary format
        kwargs: passed to variance generating function
            [Possible args include var, nctrl, nexp, use_var, alpha, beta]

    Returns:
        np array of size N_RUNS x len(DEFAULT_LABELS) x 3
        contains AUC, pAUC, and PRC for each run and metric
        arr[i][j] = (AUC, PRC, pAUC)
    """

    # TODO REMOVE ME
    start = time.time()
   
    # Hardcoded params
    N_PEPS = 10000

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid background specification")

    res = np.zeros((n_runs, len(labels), 3), dtype=np.float32)

    for i in xrange(n_runs):
        if i % 50 == 0:
            print "At iteration %d" % i
        try:
            with suppress_stdout():
                ctrl, exp, is_changed = sampler(
                    N_PEPS,
                    num_to_change,
                    fold_changes,
                    **kwargs)
                p_vals = do_stat_tests(ctrl, exp, run_modT_2samp)
                res[i,:,0], res[i,:,1], res[i,:,2] = roc_prc_scores(
                    is_changed, p_vals, fdr=0.05)
        except rpy2.rinterface.RRuntimeError:
            print "R error!"
            res[i,:,:] = np.nan

    end = time.time()
    print end - start

    return res


def err_bars_fold_change(fold_changes, num_to_change, background = "U", breaks=None, n_runs=500, threshold=0, n_peps = 10000, **kwargs):
    """ Runs multiple rounds of simulation using the given sampler
        ROC scores are calculated for each fold change subset of the original
        Summarizes ROC scores split by fold change
   
    Args:
        fold_changes: non-negative float or list of floats
            OR a function which takes one argument, the number of iterations
            and returns a vector of fold changes
            If fold_change is a function, then breaks MUST be specified
        num_to_change: int or list of ints
        background: "U" or "G" for uniform or inverse gamma variance
        breaks: optional, default None
            How to group fold_changes in summary.
            Array of breakpoints: fold changes will be rounded down to nearest breakpoint.
            MUST be sorted in ascending order
            FCs outside of breakpoints will not be included in ROC calculation
        n_runs: how many simulations to run, default 500
        threshold: Treat all fc < threshold to be control
        kwargs: passed to variance generating function
            [Possible args include var, nctrl, nexp, use_var, alpha, beta]

    Returns:
        np array of size n_runs x len(fold_changes) x len(DEFAULT_LABELS) x 3
        contains AUC, pAUC, and PRC for each run and metric
        arr[i][j][k] = (AUC, PRC, pAUC)
            i = run index
            j = fold_change index
            k = metric index (i.e. cyberT, modT, etc)
    """

    # TODO REMOVE ME
    start = time.time()

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid background specification")
   
    if breaks is None:
        unique_fcs = np.unique(fold_changes)
    else:
        unique_fcs = breaks

    res = np.zeros(
            (n_runs, 
            len(unique_fcs) if breaks is None else len(unique_fcs)-1 ,
            len(DEFAULT_LABELS_MODT_2SAMP), 
            3), dtype=np.float32)

    call_fcs = callable(fold_changes)
    gen_fcs = (fold_changes(i) for i in xrange(n_runs))

    for i in xrange(n_runs):
        if i % 50 == 0:
            print "At iteration %d" % i
        try:
            n_fold_changes = next(gen_fcs) if call_fcs else fold_changes
            with suppress_stdout():
                ctrl, exp, fcs = sampler(
                    n_peps,
                    num_to_change,
                    n_fold_changes,
                    binary_labels=False,
                    **kwargs)
                p_vals = do_stat_tests(ctrl, exp)

            # Create binary labels from scalar
            is_changed = (np.abs(fcs) >= threshold).astype(int)
            # Cast p_vals as 2D array for easier subsetting
            p_vals = np.array(p_vals)

            # Calculate ROC score for EACH fc seperately
            for (j, fc) in enumerate(unique_fcs[:-1]):
                # Subset of peptides which have fcs = 0 or in the bucket
                if breaks is None:
                    idx = np.logical_not(is_changed) | (fcs == fc)  ## Can use equality
                else:
                    idx = np.logical_not(is_changed) | (
                            (fcs >= fc) & (fcs < unique_fcs[j+1]))

                if np.sum(is_changed[idx]) == 0:
                    res[i,j,:,:] = 0  # No positive examples
                else:
                    res[i,j,:,0], res[i,j,:,1], res[i,j,:,2] = roc_prc_scores(
                            is_changed[idx], p_vals[:,idx], fdr=0.05)

        except rpy2.rinterface.RRuntimeError:
            print "R error!"
            res[i,:,:,:] = np.nan

    end = time.time()
    print end - start

    return res

TIME_FORMAT = "%Y-%m-%d_%H:%M"


def simulate_compare_one_two_sided():
    """ Compare one-sided fold change with two-sided fold change

    Compares single fc = 2**(0.5) with
             single fc = 2**(-0.5) with
             two fc = 2**([0.5, -0.5])
             for both uniform and inv gam variance models
    """
    start = time.strftime(TIME_FORMAT)
    fc1 = 0.5
    fc2 = -0.5

    res = {}
    res['hi_1_tail_uni'] = err_bars_peptide(fc1, 1000, "U", n_runs=500)
    res['hi_1_tail_gam'] = err_bars_peptide(fc1, 1000, "G", n_runs=500)
    res['lo_1_tail_uni'] = err_bars_peptide(fc2, 1000, "U", n_runs=500)
    res['lo_1_tail_gam'] = err_bars_peptide(fc2, 1000, "G", n_runs=500)
    res['2_tail_uni'] = err_bars_peptide([fc1,fc2], 500, "U", n_runs=500)
    res['2_tail_gam'] = err_bars_peptide([fc1,fc2], 500, "G", n_runs=500)

    np.save("tmp_%s.npy" % start, res) 
    return res

def simulate_multiple_fc_normal(background="G", n_runs=250, filename=None):
    """ Generate partial ROCs for normally distruited FCs (bucketed!)
        Compare to uniformly distributed FCs

        Buckets into buckets of size 0.1, truncates fold changes to range
        -1 . . . 1
    """
    start = time.strftime(TIME_FORMAT)
    N_TO_CHANGE = 5000
    if filename is None:
        filename = "tmp_%s.npy" % start
    filename_csv = filename[:-4] + ".csv"

    # Set up buckets and weight appropriately
    F_MIN, F_MAX = -8, 8
    norm = lambda x: np.exp(-x**2./(2*STD**2))/(2*pi*STD)**0.5  # Normal pdf
    buckets_all = np.setdiff1d(np.arange(F_MIN, F_MAX+1, 1), [0]) / 10.
    buckets = np.setdiff1d(np.arange(F_MIN, F_MAX+1, 1), np.arange(-4, 5, 1)) / 10.
    n_density_all = norm(buckets_all)
    n_density_all /= sum(n_density_all)
    n_density = norm(buckets)
    n_density /= sum(n_density)

    # This is for fixed bucketed fold changes
    """
    breaks, breaks_all = None, None
    fold_changes = np.repeat(
            buckets_all,
            int(N_TO_CHANGE / len(buckets_all)))
    fold_changes_n = np.repeat(
            buckets,
            np.round(N_TO_CHANGE * n_density).astype(int))
    fold_changes_n_all = np.repeat(
            buckets_all,
            np.round(N_TO_CHANGE * n_density_all).astype(int))
    """

    # This is for continuous fold changes
    breaks = buckets
    breaks_all = buckets_all
    fold_changes = np.arange(F_MIN/10., F_MAX/10., (F_MAX - F_MIN)/(10.*N_TO_CHANGE))
    ## Fixed normal FC
    # fold_changes_n_all = np.random.normal(scale=STD, size=N_TO_CHANGE)
    # fold_changes_n = fold_changes_n_all[np.abs(fold_changes_n_all) >= 0.5]
    ## Randomly generated normal FC at every iteration
    fold_changes_n_all = lambda x: np.random.normal(scale=STD, size=N_TO_CHANGE)
    def fold_changes_n(i):
        tmp = np.random.normal(scale=STD, size=N_TO_CHANGE)
        return tmp[np.abs(tmp) >= 0.5]

    res_uniform = err_bars_fold_change(fold_changes, 1, background,
            threshold=0.2, breaks=breaks_all, n_runs=n_runs)
    res_norm = err_bars_fold_change(fold_changes_n, 1, background,
            threshold=0.5, breaks=breaks, n_runs=n_runs)
    res_norm_all = err_bars_fold_change(fold_changes_n_all, 1, background,
            threshold=0.2, breaks=breaks_all, n_runs=n_runs)

    # Drop the last two labels (t-test and one-tailed t-test)
    keep_pvals = [0,1,4,5]
    res_uniform = res_uniform[:,:,keep_pvals,:]
    res_norm = res_norm[:,:,keep_pvals,:]
    res_norm_all = res_norm_all[:,:,keep_pvals,:]

    # Rearrange matrix so conditions are grouped by fold change 
    # First pad the incomplete bucket list
    missing = np.setdiff1d(buckets_all, buckets)
    idx = np.where(buckets >= missing[0])[0][0] - 1  ## Idx at which to pad
    shape = res_norm.shape
    res_norm_padded = np.concatenate((
            res_norm[:,:idx,:,:],
            np.zeros((shape[0], len(missing), shape[2], shape[3]), dtype=float),
            res_norm[:,idx:,:,:]),
            axis=1)
    # Concat along second to last dimension to squash conditions into setting
    res_arr = np.concatenate(
            (res_uniform, res_norm_padded, res_norm_all),
            axis=2)
    ## For discrete
    # res = {fc: res_arr[:,i,:,:] for i,fc in enumerate(buckets_all)}
    ## For continuous
    res = {"%.1f<=fc<%.1fd" % (buckets_all[i], buckets_all[i+1]): res_arr[:,i,:,:] 
            for i in xrange(len(buckets_all)-1)}
    labels = ["%s: %s" % (lbl, ms)
              for lbl in ["Uniform, all", "Normal", "Normal, all"]
              for ms in ["ModT", "CyberT", "fold_change"]]
    print labels
    
    # Here we save as both csv and npy file
    np.save(filename, res) 
    write_result_dict_to_df(res, labels, filename_csv)  ## TODO make less jank
    return res
    

def simulate_multiple_fc(background="G", n_runs=200, filename=None):
    """ Generate partial ROCs for uniformly distributed FCs
    """
    start = time.strftime(TIME_FORMAT)
    if filename is None:
        filename = "tmp_%s.npy" % start
    fold_changes = np.arange(0, 1, 1./100) + 0.01

    res = err_bars_fold_change(fold_changes, 100, background, n_runs=200)
    np.save(filename, res) 
    return res


def simulate_random_fc():
    """ Creates simulated datasets with random fc, normally distributed
    """
    start = time.strftime(TIME_FORMAT)
    fc = np.random.normal(1.5, 0.15, 1000)
    
    res = err_bars_peptide(fc, 1, "G")
    np.save("tmp_%s.npy" % start, res)
    return res


def simulate_fold_change_range(fold_changes = DEF_FOLD_CHANGES, **kwargs):
    """Creates simulated datasets with error bars
    """
    start = time.strftime(TIME_FORMAT)
    res = {}
    
    for f in fold_changes:
        res[f] = err_bars_peptide(f, 1000, "G", **kwargs)
        np.save("tmp_%s.npy" % start, res)
    return res


def simulate_number_experiments():
    """Runs simulated dataset with different number of samples
    """
    start = time.strftime(TIME_FORMAT)
    res = {}
    f = 0.5
    
    for n in xrange(2, 11):
        res[n] = err_bars_peptide(f, 1000, "G", nexp=n, nctrl=n)
        np.save("tmp_%s.npy" % start, res)
    return res


def simulate_number_channels_imbalanced(n_runs=250, filename=None):
    """Compared balanced and imbalanced number of channels
    """
    start = time.strftime(TIME_FORMAT)
    N_TO_CHANGE = 1000
    if filename is None:
        filename = "tmp_%s.npy" % start
    res = {}
    f = np.random.normal(scale=STD, size=2000)

    trials = [(3,3),
              (2,4),
              (1,5)]
    trials = [(5,5),  ## For n=10
              (4,6),
              (3,7),
              (2,8),
              (1,9)]
    for nctrl, nexp in trials:
        key = "nctrl=%d, nexp=%d" % (nctrl, nexp)
        res[key] = err_bars_peptide(f, 1, "G", nexp=nexp, nctrl=nctrl, n_runs=n_runs)
        np.save(filename, res)
    return res


def simulate_variance_range():                                                  
    u_std = [0.02, 0.06, 0.18]                                                  
    g_beta = [0.05, 0.1, 0.2]                                                   
    fc = 0.5                                                                 
                                                                                
    start = time.strftime(TIME_FORMAT)
    res = {}                                                                    
                                                                                
    for v in u_std: 
        res["uniform_%.2f" % v] = err_bars_peptide(fc, 1000, "U", var=v)         
        res["uniform_lap_%.2f" % v] = err_bars_peptide(
                fc, 1000, "U", var=v*(0.5**0.5), use_var=np.random.laplace)         
        np.save("tmp_%s.npy" % start, res)
    for b in g_beta:                                                            
        res["inv_gamma_%.2f" % b] = err_bars_peptide(fc, 1000, "G", beta=b)      
        res["inv_gamma_lap_%.2f" % v] = err_bars_peptide(
                fc, 1000, "G", var=v*(0.5**0.5), use_var=np.random.laplace)         
        np.save("tmp_%s.npy" % start, res)
                                                                                
    return res  

def simulate_with_gamma(n, alpha=3, beta=0.1, nctrl=3, use_var=np.random.normal):
    variances = 1 / np.random.gamma(alpha, 1./beta, n)
    noise = np.array([
        use_var(loc=0, scale=v**0.5, size=nctrl)
        for v in variances
    ])

    return noise

### Functions for simulating multiple trials on protein level ###

def err_bars_protein(m, num_to_change, fold_changes, peps_per_prot, n_runs=500,**kwargs):
    """ Runs multiple rounds of simulations using given sampler 
        Summarizes overall ROC scores
   
    Args:
        m: number of proteins
        num_to_change: int or list of ints
        fold_changes: non-negative float or list of floats
        background: "U" or "G" for uniform or inverse gamma variance
        filename: optional, if present will save results in binary format
        kwargs: passed to variance generating function
            [Possible args include var, nctrl, nexp, use_var, alpha, beta, background]

    Returns:
        np array of size N_RUNS x len(protein_pval_labels()) x 3
        contains AUC, pAUC, and PRC for each run and metric
        arr[i][j] = (AUC, PRC, pAUC)
    """

    # TODO REMOVE ME
    start = time.time()
    # TODO this is REALLY JANKY
    with suppress_stdout():
       labels = protein_pval_labels()
    # labels = protein_pval_labels()
   
    res = np.zeros((n_runs, len(labels) - 1, 3), dtype=np.float32)

    for i in xrange(n_runs):
        if i % 50 == 0:
            print "At iteration %d" % i
        try:
            with suppress_stdout():
                ctrl, exp, is_changed, protein = sample_proteins(
                    m,
                    num_to_change,
                    fold_changes,
                    peps_per_prot,
                    **kwargs)
                p_vals = do_stat_tests_protein(ctrl, exp, protein)
                # Format is_changed
                is_changed_final = extract_y_act_protein(
                        p_vals, protein, is_changed)
                # Now drop acc_num columns
                p_vals.drop('accession_number', axis=1, inplace=True)
                # Invert fold change columns
                for c in p_vals.columns:
                    if "fold_change" in c:
                        p_vals[c] = np.max(p_vals[c]) - p_vals[c] + 0.01

                res[i,:,0], res[i,:,1], res[i,:,2] = roc_prc_scores(
                    is_changed_final, [p_vals[c] for c in p_vals], fdr=0.05)
        except (rpy2.rinterface.RRuntimeError, ValueError) as e:
            # TODO determine why value errors appear
            # Infinite or NaN input? but where? p-val should be finite . . .
            print "R error"
            print e
            res[i,:,:] = np.nan

    end = time.time()
    print end - start

    return res


def simulate_protein_fold_change_range(
        fold_changes=np.arange(0.05, 0.55, 0.1),
        n_runs=150,
        **kwargs):
    """Creates simulated datasets with error bars
    """
    start = time.strftime(TIME_FORMAT)
    res = {}
    
    for f in fold_changes:
        res[f] = err_bars_protein(5000, 500, f, 2, n_runs=n_runs, **kwargs)
        np.save("tmp_protein_fc_range_%s.npy" % start, res)
    return res


def simulate_protein_num_peps(**kwargs):
    """Creates simulated datasets with error bars
    """
    start = time.strftime(TIME_FORMAT)
    res = {}
    num_peps = [1,2,4,10]  # TODO CHANGE ME 
    N_PEPS = 10000
    
    for n_p in num_peps:
        m = N_PEPS / n_p
        tc = m / 10
        res["u_%d" % n_p] = err_bars_protein(m, tc, 2**(0.3), n_p, **kwargs)
        res["g_%d" % n_p] = err_bars_protein(m, tc, 2**(0.3), n_p, background="G", **kwargs)
        np.save("tmp_protein_num_peps_%s.npy" % start, res)
    return res


def simulate_compare_with_without_central_fc_peptides(
        background="G", n_runs=500, filename=None):
    """ Compares normally distributed fold changes with and without
        including central tendency
    """
    start = time.time()  # Configure filename for results
    if filename is None:
        filename = "tmp_%s.npy" % start
   
    # Hardcoded params
    N_PEPS = 10000
    NUM_FC = 1000

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid background specification")

    res = np.zeros(  # Results array including center
            (n_runs, 
            len(DEFAULT_LABELS), 
            3), dtype=np.float32)
    res_t = res.copy()   # Results array not including center
    res_u = res.copy()   # Results array uniform fc

    for i in xrange(n_runs):
        if i % 50 == 0:
            print "At iteration %d" % i
        try:
            # Generate fold changes from normal distribution
            fold_changes = np.random.normal(0, STD, NUM_FC)
            fold_changes_t = fold_changes[np.abs(fold_changes) >= THRESH]
            fold_changes_u = np.arange(-1., 1., 2. / NUM_FC)
            with suppress_stdout():
            #     # Pvals including all FCs
            #     ctrl, exp, fcs = sampler(
            #         N_PEPS,
            #         1,
            #         fold_changes,
            #         binary_labels=False)
            #     p_vals = do_stat_tests(ctrl, exp)
            #     # Pvals when all central FCs are set to zero
            #     ctrl_t, exp_t, is_changed_t = sampler(
            #         N_PEPS,
            #         1,
            #         fold_changes_t,
            #         binary_labels=True)
            #     p_vals_t = do_stat_tests(ctrl_t, exp_t)
                  # Pvals with uniform fcs
                ctrl_u, exp_u, is_changed_u = sampler(
                    N_PEPS,
                    1,
                    fold_changes_u,
                    binary_labels=True)
                p_vals_u = do_stat_tests(ctrl_u, exp_u)

            # is_changed = (np.abs(fcs) >= THRESH).astype(int)
            # res[i,:,0], res[i,:,1], res[i,:,2] = roc_prc_scores(
            #     is_changed, p_vals, fdr=0.05)
            # res_t[i,:,0], res_t[i,:,1], res_t[i,:,2] = roc_prc_scores(
            #     is_changed_t, p_vals_t, fdr=0.05)
            res_u[i,:,0], res_u[i,:,1], res_u[i,:,2] = roc_prc_scores(
                is_changed_u, p_vals_u, fdr=0.05)

        except rpy2.rinterface.RRuntimeError:
            print "R error!"
            res[i,:,:] = np.nan

    # final = {'include_all_fcs': res,
    #          'threshold_fcs': res_t}
    final = {'uniform': res_u}
    np.save(filename, final)
    end = time.time()
    print end - start

    return final


def DEP_simulate_fold_change_range(
        fold_changes = DEF_FOLD_CHANGES,
        var = 0.06,
        nctrl=3,
        nexp=3):
    """Creates simulated datasets for multiple fold changes

    Uses a SINGLE fold change per iteration
    """
    res = {}
    
    for f in fold_changes:
        ctrl, exp, is_changed = sample_no_ctrl_uniform(
                10000,
                1000,
                f,
                var,
                nctrl=nctrl,
                nexp=nexp)
        res[f] = (
            do_stat_tests(ctrl, exp),
            is_changed
        )
    return res


def DEP_simulate_fold_change_vars(
        fold_changes = DEF_FOLD_CHANGES,
        variances = [0.01, 0.05, 0.1, 0.2],
        nctrl=3,
        nexp=3):
    """Create simulated datasets for fold changes AND variances"""

    res = {}
    
    for f in fold_changes:
        for v in variances:
            ctrl, exp, is_changed = sample_no_ctrl_uniform(
                    10000,
                    1000,
                    f,
                    v,
                    nctrl=nctrl,
                    nexp=nexp)
            res[(f,v)] = (
                do_stat_tests(ctrl, exp),
                is_changed
            )

    return res
