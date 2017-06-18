"""
Scratchpad for random functions and fragments which are useful
in the IPython notebook

IMPORTANT
Some imports are omitted because functions are ONLY for interactive shell
Run from IPython with -i flag
"""
from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout


DEF_FOLD_CHANGES = [2**i for i in np.arange(0.1, 1.1, 0.1)]


def err_bars_peptide(fold_changes, num_to_change, background = "U", **kwargs):
    """ Runs multiple rounds of simulations using given sampler 
   
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
    N_RUNS = 2
    N_PEPS = 10000

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid background specification")

    res = np.zeros((N_RUNS, len(DEFAULT_LABELS), 3), dtype=np.float32)

    for i in xrange(N_RUNS):
        if i % 50 == 0:
            print "At iteration %d" % i
        try:
            with suppress_stdout():
                ctrl, exp, is_changed = sampler(
                    N_PEPS,
                    num_to_change,
                    fold_changes,
                    **kwargs)
                p_vals = do_stat_tests(ctrl, exp)
                res[i,:,0], res[i,:,1], res[i,:,2] = roc_prc_scores(
                    is_changed, p_vals, fdr=0.05)
        except rpy2.rinterface.RRuntimeError:
            print "R error!"
            res[i,:,:] = np.nan

    # TODO filename
    # Runs 500 rounds of simulations and finds error bars on pAUC results
    # Sav0es results to binary file if filename is not None
    end = time.time()

    print end - start

    return res


def simulate_fold_change_range(fold_changes = DEF_FOLD_CHANGES, **kwargs):
    """Creates simulated datasets with error bars
    """
    res = {}
    kwargs.setdefault("background", "U")
    
    for f in fold_changes:
        res[f] = err_bars_peptide(f, 1000, **kwargs)
    return res


def simulate_number_experiments():
    """Runs simulated dataset with different number of samples
    """
    res = {}
    f = 2**0.5
    
    for n in xrange(2, 11):
        res[n] = err_bars_peptide(f, 1000, "G", nexp=n, nctrl=n)
    return res


def simulate_variance_range():                                                  
    u_std = [0.02, 0.06, 0.18]                                                  
    g_beta = [0.05, 0.1, 0.2]                                                   
    fc = 2**0.5                                                                 
                                                                                
    res = {}                                                                    
                                                                                
    for v in u_std:                                                             
        res["uniform_%.2f" % v] = err_bars_peptide(fc, 1000, "U", var=v)        
        res["uniform_lap_%.2f" % v] = err_bars_peptide(                         
                fc, 1000, "U", var=v*(0.5**0.5), use_var=np.random.laplace)     
    for b in g_beta:                                                            
        res["inv_gamma_%.2f" % b] = err_bars_peptide(fc, 1000, "G", beta=b)     
        res["inv_gamma_lap_%.2f" % v] = err_bars_peptide(                       
                fc, 1000, "G", beta=b, use_var=np.random.laplace) 
                                                                                
    return res  

def simulate_with_gamma(n, alpha=3, beta=0.1, nctrl=3, use_var=np.random.normal):
    variances = 1 / np.random.gamma(alpha, 1./beta, n)
    noise = np.array([
        use_var(loc=0, scale=v**0.5, size=nctrl)
        for v in variances
    ])

    return noise


### DEPRECATED ###
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
