import time
import numpy as np
from stat_tests import do_all_tests, TESTS
from sampler import sample
from constants import N_RUNS
from eval_results import roc_prc_scores, power_analysis

def simulate_fold_change_range(fold_changes, var_type):
    """
    creates simulated data sets and analyses it for different fold changes

    :param fold_changes:   list of fold changes
    :param var_type:       either "uniform", "gamma" or "trend". describes how variance is sampled
    :return:               dictionary of containing a numpy.ndarray (N_RUNS x number of tests x 3),
                           with arr[i][j] = (AUC, PRC, pAUC) for each fold change
    """
    return {fc : err_bars_peptide_roc(var_type, fold_change=fc) for fc in fold_changes}


def simulate_number_experiments(num_ctrls, var_type):
    """
    creates simulated dataset with different number of samples
    control and experiment channel numbers are balanced

    :param n_ctrl:    list of number of control channels (i.e. number of experiment channels)
    :param var_type:  either "uniform", "gamma" or "trend". describes how variance is sampled
    :return:          dictionary of containing a numpy.ndarray (N_RUNS x number of tests x 3),
                           with arr[i][j] = (AUC, PRC, pAUC) for each number of channels
    """
    return {n : err_bars_peptide_roc(var_type, num_ctrl=n, num_exp=n) for n in num_ctrls}



def simulate_number_channels_imbalanced(trials, var_type):
    """
    Compared balanced and imbalanced number of channels

    :param trials:    list of pairs (nctrl, nexp)
    :param var_type:  either "uniform", "gamma" or "trend". describes how variance is sampled
    :return dictionary of containing a numpy.ndarray (N_RUNS x number of tests x 3),
                           with arr[i][j] = (AUC, PRC, pAUC) for each pair
    """

    res = {}
    for nctrl, nexp in trials:
        res[str((nctrl, nexp))] = err_bars_peptide_roc(var_type, num_ctrl=nctrl, num_exp=nexp)
    return res


def simulate_variance_range(vars, betas):
    """
    creates simulated dataset with different variance (uniform, gamma different params)
    and noise distributions (normal, laplace, t)

    :param vars:   parameters for uniform distribution
    :param betas:  parameters for inverse gamma distribution
    :return:       dictionary containing a numpy.ndarray (N_RUNS x number of tests x 3) for each option
    """

    DF = 3

    def t_dist(loc, scale, size=1):
        return np.random.standard_t(DF, size=size) * scale + loc

    res = {}
    for v in vars:
        res["uniform_norm_%.2f" % v] = err_bars_peptide_roc("uniform", pep_var=v)
        res["uniform_lap_%.2f" % v] = err_bars_peptide_roc("uniform", pep_var=v, use_var=np.random.laplace)
        res["uniform_t_%.2f" % v] = err_bars_peptide_roc("uniform", pep_var=v, use_var=t_dist)
    for b in betas:
        res["invgam_norm_%.2f" % b] = err_bars_peptide_roc("gamma", beta=b)
        res["invgam_lap_%.2f" % b] = err_bars_peptide_roc("gamma", beta=b, use_var=np.random.laplace)
        res["invgam_t_%.2f" % b] = err_bars_peptide_roc("gamma", beta=b, use_var=t_dist)

    return res


def simulate_fdr_fc_range(fold_changes, var_type):
    """
    creates simulated data sets and analyses it for different fold changes

    :param fold_changes:   list of fold changes
    :param var_type:       either "uniform", "gamma" or "trend". describes how variance is sampled
    :return:               dictionary of containing a numpy.ndarray (N_RUNS x number of tests x 3),
                           with arr[i][j] = (AUC, PRC, pAUC) for each fold change
    """
    return {fc: err_bars_peptide_fdr(var_type, fold_change=fc) for fc in fold_changes}


def simulate_size_dataset(num_peps, perc_changed, var_type):
    """
    Vary size of peptide dataset and percentage of peptides with a fold change

    :param num_peps:      number of peptides in the simulated data set
    :param perc_changed:  percentage of peptides that have a fold change
    :param var_type:       either "uniform", "gamma" or "trend". describes how variance is sampled
    :return               dictionary containing a numpy.ndarray (N_RUNS x number of tests x 3) for each option
    """

    res = {}
    for n in num_peps:
        for m in perc_changed:
            res[str(n) + "_" + str(int(n*m))] = err_bars_peptide_roc(var_type, num_changed=int(n*m),  num_peps=n)

    return res


def err_bars_peptide_roc(var_type, **kwargs):
    """
    Runs multiple rounds of simulations
    Summarizes overall ROC scores

    :param var_type:  either "uniform", "gamma" or "trend". describes how variance is sampled
    :param kwargs:    optional arguments to be passed to sampler
    :return:          numpy.ndarray (N_RUNS x number of tests x 3), arr[i][j] = (AUC, PRC, pAUC)
    """

    start = time.time()
    res = np.zeros((N_RUNS, len(TESTS), 3), dtype=np.float32)

    for i in xrange(N_RUNS):
        if i % 50 == 0:
            print "At iteration %d" % i
        ctrl, exp, is_changed = sample(var_type, **kwargs)
        p_vals = do_all_tests(ctrl, exp)
        res[i, :, 0], res[i, :, 1], res[i, :, 2] = roc_prc_scores(is_changed, p_vals)

    end = time.time()
    print end - start

    return res

def err_bars_peptide_fdr(var_type, **kwargs):
    """
    Runs multiple rounds of simulations
    Summarizes overall FDR scores

    :param var_type:  either "uniform", "gamma" or "trend". describes how variance is sampled
    :param kwargs:    optional arguments to be passed to sampler
    :return:          numpy.ndarray (N_RUNS x number of tests x 4), arr[i][j] = (FP_raw, TP_raw, FP_adj, TP_adj)
    """

    start = time.time()
    res = np.zeros((N_RUNS, len(TESTS), 4), dtype=np.float32)

    for i in xrange(N_RUNS):
        if i % 50 == 0:
            print "At iteration %d" % i
        ctrl, exp, is_changed = sample(var_type, **kwargs)
        p_vals = do_all_tests(ctrl, exp)
        res[i, :, 0], res[i, :, 1], res[i, :, 2], res[i, :, 3] = power_analysis(is_changed, p_vals)
        # avoid memory filling up with unused variables
        del ctrl, exp, is_changed

    end = time.time()
    print end - start

    return res
