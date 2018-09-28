import time
import numpy as np
import pandas as pd
from stat_tests import do_all_tests, TESTS
from sampler import sample
from constants import N_RUNS, FDR, ALPHA
from eval_results import roc_prc_scores, power_analysis
from sklearn.metrics import roc_curve, auc


def simulate_FC_SD(var_type, test_labels, stdevs, pauroc):
    """
    finds the smallest fold change, given the standard deviation of
    the sample noise, that produces the desired partial AUROC score
    N_RUNS times for each sd provided

    :param var_type:       either "uniform", "gamma" or "trend". describes how variance is sampled
    :param tests:          statistical tests to use
    :param stdevs:         standard deviations to calculate fold changes for, should be in ascending order
    :param pauroc:         desired partial AUROC value
    :return:               pandas.DataFrame containing the fold changes, tests and settings as sd
    """

    res = pd.DataFrame(columns=["FC", "labels", "setting"], index=xrange(N_RUNS * len(test_labels) * len(stdevs)))
    i = 0
    for label in test_labels:
        guess = 0.1
        for sd in stdevs:
            for j in xrange(N_RUNS):
                guess = _find_FC(var_type, sd, pauroc, guess, TESTS[label])
                res.at[i, "FC"] = guess
                res.at[i, "label"] = label
                res.at[i, "setting"] = sd
                i += 1

    return res


def _find_FC(var_type, sd, pauroc_goal, start_fc, test, max_iter=100, precision=0.003, scale=0.05):
    """
    gradient descent based method for finding the smallest fold change, given the standard deviation of
    the sample noise, that produces the desired partial AUROC score

    :param var_type:       either "uniform", "gamma" or "trend". describes how variance is sampled
    :param sd:             mean standard deviation of sample noise
    :param pauroc_goal:    desired partial AUROC value
    :param start_fc:       where to start looking
    :param test:           statistical tests that has the shape test(ctrl, exp) and returns an array of p-values
    :param max_iter:       maximum number of iterations
    :param precision:      how close the pAUROC with the estimated fold change should be to the desired
                           pAUROC
    :param scale:          multiply the step size by
    :return:               approximate fold change that produces the the desired pAUROC score given the SD
    """

    fc = start_fc
    if var_type == "uniform":
        simulator = lambda fold_change: sample(var_type, pep_var=sd**2, fold_change=fold_change)
    else:
        simulator = lambda fold_change: sample(var_type, beta=(ALPHA-1) * sd**2, fold_change=fold_change)


    for x in range(max_iter):
        # simulate data
        ctrl, exp, is_changed = simulator(fc)
        p_vals = test(ctrl, exp)

        # calculate pauroc
        predicted = - np.log(p_vals)
        fpr, tpr, _ = roc_curve(is_changed, predicted)
        try:
            idx = next(i for i, v in enumerate(fpr) if v > FDR)
        except StopIteration:
            idx = len(fpr) - 1

        t_fpr, t_tpr = fpr[:idx + 1], tpr[:idx + 1]
        t_fpr[-1] = FDR

        pauroc = auc(t_fpr, t_tpr) / FDR
        print "pauroc:", pauroc, "fc:", fc

        if (abs(pauroc_goal - pauroc) < precision):
            print "\n\n\n"
            return fc

        # update fold change
        fc += scale * (0.75 - pauroc)





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
