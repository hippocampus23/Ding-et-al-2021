# Functions for simulating proteomics data

import numpy as np
import pandas as pd
import math
from constants import NUM_CHANGED, NUM_PEPTIDES, NUM_CHANNELS, VAR_COEFFICIENTS, MEAN_PEPTIDE_INTENSITY, \
                      PEPTIDE_VAR, ALPHA, BETA, OVERALL_VAR, FOLD_CHANGE


def sample(pep_var_type, num_peps = NUM_PEPTIDES, num_ctrl=NUM_CHANNELS/2, num_exp=NUM_CHANNELS/2,
           num_changed = NUM_CHANGED, fold_change=FOLD_CHANGE, use_var=np.random.normal, alpha=ALPHA,
           beta=BETA, pep_var=PEPTIDE_VAR):
    """
    Simulate a MS data set with experiment and control data

    :param pep_var_type:     type of distribution to use to sample peptide variance, either "uniform", "gamma" or "trend"
    :param num_peps:         number of peptides (i.e. number of rows of the data set)
    :param num_ctrl:         number of control channels (i.e. columns of data)
    :param num_exp:          number of experiment channels (i.e. columns of data)
    :param num_changed:      number of peptides that get a fold change
    :param fold_change:      fold change to add to perturbed peptides
    :param use_var:          function that takes arguments loc, scale and size for generating differences
                             between peptides
    :param alpha:            parameter for inverse gamma distribution
    :param beta              parameter for inverse gamma distribution
    :param pep_var               uniform variance if pep_var_type is "uniform"
    :return:
                 ctrl:       pandas.DataFrame with log 2 control data
                 exp:        pandas.DataFrame with log 2 experiment data
                 is_changed: numpy.ndarray with 0 for no fold change and 1 for fold change
    """
    # helpers for sampling variance
    def var_inv_gamma_sampler():
        return 1 / np.random.gamma(alpha, 1. / beta)

    def var_trend_sampler(pep_mean):
        if pep_mean < 10:
            return VAR_COEFFICIENTS[0] * var_inv_gamma_sampler()
        elif pep_mean > 18:
            return VAR_COEFFICIENTS[1] * var_inv_gamma_sampler()
        else:
            return VAR_COEFFICIENTS[math.ceil(pep_mean * 2)] * var_inv_gamma_sampler()

    # true peptide intensities in control data
    pep_means = np.random.normal(MEAN_PEPTIDE_INTENSITY, math.sqrt(OVERALL_VAR), num_peps)

    if pep_var_type == "uniform":
        var_sampler = lambda _: pep_var
    elif pep_var_type == "gamma":
        var_sampler = lambda _: var_inv_gamma_sampler()
    elif pep_var_type == "trend":
        var_sampler = lambda pep_mean: var_trend_sampler(pep_mean)
    else:
        raise ValueError("pep_var_type must be uniform, gamma, or trend")

    # generate control data
    var_ctrl = [var_sampler(mean) for mean in pep_means]
    ctrl = np.array([
        use_var(loc=0, scale=var_ctrl[i] ** 0.5, size=num_ctrl) + pep_means[i]
        for i in range(num_peps)
    ])

    # generate experiment data
    # For convenience, always perturb the first num_to_change peptides
    # This should not affect the analysis: can always randomize order
    pep_means[:num_changed] += fold_change
    var_exp = [var_sampler(mean) for mean in pep_means]
    exp = np.array([
        use_var(loc=0, scale=var_exp[i] ** 0.5, size=num_exp) + pep_means[i]
        for i in range(num_peps)
    ])

    is_changed = np.concatenate([np.ones(num_changed, dtype=bool), np.zeros(num_peps-num_changed, dtype=bool)])

    return pd.DataFrame(ctrl), pd.DataFrame(exp), is_changed
