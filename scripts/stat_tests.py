# Statistical tests for identifying significant fold changes

import pandas as pd
import numpy as np
from collections import OrderedDict
from scipy import stats
from rpy2.robjects import r, pandas2ri

# Setup R environment
pandas2ri.activate()
r["source"]("modT.R")
r["source"]("bayesreg.R")


def regT(ctrl, exp, **kwargs):
    """
    regularized t-test for comparing the peptide intensities of control and experiment data.

    :param ctrl:    pandas.DataFrame where each column represents a control channel
    :param exp:     pandas.DataFrame where each column represents an experiment channel
    :param kwargs:  optional arguments to be passed to bayesT function
    :return:        numpy.ndarray containing the p-values for each peptide
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError("Length of exp and ctrl data frames not identical")
    ctrl.reset_index(drop=True, inplace=True)
    exp.reset_index(drop=True, inplace=True)

    df = pd.concat([ctrl, exp], axis=1, ignore_index=True)
    df.columns = (["C%d" % (i+1) for i in xrange(ctrl.shape[1])] +
            ["E%d" % (i+1) for i in xrange(exp.shape[1])])
    res = r["bayesT"](df, numC = ctrl.shape[1], numE = exp.shape[1], ppde=False, doMulttest=True, **kwargs)
    res = pandas2ri.ri2py(res)["pVal"].values

    return res

def modT(ctrl, exp, run_2sample=False, trend=False, robust=False):
    """
    Moderated t-test for comparing the peptide intensities of control and experiment data.
    NOTE: only defined if ncol(ctrl) == 1 or ncol(ctrl) = ncol(exp)

    :param ctrl:        pandas.DataFrame where each column represents a control channel
    :param exp:         pandas.DataFrame where each column represents an experiment channel
    :param run_2sample: if True, 2-sample moderated t-test is used, else 1-sample
    :param trend:       if True, moderated t-test incorporates a mean-variance trend
    :param robust       if True, moderated t-test method is more robust against outliers (takes much more time!!!)
    :return:            numpy.ndarray containing the p-values for each peptide
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError("Length of exp and ctrl data frames not identical")

    if not run_2sample:
        if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
            data = pd.DataFrame(exp.values - ctrl.values)
        else:
            raise ValueError("Imbalanced exp and ctrl channels, 1-sample test not possible")
        # Also create id col for convenience
        data.columns = ["Log2_Ratio_%d" % i for i in xrange(data.shape[1])]
        data_cols = data.columns
        data["id"] = np.arange(len(data))
        res = r["modT_test"](data, "placeholder", id_col="id", data_col=data_cols, dframe=True, trend=trend,
                             robust=robust)
    else:
        data = pd.DataFrame(np.concatenate((ctrl.values, exp.values), axis=1))
        data.columns = (["C%d" % (i + 1) for i in xrange(ctrl.shape[1])] +
                        ["E%d" % (i + 1) for i in xrange(exp.shape[1])])
        design = np.array(([-1] * ctrl.shape[1]) + ([1] * exp.shape[1]), dtype=int)
        data_cols = data.columns
        data["id"] = np.arange(len(data))
        res = r["modT_test"](data, "placeholder", id_col="id", data_col=data_cols, dframe=True, trend=trend,
                             design=design, robust=robust)

    return pandas2ri.ri2py(res)["P.Value"].values


def t_test(ctrl, exp, ratio_test=False):
    """
    Standard t-test for comparing the peptide intensities of control and experiment data.

    :param ctrl:       pandas.DataFrame where each column represents a control channel
    :param exp:        pandas.DataFrame where each column represents an experiment channel
    :param ratio_test: if True, 1-sample t-test is performed, else 2-sample t-test
    :return:           numpy.ndarray containing the p-values for each peptide
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError("Length of exp and ctrl data frames not identical")

    if ratio_test:
        if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
            data = np.array(exp.values - ctrl.values)
            _, pvals = stats.ttest_1samp(data, 0, axis=1)
        else:
            raise ValueError("Imbalanced exp and ctrl channels, 1-sample test not possible")
    else:
        _, pvals = stats.ttest_ind(ctrl, exp, axis=1)

    return pvals

def rank_by_fold_change(ctrl, exp):
    """
    Rank peptides by absolute fold change between experiment and control.

    :param ctrl:       pandas.DataFrame where each column represents a control channel
    :param exp:        pandas.DataFrame where each column represents an experiment channel
    :return:           numpy.ndarray will be used instead of p-values
    """
    fold_change = np.abs(np.mean(ctrl, axis=1) - np.mean(exp, axis=1))
    return np.max(fold_change) + 0.01 - fold_change


# makes it easier to apply all tests in order
TESTS = OrderedDict([
    ("fold change",  rank_by_fold_change),
    ("t-test-1",     lambda ctrl, exp: t_test(ctrl, exp, True)),
    ("t-test-2",     lambda ctrl, exp: t_test(ctrl, exp, False)),
    ("modT-1",       lambda ctrl, exp: modT(ctrl, exp, False, False)),
    ("modT-2",       lambda ctrl, exp: modT(ctrl, exp, True, False)),
    ("modT-1 trend", lambda ctrl, exp: modT(ctrl, exp, False, True)),
    ("modT-2 trend", lambda ctrl, exp: modT(ctrl, exp, True, True)),
    ("RegT",         regT)
])


def do_all_tests(ctrl, exp):
    """
    Applies all statistical tests listed above
    Note: if the dimensions are unsuitable for a ratio 1-sample t-test,
          1-sample test are not applied and the p-value defaults to NaN

    :param ctrl: pandas.DataFrame where each column represents a control channel
    :param exp:  pandas.DataFrame where each column represents an experiment channel
    :return:     pandas.DataFrame with the p-values for each test

    """

    skip_tests = {}
    if ctrl.shape[1] != 1 and ctrl.shape[1] != exp.shape[1]:
        skip_tests = ["t-test-1", "modT-1", "modT-1 trend"]

    return pd.DataFrame(OrderedDict([(label, np.nan if label in skip_tests else test(ctrl, exp))
                                     for label, test in TESTS.iteritems()]))
