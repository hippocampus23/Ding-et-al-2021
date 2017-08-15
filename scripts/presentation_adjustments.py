import matplotlib

from roc import *
from sample import *

"""
Color mapping from R for peptides
"""
colors = {'CyberT': "#8DD3C7",
          'Moderated T (1-sample)': "#FFFFB3", 
          'Moderated T (2-sample)': "#BEBADA",
          'Absolute fold change': "#FB8072",              
          't-test (2-sample)': "#80B1D3",
          't-test (1-sample)': "#FDB462",
}

def plot_example_roc_curves():
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
    matplotlib.rc('font', **font)

    ctrl, exp, is_changed = sample_no_ctrl_gamma(10000, 1000, 0.4, nctrl=4, nexp=4)
    pvals = do_stat_tests(ctrl, exp, True)
    plot_both(is_changed, pvals.values.transpose(), list(pvals.columns))


"""
Transform results of variance AUC dictionary to more readable/usable ones
"""
def transform_keys(frm):
    # Replace inv_gamma with invgam
    x1 = ['invgam' + k[9:] if k.startswith('inv_gamma') else k for k in frm]
    
    x2 = []
    # Now add norm in if necessary
    for k in x1:
        toks = k.split('_')
        if len(toks) == 2:
            x2.append('_'.join((toks[0], 'norm', toks[1])))
        else:
            x2.append(k)
    
    return x2
