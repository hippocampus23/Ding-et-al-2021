import matplotlib

from roc import *
from sample import *

"""
Color mapping from R for peptides
"""
COLORS = {'CyberT': "#b79f00",
          'Moderated T (1-sample)': "#00ba38",
          'Moderated T (2-sample)': "#00bfc4",
          'Absolute fold change': "#f8766d",
          't-test (2-sample)': "#f564e3",
          't-test (1-sample)': "#619cff",
}

def plot_example_roc_curves():
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 20}
    matplotlib.rc('font', **font)
    colors = [
        COLORS['CyberT'],
        COLORS['Absolute fold change'],
        COLORS['Moderated T (1-sample)'],
        COLORS['t-test (1-sample)'],
        COLORS['t-test (2-sample)'],
        COLORS['Moderated T (2-sample)'],
    ]

    ctrl, exp, is_changed = sample_no_ctrl_gamma(10000, 1000, 0.4, nctrl=4, nexp=4)
    pvals = do_stat_tests(ctrl, exp, True)
    plot_both(is_changed, pvals.values.transpose(), list(pvals.columns), colors=colors)


"""
Transform results of variance AUC dictionary to more readable/usable ones
"""
def transform_keys(frm):
    # Replace inv_gamma with invgam
    x1 = ['invgam' + k[9:] if k.startswith('inv_gamma') else k for k in frm]
    
    x2 = []
    # Now add 'norm' for Gaussian noise
    for k in x1:
        toks = k.split('_')
        if len(toks) == 2:
            x2.append('_'.join((toks[0], 'norm', toks[1])))
        else:
            x2.append(k)
    
    return x2

"""
Format the dataframe result of running continuous fold changes
"""
def format_multiple_fc_cont(df):
    labels = df['labels']
    true_label = pd.Series([x.split(': ')[1] for x in labels])
    model = pd.Series([x.split(': ')[0] for x in labels])
    
    model[model=='Uniform, all'] = 'Uniform'
    model[model=='Normal, all'] = 'Normal'

    # TODO
    # Drop everything between -0.2 and 0.2
     
    setting = df['setting']
    true_setting = [float(x.split('<')[0] for x in setting]
    pass

