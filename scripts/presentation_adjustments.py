import sys
sys.dont_write_bytecode = True  # Avoid caching problems

import matplotlib
import matplotlib.pyplot as plt

from constants import COLORS
from format_results import *
from roc import *
from sample import *


def plot_example_roc_curves():
    # Set up font
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 20}
    matplotlib.rc('font', **font)
    """ Generate panel 2 of figure 1 """
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


def pvalue_multipanel():
    """ Generate panel comparing p-value distributions 
    Compare uniform and inverse gamma """
    f, axarr = plt.subplots(2, 5, sharex='col', sharey='row')

    ctrl_u, exp_u, _ = sample_no_ctrl_gamma(10000, 0, 0)
    pvals_u = do_stat_tests(ctrl_u, exp_u, True)
    
    ctrl_g, exp_g, _ = sample_no_ctrl_uniform(10000, 0, 0)
    pvals_g = do_stat_tests(ctrl_g, exp_g, True)

    plot_pvalue_dist(pvals_u, axarr[0])
    plot_pvalue_dist(pvals_g, axarr[1])

    for ax in axarr[0]:
        ax.set_xlabel('')
    for ax in axarr[1]:
        ax.set_title('')

    return f


def pvalue_multipanel_noise():
    """ Compare p-value distributions different noise distributions """ 
    f, axarr = plt.subplots(3, 5, sharex='col', sharey='row')
    
    DF = 3
    def t_dist(loc, scale, size=1):
        return np.random.standard_t(DF, size=size)*scale
   
    ctrl_n, exp_n, _ = sample_no_ctrl_gamma(10000, 0, 0, use_var=np.random.normal)
    plot_pvalue_dist(do_stat_tests(ctrl_n, exp_n, True), axarr[0])
    ctrl_l, exp_l, _ = sample_no_ctrl_gamma(10000, 0, 0, use_var=np.random.laplace)
    plot_pvalue_dist(do_stat_tests(ctrl_l, exp_l, True), axarr[1])
    ctrl_t, exp_t, _ = sample_no_ctrl_gamma(10000, 0, 0, use_var=t_dist)
    plot_pvalue_dist(do_stat_tests(ctrl_t, exp_t, True), axarr[2])

    for ax in axarr[0]:
        ax.set_xlabel('')
    for ax in axarr[1]:
        ax.set_xlabel('')
        ax.set_title('')
    for ax in axarr[2]:
        ax.set_title('')

    return f

def volcano_multipanel(background="U"):
    """ Generate panel comparing volcano plots
    Compare uniform and inverse gamma
    """
    f, axarr = plt.subplots(2, 5, sharex='col', sharey='row')

    if background == "U":
        sampler = sample_no_ctrl_uniform
    elif background == "G":
        sampler = sample_no_ctrl_gamma
    else:
        raise ValueError("Invalid specification for background")

    ctrl_u, exp_u, is_changed_u = sampler(10000, 1000, 0.5)
    pvals_u = do_stat_tests(ctrl_u, exp_u, True)
    
    pvals_c = pd.DataFrame.from_items([
        (col, multipletests(pvals_u[col], 0.05, method='fdr_bh')[1])
        for col in pvals_u.columns
    ])

    volcano_plots(pvals_u, ctrl_u, exp_u, is_changed_u, axarr[0])
    volcano_plots(pvals_c, ctrl_u, exp_u, is_changed_u, axarr[1])

    # Remove unnecessary labels
    for ax in axarr[0]:
        ax.set_xlabel('')
    for ax in axarr[1]:
        ax.set_title('')
    axarr[1][0].set_ylabel('$-\log_{10}$(Adjusted p-value)')

    return f

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
    true_setting = [float(x.split('<')[0]) for x in setting]
    pass

def regenerate_dataframes():
    """
    Regenerate dataframes
    """
    fc_range_gam = np.load('peptide_fc_range_gam_FINAL.npy')[()]
    fc_range_uni = np.load('peptide_fc_range_uni_FINAL.npy')[()]
    nexp = np.load('peptide_nexp_modtfix_FINAL.npy')[()]
    nexp_imba = np.load('peptide_nexp_imba_modtfix_FINAL.npy')[()]
    var = np.load('peptide_variances_FINAL.npy')[()]
    ds_size = np.load('peptide_ds_size_FINAL.npy')[()]

    # fix nexp_imba keys
    nexp_imba = {(('(%d,%d)' % k) if k != '_labels' else k):v for
            k,v in nexp_imba.iteritems()}
    # Fix ds_size keys
    ds_size = {(('%d: %d' % k) if k != '_labels' else k):v for 
            k,v in ds_size.iteritems()}
    
    write_result_dict_to_df(fc_range_gam, None).to_csv(
            'df_peptide_fc_range_gam_FINAL.csv')
    write_result_dict_to_df(fc_range_uni, None).to_csv(
            'df_peptide_fc_range_uni_FINAL.csv')
    write_result_dict_to_df(nexp, None).to_csv(
            'df_peptide_nexp_modtfix_FINAL.csv')
    write_result_dict_to_df(nexp_imba, None).to_csv(
            'df_peptide_nexp_imba_modtfix_FINAL.csv')
    write_result_dict_to_df(var, None).to_csv(
            'df_peptide_variances_FINAL.csv')
    write_result_dict_to_df(ds_size, None).to_csv(
            'df_peptide_dataset_size_FINAL.csv')



