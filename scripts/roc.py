from scipy.stats import gaussian_kde
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
from statsmodels.sandbox.stats.multicomp import multipletests

import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import LABEL_MAPPING

def plot_roc(y_act, pred, ax = None, is_pval=True, label='Area ', color=None):
    """
    Returns ax, figure
    Figure is none if ax was provided
    """
    if ax is None:
        f, ax = plt.subplots(1,1)
    else:
        f = None

    if is_pval:
        y_pred = - np.log(pred)
    else:
        y_pred = pred
    fpr, tpr, _ = roc_curve(y_act, y_pred)
    AUC = roc_auc_score(y_act, y_pred)

    lw = 3
    roc_ln, = ax.plot(fpr, tpr, lw=lw, label=label + ' - %.3f' % AUC, color=color)
    diag_ln, = ax.plot([0,1.01], [0,1.01], color='grey', lw=2, linestyle='--')
    
    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC curve')
    ax.legend(loc='lower right', fontsize='medium', title='AUC')
    # return ax, f

    return AUC, fpr, tpr


def plot_prc(y_act, pred, ax = None, is_pval=True, label='Area ', color=None):
    """
    Returns ax, figure
    Figure is none if ax was provided
    """
    if ax is None:
        f, ax = plt.subplots(1,1)
    else:
        f = None

    if is_pval:
        y_pred = - np.log(pred)
    else:
        y_pred = pred
    prec, rec, _ = precision_recall_curve(y_act, y_pred)
    AUC = auc(rec, prec)

    lw = 3
    frac_ones = np.mean(y_act)
    roc_ln, = ax.plot(rec, prec, lw=lw, label=label + ' - %.3f' % AUC, color=color)
    base_ln, = ax.plot([0,1.01], [frac_ones,frac_ones], color='grey', lw=2, linestyle='--')

    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('PRC curve')
    ax.legend(loc='lower left', fontsize= 'medium', title='AUC')
    return AUC, ax, f


def plot_partial_auc(y_act, pred, fdr=0.05, ax=None, is_pval=True, label='Area ', color=None):
    """Calculates partial AUC up to FDR specified

    TODO determine how to regularize the AUC measurement
        Currently simply divided by total area between 0 and FDR
    """
    if ax is None:
        f, ax = plt.subplots(1,1)
    else:
        f = None

    if is_pval:
        y_pred = - np.log(pred)
    else:
        y_pred = pred

    fpr, tpr, _ = roc_curve(y_act, y_pred)
    # Truncate FPR and TPR
    try:
        idx = next(i for i,v in enumerate(fpr) if v > fdr)
    except StopIteration:
        idx = len(fpr) - 1
    t_fpr, t_tpr = fpr[:idx+1], tpr[:idx+1]
    t_fpr[-1] = fdr

    AUC = auc(t_fpr, t_tpr) / fdr

    lw = 3
    frac_ones = np.mean(y_act)
    roc_ln, = ax.plot(t_fpr, t_tpr, lw=lw, label=label + ' - %.3f' % AUC, color=color)
    base_ln, = ax.plot([0,fdr+0.01], [0,fdr+0.01], color='grey', lw=2, linestyle='--')

    ax.set_xlim(0, fdr)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.set_title('Partial ROC Curve')
    # ax.legend(loc='middle right', fontsize= 'medium', title='AUC')
    return AUC, ax, f


def plot_both(y_act, p_vals, labels, title='', is_pval=True, colors=None, **kwargs):
    """
    Plots comparison of multiple methods
    
    Args:
        y_act: actual labels (vector of ones and zeros, n x 1)
        p_vals: predicted labels (list of vectors of nonnegative pvalues, smaller more significant, k x n)
            If any p_vals are None, the entry is skipped. TODO handle this better
        labels: list of strings for legend, k x 1)
        title: optional, string for figur title
        is_pval: optional, default TRUE. Can be bool or list of bool, determines
            if each set of pvals should be inverted or not
        kwargs: passed through to plotting functions [NOT CURRENTLY USED]

    Returns:
        f, axarr
        f: figure
        axarr: (1 x 2) array of matplotlib axis objects. axarr[0] plots ROC, axarr[1] plots PRC
    """
    f, axarr = plt.subplots(1,2)
    f.suptitle(title, fontsize=16)
    kwargs.pop('ax', None)

    if isinstance(is_pval, bool):
        is_pval = [is_pval] * len(p_vals)

    for i, p_val in enumerate(p_vals):
        if p_val is None:
            continue  # TODO handle this better
        if colors is not None:
            c = colors[i]
        else:
            c = None

        plot_roc(y_act, p_val, ax=axarr[0], label=labels[i], is_pval=is_pval[i], color=c)
        plot_partial_auc(y_act, p_val, ax=axarr[1], label=labels[i], fdr=0.05, is_pval=is_pval[i], color=c)
        # Shaded block on roc plot to indicate area of pAUC plot
        axarr[0].add_patch(
                patches.Rectangle(
                    (0, 0),  # Position
                    0.05,    # Width
                    1.1,     # Height (+ buffer)
                    alpha=0.1,
                    color='grey',
                    edgecolor=None
        ))
        # plot_prc(y_act, p_val, ax=axarr[1], label=labels[i], **kwargs)

    return f, axarr

def plot_pvalue_dist(pvals, axes=None):
    """
    Pvals should be df
    """
    # Plot histogram for this.
    # log-log pvalue plot to focus on smallest p values
    # TODO add labels

    # Don't plot fold change
    valid_labels = list(pvals.columns)
    valid_labels.remove('fold change')
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey='row', squeeze=True)
    else:
        f = None

    # f, (hist_axs, log_axs) = plt.subplots(2, m, sharey='row', squeeze=False)
    
    axes[0].set_ylabel('Density')
    # log_axs[0].set_ylabel('Observed p-value')
    
    for i, name in enumerate(valid_labels):
        ax = axes[i]
        # lax = log_axs[i]

        _, _, rects = ax.hist(pvals[name], bins=20, range=(0,1), normed=True, alpha=0.5)
        ax.plot([0, 1], [1, 1], color='grey', lw=1, linestyle='--')
        ax.set_title(name)
        ax.set_xlabel('p-value')
        # exp_pvals = (np.argsort(np.argsort(pval)) + 0.5) / (len(pval) + 1)
        # lax.scatter(exp_pvals, pval)
        # lax.plot([0, 1], [0 ,1], color='grey', lw=1, linestyle='--')
        # lax.set_xlim(0, 1)
        # lax.set_ylim(0, 1)
        # lax.set_xlabel('Expected p-value')
        
        # TODO plot confidence interval

    return f


def volcano_plots(pval_df, ctrl, exp, is_changed, axes=None):
    # Don't plot fold change
    valid_labels = sorted(list(pval_df.columns))
    valid_labels.remove('fold change')
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey='row', squeeze=True)
    else:
        f = None

    axes[0].set_ylabel('$-\log_{10}$(p-value)') 

    ALPHA = - np.log10(0.05)

    # Volcano plots
    fc = np.mean(exp, 1) - np.mean(ctrl, 1)
    for i, name in enumerate(valid_labels):
        ax = axes[i]

        ax.scatter(fc.values, -np.log10(pval_df[name].values), 
                   c=is_changed, alpha=0.2)
        (l, r) = ax.get_xlim()
        ax.plot([l, r], [ALPHA, ALPHA], color='grey', lw=1, linestyle='--')
        ax.set_title(name if name not in LABEL_MAPPING else LABEL_MAPPING[name])
        ax.set_ylim(bottom=0)
        ax.set_xlim(l, r)
        ax.set_xlabel('$\log_2$(fold change)')
    

    return f

def barplot_accuracy(pval_df, is_changed, ax=None):
    valid_labels = sorted(list(pval_df.columns))
    valid_labels.remove('fold change')

    if ax is None:
        f, ax = plt.subplots()
    else:
        f = None
    bax = ax
    # Bar plot
    A = 0.05
    N = len(is_changed)
    SIG = np.sum(is_changed)
    tp = np.array([np.sum((pval_df[name] <= A) & is_changed) 
            for name in valid_labels])
    fn = np.array([(SIG - tpn) for tpn in tp])
    fp = np.array([np.sum((pval_df[name] <= A) & (1-is_changed)) 
        for name in valid_labels])
    tn = np.array([(N - SIG - fpn) for fpn in fp])

    ind = np.arange(len(valid_labels))
    WIDTH = 1.0
    bax.bar(ind, tp, WIDTH, color='black')
    bax.bar(ind, fn, WIDTH, bottom=tp, color='grey')
    bax.bar(ind, fp, WIDTH, bottom=fn + tp, color='lightgray')
    bax.bar(ind, tn, WIDTH, bottom=fp + fn + tp, color='white')

    bax.set_xticks(ind + 0.5)
    bax.set_xticklabels(valid_labels)

    # TODO 
    # https://matplotlib.org/examples/pylab_examples/broken_axis.html

    return ax, f


def extract_y_act_protein(protein_ids, is_changed):
    """Converts peptide level labels to protein labels

    Args:
        protein_ids: length-n vector of protein ids
        is_changed: length-n vector of labels: same length as protein_ids

    Returns:
        length-m vector of 0-1 labels, sorted by protein id 0...n
    """
    df = pd.DataFrame({
        'protein_id': protein_ids,
        'is_changed': np.array(is_changed, dtype=int)
        }).drop_duplicates(inplace=False)
    # protein_df['order'] = np.arange(protein_df.shape[0])

    # joint = protein_df[['protein_id', 'order']].merge(
    #         df, 
    #         on='protein_id', 
    #         how='left')
    # joint = joint.set_index('order').ix[np.arange(protein_df.shape[0]),:]
    # return joint['is_changed']

    return df.sort_values('protein_id')['is_changed']


def plot_multiple_fc(multiple_fc_result, fold_changes, labels, title=""):
    """ Plots "partial" ROC with respect to FC
    TODO for now just plots pAUC
    Maybe change API . . .

    multiple_fc_result: 4D numpy array
        [i][j][k] = (AUROC, AUPRC, pAUROC) for trial i, fold change index j,
        and label k
    fold_changes: 1D array of fold changes, length must equal fc_result.shape[1]
    labels: 1D array of string labels, length must equal fc_result.shape[2]
    """
    # TODO check invariants asserted above

    # Mean ROC
    mean = np.nanmean(multiple_fc_result, axis=0)
    # Stderr ROC
    std = np.nanstd(multiple_fc_result, axis=0)

    f, ax = plt.subplots()

    for i, label in enumerate(labels):
        lbl_mn = mean[:,i,2]
        lbl_std = std[:,i,2]
        mean_ln, = ax.plot(fold_changes, lbl_mn, lw=2, label=label)
        ax.fill_between(fold_changes, lbl_mn-lbl_std, lbl_mn+lbl_std,
            facecolor=mean_ln.get_color(), alpha=0.3, interpolate=True)
    
    ax.set_ylim(0, 1)
    ax.set_xlabel('log2(fold_change)')
    ax.set_ylabel('pAUC')
    ax.set_title('pAUC versus fold_change for multiple FC in one experiment, %s' % title)
    ax.legend(loc='bottom right', fontsize= 'medium', title='AUC')


def roc_prc_scores(y_act, p_val, is_pval=True, fdr=0.05):
    """ Calculates AUROC, AUPRC, and pAUROC statistics

    Args:
        y_act: actual labels (vector of ones and zeros, n x 1)
        p_vals: predicted labels (list of vectors of nonnegative pvalues, smaller more significant, k x n)

    Returns:
        roc_auc, prc_aucs, pauc
        If p_val is a 1D array, roc_auc and prc_auc are floats
        If p_val is a list of lists, roc_auc and prc_auc are lists of floats
    """
    pval_is_list = hasattr(p_val[0], '__iter__')
    if not pval_is_list:
        p_val = [p_val]

    roc_auc = []
    prc_auc = []
    pauc = []
    for p in p_val:
        if p is None:
            pauc.append(np.nan)
            roc_auc.append(np.nan)
            prc_auc.append(np.nan)
            continue

        if is_pval:
            y_pred = - np.log(p)
        else:
            y_pred = p

        if np.all(np.isnan(y_pred)):
            # No valid p-vals
            pauc.append(np.nan)
            roc_auc.append(np.nan)
            prc_auc.append(np.nan)
            continue
        
        # Rank measurements with nan as lowest priority
        y_pred[np.isnan(y_pred)] = 0

        fpr, tpr, _ = roc_curve(y_act, y_pred)
        # Truncate FPR and TPR for partial auc
        try:
            idx = next(i for i,v in enumerate(fpr) if v > fdr)
        except StopIteration:
            idx = len(fpr) - 1

        t_fpr, t_tpr = fpr[:idx+1], tpr[:idx+1]
        t_fpr[-1] = fdr

        pauc.append(auc(t_fpr, t_tpr) / fdr)
        roc_auc.append(auc(fpr, tpr))
        prec, rec, _ = precision_recall_curve(y_act, y_pred)
        prc_auc.append(auc(rec, prec))

    if not pval_is_list:
        roc_auc = roc_auc[0]
        prc_auc = prc_auc[0] 
        pauc = pauc[0]

    return roc_auc, prc_auc, pauc


def power_analysis(is_changed, pvals, alpha=0.05):
    """ Calculates the false detection rate and power of pvals
        at a given p-value threshold, and with adjustment
        
        Returns (fp, tp, fp_adj, tp_adj)"""

    pval_is_list = hasattr(pvals[0], '__iter__')
    if not pval_is_list:
        pvals = [p_val]

    fps = []
    tps = []
    fps_adj = []
    tps_adj = []
    for pval in pvals:
        # False positives
        fp = np.sum(np.logical_not(is_changed) & (pval <= alpha))
        # Power \approx true_positives / num_changed
        tp = np.sum(is_changed & (pval <= alpha))

        # Adjust p-values using BH and repeat
        _, pval_adj, _, _ = multipletests(pval, alpha=alpha, method='fdr_bh')
        fp_adj = np.sum(np.logical_not(is_changed) & (pval_adj <= alpha))
        tp_adj = np.sum(is_changed & (pval_adj <= alpha))

        fps.append(fp)
        tps.append(tp)
        fps_adj.append(fp_adj)
        tps_adj.append(tp_adj)

    return (fps, tps, fps_adj, tps_adj)


def count_quadrants(pval, fc, is_changed, alpha=0.05):
    """ TODO documentation
    """
    not_changed = np.logical_not(is_changed)
    n_sig = np.sum(pval <= alpha)
    fc_thres = sorted(fc)[n_sig]

    tp_sig_sig = np.sum((pval <= alpha) & (fc <= fc_thres) & is_changed)
    fp_sig_sig = np.sum((pval <= alpha) & (fc <= fc_thres) & not_changed)
    tp_sig_nsig = np.sum((pval <= alpha) & (fc > fc_thres) & is_changed)
    fp_sig_nsig = np.sum((pval <= alpha) & (fc > fc_thres) & not_changed)
    tp_nsig_sig = np.sum((pval > alpha) & (fc <= fc_thres) & is_changed)
    fp_nsig_sig = np.sum((pval > alpha) & (fc <= fc_thres) & not_changed)
    tp_nsig_nsig = np.sum((pval > alpha) & (fc > fc_thres) & is_changed)
    fp_nsig_nsig = np.sum((pval > alpha) & (fc > fc_thres) & not_changed)

    return (tp_sig_sig, fp_sig_sig, tp_sig_nsig, fp_sig_nsig,
            tp_nsig_sig, fp_nsig_sig, tp_nsig_nsig, fp_nsig_nsig)


def density_scatter(data, ax=None):
    """
    Scatterplot of variance vs log2 mean intensity
    """
    x = np.mean(data.values, axis=1)
    y = np.var(data.values, axis=1)
    xy = np.vstack([x,y])

    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    if ax is None:
        fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=10, edgecolor='')
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Variance")
    ax.set_xlabel("$log_2$ mean intensity")

    return ax


def set_fontsize(ax, title=30, axis=20, tick=15):
    """ Set fontsize of ax or axarr """ 
    if not hasattr(ax, '__iter__'):
        ax = [ax]

    for a in ax:
        a.title.set_fontsize(title)
        a.xaxis.label.set_fontsize(axis)
        a.yaxis.label.set_fontsize(axis)
        a.tick_params(axis='both', which='major', labelsize=10)

# Scatterplot of fold change vs cyberT results
"""
plt.scatter(np.log2(np.argsort(np.argsort(pvals.cyberT))[random_sample]), np.log2(np.argsort(np.argsort(pvals['fold change']))[random_sample]), c=is_changed[random_sample], edgecolor=None, cmap='rainbow')

# Count number of sig by CyberT
In [120]: plt.plot(np.log2([2400, 2400]), [-5000, 25000])

In [121]: plt.plot([-5000, 25000], np.log2([2326, 2326]))
Out[121]: [<matplotlib.lines.Line2D at 0x7f53a5bf1650>]

In [122]: plt.ylim(-5, 15)
Out[122]: (-5, 15)

In [123]: plt.xlim(-5, 15)
Out[123]: (-5, 15)

TODO label quadrants, count number in each quadrant
Permute many times
Check performance of combining methods (both significant?)

Get Adobe Illustrator (possibly Photoshop/Powerpoint?)
Making figures: Don't make the text too small. Labels and ticks should be visible
Figures should be modular
Big font, feels too big
Color combination: may want to generate own color scheme
    Consistent for each measure across figures
    Fill vs rim for complex boxplot figures
    Make sure pipeline is standardized before generating figures

Plot FDR and power (separately for each type?)
    (adj and nonadj)
Rerun comparison with and without central variance(!)
    Why is uniform worse than the normal distribution?
    Why is the summary AUROC different than the mean of the buckets???
Is there a generalizable way to look at this data?
    Estimate the power: within a reasonable variance space and statistical limitations
    Absence of evidence is not evidence of absence

Quantify how close the estimated fold changes are to true for WLS
    Correlate TRUE and EST

Protein-level: explore the parameter space
    Statistically principled regularization(?)
    Integrate ANOVA???
"""

"""
Figures:
    Generate 2-sample moderated T
    Consistent color for everything

Writing:
    Look at Molecular Cellular Proteomics for inspiration on results section

    Establish the purpose of the experiment
    Design of the simulation
    Analysis
    Conclusion 1, 2, 3, ...
    Interested in how the parameter space affects the power of the analysis
    Well accepted that parameters affect statistical power
    The extent to which this is true has not been quantified.

    Start at the factual information
    Schematic of data generation
    Place panel of schematic in every plot, highlight changing items
    Be sure the noise distribution is highlighted
    Get main message from panel, inset with information

    Try keeping 2-sample/1-sample as line instead of boxplot
    18.93.6.33

New dataset:
    First-pass: data quality
    Reference channels + normalization
        Do the reference channels correlate well?
    Check variance
    ANOVA or 2-way t-tests if necessary
    DLG4 and DLG2 levels: are they down in the proper channels?

    Do PSM and compare
    Cleaning accession numbers: convert everything to gene symbols if possible
        Remove duplicate symbols, find number of unique peptides

For THU:
    Quality control reference
    Run regular pipeline: pathway + gene list sig


FOR LATER: Show volcano plot - highlight true positive / negative
        - Real world data distribution
    Separate out: uniform, inv gamma, nominal vs adjusted p
    Also quantify the FPR and power for different fold changes

Make sure all parameters are defined in the methods section
Make panels for the figures we've already finalized
    Write up figure legend. Title, A, B, etc.

Flowchart for evaluating proteomic data
    Tables for power analysis
"""


