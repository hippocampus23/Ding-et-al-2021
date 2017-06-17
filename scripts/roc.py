from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_roc(y_act, pred, ax = None, is_pval=True, label='Area '):
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
    roc_ln, = ax.plot(fpr, tpr, lw=lw, label=label + ' - %.3f' % AUC)
    diag_ln, = ax.plot([0,1.01], [0,1.01], color='grey', lw=2, linestyle='--')
    
    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC curve')
    ax.legend(loc='lower right', fontsize='medium', title='AUC')
    # return ax, f

    return AUC, fpr, tpr


def plot_prc(y_act, pred, ax = None, is_pval=True, label='Area '):
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
    roc_ln, = ax.plot(rec, prec, lw=lw, label=label + ' - %.3f' % AUC)
    base_ln, = ax.plot([0,1.01], [frac_ones,frac_ones], color='grey', lw=2, linestyle='--')

    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('PRC curve')
    ax.legend(loc='lower left', fontsize= 'medium', title='AUC')
    return AUC, ax, f


def plot_partial_auc(y_act, pred, fdr=0.05, ax=None, is_pval=True, label='Area '):
    """Calculates partial AUC up to FDR specified

    TODO determine how to regularize the AUC measurement (!)
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
    idx = next(i for i,v in enumerate(fpr) if v > fdr)
    t_fpr, t_tpr = fpr[:idx+1], tpr[:idx+1]
    t_fpr[-1] = fdr

    AUC = auc(t_fpr, t_tpr) / fdr

    lw = 3
    frac_ones = np.mean(y_act)
    roc_ln, = ax.plot(t_fpr, t_tpr, lw=lw, label=label + ' - %.3f' % AUC)
    base_ln, = ax.plot([0,fdr+0.01], [0,fdr+0.01], color='grey', lw=2, linestyle='--')

    ax.set_xlim(0, fdr+0.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.set_title('Partial ROC Curve')
    ax.legend(loc='middle right', fontsize= 'medium', title='AUC')
    return AUC, ax, f


def plot_both(y_act, p_vals, labels, title='', **kwargs):
    """
    Plots comparison of multiple methods
    
    Args:
        y_act: actual labels (vector of ones and zeros, n x 1)
        p_vals: predicted labels (list of vectors of nonnegative pvalues, smaller more significant, k x n)
            If any p_vals are None, the entry is skipped. TODO handle this better
        labels: list of strings for legend, k x 1)
        title: optional, string for figur title
        kwargs: passed through to plotting functions. Currently supports is_pval

    Returns:
        f, axarr
        f: figure
        axarr: (1 x 2) array of matplotlib axis objects. axarr[0] plots ROC, axarr[1] plots PRC
    """
    f, axarr = plt.subplots(1,2)
    f.suptitle(title, fontsize=16)
    kwargs.pop('ax', None)

    for i, p_val in enumerate(p_vals):
        if p_val is None:
            continue  # TODO handle this better
        plot_roc(y_act, p_val, ax=axarr[0], label=labels[i], **kwargs)
        plot_partial_auc(y_act, p_val, ax=axarr[1], label=labels[i], fdr=0.05, **kwargs)
        # plot_prc(y_act, p_val, ax=axarr[1], label=labels[i], **kwargs)

    return f, axarr


def extract_y_act_protein(protein_df, protein_ids, is_changed):
    """Converts peptide level labels to protein labels

    Args:
        protein_df: Pandas df which has a column labeled 'accession_number'
        protein_ids: length-n vector of protein ids
        is_changed: length-n vector of labels: same length as protein_ids

    Returns:
        length-m vector of 0-1 labels, corresponds to order of acc_nums in protein_df
    """

    df = pd.DataFrame({
        'protein_id': protein_ids,
        'is_changed': is_changed
        }).drop_duplicates(inplace=False)

    joint = protein_df.join(
            df.set_index('protein_id'), 
            on='accession_number', 
            how='inner')
    return joint['is_changed']


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
            continue

        if is_pval:
            y_pred = - np.log(p)
        else:
            y_pred = p

        fpr, tpr, _ = roc_curve(y_act, y_pred)
        # Truncate FPR and TPR for partial auc
        idx = next(i for i,v in enumerate(fpr) if v > fdr)
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


