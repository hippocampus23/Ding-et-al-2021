from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

import matplotlib.pyplot as plt
import numpy as np


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
   

def plot_both(y_act, p_vals, labels, **kwargs):
    """
    Plots comparison of multiple methods
    
    Args:
        y_act: actual labels (vector of ones and zeros, n x 1)
        p_vals: predicted labels (list of vectors of nonnegative pvalues, smaller more significant, k x n)
            If any p_vals are None, the entry is skipped. TODO handle this better
        labels: list of strings for legend, k x 1)
        kwargs: passed through to plotting functions. Currently supports is_pval

    Returns:
        f, axarr
        f: figure
        axarr: (1 x 2) array of matplotlib axis objects. axarr[0] plots ROC, axarr[1] plots PRC
    """
    f, axarr = plt.subplots(1,2)
    kwargs.pop('ax', None)

    for i, p_val in enumerate(p_vals):
        if p_val is None:
            continue  # TODO handle this better
        plot_roc(y_act, p_val, ax=axarr[0], label=labels[i], **kwargs)
        plot_prc(y_act, p_val, ax=axarr[1], label=labels[i], **kwargs)

    return f, axarr
