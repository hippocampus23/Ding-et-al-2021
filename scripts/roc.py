from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

import matplotlib.pyplot as plt
import numpy as np


def plot_roc(y_act, pred, ax = None, is_pval=True):
    # Returns ax, figure
    # Figure is none if ax was provided
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

    lw = 5
    roc_ln, = ax.plot(fpr, tpr, lw=lw, label='ROC (area = %.2f)' % AUC)
    diag_ln, = ax.plot([0,1], [0,1], color='grey', lw=lw, linestyle='--')

    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc='lower right')
    # return ax, f

    return AUC, fpr, tpr


def plot_prc(y_act, pred, ax = None, is_pval=True, label='Area '):
    # Returns ax, figure
    # Figure is none if ax was provided
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

    lw = 5
    frac_ones = np.mean(y_act)
    roc_ln, = ax.plot(rec, prec, lw=lw, label=label + '%.2f' % AUC)
    base_ln, = ax.plot([0,1], [frac_ones,frac_ones], color='grey', lw=lw, linestyle='--')

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.legend(loc='lower left')
    return AUC, ax, f
   

def plot_both(y_act, p_vals, labels, **kwargs):
    # Note: pvals should be a list
    # Plots comparison of multiple methods
    f, axarr = plt.subplots(1,2)
    kwargs.pop('ax', None)

    for i, p_val in enumerate(p_vals):
        plot_roc(y_act, p_val, ax=axarr[0], label=labels[i], **kwargs)
        plot_prc(y_act, p_val, ax=axarr[1], label=labels[i], **kwargs)

    return f, axarr
