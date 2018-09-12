import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from statsmodels.sandbox.stats.multicomp import multipletests
from constants import FDR

def roc_prc_scores(is_changed, p_vals):
    """
    Calculates AUROC, AUPRC, and pAUROC statistics
    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, k x n)
    :return:           roc_auc, prc_aucs, pauc (all lists of floats)

    """
    roc_auc = []
    prc_auc = []
    pauc = []
    for pval in p_vals.values.transpose():
        # avoid trouble with -inf in roc_curve
        # replace 0 with an arbitrary small number
        pval = [p if p != 0 else 1e-100 for p in pval]

        predicted = - np.log(pval)

        # e.g. when values are missing because a test could not be applied
        if np.all(np.isnan(predicted)):
            # No valid p-vals
            pauc.append(np.nan)
            roc_auc.append(np.nan)
            prc_auc.append(np.nan)
            continue

        # Rank measurements with nan as lowest priority
        predicted[np.isnan(predicted)] = 0
        
        fpr, tpr, _ = roc_curve(is_changed, predicted)

        # Truncate FPR and TPR for partial auc
        try:
            idx = next(i for i, v in enumerate(fpr) if v > FDR)
        except StopIteration:
            idx = len(fpr) - 1

        t_fpr, t_tpr = fpr[:idx + 1], tpr[:idx + 1]
        t_fpr[-1] = FDR

        pauc.append(auc(t_fpr, t_tpr) / FDR)
        roc_auc.append(auc(fpr, tpr))
        prec, rec, _ = precision_recall_curve(is_changed, predicted)
        prc_auc.append(auc(rec, prec))

    return roc_auc, prc_auc, pauc


def power_analysis(is_changed, pvals, alpha=0.05):
    """
    Calculates the false detection rate and power of pvals
    at a given p-value threshold, and with adjustment

    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, n x k)
    :param alpha:      FWER, family-wise error rate, 0.05 by default in multipletests
    :return:           numpy.ndarray (N_RUNS x number of tests x 4), arr[i][j] = (FP_raw, TP_raw, FP_adj, TP_adj)
    """
    fps = []
    tps = []
    fps_adj = []
    tps_adj = []
    for pval in pvals.values.transpose():
        #  e.g. when values are missing because the a test could not be applied
        if np.all(np.isnan(pval)):
            fp, tp, fp_adj, tp_adj = 0, 0, 0, 0

        else:
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

    return fps, tps, fps_adj, tps_adj
