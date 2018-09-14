import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc

from constants import LABEL_MAPPING, COLORS


def plot_density_scatter(mean_data, var_data, ax):
    """
    Generate a mean variance density scatter plot

    :param mean_data:   averaged peptide intensities
    :param var_data:    peptide variance
    :param ax:          where to plot
    """

    x = mean_data
    y = var_data
    xy = np.vstack([x,y])

    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, s=10, edgecolor="")
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Peptide variance")
    ax.set_xlabel("Mean $log2$ peptide intensity")


def make_colorbar(ax):
    """
    Make a legend for the density scatter
    Note: for best results, call this after tight_layout, and make space for it"s labels
          if necessary using subplots_adjust

    :param ax:    where to plot
    """

    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["Low", "High"])
    cbar.set_label("Arbitrary Density", labelpad=-16)


def plot_roc_curves(is_changed, pvals):
    """
    Plot roc, pauc, and prc curves side by side

    :param is_changed:   binary labels for each peptide describing whether or not it was perturbed
    :param pvals:        results from statistical tests
    :return:             plot
    """

    font = {"weight": "normal",
            "size": 8}
    matplotlib.rc("font", **font)
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4))

    # Edit column titles for pretty printing
    pvals.columns = [LABEL_MAPPING[l] if l in LABEL_MAPPING else l
                     for l in list(pvals.columns)]
    labels = sorted(list(pvals.columns))

    colors = COLORS.values()
    for i, lbl in enumerate(labels):
        c = colors[i]
        p_val = pvals[lbl]
        plot_roc(is_changed, p_val, ax=ax1, label=lbl, color=c)
        plot_partial_auc(is_changed, p_val, ax=ax2, label=lbl, fdr=0.05, color=c)
        plot_prc(is_changed, p_val, ax=ax3, label=lbl, color=c)

    # Shaded block on roc plot to indicate area of pAUC plot
    ax1.add_patch(
        matplotlib.patches.Rectangle(
            (0, 0),  # Position
            0.05,  # Width
            1.1,  # Height (+ buffer)
            alpha=0.1,
            color="grey",
            edgecolor=None
        ))
    ax1.legend(loc="lower right", fontsize="medium")
    f.tight_layout()
    return f


def volcano_multipanel(pvals, ctrl, exp, is_changed):
    """
    Volcano plots for each test using adjusted and raw p-values

    :param pvals:       results from statistical tests
    :param ctrl:        control data
    :param exp:         experiment data
    :param is_changed:  binary labels for each peptide describing whether or not it was perturbed
    :return:            plot
    """

    pvals.drop("fold change", axis=1, inplace=True)

    # Edit column titles for pretty printing
    pvals.columns = [LABEL_MAPPING[l] if l in LABEL_MAPPING else l
                       for l in list(pvals.columns)]

    matplotlib.rc("font", size=10)
    nc = len(list(pvals.columns))
    f, axarr = plt.subplots(2, nc, sharex="col", sharey="row", figsize=(20, 10))

    pvals_corrected = pd.DataFrame.from_items([
        (col, multipletests(pvals[col], 0.05, method="fdr_bh")[1])
        for col in pvals.columns
    ])

    volcano_plots(pvals, ctrl, exp, is_changed, axarr[0])
    volcano_plots(pvals_corrected, ctrl, exp, is_changed, axarr[1])

    # Remove unnecessary labels
    for ax in axarr[0]:
        ax.set_xlabel("")
    for ax in axarr[1]:
        ax.set_title("")
    axarr[1][0].set_ylabel("$-\log_{10}$(Adjusted p-value)")

    f.tight_layout()

    return f


def plot_pvalue_dist(pvals, axes=None):
    """
    Plot the p-value distribution for each test 

    :param pvals:   results from statistical tests, should be a pandas.DataFrame
    :param axes:    optional, where to plot (subplots)
    :return:        plot unless axes was specified, then None
    """

    valid_labels = sorted(list(pvals.columns))
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey="row", squeeze=True)
    else:
        f = None

    axes[0].set_ylabel("Density")

    for i, name in enumerate(valid_labels):
        ax = axes[i]

        _, _, rects = ax.hist(pvals[name], bins=20, range=(0, 1), density=True, alpha=0.5)
        ax.plot([0, 1], [1, 1], color="grey", lw=1, linestyle="--")
        ax.set_title(name if name not in LABEL_MAPPING else LABEL_MAPPING[name])
        ax.set_xlabel("p-value")

        # TODO plot confidence interval

    return f


def plot_roc(is_changed, pvals, ax=None, label="Area ", color=None):
    """
    Plot the auc curve for each test

    :param is_changed: binary labels for each peptide describing whether or not it was perturbed
    :param pvals:      results from statistical tests
    :param ax:         optional, where to plot
    :param label:
    :param color:
    :return:           AUC, ax, f
                       ROC AUC data, plot, figure is None if ax was specified
    """
    if ax is None:
        f, ax = plt.subplots(1, 1)
    else:
        f = None

    predicted = - np.log(pvals)

    fpr, tpr, _ = roc_curve(is_changed, predicted)
    AUC = roc_auc_score(is_changed, predicted)

    lw = 3
    roc_ln, = ax.plot(fpr, tpr, lw=lw, label=label, color=color)
    diag_ln, = ax.plot([0, 1.01], [0, 1.01], color="grey", lw=2, linestyle="--")

    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC curve")

    return AUC, ax, f


def plot_prc(is_changed, pvals, ax=None, label="Area ", color=None):
    """
    Plot the precision recall curve for each test

    :param is_changed: binary labels for each peptide describing whether or not it was perturbed
    :param pvals:      results from statistical tests
    :param ax:         optional, where to plot
    :param label:
    :param color:
    :return:           AUC, ax, f
                       PRC data, plot, figure is None if ax was specified
    """

    if ax is None:
        f, ax = plt.subplots(1, 1)
    else:
        f = None

    y_pred = - np.log(pvals)

    prec, rec, _ = precision_recall_curve(is_changed, y_pred)
    AUC = auc(rec, prec)

    lw = 3
    frac_ones = np.mean(is_changed)
    roc_ln, = ax.plot(rec, prec, lw=lw, label=label + " - %.3f" % AUC, color=color)
    base_ln, = ax.plot([0, 1.01], [frac_ones, frac_ones], color="grey", lw=2, linestyle="--")

    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("PRC curve")
    return AUC, ax, f


def plot_partial_auc(is_changed, pvals, fdr=0.05, ax=None, label="Area ", color=None):
    """
    Calculates partial AUC up to FDR specified
    Plots the partial AUC curve for each test

    :param is_changed: binary labels for each peptide describing whether or not it was perturbed
    :param pvals:      results from statistical tests
    :param fdr:        false discovery rate
    :param ax:         optional, where to plot
    :param label:
    :param color:
    :return:           AUC, ax, f
                       partial AUC data, plot, figure is None if ax was specified
    """

    # TODO: determine how to regularize the AUC measurement
    # Currently simply divided by total area between 0 and FDR

    if ax is None:
        f, ax = plt.subplots(1, 1)
    else:
        f = None

    y_pred = - np.log(pvals)

    fpr, tpr, _ = roc_curve(is_changed, y_pred)
    # Truncate FPR and TPR
    try:
        idx = next(i for i, v in enumerate(fpr) if v > fdr)
    except StopIteration:
        idx = len(fpr) - 1
    t_fpr, t_tpr = fpr[:idx + 1], tpr[:idx + 1]
    t_fpr[-1] = fdr

    AUC = auc(t_fpr, t_tpr) / fdr

    lw = 3
    frac_ones = np.mean(is_changed)
    roc_ln, = ax.plot(t_fpr, t_tpr, lw=lw, label=label + " - %.3f" % AUC, color=color)
    base_ln, = ax.plot([0, fdr + 0.01], [0, fdr + 0.01], color="grey", lw=2, linestyle="--")

    ax.set_xlim(0, fdr)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Partial ROC Curve")
    return AUC, ax, f


def volcano_plots(pval_df, ctrl, exp, is_changed, axes=None):
    """
    Volcano plots for all tests

    :param pval_df:     results from statistical tests (maybe adjusted), pandas.DataFrame
    :param ctrl:        control data
    :param exp:         experiment data
    :param is_changed:  binary labels for each peptide describing whether or not it was perturbed
    :param axes:        optional, where to plot
    :return:            figure, None if axes specified
    """
    valid_labels = sorted(list(pval_df.columns))
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey="row", squeeze=True)
    else:
        f = None

    axes[0].set_ylabel("$-\log_{10}$(p-value)")

    ALPHA = - np.log10(0.05)

    order = np.argsort(is_changed)
    pval_df = pval_df.iloc[order]
    ctrl = ctrl.iloc[order]
    exp = exp.iloc[order]
    is_changed = is_changed[order]

    # Volcano plots
    fc = np.mean(exp, 1) - np.mean(ctrl, 1)
    for i, name in enumerate(valid_labels):
        ax = axes[i]

        ax.scatter(fc.values, -np.log10(pval_df[name].values),
                   c=is_changed, alpha=0.25, cmap=matplotlib.cm.coolwarm, lw=0, s=3)
        (l, r) = ax.get_xlim()
        ax.plot([l, r], [ALPHA, ALPHA], color="grey", lw=1, linestyle="--")
        ax.set_title(name if name not in LABEL_MAPPING else LABEL_MAPPING[name])
        ax.set_ylim(bottom=0)
        ax.set_xlim(l, r)
        ax.set_xlabel("$\log_2$(fold change)")

    return f


def barplot_accuracy(pval_df, is_changed, ax1=None, ax2=None):
    """
    Creates stacked barplots comparing the TP, FP, TN, and FN
    of each method in the pval_df

    :param pval_df:    results from statistical tests, pandas.DataFrame
    :param is_changed: binary labels for each peptide describing whether or not it was perturbed
    :param ax1:        optional, where to plot (there are 2 axis above each other for a broken axis plot
    :param ax2:        optional, where to plot
    :return:           figure, None if ax specified
    """

    valid_labels = sorted(list(pval_df.columns))

    if ax1 is None:
        f, (ax1, ax2) = plt.subplots()
    else:
        f = None

    # Bar plot
    A = 0.05
    N = len(is_changed)
    SIG = np.sum(is_changed)
    tp = np.array([np.sum((pval_df[name] <= A) & is_changed)
                   for name in valid_labels])
    fn = np.array([(SIG - tpn) for tpn in tp])
    fp = np.array([np.sum((pval_df[name] <= A) & (1 - is_changed))
                   for name in valid_labels])
    tn = np.array([(N - SIG - fpn) for fpn in fp])

    ind = np.arange(len(valid_labels))
    WIDTH = 1.0
    for bax in (ax1, ax2):
        bax.bar(ind, tp, WIDTH, color="black", label="True Positives")
        bax.bar(ind, fn, WIDTH, bottom=tp, color="grey", label="False Negatives")
        bax.bar(ind, fp, WIDTH, bottom=fn + tp, color="lightgray", label="False Positives")
        bax.bar(ind, tn, WIDTH, bottom=fp + fn + tp, color="white", label="True Negatives")

        bax.set_xticklabels(valid_labels)
        bax.set_xticks(ind)

        bax.set_xticklabels(bax.xaxis.get_majorticklabels(), rotation=90)

    # broken axis for better readability
    ax2.set_ylim(0, 2000)
    ax1.set_ylim(9000, 10000)
    ax1.spines["bottom"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off')

    d = .015
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)              # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)        # top-right diagonal

    kwargs.update(transform=ax2.transAxes)              # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)        # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    return ax1, ax2, f
