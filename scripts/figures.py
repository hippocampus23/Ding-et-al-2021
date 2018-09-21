# TODO: don"t do this!!! it is unsafe but avoids omp error
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

import numpy as np
import pandas as pd
import matplotlib
# to avoid tkinter error on linux VM
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import t
from collections import OrderedDict
from plotting import plot_density_scatter, make_colorbar, plot_roc_curves, plot_pvalue_dist, \
                     volcano_multipanel, barplot_accuracy
from statsmodels.sandbox.stats.multicomp import multipletests
from constants import NUM_PEPTIDES, LABEL_MAPPING, PEPTIDE_VAR, MEAN_PEPTIDE_INTENSITY
from sampler import sample
from stat_tests import do_all_tests, modT, TESTS


def density_scatter(ctrl_data):
    """
    Scatterplot of variance vs log2 mean intensity
    used for figure 1BCDE

    :param ctrl_data: pd.DataFrame containing control data
    :return:          plot
    """

    matplotlib.rc("font", size=10)
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15, 5))

    def do_plot(data, ax):
        means = np.mean(data.values, axis=1)
        vars = np.var(data.values, axis=1)
        plot_density_scatter(means, vars, ax)

    # Real
    rs = np.random.choice(ctrl_data.shape[0], NUM_PEPTIDES, replace=False)
    slct_data = ctrl_data.iloc[rs]

    # Simulated
    # only NUM_ROWS avg_ctrl data provided -> will be based on real data (reordered but none omitted)
    ctrl_u, _, _ = sample("uniform", num_changed=0)
    ctrl_g, _, _ = sample("gamma", num_changed=0)
    ctrl_t, _, _ = sample("trend", num_changed=0)

    do_plot(ctrl_u, ax1)
    do_plot(ctrl_g, ax2)
    do_plot(ctrl_t, ax3)
    do_plot(slct_data, ax4)

    ax1.set_ylabel("Peptide variance")
    ax1.set_title("Uniform")
    ax2.set_title("Inverse Gamma")
    ax3.set_title("Variance Trend")
    ax4.set_title("Empirical Data")
    ax1.set_xlabel("Mean $log2$ peptide intensity")
    ax2.set_xlabel("Mean $log2$ peptide intensity")
    ax3.set_xlabel("Mean $log2$ peptide intensity")
    ax4.set_xlabel("Mean $log2$ peptide intensity")

    f.tight_layout()
    f.subplots_adjust(right=0.94)

    # Legend
    make_colorbar(ax4)

    return f


def plot_example_roc_curves(var_type):
    """
    Figure 1EFG, example curves

    :param var_type:  type of distribution to use to sample peptide variance, either "uniform", "gamma" or "trend"
    :return:          plot
    """

    ctrl, exp, is_changed = sample(var_type)
    pvals = do_all_tests(ctrl, exp)
    return plot_roc_curves(is_changed, pvals)


def pvalue_multipanel():
    """
    Generate panel comparing p-value distributions
    Compare uniform, inverse gamma and trend
    Corresponds to Figure 2BDF of manuscript

    :return:  plot
    """
    matplotlib.rc("font", size=10)
    f, axarr = plt.subplots(3, len(TESTS) + 1, sharex="col", sharey="row", figsize=(20, 15))

    for i, var_type in enumerate(["uniform", "gamma", "trend"]):
        ctrl, exp, _ = sample(var_type, num_changed=0)
        pvals = do_all_tests(ctrl, exp)
        pvals.drop("fold change", axis=1, inplace=True)

        # Edit column titles for pretty printing
        pvals.columns = [LABEL_MAPPING[l] if l in LABEL_MAPPING else l
                for l in list(pvals.columns)]
        pvals["Moderated T \n(2-sample, robust)"] = modT(ctrl, exp, robust=True, run_2sample=True)
        pvals["Moderated T \n(1-sample, robust)"] = modT(ctrl, exp, robust=True)
        plot_pvalue_dist(pvals, axarr[i])
        for ax in axarr[i]:
            ax.set_xlabel("" if i < 2 else "p-value")

    f.tight_layout()
    plt.subplots_adjust(hspace=0.4)

    return f


def volcano_multipanel_example(var_type):
    """
    Generate panel comparing volcano plots
    Compare uniform and inverse gamma
    Corresponds to Figure 4B/G of manuscript

    :param var_type:  type of distribution to use to sample peptide variance, either "uniform", "gamma" or "trend"
    :return:          plot
    """

    ctrl, exp, is_changed = sample(var_type)
    pvals = do_all_tests(ctrl, exp)

    return volcano_multipanel(pvals, ctrl, exp, is_changed)


def barplot_multipanel(var_type):
    """
    Compute FP, TP, FN, TN for raw and adj p-values
    Corresponds to Figure 4CD, 4HI in manuscript

    :param var_type:  type of distribution to use to sample peptide variance, either "uniform", "gamma" or "trend"
    :return:          plot
    """
    matplotlib.rc("font", size=8)
    f, axarr = plt.subplots(2, 2, sharey="row", sharex=True, figsize=(5, 7))

    ctrl, exp, is_changed = sample(var_type)
    pvals = do_all_tests(ctrl, exp)
    pvals.drop("fold change", axis=1, inplace=True)

    # Edit column titles for pretty printing
    pvals.columns = [LABEL_MAPPING[l] if l in LABEL_MAPPING else l
            for l in list(pvals.columns)]

    # Adjust pvals
    pvals_a = pd.DataFrame.from_dict(OrderedDict([
        (col, multipletests(pvals[col], 0.05, method="fdr_bh")[1])
        for col in pvals.columns
    ]))

    # Plot
    barplot_accuracy(pvals, is_changed, axarr[0][0], axarr[1][0])
    barplot_accuracy(pvals_a, is_changed, axarr[0][1], axarr[1][1])
    axarr[1][0].set_ylabel("Count")
    axarr[0][0].set_title("Raw p-values")
    axarr[0][1].set_title("BH adjusted")

    f.tight_layout()

    # Resize ax
    for ax in axarr[0]:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.5])

    # Add legend
    handles, labels = axarr[0][1].get_legend_handles_labels()
    plt.figlegend(handles, labels, loc="upper center", ncol=2)

    return f


def fig_1BCDE():
    data = pd.read_csv("../data_empirical/data.csv", sep=";", usecols=["C1", "C2", "C3"])
    plot = density_scatter(data)
    plot.savefig("../figures/1BCDE.eps")
    print "saved 1BCDE.eps"


def fig_1FGH(var_type):
    plot = plot_example_roc_curves(var_type)
    plot.savefig("../figures/1FGH_"+var_type+".eps")
    print "saved 1FGH_"+var_type+".eps"


def fig_2BDF():
    plot = pvalue_multipanel()
    plot.savefig("../figures/2BDF.eps")
    print "saved 2BDF.eps"


def fig_4B():
    plot = volcano_multipanel_example("uniform")
    plot.savefig("../figures/4B.eps")
    print "saved 4B.eps"


def fig_4G():
    plot = volcano_multipanel_example("gamma")
    plot.savefig("../figures/4G.eps")
    print "saved 4G.eps"


def fig_4L():
    plot = volcano_multipanel_example("trend")
    plot.savefig("../figures/4L.eps")
    print "saved 4L.eps"


def fig_4CD():
    plot = barplot_multipanel("uniform")
    plot.savefig("../figures/4CD.eps")
    print "saved 4CD.eps"


def fig_4HI():
    plot = barplot_multipanel("gamma")
    plot.savefig("../figures/4HI.eps")
    print "saved 4HI.eps"


def fig_4MN():
    plot = barplot_multipanel("trend")
    plot.savefig("../figures/4MN.eps")
    print "saved 4MN.eps"


# to be added to 6C (plotted in R)
# shapes of the noise distributions
def fig_6C_dists():
    matplotlib.rc("font", size=12)
    fig, (a1, a2, a3) = plt.subplots(3, 1, sharex="col", figsize=(6, 13))
    x = np.linspace(-2, 2, 200)
    normal = 1/np.sqrt(2 * np.pi * PEPTIDE_VAR) * np.exp(- x**2 / (2 * PEPTIDE_VAR))
    a1.plot(x, normal)
    a1.fill_between(x, normal)
    a1.set_title("Normal distribution,\n$\mu = 0$, $\sigma^2 = %0.2f$" % PEPTIDE_VAR)
    laplacian = np.exp(-abs(x) /PEPTIDE_VAR**0.5)/(2.*PEPTIDE_VAR**0.5)
    a2.plot(x, laplacian)
    a2.fill_between(x, laplacian)
    a2.set_title("Laplacian distribution,\n$\mu = 0$, $\sigma^2 = %0.2f$" % PEPTIDE_VAR)
    scaled_t = t.pdf(x / PEPTIDE_VAR**0.5, df=3) / PEPTIDE_VAR**0.5
    a3.plot(x, scaled_t)
    a3.fill_between(x, scaled_t)
    a3.set_title("Scaled t distribution,\n$df = 3$, scaled using $\sigma^2 = %0.2f$" % PEPTIDE_VAR)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.6)
    fig.savefig("../figures/6C_dists.eps")
    print "saved 6C_dists.eps"
