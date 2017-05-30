""" Module for visualizing predicted and heavy tails of distributions
"""

# TODO easier to just plot regular CDF? No, doesn't emphasize very well
# TODO simulate the empirical abs quantile change from normal and laplace
#   GIVEN normalization and difference variances

import matplotlib.pyplot as plt
import numpy as np

from scipy import stats


def absolute_scaled_deviations_from_mean(data):
    """
    Converts data distribution into CDF of absolute deviations from mean

    Returns: quantiles, cdf
        quantiles : n X 1 vector of absolute deviations from mean, sorted increasing
        cdf: n x 1 vector where i corresponds to CDF of 
            absolute deviation at ith entry of quantiles
    """

    data = data.flatten()
    n = len(data)

    mean = np.mean(data)

    sorted_data = np.sort(np.abs(data - mean))
    cdf = np.arange(0, 1.0, 1./n) 

    return sorted_data, cdf


def theoretical_cdf_normal(quantiles):
    """
    Calculates the theoretical CDF of absolute deviation
        from mean of normal distribution with var=1
    """

    var = np.sum(quantiles**2) / len(quantiles)
    print "Variance: %f" % var

    return 2 * (stats.norm.cdf(quantiles, scale=var**0.5) - 0.5)


def theoretical_cdf_laplace(quantiles):
    """
    Calculates the theoretical CDF of absolute deviation
        from mean of laplace distribution with var=1
    """
    dev = np.sum(quantiles) / len(quantiles)
    print "Sum of absolute deviations: %f" % dev

    return 2 * (stats.laplace.cdf(quantiles, scale=dev) - 0.5)


def theoretical_cdf_t(quantiles, df=3):
    """
    Calculates the theoretical CDF of absolute deviation
        from mean of t distribution with given df
    """
    
    var = np.sum(quantiles**2) / len(quantiles)

    return 2*(stats.t.cdf(quantiles, df, scale=var**0.5) - 0.5)


def _mean_recenter(data):
    """
    Recenters each row of 2D matrix around 0
    Subtracts mean of every row
    """

    return data - np.mean(data,1)[:,None]


def compare_theoretical_actual_quantiles(data, mean_recenter=False, ax=None):
    """
    Plots the actual vs theoretical expected quantiles for data
    """

    if ax is None:
        f, ax = plt.subplots()
    else:
        f = None

    if mean_recenter:
        data = _mean_recenter(data)

    data = data.flatten()
    quantiles, empirical_cdf = absolute_scaled_deviations_from_mean(data)

    normal_cdf = theoretical_cdf_normal(quantiles)
    laplace_cdf = theoretical_cdf_laplace(quantiles)
    t_cdf = theoretical_cdf_t(quantiles)

    # TODO truncate vertical axis at lowest possible empirical
    ax.set_ylim(np.log10(1 - empirical_cdf[-1])-1, 0)

    lw = 2
    ax.plot(quantiles, np.log10(1 - empirical_cdf), lw=lw, label="Empirical CDF")
    ax.plot(quantiles, np.log10(1 - normal_cdf), lw=lw, label="Theoretical Normal CDF")
    ax.plot(quantiles, np.log10(1 - laplace_cdf), lw=lw, label="Theoretical Laplacian CDF")
    ax.plot(quantiles, np.log10(1 - t_cdf), lw=lw, label="Theoretical t-distribution CDF")
    ax.legend(loc='lower left', fontsize='medium')

    ax.set_xlabel('Absolute deviation from mean')
    ax.set_ylabel('$\log_{10}(1 - CDF)$')

    return ax, f


def visualize_fits(data, mean_recenter=False, ax=None):
    """
    Compares normal and laplacian fits to given data

    If necessary, mean centers the data before comparison
    """

    if ax is None:
        f, ax = plt.subplots()
    else:
        f = None

    if mean_recenter:
        data = _mean_recenter(data)

    data = data.flatten()

    _, bins, _ = ax.hist(data, bins=100, normed=True)

    mean = np.mean(data)
    var = np.var(data)
    med = np.median(data)
    dev = (1./len(data)) * np.sum(np.abs(data - med))
    print "Variance: %f" % var
    print "Sum of absolute deviations: %f" % dev

    fit_norm = stats.norm.pdf(bins, loc=mean, scale=var**0.5)
    fit_laplace = stats.laplace.pdf(bins, loc=med, scale=dev)
    norm_ll = np.sum(np.log(stats.norm.pdf(data, loc=mean, scale=var**0.5)))
    lap_ll = np.sum(np.log(stats.laplace.pdf(data, loc=med, scale=dev)))

    ax.plot(bins, fit_norm, label="MLE normal fit: LL = %.2f" % norm_ll)
    ax.plot(bins, fit_laplace, label="MLE laplacian fit: LL = %.2f" % lap_ll)
    ax.legend()
    ax.set_ylabel('Density')
    ax.set_title('Empirical and fitted distributions of data')

    return ax, f


def visualize_dist_and_tails(data, mean_recenter=False, title=''):

    f, (ax1, ax2) = plt.subplots(2,1)
    f.suptitle(title, fontsize=16)

    visualize_fits(data, mean_recenter=mean_recenter, ax=ax1)
    compare_theoretical_actual_quantiles(data, mean_recenter=mean_recenter, ax=ax2)
