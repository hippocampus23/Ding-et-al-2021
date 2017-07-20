import numpy as np
from scipy import optimize


def run_test(n, means=None, sigs=None):
    """ n is number of samples to optimize over
    """
    means = np.array(means, dtype=float)
    sigs = np.array(sigs, dtype=float)
    if means is None:
        means = np.random.normal(size=n)
    if sigs is None:
        sigs = 1. / np.random.gamma(3, 1, size=n)

    def obj(x):
        u, sig = x
        adj_vars = sigs + sig

        return(np.sum(adj_vars) +
               np.sum((means - u)**2 / adj_vars))
    
    def grad(x):
        u, sig = x
        adj_vars = sigs + sig
        
        du = np.sum((u - means) / adj_vars)
        dsig = 0.5*np.sum(adj_vars) - 0.5*np.sum(((u-means) / adj_vars)**2)
        return np.array([du, dsig])

    x_0 = np.array([np.mean(means), np.var(means)])
    out = optimize.minimize(obj, x_0, bounds=[(None, None), (0, None)], method='TNC')
    return means, sigs, out
