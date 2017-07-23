import numpy as np
import statsmodels.api as sm
import sys
import time

sys.dont_write_bytecode = True
from rpy2.robjects import r, numpy2ri, pandas2ri
from scipy import optimize
from statsmodels.tools import add_constant

from sample import *
from roc import *

r['source']('wls.test.R')
numpy2ri.activate()

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


def time_run_wls_statsmodels(n_runs, n_samps, onesample=False):
    start = time.time() 
   
    keep_cols = ['meanC', 'meanE', 'stdC', 'stdE', 'bayesSDC', 'bayesSDE', 'nC', 'nE', 'pVal']
    #   tt <- -(m1 - m2)/sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/
    #      (n1 + n2 - 2)) * ((n1 + n2)/(n1 * n2)))
    
    def wls_degenerate(data):
        data = data.ix[0]
        nC, nE = data['nC'], data['nE']
        if use_bayes:
            stdC, stdE = data['bayesSDC'], data['bayesSDE']
        else:
            stdC, stdE = data['stdC'], data['stdE']
        return (data['meanE'] - data['meanC'],
                (((nC-1)*(stdC**2) + (nE-1)*(stdE**2)) / 
                    (nC+nE-2)) * ((nC+nE) / (nC*nE)),
                data['pVal'])

    def wls(data, use_bayes=False):
        """ Weighted least squares for peptides in protein.
            Operates on sub data frames """
        # Degenerate case, only one peptide
        if data.shape[0] == 1:
            return wls_degenerate(data)

        y = data['meanC'].append(data['meanE']).values
        if use_bayes:
            w = data['bayesSDC'].append(data['bayesSDE']).values**2
        else:
            w = data['stdC'].append(data['stdE']).values**2
        x = np.ones(data.shape[0]*2)
        x[:data.shape[0]] = 0
        
        mod_wls = sm.WLS(y, add_constant(x, prepend=False), weights=1./w)
        res_wls = mod_wls.fit()
        return (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0])

    def wls_onesample(data, use_bayes=False):
        """ Weighted least squares for one-sample peptides in protein
        """
        # Degenerate case, only one peptide
        if data.shape[0] == 1:
            return wls_degenerate(data)

        y = data['meanE'].values - data['meanC'].values
        x = np.ones(data.shape[0]) 
        nC, nE = data['nC'].values.astype(float), data['nE'].values.astype(float)
        if use_bayes:
            stdC, stdE = data['bayesSDC'].values, data['bayesSDE'].values
        else:
            stdC, stdE = data['stdC'].values, data['stdE'].values
        w = (((nC-1)*(stdC**2) + (nE-1)*(stdE**2)) / (nC+nE-2)) * ((nC+nE) / (nC*nE))
        # Weighted least squares with only constant term
        mod_wls = sm.WLS(y, x, weights=1./w)
        res_wls = mod_wls.fit()
        return (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0])

    rocs = np.zeros((n_runs, 3), dtype=float)

    for i in xrange(n_runs):
        # Generate proteins and do tests on peptides
        ctrl, exp, is_changed, proteins = sample_proteins(5000, 500, 1.5, n_samps)
        cyberT_peps = cyberT(ctrl, exp)[keep_cols]
        cyberT_peps['protein'] = proteins

        # Group by protein
        grouped = cyberT_peps.groupby('protein', sort=False)
        res = pd.DataFrame(columns=['protein_id','x','std', 'pval'], 
                           index=np.arange(len(grouped)))
        res['protein_id'] = res['protein_id'].astype(str)
        res['x'] = res['x'].astype(float)
        res['std'] = res['std'].astype(float)
        res['pval'] = res['pval'].astype(float)
        for j, (name, group) in enumerate(grouped):
            # Get the columns corresponding to meanC and meanE
            if onesample:
                x, std, pval = wls_onesample(group)
            else:
                x, std, pval = wls(group)
            res.ix[j] = [name, x, std, pval]

        rocs[i,:] = roc_prc_scores(
                extract_y_act_protein(res, proteins, is_changed).values,
                res['pval'].values)
        

    end = time.time()
    print "Total: %.0f, average %.2f" % (end-start, (end-start)/n_runs)
    return rocs

def time_run_wls_r(n_runs, n_samps):
    pass
