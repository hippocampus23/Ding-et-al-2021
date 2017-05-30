"""
Scratchpad for random functions and fragments which are useful
in the IPython notebook

IMPORTANT
No imports because these functions are NOT meant to be used in a script
"""
DEF_FOLD_CHANGES = [2**i for i in np.arange(0.1, 1.1, 0.1)]

def simulate_fold_change_range(
        fold_changes = DEF_FOLD_CHANGES,
        var = 0.06,
        nctrl=3,
        nexp=3):
    """Creates simulated datasets for multiple fold changes

    Uses a SINGLE fold change per iteration
    """
    res = {}
    
    for f in fold_changes:
        ctrl, exp, is_changed = sample_no_ctrl_uniform(
                10000,
                1000,
                f,
                var,
                nctrl=nctrl,
                nexp=nexp)
        res[f] = (
            do_stat_tests(ctrl, exp),
            is_changed
        )
    return res


def simulate_fold_change_vars(
        fold_changes = DEF_FOLD_CHANGES,
        variances = [0.01, 0.05, 0.1, 0.2],
        nctrl=3,
        nexp=3):
    """Create simulated datasets for fold changes AND variances"""

    res = {}
    
    for f in fold_changes:
        for v in variances:
            ctrl, exp, is_changed = sample_no_ctrl_uniform(
                    10000,
                    1000,
                    f,
                    v,
                    nctrl=nctrl,
                    nexp=nexp)
            res[(f,v)] = (
                do_stat_tests(ctrl, exp),
                is_changed
            )

    return res


def simulate_with_gamma(n, alpha=3, beta=0.1, nctrl=3, use_var=np.random.normal):
    variances = 1 / np.random.gamma(alpha, 1./beta, n)
    noise = np.array([
        use_var(loc=0, scale=v**0.5, size=nctrl)
        for v in variances
    ])

    return noise
