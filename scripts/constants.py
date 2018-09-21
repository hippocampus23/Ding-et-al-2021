from collections import OrderedDict

# for multiple simulations
N_RUNS = 500

# constants for simulation
MEAN_PEPTIDE_INTENSITY = 16.0
OVERALL_VAR = 4.0
NUM_PEPTIDES = 10000
NUM_CHANGED = 1000
FOLD_CHANGE = 0.5
NUM_CHANNELS = 6
# for uniform variance distribution
PEPTIDE_VAR = 0.06
# for inverse gamma variance distribution
ALPHA = 3
BETA = 0.1

COLORS = OrderedDict([
    ("Absolute fold change",                          "#C14242"),
    ("t-test (1-sample)",                             "#FDA51F"),
    ("t-test (2-sample)",                             "#0882FB"),
    ("CyberT",                                        "#20E6FC"),
    ("Moderated T (1-sample)",                        "#0BFD33"),
    ("Moderated T (2-sample)",                        "#10B571"),
    ("Moderated T with \nvariance trend (1-sample)",  "#F826FB"),
    ("Moderated T with \nvariance trend (2-sample)",  "#4D03EC")
])


# Map header columns to well-formatted legend labels
# NOTE order is important here
LABEL_MAPPING = OrderedDict([
    ("fold change",  "Absolute fold change"),
    ("t-test-1",     "t-test (1-sample)"),
    ("t-test-2",     "t-test (2-sample)"),
    ("modT-1",       "Moderated T (1-sample)"),
    ("modT-2",       "Moderated T (2-sample)"),
    ("modT-1 trend", "Moderated T with \nvariance trend (1-sample)"),
    ("modT-2 trend", "Moderated T with \nvariance trend (2-sample)"),
    ("cyberT",       "Cyber-T")
])


# TODO fit a function instead
# key is doubled intensity upper bound for easier lookup
VAR_COEFFICIENTS = {
    0  : 5.107914,
    21 : 3.899145,
    22 : 20.109211,
    23 : 18.037317,
    24 : 10.833790,
    25 : 3.693882,
    26 : 3.037527,
    27 : 2.231840,
    28 : 1.627041,
    29 : 1.316483,
    30 : 1.143547,
    31 : 0.968020,
    32 : 0.916222,
    33 : 0.838892,
    34 : 0.834427,
    35 : 0.814620,
    36 : 0.780038,
    1  : 0.727722
}

FDR = 0.05
