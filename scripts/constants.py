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
    ("Abs. FC",                        "#C14242"),
    ("t-test\n(1-sample)",            "#FDA51F"),
    ("t-test\n(2-sample)",            "#0882FB"),
    ("ModT\n(1-sample)",              "#0BFD33"),
    ("ModT\n(2-sample)",              "#10B571"),
    ("ModT\n(1-s., trend)",           "#F826FB"),
    ("ModT\n(2-s., trend)",           "#4E03EC"),
    ("RegT",                         "#20E6FC"),
])

# Map header columns to well-formatted legend labels
LABEL_MAPPING = OrderedDict([
    ("fold change",  "Abs. FC"),
    ("t-test-1",     "t-test\n(1-sample)"),
    ("t-test-2",     "t-test\n(2-sample)"),
    ("modT-1",       "ModT\n(1-sample)"),
    ("modT-2",       "ModT\n(2-sample)"),
    ("modT-1 trend", "ModT\n(1-s., trend)"),
    ("modT-2 trend", "ModT\n(2-s., trend)"),
    ("RegT",         "RegT")
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
