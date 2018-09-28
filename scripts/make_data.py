import numpy as np
from simulation import simulate_fold_change_range, simulate_number_experiments, \
                       simulate_fdr_fc_range, simulate_variance_range, simulate_size_dataset, \
                       simulate_number_channels_imbalanced, simulate_FC_SD
from format_utils import write_result_dict_to_df
from stat_tests import TESTS, regT, modT

# regenerate data for the figures that require multiple simulations

FOLD_CHANGES = np.arange(0.1, 1.0, 0.1)


def data_3B():
    res_3B = simulate_fold_change_range(FOLD_CHANGES, "uniform")
    write_result_dict_to_df(res_3B, TESTS.keys(), filename="../data_simulated/3B.csv", ana_type="roc")
    print "saved 3B.csv"


def data_3D():
    res_3D = simulate_fold_change_range(FOLD_CHANGES, "gamma")
    write_result_dict_to_df(res_3D, TESTS.keys(), filename="../data_simulated/3D.csv", ana_type="roc")
    print "saved 3D.csv"


def data_3F():
    res_3F = simulate_fold_change_range(FOLD_CHANGES, "trend")
    write_result_dict_to_df(res_3F, TESTS.keys(), filename="../data_simulated/3F.csv", ana_type="roc")
    print "saved 3F.csv"


def data_4D():
    res_4D = simulate_fdr_fc_range(FOLD_CHANGES, "uniform")
    write_result_dict_to_df(res_4D, TESTS.keys(), filename="../data_simulated/4D.csv", ana_type="fdr")
    print "saved 4D.csv"


def data_4H():
    res_4H = simulate_fdr_fc_range(FOLD_CHANGES, "gamma")
    write_result_dict_to_df(res_4H, TESTS.keys(), filename="../data_simulated/4H.csv", ana_type="fdr")
    print "saved 4H.csv"


def data_4L():
    res_4L = simulate_fdr_fc_range(FOLD_CHANGES, "trend")
    write_result_dict_to_df(res_4L, TESTS.keys(), filename="../data_simulated/4L.csv", ana_type="fdr")
    print "saved 4L.csv"


def data_5B():
    res_5B = simulate_number_experiments(range(2,11), var_type="trend")
    write_result_dict_to_df(res_5B, TESTS.keys(), filename="../data_simulated/5B.csv", ana_type="roc")
    print "saved 5B.csv"


def data_5D():
    res_5D = simulate_number_channels_imbalanced(trials=[(5,5), (4,6), (3,7), (2,8),(1,9)], var_type="trend")
    write_result_dict_to_df(res_5D, TESTS.keys(), filename="../data_simulated/5D.csv", ana_type="roc")
    print "saved 5D.csv"


def data_6C():
    res_6C = simulate_variance_range(vars=[0.02, 0.06, 0.18], betas=[0.05, 0.1, 0.2])
    write_result_dict_to_df(res_6C, TESTS.keys(), filename="../data_simulated/6C.csv", ana_type="roc")
    print "saved 6C.csv"


def data_6B():
    res_6B = simulate_FC_SD("gamma", ["t-test-2"], np.arange(0.1, 1.1, 0.1), 0.75)  #, "modT-2 trend", "RegT"
    res_6B.to_csv("../data_simulated/6B.csv")


def data_S1B():
    res_S1B = simulate_size_dataset(num_peps=[1000, 10000, 100000], perc_changed=[0.04, 0.1, 0.25], var_type="gamma")
    write_result_dict_to_df(res_S1B, TESTS.keys(), filename="../data_simulated/S1B.csv", ana_type="roc")
    print "saved S1B.csv"
