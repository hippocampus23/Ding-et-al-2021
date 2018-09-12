import numpy as np
from simulation import simulate_fold_change_range, simulate_number_experiments, \
                       simulate_fdr_fc_range, simulate_variance_range, simulate_size_dataset, \
                       simulate_number_channels_imbalanced
from formatting import write_result_dict_to_df
from stat_tests import TESTS

# script to run on Google cloud
# generates data for figures 3B, 3D, 5B, 6B, 4E, 4J
# figures can be created using the plots$make_figures() in figures.R
# TODO: 5D, supp

# make data for figures 3B, 3D and with variance trend
FOLD_CHANGES = np.arange(0.1, 1.1, 0.1)


def data_3B():
    res_3B = simulate_fold_change_range(FOLD_CHANGES, "uniform")
    write_result_dict_to_df(res_3B, TESTS.keys(), filename="../data_simulated/3B.csv", ana_type="roc")
    print "saved 3B.csv"


def data_3D():
    res_3D = simulate_fold_change_range(FOLD_CHANGES, "gamma")
    write_result_dict_to_df(res_3D, TESTS.keys(), filename="../data_simulated/3D.csv", ana_type="roc")
    print "saved 3D.csv"


def data_3B_trend():
    res_3B_trend = simulate_fold_change_range(FOLD_CHANGES, "trend")
    write_result_dict_to_df(res_3B_trend, TESTS.keys(), filename="../data_simulated/3B_trend.csv", ana_type="roc")
    print "saved 3B_trend.csv"


def data_4E_trend():
    res_4E = simulate_fdr_fc_range(FOLD_CHANGES, "uniform")
    write_result_dict_to_df(res_4E, TESTS.keys(), filename="../data_simulated/4E.csv", ana_type="fdr")
    print "saved 4E.csv"


def data_4J():
    res_4J = simulate_fdr_fc_range(FOLD_CHANGES, "gamma")
    write_result_dict_to_df(res_4J, TESTS.keys(), filename="../data_simulated/4J.csv", ana_type="fdr")
    print "saved 4J.csv"


def data_4E_trend():
    res_4E_trend = simulate_fdr_fc_range(FOLD_CHANGES, "trend")
    write_result_dict_to_df(res_4E_trend, TESTS.keys(), filename="../data_simulated/4E_trend.csv", ana_type="fdr")
    print "saved 4E_trend.csv"


def data_5B():
    res_5B = simulate_number_experiments(range(2,11), var_type="gamma")
    write_result_dict_to_df(res_5B, TESTS.keys(), filename="../data_simulated/5B.csv", ana_type="roc")
    print "saved 5B.csv"


def data_5D():
    res_5D = simulate_number_channels_imbalanced(trials=[(5,5), (4,6), (3,7), (2,8),(1,9)], var_type="gamma")
    write_result_dict_to_df(res_5D, TESTS.keys(), filename="../data_simulated/5D.csv", ana_type="roc")
    print "saved 5D.csv"


def data_6B():
    res_6B = simulate_variance_range(vars=[0.02, 0.06, 0.18], betas=[0.05, 0.1, 0.2])
    write_result_dict_to_df(res_6B, TESTS.keys(), filename="../data_simulated/6B.csv", ana_type="roc")
    print "saved 6B.csv"


def data_S1B():
    res_S1B = simulate_size_dataset(num_peps=[1000, 10000, 100000], perc_changed=[0.04, 0.1, 0.25], var_type="gamma")
    write_result_dict_to_df(res_S1B, TESTS.keys(), filename="../data_simulated/S1B.csv", ana_type="roc")
    print "saved S1B.csv"