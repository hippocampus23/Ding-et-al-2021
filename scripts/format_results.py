import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from roc import plot_both

HEADER_TEMPLATE =  "| {:>20} | {:^16} | {:^16s} | {:^16s} |"
DEF_FIELDS = ("AUROC", "AUPRC", "pAUC")
BRK = "|" + ("-" * 22) + "|" + ((("-" * 18) +"|") *3)


## Following two functions summarize results dictionary into table ##
def _subtable(res, labels):
    """ Summarizes results with mean and standard deviation for each measure
    """
    mns = np.nanmean(res, axis=0)
    stds =  np.nanstd(res, axis=0)
    lines = []

    for (i, l) in enumerate(labels):
        lines.append(
            "| %20s | %.5f+-%.5f | %.5f+-%.5f | %.5f+-%.5f |" %
            (l, mns[i,0], stds[i,0], mns[i,1], stds[i,1], mns[i,2], stds[i,2])
        )
    
    return lines


def summarize_result_dictionary(res, labels, title = "", desc="{}"):
    """ Prints entire table to summarize results
    """
    lines = [
        title,
        HEADER_TEMPLATE.format(" ", *DEF_FIELDS),
    ]
    
    for i in sorted(res.keys()):
        lines += [BRK, HEADER_TEMPLATE.format(desc.format(i), "", "", ""), BRK]
        lines += _subtable(res[i], labels)

    lines += [BRK]

    print "\n".join(lines)
    return lines


## Writes result dictionary to df. Switch to R for boxplot functions ##
def write_result_dict_to_df(res, labels, filename=None):
    """ Converts result dictionary to pandas dataframe
    NOTE: every df in res MUST have same length (n x |labels| x 3)
    Labels MUST be the same length as the second to last dimenstion of res"""
    
    len_label_dim = np.array([v.shape[-2] for v in res.itervalues()])
    if not np.all(len_label_dim == len(labels)):
        raise ValueError("Length of labels does not match second-to-last \
                input dimension for all dataframes in res")

    skeys = sorted(res.keys())
    sorted_res = [res[i].reshape(res[i].shape[0] * res[i].shape[1], -1) 
                  for i in skeys]
    concat = np.concatenate(sorted_res, axis=0)
    label_col = labels * (concat.shape[0] / len(labels))
    setting = np.repeat(skeys, [res[i].shape[0] * res[i].shape[1]
                                for i in skeys])

    print len(setting)

    out = pd.DataFrame(concat)
    out.columns = ["AUROC", "AUPRC", "pAUROC"]
    out['labels'] = label_col
    out['setting'] = setting

    if filename is not None:
        out.to_csv("../data_simulated/" + filename)
    return out


## Transform protein df results for plotting ##
def plot_protein_result_df(df, is_changed, protein, exclude_cols = ['accession_number'], **kwargs):
    """ Plots protein results on ROC curve
        
    Args: df - output of do_stats_tests_protein
               Each column is a list of pvals
          is_changed - 0/1 list of true labels
          exclude_cols - default ['accession_number'], list of cols to exclude
          kwargs - passed to plotting code
            Include title, is_pval
    """
    # Rollup changes to the protein level
    changed_df = pd.DataFrame({"p":protein, "ic":is_changed})
    ic_final = changed_df.groupby("p").first()['ic']

    df = df.drop(exclude_cols, axis=1, inplace=False)
    labels = df.columns
    plot_both(ic_final, [df[c] for c in df.columns], labels, **kwargs)
