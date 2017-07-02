import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HEADER_TEMPLATE =  "| {:>20} | {:^16} | {:^16s} | {:^16s} |"
DEF_FIELDS = ("AUROC", "AUPRC", "pAUC")
BRK = "|" + ("-" * 22) + "|" + ((("-" * 18) +"|") *3)


def subtable(res, labels):
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
        lines += subtable(res[i], labels)

    lines += [BRK]

    print "\n".join(lines)
    return lines


def write_result_dict_to_df(res, labels):
    """ Converts result dictionary to pandas dataframe
    NOTE: every df in res MUST have same length """
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
    return out
    

def _boxplot(res, labels, ax):
    pass

def boxplot_result_dictionary(res, labels, title="", desc="{}"):
    pass
