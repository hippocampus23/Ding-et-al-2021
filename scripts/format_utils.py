import numpy as np
import pandas as pd

def write_result_dict_to_df(res, labels, filename, ana_type):
    """
    Writes result dictionary to df. Switch to R for boxplot functions

    :param res:       dictionary with data for plotting
                      If data_type=="fdr", then dims must be (n x |labels| x 4)
                      else (n x |labels| x 3)
    :param labels:    labels MUST be the same length as the second to last dimenstion of res
    :param filename:  where to save data
    :param ana_type:  either "roc" or "fdr" will determine column titles
    :return:          pandas.DataFrame that can be passed to R
    """

    len_label_dim = np.array([v.shape[-2] for v in res.itervalues()])
    if not np.all(len_label_dim == len(labels)):
        raise ValueError("Length of labels does not match second-to-last \
                input dimension for all dataframes in res")

    len_label_dim = np.array([v.shape[-2] for v in res.itervalues()])
    if not np.all(len_label_dim == len(labels)):
        raise ValueError("Length of labels does not match second-to-last \
                input dimension for all dataframes in res")

    skeys = sorted(res.keys())
    sorted_res = [res[i].reshape(res[i].shape[0] * res[i].shape[1], -1)
                  for i in skeys]
    concat = np.concatenate(sorted_res, axis=0)
    label_col = labels * (concat.shape[0] / len(labels))
    temp = [res[i].shape[0] * res[i].shape[1] for i in skeys]
    setting = np.repeat(skeys, temp)

    print len(setting)

    out = pd.DataFrame(concat)
    if ana_type== "fdr":
        out.columns = ["FP_raw", "TP_raw", "FP_adj", "TP_adj"]
    elif ana_type== "roc":
        out.columns = ["AUROC", "AUPRC", "pAUROC"]

    out['labels'] = label_col
    out['setting'] = setting

    out.to_csv(filename)

    return out
