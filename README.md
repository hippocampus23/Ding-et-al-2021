# Proteomics Simulations

Comparing the performance of several statistical tests for peptide-level analysis using simulated MS data.

**Note**: All simulated and empirical data is log transformed with base 2.

## Figures
https://www.dropbox.com/sh/o1ydax3y9tql6cy/AACLABgJFt9zDSzhGd4rR7U9a?dl=0 


## Requirements

**Python (2.7+)**
- python-tk
- NumPy (1.11.1+)
- Pandas (0.8.1+)
- rpy2 (2.8.5+)
- SciPy (0.18.0+)
- h5py (2.7.0+)
- Scikit-Learn (0.17.1+)
- Matplotlib (1.5.2+)

**R (3.4+)**
- Bioconductor
- limma
- multtest
- ggplot2

## Setting up all requirements and tools needed on Linux (Debian)
### Git
Open a terminal and type the following commands to install git and clone the repository:
```bash
sudo apt-get install git
git clone https://github.com/hmsch/proteomics.git
```

### Python
*Python 2 is usually preinstalled, else install it.*

Open a terminal and type the following commands to install the required packages and
the package manager pip:
```bash
sudo apt-get install python-tk
sudo apt install python-pip
sudo pip install -r proteomcis/scripts/requirements.txt
```
Alternatively, use Conda.

### R
Open a terminal and type the following command to install R and start it:
```bash
sudo apt-get install r-base r-base-dev
R
```
Then install the required R packages:
``` R
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
biocLite("multtest")
install.packages("ggplot2")
quit()
```

## Regenerating the Data and Figures
Go into the proteomics folder: ```cd proteomics```

And explore the rules in the Makefile: ```cat Makefile```
Any data you create will be saved in the directory "data_simulated" and
figures will be saved in the directory "figures". So, feel free to copy the
existing data and figures elsewhere if you do not want to overwrite them.
To regenerate a figure or data, simply type e.g.
``` bash
make 1FGH.eps
```
Replace "1FGH.eps" with you desired target. Note, when regenerating a figure
that requires multiple rounds of simulations, the corresponding data file 
must already be in the directory "data_simulated". Either make it or use the one that is provided
in the git repository.

When making multiple targets at the same time (e.g. all data files), it is usefull to use options like
```make data -j4``` on a machine with multiple cores / CPUs to make use of all of them.

**Note**: It is recommended to create the data for the figures that require many rounds of simulations 
in a VM e.g. on the Google Cloud Platform, as they are memory intensive and have a long runtime. 
Alternatively, change N_RUNS in scripts/constants.py from 500 to a smaller number.
Creating the figures should work just fine on a Laptop. 

Screen is a useful tool to keep the simulations running even if you close the ssh connection to the VM.

## Running Custom Simulations
The ```sample``` function in ```sampler.py``` can be used to generate simulated data for a single MS run. Are variety of parameters including the total number of peptides, fold changes, and the variance distribution can be specified (see details below).

```python
def sample(pep_var_type, num_peps = NUM_PEPTIDES, num_ctrl=NUM_CHANNELS/2, num_exp=NUM_CHANNELS/2,
           num_changed=NUM_CHANGED, fold_change=FOLD_CHANGE, use_var=np.random.normal, alpha=ALPHA,
           beta=BETA, pep_var=PEPTIDE_VAR):
    """
    Simulate a MS data set with experiment and control data

    :param pep_var_type:     type of distribution to use to sample peptide variance, either "uniform", "gamma" or "trend"
    :param num_peps:         number of peptides (i.e. number of rows of the data set)
    :param num_ctrl:         number of control channels (i.e. columns of data)
    :param num_exp:          number of experiment channels (i.e. columns of data)
    :param num_changed:      number of peptides that get a fold change
    :param fold_change:      fold change to add to perturbed peptides
    :param use_var:          function that takes arguments loc, scale and size for generating differences
                             between peptides
    :param alpha:            parameter for inverse gamma distribution
    :param beta              parameter for inverse gamma distribution
    :param pep_var           uniform variance if pep_var_type is "uniform"
    :return:
                 ctrl:       pandas.DataFrame with log 2 control data
                 exp:        pandas.DataFrame with log 2 experiment data
                 is_changed: numpy.ndarray with 0 for no fold change and 1 for fold change
    """
```

## Calculating p-values
A number of statistical tests can be used to calculate the p-values for differential expression of the peptides for either simulated of empirical data. ```stat_tests.py``` contains functions for **ranking by fold change**, standard **t-tests** as well as the **regularized t-test** and **moderated t-test**. All tests take a pandas DataFrame with the control data and one with the experiment data respectively and return a numpy array of p-values. Some tests have a one-sample and two-sample option and the moderated t-test can be run with or without a mean-variance trend. Using ```do_all_tests(ctrl, exp)``` all provided tests can be run.

When using empirical data, the data can be imported into a pandas DataFrame from a csv file using
```python
ctrl = pandas.read_csv("./path/to/the/data/file.csv", sep=";", usecols=["C1", "C2", "C3"])
exp = pandas.read_csv("./path/to/the/data/file.csv", sep=";", usecols=["E1", "E2", "E3"])
```
using ```usecols``` to specify the columns corresponding to the case.

After running a test, the p-values can be BH-adjusted using
```python
_, pvals_adj, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(pvals, alpha=0.05, method='fdr_bh')
```
if desired. Alternatively, when using simulated data the roc and prc scores can be calculated or a power-analysis can be performed using functions from ```eval_results.py```.

The prc, roc or power-analysis results can be saved to a csv file using ```write_result_dict_to_df(...)``` from ```format_utils.py```. There are also various plotting functions available to visualize the results.

The p-values from ```do_all_tests(...)``` are already in a pandas DataFrame, otherwise cast them to a pandas DataFrame to save them to a csv file using
```python 
pvals.to_csv("./path/to/file/name.csv", sep=";")
```

