# Proteomics Simulations

Comparing the performance of several statistical tests for peptide-level analysis using simulated MS data.

**Note**: All simulated and empirical data is log transformed with base 2.


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
make 1EFG.eps
```
Replace "1EFG.eps" with you desired target. Note, when regenerating a figure
that requires multiple rounds of simulations, the corresponding data file 
must already be in the directory "data_simulated". Either make it or use the one that is provided
in the git repository.

**Note**: It is recommended to create the data for the figures that require many rounds of simulations 
in a VM e.g. on the Google Cloud Platform, as they are memory intensive and have a long runtime. 
Alternatively, change N_RUNS in scripts/constants.py from 500 to a smaller number.
Creating the figures should work just fine on a Laptop. 

Screen is a useful tool to keep the simulations running even if you close the ssh connection to the VM.
