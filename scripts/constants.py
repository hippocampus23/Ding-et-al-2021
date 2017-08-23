from collections import OrderedDict

# Color mapping from R for peptides                                               
COLORS = {'Moderated T (1-sample)': "#b79f00",
          'Moderated T (2-sample)': "#00ba38",                                  
          'Absolute fold change': "#00bfc4",                                  
          'CyberT': "#f8766d",                                    
          't-test (1-sample)': "#f564e3",                                       
          't-test (2-sample)': "#619cff",                                       
}       

# Map header columns to well-formatted legend labels
# NOTE order is important here
LABEL_MAPPING = OrderedDict([
        # Peptide labels
        ('cyberT', 'CyberT'),
        ('modT', 'Moderated T (1-sample)'),
        ('fold change', 'Absolute fold change'),
        ('fold_change', 'Absolute fold change'),
        ('t-test', 't-test (2-sample)'),
        ('t-test (1-sample)', 't-test (1-sample)'),
        ('modT (2-sample)', 'Moderated T (2-sample)'),
        # Protein labels
        ('cyberT_PVal_med', 'Median intensity CyberT'),
        ('modT_PVal_med', 'Median intensity Moderated T'),
        ('fold_change_med', 'Median intensity absolute fold change'),
        ('ttest_PVal_med', 'Median intensity t-test'),
        ('wls', 'Weighted least squares'),
        ('cyberT_bypep', 'CyberT by peptide'),
        ('ttest_bypep', 't-test by peptide'),
])
