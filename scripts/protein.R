extract_accession <- function(info_col) {
  return(unlist(lapply(info_col, function(x) strsplit(x,"_")[[1]][1])))
}

strip_accession_isoform <- function(accession_numbers) {
  return(unlist(lapply(accession_numbers, function(x) strsplit(x,"-")[[1]][1])))
}

map_to_psd <- function(df, by='accession_number_noiso') {
  # Generates a boolean column in the df that indicates if protein is part of psd
  psd <- read.csv('data/psd.tsv', sep="\t")$Accession
  df$is_psd <- (df[[by]] %in% psd)
  return(df)
}


find_protein_medians <- function(pepdf, use_isoform=TRUE) {
  # Finds the median peptide of each protein (as judged by acc num) for each channel
  # AND the overall median for the control and experimental conditions
  # Requires: fold_change column, C1, C2, C3, E1, E2, E3
  # All columns which end with PVal are retained

  print('Rolling up proteins')

  pValCols = colnames(pepdf)[grep('PVal$', colnames(pepdf))]

  tmp_df <- pepdf[,c('accession_number',
                     pValCols)]

  # Per-peptide pvalue reports
  pval_df <- as.data.frame(as.list(aggregate(
      . ~ accession_number,
      data = tmp_df,
      FUN = function(x) c(med=median(x), min=min(x))
  )))

  fold_change_df <- as.data.frame(as.list(aggregate(
      . ~ accession_number,
      data = pepdf[, c('accession_number', 'fold_change')],
      FUN = function(x) c(med=median(x), max=max(abs(x)))
  )))
  out <- merge(pval_df, fold_change_df, by='accession_number')
  #     out <- merge(out, fold_change_cv, by='accession_number')
  return(out)
}

