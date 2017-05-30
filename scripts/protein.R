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
    # Requires: columns
    # C1, C2, C3, E1, E2, E3, accession_number, meanE, meanC

    print('Obtaining medians')

    # acc_nums <- if (use_isoform) 'accession_number' else 'accession_number_noiso'

    # Obtain medians of each channel
    channel_medians <- aggregate(pepdf[,c('C1', 'C2', 'C3', 'E1', 'E2', 'E3')], by=pepdf['accession_number'], FUN=median)
    # Obtain summary statistics considering each channel + obs separately
    overall_medians <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data=data.frame(
            'accession_number' = rep(pepdf$accession_number, 3),
            'C' = c(pepdf$C1, pepdf$C2, pepdf$C3),
            'E' = c(pepdf$E1, pepdf$E2, pepdf$E3)),
        FUN=function(x) c(med=median(x), mn=mean(x), n=length(x), sd=sd(x), cv=sd(x)/mean(x)))))
    # Obtain fold changes from each peptide
    fold_change = data.frame(
        accession_number = pepdf$accession_number,
        # TODO change this
        # Use MEAN OF RATIOS instead of RATIO OF MEANS
        # Don't use VSN (!)
        fold_change = pepdf$meanE - pepdf$meanC)
    fold_change_cv <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data = fold_change,
        FUN=function(x) c(med=median(x), max=max(abs(x)))
    )))
    # Use bayesian variance estimates
    bayes_var <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data = pepdf[,c('accession_number', 'bayesSDC', 'bayesSDE')],
        # Calculates variance assuming that peptides are cond. independent
        FUN=function(x) c(bayesVar=sqrt(sum(x^2) / length(x)))
    )))
    p_vals <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data = pepdf[,c('accession_number', 'pVal', 'BH')],
        FUN=function(x) c(min=min(x), med=median(x))
    )))

    out <- merge(channel_medians, overall_medians, by='accession_number')
    out <- merge(out, fold_change_cv, by='accession_number')
    out <- merge(out, bayes_var, by='accession_number')
    out <- merge(out, p_vals, by='accession_number')
    # return(map_to_psd(out))
    return(out)
}

