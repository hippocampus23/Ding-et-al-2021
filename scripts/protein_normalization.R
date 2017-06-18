source('scripts/bayesreg.R')
source('scripts/multiplot.R')
library(ggplot2)
library(plyr)

# Note on isoforms
#   Isoforms are preserved in all peptide-level manipulations (ex, normalization)
#       <This may change in future, just for more easy comparison>
#       <Anything on protein level could drop isoforms?>
#   Isoforms are DROPPED when comparing to modT results


# Other possible ideas
# http://www.nature.com/nbt/journal/v28/n1/full/nbt.1592.html

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

#################
# PREPROCESSING #
#################

parse_psm_report <- function(peptidesdf) {
    # Maps the raw psm report to the appropriate proteins
    psm = read.table("data/Neurogranin_KD_PSM_021517.ssv", sep=';', quote='', header=TRUE)

    reduced_psm <- log2(psm[, c('TMT_126', 'TMT_127', 'TMT_128', 'TMT_129', 'TMT_130', 'TMT_131')])
    colnames(reduced_psm) <- c('E1', 'E2', 'E3', 'C1', 'C2', 'C3')
    reduced_psm$sequence <- psm$sequence
    # Drop any rows with NAs in intensities
    reduced_psm <- reduced_psm[complete.cases(reduced_psm),]

    reduced_peptides <- peptidesdf[c('id', 'accession_number','entry_name', 'species', 'is_psd')]
    reduced_peptides$sequence <- extract_accession(as.character(reduced_peptides$id))
    
    out <-merge(reduced_psm, reduced_peptides, by=c('sequence'))
    out$id <- make.names(out$id, unique=TRUE)
    out$accession_number_noiso <- strip_accession_isoform(out$accession_number)
    return(out)
}


map_peptides_to_proteins <- function() {
    pep = read.csv(file = "data/NgKD_Protein_Peptide_Comparison_091416.csv")
    c.input <- read.delim(file = "data/NgKD_CyberT_Peptide_Input_061316.txt")
    c.Key <- read.delim(file = "data/NgKD_CyberT_Peptide_Input_Key_061316.txt")

    # Normalize on a per-channel basis to the sum of intensities in that column
    # To make the numbers more interpretable I have normalized by the means
    # This should be fine as long as there are no NAs
    # NOTE Can also simply use the VSN applied raw intensities from pep instead!
    peptides <- c.input[,2:7]
    # peptides <- sweep(c.input[,2:7], 2, colMeans(c.input[,2:7]), '/')
    peptides$id <- c.input$id

    print('Merging peptides to acc_num')
    # Merge peptides to protein names and acc numbers
    # Then merge accession numbers to peptide table to filter to mouse proteins
    keys = merge(x=c.Key, y=unique(pep[which(pep$species == 'MOUSE'), c('accession_number', 'species')]), by='accession_number')
    peptides <- merge(x=peptides, y=keys, by='id')
    peptides$accession_number <- as.character(peptides$accession_number)
    peptides$accession_number_noiso <- strip_accession_isoform(peptides$accession_number)

    #Find duplicate peptides if any
    dupes <- duplicated(peptides$id)
    if (any(dupes)) {
        warning(paste(str(sum(dupes)), ' duplicate ids were found\nThey will be discarded'))
        peptides <- peptides[!dupes,]
    }
    return(map_to_psd(peptides))
}

map_phosphosites_to_proteins <- function() {
    phosphosites <- read.table("data/NgKD_CyberT_Phosphosite_input_012517.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
    
    phosphosites$accession_number <- extract_accession(phosphosites$id)
    phosphosites$accession_number_noiso <- strip_accession_isoform(phosphosites$accession_number)
    return(map_to_psd(phosphosites))
}


find_protein_medians <- function(pepdf, use_isoform=TRUE) {
    # Finds the median peptide of each protein (as judged by acc num) for each channel
    # AND the overall median for the control and experimental conditions

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
        fold_change = pepdf$meanE / pepdf$meanC)
    fold_change_cv <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data = fold_change,
        FUN=function(x) c(mn=mean(x), cv=sd(x)/mean(x))
    )))
    # Use bayesian variance estimates
    bayes_var <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data = pepdf[,c('accession_number', 'bayesSDC', 'bayesSDE')],
        # Calculates variance assuming that peptides are cond. independent
        FUN=function(x) c(bayesVar=sqrt(sum(x^2) / length(x)))
    )))

    out <- merge(channel_medians, overall_medians, by='accession_number')
    out <- merge(out, fold_change_cv, by='accession_number')
    out <- merge(out, bayes_var, by='accession_number')
    out$accession_number_noiso <- strip_accession_isoform(out$accession_number)
    return(map_to_psd(out))
}


normalize_peptide_to_protein <- function(pepdf, protdf) {
    # DON'T USE THIS EXCEPT FOR PHOSPHOPROTEOME DATA
    # Uses the MEDIAN protein intensity in each condition to normalize
    #   the peptide intensity measurements in preparation for 
    # TODO determine if we should normalize to ONE quantity (med in condition)
    #   or normalize to each channel

    peptide_merged <- merge(
        pepdf[,c('C1', 'C2', 'C3', 'E1', 'E2', 'E3', 'accession_number', 'accession_number_noiso', 'id')],
        protdf[,c('C.med', 'E.med', 'accession_number', 'is_psd')],
        by='accession_number')

    return(data.frame(
        accession_number=peptide_merged$accession_number,
        accession_number_noiso=peptide_merged$accession_number_noiso, 
        id=peptide_merged$id,
        'C1'=peptide_merged$C1/peptide_merged$C.med,
        'C2'=peptide_merged$C2/peptide_merged$C.med,
        'C3'=peptide_merged$C3/peptide_merged$C.med,
        'E1'=peptide_merged$E1/peptide_merged$E.med,
        'E2'=peptide_merged$E2/peptide_merged$E.med,
        'E3'=peptide_merged$E3/peptide_merged$E.med,
        is_psd=peptide_merged$is_psd
    )) 
}


do_cyberT <- function(pepdf, doVsn=FALSE, bayesInt=FALSE) {
    # Runs the cyberT t-test (with optional VSN normalization)
    # Uses default settings for cyberT
    # NOTE will throw error if duplicate names are detected

    if (sum(duplicated(pepdf$id)) > 0) {
        stop('Not all peptides have unique IDs!')
        return()
    }
    
    rownames(pepdf) <- pepdf$id
    input <- pepdf[,c('C1', 'C2', 'C3', 'E1', 'E2', 'E3')]
    if (doVsn) {
        input <- runVsn(input)
        input <- input - min(input, na.rm = TRUE) + 1.0
    }

    out <- bayesT(input, 3, 3, doMulttest=TRUE)
    out$fold_change <- out$meanE / out$meanC
    out <- merge(out, pepdf[, c('accession_number', 'accession_number_noiso', 'is_psd')], by=0)
    out$id <- out$Row.names
    out$Row.names <- NULL
    return(out)
}


peptide_to_prot_signif <- function(pepdf, protdf, sig=0.05) {
    # Calculate signif thresholds by peptide level p-vals into protein signif calculations
    # Uses Benjamini Hochberg adjusted p-vals, 'BH' column
    # WARNING: sig parameter does not work at this time because of ddply issues
    # manually change the significance threshold inside this function
    protein = ddply(pepdf, .(accession_number), summarize, count=length(accession_number), signif=sum(BH <= 0.05), signif_freq=sum(BH <= 0.05)/length(accession_number))
    protein.final = merge(protdf, protein, by='accession_number')
    protein.final$fold_change <- protein.final$E.med / protein.final$C.med
    return(protein.final)
}



protein_cutoffs <- function(protdf) {
    modt_total <- read.table("data/NgKD_proteome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE,quote = "")

    alpha = 0.05
    signif <- modt_total$Accession.Number[modt_total$adj.P.Val <= alpha]
    step = 0.1
    cutoffs <- seq(0.0, 1.0, step)

    find_overlap <- function(acc_nums) {length(intersect(acc_nums, signif))} 
    avg_pvals <- function(acc_nums) {
        mean(modt_total$adj.P.Val[modt_total$Accession.Number %in% acc_nums])
    }

    # Make into a data frame
    df = do.call(
        rbind,
        lapply(cutoffs, function(cutoff) {
                cum = protdf[protdf$signif_freq >= cutoff,]
                sub = cum[cum$signif_freq <= cutoff + step,]
                c(cutoff=cutoff,
                  cum_number=length(cum$accession_number),
                  cum_overlap=find_overlap(cum$accession_number),
                  cum_avg_pval=avg_pvals(cum$accession_number),
                  cum_avg_fold_change=mean(cum$fold_change),
                  sub_number=length(sub$accession_number),
                  sub_overlap=find_overlap(sub$accession_number),
                  sub_avg_pval=avg_pvals(sub$accession_number),
                  sub_avg_fold_change=mean(sub$fold_change)
                  )
            }
        )
    )
    return(as.data.frame(df, stringsAsFactors=FALSE))
}


### Calculate variance for each protein ###
### Playing around with pipes ###
protein_variance <- function(psm_data) {
    
}

#############################
# ACTUAL SCRIPT STARTS HERE #
#############################

# peptides <- map_peptides_to_proteins()
psm <- parse_psm_report(peptides)
# peptides_cyberT <- do_cyberT(peptides)
psm_cyberT <- do_cyberT(psm)
# proteins <- find_protein_medians(peptides_cyberT)
# proteins_final <- peptide_to_prot_signif(peptides_cyberT, proteins)
proteins_psm <- find_protein_medians(psm_cyberT)
proteins_psm_final <- peptide_to_prot_signif(psm_cyberT, proteins_psm)
# phosphos <- map_phosphosites_to_proteins()
# phosphos_cyberT <- do_cyberT(phosphos)
# phosphos_norm <- normalize_peptide_to_protein(phosphos, proteins)
# phosphos_norm_cyberT <- do_cyberT(phosphos_norm)


extract_modT_fp <- function(pepdf, proteins_final) {
    bad_proteins <- proteins_final[proteins_final$adj.P.Val <= 0.05 & proteins_final$signif_frac <= 0.01,]
}

##############################
# FIGURES AND VISUALIZATIONS #
##############################


compare_psm_peptide <- function(protdf, psm_protdf, cutoff=20, geom='box') {
    # Plots boxplot, scatterplot, or contour plot of how counts changed for psm vs peptide report
    combined <- merge(protdf, psm_protdf, by='accession_number')
    combined <- combined[combined$count.x <= cutoff,]

    if (geom == 'box') {
        p = qplot(factor(count.x), count.y, data=combined, geom='boxplot')
        p = p + labs(title='Closeup of low-count peptides (original vs PSM)')   
    } else if (geom == 'scatter') {
        p = qplot(count.x, count.y, data=combined, alpha=0.01, stroke=0.5, color='blue')
        p = p + labs(title='Scatterplot of peptide counts (original vs PSM')
    }

    p = p + geom_abline() + labs(x='Peptide count (original)', y='Peptide count (PSM)')
    return(p)
}

threshold_cvs <- function(protdf, cutoffs = c(0.005, 0.01, 0.015, 0.02, 0.025), band=FALSE) {
    # If banded, then we have to actually create cutoff bands
    # TODO
    if (band) {
        cutoffs_df = cbind(c(0, cutoffs)[1:length(cutoffs)], cutoffs)
    } else {
        cutoffs_df = cbind(0, cutoffs)
    }

    # Create one volcano plot for each cutoff
    # Set consistent xlims and ylims beforehand!
    x_min = min(log2(protdf$fold_change.mn))
    x_max = max(log2(protdf$fold_change.mn)) 
    volcanoplots <- as.list(apply(cutoffs_df, 1, function(cutoff) {
        p = qplot(log2(fold_change.mn), 
                  signif_freq,
                  data=proteins_final[proteins_final$fold_change.cv > cutoff[1] & 
                                      proteins_final$fold_change.cv <= cutoff[2],],
                  alpha=count, 
                  stroke=0)
        p = p + labs(x='Log2 fold change', y='Proportion of peptides significant', title=paste('CV at or below', cutoff[2]))
        if (band) {
            p = p + labs(title=paste("CV between", cutoff[1], "and", cutoff[2]))
        }
        p = p + xlim(x_min, x_max)
        p = p + theme(legend.position = 'none')
        return(p)
    }))

    # Create p-value comparison for each cutoff
    # TODO this is a pretty sketchy dependence on the plots function.
    merged_df <- plots$merge_proteins(protdf[,
                     c('accession_number', 'signif_freq', 'fold_change.cv', 'count')])
    pvals <- as.list(apply(cutoffs_df, 1, function(cutoff) {
        p = qplot(
            -log10(adj.P.Val),
            signif_freq,
            data = merged_df[merged_df$fold_change.cv > cutoff[1] & 
                             merged_df$fold_change.cv <= cutoff[2], ],
            alpha = count,
            stroke=0
        )
        p = p + labs(x='ModT adjusted p value', y='Proportion of peptides significant')
        p = p + theme(legend.position = 'none')
        return (p)
    }))

    multiplot(plotlist=append(volcanoplots, pvals), cols=2)
}

setup_plots <- function(protdf, phosphodf) {
    # Read in modT comparison files
    modt_total <- read.table("data/NgKD_proteome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE,quote = "")
    modt_phospho <- read.table("data/ngkd_joint_complete.csv", sep=',', header=TRUE, stringsAsFactors=FALSE, quote="\"")

    plots <- list()

    plots$merge_proteins <- function(proteins) {
        return(merge(
            modt_total[,c('Accession.Number', 'P.Value', 'adj.P.Val', 'Average.Log2.Expression', 'geneSymbol')], 
            # proteins[,c('accession_number', 'fold_change', 'signif_freq', 'count', 'is_psd')],
            proteins,
            by.x='Accession.Number',
            by.y='accession_number'
        ))
    }
    plots$merge_phospho <- function(phospho) {
        return(merge(
            modt_phospho[,c('accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA', 'P.Value', 'adj.P.Val', 'Average.Log2.Expression', 'Norm.P.Value', 'Norm.adj.P.Val', 'Norm.logFC','Gene.Symbol')],
            phospho[,c('id', 'pVal', 'BH', 'fold_change', 'is_psd')],
            by.x='accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA',
            by.y='id'
        ))
    }

    # Merge on accession_number
    plots$proteins_merged <- plots$merge_proteins(protdf)

    # Phosphopeptides
    plots$phospho_merged <- plots$merge_phospho(phosphodf)

    plots$plot_protein_pvals <- function(protdf=NULL) {
        if (protdf == NULL) {
            protdf = plots$proteins_merged
        } else {
            prodf = plots$merge_proteins(protdf)
        }
        # Compare the calculated p-values for the modT vs the cyberT proteins
        # Keeps all isoform information
        qplot(
            -log10(adj.P.Val),
            signif_freq,
            data = protdf,
            alpha = count,
            color= is_psd
        )
    }

    plots$plot_protein_foldchange <- function(protdf=plots$proteins_merged) {
        # Compare the calculated fold changes for the modT vs the cyberT proteins
        # Keeps all isoform information 
        qplot(
            log2(fold_change),
            Average.Log2.Expression,
            data = protdf,
            alpha = count,
            color = is_psd
        )
    }

    plots$plot_phospho_pvals <- function(phosphodf=plots$phospho_merged) {
        qplot(
            -log2(Norm.adj.P.Val),
            -log2(BH),
            data=phosphodf,
            alpha=0.0001,
            color=is_psd
        )
    }

    plots$plot_phospho_foldchange <- function(phosphodf=plots$phospho_merged) {
        qplot(
              log2(fold_change),
              Norm.logFC,
              data=phosphodf,
              alpha=0.1,
              color=is_psd
        )
    }
    
    return(plots)
}

# Example of how to compare the protein_final dataframes
# p = qplot(C.mn, count, data=proteins_psm_final, alpha=count, stroke=0)+ geom_density2d(aes(alpha=256))
# p = p + labs(x='count', y='Standard deviation of control intensities', title='Variance vs mean and peptide count (Controls, PSM)')

# plots <- setup_plots(proteins_final, phosphos_norm_cyberT)
