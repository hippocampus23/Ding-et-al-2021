source('scripts/bayesreg.R')
library(ggplot2)

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

#################
# PREPROCESSING #
#################

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
    return(peptides)
}

map_phosphosites_to_proteins <- function() {
    phosphosites <- read.table("data/NgKD_CyberT_Phosphosite_input_012517.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
    
    phosphosites$accession_number <- extract_accession(phosphosites$id)
    phosphosites$accession_number_noiso <- strip_accession_isoform(phosphosites$accession_number)
    return(phosphosites)
}


find_protein_medians <- function(pepdf, use_isoform=TRUE) {
    # Finds the median peptide of each protein (as judged by acc num) for each channel
    # AND the overall median for the control and experimental conditions

    print('Obtaining medians')

    # acc_nums <- if (use_isoform) 'accession_number' else 'accession_number_noiso'

    channel_medians <- aggregate(pepdf[,c('C1', 'C2', 'C3', 'E1', 'E2', 'E3')], by=pepdf['accession_number'], FUN=median)
    overall_medians <- as.data.frame(as.list(aggregate(
        . ~ accession_number,
        data=data.frame(
            'accession_number' = rep(pepdf$accession_number, 3),
            'C' = c(pepdf$C1, pepdf$C2, pepdf$C3),
            'E' = c(pepdf$E1, pepdf$E2, pepdf$E3)),
        FUN=function(x) c(med=median(x), mn=mean(x), n=length(x), sd=sd(x)))))

    return(merge(channel_medians, overall_medians, by='accession_number'))
}


normalize_peptide_to_protein <- function(pepdf, protdf) {
    # DON'T USE THIS EXCEPT FOR PHOSPHOPROTEOME DATA
    # Uses the MEDIAN protein intensity in each condition to normalize
    #   the peptide intensity measurements in preparation for 
    # TODO determine if we should normalize to ONE quantity (med in condition)
    #   or normalize to each channel

    peptide_merged <- merge(
        pepdf[,c('C1', 'C2', 'C3', 'E1', 'E2', 'E3', 'accession_number', 'id')],
        protdf[,c('C.med', 'E.med', 'accession_number')],
        by='accession_number')

    return(data.frame(
        accession_number=peptide_merged$accession_number,
        id=peptide_merged$id,
        'C1.norm'=peptide_merged$C1/peptide_merged$C.med,
        'C2.norm'=peptide_merged$C2/peptide_merged$C.med,
        'C3.norm'=peptide_merged$C3/peptide_merged$C.med,
        'E1.norm'=peptide_merged$E1/peptide_merged$E.med,
        'E2.norm'=peptide_merged$E2/peptide_merged$E.med,
        'E3.norm'=peptide_merged$E3/peptide_merged$E.med
    )) 
}


do_cyberT <- function(pepdf, doVsn=TRUE) {
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
    }

    out <- bayesT(input, 3, 3, doMulttest=TRUE)
    out <- merge(out, pepdf[, c('accession_number', 'accession_number_noiso')], by=0)
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
    cutoffs <- seq(0.1, 1.0, 0.1)

    find_overlap <- function(acc_nums) {length(intersect(acc_nums, signif))} 
    avg_pvals <- function(acc_nums) {
        mean(modt_total$adj.P.Val[modt_total$Accession.Number %in% acc_nums])
    }

    # Make into a data frame
    df = do.call(
        rbind.data.frame,
        lapply(cutoffs, function(cutoff) {
                subset = protdf$accession_number[protdf$signif_freq >= cutoff]
                c(number=length(subset),
                  overlap=find_overlap(subset),
                  avg_pval=avg_pvals(subset))
            }
        )
    )
    return(df)
}


############################
# ACTUAL STUFF STARTS HERE #
############################

# peptides <- map_peptides_to_proteins()
# proteins <- find_protein_medians(peptides)
# peptides_cyberT <- do_cyberT(peptides)
# proteins_final <- peptide_to_prot_signif(peptides_cyberT, proteins)
# phosphos <- map_phosphosites_to_proteins()
# phosphos_norm <- normalize_peptide_to_protein(phosphos, proteins)
# phosphos_cyberT <- do_cyberT(phosphos)


##############################
# FIGURES AND VISUALIZATIONS #
##############################


setup_plots <- function(protdf, phosphodf) {
    # Read in modT comparison files
    modt_total <- read.table("data/NgKD_proteome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE,quote = "")
    modt_phospho <- read.table("data/ngkd_joint_complete.csv", sep=',', header=TRUE, stringsAsFactors=FALSE, quote="\"")

    plots <- list()

    # Merge on accession_number
    plots$proteins_merged <- merge(
        modt_total[,c('Accession.Number', 'P.Value', 'adj.P.Val', 'Average.Log2.Expression')], 
        protdf[,c('accession_number', 'fold_change', 'signif_freq', 'count')],
        by.x='Accession.Number',
        by.y='accession_number'
    )

    # Phosphopeptides
    plots$phospho_merged <- merge(
        modt_phospho[,c('accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA', 'P.Value', 'adj.P.Val', 'Average.Log2.Expression', 'Norm.P.Value', 'Norm.adj.P.Val', 'Norm.logFC')],
        phosphodf[,c('id', 'pVal', 'BH', 'fold')],
        by.x='accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA',
        by.y='id'
    )

    plots$plot_protein_pvals <- function() {
        # Compare the calculated p-values for the modT vs the cyberT proteins
        # Keeps all isoform information
        qplot(
            -log2(adj.P.Val),
            signif_freq,
            data = plots$proteins_merged,
            alpha = count
        )
    }

    plots$plot_protein_foldchange <- function() {
        # Compare the calculated fold changes for the modT vs the cyberT proteins
        # Keeps all isoform information 
        qplot(
            log2(fold_change),
            Average.Log2.Expression,
            data = plots$proteins_merged,
            alpha = proteins_merged$count
        )
    }

    plots$plot_phospho_pvals <- function() {
        qplot(
            -log2(adj.P.Val),
            -log2(BH),
            data=plots$phospho_merged,
            alpha=0.1
        )
    }

    plots$plot_phospho_foldchange <- function() {
        qplot(
              Average.Log2.Expression,
              fold,
              data=plots$phospho_merged,
              alpha=0.1
        )
    }
    
    return(plots)
}

plots <- setup_plots(proteins_final, phosphos_cyberT)
