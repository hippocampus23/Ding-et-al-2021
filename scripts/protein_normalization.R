# TODO correlation the intensity scores for different peptides of a single protein
# BUT correlation across which measures?
# TODO variance standardizing normalization(?)

# http://www.nature.com/nbt/journal/v28/n1/full/nbt.1592.html


#################
# PREPROCESSING #
#################

map_peptides_to_proteins <- function() {
    # Normalize on a per-channel basis to the sum of intensities in that column
    # NOTE to make the numbers more interpretable I have normalized by the means
    # This should be fine as long as there are no NAs
    # NOTE Can also simply use the VSN applied raw intensities from pep instead!
    peptides <- sweep(c.input[,2:7], 2, colMeans(c.input[,2:7]), '/')
    peptides$id <- c.input$id

    print('Merging peptides to acc_num')
    # Merge peptides to protein names and acc numbers
    # Then merge accession numbers to peptide table to filter to mouse proteins
    keys = merge(x=c.Key, y=unique(pep[which(pep$species == 'MOUSE'), c('accession_number', 'species')]), by='accession_number')
    peptides <- merge(x=peptides, y=keys, by='id')
    return(peptides)
}


find_protein_medians <- function(pepdf) {
    # Finds the median peptide of each protein (as judged by acc num) for each channel
    # AND the overall median for the control and experimental conditions

    print('Obtaining medians')
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

    print(colnames(peptide_merged))
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

############################
# ACTUAL STUFF STARTS HERE #
############################

# peptides <- map_peptides_to_proteins
# proteins <- find_protein_medians(peptides)

