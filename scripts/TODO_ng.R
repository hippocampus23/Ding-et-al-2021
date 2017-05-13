
normalize_protein_changes_kd <- function(phospho, total) {
    # Uses columns (Log2.Med.Rep1,2,3) from phospho
    # columns (Tot.Log2.Med.Rep1,3,3) from total
    # Merges on Accession.Number.NoIso column(!)

    #TODO make less brittle: don't specify exact column names!
    reduced_total <- aggregate(total[,c('Tot.Log2.Med.Rep1',  'Tot.Log2.Med.Rep2', 'Tot.Log2.Med.Rep3')],
                  by=list(Accession.Number.NoIso=total$Accession.Number.NoIso), FUN=mean)
    joint <- merge(phospho, reduced_total, by="Accession.Number.NoIso")

    joint$Log2.Norm.Rep1 <- joint$Log2.Med.Rep1 - joint$Tot.Log2.Med.Rep1
    joint$Log2.Norm.Rep2 <- joint$Log2.Med.Rep2 - joint$Tot.Log2.Med.Rep2
    joint$Log2.Norm.Rep3 <- joint$Log2.Med.Rep3 - joint$Tot.Log2.Med.Rep3
    joint$Log2.Norm.Avg <- (joint$Log2.Norm.Rep1 + joint$Log2.Norm.Rep2 + joint$Log2.Norm.Rep3) / 3

    return(joint)
}

normalize_protein_changes_oe <- function(phospho, total) {
    # Uses columns (Log2.Med.Rep1,2) from phospho
    # columns (Tot.Log2.Med.Rep1,3) from total
    # Merges on Accession.Number.NoIso column(!)

    # Average duplicate protein entries (different isoforms of same protein)
    reduced_total <- aggregate(total[,c('Tot.Log2.Med.Rep1',  'Tot.Log2.Med.Rep2')],
                  by=list(Accession.Number.NoIso=total$Accession.Number.NoIso), FUN=mean)
    joint <- merge(phospho, reduced_total, by="Accession.Number.NoIso")

    joint$Log2.Norm.Rep1 <- joint$Log2.Med.Rep1 - joint$Tot.Log2.Med.Rep1
    joint$Log2.Norm.Rep2 <- joint$Log2.Med.Rep2 - joint$Tot.Log2.Med.Rep2
    joint$Log2.Norm.Avg <- (joint$Log2.Norm.Rep1 + joint$Log2.Norm.Rep2) / 2

    return(joint)
}

#Utils

extract_accession <- function(info_col) {
    return(unlist(lapply(info_col, function(x) strsplit(x,"_")[[1]][1])))
}

strip_accession_isoform <- function(accession_numbers) {
    return(unlist(lapply(as.character(accession_numbers), function(x) strsplit(x,"-")[[1]][1])))
}

strip_whitespace_all_cols <- function(df) {
    i <- sapply(df, is.character)
    if (sum(i) == 0) {
        print("No columns are strings! Did you remember to convert from factors?")
    }
    df[i] <- lapply(df[i], function(x) unlist(lapply(x, function(y) str_trim(y))))
}

###Actual stuff starts happening here###

ngkd_phospho <- read.table("raw_data/NgKD_phospho.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
ngkd_total <- read.table("raw_data/NgKD_proteome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
ngoe_phospho <- read.table("raw_data/NgOE_phospho.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")
ngoe_total <- read.table("raw_data/NgOE_proteome.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE,quote = "")

# Fix col2 median normalization
ngkd_phospho$Log2.Med.Rep2 <- ngkd_phospho$Log2.Med.Rep2 - median(ngkd_phospho$Log2.Med.Rep2)

# Setup useful accession number cols
ngoe_phospho$Accession.Number <- extract_accession(ngoe_phospho[,1])
ngoe_phospho$Accession.Number.NoIso <- strip_accession_isoform(ngoe_phospho$Accession.Number)
ngoe_total$Accession.Number.NoIso <- strip_accession_isoform(ngoe_total$Accession.Number)
ngkd_phospho$Accession.Number.NoIso <- strip_accession_isoform(ngkd_phospho$Accession.Number)
ngkd_total$Accession.Number.NoIso <- strip_accession_isoform(ngkd_total$Accession.Number)

#Normalize joint stuff
ngkd_joint <- normalize_protein_changes_kd(ngkd_phospho, ngkd_total)
ngoe_joint <- normalize_protein_changes_oe(ngoe_phospho, ngoe_total)
# Find phosphosites not in total proteome
ngkd_no_total <- ngkd_phospho[!(ngkd_phospho$Accession.Number.NoIso %in% ngkd_total$Accession.Number.NoIso),]
ngoe_no_total <- ngoe_phospho[!(ngoe_phospho$Accession.Number.NoIso %in% ngoe_total$Accession.Number.NoIso),]
# Duplicate Log Change column for OE just for convenience
ngoe_no_total$Average.Log2.Expression <- ngoe_no_total$logFC

#Perform moderated t-test
source('modT.r')

modt_kd <- function(output_name) {
    modT.test(ngkd_joint, paste("modT_output/", output_name="kd", sep=""), dframe=TRUE,
          data.col=c("Log2.Norm.Rep1", "Log2.Norm.Rep2", "Log2.Norm.Rep3"))
}

modt_oe <- function(output_name) {
    modT.test(ngoe_joint, paste("modT_output/", output_name="oe", sep=""), dframe=TRUE,
          data.col=c("Log2.Norm.Rep1", "Log2.Norm.Rep2"))
}

print("Call modt_kd() and modt_oe() to run modT tests!")


# Extra one-liners for finding overlaps
# sum(ngkd_no_total$adj.P.Val < 0.05 & ngkd_no_total$Average.Log2.Expression < 0)
# sum(ngkd_no_total$adj.P.Val < 0.05 & ngkd_no_total$Average.Log2.Expression > 0)
# sum(ngkd_no_total$adj.P.Val < 0.05 & ngkd_no_total$Average.Log2.Expression > 0 & (ngkd_no_total$Accession.Number.NoIso %in% psd$Accession))
# sum(ngkd_no_total$adj.P.Val < 0.05 & ngkd_no_total$Average.Log2.Expression < 0 & (ngkd_no_total$Accession.Number.NoIso %in% psd$Accession))
# sum(!duplicated(ngkd_no_total$Accession.Number.NoIso[ngkd_no_total$Accession.Number.NoIso %in% psd$Accession]))
# sum(ngkd_no_total$Accession.Number.NoIso %in% psd$Accession)
# sum(!duplicated(ngoe_no_total$Accession.Number.NoIso))
# sum(ngoe_no_total$adj.P.Val < 0.05 & ngoe_no_total$Average.Log2.Expression < 0)
