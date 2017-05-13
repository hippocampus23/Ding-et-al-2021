do_stuff <- function(df, print=TRUE) {
    if (print) {
        print("TESTING")
    }
    summary <- colSums(df)
    return(summary)
}
