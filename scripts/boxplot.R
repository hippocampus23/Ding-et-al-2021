library(ggplot2)

boxplot_results_auroc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUROC)) + 
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    scale_y_continuous(limits = quantile(df$AUROC, c(0.01, 1.0))) +
    labs(title=paste(title, "AUROC"), x=xlab, y="AUROC") +
    coord_flip()
  return(p1)
}

boxplot_results_auprc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUPRC)) + 
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    scale_y_continuous(limits = quantile(df$AUPRC, c(0.01, 1.0))) +
    labs(title=paste(title, "AUPRC"), x=xlab, y="AUPRC") +
    coord_flip()
  return(p1)
}

boxplot_results_pauc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = pAUROC)) + 
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    scale_y_continuous(limits = quantile(df$pAUROC, c(0.01, 1.0))) +
    labs(title=paste(title, "pAUROC"), x=xlab, y="pAUROC") +
    coord_flip()
  return(p1)
}
