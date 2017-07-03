library(ggplot2)

boxplot_results_auroc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUROC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="white") +
    scale_y_continuous(limits = quantile(df$AUROC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="white" ) 
    ) + 
    labs(title=paste(title, "AUROC"), x=xlab, y="AUROC")
  return(p1)
}

boxplot_results_auprc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUPRC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="white") +
    scale_y_continuous(limits = quantile(df$AUPRC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="white" ) 
    ) + 
    labs(title=paste(title, "AUPRC"), x=xlab, y="AUPRC")
  return(p1)
}

boxplot_results_pauc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = pAUROC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="white") +
    scale_y_continuous(limits = quantile(df$pAUROC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="white" ) 
    ) + 
    labs(title=paste(title, "pAUROC"), x=xlab, y="pAUROC")
  return(p1)
}

read_data_and_plot <- function(df_name, title, xlab) {
  png(paste("AUTO_", strsplit(df_name,"\\.")[[1]][1], "_pauc.png", sep=""),
      width=1280,
      height=720,
      units="px")
  df = read.csv(df_name)
  boxplot_results_pauc(df, title, xlab)
  dev.off() 
}

plots <- list(
    c("df_inv_gam_lap.csv", "Fold change varies, inverse gamma background, Laplacian noise", "log2(fold change)"),
    c("df_inv_gam_norm.csv", "Fold change varies, inverse gamma background, normal noise", "log2(fold change)"),
    c("df_uni_lap.csv", "Fold change varies, uniform background, Laplacian noise", "log2(fold change)"),
    c("df_nexp.csv", "Number of channels varies", "Number of channels"),
    c("df_random_fc.csv", "Fold change is randomly distributed", "Variance model"),
    c("df_var.csv", "Fixed fold change", "Variance model")
)
