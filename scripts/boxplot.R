library(ggplot2)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


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


lineplot_results_pauc <- function(df, title="", xlab="Setting") {
  data <- summarySE(df, measurevar='pAUROC', groupvars=c('setting', 'labels'))
  print(head(data))
  p1 <- ggplot(data, aes(x = setting, group=labels, colour=labels, y= 10^pAUROC)) +
        geom_errorbar(aes(ymin=pAUROC-se, ymax=pAUROC+se), width=.1) +
        geom_line() +
        geom_point()
  return(p1)
}

## Main interface for creating boxplots
## df_name: string OR dataframe
##      If string, will attempt to read csv from file given by df_name
## title: Title of plot
## xlab: x label of plot
## plot: one of 'AUROC', 'AUPRC', or 'pAUC'
## order_x: default NULL, function which orders [x-]settings for pretty display
## order_lab: default NULL, function which orders labels WITHIN each x-setting
read_data_and_plot <- function(df_name, title, xlab, plot="pAUC", order_x=NULL, order_lab=NULL, filename="") {
  if (is.character(df_name) & length(df_name) == 1) {
    # Read dataframe
    df <- read.csv(paste("data_simulated/", df_name, sep=""))
    if (is.null(filename)) {
        filename <- strsplit(df_name,"\\.")[[1]][1]
    }
  } else {
    # Try to coerce to dataframe 
    df <- data.frame(df_name)
    if (is.null(filename)) {
      stop("Must provide filename if df_name is not character vector")
    }
  }
  df <- df[complete.cases(df),]
  if (!is.null(order_x)) {
    fac <- as.factor(df$setting)
    df$setting <- factor(fac, levels=levels(fac)[order_x(levels(fac))])
  }
  if (!is.null(order_lab)) {
    fac <- as.factor(df$labels)
    df$labels <- factor(fac, levels=levels(fac)[order_lab(levels(fac))])
  }
  if (plot == "AUROC") {
    p <- boxplot_results_auroc(df, title, xlab)
  } else if (plot == "AUPRC") {
    p <- boxplot_results_auprc(df, title, xlab)
  } else if (plot == "pAUC") {
    p <- boxplot_results_pauc(df, title, xlab)
  } else if (plot == 'line') {
    p <- lineplot_results_pauc(df, title, xlab)
  } else {
    stop("Invalid plot specification. Must be one of 'AUROC', 'AUPRC', 'pAUC')")
  }
  
  # Name of file to save plot
  ggsave(paste("simulated_plots/boxplots/AUTO", filename,
               paste(plot,".png", sep=""), sep="_"),
      p,
      width=12.80,
      height=7.20,
      dpi=100)
  return(p)
}

replot_all <- function(plot_names) {
  lapply(plot_names, FUN = function(x) read_data_and_plot(x[1], x[2], x[3]))
}

plot_names <- list(
    c("df_inv_gam_lap.csv", "Fold change varies, inverse gamma background, Laplacian noise", "log2(fold change)"),
    c("df_inv_gam_norm.csv", "Fold change varies, inverse gamma background, normal noise", "log2(fold change)"),
    c("df_uni_lap.csv", "Fold change varies, uniform background, Laplacian noise", "log2(fold change)"),
    c("df_nexp.csv", "Number of channels varies", "Number of channels"),
    c("df_random_fc.csv", "Fold change is randomly distributed", "Variance model"),
    c("df_var.csv", "Fixed fold change", "Variance model"),
    c("df_compare_one_two_sided.csv", "One vs two sided fold change", "Fold change"),
    c("df_with_without_central_fc.csv", "Normally distributed FC, with and without central variation")
)

plot_protein_names <- list(
    c("df_prot_fc_range_uni_2.csv", "Fold change varies, 2 peptides, uniform background", "log2(fold_change)"),
    c("df_prot_fc_range_uni_4.csv", "Fold change varies, 4 peptides, uniform background", "log2(fold_change)"),
    c("df_prot_fc_range_gam_2.csv", "Fold change varies, 2 peptides, Inv Gam background", "log2(fold_change)"),
    c("df_prot_fc_range_gam_4.csv", "Fold change varies, 4 peptides, Inv Gam background", "log2(fold_change)"),
    c("df_prot_num_peps.csv", "Number of peptides varies, log2(FC)=0.3", "Variance and peptide count")
)

# Ordering function
# order(as.double(lapply(strsplit(as.character(x), '<'), function(l) l[1])))
