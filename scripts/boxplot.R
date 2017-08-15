library(ggplot2)
library(RColorBrewer)

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


# Mapping for settings
# See format_results.py for original
LABEL_MAPPING = data.frame(label=c('CyberT', 'Moderated T (1-sample)', 'Moderated T (2-sample)', 'Absolute fold change', 'Absolute fold change', 't-test (2-sample)', 't-test (1-sample)', 'Median intensity absolute fold change', 'Median intensity Moderated T', 'Median intensity CyberT', 'Median intensity t-test', 'Weighted least squares', 'CyberT by peptide', 't-test by peptide'))
rownames(LABEL_MAPPING) <- c('cyberT', 'modT', 'modT (2-sample)', 'fold change', 'fold_change', 't-test', 't-test (1-sample)', 'fold_change_med', 'modT_PVal_med', 'cyberT_PVal_med', 'ttest_PVal_med', 'wls', 'cyberT_bypep', 'ttest_bypep')
# Set up custom color map
# TODO make this more readable by manually defining a color sequence
colors_ = brewer.pal(12, "Set3")
valid_labels = unique(LABEL_MAPPING$label)
N_PEP = 6
N_PROT_PARALLEL = 4
myColorsPeptide <- colors_[1:N_PEP]
names(myColorsPeptide) <- as.character(valid_labels[1:N_PEP])
myColorsProtein <- c(
        colors_[1:N_PROT_PARALLEL],
        tail(colors_, length(valid_labels) - N_PEP - N_PROT_PARALLEL)
)
names(myColorsProtein) <- as.character(valid_labels[(N_PEP+1):length(valid_labels)])
colScale <- scale_colour_manual(c(myColorsPeptide, myColorsProtein))


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
## order_x: default NULL, order of x settings for pretty display
## order_lab: default NULL, order of labels for pretty display 
## save: TRUE to save plot, FALSE otherwise
## colors: color scale for levels. NULL for default colors
read_data_and_plot <- function(df_name, title, xlab, plot="pAUC", order_x=NULL, order_lab=NULL, filename="", save=TRUE, colors=colScale) {
  if (is.character(df_name) & length(df_name) == 1) {
    # Read dataframe
    # df <- read.csv(paste("data_simulated/", df_name, sep=""))
    df <- read.csv(df_name)
    if (is.null(filename)) {
        filename <- strsplit(df_name,"\\.")[[1]][1]
    }
  } else {
    # Try to coerce to dataframe 
    df <- data.frame(df_name)
    if (is.null(filename) && save) {
      stop("Must provide filename if df_name is not character vector")
    }
  }
  df <- df[complete.cases(df),]
  # Map labels to pretty print
  df$labels[df$labels %in% rownames(LABEL_MAPPING)] <- 
    LABEL_MAPPING[df$labels[df$labels %in% rownames(LABEL_MAPPING)], 1]
  if (!is.null(order_x)) {
    fac <- as.factor(df$setting)
    df$setting <- factor(fac, levels=levels(fac)[order_x])
  }
  if (!is.null(order_lab)) {
    fac <- as.factor(df$labels)
    df$labels <- factor(fac, levels=levels(fac)[order_lab])
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

  if (!is.null(colors)) {
    p <- p + colors
  }
  
  # Name of file to save plot
  if (save) {
    ggsave(paste("simulated_plots/boxplots/AUTO", filename,
                 paste(plot,".png", sep=""), sep="_"),
        p,
        width=12.80,
        height=7.20,
        dpi=100)
  }
  return(p)
}

replot_all <- function(plot_names) {
  lapply(plot_names, FUN = function(x) read_data_and_plot(x[1], x[2], x[3]))
}


# Code for regenerating the final set of plots for the paper
plot_final <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=18))
  # FC range for uniform and gamma distributions 
  fc_range_gam <- read_data_and_plot(
      "FINAL_DATA/df_peptide_fc_range_gam_FINAL.csv",
      "", "log2(FC)", plot="pAUC", filename="", save=FALSE, colors=colScale)
  fc_range_uni <- read_data_and_plot(
      "FINAL_DATA/df_peptide_fc_range_uni_FINAL.csv",
      "", "log2(FC)", plot="pAUC", filename="", save=FALSE, colors=colScale)
  nexp <- read_data_and_plot(
      "FINAL_DATA/df_peptide_nexp_FINAL.csv",
      "", "Number of channels", plot="pAUC", save=FALSE, colors=colScale)
  nexp_imba <- read_data_and_plot(
      "FINAL_DATA/df_peptide_nexp_imba_FINAL.csv",
      "", "Number of channels", plot="pAUC", save=FALSE, colors=colScale)

  # Plot variance results
  # This needs subsetting + faceting to be maximally useful
  var_df <- read.csv('FINAL_DATA/df_peptide_variances_FINAL.csv')
  # Split out setting into background, noise, number
  setting_df <- as.data.frame(do.call(
        'rbind', strsplit(as.character(var_df$setting), '_')))
  colnames(setting_df) <- c('background', 'noise', 'setting')
  # Pretty print settings
  background_map <- setNames(c('Inverse ~ Gamma', 'Uniform'),c('invgam', 'uniform'))
  pretty_label_map <- setNames(c('beta', 'sigma^2'),c('invgam', 'uniform'))
  noise_map <- setNames(c('Gaussian', 'Laplacian', 'Scaled ~ t'), c('norm', 'lap', 't'))
  setting_df$pretty_label <- pretty_label_map[setting_df$background]
  setting_df$background <- background_map[setting_df$background]
  setting_df$noise <- noise_map[setting_df$noise]

  var_df$setting <- NULL
  plot_df <- cbind(var_df, setting_df)
  print(head(plot_df))
  # Now plot treating the factor as the number
  var <- read_data_and_plot(plot_df, "", "Variance", save=FALSE, colors=colScale)
  var <- var + facet_grid(noise ~ background + pretty_label, scales='free',
                          switch='x', labeller=label_parsed)

  ggsave("FINAL_PLOTS/peptide_fc_range_gam_FINAL.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/peptide_fc_range_uni_FINAL.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/peptide_nexp_FINAL.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/peptide_nexp_imba_FINAL.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/peptide_var_FINAL.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
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
