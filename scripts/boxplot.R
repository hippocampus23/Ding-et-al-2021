library(ggplot2)
# library(RColorBrewer)
library(reshape2)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data, measurevar, groupvars, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    if (FALSE) { 
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )
    }

    if (!"variable" %in% colnames(data)) {
      data$variable <- 1
    }
    datac <- as.data.frame(as.list(aggregate(
        #x = data[,measurevar],
        #by = as.list(data$setting),
        . ~ setting + labels + variable,
        data = data[,c(measurevar, groupvars)],
        FUN = function(x) c(
            mn=mean(x),
            sd=sd(x),
            ci=(sd(x)/sqrt(length2(x))) * qt(conf.interval/2 + .5, length2(x)-1)
        )
    )))

    # Rename the "mean" column    
    # datac <- rename(datac, c("mean" = measurevar))

    # datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    # ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    # datac$ci <- datac$se * ciMult

    return(datac)
}


# Mapping for settings
# See format_results.py for original
LABEL_MAPPING = data.frame(label=c('CyberT', 'ModT (1-sample)', 'ModT (2-sample)', 'Absolute Fold Change', 'Absolute Fold Change', 't-test (2-sample)', 't-test (1-sample)', 't-test (1-sample)', 'Median intensity absolute fold change', 'Median intensity ModT', 'Median intensity CyberT', 'Median intensity t-test', 'Weighted least squares', 'CyberT by peptide', 't-test by peptide'), stringsAsFactors=FALSE)
rownames(LABEL_MAPPING) <- c('cyberT', 'modT', 'modT (2-sample)', 'fold change', 'fold_change', 't-test', 't-test (1-sample)', 't test (1-sample)', 'fold_change_med', 'modT_PVal_med', 'cyberT_PVal_med', 'ttest_PVal_med', 'wls', 'cyberT_bypep', 'ttest_bypep')
# Set up custom color map
# TODO make this more readable by manually defining a color sequence
  valid_labels = unique(LABEL_MAPPING$label)
  N_PEP = 6
  N_PROT_PARALLEL = 4
if (FALSE) {
  colors_ = brewer.pal(9, "Set1")
  myColorsPeptide <- colors_[1:N_PEP]
  names(myColorsPeptide) <- as.character(valid_labels[1:N_PEP])
  myColorsProtein <- c(
          colors_[1:N_PROT_PARALLEL],
          tail(colors_, length(valid_labels) - N_PEP - N_PROT_PARALLEL)
  )
  names(myColorsProtein) <- as.character(valid_labels[(N_PEP+1):length(valid_labels)])
  colScale <- c(myColorsPeptide, myColorsProtein)
}
colScale <- hcl(h = seq(15, 375, length = N_PEP + 1), l = 65, c = 100)[1:N_PEP]
names(colScale) <- valid_labels[1:N_PEP]


boxplot_results_auroc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUROC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="grey") +
    scale_y_continuous(limits = quantile(df$AUROC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="grey" ),
           # Set legend size
           legend.key.size = unit(1.5, 'lines'),
           legend.position = 'top',
    ) + 
    labs(title=title, x=xlab, y="AUROC", fill='')
  return(p1)
}

boxplot_results_auprc <- function(df, title="", xlab="Setting"){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = AUPRC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="grey") +
    scale_y_continuous(limits = quantile(df$AUPRC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="grey" ),
           # Set legend size
           legend.key.size = unit(1.5, 'lines'),
           legend.position = 'top',
    ) + 
    labs(title=title, x=xlab, y="AUPRC", fill='')
  return(p1)
}

boxplot_results_pauc <- function(df, title="", xlab="Setting", fill=FALSE){
  # Requires setting and labels column
  p1 <- ggplot(df, aes(x = as.factor(setting), fill=labels, y = pAUROC)) + 
    # Boxplot
    geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
    # Add lines between each group
    geom_vline(xintercept=seq(1.5, length(unique(df$setting))-0.5, 1),
               lwd=0.5, color="grey") +
    scale_y_continuous(limits = quantile(df$pAUROC, c(0.01, 1.0))) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="grey" ),
           # Set legend size
           legend.key.size = unit(1.5, 'lines'),
           legend.position = 'top',
    ) + 
    labs(title=title, x=xlab, y="pAUROC", fill='')
  return(p1)
}

lineplot_results_pauc <- function(df, title="", xlab="Setting") {
  data <- summarySE(df, measurevar='pAUROC', groupvars=c('setting', 'labels'))
  p1 <- ggplot(data, aes(x = setting, group=labels, colour=labels, y= 10^pAUROC)) +
        geom_errorbar(aes(ymin=pAUROC-se, ymax=pAUROC+se), width=.1) +
        geom_line() +
        geom_point()
  return(p1)
}

plot_fdr <- function(df, title="", xlab="Setting") {
  N_SIG <- 1000
  # Drop t-test (1-sample) and modT (1-sample)
  df <- df[(as.character(df$labels) != 't-test (1-sample)')
           & (as.character(df$labels) !=  'ModT (1-sample)'),]
  df <- droplevels(df)
  fdr_raw <- df$FP_raw / (df$FP_raw + df$TP_raw)
  fdr_adj <- df$FP_adj / (df$FP_adj + df$TP_adj)
  # Melt dataframe
  data <- melt(data.frame(
            setting = df$setting,
            labels = df$labels,
            FDR_raw = fdr_raw,
            FDR_adj = fdr_adj,
            TPR_raw = df$TP_raw / N_SIG,
            TPR_adj = df$TP_adj / N_SIG), 
        id.vars=c('setting', 'labels')) 
  data <- replace(data, is.na(data), 0)

  data.summary <- summarySE(data, 'value', c('setting', 'labels', 'variable'))
  # Reorder levels for more useful line types
  data.summary$variable <- factor(data.summary$variable,
          levels(data.summary$variable)[c(4,2,3,1)])
  # Set alpha for raw p-values lower for emphasis
  setting_df <- as.data.frame(do.call(
        'rbind', strsplit(as.character(data.summary$variable), '_')))
  colnames(setting_df) <- c('fp_or_tp', 'adj')
  data.summary$alpha <- 1
  data.summary$alpha[setting_df$adj == 'raw'] <- 0.9

  p1 <- ggplot(data.summary, 
               aes(x = as.numeric(as.character(setting)),
                   y = value.mn,
                   color=labels,
                   linetype=variable,
                   shape=variable,
                   alpha=alpha)) + 
        geom_line(size=0.8) + geom_point(size=3) + 
        # TODO scale the size of the line and point for raw vs adj
        scale_alpha(guide = FALSE, range = c(0.5, 1)) +
        scale_linetype_manual(values=c(1, 2, 1, 2)) +
        scale_shape_manual(values=c(16, 16, 15, 15)) +
        guides(linetype=guide_legend(keywidth=4)) + 
        labs(title=title, x=xlab, y='TPR/FDR', color='Test', linetype='', shape='')

  return(p1)
}

## Main interface for creating boxplots
## df_name: string OR dataframe
##      If string, will attempt to read csv from file given by df_name
## title: Title of plot
## xlab: x label of plot
## plot: one of 'AUROC', 'AUPRC', 'pAUC', 'line', or 'fdr'
## order_x: default NULL, order of x settings for pretty display
## order_lab: default NULL, order of labels for pretty display 
## save: TRUE to save plot, FALSE otherwise
## colors: color scale for levels. NULL for default colors
read_data_and_plot <- function(df_name, title, xlab, plot="pAUC", order_x=NULL, order_lab=NULL, filename="", save=TRUE, colors=colScale, fill=colScale) {
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
  df$labels <- as.character(df$labels)
  df$labels[df$labels %in% rownames(LABEL_MAPPING)] <- 
    LABEL_MAPPING[df$labels[df$labels %in% rownames(LABEL_MAPPING)], 1]
  df$labels <- as.factor(df$labels)

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
  } else if (plot == 'fdr') {
    p <- plot_fdr(df, title, xlab)  
  } else {
    stop("Invalid plot specification. Must be one of 'AUROC', 'AUPRC', 'pAUC')")
  }

  if (!is.null(colors)) {
    p <- p + scale_colour_manual(values=colors)
  }
  if (!is.null(fill)) {
    p <- p + scale_fill_manual(values=colors)
  }
  
  # Name of file to save plot
  if (save) {
    ggsave(
        # paste("simulated_plots/boxplots/AUTO", filename,
        #          paste(plot,".png", sep=""), sep="_"),
        filename,
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


### GENERATE FINAL FIGURES ###

plots <- list()

plots$transform_var <- function(var_df, plot='pAUC') {
  # Split out setting into background, noise, number
  setting_df <- as.data.frame(do.call(
        'rbind', strsplit(as.character(var_df$setting), '_')))
  colnames(setting_df) <- c('background', 'noise', 'setting')
  # Pretty print settings, map abbreviations to full description
  background_map <- setNames(c('Inverse ~ Gamma', 'Uniform'),c('invgam', 'uniform'))
  pretty_label_map <- setNames(c('beta', 'sigma^2'),c('invgam', 'uniform'))
  noise_map <- setNames(c('Gaussian', 'Laplacian', 'Scaled ~ t'), c('norm', 'lap', 't'))
  setting_df$pretty_label <- pretty_label_map[setting_df$background]
  setting_df$background <- background_map[setting_df$background]
  setting_df$noise <- noise_map[setting_df$noise]
  # Add columns to original DF
  var_df$setting <- NULL
  var_plot_df <- cbind(var_df, setting_df)
  # Now plot, placing the numerical variance value on the x-axis
  var <- read_data_and_plot(var_plot_df, "", "Variance", save=FALSE, colors=colScale, plot=plot)
  var <- var + facet_grid(noise ~ background + pretty_label, scales='free',
                          switch='x', labeller=label_parsed)
  return(var)
}


plots$transform_ds_size <- function(ds_size_df, plot='pAUC') {
  # Split out setting into size and fraction
  setting_df <- as.data.frame(do.call(
        'rbind', strsplit(ds_size_df$setting, ': ')))
  colnames(setting_df) <- c('size', 'setting')
  setting_df$setting <- as.numeric(setting_df$setting)
  # setting_df$setting <- as.numeric(unlist(lapply(
  #       as.character(setting_df$setting), function(x) strsplit(x," ")[[1]][1])))
  ds_size_df$setting <- NULL
  ds_size_plot_df <- cbind(ds_size_df, setting_df)
  ds_size <- read_data_and_plot(ds_size_plot_df, "", "Number of peptides changed",
                                save=FALSE, colors=colScale, plot=plot)
  ds_size <- ds_size + facet_grid(. ~ size, scales='free_x') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(ds_size)
}

# Code for regenerating the final set of plots for the paper
plots$plot_final <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=26))
  # FC range for uniform and gamma distributions 
  fc_range_gam <- read_data_and_plot(
      "FINAL_DATA/df_peptide_fc_range_gam_FINAL.csv",
      "", "log2(FC)", plot="pAUC", filename="", save=FALSE, colors=colScale)
  fc_range_uni <- read_data_and_plot(
      "FINAL_DATA/df_peptide_fc_range_uni_FINAL.csv",
      "", "log2(FC)", plot="pAUC", filename="", save=FALSE, colors=colScale)
  nexp_fix <- read_data_and_plot(
      "FINAL_DATA/df_peptide_nexp_modtfix_FINAL.csv",
      "", "Number of channels", plot="pAUC", save=FALSE, colors=colScale)
  nexp_imba_fix <- read_data_and_plot(
      "FINAL_DATA/df_peptide_nexp_imba_modtfix_FINAL.csv",
      "", "(Number control, number experimental) channels", plot="pAUC", save=FALSE, colors=colScale)

  fdr_fc_gam_df <- read.csv('FINAL_DATA/df_peptide_fdr_fc_gam_FINAL.csv')
  fdr_fc_gam <- read_data_and_plot(
      fdr_fc_gam_df[as.numeric(fdr_fc_gam_df$setting) <= 1.0,],
      "", "log2(FC)", plot='fdr', filename="", save=FALSE, colors=colScale)
  fdr_fc_uni_df <- read.csv('FINAL_DATA/df_peptide_fdr_fc_uni_FINAL.csv')
  fdr_fc_uni <- read_data_and_plot(
      fdr_fc_uni_df[as.numeric(fdr_fc_uni_df$setting) <= 1.0,],
      "", "log2(FC)", plot='fdr', filename="", save=FALSE, colors=colScale)

  # Plot variance results
  # This needs subsetting + faceting to be maximally useful
  var <- plots$transform_var(
      read.csv('FINAL_DATA/df_peptide_variances_FINAL.csv'),
      plot='pAUC')

  # Save plots
  ggsave("FINAL_PLOTS/3B.png", fc_range_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/3D.png", fc_range_uni,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/5B.png", nexp_fix,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/5C.png", nexp_imba_fix,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/6B.png", var,
        width=12.80, height=12.80, dpi=100)

  ggsave("FINAL_PLOTS/4B.png", fdr_fc_gam,
        width=12.80, height=7.20, dpi=100)
  ggsave("FINAL_PLOTS/4D.png", fdr_fc_uni,
        width=12.80, height=7.20, dpi=100)
}

plots$plot_supplementary <- function() {
  # S3
  # Fold change range AUROC and AUPRC
  fc_range_uni_df <- read.csv('FINAL_DATA/df_peptide_fc_range_uni_FINAL.csv')
  fc_range_gam_df <- read.csv('FINAL_DATA/df_peptide_fc_range_gam_FINAL.csv')
  fc_range_uni <- read_data_and_plot(
      fc_range_uni_df, "S3A", "log2(FC)", plot="AUROC",
      filename="FINAL_PLOTS/supp/S3A.eps", save=TRUE, colors=colScale)
  fc_range_uni <- read_data_and_plot(
      fc_range_uni_df, "S3B", "log2(FC)", plot="AUPRC",
      filename="FINAL_PLOTS/supp/S3B.eps", save=TRUE, colors=colScale)
  fc_range_gam <- read_data_and_plot(
      fc_range_gam_df, "S3C", "log2(FC)", plot="AUROC",
      filename="FINAL_PLOTS/supp/S3C.eps", save=TRUE, colors=colScale)
  fc_range_gam <- read_data_and_plot(
      fc_range_gam_df, "S3D", "log2(FC)", plot="AUPRC",
      filename="FINAL_PLOTS/supp/S3D.eps", save=TRUE, colors=colScale)

  # S4
  # Fold change range AUROC and AUPRC
  fc_nexp <- read.csv('FINAL_DATA/df_peptide_nexp_modtfix_FINAL.csv')
  fc_nexp_imba <- read.csv('FINAL_DATA/df_peptide_nexp_imba_modtfix_FINAL.csv')
  tmp <- read_data_and_plot(
      fc_nexp, 'S4A', 'Number of Channels', plot='AUROC',
      filename='FINAL_PLOTS/supp/S4A.eps', save=TRUE, colors=colScale)
  tmp <- read_data_and_plot(
      fc_nexp, 'S4B', '(Number control, number experimental) channels', plot='AUPRC',
      filename='FINAL_PLOTS/supp/S4B.eps', save=TRUE, colors=colScale)
  tmp <- read_data_and_plot(
      fc_nexp_imba, 'S4C', 'Number of Channels', plot='AUROC',
      filename='FINAL_PLOTS/supp/S4C.eps', save=TRUE, colors=colScale)
  tmp <- read_data_and_plot(
      fc_nexp_imba, 'S4D', '(Number control, number experimental) channels', plot='AUPRC',
      filename='FINAL_PLOTS/supp/S4D.eps', save=TRUE, colors=colScale)

  var_df <- read.csv('FINAL_DATA/df_peptide_variances_FINAL.csv')
  ggsave("FINAL_PLOTS/supp/S5A.eps", plots$transform_var(var_df, plot='AUROC'),
        width=12.80, height=12.80, dpi=100)
  ggsave("FINAL_PLOTS/supp/S5B.eps", plots$transform_var(var_df, plot='AUPRC'),
        width=12.80, height=12.80, dpi=100)

  # Plot size of dataset results
  # This also needs facetting
  ds_size <- plots$transform_ds_size(
      read.csv('FINAL_DATA/df_peptide_dataset_size_FINAL.csv',
          stringsAsFactors=FALSE),
      plot='pAUC')
  ggsave("FINAL_PLOTS/supp/S2B.png", ds_size,
        width=12.80, height=7.20, dpi=100)
}

# Ordering function
# order(as.double(lapply(strsplit(as.character(x), '<'), function(l) l[1])))
