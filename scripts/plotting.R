library(ggplot2)
library(reshape2)

# Mapping for settings
# See format_results.py for original
LABEL_MAPPING = data.frame(label=c("Absolute fold change", "t-test (1-sample)", "t-test (2-sample)",
                                   "Moderated T (1-sample)", "Moderated T (2-sample)",
                                   "Mod T, variance trend (1-sample)",
                                   "Mod T, variance trend (2-sample)", "CyberT"), stringsAsFactors=FALSE)

rownames(LABEL_MAPPING) <- c("fold change", "t-test-1", "t-test-2",
                             "modT-1", "modT-2", "modT-1 trend", "modT-2 trend", "cyberT")
# Set up custom color map
valid_labels = unique(LABEL_MAPPING$label)
N_PEP = 8

colScale <- c("#C14242", "#FDA51F", "#0882FB", "#20E6FC", "#0BFD33", "#10B571", "#F826FB", "#4D03EC")
names(colScale) <- as.character(valid_labels[1:N_PEP])

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA"s
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data, measurevar, groupvars, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA"s: if na.rm==T, don"t count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group"s data frame, return a vector with
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
        . ~ setting + labels + variable,
        data = data[,c(measurevar, groupvars)],
        FUN = function(x) c(
            mn=mean(x),
            sd=sd(x),
            ci=(sd(x)/sqrt(length2(x))) * qt(conf.interval/2 + .5, length2(x)-1)
        )
    )))

    return(datac)
}

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
           legend.key.size = unit(1.5, "lines"),
           legend.position = "top",
    ) + 
    labs(title=title, x=xlab, y="AUROC", fill="")
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
           legend.key.size = unit(1.5, "lines"),
           legend.position = "top",
    ) + 
    labs(title=title, x=xlab, y="AUPRC", fill="")
  return(p1)
}

boxplot_results_pauc <- function(df, title="", xlab="Setting", fill=FALSE, skip_labels=c()){
  # Requires setting and labels column
  # Drop skipped tests
  df <- df[!(as.character(df$labels) %in% skip_labels),]
  df <- droplevels(df)
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
           legend.key.size = unit(1.5, "lines"),
           legend.position = "top",
    ) +
    labs(title=title, x=xlab, y="pAUROC", fill="")
  return(p1)
}

lineplot_results_pauc <- function(df, title="", xlab="Setting") {
  data <- summarySE(df, measurevar="pAUROC", groupvars=c("setting", "labels"))
  p1 <- ggplot(data, aes(x = setting, group=labels, colour=labels, y= 10^pAUROC)) +
        geom_errorbar(aes(ymin=pAUROC-se, ymax=pAUROC+se), width=.1) +
        geom_line() +
        geom_point()
  return(p1)
}

plot_fdr <- function(df, title="", xlab="Setting", skip_labels=c("t-test (1-sample)", "Moderated T (1-sample)", "Absolute fold change", "Mod T, variance trend (1-sample)")) {
  N_SIG <- 1000
  # Drop weaker tests
  df <- df[!(as.character(df$labels) %in% skip_labels),]
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
        id.vars=c("setting", "labels")) 
  data <- replace(data, is.na(data), 0)

  data.summary <- summarySE(data, "value", c("setting", "labels", "variable"))
  # Reorder levels for more useful line types
  data.summary$variable <- factor(data.summary$variable,
          levels(data.summary$variable)[c(4,2,3,1)])
  # Set alpha for raw p-values lower for emphasis
  setting_df <- as.data.frame(do.call(
        "rbind", strsplit(as.character(data.summary$variable), "_")))
  colnames(setting_df) <- c("fp_or_tp", "adj")
  data.summary$alpha <- 1
  data.summary$alpha[setting_df$adj == "raw"] <- 0.9

  p1 <- ggplot(data.summary, 
               aes(x = as.numeric(as.character(setting)),
                   y = value.mn,
                   color=labels,
                   fill=labels,
                   linetype=variable,
                   shape=variable,
                   size=variable)) + 
        geom_line() + geom_point(aes(stroke=2)) + 
        # TODO scale the size of the line and point for raw vs adj
        scale_size_manual(name="", values=c(1.2, 1.2, 0.6, 0.6)) +
        scale_linetype_manual(name="", values=c(1, 2, 1, 2)) +
        scale_shape_manual(name="", values=c(22, 22, 21, 21)) +
        # scale_stroke_manual(values=c(2, 2, 1, 1)) +
        guides(size=guide_legend(keywidth=4), fill=FALSE) + 
        labs(title=title, x=xlab, y="TPR or FDR", color="Test", linetype="", shape="")

  return(p1)
}
