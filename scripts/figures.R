source("plotting.R")

## Main interface for creating boxplots
## df_name: string OR dataframe
##      If string, will attempt to read csv from file given by df_name
## title: Title of plot
## xlab: x label of plot
## plot: one of "AUROC", "AUPRC", "pAUC", "line", or "fdr"
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
  } else if (plot == "line") {
    p <- lineplot_results_pauc(df, title, xlab)
  } else if (plot == "fdr") {
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
        filename,
        p,
        width=12.80,
        height=7.20,
        dpi=100,
        bg = "transparent"
    )
  }
  return(p)
}

transform_var <- function(var_df, plot="pAUC") {
  # Split out setting into background, noise, number
  setting_df <- as.data.frame(do.call(
        "rbind", strsplit(as.character(var_df$setting), "_")))
  colnames(setting_df) <- c("background", "noise", "setting")
  # Pretty print settings, map abbreviations to full description
  background_map <- setNames(c("Inverse ~ Gamma", "Uniform"),c("invgam", "uniform"))
  pretty_label_map <- setNames(c("beta", "sigma^2"),c("invgam", "uniform"))
  noise_map <- setNames(c("Gaussian", "Laplacian", "Scaled ~ t"), c("norm", "lap", "t"))
  setting_df$pretty_label <- pretty_label_map[setting_df$background]
  setting_df$background <- background_map[setting_df$background]
  setting_df$noise <- noise_map[setting_df$noise]
  # Add columns to original DF
  var_df$setting <- NULL
  var_plot_df <- cbind(var_df, setting_df)
  # Now plot, placing the numerical variance value on the x-axis
  var <- read_data_and_plot(var_plot_df, "", "Variance", save=FALSE, plot=plot)
  var <- var + facet_grid(noise ~ background + pretty_label, scales="free",
                          switch="x", labeller=label_parsed)
  return(var)
}


transform_ds_size <- function(ds_size_df, plot="pAUC") {
  # Split out setting into size and fraction
  setting_df <- as.data.frame(do.call(
        "rbind", strsplit(as.character(ds_size_df$setting), "_")))
  colnames(setting_df) <- c("size", "setting")
  setting_df$setting <- strtoi(setting_df$setting)
  ds_size_df$setting <- setting_df$setting
  ds_size_plot_df <- cbind(ds_size_df, setting_df)
  ds_size <- read_data_and_plot(ds_size_plot_df, "", "Number of peptides changed", save=FALSE, plot=plot)
  ds_size <- ds_size + facet_grid(. ~ size, scales="free_x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(ds_size)
}


# Code for regenerating the plots for the paper
fig_3B <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot("../data_simulated/3B.csv", "", "log2(FC)", plot="pAUC",
		     save=TRUE, filename="../figures/3B.eps")
  print("saved 3B.eps")
}


fig_3D <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot("../data_simulated/3D.csv", "", "log2(FC)", plot="pAUC",
		     save=TRUE, filename="../figures/3D.eps")
  print("saved 3D.eps")
}


fig_3F <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot("../data_simulated/3F.csv",
      "", "log2(FC)", plot="pAUC", save=TRUE, filename="../figures/3F.eps")
  print("saved 3F.eps")
}


fig_4D <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  fdr_fc_uni_df <- read.csv("../data_simulated/4D.csv")
  read_data_and_plot(
      fdr_fc_uni_df[as.numeric(fdr_fc_uni_df$setting) <= 1.0,],
      "", "log2(FC)", plot="fdr", save=TRUE, filename="../figures/4D.eps")
  print("saved 4D.eps")
}


fig_4H <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  fdr_fc_gam_df <- read.csv("../data_simulated/4H.csv")
  read_data_and_plot(
      fdr_fc_gam_df[as.numeric(fdr_fc_gam_df$setting) <= 1.0,],
      "", "log2(FC)", plot="fdr", save=TRUE, filename="../figures/4H.eps")
  print("saved 4H.eps")
}


fig_4L <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  fdr_fc_trend_df <- read.csv("../data_simulated/4L.csv")
  read_data_and_plot(
      fdr_fc_trend_df[as.numeric(fdr_fc_trend_df$setting) <= 1.0,],
      "", "log2(FC)", plot="fdr", save=TRUE, filename="../figures/4L.eps")
  print("saved 4L.eps")
}


fig_5B <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot("../data_simulated/5B.csv",
      "", "Number of channels", plot="pAUC", save=TRUE, filename="../figures/5B.eps")
  print("saved 5B.eps")
}


fig_5D <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot("../data_simulated/5D.csv",
      "", "Number of channels", plot="pAUC", save=TRUE, filename="../figures/5D.eps")
  print("saved 5D.eps")
}


fig_6C <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  # Plot variance results
  # This needs subsetting + faceting to be maximally useful
  var <- transform_var(read.csv("../data_simulated/6C.csv"), plot="pAUC")
  ggsave("../figures/6C.eps", var, width=12.80, height=12.80, dpi=100, bg = "transparent")
  print("saved 6C.eps")
}


# supplementary plots

fig_S2A <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3B.csv", "Uniform peptide variance", "log2(FC)", plot="AUROC",
      filename="../figures/S2A.eps", save=TRUE)
  print("saved S2A.eps")
}


fig_S2E <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3F.csv", "Intensity dependent inverse gamma variance", "log2(FC)", plot="AUROC",
      filename="../figures/S2E.eps", save=TRUE)
  print("saved S2E.eps")
}


fig_S2B <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3B.csv", "Uniform peptide variance", "log2(FC)", plot="AUPRC",
      filename="../figures/S2B.eps", save=TRUE)
  print("saved S2B.eps")
}


fig_S2F <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3F.csv", "Intensity dependent inverse gamma variance", "log2(FC)", plot="AUPRC",
      filename="../figures/S2F.eps", save=TRUE)
  print("saved S2F.eps")
}


fig_S2C <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3D.csv", "Inverse gamma variance", "log2(FC)", plot="AUROC",
      filename="../figures/S2C.eps", save=TRUE)
  print("saved S2C.eps")
}


fig_S2D <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/3D.csv", "Inverse gamma variance", "log2(FC)", plot="AUPRC",
      filename="../figures/S2D.eps", save=TRUE)
  print("saved S2D.eps")
}


fig_S3A <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/5B.csv", "Intensity dependent inverse gamma variance", "Number of channels", plot="AUROC",
      filename="../figures/S3A.eps", save=TRUE)
  print("saved S3A.eps")
}


fig_S3B <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/5B.csv", "Intensity dependent inverse gamma variance", "Number of channels", plot="AUPRC",
      filename="../figures/S3B.eps", save=TRUE)
  print("saved S3B.eps")
}


fig_S3C <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/5D.csv", "Intensity dependent inverse gamma variance", "Number of channels", plot="AUROC",
      filename="../figures/S3C.eps", save=TRUE)
  print("saved S3C.eps")
}


fig_S3D <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  read_data_and_plot(
      "../data_simulated/5D.csv", "Intensity dependent inverse gamma variance", "Number of channels", plot="AUPRC",
      filename="../figures/S3D.eps", save=TRUE)
  print("saved S3D.eps")
}


fig_S4A <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  # Plot variance results
  # This needs subsetting + faceting to be maximally useful
  var <- transform_var(read.csv("../data_simulated/6C.csv"), plot="AUROC")
  ggsave("../figures/S4A.eps", var, width=12.80, height=12.80, dpi=100, bg = "transparent")
  print("saved S4A.eps")
}


fig_S4B <- function() {
  # Set font size and theme
  theme_set(theme_bw(base_size=23))
  # Plot variance results
  # This needs subsetting + faceting to be maximally useful
  var <- transform_var(read.csv("../data_simulated/6C.csv"), plot="AUPRC")
  ggsave("../figures/S4B.eps", var, width=12.80, height=12.80, dpi=100, bg = "transparent")
  print("saved S4B.eps")
}


fig_S1B <- function() {
  # Plot size of dataset results
  # This also needs facetting
  theme_set(theme_bw(base_size=23))
  ds_size<- transform_ds_size(
      read.csv("../data_simulated/S1B.csv", stringsAsFactors=FALSE), plot='pAUC')
  ggsave("../figures/S1B.eps", ds_size, width=12.80, height=7.20, dpi=100, bg = "transparent")
  print("saved S1B.eps")
}
