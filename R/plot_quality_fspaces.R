# Function to illustrate the quality of multidimensional functional spaces
#
# Authors: Sébastien Villéger & Camille Magneville
#
#------------------------------------------------------------------------------
#
#
#' Plot functional space quality with a chosen quality metric
#'
#' @param fspaces_quality output from the \code{\link{quality.fspaces}} function,
#'   that is a list with all data needed to illustrate quality of functional
#'   spaces based on deviation between species trait-based distance and distance
#'   in functional spaces built using PCoA (and dendrogram).
#'
#' @param fspaces_plot a \strong{vector} with names of functional spaces to consider.
#'   Should be a subset of the row names of
#'   \code{fspaces_quality$quality_fspaces}. Maximum of 10 spaces allowed to keep decent
#'   plot size.
#'
#' @param quality_metric a \strong{character string} with the name of the quality metric
#'   to illustrate. Should be one of the column names of
#'   \code{fspaces_quality$quality_fspaces}. See help of \code{\link{quality.fspaces}}
#'   for the meaning of these names regarding type of deviation and scaling of
#'   distance in functional space. Default: 'mad' (Mean absolute deviation).
#'
#' @param name_file a \strong{character string} with name of file to save the figure
#'   (without extension). Default: NULL which means plot is displayed.
#'
#' @param range_dist a \strong{vector} with minimum and maximum values to display for
#'   species pairwise distances (x-axis for all panels and y-axes of top panel).
#'   Default: NULL, which means range is 0 to maximum distance among all the
#'   functional spaces to plot.
#'
#' @param range_dev a \strong{vector} with minimum and maximum values to display for
#'   deviation to trait-based distance (y-axis of middle panel). Default: NULL,
#'   which means range is set to range of deviation among all the functional
#'   spaces to plot.
#'
#' @param range_qdev a \strong{vector} with minimum and maximum values to display for
#'   deviation to trait-based distance (y-axis of bottom panel). Default:NULL,
#'   which means range is from 0 to the maximum of (transformed) deviation among
#'   all the functional spaces to plot.
#'
#' @param gradient_deviation a \strong{vector} of 3 colors for illustrating raw deviation
#'   with \code{\link[ggplot2]{scale_colour_gradient2}}. The first value ('neg') is for the
#'   lowest negative deviation, the second value ('nul')is for null deviation and
#'   the third value ('pos') is for the highest positive deviation. Default gradient
#'   is from darkblue to grey to red.
#'
#' @param gradient_deviation_quality 2 colors (named 'low' and 'high') for
#'   illustrating transformed deviation used to compute quality metric with
#'   \code{\link[ggplot2]{scale_colour_gradient2}} (default gradient is from yellow to red).
#'
#' @param x_lab a \strong{character string} with title to display below X axis. Default
#'   is 'Trait-based distance'.
#'
#' @return a png file (resolution 300dpi) saved in the current working
#'   directory. Quality of each functional space is illustrated with three
#'   panels : - top row shows trait-based distance between species vs.
#'   space-based distance. - middle row shows trait-based distance vs. deviation
#'   between space-based and trait-based distances - bottom row shows
#'   trait-based distance between species vs. transformed deviation used to
#'   compute the quality metric All plots have the same X axis. All plots on a given
#'   row have the same Y axis and color palette. Type of distance in functional
#'   space (Euclidean in PCoA, Cophenetic on tree) are abbreviated, as well as
#'   type of transformation of distance (scaling) and of deviation (Absolute or
#'   Squared)
#'
#' @examples
#' # Load Species*Traits dataframe:
#' data("sp_tr_fruits", package = "mFD")
#' # Load Assemblages*Species dataframe:      
#' data("asb_sp_w_fruits", package = "mFD")
#' # Load Traits categories dataframe:
#' data("tr_cat_fruits", package = "mFD")   
#' # Compute functional distance 
#' sp_dist_fruits <- mFD::funct.dist(sp_tr = sp_tr_fruits,         
#'  tr_cat       = tr_cat_fruits,   
#'  dist_metric  = "kgower",         
#'  scaling      = "scaledBYrange",  
#'  stop_if_NA   = TRUE)
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#' fspaces_quality_fruits <- mFD::quality.fspaces(sp_dist = sp_dist_fruits, 
#'  maxdim_pcoa         = 10,
#'  deviation_weighting = "absolute",
#'  fdist_scaling       = FALSE,
#'  fdendro             = "average")
#'  # Illustarte the quality of functional spaces:
#'  mFD::quality.fspaces.plot(fspaces_quality = fspaces_quality_fruits, 
#'   quality_metric = "mad",
#'   fspaces_plot = c("tree_average", "pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d"),
#'   name_file = NULL, range_dist = NULL, range_dev = NULL, range_qdev = NULL,
#'   gradient_deviation  = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
#'   gradient_deviation_quality = c(low ="yellow", high = "red"),
#'   x_lab = "Trait-based distance")
#'
#' @export


quality.fspaces.plot <- function(
  fspaces_quality, quality_metric, fspaces_plot,
  name_file = NULL,
  range_dist = NULL, range_dev = NULL, range_qdev = NULL,
  gradient_deviation = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low ="yellow", high = "red"),
  x_lab = "Trait-based distance") {
  
  
  #### check inputs  ####
  
  
  # check core input with data from quality_fspaces:
  if(any(! names(fspaces_quality) %in% c("quality_fspaces", "details_trdist",
                                         "details_fspaces", "details_deviation"))) {
    stop("Error: input 'fspaces_quality' should be the output from function
       'mFD::quality_fspaces'.")
  }
  
  # check names and number of functional spaces:
  if (any(! fspaces_plot %in% row.names(fspaces_quality$quality_fspaces))) {
    stop("Error: input 'fspaces_plot' should be a subset of the row names of
      'fspaces_quality$quality_fspaces'.")
  }
  if (length(fspaces_plot) > 10){
    stop("Error: input 'fspaces_plot' should contain no more than 5 names of
     functional spaces.")
  }
  
  # check type of quality metric:
  if (any(! quality_metric %in% names(fspaces_quality$quality_fspaces))) {
    stop("Error: input 'quality_metrics' should be one of the column names of
     'fspaces_quality$quality_fspaces'.")
  }
  
  
  #### setting parameters for all plots ####
  
  # detailed name of functional spaces ----
  fspaces_nm_plot <- gsub(fspaces_plot, pattern = "_", replacement = " ")
  fspaces_nm_plot <- gsub(fspaces_nm_plot, pattern = "pcoa", replacement = "PCoA")
  fspaces_nm_plot <- gsub(fspaces_nm_plot, pattern = "tree", replacement = "Tree")
  fspaces_nm_plot <- gsub(fspaces_nm_plot, pattern = "d", replacement = "D")
  fspaces_nm_plot <- substr(fspaces_nm_plot, 1, 12)
  names(fspaces_nm_plot) <- fspaces_plot
  
  # detailed names of deviation used for quality metrics ----
  nm_dev_qual <- c("Abs. Dev. of", "Squ. Dev. of",
                   "Abs. Dev. of scaled", "Squ. Dev. of scaled")
  names(nm_dev_qual) <- c("mad","rmsd", "mad_scaled", "rmsd_scaled")
  
  # detailed names of quality metrics ----
  nm_qual_metrics <- c("Mean Abs. Dev.", "Root Mean Squ. Dev.",
                       "Mean Abs. Dev. scld", "Root Mean Squ. Dev. scld")
  names(nm_qual_metrics) <- c("mad","rmsd", "mad_scaled", "rmsd_scaled")
  
  # parameters for ggplot2 functions ----
  # scaling parameter for font size of axes and plot titles:
  scaling_text <- 7.5 + (length(fspaces_nm_plot) - 1) / 5
  
  # point size
  point_size <- 0.1
  
  
  #### extracting data for all plots ####
  
  # pairwise distances based on traits and in functional spaces (raw) ----
  df_dist <- fspaces_quality$details_fspaces$
    pairsp_fspaces_dist[, c("tr", fspaces_plot)]
  
  # raw deviation ----
  list_dev <- fspaces_quality$details_deviation
  df_dev_dist <- data.frame(list_dev$dev_distsp[, fspaces_plot])
  names(df_dev_dist) <- fspaces_plot
  
  # transformed deviation for quality metric ----
  if (quality_metric == "mad") {
    df_qdev_dist <- list_dev$abs_dev_distsp[, fspaces_plot]
  }
  if (quality_metric == "rmsd") {
    df_qdev_dist<-list_dev$sqr_dev_distsp[, fspaces_plot]
  }
  
  if (quality_metric == "mad_scaled") {
    df_qdev_dist <- list_dev$abs_dev_distsp_scaled[, fspaces_plot]
  }
  
  if (quality_metric == "rmsd_scaled") {
    df_qdev_dist <- list_dev$sqr_dev_distsp_scaled[, fspaces_plot]
  }
  
  df_qdev_dist <- data.frame(df_qdev_dist)
  names(df_qdev_dist) <- fspaces_plot
  
  
  # computing ranges for axes if not provided as input ----
  
  # ranges of distances (trait-based & in functional spaces)
  range_dist <- range_dist
  if (is.null(range_dist)) {
    range_dist <- c(0, max(df_dist))
  }
  
  # range of deviation among all spaces to plot
  range_dev <- range_dev
  if (is.null(range_dev)) {
    range_dev <- range(df_dev_dist)
  }
  
  # range of deviation among all spaces to plot
  range_qdev <- range_qdev
  if (is.null(range_qdev)) {
    range_qdev <- c(0, max(df_qdev_dist))
  }
  
  
  #### plotting quality of functional spaces ####
  
  # loop on functional spaces
  for (pos_k in 1:length(fspaces_plot)) {
    
    # name of space # k= "pcoa_3d"
    k <- fspaces_plot[pos_k]
    
    # dataframe with data to plot ----
    d_tr <- NULL
    d_sp_k <- NULL
    dev_k <- NULL
    qdev_k <- NULL
    df_plot_k <- data.frame(d_tr = df_dist[, "tr"],
                            d_sp_k = df_dist[, k],
                            dev_k = df_dev_dist[, k],
                            qdev_k = df_qdev_dist[, k])
    
    # names of Y axes and titles ----
    
    # name of Y axes depend on type of functional space: PCoa or dendrogram
    # of scaling of distance and type of deviation
    # printed only on left column = 1st space
    y_lab_dist_k <- ""
    y_lab_dev_k <- ""
    y_lab_qdev_k <- ""
    if (pos_k == 1){
      
      if (substr(k, 1, 4) == "pcoa"){
        y_lab_dist_k <- "Eucl. dist."
      }
      
      if (substr(k, 1, 4) == "tree") {
        y_lab_dist_k <- "Coph. dist."
      }
      
      y_lab_dev_k <- paste0("Dev. of ", y_lab_dist_k)
      y_lab_qdev_k <- paste0(nm_dev_qual[quality_metric], " ", y_lab_dist_k)
    }
    
    # title of plot = name of space
    # subtitle = name of metric + rounded value
    tit_k <- fspaces_nm_plot[k]
    subtit_k <- paste0(nm_qual_metrics[quality_metric], " = ",
                       round(as.numeric(fspaces_quality$quality_fspaces[k, quality_metric]), 3))
    
    # plotting trait-based distance versus raw distance in functional spaces ----
    plot_dist_k <- ggplot2::ggplot(data = df_plot_k, ggplot2::aes(x = d_tr, y = d_sp_k ) )+
      ggplot2::labs(x = NULL, y = y_lab_dist_k,
                    title = tit_k, subtitle= subtit_k ) +
      ggplot2::scale_x_continuous(limits = range_dist, expand = c(0,0) ) +
      ggplot2::scale_y_continuous(limits = range_dist, expand = c(0,0) ) +
      ggplot2::theme_bw(base_size = scaling_text) +
      ggplot2::theme(aspect.ratio = 1, plot.title = ggplot2::element_text(face = "bold"),
                     plot.margin = ggplot2::margin(2, 8, 2, 2, "pt") ) +
      ggplot2::geom_abline(ggplot2::aes(intercept=0, slope=1)) +
      
      ggplot2::geom_point(size = point_size, shape = 16, color="grey30", show.legend = FALSE)+
      ggplot2::guides(colour = "none")
    
    
    # plotting  plot for looking at raw deviation along values of distances ----
    plot_dev_k <- ggplot2::ggplot(data = df_plot_k, ggplot2::aes(x = d_tr, y = dev_k )) +
      ggplot2::labs(x = NULL, y = y_lab_dev_k) +
      ggplot2::scale_x_continuous(limits = range_dist, expand = c(0,0)) +
      ggplot2::scale_y_continuous(limits = range_dev,
                                  expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
      ggplot2::theme_bw(base_size = scaling_text) +
      ggplot2::theme(aspect.ratio = 1, plot.margin = ggplot2::margin(2, 8, 2, 2, "pt")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
      
      ggplot2::geom_point(size=point_size, shape = 16,
                          ggplot2::aes( colour = dev_k), alpha = 1)  +
      ggplot2::scale_colour_gradient2(low = gradient_deviation["neg"],
                                      high = gradient_deviation["pos"],
                                      mid = gradient_deviation["nul"],
                                      midpoint = 0,
                                      limits = range_dev,
                                      name = "Dev. raw" ) +
      ggplot2::guides(colour = "none")
    
    # plotting trait-based distance along transformed deviation ----
    # i.e. the one used for quality metric
    
    plot_qdev_k <- ggplot2::ggplot(data = df_plot_k, ggplot2::aes(x = d_tr, y = qdev_k )) +
      ggplot2::labs(x = x_lab,  y = y_lab_qdev_k) +
      ggplot2::scale_x_continuous(limits = range_dist, expand = c(0 ,0)) +
      ggplot2::scale_y_continuous(limits = range_qdev,
                                  expand = ggplot2::expansion(mult=c(0,0.03))) +
      ggplot2::theme_bw(base_size = scaling_text) +
      ggplot2::theme(aspect.ratio = 1, plot.margin = ggplot2::margin(2, 8, 2, 2, "pt")) +
      ggplot2::geom_point(size=point_size, shape = 16, ggplot2::aes( colour = qdev_k ))  +
      ggplot2::scale_colour_gradient(low = gradient_deviation_quality["low"],
                                     high = gradient_deviation_quality["high"],
                                     limits = range_qdev ,
                                     name = paste0("Dev. ", quality_metric )) +
      ggplot2::guides(colour = "none")
    
    
    # arranging panels ----
    
    # adding legends if last panel
    if (pos_k == length(fspaces_plot)) {
      plot_dev_k <- plot_dev_k + ggplot2::guides(colour = "colorbar")
      plot_qdev_k <- plot_qdev_k + ggplot2::guides(colour = "colorbar")
    }
    
    # merging the 3 plots in a single column
    col_plot_k <- (plot_dist_k / plot_dev_k / plot_qdev_k )
    
    # creating patchwork or merging with previous space(s) ----
    if (pos_k == 1) {
      patchwork_plots <- col_plot_k
    } else {
      patchwork_plots <- patchwork_plots | col_plot_k
    }
    
    
  }# and of loop on functional spaces ###
  
  
  # caption
  patchwork_plots_all <- patchwork_plots +
    patchwork::plot_annotation(caption = 'Made with mfd')
  
  # resolution and type of file
  device_file = "png"
  res_file = 300
  
  # displaying or saving
  if (is.null(name_file) == TRUE)  {
    patchwork_plots_all
  } else  {
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file),
                    plot = patchwork_plots_all,
                    device = device_file,
                    width = (0.25 + length(fspaces_plot)) * 700 / res_file,
                    height = (0.2 + 3) * 700 / res_file,
                    units = "in",
                    dpi = res_file)
  }
  
} # function end
