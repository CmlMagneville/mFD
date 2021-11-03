#' Illustrate Functional beta-Diversity indices for pairs of assemblages in a
#' multidimensional space
#'
#' Illustrate overlap between convex hulls shaping species assemblages in a
#' multidimensional functional space.\strong{Before plotting beta functional
#' diversity indices should have been computed using the}
#' \code{\link{beta.fd.multidim}} \strong{function}.
#'
#' @param output_beta_fd_multidim the list returned by
#'   \code{\link{beta.fd.multidim}} when `details_returned = TRUE`.
#'   Thus, even if this function will illustrate functional beta-diversity for 
#'   a single pair of assemblages, plots will be scaled according to all
#'   assemblages for which indices were computed.
#'
#' @param plot_asb_nm a vector with names of the 2 assemblages for which
#'  functional beta-diversity will be illustrated.
#'
#' @param beta_family a character string for the type of beta-diversity index
#'   for which values will be printed, `'Jaccard'` (default) and/or
#'   `'Sorensen'`.
#'
#' @param faxes a vector with names of axes to plot (as columns names in
#'  \code{output_beta_fd_multidim$details$input$sp_faxes_coord} ). \strong{You
#'  can only plot from 2 to 4 axes for graphical reasons}. 
#'  Default: `faxes = NULL` (the four first axes will be plotted).
#'
#' @param plot_sp_nm a vector containing species names that are to be plotted.
#'  Default: `plot_nm_sp = NULL` (no name plotted).
#'
#' @param name_file a character string with name of file to save the figure
#'   (without extension). Default: `name_file = NULL` which means plot is
#'   displayed.
#'
#' @param faxes_nm a vector with axes labels for figure. Default: as
#'  \code{faxes}).
#'
#' @param range_faxes a vector with minimum and maximum values of axes. Note
#'   that to have a fair representation of position of species in all plots,
#'   they should have the same range. Default: `faxes_lim = c(NA, NA)` (the
#'   range is computed according to the range of values among all axes).
#'
#' @param color_bg a R color name or an hexadecimal code used to fill
#'   plot background. Default: `color_bg = "grey95"`.
#'
#'@param shape_sp a vector with 3 numeric values referring to the shape of
#'  symbol used for species from the 'pool' absent from the 2 assemblages, and
#'  for species present in the 2 assemblages ('asb1', and 'asb2'), 
#'  respectively. Default: `shape_sp = c(pool = 3, asb1 = 22, asb2 = 21)` so 
#'  cross, square and circle.
#'
#' @param size_sp a numeric value referring to the size of symbols for
#'  species. Default: `is size_sp = c(pool = 0.8, asb1 = 1, asb2 = 1)`.
#'
#' @param color_sp a vector with 3 names or hexadecimal codes referring to the
#'  colour of symbol for species. Default is:
#'  `color_sp = c(pool = "grey50", asb1 = "blue", asb2 = "red")`.
#'
#' @param fill_sp a vector with 3 names or hexadecimal codes referring to the
#'   color to fill symbol (if \code{shape_sp} > 20) for species of the pool and
#'   of the 2 assemblages. Default is:
#'   `fill_sp = c(pool = NA, asb1 = "white", asb2 = "white")`.
#'
#' @param fill_vert a vector with 3 names or hexadecimal codes
#'   referring to the colour to fill symbol (if \code{shape_sp} > 20) for
#'   species being vertices of the convex hulls of the pool of species and of
#'   the 2 assemblages. Default is:
#'   `fill_vert = c(pool = NA, asb1 = "blue", asb2 = "red")`.
#'
#' @param color_ch a vector with 3 names or hexadecimal codes referring to the
#'  border of the convex hulls of the pool of species and by the 2 assemblages.
#'  Default is: `color_ch = c(pool = NA, asb1 = "blue", asb2 = "red")`.
#'
#' @param fill_ch a vector with 3 names or hexadecimal codes referring to the
#'  filling of the convex hull of the pool of species and of the 2 assemblages.
#'  Default is `fill_ch = c(pool = "white", asb1 = "blue", asb2 = "red")`.
#'
#' @param alpha_ch a vector with 3 numeric value for transparency of
#'   the filling of the convex hulls (0 = high transparency, 1 = no
#'   transparency). Default is:
#'   `alpha_ch = c(pool = 1, asb1 = 0.3, asb2 = 0.3)`.
#'
#' @param nm_size a numeric value for size of species label. Default is `3`
#'   (in points).
#'
#' @param nm_color a R color name or an hexadecimal code referring to the color
#'  of species label. Default is `black.`
#'
#' @param nm_fontface a character string for font of species labels (e.g.
#'   "italic", "bold"). Default is `'plain'`.
#'
#' @param check_input a logical value indicating whether key features the 
#'   inputs are checked (e.g. class and/or mode of objects, names of rows 
#'   and/or columns, missing values). If an error is detected, a detailed 
#'   message is returned. Default: `check.input = TRUE`.
#'
#' @return If \code{name_file} is \code{NULL}, it returns a \code{patchwork} 
#' figure with overlap between convex hulls projected in 2-dimensional spaces
#' for the given pair of assemblages. Values of functional beta-diversity 
#' indices are shown on top-right corner of the figure. If \code{name_file} is 
#' not \code{NULL}, the plot is saved locally.
#'
#' @author Sebastien Villeger and Camille Magneville
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load Species*Traits dataframe:
#'  data("fruits_traits", package = "mFD")
#'
#' # Load Assemblages*Species dataframe:
#'  data("baskets_fruits_weights", package = "mFD")
#'
#' # Load Traits categories dataframe:
#'  data("fruits_traits_cat", package = "mFD")
#'
#' # Compute functional distance
#'  sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
#'                                    tr_cat        = fruits_traits_cat,
#'                                    metric        = "gower",
#'                                    scale_euclid  = "scale_center",
#'                                    ordinal_var   = "classic",
#'                                    weight_type   = "equal",
#'                                    stop_if_NA    = TRUE)
#'
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#'  fspaces_quality_fruits <- mFD::quality.fspaces(
#'                                   sp_dist             = sp_dist_fruits,
#'                                   maxdim_pcoa         = 10,
#'                                   deviation_weighting = "absolute",
#'                                   fdist_scaling       = FALSE,
#'                                   fdendro             = "average")
#'
#' # Retrieve species coordinates matrix:
#'  sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#'
#' # Get the occurrence dataframe:
#'  asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights)
#'  asb_sp_fruits_occ <- asb_sp_fruits_summ$"asb_sp_occ"
#'
#' # Compute beta diversity indices:
#'  beta_fd_fruits <- mFD::beta.fd.multidim(
#'   sp_faxes_coord   = sp_faxes_coord_fruits[, c("PC1", "PC2", "PC3", "PC4")],
#'   asb_sp_occ       = asb_sp_fruits_occ,
#'   check_input      = TRUE,
#'   beta_family      = c("Jaccard"),
#'   details_returned = TRUE)
#'
#' # Compute beta fd plots:
#'  beta.multidim.plot(
#'    output_beta_fd_multidim = beta_fd_fruits,
#'    plot_asb_nm             = c("basket_1", "basket_6"),
#'    beta_family             = c("Jaccard"),
#'    plot_sp_nm              = c("apple", "cherry", "lemon"),
#'    faxes                   = paste0("PC", 1:4),
#'    name_file               = NULL,
#'    faxes_nm                = NULL,
#'    range_faxes             = c(NA, NA),
#'    color_bg                = "grey95",
#'    shape_sp                = c(pool = 3, asb1 = 22, asb2 = 21),
#'    size_sp                 = c(pool = 0.8, asb1 = 1, asb2 = 1),
#'    color_sp                = c(pool = "grey50", asb1 = "blue",
#'                                asb2 = "red"),
#'    fill_sp                 = c(pool = NA, asb1 = "white", asb2 = "white"),
#'    fill_vert               = c(pool = NA, asb1 = "blue", asb2 = "red"),
#'    color_ch                = c(pool = NA, asb1 = "blue", asb2 = "red"),
#'    fill_ch                 = c(pool = "white", asb1 = "blue",
#'                                asb2 = "red"),
#'    alpha_ch                = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
#'    nm_size                 = 3,
#'    nm_color                = "black",
#'    nm_fontface             = "plain",
#'    check_input             = TRUE)
#'}


beta.multidim.plot <- function(output_beta_fd_multidim,
                               plot_asb_nm,
                               beta_family,
                               plot_sp_nm = NULL,
                               faxes = NULL,
                               name_file = NULL,
                               faxes_nm = NULL, range_faxes = c(NA, NA),
                               color_bg = "grey95",
                               shape_sp = c(pool = 3, asb1 = 22, 
                                            asb2 = 21),
                               size_sp = c(pool = 0.7, asb1 = 1.2, 
                                           asb2 = 1),
                               color_sp = c(pool = "grey50", 
                                            asb1 = "blue",
                                            asb2 = "red"),
                               fill_sp = c(pool = NA, asb1 = "white",
                                           asb2 = "white"),
                               fill_vert = c(pool = NA, asb1 = "blue",
                                             asb2 = "red"),
                               color_ch = c(pool = NA, asb1 = "blue",
                                            asb2 = "red"),
                               fill_ch = c(pool = "white", asb1 = "blue",
                                           asb2 = "red"),
                               alpha_ch = c(pool = 1, asb1 = 0.3, 
                                            asb2 = 0.3),
                               nm_size = 3, nm_color = "black",
                               nm_fontface = "plain",
                               check_input = TRUE) {
  
  
  # extract dataset from inputs ####
  
  # basic check of the core input:
  if (any(names(output_beta_fd_multidim) != c("pairasb_fbd_indices",
                                              "details"))) {
    stop("Argument 'output_beta_fd_multidim' does not have elements of an ", 
         "output from 'beta.fd.multidim()' function. Please check.")
  }
  
  # get species occurrences and position in the functional space:
  sp_faxes_coord <- output_beta_fd_multidim$details$inputs$sp_faxes_coord
  asb_sp_occ <- output_beta_fd_multidim$details$inputs$asb_sp_occ
  
  # get indices values:
  beta_fd <- output_beta_fd_multidim$pairasb_fbd_indices
  beta_nm <- names(beta_fd)
  beta_fd_df <- dist.to.df(beta_fd)
  
  
  # check_inputs if asked: ####
  if (check_input) {
    
    if (length(plot_asb_nm) != 2) {
      stop("There should be 2 assemblages names in 'plot_asb_nm'. Please ",
           "check names of assemblages you want to plot.")
    }
    
    if (any(!plot_sp_nm %in% rownames(sp_faxes_coord))) {
      stop("Species names in 'plot_sp_nm' can not be found in ",
           "'sp_faxes_coord' row names. Please check names of species you ",
           "want to plot.")
    }
    
    if (! is.null(faxes)) {
      
      if (length(faxes) > 4) {
        stop("Number of functional axes should be less than 4. Please change ", 
             "the number of functional axes to plot.")
      }
      
      if (any(! faxes %in% colnames(sp_faxes_coord))) {
        stop("Names of axes to plot can not be found in 'sp_faxes_coord' ",
             "columns names. Please check names of axes you want to plot.")
      }
    }
    
    if ((! is.null(faxes_nm)) && (length(faxes_nm) != length(faxes))) {
      stop("Length of 'faxes_nm' should be equal to length of 'faxes'. ",
           "Please check congruence between these inputs.")
    }
    
    if (any(!plot_asb_nm %in% row.names(asb_sp_occ))) {
      stop("Assemblages names in 'plot_asb_nm' can not be found in ",
           "'output_beta_fd_multidim'. Please check multidimensional ",
           "functional beta-diversity has been computed on these assemblages.")
    }
    
    if (any(!tolower(substr(beta_family, 1, 3)) %in% substr(beta_nm, 1, 3))) {
      stop("Some functional beta-diversity indices are not present in ",
           "'output_beta_fd_multidim' input. Please check indices family.")
    }
  } # end of check inputs
  
  # get names and number of axes to plot: ####
  
  # if no axes selected take up to the 4 first axes in the coordinates table:
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(4, ncol(sp_faxes_coord))]
  }
  
  # if no faxes_nm provided, default names as in coordinates table:
  if (!is.null(faxes) && is.null(faxes_nm)) {
    faxes_nm <- faxes
  }
  
  # give faxes names if faxes set to NULL:
  if (is.null(faxes_nm)) {
    faxes_nm <- faxes
  }
  
  names(faxes_nm) <- faxes
  
  # get number of axes:
  nb_faxes <- length(faxes)
  
  # retrieve combinations of axes on plot: 
  axes_plot <- utils::combn(faxes, 2)
  plot_nb   <- ncol(axes_plot)
  
  
  # set graphical parameters #####
  
  # range of axes:
  user_range <- "ok"
  range_sp_coord  <- range(sp_faxes_coord)
  
  if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
    user_range <- NA
    range_faxes <- range_sp_coord +
      c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
  }
  
  # check that the range is ok if the user chose it and convex hull:
  if (!is.na(user_range)) {
    
    if (range_faxes[1] > range_sp_coord[1]) {
      stop("The first value of 'range_faxes', is higher than minimum value ", 
           "of 'sp_faxes_coord' so the convex hull can not be plotted. ", 
           "Please change the minimal value of 'range_faxes' or set ", 
           "'plot_ch' to FALSE.")
    }
    
    if (range_faxes[2] < range_sp_coord[2]) {
      stop("The second value of 'range_faxes', is lower than maximum value ", 
           "of 'sp_faxes_coord' so the convex hull can not be plotted. ", 
           "Please change the maximal value of 'range_faxes' or set ", 
           "'plot_ch' to FALSE.")
    }
    
  } 

  # build dataframe with species coordinates and option (vertices + label) ...
  # ... if required:
  sp_faxes_coord_plot <- data.frame(sp_faxes_coord,
                                    label = "")
  
  # if some species names to be plotted, adding a character variable to ...
  # ... sp_faxes_coord:
  if (!is.null(plot_sp_nm)) {
    sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
  }
  
  
  # compute convex hull of the species pool and its vertices:
  vert_pool<- output_beta_fd_multidim$"details"$"pool_vertices"  
  
  # get names of species present in each assemblage:
  sp_asb1 <- names(which(asb_sp_occ[plot_asb_nm[1], ] == 1))
  sp_asb2 <- names(which(asb_sp_occ[plot_asb_nm[2], ] == 1))
  
  # retrieve vertices of each assemblage in the n-dimensional space:
  vert_asb1 <- output_beta_fd_multidim$details$asb_vertices[[plot_asb_nm[1]]]
  vert_asb2 <- output_beta_fd_multidim$details$asb_vertices[[plot_asb_nm[2]]]
  
  
  # plot for each pair of axes ####
  
  # list to store panels:
  panels <- list()
  
  # loop on combinations:
  for (k in (1:plot_nb)) {
    
    # get species coordinates along the 2 axes:
    sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, axes_plot[1:2, k] ])
    colnames(sp_coord_xy) <- c("x", "y")
    
    # get a list with dataframe of data to plot:
    asb_sp_coord2D_k <- list(asb1 = sp_coord_xy[sp_asb1, ],
                             asb2 = sp_coord_xy[sp_asb2, ])
    vertices_nD_k <- list(asb1 = vert_asb1, asb2 = vert_asb2)
    
    
    # plot background = axes defined by range of values and names as specified:
    plot_k <- background.plot(range_faxes, faxes_nm, color_bg)
    
    # species pool
    plot_k <- pool.plot(ggplot_bg = plot_k,
                        sp_coord2D = sp_coord_xy,
                        vertices_nD = vert_pool,
                        plot_pool = TRUE,
                        color_ch = color_ch["pool"],
                        fill_ch = fill_ch["pool"],
                        alpha_ch = alpha_ch["pool"],
                        shape_pool = shape_sp["pool"],
                        size_pool = size_sp["pool"],
                        color_pool = color_sp["pool"],
                        fill_pool = fill_sp["pool"],
                        shape_vert = shape_sp["pool"],
                        size_vert = size_sp["pool"],
                        color_vert = color_sp["pool"],
                        fill_vert = fill_sp["pool"])
    
    # plot 2D convex hulls and points for the 2 assemblages:
    plot_k <- fric.plot(ggplot_bg = plot_k,
                        asb_sp_coord2D = asb_sp_coord2D_k,
                        asb_vertices_nD = vertices_nD_k,
                        plot_sp = TRUE,
                        color_ch = color_ch[c("asb1", "asb2")],
                        fill_ch = fill_ch[c("asb1", "asb2")],
                        alpha_ch = alpha_ch[c("asb1", "asb2")],
                        shape_sp = shape_sp[c("asb1", "asb2")],
                        size_sp = size_sp[c("asb1", "asb2")],
                        color_sp = color_sp[c("asb1", "asb2")],
                        fill_sp = fill_sp[c("asb1", "asb2")],
                        shape_vert = shape_sp[c("asb1", "asb2")],
                        size_vert = size_sp[c("asb1", "asb2")],
                        color_vert = color_sp[c("asb1", "asb2")],
                        fill_vert = fill_vert[c("asb1", "asb2")])
    
    
    # add species names if needed:
    if (! is.null(plot_sp_nm)) {
      x <- NULL
      y <- NULL
      plot_k <- plot_k +
        ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                 ggplot2::aes_string(x = axes_plot[1, k],
                                                     y = axes_plot[2, k],
                                                     label = "label"),
                                 size = nm_size, colour = nm_color,
                                 fontface = nm_fontface,
                                 max.overlaps = Inf,
                                 box.padding = grid::unit(2, 'lines'),
                                 force = 5,
                                 arrow = grid::arrow(
                                   length = grid::unit(0.02, 'npc')),
                                 segment.color = nm_color)
    }
    
    # save plot in a list:
    panels[[k]] <- plot_k
    
  } # end of k
  
  
  # plot for indices values #####
  
  # customize position of texts in the plot:
  spread_faxes <- (range_faxes[2] - range_faxes[1])
  hh <- c(1, 2.5, 4, 5.5)
  vv <- 0.3
  
  # plot window:
  x <- NULL
  y <- NULL
  plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                  ggplot2::aes(x = x, y = y)) +
    ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
    ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                       ymin = range_faxes[1], ymax = range_faxes[2],
                       fill = "white", colour ="black")
  
  # plot family and names of indices:
  h   <- NULL
  v   <- NULL
  top <- NULL
  plot_caption <- plot_caption +
    ggplot2::geom_text(data = data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh,
      v = range_faxes[2] - spread_faxes * rep(0.2, 4),
      top = c("Family", "Dissimilarity =", "Turnover +", " Nestedness-res.")),
      ggplot2::aes(x = h, y = v, label = top),
      size = 3, hjust = 0.5, fontface = "bold")
  
  # plot values of Jaccard indices if any:
  if ("Jaccard" %in% beta_family) {
    
    betajac_asb1_asb2 <- beta_fd_df[
      which(beta_fd_df$x1 == plot_asb_nm[1] & 
              beta_fd_df$x2 == plot_asb_nm[2]), ]
    betajac_asb1_asb2 <- as.numeric(round(betajac_asb1_asb2[, -c(1, 2)], 4))
    
    values_jac <- NULL
    values_sor <- NULL
    label <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh,
        v = range_faxes[2] - spread_faxes*rep(vv, 4),
        values_jac = c("Jaccard", betajac_asb1_asb2)),
        ggplot2::aes(x = h, y = v, label = values_jac),
        size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
  }
  
  # plot values of Sorensen indices if any:
  if ("Sorensen" %in% beta_family) {
    
    betasor_asb1_asb2 <- beta_fd_df[
      which(beta_fd_df$x1 == plot_asb_nm[1] & 
              beta_fd_df$x2 == plot_asb_nm[2]), ]
    betasor_asb1_asb2 <- as.numeric(round(betasor_asb1_asb2[, -c(1, 2)], 4))
    
    label <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes*0.15*hh,
        v = range_faxes[2] - spread_faxes*rep(vv, 4),
        values_sor = c("Sorensen", betasor_asb1_asb2)),
        ggplot2::aes(x = h, y = v, label = values_sor),
        size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.2
  }
  
  # add text about dimensionality:
  nb <- NULL
  label <- NULL
  x <- NULL
  y <- NULL
  plot_caption <- plot_caption +
    ggplot2::geom_text(data = data.frame(
      h = range_faxes[1] + spread_faxes * 0.1,
      v = range_faxes[2] - spread_faxes * vv,
      nb = paste0("NB: Indices were computed in a ",
                  ncol(sp_faxes_coord),"-dimensional space")),
      ggplot2::aes(x = h, y = v, label = nb),
      size = 3, hjust = 0, fontface = "italic")
  
  # add legend (convex hull, asb species and pool species):
  
  ## plot legend:
  values_lab <- NULL
  
  ### for 1st asb:
  plot_caption <- plot_caption +
    ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
                       xmax = range_faxes[1] + spread_faxes*0.15,
                       ymin = range_faxes[2] - spread_faxes*0.51,
                       ymax = range_faxes[2] - spread_faxes*0.55,
                       fill = color_sp[["asb1"]], alpha = alpha_ch[["asb1"]]) + 
    
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.45,
                       y = range_faxes[2] - spread_faxes*0.525,
                       label = paste0("convex hull of", sep = " ", 
                                      plot_asb_nm[1]),
                       colour = color_sp[["asb1"]], size = 3) + 
    
    ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                        y = range_faxes[2] - spread_faxes*0.58,
                        size = size_sp[["asb1"]], shape = shape_sp[["asb1"]],
                        color = color_sp[["asb1"]], fill = fill_sp[["asb1"]]) + 
    
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.45,
                       y = range_faxes[2] - spread_faxes*0.58,
                       label = paste0("shape of species from", sep = " ", 
                                      plot_asb_nm[1]),
                       colour = color_sp[["asb1"]], size = 3) 
  
  
  ### asb2:
  
  plot_caption <- plot_caption +
    ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
                       xmax = range_faxes[1] + spread_faxes*0.15,
                       ymin = range_faxes[2] - spread_faxes*0.64,
                       ymax = range_faxes[2] - spread_faxes*0.68,
                       fill = color_sp[["asb2"]], 
                       alpha = alpha_ch[["asb2"]]) + 
    
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.45,
                       y = range_faxes[2] - spread_faxes*0.665,
                       label = paste0("convex hull of", sep = " ", 
                                      plot_asb_nm[2]),
                       colour = color_sp[["asb2"]], size = 3) + 
    
    ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                        y = range_faxes[2] - spread_faxes*0.71,
                        size = size_sp[["asb2"]], 
                        shape = shape_sp[["asb2"]],
                        color = color_sp[["asb2"]], 
                        fill = fill_sp[["asb2"]]) + 
    
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.45,
                       y = range_faxes[2] - spread_faxes*0.71,
                       label = paste0("shape of species from", sep = " ", 
                                      plot_asb_nm[2]),
                       colour = color_sp[["asb2"]], size = 3)
  
  ### for global pool:
  
  plot_caption <- plot_caption +
    
    ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                        y = range_faxes[2] - spread_faxes*0.77,
                        size = size_sp[["pool"]], 
                        shape = shape_sp[["pool"]],
                        color = color_sp[["pool"]], 
                        fill = fill_sp[["pool"]]) + 
    
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.45,
                       y = range_faxes[2] - spread_faxes*0.77,
                       label = "shape of species from the global pool",
                       colour = color_sp[["pool"]], size = 3)
  
  
  # arrange panels #######
  
  # merge panels and caption:
  patchwork_plots_all <- panels.to.patchwork(panels, plot_caption)
  
  # add title and caption:
  patchwork_plots_all <- patchwork_plots_all +
    patchwork::plot_annotation(title = paste0(
      "Functional beta-diversity between '", plot_asb_nm[1], "' and '",
      plot_asb_nm[2], "'"),
      caption = "made with mFD package")
  
  
  
  ##return output ####
  
  # type, resolution and dimensions of file if to be saved
  device_file <- "png"
  res_file <- 300
  height_file <- 4 * c(1, 2, 3)
  names(height_file) <- c("1", "3", "6")
  width_file  <- 4 * c(2, 2, 3)
  names(width_file) <- c("1", "3", "6")
  
  # save in a file or returning list of ggplots:
  if (! is.null(name_file))  {
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file) ,
                    plot = patchwork_plots_all,
                    device = device_file,
                    scale = 1,
                    height= height_file[as.character(plot_nb)],
                    width = width_file[as.character(plot_nb)],
                    units= "in",
                    dpi = res_file)
  } else {
    
    # output = list of panels and caption + patchwork:
    output <- NULL
    output <- panels
    names(output) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
    output[["caption"]] <- plot_caption
    output[["patchwork"]] <- patchwork_plots_all
    
    return(output)
  }
}
