# Function to illustrate functional beta-diversity indices in a multidimensional space
#
# Authors:  Sébastien Villéger, Camille Magneville
#
#
# ------------------------------------------------------------------------------


#'Illustrate functional beta-diversity indices for pairs of assemblages in a
#'multidimensional space
#'
#'Illustrate overlap between convex hulls shaping species assemblages in a
#'multidimensional functional space.\strong{Before plotting beta functional
#'diversity indices should have been computed using the} \code{\link{beta.fd.multidim}}
#'\strong{function}.
#'
#'@param output_beta.fd.multidim the list returned by \code{\link{beta.fd.multidim}}
#'  when its input 'store_details ' is TRUE'. Thus, even if this function will
#'  illustrate functional beta-diversity for a single pair of assemblages, plots
#'  will be scaled according to all assemblages for which indices were computed.
#'
#'@param plot_asb_nm a \strong{vector} with names of the 2 assemblages for which
#'  functional beta-diversity will be illustrated.
#'
#'@param beta.family a \strong{character string} for the type of beta-diversity index for
#'  which values will be printed, 'Jaccard' (default) and/or 'Sorensen'.
#'
#'@param faxes a \strong{vector} with names of axes to plot (as columns names in
#'  \code{output_beta.fd.multidim$details$input$sp_faxes_coord} ). \strong{You
#'  can only plot from 2 to 4 axes for graphical reasons}. Default: faxes = NULL
#'  (the four first axes will be plotted).
#'
#'@param plot_sp_nm a \strong{vector} containing species names that are to be plotted.
#'  Default: plot_nm_sp = NULL (no name plotted).
#'
#'@param name_file a \strong{character string} with name of file to save the figure
#'  (without extension). Default is 'NULL' which means plot is displayed.
#'
#'@param faxes_nm a \strong{vector} with axes labels for figure. Default: as
#'  \code{faxes}).
#'
#'@param range_faxes a \strong{vector} with minimum and maximum values of axes. Note that
#'  to have a fair representation of position of species in all plots, they
#'  should have the same range. Default: faxes_lim = c(NA, NA) (the range is
#'  computed according to the range of values among all axes).
#'
#'@param color_bg a \strong{R color name or an hexadecimal cod} used to fill plot
#'  background. Default: color_bg = "grey95".
#'
#'@param shape_sp a \strong{vector} with 3 numeric values referring to the shape of
#'  symbol used for species from the 'pool' absent from the 2 assemblages, and
#'  for species present in the 2 assemblages ('asb1', and 'asb2'), respectively.
#'  Default: shape_sp = c("pool"=3, asb1=22, asb2=21) so cross, square and
#'  circle .
#'
#'@param size_sp a \strong{numeric value} referring to the size of symbols for species.
#'  Default: is size_sp = c("pool"=0.8, asb1=1, asb2=1).
#'
#'@param color_sp a \strong{vector} with 3 names or hexadecimal codes referring to the
#'  colour of symbol for species. Default is color_sp = c("pool"="grey50",
#'  asb1="blue", asb2= "red").
#'
#'@param fill_sp a \strong{vector} with 3 names or hexadecimal codes referring to the
#'  colour to fill symbol (if \code{shape_sp} >20) for species of the pool and
#'  of the 2 assemblages. Default is fill_sp = c("pool"=NA, asb1="white", asb2=
#'  "white").
#'
#'@param fill_vert a \strong{vector} with 3 names or hexadecimal codes referring to the
#'  colour to fill symbol (if \code{shape_sp} >20) for species being vertices of
#'  the convex hulls of the pool of species and of the 2 assemblages. Default is
#'  fill_vert = c("pool"=NA, asb1="blue", asb2= "red").
#'
#'@param color_ch a \strong{vector} with 3 names or hexadecimal codes referring to the
#'  border of the convex hulls of the pool of species and by the 2 assemblages.
#'  Default is color_ch = c("pool"=NA, asb1="blue", asb2= "red").
#'
#'@param fill_ch a \strong{vector} with 3 names or hexadecimal codes referring to the
#'  filling of the convex hull of the pool of species and of the 2 assemblages.
#'  Default is fill_ch = c("pool"="white", asb1="blue", asb2= "red").
#'
#'@param alpha_ch a \strong{vector} with 3 numeric value for transparency of the filling
#'  of the convex hulls (0 = high transparency, 1 = no transparency). Default is
#'  alpha_ch = c("pool"=1, asb1=0.3, asb2=0.3).
#'
#'@param nm_size a \strong{numeric value} for size of species label. Default is 3 points.
#'
#'@param nm_color a \strong{R color name or an hexadecimal code} referring to the colour
#'  of species label. Default is black.
#'
#'@param nm_fontface a \strong{character string} for font of species labels (e.g.
#'  "italic", "bold"). Default is 'plain'.
#'
#'@param check.input a \strong{logical value} defining whether key inputs (i.e. not
#'  aesthetics settings) are checked before plotting computation of indices.
#'  Possible error messages will thus may be more understandable for the user
#'  than R error messages. Default: check.input = TRUE.
#'
#'@return for the given pair of assemblages, returns a \code{patchwork} figure
#'  with overlap between convex hulls projected in 2-dimensional spaces. Values
#'  of functional beta-diversity indices are shown on top-right corner of the
#'  figure.
#'  
#'@examples
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
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$"details_fspaces"$"sp_pc_coord"
#'  # Get the occurrence dataframe:
#' asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = asb_sp_w_fruits) 
#' asb_sp_fruits_occ <- asb_sp_fruits_summ$"asb_sp_occ"
#' # Compute beta diversity indices:
#' beta_fd_fruits <- mFD::beta.fd.multidim(sp_faxes_coord_fruits[, 
#'  c("PC1", "PC2", "PC3", "PC4")], asb_sp_occ = asb_sp_fruits_occ,
#'  check.input = TRUE,
#'  beta.family = c("Jaccard"),
#'  store_details = TRUE)
#' # Compute beta fd plots:
#' beta.multidim.plot(output_beta.fd.multidim = beta_fd_fruits,
#'  plot_asb_nm = c("basket_1", "basket_3"),
#'  beta.family = c("Jaccard"),
#'  plot_sp_nm = c("apple", "cherry", "lemon"),
#'  faxes = paste0("PC", 1:4),
#'  name_file = NULL,
#'  faxes_nm = NULL, range_faxes = c(NA, NA),
#'  color_bg = "grey95",
#'  shape_sp = c("pool" = 3, asb1 = 22, asb2 = 21),
#'  size_sp = c("pool" = 0.8, asb1 = 1, asb2 = 1),
#'  color_sp = c("pool" = "grey50", asb1 = "blue", asb2 = "red"),
#'  fill_sp = c("pool" = NA, asb1 = "white", asb2 = "white"),
#'  fill_vert = c("pool" = NA, asb1 = "blue", asb2 = "red"),
#'  color_ch = c("pool" = NA, asb1 = "blue", asb2 = "red"),
#'  fill_ch = c("pool" = "white", asb1 = "blue", asb2 = "red"),
#'  alpha_ch = c("pool" = 1, asb1 = 0.3, asb2 = 0.3),
#'  nm_size = 3, nm_color = "black", nm_fontface = "plain",
#'  check.input = TRUE) 
#'
#'@export


beta.multidim.plot <- function(output_beta.fd.multidim,
                               plot_asb_nm,
                               beta.family,
                               plot_sp_nm = NULL,
                               faxes = NULL,
                               name_file = NULL,
                               faxes_nm = NULL, range_faxes = c(NA, NA),
                               color_bg = "grey95",
                               shape_sp = c("pool"=3, asb1=22, asb2=21),
                               size_sp = c("pool"=0.8, asb1=1, asb2=1),
                               color_sp = c("pool"="grey50", asb1="blue", asb2= "red"),
                               fill_sp = c("pool"= NA, asb1 ="white", asb2 = "white"),
                               fill_vert = c("pool" = NA, asb1 = "blue", asb2 = "red"),
                               color_ch = c("pool" = NA, asb1 = "blue", asb2 = "red"),
                               fill_ch = c("pool" = "white", asb1 = "blue", asb2 = "red"),
                               alpha_ch = c("pool" = 1, asb1 = 0.3, asb2 = 0.3),
                               nm_size = 3, nm_color = "black", nm_fontface = "plain",
                               check.input = TRUE) {
  
  
  ## extracting dataset from inputs ####
  
  # basic check of the core input
  if(any (names(output_beta.fd.multidim) != c("pairasb_fbd_indices", "details"))) {
    stop("Error: 'output_beta.fd.multidim' does not have elements of an output
      from 'beta.fd.multidim' function. Please check.")
  }
  
  # species occurrences and position in the functional space
  sp_faxes_coord <-output_beta.fd.multidim$details$inputs$sp_faxes_coord
  asb_sp_occ <- output_beta.fd.multidim$details$inputs$asb_sp_occ
  
  # indices values
  beta_fd <- output_beta.fd.multidim$pairasb_fbd_indices
  beta_nm <- names(beta_fd)
  
  
  ## check inputs if asked: ####
  if (check.input == TRUE) {
    
    
    if (length(plot_asb_nm) != 2) {
      stop("Error: there should be 2 assemblages names in 'plot_asb_nm'.
      Please check names of assemblages you want to plot.")
    }
    
    if (any(! plot_sp_nm %in% rownames(sp_faxes_coord))) {
      stop("Error: species names in 'plot_sp_nm' can not be found in
           'sp_faxes_coord' row names. Please check names of species you want to plot.")
    }
    
    if (! is.null(faxes)) {
      
      if (length(faxes) > 4) {
        stop("Error: Number of functional axes should be less than 4.
          Please change the number of functional axes to plot.")
      }
      
      if (! any(faxes %in% colnames(sp_faxes_coord))) {
        stop("Error: names of axes to plot can not be found in 'sp_faxes_coord'
               columns names. Please check names of axes you want to plot. ")
      }
      
    }
    
    if ((! is.null(faxes_nm)) & (length(faxes_nm) != length(faxes))) {
      stop("Error: Length of 'faxes_nm' should be equal to length of 'faxes'.
        Please check congruence between these inputs")
    }
    
    if (any(! plot_asb_nm %in% row.names(asb_sp_occ))) {
      stop("Error: assemblages names in 'plot_asb_nm' can not be found in
           'output_beta.fd.multidim'. Please check multidimensional functional
            beta-diversity has been computed on these assemblages.")
    }
    
    if (any(! tolower(substr(beta.family,1,3)) %in% substr(beta_nm,1,3))) {
      stop("Error: some functional beta-diversity indices are not present in
           'output_beta.fd.multidim' input. Please check indices family.")
    }
    
  }# end of checking inputs
  
  ## names and number of axes to plot ####
  
  asb_sp_occ <- as.data.frame(asb_sp_occ)
  
  # if no axes selected take up to the 4 first axes in the coordinates table
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(4, ncol(sp_faxes_coord))]
  }
  
  # if no faxes_nm provided, default names as in coordinates table
  if (! is.null(faxes) & is.null(faxes_nm)) {
    faxes_nm <- faxes
  }
  names(faxes_nm) <- faxes
  
  # number of axes
  nb_faxes <- length(faxes)
  
  # combinations of axes on plot
  axes_plot <- utils::combn(faxes, 2)
  plot_nb <- ncol(axes_plot)
  
  #### setting graphical parameters #####
  
  # range of axes
  if (is.na(range_faxes[1]) & is.na(range_faxes[2])) {
    range_sp_faxes_coord <- range(sp_faxes_coord)
    range_faxes <- range_sp_faxes_coord +
      c(-1, 1) * (range_sp_faxes_coord[2] - range_sp_faxes_coord[1]) * 0.05
  }
  
  # dataframe with species coordinates and option (vertices + label) if required
  sp_faxes_coord_plot <- data.frame(sp_faxes_coord,
                                    label = rep("", nrow(sp_faxes_coord)))
  
  # if some species names to be plotted, adding a binary variable to sp_faxes_coord
  if(! is.null(plot_sp_nm)) {
    sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
  }
  
  
  # computing convex hull of the species pool and its vertices if required
  if (! is.null(fill_ch["pool"])) {
    vert_pool_nm <- vertices(sp_faxes_coord, check.input = TRUE)
    vert_pool_coord <- sp_faxes_coord[vert_pool_nm, ]
  }
  
  # names of species present in species pool but not in assemblage
  sp_absent <- names(which(apply(asb_sp_occ[plot_asb_nm, ], 2, sum) == 0))
  
  # names of species present in each assemblage
  sp_asb1 <- names(asb_sp_occ[, which(asb_sp_occ[plot_asb_nm[1], ] == 1)])
  sp_asb2 <- names(asb_sp_occ[, which(asb_sp_occ[plot_asb_nm[2], ] == 1)])
  
  # vertices of each assemblage in the n-dimensional space
  vert_asb1 <- output_beta.fd.multidim$details$asb_vertices[[plot_asb_nm[1]]]
  vert_asb2 <- output_beta.fd.multidim$details$asb_vertices[[plot_asb_nm[2]]]
  
  
  
  ## plot for each pair of axes ####
  
  # list to store panels
  panels <- list()
  
  # loop on combinations
  for (k in (1:plot_nb)) {
    
    # names of axes as in input dataframe
    x <- axes_plot[1, k]
    y <- axes_plot[2, k]
    
    # background = axes defined by range of values and names as specified  ----
    plot_k <- ggplot2::ggplot(sp_faxes_coord_plot,
                              ggplot2::aes_string(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand=c(0,0)) +
      ggplot2::xlab(faxes_nm[x]) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand=c(0,0)) +
      ggplot2::ylab(faxes_nm[y]) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = color_bg)) +
      ggplot2::coord_fixed()
    
    # if required adding convex hull of pool projected in 2D ----
    if (! is.null(fill_ch["pool"])) {
      
      # vertices in 2D ordered for plotting
      vert_pool_xy <- vertices(sp_faxes_coord_plot[, c(x,y)], order_2D=TRUE,
                               check.input = TRUE )
      
      # plotting
      plot_k <- plot_k +
        ggplot2::geom_polygon(data = sp_faxes_coord_plot[vert_pool_xy, ],
                              ggplot2::aes_string(x = x, y = y ),
                              colour = color_ch["pool"],
                              fill = fill_ch["pool"], alpha = alpha_ch["pool"])
      
    }
    
    
    # points for species absent from the 2 assemblages
    plot_k <- plot_k +
      ggplot2::geom_point(data = sp_faxes_coord_plot[sp_absent, ],
                          ggplot2::aes_string(x = x, y = y),
                          colour = color_sp["pool"], fill = fill_sp["pool"],
                          shape = shape_sp["pool"], size = size_sp["pool"])
    
    
    # computing 2D convex hulls for the 2 assemblages
    vert_asb1_xy <- vertices(sp_faxes_coord_plot[sp_asb1, c(x,y)], order_2D = TRUE,
                             check.input = TRUE )
    vert_asb2_xy <- vertices(sp_faxes_coord_plot[sp_asb2, c(x,y)], order_2D = TRUE,
                             check.input = TRUE )
    
    # plotting convex hulls
    plot_k <- plot_k +
      ggplot2::geom_polygon(data = sp_faxes_coord_plot[vert_asb1_xy, ],
                            ggplot2::aes_string(x = x, y = y ),
                            colour = color_ch["asb1"],
                            fill = fill_ch["asb1"], alpha = alpha_ch["asb1"])+
      ggplot2::geom_polygon(data = sp_faxes_coord_plot[vert_asb2_xy, ],
                            ggplot2::aes_string(x = x, y = y),
                            colour = color_ch["asb2"],
                            fill = fill_ch["asb2"], alpha = alpha_ch["asb2"])
    
    # all species present in assemblages
    plot_k <- plot_k +
      ggplot2::geom_point(data = sp_faxes_coord_plot[sp_asb1, ],
                          ggplot2::aes_string(x = x, y = y),
                          colour = color_sp["asb1"], fill = fill_sp["asb1"],
                          shape = shape_sp["asb1"], size = size_sp["asb1"]) +
      ggplot2::geom_point(data = sp_faxes_coord_plot[sp_asb2, ],
                          ggplot2::aes_string(x = x, y = y ),
                          colour = color_sp["asb2"], fill = fill_sp["asb2"],
                          shape = shape_sp["asb2"], size = size_sp["asb2"])
    
    
    # filling points for species being vertices in the n-dimensional space
    plot_k <- plot_k +
      ggplot2::geom_point(data = sp_faxes_coord_plot[vert_asb1, ],
                          ggplot2::aes_string(x = x, y = y),
                          colour = color_sp["asb1"], fill = fill_vert["asb1"],
                          shape = shape_sp["asb1"], size = size_sp["asb1"]) +
      ggplot2::geom_point( data = sp_faxes_coord_plot[vert_asb2,],
                           ggplot2::aes_string(x = x, y = y ),
                           colour = color_sp["asb2"], fill = fill_vert["asb2"])
    
    
    # adding species names if needed ----
    if(! is.null(plot_sp_nm)) {
      plot_k <- plot_k +
        ggrepel::geom_text_repel(x = sp_faxes_coord_plot[, x], y = sp_faxes_coord_plot[, y],
                                 label = sp_faxes_coord_plot[, "label"],
                                 size = nm_size, colour= nm_color, fontface = nm_fontface,
                                 box.padding = grid::unit(2, 'lines'),
                                 force = 5,
                                 arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
                                 segment.color = nm_color)
    }
    
    # saving plot in a list
    panels[[k]] <- plot_k
    
  } # end of k
  
  
  ## plot for indices values #####
  
  # customizing position of texts in the plot
  spread_faxes <- (range_faxes[2]- range_faxes[1])
  hh <- c(1, 2.5, 4, 5.5)
  vv <- 0.3
  
  # plotting window
  plot_caption <- ggplot2::ggplot(sp_faxes_coord_plot) +
    ggplot2::scale_x_continuous(limits = range_faxes, expand=c(0, 0)) +
    ggplot2::scale_y_continuous(limits = range_faxes, expand=c(0, 0)) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
    ggplot2::geom_rect(xmin=range_faxes[1], xmax = range_faxes[2],
                       ymin = range_faxes[1], ymax = range_faxes[2],
                       fill = "white", colour ="black")
  
  # family and names of indices
  h <- NULL
  v <- NULL
  top <- NULL
  plot_caption <- plot_caption +
    ggplot2::geom_text(data = data.frame(
      h = range_faxes[1] + spread_faxes*0.15*hh,
      v = range_faxes[2] - spread_faxes*rep(0.2, 4),
      top = c("Family", "Dissimilarity =", "Turnover +", " Nestedness-res.")),
      ggplot2::aes(x = h, y = v, label = top ),
      size = 3, hjust = 0.5, fontface = "bold")
  
  # values of Jaccard indices if any
  if ("Jaccard" %in% beta.family ) {
    
    betajac_asb1_asb2 <- as.numeric(signif(beta_fd[
      which(beta_fd$asb.1 == plot_asb_nm[1] & beta_fd$asb.2 == plot_asb_nm[2]),
      which("jac" == substr(beta_nm, 1, 3))], 4))
    
    values_jac <- NULL
    values_sor <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text( data = data.frame(
        h = range_faxes[1] + spread_faxes*0.15*hh,
        v = range_faxes[2] - spread_faxes*rep(vv, 4),
        values_jac = c("Jaccard", betajac_asb1_asb2)),
        ggplot2::aes(x = h, y = v, label = values_jac ),
        size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
  }
  
  # values of Sorensen indices if any
  if ("Sorensen" %in% beta.family) {
    
    betasor_asb1_asb2 <- as.numeric( signif(beta_fd[
      which(beta_fd$asb.1 == plot_asb_nm[1] & beta_fd$asb.2 == plot_asb_nm[2]),
      which("sor" == substr(beta_nm, 1, 3))], 4))
    
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes*0.15*hh,
        v = range_faxes[2] - spread_faxes*rep(vv, 4),
        values_sor = c("Sorensen", betasor_asb1_asb2)),
        ggplot2::aes(x = h, y = v, label = values_sor),
        size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.2
  }
  
  # text about dimensionality
  nb <- NULL
  plot_caption<- plot_caption +
    ggplot2::geom_text(data = data.frame(
      h = range_faxes[1] + spread_faxes*0.1,
      v = range_faxes[2] - spread_faxes*vv,
      nb = paste0("NB: Indices were computed in a ",
                  ncol(sp_faxes_coord),"-dimensional space")),
      ggplot2::aes(x = h, y = v, label = nb),
      size = 3, hjust = 0, fontface = "italic")
  
  
  ### arranging panels #######
  
  # if 2 axes = 1 plot + caption
  if(plot_nb == 1) {
    patchwork_plots_all <- panels[[1]] + plot_caption +
      patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                             ncol = 2, nrow = 1, guides = "collect")
  }
  
  # if 3 axes = 3 plots + caption in a 2*2 layout
  if(plot_nb == 3) {
    patchwork_plots_all <- (panels[[1]] + plot_caption +
                              panels[[2]] + panels[[3]]) +
      patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                             ncol =2, nrow = 2, guides = "collect")
  }
  
  # if 4 axes = 6 plots + caption in a 3*3 layout with 2 empty cases
  if (plot_nb == 6) {
    patchwork_plots_all <- (panels[[1]] + patchwork::plot_spacer() + plot_caption +
                              panels[[2]] + panels[[3]] + patchwork::plot_spacer() +
                              panels[[4]] + panels[[5]] + panels[[6]]) +
      patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), widths = rep(1, 3),
                             ncol = 3, nrow = 3, guides = "collect")
  }
  
  # title and caption
  patchwork_plots_all <- patchwork_plots_all +
    patchwork::plot_annotation(title = paste0("Functional beta-diversity between '",
                                              plot_asb_nm[1],"' and '", plot_asb_nm[2], "'"),
                               caption = "made with mFD package")
  
  
  
  ## returning output ####
  
  # type, resolution and dimensions of file if to be saved
  device_file = "png"
  res_file = 300
  height_file <- 4*c(1, 2, 3) ; names(height_file) <- c("1", "3", "6")
  width_file <- 4*c(2, 2, 3) ; names(width_file) <- c("1", "3", "6")
  
  # displaying or saving
  if (is.null(name_file) == TRUE)  {
    patchwork_plots_all
  }
  else {
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file) ,
                    plot = patchwork_plots_all,
                    device = device_file,
                    scale = 1,
                    height= height_file[as.character(plot_nb)],
                    width = width_file[as.character(plot_nb)],
                    units= "in",
                    dpi = res_file)
  }
  
  
  # output = list of panels
  return(patchwork_plots_all)
  
}   # end of function ####
