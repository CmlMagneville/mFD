#' Plot Functional Space and Chosen Functional Indices
#'
#' Compute a graphical representation of functional indices. \strong{To plot
#' functional indices, functional indices values must have been retrieve through
#' the use of the} \code{\link{alpha.fd.multidim}} \strong{function}.
#'
#' @param output_alpha_fd_multidim a list of objects retrieved through the 
#' \code{\link{alpha.fd.multidim}} function.
#' 
#' @param plot_asb_nm a vector containing name(s) of assemblage(s) to plot 
#'
#' @param ind_nm a vector of character string of the name of functional
#'   indices to plot. \strong{Indices names must be written in lower case
#'   letters}. Possible indices to compute are: "fdis", "feve", "fric", "fdiv",
#'   "fori" and "fspe". Default: all the indices are computed.
#'
#' @param faxes a vector with names of axes to plot. \strong{You can only plot
#'  from 2 to 4 axes for graphical reasons: vector length should be between 2
#'  and 4}. Default: faxes = NULL (the four first axes will be plotted).
#'
#' @param faxes_nm a vector with axes labels if the user want different axes
#'  labels than \code{faxes} ones. Default: faxes_nm = faxes (labels will the
#'  the same that \code{faxes} ones).
#'
#' @param range_faxes a vector with minimum and maximum for values for axes.
#'  Note that to have a fair representation of position of species in all plots,
#'  all axes must have the same range. Default: faxes_lim = c(NA, NA) (the range
#'  is computed according to the range of values among all axes, all axes having
#'  the same range).
#'
#' @param color_bg a R color name  or an hexadecimal code used to fill plot
#'  background. Default: `color_bg = "grey95"`.
#'
#' @param size_sp a vector gathering numeric values referring to the size of 
#' species belonging to the global pool and the plotted assemblage(s).
#' It should be written  as c(pool = "...", asb1 = "...", ...). 
#' 
#' @param size_sp_nm a numeric value referring to the size of species names 
#' if plotted. 
#'
#' @param color_sp a vector gathering R color names or hexadecimal codes 
#' referring to the color of species from the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#' 
#' @param color_vert a vector gathering R color names or hexadecimal codes 
#' referring to the color of vertices from the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#' 
#' @param color_centroid_fspe a vector gathering R color name or 
#' hexadecimal code used to draw FSpe centroid (i.e. center of the 
#' functional space) color. 
#'
#' @param color_ch a vector gathering R color names or hexadecimal codes 
#' referring to the color of the convex pool of the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#' 
#' @param color_sp_nm a R color name or hexadecimal code referring to the 
#' color of names of species if plotted.
#'
#' @param fill_sp a vector gathering R color names or hexadecimal codes 
#' referring to the filled color of species from the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#' 
#' @param fill_vert a vector gathering R color names or hexadecimal codes 
#' referring to the filled color of vertices from the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#'
#' @param fill_ch a vector gathering R color names or hexadecimal codes 
#' referring to the color to fill the convex pool of the global pool and studied 
#' assemblage(s). It should be written  as c(pool = "...", asb1 = "...", ...). 
#'
#' @param alpha_ch a vector gathering numeric values referring to the opacity of 
#' convex hulls of the global pool and the plotted assemblage(s).
#' It should be written  as c(pool = "...", asb1 = "...", ...).
#' (0 = high transparency, 1 = no transparency). 
#'
#' @param shape_sp a vector gathering numeric values referring to the symbol 
#' used to draw species from the global pool and the plotted assemblage(s).
#' It should be written  as c(pool = "...", asb1 = "...", ...).
#' (0 = high transparency, 1 = no transparency). 
#'
#' @param shape_centroid_fdis a vector gathering numeric value(s) used to draw 
#' FDis centroid size. 
#'  
#' @param shape_centroid_fdiv a vector gathering numeric value(s) used to draw 
#' FDiv centroid size. 
#' 
#' @param shape_centroid_fspe a vector gathering numeric value used to draw 
#' FSpe centroid (i.e. center of the functional space) size. 
#'
#' @param plot_sp_nm a vector containing species names that are to be plotted.
#'  Default: `plot_nm_sp = NULL` (no name plotted).
#'
#' @param fontface_sp_nm a character string for font of species labels (e.g.
#'  "italic", "bold"). Default: `fontface_sp_nm = 'plain'`.
#'
#' @param save_file a logical value telling if plots should be locally 
#' saved or not.
#'
#' @param check_input a logical value indicating whether key features the inputs
#'   are checked (e.g. class and/or mode of objects, names of rows and/or
#'   columns, missing values). If an error is detected, a detailed message is
#'   returned. Default: `check.input = TRUE`.
#'
#' @return For the given assemblage, return a list of one \code{patchwork}
#'   figure per functional indice containing plots for combinations of up to
#'   four axes.
#'
#' @export
#'
#' @author Camille Magneville and Sébastien Villéger
#'
#' @examples
#' \dontrun{
#' # Load Species*Traits dataframe:
#' data("fruits_traits", package = "mFD")
#'
#' # Load Assemblages*Species dataframe:
#' data("baskets_fruits_weights", package = "mFD")
#'
#' # Load Traits categories dataframe:
#' data("fruits_traits_cat", package = "mFD")
#'
#' # Compute functional distance
#' sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
#'                                   tr_cat        = fruits_traits_cat,
#'                                   metric        = "gower",
#'                                   scale_euclid  = "scale_center",
#'                                   ordinal_var   = "classic",
#'                                   weight_type   = "equal",
#'                                   stop_if_NA    = TRUE)
#'
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#' fspaces_quality_fruits <- mFD::quality.fspaces(sp_dist = sp_dist_fruits,
#'  maxdim_pcoa         = 10,
#'  deviation_weighting = "absolute",
#'  fdist_scaling       = FALSE,
#'  fdendro             = "average")
#'
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#'
#' # Compute alpha diversity indices:
#' alpha_fd_indices_fruits <- mFD::alpha.fd.multidim(
#'   sp_faxes_coord   = sp_faxes_coord_fruits[ , c("PC1", "PC2", "PC3", "PC4")],
#'   asb_sp_w         = baskets_fruits_weights,
#'   ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv",
#'                        "fori", "fspe"),
#'   scaling          = TRUE,
#'   check_input      = TRUE,
#'   details_returned = TRUE)
#'
#' # Retrieve alpha diversity indices table:
#' fd_ind_values_fruits <- alpha_fd_indices_fruits$functional_diversity_indices
#' fd_ind_values_fruits
#' }

alpha.multidim.plot <- function(output_alpha_fd_multidim,
                                plot_asb_nm,
                                ind_nm              = c("fdis", "fnnd", "feve", 
                                                        "fric", "fdiv", "fori", 
                                                        "fspe"),
                                faxes               = NULL,
                                faxes_nm            = NULL,
                                range_faxes         = c(NA, NA),
                                color_bg            = "grey95",
                                shape_sp            = c(pool = 3, asb1 = 21, 
                                                        asb2 = 21),
                                size_sp             = c(pool = 0.7, asb1 = 1, 
                                                        asb2 = 1),
                                color_sp            = c(pool = "grey50", 
                                                        asb1 = "#0072B2",
                                                        asb2 = "#D55E00"),
                                color_vert          = c(pool = "grey50", 
                                                        asb1 = "#0072B2",
                                                        asb2 = "#D55E00"),
                                fill_sp             = c(pool =  NA, 
                                                        asb1 = "#FFFFFF30",
                                                        asb2 = "#FFFFFF30"),
                                fill_vert           = c(pool = NA, 
                                                        asb1 = "#0072B2",
                                                        asb2 = "#D55E00"),
                                color_ch            = c(pool = NA, 
                                                        asb1 = "#0072B2",
                                                        asb2 = "#D55E00"),
                                fill_ch             = c(pool = "white", 
                                                        asb1 = "#0072B2",
                                                        asb2 = "#D55E00"),
                                alpha_ch            = c(pool = 1, asb1 = 0.3, 
                                                        asb2 = 0.3),
                                shape_centroid_fdis = c(asb1 = 22,  asb2 = 22),
                                shape_centroid_fdiv = c(asb1 = 24,  asb2 = 25),
                                shape_centroid_fspe = 23,
                                color_centroid_fspe = "black",
                                size_sp_nm          = 3, 
                                color_sp_nm         = "black",
                                plot_sp_nm          = NULL,
                                fontface_sp_nm      = "plain",
                                save_file           = FALSE,
                                check_input         = TRUE) {
  
  
  
  
  # compulsory check that indices values and details available in main input:
  if (! identical(names(output_alpha_fd_multidim),
                   c("functional_diversity_indices","details"))) {
    
    stop("Input 'output_alpha_fd_multidim' should be the output of the
           'alpha/multidim' function run with 'details=TRUE'.")
    
  }
  
  # shorten names of main input: ####
  asb_fd_ind <- output_alpha_fd_multidim$functional_diversity_indices
  fd_details <- output_alpha_fd_multidim$details
  
  
  # Check inputs relative to this function ...
  # ... (funct.space.plot, already done) ####
  
  if (check_input) {
    
    # check that all functional diversity indices have the right names:
    if (any(! ind_nm %in%
        c("fide", "fdis", "fnnd", "feve", "fric", "fdiv", "fori", "fspe"))) {
      
      stop("Names of functional diversity indices in 'ind_nm' are not well written.",
           "Please re-write them. Be careful, they should all be written in ",
           "lowercase letters.")
      
    }
    
    # check that indices to plot are contained in the fd_ind_value table:
    if (any(! ind_nm %in% colnames(asb_fd_ind) ) == TRUE) {
      stop("Error: Functional diversity indices to plot must be contained in
           'fd_ind_values' columns")
    }
    
    
    # check that good number of assemblage(s) to plot:
    if (!length(plot_asb_nm) %in% c(1, 2)) {
      stop("This function can only plot one or two assemblages. Please chose ",
           "two or less assemblages to plot.")
    }
    
    # check that assemblage(s) to plot has(ve) the right name(s):
    if ( any( ! plot_asb_nm %in% row.names(asb_fd_ind)) ) {
      stop("Name(s) of assemblage(s) to plot is(are) provided in 'plot_asb_nm' do not match those in 'output_alpha_fd_multidim$asb_fd_ind. Please re-write.")
    }
    
  } # end of check input
  
  
  # Prepare data for plotting ####
  
  # get coordinates of species:
  sp_faxes_coord <- fd_details$sp_faxes_coord
  
  # create a list ot store outputs:
  list_panels <-list()
  
  # get names of assemblages:
  pool <- "pool"
  asb1 <- plot_asb_nm[1]
  nm_asb <- asb1
  two_asb <- FALSE
  
  if (length(plot_asb_nm) == 2) {
    two_asb <- TRUE
    asb2 <- plot_asb_nm[2]
    nm_asb <- paste(nm_asb, asb2, sep = "_")
  }
  
  
  # set type, resolution and dimensions of file if to be saved:
  device_file <- "png"
  res_file <- 300
  height_file <- 8 * c(2, 2, 3)
  names(height_file) <- c("1", "3", "6")
  width_file  <- 10 * c(2, 2, 3)
  names(width_file) <- c("1", "3", "6")
  
  # get number of dimensions in input:
  nb_dim <- ncol(sp_faxes_coord)
  
  # give faxes identity if faxes set to NULL:
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
  }
  
  # give faxes names if faxes set to NULL:
  if (is.null(faxes_nm)) {
    faxes_nm <- faxes
  }
  names(faxes_nm) <- faxes
  
  # get number of axes:
  nb_faxes <- length(faxes)
  
  # get combinations of axes on plot:
  axes_plot <- utils::combn(faxes, 2)
  plot_nb   <- ncol(axes_plot)
  
  # set range of axes if c(NA, NA):
  if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
    range_sp_coord  <- range(sp_faxes_coord)
    range_faxes <- range_sp_coord +
      c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
  }

  # create a dataframe with species coordinates and option (vertices + label)...
  # ... if required:
  sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")
  
  # if some species names to be plotted, adding a character variable to ...
  # ... sp_faxes_coord:
  if (! is.null(plot_sp_nm)) {
    sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
  }
  
  # get vertices of the convex hull of the species pool:
  vert_pool <- fd_details$pool_vert_nm
  
  # retrieve names and weights of species present in each assemblage:
  sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
  if (two_asb){
    sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))
  }
  
  
  # plot FRic if required: ####
  if ("fric" %in% ind_nm) {
    
    # list to store ggplot
    panels_fric <- list()
    
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # names of axes
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      vertices_nD_k <- list()
      vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
      
      if (two_asb) {
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool:
      plot_k <- pool.plot(ggplot_bg = plot_k,
                          sp_coord2D = sp_coord_xy,
                          vertices_nD = vert_pool,
                          plot_pool = TRUE,
                          color_ch = color_ch[["pool"]],
                          fill_ch = fill_ch[["pool"]],
                          alpha_ch = alpha_ch[["pool"]],
                          shape_pool = shape_sp[["pool"]],
                          size_pool = size_sp[["pool"]],
                          color_pool = color_sp[["pool"]],
                          fill_pool = fill_sp[["pool"]],
                          shape_vert = shape_sp[["pool"]],
                          size_vert = size_sp[["pool"]],
                          color_vert = color_vert[["pool"]],
                          fill_vert = fill_vert[["pool"]])
      
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
                          color_vert = color_vert[c("asb1", "asb2")],
                          fill_vert = fill_vert[c("asb1", "asb2")])
      
      # add species names if needed:
      if (! is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                      y = xy_k[2],
                                                      label = "label"),
                                   size = size_sp_nm, colour= color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fric[[k]] <- plot_k
      
    } # end of k
    
    
    # plot for indices values: #####
    
    # retrieve values to plot:
    top_fric <- c("Functional richness", asb1, "")
    values_fric <- c(round(asb_fd_ind[asb1, "fric"], 3), "")
    if (two_asb) {
      top_fric[3] <- asb2
      values_fric[2] <- round(asb_fd_ind[asb2,"fric"], 3)
    }
    
    # customize position of texts in the plot
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    
    # plot window:
    x <- NULL
    y <- NULL
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, 
                                               y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fric),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FRic values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FRic", values_fric))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
        ggplot2::aes(x = h, y = v, label = values_lab),
        size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
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

      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.525,
                         label = paste0("convex hull of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.58,
                          size = size_sp[["asb1"]], shape = shape_sp[["asb1"]],
                          color = color_sp[["asb1"]], fill = fill_sp[["asb1"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("shape of species from", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
                           xmax = range_faxes[1] + spread_faxes*0.15,
                           ymin = range_faxes[2] - spread_faxes*0.64,
                           ymax = range_faxes[2] - spread_faxes*0.68,
                           fill = color_sp[["asb2"]], 
                           alpha = alpha_ch[["asb2"]]) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("convex hull of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.71,
                            size = size_sp[["asb2"]], 
                            shape = shape_sp[["asb2"]],
                            color = color_sp[["asb2"]], 
                            fill = fill_sp[["asb2"]]) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.71,
                           label = paste0("shape of species from", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }

    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.77,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.77,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fric <- panels.to.patchwork(panels_fric, plot_caption)
    
    # title and caption
    tit_fric <- paste0("Functional Richness of '", asb1, "'")
    if (two_asb){
      tit_fric<-paste0(tit_fric, " and '", asb2,"'")
    }
    
    # create final patchwork object: 
    patchwork_fric <- patchwork_fric +
      patchwork::plot_annotation(title = tit_fric,
                                 caption = "made with mFD package")
    
    
    # saving as file or list
    if (save_file == TRUE) {
      
      # name of file built with assemblage names and number of dimensions
      file_fric <- paste0(nm_asb, "_", "FRic_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fric ,
                      plot = patchwork_fric,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units= "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fric) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fric[["caption"]] <- plot_caption
      panels_fric[["patchwork"]] <- patchwork_fric
      list_panels[["fric"]] <- panels_fric
      
    }
    
  } # end of plot FRic
  
  
  ##############################################################################
  
  
  # FDiv plot if required ####
  if ("fdiv" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fdiv <- list()
    
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes
      xy_k <- axes_plot[1:2, k]
      
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy)<-c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      vertices_nD_k <- list()
      vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
      asb_vertG_coord2D_k <- list()
      asb_vertG_coord2D_k[["asb1"]] <- fd_details$asb_G_coord[[asb1]][xy_k]
      asb_meanDtoG_k <- list()
      asb_meanDtoG_k[["asb1"]] <- fd_details$asb_mean_dist_G[[asb1]]
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
        asb_vertG_coord2D_k[["asb2"]] <- fd_details$asb_G_coord[[asb2]][xy_k]
        asb_meanDtoG_k[["asb2"]] <- fd_details$asb_mean_dist_G[[asb2]]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool:
      plot_k <- pool.plot(ggplot_bg = plot_k,
                          sp_coord2D = sp_coord_xy,
                          vertices_nD = vert_pool,
                          plot_pool = TRUE,
                          color_pool = color_ch["pool"],
                          fill_pool = fill_ch["pool"],
                          alpha_ch = alpha_ch["pool"],
                          shape_pool = shape_sp["pool"],
                          size_pool = size_sp["pool"],
                          shape_vert = shape_sp["pool"],
                          size_vert = size_sp["pool"],
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      # plot 2D convex hulls and points for the 2 assemblages:
      plot_k <- fdiv.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          asb_vertices_nD = vertices_nD_k,
                          asb_vertG_coord2D = asb_vertG_coord2D_k,
                          asb_meanDtoG = asb_meanDtoG_k,
                          plot_sp = TRUE,
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = fill_sp[c("asb1", "asb2")],
                          shape_vert = shape_sp[c("asb1", "asb2")],
                          color_vert = color_vert[c("asb1", "asb2")],
                          fill_vert = fill_vert[c("asb1", "asb2")],
                          shape_vertG = c(asb1 = 23, asb2 = 24),
                          size_vertG = c(asb1 = 3, asb2 = 3),
                          color_vertG = color_sp[c("asb1", "asb2")],
                          fill_vertG = color_sp[c("asb1", "asb2")],
                          color_meanD = color_sp[c("asb1", "asb2")])

      
      # add species names if needed:
      if (! is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                      y = xy_k[2],
                                                      label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # saving plot in a list
      panels_fdiv[[k]] <- plot_k
      
    } # end of k
    
    
    
    # plot caption ####
    
    # retrieve values to plot:
    top_fdiv <- c("Functional Divergence", asb1, "")
    values_fdiv <- c(round(asb_fd_ind[asb1, "fdiv"], 3), "")
    if (two_asb) {
      top_fdiv[3] <- asb2
      values_fdiv[2] <- round(asb_fd_ind[asb2,"fdiv"], 3) 
    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    
    # plot window:
    x <- NULL
    y <- NULL
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes,
                                               y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour = "black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fdiv),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FDiv values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FDiv", values_fdiv))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                         fill = color_sp[["asb1"]], 
                         color = color_sp[["asb1"]],
                         shape = shape_sp[["asb1"]],
                         size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.58,
                          size = 5, shape = 21,
                          color = color_sp[["asb1"]], fill = NA) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.50,
                         y = range_faxes[2] - spread_faxes*0.58,
        label = paste0("mean distance to the gravity center of species from", 
                       sep = " ", asb1), colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                           y = range_faxes[2] - spread_faxes*0.665,
                           fill = color_sp[["asb2"]], 
                           color = color_sp[["asb2"]],
                           shape = shape_sp[["asb2"]],
                           size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.735,
                            size = 5, 
                            shape = 21,
                            color = color_sp[["asb2"]], 
                            fill = fill_sp[["asb2"]]) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.50,
                           y = range_faxes[2] - spread_faxes*0.735,
          label = paste0("mean distance to the gravity center of species from", 
                         sep = " ", asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.40,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)

    # arrange panels ####
    
    # marge panels and caption:
    patchwork_fdiv <- panels.to.patchwork(panels_fdiv, plot_caption)
    
    # retrieve title and caption:
    tit_fdiv <- paste0( "Functional Divergence of '", asb1, "'")
    if (two_asb) {
      tit_fdiv <- paste0(tit_fdiv, " and '", asb2,"'")
    }

    patchwork_fdiv <- patchwork_fdiv +
    patchwork::plot_annotation(title = tit_fdiv,
                               caption = "made with mFD package")
    
    
    # save as file or in list:
    if (save_file == TRUE) {
      
      # name of file built with assemblage names and number of dimensions:
      file_fdiv <- paste0(nm_asb, "_", "FDiv_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fdiv ,
                      plot = patchwork_fdiv,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units= "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption:
      names(panels_fdiv) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fdiv[["caption"]] <- plot_caption
      panels_fdiv[["patchwork"]] <- patchwork_fdiv
      list_panels[["fdiv"]] <- panels_fdiv
      
    }

  }# end of plotting FDiv
  
  
  ##############################################################################
  
  
  # FSpe plot if required ####
  if ("fspe" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fspe <- list()

    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # retrieve coordinates of center of space:
      pool_0_xy <- fd_details$pool_O_coord[xy_k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # get a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2,]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <-fspe.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          center_coord2D = pool_0_xy,
                          pool_coord2D = sp_coord_xy,
                          plot_pool = FALSE,
                          plot_sp = TRUE,
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = color_sp[c("asb1", "asb2")],
                          color_center = color_centroid_fspe,
                          fill_center = color_centroid_fspe,
                          shape_center = shape_centroid_fspe,
                          size_center = 3,
                          color_segment = color_sp[c("asb1", "asb2")],
                          width_segment = c(asb1 = 1, asb2 = 1),
                          linetype_segment = c(asb1 = 1, asb2 = 1))
      
      
      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                      y = xy_k[2],
                                                      label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fspe[[k]] <- plot_k
      
    } # end of k
    
    # plot caption:
    
    # retrieve values to plot:
    top_fspe <- c("Functional Specialisation", asb1, "")
    values_fspe <- c(round(asb_fd_ind[asb1, "fspe"], 3), "")
    if (two_asb) {
      top_fspe[3] <- asb2
      values_fspe[2] <- round(asb_fd_ind[asb2,"fspe"], 3) 
    }
    
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
    
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fspe),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FSpe values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FSpe", values_fspe))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = color_sp[["asb1"]], 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3)
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.58,
                            fill = color_sp[["asb2"]], 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.58,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for functional gravity center:
    
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.65,
                          size = 3, 
                          shape = shape_centroid_fspe,
                          color = color_centroid_fspe, 
                          fill = color_centroid_fspe) +
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.36,
                         y = range_faxes[2] - spread_faxes*0.65,
                         label = "gravity center of functional space",
                         colour = color_centroid_fspe, size = 3)
      
    
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.72,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.36,
                         y = range_faxes[2] - spread_faxes*0.72,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fspe <- panels.to.patchwork(panels_fspe, plot_caption)
    
    # add title and caption:
    tit_fspe <- paste0( "Functional Specialization of '", asb1, "'")
    if (two_asb){
      tit_fspe <- paste0(tit_fspe, " and '", asb2, "'")
    }
    
    
    patchwork_fspe <- patchwork_fspe +
      patchwork::plot_annotation(title = tit_fspe,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      
      file_fspe <- paste0(nm_asb, "_", "fspe_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fspe ,
                      plot = patchwork_fspe,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fspe) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fspe[["caption"]] <- plot_caption
      panels_fspe[["patchwork"]] <- patchwork_fspe
      list_panels[["fspe"]] <- panels_fspe
      
    }
    

  } # end of plotting FSpe
  
  
  ##############################################################################
  

  # FDis plot if required ####
  if ("fdis" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fdis <- list()
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      asb_fide_coord2D <- list()
      asb_fide_coord2D[["asb1"]] <- asb_fd_ind[asb1, paste0("fide", sep = "_", 
                                                              xy_k)]
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        asb_fide_coord2D[["asb2"]] <- asb_fd_ind[asb2, paste0("fide", sep = "_", 
                                                                xy_k)]
      }# end of if 2 assemblages
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <- fdis.plot(ggplot_bg = plot_k,
                         asb_sp_coord2D = asb_sp_coord2D_k,
                         asb_sp_relatw = asb_sp_relw_k,
                         asb_fide_coord2D = asb_fide_coord2D,
                         plot_sp = TRUE,
                         shape_sp = shape_sp[c("asb1", "asb2")],
                         color_sp = color_sp[c("asb1", "asb2")],
                         fill_sp = fill_sp[c("asb1", "asb2")],
                         color_fide = color_sp[c("asb1", "asb2")],
                         fill_fide = fill_sp[c("asb1", "asb2")],
                         shape_fide = shape_centroid_fdis,
                         size_fide = c(asb1 = 2, asb2 = 2),
                         color_segment = color_sp[c("asb1", "asb2")],
                         width_segment = c(asb1 = 1, asb2 = 1),
                         linetype_segment = c(asb1 = 1, asb2 = 1))
      
      
      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                       y = xy_k[2],
                                                       label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fdis[[k]] <- plot_k
      
    } # end of k
    
    # plot caption:
    
    # retrieve values to plot:
    top_fdis <- c("Functional Dispersion", asb1, "")
    values_fdis <- c(round(asb_fd_ind[asb1, "fdis"], 3), "")
    if (two_asb) {
      top_fdis[3] <- asb2
      values_fdis[2] <- round(asb_fd_ind[asb2,"fdis"], 3) 
    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    x <- NULL
    y <- NULL
    
    # plot window:
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fdis),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FDis values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FDis", values_fdis))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = color_sp[["asb1"]], 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                        y = range_faxes[2] - spread_faxes*0.58,
                        size = 2, shape = shape_centroid_fdis[["asb1"]],
                        color = color_sp[["asb1"]], fill = fill_sp[["asb1"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("gravity center of", 
                                        sep = " ", asb1), 
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.665,
                            fill = color_sp[["asb2"]], 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.735,
                            size = 2, 
                            shape = shape_centroid_fdis[["asb2"]],
                            color = color_sp[["asb2"]], 
                            fill = fill_sp[["asb2"]]) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.735,
                           label = paste0("gravity center of", 
                                          sep = " ", asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.38,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fdis <- panels.to.patchwork(panels_fdis, plot_caption)
    
    # add title and caption:
    tit_fdis <- paste0( "Functional Dispersion of '", asb1, "'")
    if (two_asb) {
      tit_fdis <- paste0(tit_fdis, " and '", asb2, "'")
    }
    
    
    patchwork_fdis <- patchwork_fdis +
      patchwork::plot_annotation(title = tit_fdis,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      file_fdis <- paste0(nm_asb, "_", "fdis_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fdis,
                      plot = patchwork_fdis,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fdis) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fdis[["caption"]] <- plot_caption
      panels_fdis[["patchwork"]] <- patchwork_fdis
      list_panels[["fdis"]] <- panels_fdis
      
    }

  } # end of plotting FDis
  
  
  ##############################################################################
  
  # FIde plot if required ####
  if ("fide" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fide <- list()
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      asb_fide_coord2D <- list()
      asb_fide_coord2D[["asb1"]] <- asb_fd_ind[asb1, paste0("fide", sep = "_", 
                                                              xy_k)]
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        asb_fide_coord2D[["asb2"]] <- asb_fd_ind[asb2, paste0("fide", sep = "_", 
                                                              xy_k)]
      }# end of if 2 assemblages
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <- fide.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          asb_fide_coord2D = asb_fide_coord2D,
                          plot_sp = TRUE,
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = fill_sp[c("asb1", "asb2")],
                          color_fide = color_sp[c("asb1", "asb2")],
                          fill_fide = fill_sp[c("asb1", "asb2")],
                          shape_fide = shape_centroid_fdis,
                          size_fide = size_sp[c("asb1", "asb2")],
                          color_segment = color_sp[c("asb1", "asb2")],
                          width_segment = c(asb1 = 1, asb2 = 1),
                          linetype_segment = c(asb1 = 1, asb2 = 1))
      
      
      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                       y = xy_k[2],
                                                       label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fide[[k]] <- plot_k 
      
    } # end of k
    
    
    # plot caption:
    
    # retrieve values to plot:
    top_fide <- c("Functional Identity", asb1, "")
    if (ncol(sp_faxes_coord) == 2) {
      values_fide_PC1 <- c(round(asb_fd_ind[asb1,"fide_PC1"], 3), "")
      values_fide_PC2 <- c(round(asb_fd_ind[asb1,"fide_PC2"], 3), "")
      values_fide_PC3 <- NULL
      values_fide_PC4 <- NULL
    }
    if (ncol(sp_faxes_coord) == 3) {
      values_fide_PC1 <- c(round(asb_fd_ind[asb1,"fide_PC1"], 3), "")
      values_fide_PC2 <- c(round(asb_fd_ind[asb1,"fide_PC2"], 3), "")
      values_fide_PC3 <- c(round(asb_fd_ind[asb1,"fide_PC3"], 3), "")
      values_fide_PC4 <- NULL
    }
    if (ncol(sp_faxes_coord) == 4) {
      values_fide_PC1 <- c(round(asb_fd_ind[asb1,"fide_PC1"], 3), "")
      values_fide_PC2 <- c(round(asb_fd_ind[asb1,"fide_PC2"], 3), "")
      values_fide_PC3 <- c(round(asb_fd_ind[asb1,"fide_PC3"], 3), "")
      values_fide_PC4 <- c(round(asb_fd_ind[asb1,"fide_PC4"], 3), "")
    }
    
    if (two_asb) {
      top_fide[3] <- asb2
      if (ncol(sp_faxes_coord) == 2) {
        values_fide_PC1[2] <- round(asb_fd_ind[asb2,"fide_PC1"], 3)
        values_fide_PC2[2] <- round(asb_fd_ind[asb2,"fide_PC2"], 3)
        values_fide_PC3[2] <- NULL
        values_fide_PC4[2] <- NULL
      }
      if (ncol(sp_faxes_coord) == 3) {
        values_fide_PC1[2] <- round(asb_fd_ind[asb2,"fide_PC1"], 3)
        values_fide_PC2[2] <- round(asb_fd_ind[asb2,"fide_PC2"], 3)
        values_fide_PC3[2] <- round(asb_fd_ind[asb2,"fide_PC3"], 3)
        values_fide_PC4[2] <- NULL
      }
      if (ncol(sp_faxes_coord) == 4) {
        values_fide_PC1[2] <- round(asb_fd_ind[asb2,"fide_PC1"], 3)
        values_fide_PC2[2] <- round(asb_fd_ind[asb2,"fide_PC2"], 3)
        values_fide_PC3[2] <- round(asb_fd_ind[asb2,"fide_PC3"], 3)
        values_fide_PC4[2] <- round(asb_fd_ind[asb2,"fide_PC4"], 3)
      }

    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- c(0.3, 0.6, 0.9, 1.2)
    x <- NULL
    y <- NULL
    
    # plot window:
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fide),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FDis values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes *  0.15 * vv[c(1:4)],
      values_lab = c(values_fide_PC1, values_fide_PC2,
                     values_fide_PC3, values_fide_PC4))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = color_sp[["asb1"]], 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.58,
                          size = 2, shape = shape_centroid_fdis[["asb1"]],
                          color = color_sp[["asb1"]], fill = fill_sp[["asb1"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("gravity center of", 
                                        sep = " ", asb1), 
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.665,
                            fill = color_sp[["asb2"]], 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.735,
                            size = 2, 
                            shape = shape_centroid_fdis[["asb2"]],
                            color = color_sp[["asb2"]], 
                            fill = fill_sp[["asb2"]]) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.735,
                           label = paste0("gravity center of", 
                                          sep = " ", asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.38,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fide <- panels.to.patchwork(panels_fide, plot_caption)
    
    # add title and caption:
    tit_fide <- paste0( "Functional Identity of '", asb1, "'")
    if (two_asb) {
      tit_fide <- paste0(tit_fide, " and '", asb2, "'")
    }
    
    
    patchwork_fide <- patchwork_fide +
      patchwork::plot_annotation(title = tit_fide,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      file_fide <- paste0(nm_asb, "_", "fide_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fide,
                      plot = patchwork_fide,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fide) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fide[["caption"]] <- plot_caption
      panels_fide[["patchwork"]] <- patchwork_fide
      list_panels[["fide"]] <- panels_fide
      
    }
    
  } # end of plotting FIde
  
  
  ##############################################################################
  
  
  # FEve plot if required ####
  if ("feve" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_feve <- list()
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      asb_mst <- list()
      asb_mst[["asb1"]] <- fd_details$asb_mst[[asb1]]
      
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        asb_mst[["asb2"]] <- fd_details$asb_mst[[asb2]]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <- feve.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          asb_mst = asb_mst,
                          plot_sp = TRUE,
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = fill_sp[c("asb1", "asb2")],
                          color_mst = color_sp[c("asb1", "asb2")],
                          width_mst = c(asb1 = 1, asb2 = 1),
                          linetype_mst = c(asb1 = 1, asb2 = 1))
      
      
      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                       y = xy_k[2],
                                                       label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_feve[[k]] <- plot_k 
      
    } # end of k
    
    # retrieve values to plot:
    top_feve <- c("Functional Evenness", asb1, "")
    values_feve <- c(round(asb_fd_ind[asb1, "feve"], 3), "")
    if (two_asb) {
      top_feve[3] <- asb2
      values_feve[2] <- round(asb_fd_ind[asb2,"feve"], 3) 
    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    x <- NULL
    y <- NULL
    
    # plot window:
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_feve),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FEve values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FEve", values_feve))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = NA, 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.118,
                            xend = range_faxes[1] + spread_faxes*0.132,
                            y = range_faxes[2] - spread_faxes*0.58,
                            yend = range_faxes[2] - spread_faxes*0.58,
                            size = 1,
                            color = color_sp[["asb1"]],
                            linetype = 1) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("mst of", 
                                        sep = " ", asb1), 
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.665,
                            fill = NA, 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.118,
                            xend = range_faxes[1] + spread_faxes*0.132,
                            y = range_faxes[2] - spread_faxes*0.735,
                            yend = range_faxes[2] - spread_faxes*0.735,
                            size = 1, 
                            color = color_sp[["asb2"]], 
                            linetype = 1) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.735,
                           label = paste0("mst of", 
                                          sep = " ", asb2),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.38,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_feve <- panels.to.patchwork(panels_feve, plot_caption)
    
    # add title and caption:
    tit_feve <- paste0( "Functional Evenness of '", asb1, "'")
    if (two_asb) {
      tit_feve <- paste0(tit_feve, " and '", asb2, "'")
    }
    
    
    patchwork_feve <- patchwork_feve +
      patchwork::plot_annotation(title = tit_feve,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      file_feve <- paste0(nm_asb, "_", "feve_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_feve,
                      plot = patchwork_feve,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_feve) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_feve[["caption"]] <- plot_caption
      panels_feve[["patchwork"]] <- patchwork_feve
      list_panels[["feve"]] <- panels_feve
      
    }
    
  } # end of plotting FEve
  
  
  ##############################################################################
  
  
  # FOri plot if required ####
  if ("fori" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fori <- list()
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      asb_nn_pool <- list()
      asb_nn_pool[["asb1"]] <- fd_details$asb_nm_nn_pool[[asb1]]
      pool_coord2D <- sp_faxes_coord[, xy_k]
      
      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        asb_nn_pool[["asb2"]] <- fd_details$asb_nm_nn_pool[[asb2]]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <- fori.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          asb_nn_pool = asb_nn_pool,
                          pool_coord2D = pool_coord2D,
                          plot_pool = TRUE,
                          plot_sp = TRUE,
                          shape_pool = shape_sp[c("pool")],
                          size_pool = size_sp[c("pool")],
                          color_pool = color_sp[c("pool")],
                          fill_pool = fill_sp[c("pool")],
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = fill_sp[c("asb1", "asb2")],
                          color_segment = color_sp[c("asb1", "asb2")],
                          width_segment = c(asb1 = 1, asb2 = 1),
                          linetype_segment = c(asb1 = 1, asb2 = 1))
                          

      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                       y = xy_k[2],
                                                       label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fori[[k]] <- plot_k 
      
    } # end of k
    
    # retrieve values to plot:
    top_fori <- c("Functional Originality", asb1, "")
    values_fori <- c(round(asb_fd_ind[asb1, "fori"], 3), "")
    if (two_asb) {
      top_fori[3] <- asb2
      values_fori[2] <- round(asb_fd_ind[asb2,"fori"], 3) 
    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    x <- NULL
    y <- NULL
    
    # plot window:
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fori),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FOri values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FOri", values_fori))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = NA, 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.115,
                            xend = range_faxes[1] + spread_faxes*0.135,
                            y = range_faxes[2] - spread_faxes*0.58,
                            yend = range_faxes[2] - spread_faxes*0.58,
                            arrow = grid::arrow(length = grid::unit(0.07,
                                                                    "inches"),
                                                ends = "last",
                                                type = "open"),
                            size = 1,
                            color = color_sp[["asb1"]],
                            linetype = 1) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("nearest neighbour for each species of", 
                                        sep = " ", asb1, " in the global pool"), 
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.665,
                            fill = NA, 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.115,
                              xend = range_faxes[1] + spread_faxes*0.135,
                              y = range_faxes[2] - spread_faxes*0.735,
                              yend = range_faxes[2] - spread_faxes*0.735,
                              arrow = grid::arrow(length = grid::unit(0.07,
                                                                      "inches"),
                                                  ends = "last",
                                                  type = "open"),
                              size = 1, 
                              color = color_sp[["asb2"]], 
                              linetype = 1) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.735,
                        label = paste0("nearest neighbour for each species of", 
                                        sep = " ", asb2, " in the global pool"),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.38,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fori <- panels.to.patchwork(panels_fori, plot_caption)
    
    # add title and caption:
    tit_fori <- paste0( "Functional Originality of '", asb1, "'")
    if (two_asb) {
      tit_fori <- paste0(tit_fori, " and '", asb2, "'")
    }
    
    
    patchwork_fori <- patchwork_fori +
      patchwork::plot_annotation(title = tit_fori,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      file_fori <- paste0(nm_asb, "_", "fori_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fori,
                      plot = patchwork_fori,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fori) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fori[["caption"]] <- plot_caption
      panels_fori[["patchwork"]] <- patchwork_fori
      list_panels[["fori"]] <- panels_fori
      
    }
    
  } # end of plotting FOri
  
  ##############################################################################
  
  
  
  # FNND plot if required ####
  if ("fnnd" %in% ind_nm) {
    
    # create a list to store ggplot:
    panels_fnnd <- list()
    
    # loop on combinations:
    for (k in (1:plot_nb)) {
      
      # get names of axes:
      xy_k <- axes_plot[1:2, k]
      
      # get species coordinates along the 2 axes:
      sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
      colnames(sp_coord_xy) <- c("x", "y")
      
      # create a list with dataframes for plot:
      asb_sp_coord2D_k <- list()
      asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
      asb_sp_relw_k <- list()
      asb_sp_relw_k[["asb1"]] <- fd_details$asb_sp_relatw[asb1, sp_asb1]
      asb_nn_asb <- list()
      asb_nn_asb[["asb1"]] <- fd_details$asb_nm_nn_asb[[asb1]]

      if (two_asb){
        asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
        asb_sp_relw_k[["asb2"]] <- fd_details$asb_sp_relatw[asb2, sp_asb2]
        asb_nn_asb[["asb2"]] <- fd_details$asb_nm_nn_asb[[asb2]]
      }
      
      
      # background = axes defined by range of values and names as specified:
      plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
      
      # add species pool
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
                          color_vert = color_vert["pool"],
                          fill_vert = fill_vert["pool"])
      
      plot_k <- fnnd.plot(ggplot_bg = plot_k,
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_sp_relatw = asb_sp_relw_k,
                          asb_nn_asb = asb_nn_asb,
                          plot_sp = TRUE,
                          shape_sp = shape_sp[c("asb1", "asb2")],
                          color_sp = color_sp[c("asb1", "asb2")],
                          fill_sp = fill_sp[c("asb1", "asb2")],
                          color_segment = color_sp[c("asb1", "asb2")],
                          width_segment = c(asb1 = 1, asb2 = 1),
                          linetype_segment = c(asb1 = 1, asb2 = 1))
      
      
      # add species names if needed:
      if (!is.null(plot_sp_nm)) {
        x <- NULL
        y <- NULL
        plot_k <- plot_k +
          ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                                   ggplot2::aes_string(x = xy_k[1],
                                                       y = xy_k[2],
                                                       label = "label"),
                                   size = size_sp_nm, colour = color_sp_nm,
                                   fontface = fontface_sp_nm,
                                   max.overlaps = Inf,
                                   box.padding = grid::unit(2, 'lines'),
                                   force = 5,
                                   arrow = grid::arrow(length = grid::unit(0.02,
                                                                        'npc')),
                                   segment.color = color_sp_nm)
      }
      
      # save plot in a list:
      panels_fnnd[[k]] <- plot_k 
      
    } # end of k
    
    
    # retrieve values to plot:
    top_fnnd <- c("Functional Nearest Neighbour Distance", asb1, "")
    values_fnnd <- c(round(asb_fd_ind[asb1, "fnnd"], 3), "")
    if (two_asb) {
      top_fnnd[3] <- asb2
      values_fnnd[2] <- round(asb_fd_ind[asb2,"fnnd"], 3) 
    }
    
    # customize position of texts in the plot:
    spread_faxes <- (range_faxes[2] - range_faxes[1])
    hh <- c(1, 2.5, 4, 5.5)
    vv <- 0.3
    x <- NULL
    y <- NULL
    
    # plot window:
    plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                    ggplot2::aes(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
                         ymin = range_faxes[1], ymax = range_faxes[2],
                         fill = "white", colour ="black")
    
    # plot names of index and of assemblages:
    h   <- NULL
    v   <- NULL
    top <- NULL
    x <- NULL
    y <- NULL
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:4)],
        v = range_faxes[2] - spread_faxes * rep(0.2, 3),
        top = top_fnnd),
        ggplot2::aes(x = h, y = v, label = top),
        size = 3, hjust = 0.5, fontface = "bold")
    
    # plot FNND values:
    values_lab <- NULL
    data_caption <- data.frame(
      h = range_faxes[1] + spread_faxes * 0.15 * hh[2:4],
      v = range_faxes[2] - spread_faxes*rep(vv, 3),
      values_lab = c("FNND", values_fnnd))
    plot_caption <- plot_caption +
      ggplot2::geom_text(data = data_caption,
                         ggplot2::aes(x = h, y = v, label = values_lab),
                         size = 3, hjust = 0.5, fontface = "plain")
    vv <- vv + 0.1
    
    # add text about dimensionality:
    nb <- NULL
    plot_caption<- plot_caption +
      ggplot2::geom_text(data = data.frame(
        h = range_faxes[1] + spread_faxes * 0.1,
        v = range_faxes[2] - spread_faxes * vv,
        nb = paste0("NB: Indices were computed in a ",
                    nb_dim,"-dimensional space")),
        ggplot2::aes(x = h, y = v, label = nb),
        size = 3, hjust = 0, fontface = "italic")
    
    # add legend (convex hull, asb species and pool species):
    
    ## plot legend:
    values_lab <- NULL
    
    ### for 1st asb:
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.51,
                          fill = NA, 
                          color = color_sp[["asb1"]],
                          shape = shape_sp[["asb1"]],
                          size = 3) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.51,
                         label = paste0("relative weight of", sep = " ", 
                                        asb1),
                         colour = color_sp[["asb1"]], size = 3) + 
      
      ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.115,
                            xend = range_faxes[1] + spread_faxes*0.135,
                            y = range_faxes[2] - spread_faxes*0.58,
                            yend = range_faxes[2] - spread_faxes*0.58,
                            arrow = grid::arrow(length = grid::unit(0.07,
                                                                    "inches"),
                                                ends = "last",
                                                type = "open"),
                            size = 1,
                            color = color_sp[["asb1"]],
                            linetype = 1) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                         y = range_faxes[2] - spread_faxes*0.58,
                         label = paste0("nearest neighbour for each species of", 
                                        sep = " ", asb1, " in the assemblage"), 
                         colour = color_sp[["asb1"]], size = 3) 
    
    ### if 2nd assemblage:
    if (two_asb) {
      
      plot_caption <- plot_caption +
        ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                            y = range_faxes[2] - spread_faxes*0.665,
                            fill = NA, 
                            color = color_sp[["asb2"]],
                            shape = shape_sp[["asb2"]],
                            size = 3) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.665,
                           label = paste0("relative weight of", sep = " ", 
                                          asb2),
                           colour = color_sp[["asb2"]], size = 3) + 
        
        ggplot2::geom_segment(x = range_faxes[1] + spread_faxes*0.115,
                              xend = range_faxes[1] + spread_faxes*0.135,
                              y = range_faxes[2] - spread_faxes*0.735,
                              yend = range_faxes[2] - spread_faxes*0.735,
                              arrow = grid::arrow(length = grid::unit(0.07,
                                                                      "inches"),
                                                  ends = "last",
                                                  type = "open"),
                              size = 1, 
                              color = color_sp[["asb2"]], 
                              linetype = 1) + 
        
        ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
                           y = range_faxes[2] - spread_faxes*0.735,
                           label = paste0("nearest neighbour for each species of", 
                                          sep = " ", asb2, " in the assemblage"),
                           colour = color_sp[["asb2"]], size = 3)
      
    }
    
    ### for global pool:
    
    plot_caption <- plot_caption +
      
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
                          y = range_faxes[2] - spread_faxes*0.80,
                          size = size_sp[["pool"]], 
                          shape = shape_sp[["pool"]],
                          color = color_sp[["pool"]], 
                          fill = fill_sp[["pool"]]) + 
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.38,
                         y = range_faxes[2] - spread_faxes*0.80,
                         label = "shape of species from the global pool",
                         colour = color_sp[["pool"]], size = 3)
    
    # arrange panels ####
    
    # merge panels and caption:
    patchwork_fnnd <- panels.to.patchwork(panels_fnnd, plot_caption)
    
    # add title and caption:
    tit_fnnd <- paste0( "Functional Nearest Neighbor Distance of '", asb1, "'")
    if (two_asb) {
      tit_fnnd <- paste0(tit_fnnd, " and '", asb2, "'")
    }
    
    
    patchwork_fnnd <- patchwork_fnnd +
      patchwork::plot_annotation(title = tit_fnnd,
                                 caption = "made with mFD package")
    
    
    # save as file or in a list:
    if (save_file == TRUE) {
      
      # get name of file built with assemblage names and number of dimensions:
      file_fnnd <- paste0(nm_asb, "_", "fnnd_", nb_dim, "D" , ".", device_file)
      
      ggplot2::ggsave(filename = file_fnnd,
                      plot = patchwork_fnnd,
                      device = device_file,
                      scale = 1,
                      height= height_file[as.character(plot_nb)],
                      width = width_file[as.character(plot_nb)],
                      units = "in",
                      dpi = res_file)
    } else {
      
      # output = patchwork + list of panels and caption
      names(panels_fnnd) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
      panels_fnnd[["caption"]] <- plot_caption
      panels_fnnd[["patchwork"]] <- patchwork_fnnd
      list_panels[["fnnd"]] <- panels_fnnd
      
    }
    
  } # end of plotting FNND
  

  ##############################################################################
  

  # return complete output  = list of list of panels: ####
  
  return(list_panels)
  
}   # end of function
