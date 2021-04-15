#' Plot Species Position in a Functional Space 
#'
#' This function illustrates the position of species along pairs of axes of a 
#' functional space
#'
#' @param sp_faxes_coord a matrix of species coordinates in a
#'  multidimensional functional space. Species coordinates have been retrieved
#'  thanks to \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param faxes a vector with names of axes to plot (as columns names in
#'  \code{sp_faxes_coord}). \strong{You can only plot from 2 to 4 axes for
#'  graphical reasons}. Default: `faxes = NULL` (the four first axes will be
#'  plotted).
#'
#' @param name_file a character string with name of file to save the
#'   figure (without extension). Default: `name_file = NULL` which means plot is
#'   displayed.
#'
#' @param faxes_nm a vector with axes labels for figure. Default: as
#'  \code{faxes}).
#'
#' @param range_faxes a vector with minimum and maximum values of axes 
#' used for all plots to have a fair representation of position of species. 
#' Default: `range_faxes = c(NA, NA)` (the range is computed
#'  according to the range of values among all axes). If at least one of the 
#'  value provided is within the range of coordinates, then convex hull could 
#'  not be plotted so `plot_ch` should be `FALSE`.
#'
#' @param color_bg a R color name or an hexadecimal code used to fill plot
#'   background. Default: `color_bg = "grey95"`.
#'
#' @param color_pool a R color name or an hexadecimal code referring to the color
#'   of symbol for species. Default: `color_pool = 'darkgreen'`.
#'
#' @param fill_pool a R color name or an hexadecimal code referring to the color
#'   to fill species symbol (if \code{shape_pool} >20). Default:
#'   `fill_pool = 'white'`.
#'
#' @param shape_pool a numeric value referring to the shape of symbol used for
#'   species. Default: `shape_pool = 21` (filled circle).
#'
#' @param size_pool a numeric value referring to the size of symbol for species.
#'   Default: `size_pool = 1`.
#'
#' @param plot_ch a logical value indicating whether the convex hull shaping the
#'   pool of species should be illustrated. If `plot_ch = TRUE`, convex hull of
#'   all species in the multidimensional space described in
#'   \code{sp_faxes_coord} is computed and its projection in 2D spaces are drawn
#'   as polygons. Default: `plot_ch = TRUE`.
#'
#' @param color_ch a R color name or an hexadecimal code referring to the border
#'   of the convex hull filled by the pool of species. Default:
#'   `color_ch = "darkblue"`.
#'
#' @param fill_ch a R color name or an hexadecimal code referring to the filling
#'   of the convex hull filled by the pool of species. Default:
#'   `fill_ch = "white"`.
#'
#' @param alpha_ch a numeric value for transparency of the filling of the convex
#'   hull (0 = high transparency, 1 = no transparency). Default:
#'   `alpha_ch = 1`.
#'
#' @param plot_vertices a logical value defining whether vertices of the convex
#'   hull shaping the pool of species should be illustrated. If
#'   `plot_vertices = TRUE`, vertices of convex hull computed in the
#'   multidimensional space from `sp_faxes_coord` and are plotted with
#'   aesthetics listed below '..._vert' (species not being vertices are plotted
#'   with aesthetics described above for '.._sp'. Default:
#'   `plot_vertices = TRUE`.
#'
#' @param color_vert a character value referring to the color of symbol for
#'   vertices if `plot_vertices = TRUE`. Default: `color_vert = 'darkturquoise'`.
#'
#' @param fill_vert a character value referring to the color for filling symbol
#'  for vertices  (if \code{shape_vert} >20). Default:
#'  `fill_vert = 'darkturquoise'`.
#'
#' @param shape_vert a numeric value referring to the symbol used to show
#'   vertices position if `plot_vertices = TRUE`. Default: `shape_vert = 23`
#'   (filled diamond).
#'
#' @param size_vert a numeric value referring to the size of symbol for
#'   vertices Default: `size_vert = 1`.
#'
#' @param plot_sp_nm a vector containing species names that are to be printed
#'   near their position. Default: `plot_nm_sp = NULL` (no name plotted).
#'
#' @param nm_size a numeric value for size of species label. Default is `3`
#'   (in points).
#'
#' @param nm_color a R color name or an hexadecimal code referring to the color
#'   of species label. Default: `nm_color = 'black'`.
#'
#' @param nm_fontface a character string for font of species labels (e.g.
#'   "italic", "bold"). Default: `nm_fontface = 'plain'`.
#'
#' @param check_input a logical value indicating whether key features the inputs
#'   are checked (e.g. class and/or mode of objects, names of rows and/or
#'   columns, missing values). If an error is detected, a detailed message is
#'   returned. Default: `check.input = TRUE`.
#'
#' @return If \code{name_file} is \code{NULL}, a list containing \code{ggplot2}
#'   objects is returned: plots of functional space along all pairs of axes
#'   (named according to axes names, e.g. "PC1_PC2"), figure 'caption', and the
#'   full figure 'patchwork' built using the library \code{patchwork}.
#'   If \code{name_file} is not \code{NULL} a 300dpi png file
#'   is saved in the working directory. Ranges of axes are the same for all
#'   panels and if required projection of the convex hull computed in the
#'   multidimensional space provided as input \code{sp_faxes_coord} is
#'   illustrated with a polygon. Species being vertices of this convex hull are
#'   shown with aesthetics provided as inputs \code{..._vert}. Labels for
#'   species listed in \code{plot_sp_nm} are added with if required arrows using
#'   \code{ggrepel}. Summary about species and dimensionality are printed on
#'   top-right corner of the figure.
#'
#' @author Camille Magneville and Sebastien Villeger
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#'                                   tr_cat         = fruits_traits_cat,
#'                                   metric         = "gower",
#'                                   scale_euclid   = "scale_center",
#'                                   ordinal_var    = "classic",
#'                                   weight_type    = "equal",
#'                                   stop_if_NA     = TRUE)
#'
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#'  fspaces_quality_fruits <- mFD::quality.fspaces(
#'   sp_dist             = sp_dist_fruits,
#'   maxdim_pcoa         = 10,
#'   deviation_weighting = "absolute",
#'   fdist_scaling       = FALSE,
#'   fdendro             = "average")
#'
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#'
#' # Plot functional spaces:
#'  mFD::funct.space.plot(
#'   sp_faxes_coord    = sp_faxes_coord_fruits[, c("PC1", "PC2", "PC3", "PC4")],
#'   faxes             = NULL,
#'   name_file         = NULL,
#'   faxes_nm          = NULL,
#'   range_faxes       = c(NA, NA),
#'   color_bg          = "grey95",
#'   color_pool          = "darkturquoise",
#'   fill_pool           = "white",
#'   shape_pool          = 21,
#'   size_pool           = 1,
#'   plot_ch           = TRUE,
#'   color_ch          = "darkblue",
#'   fill_ch           = "white",
#'   alpha_ch          = 1,
#'   plot_vertices     = TRUE,
#'   color_vert        = "darkturquoise",
#'   fill_vert         = "darkturquoise",
#'   shape_vert        = 22,
#'   size_vert         = 1,
#'  plot_sp_nm         = NULL,
#'  nm_size            = 3,
#'  nm_color           = "black",
#'  nm_fontface        = "plain",
#'  check_input        = TRUE)
#' }

funct.space.plot <- function(sp_faxes_coord, faxes = NULL, 
                             name_file = NULL,
                             faxes_nm = NULL, 
                             range_faxes = c(NA, NA),
                             color_bg = "grey95",
                             color_pool = "darkturquoise", fill_pool = "white",
                             shape_pool = 21, size_pool = 1,
                             plot_ch = TRUE, color_ch = "darkblue",
                             fill_ch = "white", alpha_ch = 1,
                             plot_vertices = TRUE, 
                             color_vert = "darkturquoise",
                             fill_vert = "darkturquoise", shape_vert = 22,
                             size_vert = 1,
                             plot_sp_nm = NULL, nm_size = 3, 
                             nm_color = "black",
                             nm_fontface = "plain",
                             check_input = TRUE) {
  
  # check_inputs if asked: ####
  if (check_input) {
    
    check.sp.faxes.coord(sp_faxes_coord)
    
    if ((!is.null(plot_sp_nm)) && any(!plot_sp_nm %in%
                                      rownames(sp_faxes_coord))) {
      stop("Species names in 'plot_sp_nm' can not be found in ",
           "'sp_faxes_coord' row names. Please check names of species you ",
           "want to plot.")
    }
    
    if (! is.null(faxes)) {
      
      if(length(faxes) == 1) {
        stop("Number of functional axes should be more than one. Please change",
             " the number of functional axes to plot.")
      }
      
      if (length(faxes) > 4) {
        stop("Number of functional axes should be less than 4. Please change ",
             "the number of functional axes to plot.")
      }
      
      if (any(! faxes %in% colnames(sp_faxes_coord))) {
        stop("Names of axes to plot can not be found in 'sp_faxes_coord' ",
             "columns names. Please check names of axes you want to plot.")
      }
    }
    
    if ((!is.null(faxes_nm)) && (length(faxes_nm) != length(faxes))) {
      stop("Length of 'faxes_nm' should be equal to length of 'faxes'. Please ",
           "check congruence between these inputs.")
    }
    
  }# end of checking inputs
  
  
  # retrieve names and number of axes to plot ####
  
  # if no axes selected take up to the 4 first axes in the coordinates table:
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(c(4, ncol(sp_faxes_coord)))]
  }
  
  # if no faxes_nm provided, default names as in coordinates table:
  if (!is.null(faxes) && is.null(faxes_nm)) {
    faxes_nm <- faxes
    names(faxes_nm) <- faxes
  } 
  
  # number of axes:
  nb_faxes <- length(faxes)
  
  # combinations of axes on plot:
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
  
  # vertices
  vert_nm <- NULL
  
  # if convex hull or vertices to be plotted
  if (plot_ch || plot_vertices) {
    
    # names of vertices of the convex hull of all species
    vert_nm <- vertices(sp_faxes_coord)
    
  }
  
  
  # check that the range is ok if the user chose it and convex hull:
  if (! is.na(user_range) && plot_ch) {
    
    if (range_faxes[1] > range_sp_coord[1]) {
      stop("Error: The first value of 'range_faxes', is higher than minimum 
        value of 'sp_faxes_coord' so the convex hull can not be plotted. Please 
        change the minimal value of 'range_faxes' or set 'plot_ch' to FALSE.")
    }
    
    if (range_faxes[2] < range_sp_coord[2]) {
      stop("Error: The second value of 'range_faxes', is lower than maximum 
        value of 'sp_faxes_coord' so the convex hull can not be plotted. Please 
        change the maximal value of 'range_faxes' or set 'plot_ch' to FALSE.")
    }
    
  } 
  
  # dataframe with species coordinates and column for label
  label  <- NULL
  sp_coord_plot <- data.frame(sp_faxes_coord,
                              label = rep("", nrow(sp_faxes_coord)))
  
  
  # if some species names to be plotted, adding a binary variable to ...
  # ... sp_faxes_coord
  if(!is.null(plot_sp_nm)) {
    sp_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
  }
  
  
  # plot for each pair of axes ####
  
  # list to store panels
  panels <- list()
  
  # loop on combinations
  for (k in (1:plot_nb)) {
    
    # species coordinates along the 2 axes
    sp_coord_xy <- as.matrix(sp_coord_plot[, axes_plot[1:2, k]])
    colnames(sp_coord_xy) <- c("x", "y")
    
    # names of axes
    xy_k <- axes_plot[1:2, k]
    
    # background = axes defined by range of values and names as specified  ----
    plot_k <- background.plot(range_faxes, faxes_nm = xy_k, color_bg)
    
    # if required adding convex hull and/or vertices projected in 2D ----
    if (! plot_ch) {
      color_ch <- color_bg
      fill_ch <- color_bg
      alpha_ch <- 1
    }
    
    if  (! plot_vertices) {
      shape_vert <- shape_pool
      size_vert <- size_pool
      color_vert <- color_pool
      fill_vert <- fill_pool
    }
    
    # if size set to 0, species are not drawn:
    if (size_pool == 0 & size_vert == 0) {
      plot_pool <- FALSE 
    } else {
      plot_pool <- TRUE 
    } 
    
    
    plot_k <- pool.plot(ggplot_bg = plot_k,
                        sp_coord2D = sp_coord_xy,
                        vertices_nD = vert_nm,
                        plot_pool = plot_pool,
                        shape_pool = shape_pool, size_pool = size_pool,
                        color_pool = color_pool, fill_pool = fill_pool,
                        color_ch = color_ch, fill_ch = fill_ch, 
                        alpha_ch = alpha_ch,
                        shape_vert = shape_vert, size_vert = size_vert,
                        color_vert = color_vert, fill_vert = fill_vert)
    
    
    # adding species names if needed ----
    if (!is.null(plot_sp_nm)) {
      x <- NULL
      y <- NULL
      plot_k <- plot_k +
        ggrepel::geom_text_repel(data = sp_coord_plot,
                                 ggplot2::aes_string( x = axes_plot[1, k],
                                                      y = axes_plot[2, k],
                                                      label = "label"),
                                 size = nm_size, colour = nm_color,
                                 fontface = nm_fontface,
                                 box.padding = grid::unit(2, 'lines'),
                                 force = 5,
                                 arrow = grid::arrow(length = grid::unit(0.02,
                                                                         'npc')),
                                 segment.color = nm_color)
    }
    
    # saving plot in a list
    panels[[k]] <- plot_k
    
  } # end of k
  
  ## plot for legend and summary values #####
  spread_faxes <- (range_faxes[2] - range_faxes[1])
  x <- NULL
  y <- NULL
  plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes, y = range_faxes),
                                  ggplot2::aes(x = x, y = y)) +
    ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  
  
  if (! plot_vertices) {
    
    plot_caption <- plot_caption +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes * 0.1,
                          y = range_faxes[2] - spread_faxes * 0.5,
                          colour = color_pool, fill = fill_pool, shape = shape_pool,
                          size = size_pool) +
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                         y = range_faxes[2] - spread_faxes * 0.5,
                         label = paste0( nrow(sp_faxes_coord), " species"),
                         colour = color_pool, hjust = 0)
    ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                       y = range_faxes[2] - spread_faxes * 0.65,
                       label = paste0("plotted along ", nb_faxes,
                                      "axes from the ", ncol(sp_faxes_coord),
                                      "-dimensional space"),  hjust = 0)
  }
  
  if (plot_vertices) {
    
    plot_caption <- plot_caption +
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                         y = range_faxes[2] - spread_faxes * 0.15,
                         label = paste0( nrow(sp_faxes_coord), " species"),
                         hjust = 0, fontface = "bold") +
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                         y = range_faxes[2] - spread_faxes * 0.5,
                         label = paste0(nrow(sp_faxes_coord) - length(vert_nm),
                                        " species not vertices"),
                         colour = color_pool, hjust = 0) +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes * 0.1,
                          y = range_faxes[2] - spread_faxes * 0.5,
                          colour = color_pool, fill = fill_pool,
                          shape = shape_pool, size = size_pool) +
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                         y = range_faxes[2] - spread_faxes * 0.35,
                         label = paste0(length(vert_nm), " species vertices"),
                         colour = color_vert, hjust = 0) +
      ggplot2::geom_point(x = range_faxes[1] + spread_faxes * 0.1,
                          y = range_faxes[2] - spread_faxes * 0.35,
                          colour = color_vert, fill = fill_vert,
                          shape = shape_vert, size = size_vert) +
      
      ggplot2::geom_text(x = range_faxes[1] + spread_faxes * 0.15,
                         y = range_faxes[2] - spread_faxes * 0.65,
                         label = paste0("in the ", ncol(sp_faxes_coord),
                                        "-dimensional space"), hjust = 0)
    
  }
  
  
  #### arranging panels ####
  
  # merging panels and caption
  patchwork_plots_all <- panels.to.patchwork(panels, plot_caption)
  
  
  #### returning output ####
  
  # type, resolution and dimensions of file if to be saved
  device_file <- "png"
  res_file <- 300
  height_file <- 4 * c(1, 2, 3)
  names(height_file) <- c("1", "3", "6")
  width_file  <- 4 * c(2, 2, 3)
  names(width_file)  <- c("1", "3", "6")
  
  # saving in a file or returning list of ggplots
  if (! is.null(name_file))  {
    
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file) ,
                    plot = patchwork_plots_all,
                    device = device_file,
                    scale = 1,
                    height = height_file[as.character(plot_nb)],
                    width = width_file[as.character(plot_nb)],
                    units = "in",
                    dpi = res_file)
  } else {
    
    # output = list of panels and caption + patchwork
    output <- NULL
    output <- panels
    names(output) <- paste(axes_plot[1, ], axes_plot[2, ], sep = "_")
    output[["caption"]] <- plot_caption
    output[["patchwork"]] <- patchwork_plots_all
    
    return(output)
  }
  
} 
