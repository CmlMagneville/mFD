# Function to plot species in functional space
#
# Authors: Camille Magneville and Sébastien Villéger
#
#

# ------------------------------------------------------------------------------

#'Plot functional space with different functional axes pairs
#'
#'This function computes plots of functional spaces showing the position of
#'species along pairs of axes.
#'
#'@param sp_faxes_coord a \strong{matrix} of species coordinates in a
#'  multidimensional functional space. Species coordinates have been retrieved
#'  thanks to \code{fspace.conttr} or \code{qual_funct_space}.
#'
#'@param faxes a \strong{vector} with names of axes to plot (as columns names in
#'  \code{sp_faxes_coord} ). \strong{You can only plot from 2 to 4 axes for
#'  graphical reasons}. Default: faxes = NULL (the four first axes will be
#'  plotted).
#'
#'@param name_file a \strong{character string} with name of file to save the figure
#'  (without extension). Default is 'NULL' which means plot is displayed.
#'
#'@param faxes_nm a \strong{vector} with axes labels for figure. Default: as
#'  \code{faxes}).
#'
#'@param range_faxes_lim a \strong{vector} with minimum and maximum values of axes. Note
#'  that to have a fair representation of position of species in all plots, have
#'  the same range. Default: faxes_lim = c(NA, NA) (the range is computed
#'  according to the range of values among all axes).
#'
#'@param color_bg a \strong{R color name or an hexadecimal code} used to fill plot
#'  background. Default: color_bg = "grey95".
#'
#'@param color_sp a \strong{R color name or an hexadecimal code} referring to the colour
#'  of symbol for species. Default: color_sp = 'darkgreen'
#'
#'@param fill_sp a \strong{R color name or an hexadecimal code} referring to the colour
#'  to fill species symbol (if \code{shape_sp} >20). Default: fill_sp = 'white'
#'
#'@param shape_sp a \strong{numeric value} referring to the shape of symbol used for
#'  species. Default: shape_sp = 21 (filled circle).
#'
#'@param size_sp a \strong{numeric value} referring to the size of symbol for species.
#'  Default: size_sp = 1.
#'
#'@param plot_ch a \strong{logical value} indicating whether the convex hull shaping the
#'  pool of species should be illustrated. If plot_ch = TRUE, convex hull of all
#'  species in the multidimensional space described in \code{sp_faxes_coord} is
#'  computed and its projection in 2D spaces are drawn as polygons. Default:
#'  plot_ch = TRUE.
#'
#'@param color_ch a \strong{R color name or an hexadecimal code} referring to the border
#'  of the convex hull filled by the pool of species. Default: color_ch = "black".
#'
#'@param fill_ch a \strong{R color name or an hexadecimal code} referring to the filling
#'  of the convex hull filled by the pool of species. Default: fill_ch = "white".
#'
#'@param alpha_ch a \strong{numeric value} for transparency of the filling of the convex
#'  hull (0 = high transparency, 1 = no transparency). Default: alpha_ch = 0.3.
#'
#'@param plot_vertices a \strong{logical value} defining whether vertices of the convex
#'  hull shaping the pool of species should be illustrated. If plot_vertices =
#'  TRUE, vertices of convex hull computed in the multidimensional space from
#'  \code{sp_faxes_coord} and are plotted with aesthetics listed below
#'  '..._vert' (species not being vertices are plotted with aesthetics described
#'  above for '.._sp'. Default: plot_vertices = TRUE.
#'
#'@param color_vert a \strong{character value} referring to the color of symbol for
#'  vertices if \code{plot_vertices} = TRUE. Default: color_vert = 'blueviolet'.
#'
#'@param fill_vert a \strong{character value} referring to the color for filling symbol
#'  for vertices  (if \code{shape_vert} >20). Default: fill_vert = 'blueviolet'.
#'
#'@param shape_vert a \strong{numeric value} referring to the symbol used to show
#'  vertices position if \code{plot_vertices} = TRUE. Default: shape_vert = 23 (filled
#'  diamond).
#'
#'@param size_vert a \strong{numeric value} referring to the size of symbol for vertices
#'  Default: size_vert = 1.
#'
#'@param plot_sp_nm a \strong{vector containing} species names that are to be printed
#'  near their position. Default: plot_nm_sp = NULL (no name plotted).
#'
#'@param nm_size a \strong{numeric value} for size of species label. Default is 3 points.
#'
#'@param nm_color a \strong{R color name or an hexadecimal code} referring to the colour
#'  of species label. Default: nm_color = 'black'.
#'
#'@param nm_fontface a \strong{character string} for font of species labels (e.g.
#'  "italic", "bold"). Default: nm_fontface = 'plain'.
#'
#'@param check.input a \strong{logical value} defining whether inputs are checked before
#'  computation of indices. Possible error messages will thus may be more
#'  understandable for the user than R error messages. Default: check.input =
#'  TRUE.
#'
#'@return a list containing \code{ggplot2} objects that were built before
#'  assembling them in the figure using the library \code{patchwork}. If
#'  \code{name_file} is not \code{NULL} a 300dpi png file is saved in the
#'  working directory. Ranges of axes are the same for all panels and if
#'  required projection of the convex hull computed in the multidimensional
#'  space provided as input \code{sp_faxes_coord} is illustrated with a polygon.
#'  Species being vertices of this convex hull are shown with aesthetics
#'  provided as inputs \code{..._vert}. Labels for species listed in
#'  \code{plot_sp_nm} are added with if required arrows using \code{ggrepel}.
#'  Summary about species and dimensionality are printed on top-right corner of
#'  the figure.
#'
#'@examples
#' # Load Species*Traits dataframe:
#' data("sp_tr_fruits", package = "mFD")
#' # Load Assemblages*Species dataframe:      
#' data("asb_sp_w_fruits", package = "mFD")   
#' Load Traits categories dataframe:
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
#' # Plot functional spaces:
#' mFD::funct.space.plot(sp_faxes_coord_fruits[, c("PC1", "PC2", "PC3", "PC4")], 
#'  faxes = NULL, name_file = NULL,
#'  faxes_nm = NULL, range_faxes_lim = c(NA, NA),
#'  color_bg = "grey95",
#'  color_sp = "darkgreen", fill_sp = "white",  shape_sp = 21, size_sp = 1,
#'  plot_ch = TRUE,  color_ch = "black", fill_ch = "white", alpha_ch = 0.5,
#'  plot_vertices = TRUE, color_vert = "blueviolet", fill_vert = "blueviolet",
#'   shape_vert = 23, size_vert = 1,
#'  plot_sp_nm = NULL, nm_size = 3, nm_color = "black", nm_fontface = "plain",
#'  check.input = TRUE)
#'
#'@export


funct.space.plot <- function(sp_faxes_coord, faxes = NULL, name_file = NULL,
                             faxes_nm = NULL, range_faxes_lim = c(NA, NA),
                             color_bg = "grey95",
                             color_sp = "darkgreen", fill_sp = "white",  shape_sp = 21, size_sp = 1,
                             plot_ch = TRUE,  color_ch = "black", fill_ch = "white", alpha_ch = 0.3,
                             plot_vertices = TRUE, color_vert = "blueviolet", fill_vert = "blueviolet", shape_vert = 23, size_vert = 1,
                             plot_sp_nm = NULL, nm_size = 3, nm_color = "black", nm_fontface = "plain",
                             check.input = TRUE) {
  
  ## check inputs if asked: ####
  if (check.input == TRUE) {
    
    if(! is.matrix(sp_faxes_coord)){
      stop("Error: 'sp_faxes_coord' must be a matrix. Please change its class.")
    }
    sp_faxes_coord <- as.data.frame(sp_faxes_coord)
    
    if (any(is.na(sp_faxes_coord))) {
      stop("Error: The species*coordinates dataframe contains NA. Please check.")
    }
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates dataframe.
             Please add species names as row names.")
    }
    if (is.null(colnames(sp_faxes_coord))) {
      stop("Error: No column names provided in species*coordinates dataframe.
             Please add dimensions labels as column names.")
    }
    
    if ( (! is.null(plot_sp_nm)) &
         any( ! plot_sp_nm %in% rownames(sp_faxes_coord))) {
      stop("Error: species names in 'plot_sp_nm' can not be found in
           'sp_faxes_coord'row names. Please check names of species you want to plot.")
    }
    
    if (! is.null(faxes)) {
      
      if (length(faxes) > 4) {
        stop("Error: Number of functional axes should be less than 4.
          Please change the number of functional axes to plot.")
      }
      
      if (! any(faxes %in% colnames(sp_faxes_coord))) {
        stop("Error: names of axes to plot can not be found in 'sp_faxes_coord'
               columns names. Please check names of axes you want to plot.")
      }
      
    }
    
    
    if ((! is.null(faxes_nm)) & (length(faxes_nm) != length(faxes))) {
      stop("Error: Length of 'faxes_nm' should be equal to length of 'faxes'.
        Please check congruence between these inputs")
    }
  }# end of checking inputs
  
  
  #### names and number of axes to plot ####
  
  # if no axes selected take up to the 4 first axes in the coordinates table
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(c(4, ncol(sp_faxes_coord)))]
  }
  
  # if no faxes_nm provided, default names as in coordinates table
  if (!is.null(faxes) & is.null(faxes_nm)) {
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
  if (is.na(range_faxes_lim[1]) & is.na(range_faxes_lim[2])) {
    range_sp_coord <- range(sp_faxes_coord)
    range_faxes_lim <- range_sp_coord +
      c(-1,1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
  }
  
  # dataframe with species coordinates and option (vertices + label) if required
  vertex <- NULL
  label <- NULL
  sp_coord_plot <- data.frame(sp_faxes_coord,
                              vertex = rep(FALSE, nrow(sp_faxes_coord)),
                              label = rep("", nrow(sp_faxes_coord)))
  
  
  # if convex hull to be plotted
  if (plot_ch == TRUE | plot_vertices == TRUE) {
    
    # names and coordinates of vertices of the convex hull of all species
    vert_nm <- vertices(sp_faxes_coord, check.input = TRUE)
    vert_coord <- sp_faxes_coord[vert_nm, ]
    
    # if vertices to be shown, adding a binary variable to sp_faxes_coord
    if (plot_vertices == TRUE) {
      sp_coord_plot[vert_nm, "vertex"] <- TRUE
    }
    
  }
  
  
  # if some species names to be plotted, adding a binary variable to sp_faxes_coord
  if(! is.null(plot_sp_nm)) {
    sp_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
  }
  
  
  
  #### plot for each pair of axes ####
  
  # list to store panels
  panels <- list()
  
  # loop on combinations
  for (k in (1:plot_nb)) {
    
    # names of axes as in input dataframe
    x <- axes_plot[1, k]
    y <- axes_plot[2, k]
    
    # background = axes defined by range of values and names as specified  ----
    plot_k <- ggplot2::ggplot(sp_coord_plot,
                              ggplot2::aes_string(x = x, y = y)) +
      ggplot2::scale_x_continuous(limits = range_faxes_lim, expand=c(0,0)) +
      ggplot2::xlab(faxes_nm[x]) +
      ggplot2::scale_y_continuous(limits = range_faxes_lim, expand=c(0,0)) +
      ggplot2::ylab(faxes_nm[y]) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = color_bg)) +
      ggplot2::coord_fixed()
    
    # if required adding convex hull of pool projected in 2D ----
    if (plot_ch == TRUE) {
      
      # vertices in 2D(= subset of the vertices in n dimensions)
      vert_nm_xy <- vertices(sp_coord_plot[, c(x,y)], check.input = TRUE)
      vert_coord_xy <- sp_coord_plot[vert_nm_xy, c(x, y)]
      
      # sorting vertices clockwise for plotting as a polygon later
      vert_coord_xy_sorted <- vert_coord_xy[order(-1 * atan2(
        vert_coord_xy[, x] - mean(range(vert_coord_xy[, x])),
        vert_coord_xy[, y] - mean(range(vert_coord_xy[, y])) )
      ), ]
      
      # plotting
      plot_k <- plot_k +
        ggplot2::geom_polygon(data = vert_coord_xy_sorted,
                              colour = color_ch, fill = fill_ch, alpha = alpha_ch)
      
    }
    
    # points for species with aesthetics according to vertex status if needed  ----
    plot_k <- plot_k +
      ggplot2::geom_point(ggplot2::aes(colour = vertex, fill = vertex,
                                       shape = vertex, size = vertex), size = size_sp,
                          show.legend= FALSE) +
      ggplot2::scale_color_manual(name = "vertex",
                                  values = c("TRUE" = color_vert , "FALSE" =  color_sp)) +
      ggplot2::scale_fill_manual(name = "vertex",
                                 values = c("TRUE" = fill_vert , "FALSE" =  fill_sp)) +
      ggplot2::scale_shape_manual(name = "vertex",
                                  values = c("TRUE" = shape_vert , "FALSE" =  shape_sp)) +
      ggplot2::scale_size_manual(name = "vertex",
                                 values = c("TRUE" = size_vert , "FALSE" =  size_sp))
    
    
    
    # adding species names if needed ----
    if(! is.null(plot_sp_nm)) {
      plot_k <- plot_k +
        ggrepel::geom_text_repel(x = sp_coord_plot[, x], y = sp_coord_plot[, y],
                                 label = sp_coord_plot[,"label"],
                                 size = nm_size, colour= nm_color, fontface = nm_fontface,
                                 box.padding   = grid::unit(2, 'lines'),
                                 force = 5,
                                 arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
                                 segment.color = nm_color)
    }
    
    # saving plot in a list
    panels[[k]] <- plot_k
    
  } # end of k
  
  ## plot for legend and summary values #####
  spread_faxes <- (range_faxes_lim[2]- range_faxes_lim[1])
  
  plot_caption <- ggplot2::ggplot(sp_faxes_coord) +
    ggplot2::scale_x_continuous(limits = range_faxes_lim, expand=c(0,0)) +
    ggplot2::scale_y_continuous(limits = range_faxes_lim, expand=c(0,0)) +
    ggplot2::theme_void() + ggplot2::theme(legend.position="none") +
    ggplot2::geom_rect( xmin = range_faxes_lim[1], xmax = range_faxes_lim[2],
                        ymin = range_faxes_lim[1], ymax = range_faxes_lim[2],
                        fill = "white", colour = "black")
  
  
  if (plot_vertices == FALSE) {
    plot_caption <- plot_caption +
      ggplot2::geom_point (x = range_faxes_lim[1] + spread_faxes*0.1,
                           y = range_faxes_lim[2] - spread_faxes*0.5,
                           colour = color_sp, fill = fill_sp, shape = shape_sp, size = size_sp) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.5,
                         label = paste0( nrow(sp_faxes_coord), " species"), colour = color_sp, hjust = 0)
    ggplot2::geom_text( x = range_faxes_lim[1] + spread_faxes*0.15,
                        y = range_faxes_lim[2] - spread_faxes*0.65,
                        label = paste0("plotted along ", nb_faxes,"axes from the ", ncol(sp_faxes_coord), "-dimensional space"),  hjust = 0)
  }
  
  if (plot_vertices == TRUE) {
    
    plot_caption <- plot_caption +
      ggplot2::geom_text( x = range_faxes_lim[1] + spread_faxes*0.15,
                          y = range_faxes_lim[2] - spread_faxes*0.15,
                          label = paste0( nrow(sp_faxes_coord), " species"), hjust = 0, fontface = "bold") +
      
      ggplot2::geom_text( x = range_faxes_lim[1] + spread_faxes*0.15,
                          y = range_faxes_lim[2] - spread_faxes*0.5,
                          label = paste0(nrow(sp_faxes_coord) - length(vert_nm), " species not vertices"), colour = color_sp, hjust = 0) +
      ggplot2::geom_point (x = range_faxes_lim[1] + spread_faxes*0.1,
                           y = range_faxes_lim[2] - spread_faxes*0.5,
                           colour = color_sp, fill = fill_sp, shape = shape_sp, size = size_sp) +
      
      ggplot2::geom_text( x = range_faxes_lim[1] + spread_faxes*0.15,
                          y = range_faxes_lim[2] - spread_faxes*0.35,
                          label = paste0(length(vert_nm), " species vertices"), colour = color_vert, hjust = 0) +
      ggplot2::geom_point (x= range_faxes_lim[1] + spread_faxes*0.1,
                           y= range_faxes_lim[2] - spread_faxes*0.35,
                           colour = color_vert, fill = fill_vert, shape = shape_vert, size = size_vert) +
      
      ggplot2::geom_text( x = range_faxes_lim[1] + spread_faxes*0.15,
                          y = range_faxes_lim[2] - spread_faxes*0.65,
                          label = paste0("in the ", ncol(sp_faxes_coord), "-dimensional space"), hjust = 0)
    
  }
  
  
  #### arranging panels ####
  
  # if 2 axes = 1 plot + caption
  if(plot_nb == 1) {
    patchwork_plots_all <- panels[[1]] + plot_caption +
      patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1, 1),
                             ncol = 2, nrow = 1, guides = "collect")
    patchwork_plots_all <- patchwork_plots_all +
      patchwork::plot_annotation(title = "Position of species along functional axes",
                                 caption = "made with mFD package")
    panels[[2]] <- patchwork_plots_all
  }
  
  # if 3 axes = 3 plots + caption in a 2*2 layout
  if(plot_nb == 3) {
    patchwork_plots_all <- (panels[[1]] + plot_caption +
                              panels[[2]] + panels[[3]]) +
      patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1, 1),
                             ncol = 2, nrow = 2, guides = "collect")
    patchwork_plots_all <- patchwork_plots_all +
      patchwork::plot_annotation(title = "Position of species along functional axes",
                                 caption = "made with mFD package")
    panels[[4]] <- patchwork_plots_all
  }
  
  # if 4 axes = 6 plots + caption in a 3*3 layout with 2 empty cases
  if (plot_nb == 6) {
    patchwork_plots_all <- (panels[[1]] + patchwork::plot_spacer() + plot_caption +
                              panels[[2]] + panels[[3]] + patchwork::plot_spacer() +
                              panels[[4]] + panels[[5]] + panels[[6]]) +
      patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), widths = rep(1, 3),
                             ncol = 3, nrow = 3, guides = "collect")
    patchwork_plots_all <- patchwork_plots_all +
      patchwork::plot_annotation(title = "Position of species along functional axes",
                                 caption = "made with mFD package")
    panels[[7]] <- patchwork_plots_all
  }
  
  
  
  
  #### returning output ####
  
  # type, resolution and dimensions of file if to be saved
  device_file = "png"
  res_file = 300
  height_file <- 4*c(1,2,3) ; names(height_file) <- c("1", "3", "6")
  width_file <- 4*c(2,2,3) ; names(width_file) <- c("1", "3", "6")
  
  # displaying or saving
  if (is.null(name_file) == TRUE)  {
    patchwork_plots_all
    
  } else  {
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file) ,
                    plot = patchwork_plots_all,
                    device = device_file,
                    scale = 1,
                    height = height_file[as.character(plot_nb)],
                    width = width_file[as.character(plot_nb)],
                    units = "in",
                    dpi = res_file)
  }
  
  
  # output = list of panels and patchwork object:
  return(panels)
  
} #### end of function ####



