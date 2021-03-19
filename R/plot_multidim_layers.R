#' Plot background of multidimensional plots
#'
#' This function creates a ggplot object with customized axes range 
#' (same for both), names and background
#'
#' @param range_faxes a vector with minimum and maximum values of axes. Note
#'   that to have a fair representation of position of species in all plots,
#'   they should have the same range. Default: `range_faxes = c(NA, NA)` (the 
#'   range is computed according to the range of values among all axes).
#' 
#' @param faxes_nm a vector with axes labels for figure.
#' 
#' @param color_bg a R color name  or an hexadecimal code used to fill plot
#'  background. Default: `color_bg = "grey95"`.
#'
#' @return a ggplot object plotting background of multidimensional graphs
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'  background <- background.plot(range_faxes = c(-1, 2), 
#'                                faxes_nm    = c("PC 1", "PC 2"), 
#'                                color_bg    = "grey90") 
#' }
#'



background.plot <- function(range_faxes, faxes_nm, color_bg) {
  
    ggplot_bg <- ggplot2::ggplot( ) +
      ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::xlab(faxes_nm[1]) +
      ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
      ggplot2::ylab(faxes_nm[2]) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = color_bg))+
      ggplot2::coord_fixed()
    
    return(ggplot_bg)
    
  
} 




#' Plot species from the pool
#' 
#' Plot all species from the study case and associated convex hull
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in the pool for a given pair of axes
#' 
#' @param vertices_nD a list (with names as in sp_coord2D) of vectors with 
#' names of species being vertices in n dimensions.
#' 
#' @param plot_pool a logical value indicating whether species of each 
#' assemblage should be plotted or not. Default: plot_pool = TRUE.
#' 
#' @param color_ch a R color name or an hexadecimal code referring to the border
#'  of the convex hull filled by the pool of species. Default: 
#'  `color_ch = "black"`. If several assemblages it should be a vector with 
#'  names as in `sp_coord2D`.
#' 
#' @param color_pool a R color name or an hexadecimal code referring to the 
#' color of the pool.  This color is also used for FRic convex hull color. 
#' Default: `color_pool = "#0072B2"`.
#' 
#' @param color_vert a R color name or an hexadecimal code referring to the
#'   color of vertices if plotted. If color_vert = NA, vertices are not plotted
#'   (for shapes only defined by color, ie shape inferior to 20. Otherwise fill
#'   must also be set to NA). Default: `color_vert =  NA`. 
#'   
#' @param fill_ch a R color name or an hexadecimal code referring to the filling
#'  of the convex hull filled by the pool of species. Default is: 
#'  `fill_ch = "white"`. If several assemblages it should be a vector with 
#'  names as in `sp_coord2D`.
#' 
#' @param fill_pool a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20) and the assemblage convex
#'  hull. Default: `fill_pool = '#0072B2'`.
#' 
#' @param fill_vert a character value referring to the color for filling symbol
#'   for vertices (if \code{shape_vert} >20). If `fill = NA` and `color = NA`,
#'   vertices are not plotted (if \code{shape_vert} superior to 20. Otherwise
#'   `color_vert = NULL` is enough). Default is `NA`. If several assemblages 
#'   it should be a vector with names as in `sp_coord2D`.
#' 
#' @param shape_pool a numeric value referring to the shape used to plot species
#'  pool. Default: `shape_pool = 16`(filled circle). 
#' 
#' @param shape_vert a numeric value referring to the shape used to plot
#'   vertices if vertices should be plotted in a different way than other
#'   species. If `shape_vert = NA`, no vertices plotted. Default: 
#'   `shape_vert = NA`. If several assemblages it should be a vector with 
#'   names as in `sp_coord2D`.
#' 
#' @param size_pool a numeric value referring to the size of species belonging 
#' to the global pool. Default: `size_pool = 1`.
#' 
#' @param size_vert a numeric value referring to the size of symbol for vertices
#'  Default: `size_vert = 1`.  If several assemblages it should be a vector 
#'  with names as in `sp_coord2D`.
#' 
#' @param alpha_ch a numeric value for transparency of the filling of the convex
#'  hull (0 = high transparency, 1 = no transparency). Default: 
#'  `alpha_ch = 0.3`.
#'
#' @return a ggplot object plotting background of multidimensional graphs and
#' species from the global pool (associated convex hull if asked)
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
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
#' sp_faxes_coord_fruits_2D <- 
#'   fspaces_quality_fruits$details_fspaces$sp_pc_coord[, c("PC1", "PC2")]
#' 
#' # Set faxes limits:
#' # set range of axes if c(NA, NA):
#'  range_sp_coord_fruits  <- range(sp_faxes_coord_fruits_2D)
#'  range_faxes_lim <- range_sp_coord_fruits + 
#'    c(-1, 1)*(range_sp_coord_fruits[2] - 
#'    range_sp_coord_fruits[1]) * 0.05
#'  
#'  # Retrieve the background plot:
#'  ggplot_bg_fruits <- mFD::background.plot(
#'                                range_faxes = ranges_faxes_lim, 
#'                                faxes_nm    = c("PC 1", "PC 2"), 
#'                                color_bg    = "grey90") 
#'                                
#'  # Retrieve vertices names:
#'  vert_nm_fruits <- vertices(sp_faxes_coord_fruits, 
#'   order_2D = TRUE, check_input = TRUE)
#'   
#'  # Plot the pool:
#'  plot_pool_fruits <- pool.plot(ggplot_bg   = ggplot_bg_fruits,
#'            sp_coord2D    = sp_coord_fruits_2D,
#'            vertices_nD   = vert_nm_fruits,
#'            plot_pool     = TRUE,
#'            shape_pool    = 3, 
#'            size_pool     = 0.8, 
#'            color_pool    = "grey95", 
#'            fill_pool     = NA, 
#'            color_ch      = NA, 
#'            fill_ch       = "white", 
#'            alpha_ch      = 1, 
#'            shape_vert    = 3, 
#'            size_vert     = 1, 
#'            color_vert    = "black", 
#'            fill_vert     = NA) 
#'  plot_pool_fruits
#'                                
#' }

pool.plot <-function(ggplot_bg,
                    sp_coord2D,
                    vertices_nD,
                    plot_pool = TRUE,
                    shape_pool = 3, size_pool = 0.8, color_pool = "grey95", 
                    fill_pool = NA, color_ch = NA, fill_ch = "white", 
                    alpha_ch = 1, shape_vert = 3, size_vert = 1, 
                    color_vert = "black", fill_vert = NA) {
  
  
  # prepare data for plotting ####
  
  # dataframe with species coordinates and vertex status in nD:
  sp_xyv <- data.frame(x = sp_coord2D[, 1], y = sp_coord2D[, 2], vert = "sp")
  
  
  # if required vertices and convex hull:
  if (! is.null(vertices_nD)) {
    sp_xyv[vertices_nD, "vert"] <- "vert"
    
    # coordinates of species vertices in nD:
    vertnD_xy <- sp_xyv[sp_xyv$vert == "vert", c("x", "y")]
    
    #  species being vertices of the convex hull in 2D:
    vert2D <- vertices(as.matrix(vertnD_xy), order_2D = TRUE)
    
  } # end of if vertices
  
  
  # plotting layers ####
  
  # default plot is background:
  ggplot_pool <- ggplot_bg
  
  # plotting convex hull if required:
  x <- NULL
  y <- NULL
  if (! is.null(vertices_nD)) {
    ggplot_pool <- ggplot_pool +
      ggplot2::geom_polygon(data = vertnD_xy[vert2D, ],
                            ggplot2::aes(x = x, y = y) ,
                            colour = color_ch, fill = fill_ch, 
                            alpha = alpha_ch)
  } # end of if vertices
  
  # plotting species if required:
  if (plot_pool) {
    
    # plotting species not being vertices
    x <- NULL
    y <- NULL
    ggplot_pool <- ggplot_pool +
      ggplot2::geom_point(data = sp_xyv[sp_xyv$vert == "sp", ],
                          ggplot2::aes(x = x, y = y),
                          size = size_pool, shape = shape_pool,
                          colour = color_pool, fill = fill_pool)
    
    # plotting species being vertices
    if (! is.null(vertices_nD)) {
      ggplot_pool <- ggplot_pool +
        ggplot2::geom_point(data = sp_xyv[sp_xyv$vert == "vert", ],
                            ggplot2::aes(x = x, y = y),
                            size = size_vert, shape = shape_vert,
                            colour = color_vert, fill = fill_vert)
    } # end of if vertices
    
  } # end of if plot species
  
  return(ggplot_pool)
  
}




#' Plot species
#' 
#' Plot species for one to n assemblages
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of axes.
#' 
#' @param asb_vertices_nD a list (with names as in asb_sp_coord2D) of vectors 
#' with names of species being vertices in n dimensions.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param size_sp a numeric value referring to the size of species belonging to
#'  the plotted assemblage. Default: `size_sp = 1`. If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param size_vert a numeric value referring to the size of symbol for vertices
#'  Default: `size_vert = 1`.  If several assemblages it should be a vector 
#'  with names as in `asb_sp_coord2D`.
#' 
#' @param shape_sp a numeric value referring to the shape of species 
#' belonging to the plotted assemblage. 
#' Default: `shape_sp = 16`. If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param shape_vert a numeric value referring to the shape of symbol for 
#' vertices Default: `shape_vert = NA`.  If several assemblages it should be 
#' a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage.  This color is also used for FRic
#'  convex hull color. Default: `color_sp = "#0072B2"`. If several assemblages 
#'  it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_vert a R color name or an hexadecimal code referring to the
#'   color of vertices if plotted. If color_vert = NA, vertices are not plotted
#'   (for shapes only defined by color, ie shape inferior to 20. Otherwise fill
#'   must also be set to NA). Default: `color_vert =  NA`. If several 
#'   assemblages it should be a vector with names as in `asb_sp_coord2D`.
#'
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20) and the assemblage convex
#'  hull. Default: `fill_sp = '#0072B2'`. If several assemblages it should 
#'  be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param fill_vert a character value referring to the color for filling symbol
#'   for vertices (if \code{shape_vert} >20). If `fill = NA` and `color = NA`,
#'   vertices are not plotted (if \code{shape_vert} superior to 20. Otherwise
#'   `color_vert = NULL` is enough). Default is `NA`. If several assemblages 
#'   it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param limits_relatw a vector of two numbers giving the limits to set for 
#' the scale of species relative weights.
#' 
#' @param range_size_relatw a vector of two numbers specifying the minimum and 
#' the maximum size of the plotting symbol for relative weights.
#' 
#' @return a ggplot object with species plotted on the background plot
#' 
#' @note if \code{asb_vertices_nD = NULL}, all arguments for vertices are 
#' ignored and if \code{asb_sp_relatw != NULL}, \code{size_sp} and 
#' \code{size_vert} are ignored.
#' 
#' 

sp.plot <-function(ggplot_bg,
                  asb_sp_coord2D,
                  asb_vertices_nD = NULL,
                  asb_sp_relatw = NULL,
                  size_sp, shape_sp, color_sp, fill_sp,
                  size_vert, shape_vert, color_vert, fill_vert,
                  limits_relatw = c(0, 1),
                  range_size_relatw = c(1, 10)) {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  # dataframe with variables needed for plotting ####
  asb_sp_xywv <- NULL
  
  # loop on assemblages:
  for (z in asb_nm) {
    
    x <- NULL
    y <- NULL
    w <- NULL
    vert <- NULL
    asb <- NULL
    # species coordinates and name of assemblage: 
    data_z <- data.frame(x = asb_sp_coord2D[[z]][, 1],
                         y = asb_sp_coord2D[[z]][, 2],
                         w = 0,
                         vert = "no",
                         asb = z)
    
    
    # adding species weight:
    if (! is.null (asb_sp_relatw))  {
      data_z[, "w"] <- asb_sp_relatw[[z]][row.names(data_z)] 
    } else { 
      data_z[, "w"] <- size_sp[[z]]/100
    }

    # if needed, adding vertex status and size if fixed
    if (! is.null (asb_vertices_nD)) {
      data_z[asb_vertices_nD[[z]], "vert"] <- "vert"
      
      if(is.null (asb_sp_relatw)) {
        data_z[asb_vertices_nD[[z]], "w"] <- size_vert[[z]]/100
      }
    }
    
    # row binding:
    asb_sp_xywv <- rbind(asb_sp_xywv, data_z)
    
  } # end of z
  
  
  # define aesthetics ####
  
  # add new variable mixing assemblage and vertex status:
  asb_vert <- NULL
  asb_sp_xywv$asb_vert <- paste(asb_sp_xywv$asb, asb_sp_xywv$vert, sep = "_")
  
  # define aesthetics for species according to asb and vertex status:
  lev_asb_vert <- c(paste(asb_nm, "no", sep = "_"),
                    paste(asb_nm, "vert", sep = "_"))
  
  # if no vertices, setting aesthetics as for species:
  if (is.null(asb_vertices_nD)) {
    shape_vert <- shape_sp
    color_vert <- color_sp
    fill_vert <- fill_sp
  }
  
  shape_asb_vert <- c(shape_sp, shape_vert)
  names(shape_asb_vert) <- lev_asb_vert
  
  color_asb_vert <- c(color_sp, color_vert)
  names(color_asb_vert) <- lev_asb_vert
  
  fill_asb_vert <- c(fill_sp, fill_vert)
  names(fill_asb_vert) <- lev_asb_vert
  
  # reorder species according to decreasing weight:
  asb_sp_xywv <- asb_sp_xywv[order(asb_sp_xywv$w, decreasing = TRUE), ]
  
  
  ## plotting ####
  
  # default plot is background and set size of points:
  ggplot_sp <- ggplot_bg +
    ggplot2::scale_size(limits = limits_relatw, range = range_size_relatw)
  
  
  # plot species as points with chosen shape, size and colors:
  ggplot_sp <- ggplot_sp +
    ggplot2::geom_point(data = asb_sp_xywv,
                        ggplot2::aes(x = x, y = y, size = w,
                            shape = asb_vert,
                            colour = asb_vert, fill = asb_vert),
                        show.legend = FALSE) +
    ggplot2::scale_shape_manual(name = "asb_vert", values = shape_asb_vert) +
    ggplot2::scale_colour_manual(name = "asb_vert", values = color_asb_vert) +
    ggplot2::scale_fill_manual(name = "asb_vert", values = fill_asb_vert)
  
  return(ggplot_sp)
  
}





#' Plot FRic index
#'
#' This function plots FRic index for a given pair of functional axes and one 
#' or several assemblages. It adds convex hull(s), points and vertices of 
#' 1:N assemblages on the background plot
#'
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of axes for a given pair 
#' of functional axes.
#' 
#' @param asb_vertices_nD a list (with names as in asb_sp_coord2D) of vectors 
#' with names of species being vertices in n dimensions.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: plot_sp = TRUE.
#' 
#' @param color_ch a R color name or an hexadecimal code referring to the border
#'  of the convex hull filled by the pool of species. Default: 
#'  `color_ch = "black"`. If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage.  This color is also used for FRic
#'  convex hull color. Default: `color_sp = "#0072B2"`. If several assemblages 
#'  it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_vert a R color name or an hexadecimal code referring to the
#'   color of vertices if plotted. If color_vert = NA, vertices are not plotted
#'   (for shapes only defined by color, ie shape inferior to 20. Otherwise fill
#'   must also be set to NA). Default: `color_vert =  NA`. If several 
#'   assemblages it should be a vector with names as in `asb_sp_coord2D`.
#'
#' @param fill_ch a R color name or an hexadecimal code referring to the filling
#'  of the convex hull filled by the pool of species. Default is: 
#'  `fill_ch = "white"`. If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20) and the assemblage convex
#'  hull. Default: `fill_sp = '#0072B2'`. If several assemblages it should 
#'  be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param fill_vert a character value referring to the color for filling symbol
#'   for vertices (if \code{shape_vert} >20). If `fill = NA` and `color = NA`,
#'   vertices are not plotted (if \code{shape_vert} superior to 20. Otherwise
#'   `color_vert = NULL` is enough). Default is `NA`. If several assemblages 
#'   it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param shape_vert a numeric value referring to the shape used to plot
#'   vertices if vertices should be plotted in a different way than other
#'   species. If `shape_vert = NA`, no vertices plotted. Default: 
#'   `shape_vert = NA`. If several assemblages it should be a vector with 
#'   names as in `asb_sp_coord2D`.
#' 
#' @param size_sp a numeric value referring to the size of species belonging to
#'  the plotted assemblage. Default: `size_sp = 1`. If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param size_vert a numeric value referring to the size of symbol for vertices
#'  Default: `size_vert = 1`.  If several assemblages it should be a vector 
#'  with names as in `asb_sp_coord2D`.
#' 
#' @param alpha_ch a numeric value for transparency of the filling of the convex
#'  hull (0 = high transparency, 1 = no transparency). Default: 
#'  `alpha_ch = 0.3`.
#'
#' @return a ggplot object plotting background of multidimensional graphs and
#' FRic convex hulls
#' 
#' @note first convex hull (2D projection) of all assemblages are added to 
#' the background ggplot object (`ggplot_bg`), then species not being vertices 
#' in n Dimensions-space and lastly species being vertices are plotted
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
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
#'   sp_faxes_coord_fruits <- 
#'     fspaces_quality_fruits$details_fspaces$sp_pc_coord
#'     
#'  # Retrieve species coordinates matrix for the assemblage "basket_1":
#'   sp_filter <- mFD::sp.filter(asb_nm          = c("basket_1"), 
#'                               sp_faxes_coord = sp_faxes_coord_fruits, 
#'                               asb_sp_w       = baskets_fruits_weight)
#'   sp_faxes_coord_fruits_b1 <- sp_filter$`species coordinates`
#'  
#'  # Reduce it to the two studid axes: PC1 and PC2:
#'  sp_faxes_coord_fruits_b1_2D <- sp_faxes_coord_fruits_b1[, c("PC1", "PC2")]
#'
#' # Set faxes limits:
#' # set range of axes if c(NA, NA):
#'  range_sp_coord  <- range(sp_faxes_coord_fruits)
#'  range_faxes_lim <- range_sp_coord + c(-1, 1)*(range_sp_coord_fruits[2] - 
#'  range_sp_coord_fruits[1]) * 0.05
#'  
#'  # Retrieve the background plot:
#'  ggplot_bg_fruits <- mFD::background.plot(
#'                                range_faxes = ranges_faxes_lim, 
#'                                faxes_nm    = c("PC 1", "PC 2"), 
#'                                color_bg    = "grey90") 
#'                                
#'  # Retrieve vertices names:
#'  vert_nm_fruits <- vertices(sp_faxes_coord_fruits_b1, 
#'   order_2D = TRUE, check_input = TRUE)
#'                                
#'  # Plot in white the convex hull of all fruits species:
#'  ggplot_fric <- mFD::fric.plot(
#'            ggplot_bg       = ggplot_bg_fruits,
#'            asb_sp_coord2D  = list(basket_1 = sp_faxes_coord_fruits_b1_2D),
#'            asb_vertices_nD = vert_nm_fruits,
#'            plot_sp         = TRUE,
#'            color_ch        = "black", 
#'            fill_ch         = "white", 
#'            alpha_ch        = 0.3,
#'            size_sp = 1,
#'            shape_sp = 16,
#'            color_sp = "red",
#'            fill_sp = "red",
#'            size_vert = 1,
#'            color_vert = "red",
#'            fill_vert = "red",
#'            shape_vert = 16)
#'  ggplot_fric
#'
#' }
#'
#'

fric.plot <- function(ggplot_bg,
                      asb_sp_coord2D,
                      asb_vertices_nD,
                      plot_sp = TRUE,
                      color_ch, fill_ch, alpha_ch,
                      shape_sp, size_sp, color_sp, fill_sp,
                      shape_vert, size_vert, color_vert, fill_vert) {
  
  # names of assemblages
  asb_nm <- names(asb_sp_coord2D)
  
  ## plotting layers ####
  
  # default plot is background
  ggplot_fric <- ggplot_bg
  
  # computing and plotting 2D convex hull(s)
  for (k in asb_nm) {
    # coordinates of species vertices in nD
    k_vertnD_xy <- asb_sp_coord2D[[k]][asb_vertices_nD[[k]], ]
    
    # species being vertices of the convex hull in 2D
    k_vert2D <- vertices(k_vertnD_xy, order_2D = TRUE)
    
    # dataframe with coordinates of vertices
    x <- NULL
    y <- NULL
    k_vertnD_xy <- data.frame(x = k_vertnD_xy[k_vert2D, 1],
                              y = k_vertnD_xy[k_vert2D, 2])
    
    # plotting convex hulls
    ggplot_fric <- ggplot_fric +
      ggplot2::geom_polygon(data = k_vertnD_xy,
                            ggplot2::aes(x = x, y = y),
                            colour = color_ch[k],
                            fill = fill_ch[k],
                            alpha = alpha_ch[k])
  }
  
  
  # plotting species if required ----
  if(plot_sp) {
    ggplot_fric <- sp.plot(ggplot_bg = ggplot_fric,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = asb_vertices_nD,
                           asb_sp_relatw = NULL,
                           shape_sp = shape_sp,
                           size_sp = size_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp,
                           shape_vert = shape_vert,
                           size_vert = size_vert,
                           color_vert = color_vert,
                           fill_vert = fill_vert)
    
  }# end of if plot_sp

  return(ggplot_fric)
  
}




#' Plot FDiv indice
#'
#' This plot fDiv indice for a given pair of functional axes and one or several 
#' assemblages. This function adds mean distance to center of gravity of 
#' vertices, points and vertices of 1:N assemblages on the background plot
#'
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_vertices_nD a list (with names as in asb_sp_coord2D) of vectors
#'   with names of species being vertices in n dimensions.
#' 
#' @param asb_vertG_coord2D a list (with names as in asb_sp_coord2D) containing 
#' for each assemblage the coordinates of center of gravity of vertices for a 
#' given pair of axes
#' 
#' @param asb_meanDtoG a list (with names as in asb_sp_coord2D) containing 
#' for each assemblage the mean distance to the center of gravity of vertices
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: plot_sp = TRUE.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#' 
#' @param color_vert a R color name or an hexadecimal code referring to the
#'   color of vertices if plotted. If color_vert = NA, vertices are not plotted
#'   (for shapes only defined by color, ie shape inferior to 20. Otherwise fill
#'   must also be set to NA). Default: `color_vert =  NA`. If several 
#'   assemblages it should be a vector with names as in `asb_sp_coord2D`.
#'   
#' @param color_vertG a R color name or an hexadecimal code referring to the
#'   color of the center of gravity of vertices. If several assemblages it 
#'   should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_meanD a R color name or an hexadecimal code referring to the
#'   color of the circle representing the mean distance of species to
#'   the center of gravity of the vertices. If several assemblages it should be 
#'   a vector with names as in `asb_sp_coord2D`.
#' 
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should 
#'  be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param fill_vert a R color name or an hexadecimal code referring to the 
#'   color for filling symbol
#'   for vertices (if \code{shape_vert} >20). If `fill = NA` and `color = NA`,
#'   vertices are not plotted (if \code{shape_vert} superior to 20. Otherwise
#'   `color_vert = NULL` is enough). Default is `NA`. If several assemblages 
#'   it should be a vector with names as in `asb_sp_coord2D`.
#'   
#' @param fill_vertG a R color name or an hexadecimal code 
#' referring to the color to fill the center of gravity of vertices. If several 
#' assemblages it should be a vector with names as in `asb_sp_coord2D`.
#'  
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param shape_vert a numeric value referring to the shape used to plot
#'   vertices if vertices should be plotted in a different way than other
#'   species. If `shape_vert = NA`, no vertices plotted. Default: 
#'   `shape_vert = NA`. If several assemblages it should be a vector with 
#'   names as in `asb_sp_coord2D`.
#'   
#' @param shape_vertG a numeric value referring to the shape to use to
#' plot the center of gravity of vertices. If several assemblages it should be 
#' a vector with names as in `asb_sp_coord2D`.
#'  
#' @param size_vertG a numeric value referring to the size to use to
#' plot the center of gravity of vertices. If several assemblages it should be 
#' a vector with names as in `asb_sp_coord2D`.
#'  
#'
#' @return a ggplot object plotting background of multidimensional graphs and
#' FDiv indice
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
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
#' # Set faxes limits:
#' # set range of axes if c(NA, NA):
#'  range_sp_coord_fruits  <- range(sp_faxes_coord_fruits)
#'  range_faxes_lim <- range_sp_coord_fruits + 
#'  c(-1, 1)*(range_sp_coord_fruits[2] - 
#'  range_sp_coord_fruits[1]) * 0.05
#'  
#'  # Retrieve the background plot:
#'  ggplot_bg_fruits <- mFD::background.plot(
#'                                range_faxes = range_faxes_lim, 
#'                                faxes_nm    = c("PC 1", "PC 2"), 
#'                                color_bg    = "grey90") 
#'                                
#'  # Retrieve the matrix of species coordinates for "basket_1" and PC1 and PC2:
#'  sp_filter <- mFD::sp.filter(asb_nm          = "basket_1", 
#'                              sp_faxes_coord = sp_faxes_coord_fruits, 
#'                              asb_sp_w       = baskets_fruits_weights)
#'  fruits_asb_sp_coord_b1 <- sp_filter$`species coordinates`
#'  fruits_asb_sp_coord2D_b1 <- fruits_asb_sp_coord_b1[, c("PC1", "PC2")]
#'                                
#'  # Use alpha.fd.multidim() function to get inputs to plot FDiv:
#'  alpha_fd_indices_fruits <- mFD::alpha.fd.multidim(
#'   sp_faxes_coord    = sp_faxes_coord_fruits[, c("PC1", "PC2", "PC3", "PC4")],
#'   asb_sp_w         = baskets_fruits_weights,
#'   ind_vect         = c("fdiv"),
#'   scaling          = TRUE,
#'   check_input      = TRUE,
#'   details_returned = TRUE)
#'   
#'  # Retrieve inputs of the fdiv.plot() function for "basket_1" and PC1, PC2...
#'  # ... through alpha.fd.multidim outputs:
#'  fruits_asb_sp_relatw_b1 <- 
#'          alpha_fd_indices_fruits$details$asb_sp_relatw["basket_1", ]
#'  fruits_asb_vertices_nD_b1_2D <- 
#'                       alpha_fd_indices_fruits$details$asb_vert_nm["basket_1"]
#'  fruits_asb_vertG_coord_b1 <- 
#'                       alpha_fd_indices_fruits$details$asb_G_coord["basket_1"]
#'  fruits_asb_vertG_coord_b1_2D <- 
#'              fruits_asb_vertG_coord_b1[[1]][c("PC1", "PC2")]
#'  fruits_asb_meanDtoG_b1 <- 
#'                   alpha_fd_indices_fruits$details$asb_mean_dist_G["basket_1"]
#'  
#'  # Retrieve FDiv plot:
#'  fdiv_plot <- fdiv.plot(
#'            ggplot_bg         = ggplot_bg_fruits,
#'            asb_sp_coord2D    = list(basket_1 = fruits_asb_sp_coord2D_b1),
#'            asb_sp_relatw     = list(basket_1 = fruits_asb_sp_relatw_b1),
#'            asb_vertices_nD   = fruits_asb_vertices_nD,
#'            asb_vertG_coord2D = list(basket_1 = fruits_asb_vertG_coord_b1_2D),
#'            asb_meanDtoG      = fruits_asb_meanDtoG,
#'            plot_sp           = TRUE,
#'            shape_sp          = 16,
#'            color_sp          = "red",
#'            fill_sp           = "red",
#'            color_vert        = "red",
#'            fill_vert         = "red",
#'            shape_vert        = 16,
#'            shape_vertG       = list(basket_1 = 18),
#'            size_vertG        = list(basket_1 = 2),
#'            color_vertG       = list(basket_1 = "blue"),
#'            fill_vertG        = list(basket_1 = "blue"),
#'            color_meanD       = list(basket_1 = "red"))
#'  fdiv_plot
#'  
#'
#' }
#'
#'


fdiv.plot <- function(ggplot_bg,
                      asb_sp_coord2D,
                      asb_sp_relatw,
                      asb_vertices_nD,
                      asb_vertG_coord2D,
                      asb_meanDtoG,
                      plot_sp = TRUE,
                      shape_sp, color_sp, fill_sp,
                      shape_vert = NA, color_vert = NA, fill_vert = NA,
                      shape_vertG, size_vertG,
                      color_vertG, fill_vertG,
                      color_meanD) {
  
  # get the names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  # prepare data for plotting ####
  
  # list of dataframe with coordinates of center of gravity of vertices...
  # ... and mean distance to it:
  asb_vertG_xyr <- list()
  x <- NULL
  y <- NULL
  r <- NULL
  
  for (z in asb_nm) {
    asb_vertG_xyr[[z]] <- data.frame(x = asb_vertG_coord2D[[z]][1],
                                     y = asb_vertG_coord2D[[z]][2],
                                     r = asb_meanDtoG[[z]])
    rownames(asb_vertG_xyr[[z]]) <- z
  }
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_fdiv <- ggplot_bg
  
  # plotting all species if required:
  if (plot_sp) {
    ggplot_fdiv <- sp.plot(ggplot_bg = ggplot_fdiv,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = asb_vertices_nD,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp,
                           shape_vert = shape_vert,
                           color_vert = color_vert,
                           fill_vert = fill_vert)
  }
  
  # plotting center of gravity of vertices in nD and mean distance to it:
  for (k in asb_nm) {
    x <- NULL
    y <- NULL
    ggplot_fdiv <- ggplot_fdiv +
      ggplot2::geom_point(data = asb_vertG_xyr[[k]], 
                          ggplot2::aes(x = x , y = y),
                          colour = color_vertG[[k]], fill = fill_vertG[[k]],
                          shape = shape_vertG[[k]], size = size_vertG[[k]]) +
      ggforce::geom_circle(data = asb_vertG_xyr[[k]], ggplot2::aes(x0 = x, 
                                                                   y0 = y, 
                                                                    r = r),
                           colour = color_meanD[[k]],
                           show.legend = FALSE, inherit.aes = FALSE)
  }
  
  return(ggplot_fdiv)
  
}




#' Plot FIde index
#' 
#' Plot FIde index given a pair of functional axes for one or several 
#' assemblages. It adds the centroid of species for each assemblage and 
#' segments showing centroid coordinates on functional axes on the background 
#' plot.
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_fide_coord2D a list (with names as in asb_sp_coord2D) of 
#' vectors with coordinates of the centroid of species for each assemblage
#' for a given pair of axes.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#'  
#' @param color_centroid a R color name or an hexadecimal code referring to the 
#' color of the species centroid from the studied assemblage. If several 
#' assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_segment a R color name or an hexadecimal code referring to the 
#' color of the segment linking axes and centroid from the studied assemblage. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param fill_fide a R color name or an hexadecimal code referring to 
#' the colour to fill assemblage centroid symbol (if \code{shape_sp} > 20).
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param shape_fide a numeric value referring to the shape used to 
#' plot fide centroid of the studied assemblage. Default: `shape_centroid = 18` 
#'  (filled diamond). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param size_fide a numeric value referring to the size of species fide
#' centroid but not the plotted assemblage. Default: `size_sp = 1`. If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param width_segment a numeric value referring to the size of the segment 
#' linking fide centroid and functional axes. Default: `width_segment = 1`.  If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param linetype_segment a character string referring to the linetype used to 
#' draw the segment linking fide centroid and functional axes. 
#' Default: `linetype_segment = "dashed"`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object with FIde index, species and background
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' 
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
#' # Set faxes limits:
#' # set range of axes if c(NA, NA):
#'  range_sp_coord_fruits  <- range(sp_faxes_coord_fruits)
#'  range_faxes_lim <- range_sp_coord_fruits + 
#'  c(-1, 1)*(range_sp_coord_fruits[2] - 
#'  range_sp_coord_fruits[1]) * 0.05
#'  
#'  # Retrieve the background plot:
#'  ggplot_bg_fruits <- mFD::background.plot(
#'                                range_faxes = range_faxes_lim, 
#'                                faxes_nm    = c("PC 1", "PC 2"), 
#'                                color_bg    = "grey90") 
#'                                
#'  # Retrieve the matrix of species coordinates for "basket_1" and PC1 and PC2:
#'  sp_filter <- mFD::sp.filter(asb_nm          = "basket_1", 
#'                              sp_faxes_coord = sp_faxes_coord_fruits, 
#'                              asb_sp_w       = baskets_fruits_weights)
#'  fruits_asb_sp_coord_b1 <- sp_filter$`species coordinates`
#'  fruits_asb_sp_coord2D_b1 <- fruits_asb_sp_coord_b1[, c("PC1", "PC2")]
#'                                
#'  # Use alpha.fd.multidim() function to get inputs to plot FIde:
#'  alpha_fd_indices_fruits <- mFD::alpha.fd.multidim(
#'   sp_faxes_coord    = sp_faxes_coord_fruits[, c("PC1", "PC2", "PC3", "PC4")],
#'   asb_sp_w         = baskets_fruits_weights,
#'   ind_vect         = c("fide"),
#'   scaling          = TRUE,
#'   check_input      = TRUE,
#'   details_returned = TRUE)
#'   
#'  # Retrieve fide values through alpha.fd.multidim outputs:
#'  fruits_asb_fide_coord2D <- 
#'   alpha_fd_indices_fruits$functional_diversity_indices[c("fide_PC1", 
#'                                                          "fide_PC2")]
#'  fruits_asb_fide_coord2D_b1 <- fruits_asb_fide_coord2D[c("basket_1"), ]
#'  fruits_asb_sp_relatw_b1 <- 
#'          alpha_fd_indices_fruits$details$asb_sp_relatw["basket_1", ]
#'
#'
#'  
#'  # Retrieve FIde plot:
#'  fide_plot <- fide.plot(ggplot_bg = ggplot_bg_fruits,
#'            asb_sp_coord2D = list(basket_1 = fruits_asb_sp_coord2D_b1),
#'            asb_sp_relatw = list(basket_1 = fruits_asb_sp_relatw_b1),
#'            asb_fide_coord2D = list(basket_1 = fruits_asb_fide_coord2D_b1),
#'            plot_sp = TRUE,
#'            shape_sp = 16,
#'            color_sp = "red",
#'            fill_sp = "red",
#'            shape_fide = list(basket_1 = 18),
#'            size_fide = list(basket_1 = 5),
#'            color_fide = list(basket_1 = "blue"),
#'            fill_fide = list(basket_1 = "blue"),
#'            color_segment = list(basket_1 = "red"),
#'            width_segment = list(basket_1 = 1),
#'            linetype_segment = list(basket_1 = "dashed"))
#'            
#'  fide_plot
#'  
#' }


fide.plot <-function(ggplot_bg,
                    asb_sp_coord2D,
                    asb_sp_relatw,
                    asb_fide_coord2D,
                    plot_sp = TRUE,
                    shape_sp, color_sp, fill_sp,
                    shape_fide, size_fide,
                    color_fide, fill_fide,
                    color_segment, width_segment, linetype_segment) {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  
  ## prepare data for plotting:
  
  # actual limits of plot:
  ggplot_bg_lim <- 
    ggplot2::ggplot_build(ggplot_bg)$layout$panel_params[[1]]$x.range
  
  # list of dataframe with coordinates of center of gravity and links to axes:
  asb_centroid_xy <- list()
  asb_centroidtoaxes_xyxy <- list()
  
  for (z in asb_nm){
    
    cent_z <- asb_fide_coord2D[[z]]
    
    x <- NULL
    y <- NULL
    asb_centroid_xy[[z]] <- data.frame(x = cent_z[1, 1], y = cent_z[1, 2])
    
    asb_centroidtoaxes_xyxy[[z]] <- data.frame(
                          x = c(ggplot_bg_lim[1], cent_z[1, 1]),
                          y = c(cent_z[1, 2], ggplot_bg_lim[1]),
                          xend = rep(cent_z[1, 1], 2),
                          yend = rep(cent_z[1, 2], 2))
  }
  
  
  # plotting layers ####
  
  # default plot is background:
  ggplot_fide <- ggplot_bg
  
  # plot species if required:
  if (plot_sp) {
    ggplot_fide <- sp.plot(ggplot_bg = ggplot_fide,
                          asb_sp_coord2D = asb_sp_coord2D,
                          asb_vertices_nD = NULL,
                          asb_sp_relatw = asb_sp_relatw,
                          shape_sp = shape_sp,
                          color_sp = color_sp,
                          fill_sp = fill_sp)
  }
  
  
  # plotting projection of mean along axes:
  for (k in asb_nm) {
    xend <- NULL
    yend <- NULL
    ggplot_fide <- ggplot_fide +
      ggplot2::geom_segment(data = asb_centroidtoaxes_xyxy[[k]],
                            ggplot2::aes(x = x , 
                                         y = y, xend = xend, yend = yend),
                            colour = color_segment[[k]],
                            size = width_segment[[k]],
                            linetype = linetype_segment[[k]])
  }
  
  # plotting center of gravity of species:
  for (k in asb_nm) {
    ggplot_fide <- ggplot_fide +
      ggplot2::geom_point(data = asb_centroid_xy[[k]], 
                          ggplot2::aes(x = x , y = y),
                          colour = color_fide[[k]], fill = fill_fide[[k]],
                          shape = shape_fide[[k]], size = size_fide[[k]])
  }

  return(ggplot_fide)
  
}





#' Plot FDis index
#' 
#' This function plots FDis index for a given pair of functional axes and for
#' one or several assemblages. It adds segments between species relative weights
#' and assemblage cenntroid on the background plot.
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_centroid_coord2D a list (with names as in asb_sp_coord2D) of 
#' vectors with coordinates of the centroid of species for each assemblage
#' for a given pair of axes.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param shape_centroid a numeric value referring to the shape used to 
#' plot centroid of the studied assemblage. Default: `shape_centroid = 18` 
#'  (filled diamond). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#' 
#' @param color_centroid a R color name or an hexadecimal code referring to the 
#' color of the species centroid from the studied assemblage. If several 
#' assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param color_segment a R color name or an hexadecimal code referring to the 
#' color of the segment linking axes and centroid from the studied assemblage. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param fill_centroid a R color name or an hexadecimal code referring to 
#' the colour to fill assemblage centroid symbol (if \code{shape_sp} > 20).
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param size_centroid a numeric value referring to the size of species 
#' centroid but not the plotted assemblage. Default: `size_sp = 1`. If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param width_segment a numeric value referring to the size of the segment 
#' linking centroid and functional axes. Default: `width_segment = 1`.  If 
#'  several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param linetype_segment a character string referring to the linetype used to 
#' draw the segment linking centroid and functional axes. 
#' Default: `linetype_segment = "dashed`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object showing FDis index for one or several assemblage(s) 
#' and a given pair of functional axes.
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' }
#' 

fdis.plot <- function(ggplot_bg,
                    asb_sp_coord2D,
                    asb_sp_relatw,
                    asb_centroid_coord2D,
                    plot_sp = TRUE,
                    shape_sp, color_sp, fill_sp, 
                    shape_centroid, size_centroid,
                    color_centroid, fill_centroid,
                    color_segment, width_segment, linetype_segment) {
  
  # get names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  
  ## retrieve data for plotting:
  
  # list of dataframe with coordinates of center of gravity and links...
  # ... to species:
  
  asb_centroid_xy <- list()
  asb_sptocentroid_xyxy <- list()
  
  for (z in asb_nm) {
    
    sp_z <- asb_sp_coord2D[[z]]
    cent_z <- asb_centroid_coord2D[[z]]
    
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    
    asb_centroid_xy[[z]] <- data.frame(x = cent_z[1, 1], y = cent_z[1, 2])
    
    asb_sptocentroid_xyxy[[z]] <- data.frame(x = sp_z[, 1],
                                             y = sp_z[, 2],
                                             xend = cent_z[1, 1],
                                             yend = cent_z[1, 2])
  }
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_fdis <- ggplot_bg
  
  # plotting species if required:
  
  if(plot_sp) {
    ggplot_fdis <- sp.plot(ggplot_bg = ggplot_fdis,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = NULL,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp)
  }
  
  
  # plotting distance to center of gravity of species:
  for (k in asb_nm) {
    ggplot_fdis <- ggplot_fdis +
      ggplot2::geom_segment(data = asb_sptocentroid_xyxy[[k]],
                            ggplot2::aes(x = x , y = y, xend= xend, 
                                         yend = yend),
                            colour = color_segment[k],
                            size = width_segment[k],
                            linetype = linetype_segment[k])
  }
  
  # plotting center of gravity of species ----
  for (k in asb_nm) {
    ggplot_fdis <- ggplot_fdis +
      ggplot2::geom_point(data = asb_centroid_xy[[k]], 
                          ggplot2::aes(x = x , y = y),
                          colour = color_centroid[k], fill = fill_centroid[k],
                          shape = shape_centroid[k], size = size_centroid[k])
  }
  
  
  
  return(ggplot_fdis)
  
}




#' Plot FEve index
#' 
#' This function plots FEve index for a given pair of functional axes and one 
#' or several assemblages. It adds Minimum Spanning Tree (MST) of a given 
#' assemblage on the background plot.
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_mst a list (with names as in `asb_sp_coord2D`) of 
#' vectors with names of species linked in the MST of the 
#' studied assemblage. If several assemblages it should be a vector with names 
#' as in `asb_sp_coord2D`.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#'  
#' @param color_mst a R color name or an hexadecimal code referring to the color
#'  of the MST from the studied assemblage. Default: `color_MST = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#' 
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param width_mst a numeric value referring to the size of the line 
#' representing MST. Default: `width_mst = 1`.  If 
#' several assemblages it should be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param linetype_mst a character string referring to the linetype used to 
#' draw the MST. Default: `linetype_mst = "dashed`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object showing FEve index on the backrgound plot
#' 
#' @author Camille Magneville and sébastien Villéger
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' }
#'  
#'  

feve.plot <- function(ggplot_bg,
                      asb_sp_coord2D,
                      asb_sp_relatw,
                      asb_mst,
                      plot_sp = TRUE,
                      shape_sp, color_sp, fill_sp,
                      color_mst, width_mst, linetype_mst) {
  
  # get names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_feve <- ggplot_bg
   
  # plotting species if required:
  if (plot_sp) {
    ggplot_feve <- sp.plot(ggplot_bg = ggplot_feve,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = NULL,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp)
  }
  
  # plotting minimum spanning tree:
  for (k in asb_nm) {
    
    # branches of MST in 2D:
    mst_k <- dendextend::dist_long(asb_mst[[k]])
    k_branches <- mst_k[which(mst_k$distance == 1), ]
    
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    
    k_branches_xyxy <- data.frame(k_branches,
                                x = asb_sp_coord2D[[k]][k_branches$rows, 1],
                                y = asb_sp_coord2D[[k]][k_branches$rows, 2],
                                xend = asb_sp_coord2D[[k]][k_branches$cols, 1],
                                yend = asb_sp_coord2D[[k]][k_branches$cols, 2])
    
    ggplot_feve <- ggplot_feve +
      ggplot2::geom_segment(data = k_branches_xyxy,
                          ggplot2::aes(x = x , y = y, xend = xend, yend = yend),
                          colour = color_mst[k],
                          size = width_mst[k],
                          linetype = linetype_mst[k])
  }
  
  
  return(ggplot_feve)
  
}# end of function




#' Plot FNND index
#'
#' This function plots FNND index for a given pair of functional axes and for 
#' one or several assemblages
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_nn_asb a list gathering for each species of a studied assemblage 
#' its nearest neighbour(s) in the assemblage. If several assemblages it 
#' should  be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#'  should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#'  If several assemblages it should be a vector with names as in 
#'  `asb_sp_coord2D`.
#'  
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param color_segment a R color name or an hexadecimal code referring to the 
#' color of the segment linking nearest neighbors in the studied assemblage. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param width_segment a numeric value referring to the size of the segment 
#' linking nearest neighbors in the studied assemblage. Default: 
#' `width_segment = 1`. If several assemblages it should be a vector with 
#' names as in `asb_sp_coord2D`.
#' 
#' @param linetype_segment a character string referring to the linetype used to 
#' link nearest neighbors in the studied assemblages. 
#' Default: `linetype_segment = "dashed`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object with FNND index
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' }


fnnd.plot <- function(ggplot_bg,
                    asb_sp_coord2D,
                    asb_sp_relatw,
                    asb_nn_asb,
                    plot_sp = TRUE,
                    shape_sp, color_sp, fill_sp,
                    color_segment, width_segment, linetype_segment) {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_fnnd <- ggplot_bg
  
  # plotting species if required:
  if (plot_sp) {
    ggplot_fnnd <- sp.plot(ggplot_bg = ggplot_fnnd,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = NULL,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp)
  }
  
  # plotting nearest neighbours:
  for (k in asb_nm) {
    
    # nearest neighbours (could be > 1 per species)
    k_sp_xy <- asb_sp_coord2D[[k]]
    k_nn <- asb_nn_asb[[k]]
    
    k_segments_xyxy <- data.frame(sp_f = NULL, sp_nn = NULL, x = NULL, y = NULL,
                                xend = NULL, yend = NULL)
    
    for (i in names(k_nn)) {
      
      i_nn <- k_nn[[i]]
      
      k_segments_xyxy <- rbind(k_segments_xyxy,
                             data.frame(sp_f = rep(i, length(i_nn)),
                                        sp_nn = i_nn,
                                        x = k_sp_xy[i, 1],
                                        y = k_sp_xy[i, 2],
                                        xend = k_sp_xy[i_nn, 1],
                                        yend = k_sp_xy[i_nn, 2]))
    }
    
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    ggplot_fnnd <- ggplot_fnnd +
      ggplot2::geom_segment(data = k_segments_xyxy,
                            ggplot2::aes(x = x, y = y, xend = xend, 
                                         yend = yend),
                            colour = color_segment[k],
                            size = width_segment[k],
                            linetype = linetype_segment[k],
                            arrow = grid::arrow(length = grid::unit(0.10,
                                                                    "inches"),
                                                ends = "last",
                                                type = "open"))
  }

  return(ggplot_fnnd)
  
}




#' Plot Fori
#' 
#' This function plots FOri index for a given pair of functional axes and for 
#' one or several assemblages. It adds the distance of species from the studied
#' to their nearest neighbour(s) from the global pool.
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param asb_nn_pool a list gathering for each species of a studied assemblage 
#' its nearest neighbour(s) in the global pool. If several assemblages it 
#' should  be a vector with names as in `asb_sp_coord2D`.
#' 
#' @param pool_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in the global pool for a given pair of functional axes.
#' 
#' @param plot_pool a logical value indicating whether species of each 
#' assemblage should be plotted or not. Default: `plot_pool = TRUE`.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#' should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param color_pool a R color name or an hexadecimal code referring to the 
#' color of the pool.  This color is also used for FRic convex hull color. 
#' Default: `color_pool = "#0072B2"`.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#' of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param fill_pool a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20) and the assemblage convex
#'  hull. Default: `fill_pool = '#0072B2'`.
#'  
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param shape_pool a numeric value referring to the shape used to plot species
#'  pool. Default: `shape_pool = 16`(filled circle). 
#'  
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#' 
#' @param size_pool a numeric value referring to the size of species belonging 
#' to the global pool. Default: `size_pool = 1`.
#' 
#' @param color_segment a R color name or an hexadecimal code referring to the 
#' color of the segment linking nearest neighbors in the global pool. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param width_segment a numeric value referring to the size of the segment 
#' linking nearest neighbors in the global pool. Default: 
#' `width_segment = 1`. If several assemblages it should be a vector with 
#' names as in `asb_sp_coord2D`.
#' 
#' @param linetype_segment a character string referring to the linetype used to 
#' link nearest neighbors in the global pool. 
#' Default: `linetype_segment = "dashed`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object with FOri index
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' }
#' 

fori.plot <- function(ggplot_bg,
                    asb_sp_coord2D,
                    asb_sp_relatw,
                    asb_nn_pool,
                    pool_coord2D,
                    plot_pool = TRUE,
                    plot_sp = TRUE,
                    shape_pool, size_pool, color_pool, fill_pool,
                    shape_sp, color_sp, fill_sp,
                    color_segment, width_segment, linetype_segment) {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_fori <- ggplot_bg
  
  # plotting species of the pool if required:
  if (plot_pool) {
    ggplot_fori <- pool.plot(ggplot_bg = ggplot_fori,
                            sp_coord2D = pool_coord2D,
                            vertices_nD = NULL,
                            shape_pool = shape_pool, size_pool = size_pool,
                            color_pool = color_pool, fill_pool = fill_pool)
  }
  
  # plotting species if required:
  if (plot_sp) {
    ggplot_fori <- plot_sp(ggplot_bg = ggplot_fori,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = NULL,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp)
  }
  
  
  # plotting nearest neighbour(s) from the pool (could be >1 per species):
  for (k in asb_nm) {
    
    k_nn <- asb_nn_pool[[k]]
    
    k_segments_xyxy <- data.frame(sp_f = NULL, sp_nn = NULL, x = NULL, y= NULL,
                                  xend = NULL,  yend = NULL)
    
    for (i in names(k_nn)) {
      i_nn <- k_nn[[i]]
      k_segments_xyxy <- rbind(k_segments_xyxy,
                             data.frame(sp_f = rep(i, length(i_nn)),
                                        sp_nn = i_nn,
                                        x = pool_coord2D[i, 1],
                                        y = pool_coord2D[i, 2],
                                        xend = pool_coord2D[i_nn, 1],
                                        yend = pool_coord2D[i_nn, 2]))
    }
    
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    ggplot_fori <- ggplot_fori +
      ggplot2::geom_segment(data = k_segments_xyxy,
                            ggplot2::aes(x = x , y = y, xend = xend, 
                                         yend = yend),
                            colour = color_segment[k],
                            size = width_segment[k],
                            linetype = linetype_segment[k],
                            arrow = grid::arrow(length = grid::unit(0.10,
                                                                    "inches"),
                                                ends = "last",
                                                type = "open"))
  }
  
  return(ggplot_fori)
  
}




#' Plot FSpe
#' 
#' This function plots FSpe index for a given pair of functional axes and for
#' one or several assemblages. It adds the mean position of species from the 
#' global pool and the distance of each species from the studied assemblage(s) 
#' on the background plot
#' 
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage for a given pair of functional axes.
#' 
#' @param asb_sp_relatw a list of vector gathering species relative weight in
#'   each assemblage. It can be retrieved through the
#'   \code{\link{alpha.fd.multidim}}. If several assemblages it should be a 
#'   vector with names as in `asb_sp_coord2D`.
#' 
#' @param center_coord2D a list containing the coordinates of the center of the
#' global pool for two given functional axes 
#' 
#' @param pool_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in the global pool for a given pair of functional axes.
#' 
#' @param plot_pool a logical value indicating whether species of each 
#' assemblage should be plotted or not. Default: `plot_pool = TRUE`.
#' 
#' @param plot_sp a logical value indicating whether species of each assemblage 
#' should be plotted or not. Default: `plot_sp = TRUE`
#' 
#' @param color_pool a R color name or an hexadecimal code referring to the 
#' color of the pool.  This color is also used for FRic convex hull color. 
#' Default: `color_pool = "#0072B2"`.
#' 
#' @param color_sp a R color name or an hexadecimal code referring to the color
#' of species from the studied assemblage. Default: `color_sp = "#0072B2"`. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param color_center a R color name or an hexadecimal code referring to the 
#' color of the center of the global pool. Default: `color_center = "#0072B2"`. 
#' 
#' @param color_segment a R color name or an hexadecimal code referring to the 
#' color of the segment linking nearest neighbors in the global pool. 
#' If several assemblages it should be a vector with names as in 
#' `asb_sp_coord2D`.
#' 
#' @param fill_pool a R color name or an hexadecimal code referring to the 
#' colour to fill species symbol (if \code{shape_sp} > 20) and the assemblage 
#' convex hull. Default: `fill_pool = '#0072B2'`.
#'  
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20). Default: 
#'  `fill_sp = '#0072B2'`. If several assemblages it should  be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param fill_center a R color name or an hexadecimal code referring to the 
#' colour to fill the center of the global pool (if \code{shape_sp} > 20). 
#' Default: `fill_center = '#0072B2'`.
#' 
#' @param shape_pool a numeric value referring to the shape used to plot species
#'  pool. Default: `shape_pool = 16`(filled circle). 
#'  
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: `shape_sp = 16` 
#'  (filled circle). If several assemblages it should be a vector with 
#'  names as in `asb_sp_coord2D`.
#'  
#' @param shape_center a numeric value referring to the shape used to plot the 
#' center of the global pool.Default: `shape_center = 16`(filled circle).
#' 
#' @param size_pool a numeric value referring to the size of species belonging 
#' to the global pool. Default: `size_pool = 1`.
#' 
#' @param size_center a numeric value referring to the size of the center of the
#' global pool. Default: `size_center = 1`.
#' 
#' @param width_segment a numeric value referring to the size of the segment 
#' linking nearest neighbors in the global pool. Default: 
#' `width_segment = 1`. If several assemblages it should be a vector with 
#' names as in `asb_sp_coord2D`.
#' 
#' @param linetype_segment a character string referring to the linetype used to 
#' link nearest neighbors in the global pool. 
#' Default: `linetype_segment = "dashed`.  If several assemblages it should be a 
#' vector with names as in `asb_sp_coord2D`.
#' 
#' @return a ggplot object with FSpe index on the background plot
#' 
#' @author Camille Magneville and Sébastien Villéger
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' }


fspe.plot <- function(ggplot_bg,
                    asb_sp_coord2D,
                    asb_sp_relatw,
                    center_coord2D,
                    pool_coord2D,
                    plot_pool = TRUE,
                    plot_sp = TRUE,
                    shape_pool, size_pool, color_pool, fill_pool,
                    shape_sp, color_sp, fill_sp,
                    color_center, fill_center,
                    shape_center, size_center,
                    color_segment, width_segment, linetype_segment) {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  ## plotting layers ####
  
  # default plot is background:
  ggplot_fspe <- ggplot_bg
  
  # plotting species of the pool if required:
  if (plot_pool) {
    ggplot_fspe <- pool.plot(ggplot_bg = ggplot_fspe,
                            sp_coord2D = pool_coord2D,
                            vertices_nD = NULL,
                            shape_pool = shape_pool, size_pool = size_pool,
                            color_pool = color_pool, fill_pool = fill_pool)
  }
  
  # plotting species if required:
  if(plot_sp) {
    ggplot_fspe <- sp.plot(ggplot_bg = ggplot_fspe,
                           asb_sp_coord2D = asb_sp_coord2D,
                           asb_vertices_nD = NULL,
                           asb_sp_relatw = asb_sp_relatw,
                           shape_sp = shape_sp,
                           color_sp = color_sp,
                           fill_sp = fill_sp)
  }
  
  # coordinates of center of functional space:
  center_xy <- data.frame(x = center_coord2D[1], y = center_coord2D[2],
                         row.names = NULL )
  
  
  # plotting distance to center of space:
  for (k in asb_nm) {
    
    # coordinates of segments to center of functional space:
    k_sptocenter_xyxy <- data.frame(x = asb_sp_coord2D[[k]][, 1],
                                  y = asb_sp_coord2D[[k]][, 2],
                                  xend = center_xy[1, 1],
                                  yend = center_xy[1, 2])
    # plotting segments:
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    ggplot_fspe <- ggplot_fspe +
      ggplot2::geom_segment(data = k_sptocenter_xyxy,
                            ggplot2::aes(x = x , y = y, xend = xend, 
                                          yend = yend),
                            colour = color_segment[k],
                            size = width_segment[k],
                            linetype = linetype_segment[k])
  }
  
  # plotting center of space:
  ggplot_fspe <- ggplot_fspe +
    ggplot2::geom_point(data = center_xy, ggplot2::aes(x = x , y = y),
                        colour = color_center, fill = fill_center,
                        shape = shape_center, size = size_center)
  
  return(ggplot_fspe)
  
}
