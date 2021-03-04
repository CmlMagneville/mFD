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
  
  x1 <- NULL
  y1 <- NULL
  ggplot_bg <- ggplot2::ggplot(data.frame(x1 = range_faxes, y1 = range_faxes),
                             ggplot2::aes(x = x1, y = y1)) +
    ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
    ggplot2::xlab(faxes_nm[1]) +
    ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0,0)) +
    ggplot2::ylab(faxes_nm[2]) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = color_bg)) +
    ggplot2::coord_fixed()
  
  return(ggplot_bg)
  
} # end of plot_background






#' Plot FRic indice
#'
#' This function adds convex hull(s), points and vertices of 1:N assemblages on
#' the backrgound plot
#'
#' @param ggplot_bg a ggplot object of the plot background retrieved through
#' the \code{\link{background.plot}} function.
#' 
#' @param asb_sp_coord2D a list of matrix (ncol = 2) with coordinates of 
#' species present in each assemblage
#' 
#' @param vertices_nD a list (with names as in asb_sp_coord2D) of vectors with 
#' names of species being vertices in n dimensions.
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
#'  the global pool but not the plotted assemblage. Default: `size_sp = 1`. If 
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
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
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
#'  vert_nm_fruits <- vertices(sp_faxes_coord_fruits, 
#'   order_2D = TRUE, check_input = TRUE)
#'                                
#'  # Plot in white the convex hull of all fruits species:
#'  ggplot_fric <- mFD::fric.plot(
#'            ggplot_bg      = ggplot_bg_fruits,
#'            asb_sp_coord2D = list(sp_faxes_coord_fruits),
#'            vertices_nD    = vert_nm_fruits,
#'              color_ch     = "black", 
#'              fill_ch      = "white", 
#'              alpha_ch     = 0.3,
#'              shape_sp     = 16, 
#'              size_sp      = 1, 
#'              color_sp     = "#0072B2", 
#'              fill_sp      = "#0072B2",
#'              shape_vert   = 15, 
#'              size_vert    = 1, 
#'              color_vert   = "#0072B2", 
#'              fill_vert    = "#0072B2")
#'  ggplot_fric
#'
#' }
#'
#'

fric.plot <- function(ggplot_bg,
                    asb_sp_coord2D,
                    vertices_nD,
                    color_ch, fill_ch, alpha_ch,
                    shape_sp, size_sp, color_sp, fill_sp,
                    shape_vert, size_vert, color_vert, fill_vert)  {
  
  # names of assemblages:
  asb_nm <- names(asb_sp_coord2D)
  
  # default plot is background:
  ggplot_fric <- ggplot_bg
  
  
  # computing and plotting 2D convex hull(s) ----
  
  for (k in asb_nm) {
    
    # vertices of the convex hull in 2D:
    vert2D_k <- vertices(asb_sp_coord2D[[k]], order_2D = TRUE)
    
    # plotting convex hulls:
    ggplot_fric <- ggplot_fric +
      ggplot2::geom_polygon(data = asb_sp_coord2D[[k]][vert2D_k, ],
                            colour = color_ch[k],
                            fill = fill_ch[k],
                            alpha = alpha_ch[k])
  }# end of k
  
  
  # plotting species not vertices in nD ----
  for (k in asb_nm) {
    
    sp_k <- rownames(asb_sp_coord2D[[k]])
    not_vertices_nD_k <- sp_k[which(! (sp_k %in% vertices_nD[[k]]))]
    ggplot_fric <- ggplot_fric +
      ggplot2::geom_point(data = asb_sp_coord2D[[k]][not_vertices_nD_k, ],
                          colour = color_sp[k], fill = fill_sp[k],
                          shape = shape_sp[k], size = size_sp[k])
  }# end of k
  
  
  # plotting species vertices in nD ----
  for (k in asb_nm) {    
    ggplot_fric <- ggplot_fric +
    ggplot2::geom_point(data = asb_sp_coord2D[[k]][vertices_nD[[k]], ],
                        colour = color_vert[k], fill = fill_vert[k],
                        shape = shape_vert[k], size = size_vert[k])
  }# end of k
  
  return(ggplot_fric)
  
}# end of function



