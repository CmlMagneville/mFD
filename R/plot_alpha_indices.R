#' Plot functional space and chosen functional indices
#'
#' Compute a graphical representation of functional indices. \strong{To plot
#' functional indices, functional indices values must have been retrieve through
#' the use of the} \code{\link{alpha.fd.multidim}} \strong{function}.
#'
#' @param sp_faxes_coord a matrix of species coordinates in a chosen
#'   functional space. Species coordinates have been retrieved thanks to
#'   \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param asb_sp_w a matrix linking weight of species (columns) and a set of
#'  assemblages (rows)
#'
#' @param ind_vect a vector of character string of the name of functional
#'   indices to plot. \strong{Indices names must be written in lower case
#'   letters}. Possible indices to compute are: "fdis", "feve", "fric", "fdiv",
#'   "fori" and "fspe". Default: all the indices are computed.
#'
#' @param details_list a list retrieved with the
#'  \code{\link{alpha.fd.multidim}} function containing information about
#'  features used to compute indices.
#'
#' @param asb_vect a vector containing the name of the assemblages to represent.
#'  Two assemblages can be represented in each plot if asked by the user.
#'
#' @param fd_ind_values a dataframe containing the numeric value of each
#'   computed indices. This dataframe is an output of the
#'   \code{\link{alpha.fd.multidim}} function.
#'
#' @param faxes a vector with names of axes to plot. \strong{You can only plot
#'  from 2 to 4 axes for graphical reasons: vector length should be between 2
#'  and 4}. Default: faxes = NULL (the four first axes will be plotted).
#'
#' @param faxes_nm a vector with axes labels if the user want different axes
#'  labels than \code{faxes} ones. Default: faxes_nm = faxes (labels will the
#'  the same that \code{faxes} ones).
#'
#' @param range_faxes_lim a vector with minimum and maximum for values for axes.
#'  Note that to have a fair representation of position of species in all plots,
#'  all axes must have the same range. Default: faxes_lim = c(NA, NA) (the range
#'  is computed according to the range of values among all axes, all axes having
#'  the same range).
#'
#' @param color_bg a R color name  or an hexadecimal code used to fill plot
#'  background. Default: color_bg = "grey95".
#'
#' @param size_sp a numeric value referring to the size of species belonging to
#'  the global pool but not the plotted assemblage. Default: size_sp = 1.
#'
#' @param size_centroid a numeric value referring to the size of the centroid
#'  point. Used for FDiv, FSpe, FDis plotting. Default: size_centroid = 1.
#'
#' @param size_centroid_asb2 a numeric value referring to the size of the
#'  centroid point for the second assemblage to plot if there is one. Used for
#'  FDiv, FSpe, FDis plotting. Default: size_centroid = 1.
#'
#' @param size_vert a numeric value referring to the size of symbol for vertices
#'  Default: size_vert = 1.
#'
#' @param size_sp_nm a numeric value for size of species label. Default:
#'  size_sp_nm = 3.
#'
#' @param color_sp a R color name or an hexadecimal code referring to the color
#'  of species from the studied assemblage.  This color is also used for FRic
#'  convex hull color. Default: color_sp = "#0072B2".
#'
#' @param color_sp_asb2 a R color name or an hexadecimal code referring to the
#'  colour of species from the second assemblage to plot if there is one. This
#'  color is also used for FRic convex hull color. Default: color_sp_asb2 =
#'  "#D55E00".
#'
#' @param color_sp_gp a R color name or an hexadecimal code to give to species
#'  that do not belong to the studied assemblages. Default: color_sp_global_pool
#'  = "gray80".
#'
#' @param color_segment a R color name or an hexadecimal code referring to the
#'  color of segments linking species of the studied assemblage and/or centroid
#'  for FDiv, FEve, FSpe, FOri, FNND, FDis computation. Defaut is the same than
#'  \code{color_sp} ; Default: color_segment = color_segment = "#0072B2".
#'
#' @param color_segment_asb2 a R color name or an hexadecimal code referring to
#'  the color of segments linking species and/or centroid of the second
#'  assemblage to plot if there is one for FDiv, FEve, FSpe, FOri, FNND, FDis
#'  computation. Defaut is the same than \code{color_sp_asb2} ; Default:
#'  color_segment = color_segment = "#0072B2".
#'
#' @param color_centroid a R color name or an hexadecimal code referring to the
#'  color of the centroid point. Used for FDiv, FSpe, FDis plotting. Default:
#'  color_centroid = 1. Defaut is the same than \code{color_sp} ; Default:
#'  color_segment = color_segment = "#0072B2".
#'
#' @param color_centroid_asb2 a R color name or an hexadecimal code referring to
#'  the color of the centroid point for the second assemblage to plot if there
#'  is one. Used for FDiv, FSpe, FDis plotting. Default: color_centroid = 1.
#'  Defaut is the same than \code{color_sp_asb2} ; Default: color_segment =
#'  color_segment = "#0072B2".
#'
#' @param color_vert a R color name or an hexadecimal code referring to the
#'   color of vertices if plotted. If color_vert = NA, vertices are not plotted
#'   (for shapes only defined by color, ie shape inferior to 20. Otherwise fill
#'   must also be set to NA). Default: color_vert =  NA.
#'
#' @param color_vert_asb2 a R color name or an hexadecimal code referring to the
#'  color of vertices if plotted for the second assemblage. If color_vert = NA,
#'  vertices are not plotted (for shapes only defined by color, ie shape < 20.
#'  Otherwise fill must also be set to NA). Default: color_vert_asb2 =  NA.
#'
#' @param color_ch a R color name or an hexadecimal code referring to the border
#'  of the convex hull filled by the pool of species. Default: color_ch =
#'  "black".
#'
#' @param color_sp_nm a R color name or an hexadecimal code referring to the
#'  colour of species label. Default: color_sp_nm = "black".
#'
#' @param fill_sp a R color name or an hexadecimal code referring to the colour
#'  to fill species symbol (if \code{shape_sp} > 20) and the assemblage convex
#'  hull. Default: fill_sp = '#0072B2'
#'
#' @param fill_sp_gp a R color name or an hexadecimal code referring to the
#'  colour to fill symbol for species from the global pool (if \code{shape_sp}
#'   superior to 20). Default: fill_sp_gp = "grey80".
#'
#' @param fill_sp_asb2 a R color name or an hexadecimal code referring to the
#'   colour to fill species symbol for the second assemblage (if \code{shape_sp}
#'   superior to 20) and its associated convex hull. Default: fill_sp_asb2 =
#'   'white'.
#'
#' @param fill_vert a character value referring to the color for filling symbol
#'   for vertices (if \code{shape_vert} >20). If fill = NA and color = NA,
#'   vertices are not plotted (if \code{shape_vert} superior to 20. Otherwise
#'   color_vert = NULL is enough). Default is NA.
#'
#' @param fill_vert_asb2 a character value referring to the color for filling
#'   symbol for vertices (if \code{shape_vert} superior to 20) for the second
#'   assemblage. If fill = NA and color = NA, vertices are not plotted (if
#'   \code{shape_vert} superior to 20. Otherwise color_vert = NA is enough).
#'   Default is NA.
#'
#' @param fill_ch a R color name or an hexadecimal code referring to the filling
#'  of the convex hull filled by the pool of species. Default is "white".
#'
#' @param fill_centroid a R color name or an hexadecimal code referring to the
#'  filling of the centroid (if \code{shape_centroid} superior to 20). Default:
#'  fill_centroid = '#0072B2'.
#'
#' @param fill_centroid_asb2 a R color name or an hexadecimal code referring to
#'  the filling of the centroid of the second assemblage (if
#'  \code{shape_centroid} superior to 20). Default: fill_centroid = "#D55E00".
#'
#' @param alpha_ch a numeric value for transparency of the filling of the convex
#'  hull (0 = high transparency, 1 = no transparency). Default is 0.3.
#'
#' @param shape_sp_gp a numeric value referring to the shape used to plot
#'   species belonging to the global pool but not to the studied assemblage(s).
#'   Default: shape_sp_global_pool = 3 (horizontal cross).
#'
#' @param shape_sp a numeric value referring to the shape used to plot species
#'  belonging to the studied assemblage. Default: shape_sp = 16 (filled circle).
#'
#' @param shape_sp_asb2 a numeric value referring to the shape used to plot
#'  species belonging to the second assemblage to plot if there is one. Default:
#'  shape_sp_asb2 = 15 (filled square)
#'
#' @param shape_vert a numeric value referring to the shape used to plot
#'   vertices if vertices should be plotted in a different way than other
#'   species. If shape_vert = NA, no vertices plotted. Default: shape_vert = NA.
#'
#' @param shape_vert_asb2 a numeric value referring to the shape used to plot
#'  vertices of the second assemblage if vertices should be plotted in a
#'  different way than other species. If shape_vert = NA, no vertices plotted.
#'  Default: shape_vert_asb2 = NA.
#'
#' @param shape_centroid a numeric value referring to the shape used to plot
#'  centroid for the given assemblage. Default: shape_centroid = 10 (horizontal
#'  cross inside an empty circle).
#'
#' @param shape_centroid_asb2 a numeric value referring to the shape used to
#'   plot centroid for the second assemblage to plot if there is one. Default:
#'   shape_centroid_asb2 = 12 (horizontal cross inside an empty square).
#'
#' @param shape_vert a numeric value referring to the symbol used to show
#'  vertices position if \code{plot_vertices} = TRUE. Default is 23 (filled
#'  diamond).
#'
#' @param segment_size a numeric value referring to the size of the segment used
#'  to link species of a given assemblage and centroid. Default: segment_size =
#'  1
#'
#' @param segment_size_asb2 a numeric value referring to the size of the segment
#'  used to link species of the second assemblage and centroid. Default:
#'  segment_size = 0.5
#'
#' @param linetype_segment a character string or the associated numeric value
#'  referring the type of line used for segments linking species of the studied
#'  assemblage and/or centroid for FDiv, FEve, FSpe, FOri, FNND, FDis
#'  computation. Default: linetype_segment = "solid".
#'
#' @param linetype_segment_asb2 a character string or the associated numeric
#'  value referring the type of line used for segments linking species and/or
#'  centroid of the second assemblage to plot if there is one. Used for FDiv,
#'  FEve, FSpe, FOri, FNND, FDis computation. Default: linetype_segment =
#'  "dashed".
#'
#' @param scale_inf a numeric value referring to the minimal size of a point in
#'   a plot according to each species's relative weight. Default: scale_inf = 1.
#'
#' @param scale_sup a numeric value referring to the minimal size of a point in
#'   a plot according to each species's relative weight. Default: scale_inf = 3.
#'
#' @param plot_sp_nm a vector containing species names that are to be plotted.
#'  Default: plot_nm_sp = NULL (no name plotted).
#'
#' @param plot_ch a logical value indicating whether the convex hull shaping the
#'  pool of species should be illustrated. If plot_ch = TRUE, convex hull of all
#'  species in the multidimensional space described in \code{sp_faxes_coord} is
#'  computed and its projection in 2D spaces are drawn as polygons. Default:
#'  plot_ch = TRUE.
#'
#' @param fontface_nm a character string for font of species labels (e.g.
#'  "italic", "bold"). Default: fontface_nm = 'plain'.
#'
#' @param name_file a character string with name of file to save the figure
#'  (without extension). Default is 'NULL' which means plot is displayed.
#'  If several plots are to be saved (for several indices), then files are
#'  named as follow "name_file1", "name_file2"...
#'
#' @param check_input a logical value defining whether inputs are checked before
#'  computation of indices. Possible error messages will thus may be more
#'  understandable for the user than R error messages. Default: check_input =
#'  TRUE.
#'
#' @return for the given assemblage, return a list of one \code{patchwork}
#'   figure per functional indice containing plots for combinations of up to
#'   four axes.
#'
#' @examples
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
#' # Retrieve a list of details needed for plots:
#' details_list_fruits <- alpha_fd_indices_fruits$details
#' 
#' # Compute the plots:
#' mFD::alpha.multidim.plot(sp_faxes_coord = sp_faxes_coord_fruits, 
#'  asb_sp_w = baskets_fruits_weights,
#'  ind_vect = c("fdis",
#'              "fnnd", "feve",
#'              "fric", "fdiv",
#'              "fori", "fspe"),
#'  details_list = details_list_fruits,
#'  asb_vect = c("basket_1"),
#'  fd_ind_values = fd_ind_values_fruits,
#'  faxes = NULL, faxes_nm = NULL,
#'  range_faxes_lim = c(NA, NA),
#'  color_bg = "grey95",
#'  size_sp = 1,
#'  size_centroid = 1,
#'  size_centroid_asb2 = 1,
#'  size_vert = 1,
#'  size_sp_nm = 3,
#'  color_sp = "#0072B2", color_sp_asb2 = "#D55E00",
#'  color_sp_gp = "gray80",
#'  color_segment = "#0072B2",
#'  color_segment_asb2 = "#CC79A7",
#'  color_centroid = '#0072B2',
#'  color_centroid_asb2 = "#D55E00",
#'  color_vert = NA, color_vert_asb2 = NA,
#'  color_ch = "black",
#'  color_sp_nm = "black",
#'  fill_sp = "white",
#'  fill_sp_asb2 = "white",
#'  fill_sp_gp = "gray80",
#'  fill_vert = NA,
#'  fill_vert_asb2 = NA,
#'  fill_ch = "white",
#'  fill_centroid = '#0072B2',
#'  fill_centroid_asb2 = "#D55E00",
#'  alpha_ch = 0.3,
#'  shape_sp_gp = 3,
#'  shape_sp = 16,
#'  shape_sp_asb2 = 15,
#'  shape_vert = NA,
#'  shape_vert_asb2 = NA,
#'  shape_centroid = 10,
#'  shape_centroid_asb2 = 12,
#'  segment_size = 1,
#'  segment_size_asb2 = 0.5,
#'  linetype_segment = "solid",
#'  linetype_segment_asb2 = "dashed",
#'  scale_inf = 1, scale_sup = 3,
#'  plot_sp_nm = NULL,
#'  plot_ch = TRUE,
#'  fontface_nm = "plain",
#'  name_file = NULL,
#'  check_input = TRUE)
#'  
#' @author Camille Magneville and Sébastien Villéger
#'
#' @importFrom ggplot2 aes_ geom_hline geom_vline geom_segment geom_rect
#' @importFrom ggplot2 geom_point geom_polygon scale_size scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous theme theme_void ggplot ggsave
#' @importFrom grid arrow unit
#' @importFrom utils tail
#' @importFrom rlist list.append
#' @importFrom patchwork plot_layout plot_annotation plot_spacer
#'
#' @export


alpha.multidim.plot <- function(sp_faxes_coord, asb_sp_w,
                                ind_vect = c("fdis",
                                             "fnnd", "feve",
                                             "fric", "fdiv",
                                             "fori", "fspe"),
                                details_list,
                                asb_vect,
                                fd_ind_values,
                                faxes = NULL, faxes_nm = NULL,
                                range_faxes_lim = c(NA, NA),
                                color_bg = "grey95",
                                size_sp = 1,
                                size_centroid = 1,
                                size_centroid_asb2 = 1,
                                size_vert = 1,
                                size_sp_nm = 3,
                                color_sp = "#0072B2", color_sp_asb2 = "#D55E00",
                                color_sp_gp = "gray80",
                                color_segment = "#0072B2",
                                color_segment_asb2 = "#CC79A7",
                                color_centroid = '#0072B2',
                                color_centroid_asb2 = "#D55E00",
                                color_vert = NA, color_vert_asb2 = NA,
                                color_ch = "black",
                                color_sp_nm = "black",
                                fill_sp = "white",
                                fill_sp_asb2 = "white",
                                fill_sp_gp = "gray80",
                                fill_vert = NA,
                                fill_vert_asb2 = NA,
                                fill_ch = "white",
                                fill_centroid = '#0072B2',
                                fill_centroid_asb2 = "#D55E00",
                                alpha_ch = 0.3,
                                shape_sp_gp = 3,
                                shape_sp = 16,
                                shape_sp_asb2 = 15,
                                shape_vert = NA,
                                shape_vert_asb2 = NA,
                                shape_centroid = 10,
                                shape_centroid_asb2 = 12,
                                segment_size = 1,
                                segment_size_asb2 = 0.5,
                                linetype_segment = "solid",
                                linetype_segment_asb2 = "dashed",
                                scale_inf = 1, scale_sup = 3,
                                plot_sp_nm = NULL,
                                plot_ch = TRUE,
                                fontface_nm = "plain",
                                name_file = NULL,
                                check_input = TRUE) {
  
  
  #### Retrieve all the needed elements from details_list ####
  
  asb_sp_relatw                <- details_list$asb_sp_relatw
  mst_list                    <- details_list$mst_list
  grav_center_vert_coord_list <- details_list$grav_center_vert_coord_list
  grav_center_global_pool        <- details_list$grav_center_global_pool
  nm_nn_global_pool_list              <- details_list$nm_nn_global_pool_list
  nm_nn_asb_list                  <- details_list$nm_nn_asb_list
  
  
  # Check inputs relative to this function (funct.space.plot, already done) ####
  
  if (check_input == TRUE) {
    
    # check that functional indices have the right names:
    for (ind in ind_vect) {
      if (! ind %in% c("fdis","fnnd", "feve", "fric", "fdiv", "fori", "fspe")) {
        stop("Error: Provided names of functional indices are not well written.
             Please re-write them. Be careful, they should all be written in 
             lowercase letters.")
      }
    }
    
    # check that good number of assemblage(s) to plot:
    if (! length(asb_vect) %in% c(1, 2)) {
      stop("Error: This function can only plot one or two assemblages.
           Please chose two or less assemblages to plot.")
    }
    
    # check that assemblage(s) to plot has(ve) the right name(s):
    for (asb in asb_vect) {
      if (! asb %in% colnames(asb_sp_relatw)) {
        stop("Error: Provided name(s) of assemblage(s) to plot is(are) not well 
        written.
             Please re-write.")
      }
    }
    
    if (any(is.na(sp_faxes_coord))) {
      stop("Error: The species*coordinates matrix contains NA. Please check.")
    }
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    
    if (is.matrix(asb_sp_w) == FALSE) {
      stop("Error: 'asb_sp_w' must be a matrix")
    }
    
    if (any(is.na(asb_sp_w))) {
      stop("Error: The species*weights matrix contains NA. Please check.")
    }
    if (is.null(rownames(asb_sp_w))) {
      stop("Error: No row names provided in species*weights matrix.
             Please add assemblages names as row names.")
    }
    if (is.null(colnames(asb_sp_w))) {
      stop("Error: No column names provided in species*weight matrix.
             Please add species names as column names.")
    }
    
    # Add a stop if there is a negative value in the occurrence dataframe:
    if (any(asb_sp_w < 0)) {
      stop("Error: The species*weight matrix should not contain negative values.
           Please check.")
    }
    
    asb_sp_w <- as.data.frame(asb_sp_w)
    isnum_vect <- sapply(asb_sp_w, is.numeric)
    
    if (FALSE %in% isnum_vect) {
      stop("Error: The 'asp_sp_w' matrix must only contain numeric values. 
           Please convert values")
    }
    
    if (any(!(colnames(asb_sp_w) %in% rownames(sp_faxes_coord)))) {
      stop(paste("Error: Mismatch between names in species*weight and
                   species*coordinates matrix. Please check."))
    }
    
    # Add a stop if some species do not belong to any assemblage:
    if (min(apply(asb_sp_w, 2, sum)) == 0){
      stop("Error: Some species are absent from all assemblages.")
    }
    # Add a stop if some asb do not contain species:
    if (min(apply(asb_sp_w, 1, sum)) == 0){
      stop("Error: Some assemblages do not contain species.")
    }
  }
  
  
  #### Diverse data manipulation ####
  
  # change sp_faxes_coord format to use colnames:
  sp_faxes_coord <- as.data.frame(sp_faxes_coord)
  
  # give faxes identity if faxes set to NULL:
  if (is.null(faxes)) {
    faxes <- colnames(sp_faxes_coord)[1:min(c(4, ncol(sp_faxes_coord)))]
  }
  
  # give faxes names if faxes set to NULL:
  if (is.null(faxes_nm)) {
    faxes_nm <- faxes
  }
  
  # give a value to plot_vertices argument of funct.space.plot:
  if (is.null(color_vert)) {
    plot_vertices <- FALSE
  } else {
    plot_vertices <- TRUE
  }
  
  
  #### Plot functional space ####
  
  # change sp_faxes_coord format to use in funct.space.plot:
  sp_faxes_coord <- as.matrix(sp_faxes_coord)
  
  # plot functional space using function plot.funct.space:
  funct_space_output <- funct.space.plot(sp_faxes_coord, faxes = faxes, 
                                         name_file = name_file,
                                         faxes_nm = faxes_nm, 
                                         range_faxes_lim = c(NA, NA),
                                         color_bg = color_bg,
                                         color_sp = NA, fill_sp = NA,  
                                         shape_sp = shape_sp_gp, 
                                         size_sp = size_sp,
                                         plot_ch = plot_ch,  
                                         color_ch = color_ch, 
                                         fill_ch = fill_ch, 
                                         alpha_ch = alpha_ch,
                                         plot_vertices = FALSE, 
                                         color_vert = color_vert, 
                                         fill_vert = fill_vert,
                                         shape_vert = shape_vert,
                                         size_vert = size_vert,
                                         plot_sp_nm = plot_sp_nm, 
                                         nm_size = size_sp_nm, 
                                         nm_color = color_sp_nm, 
                                         nm_fontface = fontface_nm,
                                         check_input = check_input)
  
  
  #### Retrieve information, create df and format data for the 1st asb ####
  
  # retrieve the number of plots done of axes combination...
  # ... one plot is the summary of all the combination of plots, does not count:
  plot_nb   <- length(funct_space_output)
  
  # retrieve assemblage name (1st one if 2 to plot):
  asb_k <- asb_vect[1]
  
  # retrieve coordinates of species belonging to asb_k:
  sp_filter <- sp.filter(asb_k, sp_faxes_coord, asb_sp_w)
  sp_faxes_coord_k <- sp_filter$`species coordinates`
  sp_faxes_coord_k <- data.matrix(sp_faxes_coord_k)
  asb_sp_relatw_k <- sp_filter$`species relative weight`
  asb_sp_relatw_k <- data.matrix(asb_sp_relatw_k)
  
  # save species coord of asb 1 (in case another one if asked to be plotted):
  sp_coord_asb_1 <- sp_faxes_coord_k
  
  # format "asb_sp_relatw_k" so easy to plot:
  sp_faxes_coord_k2 <- sp_faxes_coord_k
  sp_faxes_coord_k2 <- as.data.frame(sp_faxes_coord_k2)
  weight <- c()
  for (n in (1:ncol(asb_sp_relatw_k))) {
    weight <- append(weight, asb_sp_relatw_k[1, n])
  }
  sp_faxes_coord_k2$w <- weight
  sp_faxes_coord_k2 <- as.data.frame(sp_faxes_coord_k2)
  
  # keep species weight of asb 1 for plotting:
  sp_coord_asb_w_1 <- sp_faxes_coord_k2$w
  
  # change again sp_faxes_coord format to use it easily:
  sp_faxes_coord <- as.data.frame(sp_faxes_coord)
  
  # create a table that will sum up coordinates of all species to plot...
  # ... (global pool) and add a column to know if there are...
  # ... present on the assemblage to plot:
  sp_coord2 <- sp_faxes_coord
  sp_coord2$asb <- rownames(sp_coord2)
  sp_coord2$asb[which(! sp_coord2$asb %in% rownames(sp_faxes_coord_k))] <- "no"
  sp_coord2$asb[which(sp_coord2$asb %in% rownames(sp_faxes_coord_k))] <- asb_k
  sp_coord2$asb <- as.factor(sp_coord2$asb)
  # count the number of column of this dataframe used to retrieve last ...
  # ... column when plotting:
  col_nb <- ncol(sp_coord2)
  
  # set range of axes if c(NA, NA):
  if (is.na(range_faxes_lim[1]) & is.na(range_faxes_lim[2])) {
    range_sp_coord <- range(sp_faxes_coord)
    range_faxes_lim <- range_sp_coord +
      c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
  }
  
  # create a return list for each index:
  return_fric_list <- list()
  return_fdiv_list <- list()
  return_feve_list <- list()
  return_fspe_list <- list()
  return_fdis_list <- list()
  return_fori_list <- list()
  return_fnnd_list <- list()
  
  
  # Do the plot for each index and along asked dimensions for the 1st asb ####
  
  
  # span across the plots (up to 6 plots for 4 dimensions):
  for (i in (1:(plot_nb - 1))) {
    
    ## retrieve information on the graph to plot:
    
    # retrieve the plot to work on:
    assign(paste0("plot_funct", sep = "_", i), funct_space_output[[i]])
    
    # retrieve the combination of dimensions of the plot:
    faxes_nm <- c(eval(parse(text = paste0("plot_funct", sep = "_", i,
                                           "$labels$x"))), 
                  eval(parse(text = (paste0("plot_funct", sep = "_", i, 
                                            "$labels$y")))))
    
    # retrieve vertices names of asb_k along the 2 dimensions to be plotted ...
    # ... before all indices and not just for FRic because needed if 
    # ... vertices number  is to show in graphs caption for all indices:
    vert_nm_asb_k <- vertices(sp_faxes_coord_k[, c(faxes_nm[[1]], 
                                                   faxes_nm[[2]])], 
                              check_input = TRUE)
    
    # retrieve coordinates of species of asb_k along the 2 plotted dimensions:
    vert_sp_faxes_coord_k <- sp_faxes_coord_k[which(rownames(
                                        sp_faxes_coord_k) %in% vert_nm_asb_k),
                                              c(faxes_nm[[1]], faxes_nm[[2]])]
    
    # then vertices must be ordered so that the convex hull...
    #... can be printed (outside path): sort them clockwise:
    # (https://stackoverflow.com/questions/48249540/plot-convex-hull-given-
    # ... by-quickhull-algorithm-in-r-convhulln-function)
    # find the gravity center of present species and then compute the ...
    # ... angle value of each species to this points, then order angle values:
    vert_sp_faxes_coord_k <- vert_sp_faxes_coord_k[order(-1 * atan2(
      vert_sp_faxes_coord_k[, 
                            faxes_nm[[2]]] - mean(range(vert_sp_faxes_coord_k[, 
                                                          faxes_nm[[2]]])),
      vert_sp_faxes_coord_k[, 
                            faxes_nm[[1]]] - mean(range(vert_sp_faxes_coord_k[,
                                                          faxes_nm[[1]]])))), ]
    
    # convert the format so that it can be used with ggplot2:
    vert_sp_faxes_coord_k <- as.data.frame(vert_sp_faxes_coord_k)
    sp_faxes_coord_k <- as.data.frame(sp_faxes_coord_k)
    
    # save the vert_sp_coord for caption (asb1)
    vert_sp_coord_asb1 <- vert_sp_faxes_coord_k
    
    ## plot each index:
    
    # FRic index:
    
    if ("fric" %in% ind_vect) {
      
      # check that FRic can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fric"])) {
        stop(paste0("Error: FRic value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'), 
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            ggplot2::aes_(x = sp_coord2[which(
                                            sp_coord2[, col_nb] == 'no'),
                                                        c(faxes_nm[[1]])],
                                          y = sp_coord2[which(
                                            sp_coord2[, col_nb] == 'no'),
                                                        c(faxes_nm[[2]])]),
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_polygon(data = vert_sp_faxes_coord_k,
                              ggplot2::aes_(x = vert_sp_faxes_coord_k[, 
                                                              faxes_nm[[1]]],
                                            y = vert_sp_faxes_coord_k[, 
                                                              faxes_nm[[2]]]),
                              fill = color_sp, alpha = alpha_ch) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k, 
                          ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                        y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                          colour = color_sp, shape = shape_sp, fill = fill_sp) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                          ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                          y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                          shape = shape_vert, size = size_vert) +
        ggplot2::theme(legend.position = "none")
      
      return_fric_list[[i]] <- plot_k
      
    }
    
    # FDiv index:
    
    if ("fdiv" %in% ind_vect) {
      
      # check that FDiv can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fdiv"])) {
        stop(paste0("Error: FDiv value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve gravity center of vertices for the studied assemblage:
      grav_center_vert_coord_list <- as.data.frame(grav_center_vert_coord_list)
      grav_center_vert_asb <- grav_center_vert_coord_list
      
      # plot fdiv:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            ggplot2::aes_(x = sp_coord2[which(
                              sp_coord2[, col_nb] == 'no'),
                                                        c(faxes_nm[[1]])],
                                          y = sp_coord2[which(
                                            sp_coord2[, col_nb] == 'no'),
                                                        c(faxes_nm[[2]])]),
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2, 
                        ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                      y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                      size = sp_faxes_coord_k2$w),
                          colour = color_sp, shape = shape_sp, fill = fill_sp) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::geom_segment(data = sp_faxes_coord_k, 
                           ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                         y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                              xend = grav_center_vert_asb[faxes_nm[[1]],
                                             paste0("grav_center_vert_coord",
                                                    sep = "_", asb_k)],
                              yend = grav_center_vert_asb[faxes_nm[[2]],
                                             paste0("grav_center_vert_coord",
                                                     sep = "_", asb_k)],
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                        ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                        shape = shape_vert, size = size_vert) +
        
        ggplot2::geom_point(data = grav_center_vert_asb,
                        ggplot2::aes_(x = grav_center_vert_asb[faxes_nm[[1]], 
                          paste0("grav_center_vert_coord", sep = "_", asb_k)],
                                        y = grav_center_vert_asb[faxes_nm[[2]],
                          paste0("grav_center_vert_coord", sep = "_", asb_k)]),
                          colour = color_centroid, shape = shape_centroid, 
                          size = size_centroid, fill = fill_centroid) +
        
        ggplot2::theme(legend.position = "none")
      
      return_fdiv_list[[i]] <- plot_k
      
    }
    
    # FEve index:
    
    if ("feve" %in% ind_vect) {
      
      # check that FEve can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "feve"])) {
        stop(paste0("Error: FEve value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve mst information for the studied assemblage:
      mst_asb_k <- eval(parse(text = paste0("mst_list[['mst", sep = '_', asb_k, 
                                            "']]")))
      mst_asb_k <- as.matrix(mst_asb_k)
      mst_asb_k <- as.data.frame(mst_asb_k)
      
      # gather species coordinates into a big dataframe for plotting:
      segment_coord <- data.frame(sp_start = NA, sp_stop = NA)
      j <- 1
      for (m1 in 1:(ncol(mst_asb_k) - 1)) {
        for (m2 in (m1 + 1):nrow(mst_asb_k)) {
          if (mst_asb_k[m2, m1] == 1) {
            segment_coord[j, "sp_start"] <- rownames(mst_asb_k)[m2]
            segment_coord[j, "sp_stop"]  <- colnames(mst_asb_k)[m1]
            j <- j + 1
          }
        }
      }
      
      for (n in (1:nrow(segment_coord))) {
        segment_coord[n, paste0(faxes_nm[[1]], 
            sep = "_", "start")] <- sp_faxes_coord_k[segment_coord$sp_start[n], 
                                                     faxes_nm[[1]]]
        segment_coord[n, paste0(faxes_nm[[2]], 
            sep = "_", "start")] <- sp_faxes_coord_k[segment_coord$sp_start[n], 
                                                      faxes_nm[[2]]]
        segment_coord[n, paste0(faxes_nm[[1]], 
            sep = "_", "stop")] <- sp_faxes_coord_k[segment_coord$sp_stop[n], 
                                                      faxes_nm[[1]]]
        segment_coord[n, paste0(faxes_nm[[2]], 
            sep = "_", "stop")] <- sp_faxes_coord_k[segment_coord$sp_stop[n],
                                                    faxes_nm[[2]]]
      }
      
      # plot FEve:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_segment(data = segment_coord, 
                              ggplot2::aes_(x = segment_coord[, 3],
                                            y = segment_coord[, 4]),
                              xend =  segment_coord[, 5],
                              yend =  segment_coord[, 6],
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2,
                            ggplot2::aes_(x = sp_faxes_coord_k2[, 
                                                                faxes_nm[[1]]],
                                          size = sp_faxes_coord_k2$w), 
                            colour = color_sp,
                            shape = shape_sp, fill = fill_sp) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                        ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                            shape = shape_vert, size = size_vert) +
        ggplot2::theme(legend.position = 'none')
      
      return_feve_list[[i]] <- plot_k
    
    }
    
    
    # FSpe index:
    
    if ("fspe" %in% ind_vect) {
      
      # check that FSpe can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fspe"])) {
        stop(paste0("Error: FSpe value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve coordinates of gravity center of the global pool:
      grav_center_global_pool <- as.data.frame(grav_center_global_pool)
      
      # plot fspe:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            colour = color_sp_gp, shape = shape_sp_gp, fill = fill_sp_gp) +
        
        ggplot2::geom_segment(data = sp_faxes_coord_k,
                              ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                            y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                              xend = grav_center_global_pool[faxes_nm[[1]], ],
                              yend = grav_center_global_pool[faxes_nm[[2]], ],
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2,
                            ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                          y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                          size = sp_faxes_coord_k2$w), colour = color_sp,
                            shape = shape_sp, fill = fill_sp) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                                                     y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, shape = shape_vert, size = size_vert) +
        
        ggplot2::geom_point(data = grav_center_global_pool,
                            ggplot2::aes_(x = grav_center_global_pool[faxes_nm[[1]], ],
                                          y = grav_center_global_pool[faxes_nm[[2]], ]),
                            colour = color_centroid, shape = shape_centroid,
                            size = size_centroid, fill = fill_centroid) +
        
        ggplot2::theme(legend.position='none')
      
      return_fspe_list[[i]] <- plot_k
    }
    
    
    # FDis index:
    
    if ("fdis" %in% ind_vect) {
      
      # check that FDis can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fdis"])) {
        stop(paste0("Error: FDis value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve fdis values
      fide_values <- fd_ind_values[asb_k, c(paste0("fide", sep = '_', 
                                                   faxes_nm[[1]]), 
                                            paste0("fide", sep = '_', 
                                                   faxes_nm[[2]]))]
      
      # plot fdis:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_segment(data = sp_faxes_coord_k,
                              ggplot2::aes_(x = sp_faxes_coord_k[, 
                                                                faxes_nm[[1]]],
                                            y = sp_faxes_coord_k[, 
                                                              faxes_nm[[2]]]),
                              xend = fide_values[, paste0("fide", sep = "_", 
                                                          faxes_nm[[1]])],
                              yend = fide_values[, paste0("fide", sep = "_", 
                                                          faxes_nm[[2]])],
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment) +
        
        ggplot2::geom_hline(yintercept = fide_values[, paste0("fide", sep = "_", 
                                                              faxes_nm[[2]])],
                            linetype = "dotted", color = "red", size = 1) +
        
        ggplot2::geom_vline(xintercept = fide_values[, paste0("fide", sep = "_", 
                                                              faxes_nm[[1]])],
                            linetype = "dotted", color = "red", size = 1) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2,
                            ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                          y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                          size = sp_faxes_coord_k2$w), 
                            colour = color_sp,
                            shape = shape_sp, fill = fill_sp) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                         ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                       y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                            shape = shape_vert, size = size_vert) +
        
        ggplot2::geom_point(data = grav_center_global_pool,
                    ggplot2::aes_(x = fide_values[, 
                                    paste0("fide", sep = "_", faxes_nm[[1]])],
                                  y = fide_values[, 
                                    paste0("fide", sep = "_", faxes_nm[[2]])]),
                            colour = color_centroid, shape = shape_centroid, 
                            size = size_centroid, fill = fill_centroid) +
        
        ggplot2::theme(legend.position = 'none')
      
      return_fdis_list[[i]] <- plot_k
    }
    
    
    # FOri index:
    
    if ('fori' %in% ind_vect) {
      
      # check that FOri can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fori"])) {
        stop(paste0("Error: FOri value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve nearest neighbor identity in the global pool...
      # ... for species present in the studied assemblage:
      nn_asb_k <- nm_nn_global_pool_list[asb_k]
      # create a dataframe that will contain nn global pool id for each...
      # ... sp of the asb:
      nn_asb_k_coord <- data.frame(asb_nm = NA, species_nm = NA, nn_nm = NA,
                                   coord_sp_1 = NA,
                                   coord_sp_2 = NA,
                                   coord_nn_1 = NA,
                                   coord_nn_2 = NA)
      
      # complete the df:
      n <- 1
      for (m in (1:length(nn_asb_k[[1]]))) {
        # complete the df for sp_nm and nn_nm:
        sp_nm <- names(nn_asb_k[[1]][m])
        nn_nm <- nn_asb_k[[1]][[m]]
        # for all the nn of each species (can have several nn):
        for(nn in 1:length(nn_nm)) {
          nn_asb_k_coord[n, "asb_nm"] <- asb_k
          nn_asb_k_coord[n, "species_nm"] <- as.character(sp_nm)
          nn_asb_k_coord[n, "nn_nm"] <- as.character(nn_nm[nn])
          # and add coord:
          nn_asb_k_coord[n, "coord_sp_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[1]]]
          nn_asb_k_coord[n, "coord_sp_2"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[2]]]
          nn_asb_k_coord[n, "coord_nn_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[1]]]
          nn_asb_k_coord[n, "coord_nn_2"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[2]]]
          #Î n
          n <- n + 1
        }
      }
      
      # plot fori:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_segment(data = nn_asb_k_coord, 
                              ggplot2::aes_(x = nn_asb_k_coord[, 4],
                                            y = nn_asb_k_coord[, 5],
                                            xend =  nn_asb_k_coord[, 6],
                                            yend =  nn_asb_k_coord[, 7]),
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment,
                              arrow = grid::arrow(length = grid::unit(0.10, 
                                                                      "inches"),
                                                  ends = "last", 
                                                  type = "open")) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2,
                          ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                        y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                        size = sp_faxes_coord_k2$w),
                            colour = color_sp, shape = shape_sp, 
                            fill = fill_sp) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                        ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                            shape = shape_vert, size = size_vert) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::theme(legend.position = 'none')
      
      return_fori_list[[i]] <- plot_k
    }
    
    
    # FNND index:
    
    if ("fnnd" %in% ind_vect) {
      
      # check that FNND can be computed for this asb:
      if (is.na(fd_ind_values[asb_k, "fnnd"])) {
        stop(paste0("Error: FNND value can not be computed for",
                    sep = " ", asb_k, ".", sep = " ",
                    "The associated figure can not be computed."))
      }
      
      # retrieve nearest neighbor identity in the asb...
      # ... for species present in the studied assemblage:
      nn_asb_k <- nm_nn_asb_list[asb_k]
      # create a dataframe that will contain nn global pool id for each ...
      # ... sp of the asb:
      nn_asb_k_coord <- data.frame(asb_nm = NA, species_nm = NA, nn_nm = NA,
                                   coord_sp_1 = NA,
                                   coord_sp_2 = NA,
                                   coord_nn_1 = NA,
                                   coord_nn_2 = NA)
      
      # complete the df:
      n <- 1
      for (m in (1:length(nn_asb_k[[1]]))) {
        # complete the df for sp_nm and nn_nm:
        sp_nm <- names(nn_asb_k[[1]][m])
        nn_nm <- nn_asb_k[[1]][[m]]
        # for all the nn of each species (can have several nn):
        for(nn in 1:length(nn_nm)) {
          nn_asb_k_coord[n, "asb_nm"] <- asb_k
          nn_asb_k_coord[n, "species_nm"] <- as.character(sp_nm)
          nn_asb_k_coord[n, "nn_nm"] <- as.character(nn_nm[nn])
          # and add coord:
          nn_asb_k_coord[n, "coord_sp_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[1]]]
          nn_asb_k_coord[n, "coord_sp_2"] <- sp_faxes_coord[nn_asb_k_coord[n,
                                                  "species_nm"], faxes_nm[[2]]]
          nn_asb_k_coord[n, "coord_nn_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[1]]]
          nn_asb_k_coord[n, "coord_nn_2"] <- sp_faxes_coord[nn_asb_k_coord[n,
                                                  "nn_nm"], faxes_nm[[2]]]
          # n
          n <- n + 1
        }
      }
      
      
      # plot fnnd:
      plot_k <- get(paste0("plot_funct", sep = '_', i)) +
        
        ggplot2::geom_point(data = sp_coord2[which(sp_coord2[, col_nb] == 'no'),
                                             c(faxes_nm[[1]], faxes_nm[[2]])],
                            colour = color_sp_gp, shape = shape_sp_gp, 
                            fill = fill_sp_gp) +
        
        ggplot2::geom_segment(data = nn_asb_k_coord, 
                              ggplot2::aes_(x = nn_asb_k_coord[, 4],
                                            y = nn_asb_k_coord[, 5],
                                            xend =  nn_asb_k_coord[, 6],
                                            yend =  nn_asb_k_coord[, 7]),
                              colour = color_segment, size = segment_size,
                              linetype = linetype_segment,
                              arrow = grid::arrow(length = grid::unit(0.10, 
                                                                      "inches"),
                                                  ends = "last", 
                                                  type = "open")) +
        
        ggplot2::geom_point(data = sp_faxes_coord_k2,
                          ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                        y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                        size = sp_faxes_coord_k2$w),
                            colour = color_sp, shape = shape_sp, 
                            fill = fill_sp) +
        
        ggplot2::geom_point(data = vert_sp_coord_asb1, 
                        ggplot2::aes_(x = vert_sp_coord_asb1[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb1[, faxes_nm[[2]]]),
                            color = color_vert, fill = fill_vert, 
                            shape = shape_vert, size = size_vert) +
        
        ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
        
        ggplot2::theme(legend.position = 'none')
      
      return_fnnd_list[[i]] <- plot_k
    }
    
  } # end of the loop for the 1st asb
  
  #### Retrieve information, create df and format data for the 2nd asb ####
  
  
  # if there is a 2nd asb to plot:
  if (length(asb_vect) == 2) {
    
    # retrieve assemblage name:
    asb_k <- asb_vect[2]
    
    # retrieve coordinates of species belonging to asb_k:
    sp_filter <- sp.filter(asb_k, sp_faxes_coord, asb_sp_w)
    sp_faxes_coord_k <- sp_filter$`species coordinates`
    asb_sp_relatw_k <- sp_filter$`species relative weight`
    
    # keep species coordinates of assemblage 2 for plotting caption:
    sp_coord_asb_2 <- sp_faxes_coord_k
    
    # format "asb_sp_relatw_k" so easy to plot:
    sp_faxes_coord_k2 <- sp_faxes_coord_k
    sp_faxes_coord_k2 <- as.data.frame(sp_faxes_coord_k2)
    weight <- c()
    for (n in (1:ncol(asb_sp_relatw_k))) {
      weight <- append(weight, asb_sp_relatw_k[1, n])
    }
    sp_faxes_coord_k2$w <- weight
    sp_faxes_coord_k2 <- as.data.frame(sp_faxes_coord_k2)
    
    # keep species weight of assemblage 2 for plotting caption:
    sp_coord_asb_w_2 <- sp_faxes_coord_k2$w
    
    # create a table that will sum up coordinates of all species to plot...
    # ... (global pool) and add a column to know if there are...
    # ... present on the assemblage to plot:
    sp_coord2 <- sp_faxes_coord
    sp_coord2$asb <- rownames(sp_coord2)
    sp_coord2$asb[which(! sp_coord2$asb %in% rownames(sp_faxes_coord_k))] <- "no"
    sp_coord2$asb[which(sp_coord2$asb %in% rownames(sp_faxes_coord_k))] <- asb_k
    sp_coord2$asb <- as.factor(sp_coord2$asb)
    
    # count the number of column of this dataframe used to retrieve last ...
    # ... column when plotting:
    col_nb <- ncol(sp_coord2)
    
    
    # Do the plot for each index and along asked dimensions for the 2nd asb ####
    
    
    # span across the plots (up to 6 plots for 4 dimensions):
    for (i in (1:(plot_nb - 1))) {
      
      ## retrieve information on the graph to plot:
      
      # retrieve the plot to work on:
      assign(paste0("plot_funct", sep = "_", i), funct_space_output[[i]])
      
      # retrieve the combination of dimensions of the plot:
      faxes_nm <- c(eval(parse(text = paste0("plot_funct", 
                                             sep = "_", i, 
                                             "$labels$x"))), 
                    eval(parse(text = (paste0("plot_funct", sep = "_", i, 
                                              "$labels$y")))))
      
      # retrieve vertices names of asb_k along the 2 dimensions to be ...
      # ... plotted before all indices and not just for FRic because needed ...
      # ... if vertices number is to show in graphs caption for all indices:
      vert_nm_asb_k <- vertices(sp_faxes_coord_k[, c(faxes_nm[[1]], 
                                                     faxes_nm[[2]])], 
                                check_input = TRUE)
      
      # retrieve coordinates of species of asb_k along the 2 plotted dimensions:
      vert_sp_faxes_coord_k <- sp_faxes_coord_k[which(rownames(sp_faxes_coord_k) 
                                                      %in% vert_nm_asb_k),
                                                c(faxes_nm[[1]], faxes_nm[[2]])]
      
      # then vertices must be ordered so that the convex hull...
      #... can be printed (outside path): sort them clockwise:
      # (https://stackoverflow.com/questions/48249540/plot-convex-hull-given-
      # ... by-quickhull-algorithm-in-r-convhulln-function)
      # find the gravity center of present species and then compute the ...
      # ... angle value of each species to this points, then order angle values:
      vert_sp_faxes_coord_k <- vert_sp_faxes_coord_k[order(-1 * atan2(
        vert_sp_faxes_coord_k[, 
         faxes_nm[[2]]] - mean(range(vert_sp_faxes_coord_k[, 
                                                           faxes_nm[[2]]])),
        vert_sp_faxes_coord_k[, 
         faxes_nm[[1]]] - mean(range(vert_sp_faxes_coord_k[, 
                                                          faxes_nm[[1]]])))), ]
      
      # convert the format so that it can be used with ggplot2:
      vert_sp_faxes_coord_k <- as.data.frame(vert_sp_faxes_coord_k)
      sp_faxes_coord_k <- as.data.frame(sp_faxes_coord_k)
      
      # save the vert_sp_coord for caption (asb2)
      vert_sp_coord_asb2 <- vert_sp_faxes_coord_k
      
      
      ## plot each index:
      
      # FRic index:
      
      if ("fric" %in% ind_vect) {
        
        # check that FRic can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fric"])) {
          stop(paste0("Error: FRic value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        plot_k <- get("return_fric_list")[[i]] +
          
          ggplot2::geom_polygon(data = vert_sp_faxes_coord_k,
                      ggplot2::aes_(x = vert_sp_faxes_coord_k[, faxes_nm[[1]]],
                                    y = vert_sp_faxes_coord_k[, faxes_nm[[2]]]),
                                fill = color_sp_asb2, alpha = alpha_ch) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k, 
                      ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                   y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                              colour = color_sp_asb2, shape = shape_sp_asb2, 
                              fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                        ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                        shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::theme(legend.position = "none")
        
        return_fric_list[[i]] <- plot_k
        
      }
      
      # FDiv index:
      
      if ("fdiv" %in% ind_vect) {
        
        # check that FDiv can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fdiv"])) {
          stop(paste0("Error: FDiv value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve gravity center of vertices for the studied assemblage:
        grav_center_vert_coord_list <- 
          as.data.frame(grav_center_vert_coord_list)
        grav_center_vert_asb <- grav_center_vert_coord_list
        
        # plot fdiv:
        plot_k <- get("return_fdiv_list")[[i]] +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2, 
                    ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                  y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                  size = sp_faxes_coord_k2$w),
                              colour = color_sp_asb2, shape = shape_sp_asb2, 
                              fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          ggplot2::geom_segment(data = sp_faxes_coord_k, 
                    ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                  y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                                xend = grav_center_vert_asb[faxes_nm[[1]],
                                             paste0("grav_center_vert_coord",
                                             sep = "_", asb_k)],
                                yend = grav_center_vert_asb[faxes_nm[[2]],
                                             paste0("grav_center_vert_coord",
                                             sep = "_", asb_k)],
                                colour = color_segment_asb2, 
                                size = segment_size_asb2,
                                linetype = linetype_segment_asb2) +
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                              ggplot2::aes_(x = vert_sp_coord_asb2[, 
                                                                faxes_nm[[1]]],
                                            y = vert_sp_coord_asb2[, 
                                                                faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                              shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::geom_point(data = grav_center_vert_asb,
                        ggplot2::aes_(x = grav_center_vert_asb[faxes_nm[[1]], 
                          paste0("grav_center_vert_coord", sep = "_", asb_k)],
                                      y = grav_center_vert_asb[faxes_nm[[2]], 
                          paste0("grav_center_vert_coord", sep = "_", asb_k)]),
                          colour = color_centroid_asb2, 
                          shape = shape_centroid_asb2, size = size_centroid_asb2,
                          fill = fill_centroid_asb2) +
          
          ggplot2::theme(legend.position = "none")
        
        return_fdiv_list[[i]] <- plot_k
        
      }
      
      # FEve index:
      
      if ("feve" %in% ind_vect) {
        
        # check that FEve can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "feve"])) {
          stop(paste0("Error: FEve value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve mst information for the studied assemblage:
        mst_asb_k <- eval(parse(text = paste0('mst_list$mst', sep = '_', 
                                              asb_k)))
        mst_asb_k <- as.matrix(mst_asb_k)
        mst_asb_k <- as.data.frame(mst_asb_k)
        
        # gather species coordinates into a big dataframe for plotting:
        segment_coord <- data.frame(sp_start = NA, sp_stop = NA)
        j <- 1
        for (m1 in 1:(ncol(mst_asb_k) - 1)) {
          for (m2 in (m1 + 1):nrow(mst_asb_k)) {
            if (mst_asb_k[m2, m1] == 1) {
              segment_coord[j, "sp_start"] <- rownames(mst_asb_k)[m2]
              segment_coord[j, "sp_stop"]  <- colnames(mst_asb_k)[m1]
              j <- j + 1
            }
          }
        }
        
        for (n in (1:nrow(segment_coord))) {
          segment_coord[n, 
               paste0(faxes_nm[[1]], sep = '_', 
                      "start")] <- sp_faxes_coord_k[segment_coord$sp_start[n], 
                                                    faxes_nm[[1]]]
          segment_coord[n,
               paste0(faxes_nm[[2]], sep = '_', 
                      "start")] <- sp_faxes_coord_k[segment_coord$sp_start[n], 
                                                    faxes_nm[[2]]]
          segment_coord[n,
                paste0(faxes_nm[[1]], sep = '_', 
                       "stop")] <- sp_faxes_coord_k[segment_coord$sp_stop[n], 
                                                    faxes_nm[[1]]]
          segment_coord[n, 
                paste0(faxes_nm[[2]], sep = '_', 
                       "stop")] <- sp_faxes_coord_k[segment_coord$sp_stop[n], 
                                                    faxes_nm[[2]]]
        }
        
        # plot FEve:
        plot_k <- get("return_feve_list")[[i]] +
          
          ggplot2::geom_segment(data = segment_coord, 
                                ggplot2::aes_(x = segment_coord[, 3],
                                              y = segment_coord[, 4]),
                                xend =  segment_coord[, 5],
                                yend =  segment_coord[, 6],
                                colour = color_segment_asb2, 
                                size = segment_size_asb2,
                                linetype = linetype_segment_asb2) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2,
                           ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                        size = sp_faxes_coord_k2$w), 
                              colour = color_sp_asb2,
                              shape = shape_sp_asb2, fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                        ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                              shape = shape_vert_asb2, size = size_vert) +
          ggplot2::theme(legend.position = 'none')
        
        return_feve_list[[i]] <- plot_k
        
        
      }
      
      
      # FSpe index:
      
      if ("fspe" %in% ind_vect) {
        
        # check that FSpe can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fspe"])) {
          stop(paste0("Error: FSpe value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve coordinates of gravity center of the global pool:
        grav_center_global_pool <- as.data.frame(grav_center_global_pool)
        
        # plot fspe:
        plot_k <- get("return_fspe_list")[[i]] +
          
          ggplot2::geom_segment(data = sp_faxes_coord_k,
                          ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                        y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                                xend = grav_center_global_pool[faxes_nm[[1]], ],
                                yend = grav_center_global_pool[faxes_nm[[2]], ],
                                colour = color_segment_asb2, 
                                size = segment_size_asb2,
                                linetype = linetype_segment_asb2) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2,
                        ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                      y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                      size = sp_faxes_coord_k2$w), 
                                      colour = color_sp_asb2,
                                      shape = shape_sp_asb2, 
                                      fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                        ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                              shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::geom_point(data = grav_center_global_pool,
                  ggplot2::aes_(x = grav_center_global_pool[faxes_nm[[1]], ],
                                y = grav_center_global_pool[faxes_nm[[2]], ]),
                              colour = color_centroid_asb2, 
                              shape = shape_centroid_asb2,
                              size = size_centroid_asb2, 
                              fill = fill_centroid_asb2) +
          
          ggplot2::theme(legend.position='none')
        
        return_fspe_list[[i]] <- plot_k
      }
      
      
      # FDis index:
      
      if ("fdis" %in% ind_vect) {
        
        # check that FDis can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fdis"])) {
          stop(paste0("Error: FDis value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve fdis values
        fide_values <- fd_ind_values[asb_k, 
                                     c(paste0("fide", sep = '_', faxes_nm[[1]]), 
                                      paste0("fide", sep = '_', faxes_nm[[2]]))]
        
        # plot fdis:
        plot_k <- get("return_fdis_list")[[i]] +
          
          ggplot2::geom_segment(data = sp_faxes_coord_k,
                           ggplot2::aes_(x = sp_faxes_coord_k[, faxes_nm[[1]]],
                                         y = sp_faxes_coord_k[, faxes_nm[[2]]]),
                xend = fide_values[, paste0("fide", sep = "_", faxes_nm[[1]])],
                yend = fide_values[, paste0("fide", sep = "_", faxes_nm[[2]])],
                colour = color_segment_asb2, size = segment_size_asb2,
                linetype = linetype_segment_asb2) +
          
          ggplot2::geom_hline(yintercept = fide_values[, paste0("fide", 
                                                                sep = "_", 
                                                                faxes_nm[[2]])],
                              linetype = "dotted", color = "red", size = 1) +
          
          ggplot2::geom_vline(xintercept = fide_values[, paste0("fide", 
                                                                sep = "_", 
                                                                faxes_nm[[1]])],
                              linetype = "dotted", color = "red", size = 1) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2,
                        ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                      y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                      size = sp_faxes_coord_k2$w), 
                                      colour = color_sp_asb2,
                                      shape = shape_sp_asb2, 
                                      fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                        ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                              shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::geom_point(data = grav_center_global_pool,
                  ggplot2::aes_(x = fide_values[, 
                                                paste0("fide", sep = "_", 
                                                       faxes_nm[[1]])],
                                y = fide_values[, paste0("fide", sep = "_",
                                                         faxes_nm[[2]])]),
                              colour = color_centroid_asb2, 
                              shape = shape_centroid_asb2, 
                              size = size_centroid_asb2,
                              fill = fill_centroid_asb2) +
          
          ggplot2::theme(legend.position='none')
        
        return_fdis_list[[i]] <- plot_k
      }
      
      
      # FOri index:
      
      if ('fori' %in% ind_vect) {
        
        # check that FOri can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fori"])) {
          stop(paste0("Error: FOri value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve nearest neighbor identity in the global pool...
        # ... for species present in the studied assemblage:
        nn_asb_k <- nm_nn_global_pool_list[asb_k]
        # create a dataframe that will contain nn global pool id for ...
        # ... each sp of the asb:
        nn_asb_k_coord <- data.frame(asb_nm = NA, species_nm = NA, nn_nm = NA,
                                     coord_sp_1 = NA,
                                     coord_sp_2 = NA,
                                     coord_nn_1 = NA,
                                     coord_nn_2 = NA)
        
        # complete the df:
        n <- 1
        for (m in (1:length(nn_asb_k[[1]]))) {
          # complete the df for sp_nm and nn_nm:
          sp_nm <- names(nn_asb_k[[1]][m])
          nn_nm <- nn_asb_k[[1]][[m]]
          # for all the nn of each species (can have several nn):
          for(nn in 1:length(nn_nm)) {
            nn_asb_k_coord[n, "asb_nm"] <- asb_k
            nn_asb_k_coord[n, "species_nm"] <- as.character(sp_nm)
            nn_asb_k_coord[n, "nn_nm"] <- as.character(nn_nm[nn])
            # and add coord:
            nn_asb_k_coord[n, "coord_sp_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[1]]]
            nn_asb_k_coord[n, "coord_sp_2"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[2]]]
            nn_asb_k_coord[n, "coord_nn_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[1]]]
            nn_asb_k_coord[n, "coord_nn_2"] <- sp_faxes_coord[nn_asb_k_coord[n,
                                                  "nn_nm"], faxes_nm[[2]]]
            #Î n
            n <- n + 1
          }
        }
        
        
        # plot fori:
        plot_k <- get("return_fori_list")[[i]] +
          
          ggplot2::geom_segment(data = nn_asb_k_coord, 
                                ggplot2::aes_(x = nn_asb_k_coord[, 4],
                                              y = nn_asb_k_coord[, 5],
                                              xend =  nn_asb_k_coord[, 6],
                                              yend =  nn_asb_k_coord[, 7]),
                                colour = color_segment_asb2, 
                                size = segment_size_asb2,
                                linetype = linetype_segment_asb2,
                                arrow = grid::arrow(length = grid::unit(0.10, 
                                                                        "inches"),
                                                    ends = "last", 
                                                    type = "open")) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2,
                          ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                        y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                            size = sp_faxes_coord_k2$w),
                              colour = color_sp_asb2, shape = shape_sp_asb2, 
                              fill = fill_sp_asb2) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                         ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                       y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, 
                              fill = fill_vert_asb2, 
                              shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::theme(legend.position='none')
        
        return_fori_list[[i]] <- plot_k
      }
      
      
      # FNND index:
      
      if ("fnnd" %in% ind_vect) {
        
        # check that FNND can be computed for this asb:
        if (is.na(fd_ind_values[asb_k, "fnnd"])) {
          stop(paste0("Error: FNND value can not be computed for",
                      sep = " ", asb_k, ".", sep = " ",
                      "The associated figure can not be computed."))
        }
        
        # retrieve nearest neighbor identity in the asb...
        # ... for species present in the studied assemblage:
        nn_asb_k <- nm_nn_asb_list[asb_k]
        # create a dataframe that will contain nn global pool id for each ...
        # ... sp of the asb:
        nn_asb_k_coord <- data.frame(asb_nm = NA, species_nm = NA, nn_nm = NA,
                                     coord_sp_1 = NA,
                                     coord_sp_2 = NA,
                                     coord_nn_1 = NA,
                                     coord_nn_2 = NA)
        
        # complete the df:
        n <- 1
        for (m in (1:length(nn_asb_k[[1]]))) {
          # complete the df for sp_nm and nn_nm:
          sp_nm <- names(nn_asb_k[[1]][m])
          nn_nm <- nn_asb_k[[1]][[m]]
          # for all the nn of each species (can have several nn):
          for(nn in 1:length(nn_nm)) {
            nn_asb_k_coord[n, "asb_nm"] <- asb_k
            nn_asb_k_coord[n, "species_nm"] <- as.character(sp_nm)
            nn_asb_k_coord[n, "nn_nm"] <- as.character(nn_nm[nn])
            # and add coord:
            nn_asb_k_coord[n, "coord_sp_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[1]]]
            nn_asb_k_coord[n, "coord_sp_2"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "species_nm"], faxes_nm[[2]]]
            nn_asb_k_coord[n, "coord_nn_1"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[1]]]
            nn_asb_k_coord[n, "coord_nn_2"] <- sp_faxes_coord[nn_asb_k_coord[n, 
                                                  "nn_nm"], faxes_nm[[2]]]
            n <- n + 1
          }
        }
        
        # plot fnnd:
        plot_k <- get("return_fnnd_list")[[i]] +
          
          ggplot2::geom_segment(data = nn_asb_k_coord, 
                                ggplot2::aes_(x = nn_asb_k_coord[, 4],
                                              y = nn_asb_k_coord[, 5],
                                              xend =  nn_asb_k_coord[, 6],
                                              yend =  nn_asb_k_coord[, 7]),
                                colour = color_segment_asb2, 
                                size = segment_size_asb2,
                                linetype = linetype_segment_asb2,
                                arrow = grid::arrow(length = grid::unit(0.10, 
                                                                        "inches"),
                                                    ends = "last", 
                                                    type = "open")) +
          
          ggplot2::geom_point(data = sp_faxes_coord_k2,
                         ggplot2::aes_(x = sp_faxes_coord_k2[, faxes_nm[[1]]],
                                       y = sp_faxes_coord_k2[, faxes_nm[[2]]],
                                       size = sp_faxes_coord_k2$w),
                              colour = color_sp_asb2, shape = shape_sp_asb2, 
                         fill = fill_sp_asb2) +
          
          ggplot2::geom_point(data = vert_sp_coord_asb2, 
                        ggplot2::aes_(x = vert_sp_coord_asb2[, faxes_nm[[1]]],
                                      y = vert_sp_coord_asb2[, faxes_nm[[2]]]),
                              color = color_vert_asb2, fill = fill_vert_asb2, 
                        shape = shape_vert_asb2, size = size_vert) +
          
          ggplot2::scale_size(range = c(scale_inf, scale_sup)) +
          
          ggplot2::theme(legend.position='none')
        
        return_fnnd_list[[i]] <- plot_k
      }
      
      
    } # end of the loop on functional plots
    
  } # end of the loop if 2nd asb
  
  
  #### Plot for caption (legend and summary) ####
  
  
  ## If 1 asb to plot:
  
  if (length(asb_vect) == 1) {
    
    spread_faxes <- (range_faxes_lim[2] - range_faxes_lim[1])
    
    plot_caption0 <- ggplot2::ggplot(sp_faxes_coord) +
      ggplot2::scale_x_continuous(limits = range_faxes_lim, expand = c(0,0)) +
      ggplot2::scale_y_continuous(limits = range_faxes_lim, expand = c(0,0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect(xmin = range_faxes_lim[1], xmax = range_faxes_lim[2],
                         ymin = range_faxes_lim[1], ymax = range_faxes_lim[2],
                         fill = "white", colour = "black")
    
    # if indices to plot are not reduced to FRic:
    plot_caption <- plot_caption0 +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.10,
                          colour = color_sp, fill = fill_sp, shape = shape_sp,
                          size = size_sp) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.10,
                         label = paste0(nrow(sp_coord_asb_1), " species for ", 
                                        asb_vect[1]),
                         colour = color_sp, hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.25,
                          colour = color_vert, shape = shape_vert,
                          size = size_vert, fill = fill_vert) +
      
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.25,
                         label = c("vertices"),
                         colour = color_vert, hjust = 0) +
      
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.40,
                         label = paste0("plotted along ", length(faxes),
                                        "axes from the ",
                                        ncol(sp_coord_asb_1), 
                                        "-dimensional space"),  
                         hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.55,
                          colour = color_sp, shape = shape_sp,
                          size = scale_inf, fill = fill_sp) +
      
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.55,
                         label = paste0("relative weight", sep = " ", 
                                        min(sp_coord_asb_w_1)),
                         colour = color_sp, hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.70,
                          colour = color_sp, fill = fill_sp, shape = shape_sp,
                          size = scale_sup) +
      
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.70,
                         label = paste0("relative weight", sep = " ", 
                                        max(sp_coord_asb_w_1)),
                         colour = color_sp, hjust = 0)
    
    # if FRic index to plot (caption different from other indices):
    if ("fric" %in% ind_vect) {
      
      plot_caption_fric <- plot_caption0 +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.25,
                            colour = color_sp, fill = fill_sp, shape = shape_sp,
                            size = size_sp) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.25,
                           label = paste0(nrow(sp_coord_asb_1), " species for ", 
                                          asb_vect[1]),
                           colour = color_sp, hjust = 0) +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.45,
                            colour = color_vert, shape = shape_vert,
                            size = size_vert, fill = fill_vert) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.45,
                           label = c("vertices"),
                           colour = color_vert, hjust = 0) +
        
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.65,
                           label = paste0("plotted along ", length(faxes),
                                          "axes from the ",
                                          ncol(sp_coord_asb_1), 
                                          "-dimensional space"),
                           hjust = 0) +
        
        ggplot2::geom_rect(xmin = range_faxes_lim[1] + spread_faxes*0.05,
                           xmax = range_faxes_lim[1] + spread_faxes*0.15,
                           ymin = range_faxes_lim[2] - spread_faxes*0.88,
                           ymax = range_faxes_lim[2] - spread_faxes*0.92,
                           fill = color_sp, alpha = alpha_ch) +
        
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.20,
                           y = range_faxes_lim[2] - spread_faxes*0.90,
                           label = paste0("convex hull of", sep = " ", 
                                          asb_vect[1]),
                           colour = color_sp, hjust = 0)
    }
    
  } # end loop if one asb
  
  
  ## If 2 asb to plot:
  
  if (length(asb_vect) == 2) {
    spread_faxes <- (range_faxes_lim[2]- range_faxes_lim[1])
    
    plot_caption0 <- ggplot2::ggplot(sp_faxes_coord) +
      ggplot2::scale_x_continuous(limits = range_faxes_lim, expand = c(0,0)) +
      ggplot2::scale_y_continuous(limits = range_faxes_lim, expand = c(0,0)) +
      ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_rect( xmin = range_faxes_lim[1], xmax = range_faxes_lim[2],
                          ymin = range_faxes_lim[1], ymax = range_faxes_lim[2],
                          fill = "white", colour = "black")
    
    # if indices to plot are not reduced to FRic:
    
    plot_caption <- plot_caption0 +
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.10,
                          colour = color_sp, fill = fill_sp, shape = shape_sp,
                          size = size_sp) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.10,
                         label = paste0(nrow(sp_coord_asb_1), " species for ", 
                                        asb_vect[1]),
                         colour = color_sp, hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.20,
                          colour = color_sp_asb2, fill = fill_sp_asb2,
                          shape = shape_sp_asb2, size = size_sp) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.20,
                         label = paste0(nrow(sp_coord_asb_2), " species for ", 
                                        asb_vect[2]),
                         colour = color_sp_asb2, hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.30,
                          colour = color_vert, shape = shape_vert,
                          size = size_vert, fill = fill_vert) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.30,
                         label = paste0("vertices of", sep = " ", asb_vect[1]),
                         colour = color_vert, hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.40,
                          colour = color_vert_asb2, shape = shape_vert_asb2,
                          size = size_vert, fill = fill_vert_asb2) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.40,
                         label = paste0("vertices of", sep = " ", asb_vect[2]),
                         colour = color_vert_asb2, hjust = 0) +
      
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.60,
                         label = paste0("plotted along ", length(faxes),
                                        "axes from the ",
                                        ncol(sp_faxes_coord), 
                                        "-dimensional space"),
                         hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.80,
                          colour = "black", fill = "black", shape = shape_sp,
                          size = scale_inf) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.80,
                         label = paste0("relative weight", sep = " ",
                                        min(min(sp_coord_asb_w_1), 
                                            min(sp_coord_asb_w_2))),
                         colour = "black", hjust = 0) +
      
      ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                          y = range_faxes_lim[2] - spread_faxes*0.90,
                          colour = "black", fill = "black", shape = shape_sp,
                          size = scale_sup) +
      ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                         y = range_faxes_lim[2] - spread_faxes*0.90,
                         label = paste0("relative weight", sep = " ",
                                        max(max(sp_coord_asb_w_1), 
                                            max(sp_coord_asb_w_2))),
                         colour = "black", hjust = 0)
    
    # if FRic index to plot (caption different from other indices):
    if ("fric" %in% ind_vect) {
      
      plot_caption_fric <- plot_caption0 +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.10,
                            colour = color_sp, fill = fill_sp, shape = shape_sp,
                            size = size_sp) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.10,
                           label = paste0(nrow(sp_coord_asb_1), " species"),
                           colour = color_sp, hjust = 0) +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.20,
                            colour = color_sp_asb2, fill = fill_sp_asb2,
                            shape = shape_sp_asb2,
                            size = size_sp) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.20,
                           label = paste0(nrow(sp_coord_asb_2), " species"),
                           colour = color_sp_asb2, hjust = 0) +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.30,
                            colour = color_vert, shape = shape_vert,
                            size = size_vert, fill = fill_vert) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.30,
                           label = paste0("vertices of", sep = " ", 
                                          asb_vect[1]),
                           colour = color_vert, hjust = 0) +
        
        ggplot2::geom_point(x = range_faxes_lim[1] + spread_faxes*0.1,
                            y = range_faxes_lim[2] - spread_faxes*0.40,
                            colour = color_vert_asb2, shape = shape_vert_asb2,
                            size = size_vert, fill = fill_vert_asb2) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.40,
                           label = paste0("vertices of", sep = " ", 
                                          asb_vect[2]),
                           colour = color_vert_asb2, hjust = 0) +
        
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.15,
                           y = range_faxes_lim[2] - spread_faxes*0.50,
                           label = paste0("plotted along ", length(faxes),
                                          "axes from the ",
                                          ncol(sp_faxes_coord), 
                                          "-dimensional space"),
                           hjust = 0) +
        
        ggplot2::geom_rect(xmin = range_faxes_lim[1] + spread_faxes*0.05,
                           xmax = range_faxes_lim[1] + spread_faxes*0.15,
                           ymin = range_faxes_lim[2] - spread_faxes*0.63,
                           ymax = range_faxes_lim[2] - spread_faxes*0.67,
                           fill = color_sp, alpha = alpha_ch) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.20,
                           y = range_faxes_lim[2] - spread_faxes*0.65,
                           label = paste0("convex hull of", sep = " ", 
                                          asb_vect[1]),
                           colour = color_sp, hjust = 0) +
        
        ggplot2::geom_rect(xmin = range_faxes_lim[1] + spread_faxes*0.05,
                           xmax = range_faxes_lim[1] + spread_faxes*0.15,
                           ymin = range_faxes_lim[2] - spread_faxes*0.83,
                           ymax = range_faxes_lim[2] - spread_faxes*0.87,
                           fill = color_sp_asb2, alpha = alpha_ch) +
        ggplot2::geom_text(x = range_faxes_lim[1] + spread_faxes*0.20,
                           y = range_faxes_lim[2] - spread_faxes*0.85,
                           label = paste0("convex hull of", sep = " ", 
                                          asb_vect[2]),
                           colour = color_sp_asb2, hjust = 0)
    }
  } # end of if 2 asb
  
  #### Create patchwork objects ####
  
  ## create a list that will contain all plots for given indices:
  return_panels_list <- list()
  
  
  # if fric asked:
  if ("fric" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fric <- return_fric_list[[1]] + plot_caption_fric +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FRic representation for",
                                                  sep = " ", asb_k,
                                                  sep = " ","FRic =",
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fric"], 
                                                        4)),
                                   subtitle = 
                                     "2D projection of the multidimensional 
                                   convex-hull",
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fric <- (return_fric_list[[1]] + plot_caption_fric + 
                             return_fric_list[[2]] + return_fric_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FRic representation for",
                                                  sep = " ", asb_k,
                                                  sep = " ","FRic =",
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fric"], 
                                                        4)),
                                   subtitle = 
                                     "2D projection of the multidimensional 
                                   convex-hull",
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fric <- (return_fric_list[[1]] + plot_caption_fric +
                             patchwork::plot_spacer() + return_fric_list[[2]] +
                             return_fric_list[[3]] + patchwork::plot_spacer() +
                             return_fric_list[[4]] + return_fric_list[[5]] +
                             return_fric_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FRic representation for",
                                                  sep = " ", asb_k,
                                                  sep = " ","FRic =",
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fric"], 
                                                        4)),
                                   subtitle = "2D projection of the 
                                   multidimensional convex-hull",
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list,  return_plot_fric)
    
  }
  
  
  # if fdiv asked:
  if ("fdiv" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fdiv <- return_fdiv_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FDiv representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FDiv =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdiv"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   of the vertices of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fdiv <- (return_fdiv_list[[1]] + plot_caption +
                             return_fdiv_list[[2]] + return_fdiv_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FDiv representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FDiv =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdiv"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   of the vertices of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fdiv <- (return_fdiv_list[[1]] + plot_caption + 
                             patchwork::plot_spacer() +
                             return_fdiv_list[[2]] + return_fdiv_list[[3]] +
                             patchwork::plot_spacer() + return_fdiv_list[[4]] +
                             return_fdiv_list[[5]] +  return_fdiv_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FDiv representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FDiv =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdiv"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   of the vertices of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_fdiv)
    
  }
  
  
  # if feve asked:
  if ("feve" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_feve <- return_feve_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FEve representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FEve =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "feve"], 
                                                        4)),
                                   subtitle = "Minimum Spanning Tree linking 
                                   species of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_feve <- (return_feve_list[[1]] + plot_caption +
                             return_feve_list[[2]] + return_feve_list[[3]])+
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FEve representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FEve =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "feve"], 
                                                        4)),
                                   subtitle = "Minimum Spanning Tree linking 
                                   species of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_feve <- (return_feve_list[[1]] + plot_caption +
                             patchwork::plot_spacer() + return_feve_list[[2]] +
                             return_feve_list[[3]] +
                             patchwork::plot_spacer() + return_feve_list[[4]] +
                             return_feve_list[[5]] + return_feve_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FEve representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FEve =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "feve"], 
                                                        4)),
                                   subtitle = "Minimum Spanning Tree linking 
                                   species of the assemblage", 
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_feve)
    
  }
  
  
  # if fspe asked:
  if ("fspe" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fspe <- return_fspe_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FSpe representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FSpe =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fspe"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   OF the vertices of the global pool", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fspe <- (return_fspe_list[[1]] + plot_caption +
                             return_fspe_list[[2]] + return_fspe_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FSpe representation for", 
                                                  sep = " ",
                                                  asb_k, sep = " ", "FSpe =",
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fspe"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   of the vertices of the global pool",
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fspe <- (return_fspe_list[[1]] + plot_caption + 
                             patchwork::plot_spacer() +
                             return_fspe_list[[2]] + return_fspe_list[[3]] + 
                             patchwork::plot_spacer() +
                             return_fspe_list[[4]] + return_fspe_list[[5]] + 
                             return_fspe_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FSpe representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FSpe =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fspe"], 
                                                        4)),
                                   subtitle = "Distances to the gravity center 
                                   of the vertices of the global pool", 
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_fspe)
    
  }
  
  
  # if fdis asked:
  if ("fdis" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fdis <- return_fdis_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FDis representation for", 
                                                  sep = " ", asb_k, sep = " ", 
                                                  "FDis =", sep = " ", 
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdis"], 
                                                        4)),
                                   subtitle = "Distances to the average position 
                                   of species present in the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fdis <- (return_fdis_list[[1]] + plot_caption + 
                             return_fdis_list[[2]] +
                             return_fdis_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FDis representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FDis =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdis"], 
                                                        4)),
                                   subtitle = "Distances to the average position 
                                   of species present in the assemblage",
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fdis <- (return_fdis_list[[1]] + plot_caption + 
                             patchwork::plot_spacer() +
                             return_fdis_list[[2]] + return_fdis_list[[3]] + 
                             patchwork::plot_spacer() +
                             return_fdis_list[[4]] + return_fdis_list[[5]] + 
                             return_fdis_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FDis representation for", 
                                                  sep = " ",
                                                  asb_k, sep = " ", "FDis =",
                                                  sep = " ", 
                                                  round(fd_ind_values[asb_k, 
                                                                      "fdis"], 
                                                        4)),
                                   subtitle = "Distances to the average position 
                                   of species present in the assemblage", 
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_fdis)
    
  }
  
  
  # if fori asked:
  if ("fori" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fori <- return_fori_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FOri representation for", 
                                                  sep = " ", asb_k, sep = " ", 
                                                  "FOri =", sep = " ", 
                                                  round(fd_ind_values[asb_k, 
                                                                      "fori"], 
                                                        4)), 
                                   subtitle = "Distances to nearest neighbour 
                                   from the global pool", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fori <- (return_fori_list[[1]] + plot_caption + 
                             return_fori_list[[2]] +
                             return_fori_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FOri representation for", 
                                                  sep = " ",
                                                  asb_k, sep = " ", "FOri =", 
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fori"], 
                                                        4)),
                                   subtitle = "Distances to nearest neighbour 
                                   from the global pool",
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fori <- (return_fori_list[[1]] + plot_caption + 
                             patchwork::plot_spacer() +
                             return_fori_list[[2]] + return_fori_list[[3]] + 
                             patchwork::plot_spacer() +
                             return_fori_list[[4]] + return_fori_list[[5]] + 
                             return_fori_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FOri representation for",
                                                  sep = " ", asb_k, sep = " ",
                                                  "FOri =", sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fori"], 
                                                        4)),
                                   subtitle = "Distances to nearest neighbour 
                                   from the global pool",
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_fori)
    
  }
  
  # if fnnd asked:
  if ("fnnd" %in% ind_vect) {
    
    # if 1 plot to return (combination of 2 dimensions only):
    if (plot_nb == 2) {
      return_plot_fnnd <- return_fnnd_list[[1]] + plot_caption +
        patchwork::plot_layout(byrow = TRUE, heights = c(1), widths = c(1,1),
                               ncol = 2, nrow = 1 ) +
        patchwork::plot_annotation(title = paste0("FNND representation for", 
                                                  sep = " ", asb_k, sep = " ", 
                                                  "FNND =", sep = " ", 
                                                  round(fd_ind_values[asb_k, 
                                                                      "fnnd"], 
                                                        4)), 
                                   subtitle = "Distances to nearest neighbour 
                                   from the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 3 plots to return (combination of up to 3 dimensions):
    if (plot_nb == 4) {
      return_plot_fnnd <- (return_fnnd_list[[1]] + plot_caption + 
                             return_fnnd_list[[2]] +
                             return_fnnd_list[[3]]) +
        patchwork::plot_layout(byrow = TRUE, heights = c(1, 1), widths = c(1,1),
                               ncol = 2, nrow = 2 ) +
        patchwork::plot_annotation(title = paste0("FNND representation for", 
                                                  sep = " ",
                                                  asb_k, sep = " ", "FNND =", 
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fnnd"], 
                                                        4)),
                                   subtitle = "Distances to nearest neighbour 
                                   from the assemblage", 
                                   caption = "made with mFD package")
    }
    
    # if 6 plots to return (combination of up to 4 dimensions):
    if (plot_nb == 7) {
      return_plot_fnnd <- (return_fnnd_list[[1]] + plot_caption + 
                             patchwork::plot_spacer() +
                             return_fnnd_list[[2]] + return_fnnd_list[[3]] + 
                             patchwork::plot_spacer() +
                             return_fnnd_list[[4]] + return_fnnd_list[[5]] + 
                             return_fnnd_list[[6]]) +
        patchwork::plot_layout(byrow = TRUE, heights = rep(1, 3), 
                               widths = rep(1, 3),
                               ncol = 3, nrow = 3 ) +
        patchwork::plot_annotation(title = paste0("FNND representation for", 
                                                  sep = " ",
                                                  asb_k, sep = " ", "FNND =", 
                                                  sep = " ",
                                                  round(fd_ind_values[asb_k, 
                                                                      "fnnd"], 
                                                        4)),
                                   subtitle = "Distances to nearest neighbour 
                                   from the assemblage",
                                   caption = "made with mFD package")
    }
    
    return_panels_list <- rlist::list.append(return_panels_list, 
                                             return_plot_fnnd)
    
  }
  
  ## returning output ####
  
  # only keep the summed up plots in return_panels_list ie the last 
  # "nb of ind" ones:
  return_panels_list <- utils::tail(return_panels_list, length(ind_vect))
  
  # type, resolution and dimensions of file if to be saved
  device_file = "png"
  res_file= 300
  
  # displaying or saving
  if (is.null(name_file) == TRUE )  {
    return(return_panels_list)
    
  } else  {
    for (j in (1:length(return_panels_list))) {
      ggplot2::ggsave(filename = paste0(name_file, j, ".", device_file),
                      plot = return_panels_list[j],
                      device = device_file,
                      scale = 1,
                      height= 4,
                      width = 5,
                      units= "in",
                      dpi = res_file)
    }
  }
}