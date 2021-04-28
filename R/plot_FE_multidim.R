#' Illustrate Functional Diversity indices based on Functional Entities
#' 
#' Graphical representation of distribution of species in Functional
#' Entities (FE) and of indices from Mouillot _et al._ (2014). \strong{To plot
#' functional indices, functional indices values must have been computed first
#' through the use of the} \code{\link{alpha.fd.fe}} function.
#'
#' @param alpha_fd_fe output from the function \code{\link{alpha.fd.fe}} 
#'   applied on assemblage of interest with \code{details_returned = TRUE}.
#'
#' @param plot_asb_nm a vector containing the name of the assemblage to plot.
#'
#' @param name_file a character string with name of file to save the figure
#'   (without extension). Default: `name_file = NULL` which means plot is
#'   displayed.
#'
#' @param plot_ind_nm a vector containing the names of the indices to plot. It
#'  can be \code{'fred'} to plot functional redundancy (FRed), \code{'fored'} to
#'  plot functional over-redundancy (FOred) and/or \code{'fvuln'} to plot
#'  functional vulnerability (FVuln). Default is all `3` indices.
#'
#' @param color_fill_fored a R color name or an hexadecimal code referring to
#'   the color used to fill the part of barplots that contain species in excess
#'   in species-rich FEs. It refers to the FORed value. Default:
#'   `color_fill_fored = "darkolivegreen2"`.
#'
#' @param color_line_fred a R color name or an hexadecimal code referring to 
#'   the color used to draw the horizontal line referring to the FRed value.
#'   Default: `color_line_fred = "darkolivegreen4"`.
#'
#' @param color_fill_bar a R color name or an hexadecimal code referring to the
#'  color used to draw barplots. Default: `color_fill_bar = "grey80"`.
#'
#' @param color_fill_fvuln a R color name or an hexadecimal code referring to
#'   the color used to fill barplot containing only one species for 
#'   illustrating FVuln. Default: `color_fill_fvuln = "lightcoral"`.
#'
#' @param color_arrow_fvuln  a R color name or an hexadecimal code referring to
#'  the color used to draw the horizontal arrow showing the proportion of FEs
#'  containing only one species for illustrating FVuln. If there is only
#'  one FE containing one species, the arrow will be a point. Default:
#'  `color_arrow_fvuln = "indianred4"`.
#'
#' @param size_line_fred a numeric value referring to the size of the 
#'   horizontal line illustrating FRed. Default: `size_line_fred = 1.5.`
#'
#' @param size_arrow_fvuln a numeric value referring to the size of the arrow
#'  showing the proportion of FEs containing only one species. Default:
#'  `size_arrow_fvuln = 1`.
#'
#' @param check_input a logical value indicating whether key features the 
#'   inputs are checked (e.g. class and/or mode of objects, names of rows 
#'   and/or columns, missing values). If an error is detected, a detailed 
#'   message is returned. Default: `check_input = TRUE`.
#'
#' @return A \code{patchwork} object with a barplot of number of species per 
#'   FE. Indices names provided in 'plot_ind_nm' are illustrated. Functional
#'   Redundancy (average number of species per FE) is illustrated with a
#'   horizontal line. Functional Over-redundancy (proportion of species in 
#'   excess in FE richer than average) is illustrated with top part of these 
#'   bars filled with 'color_fill_fored'. Functional Vulnerability (proportion 
#'   of FE with a single species) is illustrated with bars of these vulnerable 
#'   FE filled with 'color_fill_fvuln' and the double-head arrow highlighting 
#'   their number. FE-based indices values on top of the plot. if 
#'   \code{name_file} is provided, plot saved as a 300dpi png file in the 
#'   working directory.
#'
#' @references
#' Mouillot _et al._ (2014) Functional over-redundancy and high functional 
#'   vulnerability in global fish faunas on tropical reefs. _PNAS_, **111**, 
#'   13757-13762.
#'
#' @author Camille Magneville and Sebastien Villeger
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load Species*Traits dataframe
#' data("fruits_traits", package = "mFD")
#' 
#' # Load Traits categories dataframe
#' data("fruits_traits_cat", package = "mFD")
#' 
#' # Load Assemblages*Species matrix
#' data("baskets_fruits_weights", package = "mFD")
#' 
#' # Remove continuous trait
#' fruits_traits  <- fruits_traits[ , -5]
#' fruits_traits_cat <- fruits_traits_cat[-5, ]
#' 
#' # Compute gathering species into FEs
#' sp_to_fe_fruits <- mFD::sp.to.fe(
#'   sp_tr       = fruits_traits, 
#'   tr_cat      = fruits_traits_cat, 
#'   fe_nm_type  = "fe_rank", 
#'   check_input = TRUE)
#'  
#' # Get the occurrence matrix
#' asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights) 
#' asb_sp_fruits_occ <- asb_sp_fruits_summ$"asb_sp_occ"
#' 
#' # Compute alpha fd indices
#' alpha_fd_fe_fruits <- mFD::alpha.fd.fe(
#'   asb_sp_occ       = asb_sp_fruits_occ, 
#'   sp_to_fe         = sp_to_fe_fruits,
#'   ind_nm           = c("fred", "fored", "fvuln"),
#'   check_input      = TRUE, 
#'   details_returned = TRUE)
#'   
#' # Plot fd fe indices
#' mFD::alpha.fd.fe.plot(
#'   alpha_fd_fe       = alpha_fd_fe_fruits, 
#'   plot_asb_nm       = c("basket_1"), 
#'   plot_ind_nm       = c("fred", "fored", "fvuln"),
#'   name_file         = NULL,
#'   color_fill_fored  = "darkolivegreen2",
#'   color_line_fred   = "darkolivegreen4",
#'   color_fill_bar    = "grey80",
#'   color_fill_fvuln  = "lightcoral",
#'   color_arrow_fvuln = "indianred4",
#'   size_line_fred    = 1.5,
#'   size_arrow_fvuln  = 1,
#'   check_input       = TRUE) 
#' }

alpha.fd.fe.plot <- function(alpha_fd_fe,
                             plot_asb_nm,
                             plot_ind_nm = c("fred", "fored", "fvuln"),
                             name_file = NULL,
                             color_fill_fored = "darkolivegreen2",
                             color_line_fred = "darkolivegreen4",
                             color_fill_bar = "grey80",
                             color_fill_fvuln = "lightcoral",
                             color_arrow_fvuln = "indianred4",
                             size_line_fred = 1.5,
                             size_arrow_fvuln = 1,
                             check_input = TRUE) {
  
  
  # if asked to check_inputs:
  if (check_input) {
    
    if (any(!names(alpha_fd_fe) %in% c("asb_fdfe", "details_fdfe"))) {
      stop("Input 'alpha_fd_fe' should be output from function 'alpha.fd.fe' ",
           "applied with details_returned = TRUE. Please check.")
    }
    
    
    # check name(s) of indice(s) to plot:
    if (!any(plot_ind_nm %in% c("fred", "fored", "fvuln"))) {
      stop("Names of functional diversity indices to plot are not within ", 
           "names allowed. Be careful, they should all be written in ",
           "lowercase letters.")
    }
    
    # check indice(s) were computed:
    if (!any(plot_ind_nm %in% colnames(alpha_fd_fe$asb_fdfe))) {
      stop("Names of functional diversity indices to plot should have been ",
           "computed with 'alpha.fd.fe'. Please check.")
    }
    
    # check name of the assemblage to plot:
    if (!plot_asb_nm %in% rownames(alpha_fd_fe$asb_fdfe)){
      stop("Name of assemblage for which FE-based indices should be plotted ",
           "is not  present in input 'alpha_fd_fe'. Please check function ",
           "'alpha.fd.fe' was applied on this assemblage.")
    }
  }
  
  
  # number of species per FE present in assemblage to plot
  fe_nbsp_k <- alpha_fd_fe$details_fdfe$asb_fe_nbsp[plot_asb_nm, ]
  fe_nbsp_k <- fe_nbsp_k[which(fe_nbsp_k > 0)]
  
  # FE sorted according to decreasing number of species
  fe_nbsp_k <- sort(fe_nbsp_k, decreasing = TRUE)
  
  # names of FE after sorting according to nb of species
  fe_nm   <- NULL
  fe_nbsp <- NULL
  data_k  <- data.frame(fe_nm = names(fe_nbsp_k), fe_nbsp = fe_nbsp_k)
  
  breaks_max <- max(fe_nbsp_k)
  breaks <- seq(0, breaks_max, by = 1)
  
  # plotting number of species for all FE
  plot_k <- ggplot2::ggplot(data = data_k, 
                            ggplot2::aes(
                            x = stats::reorder(fe_nm, -fe_nbsp), 
                            y = fe_nbsp)) +
    ggplot2::geom_bar(stat = "identity", fill = color_fill_bar) +
    ggplot2::scale_y_continuous(breaks = breaks) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   panel.background = ggplot2::element_rect(fill = "white", 
                                                            colour = "grey50", 
                                                            size = 1),
                   panel.grid.major = ggplot2::element_line(size = 0.5,
                                                            linetype = 'solid',
                                                            colour = "grey85"),
                   panel.grid.minor = ggplot2::element_line(size = 0.5, 
                                                          linetype = 'solid', 
                                                          colour = "grey85")) +
    ggplot2::xlab("Functional Entities") +
    ggplot2::ylab("Number of species")
  
  # title and subtitle
  tit_k <- paste0("FD based on FE for ", plot_asb_nm, " : ")
  sub_k <- ""
  
  # adding illustration of funct vulnerability if needed
  if ("fvuln" %in% plot_ind_nm) {
    
    # names of vulnerable FE (sorted alphabetical order)
    fe_vuln_k<-sort(data_k[which(data_k$fe_nbsp==1), "fe_nm"])
    
    if (length(fe_vuln_k) > 0) {
      plot_k <- plot_k +
        ggplot2::geom_bar(data = data_k[fe_vuln_k, ], 
                          ggplot2::aes(x = fe_nm, y = fe_nbsp),
                          stat = "identity", fill = color_fill_fvuln)
      if (length(fe_vuln_k) > 1) {
        plot_k <- plot_k +
          ggplot2::geom_segment(ggplot2::aes(x = fe_vuln_k[1] , y = 0.5,
                                            xend = fe_vuln_k[
                                              length(fe_vuln_k)], 
                                            yend = 0.5),
                                arrow = grid::arrow(ends = "both", 
                                                    length = grid::unit(0.25, 
                                                                        "cm")),
                                color = color_arrow_fvuln, 
                                size = size_arrow_fvuln)
      }# end of if >1 vulnerable
    } # end of if >=1 vulnerable FE
    sub_k<-paste0(sub_k, "  FVuln=", round(alpha_fd_fe$asb_fdfe[plot_asb_nm, 
                                                                'fvuln'], 3))
  }
  
  
  # adding illustration of funct over-redundancy if needed
  if ("fored" %in% plot_ind_nm) {
    fred_k <- alpha_fd_fe$asb_fdfe[plot_asb_nm, 'fred']
    
    # names of over-redundant FE (sorted alphabetical order)
    fe_ored_k <- sort(data_k[which(data_k$fe_nbsp > fred_k), "fe_nm"])
    
    if (length(fe_ored_k)>0) {
      plot_k <- plot_k +
        ggplot2::geom_bar(data = data_k[fe_ored_k, ], 
                          ggplot2::aes(x = fe_nm, y = fe_nbsp),
                          stat = "identity", fill = color_fill_fored) +
        ggplot2::geom_bar(data = data_k[fe_ored_k, ], ggplot2::aes(x = fe_nm, 
                                                                   y = fred_k),
                          stat = "identity", fill = color_fill_bar)
    }#edn of if FE over-redundant
    sub_k <- paste0("  FORed=", round(alpha_fd_fe$asb_fdfe[plot_asb_nm, 
                                                           'fored'], 3), sub_k)
  }
  
  # adding illustration of funct redundancy if needed
  if ("fred" %in% plot_ind_nm) {
    plot_k <- plot_k +
      ggplot2::geom_hline(yintercept = alpha_fd_fe$asb_fdfe[plot_asb_nm, 
                                                            'fred'],
                          color = color_line_fred, size = size_line_fred)
    sub_k <- paste0("  FRed=", round(alpha_fd_fe$asb_fdfe[plot_asb_nm, 'fred'], 
                                     3), sub_k)
  }
  
  # create patchwork object to return:
  return_plot <- plot_k +
    patchwork::plot_annotation(title = tit_k,
                               subtitle = sub_k,
                               caption = "made with mFD package")
  
  ## returning output ####
  
  # type, resolution and dimensions of file if to be saved
  device_file <- "png"
  res_file <- 300
  
  # displaying or saving
  if (is.null(name_file)) {
    
    return(return_plot)
    
  } else { 
    
    ggplot2::ggsave(filename = paste0(name_file, ".", device_file),
                    plot = return_plot,
                    device = device_file,
                    scale = 1,
                    height = 4,
                    width = 5,
                    units = "in",
                    dpi = res_file)
  }
}
