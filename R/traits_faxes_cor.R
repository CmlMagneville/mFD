#' Correlation between Traits and Axes
#' 
#' Compute relationship between all traits and all axes of the functional space. 
#' For continuous trait a linear model is computed and r2 and p-value are 
#' returned. For other types of traits, a Kruskal-Wallis test is computed and  
#' eta2 statistics is returned.
#' 
#' Option allows to plot trait-axis relationships with scatterplot and boxplot 
#' for continuous and non-continuous traits, respectively.
#'
#' @param sp_tr a data frame containing species as rows and traits as columns.
#'
#' @param sp_faxes_coord a matrix of species coordinates in a multidimensional
#'   functional space. Species coordinates have been retrieved 
#'   thanks to \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param tr_nm a vector gathering the names of traits (as in `sp_tr`) to 
#'   consider. If `NULL` all traits are considered.
#'
#' @param faxes_nm a vector gathering the names of PCoA axes (as in 
#'   `sp_faxes_coord`) to consider.
#'
#' @param plot a logical value indicating whether plot illustrating relations
#'   between trait and axes should be drawn. **You can only plot relationships 
#'   for up to 10 traits and/or 10 axes**.
#'
#' @param name_file the file name (without extension) to save the plot as a 300 
#'   dpi JPEG file. Default is `NULL` which means plot is only displayed. If 
#'   `plot = FALSE` this argument is ignored.
#'
#' @param color_signif an R color name or an hexadecimal code referring to 
#'   the color of points when relationships between the trait and the axis is
#'   significant. Default is `dark blue`.
#'
#' @return A data frame with for each combination of trait and axis (rows), the  
#'   name of the test performed, and the corresponding statistics and p-value. 
#'   If `plot = TRUE` a multi-panel figure with traits as columns and axes as 
#'   rows is also plotted. When relationships between trait and axis is 
#'   significant the points are colored, else they remain grayish.
#'
#' @author Nicolas Loiseau & Sebastien Villeger
#'
#' @export
#' 
#' @importFrom ggplot2 ggplot aes xlab ylab theme_bw geom_point geom_boxplot 
#'   geom_jitter ggsave
#' @importFrom patchwork plot_annotation
#' @importFrom rstatix kruskal_effsize
#' @importFrom stats lm summary.lm kruskal.test
#' 
#' @examples
#' # Load Species x Traits Data
#' data("sp_tr_fruits", package = "mFD")
#' 
#' # Load Traits categories dataframe
#' data("tr_cat_fruits", package = "mFD")
#' 
#' # Compute Functional Distance
#' sp_dist_fruits <- mFD::funct.dist(
#'   sp_tr       = sp_tr_fruits,         
#'   tr_cat      = tr_cat_fruits,   
#'   dist_metric = "kgower",         
#'   scaling     = "scaledBYrange",  
#'   stop_if_NA  = TRUE)
#'   
#' # Compute Functional Spaces Quality (to retrieve species coordinates)
#' fspaces_quality_fruits <- mFD::quality.fspaces(
#'   sp_dist             = sp_dist_fruits, 
#'   maxdim_pcoa         = 10,
#'   deviation_weighting = "absolute",
#'   fdist_scaling       = FALSE,
#'   fdendro             = "average")
#'   
#' # Retrieve Species Coordinates
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#' 
#' # Compute Correlation between Traits and Functional Axes
#' mFD::traits.faxes.cor(
#'   sp_tr          = sp_tr_fruits, 
#'   sp_faxes_coord = sp_faxes_coord_fruits, 
#'   tr_nm          = NULL, 
#'   faxes_nm       = NULL,
#'   name_file      = NULL, 
#'   color_signif   = "darkblue")


traits.faxes.cor <- function(sp_tr, sp_faxes_coord, tr_nm = NULL, 
                             faxes_nm = NULL, plot = FALSE, name_file = NULL,
                             color_signif = "darkblue") {
  
  
  ## Check inputs ----
  
  if (missing(sp_tr)) {
    stop("Argument 'sp_tr' is mandatory.")
  }
  
  if (missing(sp_faxes_coord)) {
    stop("Argument 'sp_faxes_coord' is mandatory.")
  }
  
  check.sp.tr(sp_tr)
  
  check.sp.faxes.coord(sp_faxes_coord)

  
  if (!identical(sort(rownames(sp_tr)), sort(rownames(sp_faxes_coord)))) {
    stop("Species names mismatch between 'sp_tr' and 'sp_faxes_coord'.")
  }
  
  # If not provided, getting 'tr_nm' and 'faxes_nm' else checking
  
  if (is.null(tr_nm)) {
    tr_nm <- names(sp_tr)
  } else {
    if (any(!(tr_nm %in% names(sp_tr)))) {
      stop("Trait names should be as in 'sp_tr'.")
    }
  }
  
  if (is.null(faxes_nm)){
    faxes_nm <- colnames(sp_faxes_coord)
  } else {
    if (any(!(faxes_nm %in% colnames(sp_faxes_coord)))) {
      stop("Axes names should be as in 'sp_tr'.")
    }
  }
  
  
  # If needed, checking number of plots to draw
  if (plot) {
    
    if (length(faxes_nm) > 10) {
      stop("Number of axes to plot should be < 11.")
    }
    
    if (length(tr_nm) > 10) {
      stop("Number of traits to plot should be < 11.")
    }
  }
  
  
  ## Preparing ----
  
  # combinations trait*axes:
  res <- as.data.frame(matrix(NA, length(tr_nm) * length(faxes_nm), 6 ,
                              dimnames = list(NULL, c("trait", "axis", "test", 
                                                      "stat", "value", "p.value"
                                                      ))))
  
  
  # flag for moving down combinations of traits and axes:
  flag <- 0
  
  # loop on traits and axes ####
  for (i in tr_nm) {
    for (j in faxes_nm) {
      
      # moving through combinations of traits and axes:
      flag <- flag + 1
      
      # data:
      trait   <- NULL
      axis    <- NULL
      data_ij <- data.frame("trait" = sp_tr[ , i], 
                            "axis"  = sp_faxes_coord[ , j])
      
      # if trait is continuous Linear Model ----
      if (is.numeric(data_ij$trait)) {
        
        # test
        test_ij <- "Linear Model"
        
        lm_ij <- summary(stats::lm(trait ~ axis, data = data_ij))
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, j, test_ij, "r2")
        res[flag, c("value", "p.value")] <- c(round(lm_ij$r.squared, 3),
                                              round(lm_ij$coefficients[2, 4], 
                                                    4))
      } else {
        
        # if trait is not continuous Kruskal Wallis test  ----
        
        test_ij <- "Kruskal-Wallis"
        
        kw_ij <- stats::kruskal.test(axis ~ trait, data = data_ij)
        
        kw_ij_eta2 <- rstatix::kruskal_effsize(data = data_ij, axis ~ trait)
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, j, test_ij, 
                                                           "eta2")
        res[flag, c("value", "p.value")] <- c(round(kw_ij_eta2$effsize, 3),
                                              round(kw_ij$p.value, 4))
      }
      
      # plotting if needed ----
      
      if (plot) {
        
        # X axis has trait names only for last Faxes (bottom row)
        x_lab_ij <- NULL
        if (j == faxes_nm[length(faxes_nm)]) {
          x_lab_ij <- i
        }
        
        # Y axis will have Faxes names only for first trait (left column)
        y_lab_ij <- NULL
        if (i == tr_nm[1]) {
          y_lab_ij <- j
        }
        
        # Main color according to significance of relationship:
        if (res[flag, "p.value"] < 0.05) {
          col_cor <- color_signif
        } else {
          col_cor <- "gray90"
        }
        
        # Empty plot
        gg_ij <- ggplot2::ggplot(data_ij, ggplot2::aes(trait, axis)) +
          ggplot2::xlab(x_lab_ij) + 
          ggplot2::ylab(y_lab_ij) +
          ggplot2::theme_bw()
        
        # if continuous trait, adding points:
        if (res[flag, "stat"] == "r2") {
          
          gg_ij <- gg_ij + ggplot2::geom_point(size = 2, col = col_cor)
          
        } else {
          
          # if ordinal or categorical trait violin plot + jittered points:
          gg_ij <- gg_ij + ggplot2::geom_boxplot(colour = col_cor) +
            ggplot2::geom_jitter(colour = col_cor, width = 0.3, size = 1.5)
        }
        
        # if not first panel merging down with previous one(s):
        if (j == faxes_nm[1]) {
          col_j <- gg_ij
        } else {
          col_j <- col_j / gg_ij
        }
      } # end of if plot
    } # end of j
    
    # merging plots if needed:
    if (plot) {
      # if not first trait, merging with other side by side:
      if (i == tr_nm[1]) {
        tr_faxes_plot <- col_j
      } else {
        tr_faxes_plot <- tr_faxes_plot | col_j
      }
    } # end of if plot
  } # end of i
  
  
  ## Function Return ----
  
  if (!plot) {
    
    return("tr_faxes_stat" = res)    # default
  
  } else {
    
    # if plot to be returned:
    # adding title and legend to plot:
    tr_faxes_plot <- tr_faxes_plot + 
      patchwork::plot_annotation(title = "Relation between traits and PCoA axes",
                                 caption = "Made with mFD package")
    
    # displaying or saving:
    if (is.null(name_file)){
      
      # table with statistics and plot:
      return(list("tr_faxes_stat" = res, "tr_faxes_plot" = tr_faxes_plot))
    
    } else {
      
      # saving plot as a 300dpi jepg file:
      ggplot2::ggsave(filename = paste0(name_file, ".jpeg"),
                      plot     = tr_faxes_plot,
                      device   = "jpeg",
                      width    = length(tr_nm) * 3,
                      height   = length(faxes_nm) * 2.33,
                      units    = "in",
                      dpi      = 300)
      
      # and returning table with statistics
      return("tr_faxes_stat" = res)
      
    } # end of if jpeg
  } # end of if plot
} # end of function
