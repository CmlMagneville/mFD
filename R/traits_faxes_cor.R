# Function to test relationship between traits and axes of a functional space
#
# Authors: Nicolas Loiseau and Sébastien Villéger
#
# ------------------------------------------------------------------------------
#' Compute relationship between all traits and all axes of the functional space.
#' 
#' For continuous trait a linear model is computed and r2 and p-value are
#' returned. For other types of traits, a Kruskal-Wallis test is computed and
#' eta2 statistics is returned.
#' Option allows to plot trait-axis relationships with scatterplot and boxplot
#' for continuous and non-continuous traits, respectively.
#'
#'@param sp_tr a \strong{dataframe} containing species as rows and traits as columns.
#'
#'@param sp_faxes_coord  a \strong{matrix} of species coordinates in a multidimensional
#'  functional space. Species coordinates have been retrieved thanks to
#'  \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#'@param tr_nm a \strong{vector} gathering the names of traits (as in \code{sp_tr}) to consider.
#'
#'@param faxes_nm a \strong{vector} gathering the names of PCoA axes (as in \code{sp_faxes_coord}) to
#' consider.
#'
#'@param plot a \strong{logical value} indicating whether plot illustrating relations
#' between trait and axes should be drawn. \strong{You can only plot
#' relationships for up to 10 traits and/or 10 axes}.
#'
#'@param name_file a \strong{character string} with name of file to save the plot
#' as a 300dpi jpeg file. Default is 'NULL' which means plot is displayed. If
#'  plot is FALSE, this input is ignored.
#'
#'@param color_signif a \strong{R color name or an hexadecimal code} referring to the
#' colour of points when relationships between the trait and the axis is
#' significant. Default is dark blue.
#'
#'@return a table with for each combination of trait and axis (rows),
#' name of test performed, and corresponding statistics and p-value.
#' If plot == TRUE a multi-panel figure with traits as columns and axes as rows.
#' When relationships between trait and axis is significant the points are
#'  colored, else they remain grayish.
#'
#'@examples
#' load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#' load(system.file("extdata", "sp_faxes_coord_fruits", package = "mFD"))
#' traits.faxes.cor(sp_tr, sp_faxes_coord, tr_nm = NULL, faxes_nm = NULL,
#'  plot = FALSE, name_file = NULL,
#'  color_signif = "darkblue")
#'
#'@export

traits.faxes.cor <- function(sp_tr, sp_faxes_coord,
                             tr_nm = NULL, faxes_nm = NULL,
                             plot = FALSE, name_file = NULL,
                             color_signif = "darkblue") {
  
  
  ## Check inputs ####
  
  if (any(is.na(sp_tr))) {
    stop("Error: The species*traits dataframe contains NA. Please check.")
  }
  
  if (any(is.na(sp_faxes_coord))) {
    stop("Error: The species*coordinates matrix contains NA. Please check.")
  }
  
  if (is.null(rownames(sp_tr))) {
    stop("Error: No row names provided in species*traits dataframe.
           Please add species names as row names.")
  }
  
  if (is.null(colnames(sp_tr))) {
    stop("Error: No column names provided in species*traits dataframe.
           Please add traits names as column names.")
  }
  
  if (is.null(rownames(sp_faxes_coord))) {
    stop("Error: No row names provided in species*coordinates matrix.
           Please add species names as row names.")
  }
  
  if (is.null(colnames(sp_faxes_coord))) {
    stop("Error: No column names provided in species*coordinates matrix.
           Please add axes names as column names.")
  }
  
  if (!identical( rownames(sp_tr), rownames(sp_faxes_coord) ) ) {
    stop(paste("Error: Mismatch between species names in 'sp_tr' and
                 'sp_faxes_coord' ."))
  }
  
  # if not provided, getting tr_nm and faxes_nm else checking:
  
  if (is.null(tr_nm)){
    tr_nm <- names(sp_tr)
  }
  else {
    if (any(! tr_nm %in% names(sp_tr))) {
      stop(paste("Error: Trait names should be as in 'sp_tr'"))
    }
  }
  
  if (is.null(faxes_nm)){
    faxes_nm <- colnames(sp_faxes_coord)
  } else {
    if (any(! faxes_nm %in% colnames(sp_faxes_coord))) {
      stop(paste("Error: axes names should be as in 'sp_tr'"))
    }
  }
  
  
  ## preparing ####
  
  # combinations trait*axes:
  res <- as.data.frame(matrix(NA, length(tr_nm)*length(faxes_nm), 6 ,
                              dimnames = list(NULL, c("trait", "axis",
                                                      "test", "stat",
                                                      "value", "p.value"))))
  
  # if needed, checking number of plots to draw:
  if (plot == TRUE) {
    
    # checking number of traits and axes did not exceed limits for graphics:
    if (length(faxes_nm) > 10) {
      stop(paste("Error: Number of axes to plot should be < 11"))
    }
    
    if (length(tr_nm) > 10) {
      stop(paste("Error: Number of traits to plot should be < 11"))
    }
  } # end of if plot
  
  
  # flag for moving down combinations of traits and axes:
  flag <- 0
  
  # loop on traits and axes ####
  for (i in tr_nm) {
    for (j in faxes_nm) {
      
      # moving through combinations of traits and axes:
      flag <- flag + 1
      
      # data:
      trait <- NULL
      axis <- NULL
      data_ij <- data.frame(trait = sp_tr[, i],
                            axis = sp_faxes_coord[, j])
      
      # if trait is continuous Linear Model ----
      if (is.numeric(data_ij$trait) == TRUE) {
        
        # test
        test_ij <- "Linear Model"
        
        lm_ij <- summary(stats::lm(data_ij[, "trait"] ~ data_ij[, "axis"]))
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, j, test_ij, "r2")
        res[flag, c("value", "p.value")] <- c(round(lm_ij$r.squared, 3),
                                              round(lm_ij$coefficients[2, 4], 4))
      }
      
      else {
        
        # if trait is not continuous Kruskal Wallis test  ----
        
        test_ij <- "Kruskal-Wallis"
        
        kw_ij <- stats::kruskal.test(axis ~ trait, data = data_ij)
        
        kw_ij_eta2 <- rstatix::kruskal_effsize(data = data_ij, axis ~ trait)
        
        res[flag, c("trait", "axis", "test", "stat")] <- c(i, j, test_ij, "eta2")
        res[flag, c("value", "p.value")] <- c(round(kw_ij_eta2$effsize, 3),
                                              round(kw_ij$p.value, 4))
      } # end if not numeric
      
      # plotting if needed ----
      
      if (plot == TRUE) {
        
        # X axis has trait names only for last Faxes = bottom row:
        if (j == faxes_nm[length(faxes_nm)]) {
          x_lab_ij <- i
        }
        else {
          x_lab_ij <- NULL
        }
        
        # y axis will have Faxes names only for first trait (left column):
        if (i == tr_nm[1]) {
          y_lab_ij <- j
        } else {
          y_lab_ij <- NULL
        }
        
        # main color according to significance of relationship:
        if (res[flag, c("p.value")] < 0.05) {
          col_cor <- color_signif
        }
        else {
          col_cor <- "gray90"
        }
        
        # empty plot:
        gg_ij <- ggplot2::ggplot(data_ij, ggplot2::aes(trait, axis)) +
          ggplot2::xlab(x_lab_ij) + ggplot2::ylab(y_lab_ij) +
          ggplot2::theme_bw()
        
        # if continuous trait, adding points:
        if (res[flag, c("stat")] == "r2") {
          gg_ij <- gg_ij + ggplot2::geom_point(size = 2, col = col_cor)
        }
        else {
          # if ordinal or categorical trait violin plot + jittered points:
          gg_ij <- gg_ij + ggplot2::geom_boxplot(colour = col_cor) +
            ggplot2::geom_jitter(colour = col_cor,
                                 width = 0.3, size = 1.5)
        }
        
        # if not first panel merging down with previous one(s):
        if (j == faxes_nm[1]) {
          col_j <- gg_ij
        }
        else {
          col_j <- col_j / gg_ij
        }
        
      } # end of if plot
      
      
    } # end of j
    
    # merging plots if needed:
    if(plot == TRUE) {
      # if not first trait, merging with other side by side:
      if (i == tr_nm[1]) {
        tr_faxes_plot <- col_j
      }
      else {
        tr_faxes_plot <- tr_faxes_plot | col_j
      }
    } # end of if plot
    
  } # end of i
  
  #### returning outputs ####
  
  # default = only table with statistics:
  if (plot == FALSE) {
    return(tr_faxes_stat = res)
  }
  else {
    # if plot to be returned:
    # adding title and legend to plot:
    tr_faxes_plot <- tr_faxes_plot + patchwork::plot_annotation(
      title = "Relation between Traits and PCoA axes",
      caption = "made with mFD package")
    
    # displaying or saving:
    if (is.null(name_file) == TRUE ){
      # table with statistics and plot:
      return(list(tr_faxes_stat = res, tr_faxes_plot = tr_faxes_plot))
    }
    else {
      # saving plot as a 300dpi jepg file:
      ggplot2::ggsave(filename = paste0(name_file, ".", "jpeg"),
                      plot = tr_faxes_plot,
                      device = "jpeg",
                      width = length(tr_nm)*3,
                      height = length(faxes_nm)*2.33,
                      units = "in",
                      dpi = 300 )
      
      # and returning table with statistics
      return(tr_faxes_stat = res)
      
    } # end of if jpeg
    
  } # end of if plot
  
} # end of function

