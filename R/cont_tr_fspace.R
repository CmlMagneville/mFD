#' Scale Continuous Traits
#'
#' This function standardizes continuous traits. It can be useful before 
#' computing functional space. You will have to choose which standardized method 
#' to use based on your data. For this function to work, there must be no NA in 
#' your `sp_tr` data frame.
#'
#' @param sp_tr a data frame of traits values (columns) for each species (rows).
#'   Note that species names **must be** specified in the row names and traits 
#'   must be **continuous**.
#'
#' @param std_method a character string referring to the standardization method.
#'   Possible values: 
#'   `range` (standardize by the range), 
#'   `center` (use the center transformation: \eqn{x' = x - mean(x)}), 
#'   `scale` (use the scale transformation: \eqn{x' = \frac{x}{sd(x)}}), or
#'   `scale_center` (use the scale-center transformation: 
#'   \eqn{x' = \frac{x - mean(x)}{sd(x)}}). 
#'   Default is `scale_center`.
#'
#' @return A data frame of standardized trait values (columns) for each species
#'   (rows).
#'
#' @author Camille Magneville & Sebastien Villeger
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' load(system.file('extdata', 'sp_tr_cestes_df', package = 'mFD'))
#' 
#' mFD::tr.cont.scale(sp_tr = sp_tr, std_method = 'scale_center')
#' }

tr.cont.scale <- function(sp_tr, std_method = "scale_center") {
  
  ## Check Inputs ----
  
  check.sp.tr(sp_tr)
  
  if (any(!apply(sp_tr, 2, is.numeric))) {
    stop("Species x traits data frame must contain only numerical variables.")
  }
  
  std_method <- match.arg(std_method, c("range", "center", "scale", 
                                        "scale_center"))
  
  
  ## Standardization ----
  
  if (std_method == "range") {
    sp_tr <- apply(sp_tr, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  }
  
  if (std_method == "center") {
    sp_tr <- apply(sp_tr, 2, function(x) x - mean(x))
  }
  
  if (std_method == "scale") {
    sp_tr <- apply(sp_tr, 2, function(x) x / stats::sd(x))
  }
  
  if (std_method == "scale_center") {
    sp_tr <- apply(sp_tr, 2, function(x) (x - mean(x) / stats::sd(x)))
  }
  
  return(sp_tr)
}



#' Build a Functional Space based on Continuous Traits Only
#'
#' This function computes a functional space based on continuous standardized 
#' traits or continuous raw traits matrix. User can either choose to compute 
#' functional space based on PCA analysis or using one trait for one functional 
#' axis. For PCA analysis, center and scale arguments are considered `FALSE`: if 
#' you want to center, scale or standardize by any mean your data, please use
#' \code{\link{tr.cont.scale}} function. Option makes it possible to compute 
#' correlation between traits.
#'
#' @param sp_tr a data frame of traits values (columns) for each species (rows).
#'   Note that species names **must be** specified in the row names and traits 
#'   must be **continuous** (raw or standardized).
#'
#' @param pca a logical value. If `TRUE` a PCA analysis is computed, elsewhere
#'   the functional space is computed with one trait for each dimension. 
#'   Default is `TRUE`.
#'
#' @param nb_dim an integer referring to the maximum number of dimensions for
#'   multidimensional functional spaces. Final number of dimensions depends
#'   on the number of positive eigenvalues obtained with the PCA. High value for
#'   `nb_dim` can increase computation time. Default is `nb_dim = 7`.
#'   
#' @param scaling a string value to compute (or not) scaling of traits using the
#'   \code{\link{tr.cont.scale}} function. Possible options are:
#'   `range` (standardize by the range), 
#'   `center` (use the center transformation: \eqn{x' = x - mean(x)}), 
#'   `scale` (use the scale transformation: \eqn{x' = \frac{x}{sd(x)}}),
#'   `scale_center` (use the scale-center transformation: 
#'   \eqn{x' = \frac{x - mean(x)}{sd(x)}}), or
#'   `no_scale`
#'   Default is `scale_center`.
#'
#' @param compute_corr a string value to compute Pearson correlation
#'   coefficients between traits (`compute_corr = 'pearson'`). You can choose 
#'   not to compute correlation coefficient by setting `compute_corr` to `none`.
#'
#' @return A list with a data frame containing species coordinates on each 
#'   functional axis, a species distance matrix, and a correlation coefficients 
#'   between traits (if asked).
#'
#' @author Camille Magneville & Sebastien Villeger
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' load(system.file('extdata', 'sp_tr_cestes_df', package = 'mFD'))
#' 
#' mFD::tr.cont.fspace(
#'     sp_tr        = sp_tr, 
#'     pca          = TRUE, 
#'     nb_dim       = 7, 
#'     scaling      = 'scale_center',
#'     compute_corr = 'pearson')
#' }


tr.cont.fspace <- function(sp_tr, pca = TRUE, nb_dim = 7, 
                           scaling = "scale_center", compute_corr = "pearson") {
  
  
  ## Check Inputs ----
  
  check.sp.tr(sp_tr)
  
  if (any(!apply(sp_tr, 2, is.numeric))) {
    stop("Species x traits data frame must contain only numerical variables.")
  }
  
  scaling <- match.arg(scaling, c("range", "center", 
                                  "scale", "scale_center", "no_scale"))
  
  compute_corr <- match.arg(compute_corr, c("pearson", 
                                            "none"))
  
  
  if (pca) {
    
    if (ncol(sp_tr) < 3) {
      stop("There must be at least 3 traits in 'sp_tr'.")
    }
    
    if (nrow(sp_tr) < 3) {
      stop("There must be at least 3 species in 'sp_tr'.")
    }
  }
  
  if (nb_dim < 2) {
    stop("Number of dimensions must be higher than 1.")
  }
  
  if (nb_dim > ncol(sp_tr)) {
    stop("Number of dimensions must be lower than the number of traits.")
  }
  
  
  ## Functions Definition ----
  
  compute_corr_coef <- function(sp_tr) {
    # Compute Correlation Matrix
    
    if (nrow(sp_tr) <= 4) {
      stop("If you want to compute correlation matrix you must have more than ", 
           "4 observations.")
    }
    
    Hmisc::rcorr(as.matrix(sp_tr), type = "pearson")
  }
  
  deviation_dist <- function(sp_dist_init, sp_dist_multidim) {
    # mSD and mAD
    
    deviation_dist <- sp_dist_multidim - sp_dist_init
    
    c(mAD = mean(abs(deviation_dist), 6), mSD = mean(((deviation_dist)^2), 
                                                     6))
  }
  
  
  ## Standardize Traits ----
  
  if (scaling != "no_scale") {
    sp_tr <- tr.cont.scale(sp_tr, std_method = scaling)
  }
  
  
  if (pca) {
    
    # compute functional dissimilarity matrix used for
    # computing quality of...  ... functional spaces:
    sp_dist_init <- cluster::daisy(sp_tr, metric = "euclidean")
    
    # compute PCA analysis:
    pca_analysis <- FactoMineR::PCA(sp_tr, ncp = nb_dim, 
                                    graph = FALSE, scale.unit = FALSE)
    
    sp_faxes_coord <- as.data.frame(pca_analysis$ind$coord)
    
    # restrict the number of column to nb_dim:
    if (ncol(sp_faxes_coord) > nb_dim) {
      sp_faxes_coord <- sp_faxes_coord[, 1:nb_dim]
    }
    
    
    # matrix to store quality results:
    quality_nbdim <- matrix(NA, nb_dim - 1, 2, 
                            dimnames = list(paste0(2:nb_dim, "D"), 
                                            c("mAD", "mSD")))
    
    sp_dist_multidim <- list()
    
    for (i in 2:nb_dim) {
      
      sp_dist_multidim2 <- stats::dist(sp_faxes_coord[, 
                                                      1:i], 
                                       method = "euclidean")
      quality_nbdim[paste0(i, "D"), c("mAD", 
                                      "mSD")] <- deviation_dist(sp_dist_init, 
                                                              sp_dist_multidim2)
      
      sp_dist_multidim[[i - 1]] <- sp_dist_multidim2
      names(sp_dist_multidim)[i - 1] <- paste0(i, "D")
    }
    
    
    if (compute_corr == "pearson") {
      
      corr_tr_coeff <- compute_corr_coef(sp_tr)
      
      return_list1 <- list(quality_nbdim, as.matrix(sp_faxes_coord), 
                           sp_dist_multidim, corr_tr_coeff)
      
      names(return_list1) <- c("quality_metrics", 
                               "sp_faxes_coord", 
                               "sp_dist", 
                               "tr_correl")
      return(return_list1)
      
    } else {
      
      # no Pearson correlation
      
      return_list1 <- list(quality_nbdim, sp_faxes_coord, 
                           sp_dist_multidim)
      
      names(return_list1) <- c("quality_metrics", 
                               "sp_faxes_coord", 
                               "sp_dist")
      return(return_list1)
    }
    
  } else {
    # no PCA
    
    # compute distance matrix between species for
    # computing ...  ... multidimensional space:
    sp_dist_init <- cluster::daisy(sp_tr, metric = "euclidean")
    sp_faxes_coord <- sp_tr
    
    if (compute_corr == "pearson") {
      
      corr_tr_coeff <- compute_corr_coef(sp_tr)
      return_list2 <- list(sp_faxes_coord, sp_dist_init, 
                           corr_tr_coeff)
      
      names(return_list2) <- c("sp_faxes_coord", 
                               "sp_dist", 
                               paste("tr_correl"))
      return(return_list2)
      
    } else {
      # no Pearson correlation
      
      return_list2 <- list(as.matrix(sp_faxes_coord), 
                           sp_dist_init)
      
      names(return_list2) <- c("sp_faxes_coord", 
                               "sp_dist")
      return(return_list2)
    }
  }
}
