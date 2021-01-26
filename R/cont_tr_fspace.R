# Functions to build a functional space based on continuous traits: possibility
# to scale continuous traits(1st function) and build functional space (2nd function)
#
# Authors: Camille Magneville & Sébastien Villéger
#
#

# ------------------------------------------------------------------------------


#' Scale continuous traits
#'
#' Standardize continuous traits. Can be useful before computing functional
#' space. You will have to choose which standardized method to use based on your
#' data. For this function to work, there must be no NA in your \code{sp_tr}
#' dataframe.
#'
#' @param sp_tr a \strong{dataframe} of traits values for each species. Here traits must
#'   be \strong{continuous}
#'
#' @param std_method a \strong{character string} referring to the standardization method.
#'   Possible values: \emph{range} (standardize by the range),
#'   \emph{center} (use the center transformation: \eqn{x'= x - mean(x)}),
#'   \emph{scale} (use the scale transformation: \eqn{x' = \frac{x}{sd(x)}})
#'   or \emph{scale_center} (use the scale-center transformation: \eqn{x' =
#'   \frac{x - mean(x)}{sd(x)}}). Default is \emph{scale_center}.
#'
#' @return a dataframe of standardized trait values for each species \code{sp_tr}
#'
#' @examples
#' load(system.file("extdata", "sp_tr_cestes_df", package = "mFD"))
#' mFD::tr.cont.scale(sp_tr, std_method = "scale_center")
#'
#' @export


tr.cont.scale <- function(sp_tr, std_method = "scale_center") {
  
  if (is.null(colnames(sp_tr))) {
    stop("Error: No column names provided in traits table.
         Analysis will not go through, please add traits as column names")
  }
  if (is.null(rownames(sp_tr))) {
    stop("Error: No row names provided in traits table.
         Analysis will not go through, please add species names as row names")
  }
  if (is.null(sp_tr)) {
    stop("Error: NA in traits table. Analysis will not go through.")
  }
  
  # for range standardization:
  if (std_method == "range"){
    sp_tr <- apply(sp_tr, 2, function (x) (x - min(x))/(max(x) - min(x)))
  }
  # for center standardization:
  if (std_method == "center"){
    sp_tr <- apply(sp_tr, 2, function (x) x - mean(x))
  }
  # for scale standardization:
  if (std_method == "scale"){
    sp_tr <- apply(sp_tr, 2, function (x) x / stats::sd(x))
  }
  # for scale_center standardization:
  if (std_method == "scale_center"){
    sp_tr <- apply(sp_tr, 2, function (x) (x - mean(x) / stats::sd(x)))
  }
  return(sp_tr)
}


# ------------------------------------------------------------------------------


#' Build functional space
#'
#' Compute functional space based on continuous standardized traits or
#' continuous raw traits matrix. User can either choose to compute functional
#' space based on PCA analysis or using one trait for one functional axe. For
#' PCA analysis, center and scale arguments are considered FALSE: if you want to
#' center, scale or standardize by any mean your data, please use
#' \code{scale.conttr} function. Option makes it possible to compute correlation
#' between traits.
#'
#' @param sp_tr a \strong{dataframe} of raw/standardized traits values for each species.
#'   Here traits must be \strong{continuous}
#'
#' @param pca a \strong{logical value}, \code{TRUE} to compute PCA analysis, \code{FALSE}
#'   to compute functional space with one trait for each dimension. Default:
#'   TRUE
#'
#' @param nb_dim a \strong{numerical value} referring to the maximum number of dimensions
#'   for multidimensional functional spaces. Final number of dimensions depends
#'   on the number of positive eigenvalues obtained with PCA. High value for
#'   nb_dim can increase computation time. Default: nbdim=7.
#'   
#' @param scaling a \strong{string value} to compute or not scaling of traits
#'   using the \code{\link{tr.cont.scale}} function. Possible options are
#'   standardizing by the range \code{range}, center standardization
#'   \code{center}, scale standardization \code{scale}, scale center
#'   standardization \code{scale_center} or no scaling \code{no_scale}. Default : scale = 
#'
#' @param compute_corr a \strong{string value} to compute Pearson correlation
#'   coefficients between traits \code{"pearson"}. You can choose not to compute
#'   correlation coefficient by \code{"none"}
#'
#' @return a dataframe containing species coordinates on each functional axe
#'   \code{sp_faxes_coord}; species distance matrix \code{sp_dist_init};
#'   correlation coefficient between traits if asked.
#'
#' @examples
#' load(system.file("extdata", "sp_tr_cestes_df", package = "mFD"))
#' mFD::tr.cont.fspace(sp_tr, pca = TRUE, nb_dim = 7, scaling = "scale_center",
#'  compute_corr = "pearson")
#'
#' @export


tr.cont.fspace <- function(sp_tr, pca = TRUE, nb_dim = 7, scaling = "scale_center",
                           compute_corr = "pearson") {
  
  
  readline("If you want to standardize traits values, use scale.conttr
           function before fspace.conttr one \n(Press enter to continue)")
  
  if (any(is.na(sp_tr))) {
    stop("Error: There must be no NA in traits table.")
  }
  if (pca == TRUE) {
    if (ncol(sp_tr) < 3) {
      stop("Error: There must be at least 3 traits in 'sp_tr'.")
    }
    if (nrow(sp_tr) < 3) {
      stop("Error: There must be at least 3 species in 'sp_tr'.")
    }
  }
  if (nb_dim < 2) {
    stop("Error: Number of dimensions must be higher than 1.")
  }
  if (nb_dim > ncol(sp_tr)) {
    stop("Error: Number of dimensions must be lower than the number of traits.")
  }
  if (is.null(colnames(sp_tr))) {
    stop(
      "Error: No column names provided in traits table. Analysis will not go
      through, please add traits as column names"
    )
  }
  if (is.null(rownames(sp_tr))) {
    stop(
      "Error: No row names provided in traits table. Analysis will not go
      through, please add species names as row names"
    )
  }
  
  # function to compute correlation coefficient to be used later:
  compute_corr_coef <- function(sp_tr) {
    if (nrow(sp_tr) <= 4) {
      stop("Error: If you want to compute correlation matrix,
           you must have more than 4 observations.")
    }
    # convert species traits data.frame as a matrix for computing ...
    # ... traits correlation coefficients:
    sp_tr2 <- as.matrix(sp_tr)
    corr_tr <- Hmisc::rcorr(sp_tr2, type = "pearson")
  }
  
  # scale if needed: 
  if (scaling != c("no_scale")) {
    sp_tr <- mFD::tr.cont.scale(sp_tr, std_method = scaling)
  }

  if (pca == TRUE) {
    # compute functional dissimilarity matrix used for computing quality of...
    # ... functional spaces:
    sp_dist_init <- cluster::daisy(sp_tr, metric = "euclidean")
    
    # compute PCA analysis:
    pca_analysis <- FactoMineR::PCA(sp_tr, ncp = nb_dim, graph = FALSE, scale.unit = FALSE)
    sp_faxes_coord <- pca_analysis$ind$coord
    sp_faxes_coord <- as.data.frame(sp_faxes_coord)
    # restrict the number of column to nb_dim:
    if (ncol(sp_faxes_coord) > nb_dim){
      sp_faxes_coord <- dplyr::select(sp_faxes_coord, c(1:nb_dim))
    }
    
    # matrix to store quality results:
    quality_nbdim <- matrix(NA, nb_dim - 1, 2,
                            dimnames = list(paste0(2:nb_dim, "D"),
                                            c("mAD", "mSD")))
    
    # compute mSD and mAD based on deviation_dist function:
    deviation_dist <- function(sp_dist_init, sp_dist_multidim) {
      deviation_dist <- sp_dist_multidim - sp_dist_init
      mAD <- mean(abs(deviation_dist), 6)
      mSD <- mean(((deviation_dist)^2), 6)
      return(c(mAD = mAD, mSD = mSD))
    }
    sp_dist_multidim <- list()
    for (i in 2:nb_dim) {
      sp_dist_multidim2 <- stats::dist(sp_faxes_coord[, 1:i], method = "euclidean")
      quality_nbdim[paste0(i,"D"),
                    c("mAD", "mSD")] <- deviation_dist(sp_dist_init,
                                                       sp_dist_multidim2)
      sp_dist_multidim <- rlist::list.append(sp_dist_multidim,
                                             sp_dist_multidim2)
    }
    names(sp_dist_multidim) <- c(paste(2:nb_dim, "D", sep=""))
    
    # compute coeff correlation between traits:
    if (compute_corr == "pearson") {
      corr_tr_coeff <- compute_corr_coef(sp_tr)
      return_list1 <- list(quality_nbdim, as.matrix(sp_faxes_coord), sp_dist_multidim,
                           corr_tr_coeff)
      names(return_list1) <- c("mAD and mSD for each functional space",
                               "species coordinates in functional space",
                               "species distance in functional space",
                               "correlation coefficients between traits
                               and their associated pvalue")
      return(return_list1)
    }
    if (compute_corr == "none") {
      return_list1 <- list(quality_nbdim, sp_faxes_coord, sp_dist_multidim)
      names(return_list1) <- c("mAD and mSD for each functional space",
                               "species coordinates in functional space",
                               "species distance in functional space")
      return(return_list1)
    }
  }
  
  
  if (pca == FALSE) {
    # compute distance matrix between species for computing ...
    # ... multidimensional space:
    sp_dist_init <- cluster::daisy(sp_tr, metric = "euclidean")
    sp_faxes_coord <- sp_tr
    
    # compute coeff correlation between traits:
    if (compute_corr == "pearson") {
      corr_tr_coeff <- compute_corr_coef(sp_tr)
      return_list2 <- list(sp_faxes_coord, sp_dist_init, corr_tr_coeff)
      names(return_list2) <- c("species coordinates in functional space",
                               "species distances in functional space",
                               "correlation coefficients between traits
                               and their associated pvalue")
      return(return_list2)
    }
    if (compute_corr == "none") {
      return_list2 <- list(as.matrix(sp_faxes_coord), sp_dist_init)
      names(return_list2) <- c("species coordinates in functional space",
                               "species distances in functional space")
      return(return_list2)
    }
  }
}

