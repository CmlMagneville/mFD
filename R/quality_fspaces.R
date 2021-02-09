# Function to compute the quality of multidimensional functional spaces
#
# Authors: Sébastien Villéger, Eva Maire & Camille Magneville
#
#------------------------------------------------------------------------------


#' Compute functional spaces and return their quality
#'
#' Compute a Principal Coordinates Analysis (PCoA) using functional distance
#' between species. Then the function evaluates the quality of spaces
#' built using an increasing number of Principal components. Quality is evaluated
#' as the (absolute or squared) deviation between trait-based distance (input)
#' and distance in the PCoA-based space (raw Euclidean distance or scaled
#' distance according to its maximum value and maximum of trait-based distance).
#' Option to compute a functional dendrogram and its quality.
#' This function is based on the framework presented in \emph{Maire et al 2015},
#' DOI:10.1111/geb.12299).
#'
#' @param sp_dist a \strong{dist object} with pairwise distance among all species (at
#'   least 3 species needed). Functional distance matrix from trait values can
#'   be computed using \code{\link{funct.dist}} function.
#'
#' @param maxdim_pcoa a single \strong{numeric value} with maximum number of PCoA axes to
#'   consider to build multidimensional functional spaces. Default: maxdim_pcoa
#'   = 10. See note about number of axes actually considered.
#'
#' @param deviation_weighting a \strong{character string} referring to the method(s)
#'   used to weight the differences between species pairwise distance in the
#'   functional space and trait-based distance. \code{'absolute'} (default) means
#'   absolute differences are used to compute mean absolute deviation \emph{mad}
#'   index; \code{'squared'} means squared differences are used to compute root of
#'   mean squared deviation \emph{rmsd} index.
#'   Both values could be provided to compare quality metrics.
#'
#' @param fdist_scaling a \strong{vector} with logical value(s) specifying whether
#'   distances in the functional space should be scaled before computing
#'   differences with trait-based distances. Scaling ensures that trait-based
#'   distances and distances in the functional space have the same maximum.
#'   Default: FALSE. Both values could be provided to compare quality metrics.
#'
#' @param fdendro a \strong{character string} indicating the clustering algorithm to use
#'   to compute dendrogram. Should be one of the method recognized by
#'   \code{\link[stats]{hclust}} (e.g. 'average' for UPGMA). Default: fdendro = NULL (so
#'   no dendrogram computed).
#'
#' @return a list with: \itemize{
#'
#'   \item \code{$quality_fspaces}: a dataframe with quality metric(s) for each
#'   functional space. Functional spaces are named as 'pcoa_.d'
#'   and if required 'tree_clustering method'. Quality metrics are named after
#'   deviation_weighting ('mad' for 'absolute' and and 'rmsd' for 'squared')
#'   and if fdist_scaling is TRUE with suffix '_scaled'.
#'
#'   \item \code{$details_trdist} a list with 2 elements:
#'   \code{$trdist_summary} a vector with minimum (min), maximum (max),
#'   mean (mean) and standard deviation (sd) of \code{sp_dist} ;
#'   \code{$trdist_euclidean} a logical value indicating whether
#'   \code{sp_dist} checks Euclidean properties
#'
#'   \item \code{$details_fspaces} a list with 4 elements: \code{$sp_pc_coord}
#'   a matrix with coordinates of species (rows) along Principal Components
#'   (columns) with positive eigenvalues ; \code{$pc_eigenvalues} a matrix
#'   with eigenvalues of axes from PCoA ; \code{$dendro} a hclust
#'   object with the dendrogram details (null if no dendrogram computed) ;
#'   \code{$pairsp_fspaces_dist} a dataframe containing for each pair of
#'   species (rows), their names in the 2 first columns ('sp.x' and 'sp.y'),
#'   their distance based on trait-values ('tr'), and their Euclidean (for PCoA) or
#'   cophenetic (for dendrogram if computed) distance in each of the functional
#'   space computed ('pcoa_1d', 'pcoa_2d', ... , 'tree_clust');
#'   if fdist_scaling = TRUE, \code{$pairsp_fspaces_dist_scaled} a dataframe
#'   with scaled values of distances in functional spaces.
#'
#'   \item \code{$details_deviation} a list of dataframes:
#'   \code{$dev_distsp} a dataframe containing for each space (columns) the
#'   difference for all species pairs (rows) of the distance in the functional
#'   space and trait-based distance (i.e. positive deviation indicates
#'   overestimation of actual distance) ; \code{$abs_dev_distsp} and/or
#'   \code{$sqr_dev_distsp}, dataframes with for each space (columns) and all
#'   species pairs (rows) the absolute or squared deviation of distance ; if
#'   fdist_scaling = TRUE \code{$dev_distsp_scaled}, and
#'   \code{$abs_dev_distsp_scaled} and/or \code{$sqr_dev_distsp_scaled},
#'   dataframes with deviation computed on scaled distance in functional spaces.
#'
#'   }
#'
#' @note the maximum number of dimensions considered for assessing quality of
#'  functional spaces depends on number of PC axes with positive eigenvalues
#'  (i.e. axes with negative eigenvalues are not considered); so it could be
#'  lower than \code{$maxdim_pcoa}.
#'   The quality metric obtained with deviation_weighting = 'squared' and
#'   fdist_scaling = TRUE is equivalent to the square-root of the 'mSD'
#'   originally suggested in \emph{Maire et al. 2015}.
#'
#' @examples
#' # Load Species x Traits Data
#' data("fruits_traits", package = "mFD")
#'
#' # Load Traits x Categories Data
#' data("fruits_traits_cat", package = "mFD")
#'
#' # Compute Functional Distance
#' sp_dist_fruits <- mFD::funct.dist(
#'   sp_tr       = fruits_traits,
#'   tr_cat      = fruits_traits_cat,
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
#' @importFrom dendextend is.dist dist_long
#' @importFrom stats sd dist hclust cophenetic
#' @importFrom ade4 is.euclid
#' @importFrom ape pcoa
#'
#' @export


quality.fspaces <- function(sp_dist, fdendro = NULL, maxdim_pcoa = 10,
                            deviation_weighting = "absolute",
                            fdist_scaling = FALSE) {

  # check_inputs #

  # check distance object:
  if ( ! dendextend::is.dist(sp_dist)) {
    stop("Error: Trait-based distance between species should be
             provided as a dist object." )
  }
  # check species names:
  if (is.null(labels(sp_dist))) {
    stop("Error: No species names provided in distance matrix.")
  }
  # check the number of species:
  if (length(sp_dist) < 3) {
    stop("Error: there must be at least 3 species in 'sp_dist'.")
  }
  # check the absence of NA:
  if (any(is.na(sp_dist))) {
    stop("Error: NA are not allowed in 'sp_dist'.")
  }

  # check that the maximum number of axes to keep after PCoA:
  if (maxdim_pcoa < 1) {
    stop("Error: the number of pcoa axes must be higher than 0.")
  }

  # check that the name of quality metric is correct:
  if ( any( ! deviation_weighting %in% c("absolute", "squarred") ) ) {
    stop("Error: input 'deviation_weighting' should be 'absolute' and/or 'squarred'.")
  }

  # check that the scaled_metric input is logical:
  if ( any (! is.logical(fdist_scaling) ) ) {
    stop("Error: input 'fdist_scaling' should be TRUE and/or FALSE.")
  }


  # summary of input = trait-based distance between species ####

  # species names from input distance matrix
  sp_nm <- labels(sp_dist)

  # compute min, max, mean, standard deviation of distances
  trdist_summary <- c(min = min(sp_dist), max = max(sp_dist),
                      mean = mean(sp_dist), sd = stats::sd(sp_dist))

  # checking whether it is Euclidean
  trdist_euclidean <- ade4::is.euclid(sp_dist)

  # storing  summary of inputs
  details_trdist <- list(trdist_summary = trdist_summary,
                         trdist_euclidean = trdist_euclidean)

  # Computing functional spaces  #

  # computing PCoA ----
  pcoa_trdist <- ape::pcoa(sp_dist)
  # number of dimensions to keep given the input from user and number of PC...
  # ... with positive eigenvalues:
  nbdim <- min(maxdim_pcoa, ncol(pcoa_trdist$vectors))
  # keeping species coordinates on the 'nbdim' axes and renaming axes:
  sp_pc_coord <- pcoa_trdist$vectors[, 1:nbdim]
  colnames(sp_pc_coord) <- paste("PC", 1:nbdim, sep = "")
  ## check not run: any( rownames(sp_pc_coord)!= sp_nm)

  # storing details about PCoA:
  details_fspaces<-list (sp_pc_coord = sp_pc_coord,
                         pc_eigenvalues = pcoa_trdist$values[1:nbdim, ])

  # vector with names of all spaces from PCoA:
  fspaces_nm <- paste0("pcoa_", 1:nbdim, "d")

  # turning dist object with trait-based ditance in a 3-variables dataframe ...
  # ...with 1 row for each pair of species (names in the 2 first columns):
  df_distsp <- dendextend::dist_long(sp_dist)
  names(df_distsp) <- c("sp.x", "sp.y", "tr")

  # computing distance between species on increasing number of PCoA axes....
  #... and adding these pairwise distance ot the distance dataframe:

  # loop on number of axes kept
  for (k in 1:nbdim) {
    # computing Euclidean distance between species
    dist_k_dim<- stats::dist(sp_pc_coord[, 1:k])
    # storing these distances as additional column to distance dataframe:
    df_distsp[, paste0("pcoa_", k, "d")] <-
      dendextend::dist_long(dist_k_dim)$distance
  }

  # if needed, computing quality of dendrogram ----
  # if no dendrogram, null output:
  fdendro_hclust <- NULL
  # if dendrogram computed, compute the quality:
  if ( !is.null(fdendro)) {
    # checking if the name of clustering method is correct:
    if ( ! fdendro %in% c("ward.D", "ward.D2", "single", "complete", "average",
                          "mcquitty", "median", "centroid")) {
      stop("Error: Name of the clustering algorithm does not match with those
           allowed by stats::hclust. Please check.")
    }
    # adding dendro name to list of spaces
    fspaces_nm <- c(fspaces_nm, paste0("tree_", fdendro))
    # compute dendrogram using defined clustering algorithm and storing:
    fdendro_hclust <- stats::hclust(sp_dist, method = fdendro)
    details_fspaces$dendro <- fdendro_hclust
    # compute cophenetic distance between all species along the dendrogram:
    dist_fdendro <- stats::cophenetic(fdendro_hclust)
    ## check not run: any( rownames(as.matrix(dist_fdendro))!= sp_nm)
    # add pairwise cophenetic distance to the dataframe:
    df_distsp[, paste0("tree_", fdendro)] <-
      dendextend::dist_long(dist_fdendro)$distance
  }

  # storing distances in functional spaces
  details_fspaces$pairsp_fspaces_dist <- df_distsp


  # computing quality of all functional spaces ####

  # dataframe to store quality metrics of all spaces
  q_fspaces <- data.frame(fspaces_nm, row.names = fspaces_nm)


  # if required on raw distances in the functional spaces ----
  if (FALSE %in% fdist_scaling) {

    # compute deviation between distance in each functional space and ...
    # ... trait-based distance, and storing in a list:
    dev_distsp <- data.frame (df_distsp[, c("sp.x", "sp.y")],
                              df_distsp[,fspaces_nm ] - df_distsp[, "tr"])
    details_deviation <- list(dev_distsp = dev_distsp)

    # if required based on absolute deviation
    if ("absolute" %in% deviation_weighting) {

      # compute absolute deviation and storing:
      abs_dev_distsp <- data.frame (dev_distsp[, c("sp.x", "sp.y")],
                                    abs(dev_distsp[, fspaces_nm]))
      details_deviation$abs_dev_distsp <- abs_dev_distsp

      # mean absolute deviation:
      q_fspaces[fspaces_nm, "mad"] <- apply(abs_dev_distsp[, fspaces_nm],
                                            2, mean)
    }

    # if required based on squared deviation:
    if ("squarred" %in% deviation_weighting) {

      # compute squared deviation and storing:
      sqr_dev_distsp <- data.frame(dev_distsp[, c("sp.x", "sp.y")],
                                    (dev_distsp[, fspaces_nm])^2)
      details_deviation$sqr_dev_distsp <- sqr_dev_distsp

      # root of mean squared deviation:
      q_fspaces[fspaces_nm, "rmsd"] <- sqrt(apply(sqr_dev_distsp[, fspaces_nm],
                                                 2, mean))
    }

  }

  # if required, computing metrics on scaled distances in the functional spaces----
  if (TRUE %in% fdist_scaling) {

    # scaling distance in each functional space according to maximum distance...
    # ... in the functional space and maxium trait-based distance :
    df_distsp_scaled <- data.frame (
      df_distsp[, c("sp.x", "sp.y", "tr")],
      apply(df_distsp[, fspaces_nm], 2,
            function(x) {x / max(x) * max(df_distsp[, "tr"])})
    )
    details_fspaces$pairsp_fspaces_dist_scaled <- df_distsp_scaled

    # compute deviation between scaled distance and trait-based distance:
    dev_distsp_scaled <- data.frame(df_distsp_scaled[, c("sp.x", "sp.y")],
                                     df_distsp_scaled[, fspaces_nm] - df_distsp_scaled[, "tr"])
    details_deviation$dev_distsp_scaled <-  dev_distsp_scaled

    # if required based on absolute deviation
    if ("absolute" %in% deviation_weighting ) {

      # compute absolute deviation and storing:
      abs_dev_distsp_scaled <- data.frame (dev_distsp_scaled[, c("sp.x", "sp.y")],
                                           abs(dev_distsp_scaled[, fspaces_nm]))
      details_deviation$abs_dev_distsp_scaled <- abs_dev_distsp_scaled

      # mean absolute deviation:
      q_fspaces[fspaces_nm, "mad_scaled"] <- apply(abs_dev_distsp_scaled[, fspaces_nm],
                                                  2, mean)
    }

    # if required based on squared deviation:
    if ("squarred" %in% deviation_weighting) {

      # compute squared deviation and storing:
      sqr_dev_distsp_scaled <- data.frame (dev_distsp_scaled[, c("sp.x", "sp.y")],
                                           (dev_distsp_scaled[, fspaces_nm])^2)
      details_deviation$sqr_dev_distsp_scaled <- sqr_dev_distsp_scaled

      # root of mean squared deviation:
      q_fspaces[fspaces_nm, "rmsd_scaled"] <- sqrt(apply(sqr_dev_distsp_scaled[, fspaces_nm],
                                                        2, mean))
    }

  }

  quality_fspaces <- as.matrix(q_fspaces)
  quality_fspaces <-  quality_fspaces[, -1, drop = FALSE]
  quality_fspaces <-  as.data.frame(quality_fspaces)
  for (i in (1:ncol(quality_fspaces))) {
    quality_fspaces[, i] <- as.numeric(quality_fspaces[, i])
  }
  colnames(quality_fspaces) <- colnames(q_fspaces)[-1]


  # grouping and returning results:
  return_list <- list(quality_fspaces =  quality_fspaces,
                      details_trdist = details_trdist,
                      details_fspaces = details_fspaces,
                      details_deviation = details_deviation)

  if (min(sp_dist) == 0) {
    warning("Functional distance between some species is equal to 0 (explains the Warning message 1). You can choose to gather species into Functional Entities gathering species with similar traits values")
  }

  return(return_list)
} # end of function
