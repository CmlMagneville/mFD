#' Compute Functional Identity
#'
#' This function computes the weighted average position along each axis. FIde is
#' computed using relative weight so that it is not affected by unit (e.g. g or
#' kg for biomass). In the special case where 'weight' is filled with only 0/1
#' (absence/presence), then FIde will be computed assuming that all species have
#' the same weight. The results of this function are used in FSpe, FOri and FNND
#' computation.
#'
#' @param sp_faxes_coord_k a matrix of species coordinates present in a
#'   given assemblage in a chosen functional space with only needed axes.
#'   Species coordinates have been retrieved thanks to \code{tr.cont.fspace} or
#'   \code{\link{quality.fspaces}} and filtered thanks to
#'   \code{\link{sp.filter}}.
#'
#' @param asb_sp_relatw_k a matrix containing species relative weight
#'   (columns) for a given assemblage.
#'   
#' @param k a character string referring to the assemblage studied.
#'
#' @param check_input a logical value allowing to test or not the
#'   inputs. Possible error messages will thus may be more understandable for
#'   the user than R error messages. Species coordinates matrix and
#'   species*weight data frame must not contain NA, their rownames must be filled
#'   and they must have similar names values. Default: `check_input = FALSE`.
#'
#' @return A matrix containing functional identity values for a given assemblage
#'   along the dimensions (columns). Number of dimensions is fixed to the number
#'   of dimensions in \code{sp_faxes_coord} data frame.
#'   
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fide.computation <- function(asb_sp_relatw_k, sp_faxes_coord_k, 
                             k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'relat_sp_w_asb_k' and ", 
           "'sp_faxes_coord_k'. Please check.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }

  return(asb_sp_relatw_k %*% sp_faxes_coord_k)
}



#' Compute Functional Dispersion
#'
#' This function computes the weighted deviation to center of gravity of species
#' present in a given assemblage. It is the weighted mean distance to the
#' weighted centroid. Its calculation requires FIde computation that can be
#' achieved thanks to \code{fide.computation} function. FDis value can be scaled
#' by the maximum value possible given species pool (i.e. the most distant
#' species pair have half of total weight) to standardize values.
#'
#' @param sp_faxes_coord_k a matrix of species coordinates present in a
#'   given assemblage in a chosen functional space with only needed axes.
#'   Species coordinates have been retrieved thanks to
#'   \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}} and filtered
#'   thanks to \code{\link{sp.filter}}.
#'   
#' @param asb_sp_relatw_k a matrix containing species relative weight
#'   (columns) for a given assemblage.
#'
#' @param fide_asb a matrix containing functional identity values a
#'   given assemblage along the dimensions (columns). Can be retrieved after
#'   \code{fide.computation} function and is compute if NULL. Default: fide_asb
#'   = NULL.
#'
#' @param k a character string referring to the assemblage studied.
#'
#' @param check_input a logical value allowing to test or not the inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix and species*weight
#'   dataframe must not contain NA, their rownames must be filled and they must
#'   have similar names values. Default: check_input = FALSE.
#'
#' @return A matrix containing functional dispersion for a given assemblage.
#'   
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fdis.computation <- function(asb_sp_relatw_k, sp_faxes_coord_k, 
                             fide_asb = NULL, k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'relat_sp_w_asb_k' and ", 
           "'sp_faxes_coord_k'. Please check.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  
  # compute fide if NULL:
  if (is.null(fide_asb)) {
    
    fide_asb_k <- asb_sp_relatw_k %*% sp_faxes_coord_k
    # compute distance to the centroid to compute fdis:
    dist_centr_k <- apply(sp_faxes_coord_k, 1, 
                          function(x) (sum((x - fide_asb_k) ^ 2)) ^ 0.5)
    # compute fdis value:
    fdis_asb_k <- (asb_sp_relatw_k %*% dist_centr_k)
    
    return(fdis_asb_k)
  }
  
  if (!is.null(fide_asb)) {
    # compute distance to the centroid to compute fdis:
    dist_centr_k <- apply(sp_faxes_coord_k, 1, 
                          function(x) {
                            (sum((x - fide_asb[k, 
                                       colnames(sp_faxes_coord_k)]) ^ 2)) ^ 0.5
                          })
    # compute fdis value:
    fdis_asb_k <- (asb_sp_relatw_k %*% dist_centr_k)
    
    return(fdis_asb_k)
  }
}



#' Compute Functional Richness
#'
#' This function computes the volume of functional space filled by species
#' present in a given assemblage.
#'
#' @param sp_faxes_coord_k a matrix of species coordinates present in a given
#'   assemblage in a chosen functional space with only needed axes. Species
#'   coordinates have been retrieved thanks to \code{\link{tr.cont.fspace}} or
#'   \code{qual.funct.space} and filtered thanks to \code{\link{sp.filter}}.
#'
#' @param k a character string referring to the assemblage studied.
#'
#' @param check_input a logical value allowing to test or not the inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix must not contain NA, its
#'   rownames must be filled and the number of species should strictly be higher
#'   than the number of axes to compute the convex hull. Default: check_input =
#'   FALSE.
#'
#' @return a list containing: \strong{$fric} a vector with fric value for a
#'   given assemblage and \strong{vertices_nm} a vector containing names of the
#'   species being vertices (species are ordered as in row names of input)
#'
#' @note Computation with qconvex algorithm is led using option 'Tv' so result
#'   are verified for structure, convexity, and point inclusion. FRic value is
#'   based on axes units. See \code{\link{alpha.fd.multidim}} for option to 
#'   scale values using volume by species pool.
#'
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fric.computation <- function(sp_faxes_coord_k, k, check_input = check_input) {
  
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (nrow(sp_faxes_coord_k) < ncol(sp_faxes_coord_k)) {
      stop("Number of species should strictly be higher than the number of ", 
           "axes to compute the convex hull. Problem for assemblage ", k)
    }
  }
  
  # applying convhulln function to compute
  # convexhull...  ... with options = 'FA', to
  # compute the general area of the functional hull,
  # if convex hulln can not be computed
  # (coplanearity), takes the value: NA
  conv_fa_k <- tryCatch(geometry::convhulln(sp_faxes_coord_k, 
                                            option = "FA"), 
                        error = function(err) "NA")
  
  # extracting unique names of vertices from the
  # matrix with identity of ...  ... species for each
  # facet if vertices have been computed (no
  # coplanearity pbs): ... and get the raw fric
  # value:
  if (!is.character(conv_fa_k)) {
    vert_nm_k <- row.names(sp_faxes_coord_k)[
      sort(unique(as.vector((conv_fa_k$hull))))]
    fric <- conv_fa_k$vol
  }
  
  # if vertices have not been computed (coplanearity
  # pbs):
  if (is.character(conv_fa_k)) {
    vert_nm_k <- NA
    fric <- NA
  }
  
  return_list <- list(fric = fric, vertices_nm = vert_nm_k)
  
  return(return_list)
}



#' Compute Functional Divergence (FDiv) Index for One Assemblage
#'
#' This function to compute Functional Divergence (FDiv) index for one
#' assemblage. FDiv indice accounts for deviation of biomass to the center of
#' gravity of the vertices shaping the convex hull.
#' FDiv is scaled between 0 and 1. For details about FDiv index see Villeger _et
#' al._ 2008. Use \code{\link{alpha.fd.multidim}} to compute FDiv over multiple
#' assemblages (and together with other FD indices) .
#'
#' @param sp_faxes_coord_k a matrix with species coordinates for species present
#'   in a given assemblage along functional axes.
#'
#' @param asb_sp_relatw_k a matrix containing species relative weight
#'   (columns) for a given assemblage.
#'
#' @param vert_nm a vector with names of the species being vertices of the
#'   convex hull (so should be a subset of colnames of \code{sp_faxes_coord_k}.
#'   This vector can be provided through \code{vertices}) function or
#'   retrieved after the \code{fric.computation} function. Default: vert_nm =
#'   NULL.
#'
#' @param k a character string referring to the assemblage studied.
#'
#' @param check_input a logical value allowing to test or not the inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix and species*weight
#'   dataframe must not contain NA, their rownames must be filled, they must
#'   have similar names values, species acting as vertices must be contained in
#'   the species coordinates matrix and  the number of species should strictly
#'   be higher than number of axes for computing the convex hull. Default:
#'   check_input = FALSE.
#'
#' @return A list with \strong{$fdiv} a single value vector ; \strong{$details}
#'   a list with $vertices_nm a vector containing names of the species being
#'   vertices ; \strong{$B_coord} a vector with coordinates of the center of
#'   gravity ; \strong{$mean_dtoB} a single value with average distance of
#'   species to center of gravity of vertices.
#'
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd


fdiv.computation <- function(sp_faxes_coord_k, asb_sp_relatw_k, 
                             vert_nm = NULL, k, check_input = check_input) {
  
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'asb_sp_relatw_k' and ", 
           "'sp_faxes_coord_k'.")
    }
    
    if (!is.null(vert_nm)) {
      if (any((vert_nm %in% row.names(sp_faxes_coord_k) == FALSE))) {
        stop("Names of the vertices are not all present in species ", 
             "coordinates matrix. Please check.")
      }
    }
    
    if (is.null(vert_nm)) {
      if (nrow(sp_faxes_coord_k) <= ncol(sp_faxes_coord_k)) {
        stop("Number of species should strictly be higher than number of axes ", 
             "for computing the convex hull. Please check.")
      }
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  
  # if vertices names are not provided, compute
  # vertices:
  if (is.null(vert_nm)) {
    # computes vertices names:
    vert_nm <- vertices(sp_faxes_coord_k, check_input = FALSE)
  }
  # compute fdiv values if vertices have been
  # computed (no coplanearity pbs):
  
  if (is.character(vert_nm)) {
    # get the coordinates of the vertices center of
    # gravity (named B):
    B_coord <- apply(sp_faxes_coord_k[vert_nm, ], 2, mean)
    # compute the euclidean distance of all species to
    # B:
    dtoB <- apply(sp_faxes_coord_k, 1, function(x) {
      (sum((x - B_coord)^2))^0.5
    })
    # compute mean of distances to B:
    mean_dtoB <- mean(dtoB)
    # compute deviation of species distances to B to
    # their mean:
    dev_dtoB <- dtoB - mean_dtoB
    # computes the weighted mean of raw deviations:
    ab_dev <- asb_sp_relatw_k * dev_dtoB
    # compute the weighted mean of absolute deviations:
    ab_absdev <- asb_sp_relatw_k * abs(dev_dtoB)
    # computing fdiv index:
    fdiv_asb_k <- (sum(ab_dev) + mean_dtoB) / (sum(ab_absdev) + mean_dtoB)
  }
  
  # if vertices have not been computed (coplanearity
  # pb):
  if (!is.character(vert_nm)) {
    fdiv_asb_k <- NA
    B_coord <- NA
    mean_dtoB <- NA
  }
  
  return_list <- list(fdiv = fdiv_asb_k, details = list(vertices_nm = vert_nm, 
                                                        B_coord = B_coord, 
                                                        mean_dtoB = mean_dtoB))
  
  return(return_list)
}



#' Compute Functional Evenness (FEve)
#'
#' This function computes the regularity of distribution of species weights in
#' the functional space.
#'
#' @param sp_faxes_coord_k a matrix of species coordinates present in a
#'   given assemblage in a chosen functional space with only needed axes.
#'   Species coordinates have been retrieved thanks to 
#'   \code{\link{tr.cont.fspace}} or
#'   \code{\link{quality.fspaces}} and filtered thanks to 
#'   \code{\link{sp.filter}}.
#'
#' @param asb_sp_relatw_k a \strong{matrix} containing species relative weight 
#' (columns) for a given assemblage.
#'
#' @param k a \strong{character string} referring to the assemblage studied.
#'
#' @param check_input a \strong{logical value} allowing to test or not the 
#' inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix and species*weight
#'   dataframe must not contain NA, their rownames must be filled, they must
#'   have similar names values, and the number of species in the assemblage must
#'   be higher than three to compute feve. Default: check_input = FALSE.
#'
#' @return a matrix containing functional evenness for a given assemblage.
#' 
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

feve.computation <- function(asb_sp_relatw_k, sp_faxes_coord_k, 
                             k, check_input = check_input) {
  
  # get the number of species present in the
  # assemblage:
  sp_nb_asb_k <- nrow(sp_faxes_coord_k)
  
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'relat_sp_w_asb_k' and ", 
           "'sp_faxes_coord_k'.")
    }
    
    if (sp_nb_asb_k < 3) {
      stop("There must be at least 3 species in the assemblage to compute ",
           "feve. The assemblage ", k, " contains less than 3 species.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  
  # compute distances between species to calculate
  # weighted evenness indice:
  sp_dist_asb_k <- stats::dist(sp_faxes_coord_k, method = "euclidian")
  
  # compute a dist object summarizing the mst for
  # species:
  mst_asb_k <- mst.computation(sp_faxes_coord_k)
  
  # compute a dist_long object that can be useful to
  # return:
  mst_asb_k_2 <- dendextend::dist_long(sp_dist_asb_k)
  names(mst_asb_k_2) <- c("sp.x", "sp.y", "feve_mst")
  mst_asb_k_2 <- mst_asb_k_2[which(mst_asb_k_2$feve_mst != 0), ]
  
  # compute the pairwise cumulative abundances:
  cum_ab_asb_k <- matrix(0, nrow = sp_nb_asb_k, ncol = sp_nb_asb_k)
  for (i in (1:sp_nb_asb_k)) {
    for (j in (1:sp_nb_asb_k)) {
      cum_ab_asb_k[i, j] <- asb_sp_relatw_k[i] + asb_sp_relatw_k[j]
    }
  }
  cum_ab_asb_k <- stats::as.dist(cum_ab_asb_k)
  
  # compute the weighted evenness index for the
  # (number of species - 1)...  ... segments linking
  # species:
  ew_asb_k <- rep(0, sp_nb_asb_k - 1)
  ind <- 1
  for (m in (1:((sp_nb_asb_k - 1) * sp_nb_asb_k / 2))) {
    if (mst_asb_k[m] != 0) {
      ew_asb_k[ind] <- sp_dist_asb_k[m] / (cum_ab_asb_k[m])
      ind <- ind + 1
    }
  }
  
  # compute the minimum between partial weighted
  # evenness (pew) index and ...  ... 1/(number of
  # species - 1):
  min_pew_asb_k <- rep(0, sp_nb_asb_k - 1)
  comp_value <- 1 / (sp_nb_asb_k - 1)
  for (l in (1:(sp_nb_asb_k - 1))) {
    min_pew_asb_k[l] <- min((ew_asb_k[l] / sum(ew_asb_k)), comp_value)
  }
  
  # compute feve value:
  feve_asb_k <- round(((sum(min_pew_asb_k)) - comp_value)/(1 - comp_value), 6)
  
  return_list <- list(feve = feve_asb_k, mst = mst_asb_k, 
                      mst_2 = mst_asb_k_2)
  
  return(return_list)
}



#' Compute Functional mMan Pairwise Distance (FMPD)
#'
#' This function computes the mean weighted distance between all pairs of
#' species. FMPD value can be scaled by the maximum value possible given species
#' pool (i.e. the most distant species pair have total weight) to standardize
#' values.
#'
#' @param sp_faxes_coord_k a \strong{matrix} of species coordinates present in 
#'   a given
#'   assemblage in a chosen functional space with only needed axes. Species
#'   coordinates have been retrieved thanks to \code{\link{tr.cont.fspace}} or
#'   \code{\link{quality.fspaces}} and filtered thanks to \code{sp.filter}.
#'
#' @param asb_sp_relatw_k a \strong{matrix} containing species relative weight 
#'   (columns) for a given assemblage.
#'
#' @param k a \strong{character string} referring to the assemblage studied.
#'
#' @param check_input a \strong{logical value} allowing to test or not the 
#'   inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix and species*weight
#'   dataframe must not contain NA, their rownames must be filled and they must
#'   have similar names values. Default: check_input = FALSE.
#'
#' @return A matrix containing functional mean pairwise distance for a given
#'   assemblage.
#'   
#' @author Camille Magneville and Sebastien Villeger
#'
#' @noRd

fmpd.computation <- function(asb_sp_relatw_k, sp_faxes_coord_k, 
                             k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'asb_sp_relatw_k' and ", 
           "'sp_faxes_coord_k'.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  
  # compute mean distance between species of a given
  # community to...  ... compute fmpd value:
  dist_sp_asb_k <- as.matrix(stats::dist(sp_faxes_coord_k, 
                                         method = "euclidean"))
  dist_sp_asb_k[which(dist_sp_asb_k == 0)] <- NA
  mean_dist_sp_asb_k <- apply(dist_sp_asb_k, 1, mean, na.rm = TRUE)
  # compute fmpd value for a given assemblage:
  fmpd_asb_k <- (asb_sp_relatw_k %*% mean_dist_sp_asb_k)
  
  return(fmpd_asb_k)
}



#' Compute Functional Mean Nearest Neighbor Distance (FNND)
#'
#' This function computes the weighted mean distance to nearest neighbor. It
#' uses /code{dist.nearneighb} function that computes distance to the nearest
#' neighbor for each species.FNND value can be scaled by the maximum value
#' possible given species pool (i.e. the most distant species pair have total
#' weight) to standardize values.
#'
#' @param sp_faxes_coord_k a \strong{matrix} of species coordinates present in 
#'   a given
#'   assemblage in a chosen functional space with only needed axes. Species
#'   coordinates have been retrieved thanks to \code{\link{tr.cont.fspace}} or
#'   \code{\link{quality.fspaces}} and filtered thanks to
#'   \code{\link{sp.filter}}.
#'
#' @param asb_sp_relatw_k a \strong{matrix} containing species relative weight 
#'   (columns) for a given assemblage.
#'
#' @param k a \strong{character string} referring to the assemblage studied.
#'
#' @param check_input a \strong{logical value} allowing to test or not the 
#'   inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species coordinates matrix and species*weight
#'   dataframe must not contain NA, their rownames must be filled and they must
#'   have similar names values. Default: check_input = FALSE.
#'
#' @return A matrix containing functional mean nearest neighbor distance for a
#'   given assemblage.
#'
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fnnd.computation <- function(asb_sp_relatw_k, sp_faxes_coord_k, 
                             k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(sp_faxes_coord_k))) {
      stop("The species*coordinates matrix contains NA. Please check.")
    }
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord_k))) {
      stop("No row names provided in species*coordinates matrix. Please add ", 
           "species names as row names.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (sum(colnames(asb_sp_relatw_k) %in% rownames(sp_faxes_coord_k)) != 
        nrow(sp_faxes_coord_k)) {
      stop("Mismatch between names in 'asb_sp_relatw_k' and ", 
           "'sp_faxes_coord_k'.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  
  # create a list to store distance to the nn for
  # each species of...  ...the assemblage:
  dist_nn_k <- list()
  
  # create a list to store the identity of the nn for
  # each species of...  ...the assemblage:
  nm_nn_k <- list()
  
  # compute distance to the nearest neighbor in a
  # given assemblage to...  ... compute fnnd value
  # and compute name of the nearest neighbor to
  # return:
  for (i in (1:nrow(sp_faxes_coord_k))) {
    ref_sp <- rownames(sp_faxes_coord_k)[i]
    dist_nn_sp_asb_k <- dist.nearneighb(sp_faxes_coord_k, ref_sp)
    dist_nn_k[ref_sp] <- 
  dist_nn_sp_asb_k$`distance of the reference species to its nearest neighbour`
    nms <- list(dist_nn_sp_asb_k$`nearest neighbour identity`)
    names(nms) <- ref_sp
    nm_nn_k[ref_sp] <- nms
  }
  
  # compute fnnd:
  fnnd_asb_k <- (asb_sp_relatw_k %*% unlist(dist_nn_k))
  
  # get the return list of outputs:
  return_list <- list(fnnd = fnnd_asb_k, details = list(nm_nn_k = nm_nn_k, 
                                                        dist_nn_k = dist_nn_k))
  
  return(return_list)
}



#' Compute Functional Originality
#'
#' This function computes the weighted mean distance to nearest species from the
#' species pool. FOri value can be scaled by the maximum distance to the nearest
#' neighbour possible in the global species pool (i.e. an assemblage hosting
#' only the most original species).
#'
#' @param dist_nn_global_pool a \strong{vector} containing the minimal distance 
#'   to the nearest neighbor for each species of the global pool of species.
#'
#' @param asb_sp_relatw_k a \strong{matrix} containing species relative weight
#'   (columns) for a given assemblage.
#'
#' @param k a \strong{character string} referring to the assemblage studied.
#'
#' @param check_input a \strong{logical value} allowing to test or not the 
#'   inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species*weight dataframe must not contain NA and its
#'   rownames must be filled. Default: check_input = FALSE.
#'
#' @return a matrix containing functional originality for a given assemblage.
#'
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fori.computation <- function(dist_nn_global_pool, asb_sp_relatw_k, 
                             k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  nm_sp_asb_k <- colnames(asb_sp_relatw_k)
  fori_asb_k <- (asb_sp_relatw_k %*% as.matrix(dist_nn_global_pool[
    nm_sp_asb_k]))
  
  return(fori_asb_k)
}



#' Compute Functional Specialization
#'
#' This function computes the weighted mean distance to the centroid of the
#' global species pool. It computes the average position of all the species
#' present given assemblage. FSpe value can be scaled by the maximum distance to
#' the global pool centroid (i.e. an assemblage hosting only the most
#' specialized species).
#'
#' @param special_sp_global_pool a \strong{vector} containing the distance to 
#'   the centroid of the global pool of species for each species of the global 
#'   pool (it is called specialization).
#'
#' @param asb_sp_relatw_k a \strong{matrix} containing species relative weight
#'   (columns) for a given assemblage.
#'
#' @param k a \strong{character string} referring to the assemblage studied.
#'
#' @param check_input a \strong{logical value} allowing to test or not the 
#'   inputs.
#'   Possible error messages will thus may be more understandable for the user
#'   than R error messages. Species*weight dataframe must not contain NA and its
#'   rownames must be filled. Default: check_input = FALSE.
#'
#' @return A matrix containing functional specialization for a given assemblage.
#'
#' @author Camille Magneville and Sebastien Villeger
#' 
#' @noRd

fspe.computation <- function(asb_sp_relatw_k, special_sp_global_pool, 
                             k, check_input = check_input) {
  # check_inputs if required:
  if (check_input) {
    
    if (any(is.na(asb_sp_relatw_k))) {
      stop("The species*weights dataframe contains NA. Please check.")
    }
    
    if (is.null(rownames(asb_sp_relatw_k))) {
      stop("No row names provided in species*weights dataframe. Please add ", 
           "assemblages names as row names.")
    }
    
    if (round(sum(asb_sp_relatw_k), 10) != 1) {
      stop("The sum of relative weights is not equal to one for ", k)
    }
  }
  nm_sp_asb_k <- colnames(asb_sp_relatw_k)
  special_sp_asb_k <- special_sp_global_pool[nm_sp_asb_k]
  fspe_asb_k <- asb_sp_relatw_k %*% special_sp_asb_k
  
  return(fspe_asb_k)
}
