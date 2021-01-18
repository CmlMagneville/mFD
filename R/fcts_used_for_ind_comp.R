# Four functions used to compute functional indices
#
# Authors: Camille Magneville & Sébastien Villéger
#
#

# ------------------------------------------------------------------------------


#' Compute distances of all points to a given point in the functional space
#'
#' This function computes the distances of all species to a reference species.
#' It is used in FSpe, FOri and FNND computation.
#'
#' @param sp_faxes_coord a \bold{matrix} of species coordinates in a chosen functional
#'   space. Species coordinates have been retrieved thanks to
#'   \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param ref_sp a \bold{character string} referring to the name of the reference
#'   species
#'
#' @return a vector of species distances to the reference species \code{dist_refsp}
#'

dist.point <- function(sp_faxes_coord, ref_sp) {
  # build the matrix of distances between species to be used to compute ...
  # ... distances to a given species (a reference species):
  sp_faxes_coord <- as.data.frame(sp_faxes_coord)
  dist_sp <- as.matrix(stats::dist(sp_faxes_coord, method = "euclidean"))
  dist_refsp <- dist_sp[which(rownames(dist_sp) == ref_sp), ]
  return(dist_refsp)
}


# ------------------------------------------------------------------------------


#' Compute distance of a given point to its nearest neighbor in the functional
#' space and the identity of the nearest neighbor
#'
#' This function is used in functional indices computation.
#'
#' @param sp_faxes_coord a \bold{matrix} of species coordinates in a chosen
#'   functional space. Species coordinates have been retrieved thanks to
#'   \code{fspace.conttr} or \code{qual.funct.space}.
#'
#' @param ref_sp a \bold{character string} referring to the name of the reference
#'   species
#'
#' @return a list containing the nearest neighbor identity \code{nn_id} and a
#'   list of the distance of the reference point to its nearest neighbor
#'   \code{nn_ref_sp_dist}.


dist.nearneighb <- function(sp_faxes_coord, ref_sp) {
  dist_sp <- as.matrix(stats::dist(sp_faxes_coord, method = "euclidean"))
  dist_sp[which(dist_sp == 0)] <- NA
  # search the nearest neighbor identity of each point:
  nn_id <- dist_sp - apply(dist_sp, 1, min, na.rm = TRUE)
  nn_id[which(nn_id != 0)] <- NA
  nn_id[which(nn_id == 0)] <- 1
  # keep information only for the reference species:
  nn_refsp_id <- nn_id[which(rownames(nn_id) == ref_sp), ]
  nn_refsp_id <- nn_refsp_id[which(nn_refsp_id == 1)]
  nn_refsp_id <- as.data.frame(nn_refsp_id)
  nn_id <- rownames(nn_refsp_id)
  nn_ref_sp_dist <- dist_sp[
    which(rownames(dist_sp) == ref_sp),
    which(colnames(dist_sp) %in% as.vector(nn_id[1]))]
  return_list <- list(nn_id, nn_ref_sp_dist)
  names(return_list) <- c(
    "nearest neighbour identity",
    "distance of the reference species to its nearest neighbour")
  return(return_list)
}


# ------------------------------------------------------------------------------


#' Compute the minimum spanning tree (MST) linking species of a given assemblage
#'
#' This function computes the MST linking species of a given assemblage and is
#' used to compute FEve index.
#'
#' @param sp_faxes_coord_k a \bold{matrix} relating species coordinates for species
#'   present in a given assemblage
#'
#' @return a dist object summarizing the MST for all species of a given
#'   assemblage \code{mst_asb_k}
#'

mst.computation <- function(sp_faxes_coord_k) {
  sp_dist_asb_k <- stats::dist(sp_faxes_coord_k, method = "euclidian")
  mst_asb_k <- ape::mst(sp_dist_asb_k)
  mst_asb_k <- stats::as.dist(mst_asb_k)
  return(mst_asb_k)
}


# ------------------------------------------------------------------------------


#' Compute vertices of the minimal convex hull shaping species from a single
#' assemblage in a multidimensional functional space.
#'
#' This function identifies species that are vertices of the minimal convex hull
#' enclosing a community in a multidimensional functional space. This function is
#' using the Qhull (http://www.qhull.org/) library through the
#' \code{geometry::convhulln} function.
#'
#' @param sp_faxes_coord a \bold{matrix} of species coordinates in a chosen functional
#'   space. Species coordinates have been retrieved thanks to
#'   \code{fspace.conttr} or \code{qual.funct.space}.
#'
#' @param order_2D a \bold{logical value} defining whether vertices names are reordered
#' so that they define a convex polygon in 2D which is convenient for plotting.
#' Default is FALSE, vertices ordered as in row names of 'sp_faxes_coord'
#'
#' @param check.input a \bold{logical value} defining whether inputs are checked before
#'   computation: species names must be put as row.names, there must be no NA
#'   and species number must be superior to (axes number + 1). Default:
#'   check.input = TRUE.
#'
#' @return a vector containing names of species being vertices \code{vert_nm}
#'
#'

vertices <- function(sp_faxes_coord, order_2D = FALSE, check.input = FALSE) {
  
  # checking input if required:
  if (check.input == TRUE) {
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
           Please add species names as row names.")
    }
    if (any(is.na(sp_faxes_coord))) {
      stop("Error: the species*coordinates matrix must not contain NA.
           Please check")
    }
    if (nrow(sp_faxes_coord) <= ncol(sp_faxes_coord)) {
      stop("Error: Number of species should be strictly higher than number of
           axes for computing the convex hull.")
    }
    
    if (order_2D == TRUE & ncol(sp_faxes_coord)!=2) {
      stop("Error: Ordering vertices could be done only if there are 2 axes.")
    }
  }
  
  # applying convhulln function to compute convexhull...
  # ... with options = 'Fx', to only compute vertices:
  # if vertices can be computed (coplanearity), convFx == NA:
  conv_Fx <- tryCatch(geometry::convhulln(sp_faxes_coord, options = "Fx"),
                      error = function(err){"NA"})
  
  # extracting unique names of vertices from the matrix with identity ...
  # ...of species for each facet:
  
  # if vertices have been computed (no coplanearity)
  if (! is.character(unlist(conv_Fx))) {
    vert_nm <- row.names(sp_faxes_coord)[sort(unique(as.vector((conv_Fx))))]
    
    # if required and only 2 dimensions, reordering vertices names clockwise:
    if (order_2D == TRUE & ncol(sp_faxes_coord) == 2) {
      
      vert_nm <- vert_nm[ order( -1 * atan2(
        sp_faxes_coord[vert_nm, 1] - mean(range(sp_faxes_coord[vert_nm, 1])),
        sp_faxes_coord[vert_nm, 2] - mean(range(sp_faxes_coord[vert_nm, 2])))) ]
    }
  }
  
  # if vertices have not been computed (coplanearity problem):
  if (is.character(unlist(conv_Fx))) {
    vert_nm <- NA
  }
  return(vert_nm)
}
