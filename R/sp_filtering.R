# Function to filter species by community
#
# Authors: Camille Magneville & Sébastien Villéger
#
#

# ------------------------------------------------------------------------------


#' Retrieve information about species in a given assemblage
#'
#' This function computes names of species present in an given assemblage, their
#' coordinates in the functional space and their weights. It is used in the
#' \code{alpha_FD_multidim} function to filter species and compute each
#' functional indices for each community.
#'
#' @param asb_nm a \strong{string} object referring to the name of a given
#' community.
#'
#' @param sp_faxes_coord a \strong{matrix} of species coordinates in a chosen functional
#'   space. Species coordinates have been retrieved thanks to
#'   \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param asb_sp_w a \strong{matrix} linking weight of species (columns) and a set of
#'   assemblages (rows).
#'
#' @return a vector containing names of species present in a given assemblage \code{sp_name_asb_k}, a
#'   matrix containing coordinates of species present in a given assemblage \code{sp_faxes_coord_k}, a
#'   matrix containing weight of species present in a given assemblage \code{asb_sp_w_k}, a matrix
#'   containing relative weight of species present in a given assemblage \code{asb_sp_relatw_k}.
#'
#' @examples
#' # Load Species*Traits dataframe:
#' data("sp_tr_fruits", package = "mFD")
#' # Load Assemblages*Species dataframe:      
#' data("asb_sp_w_fruits", package = "mFD") 
#' # Load Traits categories dataframe:
#' data("tr_cat_fruits", package = "mFD")  
#' # Compute functional distance 
#' sp_dist_fruits <- mFD::funct.dist(sp_tr = sp_tr_fruits,         
#'  tr_cat       = tr_cat_fruits,   
#'  dist_metric  = "kgower",         
#'  scaling      = "scaledBYrange",  
#'  stop_if_NA   = TRUE)
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#' fspaces_quality_fruits <- mFD::quality.fspaces(sp_dist = sp_dist_fruits, 
#'  maxdim_pcoa         = 10,
#'  deviation_weighting = "absolute",
#'  fdist_scaling       = FALSE,
#'  fdendro             = "average")
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$"details_fspaces"$"sp_pc_coord"
#' # Filter species of basket_1 assemblage:
#' sp.filter(asb_nm = "basket_1", sp_faxes_coord_fruits, asb_sp_w_fruits)
#'
#'@export

sp.filter <- function(asb_nm, sp_faxes_coord, asb_sp_w) {
  asb_sp_w <- as.data.frame(asb_sp_w)
  asb_sp_w2 <- asb_sp_w[asb_nm, ]
  sp_name_asb_k <- colnames(asb_sp_w2)[which(asb_sp_w2 != 0)]
  sp_faxes_coord_k <- sp_faxes_coord[sp_name_asb_k, , drop = FALSE]
  asb_sp_w_k <- asb_sp_w2[asb_nm, sp_name_asb_k, drop = FALSE]
  asb_sp_relatw_k <- asb_sp_w_k / sum(asb_sp_w_k)
  
  return_list <- list(sp_name_asb_k, sp_faxes_coord_k, asb_sp_w_k,
                      asb_sp_relatw_k)
  names(return_list) <- c("species names", "species coordinates",
                          "species weight", "species relative weight")
  return(return_list)
}
