# Function to compute a set of functional diversity indices based on...
# ... functional entities
#
# Authors: Camille Magneville, Sébastien Villéger
#
#

# ------------------------------------------------------------------------------

#' Compute the set of indices based on number of species in functional entities (FEs)
#' following Mouillot et al (2014)
#'
#'@param asb_sp_occ a \strong{matrix} linking occurrences (coded as 0/1) of species
#'  (columns) in a set of assemblages (rows). Warning: \strong{An assemblage
#'  must contain at least one species}.
#'
#'@param sp_to_fe a \strong{list} with details of species clustering into FE from
#'  \code{sp.to.fe}.
#'
#'@param ind_nm a \strong{vector} of character strings with the names of functional
#'  diversity indices to compute among "fred", "fored" and "fvuln".
#'  \bold{Indices names must be written in lower case letters}. Default: all
#'  the indices are computed.
#'
#'@param check_input a \strong{logical value} allowing to test or not the inputs.
#'  Possible error messages will thus may be more understandable for the user
#'  than R error messages. Default: check_input = TRUE.
#'
#'@param details_returned a \strong{logical value} indicating whether the user wants to
#'  details_returned. Details are used in graphical functions and thus must be kept
#'  if the user want to have graphical outputs for the computed indices.
#'
#'@return  \itemize{ \item \emph{asb_fdfe} a matrix containing for each
#'  assemblage (rows), values of functional diversity indices (same names than
#'  in 'ind_nm') as well as the number of species ('nb_sp') and the number of FE
#'  (nb_fe). \item if \emph{details_returned} is TRUE, \item \emph{details_fdfe} a
#'  list with \emph{asb_fe_nbsp} a matrix with number of species per FE in each
#'  assemblage.}
#'
#'@examples
#' # Load Species*Traits dataframe:
#' data("fruits_traits", package = "mFD")
#' # Load Traits categories dataframe:
#' data("fruits_traits_cat", package = "mFD")
#' # Load Assemblages*Species matrix:
#' data("baskets_fruits_weights", package = "mFD")
#' # Remove continuous trait:
#' fruits_traits <- fruits_traits[, -5]
#' fruits_traits_cat <- fruits_traits_cat[-5, ]
#' # Compute gathering species into FEs:
#' sp_to_fe_fruits <- mFD::sp.to.fe(sp_tr = fruits_traits, tr_cat = fruits_traits_cat, 
#'  fe_nm_type = "fe_rank", check_input = TRUE)
#' # Get the occurrence dataframe:
#' asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights) 
#' asb_sp_fruits_occ <- asb_sp_fruits_summ$"asb_sp_occ"
#' # Compute alpha fd indices:
#' alpha.fd.fe(asb_sp_occ = asb_sp_fruits_occ, sp_to_fe = sp_to_fe_fruits,
#'  ind_nm = c("fred", "fored", "fvuln"),
#'  check_input = TRUE, details_returned = TRUE)
#'
#'
#'@export


alpha.fd.fe <- function(asb_sp_occ, sp_to_fe,
                        ind_nm = c("fred", "fored", "fvuln"),
                        check_input = TRUE, details_returned = TRUE) {
  
  
  # check_inputs if asked:
  if(check_input == TRUE) {
    if (is.null(rownames(asb_sp_occ))) {
      stop("Error: No row names provided in species*assemblage matrix.
             Please add assemblage names as column names.")
    }
    if (is.null(colnames(asb_sp_occ))) {
      stop("Error: No column names provided in species*assemblage matrix.
             Please add species names as column names.")
    }
    if (any(is.na(asb_sp_occ))) {
      stop("Error: The species*assemblage matrix contains NA. Please check.")
    }
    
    if (any((asb_sp_occ != 0 & asb_sp_occ != 1))) {
      stop("Error: The species*assemblage occurrence matrix should only contain 0 or 1. Please check.")
    }
    
    if (any(! ind_nm %in% c("fred", "fored", "fvuln"))) {
      stop("Error: Name of the functional indice to compute does not match with
         those allowed. Please check.")
    }
    
    if (!("sp_fe" %in% names(sp_to_fe))) {
      stop("Error: Input 'sp_to_fe' should contain a vector named 'sp_fe' (as returned by function 'sp.to.fe').
           Please check.")
    }
    
    if (any(!(colnames(asb_sp_occ) %in% names(sp_to_fe$sp_fe)))) {
      stop("Error: All species in species*assemblage occurence matrix  should be in the input 'sp_to_fe'.
          Please check that the function 'sp.to.fe' was applied on trait values of the same set of species.")
    }
    
    # Add a stop if some species do not belong to any assemblage:
    if (min(apply(asb_sp_occ, 2, sum)) == 0){
      stop("Error: Some species are absent from all assemblages.")
    }
    # Add a stop if some asb do not contain species:
    if (min(apply(asb_sp_occ, 1, sum)) == 0){
      stop("Error: Some assemblages do not contain species.")
    }
    
  }
  
  # change occ matrix as a df:
  asb_sp_occ <- as.data.frame(asb_sp_occ)
  
  # names of assemblages and of FE in the pool
  nmasb <- rownames(asb_sp_occ)
  nmfe_pool <- unique(sp_to_fe$sp_fe)
  
  # matrix to store indices values for each assemblage:
  asb_fdfe <- matrix(NA, length(nmasb), length(ind_nm) + 2,
                     dimnames = list(nmasb, c("nb_sp", "nb_fe", ind_nm)))
  
  # matrix to store number of species per FE for each assemblage:
  asb_fe_nbsp <- matrix(0, length(nmasb), length(nmfe_pool),
                        dimnames = list(nmasb, nmfe_pool))
  
  
  # for each assemblage, compute indices:
  for (k in nmasb) {
    
    # retrieve names of species present in assemblage k:
    asb_sp_occ_k <- asb_sp_occ[k, , drop = FALSE]
    to_remove <- colnames(asb_sp_occ_k[, asb_sp_occ_k == 0])
    nmsp_k <- colnames(asb_sp_occ[which(! (colnames(asb_sp_occ) %in% to_remove))])
    
    # number of species
    nbsp_k <- length(nmsp_k)
    asb_fdfe[k, "nb_sp"] <- nbsp_k
    
    # FEs to which these species belong
    fe_k <- sp_to_fe$sp_fe[nmsp_k]
    
    # names of FEs present
    nmfe_k <- unique(fe_k)
    
    # number of FEs
    nbfe_k <- length(nmfe_k)
    asb_fdfe[k, "nb_fe"] <- nbfe_k
    
    # ration between number of species and number of FE (i.e. Funct redundancy)
    fred_k <- nbsp_k / nbfe_k
    
    # number of species in each FE
    fe_nbsp_k <- unlist(lapply(nmfe_k, function(x) {length(which(fe_k == x))}))
    names(fe_nbsp_k) <- nmfe_k
    asb_fe_nbsp[k, nmfe_k] <- fe_nbsp_k
    
    
    # compute functional redundancy for assemblage k:
    if ("fred" %in% ind_nm) {
      asb_fdfe[k, "fred"] <- fred_k # or mean(fe_nbsp_k)
    }
    
    # compute functional over-redundancy for assemblage k:
    if ("fored" %in% ind_nm) {
      fored_k <- sum(sapply(fe_nbsp_k,
                            function(x) {max(c(x, fred_k) - fred_k) })) / nbsp_k
      asb_fdfe[k, "fored"] <- fored_k
    }
    
    # compute functional vulnerability for assemblage k:
    if ("fvuln" %in% ind_nm) {
      fvuln_k <- length(which(fe_nbsp_k == 1)) / nbfe_k
      asb_fdfe[k, "fvuln"] <- fvuln_k
    }
    
  }
  
  # outputs
  if (details_returned == TRUE) {
    res <- list(asb_fdfe = asb_fdfe,
                details_fdfe = list(asb_fe_nbsp = asb_fe_nbsp))
  } else {
    res <- asb_fdfe
  }
  
  return(res)
  
}
