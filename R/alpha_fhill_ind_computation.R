# ------------------------------------------------------------------------------
# Function to compute  functional alpha diversity indices based on Hill numbers
# applied to distance between species
#
# Authors: Sébastien Villéger and Camille Magneville
#
# ------------------------------------------------------------------------------


#'Compute functional alpha diversity indices based on Hill numbers
#'
#' Compute functional alpha diversity applied to distance between species
#' following the framework from Chao et al.2019, Ecological Monographs
#' (89:e01343), DOI: 10.1002/ecm.1343). FD is computed applying the special case
#' where function 'f' in equation 3c is linear:f(dij(tau)) = dij(tau)/tau, hence
#' f(0) = 0 and f(tau) = 1.
#'
#'@param asb_sp_w a \strong{matrix} with weight of species (columns) in a set of
#'  assemblages (rows). Rows and columns should have names. NA are not allowed.
#'
#'@param sp_dist a \strong{matrix or dist object} with distance between species. Species
#' names should be provided and match those in 'asb_sp_w'. NA are not allowed.
#'
#'@param q a \strong{vector} containing values referring to the order of diversity to use
#'
#'@param tau a \strong{character string} with name of function to apply to distance
#' matrix (i.e. among all pairs of species) to get the threshold used to define
#' 'functionally indistinct set of species'. Could be qet to 'mean' (default),
#' 'min' or 'max'.
#'
#'@param check_input a \strong{logical value} defining whether inputs are checked before
#'  computation of indices. Possible error messages will thus may be more
#'  understandable for the user than R error messages. Default: check_input =
#'  TRUE.
#'
#'@param details_returned a \strong{logical value} indicating whether the user want to store
#' values used for computing indices (see list below)
#'
#'@return a list with: \itemize{
#'
#'  \item \emph{asb_FD_Hill} a matrix containing indices values for each level of
#'  q (columns, named as 'FD_qx') for each assemblage (rows, named as in
#'  \strong{asb_sp_w})
#'  \item \emph{tau_dist} the threshold value applied to distance between
#'  species to compute diversity according to function provided in \strong{tau}
#'
#'  \item if \strong{details_returned} turned to TRUE a list \emph{details} with
#'  \itemize{
#'  \item \emph{asb_totw} a vector with total weight of each assemblage
#'  \item \emph{asb_sp_relw} a matrix with relative weight of species in
#'  assemblages
#'  }
#'  }
#'
#'@note FD computed with q=2 and tau = 'max' is equivalent to the Rao's quadratic
#' entropy from Ricotta & Szeidl (2009, J Theor Biol).
#' FD computed with tau = 'min' is equivalent to Hill number taxonomic diversity,
#' thus with q=0 it is species richness (S), with q = 1 it is exponential of Shannon
#' entropy (H) and with q = 2 it is 1/(1-D) where D is Simpson diversity
#' FD is computed applying the special case where function 'f' in equation 3c
#' is linear:f(dij(tau)) = dij(tau)/tau, hence f(0)=0 and f(tau)=1.
#'
#'@examples
#' # Load Species*Traits dataframe:
#' data("fruits_traits", package = "mFD")
#' # Load Assemblages*Species dataframe:      
#' data("baskets_fruits_weights", package = "mFD")   
#' # Compute functional distance 
#' sp_dist_fruits <- mFD::funct.dist(sp_tr = fruits_traits,         
#'  tr_cat       = fruits_traits_cat,   
#'  dist_metric  = "kgower",         
#'  scaling      = "scaledBYrange",  
#'  stop_if_NA   = TRUE)
#' # Compute alpha fd hill indices:
#' alpha.fd.hill(asb_sp_w = baskets_fruits_weights, sp_dist = sp_dist_fruits, q = c(0, 1, 2),
#'  tau = "mean", check_input = TRUE, details_returned = TRUE)
#'
#'@export

alpha.fd.hill <- function(asb_sp_w,
                          sp_dist,
                          q = c(0, 1, 2),
                          tau = "mean",
                          check_input = TRUE,
                          details_returned = TRUE) {
  
  
  ####  distance between species stored in a matrix  ####
  
  sp_sp_dist <- sp_dist
  
  if (is.matrix(sp_sp_dist) == FALSE) {
    sp_sp_dist <- as.matrix(sp_sp_dist)
  }
  
  ## check_inputs if required #####
  if (check_input == TRUE) {
    
    if (is.matrix(asb_sp_w) == FALSE) {
      stop("Error: 'asb_sp_w' must be a matrix")
    }
    
    if (any(is.na(sp_dist))) {
      stop("Error: The species distances matrix contains NA. Please check.")
    }
    if (is.null(rownames(sp_sp_dist))) {
      stop("Error: No row names provided in species distances matrix.
             Please add species names as row names.")
    }
    if (any(is.na(asb_sp_w))) {
      stop("Error: The species*weights dataframe contains NA. Please check.")
    }
    if (is.null(rownames(asb_sp_w))) {
      stop("Error: No row names provided in species*weights dataframe.
             Please add assemblages names as row names.")
    }
    if (is.null(colnames(asb_sp_w))) {
      stop("Error: No column names provided in species*assemblage dataframe.
             Please add species names as column names.")
    }
    
    isnum_vect <- sapply(asb_sp_w, is.numeric)
    
    if (FALSE %in% isnum_vect) {
      stop("Error: The 'asp_sp_w' dataframe must only contain numeric values. Please convert values")
    }
    
    if (any(! (colnames(asb_sp_w) %in% rownames(sp_sp_dist) ) ) ) {
      stop(paste("Error: Mismatch between names in species*weight and
                   species distances matrix. Please check."))
    }
    
    if(any(! q %in% c(0,1,2) ) ) {
      stop(paste("Error: q should be 0, 1 and/or 2.
                  Please check."))
    }
    
    if(any(! tau %in% c("min", "mean", "max") ) ) {
      stop(paste("Error: tau should be 'min', 'mean' or 'max'. Please check."))
    }
    
    # Add a stop if some species do not belong to any assemblage:
    if (min(apply(asb_sp_w, 2, sum)) == 0){
      stop("Error: Some species are absent from all assemblages.")
    }
    # Add a stop if some asb do not contain species:
    if (min(apply(asb_sp_w, 1, sum)) == 0){
      stop("Error: Some assemblages do not contain species.")
    }
    
    # Add a stop if there is a negative value in the occurrence dataframe:
    if (any(asb_sp_w < 0)) {
      stop("Error: The species*weight dataframe should not contain negative values.
           Please check.")
    }
    
  }
  
  ####  preliminary operations ####
  
  # ensuring species are in the same order in both matrices:
  asb_sp_w <- as.matrix(asb_sp_w)
  sp_sp_dist <- sp_sp_dist[colnames(asb_sp_w), colnames(asb_sp_w)]
  
  
  # computing total weight per assemblage and relative weights of species ----
  asb_totw <- apply(asb_sp_w, 1, sum)
  asb_sp_relw <- asb_sp_w/asb_totw
  
  # computing tau as mean or max on distances ----
  tau_dist <- NULL
  
  if( tau == "min") {
    tau_dist <- min(sp_dist)
    
    # special case of null distance outside diagonal
    if (tau_dist == 0) {
      tau_dist<- min(sp_dist[sp_dist != 0])
      cat("Warning: some species has null functional distance,
          'tau' was set to the minimum non-null distance")
    }
  }
  
  if( tau == "mean") {
    tau_dist <- mean(sp_dist)
  }
  
  if( tau == "max") {
    tau_dist <- max(sp_dist)
  }
  
  # applying tau threshold to distance matrix
  dij_tau <- sp_sp_dist
  dij_tau[which(dij_tau > tau_dist, arr.ind = T)] <- tau_dist
  
  
  ####  computing diversity of assemblages ####
  
  # empty matrix to store outputs
  asb_FD_Hill <- matrix(NA, nrow(asb_sp_w), length(q),
                        dimnames = list(row.names(asb_sp_w), paste0("FD_q", q)))
  
  # loop on assemblages (equations id refers to those in Chao et al 2019
  for (k in row.names(asb_sp_w)) {
    
    # total weight (n+ in eq 4a)
    nplus_k <- asb_totw[k]
    
    # species names present in assemblage k
    sp_k <- colnames(asb_sp_w)[which(asb_sp_w[k, ] > 0)]
    
    # f(dij(tau)) with f being linear so dij(tau)/tau (see bottom right of p7)
    f_dij_tau <- dij_tau[sp_k, sp_k]/tau_dist
    
    # 'abundance of species given tau' (eq 3c)
    a_k <- (1 - f_dij_tau) %*% asb_sp_w[k, sp_k]
    
    # attribute contribution of species given tau (eq 3d)
    v_k <- asb_sp_w[k, sp_k]/a_k[, 1]
    
    # diversity of order 0 (eq 4c)
    if (0 %in% q) {
      asb_FD_Hill[k, "FD_q0"] <- sum(v_k)
    }
    
    # diversity of order 1 (eq 4d)
    if (1 %in% q) {
      asb_FD_Hill[k, "FD_q1"] <- exp(sum(- v_k*a_k/nplus_k*log(a_k/nplus_k)))
    }
    
    # diversity of order 2 (eq 4e)
    if (2 %in% q) {
      asb_FD_Hill[k, "FD_q2"] <- 1 / (sum(v_k*(a_k/nplus_k)^2))
    }
    
    
  }# end of k
  
  
  #### outputs ####
  
  # indices values
  res <- asb_FD_Hill
  
  # details if required
  if (details_returned == TRUE) {
    res <- list(asb_FD_Hill = asb_FD_Hill,
                tau_dist = tau_dist,
                details= list( asb_totw = asb_totw,
                               asb_sp_relw = asb_sp_relw)
    )
  }
  
  # returning
  return(res)
  
  
} # end of function
