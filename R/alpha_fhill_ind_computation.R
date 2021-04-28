#' Compute Functional alpha-Diversity indices based on Hill Numbers
#'
#' Compute functional alpha diversity applied to distance between species
#' following the framework from Chao _et al._(2019).
#'
#' @param asb_sp_w a matrix with weight of species (columns) in a set
#'   of assemblages (rows). Rows and columns should have names. NA are not
#'   allowed.
#'
#' @param sp_dist a matrix or dist object with distance between
#'   species. Species names should be provided and match those in 'asb_sp_w'. 
#'   NA are not allowed.
#'
#' @param q a vector containing values referring to the order of
#'   diversity to consider, could be 0, 1 and/or 2.
#'
#' @param tau a character string with name of function to apply to
#'   distance matrix (i.e. among all pairs of species) to get the threshold 
#'   used to define 'functionally indistinct set of species'. Could be 'mean'
#'   (default), 'min' or 'max'. If tau = 'min" and there are null values in
#'   \code{sp_dist}, the threshold is the lowest strictly positive value and a
#'   warning message is displayed.
#'
#' @param check_input a logical value indicating whether key features the 
#'   inputs are checked (e.g. class and/or mode of objects, names of rows 
#'   and/or columns, missing values). If an error is detected, a detailed 
#'   message is returned. Default: `check.input = TRUE`.
#'
#' @param details_returned a logical value indicating whether the user
#'   want to store values used for computing indices (see list below)
#'
#' @return A list with: \itemize{
#'
#'  \item \emph{asb_FD_Hill} a matrix containing indices values for each level
#'  of q (columns, named as 'FD_qx') for each assemblage (rows, named as in
#'  \strong{asb_sp_w})
#'  \item \emph{tau_dist} the threshold value applied to distance between
#'  species to compute diversity according to function provided in \strong{tau}
#'
#'  \item if \strong{details_returned} turned to TRUE a list \emph{details} 
#'  with
#'  \itemize{
#'  \item \emph{asb_totw} a vector with total weight of each assemblage
#'  \item \emph{asb_sp_relw} a matrix with relative weight of species in
#'  assemblages
#'  }
#'  }
#'
#' @note FD is computed applying the special case where function 'f' in 
#'   equation 3c is linear:f(dij(tau)) = dij(tau)/tau, hence f(0) = 0 
#'   and f(tau) = 1. FD computed with q=2 and tau = 'max' is equivalent to 
#'   the Rao's quadratic entropy from Ricotta & Szeidl (2009, J Theor Biol). 
#'   FD computed with tau = 'min' is equivalent to Hill number taxonomic 
#'   diversity, thus with q=0 it is species richness (S), with q = 1 it is
#'   exponential of Shannon entropy (H) and with q = 2 it is 1/(1-D) where D is 
#'   Simpson diversity.  Note that even when q=0, weights of species are 
#'   accounted for in FD. Hence to compute FD based only on distance between 
#'   species present in an assemblage (i.e. a richness-like index) , asb_sp_w 
#'   has to contain only species presence/absence coded as 0/1 with q=0 and 
#'   tau=”mean”. If asb_sp_w contains only 0/1 and q>0, it means that all 
#'   species have the same contribution to FD.
#'   
#' @references 
#' Chao _et al._ (2019) An attribute‐diversity approach to functional
#'   diversity, functional beta diversity, and related (dis)similarity 
#'   measures. _Ecological Monographs_, **89**, e01343.
#'
#' @author Sebastien Villeger and Camille Magneville
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load Species*Traits dataframe:
#' data('fruits_traits', package = 'mFD')
#' 
#' # Load Assemblages*Species dataframe:      
#' data('baskets_fruits_weights', package = 'mFD') 
#'   
#' # Compute functional distance 
#' sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
#'                                   tr_cat        = fruits_traits_cat,
#'                                   metric        = "gower",
#'                                   scale_euclid  = "scale_center",
#'                                   ordinal_var   = "classic",
#'                                   weight_type   = "equal",
#'                                   stop_if_NA    = TRUE)
#' 
#' # Compute alpha fd hill indices:
#' alpha.fd.hill(
#'    asb_sp_w         = baskets_fruits_weights, 
#'    sp_dist          = sp_dist_fruits, 
#'    q                = c(0, 1, 2),
#'    tau              = 'mean', 
#'    check_input      = TRUE, 
#'    details_returned = TRUE)
#' }

alpha.fd.hill <- function(asb_sp_w, sp_dist, 
                          q = c(0, 1, 2), tau = "mean", check_input = TRUE, 
                          details_returned = TRUE) {
  
  
  #### distance between species stored in a matrix ####
  
  sp_sp_dist <- sp_dist
  
  if (!is.matrix(sp_sp_dist)) {
    sp_sp_dist <- as.matrix(sp_sp_dist)
  }
  
  ## check_inputs if required #####
  if (check_input) {
    
    check.asb.sp.w(asb_sp_w)
    
    if (any(is.na(sp_dist))) {
      stop("The species distances matrix contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_sp_dist))) {
      stop("No row names provided in species distances matrix. Please add ", 
           "species names as row names.")
    }
    
    if (any(!(colnames(asb_sp_w) %in% rownames(sp_sp_dist)))) {
      stop("Mismatch between names in species*weight and species distances ", 
           "matrix. Please check.")
    }
    
    if (any(!q %in% c(0, 1, 2))) {
      stop("q should be 0, 1 and/or 2. Please check.")
    }
    
    if (any(!tau %in% c("min", "mean", "max"))) {
      stop("tau should be 'min', 'mean' or 'max'. Please check.")
    }
    
  }
  
  #### preliminary operations ####
  
  # ensuring species are in the same order in both
  # matrices:
  asb_sp_w <- as.matrix(asb_sp_w)
  sp_sp_dist <- sp_sp_dist[colnames(asb_sp_w), colnames(asb_sp_w)]
  
  
  # computing total weight per assemblage and
  # relative weights of species ----
  asb_totw <- apply(asb_sp_w, 1, sum)
  asb_sp_relw <- asb_sp_w / asb_totw
  
  # computing tau as mean or max on distances ----
  tau_dist <- NULL
  
  if (tau == "min") {
    tau_dist <- min(sp_dist)
    
    # special case of null distance outside diagonal
    if (tau_dist == 0) {
      tau_dist <- min(sp_dist[sp_dist != 0])
      cat("Warning: some species has null functional distance,
          'tau' was set to the minimum non-null distance")
    }
  }
  
  if (tau == "mean") {
    tau_dist <- mean(sp_dist)
  }
  
  if (tau == "max") {
    tau_dist <- max(sp_dist)
  }
  
  # applying tau threshold to distance matrix
  dij_tau <- sp_sp_dist
  dij_tau[which(dij_tau > tau_dist, arr.ind = TRUE)] <- tau_dist
  
  
  #### computing diversity of assemblages ####
  
  # empty matrix to store outputs
  asb_FD_Hill <- matrix(NA, nrow(asb_sp_w), length(q), 
                        dimnames = list(row.names(asb_sp_w), 
                                        paste0("FD_q", q)))
  
  # loop on assemblages (equations id refers to those
  # in Chao et al 2019
  for (k in row.names(asb_sp_w)) {
    
    # total weight (n+ in eq 4a)
    nplus_k <- asb_totw[k]
    
    # species names present in assemblage k
    sp_k <- colnames(asb_sp_w)[which(asb_sp_w[k, ] > 0)]
    
    # f(dij(tau)) with f being linear so dij(tau)/tau
    # (see bottom right of p7)
    f_dij_tau <- dij_tau[sp_k, sp_k] / tau_dist
    
    # 'abundance of species given tau' (eq 3c)
    a_k <- (1 - f_dij_tau) %*% asb_sp_w[k, sp_k]
    
    # attribute contribution of species given tau (eq
    # 3d)
    v_k <- asb_sp_w[k, sp_k] / a_k[ , 1]
    
    # diversity of order 0 (eq 4c)
    if (0 %in% q) {
      asb_FD_Hill[k, "FD_q0"] <- sum(v_k)
    }
    
    # diversity of order 1 (eq 4d)
    if (1 %in% q) {
      asb_FD_Hill[k, "FD_q1"] <- exp(sum(-v_k * 
                                           a_k / nplus_k * log(a_k / nplus_k)))
    }
    
    # diversity of order 2 (eq 4e)
    if (2 %in% q) {
      asb_FD_Hill[k, "FD_q2"] <- 1 / (sum(v_k * (a_k / nplus_k) ^ 2))
    }
  }  # end of k
  
  
  #### outputs ####
  
  # indices values
  res <- asb_FD_Hill
  
  # details if required
  if (details_returned) {
    res <- list(asb_FD_Hill = asb_FD_Hill, tau_dist = tau_dist, 
                details = list(asb_totw = asb_totw, asb_sp_relw = asb_sp_relw))
  }
  
  # returning
  return(res)
}  # end of function
