# ------------------------------------------------------------------------------
# Function to compute  functional beta-diversity indices based on Hill numbers
# applied to distance between species
#
# Authors: Sébastien Villéger and Camille Magneville
#
# ------------------------------------------------------------------------------


#'Compute functional beta-diversity indices based on Hill numbers
#' applied to distance between species following the framework from Chao et al.
#' 2019, Ecological Monographs (89:e01343), DOI: 10.1002/ecm.1343)
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
#'@param beta_type a \strong{character string} with name of framework used for computing
#' beta-diversity, either 'Jaccard' (default) or 'Sorensen'.
#'
#'@param check.input a \strong{logical value} defining whether inputs are checked before
#'  computation of indices. Possible error messages will thus may be more
#'  understandable for the user than R error messages. Default: check.input =
#'  TRUE.
#'
#'@param store.details a \strong{logical value} indicating whether the user want to store
#' values used for computing indices (see list below)
#'
#'@return a list with: \itemize{
#'
#'  \item \emph{asb_FDbeta} a dataframe with a row for each pair of assemblages
#'  (names in 2 first columns, as in \strong{asb_sp_w}) and beta-diversity
#'  for each value of q in other column(s)
#'
#'  \item if \strong{store.details} turned to TRUE a list \emph{details} with
#'  \itemize{
#'  \item \emph{asb_FDalpha} a dataframe with mean alpha diversity of each pair
#'  of assemblages (rows) and values of q (columns)
#'  \item \emph{asb_FDgamma} a dataframe with gamma diversity of each pair
#'  of assemblages (rows) and values of q (columns)
#'  }
#'  }
#
#'@note when q=1 Jaccard-like and Sorensen-like beta-diversity are identical.
#' FD computed with tau='min' is equivalent to Hill number taxonomic beta
#' diversity.
#'
#'@examples
#' load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#' sp_tr <- sp_tr[, -c(6:8)]
#' load(system.file("extdata", "asb_sp_w_fruits", package = "mFD"))
#' asb_sp_w <- as.matrix(asb_sp_w)
#' sp_dist <- cluster::daisy(sp_tr, metric = "gower")
#'  beta.fd.hill(asb_sp_w, sp_dist, q = c(0,1,2), tau = "mean",
#'  beta_type = "Jaccard", check.input = TRUE, store.details = TRUE)
#'
#'@export

beta.fd.hill <- function(asb_sp_w,
                         sp_dist,
                         q = c(0,1,2),
                         tau = "mean",
                         beta_type = "Jaccard",
                         check.input = TRUE,
                         store.details = TRUE) {
  
  
  
  #  distance between species stored in a matrix  ####
  sp_sp_dist <- sp_dist
  
  if (is.matrix(sp_sp_dist) == FALSE) {
    sp_sp_dist <- as.matrix(sp_sp_dist)
  }
  
  
  ## check inputs if required #####
  if (check.input == TRUE) {
    
    if (any(is.na(sp_dist))) {
      stop("Error: The species distances matrix contains NA. Please check.")
    }
    if (is.null(rownames(sp_sp_dist))) {
      stop("Error: No row names provided in species distance matrix.
             Please add species names as row names.")
    }
    if (any(is.na(asb_sp_w))) {
      stop("Error: The species*weights matrix contains NA. Please check.")
    }
    if (is.null(rownames(asb_sp_w))) {
      stop("Error: No row names provided in species*weights dataframe.
             Please add assemblages names as row names.")
    }
    if (is.null(colnames(asb_sp_w))) {
      stop("Error: No column names provided in species*weights dataframe.
             Please add species names as column names.")
    }
    if (any(! (colnames(asb_sp_w) %in% rownames(sp_sp_dist)))) {
      stop(paste("Error: Mismatch between names in species*weight and
                   species distances matrix. Please check."))
    }
    
    if (any(! q %in% c(0, 1, 2))) {
      stop(paste("Error: q should be 0, 1 and/or 2.
                  Please check."))
    }
    
    if (any(! tau %in% c("min", "mean", "max"))) {
      stop(paste("Error: tau should be 'mean' or 'max'. Please check."))
    }
    
    if (any(! beta_type %in% c("Jaccard", "Sorensen"))) {
      stop(paste("Error: beta_type should be 'Jaccard' or 'Sorensen'. Please check."))
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
    
    isnum_vect <- sapply(asb_sp_w, is.numeric)
    
    if (FALSE %in% isnum_vect) {
      stop("Error: The 'asp_sp_w' dataframe must only contain numeric values. Please convert values")
    }
    
  }# end of checking inputs
  
  #  preliminary operations ####
  
  # ensuring species are in the same order in  (asb_sp_w)]
  
  # names and number of assemblages
  asb_sp_w <- as.matrix(asb_sp_w)
  asb_nm <- row.names(asb_sp_w)
  asb_nb <-length(asb_nm)
  
  # computing total weight per assemblage  ----
  asb_totw <- apply(asb_sp_w, 1, sum)
  if(any(asb_totw == 0)) {
    stop(paste("Error: all assemblages should contain at least one species.
               Please check."))
  }
  
  
  # computing tau as mean or max on distances ----
  tau_dist <- NULL
  
  if(tau == "min") {
    tau_dist <- min(sp_dist)
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
  
  
  
  # dissimilarity between assemblages ####
  
  # combinations of assemblages
  asb_pairs <- t (utils::combn(asb_nm, 2))
  asb_pairs_nb <- nrow(asb_pairs)
  colnames(asb_pairs) <- paste0("asb.", 1:2)
  
  # dataframe to store diversity values of order q
  asb_FDgamma <- matrix(NA, asb_pairs_nb,length(q),
                        dimnames = list(NULL, paste0("q", q)))
  asb_FDalpha <- matrix(NA, asb_pairs_nb,length(q),
                        dimnames = list(NULL, paste0("q", q)))
  asb_FDbeta <- matrix(NA, asb_pairs_nb,length(q),
                       dimnames = list(NULL, paste0("q", q)))
  
  # loop on pairs of assemblages
  for (x in 1:asb_pairs_nb)
  {
    
    # computing core variables for the pair of assemblages ----
    # notations as in Chao et al 2019, page 16, bottom right (with p for +)
    
    # weights of species (rows) in the 2 assemblages (columns)
    # (nik, bottom right p16)
    x_nik <- t(asb_sp_w[asb_pairs[x, ], ])
    
    # total weight of species in the 2 assemblages
    x_npp <- sum(x_nik)
    
    # total weight of each species among the 2 assemblages
    x_nip <- apply(x_nik, 1, sum)
    
    # keeping only weight and distance of species present in pair of assemblages
    x_sp <- names(which(x_nip > 0))
    x_nip <- x_nip[x_sp]
    x_nik <- x_nik[x_sp, ]
    x_sp_dist <- dij_tau[x_sp, x_sp]
    
    # weight of functionally distinct group of species (aik, ai+ and vi+)
    x_sp_aik <- (1 - x_sp_dist / tau_dist) %*% x_nik
    x_sp_aip <- apply(x_sp_aik, 1, sum)
    x_sp_vip <- x_nip / x_sp_aip
    
    # species occurrences
    x_sp_01 <- x_sp_aik
    x_sp_01[which(x_sp_01 > 0)] <- 1
    
    # computing alpha, gamma and beta diversity according to levels of q ----
    
    # q=0 ----
    if (0 %in% q)
    {
      # alpha diversity (eq 7a) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      asb_FDalpha[x, "q0"] <- sum(x_sp_vip * x_sp_01) / 2
      
      # gamma diversity (eq 6a)
      asb_FDgamma[x, "q0"] <- sum(x_sp_vip)
      
      # beta Jaccard or Sorensen
      if (beta_type == "Sorensen") {
        asb_FDbeta[x, "q0"] <- (asb_FDgamma[x, "q0"] / asb_FDalpha[x, "q0"]) - 1
      }
      if (beta_type == "Jaccard") {
        asb_FDbeta[x, "q0"] <- (1 - (asb_FDalpha[x, "q0"] / asb_FDgamma[x, "q0"])) / (0.5)
      }
    } # end of q=0
    
    # q=1 -----
    if (1 %in% q)
    {
      # alpha diversity (eq 7b) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      asb_FDalpha[x, "q1"] <- 0.5 * exp((-1) * sum(x_sp_vip * (x_sp_aik / x_npp) *
                                                     log(x_sp_aik / x_npp), na.rm = T))
      
      # gamma diversity (eq 6b)
      asb_FDgamma[x, "q1"] <- exp((-1) * sum(x_sp_vip * (x_sp_aip / x_npp) *
                                               log(x_sp_aip / x_npp)))
      
      
      # beta Jaccard or Sorensen are identical
      asb_FDbeta[x, "q1"] <- log(asb_FDgamma[x, "q1"] / asb_FDalpha[x, "q1"]) / log(2)
      
    } # end of q=1
    
    
    
    # q=2 ----
    if (2 %in% q)
    {
      # alpha diversity (eq 7a) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      asb_FDalpha[x, "q2"] <- 0.5 / ( sum(x_sp_vip*((x_sp_aik/x_npp)^2)))
      
      # gamma diversity (eq 6a)
      asb_FDgamma[x, "q2"] <- 1 / ( sum(x_sp_vip*((x_sp_aip/x_npp)^2)))
      
      # beta Jaccard or Sorensen
      if (beta_type == "Sorensen") {
        asb_FDbeta[x, "q2"] <- (1 - (asb_FDalpha[x, "q2"]/asb_FDgamma[x, "q2"]))/(0.5)
      }
      if (beta_type == "Jaccard") {
        asb_FDbeta[x,"q2"] <- (asb_FDgamma[x, "q2"]/asb_FDalpha[x, "q2"]) - 1
      }
    } # end of q=2
    
    
  } # end of loop on pairs
  
  # returning outputs
  asb_FDbeta <- data.frame(asb_pairs, asb_FDbeta)
  res <- asb_FDbeta
  
  if (store.details == TRUE)
  {
    res <- list(asb_FDbeta = asb_FDbeta, details = list(
      asb_FDalpha = data.frame(asb_pairs, asb_FDalpha),
      asb_FDgamma = data.frame(asb_pairs, asb_FDgamma)
    )
    )
  }
  return(res)
  
  
} # end of function
