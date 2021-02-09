# ------------------------------------------------------------------------------
# Function to compute  functional beta-diversity indices based on Hill numbers
# applied to distance between species
#
# Authors: Sébastien Villéger and Camille Magneville
#
# ------------------------------------------------------------------------------


#'Compute functional beta-diversity indices based on Hill numbers 
#'
#'Compute functional beta-diversity indices based on Hill numbers applied to
#'distance between species following the framework from Chao et al. 2019,
#'Ecological Monographs (89:e01343), DOI: 10.1002/ecm.1343). FD is computed
#'applying the special case where function 'f' in equation 3c is
#'linear: f(dij(tau)) = dij(tau)/tau, hence f(0) = 0 and f(tau) = 1.
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
#'  \item \emph{asb_FDbeta} a dataframe with a row for each pair of assemblages
#'  (names in 2 first columns, as in \strong{asb_sp_w}) and beta-diversity
#'  for each value of q in other column(s)
#'
#'  \item if \strong{details_returned} turned to TRUE a list \emph{details} with
#'  \itemize{
#'  \item \emph{asb_FDalpha} a dataframe with mean alpha diversity of each pair
#'  of assemblages (rows) and values of q (columns)
#'  \item \emph{asb_FDgamma} a dataframe with gamma diversity of each pair
#'  of assemblages (rows) and values of q (columns)
#'  }
#'  }
#'
#'@note when q=1 Jaccard-like and Sorensen-like beta-diversity are identical.
#' FD computed with tau='min' is equivalent to Hill number taxonomic beta
#' diversity. If tau='min' and there are species with null distance, tau is
#' set to the minimum non-null value and a warning message is displayed.
#' Indices values are stored as \emph{dist} objects to optimize memory.
#'
#'@examples
#' # Load Species*Traits dataframe:
#' data("sp_tr_fruits", package = "mFD")
#' # Load Assemblages*Species dataframe:      
#' data("asb_sp_w_fruits", package = "mFD")   
#' # Compute functional distance 
#' sp_dist_fruits <- mFD::funct.dist(sp_tr = sp_tr_fruits,         
#'  tr_cat       = tr_cat_fruits,   
#'  dist_metric  = "kgower",         
#'  scaling      = "scaledBYrange",  
#'  stop_if_NA   = TRUE)
#' # Compute beta functional hill indices:
#' beta.fd.hill(asb_sp_w = asb_sp_w_fruits, sp_dist = sp_dist_fruits, 
#'  q = c(0,1,2), tau = "mean",
#'  beta_type = "Jaccard", check_input = TRUE, details_returned = TRUE)
#'
#'@export

beta.fd.hill <- function(asb_sp_w,
                         sp_dist,
                         q = c(0,1,2),
                         tau = "mean",
                         beta_type = "Jaccard",
                         check_input = TRUE,
                         details_returned = TRUE) {



  #  distance between species stored in a matrix  ####
  sp_sp_dist <- sp_dist

  if (is.matrix(sp_sp_dist) == FALSE) {
    sp_sp_dist <- as.matrix(sp_sp_dist)
  }


  ## check_inputs if required #####
  if (check_input == TRUE) {

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

    # special case of null distance outside diagonal:
    if (tau_dist == 0) {
      tau_dist <- min(sp_dist[sp_dist != 0])
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


  # dissimilarity between assemblages ####

  # list to store diversity values
  beta_fd_q <- list()
  malpha_fd_q <- list()
  gamma_fd_q <- list()

  # matrices to store diversity values of order q
  mat_res <- matrix(NA, asb_nb, asb_nb, dimnames = list(asb_nm, asb_nm))
  if (0 %in% q) {
    beta_fd_q$q0 <- mat_res
    gamma_fd_q$q0 <- mat_res
    malpha_fd_q$q0 <- mat_res
  }
  if (1 %in% q) {
    beta_fd_q$q1 <- mat_res
    gamma_fd_q$q1 <- mat_res
    malpha_fd_q$q1 <- mat_res
  }
  if (2 %in% q) {
    beta_fd_q$q2 <- mat_res
    gamma_fd_q$q2 <- mat_res
    malpha_fd_q$q2 <- mat_res
  }


  # combinations of assemblages
  asb_pairs <- t (utils::combn(asb_nm, 2))
  asb_pairs_nb <- nrow(asb_pairs)
  colnames(asb_pairs) <- paste0("asb.", 1:2)

  # loop on pairs of assemblages
  for (x in 1:asb_pairs_nb) {

    # names of assemblages in the pair x:
    asb_nm_x<-asb_pairs[x,]


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
    if (0 %in% q) {

      # hence sum of species attribute contribution depends on their occurrence:
      x_malpha_q0 <- sum(x_sp_vip*x_sp_01)/2

      # gamma diversity (eq 6a)
      x_gamma_q0 <- sum(x_sp_vip)

      # beta Jaccard or Sorensen
      if (beta_type == "Sorensen") {
        x_beta_q0 <- (x_gamma_q0/x_malpha_q0) - 1
      }
      if (beta_type == "Jaccard") {
        x_beta_q0 <- (1 - (x_malpha_q0/x_gamma_q0))/(0.5)
      }

      # storing values
      malpha_fd_q$q0[asb_nm_x[2], asb_nm_x[1]] <- x_malpha_q0
      gamma_fd_q$q0[asb_nm_x[2], asb_nm_x[1]] <- x_gamma_q0
      beta_fd_q$q0[asb_nm_x[2], asb_nm_x[1]] <- x_beta_q0
    } # end of q=0

    # q=1 -----
    if (1 %in% q) {

      # alpha diversity (eq 7b) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      x_malpha_q1 <- 0.5 * exp((-1)* sum(x_sp_vip*(x_sp_aik/x_npp)*
                                          log(x_sp_aik/x_npp), na.rm = T))

      # gamma diversity (eq 6b)
      x_gamma_q1 <- exp((-1) * sum(x_sp_vip*(x_sp_aip/x_npp)*
                                    log(x_sp_aip/x_npp)))

      # beta Jaccard or Sorensen are identical
      x_beta_q1 <- log(x_gamma_q1/x_malpha_q1)/log(2)

      # storing values
      malpha_fd_q$q1[asb_nm_x[2], asb_nm_x[1]] <- x_malpha_q1
      gamma_fd_q$q1[asb_nm_x[2], asb_nm_x[1]] <- x_gamma_q1
      beta_fd_q$q1[asb_nm_x[2], asb_nm_x[1]] <- x_beta_q1

    } # end of q=1



    # q=2 ----
    if (2 %in% q)
    {
      # alpha diversity (eq 7a) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      x_malpha_q2 <- 0.5 / (sum(x_sp_vip * ((x_sp_aik / x_npp)^2)))

      # gamma diversity (eq 6a)
      x_gamma_q2 <-1 / (sum(x_sp_vip * ((x_sp_aip / x_npp)^2)))

      # beta Jaccard or Sorensen
      if (beta_type == "Sorensen") {
        x_beta_q2 <- (1 - (x_malpha_q2 / x_gamma_q2)) / (0.5)
      }
      if (beta_type == "Jaccard") {
        x_beta_q2 <- (x_gamma_q2 / x_malpha_q2) - 1
      }

      # storing values
      malpha_fd_q$q2[asb_nm_x[2], asb_nm_x[1]] <- x_malpha_q2
      gamma_fd_q$q2[asb_nm_x[2], asb_nm_x[1]] <- x_gamma_q2
      beta_fd_q$q2[asb_nm_x[2], asb_nm_x[1]] <- x_beta_q2
    } # end of q=2


  } # end of loop on pairs

  # matrix with indices values as dist objects
  malpha_fd_q <- lapply(malpha_fd_q, stats::as.dist)
  gamma_fd_q <- lapply(gamma_fd_q, stats::as.dist)
  beta_fd_q <- lapply(beta_fd_q, stats::as.dist)

  # returning outputs
  res <- beta_fd_q

  if (details_returned == TRUE) {
    res <- list(beta_fd_q = beta_fd_q,
              details = list(malpha_fd_q = malpha_fd_q, gamma_fd_q = gamma_fd_q))
  }
  return(res)


} # end of function
