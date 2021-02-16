#' Compute a set of alpha functional indices for a set of assemblages
#'
#' This function computes a set of multidimensional space based indices of alpha
#' functional diversity. The user can choose which functional indices to
#' compute.
#'
#' @param sp_faxes_coord a \strong{matrix} of species coordinates in a chosen
#'   functional space. Species coordinates have been retrieved thanks to
#'   \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param asb_sp_w a \strong{matrix} linking weight of species (columns) and a
#'   set of assemblages (rows).
#'
#' @param ind_vect a \strong{vector} of character string of the name of
#'   functional indices to compute. \strong{Indices names must be written in
#'   lower case letters}. Possible indices to compute are: 'fide', fdis',
#'   'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 'fori' and 'fspe'. Default: all the
#'   indices are computed.
#'
#' @param scaling a \strong{logical value} indicating if scaling is to be done
#'   (TRUE) or not (FALSE) on functional indices. Scaling is used to be able to
#'   compare indices values between assemblages. Default: scaling = TRUE.
#'
#' @param check_input a \strong{logical value} defining whether inputs are
#'   checked before computation of indices. Possible error messages will thus
#'   may be more understandable for the user than R error messages. Default:
#'   check_input = TRUE.
#'
#' @param details_returned a \strong{logical value} indicating whether the user
#'   want to store details. Details are used in graphical functions and thus
#'   must be kept if the user want to have graphical outputs for the computed
#'   indices.
#'
#' @return a list with: \itemize{
#'
#'  \item \emph{functional_diversity_indices} matrix containing indices values
#'  (columns) for each assemblage (rows)
#'
#'  \item \emph{details} list: a \strong{asb_sp_occ} data.frame of species
#'  occurrences in each assemblage ; a \strong{asb_sp_relatw} matrix of
#'  relative weight of species in each assemblage ; a \strong{sp_coord_all_asb}
#'  list of matrices of species coordinates along functional axes for species
#'  present in each assemblage ; a \strong{vert_nm_all_asb} list of vectors of
#'  species names being vertices of the convex hull for each assemblage ; a
#'  \strong{mst_all_asb} list of data.frames summarizing link between species in
#'  the minimum spanning tree of each assemblage ; a
#'  \strong{grav_center_vert_coord_all_asb} list of vectors of coordinates of
#'  the vertices gravity center for each assemblage ; a
#'  \strong{mean_dtogravcenter_all_asb} list of vectors containing mean distance
#'  to the species gravity center for each assemblage ; a
#'  \strong{dist_gravcenter_global_pool} vector containing the distance of each
#'  species to the gravity center of all species from the global pool ; a
#'  \strong{dist_nn_global_pool} data.frame showing the distances of each
#'  species from the global pool to its nearest neighbor ; a
#'  \strong{nm_nn_all_asb} data.frame containing the name of each nearest
#'  neighbor of each species present in a given assemblage ; a
#'  \strong{dist_nn_all_asb} data.frame containing distance of each species
#'  present in a given assemblage to its nearest neighbor.}
#'  
#' @author Camille Magneville & Sebastien Villeger
#'
#' @examples
#' # Load Species*Traits dataframe:
#' data('fruits_traits', package = 'mFD')
#' 
#' # Load Assemblages*Species dataframe:      
#' data('baskets_fruits_weights', package = 'mFD')
#' 
#' # Load Traits categories dataframe:
#' data('fruits_traits_cat', package = 'mFD')  
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
#' # Compute functional spaces quality to retrieve species coordinates matrix:
#' fspaces_quality_fruits <- mFD::quality.fspaces(
#'   sp_dist             = sp_dist_fruits, 
#'   maxdim_pcoa         = 10,
#'   deviation_weighting = 'absolute',
#'   fdist_scaling       = FALSE,
#'   fdendro             = 'average')
#'   
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#' 
#' # Compute alpha diversity indices
#' alpha_fd_indices_fruits <- mFD::alpha.fd.multidim(
#'   sp_faxes_coord   = sp_faxes_coord_fruits[ , c('PC1', 'PC2', 'PC3', 'PC4')],
#'   asb_sp_w         = baskets_fruits_weights, 
#'   ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
#'                        'fori', 'fspe'),
#'   scaling          = TRUE, 
#'   check_input      = TRUE, 
#'   details_returned = TRUE)
#'   
#' # Retrieve alpha diversity indices table
#' fd_ind_values_fruits <- alpha_fd_indices_fruits$functional_diversity_indices
#' fd_ind_values_fruits
#'
#' @importFrom stats dist
#' @importFrom geometry convhulln
#' 
#' @export


alpha.fd.multidim <- function(sp_faxes_coord, asb_sp_w, 
                              ind_vect = c("fide", "fdis", "fmpd", "fnnd", 
                                           "feve", "fric", "fdiv", "fori", 
                                           "fspe"), scaling = TRUE, 
                              check_input = TRUE, details_returned = TRUE) {
  
  
  ## check_input if asked:
  if (check_input == TRUE) {
    
    if (any(!ind_vect %in% c("fide", "fdis", "fmpd", 
                             "fnnd", "feve", "fric", "fdiv", "fori", 
                             "fspe"))) {
      stop("Error: Indices names are wrong.
           Please check and rewright them following the function help.")
    }
    
    if (!is.matrix(sp_faxes_coord)) {
      stop("Error: species coordinates on functional axes should be provided as
      a matrix. Please check.")
    }
    
    if (any(is.na(sp_faxes_coord))) {
      stop("Error: The species*coordinates matrix contains NA. Please check.")
    }
    
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    
    if (is.matrix(asb_sp_w) == FALSE) {
      stop("Error: 'asb_sp_w' must be a matrix")
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
      stop("Error: The 'asp_sp_w' dataframe must only contain numeric values. 
           Please convert values")
    }
    
    if (any(!(colnames(asb_sp_w) %in% rownames(sp_faxes_coord)))) {
      stop(paste("Error: Mismatch between names in species*weight and
                   species*coordinates matrix. Please check."))
    }
    
    # Add a stop if some species do not belong to any
    # assemblage:
    if (any(apply(asb_sp_w, 2, sum)) == 0) {
      stop("Error: Some species are absent from all assemblages.")
    }
    # Add a stop if some asb do not contain species:
    if (any(apply(asb_sp_w, 1, sum)) == 0) {
      stop("Error: Some assemblages do not contain species.")
    }
    
    # Add a stop if there is a negative value in the
    # occurrence dataframe:
    if (any(asb_sp_w < 0)) {
      stop("Error: The species*weight dataframe should not contain negative 
          values.Please check.")
    }
    
  }
  
  ## create outputs with NULL objects so work even if
  ## not all ind computed:
  asb_sp_occ <- NULL
  asb_sp_relatw <- NULL
  sp_faxes_coord_list <- NULL
  vert_nm_list <- NULL
  mst_list <- NULL
  grav_center_vert_coord_list <- NULL
  mean_dtogravcenter_list <- NULL
  dtogravcenter_global_pool_list <- NULL
  nm_nn_global_pool_list <- NULL
  dist_nn_global_pool_list <- NULL
  nm_nn_all_asb <- NULL
  dist_nn_asb_list <- NULL
  grav_center_global_pool <- NULL
  centroid_global_pool <- NULL
  
  ## build matrices and vectors to store values to
  ## compute indices or ...  ... get outputs (used in
  ## graphical functions):
  
  sp_faxes_coord <- as.data.frame(sp_faxes_coord)
  # If asb_sp_w is not a df, convert it:
  asb_sp_w <- as.data.frame(asb_sp_w)
  
  # build a matrix to store indices values for each
  # assemblage:
  if ("fide" %in% ind_vect) {
    asb_ind_values <- matrix(NA, nrow(asb_sp_w), 
                             (length(ind_vect) - 1))
  } else {
    asb_ind_values <- matrix(NA, nrow(asb_sp_w), 
                             length(ind_vect))
  }
  
  rownames(asb_ind_values) <- rownames(asb_sp_w)
  
  if ("fide" %in% ind_vect) {
    colnames(asb_ind_values) <- ind_vect[!ind_vect %in% 
                                           c("fide")]
  } else {
    colnames(asb_ind_values) <- ind_vect
  }
  # create a list to store vertices identities if
  # fric or fdiv computed:
  vert_nm_all_asb <- list()
  
  # create a list to store coordinates of species per
  # assemblage to use it...  ... in output so it can
  # be then used in graphical functions:
  sp_coord_all_asb <- list()
  
  # create a list to store mst per assemblage if feve
  # computed to use it...  ... in output so it can be
  # then used in graphical functions...  ... created
  # in general loop for graphical issue:
  mst_all_asb <- list()
  
  
  # create a list to store a vector for each
  # assemblage with coordinates ...  ...of the
  # gravity center of vertices if fdiv computed...
  # ... created in general loop for graphical issue:
  grav_center_vert_coord <- list()
  
  
  # create a list to store the mean distance of
  # species of an assemblage to...  ... the gravity
  # center of species present in the assemblage...
  # ... created in general loop for graphical issue:
  mean_dtogravcenter_all_asb <- list()
  
  
  # create a list to store the identity of the
  # nearest neighbor (fnnd) for ...  ... each species
  # in a given assemblage if fnnd is computed...  ...
  # created in general loop for graphical issue:
  nm_nn_all_asb <- list()
  
  # create a matrix to store the distance to the
  # nearest neighbour (fnnd) for ...  ... each
  # species in a given assemblage if fnnd is
  # computed...  ... created in general loop for
  # graphical issue:
  dist_nn_all_asb <- matrix(NA, nrow(asb_sp_w), nrow(sp_faxes_coord))
  rownames(dist_nn_all_asb) <- rownames(asb_sp_w)
  colnames(dist_nn_all_asb) <- rownames(sp_faxes_coord)
  dist_nn_all_asb <- as.data.frame(dist_nn_all_asb)
  
  
  # create a list to store names of nearest neighbor
  # in the global pool ...  ... for each species of a
  # given assemblage:
  nm_nn_global_pool <- list()
  
  # create an empty dataframe that will contain the
  # distance to the nearest...  ... neighbor of the
  # global pool for each species of the global
  # pool... ... created in general loop for graphical
  # issue:
  dist_nn_global_pool <- data.frame()
  
  # create an empty dataframe that will contain the
  # distance to the...  centroid of the global pool
  # for each species of the global pool...  ...
  # created in general loop for graphical issue:
  special_sp_global_pool <- data.frame()
  
  # create of a matrix to store fide values used to
  # compute fdis values:
  fide_asb <- matrix(NA, nrow(asb_sp_w), ncol(sp_faxes_coord))
  rownames(fide_asb) <- rownames(asb_sp_w)
  colnames(fide_asb) <- paste0("fide", sep = "_", 
                               colnames(sp_faxes_coord))
  
  # compute distance between species in the global
  # pool for scaling and ...  ... fori computation
  dist_sp <- as.matrix(stats::dist(sp_faxes_coord, 
                                   method = "euclidean"))
  dist_sp[which(dist_sp == 0)] <- NA
  
  # retrieve a matrix with relative weight of species
  # per assemblage to use it...  ... in output so it
  # can be then used in graphical functions:
  relat_w_all_asb <- asb_sp_w
  relat_w_all_asb <- apply(relat_w_all_asb, 1, function(x) {
    x/sum(x)
  })
  
  # retrieve a matrix of presence/absence of species
  # per assemblage:
  asb_sp_w <- as.matrix(asb_sp_w)
  asb_sp_summary <- asb.sp.summary(asb_sp_w)
  asb_sp_occ <- asb_sp_summary$asb_sp_occ
  
  # retrieve species richness per assemblage and add
  # it to the indices matrix:
  asb_sprichn <- apply(asb_sp_occ, 1, sum)
  asb_sprichn <- as.matrix(asb_sprichn)
  rownames(asb_sprichn) <- rownames(asb_sp_w)
  colnames(asb_sprichn) <- "sp_richn"
  asb_ind_values <- as.data.frame(asb_ind_values)
  asb_sprichn <- as.data.frame(asb_sprichn)
  asb_ind_values <- cbind(asb_sprichn, asb_ind_values)
  
  
  ## compute indices for each assemblage: loop on
  ## assemblages to compute if asked fdis, fmpd, fnnd,
  ## feve, fric, fdiv:
  for (n in (1:nrow(asb_sp_w))) {
    k <- rownames(asb_sp_w)[n]
    
    # retrieve information on species belonging to
    # assemblage named k to... ... compute indices:
    sp_filter <- sp.filter(k, sp_faxes_coord, asb_sp_w)
    sp_faxes_coord_k <- sp_filter$`species coordinates`
    sp_faxes_coord_k <- data.matrix(sp_faxes_coord_k)
    asb_sp_w_k <- sp_filter$`species weight`
    asb_sp_w_k <- data.matrix(asb_sp_w_k)
    asb_sp_relatw_k <- sp_filter$`species relative weight`
    asb_sp_relatw_k <- data.matrix(asb_sp_relatw_k)
    
    
    # fill the matrix of species coordinates for all
    # assemblages:
    sp_coord_all_asb[[paste0("sp_coord", sep = "_", 
                             k)]] <- sp_faxes_coord_k
    
    # fide:
    if ("fide" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal to 
                      one for", 
                      k, sep = ""))
        }
      }
      # compute fide value and store in the fide matrix:
      fide <- fide.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, k, check_input = check_input)
      fide_asb[k, ] <- fide
    }
    
    # fide and fdis:
    if ("fdis" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(round(sum(asb_sp_relatw_k), 
                        10), 10) != 1) {
          stop(paste0("Error: the sum of relative weights is not equal to 
                      one for", 
                      sep = " ", k))
        }
      }
      # compute fide value and store in the fide matrix:
      fide <- fide.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, k, check_input = check_input)
      fide_asb[k, ] <- fide
      fdis <- fdis.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, fide_asb = NULL, 
                               k, check_input = check_input)
      # fill the matrix of indices if no scaling: scale
      # and fill the matrix of indices if scaling:
      if (scaling == TRUE) {
        # compute distance between species in the global
        # pool for scaling:
        fdis <- fdis/(max(dist_sp, na.rm = TRUE)/2)
        asb_ind_values[k, "fdis"] <- fdis
      }
      asb_ind_values[k, "fdis"] <- fdis
    }
    
    # fric:
    if ("fric" %in% ind_vect) {
      if (check_input == TRUE) {
        if (nrow(sp_faxes_coord_k) < ncol(sp_faxes_coord_k)) {
          stop(paste0("Error: Number of species should strictly be higher than
          the number of axes to compute the convex hull.
          It is not the case for", 
                      sep = " ", k, sep = ".", "Remove this assemblage or
                      decrease the number of functional axes.
          FRic can not be computed here."))
        }
      }
      fric <- fric.computation(sp_faxes_coord_k, 
                               k, check_input = check_input)
      # scale fric value and fill the matrix of indices
      # if scaling:
      if (scaling == TRUE) 
      {
        # applying convhulln function to compute convexhull
        # for all species...  ... with options = 'FA', to
        # compute the general area of the ...  ...
        # functional hull:
        fric_value <- fric$fric
        
        # if FRic value for this asb can be computed (no
        # coplanearity):
        if (!is.na(fric_value)) {
          conv_fa_all <- tryCatch(geometry::convhulln(sp_faxes_coord, 
                                                      option = "FA"), 
                                  error = function(err) {
                                                        "NA"
                                                         })
          fric_value <- fric$fric
          # if convex hull of the gp can be computed:
          if (!is.character(conv_fa_all)) {
            fric_value <- fric_value/conv_fa_all$vol
          }
        }
        
        # if Fric value for this asb can not be computed
        # (coplanearity):
        if (is.na(fric_value)) {
          fric_value <- NA
        }
        
        # if convex hull of the gp can not be computed:
        if (is.character(conv_fa_all)) {
          fric_value <- NA
        }
        
      }  # end if scaling == TRUE
      
      asb_ind_values[k, "fric"] <- fric_value
      vert_nm_all_asb[[paste0("vert_nm", sep = "_", 
                              k)]] <- fric$vertices_nm
    }
    
    # fdiv:
    if ("fdiv" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal to 
                      one for", 
                      k, sep = " "))
        }
      }
      # retrieve vert_nm values if fric computed before:
      if ("fric" %in% ind_vect) {
        vert_nm <- fric$vertices_nm
      } else {
        vert_nm <- NULL
      }
      # compute fdiv value:
      if (is.character(vert_nm)) {
        fdiv <- fdiv.computation(sp_faxes_coord_k, 
                                 asb_sp_relatw_k, vert_nm = vert_nm, 
                                 k, check_input = check_input)
        # scale and fill the matrix if scaling:
        fdiv_value <- fdiv$fdiv
        # fill the matrix:
        asb_ind_values[k, "fdiv"] <- fdiv_value
        # fill the vertices list if no fric computed:
        if (!"fric" %in% ind_vect) {
          vert_nm_all_asb[[paste0("vert_nm", 
                                  sep = "_", k)]] <- vert_nm
        }
        # fill the list of gravity center coordinates for
        # each assemblage:
        grav_center_vert_coord[[paste0("grav_center_vert_coord", 
                                       sep = "_", k)]] <- fdiv$details$B_coord
        # fill the list of mean distance to gravity center
        # of species present...  ... in each assemblage:
        mean_dtogravcenter_all_asb[[paste0("mean_dist_to_sp_gravcent_asb", 
                                           sep = "_", 
                                           k)]] <- fdiv$details$mean_dtoB
      }
    }
    if (!is.character(vert_nm)) {
      fdiv_value <- NA
      asb_ind_values[k, "fdiv"] <- fdiv_value
      # fill the list of gravity center coordinates for
      # each assemblage:
      grav_center_vert_coord[[paste0("grav_center_vert_coord", 
                                     sep = "_", k)]] <- NA
      # fill the list of mean distance to gravity center
      # of species present...  ... in each assemblage:
      mean_dtogravcenter_all_asb[[paste0("mean_dist_to_sp_gravcent_asb", 
                                         sep = "_", k)]] <- NA
    }
    # feve:
    if ("feve" %in% ind_vect) {
      if (check_input == TRUE) {
        if (nrow(sp_faxes_coord_k) < 3) {
          stop("Error: there must be at least 3 species to compute feve.")
        }
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal 
                      to one for", 
                      k, sep = " "))
        }
      }
      # compute feve value and fill the matrix of indices
      # if no scaling:
      feve <- feve.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, k, check_input = check_input)
      feve_value <- feve$feve
      asb_ind_values[k, "feve"] <- feve_value
      # get the mst for the assemblage and fill the mst
      # list of all assemblages:
      feve_mst <- feve$mst
      mst_all_asb[[paste0("mst", sep = "_", k)]] <- feve_mst
    }
    
    # fmpd:
    if ("fmpd" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal 
                      to one for", 
                      k, sep = ""))
        }
      }
      # compute fmpd value and fill the matrix of indices
      # if no scaling:
      fmpd <- fmpd.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, k, check_input = check_input)
      # scale fmpd value and fill the matrix of indices
      # if scaling:
      if (scaling == TRUE) {
        mean_dist_sp <- apply(dist_sp, 1, mean, 
                              na.rm = TRUE)
        fmpd <- fmpd/(max(mean_dist_sp, na.rm = TRUE))
        asb_ind_values[k, "fmpd"] <- fmpd
      }
      asb_ind_values[k, "fmpd"] <- fmpd
    }
    
    # fnnd:
    if ("fnnd" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 1) {
          stop(paste0("Error: the sum of relative weights is not equal 
                      to one for", 
                      k, sep = ""))
        }
      }
      # compute fmpd value and fill the matrix of indices
      # if no scaling:
      fnnd <- fnnd.computation(asb_sp_relatw_k, 
                               sp_faxes_coord_k, k, check_input = check_input)
      # scale fnnd value and fill the matrix of indices
      # if scaling:
      if (scaling == TRUE) {
        dist_nn <- list()
        for (i in (1:nrow(sp_faxes_coord))) {
          ref_sp <- rownames(sp_faxes_coord)[i]
          dist_nn_sp <- dist.nearneighb(sp_faxes_coord, 
                                        ref_sp)
          dist_nn[ref_sp] <- dist_nn_sp$`distance of the reference species to its nearest neighbour`
        }
        fnnd_value <- fnnd$fnnd
        fnnd_value <- fnnd_value/(max(unlist(dist_nn), 
                                      na.rm = TRUE))
        asb_ind_values[k, "fnnd"] <- fnnd_value
      }
      # fill the matrix with the fnnd value if no
      # scaling:
      fnnd_value <- fnnd$fnnd
      asb_ind_values[k, "fnnd"] <- fnnd_value
      
      # retrieve data to use for outputs: a list for nms
      # (because can have several nn) and a df for
      # distances (because only one nn dist):
      nms_all <- list(fnnd$details$nm_nn_k)
      names(nms_all) <- k
      nm_nn_all_asb[k] <- nms_all
      dist_nn_all_asb[k, 
                      c(rownames(sp_faxes_coord_k))] <- fnnd$details$dist_nn_k
    }
    
    # fori:
    if ("fori" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal to one for", 
                      k, sep = ""))
        }
      }
      # compute the minimal distance to the nearest
      # neighbor for each...  ... species in the global
      # pool:
      dist_nn_global_pool <- apply(dist_sp, 1, 
                                   min, na.rm = TRUE)
      fori <- fori.computation(dist_nn_global_pool, 
                               asb_sp_relatw_k, k, check_input = check_input)
      
      # create a list to store the identity of the nn for
      # each species of...  ...the assemblage:
      nm_nn_gp <- list()
      
      # compute identity to the nearest neighbor in the
      # gp for plots...
      for (i in (1:nrow(sp_faxes_coord_k))) {
        ref_sp <- rownames(sp_faxes_coord_k)[i]
        dist_nn_sp_gp <- dist.nearneighb(sp_faxes_coord, 
                                         ref_sp)
        nms <- list(dist_nn_sp_gp$`nearest neighbour identity`)
        names(nms) <- ref_sp
        nm_nn_gp[ref_sp] <- nms
      }
      
      # retrieve data to use for outputs: a list for nms
      # (because can have several nn) and a df for
      # distances (because only one nn dist):
      nms_all_gp <- list(nm_nn_gp)
      names(nms_all_gp) <- k
      nm_nn_global_pool[k] <- nms_all_gp
      dist_nn_global_pool <- as.data.frame(dist_nn_global_pool)
      
      # scale fori value and fill the matrix of indices
      # if scaling:
      if (scaling == TRUE) {
        fori <- fori/max(dist_nn_global_pool)
        asb_ind_values[k, "fori"] <- fori
      }
      asb_ind_values[k, "fori"] <- fori
    }
    
    # fspe:
    if ("fspe" %in% ind_vect) {
      # check relative weights sum equals to 1:
      if (check_input == TRUE) {
        if (round(sum(asb_sp_relatw_k), 10) != 
            1) {
          stop(paste0("Error: the sum of relative weights is not equal 
                      to one for", 
                      k, sep = ""))
        }
      }
      # compute specialization of each species in the
      # global pool i.e. ...  ... distances to the
      # centroid of the global pool of species:
      # coordinates of the gravity center of the
      # vertices:
      centroid_global_pool <- apply(sp_faxes_coord, 
                                    2, mean)
      special_sp_global_pool <- apply(sp_faxes_coord, 
                                      1, function(x) {
                                        (sum((x - centroid_global_pool)^2))^0.5
                                      })
      # compute fspe value:
      fspe <- fspe.computation(asb_sp_relatw_k, 
                               special_sp_global_pool, k, 
                               check_input = check_input)
      if (scaling == TRUE) {
        fspe <- fspe/max(special_sp_global_pool)
        asb_ind_values[k, "fspe"] <- fspe
      }
      asb_ind_values[k, "fspe"] <- fspe
    }
    print(paste0(k, sep = " ", "done"))
    print(paste0(k, sep = " ", "done ", round((n/nrow(asb_sp_w)) * 
                                                100, digits = 1), "%"))
  }  # end loop on assemblages
  
  
  ## create a matrix linking fide values and other
  ## indices:
  asb_ind_values_all <- cbind(asb_ind_values, fide_asb)
  
  ## construct the return list:
  if (details_returned == TRUE) {
    return_list <- list(functional_diversity_indices = asb_ind_values_all, 
          details = list(asb_sp_occ = asb_sp_occ, 
            asb_sp_relatw = relat_w_all_asb, 
            sp_faxes_coord_list = sp_coord_all_asb, 
            vert_nm_list = vert_nm_all_asb, 
            mst_list = mst_all_asb, 
            grav_center_vert_coord_list = grav_center_vert_coord, 
            mean_dtogravcenter_list = mean_dtogravcenter_all_asb, 
            dtogravcenter_global_pool_list = special_sp_global_pool, 
            nm_nn_global_pool_list = nm_nn_global_pool, 
            dist_nn_global_pool_list = dist_nn_global_pool, 
            nm_nn_asb_list = nm_nn_all_asb, dist_nn_asb_list = dist_nn_all_asb, 
            grav_center_global_pool = centroid_global_pool))
  } else {
    return_list <- (functional_diversity_indices = asb_ind_values_all)
  }
  
  return(return_list)
  
}  # end function