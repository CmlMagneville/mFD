#' Compute Functional beta-Diversity indices for pairs of assemblages in a
#' multidimensional space
#'
#' Computes a set of indices of pairwise functional beta-diversity
#' (dissimilarity and its turnover and nestedness-resultant components) based 
#' on overlap between convex hulls in a multidimensional space. For details 
#' about indices formulas see Villeger _et al._ (2013). This functions stands 
#' on \code{\link[betapart]{functional.betapart.core.pairwise}}.
#'
#' @param sp_faxes_coord a matrix with coordinates of species (rows) on
#'  functional axes (columns). Species coordinates have been retrieved thanks 
#'  to \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#' @param asb_sp_occ a matrix with presence/absence (coded as 0/1)
#'  of species (columns) in a set of assemblages (rows).
#'
#' @param check_input a logical value defining whether inputs are
#'   checked before computation of indices. Possible error messages will thus
#'   may be more understandable for the user than R error messages. Default:
#'   `check_input = TRUE`.
#'   
#' @param beta_family a character string for the type of beta-diversity
#'   index to use, 'Jaccard' (default) and/or 'Sorensen'.
#'
#' @param details_returned a logical value indicating whether the user
#'   wants to details_returned. Details are used in the graphical function
#'   \code{beta.multidim.plot} and thus must be kept if the user want to have
#'   graphical outputs for the computed indices.
#'   
#' @param betapart_step a logical value indicating whether the
#'   computation progress should be displayed in the R console. Default: 
#'   `betapart_step = TRUE`.
#'
#' @param betapart_chullopt a A list of two named vectors of character conv1 
#'   and conv2 defining the options that will be used to compute the convex 
#'   hulls (through the options of geometry::convhulln function). For further 
#'   details  see help of 
#'   \code{\link[betapart]{functional.betapart.core.pairwise}}.
#'   Default: `betapart_chullopt = c(conv1 = 'Qt', conv2 = 'QJ')`.
#'
#' @param betapart_para a logical value indicating whether internal
#'  parallelization should be used to compute pairwise dissimilarities. 
#'  Default: `betapart_para = FALSE`.
#'
#' @param betapart_para_opt a list with details about parallelization.
#'   Default value means those parameters are set according to computer
#'   specifications. \code{nc} is the number of cores (default = 4), 
#'   \code{type} is a character string with code of method used 
#'   (default PSOCK), \code{LB} is a boolean specifying whether load-balancing 
#'   is applied (default is TRUE) and \code{size} is a numeric value for number 
#'   of tasks performed at each time (default is 1). See help of
#'   \code{\link[betapart]{functional.betapart.core.pairwise}} for more 
#'   details.
#'
#' @return A list with: \itemize{ \item \emph{pairasb_fbd_indices} a list with 
#' for each
#'   index a \emph{dist} object with values for all pairs of assemblages.
#'   Indices names start with the abbreviation of the type of dissimilarity
#'   index ('jac' for Jaccard-like and 'sor' for Sorensen-like dissimilarity)
#'   and end with abbreviation of component ('diss': dissimilarity, 'turn' its
#'   turnover component, and 'nest' its nestedness-resultant component).
#'   \item if \emph{store_details} is TRUE, \item \emph{details_beta} list:
#'  \strong{inputs} a list with \emph{sp_faxes_coord} and \emph{asb_sp_occ} 
#'  on which indices were computed (required for drawing graphics), 
#'  \strong{pool_vertices} a list of vectors (1 per assemblage) with names of
#'  species being vertices of the convex hull shaping all species; 
#'  \strong{asb_FRic} a vector with volume of the convex hull shaping each
#'  assemblage (relative to volume filled by all species) ;
#'  \strong{asb_vertices} a list of vectors (1 per assemblage) with names of
#'  species being vertices of the convex hull}
#'
#' @note All assemblages should have a number of species strictly higher than
#' the number of functional axes.
#' Computing intersection of convex hulls in space of >5 dimensions is yet
#' impossible with most computers.
#' This function uses R libraries 'betapart' (> =1.5.4) for indices 
#' computation.
#' Indices values are stored as \emph{dist} objects to optimize memory.
#' See below example of how merging distance values in a \emph{dataframe} with
#' \code{dist.to.df}.
#'
#' @references 
#' Villeger _et al._ (2013) Decomposing functional beta-diversity reveals that
#'   low functional beta-diversity is driven by low functional turnover in 
#'   European fish assemblages. _Global Ecology and Biogeography_, **22**, 
#'   671-681.
#'   
#' @author Sebastien Villeger and Camille Magneville
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  # Load Species*Traits dataframe:
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
#'  sp_dist             = sp_dist_fruits, 
#'  maxdim_pcoa         = 10,
#'  deviation_weighting = 'absolute',
#'  fdist_scaling       = FALSE,
#'  fdendro             = 'average')
#'  
#' # Retrieve species coordinates matrix:
#' sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#' 
#' # Get the occurrence dataframe:
#' asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights) 
#' asb_sp_fruits_occ <- asb_sp_fruits_summ$'asb_sp_occ'
#' 
#' # Compute beta diversity indices:
#' beta_fd_fruits <- mFD::beta.fd.multidim(
#'   sp_faxes_coord   = sp_faxes_coord_fruits[, c('PC1', 'PC2', 'PC3', 'PC4')], 
#'   asb_sp_occ       = asb_sp_fruits_occ,
#'   check_input      = TRUE,
#'   beta_family      = c('Jaccard'),
#'   details_returned = TRUE)
#' 
#' # merging pairwise beta-diversity indices in a data.frame
#' dist.to.df(beta_fd_fruits$pairasb_fbd_indices)
#' }


beta.fd.multidim <- function(sp_faxes_coord, 
                             asb_sp_occ, 
                             check_input = TRUE, 
                             beta_family = "Jaccard", 
                             details_returned = TRUE, 
                             betapart_step = TRUE, 
                             betapart_chullopt = list(conv1 = 'Qt', 
                                                      conv2 = 'QJ'),
                             betapart_para = FALSE, 
                             betapart_para_opt = betapart::beta.para.control())
  {
  
  
  ## check_input if asked:
  if (check_input) {
    
    check.sp.faxes.coord(sp_faxes_coord)
    
    check.asb.sp.w.occ(asb_sp_occ)
    
    if (any(!(colnames(asb_sp_occ) %in% rownames(sp_faxes_coord)))) {
      stop("Mismatch in species names between species occurence and species ", 
           "coordinates matrices. Please check.")
    }
    
    if (any(apply(asb_sp_occ, 1, sum) <= ncol(sp_faxes_coord))) {
      stop("All assemblages should have more species than number of ", 
           "functional axes")
    }
    
    
    if (ncol(sp_faxes_coord) > 5) {
      stop("Computing beta functional diversity in a", ncol(sp_faxes_coord), 
           "-dimensions space could exceed computing power of your computer. ", 
           "Consider keeping only five dimensions.")
    }
    
    if (any(!beta_family %in% c("Jaccard", "Sorensen"))) {
      stop("beta diversity index should be 'Jaccard' and/or 'Sorensen'. ", 
           "Please check.")
    }
  }
  
  
  # ensuring species are in the same order
  # in both matrices
  sp_faxes_coord <- sp_faxes_coord[colnames(asb_sp_occ), ]
  
  # calling functional.betapart.core
  # function to compute beta functional
  # diversity = computing convex hulls and
  # their intersections for all pairs of
  # assemblages 2 last parameters are for
  # parallelization
  F_betapart_core <- betapart::functional.betapart.core.pairwise(
    x              = asb_sp_occ,
    traits         = sp_faxes_coord,
    return.details = TRUE,
    convhull.opt   = betapart_chullopt,
    parallel       = betapart_para,
    opt.parallel   = betapart_para_opt,
    progress       = betapart_step)
  
  # list to store dist objects with beta-diversity values
  beta_fd_ind <- list()
  
  # computing functional beta diversity
  # indices for all pairs of assemblages
  # according to the type of index
  # (indices) selected
  if ("Jaccard" %in% beta_family) {
    
    F_beta_jac <- betapart::functional.beta.pair(F_betapart_core, 
                                                 index.family = "jaccard")
    
    beta_fd_ind$jac_diss <- F_beta_jac$funct.beta.jac
    beta_fd_ind$jac_turn <- F_beta_jac$funct.beta.jtu
    beta_fd_ind$jac_nest <- F_beta_jac$funct.beta.jne
    
  }
  
  if ("Sorensen" %in% beta_family) {
    
    F_beta_sor <- betapart::functional.beta.pair(F_betapart_core, 
                                                 index.family = "sorensen")
    
    beta_fd_ind$sor_diss <- F_beta_sor$funct.beta.sor
    beta_fd_ind$sor_turn <- F_beta_sor$funct.beta.sim
    beta_fd_ind$sor_nest <- F_beta_sor$funct.beta.sne
  }
  
  
  # names of species being vertices in each assemblage
  asb_vertices <- lapply(F_betapart_core$details$CH$coord_vertices, 
                         function(x) rownames(x))
  
  # volume of convex hull filled by each assemblage (FRic)
  asb_FRic <- F_betapart_core$details$CH$FRi$"FRi"
  names(asb_FRic) <- row.names(F_betapart_core$details$CH$FRi)
  
  # compute the volume of the convex hull of the species pool
  ch_pool <- tryCatch(geometry::convhulln(sp_faxes_coord, 
                                          option = "FA"), 
                      error = function(err) "NA")
  
  # if convex hull of the species pool computed, scaling FRic values:
  if (!is.character(ch_pool)) {
    
    pool_vertices <- unique(row.names(sp_faxes_coord)[ch_pool$hull])
    asb_FRic <- asb_FRic / ch_pool$vol 
    
  } else {
    
    pool_vertices <- NA
    asb_FRic <- NA
  }
  
  
  ## results to return
  if (details_returned) {
    
    return_list <- list("pairasb_fbd_indices" = beta_fd_ind, 
                        "details" = list(
                          "inputs" = list(
                            "sp_faxes_coord" = sp_faxes_coord, 
                            "asb_sp_occ" = asb_sp_occ), 
                          "pool_vertices" = pool_vertices, 
                          "asb_FRic" = asb_FRic, 
                          "asb_vertices" = asb_vertices))
    
  } else {
    
    return_list <- list("pairasb_fbd_indices" = beta_fd_ind)
  }
  
  return(return_list)
  
}  # end function
