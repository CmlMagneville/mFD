# Function to compute functional beta-diversity indices in a multidimensional space
#
# Authors:  Sébastien Villéger, Camille Magneville
#
#
# ------------------------------------------------------------------------------


#'Computes functional beta-diversity indices for pairs of assemblages in a
#'multidimensional space.
#'
#'Computes a set of indices of pairwise functional beta-diversity (dissimilarity
#'and its turnover and nestedness-resultant components) based on overlap between
#'convex hulls in a multidimensional space. For details about indices formulas
#'see Villeger et al. 2013, Global Ecology and Biogeography (22:671-681).
#'
#'@param sp_faxes_coord a \strong{matrix} with coordinates of species (rows) on
#'  functional axes (columns). Species coordinates have been retrieved thanks to
#'  \code{\link{tr.cont.fspace}} or \code{\link{quality.fspaces}}.
#'
#'@param asb_sp_occ a \strong{matrix} with presence/absence (coded as 0/1)
#'  of species (columns) in a set of assemblages (rows).
#'
#'@param check.input a \strong{logical value} defining whether inputs are checked before
#'  computation of indices. Possible error messages will thus may be more
#'  understandable for the user than R error messages. Default: check.input =
#'  TRUE.
#'
#'@param beta.family a \strong{character string} for the type of beta-diversity index to
#'  use, 'Jaccard' (default) and/or 'Sorensen'
#'
#'@param store_details a \strong{logical value} indicating whether the user wants to
#'  store details. Details are used in the graphical function
#'  \code{beta.multidim.plot} and thus must be kept if the user want to have
#'  graphical outputs for the computed indices.
#'
#'@param betapart.step a \strong{logical value} indicating whether the computation
#'  progress tracking file "step.fbc.txt" should be created. Setting it to FALSE
#'  will speed up the function. Default is FALSE, and it is automatically turned
#'  to FALSE when 'betapart.para' is TRUE.
#'
#'@param betapart.para a \strong{logical value} indicating whether internal
#'  parallelization should be used to compute pairwise dissimilarities. Default
#'  is FALSE.
#'
#'@param betapart.para.opt a \strong{list} with details about parallelization. Default
#'  value means those parameters are set according to computer specifications.
#'  'nc' is number of cores (default = 4), 'type' is a character string with
#'  code of method used (default PSOCK), 'LB' is a boolean specifying whether
#'  load-balancing is applied (default is TRUE) and 'size' is a numeric value
#'  for number of tasks performed at each time (default is 1). See help of
#'  \code{\link[betapart]{functional.betapart.core}} for more details.
#'
#'@return a list with: \itemize{
#'  \item inputs a list with \item \emph{pairasb_fbd_indices} a dataframe with for all pairs of
#'  assemblages (rows) names of the assemblages in the 2 first columns and
#'  columns with beta functional diversity indices following the additive
#'  decomposition of Villéger et al 2013 : 'F_diss': functional dissimilarity
#'  index, 'F_turn': its turnover component, and 'F_nest': its
#'  nestedness-resultant component, i.e. F_diss = F_turn + F_nest. Those indices
#'  names are preceded by the abbreviation of the type of indices used, 'jac'
#'  for decomposition of Jaccard-like functional dissimilarity and 'sor' for
#'  Sorensen-like dissimilarity.
#'  \item \emph{details_beta} list if \emph{store_details} is TRUE:
#'  \emph{sp_faxes_coord} and \emph{asb_sp_occ} on which indices were computed (convenient for drawing
#'  graphics) ; a \strong{asb_FRic_raw} vector with volume of
#'  the convex hull shaping each assemblage ; a \strong{asb_FRic} vector with volume of
#'  the convex hull shaping each assemblage ; a \strong{asb_vertices} a list of
#'  vectors (1 per assemblage) with names of species being vertices of the
#'  convex hull }
#'
#'@section Notes: \strong{All assemblages should have a number of species strictly
#'  higher than the number of functional axes}. Computing intersection of convex
#'  hulls in space of >5 dimensions is yet impossible with most computers.
#'  This function uses R libraries 'betapart' (>=1.5.2) for indices computation
#'  and 'dendextend' for formatting values in a dataframe.
#'
#'@examples
#' load(system.file("extdata", "sp_faxes_coord_fruits", package = "mFD"))
#' load(system.file("extdata", "asb_sp_w_fruits", package = "mFD"))
#' asb_sp_w <- as.matrix(asb_sp_w)
#' asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w)
#' asb_sp_occ <- asb_sp_summ$asb_sp_occ
#' beta_fd_indices <- mFD::beta.fd.multidim(sp_faxes_coord[, 
#'  c("PC1", "PC2", "PC3", "PC4")], asb_sp_occ,
#'  check.input = TRUE,
#'  beta.family = c("Jaccard"),
#'  store_details = TRUE)
#'
#'@export


beta.fd.multidim <- function(sp_faxes_coord,
                             asb_sp_occ,
                             check.input = TRUE,
                             beta.family = "Jaccard",
                             store_details = TRUE,
                             betapart.step = FALSE,
                             betapart.para = FALSE,
                             betapart.para.opt = betapart::beta.para.control()) {
  
  
  
  ## check input if asked:
  if (check.input == TRUE) {
    
    if (!is.matrix(sp_faxes_coord)) {
      stop("Error: species coordinates on functional axes should be provided as
      a matrix. Please check.")
    }
    if (is.null(rownames(sp_faxes_coord))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    if (is.null(colnames(sp_faxes_coord))) {
      stop("Error: No column names provided in species*coordinates matrix.
             Please add axes names as column names.")
    }
    if (any(is.na(sp_faxes_coord))) {
      stop("Error: The species*coordinates matrix contains NA. Please check.")
    }
    
    if (is.null(rownames(asb_sp_occ))) {
      stop("Error: No row names provided in the species occurence matrix.
             Please add assemblages names as row names.")
    }
    if (is.null(colnames(asb_sp_occ))) {
      stop("Error: No column names provided in the species occurence matrix.
             Please add species names as column names.")
    }
    if (any(is.na(asb_sp_occ))) {
      stop("Error: The species occurence matrix contains NA. Please check.")
    }
    if (any(asb_sp_occ != 0 & asb_sp_occ != 1)) {
      stop("Error: The species occurence matrix should contain only 0 and 1.
           Please check.")
    }
    
    if (any(!(colnames(asb_sp_occ) %in% rownames(sp_faxes_coord)))) {
      stop(paste("Error: Mismatch in species names between species occurence
        and species coordinates matrices. Please check."))
    }
    
    if (any(apply(asb_sp_occ, 1, sum) <= ncol(sp_faxes_coord))) {
      stop(paste("Error: all assemblages should have more species
                   than number of functional axes"))
    }
    
    
    if (ncol(sp_faxes_coord) > 5) {
      stop(paste("Computing beta functional diversity in a",
                 ncol(sp_faxes_coord), "-dimensions space could exceed
                   computing power of your computer.
                   Consider keeping only five dimensions"))
    }
    
    if(any(! beta.family %in% c("Jaccard", "Sorensen"))) {
      stop(paste("Error: beta diversity index should be 'Jaccard' and/or
                   'Sorensen'. Please check."))
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
  
  
  # ensuring species are in the same order in both matrices
  sp_faxes_coord <- sp_faxes_coord[colnames(asb_sp_occ), ]
  
  # calling functional.betapart.core function to compute beta functional diversity
  # = computing convex hulls and their intersections for all pairs of assemblages
  # 2 last parameters are for parallelization
  F_betapart_core <- betapart::functional.betapart.core(x = asb_sp_occ,
                                                        traits = sp_faxes_coord,
                                                        multi = FALSE,
                                                        return.details = TRUE,
                                                        fbc.step = betapart.step,
                                                        parallel = betapart.para,
                                                        opt.parallel = betapart.para.opt)
  
  
  # computing functional beta diversity indices for all pairs of assemblages
  # according to the type of index (indices) selected
  if("Jaccard" %in% beta.family) {
    F_beta_jac <- betapart::functional.beta.pair(F_betapart_core, index.family = "jaccard")
    
    # indices values in a dataframe where rows are pairs of assemblages
    F_beta_jac_df <- dendextend::dist_long(F_beta_jac$funct.beta.jac)
    names(F_beta_jac_df) <- c("asb.2", "asb.1", "jac_diss")
    
    F_beta_jac_df <- data.frame(F_beta_jac_df,
                                jac_turn = dendextend::dist_long(F_beta_jac$funct.beta.jtu)$distance ,
                                jac_nest = dendextend::dist_long(F_beta_jac$funct.beta.jne)$distance)
    
    # dataframe after reordering columns and proper row names
    pairasb_fbd_indices <- F_beta_jac_df[, c("asb.1", "asb.2",
                                             "jac_diss", "jac_turn", "jac_nest")]
    row.names(pairasb_fbd_indices) <- NULL
    
  }
  
  if( "Sorensen" %in% beta.family) {
    F_beta_sor <- betapart::functional.beta.pair(F_betapart_core, index.family = "sorensen")
    
    # indices values in a dataframe where rows are pairs of assemblages
    F_beta_sor_df <- dendextend::dist_long(F_beta_sor$funct.beta.sor)
    names(F_beta_sor_df) <- c("asb.2", "asb.1", "sor_diss")
    
    F_beta_sor_df <- data.frame(F_beta_sor_df,
                                sor_turn = dendextend::dist_long(F_beta_sor$funct.beta.sim)$distance,
                                sor_nest = dendextend::dist_long(F_beta_sor$funct.beta.sne)$distance)
    
    # dataframe after reordering columns and proper row names
    pairasb_fbd_indices <- F_beta_sor_df[, c("asb.1", "asb.2",
                                             "sor_diss", "sor_turn", "sor_nest")]
    row.names(pairasb_fbd_indices) <- NULL
    
  }
  
  # if both families of indices should be returned
  if(("Sorensen" %in% beta.family) & ("Jaccard" %in% beta.family)) {
    pairasb_fbd_indices <- cbind.data.frame(
      F_beta_jac_df[, c("asb.1", "asb.2","jac_diss", "jac_turn", "jac_nest")],
      F_beta_sor_df[, c( "sor_diss", "sor_turn", "sor_nest")])
    row.names(pairasb_fbd_indices) <- NULL
    
  } #end of if both indices
  
  # compute the volume of the convex hull of the gp to compute real FRic value:
  conv_fa_all <- tryCatch(geometry::convhulln(sp_faxes_coord, option = "FA"),
                          error = function(err){"NA"})
  # if convex hull of the gp can be computed:
  if(is.character(conv_fa_all)) {
    fric <- F_betapart_core$details$CH$FRi/conv_fa_all
  }
  # if convex hull of the gp can not be computed:
  if(! is.character(conv_fa_all)) {
    fric <- NA
  }
  
  
  ## results to return
  if(store_details == TRUE) {
    return_list <- list(pairasb_fbd_indices = pairasb_fbd_indices,
                        details = list(inputs=list(sp_faxes_coord = sp_faxes_coord,
                                                   asb_sp_occ = asb_sp_occ),
                                       asb_FRic_raw = F_betapart_core$details$CH$FRi,
                                       asb_FRic = fric,
                                       asb_vertices = lapply(F_betapart_core$details$CH$coord_vertices, row.names)))
  } else {
    return_list <- (pairasb_fbd_indices = pairasb_fbd_indices )
  }
  
  return(return_list)
  
} # end function
