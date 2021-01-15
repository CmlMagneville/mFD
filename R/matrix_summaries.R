# Functions to compute summaries on species*traits dataframe and
# assemblage*species dataframe
#
# Authors: Camille Magneville & Sébastien Villéger
#
#

# ------------------------------------------------------------------------------


#' Summarize species*traits dataframe
#'
#' Computes summary data helping to choose the type of analysis you can do with
#' your data. For this function to work, there must be no NA in your
#' \code{sp_tr} dataframe.
#'
#' @param tr_cat a \bold{dataframe} containing three columns for each trait (rows):
#'   \itemize{ \item \strong{trait_name}: names of all traits as in \code{sp_tr}
#'   data.frame \item \strong{trait_type}: category codes for each trait as
#'   followed: \emph{N} for Nominal traits (factor variable), \emph{O} for
#'   Ordinal traits (ordered variable), \emph{C} for Circular traits (integer
#'   values), \emph{Q} for quantitative traits (numeric values) and
#'   \emph{F} for fuzzy traits (described with several values defined with
#'   several column). } An option is to add a fourth column with a numeric
#'   vector of length n (traits number) to specify a weight for each trait.
#'
#' @param sp_tr a \bold{dataframe} of traits values (columns) for each species (rows)
#'
#' @return a table summarizing for each trait the number of species per modality
#'   for non-continuous trait and min/max/mean/median/quartiles for continuous
#'   traits \code{tr_summary_list} ; a list containing traits type
#'   \code{tr_types} ; a list containing modalities for each non continuous
#'   trait \code{non_conttr_mod_list}
#'
#' @examples
#'  load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#'  load(system.file("extdata", "sp_tr_cat_fruits_df", package = "mFD"))
#'  mFD::sp.tr.summary(sp_tr_cat, sp_tr)
#'
#' @export


sp.tr.summary <- function(tr_cat, sp_tr) {
  
  # retrieve traits information (fuzzy coded traits are modified thereafter):
  # get traits names:
  tr_nm <- names(sp_tr)
  # get traits number:
  tr_nb <- length(tr_nm)
  # get traits type:
  tr_type <- c(tr_cat$trait_type)
  names(tr_type) <- tr_nm
  
  # retrieve species number:
  sp_nb <- nrow(sp_tr)
  
  
  # check inputs:
  
  if (! is.data.frame(sp_tr))  {
    stop("Error: your species-traits data must be gathered in a matrix")
  }
  if (is.null(colnames(sp_tr))) {
    stop(
      "Error: No column names provided in traits table. Analysis will not go
      through, please add traits as column names"
    )
  }
  if (is.null(rownames(sp_tr))) {
    stop(
      "Error: No row names provided in traits table. Analysis will not go
      through, please add species names as row names"
    )
  }
  if (any(is.na(sp_tr))) {
    stop("Error: There must be no NA in traits table.")
  }
  
  if (any(is.na(tr_cat$trait_type))) {
    stop("Error: Trait type in traits*category data.frame contains NA.
             Please check and specify type of all traits.")
  }
  if (any(tr_nm != tr_cat$trait_name)) {
    stop("Error: Trait names differ between species*traits data.frame and
             traits*category data.frame. Please check.")
  }
  if (any(is.na(tr_cat$trait_type))) {
    stop("Error: Trait type in traits*category data.frame contains NA.
             Please check and specify type of all traits.")
  }
  
  # check nominal traits:
  if ("N" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "N")]) {
      if (is.factor(sp_tr[, k]) == FALSE) {
        stop(paste0("Error: Trait '", k, "'is supposed to be nominal
            but is not described with a 'factor' variable"))
      }
    }
  }
  # check that ordinal traits have ordered categories:
  if ("O" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "O")]) {
      if (is.ordered(sp_tr[, k]) == FALSE) {
        stop(paste0("Error: Trait '", k, "'is supposed to be ordinal
            but is not described with an 'ordered' variable"))
      }
    }
  }
  # check that circular traits are coded with integer:
  if ("C" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "C")]) {
      if (is.integer(sp_tr[, k]) == FALSE) {
        stop( paste0("Error: Trait '", k, "'is supposed to be circular
          but is not described with an 'integer' variable"))
      }
    }
  }
  # check that continuous trait have non unique values:
  if ("Q" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "Q")]) {
      if (is.numeric(sp_tr[, k]) == FALSE) {
        stop(paste0("Error: Trait '", k, "'is supposed to be continuous
                  but is not described with a 'numeric' variable"))
      }
    }
  }
  # check fuzzy coded traits:
  if ("F" %in% tr_cat$trait_type) {
    # retrieve names of fuzzy-coded traits:
    nm_fuzzy <- unique(stats::na.omit(tr_cat$fuzzy_name))
    # check that fuzzy traits are described with more than one variable:
    for (k in nm_fuzzy) {
      # get the names of variables for fuzzy trait k:
      var_k <- tr_cat$trait_name[which(tr_cat$fuzzy_name == k)]
      # check that there are at least 2 columns per fuzzy coded trait:
      if(length(var_k) < 2) {
        stop(paste0("Error: Fuzzy-coded trait '", k,
                    "' is described with a single variable,
            consider changing its type to 'nominal'"))
      }
      # check that variables are continuous:
      if(any(apply(sp_tr[, var_k ], 2, is.numeric) == FALSE)) {
        stop( paste0("Error: Fuzzy-coded trait '", k,
                     "' is not described with 'numeric' variables"))
      }
    }
    
    
    # update names, type and number of traits:
    tr_nm <- c(tr_cat$trait_name[tr_cat$trait_type != "F"], nm_fuzzy)
    tr_nb <- length(tr_nm)
    tr_type <- c(tr_cat$trait_type[tr_cat$trait_type != "F"],
                 rep("F", length(nm_fuzzy)))
    names(tr_type) <- tr_nm
  }
  
  
  # If no fuzzy traits:
  if (! "F" %in% tr_type) {
    # Table with traits summary:
    tr_summary_list <- summary(sp_tr)
    # Vector containing traits types:
    tr_types <- sapply(sp_tr, class)
    # List containing modalities for non continuous traits:
    sp_non_conttr <- sp_tr[, lapply(sp_tr, is.numeric) == FALSE]
    non_conttr_modalities_list <- lapply(sp_non_conttr, unique)
    
    return_list <- list(tr_summary_list = tr_summary_list, tr_types = tr_types,
                        non_conttr_mod_list = non_conttr_modalities_list)
    return(return_list)
    
  # If fuzzy coded traits:
  } else {
    # Table with traits summary for non fuzzy traits:
    tr_summary_list <- summary(sp_tr[, tr_cat$trait_name[
      which(tr_cat$trait_type != "F")], ])
    # Vector containing traits type for non fuzzy traits:
    tr_types <- sapply(sp_tr[, tr_cat$trait_name[
      which(tr_cat$trait_type != "F")], ], class)
    # List containing modalities for non continuous traits:
    sp_non_conttr <- sp_tr[, tr_cat$trait_name[
      which(tr_cat$trait_type != "F" & tr_cat$trait_type != "Q")], ]
    mod_list <- lapply(sp_non_conttr, unique)
    # for fuzzy traits:
    tr_summary_fuzzy_list <- summary(sp_tr[, tr_cat$trait_name[
      which(tr_cat$trait_type == "F")], ])
    
    return_list <- list(tr_summary_non_fuzzy_list = tr_summary_list,
                        tr_summary_fuzzy_list = tr_summary_fuzzy_list,
                        tr_types = tr_types,
                        mod_list = mod_list)
    return(return_list)
  }
  
}


# ------------------------------------------------------------------------------


#' Summarize assemblage*species dataframe
#'
#' Computes summary helping you to picture assemblages. For this function to
#' work, there must be no NA in your \code{assemblage_species} dataframe.
#'
#' @param asb_sp_w a \bold{matrix} showing assemblages (rows) composition in
#'   species (columns).
#'
#' @return a species occurrences data.frame in each assemblage
#'   \code{asb_sp_w_occ}. Use the Data Viewer to see it ; a vector gathering
#'   species biomass/abundance per species \code{tot_ab_all_sp} ; a vector
#'   gathering total abundance/biomass per assemblage \code{tot_ab_all_asb} ; a
#'   vector gathering species richness per assemblage \code{sp_richn_all_asb} ;
#'   a list gathering the names of species of each assemblage \code{sp_nm_asb}.
#'
#' @examples
#' load(system.file("extdata", "asb_sp_w_fruits", package = "mFD"))
#' asb_sp_w <- as.matrix(asb_sp_w)
#' mFD::asb.sp.summary(asb_sp_w)
#'
#' @export

asb.sp.summary <- function(asb_sp_w) {
  
  if (is.matrix(asb_sp_w) == FALSE) {
    stop("Error: 'asb_sp_w' must be a matrix")
  }
  
  if (is.null(colnames(asb_sp_w))) {
    stop(
      "Error: No column names provided in asb_sp_w matrix.
      Analysis will not go through,
      please add assemblages names as column names"
    )
  }
  if (is.null(rownames(asb_sp_w))) {
    stop(
      "Error: No row names provided in asb_sp_w matrix.
      Analysis will not go through,
      please add species names as row names"
    )
  }
  if (is.null(asb_sp_w)) {
    stop("Error: NA in asb_sp_w matrix. Analysis will not go through.")
  }
  
  isnum_vect <- sapply(asb_sp_w, is.numeric)
  
  if (FALSE %in% isnum_vect) {
    stop("Error: The 'asp_sp_w' matrix must only contain numeric values. Please convert values")
  }
  
  # Add a stop if there is a negative value in the occurrence dataframe:
  if (any(asb_sp_w < 0)) {
    stop("Error: The species*weight matrix should not contain negative values.
           Please check.")
  }
  
  # convert asb_sp_w to df so can use dplyr fcts:
  asb_sp_w <- as.data.frame(asb_sp_w)
  
  # Species occurrence dataframe:
  asb_sp_w_occ <- replace(asb_sp_w, asb_sp_w > 1, 1)
  
  # Vector containing the number of occurrences of each species:
  nbocc_sp <- apply(asb_sp_w, 2, sum)
  
  # Vector containing total abundance/biomass per assemblage:
  asb_totab <- apply(asb_sp_w, 1, sum)
  
  # Vector containing species richness of each assemblage:
  asb_sp_wrichn <- apply(asb_sp_w_occ, 1, sum)
  
  # Construction of a list containing vectors with species names for each...
  # ... assemblage:
  L <- list()
  for (i in (1:nrow(asb_sp_w_occ))) {
    data <- dplyr::slice(asb_sp_w_occ, i)
    data2 <- data[which(apply(data, 2, max) == TRUE)]
    rownames(data2) <- rownames(asb_sp_w_occ[i, ])
    data3 <- as.vector(data2)
    L[[i]] <-  data3
  }
  
  # Create the return object:
  return_list <- list(asb_sp_occ = asb_sp_w_occ,
                      tot_ab_all_sp = nbocc_sp,
                      tot_ab_all_asb = asb_totab,
                      sp_richn_all_asb = asb_sp_wrichn,
                      sp_nm_asb = L)
  
  # Add a warning if some species do not belong to any assemblage:
  if (min(apply(asb_sp_w, 2, sum)) == 0){
    warning("Some species are absent from all assemblages")
  }
  # Add a warning if some asb do not contain species:
  if (min(apply(asb_sp_w, 1, sum)) == 0){
    warning("Some assemblages do not contain species.")
  }
  
  return(return_list)
}
