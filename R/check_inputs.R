##
## HELPERS FUNCTIONS TO CHECK {mFD} FUNCTIONS INPUTS
##


#' Traits x Categories Data Frame
#' @noRd

check.tr.cat <- function(tr_cat) {
  
  if (!is.data.frame(tr_cat)) {
    stop("Traits x category data must be gathered in a data frame.")
  }
  
  if (ncol(tr_cat) < 2) {
    stop("Traits x category data frame must have at least 2 columns.")
  }
  
  if (ncol(tr_cat) == 2) {
    valid_names  <- c("trait_name", "trait_type")
    if (!identical(colnames(tr_cat), valid_names)) {
      stop("The 2 columns of the traits x category data frame must be ",
           "'trait_name', and 'trait_type' in this exact order.")
    }
  }
  
  # Check for a third column if some traits are fuzzy:
  if (ncol(tr_cat) == 3 & any(tr_cat$"trait_type" == "F")) {
    valid_names  <- c("trait_name", "trait_type", "fuzzy_name")
    if (!identical(colnames(tr_cat), valid_names)) {
      stop("The 3 first columns of the traits x category data frame must be ",
           "'trait_name', 'trait_type', and 'fuzzy_name' in this exact order.
           Weight can not be provided if fuzzy traits are used with mFD,
           please use gawdis package instead")
    }
  } 
  
  # Check for a fourth column if trait weight is needed:
  if (ncol(tr_cat) == 3 & (is.element(c("F"), tr_cat$trait_type) == FALSE)) {
    valid_names  <- c("trait_name", "trait_type","trait_weight")
    if (!identical(colnames(tr_cat), valid_names)) {
      stop("The 3 first columns of the traits x category data frame must be ",
           "'trait_name', 'trait_type', and 'trait_weight' in this exact ", 
           "order.")
    }
  }
  
  
  if (any(is.na(tr_cat$"trait_type"))) {
    stop("Trait type in traits x category data frame contains NA. Please ",
         "check and specify type of all traits.")
  }
  
  valid_traits <- c("N", "O", "C", "Q", "F")
  
  if (any(!(tr_cat$"trait_type" %in% valid_traits))) {
    stop("Trait type in traits x category should be among 'N', 'O', 'C', ",
         "'Q', 'F'. Please check type of all traits.")
  }
  
  if (any(tr_cat$"trait_type" == "F")) {
    
    if (any(is.na(tr_cat[which(tr_cat$"trait_type" == "F"), "fuzzy_name"]))) {
      stop("Missing trait names in 'fuzzy_name' for fuzzy traits.")
    }
    
    if (any(table(tr_cat[which(tr_cat$"trait_type" == "F"), "fuzzy_name"]) < 
            2)) {
      stop("Fuzzy traits need to have at least two categories.")
    }
  }
  
  invisible(NULL)
}



#' Species x Traits Data Frame
#' @noRd

check.sp.tr <- function(sp_tr, tr_cat = NULL, stop_if_NA = TRUE) {
  
  if (!is.logical(stop_if_NA)) {
    stop("Argument 'stop_if_NA' must be a boolean.")
  }
  
  if (!is.data.frame(sp_tr)) {
    stop("Species x traits data must be gathered in a data frame.")
  }
  
  if (any(sort(rownames(sp_tr)) == 1:nrow(sp_tr))) {
    stop("No row names provided in the species x traits data frame. Analysis ",
         "will not go through. Please add species names as row names.")
  }
  
  if (stop_if_NA) {
    if (any(is.na(sp_tr))) {
      stop("Species x traits data frame contains NA. If you want to ",
           "continue with missing traits (Be careful: Functional measures ", 
           "are sensitive to missing traits), set 'stop_if_NA' parameter ",
           "to FALSE. Otherwise you can delete species with missing or ",
           "extrapolate missing traits (Johnson et al. (2020).")
    } 
  }
  
  
  if (!is.null(tr_cat)) {
    
    check.tr.cat(tr_cat)
    
    if (length(names(sp_tr)) != length(tr_cat$"trait_name")) {
      stop("Trait numbers differ between species x traits data frame and ",
           "traits x category data frame. Please check.")
    }
    
    if (any(names(sp_tr) != tr_cat$"trait_name")) {
      stop("Trait names differ between species x traits data frame and ",
           "traits x category data frame or are not in the ame order. 
           Please check.")
    }
  }
  
  invisible(NULL)
}



#' Assemblages x Species Matrix
#' @noRd

check.asb.sp.w <- function(asb_sp_w) { 
  
  if (!is.matrix(asb_sp_w)) {
    stop("The argument 'asb_sp_w' must be a matrix.")
  }
  
  if (!(is.numeric(asb_sp_w))) {
    stop("The 'asp_sp_w' matrix must only contain numeric values. Please ", 
         "convert values.")
  }
  
  if (any(sort(rownames(asb_sp_w)) == 1:nrow(asb_sp_w))) {
    stop("No row names provided in 'asb_sp_w' matrix. Analysis will not go ", 
         "through. Please add assemblages names as row names.")
  }
  
  if (any(is.na(asb_sp_w))) {
    stop("The 'asb_sp_w' matrix contains NA. Analysis will not go through.")
  }
  
  
  if (any(asb_sp_w < 0)) {
    stop("The species x weight matrix should not contain negative values. ",
         "Please check.")
  }
  
  if (min(apply(asb_sp_w, 2, sum)) == 0) {
    warning("Some species are absent from all assemblages.")
  }
  
  if (min(apply(asb_sp_w, 1, sum)) == 0) {
    warning("Some assemblages do not contain species.")
  }
  
  invisible(NULL)
}



#' Species Coordinates on Axes
#' @noRd

check.sp.faxes.coord <- function(sp_faxes_coord) {
  
  if (!is.matrix(sp_faxes_coord)) {
    stop("The species x coordinates data must be a matrix.")  
  }
  
  if (any(is.na(sp_faxes_coord))) {
    stop("The species x coordinates matrix contains NA. Please check.")
  }
  
  if (any(sort(rownames(sp_faxes_coord)) == 1:nrow(sp_faxes_coord))) {
    stop("No row names provided in the species x coordinates matrix. Please ", 
         "add species names as row names.")
  }
  
  if (is.null(colnames(sp_faxes_coord))) {
    stop("No column names provided in species*coordinates matrix. Please add ",
         "axes names as column names.")
  }
  
  invisible(NULL)
}



#' Nominal Traits
#' @noRd

check.nominal <- function(tr_cat, sp_tr) {
  
  if ("N" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "N")]) {
      if (!is.factor(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be nominal but is not described ", 
             "with a 'factor' variable.")
      }
    }
  }
  
  invisible(NULL)
}



#' Ordinal Traits
#' @noRd

check.ordinal <- function(tr_cat, sp_tr) {
  
  if ("O" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "O")]) {
      if (!is.ordered(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be ordinal but is not described ",
             "with an 'ordered' variable.")
      }
    }
  }
  
  invisible(NULL)
}



#' Circular Traits
#' @noRd

check.circular <- function(tr_cat, sp_tr) {
  
  if ("C" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "C")]) {
      if (!is.integer(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be circular but is not ",
             "described with an 'integer' variable.")
      }
    }
  }
  
  invisible(NULL)
}



#' Continuous Traits
#' @noRd

check.continuous <- function(tr_cat, sp_tr) {
  
  if ("Q" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "Q")]) {
      if (!is.numeric(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be continuous but is not ",
             "described with a 'numeric' variable.")
      }
    }
  }
  
  invisible(NULL)
}



#' Fuzzy Traits
#' @noRd

check.fuzzy <- function(tr_cat, sp_tr) {
  
  if ("F" %in% tr_cat$"trait_type") {
  
    nm_fuzzy <- unique(stats::na.omit(tr_cat$"fuzzy_name"))

    for (k in nm_fuzzy) {
      
      var_k <- tr_cat$"trait_name"[which(tr_cat$"fuzzy_name" == k)]
      
      if (length(as.character(var_k)) < 2) {
        stop("Fuzzy-coded trait '", k, "' is described with a single ",
             "variable. Consider changing its type to 'nominal'.")
      }
      
      if (!any(apply(sp_tr[ , as.character(var_k)], 2, is.numeric))) {
        stop("Fuzzy-coded trait '", k, "' is not described with 'numeric' ", 
             "variables.")
      }
    }
  }
  
  invisible(NULL)
}


#' ...
#' @noRd

check.asb.sp.w.occ <- function(asb_sp_occ) {
  
  if (is.null(rownames(asb_sp_occ))) {
    stop("No row names provided in the species occurence matrix. Please add ",
         "assemblages names as row names.")
  }
  
  if (is.null(colnames(asb_sp_occ))) {
    stop("No column names provided in the species occurence matrix. Please ",
         "add species names as column names.")
  }
  
  if (any(is.na(asb_sp_occ))) {
    stop("The species occurence matrix contains NA. Please check.")
  }
  
  if (any(asb_sp_occ != 0 & asb_sp_occ != 1)) {
    stop("The species occurence matrix should contain only 0 and 1. Please ",
         "check.")
  }
  
  # Add a stop if some species do not belong to any
  # assemblage:
  if (min(apply(asb_sp_occ, 2, sum)) == 0) {
    stop("Some species are absent from all assemblages.")
  }
  
  # Add a stop if some asb do not contain species:
  if (min(apply(asb_sp_occ, 1, sum)) == 0) {
    stop("Some assemblages do not contain species.")
  }
  
  invisible(NULL)
}
