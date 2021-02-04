# Helpers to Check Functions Inputs
# 
#' @importFrom stats na.omit
#' @keywords internal

check.tr.cat <- function(tr_cat) {
  
  if (!is.data.frame(tr_cat)) {
    stop("Traits x category data must be gathered in a data frame.")
  }
  
  if (ncol(tr_cat) < 3) {
    stop("Traits x category data frame must have 3 (or 4) columns.")
  }
  
  valid_names  <- c("trait_name", "trait_type", "fuzzy_name")
  valid_traits <- c("N", "O", "C", "Q", "F")
  
  if (!identical(colnames(tr_cat)[1:3], valid_names)) {
    stop("The 3 first columns of the traits x category data frame must be ",
         "'trait_name', 'trait_type', and 'fuzzy_name' in this exact order.")
  }
  
  if (any(is.na(tr_cat$"trait_type"))) {
    stop("Trait type in traits x category data frame contains NA. Please ",
         "check and specify type of all traits.")
  }
  
  if (any(!(tr_cat$"trait_type" %in% valid_traits))) {
    stop("Trait type in traits x category should be among 'N', 'O', 'C', 'Q', ",
         "'F'. Please check type of all traits.")
  }
}



check.sp.tr <- function(sp_tr, tr_cat = NULL, stop_if_NA = TRUE) {
  
  if (!is.logical(stop_if_NA)) {
    stop("Argument 'stop_if_NA' must be a boolean.")
  }
  
  if (!is.data.frame(sp_tr)) {
    stop("Species x traits data must be gathered in a data frame.")
  }
  
  if (any(sort(rownames(sp_tr)) == 1:nrow(sp_tr))) {
    stop("No row names provided in the species x traits data frame. Analysis ",
         "will not go through, please add species names as row names.")
  }
  
  if (stop_if_NA) {
    if (any(is.na(sp_tr))) {
      stop("There must be no NA in the species x traits data frame.")
    } 
  }
  
  
  if (!is.null(tr_cat)) {
    
    check.tr.cat(tr_cat)
    
    if (length(names(sp_tr)) != length(tr_cat$"trait_name")) {
      stop("Trait numbers differ between species x traits data frame and ",
           "traits x category data frame. Please check.")
    }
    
    if (any(sort(names(sp_tr)) != sort(tr_cat$"trait_name"))) {
      stop("Trait names differ between species x traits data frame and ",
           "traits x category data frame. Please check.")
    }
  }
}



check.nominal <- function(tr_cat, sp_tr) {
  
  if ("N" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "N")]) {
      if (!is.factor(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be nominal but is not described ", 
             "with a 'factor' variable.")
      }
    }
  }
}



check.ordinal <- function(tr_cat, sp_tr) {
  
  if ("O" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "O")]) {
      if (!is.ordered(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be ordinal but is not described ",
             "with an 'ordered' variable.")
      }
    }
  }
}



check.circular <- function(tr_cat, sp_tr) {
  
  if ("C" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "C")]) {
      if (!is.integer(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be circular but is not described ",
             "with an 'integer' variable.")
      }
    }
  }
}



check.continuous <- function(tr_cat, sp_tr) {
  
  if ("Q" %in% tr_cat$"trait_type") {
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "Q")]) {
      if (!is.numeric(sp_tr[ , k])) {
        stop("Trait '", k, "' is supposed to be continuous but is not ",
             "described with a 'numeric' variable.")
      }
    }
  }
}



check.fuzzy <- function(tr_cat, sp_tr) {
  
  if ("F" %in% tr_cat$"trait_type") {
  
    nm_fuzzy <- unique(stats::na.omit(tr_cat$"fuzzy_name"))

    for (k in nm_fuzzy) {
      
      var_k <- tr_cat$"trait_name"[which(tr_cat$"fuzzy_name" == k)]
      
      if (length(var_k) < 2) {
        stop("Fuzzy-coded trait '", k, "' is described with a single ",
             "variable. Consider changing its type to 'nominal'.")
      }
      
      if (!any(apply(sp_tr[ , var_k], 2, is.numeric))) {
        stop("Fuzzy-coded trait '", k, "' is not described with 'numeric' ", 
             "variables.")
      }
    }
  }
}



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
    stop("The species x weight matrix should not contain negative values.",
         "Please check.")
  }
  
  
  if (min(apply(asb_sp_w, 2, sum)) == 0) {
    warning("Some species are absent from all assemblages.")
  }
  
  if (min(apply(asb_sp_w, 1, sum)) == 0) {
    warning("Some assemblages do not contain species.")
  }
}



check.asb.sp.tr <- function() { }

check.sp.coord <- function() { }

