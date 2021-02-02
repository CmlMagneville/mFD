#' Summarize Species x Traits Data Frame
#'
#' This function computes a summary data helping to choose the type of analysis 
#' you can do with your data. For this function to work, there must be no NA in 
#' your `sp_tr` data frame.
#'
#' @param tr_cat a data frame containing three columns for each trait (rows):
#'   \itemize{
#'   \item \strong{trait_name}: the name of all traits as in \code{sp_tr} data 
#'   frame;
#'   \item \strong{trait_type}: the category code for each trait as followed:
#'   \code{N} for Nominal traits (factor variable), \code{O} for Ordinal traits
#'   (ordered variable), \code{C} for Circular traits (integer values), \code{Q}
#'   for quantitative traits (numeric values), and \code{F} for 
#'   fuzzy traits (i.e. described with several values defined with several 
#'   column);
#'   \item \strong{fuzzy_name}: name of fuzzy-coded trait to which 'sub-trait'
#'   belongs (if trait is not fuzzy, ignored so could be trait name or NA).
#'   }
#'   
#'   An option is to add a fourth column with a numeric vector of length n 
#'   (traits number) to specify a weight for each trait.
#'
#' @param sp_tr a data frame of traits values (columns) for each species (rows).
#'   Note that species names **must be** specified in the row names.
#'
#' @return 
#'   If there is no fuzzy-coded trait, a three-elements list with:
#'   \item{tr_summary_list}{a table summarizing for each trait the 
#'   number of species per modality for non-continuous trait and min, max, 
#'   mean, median, and quartiles for continuous traits.}
#'   \item{tr_types}{a list containing traits type.}
#'   \item{non_conttr_mod_list}{a list containing modalities for non-continuous 
#'   trait.}
#'   
#'   If there is fuzzy-coded trait, a four-elements list with:
#'   \item{tr_summary_non_fuzzy_list}{a table summarizing for each trait the 
#'   number of species per modality for non-continuous trait and min, max, 
#'   mean, median, and quartiles for continuous traits.}
#'   \item{tr_summary_fuzzy_list}{__MISSING DESCRIPTION__}
#'   \item{tr_types}{a list containing traits type.}
#'   \item{mod_list}{a list containing modalities for non-continuous trait.}
#'
#' @author Camille Magneville & Sebastien Villeger
#'
#' @export
#' @importFrom stats na.omit
#'
#' @examples
#' load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#' load(system.file("extdata", "sp_tr_cat_fruits_df", package = "mFD"))
#' mFD::sp.tr.summary(sp_tr_cat, sp_tr)


sp.tr.summary <- function(tr_cat, sp_tr) {
  
  
  ## Retrieve Traits Informations ----
  ## (fuzzy-coded traits are modified thereafter)
  
  sp_nb   <- nrow(sp_tr)            # species number
  tr_nm   <- names(sp_tr)           # traits names
  tr_nb   <- length(tr_nm)          # traits number
  tr_type <- c(tr_cat$trait_type)   # traits type
  
  names(tr_type) <- tr_nm

   
  ## Check Inputs ----
  
  if (!is.data.frame(sp_tr)) {
    stop("Your species x traits data must be gathered in a data frame.")
  }
  
  if (any(rownames(sp_tr) == 1:nrow(sp_tr))) {
    stop(paste("No row names provided in traits data frame. Analysis will",
               "not go through, please add species names as row names."))
  }
  
  if (any(is.na(sp_tr))) {
    stop("There must be no NA in traits data.")
  }
  
  if (any(is.na(tr_cat$trait_type))) {
    stop(paste("Trait type in traits x category data frame contains NA. Please",
               "check and specify type of all traits."))
  }
  
  if (any(sort(tr_nm) != sort(tr_cat$trait_name))) {
    stop(paste("Trait names differ between species x traits data frame and",
               "traits x category data frame. Please check."))
  }
 
   
  ## Check for Nominal Traits ----
  
  if ("N" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "N")]) {
      if (!is.factor(sp_tr[ , k])) {
        stop(paste0("Trait '", k, "'is supposed to be nominal but is not ", 
                    "described with a 'factor' variable."))
      }
    }
  }
  
  
  ## Check that Ordinal Traits have ordered categories ----
  
  if ("O" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "O")]) {
      if (!is.ordered(sp_tr[ , k])) {
        stop(paste0("Trait '", k, "'is supposed to be ordinal but is not ", 
                    "described with an 'ordered' variable."))
      }
    }
  }
  
  
  ## Check that Circular Traits are coded with integer ----
  
  if ("C" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "C")]) {
      if (!is.integer(sp_tr[ , k])) {
        stop(paste0("Trait '", k, "'is supposed to be circular but is not ", 
                    "described with an 'integer' variable."))
      }
    }
  }
  
  
  ## Check that Continuous Trait have NON-Unique Values ----
  
  if ("Q" %in% tr_cat$trait_type) {
    for (k in tr_cat$trait_name[which(tr_cat$trait_type == "Q")]) {
      if (!is.numeric(sp_tr[ , k])) {
        stop(paste0("Trait '", k, "'is supposed to be continuous but is not ", 
                    "described with a 'numeric' variable."))
      }
    }
  }
  
  
  ## Check Fuzzy-coded Traits ----
  
  if ("F" %in% tr_cat$trait_type) {
    
    # retrieve names of fuzzy-coded traits:
    nm_fuzzy <- unique(stats::na.omit(tr_cat$fuzzy_name))
    
    # Check that fuzzy traits are described with more than one variable:
    for (k in nm_fuzzy) {
      # get the names of variables for fuzzy trait k:
      var_k <- tr_cat$trait_name[which(tr_cat$fuzzy_name == k)]
      
      # Check that there are at least 2 columns per fuzzy coded trait:
      if (length(var_k) < 2) {
        stop(paste0("Fuzzy-coded trait '", k, "' is described with a single ", 
                    "variable. Consider changing its type to 'nominal'."))
      }
      
      # Check that variables are continuous:
      if (!any(apply(sp_tr[ , var_k ], 2, is.numeric))) {
        stop(paste0("Fuzzy-coded trait '", k, "' is not described with ",
                    "'numeric' variables"))
      }
    }
    
    # Update names, type and number of traits
    tr_nm   <- c(tr_cat$trait_name[tr_cat$trait_type != "F"], nm_fuzzy)
    tr_nb   <- length(tr_nm)
    tr_type <- c(tr_cat$trait_type[tr_cat$trait_type != "F"],
                 rep("F", length(nm_fuzzy)))
    names(tr_type) <- tr_nm
  }
  
  
  ## Function Return  ----
  
  if (!("F" %in% tr_type)) {          # no fuzzy traits
    
    # Table with traits summary
    tr_summary_list <- summary(sp_tr)
    
    # Vector containing traits types
    tr_types <- sapply(sp_tr, class)
    
    # List containing modalities for non continuous traits
    sp_non_conttr <- sp_tr[, lapply(sp_tr, is.numeric) == FALSE]
    non_conttr_modalities_list <- lapply(sp_non_conttr, unique)
    
    return(list("tr_summary_list"     = tr_summary_list, 
                "tr_types"            = tr_types,
                "non_conttr_mod_list" = non_conttr_modalities_list))
    
  } else {                            # fuzzy coded traits

    # Table with traits summary for non fuzzy traits
    tr_summary_list <- summary(sp_tr[ , tr_cat$trait_name[
      which(tr_cat$trait_type != "F")], ])
    
    # Vector containing traits type for non fuzzy traits
    tr_types <- sapply(sp_tr[ , tr_cat$trait_name[
      which(tr_cat$trait_type != "F")], ], class)
    
    # List containing modalities for non continuous traits
    sp_non_conttr <- sp_tr[ , tr_cat$trait_name[
      which(tr_cat$trait_type != "F" & tr_cat$trait_type != "Q")], ]
    mod_list <- lapply(sp_non_conttr, unique)
    
    # For fuzzy traits
    tr_summary_fuzzy_list <- summary(sp_tr[ , tr_cat$trait_name[
      which(tr_cat$trait_type == "F")], ])
    
    return(list("tr_summary_non_fuzzy_list" = tr_summary_list,
                "tr_summary_fuzzy_list"     = tr_summary_fuzzy_list,
                "tr_types"                  = tr_types,
                "mod_list"                  = mod_list))
  }
}



#' Summarize Assemblage x Species Data Frame
#'
#' This function computes a summary helping you to picture assemblages. For this 
#' function to work, there must be no NA in your `asb_sp_w` data frame.
#'
#' @param asb_sp_w a matrix showing assemblages (rows) composition in species 
#'   (columns). Note that species names **must be** the names of rows.
#'
#' @return A list with:
#'   \item{asb_sp_w_occ}{a data frame with species occurrences in each 
#'   assemblage.}
#'   \item{tot_ab_all_sp}{a vector gathering species biomass/abundance per 
#'   species.}
#'   \item{tot_ab_all_asb}{a vector gathering total abundance/biomass per 
#'   assemblage.}
#'   \item{sp_richn_all_asb}{a vector gathering species richness per 
#'   assemblage.}
#'   \item{sp_nm_asb}{a list gathering the names of species of each assemblage.}
#'   
#' @author Camille Magneville & Sebastien Villeger
#'
#' @export
#' 
#' @examples
#' load(system.file("extdata", "asb_sp_w_fruits", package = "mFD"))
#' asb_sp_w <- as.matrix(asb_sp_w)
#' mFD::asb.sp.summary(asb_sp_w)

asb.sp.summary <- function(asb_sp_w) {
  
  
  ## Check Input ----
  
  if (!is.matrix(asb_sp_w)) {
    stop("The argument 'asb_sp_w' must be a matrix.")
  }
  
  if (!(is.numeric(asb_sp_w))) {
    stop(paste("The 'asp_sp_w' matrix must only contain numeric values. Please",
               "convert values."))
  }
  
  if (any(rownames(asb_sp_w) == 1:nrow(asb_sp_w))) {
    stop(paste("No row names provided in 'asb_sp_w' matrix. Analysis will not", 
               "go through. Please add species names as row names."))
  }
  
  if (any(is.na(asb_sp_w))) {
    stop("The 'asb_sp_w' matrix contains NA. Analysis will not go through.")
  }

  
  if (any(asb_sp_w < 0)) {
    stop(paste("The species x weight matrix should not contain negative values.",
               "Please check."))
  }
  
  # Convert matrix to data frame
  asb_sp_w <- as.data.frame(asb_sp_w)
  
  # Species occurrence data frame
  asb_sp_w_occ <- replace(asb_sp_w, asb_sp_w != 0, 1)
  
  # Vector containing the number of occurrences of each species
  nbocc_sp <- apply(asb_sp_w, 2, sum)
  
  # Vector containing total abundance/biomass per assemblage
  asb_totab <- apply(asb_sp_w, 1, sum)
  
  # Vector containing species richness of each assemblage
  asb_sp_wrichn <- apply(asb_sp_w_occ, 1, sum)
  
  
  # Construction of a list containing vectors with species names for each...
  # ... assemblage
  L <- list()
  
  for (i in 1:nrow(asb_sp_w_occ)) {
    
    data  <- asb_sp_w_occ[i, ]
    data2 <- data[which(apply(data, 2, max) == TRUE)]
    rownames(data2) <- rownames(asb_sp_w_occ[i, ])
    L[[i]] <- as.vector(data2)
  }
  
  
  ## Function Warnings ----
  
  # Add a warning if some species do not belong to any assemblage
  if (min(apply(asb_sp_w, 2, sum)) == 0) {
    warning("Some species are absent from all assemblages.")
  }
  
  # Add a warning if some asb do not contain species
  if (min(apply(asb_sp_w, 1, sum)) == 0) {
    warning("Some assemblages do not contain species.")
  }
  
  
  ## Function Return ----
  
  list("asb_sp_occ"       = asb_sp_w_occ,
       "tot_ab_all_sp"    = nbocc_sp,
       "tot_ab_all_asb"   = asb_totab,
       "sp_richn_all_asb" = asb_sp_wrichn,
       "sp_nm_asb"        = L)
}
