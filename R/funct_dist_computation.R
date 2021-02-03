#' Compute Functional Distance between Species
#'
#' For a given combination of traits, this function returns the functional
#' distance matrix between species.
#'
#' @param sp_tr a data frame of traits values (columns) for each species (rows).
#' 
#' @param tr_cat a data frame containing three columns for each trait (rows):
#'   \itemize{
#'   \item \strong{trait_name}: the name of all traits as in \code{sp_tr} data 
#'   frame;
#'   \item \strong{trait_type}: the category code for each trait as followed:
#'   \code{N} for Nominal traits (factor variable), \code{O} for Ordinal traits
#'   (ordered variable), \code{C} for Circular traits (integer values), \code{Q}
#'   for quantitative traits (numeric values) that is allowed \strong{only} if
#'   there are at least 2 species with the same value, and \code{F} for fuzzy 
#'   traits (i.e. described with several values defined with several column);
#'   \item \strong{fuzzy_name}: name of fuzzy-coded trait to which 'sub-trait'
#'   belongs (if trait is not fuzzy, ignored so could be trait name or NA).
#'   }
#'   
#'   An option is to add a fourth column with a numeric vector of length n 
#'   (traits number) to specify a weight for each trait.
#'
#' @param dist_metric the distance to be computed:
#'   `euclidean`, the Euclidean distance, 
#'   `classical_gower`, the Classical Gower distance as defined by Gower (1971),
#'   or `kgower`, the Gower distance modified by Pavoine _et al._ (2009).
#'   
#'   Based on the \code{\link[ade4]{dist.ktab}} function: "_If `kgower` is 
#'   chosen user have several choice to scale quantitative variables A string 
#'   that can have three values: either \code{scaledBYrange} if the quantitative
#'   variables must be scaled by their range, or \code{scaledBYsd} if they must 
#'   be scaled by their standard deviation, or \code{noscale} if they should not 
#'   be scaled. This last option can be useful if the values have already been
#'   normalized by the known range of the whole population instead of the 
#'   observed range measured on the sample. If x contains data from various 
#'   types, then the option \code{scaledBYsd} is not suitable (a warning will 
#'   appear if the option selected with that condition)"_.
#' 
#' @param scaling a character string referring to the way traits must be scaled. 
#'   There are three options: 
#'   `scaledBYrange` (if traits must be scaled byrange),
#'   `scaledBYsd` (if traits must be scaled by their standard deviation), or
#'   `noscale` (if traits do not have to be scaled). 
#' 
#' @param stop_if_NA a logical value to stop or not the process if the
#'   `sp_tr` data frame contains NA. Functional measures are sensitive to
#'   missing traits. For further explanations, see the Note section.
#'   Default is `TRUE`.
#'
#' @return A `dist` object containing distance between each pair of species.
#'
#' @note If the `sp_tr` data frame contains `NA` you can either
#'   chose to compute anyway functional distances (but keep in mind that
#'   **Functional measures are sensitive to missing traits!**) or you can
#'   delete species with missing or extrapolate missing traits (see
#'   Johnson _et al._ (2020)).
#'
#' @references 
#'   Gower, J.C. (1971) A general coefficient of similarity and some of its 
#'   properties. _Biometrics_, **27**, 857-871.\cr
#'   Johnson _et al._ (2020) __\{\{ ADD THE COMPLETE REFERENCE \}\}__\cr
#'   Pavoine _et al._ (2009) __\{\{ ADD THE COMPLETE REFERENCE \}\}__
#' 
#' @author Nicolas Loiseau & Sébastien Villéger
#'
#' @export
#' @importFrom ade4 prep.fuzzy ktab.list.df dist.ktab
#' @importFrom cluster daisy
#' 
#' @examples
#' load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#' sp_tr <- sp_tr[ , -c(6:8)]
#' 
#' load(system.file("extdata", "sp_tr_cat_fruits_df", package = "mFD"))
#' sp_tr_cat <- sp_tr_cat[-c(6:8), ]
#' 
#' mFD::funct.dist(sp_tr, sp_tr_cat, dist_metric = "classical_gower", 
#'                 scaling = "noscale", stop_if_NA = TRUE)

funct.dist <- function(sp_tr, tr_cat, dist_metric, scaling, stop_if_NA = TRUE) {
  
  
  ## Check Inputs ----
  
  if (any(is.na(sp_tr)) && stop_if_NA) {
    stop("Species x traits data frame contains NA. If you want to ",
         "continue with missing traits (Be careful: Functional measures ", 
         "are sensitive to missing traits), set 'stop_if_NA' parameter ",
         "to FALSE. Otherwise you can delete species with missing or ",
         "extrapolate missing traits (Johnson et al. (2020).")
  }
  
  if (!is.data.frame(sp_tr)) {
    stop("Species x traits data must be gathered in a matrix.")
  }
  
  if (any(rownames(sp_tr) == 1:nrow(sp_tr))) {
    stop(paste("No row names provided in traits data frame. Analysis will",
               "not go through, please add species names as row names."))
  }
  
  if (any(is.na(sp_tr))) {
    stop("There must be no NA in traits data frame.")
  }
  
  if (any(is.na(tr_cat$trait_type))) {
    stop("Trait type in traits x category data frame contains NA. Please ",
         "check and specify type of all traits.")
  }
  
  tr_nm <- names(sp_tr)
  
  if (any(tr_nm != tr_cat$trait_name)) {
    stop("Trait names differ between species x traits data frame and ",
         "traits x category data frame. Please check.")
  }
  
  if (ncol(tr_cat) == 4) {
    cat("tr_cat has 4 columns, if you use classical_gower traits will be",
        "weighted.\n")
  } 
  
  if (any(!(dist_metric %in% c("euclidean", "classical_gower", "kgower")))) {
    stop("Argument 'dist_metric' should be 'euclidean', 'classical_gower', ", 
         "or 'kgower'.")
  }
  
  
  ## Compute Distances ----
  
  if (dist_metric == "classical_gower") {
    
    if (ncol(tr_cat) == 4) {
      
      ktab_dist <- cluster::daisy(sp_tr, "gower", weights = tr_cat$weight)
      
    } else {
      
      ktab_dist <- cluster::daisy(sp_tr, "gower")
    }
    
  }
  
  if (dist_metric == "euclidean") {
    
    if (!is.numeric(sp_tr)) {
      stop("Some traits are not numerical.")
    }

    ktab_dist <- stats::dist(sp_tr, method = "euclidean")
  }
    
  
  if (dist_metric == "kgower") {
    
    # Need to prepare functional trait in function of their nature
    
    
    # Quantitative traits
    
    quant_trait <- NULL
    
    if (any(tr_cat$trait_type == "Q")) {
      quant_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "Q"],
                            drop = FALSE]
    }
    
    
    # Ordinal Traits
    
    ord_trait <- NULL
    
    if (any(tr_cat$trait_type == "O")) {
      ord_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "O"],
                          drop = FALSE]
    }
    
    
    # Circular traits
    
    circ_trait <- NULL
    
    if (any(tr_cat$trait_type == "C")) {
      circ_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "C"],
                           drop = FALSE]
      
      circ_trait <- ade4::prep.circular(circ_trait, 1, 12)
    }
    
    
    # Fuzzy traits 
    # (basically several categories that are considered a single trait)
    
    fuzz_trait <- NULL
    
    if (any(tr_cat$trait_type == "F")) {
      
      # Select the fuzzy traits
      fuzz_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "F"],
                           drop = FALSE]
      
      # Count the number of fuzzy categories
      fuzz_cat <- table(tr_cat[tr_cat$trait_type == "F", ]$fuzzy_name)
      
      # Order the trait names based on the order of the categories
      fuzz_names_ordered <- unlist(lapply(names(fuzz_cat), function(x) {
        tr_cat$trait_name[tr_cat$fuzzy_name == x & !is.na(tr_cat$fuzzy_name)]
      }))
      
      # Reorder the traits according to the names
      fuzz_trait <- fuzz_trait[ , fuzz_names_ordered]
      
      fuzz_trait <- ade4::prep.fuzzy(fuzz_trait,
                                     col.blocks = as.numeric(fuzz_cat),
                                     labels     = names(fuzz_cat))
    }
    
    
    # Binary traits
    
    bin_trait <- NULL
    
    if (any(tr_cat$trait_type == "B")) {
      bin_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "B"], 
                          drop = FALSE]
      
      bin_trait <- ade4::prep.binary(bin_trait, col.blocks = ncol(bin_trait))
    }
    
    
    # Nominal traits
    
    nom_trait <- NULL
    
    if (any(tr_cat$trait_type == "N")) {
      nom_trait <- sp_tr[ , tr_cat$trait_name[tr_cat$trait_type == "N"], 
                          drop = FALSE]
    }
    
    # Combine all traits
    all_trait <- list("Q" = quant_trait, "O" = ord_trait, "F" = fuzz_trait,
                      "B" = bin_trait, "N" = nom_trait, "C" = circ_trait)
    
    # Remove NULL data frames
    not_null  <- unlist(lapply(all_trait, function(x) {
      ifelse(!is.null(x), TRUE, FALSE) 
    }))
    all_trait <- all_trait[not_null]
    
    ktab_list <- ade4::ktab.list.df(all_trait)
    
    
    # Compute the final distance using the retained categories:
    
    if (scaling == "scaledBYrange") {
      ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),
                                   option = "scaledBYrange")
    }
    
    if (scaling == "scaledBYsd") {
      ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),
                                   option = "scaledBYsd")
    }
    
    if (scaling == "noscale") {
      ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),
                                   option = "noscale")
    }
  }
  
  # if some species have a functional distance equal to 0, tell the user that it
  # could be a god idea to gather species into FEs:
  
  if (min(ktab_dist) == 0) {
    warning("Functional distance between some species is equal to 0. You can ",
            "choose to gather species into Functional Entities gathering ",
            "species with similar traits values.")
  }
  
  return(ktab_dist)
}