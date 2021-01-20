# Function to compute distance between species
#
# Authors: Nicolas Loiseau & Sébastien Villéger
#
#

# ------------------------------------------------------------------------------

#' Compute functional distance between species
#'
#' For a given combination of traits, this function returns the functional
#' distance matrix between species
#'
#' @param sp_tr a \strong{dataframe} of traits values (columns) for each species (rows)
#'
#' @param tr_cat a \strong{dataframe} containing three columns for each trait (rows):
#'  \itemize{
#'  \item \strong{trait_name}: names of all traits as in \code{sp_tr} data.frame
#'  \item \strong{trait_type}: category codes for each trait as followed:
#'  \emph{N} for Nominal traits (factor variable), \emph{O} for Ordinal traits
#'  (ordered variable), \emph{C} for Circular traits (integer values), \emph{Q}
#'  for quantitative traits (numeric values) that is allowed \strong{only} if
#'  there are at least 2 species with the same value and \emph{F} for fuzzy
#'  traits (described with several values defined with several column).
#'  }
#'  An option is to add a fourth column with a numeric vector of length n
#'  (traits number) to specify a weight for each trait. The name of the fourth
#'  column has to be 'weight'.
#'
#' @param dist_metric if euclidean, use euclidean distance,
#'                        classical_gower: use classical gower distance \emph{(J. C. (1971))}
#'                        kgower: gower distance modified by Pavoine \emph{Pavoine, S. et al (2009)}
#' Based on the \code{\link[ade4]{dist.ktab}} function: "If kgower
#' is chosen user have several choice to scale quantitative variables A string
#' that can have three values: either \code{"scaledBYrange"} if the quantitative
#' variables must be scaled by their range, or \code{"scaledBYsd"} if they must be
#' scaled by their standard deviation, or \code{"noscale"} if they should not be
#' scaled. This last option can be useful if the values have already been
#' normalized by the known range of the whole population instead of the observed
#' range measured on the sample. If x contains data from various types, then the
#' option \code{"scaledBYsd"} is not suitable (a warning will appear if the option
#' selected with that condition)."
#' 
#' @param stop_if_NA a \strong{logical value} to stop or not the process if the
#'   \code{sp_tr} data frame contains NA. Functional measures are sensitive to
#'   missing traits. For further explanations, see the Note section.
#'   Default: stop_if_NA = TRUE
#'
#' @return a dist object containing distance between each pair of species
#'
#' @section Notes: If the \code{sp_tr} data frame contains NA you can either
#'   chose to compute anyway functional distances (but keep in mind that
#'   \strong{Functional measures are sensitive to missing traits!}) or you can
#'   delete species with missing or extrapolate missing traits (see
#'   \emph{Johnson et al. (2020)})
#' 
#' @examples
#' load(system.file("extdata", "sp_tr_fruits_df", package = "mFD"))
#' sp_tr <- sp_tr[, -c(6:8)]
#' load(system.file("extdata", "sp_tr_cat_fruits_df", package = "mFD"))
#' sp_tr_cat <- sp_tr_cat[-c(6:8), ]
#' mFD::funct.dist(sp_tr, sp_tr_cat, dist_metric = "classical_gower", stop_if_NA = TRUE)
#'
#' @export


funct.dist <- function(sp_tr, tr_cat, dist_metric, stop_if_NA = TRUE) {
  
  if (any(is.na(sp_tr)) == TRUE) {
    
    if (stop_if_NA == TRUE) {
      stop("Error: Species-traits dataframe contains NA. 
      If you want to continue with missing traits (Be careful: Functional measures 
      are sensitive to missing traits), set 'stop_if_NA' parameter to FALSE. 
      Otherwise you can delete species with missing or extrapolate missing traits (Johnson et al. (2020)")
    }
    
  }
  
  if (! is.data.frame(sp_tr))  {
    stop("Error: Your species-traits data must be gathered in a matrix")
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
  
  tr_nm <- names(sp_tr)
  if (any(tr_nm != tr_cat$trait_name)) {
    stop("Error: Trait names differ between species*traits data.frame and
             traits*category data.frame. Please check.")
  }
  
  if(ncol(tr_cat) == 4)  cat("tr_cat has 4 columns, if you use classical_gower
                             traits will be weighted", "\n") # end of checking weight information
  
  if(dist_metric == "classical_gower") {
    
    if(ncol(tr_cat) == 4) {
      ktab_dist <- cluster::daisy(sp_tr, "gower", weights = tr_cat$weight)
      
    } else {
      ktab_dist <- cluster::daisy(sp_tr, "gower")
    }
    
  } # compute Gower distance based on Gower, J. C. (1971)
  
  else if (dist_metric == "euclidean") {
    
    if (!is.numeric(sp_tr)) {
      stop("Some traits are not numerical")
    }   # end of checking if variables are numerical for euclidean distance
    
    # compute Euclidean distance:
    
    ktab_dist <- stats::dist(sp_tr, method = "euclidean")
    
  } else {
    
    # compute modifed Gower distance bases on Pavoine, S. et al (2009)
    
    # Need to prepare functional trait in function of their nature:
    
    # Quantitative traits:
    
    quant_trait <- NULL
    
    if (any(tr_cat$trait_type == "Q")) {
      quant_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "Q"],
                           drop = FALSE]
    }
    
    # Ordinal Traits:
    
    ord_trait <- NULL
    
    if (any(tr_cat$trait_type == "O")) {
      ord_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "O"],
                         drop = FALSE]
    }
    
    # Circular traits:
    
    circ_trait <- NULL
    
    if (any(tr_cat$trait_type == "C")) {
      circ_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "C"],
                          drop = FALSE]
      
      circ_trait <- ade4::prep.circular(circ_trait,1,12)
    }
    
    # Fuzzy traits (basically several categories that are considered a single trait):
    
    fuzz_trait <- NULL
    
    if (any(tr_cat$trait_type == "F")) {
      
      # Select the fuzzy traits
      fuzz_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "F"],
                          drop = FALSE]
      
      # Count the number of fuzzy categories
      fuzz_cat <- table(tr_cat[tr_cat$trait_type == "F", ]$fuzzy_name)
      
      # Order the trait names based on the order of the categories
      fuzz_names_ordered <- unlist(lapply(names(fuzz_cat),
                                          function(x) {
                                            tr_cat$trait_name[tr_cat$fuzzy_name == x & !is.na(tr_cat$fuzzy_name)]
                                          }))
      
      # Reorder the traits according to the names
      fuzz_trait <- fuzz_trait[, fuzz_names_ordered]
      
      fuzz_trait <- ade4::prep.fuzzy(fuzz_trait,
                                     col.blocks = as.numeric(fuzz_cat),
                                     labels     = names(fuzz_cat))
    }
    
    # Binary traits:
    
    bin_trait <- NULL
    
    if (any(tr_cat$trait_type == "B")) {
      bin_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "B"], drop = FALSE]
      
      bin_trait <- ade4::prep.binary(bin_trait, col.blocks = ncol(bin_trait))
    }
    
    # Norminal traits:
    
    nom_trait <- NULL
    
    if (any(tr_cat$trait_type == "N")) {
      nom_trait <- sp_tr[, tr_cat$trait_name[tr_cat$trait_type == "N"], drop = FALSE]
    }
    
    # Combine all traits (remove 'NULL' data.frames using 'plyr::compact()')
    all_trait <- plyr::compact(list(Q = quant_trait,
                                    O = ord_trait,
                                    F = fuzz_trait,
                                    B = bin_trait,
                                    N = nom_trait,
                                    C = circ_trait))
    
    ktab_list <- ade4::ktab.list.df(all_trait)
    
    
    # Compute the final distance using the retained categories
    switch(utils::menu(c("Quantitative variables must be scaled by their range : scaledBYrange",
                         
                         "Quantitative variables must be scaled by their standard deviation :scaledBYsd",
                         
                         "Quantitative variables should not be scaled : noscale"),
                       
                       title="Do you want to scale quantitative varibles ?") + 1,
           cat("Nothing done\n"),
           
           ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),option = "scaledBYrange"),
           
           ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),option = "scaledBYsd"),
           
           ktab_dist <- ade4::dist.ktab(ktab_list, names(all_trait),option = "noscale")) #End of the choice of the method to scale
  }
  
  # if some species have a functional distance equal to 0, tell the user that it
  # could be a god idea to gather species into FEs:
  
  if (min(ktab_dist) == 0) {
    warning("Functional distance between some species is equal to 0. You can choose to gather species into Functional Entities gathering species with similar traits values")
  }
  
  return(ktab_dist)
  
}

