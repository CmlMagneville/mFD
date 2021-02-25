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
#'  \item \strong{trait_weight}: Optional, a numeric vector of length n (traits
#'  number) to specify a weight for each trait.
#'   }
#'   
#' @param metric the distance to be computed:
#'   `euclidean`, the Euclidean distance, 
#'   `gower`, the Classical Gower distance as defined by Gower (1971), extent by
#'   de Bello _et al._ (2021) and based on the \code{\link[gawdis]{gawdis}} 
#'   function.
#'   
#' @param scale_euclid only when computing euclidean distance a string value to
#'   compute (or not) scaling of quantitative traits using the
#'   \code{\link{tr.cont.scale}} function.
#'   Possible options are:
#'   `range` (standardize by the range: 
#'   \eqn{({x' = x - min(x) )} / (max(x) - min (x))})
#'   `center` (use the center transformation: \eqn{x' = x - mean(x)}), 
#'   `scale` (use the scale transformation: \eqn{x' = \frac{x}{sd(x)}}),
#'   `scale_center` (use the scale-center transformation: 
#'   \eqn{x' = \frac{x - mean(x)}{sd(x)}}), or
#'   `noscale` traits are not scaled
#'   Default is `scale_center`.
#'   
#' @param ordinal_var a character string specifying the method to be used for
#'   ordinal variables (i.e. ordered).
#'    `classic` simply treats ordinal variables as continuous variables;
#'    `metric` refers to Eq. 3 of Podani (1999); 
#'    `podani` refers to Eqs. 2a-b of Podani (1999), 
#'   Both options convert ordinal variables to ranks. Default is `classic`.
#'  
#' @param weight_type the type of used method to weight traits. 
#'   `user` – user defined weights in tr_cat, 
#'   `equal` – all traits having the same weight.
#'   More methods are available using \code{\link[gawdis]{gawdis}} from `gawdis` 
#'   package. To compute gower distance with fuzzy trait and traits weight 
#'   please refer to \code{\link[gawdis]{gawdis}}. Default is `equal`.
#' 
#' @param stop_if_NA a logical value to stop or not the process if the
#'   `sp_tr` data frame contains NA. Functional measures are sensitive to
#'   missing traits. For further explanations, see the Note section.
#'   Default is `TRUE`.
#'
#' @return a `dist` object containing distance between each pair of species.
#'
#' @note If the `sp_tr` data frame contains `NA` you can either
#'   chose to compute anyway functional distances (but keep in mind that
#'   **Functional measures are sensitive to missing traits!**) or you can
#'   delete species with missing or extrapolate missing traits (see
#'   Johnson _et al._ (2020)).
#'
#' @references 
#' de Bello _et al._ (2021) Towards a more balanced combination of multiple
#'   traits when computing functional differences between species. 
#'   _Method in Ecology and Evolution_, accepted, 
#'   doi: https://doi.org/10.1111/2041-210X.13537.
#' Gower (1971 ) A general coefficient of similarity and some of its 
#'   properties. _Biometrics_, **27**, 857-871.\cr
#' Johnson _et al._ (2020) Handling missing values in trait data. 
#'   _Global Ecology and Biogeography_, **30**, 51-62.\cr
#' Podani (1999) Extending Gower's general coefficient of similarity to ordinal 
#'   characters, _Taxon_, **48**, 331-340.
#'     
#' @author Nicolas Loiseau & Sebastien Villeger
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load Species x Traits data
#' data("fruits_traits", package = "mFD")
#'
#' # Load Traits x Categories data
#' data("fruits_traits_cat", package = "mFD")
#' 
#' # Remove fuzzy traits for this example and thus remove lat column:
#' fruits_traits     <- fruits_traits[ , -c(6:8)]
#' fruits_traits_cat <- fruits_traits_cat[-c(6:8), ]
#' fruits_traits_cat <- fruits_traits_cat[ , -3]
#' 
#' # Compute Functional Distance
#' sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
#'                                   tr_cat        = fruits_traits_cat,
#'                                   metric        = "gower",
#'                                   scale_euclid  = "scale_center",
#'                                   ordinal_var   = "classic",
#'                                   weight_type   = "equal",
#'                                   stop_if_NA    = TRUE)
#' }

funct.dist <- function(sp_tr,
                       tr_cat,
                       metric,
                       scale_euclid  = "scale_center",
                       ordinal_var = "classic",
                       weight_type = "equal",
                       stop_if_NA  = TRUE) {
  
  ## Check Inputs ----
  
  # with generic functions:
  check.sp.tr(sp_tr, tr_cat, stop_if_NA)
  check.tr.cat(tr_cat)
  check.nominal(tr_cat, sp_tr)
  check.ordinal(tr_cat, sp_tr)
  check.circular(tr_cat, sp_tr)
  check.continuous(tr_cat, sp_tr)
  check.fuzzy(tr_cat, sp_tr)
    
  # checks associated with funct.dist:
  
  metric      <- match.arg(metric, c("euclidean", "gower"))
  weight_type <- match.arg(weight_type, c("equal", "user"))
  ordinal_var <- match.arg(ordinal_var, c("classic", "metric", "podani"))
  
  
  if (weight_type == "user" && sum(colnames(tr_cat) %in% "trait_weight") == 0) {
    stop("A fourth colunm trait_weight must be added in tr_cat to specify ", 
         "weight for each trait.")
  }
  
  if (sum(colnames(tr_cat) %in% "trait_weight") != 0) { 
    
    weight_type <- "user"
    tr_weights  <- tr_cat$trait_weight
    
    if (any(tr_cat$trait_weight <= 0)){
      stop("Trait weight cannot be negative.")
    }
    
  } else { 
    
    tr_weights  <- NULL 
  }
  
  if (any(tr_cat$"trait_type" == "F") && !is.null(tr_weights)) {
    stop("To compute Gower distance with fuzzy traits and trait weights ", 
         "please use gawdis package.")
  }
  
  ## Compute Distances ----
  
  # . . . Euclidean distance
  
  if (metric == "euclidean") {
    
    if (any(tr_cat$"trait_type" == "N")){
      stop("At least one trait is nominal.Species x traits data frame must ", 
           "contain only numerical variables.")
    }
    
    if (any(tr_cat$"trait_type" == "F")){
      stop("At least one trait is fuzzy.Gower distance should be used to ", 
           "consider fuzzy trait.")
    }
    
    for (k in tr_cat$"trait_name"[which(tr_cat$"trait_type" == "Q")]) {
      if (!is.numeric(sp_tr[ , k])) {
        warning("To compute euclidean distance species x traits data frame ", 
                "must contain only numerical variables.")
        stop("Trait '", k, "' is supposed to be continuous but is not ",
             "described with a 'numeric' variable.")
      }
    }
    
    scale_euclid <- match.arg(scale_euclid, c("scale_center","range", "center", 
                                              "scale", "noscale"))
    
    if (scale_euclid != "noscale") {
      sp_tr <- tr.cont.scale(sp_tr, std_method = scale_euclid)
    }
    
    tab_dist <- stats::dist(sp_tr, method = "euclidean")
  }
  
  # . . . Gower distance  
  
  if (metric == "gower") {
    
    if (any(apply(sp_tr, 2, is.numeric))) {
      warning("Species x traits data frame contains only numerical variables. ", 
              "Euclidean can be used.")
      # Not a stop because user may want to weight traits
    }
    
    
    # Compute the final distance using the retained categories and weight:      
    
    #if any fuzzy traits   
    if (any(tr_cat$"trait_type" == "F")) {
      
      # Select the fuzzy traits
      fuzz_cat <- table(tr_cat[tr_cat$"trait_type" == "F", ]$"fuzzy_name")
      
      # Fill all fuzzy names of traits
      tr_cat[is.na(tr_cat$fuzzy_name), ]$fuzzy_name <- tr_cat[is.na(
                                                tr_cat$fuzzy_name), ]$trait_name 
      
      tab_dist <- gawdis::gawdis(sp_tr,
                                 W                 = tr_weights, 
                                 asym.bin          = NULL, 
                                 ord               = ordinal_var,
                                 w.type            = weight_type, 
                                 groups            = tr_cat$fuzzy_name,
                                 groups.weight     = FALSE, 
                                 fuzzy             = c(names(fuzz_cat)), 
                                 opti.getSpecDists = NULL,
                                 opti.f            = NULL,
                                 opti.min.weight   = 0.01,
                                 opti.max.weight   = 1,
                                 opti.maxiter      = 300,
                                 silent            = TRUE)  
      
    } else { 
      
      #without fuzzy traits   
      tab_dist <- gawdis::gawdis(sp_tr,
                                 W                 = tr_weights, 
                                 asym.bin          = NULL, 
                                 ord               = ordinal_var,
                                 w.type            = weight_type, 
                                 groups            = NULL,
                                 groups.weight     = FALSE, 
                                 fuzzy             = NULL, 
                                 opti.getSpecDists = NULL,
                                 opti.f            = NULL,
                                 opti.min.weight   = 0.01,
                                 opti.max.weight   = 1,
                                 opti.maxiter      = 300,
                                 silent            = TRUE) 
    }
  }    
  
  # if some species have a functional distance equal to 0, tell the user that it
  # could be a god idea to gather species into FEs:
  
  if (min(tab_dist) == 0) {
    warning("Functional distance between some species is equal to 0. You can ",
            "choose to gather species into Functional Entities gathering ",
            "species with similar traits values.")
  }
  
  return(tab_dist)
} 
