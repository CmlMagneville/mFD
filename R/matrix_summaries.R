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
#' # Load Species x Traits data
#' data("sp_tr_fruits", package = "mFD")
#'
#' # Load Traits x Categories data
#' data("tr_cat_fruits", package = "mFD")
#'
#' # Summarize Species x Traits data
#' mFD::sp.tr.summary(tr_cat = tr_cat_fruits, sp_tr = sp_tr_fruits)


sp.tr.summary <- function(tr_cat, sp_tr) {


  ## Check Inputs ----

  check.sp.tr(sp_tr, tr_cat, stop_if_NA = TRUE)


  ## Checks Traits Formats ----

  check.nominal(tr_cat, sp_tr)
  check.ordinal(tr_cat, sp_tr)
  check.circular(tr_cat, sp_tr)
  check.continuous(tr_cat, sp_tr)
  check.fuzzy(tr_cat, sp_tr)


  ## Retrieve Traits Informations ----

  tr_type        <- tr_cat$"trait_type"
  names(tr_type) <- tr_cat$"trait_name"


  ## Transform Fuzzy-coded Traits ----

  if ("F" %in% tr_cat$"trait_type") {

    # retrieve names of fuzzy-coded traits:
    nm_fuzzy <- unique(stats::na.omit(tr_cat$"fuzzy_name"))

    # Update names, type and number of traits
    tr_nm   <- c(tr_cat$"trait_name"[tr_cat$"trait_type" != "F"], nm_fuzzy)
    tr_nb   <- length(tr_nm)
    tr_type <- c(tr_cat$"trait_type"[tr_cat$"trait_type" != "F"],
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
    sp_non_conttr <- sp_tr[ , lapply(sp_tr, function(x) !is.numeric(x))]
    non_conttr_modalities_list <- lapply(sp_non_conttr, unique)

    return(list("tr_summary_list"     = tr_summary_list,
                "tr_types"            = tr_types,
                "non_conttr_mod_list" = non_conttr_modalities_list))

  } else {                            # fuzzy coded traits

    # Table with traits summary for non fuzzy traits
    tr_summary_list <- summary(sp_tr[ , tr_cat$"trait_name"[
      which(tr_cat$"trait_type" != "F")], ])

    # Vector containing traits type for non fuzzy traits
    tr_types <- sapply(sp_tr[ , tr_cat$"trait_name"[
      which(tr_cat$"trait_type" != "F")], ], class)

    # List containing modalities for non continuous traits
    sp_non_conttr <- sp_tr[ , tr_cat$"trait_name"[
      which(tr_cat$"trait_type" != "F" & tr_cat$"trait_type" != "Q")], ]
    mod_list <- lapply(sp_non_conttr, unique)

    # For fuzzy traits
    tr_summary_fuzzy_list <- summary(sp_tr[ , tr_cat$"trait_name"[
      which(tr_cat$"trait_type" == "F")], ])

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
#'   \item{asb_sp_w_occ}{a matrix with species occurrences in each
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
#' # Load Assemblages x Species Matrix
#' data("asb_sp_w_fruits", package = "mFD")
#' 
#' # Summarize Assemblages Data
#' mFD::asb.sp.summary(asb_sp_w = asb_sp_w_fruits)

asb.sp.summary <- function(asb_sp_w) {


  ## Check Input ----

  check.asb.sp.w(asb_sp_w)

  # Species occurrence data frame
  asb_sp_w_occ <- replace(asb_sp_w, asb_sp_w != 0, 1)

  # Vector containing the number of occurrences of each species
  nbocc_sp <- apply(as.data.frame(asb_sp_w), 2, sum)

  # Vector containing total abundance/biomass per assemblage
  asb_totab <- apply(as.data.frame(asb_sp_w), 1, sum)

  # Vector containing species richness of each assemblage
  asb_sp_wrichn <- apply(asb_sp_w_occ, 1, sum)


  # Construction of a list containing vectors with species names for each...
  # ... assemblage
  L <- list()

  for (i in 1:nrow(asb_sp_w_occ)) {

    data  <- asb_sp_w_occ[i, ]
    data2 <- data[which(apply(as.data.frame(data), 2, max) == TRUE)]
    rownames(data2) <- rownames(asb_sp_w_occ[i, ])
    L[[i]] <- as.vector(data2)
  }


  ## Function Return ----

  list("asb_sp_occ"       = as.matrix(asb_sp_w_occ),
       "tot_ab_all_sp"    = nbocc_sp,
       "tot_ab_all_asb"   = asb_totab,
       "sp_richn_all_asb" = asb_sp_wrichn,
       "sp_nm_asb"        = L)
}
