#' Compute Functional Entities composition based on a Species x Traits matrix
#'
#' @param sp_tr a data frame containing species as rows and traits as columns.
#'
#' @param tr_cat a data frame containing three columns for each trait (rows):
#' \itemize{
#'   \item \strong{trait_name}: names of all traits as in \code{sp_tr} data 
#'     frame;
#'   \item \strong{trait_type}: category codes for each trait as followed: 
#'     \emph{N} for Nominal traits (factor variable), 
#'     \emph{O} for Ordinal traits (ordered variable), 
#'     \emph{C} for Circular traits (integer values), 
#'     \emph{Q} for Quantitative traits (numeric values) that is allowed
#'       \strong{only} if there are at least 2 species with the same value, and
#'     \emph{F} for fuzzy-coded traits (i.e. described with several 
#'       'sub-traits');
#'   \item \strong{fuzzy_name} name of fuzzy-coded trait to which 'sub-trait'
#'     belongs (if trait is not fuzzy, ignored so could be trait name or NA).
#' }
#'
#' @param fe_nm_type a character string referring to the type of naming
#'   functional entities. Two possible values: \emph{"fe_rank"} (FE are named
#'   after their decreasing rank in term of number of species \emph{i.e.} fe_1
#'   is the one gathering most species) and \emph{"tr_val"} (FE are named after
#'   names of traits and of trait values for each FE, if possible, \emph{see
#'   details below}). Default: `fe_nm_type = "fe_rank"`.
#'
#' @param check_input a logical value indicating whether key features the 
#'   inputs are checked (e.g. class and/or mode of objects, names of rows 
#'   and/or columns, missing values). If an error is detected, a detailed 
#'   message is returned. Default: `check.input = TRUE`.
#'
#' @details `fe_nm_type = "tr_val"` is allowed \strong{only} if: 
#' \itemize{
#'   \item there are less than 7 traits;
#'   \item none of them is fuzzy-coded (so that names are not too long)
#'   \item all trait names and all trait values have different 2 first letters 
#' }
#' 
#' If these 3 conditions are met, names of Functional Entities are made as a 
#' character string of up to 2 letters for trait name in upper case font then 
#' up to 2 letters for trait value in lower case font, separated by "_" between 
#' traits. Trait names are abbreviated to a single letter whenever possible.
#' \emph{Examples:} ("TAc2_TBxx_TCyy", "TAc3_TBff_TCyy") or 
#' ("A2_Bx_Cy", "A3_Bf_Cy")
#'
#' @return A list of objects containing: 
#' \itemize{
#'   \item \strong{fe_nm}: a vector with names of all FE (following 
#'     fe_nm_type). FE are ordered according to the decreasing number of 
#'     species they gather. 
#'   \item \strong{sp_fe}: a vector containing for each species the name of the 
#'     FE it belongs to. FE order is done according to decreasing number of 
#'     species.
#'   \item \strong{fe_tr}: a data frame containing traits values (variables in 
#'     columns) for each FE (rows). FE order is done according to decreasing 
#'     number of species.
#'   \item \strong{fe_nb_sp}: a vector with species number per FE. If all FE 
#'     have only one species, a warning message is returned. FE are ordered 
#'     according to the decreasing number of species they gather. 
#'   \item \strong{details_fe}: a list containing: \emph{fe_codes} a vector 
#'     containing character referring to traits values (like a barcode) with 
#'     names as in \code{fe_nm_type} and sorted according to \code{fe_nb_sp}; 
#'     \emph{tr_uval} a list containing for each trait a vector of its unique 
#'     values or a data frame for fuzzy-coded traits; \emph{fuzzy_E} a list 
#'     with for each fuzzy-coded trait a data frame with names of entities (E) 
#'     and names of species (sp); \emph{tr_nb_uval} a vector with number of 
#'     unique values per trait (or combinations for fuzzy-coded traits); 
#'     \emph{max_nb_fe} the maximum number of FE possible given number of 
#'     unique values per trait.
#' }
#'
#' @author Sebastien Villeger, Nicolas Loiseau, and Camille Magneville
#'
#' @export
#' 
#' @examples
#' # Load species traits data:
#'  data("fruits_traits", package = "mFD")
#' 
#' # Transform species traits data:
#' # Only keep the first 4 traits to illustrate FEs:
#'  fruits_traits <- fruits_traits[ , c(1:4)]   
#'
#' # Load trait types data:
#'  data("fruits_traits_cat", package = "mFD")
#' 
#' # Transform the trait types data to only keep traits 1 - 4:
#'  fruits_traits_cat <- fruits_traits_cat[c(1:4), ]
#'
#' # Gather species into FEs:
#' ## gathering species into FEs (FEs named according to the decreasing...
#' ## ...  number of species they gather):
#'  sp_FEs <- mFD::sp.to.fe(
#'       sp_tr      = fruits_traits, 
#'       tr_cat     = fruits_traits_cat, 
#'       fe_nm_type = "fe_rank")
#' 
#' ## display FEs names:
#' sp_FEs$fe_nm
#' 
#' ## display for each species the name of the FE it belongs to:
#' sp_FEs$sp_fe
#' 
#' ## display trait values for each FE:
#' sp_FEs$fe_tr
#'  
#' ## display the number of species per FEs:
#' sp_FEs$fe_nb_sp


sp.to.fe <-  function(sp_tr, tr_cat, fe_nm_type = "fe_rank", 
                      check_input = TRUE) {
  
  
  # define key parameters used in the function:
  # define species names:
  sp_nm <- row.names(sp_tr)
  # define species number
  sp_nb <- length(sp_nm)
  
  
  # retrieve traits information (fuzzy coded traits are modified thereafter):
  # get traits names:
  tr_nm <- names(sp_tr)
  # get traits number:
  tr_nb <- length(tr_nm)
  # get traits type:
  tr_type <- c(tr_cat$trait_type)
  names(tr_type) <- tr_nm
  
  
  ## check_inputs:
  
  if (sp_nb < 2) {
    stop("There must be at least 2 species in species x traits data frame.")
  }
  
  if (length(tr_nm) < 2) {
    stop("There must be at least 2 traits in the species x traits data frame.")
  }
  
  if (any(is.na(sp_tr))) {
    stop("The species x traits data frame contains NA. Please check.")
  }
  
  if (! is.data.frame(sp_tr))  {
    stop("Your species x traits data must be gathered in a matrix.")
  }
  
  if (is.null(colnames(sp_tr))) {
    stop("No column names provided in traits table. Analysis will not go ",
         "through, please add traits as column names.")
  }
  
  if (is.null(rownames(sp_tr))) {
    stop("No row names provided in traits table. Analysis will not go ",
         "through. Please add species names as row names.")
  }
  
  if (any(is.na(tr_cat$trait_type))) {
    stop("Trait type in traits x category data frame contains NA. Please ",
         "check and specify type of all traits.")
  }
  
  if (any(tr_nm != tr_cat$trait_name)) {
    stop("Trait names differ between species x traits data frame and traits ",
         "x category data.frame. Please check.")
  }
  
  if (any(is.na(tr_cat$trait_type))) {
    stop("Trait type in traits x category data frame contains NA. Please ", 
         "check and specify type of all traits.")
  }
  
  if (any(!(tr_cat$trait_type %in% c("N", "O", "C", "Q", "F")))) {
    stop("Trait type in traits x category should be among 'N', 'O', 'C', ",
         "'Q', 'F'. Please check type of all traits.")
  }
  
  # check nominal traits:
  check.nominal(tr_cat, sp_tr)
  
  # check that ordinal traits have ordered categories:
  check.ordinal(tr_cat, sp_tr)
  
  # check that circular traits are coded with integer:
  check.circular(tr_cat, sp_tr)
  
  # check that continuous trait have non unique values:
  check.continuous(tr_cat, sp_tr)
  
  # check fuzzy coded traits:
  if ("F" %in% tr_cat$trait_type) {
    # retrieve names of fuzzy-coded traits:
    nm_fuzzy <- unique(stats::na.omit(tr_cat$fuzzy_name))
    # check that fuzzy traits are described with more than one variable:
    for (k in nm_fuzzy) {
      # get the names of variables for fuzzy trait k:
      var_k <- tr_cat$trait_name[which(tr_cat$fuzzy_name == k)]
      
      # check that there are at least 2 columns per fuzzy coded trait:
      if (length(var_k) < 2) {
        stop("Fuzzy-coded trait '", k, "' is described with a single ",
             "variable, consider changing its type to 'nominal'")
      }
      
      # check that variables are continuous:
      if (any(apply(sp_tr[, var_k ], 2, is.numeric) == FALSE)) {
        stop("Fuzzy-coded trait '", k, "' is not described with 'numeric' ",
             "variables.")
      }
    }
    
    # update names, type and number of traits:
    tr_nm <- c(tr_cat$trait_name[tr_cat$trait_type != "F"], nm_fuzzy)
    tr_nb <- length(tr_nm)
    tr_type <- c(tr_cat$trait_type[tr_cat$trait_type != "F"],
                 rep("F", length(nm_fuzzy)))
    names(tr_type) <- tr_nm
  }
  
  
  ## retrieve combination of unique values of traits:
  
  # code a function to paste name of a trait (uppercase) ...
  #... to its values (lowercase):
  paste.tr.val <- function(tr, x) {
    paste0(toupper(x), tolower(as.character(tr[, x])))
  }
  
  # create a list to store outputs:
  tr_uval <- list()
  
  # create a vector to store functional barcodes of each specis...
  # ... (combinations of trait name and values):
  sp_fbarcode <- rep("", sp_nb)
  
  # complete the list and the vector for non fuzzy traits:
  tr_notfuzzy <- names(tr_type)[tr_type != "F"]
  if (length(tr_notfuzzy) > 0) {
    for (k in tr_notfuzzy) {
      # get traits unique values:
      tr_uval[[k]] <- unique(as.character(sp_tr[, k]))
      # retrieve traits values:
      sp_fbarcode <- paste(sp_fbarcode, paste.tr.val(sp_tr, k) , sep = "_")
    }
  }
  
  # group species into a nominal-like trait for fuzzy traits:
  if ("F" %in% tr_type) {
    # create a list to store functional entities names and their species:
    fuzzy_E <- list()
    # complete fuzzy_E list for all fuzzy traits, output list and vector:
    for (k in nm_fuzzy) {
      # retrieve variables for fuzzy trait k for each species:
      fuzzy_k <- sp_tr[ , tr_cat$trait_name[which(tr_cat$fuzzy_name == k)]]
      # create a vector for first variable:
      val_k <- rep("", sp_nb)
      # loop on other variables:
      for (z in colnames(fuzzy_k)) {
        val_k <- paste(val_k, paste.tr.val(fuzzy_k, z), sep = "_")
      }
      # removing useless first "_":
      val_k <- substr(val_k, 2, nchar(val_k))
      # create a data frame with entities and species:
      val_k_df <- data.frame(E = val_k, sp = row.names(fuzzy_k))
      fuzzy_E[[k]] <- val_k_df
      # get the unique combinations of values for fuzzy trait k:
      u_val_k <- as.character(unique(val_k_df$E))
      # add code of fuzzy entities to barcode of all traits: ...
      # ... get all combinations of traits to draw FEs:
      sp_fbarcode <- paste(sp_fbarcode, u_val_k, sep = "_")
      # there are variable values for each fuzzy entity so ...
      # ... pick values from first species belonging to the fuzzy entity...
      # ... so we can get all combinations of "traits" forming ...
      # ... the fuzzy trait k:
      E_val_k <- lapply(u_val_k, function(x) {
        fuzzy_k[val_k_df$sp[which(val_k_df$E == x)[1]], ]})
      # replace names by entities names
      names(E_val_k) <- u_val_k
      # gather vector as a dataframe to store unique values:
      tr_uval[[k]] <- as.data.frame(do.call(rbind, E_val_k))
    }
  }
  
  # create a vector to store the number of unique values per trait:
  tr_nb_uval <- rep(NA, length(tr_nm))
  names(tr_nb_uval) <- tr_nm
  # store the number of unique value per trait:
  for (k in tr_nm) {
    if (tr_type[k] != "F"){
      tr_nb_uval[k] <- length(tr_uval[[k]])
    }
    if (tr_type[k] == "F") {
      tr_nb_uval[k] <- nrow(tr_uval[[k]])
    }
  }
  
  # retrieve the maximum number of possible FE given the number of ...
  # ... unique values per trait:
  max_nb_fe <- prod(tr_nb_uval)
  
  
  ## arrange functional entities:
  
  # remove useless "_" from FEs codes:
  sp_fe <- substr(sp_fbarcode, 2, nchar(sp_fbarcode))
  
  # add species names to vector of species functional barcodes:
  names(sp_fe) <- sp_nm
  
  # get the codes of unique FE:
  fe_codes <- as.character(unique(sp_fe))
  
  # get FE number:
  fe_nb <- length(fe_codes)
  
  # retrieve the number of species per FE:
  fe_nb_sp <- unlist(lapply(fe_codes, function(x) length(which(sp_fe == x))))
  names(fe_nb_sp) <- fe_codes
  
  # return a warning if all species are different (i.e. only 1 sp per FE):
  if (max(fe_nb_sp) == 1) {
    warning("All Functional Entities have a single species.")
  }
  
  # decrease order of FE according to their number of species to order fe... 
  # ... codes:
  # order fe with increasing species number
  order_fe <- order(fe_nb_sp, decreasing = TRUE)
  # reorder vectors with increasing number of species per FE:
  fe_nb_sp_ord <- fe_nb_sp[order_fe]
  # reorder fe_codes with increasing number of species per FE:
  fe_codes_ord <-fe_codes[order_fe]
  
  # get trait values for FEs:
  # pick traits values of the first species (one species) for each FE:
  fe_trait_list <- lapply(fe_codes,
                          function(x) {sp_tr[names(which(sp_fe == x)[1]), ]})
  
  # replace species names by entities names:
  names(fe_trait_list) <- fe_codes
  
  # gather vectors as a dataframe to store unique values:
  fe_trait <- as.data.frame(do.call(rbind, fe_trait_list))
  
  # reorder rows of fe*traits dataframe according to the number of ...
  # ... species per FE:
  fe_trait_ord <- fe_trait[fe_codes_ord, ]
  
  
  ## code Functional Entities names:
  
  # write FEs names based on decreasing number of species per FE if asked:
  if (fe_nm_type == "fe_rank") {
    fe_nm <- paste0("fe_", 1:fe_nb)
  }
  
  # write FEs names based on traits and values names if ...
  # ... asked (not for fuzzy traits):
  if (fe_nm_type == "tr_val") {
    # check that there is no fuzzy-coded trait:
    if ("F" %in% tr_type) {
      stop("Functional entities could not be named according to names of ",
           "traits and of trait values when there is a fuzzy-coded trait. ",
           "Set 'fe_nm_type' to 'fe_rank'.")
    }
    # check that there is less than 6 traits:
    if (tr_nb > 6) {
      stop("Functional entities could not be named according to names of ",
           "traits and of trait values becasue there are more than 6 traits. ",
           "Set 'fe_nm_type' to 'fe_rank'.")
    }
    # set traits code default value to one letter:
    traits_codes <- substr(names(fe_trait_ord), 1, 1)
    names(traits_codes) <- tr_nm
    # if some traits have the same 1st letter, up to 2 first letters are kept:
    if (length(unique(traits_codes) != tr_nb)) {
      traits_codes <- substr(names(sp_tr), 1, 2)
      names(traits_codes) <- tr_nm
    }
    # check that codes for traits are unique so traits can be identified:
    if (length(unique(traits_codes)) != tr_nb) {
      stop("Two first letters of trait names should be unique. Please change ",
           "'fe_nm_type' to 'fe_rank'.")
    }
    # build FE names after check for each trait...
    # ... that levels have unique two first letters:
    # create a vector that will contain FEs names:
    fe_nm <- rep("", fe_nb)
    # for each trait:
    for (t in tr_nm) {
      # get the levels of trait t:
      level_t <- as.character(fe_trait_ord[, t])
      # set default code to one letter:
      level_t_codes <- substr(level_t, 1, 1)
      # if some codes have the same 1st letter, up to 2 first letters are kept:
      if (length(unique(level_t_codes)) < length(unique(level_t))) {
        level_t_codes <- substr(level_t_codes, 1, 2)
      }
      # if some codes have the same two 1st letters, error message:
      if (length(unique(level_t_codes)) != length(unique(level_t))) {
        stop("Two first letters of levels of trait '", t, "' are not unique. ",
             "Please change 'fe_nm_type' to 'fe_rank'.")
      }
      # add trait name and values among FE to the names of FE:
      fe_nm <- paste0(fe_nm, "_", toupper(traits_codes[t]),
                      tolower(level_t_codes))
    }
    # remove the useless first "_":
    fe_nm <- substr(fe_nm, 2, nchar(fe_nm))
  }
  
  # replace FEs barcodes by fe_names in fe_codes, fe_nb_sp and fe_trait...
  # ... which all had the same order as fe_codes (based on fe_nb_sp):
  names(fe_nb_sp_ord) <- fe_nm
  names(fe_codes_ord) <- fe_nm
  rownames(fe_trait_ord) <- fe_nm
  
  # replace FEs codes in sp_fe by proper FE names:
  sp_fe_ord <- sp_fe
  # for each FE code:
  for (z in names(fe_codes_ord)) {
    sp_fe_ord[which(sp_fe_ord == fe_codes_ord[z])] <- z
  }
  
  
  ## return results:
  
  return_list <-list("fe_nm" = fe_nm, "sp_fe" = sp_fe_ord, 
                     "fe_tr" = fe_trait_ord, "fe_nb_sp" = fe_nb_sp_ord,
                     "details_fe" = list("fe_codes" = fe_codes_ord, 
                                         "tr_uval" = tr_uval,
                                         "tr_nb_uval" = tr_nb_uval, 
                                         "max_nb_fe" = max_nb_fe))
  
  if ("F" %in% tr_type) {
    return_list$details_fe$fuzzy_E <- fuzzy_E
  }
  
  return(return_list)
}
