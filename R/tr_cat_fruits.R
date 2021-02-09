#' Fruits Traits Informations
#' 
#' This dataset summarizes information about the 6 traits used in the 
#' `sp_tr_fruits` dataset.
#'   
#' @format A data frame with 8 rows (traits) and the following three columns:
#' \describe{
#'   \item{trait_name}{a character giving the trait name}
#'   \item{trait_type}{a character giving the trait type (**O** for Ordinal 
#'   trait, **N** for Nominal trait, **Q** for Quantitative trait, and **F** for
#'   Fuzzy-coded trait)}
#'   \item{fuzzy_name}{a character giving the name of fuzzy-coded traits 
#'   (i.e. `Use`) to  which 'sub-traits' (i.e. `raw`, `pastry`, and `jam`) 
#'   belongs}
#' }
#' 
#' @note If your dataset does not contain fuzzy trait, the column `fuzzy_name`
#' can be ignored but the first two columns are mandatory.
#' 
#' Traits in this dataset correspond to columns (traits) of the `sp_tr_fruits`
#' dataset.
#'
#' @seealso `sp_tr_fruits`, `asb_sp_w_fruits`
#' 
#' @examples 
#' # Load Traits Information
#' data("tr_cat_fruits", package = "mFD")
#' tr_cat_fruits
#' 
#' # Load Species x Traits Data Frame
#' data("sp_tr_fruits", package = "mFD")
#' 
#' # Summarize Species x Traits Data
#' mFD::sp.tr.summary(tr_cat = tr_cat_fruits, sp_tr = sp_tr_fruits)

"tr_cat_fruits"