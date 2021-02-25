#' Dataset: Fruits Traits Informations
#' 
#' This dataset summarizes information about the 6 traits used in the 
#' `fruits_traits` dataset.
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
#' Traits in this dataset correspond to columns (traits) of the `fruits_traits`
#' dataset.
#'
#' @seealso `fruits_traits`, `baskets_fruits_weights`
#' 
#' @examples 
#' \dontrun{
#' # Load Traits Information
#' data("fruits_traits_cat", package = "mFD")
#' fruits_traits_cat
#' 
#' # Load Species x Traits Data Frame
#' data("fruits_traits", package = "mFD")
#' 
#' # Summarize Species x Traits Data
#' mFD::sp.tr.summary(tr_cat = fruits_traits_cat, sp_tr = fruits_traits)
#' }

"fruits_traits_cat"