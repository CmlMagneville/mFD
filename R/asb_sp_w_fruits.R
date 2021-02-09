#' Baskets Composition in Fruits Species
#' 
#' This dataset represents the abundance of 25 fruits species in 10 baskets.
#'   
#' @format A matrix of integers with 10 rows (baskets) and 25 columns (species).
#' 
#' @seealso `sp_tr_fruits`, `tr_cat_fruits`
#' 
#' @examples 
#' # Load Assemblages x Species Matrix
#' data("asb_sp_w_fruits", package = "mFD")
#' asb_sp_w_fruits[1:5, 1:5]
#' 
#' # Summarize Assemblages Data
#' mFD::asb.sp.summary(asb_sp_w = asb_sp_w_fruits)

"asb_sp_w_fruits"