#' Dataset: Traits Values of Fruits Species
#' 
#' This dataset represents the value of 6 traits for 15 fruits species. 
#' **Important:** Species names must be specified as the data frame row names 
#' (not in an additional column).
#'   
#' @format A data frame with 25 rows (species) and the following columns 
#' (traits):
#' \describe{
#'   \item{Size}{an ordered factor describing the size of fruits species}
#'   \item{Plant}{an unordered factor describing the type of plant}
#'   \item{Climate}{an ordered factor describing the climate regions}
#'   \item{Seed}{an ordered factor describing the type of seed}
#'   \item{Sugar}{a numeric describing the quantity of sugar}
#'   \item{Use.raw}{an integer (percentage) describing the proportion of a raw 
#'   use (fuzzy trait) of the fruit}
#'   \item{Use.pastry}{an integer (percentage) describing the proportion of a 
#'   pastry use (fuzzy trait) of the fruit}
#'   \item{Use.jam}{an integer (percentage) describing the proportion of a 
#'   jam use (fuzzy trait) of the fruit}
#' }
#'
#' @seealso `fruits_traits_cat`, `baskets_fruits_weights`
#' 
#' @examples 
#' \dontrun{
#' # Load Species x Traits Data Frame
#' data("fruits_traits", package = "mFD")
#' fruits_traits
#' 
#' # Load Traits Information
#' data("fruits_traits_cat", package = "mFD")
#' fruits_traits_cat
#' 
#' # Summarize Species x Traits Data
#' mFD::sp.tr.summary(tr_cat = fruits_traits_cat, sp_tr = fruits_traits)
#' }

"fruits_traits"
