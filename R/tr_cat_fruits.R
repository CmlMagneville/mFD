#' @title Dataframe gathering traits  for use in mFD example and tutorials
#' @description A dataframe summarizing traits names with traits names on the first column, 
#'  traits_types on the second column and a third column summarizing name of
#'  fuzzy-coded trait to which 'sub-trait' belongs (if trait is not fuzzy,
#'  ignored so could be trait name or NA).
#'   
#' @format A dataframe with 8 rows and 3 columns
#' \describe{
#'   \item{tr}{8 traits types}
#'   \item{cat}{first column: traits names, second column: traits types, third column: fuzzy traits names}
#' }
#'
"tr_cat_fruits"