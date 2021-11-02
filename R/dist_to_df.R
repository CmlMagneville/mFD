#' Merge distance object(s) into a single data frame
#' 
#' This function merges distance object(s) into a single data frame which rows 
#' are pairs of elements and column(s) distance metric(s). It stands on the 
#' \code{\link[dendextend]{dist_long}} function.
#' 
#' @param list_dist a list of dist object(s). All dist objects should have a
#'   name (e.g. name of distance metric) and the same labels (i.e. names of 
#'   sets between which distance was computed).
#'
#' @return A data frame which first and second columns (names `x1` and `x2`)
#' contain names of the 2 sets involved in each pair, and with one column for
#' each dist object (named after its name in \code{list_dist}.
#' 
#' @author Sebastien Villeger
#' 
#' @export
#' 
#' @examples 
#' # Create dist objects: 
#' dist_A <- round(dist(matrix(runif(1:100), 5, 2, 
#'                       dimnames = list(letters[1:5], NULL))), 2)
#' dist_B <- round(dist(matrix(runif(1:100), 5, 2, 
#'                       dimnames = list(letters[1:5], NULL))), 2)
#' dist_C <- round(dist(matrix(runif(1:100), 5, 2, 
#'                       dimnames = list(letters[1:5], NULL))), 2)
#'
#' # First example with only 1 distance:
#' dist.to.df(list(dA = dist_A))
#'
#' # Second example with 3 distances:
#' dist.to.df(list(d1 = dist_A, d2 = dist_B, d3 = dist_C))

dist.to.df <- function(list_dist) {
  
  # names and number of dist objects:
  dist_nm <- names(list_dist)
  dist_nb <- length(dist_nm)
  
  # labels of the fist dist object:
  dist1_labels <- labels(list_dist[[1]])
  
  
  #### checking inputs #####
  
  # checking list contains only dist
  # objects:
  if (any(unlist(lapply(list_dist, class)) != "dist")) {
    stop("Input 'list_dist' should contain only 'dist' object. Correct using ", 
         "'as.dist()' if necessary.")
  }
  
  # checking input is a list
  if (!is.list(list_dist)) {
    stop("Input 'list_dist' should be a list (even when only one dist object ", 
         "is provided).")
  }
  
  
  # checking all dist objects have names:
  if (is.null(dist_nm) || any(nchar(dist_nm) == 0)) {
    stop("Some of dist objects in 'list_dist' do not all have a name. Name ", 
         "all dist objects within the list (e.g. with distance metric).")
  }
  
  # checking 1st dist objects has labels:
  if (is.null(dist1_labels) || any(nchar(dist1_labels) == 0)) {
    stop("First dist object in 'list_dist' does not have labels. Provide row ", 
         "names as character strings for the sets*variables matrix before ", 
         "computing distance.")
  }
  
  
  #### reference for pairs of sets is first
  #### object of list ####
  
  # applying dist_long:
  df_dist <- dendextend::dist_long(list_dist[[1]])
  
  # reversing and renaming order of
  # columns:
  df_dist <- data.frame(df_dist$cols, df_dist$rows, 
                        df_dist$distance)
  names(df_dist) <- c("x1", "x2", dist_nm[1])
  
  
  #### if other dist object, binding values
  #### with first ones as new column(s) ####
  
  if (dist_nb > 1) {
    
    # loop on dist objects:
    for (k in 2:length(list_dist)) {
      # name of dist objects:
      k_nm <- dist_nm[k]
      
      # names of its labels:
      k_labels <- labels(list_dist[[k_nm]])
      
      # checking same labels than first dist
      # object:
      if (is.null(k_labels) || any(dist1_labels != k_labels)) {
        stop("Element ", k_nm, " does not have the same labels than first ", 
             "dist object ", dist_nm[1], ". Check labels of inputs.")
      }
      
      # applying dist_long:
      df_k <- dendextend::dist_long(list_dist[[k_nm]])
      
      # merging with first dataframe as
      # variable with dist object name:
      df_dist[[k_nm]] <- df_k$distance
      
    }  # end of k
  }  # end of if at least 2 dist objects
  
  # output:
  return(df_dist)
}
