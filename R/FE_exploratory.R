#' Build the assemblage-FEs dataframe from the assemblages-species one 
#' 
#' This function computes an occurrence data frame with assemblages in rows and 
#' Functional Entities (FEs) in columns, from the dataframe of 
#' assemblage-species and the output of the \code{mFD::sp.to.fe()} function.
#' 
#' @param sp_fe list gathering to which FE belongs eac species. It is the output 
#' of the mFD::sp.to.fe() function - $sp_fe. 
#' @param asb_sp_w the assemblage * species data frame with assemblages being 
#' rows and species being columns
#'
#' @return an occurrence dataframe with studied FEs in columns and assemblages 
#' in rows. \strong{Be careful: It's an occurrence data frame, not an 
#' abundance one} 
#' 
#' @author Camille Magneville
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
#'  # Load Assemblages*Species matrix:
#' data('baskets_fruits_weights', package = 'mFD')
#'
#' # Gather species into FEs:
#' ## gathering species into FEs (FEs named according to the decreasing...
#' ## ...  number of species they gather):
#'  sp_FEs_fruits <- mFD::sp.to.fe(
#'       sp_tr      = fruits_traits, 
#'       tr_cat     = fruits_traits_cat, 
#'       fe_nm_type = "fe_rank")
#'       
#' # Get the list which gather to which FE belongs each species:
#' sp_fes_list <- sp_FEs_fruits$sp_fe
#'       
#' # Build the Assemblages*FEs data frame:
#' asb_fes <- mFD::from.spfe.to.feasb(sp_fe = sp_fes_list,
#'                                    asb_sp_w = baskets_fruits_weights)
#' asb_fes
#' 
#'
#' @export
#'

from.spfe.to.feasb <- function(sp_fe, asb_sp_w) {
  
  
  # Retrieve the FEs names:
  fe_nm <- unique(sp_fe)
  
  # Create a dataframe that will contain fe*asb data:
  asb_fe <- as.data.frame(matrix(0,
                                 ncol = length(fe_nm),
                                 nrow = nrow(asb_sp_w)))
  colnames(asb_fe) <- fe_nm
  rownames(asb_fe) <- rownames(asb_sp_w)
  
  
  # Fill the asb_fe df by looping on species, linking sp -> FE:
  for (i in names(sp_fe)) {
    
    # get the FE associated with the i species:
    assoc_FE <- sp_fe[which(names(sp_fe) == i)][[1]]
    
    # get the presabs values of the i species:
    values_sp <- asb_sp_w[, which(colnames(asb_sp_w) == i)]
    
    # fill the asb_fe df:
    asb_fe[, assoc_FE] <- asb_fe[, assoc_FE] + values_sp
    
  }
  
  # Build the asb_fe df that only contains occurrence of the FE in each ...
  # ... assemblage:
  asb_fe_occ <- asb_fe
  asb_fe_occ[asb_fe_occ > 1] <- 1
  
  return("asb_fe_occ" = asb_fe_occ)
  
}


#' Get a data frame linking Functional Entities names and species names
#'
#' @param sp_to_fe a list which is the output of the \code{mFD::sp.to.fe()} 
#' function 
#'
#' @return a data frame linking FEs names to species names with columns being
#' FE names and Species names.
#' 
#' @author Camille Magneville
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
#'  # Load Assemblages*Species matrix:
#' data('baskets_fruits_weights', package = 'mFD')
#'
#' # Gather species into FEs:
#' ## gathering species into FEs (FEs named according to the decreasing...
#' ## ...  number of species they gather):
#'  sp_FEs_fruits <- mFD::sp.to.fe(
#'       sp_tr      = fruits_traits, 
#'       tr_cat     = fruits_traits_cat, 
#'       fe_nm_type = "fe_rank")
#'       
#' # Create the data frame containing species and FEs names:
#' fe_sp_df <- mFD::fe.sp.df.computation(sp_to_fe = sp_FEs_fruits)
#' 
#' @export
#'

fe.sp.df.computation <- function(sp_to_fe) {
  
  
  # Create a dataframe that contain fe_nm and sp_nm:
  fe_sp_df <- as.data.frame(matrix(ncol = 2, nrow = 1))
  colnames(fe_sp_df) <- c("fe_nm", "species_nm")
  
  # loop on fes:
  for (i in (1:length(sp_to_fe$fe_nm))) {
    
    # Get the names of species belonging to the studied FE:
    sp_nm <- mFD::search.sp.nm(sp_to_fe, paste0("fe", sep = "_", i))
    
    # loop on species:
    for (j in (1:length(sp_nm))) {
      
      fe_sp_df[nrow(fe_sp_df) + 1, "fe_nm"] <- paste0("fe", sep = "_", i)
      fe_sp_df[nrow(fe_sp_df), "species_nm"] <- sp_nm[j]
      
    } # end loop on species of the studied FE
    
  } # end loop on FEs
  
  final_fe_sp_df <- fe_sp_df[-1, ]
  
  return(final_fe_sp_df)
  
}


#' Get the names of species belonging to a specific Functional Entity (FE)
#'
#' @param sp_to_fe a list which is the output of the mFD::sp.to.fe() function 
#' 
#' @param nm_fe a character string referring to the name of the FE to study
#'
#' @return a vector containing the names of the species belonging to the studied
#' FE.
#' 
#' @author Camille Magneville
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
#'  # Load Assemblages*Species matrix:
#' data('baskets_fruits_weights', package = 'mFD')
#'
#' # Gather species into FEs:
#' ## gathering species into FEs (FEs named according to the decreasing...
#' ## ...  number of species they gather):
#'  sp_FEs_fruits <- mFD::sp.to.fe(
#'       sp_tr      = fruits_traits, 
#'       tr_cat     = fruits_traits_cat, 
#'       fe_nm_type = "fe_rank")
#'       
#' # Look for te names of the species belonging to the FE called "fe_3":
#' mFD::search.sp.nm(sp_to_fe = sp_FEs_fruits,
#'                  nm_fe = "fe_3")
#' 
#' @export
#'


search.sp.nm <- function(sp_to_fe, nm_fe) {
  
  
  # Check that the fe studied is belonging to the list of existing FEs:
  if (! nm_fe %in% sp_to_fe$fe_nm){
    stop("This Functional Entity doesn't exist. Please check its spelling.")
  }
  
  sp_vect <- c()
  
  for (i in names(sp_to_fe$sp_fe)) {
    
    if (sp_to_fe$sp_fe[[i]] == nm_fe) {
      
      sp_vect <- append(sp_vect, i)
      
    }
    
  }
  
  return(sp_vect)
  
}


#' Convert the data frame of FEs coordinates to a species coordinates one
#'
#' @param fe_faxes_coord the data frame gathering FEs coordinates along
#' functional axes, with columns representing axes and rows representing FEs.
#' 
#' @param asb_sp_w the data frame gathering species*assemblages, with species
#' in columns and assemblages in rows.
#' 
#' @param sp_to_fe the output of the mFD::sp.to.fe() function.
#'
#' @return a data frame gathering species coordinates, with axes being columns
#' and species being rows. \strong{Note: Works for up to 10 PC}.
#' 
#' @author Camille Magneville
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
#'  fruits_traits_cat <- fruits_traits_cat[c(1:4), -3]
#'  
#'  # Load Assemblages*Species matrix:
#' data('baskets_fruits_weights', package = 'mFD')
#'
#' # Gather species into FEs:
#' ## gathering species into FEs (FEs named according to the decreasing...
#' ## ...  number of species they gather):
#'  sp_FEs_fruits <- mFD::sp.to.fe(
#'       sp_tr      = fruits_traits, 
#'       tr_cat     = fruits_traits_cat, 
#'       fe_nm_type = "fe_rank")
#'       
#' # Compute the functional distance between FEs:
#' fe_dist_fruits <- mFD::funct.dist(sp_tr         = sp_FEs_fruits$fe_tr,
#'                                   tr_cat        = fruits_traits_cat,
#'                                   metric        = "gower",
#'                                   scale_euclid  = "scale_center",
#'                                   ordinal_var   = "classic",
#'                                   weight_type   = "equal",
#'                                   stop_if_NA    = TRUE)
#'  
#'  # Compute the quality of functional spaces:
#'  fspaces_quality_fruits <- mFD::quality.fspaces(
#'                                           sp_dist = fe_dist_fruits,
#'                                           maxdim_pcoa = 10,
#'                                           deviation_weighting = "absolute",
#'                                           fdist_scaling       = FALSE,
#'                                           fdendro             = "average")
#'  fspaces_quality_fruits$quality_fspaces
#'  
#'  # Get coordinates of FEs in the 3D space:
#'  fe_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
#'  fe_faxes_coord_fruits <- fe_faxes_coord_fruits[, 1:3]
#'       
#' # Get the matrix of species coordinates, from the matrix of FEs coordinates:
#' sp_faxes_coord <- mFD::from.fecoord.to.spcoord(
#'                                       fe_faxes_coord = fe_faxes_coord_fruits,
#'                                       asb_sp_w = baskets_fruits_weights,
#'                                       sp_to_fe = sp_FEs_fruits)
#' 
#' @export
#'

from.fecoord.to.spcoord <- function(fe_faxes_coord, 
                                    asb_sp_w,
                                    sp_to_fe) {
  
  
  # Create a data frame with as many columns as there are PCs and ...
  # ... as many rows as they are species:
  sp_faxes_coord <- as.data.frame(matrix(0, ncol = ncol(fe_faxes_coord),
                                         nrow = ncol(asb_sp_w)))
  colnames(sp_faxes_coord) <- colnames(fe_faxes_coord)
  rownames(sp_faxes_coord) <- colnames(asb_sp_w)
  
  # Loop on each FE:
  for (i in (1:nrow(fe_faxes_coord))) {
    
    # Get the names of the species belonging to this FE:
    sp_nm <- mFD::search.sp.nm(sp_to_fe, rownames(fe_faxes_coord)[i])
    
    values <- fe_faxes_coord[i, ]
    
    # Fill the sp faxes coord:
    for (j in c(1:ncol(fe_faxes_coord))) {
      sp_faxes_coord[which(rownames(sp_faxes_coord) %in% sp_nm), j] <- values[j]
    }

  }
  
  return(sp_faxes_coord)
  
}
