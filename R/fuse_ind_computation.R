# Functions to compute FUSE Functionally Unique, Specialized and Endangered
# ... indice and related functions
#
# Authors: Fabien Leprieur & Camille Albouy
#
#------------------------------------------------------------------------------
#
#
#' Compute FUSE (Functionally Unique, Specialized and Endangered)
#'
#' @param sp_dist a dist object provided by the daisy() of the cluster
#' package or dist.ktab() of the ade4 package
#'
#' @param sp_faxes_coord a data.frame with the coordinates of the species on a
#'   multidimensional space based on a selected number of axes derived from a
#'   Principal Coordinate Analysis (PCOA). The species are in rows and the PCOA
#'   axes are in column.
#'
#' @param nb_NN a numerical value giving the number of nearest neighbor to
#'   consider. Default: nb_NN = 5
#'
#' @param GE a numerical vector giving the IUCN status rank (DD = NA, LC = 0,
#'   NT = 1, VU = 2, EN = 3, CR = 4) or the IUCN extinction probability associated with
#'   each status, see Mooers et al.(2008) for example with DD = NA, LC = 0, NT = 0.1,
#'   VU = 0.4, EN = 0.666, CR = 0.999)
#'
#' @param standGE boolean value to standardize the GE values (TRUE or FALSE)
#'
#' @return a dataframe with species in rows and the different metrics in
#'   columns. The metrics are: \itemize{
#'   \item FUSE: functionally unique, specialized and endangered (see Pimiento
#'   et al. 2020 and Griffin et. al 2020)
#'   \item FDGE: functionally distinctive and endangered
#'   \item FUn_std: functional uniqueness standardized between 0 and 1 (see
#'   Mouillot et al. 2013 and Griffin et al. 2020)
#'   \item FSp_std:functional specialization standardized between 0 and 1  (see
#'   Mouillot et al. 2013 and Griffin et al. 2020)
#'   \item FDist_std:functional distinctiveness standardized between 0 and 1
#'   (see Violle et al. 2007 and Griffin et al. 2020)
#'   }
#'
#'@examples
#' sp_tr <- read.csv(system.file("extdata", "data_traits_MMA_ursus.csv", package = "mFD"),
#' dec=",", sep=";", header=TRUE,row.names=1, na.strings="NA")
#' # Trait compilation and ordination:
#'  dimorphism <- ordered(sp_tr$dimorphism)
#'  breeding_site   <- ordered(sp_tr$breeding_site)
#'  social_behavior <- ordered(sp_tr$social_behavior)
#'  weight_max <- log(sp_tr$adult_weight_max)
#'  social_group <- log(sp_tr$social_group_mean)
#' # Trait Matrix construction:
#'  sp_tr_end <- data.frame(main_diet = sp_tr$main_diet, 
#'   foraging_water_depth = sp_tr$foraging_water_depth,
#'   foraging_location = sp_tr$foraging_location, 
#'   fasting_strategy = sp_tr$fasting_strategy,
#'   female_sexual_maturity = sp_tr$female_sexual_maturity, 
#'   weaning = sp_tr$weaning,
#'   gestation = sp_tr$gestation, inter_litter = sp_tr$inter_litter,
#'   breeding_site = sp_tr$breeding_site, social_group = sp_tr$social_group_mean,
#'   social_behavior = sp_tr$social_behavior, weight_max = sp_tr$adult_weight_max,
#'   dimorphism = sp_tr$dimorphism)
#'  rownames(sp_tr_end) <- rownames(sp_tr)
#' # Function weigthing vector:
#'  v <- c(0.25, 0.25, 0.25, 0.25, 0.20, 0.20, 0.20, 0.20, 0.20, 0.5, 0.5, 0.5, 0.5)
#' # Gower distance calculation:
#'  sp_tr_end$main_diet <- as.factor(sp_tr_end$main_diet)
#'  sp_tr_end$foraging_water_depth <- as.factor(sp_tr_end$foraging_water_depth)
#'  sp_tr_end$foraging_location <- as.factor(sp_tr_end$foraging_location)
#'  sp_tr_end$breeding_site <- as.factor(sp_tr_end$breeding_site)
#'  sp_tr_end$social_behavior <- as.factor(sp_tr_end$social_behavior)
#'  sp_dist_tr <- cluster::daisy(sp_tr_end, metric = c("gower"), 
#'   type = list(symm = c(4)), weights = v)
#' # Principal coordinate analyses
#'  Pcoa <- ade4::dudi.pco(ade4::quasieuclid(sp_dist_tr), scann = FALSE, nf = 40)
#'  sp_faxes_coord <- Pcoa$li[1:40]
#'# FUSE calculation:
#'  FUSE_res <- mFD::fuse(sp_dist = sp_dist_tr, sp_faxes_coord = as.matrix(sp_faxes_coord), 
#'   nb_NN = 5, GE = sp_tr$IUCN_num,
#'   standGE = TRUE)
#'  FUSE_res2 <- mFD::fuse(sp_dist = sp_dist_tr, sp_faxes_coord = as.matrix(sp_faxes_coord), 
#'   nb_NN = 5,GE = sp_tr$IUCN_50,
#'   standGE = TRUE)
#'  FUSE_res3 <- mFD::fuse(sp_dist = sp_dist_tr, sp_faxes_coord = as.matrix(sp_faxes_coord), 
#'   nb_NN = 5, GE = sp_tr$IUCN_100,
#'  standGE = TRUE)
#'
#'@export


fuse <- function(sp_dist, sp_faxes_coord, nb_NN = 5, GE, standGE = F) {
  
  if (! identical(row.names(as.matrix(sp_dist)), row.names(sp_faxes_coord))) {
    stop("Coords lines do not match with the distance matrix")
  }
  
  if (!is.matrix(sp_faxes_coord)) {
    stop("Error: species coordinates on functional axes should be provided as
      a matrix. Please check.")
  }
  
  if (any(is.na(sp_faxes_coord))) {
    stop("Error: The species*coordinates matrix contains NA. Please check.")
  }
  
  if (is.null(rownames(sp_faxes_coord))) {
    stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
  }
  
  if (is.null(rownames(sp_faxes_coord))) {
    stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
  }
  
  if (any(is.na(sp_dist))) {
    stop("Error: The species distances matrix contains NA. Please check.")
  }
  if (is.null(rownames(as.matrix(sp_dist)))) {
    stop("Error: No row names provided in species distances matrix.
             Please add species names as row names.")
  }
  
  nm <- rownames(sp_faxes_coord)
  
  # Specialization calculation
  O <- apply(sp_faxes_coord, 2, mean)
  spe <- apply(sp_faxes_coord, 1, function(x){
    sum((x - O)^2)^0.5})
  
  # Distinctivness calculation
  dist_sp <- as.matrix(sp_dist)
  Fdistinct <- apply(dist_sp, 1, mean)
  
  # Uniqueness calculation
  uni_res <- get_indicator(sp_dist = as.matrix(sp_dist), nb_NN = nb_NN)
  uniqu <- uni_res$Average_uniqueness[, "Mean"]
  
  if (standGE == TRUE) {
    GE <- as.vector(vegan::decostand(GE, "range", na.rm = TRUE))
  }
  
  # FUSE metrics calculation
  FUn_std <- as.vector(vegan::decostand(uniqu, "range"))
  FUGE <- log(1 + (FUn_std * GE))
  FSp_std <- as.vector(vegan::decostand(spe, "range"))
  FSGE <- log(1 + (FSp_std * GE))
  FUS <- FUn_std + FSp_std
  FDist_std <- as.vector(vegan::decostand(Fdistinct, "range"))
  FUSE <- stats::setNames(FUGE + FSGE, nm = nm)
  FUSE_alt <- log(1 + FUS * GE)
  FDGE <- log(1 + (FDist_std * GE))
  
  data.frame(cbind(FUSE, FUSE_alt, FUGE, FDGE, FSGE, FUn_std, FSp_std, FDist_std))
  
}

#------------------------------------------------------------------------------
#
#
#' Compute the nearest neighbors for all considered species for FUSE computation
#'
#' @param sp_dist  a matrix object representing the distance in an euclidean
#'   space between species based on species traits
#'
#' @param nb_NN a numerical value giving the number of nearest neighbor to
#'   consider
#'


get_indicator <- function(sp_dist, nb_NN){
  
  w <- reshape2::melt(sp_dist)
  s <- split(w, f = w[, 2])
  
  Res <- lapply(s, function(x) {
    get_dist_func(nb_NN = nb_NN, data = x)})
  Res_mean_sd <- do.call(rbind, lapply(1:length(Res), function(i){
    Res[[i]][[1]]}))
  NN <- lapply(1:length(Res), function(i){Res[[i]][[2]]})
  
  rownames(Res_mean_sd) <- names(NN) <- names(Res)
  list(Average_uniqueness = Res_mean_sd, Nearest_neighbour = NN)
  
}


#------------------------------------------------------------------------------
#
#
#' Calculate the nearest neighbors for a given species
#'
#'@param nb_NN a numerical value giving the number of nearest neighbor to consider
#'
#'@param data a dataframe considering all the distance between the considered
#'  species and all of these neighbors
#'

get_dist_func <- function(nb_NN, data){
  
  data <- data[order(data[, 3], decreasing = F), ]
  data <- data[-1, ]
  mm <- mean(data[1:nb_NN, 3])
  sd <- sd(data[1:nb_NN, 3])
  sp <- as.character(data[1:nb_NN, 1])
  list(c(Mean = mm, Sd = sd), Species = sp)
  
}