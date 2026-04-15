# Compute a set of alpha functional indices for a set of assemblages

This function computes a set of multidimensional space based indices of
alpha functional diversity. The user can choose which functional indices
to compute.

## Usage

``` r
alpha.fd.multidim(
  sp_faxes_coord,
  asb_sp_w,
  ind_vect = c("fide", "fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", "fspe"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE,
  verbose = TRUE
)
```

## Arguments

- sp_faxes_coord:

  a matrix of species coordinates in a chosen functional space. Species
  coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- asb_sp_w:

  a matrix linking weight of species (columns) and a set of assemblages
  (rows).

- ind_vect:

  a vector of character string of the name of functional indices to
  compute. **Indices names must be written in lower case letters**.
  Possible indices to compute are: 'fide', fdis', 'fmpd', 'fnnd',
  'feve', 'fric', 'fdiv', 'fori' and 'fspe'. Default: all the indices
  are computed.

- scaling:

  a logical value indicating if scaling is to be done (TRUE) or not
  (FALSE) on functional indices. Scaling is used to be able to compare
  indices values between assemblages. Default: scaling = TRUE.

- check_input:

  a logical value indicating whether key features the inputs are checked
  (e.g. class and/or mode of objects, names of rows and/or columns,
  missing values). If an error is detected, a detailed message is
  returned. Default: `check.input = TRUE`.

- details_returned:

  a logical value indicating whether the user want to store details.
  Details are used in graphical functions and thus must be kept if the
  user want to have graphical outputs for the computed indices.

- verbose:

  a logical value indicating whether progress details should be printed
  in the console. If `FALSE` does not provide percent progress when
  computing diversity indices.

## Value

The following list is returned:

- *functional_diversity_indices* matrix containing indices values
  (columns) for each assemblage (rows)

- *details* list: a **asb_sp_occ** data.frame of species occurrences in
  each assemblage ; a **asb_sp_relatw** matrix of relative weight of
  species in each assemblage ; a **sp_coord_all_asb** list of matrices
  of species coordinates along functional axes for species present in
  each assemblage ; a **vert_nm_all_asb** list of vectors of species
  names being vertices of the convex hull for each assemblage ; a
  **mst_all_asb** list of data.frames summarizing link between species
  in the minimum spanning tree of each assemblage ; a
  **grav_center_vert_coord_all_asb** list of vectors of coordinates of
  the vertices gravity center for each assemblage ; a
  **mean_dtogravcenter_all_asb** list of vectors containing mean
  distance to the species gravity center for each assemblage ; a
  **dist_gravcenter_global_pool** vector containing the distance of each
  species to the gravity center of all species from the global pool ; a
  **dist_nn_global_pool** data.frame showing the distances of each
  species from the global pool to its nearest neighbor ; a
  **nm_nn_all_asb** data.frame containing the name of each nearest
  neighbor of each species present in a given assemblage ; a
  **dist_nn_all_asb** data.frame containing distance of each species
  present in a given assemblage to its nearest neighbor.

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
# Load Species*Traits dataframe:
data('fruits_traits', package = 'mFD')

# Load Assemblages*Species dataframe:      
data('baskets_fruits_weights', package = 'mFD')

# Load Traits categories dataframe:
data('fruits_traits_cat', package = 'mFD')  
 
# Compute functional distance 
sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
                                  tr_cat        = fruits_traits_cat,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)
#> [1] "Running w.type=equal on groups=c(Size)"
#> [1] "Running w.type=equal on groups=c(Plant)"
#> [1] "Running w.type=equal on groups=c(Climate)"
#> [1] "Running w.type=equal on groups=c(Seed)"
#> [1] "Running w.type=equal on groups=c(Sugar)"
#> [1] "Running w.type=equal on groups=c(Use,Use,Use)"
 
# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fruits <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fruits, 
  maxdim_pcoa         = 10,
  deviation_weighting = 'absolute',
  fdist_scaling       = FALSE,
  fdendro             = 'average')
#> Registered S3 method overwritten by 'dendextend':
#>   method     from 
#>   rev.hclust vegan
  
# Retrieve species coordinates matrix:
sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord

# Compute alpha diversity indices
alpha_fd_indices_fruits <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fruits[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = baskets_fruits_weights, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)
#> basket_1 done 10%
#> basket_2 done 20%
#> basket_3 done 30%
#> basket_4 done 40%
#> basket_5 done 50%
#> basket_6 done 60%
#> basket_7 done 70%
#> basket_8 done 80%
#> basket_9 done 90%
#> basket_10 done 100%
  
# Retrieve alpha diversity indices table
fd_ind_values_fruits <- alpha_fd_indices_fruits$functional_diversity_indices
fd_ind_values_fruits
#>           sp_richn      fdis      fmpd      fnnd     feve        fric      fdiv
#> basket_1         8 0.4572320 0.6366416 0.5778912 0.648241 0.127034350 0.5385777
#> basket_2         8 0.6797564 0.7244031 0.8106286 0.797667 0.127034350 0.7974857
#> basket_3         8 0.7002352 0.7308197 0.8261547 0.791014 0.127034350 0.7995983
#> basket_4         8 0.2854746 0.3351787 0.3259376 0.683059 0.004136639 0.6176418
#> basket_5         8 0.3213250 0.3509955 0.3546202 0.830008 0.004136639 0.6875893
#> basket_6         8 0.7577626 0.7829217 0.8741800 0.779539 0.110378639 0.8866375
#> basket_7         8 0.7907944 0.8111233 0.8785200 0.793622 0.110378639 0.8937828
#> basket_8         8 0.4190958 0.5149757 0.4450844 0.586306 0.014059458 0.6059633
#> basket_9         8 0.5107095 0.5506783 0.5371850 0.787953 0.014059458 0.6705378
#> basket_10        8 0.4825409 0.5352984 0.5167429 0.783877 0.031100594 0.7230460
#>                fori      fspe    fide_PC1    fide_PC2     fide_PC3     fide_PC4
#> basket_1  0.3789936 0.3931170  0.01899853  0.02941357 -0.005068006 -0.018344605
#> basket_2  0.5320906 0.5687316 -0.01331858  0.05692236 -0.023750910 -0.003179800
#> basket_3  0.5447687 0.5790966 -0.03243182  0.05008021 -0.026841109  0.006822367
#> basket_4  0.2893770 0.2604065 -0.01050901 -0.01509452 -0.024893697 -0.048687531
#> basket_5  0.3035437 0.3019012 -0.01540439 -0.01115358 -0.006787842 -0.067223460
#> basket_6  0.7708788 0.7829723 -0.21168832  0.07159069 -0.056906946  0.037226970
#> basket_7  0.7237701 0.7451308 -0.11454663  0.14176088 -0.058487715  0.029260412
#> basket_8  0.3324008 0.6220539  0.20828315 -0.01234992  0.047915332  0.039808906
#> basket_9  0.3759408 0.5785916  0.14913920 -0.01500900  0.076594090  0.010513208
#> basket_10 0.3934059 0.4219834  0.02616413  0.01620247  0.002831861 -0.068301807
```
