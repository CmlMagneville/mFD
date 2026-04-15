# Compute distance of a given point to its nearest neighbor in the functional space and the identity of the nearest neighbor

This function is used in functional indices computation.

## Usage

``` r
dist.nearneighb(sp_faxes_coord, ref_sp)
```

## Arguments

- sp_faxes_coord:

  a matrix of species coordinates in a chosen functional space. Species
  coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- ref_sp:

  a character string referring to the name of the reference species.

## Value

A list containing the nearest neighbor identity `nn_id` and a list of
the distance of the reference point to its nearest neighbor
`nn_ref_sp_dist`.

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
# Load Species*Traits dataframe:
 data("fruits_traits", package = "mFD")

# Load Assemblages*Species dataframe:      
 data("baskets_fruits_weights", package = "mFD") 

# Load Traits categories dataframe:
 data("fruits_traits_cat", package = "mFD") 
 
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
                                  deviation_weighting = "absolute",
                                  fdist_scaling       = FALSE,
                                  fdendro             = "average")
 
# Retrieve species coordinates matrix:
 sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord

# Compute the distance of "pear" to its nearest neighbor(s):
 dist_nn_pear <- dist.nearneighb(sp_faxes_coord_fruits, ref_sp = "pear")
 dist_nn_pear
#> $`nearest neighbour identity`
#> [1] "apple"
#> 
#> $`distance of the reference species to its nearest neighbour`
#> [1] 0.1270335
#> 
```
