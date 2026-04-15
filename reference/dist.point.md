# Compute distances of all points to a given point in the functional space

This function computes the distances of all species to a reference
species. It is used in FSpe, FOri and FNND computation.

## Usage

``` r
dist.point(sp_faxes_coord, ref_sp)
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

A vector of species distances to the reference species.

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

# Retrieve the distances of all species to "pear":
 dist_pear <- dist.point(sp_faxes_coord_fruits, ref_sp = "pear")
 dist_pear
#>         apple       apricot        banana       currant    blackberry 
#>     0.1270335     0.2810337     0.3857251     0.3700010     0.3272434 
#>     blueberry        cherry         grape    grapefruit     kiwifruit 
#>     0.4182818     0.3786579     0.3233137     0.2405707     0.1688475 
#>         lemon          lime        litchi         mango         melon 
#>     0.2491467     0.3896421     0.4622051     0.4371372     0.2904542 
#>        orange passion_fruit         peach          pear     pineapple 
#>     0.1809690     0.4018250     0.2216462     0.0000000     0.4871492 
#>          plum     raspberry    strawberry     tangerine   water_melon 
#>     0.2515970     0.2924277     0.3746643     0.1998687     0.3290931 
```
