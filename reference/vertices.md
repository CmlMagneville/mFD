# Compute vertices of the Minimal Convex Hull shaping species from a single assemblage in a multidimensional functional space

This function identifies species that are vertices of the minimal convex
hull enclosing a community in a multidimensional functional space. This
function is using the
[`convhulln`](https://rdrr.io/pkg/geometry/man/convhulln.html) function.

## Usage

``` r
vertices(sp_faxes_coord, order_2D = FALSE, check_input = FALSE)
```

## Arguments

- sp_faxes_coord:

  a matrix of species coordinates in a chosen functional space. Species
  coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- order_2D:

  a logical value defining whether vertices names are reordered so that
  they define a convex polygon in 2D which is convenient for plotting.
  Default is `FALSE`, vertices ordered as in row names of
  'sp_faxes_coord'.

- check_input:

  a logical value defining whether inputs are checked before
  computation: species names must be put as row.names, there must be no
  NA and species number must be superior to (axes number + 1). Default:
  `check_input = TRUE`.

## Value

A vector containing names of species being vertices `vert_nm`.

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

# Compute vertices and order them clockwise:
 vert_nm <- vertices(sp_faxes_coord_fruits[ , c("PC1", "PC2")], 
  order_2D = TRUE, check_input = TRUE)
 vert_nm
#> [1] "cherry"      "currant"     "blueberry"   "water_melon" "pineapple"  
#> [6] "banana"      "mango"       "litchi"     
```
