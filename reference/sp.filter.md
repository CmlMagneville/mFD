# Retrieve information about species in a given assemblage

This function computes names of species present in an given assemblage,
their coordinates in the functional space and their weights. It is used
in the `alpha_FD_multidim` function to filter species and compute each
functional indices for each community.

## Usage

``` r
sp.filter(asb_nm, sp_faxes_coord, asb_sp_w)
```

## Arguments

- asb_nm:

  a string object referring to the name of a given community.

- sp_faxes_coord:

  a matrix of species coordinates in a chosen functional space. Species
  coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- asb_sp_w:

  a matrix linking weight of species (columns) and a set of assemblages
  (rows).

## Value

A vector containing names of species present in a given assemblage
`sp_name_asb_k`, a matrix containing coordinates of species present in a
given assemblage `sp_faxes_coord_k`, a matrix containing weight of
species present in a given assemblage `asb_sp_w_k`, a matrix containing
relative weight of species present in a given assemblage
`asb_sp_relatw_k`.

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
sp_dist_fruits <- mFD::funct.dist(
 sp_tr         = fruits_traits,
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

# Filter species of basket_1 assemblage:
sp.filter(asb_nm         = "basket_1", 
          sp_faxes_coord = sp_faxes_coord_fruits, 
          asb_sp_w       = baskets_fruits_weights)
#> $`species names`
#> [1] "apple"         "banana"        "cherry"        "lemon"        
#> [5] "melon"         "passion_fruit" "pear"          "strawberry"   
#> 
#> $`species coordinates`
#>                         PC1         PC2          PC3         PC4         PC5
#> apple          0.0265193144  0.01305413 -0.043717435 -0.02121146 -0.14821854
#> banana        -0.3346937349  0.09137141 -0.070877547  0.06988872 -0.05000085
#> cherry         0.1532479527 -0.25099690 -0.135570593  0.05621669  0.02650788
#> lemon         -0.1116123197  0.05271882  0.138487467 -0.10319375  0.02307706
#> melon          0.0475710274  0.22821444 -0.009073564 -0.13015462  0.03728241
#> passion_fruit -0.1743734750 -0.09195180  0.051587716  0.22053201 -0.14745293
#> pear          -0.0008010876  0.01968943 -0.013204612 -0.07056511 -0.07014583
#> strawberry     0.2743910893  0.09325097  0.046622086  0.09331814  0.04954010
#>                        PC6          PC7         PC8          PC9         PC10
#> apple          0.040965006  0.039965979 -0.02257727 -0.012017348  0.010675722
#> banana         0.020752971 -0.041476147 -0.08303923  0.042337245 -0.015830902
#> cherry        -0.020780891  0.004803343 -0.04262059  0.030307476  0.022668724
#> lemon          0.101128620  0.044498483 -0.03408272 -0.023135838  0.021535320
#> melon         -0.136045993  0.036709357  0.01345963 -0.022081148 -0.005067072
#> passion_fruit  0.066901090 -0.006120818  0.06645768  0.020143504  0.008008755
#> pear          -0.004973365  0.023163490 -0.05589114 -0.008962102 -0.038038792
#> strawberry     0.001191018  0.094872368  0.01891473  0.012390681  0.018860955
#> 
#> $`species weight`
#>          apple banana cherry lemon melon passion_fruit pear strawberry
#> basket_1   400    100    150   200   200           100  600        250
#> 
#> $`species relative weight`
#>          apple banana cherry lemon melon passion_fruit pear strawberry
#> basket_1   0.2   0.05  0.075   0.1   0.1          0.05  0.3      0.125
#> 
```
