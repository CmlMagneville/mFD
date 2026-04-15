# Build the assemblage-FEs dataframe from the assemblages-species one

This function computes an occurrence data frame with assemblages in rows
and Functional Entities (FEs) in columns, from the dataframe of
assemblage-species and the output of the
[`mFD::sp.to.fe()`](https://cmlmagneville.github.io/mFD/reference/sp.to.fe.md)
function.

## Usage

``` r
from.spfe.to.feasb(sp_fe, asb_sp_w)
```

## Arguments

- sp_fe:

  list gathering to which FE belongs eac species. It is the output of
  the mFD::sp.to.fe() function - \$sp_fe.

- asb_sp_w:

  the assemblage \* species data frame with assemblages being rows and
  species being columns

## Value

an occurrence dataframe with studied FEs in columns and assemblages in
rows. **Be careful: It's an occurrence data frame, not an abundance
one**

## Author

Camille Magneville

## Examples

``` r
# Load species traits data:
 data("fruits_traits", package = "mFD")

# Transform species traits data:
# Only keep the first 4 traits to illustrate FEs:
 fruits_traits <- fruits_traits[ , c(1:4)]   

# Load trait types data:
 data("fruits_traits_cat", package = "mFD")

# Transform the trait types data to only keep traits 1 - 4:
 fruits_traits_cat <- fruits_traits_cat[c(1:4), ]
 
 # Load Assemblages*Species matrix:
data('baskets_fruits_weights', package = 'mFD')

# Gather species into FEs:
## gathering species into FEs (FEs named according to the decreasing...
## ...  number of species they gather):
 sp_FEs_fruits <- mFD::sp.to.fe(
      sp_tr      = fruits_traits, 
      tr_cat     = fruits_traits_cat, 
      fe_nm_type = "fe_rank")
      
# Get the list which gather to which FE belongs each species:
sp_fes_list <- sp_FEs_fruits$sp_fe
      
# Build the Assemblages*FEs data frame:
asb_fes <- mFD::from.spfe.to.feasb(sp_fe = sp_fes_list,
                                   asb_sp_w = baskets_fruits_weights)
asb_fes
#>           fe_1 fe_2 fe_6 fe_7 fe_3 fe_8 fe_9 fe_10 fe_11 fe_12 fe_4 fe_13 fe_14
#> basket_1     1    0    1    0    0    0    1     0     0     0    1     0     0
#> basket_2     1    0    1    0    0    0    1     0     0     0    1     0     0
#> basket_3     1    0    1    0    0    0    1     0     0     0    1     0     0
#> basket_4     1    1    0    0    0    0    0     0     0     1    1     0     0
#> basket_5     1    1    0    0    0    0    0     0     0     1    1     0     0
#> basket_6     1    0    1    0    0    0    0     0     0     0    1     1     1
#> basket_7     1    0    1    0    0    0    0     0     0     0    1     1     1
#> basket_8     0    0    0    1    1    1    1     1     0     0    1     0     0
#> basket_9     0    0    0    1    1    1    1     1     0     0    1     0     0
#> basket_10    1    1    0    0    0    0    0     1     1     0    0     0     0
#>           fe_15 fe_5 fe_16 fe_17 fe_18 fe_19 fe_20
#> basket_1      0    1     1     0     0     1     0
#> basket_2      0    1     1     0     0     1     0
#> basket_3      0    1     1     0     0     1     0
#> basket_4      0    0     0     1     0     0     1
#> basket_5      0    0     0     1     0     0     1
#> basket_6      1    1     0     0     1     0     0
#> basket_7      1    1     0     0     1     0     0
#> basket_8      0    0     0     0     0     1     0
#> basket_9      0    0     0     0     0     1     0
#> basket_10     0    1     0     0     0     1     0

```
