# Get a data frame linking Functional Entities names and species names

Get a data frame linking Functional Entities names and species names

## Usage

``` r
fe.sp.df.computation(sp_to_fe)
```

## Arguments

- sp_to_fe:

  a list which is the output of the
  [`mFD::sp.to.fe()`](https://cmlmagneville.github.io/mFD/reference/sp.to.fe.md)
  function

## Value

a data frame linking FEs names to species names with columns being FE
names and Species names.

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
      
# Create the data frame containing species and FEs names:
fe_sp_df <- mFD::fe.sp.df.computation(sp_to_fe = sp_FEs_fruits)
```
