# Compute the set of indices based on number of species in Functional Entities

This function computes the set of indices based on number of species in
Functional Entities (FEs) following Mouillot *et al.* (2014).

## Usage

``` r
alpha.fd.fe(
  asb_sp_occ,
  sp_to_fe,
  ind_nm = c("fred", "fored", "fvuln"),
  check_input = TRUE,
  details_returned = TRUE
)
```

## Arguments

- asb_sp_occ:

  a matrix linking occurrences (coded as 0/1) of species (columns) in a
  set of assemblages (rows). Warning: **An assemblage must contain at
  least one species**.

- sp_to_fe:

  a list with details of species clustering into FE from
  [`sp.to.fe`](https://cmlmagneville.github.io/mFD/reference/sp.to.fe.md).

- ind_nm:

  a vector of character strings with the names of functional diversity
  indices to compute among 'fred', 'fored' and 'fvuln'. **Indices names
  must be written in lower case letters**. Default: all the indices are
  computed.

- check_input:

  a logical value indicating whether key features the inputs are checked
  (e.g. class and/or mode of objects, names of rows and/or columns,
  missing values). If an error is detected, a detailed message is
  returned. Default: `check.input = TRUE`.

- details_returned:

  a logical value indicating whether details about indices computation
  should be returned. These details are required by
  [`alpha.fd.fe.plot`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.fe.plot.md)
  to plot FEs indices.

## Value

A list with:

- *asb_fdfe* a matrix containing for each assemblage (rows), values of
  functional diversity indices (same names than in 'ind_nm') as well as
  the number of species ('nb_sp') and the number of FE (nb_fe);

- if *details_returned* is `TRUE`,

- *details_fdfe* a list with *asb_fe_nbsp* a matrix with number of
  species per FE in each assemblage.

## References

Mouillot *et al.* (2014) Functional over-redundancy and high functional
vulnerability in global fish faunas on tropical reefs. *PNAS*, **38**,
13757-13762.

## Author

Camille Magneville

## Examples

``` r
# Load Species*Traits dataframe:
data('fruits_traits', package = 'mFD')

# Load Traits categories dataframe:
data('fruits_traits_cat', package = 'mFD')

# Load Assemblages*Species matrix:
data('baskets_fruits_weights', package = 'mFD')

# Remove continuous trait:
fruits_traits <- fruits_traits[, -5]
fruits_traits_cat <- fruits_traits_cat[-5, ]

# Compute gathering species into FEs:
sp_to_fe_fruits <- mFD::sp.to.fe(sp_tr = fruits_traits, 
 tr_cat = fruits_traits_cat, 
 fe_nm_type = 'fe_rank', check_input = TRUE)
#> Warning: All Functional Entities have a single species.
 
# Get the occurrence dataframe:
asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights) 
asb_sp_fruits_occ <- asb_sp_fruits_summ$'asb_sp_occ'

# Compute alpha fd indices:
alpha.fd.fe(
   asb_sp_occ       = asb_sp_fruits_occ, 
   sp_to_fe         = sp_to_fe_fruits,
   ind_nm           = c('fred', 'fored', 'fvuln'),
   check_input      = TRUE, 
   details_returned = TRUE)
#> $asb_fdfe
#>           nb_sp nb_fe fred fored fvuln
#> basket_1      8     8    1     0     1
#> basket_2      8     8    1     0     1
#> basket_3      8     8    1     0     1
#> basket_4      8     8    1     0     1
#> basket_5      8     8    1     0     1
#> basket_6      8     8    1     0     1
#> basket_7      8     8    1     0     1
#> basket_8      8     8    1     0     1
#> basket_9      8     8    1     0     1
#> basket_10     8     8    1     0     1
#> 
#> $details_fdfe
#> $details_fdfe$asb_fe_nbsp
#>           fe_1 fe_2 fe_3 fe_4 fe_5 fe_6 fe_7 fe_8 fe_9 fe_10 fe_11 fe_12 fe_13
#> basket_1     1    0    1    0    0    0    1    0    0     0     1     0     0
#> basket_2     1    0    1    0    0    0    1    0    0     0     1     0     0
#> basket_3     1    0    1    0    0    0    1    0    0     0     1     0     0
#> basket_4     1    0    0    0    0    0    0    0    0     1     1     0     0
#> basket_5     1    0    0    0    0    0    0    0    0     1     1     0     0
#> basket_6     1    0    1    0    0    0    0    0    0     0     0     1     1
#> basket_7     1    0    1    0    0    0    0    0    0     0     0     1     1
#> basket_8     0    0    0    1    1    1    1    1    0     0     1     0     0
#> basket_9     0    0    0    1    1    1    1    1    0     0     1     0     0
#> basket_10    1    1    0    0    0    0    0    1    1     0     0     0     0
#>           fe_14 fe_15 fe_16 fe_17 fe_18 fe_19 fe_20 fe_21 fe_22 fe_23 fe_24
#> basket_1      0     1     0     1     0     1     0     0     0     1     0
#> basket_2      0     1     0     1     0     1     0     0     0     1     0
#> basket_3      0     1     0     1     0     1     0     0     0     1     0
#> basket_4      0     0     1     0     1     1     0     1     0     0     1
#> basket_5      0     0     1     0     1     1     0     1     0     0     1
#> basket_6      1     0     1     0     0     0     1     0     0     0     0
#> basket_7      1     0     1     0     0     0     1     0     0     0     0
#> basket_8      0     0     0     0     0     0     0     0     1     1     0
#> basket_9      0     0     0     0     0     0     0     0     1     1     0
#> basket_10     0     1     0     0     0     1     0     1     0     1     0
#>           fe_25
#> basket_1      0
#> basket_2      0
#> basket_3      0
#> basket_4      0
#> basket_5      0
#> basket_6      1
#> basket_7      1
#> basket_8      0
#> basket_9      0
#> basket_10     0
#> 
#> 
```
