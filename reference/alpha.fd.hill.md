# Compute Functional alpha-Diversity indices based on Hill Numbers

Compute functional alpha diversity applied to distance between species
following the framework from Chao *et al.*(2019).

## Usage

``` r
alpha.fd.hill(
  asb_sp_w,
  sp_dist,
  q = c(0, 1, 2),
  tau = "mean",
  check_input = TRUE,
  details_returned = TRUE
)
```

## Arguments

- asb_sp_w:

  a matrix with weight of species (columns) in a set of assemblages
  (rows). Rows and columns should have names. NA are not allowed.

- sp_dist:

  a matrix or dist object with distance between species. Species names
  should be provided and match those in 'asb_sp_w'. NA are not allowed.

- q:

  a vector containing values referring to the order of diversity to
  consider, could be 0, 1 and/or 2.

- tau:

  a character string with name of function to apply to distance matrix
  (i.e. among all pairs of species) to get the threshold used to define
  'functionally indistinct set of species'. Could be 'mean' (default),
  'min' or 'max'. If tau = 'min' and there are null values in `sp_dist`,
  the threshold is the lowest strictly positive value and a warning
  message is displayed.

- check_input:

  a logical value indicating whether key features the inputs are checked
  (e.g. class and/or mode of objects, names of rows and/or columns,
  missing values). If an error is detected, a detailed message is
  returned. Default: `check.input = TRUE`.

- details_returned:

  a logical value indicating whether the user want to store values used
  for computing indices (see list below)

## Value

A list with:

- *asb_FD_Hill* a matrix containing indices values for each level of q
  (columns, named as 'FD_qx') for each assemblage (rows, named as in
  **asb_sp_w**)

- *tau_dist* the threshold value applied to distance between species to
  compute diversity according to function provided in **tau**

- if **details_returned** turned to TRUE a list *details* with

  - *asb_totw* a vector with total weight of each assemblage

  - *asb_sp_relw* a matrix with relative weight of species in
    assemblages

## Note

FD is computed applying the special case where function 'f' in equation
3c is linear:f(dij(tau)) = dij(tau)/tau, hence f(0) = 0 and f(tau) = 1.
FD computed with q=2 and tau = 'max' is equivalent to the Rao's
quadratic entropy from Ricotta & Szeidl (2009, J Theor Biol). FD
computed with tau = 'min' is equivalent to Hill number taxonomic
diversity, thus with q=0 it is species richness (S), with q = 1 it is
exponential of Shannon entropy (H) and with q = 2 it is 1/(1-D) where D
is Simpson diversity. Note that even when q=0, weights of species are
accounted for in FD. Hence to compute FD based only on distance between
species present in an assemblage (i.e. a richness-like index) , asb_sp_w
has to contain only species presence/absence coded as 0/1 with q=0 and
tau = 'mean'. If asb_sp_w contains only 0/1 and q\>0, it means that all
species have the same contribution to FD.

## References

Chao *et al.* (2019) An attribute diversity approach to functional
diversity, functional beta diversity, and related (dis)similarity
measures. *Ecological Monographs*, **89**, e01343.

## Author

Sebastien Villeger and Camille Magneville

## Examples

``` r
# Load Species*Traits dataframe:
data('fruits_traits', package = 'mFD')

# Load Assemblages*Species dataframe:      
data('baskets_fruits_weights', package = 'mFD') 
  
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

# Compute alpha fd hill indices:
alpha.fd.hill(
   asb_sp_w         = baskets_fruits_weights, 
   sp_dist          = sp_dist_fruits, 
   q                = c(0, 1, 2),
   tau              = 'mean', 
   check_input      = TRUE, 
   details_returned = TRUE)
#> $asb_FD_Hill
#>              FD_q0    FD_q1    FD_q2
#> basket_1  3.912219 3.082878 2.615626
#> basket_2  4.495517 4.198786 3.972157
#> basket_3  4.493944 4.266210 4.101477
#> basket_4  1.758038 1.740013 1.722964
#> basket_5  1.891440 1.872348 1.853837
#> basket_6  4.513411 4.049255 3.730993
#> basket_7  4.474346 4.175529 3.944302
#> basket_8  3.367259 2.516060 2.156245
#> basket_9  3.393249 3.039767 2.749451
#> basket_10 3.378281 3.168822 3.006560
#> 
#> $tau_dist
#> [1] 0.3593121
#> 
#> $details
#> $details$asb_totw
#>  basket_1  basket_2  basket_3  basket_4  basket_5  basket_6  basket_7  basket_8 
#>      2000      2000      2000      2000      2000      2000      2000      2000 
#>  basket_9 basket_10 
#>      2000      2000 
#> 
#> $details$asb_sp_relw
#>           apple apricot banana currant blackberry blueberry cherry grape
#> basket_1  0.200     0.0   0.05    0.00       0.00      0.00  0.075  0.00
#> basket_2  0.100     0.0   0.20    0.00       0.00      0.00  0.125  0.00
#> basket_3  0.100     0.0   0.25    0.00       0.00      0.00  0.125  0.00
#> basket_4  0.150     0.0   0.00    0.00       0.00      0.00  0.000  0.00
#> basket_5  0.100     0.0   0.00    0.00       0.00      0.00  0.000  0.00
#> basket_6  0.050     0.0   0.10    0.00       0.00      0.00  0.000  0.00
#> basket_7  0.050     0.0   0.10    0.00       0.00      0.00  0.000  0.00
#> basket_8  0.000     0.0   0.00    0.10       0.15      0.10  0.100  0.10
#> basket_9  0.000     0.0   0.00    0.05       0.05      0.05  0.050  0.20
#> basket_10 0.175     0.1   0.00    0.00       0.00      0.00  0.000  0.15
#>           grapefruit kiwifruit lemon lime litchi mango melon orange
#> basket_1        0.00      0.00  0.10  0.0   0.00  0.00  0.10   0.00
#> basket_2        0.00      0.00  0.05  0.0   0.00  0.00  0.25   0.00
#> basket_3        0.00      0.00  0.05  0.0   0.00  0.00  0.20   0.00
#> basket_4        0.00      0.05  0.05  0.0   0.00  0.00  0.00   0.20
#> basket_5        0.00      0.15  0.15  0.0   0.00  0.00  0.00   0.15
#> basket_6        0.00      0.00  0.00  0.1   0.10  0.25  0.00   0.05
#> basket_7        0.00      0.00  0.00  0.1   0.05  0.10  0.00   0.05
#> basket_8        0.00      0.00  0.05  0.0   0.00  0.00  0.00   0.00
#> basket_9        0.00      0.00  0.15  0.0   0.00  0.00  0.00   0.00
#> basket_10       0.15      0.00  0.00  0.0   0.00  0.00  0.20   0.00
#>           passion_fruit peach pear pineapple  plum raspberry strawberry
#> basket_1           0.05  0.00 0.30      0.00 0.000      0.00      0.125
#> basket_2           0.05  0.00 0.10      0.00 0.000      0.00      0.125
#> basket_3           0.05  0.00 0.10      0.00 0.000      0.00      0.125
#> basket_4           0.00  0.15 0.20      0.00 0.100      0.00      0.000
#> basket_5           0.00  0.15 0.15      0.00 0.100      0.00      0.000
#> basket_6           0.00  0.00 0.00      0.25 0.000      0.00      0.000
#> basket_7           0.00  0.00 0.00      0.25 0.000      0.00      0.000
#> basket_8           0.00  0.00 0.00      0.00 0.000      0.20      0.200
#> basket_9           0.00  0.00 0.00      0.00 0.000      0.25      0.200
#> basket_10          0.00  0.00 0.10      0.00 0.075      0.00      0.050
#>           tangerine water_melon
#> basket_1       0.00         0.0
#> basket_2       0.00         0.0
#> basket_3       0.00         0.0
#> basket_4       0.10         0.0
#> basket_5       0.05         0.0
#> basket_6       0.00         0.1
#> basket_7       0.00         0.3
#> basket_8       0.00         0.0
#> basket_9       0.00         0.0
#> basket_10      0.00         0.0
#> 
#> 
```
