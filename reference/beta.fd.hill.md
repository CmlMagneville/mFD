# Compute Functional beta-Diversity indices based on Hill Numbers

Compute functional beta-diversity indices based on Hill numbers applied
to distance between species following the framework from Chao *et al.*
(2019).

## Usage

``` r
beta.fd.hill(
  asb_sp_w,
  sp_dist,
  q = c(0, 1, 2),
  tau = "mean",
  beta_type = "Jaccard",
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

  a vector containing values referring to the order of diversity to use

- tau:

  a character string with name of function to apply to distance matrix
  (i.e. among all pairs of species) to get the threshold used to define
  'functionally indistinct set of species'. Could be qet to 'mean'
  (default), 'min' or 'max'.

- beta_type:

  a character string with name of framework used for computing
  beta-diversity, either 'Jaccard' (default) or 'Sorensen'.

- check_input:

  a logical value indicating whether key features the inputs are checked
  (e.g. class and/or mode of objects, names of rows and/or columns,
  missing values). If an error is detected, a detailed message is
  returned. Default: `check_input = TRUE`.

- details_returned:

  a logical value indicating whether the user want to store values used
  for computing indices (see list below)

## Value

A list with:

- *asb_FDbeta* a list with for each value of q a *dist* object with beta
  functional diversity indices for all pairs of assemblages item if
  **store.details** turned to TRUE a list *details* with

  - *malpha_fd_q* a list with for each value of q a *dist* object with
    mean alpha functional diversity indices for all pairs of assemblages

  - *gamma_fd_q* a list with for each value of q a *dist* object with
    gamma functional diversity indices for all pairs of assemblages

## Note

When q=1 Jaccard-like and Sorensen-like beta-diversity are identical. FD
computed with tau='min' is equivalent to Hill number taxonomic beta
diversity. If tau='min' and there are species with null distance, tau is
set to the minimum non-null value and a warning message is displayed.
Indices values are stored as *dist* objects to optimize memory. See
below example of how merging distance values in a *dataframe* with
[`dist.to.df`](https://cmlmagneville.github.io/mFD/reference/dist.to.df.md)

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

# Load Traits types dataframe:      
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

# Compute beta functional hill indices:
baskets_beta <- beta.fd.hill(
      asb_sp_w         = baskets_fruits_weights, 
      sp_dist          = sp_dist_fruits, 
      q                = c(0,1,2), 
      tau              = 'mean',
      beta_type        = 'Jaccard', 
      check_input      = TRUE, 
      details_returned = TRUE)
 
# Then use the mFD::dist.to.df function to ease visualizing result:
## for q = 0:
mFD::dist.to.df(list_dist = list(FDq2 = baskets_beta$beta_fd_q$q0))
#>          x1        x2       FDq2
#> 1  basket_1  basket_2 0.00000000
#> 2  basket_1  basket_3 0.00000000
#> 3  basket_1  basket_4 0.00000000
#> 4  basket_1  basket_5 0.00000000
#> 5  basket_1  basket_6 0.00000000
#> 6  basket_1  basket_7 0.00000000
#> 7  basket_1  basket_8 0.12848199
#> 8  basket_1  basket_9 0.12605403
#> 9  basket_1 basket_10 0.00000000
#> 10 basket_2  basket_3 0.00000000
#> 11 basket_2  basket_4 0.00000000
#> 12 basket_2  basket_5 0.00000000
#> 13 basket_2  basket_6 0.00000000
#> 14 basket_2  basket_7 0.00000000
#> 15 basket_2  basket_8 0.10650862
#> 16 basket_2  basket_9 0.10367515
#> 17 basket_2 basket_10 0.00000000
#> 18 basket_3  basket_4 0.00000000
#> 19 basket_3  basket_5 0.00000000
#> 20 basket_3  basket_6 0.00000000
#> 21 basket_3  basket_7 0.00000000
#> 22 basket_3  basket_8 0.10044178
#> 23 basket_3  basket_9 0.09776647
#> 24 basket_3 basket_10 0.00000000
#> 25 basket_4  basket_5 0.00000000
#> 26 basket_4  basket_6 0.18524852
#> 27 basket_4  basket_7 0.18088474
#> 28 basket_4  basket_8 0.00000000
#> 29 basket_4  basket_9 0.00000000
#> 30 basket_4 basket_10 0.00000000
#> 31 basket_5  basket_6 0.17793331
#> 32 basket_5  basket_7 0.17358512
#> 33 basket_5  basket_8 0.00000000
#> 34 basket_5  basket_9 0.00000000
#> 35 basket_5 basket_10 0.00000000
#> 36 basket_6  basket_7 0.00000000
#> 37 basket_6  basket_8 0.13417736
#> 38 basket_6  basket_9 0.13382527
#> 39 basket_6 basket_10 0.00000000
#> 40 basket_7  basket_8 0.13138912
#> 41 basket_7  basket_9 0.13220347
#> 42 basket_7 basket_10 0.00000000
#> 43 basket_8  basket_9 0.00000000
#> 44 basket_8 basket_10 0.00000000
#> 45 basket_9 basket_10 0.00000000
## for q = 1:
mFD::dist.to.df(list_dist = list(FDq2 = baskets_beta$beta_fd_q$q1))
#>          x1        x2         FDq2
#> 1  basket_1  basket_2 0.0465649552
#> 2  basket_1  basket_3 0.0633934986
#> 3  basket_1  basket_4 0.0390229310
#> 4  basket_1  basket_5 0.0416185003
#> 5  basket_1  basket_6 0.4025770889
#> 6  basket_1  basket_7 0.2582156175
#> 7  basket_1  basket_8 0.2549851440
#> 8  basket_1  basket_9 0.1703571633
#> 9  basket_1 basket_10 0.0172134140
#> 10 basket_2  basket_3 0.0029886932
#> 11 basket_2  basket_4 0.1315828608
#> 12 basket_2  basket_5 0.1434349036
#> 13 basket_2  basket_6 0.2382682836
#> 14 basket_2  basket_7 0.1167463677
#> 15 basket_2  basket_8 0.2909423919
#> 16 basket_2  basket_9 0.2449369402
#> 17 basket_2 basket_10 0.0780014178
#> 18 basket_3  basket_4 0.1417318605
#> 19 basket_3  basket_5 0.1583078065
#> 20 basket_3  basket_6 0.2056050375
#> 21 basket_3  basket_7 0.1045895081
#> 22 basket_3  basket_8 0.3043051662
#> 23 basket_3  basket_9 0.2597238131
#> 24 basket_3 basket_10 0.0987942065
#> 25 basket_4  basket_5 0.0008213232
#> 26 basket_4  basket_6 0.4526919361
#> 27 basket_4  basket_7 0.3719680099
#> 28 basket_4  basket_8 0.3961862747
#> 29 basket_4  basket_9 0.2769826915
#> 30 basket_4 basket_10 0.0440610254
#> 31 basket_5  basket_6 0.4788884130
#> 32 basket_5  basket_7 0.3962066363
#> 33 basket_5  basket_8 0.4087861766
#> 34 basket_5  basket_9 0.2759382080
#> 35 basket_5 basket_10 0.0404365649
#> 36 basket_6  basket_7 0.0545121456
#> 37 basket_6  basket_8 0.6639322400
#> 38 basket_6  basket_9 0.5921222975
#> 39 basket_6 basket_10 0.3977446026
#> 40 basket_7  basket_8 0.4726732200
#> 41 basket_7  basket_9 0.4300437283
#> 42 basket_7 basket_10 0.2751488242
#> 43 basket_8  basket_9 0.0303705511
#> 44 basket_8 basket_10 0.2652284919
#> 45 basket_9 basket_10 0.1601080610
## for q = 2:
mFD::dist.to.df(list_dist = list(FDq2 = baskets_beta$beta_fd_q$q2))
#>          x1        x2        FDq2
#> 1  basket_1  basket_2 0.058982325
#> 2  basket_1  basket_3 0.078716397
#> 3  basket_1  basket_4 0.029573623
#> 4  basket_1  basket_5 0.027059789
#> 5  basket_1  basket_6 0.484115290
#> 6  basket_1  basket_7 0.292594562
#> 7  basket_1  basket_8 0.290545545
#> 8  basket_1  basket_9 0.185475113
#> 9  basket_1 basket_10 0.011136995
#> 10 basket_2  basket_3 0.004420448
#> 11 basket_2  basket_4 0.161833512
#> 12 basket_2  basket_5 0.162571972
#> 13 basket_2  basket_6 0.260541701
#> 14 basket_2  basket_7 0.097053161
#> 15 basket_2  basket_8 0.294504888
#> 16 basket_2  basket_9 0.225897615
#> 17 basket_2 basket_10 0.058298760
#> 18 basket_3  basket_4 0.172877455
#> 19 basket_3  basket_5 0.178123024
#> 20 basket_3  basket_6 0.207102590
#> 21 basket_3  basket_7 0.081951839
#> 22 basket_3  basket_8 0.308336365
#> 23 basket_3  basket_9 0.241482168
#> 24 basket_3 basket_10 0.082649928
#> 25 basket_4  basket_5 0.001049851
#> 26 basket_4  basket_6 0.511165067
#> 27 basket_4  basket_7 0.421141181
#> 28 basket_4  basket_8 0.482330219
#> 29 basket_4  basket_9 0.342926459
#> 30 basket_4 basket_10 0.050817451
#> 31 basket_5  basket_6 0.532084800
#> 32 basket_5  basket_7 0.438841544
#> 33 basket_5  basket_8 0.496554512
#> 34 basket_5  basket_9 0.336894275
#> 35 basket_5 basket_10 0.044052657
#> 36 basket_6  basket_7 0.068382884
#> 37 basket_6  basket_8 0.759422492
#> 38 basket_6  basket_9 0.680332414
#> 39 basket_6 basket_10 0.453325136
#> 40 basket_7  basket_8 0.478528108
#> 41 basket_7  basket_9 0.431531941
#> 42 basket_7 basket_10 0.265928889
#> 43 basket_8  basket_9 0.020812705
#> 44 basket_8 basket_10 0.345652088
#> 45 basket_9 basket_10 0.219780509
```
