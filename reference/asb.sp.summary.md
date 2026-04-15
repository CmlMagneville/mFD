# Summarize Assemblage x Species data frame

This function computes a summary helping you to picture assemblages. For
this function to work, there must be no NA in your `asb_sp_w` data
frame.

## Usage

``` r
asb.sp.summary(asb_sp_w)
```

## Arguments

- asb_sp_w:

  a matrix showing assemblages (rows) composition in species (columns).
  Note that species names **must be** the names of rows.

## Value

A list with:

- asb_sp_w_occ:

  a matrix with species occurrences in each assemblage.

- sp_tot_w:

  a vector gathering species biomass/abundance per species.

- asb_tot_w:

  a vector gathering total abundance/biomass per assemblage.

- asb_sp_richn:

  a vector gathering species richness per assemblage.

- asb_sp_nm:

  a list gathering the names of species of each assemblage.

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
# Load Assemblages x Species Matrix
data('baskets_fruits_weights', package = 'mFD')

# Summarize Assemblages Data
mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights)
#> $asb_sp_occ
#>           apple apricot banana currant blackberry blueberry cherry grape
#> basket_1      1       0      1       0          0         0      1     0
#> basket_2      1       0      1       0          0         0      1     0
#> basket_3      1       0      1       0          0         0      1     0
#> basket_4      1       0      0       0          0         0      0     0
#> basket_5      1       0      0       0          0         0      0     0
#> basket_6      1       0      1       0          0         0      0     0
#> basket_7      1       0      1       0          0         0      0     0
#> basket_8      0       0      0       1          1         1      1     1
#> basket_9      0       0      0       1          1         1      1     1
#> basket_10     1       1      0       0          0         0      0     1
#>           grapefruit kiwifruit lemon lime litchi mango melon orange
#> basket_1           0         0     1    0      0     0     1      0
#> basket_2           0         0     1    0      0     0     1      0
#> basket_3           0         0     1    0      0     0     1      0
#> basket_4           0         1     1    0      0     0     0      1
#> basket_5           0         1     1    0      0     0     0      1
#> basket_6           0         0     0    1      1     1     0      1
#> basket_7           0         0     0    1      1     1     0      1
#> basket_8           0         0     1    0      0     0     0      0
#> basket_9           0         0     1    0      0     0     0      0
#> basket_10          1         0     0    0      0     0     1      0
#>           passion_fruit peach pear pineapple plum raspberry strawberry
#> basket_1              1     0    1         0    0         0          1
#> basket_2              1     0    1         0    0         0          1
#> basket_3              1     0    1         0    0         0          1
#> basket_4              0     1    1         0    1         0          0
#> basket_5              0     1    1         0    1         0          0
#> basket_6              0     0    0         1    0         0          0
#> basket_7              0     0    0         1    0         0          0
#> basket_8              0     0    0         0    0         1          1
#> basket_9              0     0    0         0    0         1          1
#> basket_10             0     0    1         0    1         0          1
#>           tangerine water_melon
#> basket_1          0           0
#> basket_2          0           0
#> basket_3          0           0
#> basket_4          1           0
#> basket_5          1           0
#> basket_6          0           1
#> basket_7          0           1
#> basket_8          0           0
#> basket_9          0           0
#> basket_10         0           0
#> 
#> $sp_tot_w
#>         apple       apricot        banana       currant    blackberry 
#>          1850           200          1400           300           400 
#>     blueberry        cherry         grape    grapefruit     kiwifruit 
#>           300           950           900           300           400 
#>         lemon          lime        litchi         mango         melon 
#>          1200           400           300           700          1500 
#>        orange passion_fruit         peach          pear     pineapple 
#>           900           300           600          1900          1000 
#>          plum     raspberry    strawberry     tangerine   water_melon 
#>           550           900          1650           300           800 
#> 
#> $asb_tot_w
#>  basket_1  basket_2  basket_3  basket_4  basket_5  basket_6  basket_7  basket_8 
#>      2000      2000      2000      2000      2000      2000      2000      2000 
#>  basket_9 basket_10 
#>      2000      2000 
#> 
#> $asb_sp_richn
#>  basket_1  basket_2  basket_3  basket_4  basket_5  basket_6  basket_7  basket_8 
#>         8         8         8         8         8         8         8         8 
#>  basket_9 basket_10 
#>         8         8 
#> 
#> $asb_sp_nm
#> $asb_sp_nm$basket_1
#>         apple        banana        cherry         lemon         melon 
#>             1             1             1             1             1 
#> passion_fruit          pear    strawberry 
#>             1             1             1 
#> 
#> $asb_sp_nm$basket_2
#>         apple        banana        cherry         lemon         melon 
#>             1             1             1             1             1 
#> passion_fruit          pear    strawberry 
#>             1             1             1 
#> 
#> $asb_sp_nm$basket_3
#>         apple        banana        cherry         lemon         melon 
#>             1             1             1             1             1 
#> passion_fruit          pear    strawberry 
#>             1             1             1 
#> 
#> $asb_sp_nm$basket_4
#>     apple kiwifruit     lemon    orange     peach      pear      plum tangerine 
#>         1         1         1         1         1         1         1         1 
#> 
#> $asb_sp_nm$basket_5
#>     apple kiwifruit     lemon    orange     peach      pear      plum tangerine 
#>         1         1         1         1         1         1         1         1 
#> 
#> $asb_sp_nm$basket_6
#>       apple      banana        lime      litchi       mango      orange 
#>           1           1           1           1           1           1 
#>   pineapple water_melon 
#>           1           1 
#> 
#> $asb_sp_nm$basket_7
#>       apple      banana        lime      litchi       mango      orange 
#>           1           1           1           1           1           1 
#>   pineapple water_melon 
#>           1           1 
#> 
#> $asb_sp_nm$basket_8
#>    currant blackberry  blueberry     cherry      grape      lemon  raspberry 
#>          1          1          1          1          1          1          1 
#> strawberry 
#>          1 
#> 
#> $asb_sp_nm$basket_9
#>    currant blackberry  blueberry     cherry      grape      lemon  raspberry 
#>          1          1          1          1          1          1          1 
#> strawberry 
#>          1 
#> 
#> $asb_sp_nm$basket_10
#>      apple    apricot      grape grapefruit      melon       pear       plum 
#>          1          1          1          1          1          1          1 
#> strawberry 
#>          1 
#> 
#> 
```
