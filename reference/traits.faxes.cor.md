# Correlation between Traits and Axes

Compute relationship between all traits and all axes of the functional
space. For continuous trait a linear model is computed and r2 and
p-value are returned. For other types of traits, a Kruskal-Wallis test
is computed and eta2 statistics is returned. Option allows to plot
trait-axis relationships with scatterplot and boxplot for continuous and
non-continuous traits, respectively.

## Usage

``` r
traits.faxes.cor(
  sp_tr,
  sp_faxes_coord,
  tr_nm = NULL,
  faxes_nm = NULL,
  plot = FALSE,
  name_file = NULL,
  color_signif = "darkblue",
  color_non_signif = "gray80",
  stop_if_NA = TRUE
)
```

## Arguments

- sp_tr:

  a data frame containing species as rows and traits as columns.

- sp_faxes_coord:

  a matrix of species coordinates in a multidimensional functional
  space. Species coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- tr_nm:

  a vector gathering the names of traits (as in `sp_tr`) to consider. If
  `NULL` all traits are considered.

- faxes_nm:

  a vector gathering the names of PCoA axes (as in `sp_faxes_coord`) to
  consider.

- plot:

  a logical value indicating whether plot illustrating relations between
  trait and axes should be drawn. **You can only plot relationships for
  up to 10 traits and/or 10 axes**.

- name_file:

  the file name (without extension) to save the plot as a 300 dpi JPEG
  file. Default is `NULL` which means plot is only displayed. If
  `plot = FALSE` this argument is ignored.

- color_signif:

  an R color name or an hexadecimal code referring to the color of
  points when relationships between the trait and the axis is
  significant. Default is `darkblue`.

- color_non_signif:

  an R color name or an hexadecimal code referring to the color of
  points when relationships between the trait and the axis are not
  significant. Default is `gray80`.

- stop_if_NA:

  a logical value to stop or not the process if the `sp_tr` data frame
  contains NA. Functional measures are sensitive to missing traits. For
  further explanations, see the Note section. Default is `TRUE`.

## Value

1 data frame with for each combination of trait and axis (rows), the
name of the test performed, and the corresponding statistics and
p-value. If `plot = TRUE` a multi-panel figure with traits as columns
and axes as rows is also plotted. When relationships between trait and
axis is significant the points are colored, else they remain grayish.

## Author

Nicolas Loiseau and Sebastien Villeger

## Examples

``` r
# Load Species x Traits Data
data("fruits_traits", package = "mFD")

# Load Traits categories dataframe
data("fruits_traits_cat", package = "mFD")

# Compute Functional Distance
sp_dist_fruits <- mFD::funct.dist(sp_tr  = fruits_traits,
 tr_cat = fruits_traits_cat,
 metric = "gower",
 scale_euclid  = "scale_center",
 ordinal_var = "classic",
 weight_type = "equal",
 stop_if_NA  = TRUE)
#> [1] "Running w.type=equal on groups=c(Size)"
#> [1] "Running w.type=equal on groups=c(Plant)"
#> [1] "Running w.type=equal on groups=c(Climate)"
#> [1] "Running w.type=equal on groups=c(Seed)"
#> [1] "Running w.type=equal on groups=c(Sugar)"
#> [1] "Running w.type=equal on groups=c(Use,Use,Use)"
  
# Compute Functional Spaces Quality (to retrieve species coordinates)
fspaces_quality_fruits <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fruits, 
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")
  
# Retrieve Species Coordinates
sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord

# Compute Correlation between Traits and Functional Axes
mFD::traits.faxes.cor(
  sp_tr          = fruits_traits, 
  sp_faxes_coord = sp_faxes_coord_fruits, 
  tr_nm          = NULL, 
  faxes_nm       = NULL,
  name_file      = NULL, 
  color_signif   = "darkblue",
  color_non_signif = "gray80")
#>         trait axis           test stat  value p.value
#> 1        Size  PC1 Kruskal-Wallis eta2  0.284  0.0461
#> 2        Size  PC2 Kruskal-Wallis eta2  0.316  0.0352
#> 3        Size  PC3 Kruskal-Wallis eta2  0.008  0.3839
#> 4        Size  PC4 Kruskal-Wallis eta2  0.285  0.0459
#> 5        Size  PC5 Kruskal-Wallis eta2  0.151  0.1349
#> 6        Size  PC6 Kruskal-Wallis eta2  0.312  0.0366
#> 7        Size  PC7 Kruskal-Wallis eta2  0.100  0.1990
#> 8        Size  PC8 Kruskal-Wallis eta2  0.164  0.1216
#> 9        Size  PC9 Kruskal-Wallis eta2  0.299  0.0408
#> 10       Size PC10 Kruskal-Wallis eta2  0.132  0.1564
#> 11      Plant  PC1 Kruskal-Wallis eta2  0.315  0.0222
#> 12      Plant  PC2 Kruskal-Wallis eta2  0.413  0.0086
#> 13      Plant  PC3 Kruskal-Wallis eta2  0.197  0.0675
#> 14      Plant  PC4 Kruskal-Wallis eta2 -0.057  0.6157
#> 15      Plant  PC5 Kruskal-Wallis eta2  0.266  0.0354
#> 16      Plant  PC6 Kruskal-Wallis eta2  0.089  0.1818
#> 17      Plant  PC7 Kruskal-Wallis eta2 -0.087  0.7598
#> 18      Plant  PC8 Kruskal-Wallis eta2  0.187  0.0745
#> 19      Plant  PC9 Kruskal-Wallis eta2 -0.014  0.4378
#> 20      Plant PC10 Kruskal-Wallis eta2 -0.138  0.9927
#> 21    Climate  PC1 Kruskal-Wallis eta2  0.746  0.0001
#> 22    Climate  PC2 Kruskal-Wallis eta2 -0.048  0.6216
#> 23    Climate  PC3 Kruskal-Wallis eta2  0.030  0.2636
#> 24    Climate  PC4 Kruskal-Wallis eta2  0.194  0.0433
#> 25    Climate  PC5 Kruskal-Wallis eta2 -0.015  0.4334
#> 26    Climate  PC6 Kruskal-Wallis eta2 -0.032  0.5232
#> 27    Climate  PC7 Kruskal-Wallis eta2 -0.012  0.4199
#> 28    Climate  PC8 Kruskal-Wallis eta2 -0.022  0.4706
#> 29    Climate  PC9 Kruskal-Wallis eta2  0.116  0.1026
#> 30    Climate PC10 Kruskal-Wallis eta2  0.024  0.2821
#> 31       Seed  PC1 Kruskal-Wallis eta2  0.111  0.1082
#> 32       Seed  PC2 Kruskal-Wallis eta2  0.475  0.0020
#> 33       Seed  PC3 Kruskal-Wallis eta2  0.435  0.0031
#> 34       Seed  PC4 Kruskal-Wallis eta2  0.057  0.1957
#> 35       Seed  PC5 Kruskal-Wallis eta2  0.034  0.2524
#> 36       Seed  PC6 Kruskal-Wallis eta2 -0.074  0.8344
#> 37       Seed  PC7 Kruskal-Wallis eta2 -0.081  0.8981
#> 38       Seed  PC8 Kruskal-Wallis eta2 -0.040  0.5700
#> 39       Seed  PC9 Kruskal-Wallis eta2 -0.040  0.5700
#> 40       Seed PC10 Kruskal-Wallis eta2 -0.063  0.7397
#> 41      Sugar  PC1   Linear Model   r2  0.089  0.1486
#> 42      Sugar  PC2   Linear Model   r2  0.195  0.0273
#> 43      Sugar  PC3   Linear Model   r2  0.218  0.0187
#> 44      Sugar  PC4   Linear Model   r2  0.022  0.4814
#> 45      Sugar  PC5   Linear Model   r2  0.094  0.1357
#> 46      Sugar  PC6   Linear Model   r2  0.305  0.0042
#> 47      Sugar  PC7   Linear Model   r2  0.032  0.3948
#> 48      Sugar  PC8   Linear Model   r2  0.004  0.7668
#> 49      Sugar  PC9   Linear Model   r2  0.000  0.9473
#> 50      Sugar PC10   Linear Model   r2  0.000  0.9252
#> 51    Use.raw  PC1   Linear Model   r2  0.497  0.0001
#> 52    Use.raw  PC2   Linear Model   r2  0.004  0.7551
#> 53    Use.raw  PC3   Linear Model   r2  0.103  0.1168
#> 54    Use.raw  PC4   Linear Model   r2  0.218  0.0187
#> 55    Use.raw  PC5   Linear Model   r2  0.002  0.8154
#> 56    Use.raw  PC6   Linear Model   r2  0.117  0.0947
#> 57    Use.raw  PC7   Linear Model   r2  0.023  0.4740
#> 58    Use.raw  PC8   Linear Model   r2  0.000  0.9295
#> 59    Use.raw  PC9   Linear Model   r2  0.011  0.6154
#> 60    Use.raw PC10   Linear Model   r2  0.005  0.7474
#> 61 Use.pastry  PC1   Linear Model   r2  0.031  0.4036
#> 62 Use.pastry  PC2   Linear Model   r2  0.000  0.9630
#> 63 Use.pastry  PC3   Linear Model   r2  0.061  0.2349
#> 64 Use.pastry  PC4   Linear Model   r2  0.360  0.0015
#> 65 Use.pastry  PC5   Linear Model   r2  0.283  0.0062
#> 66 Use.pastry  PC6   Linear Model   r2  0.091  0.1432
#> 67 Use.pastry  PC7   Linear Model   r2  0.100  0.1228
#> 68 Use.pastry  PC8   Linear Model   r2  0.002  0.8436
#> 69 Use.pastry  PC9   Linear Model   r2  0.004  0.7750
#> 70 Use.pastry PC10   Linear Model   r2  0.010  0.6325
#> 71    Use.jam  PC1   Linear Model   r2  0.567  0.0000
#> 72    Use.jam  PC2   Linear Model   r2  0.006  0.7239
#> 73    Use.jam  PC3   Linear Model   r2  0.055  0.2588
#> 74    Use.jam  PC4   Linear Model   r2  0.033  0.3835
#> 75    Use.jam  PC5   Linear Model   r2  0.082  0.1655
#> 76    Use.jam  PC6   Linear Model   r2  0.050  0.2831
#> 77    Use.jam  PC7   Linear Model   r2  0.153  0.0533
#> 78    Use.jam  PC8   Linear Model   r2  0.003  0.8123
#> 79    Use.jam  PC9   Linear Model   r2  0.029  0.4193
#> 80    Use.jam PC10   Linear Model   r2  0.000  0.9327
```
