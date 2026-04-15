# Compute Functional Diversity Hill Indices

## About this tutorial

### What is this tutorial about?

  

This tutorial explains how to compute the family of indices presented in
[Chao *et al.*
(2019)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343)
using `mFD`.

  

### Let’s load data and compute functional distance

  

The data set used to illustrate this tutorial is the **fruits dataset**
based on 25 types of fruits (*i.e.* species) distributed in 10 fruits
baskets (*i.e.* assemblages). Each fruit is characterized by six traits
values summarized in the following table:

| Trait name | Trait measurement | Trait type  | Number of classes |            Classes code            | Unit |
|:----------:|:-----------------:|:-----------:|:-----------------:|:----------------------------------:|:----:|
|    Size    | Maximal diameter  |   Ordinal   |         5         |   0-1 ; 1-3 ; 3-5 ; 5-10 ; 10-20   |  cm  |
|   Plant    |    Growth form    | Categorical |         4         |      tree; shrub; vine; forb       |  NA  |
|  Climate   |  Climatic niche   |   Ordinal   |         3         | temperate ; subtropical ; tropical |  NA  |
|    Seed    |     Seed type     |   Ordinal   |         3         |          none ; pip ; pit          |  NA  |
|   Sugar    |       Sugar       | Continuous  |        NA         |                 NA                 | g/kg |
|    Use     |    Use as food    |    Fuzzy    |         3         |         raw ; pastry ; jam         |  %   |

  

We load the three objects used to compute functional framework (for more
explanations, see [mFD General
Workflow](https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html)):

  

- a **data frame** summarizing traits values for each species called
  `fruits_traits` in this tutorial:

``` r
data("fruits_traits", package = "mFD")

knitr::kable(head(fruits_traits),
      caption = "Species x traits data frame based on the **fruits** dataset")
```

|            | Size    | Plant | Climate   | Seed | Sugar | Use.raw | Use.pastry | Use.jam |
|:-----------|:--------|:------|:----------|:-----|------:|--------:|-----------:|--------:|
| apple      | 5-10cm  | tree  | temperate | pip  | 103.9 |      50 |         50 |       0 |
| apricot    | 3-5cm   | tree  | temperate | pit  |  92.4 |      40 |         10 |      50 |
| banana     | 10-20cm | tree  | tropical  | none | 122.3 |      70 |         20 |      10 |
| currant    | 0-1cm   | shrub | temperate | pip  |  73.7 |      10 |         10 |      80 |
| blackberry | 1-3cm   | shrub | temperate | pip  |  48.8 |      30 |         10 |      60 |
| blueberry  | 0-1cm   | forb  | temperate | pip  | 100.0 |      10 |         40 |      50 |

Species x traits data frame based on the **fruits** dataset

  

- a **matrix** summarizing species assemblages called
  `baskets_fruits_weights` in this tutorial. Weights in this matrix can
  be occurrence data, abundance, biomass, coverage, etc. The studied
  example works with biomass (*i.e.* grams of a fruit in a basket) and
  this matrix looks as follows:

``` r
data("baskets_fruits_weights", package = "mFD")

knitr::kable(as.data.frame(baskets_fruits_weights[1:6, 1:6]), 
             caption = "Species x assemblages matrix based on the **fruits** dataset")
```

|          | apple | apricot | banana | currant | blackberry | blueberry |
|:---------|------:|--------:|-------:|--------:|-----------:|----------:|
| basket_1 |   400 |       0 |    100 |       0 |          0 |         0 |
| basket_2 |   200 |       0 |    400 |       0 |          0 |         0 |
| basket_3 |   200 |       0 |    500 |       0 |          0 |         0 |
| basket_4 |   300 |       0 |      0 |       0 |          0 |         0 |
| basket_5 |   200 |       0 |      0 |       0 |          0 |         0 |
| basket_6 |   100 |       0 |    200 |       0 |          0 |         0 |

Species x assemblages matrix based on the **fruits** dataset

  

- a **data frame** summarizing traits types called `fruits_traits_cat`
  in this tutorial (for details, see [mFD General
  Workflow](https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html)):

``` r
data("fruits_traits_cat", package = "mFD")
knitr::kable(head(fruits_traits_cat), 
             caption = "Traits types based on **fruits & baskets** dataset")
```

| trait_name | trait_type | fuzzy_name |
|:-----------|:-----------|:-----------|
| Size       | O          | NA         |
| Plant      | N          | NA         |
| Climate    | O          | NA         |
| Seed       | O          | NA         |
| Sugar      | Q          | NA         |
| Use.raw    | F          | Use        |

Traits types based on **fruits & baskets** dataset

  

Then, we can compute functional distance using the
[`mFD::funct.dist()`](https://cmlmagneville.github.io/mFD/reference/funct.dist.md)
function as follows:

**USAGE**

    ## [1] "Running w.type=equal on groups=c(Size)"
    ## [1] "Running w.type=equal on groups=c(Plant)"
    ## [1] "Running w.type=equal on groups=c(Climate)"
    ## [1] "Running w.type=equal on groups=c(Seed)"
    ## [1] "Running w.type=equal on groups=c(Sugar)"
    ## [1] "Running w.type=equal on groups=c(Use,Use,Use)"

  

## Generalisation of Hill numbers for alpha functional diversity

  

The family of indices presented in [Chao *et al.*
(2019)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343)
allows computing FD based on pairwise distance between species and their
weights in assemblages. This generalization of Hill numbers framework is
based on two parameters:

- `q`: the importance of species weight compared to species distance.
  Values allowed in `mFD` are 0, 1, 2 (the most often used, see below).

- `tau`: the threshold of functional distinctiveness between any two
  species (*i.e.* all species with distance above this threshold are
  considered as functionally equally distinct). Values allowed in `mFD`
  are “min(imum)”, “mean(imum)” and “max(imum)”.

  

Indices are expressed as effective number of functionally equally
distinct species (or *virtual functional groups*) and could thus be
directly compared to taxonomic Hill numbers (including species
richness).

  

**NOTE** For more details about the properties of Hill numbers FD read
[Chao *et al.*
(2019)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343)
and especially its Figures 1 & 2.

  

All these indices can be computed with the function
[`mFD::alpha.fd.hill()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.hill.md).

  

Here we start by comparing the **‘classical’ Rao’s quadratic entropy
expressed in Hill numbers** following [Ricotta & Szeidl
(2009)](https://www.sciencedirect.com/science/article/pii/S0040580909001117)
which is the special case with `q = 2` and `tau = "max"`.

  

**USAGE**

``` r
baskets_FD2max <- mFD::alpha.fd.hill(
  asb_sp_w = baskets_fruits_weights, 
  sp_dist  = fruits_gower, 
  tau      = "max", 
  q        = 2)
```

  

Then, we can compute **Hill numbers FD of order 2** computed with
`tau = "mean"` and `q = 2` as recommended in [Chao *et al.*
(2019)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343)

  

**USAGE**

``` r
baskets_FD2mean <- mFD::alpha.fd.hill(
  asb_sp_w = baskets_fruits_weights, 
  sp_dist  = fruits_gower, 
  tau      = "mean", 
  q        = 2)
```

  

We can now compare these two metrics:

``` r
round(cbind(FD2max  = baskets_FD2max$"asb_FD_Hill"[ , 1], 
            FD2mean = baskets_FD2mean$"asb_FD_Hill"[ , 1]), 2)
```

    ##           FD2max FD2mean
    ## basket_1    1.50    2.62
    ## basket_2    1.83    3.97
    ## basket_3    1.86    4.10
    ## basket_4    1.27    1.72
    ## basket_5    1.30    1.85
    ## basket_6    1.74    3.73
    ## basket_7    1.82    3.94
    ## basket_8    1.40    2.16
    ## basket_9    1.53    2.75
    ## basket_10   1.53    3.01

  

We can see that FD computed with `tau = "max"` is less variable (ranging
from 1.50 to only 1.86) than FD computed with `tau = "min"` (ranging
from 1.72 to 4.10) illustrating its **higher sensitivity to functional
differences between species**.

  

**NB** **Note that even with** `q = 0`**, weights of species are still
accounted for by FD.** Hence, if the goal is to compute a richness-like
index (*i.e.* accounting only for distance between species present),
function
[`mFD::alpha.fd.hill()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.hill.md)
should be applied to **species occurrence data** (coded as 0/1,
previously computed using sp.tr.summary) so that all species have the
same weight). Species occurrence data can be retrieve through the
[`mFD::asb.sp.summary()`](https://cmlmagneville.github.io/mFD/reference/asb.sp.summary.md)
function:

  

**USAGE**

``` r
# Retrieve species occurrences data:
baskets_summary    <- mFD::asb.sp.summary(baskets_fruits_weights)
baskets_fruits_occ <- baskets_summary$"asb_sp_occ"

head(baskets_fruits_occ)
```

    ##          apple apricot banana currant blackberry blueberry cherry grape
    ## basket_1     1       0      1       0          0         0      1     0
    ## basket_2     1       0      1       0          0         0      1     0
    ## basket_3     1       0      1       0          0         0      1     0
    ## basket_4     1       0      0       0          0         0      0     0
    ## basket_5     1       0      0       0          0         0      0     0
    ## basket_6     1       0      1       0          0         0      0     0
    ##          grapefruit kiwifruit lemon lime litchi mango melon orange
    ## basket_1          0         0     1    0      0     0     1      0
    ## basket_2          0         0     1    0      0     0     1      0
    ## basket_3          0         0     1    0      0     0     1      0
    ## basket_4          0         1     1    0      0     0     0      1
    ## basket_5          0         1     1    0      0     0     0      1
    ## basket_6          0         0     0    1      1     1     0      1
    ##          passion_fruit peach pear pineapple plum raspberry strawberry tangerine
    ## basket_1             1     0    1         0    0         0          1         0
    ## basket_2             1     0    1         0    0         0          1         0
    ## basket_3             1     0    1         0    0         0          1         0
    ## basket_4             0     1    1         0    1         0          0         1
    ## basket_5             0     1    1         0    1         0          0         1
    ## basket_6             0     0    0         1    0         0          0         0
    ##          water_melon
    ## basket_1           0
    ## basket_2           0
    ## basket_3           0
    ## basket_4           0
    ## basket_5           0
    ## basket_6           1

``` r
# Compute alpha FD Hill with q = 0:
baskets_FD0mean <- mFD::alpha.fd.hill(
  asb_sp_w = baskets_fruits_occ, 
  sp_dist  = fruits_gower, 
  tau      = "mean", 
  q        = 0)

round(baskets_FD0mean$"asb_FD_Hill", 2)
```

    ##           FD_q0
    ## basket_1   4.73
    ## basket_2   4.73
    ## basket_3   4.73
    ## basket_4   1.93
    ## basket_5   1.93
    ## basket_6   4.57
    ## basket_7   4.57
    ## basket_8   3.67
    ## basket_9   3.67
    ## basket_10  3.52

  

We can see that baskets with same composition of fruits species have
same FD values (*e.g* *basket_1*, *basket_2* and *basket_3*)

  

## Generalisation of Hill numbers for beta functional diversity

  

Framework of [Chao *et al.*
(2019)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343)
also allows computing beta-diversity, with 2 framework similar to
Jaccard and Sorensen ones for taxonomic diversity. The
`mFD:beta.fd.hill()` function computes these indices.

  

**NB** Note that **total weight of assemblages is affecting computation
of functional beta-diversity**. Hence **if it is does not reflect an
ecological pattern** (*e.g*. rather difference in sampling effort), it
is recommended to **apply**
[`mFD::beta.fd.hill()`](https://cmlmagneville.github.io/mFD/reference/beta.fd.hill.md)
**to relative weight of species in assemblages**.

  

``` r
# retrieve total weight per basket:
baskets_summary$"asb_tot_w"
```

    ##  basket_1  basket_2  basket_3  basket_4  basket_5  basket_6  basket_7  basket_8 
    ##      2000      2000      2000      2000      2000      2000      2000      2000 
    ##  basket_9 basket_10 
    ##      2000      2000

``` r
# Here baskets all contain 2000g of fruits, we illustrate how to compute...
# relative weights using the output of asb.sp.summary:

baskets_fruits_relw <- baskets_fruits_weights / baskets_summary$"asb_tot_w"
apply(baskets_fruits_relw, 1, sum)
```

    ##  basket_1  basket_2  basket_3  basket_4  basket_5  basket_6  basket_7  basket_8 
    ##         1         1         1         1         1         1         1         1 
    ##  basket_9 basket_10 
    ##         1         1

  

Now we can compute functional beta-diversity of order `q = 2` (with
`tau = "mean"` for higher sensitivity) with Jaccard-type index:

  

**USAGE**

``` r
# Compute index:
baskets_betaq2 <- mFD::beta.fd.hill(
  asb_sp_w  = baskets_fruits_relw, 
  sp_dist   = fruits_gower, 
  q         = 2,
  tau       = "mean", 
  beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result
mFD::dist.to.df(list_dist = list("FDq2" = baskets_betaq2$"beta_fd_q"$"q2"))
```

    ## Registered S3 method overwritten by 'dendextend':
    ##   method     from 
    ##   rev.hclust vegan

    ##          x1        x2        FDq2
    ## 1  basket_1  basket_2 0.058982325
    ## 2  basket_1  basket_3 0.078716397
    ## 3  basket_1  basket_4 0.029573623
    ## 4  basket_1  basket_5 0.027059789
    ## 5  basket_1  basket_6 0.484115290
    ## 6  basket_1  basket_7 0.292594562
    ## 7  basket_1  basket_8 0.290545545
    ## 8  basket_1  basket_9 0.185475113
    ## 9  basket_1 basket_10 0.011136995
    ## 10 basket_2  basket_3 0.004420448
    ## 11 basket_2  basket_4 0.161833512
    ## 12 basket_2  basket_5 0.162571972
    ## 13 basket_2  basket_6 0.260541701
    ## 14 basket_2  basket_7 0.097053161
    ## 15 basket_2  basket_8 0.294504888
    ## 16 basket_2  basket_9 0.225897615
    ## 17 basket_2 basket_10 0.058298760
    ## 18 basket_3  basket_4 0.172877455
    ## 19 basket_3  basket_5 0.178123024
    ## 20 basket_3  basket_6 0.207102590
    ## 21 basket_3  basket_7 0.081951839
    ## 22 basket_3  basket_8 0.308336365
    ## 23 basket_3  basket_9 0.241482168
    ## 24 basket_3 basket_10 0.082649928
    ## 25 basket_4  basket_5 0.001049851
    ## 26 basket_4  basket_6 0.511165067
    ## 27 basket_4  basket_7 0.421141181
    ## 28 basket_4  basket_8 0.482330219
    ## 29 basket_4  basket_9 0.342926459
    ## 30 basket_4 basket_10 0.050817451
    ## 31 basket_5  basket_6 0.532084800
    ## 32 basket_5  basket_7 0.438841544
    ## 33 basket_5  basket_8 0.496554512
    ## 34 basket_5  basket_9 0.336894275
    ## 35 basket_5 basket_10 0.044052657
    ## 36 basket_6  basket_7 0.068382884
    ## 37 basket_6  basket_8 0.759422492
    ## 38 basket_6  basket_9 0.680332414
    ## 39 basket_6 basket_10 0.453325136
    ## 40 basket_7  basket_8 0.478528108
    ## 41 basket_7  basket_9 0.431531941
    ## 42 basket_7 basket_10 0.265928889
    ## 43 basket_8  basket_9 0.020812705
    ## 44 basket_8 basket_10 0.345652088
    ## 45 basket_9 basket_10 0.219780509

We can see that *basket 1* is similar (beta \< 0.1) to *baskets
2,3,4,5,10* and that it is the most dissimilar to *basket 8 (beta \>
0.5)*. *Baskets 4 and 5* are highly dissimilar (beta \> 0.8) to *basket
8*.

  

## References

  

- Chao *et al.* (2019) An attribute diversity approach to functional
  diversity, functional beta diversity, and related (dis)similarity
  measures. *Ecological Monographs*, **89**, e01343.

- Ricotta & Szeidl (2009) Diversity partitioning of Rao’s quadratic
  entropy. *Theoretical Population Biology*, **76**, 299-302.
