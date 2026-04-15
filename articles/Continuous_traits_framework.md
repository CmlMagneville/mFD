# Continuous Traits Framework

## About this tutorial

  

This tutorial explains the workflow used to compute functional space
based on continuous traits and it shows how to retrieve species
coordinates and species functional distances in the functional space.

  

**DATA** This tutorial uses a dataset from one of the 80 CESTES database
Jeliazkov & the CESTES consortium (2019)) based on \[Villeger *et al.*
2012\]
(<https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0040679>).
This data frame contains 45 fish species from the Terminos Lagoon (Gulf
of Mexico) gathered into 36 sites considered as assemblages. Each
species is described with 16 continuous morphological traits.

  

When the dataset **only** gathers continuous traits, the functional
space can be computed using one trait for one dimension or using
Principal Component Analysis (PCA: convert correlations among samples
into a 2D plot). **NB** Using a PCoA on continuous traits and euclidean
distance is the same than using a PCA (clusters made by minimizing the
linear distance (PCoA) are the same as those obtained by maximizing
linear correlations (PCA)).

  

## 1. Load dataset

  

The species traits data frame has rows corresponding to species and
columns corresponding to traits. The different traits are summed up in
the following table:

  

| Trait name |       Trait signification        |
|:----------:|:--------------------------------:|
|    logM    |            log(mass)             |
|    Ogsf    |        Oral gape surface         |
|    OgSh    |         Oral gape shape          |
|    OgPo    |        Oral gape position        |
|    GrLg    |        Gill raker length         |
|    GtLg    |            Gut length            |
|    EySz    |             Eye size             |
|    EyPo    |           Eye position           |
|    BdSh    |      Body transversal shape      |
|    BdSf    |     Body transversal surface     |
|    PfPo    |      Pectoral fin position       |
|    PfSh    | Aspect ratio of the pectoral fin |
|    CpHt    |    Caudal peduncle throttling    |
|    CfSh    |  Aspect ratio of the caudal fin  |
|    FsRt    |        Fins surface ratio        |
|    FsSf    | Fins surface to body size ratio  |

  

To work with `mFD` with only continuous traits, you must load two
objects:

  

- `sp_tr`: species x traits data frame

``` r
# load dataset:
sp_tr <- read.csv(system.file("extdata", "data_cestes_sp_tr.csv", 
                              package = "mFD"), dec = ",", sep = ":")

rownames(sp_tr) <- sp_tr$"Sp"
sp_tr <- sp_tr[ , -1]

# display the table:
knitr::kable(head(sp_tr), 
             caption = "Species x Traits data frame based on *CESTES* dataset")
```

|                             |  logM |  OgSf |  OgSh |  OgPo |  EySz |  GrLg |  GtLg |  EyPo |  BdSh |  BdSf |  PfPo |  PfSh |  CpHt |  CfSh |  FsRt |  FsSf |
|:----------------------------|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| Achirus_lineatus            | 2.187 | 0.072 | 0.947 | 1.000 | 0.151 | 0.000 | 1.782 | 1.000 | 0.143 | 2.168 | 0.000 | 0.000 | 1.121 | 0.767 | 0.000 | 1.158 |
| Anchoa_mitchilli            | 0.706 | 0.283 | 2.054 | 0.508 | 0.474 | 0.381 | 0.688 | 0.584 | 3.292 | 3.974 | 0.773 | 2.708 | 2.493 | 3.108 | 0.504 | 2.618 |
| Archosargus_probatocephalus | 2.674 | 0.082 | 0.754 | 0.221 | 0.282 | 0.035 | 2.591 | 0.652 | 3.091 | 1.679 | 0.666 | 4.740 | 2.584 | 2.393 | 1.526 | 1.504 |
| Archosargus_rhomboidalis    | 3.327 | 0.056 | 0.648 | 0.273 | 0.294 | 0.032 | 3.467 | 0.647 | 3.102 | 1.596 | 0.638 | 6.868 | 2.906 | 2.684 | 1.060 | 1.502 |
| Ariopsis_felis              | 3.110 | 0.173 | 0.513 | 0.346 | 0.263 | 0.128 | 2.078 | 0.647 | 1.021 | 1.658 | 0.806 | 3.480 | 3.592 | 4.052 | 0.781 | 1.648 |
| Bagre_marinus               | 2.170 | 0.248 | 0.519 | 0.489 | 0.357 | 0.142 | 2.154 | 0.613 | 1.016 | 2.087 | 0.662 | 3.674 | 4.205 | 3.460 | 0.625 | 1.678 |

Species x Traits data frame based on *CESTES* dataset

  

- `asb_sp_w`: species x assemblages data frame summarizing biomass
  recorded in a volume of 4500m³ per site and per species:

``` r
# load dataset:
asb_sp_w <- read.csv(system.file("extdata", "data_cestes_asb_sp_w.csv", 
                                 package = "mFD"), dec = ",", sep = ":")

rownames(asb_sp_w) <- paste0("site", sep = "_", asb_sp_w$Sites)
asb_sp_w <- asb_sp_w[ , -1]

asb_sp_w$Urobatis_jamaicensis <- as.numeric(asb_sp_w$Urobatis_jamaicensis)

# remove sites 12, 23, 35 because FRic can not be computed on it...
# ... (for a clean example):
asb_sp_w <- asb_sp_w[-c(11, 22, 33), ]

# display the table:
knitr::kable(asb_sp_w[1:7, 1:6], 
             caption = "Species x Assemblages data frame based on *CESTES* dataset for the first six species and first seven sites")
```

|        | Achirus_lineatus | Anchoa_mitchilli | Archosargus_probatocephalus | Archosargus_rhomboidalis | Ariopsis_felis | Bagre_marinus |
|:-------|-----------------:|-----------------:|----------------------------:|-------------------------:|---------------:|--------------:|
| site_1 |                0 |                0 |                           0 |                        0 |          169.8 |          66.5 |
| site_2 |                0 |                0 |                           0 |                        0 |            0.0 |          29.5 |
| site_3 |                0 |                0 |                           0 |                        0 |          592.4 |           0.0 |
| site_5 |                0 |                0 |                           0 |                        0 |            0.0 |           0.0 |
| site_6 |                0 |                0 |                           0 |                        0 |            0.0 |           0.0 |
| site_7 |                0 |                0 |                           0 |                        0 |          135.4 |           0.0 |
| site_8 |                0 |                0 |                           0 |                        0 |            0.0 |           0.0 |

Species x Assemblages data frame based on *CESTES* dataset for the first
six species and first seven sites

  

## 2. Compute the functional space

  

Based on the species-trait data frame or the species-standardized traits
data frame, `mFD` allows to build a functional space based on a PCA or
using each trait as a dimension. (**NB** Using up to the 1.0.3 version
of the `mFD` package does not allow weighting continuous traits, it will
be done in a next version of the package. You can use the
[`col.w`](https://rdrr.io/cran/FactoMineR/man/PCA.html) argument of the
PCA function of the FactomineR package.). The function used to compute
functional space with continuous traits is called
[`mFD::tr.cont.fspace()`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
and is used as follow:

  

**USAGE**

``` r
mFD::tr.cont.fspace(
  sp_tr        = sp_tr, 
  pca          = TRUE, 
  nb_dim       = 7, 
  scaling      = "scale_center",
  compute_corr = "pearson")
```

  

It takes as inputs:

- the *sp_tr* data frame

- a *pca* argument that must be set to TRUE if you want to compute a PCA
  or to FALSE if you want to use each trait as a dimension to construct
  the multidimensional space

- a *nb_dim* argument referring to the maximum number of dimensions for
  multidimensional functional spaces. Final number of dimensions depends
  on the number of positive eigenvalues obtained with PCA if *pca =
  TRUE* or the number of traits used if *pca = FALSE*. **NB** High value
  for *nb_dim* can increase computation time.

- a *scaling* argument allowing traits values to be standardized. They
  can be standardized in several ways: standardization by the range
  value of the trait, center-transformation, scale transformation or
  scale-center transformation can be used. You can also chose not to
  standardize traits values. **NOTE** Scaling ensures that trait-based
  distances and distances in the functional space have the same maximum.
  Scaling distances implies that the quality of the functional space
  accounts for congruence in distances rather than their equality

- a *compute_corr* argument which refers to a string value to compute
  Pearson correlation coefficients between traits using *“pearson”* or
  not using *“none”*.

  

In this example, we will compute a PCA based on a maximum number of 7
dimensions and get Pearson’s correlation coefficients:

  

``` r
fspace <- mFD::tr.cont.fspace(
  sp_tr        = sp_tr, 
  pca          = TRUE, 
  nb_dim       = 10, 
  scaling      = "scale_center",
  compute_corr = "pearson")
```

  

If the PCA is computed, the output contains:

- quality metrics for spaces from 2 dimensions to *nb_dim* dimensions:

**NB** **mean absolute deviation (mad)** reflects the actual magnitude
of errors that affect distances, hence FD metrics ; **mean squared
deviation (msd)** reflects the potential risk associated with a few
species pairs being strongly misplaced in the functional space ([Maire
*et al.*
(2015)](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12299)).

  

``` r
fspace$"quality_metrics"
```

    ##            mAD         mSD
    ## 2D  1.83867581 3.380728750
    ## 3D  1.22034177 1.489235074
    ## 4D  0.90241642 0.814355528
    ## 5D  0.59513109 0.354181166
    ## 6D  0.46267807 0.214071066
    ## 7D  0.33359921 0.111288894
    ## 8D  0.23662790 0.055993113
    ## 9D  0.16252824 0.026415439
    ## 10D 0.09772059 0.009549359

  

**NB** The lower the quality metric is, the better the quality of your
space is. Here, thanks to mAD and mSD value, we can see that as the
number of dimensions increases, the quality increases. However, to
decrease computation time, we can chose to work with the 6D space which
has good quality of functional space. Generally, you must keep in mind a
trade-off between the number of axes and quality of functional space.
Increasing the number of functional axes increases computation time.

  

- eigenvalues, percentage of variance explained and cumulative
  percentage of variance explained for each axis up to *nb_dim*
  dimensions:

  

``` r
fspace$"eigenvalues_percentage_var"
```

    ##      eigenvalue percentage of variance cumulative percentage of variance
    ## PC1   5.0894430              32.531951                          32.53195
    ## PC2   2.3267315              14.872574                          47.40452
    ## PC3   1.9839001              12.681179                          60.08570
    ## PC4   1.6167089              10.334077                          70.41978
    ## PC5   1.4277215               9.126061                          79.54584
    ## PC6   0.7557779               4.830967                          84.37681
    ## PC7   0.6704457               4.285519                          88.66233
    ## PC8   0.5134516               3.282006                          91.94433
    ## PC9   0.3510937               2.244207                          94.18854
    ## PC10  0.2987544               1.909652                          96.09819

  

- a matrix of species coordinates in the functional space:

  

``` r
head(fspace$"sp_faxes_coord")
```

    ##                                    PC1         PC2        PC3        PC4
    ## Achirus_lineatus            -3.1102450 -4.18756669  0.8212298 -0.7162535
    ## Anchoa_mitchilli             3.5220130  0.04037005  2.7468533  0.9280388
    ## Archosargus_probatocephalus  0.5598244  0.56101682 -1.9905123 -0.2138744
    ## Archosargus_rhomboidalis     0.7557227  1.05395885 -2.3524317 -1.0131106
    ## Ariopsis_felis               0.7945539  0.63147268 -1.9667765 -0.6958125
    ## Bagre_marinus                1.1521566  0.25626547 -0.7971292 -0.6188243
    ##                                    PC5        PC6         PC7         PC8
    ## Achirus_lineatus             0.1119252  0.8809080  0.22663165  0.52316085
    ## Anchoa_mitchilli             1.7894441  2.0677051 -0.85219538  0.56621740
    ## Archosargus_probatocephalus -0.5774869 -0.4797906 -0.33549421 -0.31168020
    ## Archosargus_rhomboidalis    -1.2508715 -0.4855985  0.02667960 -0.01979034
    ## Ariopsis_felis               0.5466197 -0.5300129 -0.03313655  1.72467762
    ## Bagre_marinus                0.8995832 -0.3424424 -0.29934540  2.00458491
    ##                                      PC9       PC10
    ## Achirus_lineatus            -0.211572742  0.5211315
    ## Anchoa_mitchilli             0.804586482 -0.8717554
    ## Archosargus_probatocephalus  0.068855803  0.5161352
    ## Archosargus_rhomboidalis     0.089848798  1.3715056
    ## Ariopsis_felis              -0.450194649 -0.5236044
    ## Bagre_marinus                0.007347954  0.2391765

  

- a dist object containing species euclidean distances in the functional
  space (here 5D space and for the first five species):

  

``` r
dist_mat <- as.matrix(fspace$sp_dist_multidim$"6D")
dist_mat[1:5, 1:5]
```

    ##                             Achirus_lineatus Anchoa_mitchilli
    ## Achirus_lineatus                    0.000000         8.514492
    ## Anchoa_mitchilli                    8.514492         0.000000
    ## Archosargus_probatocephalus         6.819349         6.699577
    ## Archosargus_rhomboidalis            7.503606         7.362973
    ## Ariopsis_felis                      6.958710         6.398155
    ##                             Archosargus_probatocephalus
    ## Achirus_lineatus                               6.819349
    ## Anchoa_mitchilli                               6.699577
    ## Archosargus_probatocephalus                    0.000000
    ## Archosargus_rhomboidalis                       1.226627
    ## Ariopsis_felis                                 1.248610
    ##                             Archosargus_rhomboidalis Ariopsis_felis
    ## Achirus_lineatus                            7.503606       6.958710
    ## Anchoa_mitchilli                            7.362973       6.398155
    ## Archosargus_probatocephalus                 1.226627       1.248610
    ## Archosargus_rhomboidalis                    0.000000       1.913729
    ## Ariopsis_felis                              1.913729       0.000000

  

- a dist object containing species distances based on traits (here for
  the first five species):

  

``` r
dist_mat <- as.matrix(fspace$sp_dist_init)
dist_mat[1:5, 1:5]
```

    ##                             Achirus_lineatus Anchoa_mitchilli
    ## Achirus_lineatus                    0.000000         8.810865
    ## Anchoa_mitchilli                    8.810865         0.000000
    ## Archosargus_probatocephalus         6.954479         7.010258
    ## Archosargus_rhomboidalis            7.643421         7.834878
    ## Ariopsis_felis                      7.165968         6.725278
    ##                             Archosargus_probatocephalus
    ## Achirus_lineatus                               6.954479
    ## Anchoa_mitchilli                               7.010258
    ## Archosargus_probatocephalus                    0.000000
    ## Archosargus_rhomboidalis                       1.814132
    ## Ariopsis_felis                                 2.780871
    ##                             Archosargus_rhomboidalis Ariopsis_felis
    ## Achirus_lineatus                            7.643421       7.165968
    ## Anchoa_mitchilli                            7.834878       6.725278
    ## Archosargus_probatocephalus                 1.814132       2.780871
    ## Archosargus_rhomboidalis                    0.000000       3.371347
    ## Ariopsis_felis                              3.371347       0.000000

  

- a correlation matrix containing correlation between traits and their
  associated pvalue:

  

``` r
fspace$"tr_correl"
```

    ##       logM  OgSf  OgSh  OgPo  EySz  GrLg  GtLg  EyPo  BdSh  BdSf  PfPo  PfSh
    ## logM  1.00 -0.05 -0.62 -0.41  0.07 -0.27  0.22  0.48 -0.48 -0.69 -0.35 -0.27
    ## OgSf -0.05  1.00  0.19  0.04  0.14  0.28 -0.22  0.03 -0.18 -0.05  0.17 -0.07
    ## OgSh -0.62  0.19  1.00  0.14  0.06  0.35 -0.19 -0.48  0.75  0.57  0.31  0.34
    ## OgPo -0.41  0.04  0.14  1.00 -0.15 -0.06 -0.11  0.22 -0.05  0.16 -0.30 -0.37
    ## EySz  0.07  0.14  0.06 -0.15  1.00  0.34 -0.21 -0.14 -0.21 -0.15  0.12  0.41
    ## GrLg -0.27  0.28  0.35 -0.06  0.34  1.00  0.13 -0.38  0.10  0.16  0.50  0.21
    ## GtLg  0.22 -0.22 -0.19 -0.11 -0.21  0.13  1.00 -0.18 -0.02 -0.21  0.24 -0.15
    ## EyPo  0.48  0.03 -0.48  0.22 -0.14 -0.38 -0.18  1.00 -0.62 -0.31 -0.86 -0.65
    ## BdSh -0.48 -0.18  0.75 -0.05 -0.21  0.10 -0.02 -0.62  1.00  0.50  0.35  0.47
    ## BdSf -0.69 -0.05  0.57  0.16 -0.15  0.16 -0.21 -0.31  0.50  1.00  0.16  0.05
    ## PfPo -0.35  0.17  0.31 -0.30  0.12  0.50  0.24 -0.86  0.35  0.16  1.00  0.50
    ## PfSh -0.27 -0.07  0.34 -0.37  0.41  0.21 -0.15 -0.65  0.47  0.05  0.50  1.00
    ## CpHt -0.46 -0.07  0.50 -0.05 -0.36  0.12  0.09 -0.61  0.75  0.48  0.41  0.35
    ## CfSh -0.44 -0.05  0.37 -0.13 -0.09  0.27  0.28 -0.73  0.53  0.35  0.61  0.40
    ## FsRt  0.08  0.25 -0.21 -0.30 -0.42 -0.06  0.11  0.00 -0.05  0.13  0.15 -0.17
    ## FsSf  0.32  0.19 -0.13 -0.58  0.29  0.00 -0.28  0.30 -0.19 -0.01 -0.25  0.15
    ##       CpHt  CfSh  FsRt  FsSf
    ## logM -0.46 -0.44  0.08  0.32
    ## OgSf -0.07 -0.05  0.25  0.19
    ## OgSh  0.50  0.37 -0.21 -0.13
    ## OgPo -0.05 -0.13 -0.30 -0.58
    ## EySz -0.36 -0.09 -0.42  0.29
    ## GrLg  0.12  0.27 -0.06  0.00
    ## GtLg  0.09  0.28  0.11 -0.28
    ## EyPo -0.61 -0.73  0.00  0.30
    ## BdSh  0.75  0.53 -0.05 -0.19
    ## BdSf  0.48  0.35  0.13 -0.01
    ## PfPo  0.41  0.61  0.15 -0.25
    ## PfSh  0.35  0.40 -0.17  0.15
    ## CpHt  1.00  0.81  0.17 -0.22
    ## CfSh  0.81  1.00 -0.03 -0.25
    ## FsRt  0.17 -0.03  1.00  0.17
    ## FsSf -0.22 -0.25  0.17  1.00
    ## 
    ## n= 45 
    ## 
    ## 
    ## P
    ##      logM   OgSf   OgSh   OgPo   EySz   GrLg   GtLg   EyPo   BdSh   BdSf  
    ## logM        0.7651 0.0000 0.0057 0.6264 0.0738 0.1394 0.0010 0.0009 0.0000
    ## OgSf 0.7651        0.2153 0.7759 0.3746 0.0632 0.1512 0.8559 0.2328 0.7626
    ## OgSh 0.0000 0.2153        0.3540 0.7084 0.0196 0.2164 0.0009 0.0000 0.0000
    ## OgPo 0.0057 0.7759 0.3540        0.3200 0.6850 0.4850 0.1492 0.7320 0.2970
    ## EySz 0.6264 0.3746 0.7084 0.3200        0.0217 0.1688 0.3552 0.1752 0.3416
    ## GrLg 0.0738 0.0632 0.0196 0.6850 0.0217        0.4014 0.0092 0.5024 0.3083
    ## GtLg 0.1394 0.1512 0.2164 0.4850 0.1688 0.4014        0.2431 0.8801 0.1633
    ## EyPo 0.0010 0.8559 0.0009 0.1492 0.3552 0.0092 0.2431        0.0000 0.0365
    ## BdSh 0.0009 0.2328 0.0000 0.7320 0.1752 0.5024 0.8801 0.0000        0.0005
    ## BdSf 0.0000 0.7626 0.0000 0.2970 0.3416 0.3083 0.1633 0.0365 0.0005       
    ## PfPo 0.0198 0.2608 0.0364 0.0465 0.4391 0.0005 0.1170 0.0000 0.0179 0.2825
    ## PfSh 0.0713 0.6297 0.0237 0.0129 0.0049 0.1765 0.3167 0.0000 0.0010 0.7222
    ## CpHt 0.0017 0.6657 0.0005 0.7423 0.0137 0.4351 0.5597 0.0000 0.0000 0.0008
    ## CfSh 0.0027 0.7412 0.0126 0.4075 0.5349 0.0752 0.0648 0.0000 0.0002 0.0168
    ## FsRt 0.6006 0.1010 0.1606 0.0483 0.0043 0.6867 0.4784 0.9796 0.7473 0.4052
    ## FsSf 0.0344 0.2210 0.3825 0.0000 0.0554 0.9783 0.0604 0.0428 0.2024 0.9353
    ##      PfPo   PfSh   CpHt   CfSh   FsRt   FsSf  
    ## logM 0.0198 0.0713 0.0017 0.0027 0.6006 0.0344
    ## OgSf 0.2608 0.6297 0.6657 0.7412 0.1010 0.2210
    ## OgSh 0.0364 0.0237 0.0005 0.0126 0.1606 0.3825
    ## OgPo 0.0465 0.0129 0.7423 0.4075 0.0483 0.0000
    ## EySz 0.4391 0.0049 0.0137 0.5349 0.0043 0.0554
    ## GrLg 0.0005 0.1765 0.4351 0.0752 0.6867 0.9783
    ## GtLg 0.1170 0.3167 0.5597 0.0648 0.4784 0.0604
    ## EyPo 0.0000 0.0000 0.0000 0.0000 0.9796 0.0428
    ## BdSh 0.0179 0.0010 0.0000 0.0002 0.7473 0.2024
    ## BdSf 0.2825 0.7222 0.0008 0.0168 0.4052 0.9353
    ## PfPo        0.0005 0.0047 0.0000 0.3207 0.1029
    ## PfSh 0.0005        0.0181 0.0070 0.2514 0.3267
    ## CpHt 0.0047 0.0181        0.0000 0.2722 0.1433
    ## CfSh 0.0000 0.0070 0.0000        0.8596 0.0914
    ## FsRt 0.3207 0.2514 0.2722 0.8596        0.2563
    ## FsSf 0.1029 0.3267 0.1433 0.0914 0.2563

  

Here we can notice that there is no strong correlation between traits.
**NB** However, if some strong correlation is to be found, then one of
the two correlated trait can be remove from the analysis.

If the PCA is not computed, outputs are the same except that mad and msd
are not computed and that only one distance object is returned.

  

## 3. Plot functional space, compute and illustrate indices

  

Then, based on the species coordinates matrix, steps are similar as
those listed in the [mFD General
Workflow](https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html),
from step 5 till the end.

  

## References

  

- Maire *et al.* (2015) How many dimensions are needed to accurately
  assess functional diversity? A pragmatic approach for assessing the
  quality of functional spaces. *Global Ecology and Biogeography*,
  **24**, 728-740.

- Villeger *et al.* (2012) Low Functional beta Diversity Despite High
  Taxonomic beta Diversity among Tropical Estuarine Fish Communities.
  *PLoS ONE*, **7**, e40679.
