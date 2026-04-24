# Compute functional spaces and their quality

Compute a Principal Coordinates Analysis (PCoA) using functional
distance between species. Then the function evaluates the quality of
spaces built using an increasing number of principal components. Quality
is evaluated as the (absolute or squared) deviation between trait-based
distance (input) and distance in the PCoA-based space (raw Euclidean
distance or scaled distance according to its maximum value and maximum
of trait-based distance). Option to compute a functional dendrogram and
its quality. This function is based on the framework presented in Maire
*et al.* (2015).

## Usage

``` r
quality.fspaces(
  sp_dist,
  fdendro = NULL,
  maxdim_pcoa = 10,
  deviation_weighting = "absolute",
  fdist_scaling = FALSE
)
```

## Arguments

- sp_dist:

  a dist object with pairwise distance among all species (at least 3
  species needed). Functional distance matrix from trait values can be
  computed using
  [`funct.dist`](https://cmlmagneville.github.io/mFD/reference/funct.dist.md)
  function.

- fdendro:

  a character string indicating the clustering algorithm to use to
  compute dendrogram. Should be one of the method recognized by
  [`hclust`](https://rdrr.io/r/stats/hclust.html) (e.g. 'average' for
  UPGMA). Default: `fdendro = NULL` (so no dendrogram computed).

- maxdim_pcoa:

  a single numeric value with maximum number of PCoA axes to consider to
  build multidimensional functional spaces. Default: `maxdim_pcoa = 10`.
  See below about number of axes actually considered.

- deviation_weighting:

  a character string referring to the method(s) used to weight the
  differences between species pairwise distance in the functional space
  and trait-based distance. `'absolute'` (default) means absolute
  differences are used to compute mean absolute deviation *mad* index;
  `'squared'` means squared differences are used to compute root of mean
  squared deviation *rmsd* index. Both values could be provided to
  compare quality metrics.

- fdist_scaling:

  a vector with logical value(s) specifying whether distances in the
  functional space should be scaled before computing differences with
  trait-based distances. Scaling ensures that trait-based distances and
  distances in the functional space have the same maximum. Default:
  `FALSE`. Both values could be provided to compare quality metrics.

## Value

A list with:

- `$quality_fspaces`: a data frame with quality metric(s) for each
  functional space. Functional spaces are named as 'pcoa\_.d' and if
  required 'tree_clustering method'. Quality metrics are named after
  deviation_weighting ('mad' for 'absolute' and and 'rmsd' for
  'squared') and if fdist_scaling is TRUE with suffix '\_scaled'.

- `$details_trdist` a list with 2 elements: `$trdist_summary` a vector
  with minimum (min), maximum (max), mean (mean) and standard deviation
  (sd) of `sp_dist`; `$trdist_euclidean` a logical value indicating
  whether `sp_dist` checks Euclidean properties.

- `$details_fspaces` a list with 4 elements: `$sp_pc_coord` a matrix
  with coordinates of species (rows) along Principal Components
  (columns) with positive eigenvalues ; `$pc_eigenvalues` a matrix with
  eigenvalues of axes from PCoA ; `$dendro` a hclust object with the
  dendrogram details (null if no dendrogram computed) ;
  `$pairsp_fspaces_dist` a dataframe containing for each pair of species
  (rows), their names in the 2 first columns ('sp.x' and 'sp.y'), their
  distance based on trait-values ('tr'), and their Euclidean (for PCoA)
  or cophenetic (for dendrogram if computed) distance in each of the
  functional space computed ('pcoa_1d', 'pcoa_2d', ... , 'tree_clust');
  if `fdist_scaling = TRUE`, `$pairsp_fspaces_dist_scaled` a data frame
  with scaled values of distances in functional spaces.

- `$details_deviation` a list of data frames: `$dev_distsp` a dataframe
  containing for each space (columns) the difference for all species
  pairs (rows) of the distance in the functional space and trait-based
  distance (i.e. positive deviation indicates overestimation of actual
  distance) ; `$abs_dev_distsp` and/or `$sqr_dev_distsp`, dataframes
  with for each space (columns) and all species pairs (rows) the
  absolute or squared deviation of distance ; if `fdist_scaling = TRUE`
  `$dev_distsp_scaled`, and `$abs_dev_distsp_scaled` and/or
  `$sqr_dev_distsp_scaled`, data frames with deviation computed on
  scaled distance in functional spaces.

## Note

The maximum number of dimensions considered for assessing quality of
functional spaces depends on number of PC axes with positive eigenvalues
(i.e. axes with negative eigenvalues are not considered); so it could be
lower than `$maxdim_pcoa`. The quality metric obtained with
deviation_weighting = 'squared' and `fdist_scaling = TRUE` is equivalent
to the square-root of the 'mSD' originally suggested in Maire *et al.*
(2015).

## References

Maire *et al.* (2015) How many dimensions are needed to accurately
assess functional diversity? A pragmatic approach for assessing the
quality of functional spaces *Global Ecology and Biogeography*, **24**,
728-740.

## Author

Sebastien Villeger, Eva Maire, and Camille Magneville

## Examples

``` r
# Load Species x Traits Data
data("fruits_traits", package = "mFD")

# Load Traits x Categories Data
data("fruits_traits_cat", package = "mFD")

# Compute Functional Distance
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

# Compute Functional Spaces Quality (to retrieve species coordinates)
fspaces_quality_fruits <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fruits,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")
 fspaces_quality_fruits
#> $quality_fspaces
#>                     mad
#> pcoa_1d      0.13350558
#> pcoa_2d      0.06475567
#> pcoa_3d      0.04261771
#> pcoa_4d      0.03267394
#> pcoa_5d      0.03376338
#> pcoa_6d      0.04017987
#> pcoa_7d      0.04464819
#> pcoa_8d      0.04841744
#> pcoa_9d      0.05097090
#> pcoa_10d     0.05237960
#> tree_average 0.07154738
#> 
#> $details_trdist
#> $details_trdist$trdist_summary
#>       min       max      mean        sd 
#> 0.0289960 0.7191683 0.3593121 0.1369976 
#> 
#> $details_trdist$trdist_euclidean
#> [1] FALSE
#> 
#> 
#> $details_fspaces
#> $details_fspaces$sp_pc_coord
#>                         PC1          PC2          PC3          PC4          PC5
#> apple          0.0265193144  0.013054128 -0.043717435 -0.021211462 -0.148218537
#> apricot        0.1032240731 -0.148535464 -0.084372792 -0.059069034  0.051197559
#> banana        -0.3346937349  0.091371413 -0.070877547  0.069888722 -0.050000845
#> currant        0.2941221208 -0.002028080  0.059840313  0.027477102  0.036184138
#> blackberry     0.2398366750  0.002921130  0.096332589  0.012767053  0.063869252
#> blueberry      0.3146904236  0.065534392 -0.061375184  0.170320608 -0.043827091
#> cherry         0.1532479527 -0.250996898 -0.135570593  0.056216685  0.026507879
#> grape          0.0320396363 -0.175761176  0.131682624 -0.063280205 -0.101388849
#> grapefruit    -0.1565337233  0.091667465  0.051101664 -0.129100040  0.021725499
#> kiwifruit     -0.0076484035 -0.001811478  0.049335316 -0.133654078 -0.149411856
#> lemon         -0.1116123197  0.052718824  0.138487467 -0.103193749  0.023077059
#> lime          -0.2106569060  0.003790913  0.219519143  0.031204379  0.104197631
#> litchi        -0.2423438344 -0.252455139  0.019104502  0.089864709  0.130046060
#> mango         -0.3558930629 -0.084500319 -0.155978719 -0.069529132  0.072613533
#> melon          0.0475710274  0.228214444 -0.009073564 -0.130154617  0.037282406
#> orange        -0.0673938482  0.026151045 -0.002061839 -0.019971288 -0.014393513
#> passion_fruit -0.1743734750 -0.091951801  0.051587716  0.220532013 -0.147452934
#> peach          0.0403146225 -0.071052366 -0.144444375 -0.097071825  0.002288305
#> pear          -0.0008010876  0.019689428 -0.013204612 -0.070565107 -0.070145826
#> pineapple     -0.2308886096  0.336194377 -0.081546943  0.178095824  0.031816021
#> plum           0.0896173933 -0.145183768 -0.102585094 -0.018724455  0.019074804
#> raspberry      0.2180001689  0.011254750  0.088794687  0.026582440  0.044477454
#> strawberry     0.2743910893  0.093250971  0.046622086  0.093318143  0.049540096
#> tangerine     -0.0989381984 -0.035898731  0.072512352  0.008770777 -0.019765660
#> water_melon    0.1582027062  0.224361939 -0.120111760 -0.069513464  0.030707416
#>                        PC6          PC7          PC8          PC9          PC10
#> apple          0.040965006  0.039965979 -0.022577267 -0.012017348 -1.067572e-02
#> apricot        0.059097707 -0.027540550  0.056358883 -0.016901596  1.356995e-02
#> banana         0.020752971 -0.041476147 -0.083039226  0.042337245  1.583090e-02
#> currant        0.005808690 -0.114447728 -0.019688896 -0.024109923 -1.297464e-02
#> blackberry     0.035900843 -0.024364564  0.003797170  0.038186173 -9.737454e-03
#> blueberry     -0.072940389 -0.011254549  0.007080864 -0.046110072  4.357324e-03
#> cherry        -0.020780891  0.004803343 -0.042620591  0.030307476 -2.266872e-02
#> grape         -0.195578235 -0.001035913 -0.003686485  0.033761747 -3.964539e-03
#> grapefruit    -0.015750045 -0.051616104 -0.001017932  0.004568718 -2.143000e-02
#> kiwifruit     -0.010962090  0.006877412  0.052415441  0.009613759 -3.049371e-05
#> lemon          0.101128620  0.044498483 -0.034082720 -0.023135838 -2.153532e-02
#> lime           0.053717374  0.022699643  0.035206514  0.015350425  2.127588e-02
#> litchi        -0.107160270  0.013020916 -0.019139655 -0.016630151 -1.110479e-02
#> mango         -0.031431302 -0.007224234  0.043541389  0.011834427  1.624666e-04
#> melon         -0.136045993  0.036709357  0.013459634 -0.022081148  5.067072e-03
#> orange         0.091930421 -0.061574888 -0.012476483 -0.021006738 -2.425390e-02
#> passion_fruit  0.066901090 -0.006120818  0.066457679  0.020143504 -8.008755e-03
#> peach          0.073541526  0.076191641 -0.010596386 -0.012269290 -3.260209e-02
#> pear          -0.004973365  0.023163490 -0.055891141 -0.008962102  3.803879e-02
#> pineapple     -0.044607834  0.025787798  0.006751938 -0.010770300 -4.320707e-03
#> plum           0.054398550  0.019423316  0.018261782 -0.013010921  5.390245e-02
#> raspberry      0.045504907  0.029550728 -0.045185074  0.033817649  3.060490e-02
#> strawberry     0.001191018  0.094872368  0.018914728  0.012390681 -1.886095e-02
#> tangerine     -0.009194489 -0.033661520 -0.006545485 -0.078234499  1.424240e-02
#> water_melon   -0.001413820 -0.057247459  0.034301318  0.052928120  5.115953e-03
#> 
#> $details_fspaces$pc_eigenvalues
#>    Eigenvalues Relative_eig Rel_corr_eig Broken_stick Cum_corr_eig
#> 1   0.91554391  0.516168083   0.21167674   0.16236050    0.2116767
#> 2   0.47206587  0.266142709   0.12209464   0.11888224    0.3337714
#> 3   0.23421119  0.132044286   0.07404825   0.09714311    0.4078196
#> 4   0.23174667  0.130654830   0.07355042   0.08265036    0.4813700
#> 5   0.13715945  0.077328167   0.05444390   0.07178079    0.5358139
#> 6   0.12085527  0.068136146   0.05115047   0.06308514    0.5869644
#> 7   0.05019228  0.028297553   0.03687662   0.05583876    0.6238410
#> 8   0.03226312  0.018189400   0.03325495   0.04962758    0.6570960
#> 9   0.02181729  0.012300218   0.03114490   0.04419280    0.6882409
#> 10  0.01049296  0.005915751   0.02885740   0.03936188    0.7170983
#>    Cumul_br_stick
#> 1       0.1623605
#> 2       0.2812427
#> 3       0.3783858
#> 4       0.4610362
#> 5       0.5328170
#> 6       0.5959021
#> 7       0.6517409
#> 8       0.7013685
#> 9       0.7455613
#> 10      0.7849232
#> 
#> $details_fspaces$dendro
#> 
#> Call:
#> stats::hclust(d = sp_dist, method = fdendro)
#> 
#> Cluster method   : average 
#> Number of objects: 25 
#> 
#> 
#> $details_fspaces$pairsp_fspaces_dist
#>              sp.x          sp.y         tr     pcoa_1d    pcoa_2d    pcoa_3d
#> 2         apricot         apple 0.24018412 0.076704759 0.17887095 0.18343302
#> 3          banana         apple 0.37434510 0.361213049 0.36960582 0.37060240
#> 4         currant         apple 0.36866051 0.267602806 0.26802749 0.28733768
#> 5      blackberry         apple 0.32115315 0.213317361 0.21355789 0.25538399
#> 6       blueberry         apple 0.33047439 0.288171109 0.29291085 0.29344261
#> 7          cherry         apple 0.28842199 0.126728638 0.29288751 0.30695292
#> 8           grape         apple 0.30697774 0.005520322 0.18889598 0.25777291
#> 9      grapefruit         apple 0.26126651 0.183053038 0.19921966 0.22063348
#> 10      kiwifruit         apple 0.15238928 0.034167718 0.03726150 0.10023589
#> 11          lemon         apple 0.25445735 0.138131634 0.14371373 0.23206090
#> 12           lime         apple 0.40893135 0.237176220 0.23735704 0.35444585
#> 13         litchi         apple 0.48974636 0.268863149 0.37786580 0.38305242
#> 14          mango         apple 0.43111819 0.382412377 0.39465947 0.41031536
#> 15          melon         apple 0.28274850 0.021051713 0.21618773 0.21894596
#> 16         orange         apple 0.18148033 0.093913163 0.09482200 0.10356833
#> 17  passion_fruit         apple 0.31053460 0.200892789 0.22668074 0.24590085
#> 18          peach         apple 0.14865135 0.013795308 0.08523035 0.13194745
#> 19           pear         apple 0.06793207 0.027320402 0.02811461 0.04149053
#> 20      pineapple         apple 0.44936314 0.257407924 0.41313250 0.41486086
#> 21           plum         apple 0.21118812 0.063098079 0.17035433 0.18023873
#> 22      raspberry         apple 0.28826382 0.191480855 0.19148931 0.23286824
#> 23     strawberry         apple 0.33417000 0.247871775 0.26052246 0.27574115
#> 24      tangerine         apple 0.20798299 0.125457513 0.13466986 0.17789135
#> 25    water_melon         apple 0.30199453 0.131683392 0.24898094 0.26043732
#> 28         banana       apricot 0.52362013 0.437917808 0.49932686 0.49950920
#> 29        currant       apricot 0.29514305 0.190898048 0.24063765 0.28054215
#> 30     blackberry       apricot 0.24763570 0.136612602 0.20396594 0.27250053
#> 31      blueberry       apricot 0.34708347 0.211466350 0.30090517 0.30178272
#> 32         cherry       apricot 0.12305056 0.050023880 0.11402076 0.12498779
#> 33          grape       apricot 0.36685884 0.071184437 0.07621328 0.22910348
#> 34     grapefruit       apricot 0.37815310 0.259757796 0.35379593 0.37884679
#> 35      kiwifruit       apricot 0.26927586 0.110872477 0.18390387 0.22737302
#> 36          lemon       apricot 0.37134393 0.214836393 0.29437725 0.36922169
#> 37           lime       apricot 0.44248460 0.313880979 0.34889052 0.46268229
#> 38         litchi       apricot 0.38296079 0.345567908 0.36085520 0.37539849
#> 39          mango       apricot 0.38948413 0.459117136 0.46356126 0.46905911
#> 40          melon       apricot 0.39963509 0.055653046 0.38083823 0.38821093
#> 41         orange       apricot 0.24382146 0.170617921 0.24418405 0.25768381
#> 42  passion_fruit       apricot 0.43304196 0.277597548 0.28330568 0.31424094
#> 43          peach       apricot 0.11099248 0.062909451 0.09980596 0.11648959
#> 44           pear       apricot 0.22073690 0.104025161 0.19778991 0.21020409
#> 45      pineapple       apricot 0.62869977 0.334112683 0.58872260 0.58872939
#> 46           plum       apricot 0.02899600 0.013606680 0.01401341 0.02297963
#> 47      raspberry       apricot 0.25717061 0.114776096 0.19673959 0.26209434
#> 48     strawberry       apricot 0.32731921 0.171167016 0.29624116 0.32391123
#> 49      tangerine       apricot 0.26786408 0.202162272 0.23142303 0.27958821
#> 50    water_melon       apricot 0.31181041 0.054978633 0.37692854 0.37861907
#> 54        currant        banana 0.65209651 0.628815856 0.63571444 0.64901464
#> 55     blackberry        banana 0.60458916 0.574530410 0.58129910 0.60487013
#> 56      blueberry        banana 0.66239524 0.649384159 0.64989794 0.64996741
#> 57         cherry        banana 0.55397589 0.487941688 0.59607311 0.59957347
#> 58          grape        banana 0.53414780 0.366733371 0.45371046 0.49687404
#> 59     grapefruit        banana 0.28056111 0.178160012 0.17816026 0.21591666
#> 60      kiwifruit        banana 0.42168387 0.327045331 0.34006132 0.36068384
#> 61          lemon        banana 0.35708528 0.223081415 0.22640526 0.30837161
#> 62           lime        banana 0.34489261 0.124036829 0.15184031 0.32769760
#> 63         litchi        banana 0.38358308 0.092349900 0.35601292 0.36720835
#> 64          mango        banana 0.22343975 0.021199328 0.17714479 0.19652605
#> 65          melon        banana 0.46870976 0.382264762 0.40602015 0.41069709
#> 66         orange        banana 0.28208805 0.267299887 0.27514165 0.28361687
#> 67  passion_fruit        banana 0.31482060 0.160320260 0.24353642 0.27259444
#> 68          peach        banana 0.45632978 0.375008357 0.40867194 0.41524069
#> 69           pear        banana 0.33722666 0.333892647 0.34150052 0.34633622
#> 70      pineapple        banana 0.15956682 0.103805125 0.26592064 0.26613460
#> 71           plum        banana 0.49462413 0.424311128 0.48579655 0.48683021
#> 72      raspberry        banana 0.57169983 0.552693904 0.55847044 0.58084805
#> 73     strawberry        banana 0.64184843 0.609084824 0.60908772 0.62031768
#> 74      tangerine        banana 0.30626110 0.235755536 0.26791484 0.30387337
#> 75    water_melon        banana 0.50209721 0.492896441 0.51052266 0.51289121
#> 80     blackberry       currant 0.10451285 0.054285446 0.05451059 0.06559795
#> 81      blueberry       currant 0.14929723 0.020568303 0.07062395 0.14028877
#> 82         cherry       currant 0.30253705 0.140874168 0.28606119 0.34643387
#> 83          grape       currant 0.41200189 0.262082485 0.31443667 0.32253951
#> 84     grapefruit       currant 0.46381813 0.450655844 0.46029289 0.46037583
#> 85      kiwifruit       currant 0.39202880 0.301770524 0.30177060 0.30195339
#> 86          lemon       currant 0.45700896 0.405734440 0.40941136 0.41689691
#> 87           lime       currant 0.52814963 0.504779027 0.50481257 0.52946488
#> 88         litchi       currant 0.59477051 0.536465955 0.59203837 0.59343815
#> 89          mango       currant 0.68462718 0.650015184 0.65522623 0.68985452
#> 90          melon       currant 0.38540695 0.246551093 0.33734116 0.34430827
#> 91         orange       currant 0.37000846 0.361515969 0.36261255 0.36785831
#> 92  passion_fruit       currant 0.54535673 0.468495596 0.47704759 0.47711896
#> 93          peach       currant 0.38667582 0.253807498 0.26302585 0.33303878
#> 94           pear       currant 0.34921329 0.294923208 0.29572174 0.30460944
#> 95      pineapple       currant 0.64606505 0.525010730 0.62452438 0.64032886
#> 96           plum       currant 0.32413906 0.204504727 0.24963120 0.29782167
#> 97      raspberry       currant 0.14839119 0.076121952 0.07727215 0.08251873
#> 98     strawberry       currant 0.20237818 0.019731032 0.09730062 0.09819436
#> 99      tangerine       currant 0.39634047 0.393060319 0.39451696 0.39472043
#> 100   water_melon       currant 0.28535284 0.135919415 0.26405781 0.31954542
#> 106     blueberry    blackberry 0.21744644 0.074853749 0.09758844 0.18545955
#> 107        cherry    blackberry 0.27119131 0.086588722 0.26827593 0.35461396
#> 108         grape    blackberry 0.36449453 0.207797039 0.27405652 0.27632698
#> 109    grapefruit    blackberry 0.41493715 0.396370398 0.40618395 0.40869456
#> 110     kiwifruit    blackberry 0.34452145 0.247485078 0.24753032 0.25195239
#> 111         lemon    blackberry 0.35249611 0.351448995 0.35495944 0.35745383
#> 112          lime    blackberry 0.42363678 0.450493581 0.45049442 0.46703335
#> 113        litchi    blackberry 0.54726315 0.482180509 0.54563274 0.55107101
#> 114         mango    blackberry 0.63711982 0.595729738 0.60210998 0.65283798
#> 115         melon    blackberry 0.33789960 0.192265648 0.29618095 0.31437814
#> 116        orange    blackberry 0.32250111 0.307230523 0.30810749 0.32343730
#> 117 passion_fruit    blackberry 0.49784937 0.414210150 0.42493637 0.42728565
#> 118         peach    blackberry 0.33916847 0.199522052 0.21279363 0.32133265
#> 119          pear    blackberry 0.30170593 0.240637763 0.24122129 0.26492661
#> 120     pineapple    blackberry 0.59855769 0.470725285 0.57676109 0.60356812
#> 121          plum    blackberry 0.27663170 0.150219282 0.21095235 0.28994678
#> 122     raspberry    blackberry 0.04387834 0.021836506 0.02337268 0.02455814
#> 123    strawberry    blackberry 0.09809427 0.034554414 0.09671343 0.10874108
#> 124     tangerine    blackberry 0.34883311 0.338774873 0.34099178 0.34182276
#> 125   water_melon    blackberry 0.23784549 0.081633969 0.23600876 0.32023162
#> 132        cherry     blueberry 0.28859335 0.161442471 0.35532482 0.36298855
#> 133         grape     blueberry 0.44654304 0.282650787 0.37163829 0.41879152
#> 134    grapefruit     blueberry 0.55856990 0.471224147 0.47194823 0.48516613
#> 135     kiwifruit     blueberry 0.44969267 0.322338827 0.32929893 0.34741128
#> 136         lemon     blueberry 0.55176074 0.426302743 0.42649533 0.47100249
#> 137          lime     blueberry 0.62290140 0.525347330 0.52896321 0.59891877
#> 138        litchi     blueberry 0.62931166 0.557034258 0.64140822 0.64643753
#> 139         mango     blueberry 0.71916833 0.670583487 0.68716274 0.69364433
#> 140         melon     blueberry 0.35782967 0.267119396 0.31275801 0.31710098
#> 141        orange     blueberry 0.42423826 0.382084272 0.38410863 0.38866118
#> 142 passion_fruit     blueberry 0.50717061 0.489063899 0.51379509 0.52606655
#> 143         peach     blueberry 0.40959110 0.274375801 0.30649310 0.31755080
#> 144          pear     blueberry 0.36523546 0.315491511 0.31880504 0.32242373
#> 145     pineapple     blueberry 0.52444431 0.545579033 0.60902653 0.60936049
#> 146          plum     blueberry 0.33626929 0.225073030 0.30831804 0.31105991
#> 147     raspberry     blueberry 0.21890054 0.096690255 0.11088411 0.18667157
#> 148    strawberry     blueberry 0.13753399 0.040299334 0.04891058 0.11855655
#> 149     tangerine     blueberry 0.43088162 0.413628622 0.42588416 0.44643386
#> 150   water_melon     blueberry 0.27000500 0.156487717 0.22296770 0.23057446
#> 158         grape        cherry 0.32461636 0.121208316 0.14265998 0.30294579
#> 159    grapefruit        cherry 0.50120366 0.309781676 0.46193458 0.49822694
#> 160     kiwifruit        cherry 0.39232642 0.160896356 0.29661593 0.34952998
#> 161         lemon        cherry 0.49439449 0.264860272 0.40298164 0.48734179
#> 162          lime        cherry 0.56553516 0.363904859 0.44423369 0.56871108
#> 163        litchi        cherry 0.34071831 0.395591787 0.39559447 0.42475802
#> 164         mango        cherry 0.43057498 0.509141016 0.53567311 0.53606173
#> 165         melon        cherry 0.52268565 0.105676925 0.49072510 0.50676683
#> 166        orange        cherry 0.36687202 0.220641801 0.35425102 0.37857413
#> 167 passion_fruit        cherry 0.44657426 0.327621428 0.36418559 0.40946230
#> 168         peach        cherry 0.18555819 0.112933330 0.21244757 0.21263282
#> 169          pear        cherry 0.30786921 0.154049040 0.31145175 0.33462759
#> 170     pineapple        cherry 0.68930028 0.384136562 0.70167976 0.70375638
#> 171          plum        cherry 0.11223637 0.063630559 0.12347172 0.12780184
#> 172     raspberry        cherry 0.27264541 0.064752216 0.27012733 0.35115318
#> 173    strawberry        cherry 0.30238997 0.121143137 0.36494144 0.40789267
#> 174     tangerine        cherry 0.36023699 0.252186151 0.33145901 0.39136120
#> 175   water_melon        cherry 0.43486097 0.004954754 0.47538466 0.47563594
#> 184    grapefruit         grape 0.36622405 0.188573360 0.32722773 0.33700338
#> 185     kiwifruit         grape 0.18663975 0.039688040 0.17841984 0.19650628
#> 186         lemon         grape 0.39981893 0.143651956 0.26988700 0.26997278
#> 187          lime         grape 0.43055556 0.242696542 0.30189496 0.31441345
#> 188        litchi         grape 0.31723138 0.274383471 0.28490043 0.30633656
#> 189         mango         grape 0.47737471 0.387932699 0.39852268 0.49149707
#> 190         melon         grape 0.38770604 0.015531391 0.40427407 0.42807690
#> 191        orange         grape 0.37532676 0.099433485 0.22506791 0.26180746
#> 192 passion_fruit         grape 0.40452395 0.206413111 0.22277878 0.23673947
#> 193         peach         grape 0.40714424 0.008274986 0.10503528 0.29542940
#> 194          pear         grape 0.25369769 0.032840724 0.19819044 0.24550308
#> 195     pineapple         grape 0.67553280 0.262928246 0.57552563 0.61375614
#> 196          plum         grape 0.36210526 0.057577757 0.06519337 0.24316978
#> 197     raspberry         grape 0.35584762 0.185960533 0.26373524 0.26719965
#> 198    strawberry         grape 0.42599623 0.242351453 0.36207977 0.37193690
#> 199     tangerine         grape 0.26566142 0.130977835 0.19161601 0.20054381
#> 200   water_melon         grape 0.48776016 0.126163070 0.41954216 0.48930158
#> 210     kiwifruit    grapefruit 0.21998834 0.148885320 0.17579861 0.17580748
#> 211         lemon    grapefruit 0.11692821 0.044921404 0.05945527 0.10569393
#> 212          lime    grapefruit 0.23099817 0.054123183 0.10320662 0.19752482
#> 213        litchi    grapefruit 0.38232601 0.085810111 0.35466004 0.35610049
#> 214         mango    grapefruit 0.25652611 0.199359340 0.26604367 0.33713724
#> 215         melon    grapefruit 0.20074023 0.204104751 0.24556837 0.25283371
#> 216        orange    grapefruit 0.13433164 0.089139875 0.11062693 0.12273824
#> 217 passion_fruit    grapefruit 0.40513445 0.017839752 0.18448385 0.18448449
#> 218         peach    grapefruit 0.31564547 0.196848346 0.25539580 0.32166017
#> 219          pear    grapefruit 0.19333444 0.155732636 0.17156192 0.18321787
#> 220     pineapple    grapefruit 0.36745893 0.074354886 0.25558181 0.28795436
#> 221          plum    grapefruit 0.38896728 0.246151117 0.34159754 0.37457776
#> 222     raspberry    grapefruit 0.40629024 0.374533892 0.38306898 0.38491896
#> 223    strawberry    grapefruit 0.47643884 0.430924813 0.43092772 0.43095100
#> 224     tangerine    grapefruit 0.14096667 0.057595525 0.13996563 0.14159377
#> 225   water_melon    grapefruit 0.30674673 0.314736430 0.34156528 0.38207444
#> 236         lemon     kiwifruit 0.21317918 0.103963916 0.11739697 0.14741152
#> 237          lime     kiwifruit 0.36765318 0.203008503 0.20308579 0.26496485
#> 238        litchi     kiwifruit 0.48051948 0.234695431 0.34337180 0.34470000
#> 239         mango     kiwifruit 0.42189130 0.348244659 0.35792707 0.41263257
#> 240         melon     kiwifruit 0.24147034 0.055219431 0.23656101 0.24366516
#> 241        orange     kiwifruit 0.19212107 0.059745445 0.06596530 0.08362469
#> 242 passion_fruit     kiwifruit 0.35181277 0.166725072 0.18953239 0.18954577
#> 243         peach     kiwifruit 0.22050450 0.047963026 0.08423035 0.21129439
#> 244          pear     kiwifruit 0.08445721 0.006847316 0.02256490 0.06648622
#> 245     pineapple     kiwifruit 0.50858170 0.223240206 0.40507302 0.42569275
#> 246          plum     kiwifruit 0.28009005 0.097265797 0.17325198 0.23042582
#> 247     raspberry     kiwifruit 0.33587454 0.225648572 0.22602656 0.22944508
#> 248    strawberry     kiwifruit 0.40602314 0.282039493 0.29762921 0.29764157
#> 249     tangerine     kiwifruit 0.19875611 0.091289795 0.09744623 0.10016458
#> 250   water_melon     kiwifruit 0.32536283 0.165851110 0.28046569 0.32767868
#> 262          lime         lemon 0.15447400 0.099044586 0.11047068 0.13700330
#> 263        litchi         lemon 0.41592088 0.130731515 0.33199680 0.35280897
#> 264         mango         lemon 0.35729271 0.244280743 0.28018240 0.40646342
#> 265         melon         lemon 0.31766844 0.159183347 0.23693470 0.27912776
#> 266        orange         lemon 0.12752248 0.044218471 0.05158605 0.14971716
#> 267 passion_fruit         lemon 0.39832529 0.062761155 0.15769766 0.18005588
#> 268         peach         lemon 0.30883630 0.151926942 0.19596199 0.34416788
#> 269          pear         lemon 0.18652528 0.110811232 0.11562902 0.19073688
#> 270     pineapple         lemon 0.44398310 0.119276290 0.30754711 0.37815389
#> 271          plum         lemon 0.38215812 0.201229713 0.28223896 0.37118030
#> 272     raspberry         lemon 0.33286020 0.329612489 0.33221027 0.33590629
#> 273    strawberry         lemon 0.41422675 0.386003409 0.38812561 0.39884926
#> 274     tangerine         lemon 0.13415751 0.012674121 0.08951930 0.11120441
#> 275   water_melon         lemon 0.40751332 0.269815026 0.31978353 0.41126034
#> 288        litchi          lime 0.27999084 0.031686928 0.25819779 0.32685185
#> 289         mango          lime 0.32085761 0.145236157 0.16996730 0.41217415
#> 290         melon          lime 0.43173840 0.258227933 0.34212218 0.41146350
#> 291        orange          lime 0.28199648 0.143263058 0.14499751 0.26480637
#> 292 passion_fruit          lime 0.30279928 0.036283431 0.10238728 0.19668279
#> 293         peach          lime 0.46331030 0.250971529 0.26189354 0.44839454
#> 294          pear          lime 0.34099928 0.209855818 0.21045719 0.31377153
#> 295     pineapple          lime 0.43179043 0.020231704 0.33301859 0.44893449
#> 296          plum          lime 0.45329878 0.300274299 0.33519861 0.46487552
#> 297     raspberry          lime 0.40400086 0.428657075 0.42872205 0.44820919
#> 298    strawberry          lime 0.48536741 0.485047995 0.49322881 0.52265481
#> 299     tangerine          lime 0.20529817 0.111718708 0.11855943 0.18885797
#> 300   water_melon          lime 0.53774489 0.368859612 0.42977784 0.54777563
#> 314         mango        litchi 0.16014333 0.113549228 0.20273690 0.26787382
#> 315         melon        litchi 0.57047466 0.289914862 0.56133223 0.56203904
#> 316        orange        litchi 0.39142871 0.174949986 0.32898162 0.32966182
#> 317 passion_fruit        litchi 0.36507035 0.067970359 0.17430230 0.17730327
#> 318         peach        litchi 0.42324620 0.282658457 0.33586124 0.37356526
#> 319          pear        litchi 0.43646631 0.241542747 0.36387575 0.36530733
#> 320     pineapple        litchi 0.52496809 0.011455225 0.58876097 0.59730243
#> 321          plum        litchi 0.37820721 0.331961228 0.34886302 0.36947769
#> 322     raspberry        litchi 0.53861624 0.460344003 0.53052757 0.53508525
#> 323    strawberry        litchi 0.60876485 0.516734924 0.62171352 0.62232219
#> 324     tangerine        litchi 0.28176338 0.143405636 0.25973420 0.26516835
#> 325   water_melon        litchi 0.67052878 0.400546541 0.62272952 0.63810126
#> 340         melon         mango 0.44467477 0.403464090 0.51046429 0.53118257
#> 341        orange         mango 0.31461871 0.288499215 0.30899113 0.34520418
#> 342 passion_fruit         mango 0.37159368 0.181519588 0.18167247 0.27584182
#> 343         peach         mango 0.34643620 0.396207685 0.39643584 0.39660360
#> 344          pear         mango 0.37783813 0.355091975 0.37006191 0.39664879
#> 345     pineapple         mango 0.38300658 0.125004453 0.43887372 0.44514068
#> 346          plum         mango 0.38473055 0.445510456 0.44962434 0.45278353
#> 347     raspberry         mango 0.62847292 0.573893232 0.58182684 0.63121826
#> 348    strawberry         mango 0.69862152 0.630284152 0.65486917 0.68549305
#> 349     tangerine         mango 0.30646853 0.256954864 0.26151083 0.34726947
#> 350   water_melon         mango 0.51038545 0.514095769 0.59974191 0.60081344
#> 366        orange         melon 0.32248030 0.114964876 0.23247912 0.23258483
#> 367 passion_fruit         melon 0.59328311 0.221944502 0.38957129 0.39426587
#> 368         peach         melon 0.33712746 0.007256405 0.29935477 0.32854001
#> 369          pear         melon 0.21481643 0.048372115 0.21406201 0.21410187
#> 370     pineapple         melon 0.33338536 0.278459637 0.29866274 0.30733016
#> 371          plum         melon 0.41044927 0.042046366 0.37575806 0.38721897
#> 372     raspberry         melon 0.32925269 0.170429142 0.27589419 0.29273844
#> 373    strawberry         melon 0.28829018 0.226820062 0.26393651 0.26974893
#> 374     tangerine         melon 0.32911533 0.146509226 0.30202768 0.31285298
#> 375   water_melon         melon 0.10600649 0.110631679 0.11069874 0.15679187
#> 392 passion_fruit        orange 0.30716644 0.106979627 0.15935157 0.16814041
#> 393         peach        orange 0.21666597 0.107708471 0.14508486 0.20327913
#> 394          pear        orange 0.14587149 0.066592761 0.06690552 0.06782706
#> 395     pineapple        orange 0.38716769 0.163494761 0.35051021 0.35940964
#> 396          plum        orange 0.25463564 0.157011242 0.23239653 0.25320559
#> 397     raspberry        orange 0.31385420 0.285394017 0.28578251 0.29987756
#> 398    strawberry        orange 0.38400280 0.341784938 0.34830926 0.35169514
#> 399     tangerine        orange 0.10966533 0.031544350 0.06960762 0.10201240
#> 400   water_melon        orange 0.30334249 0.225596554 0.30030212 0.32267189
#> 418         peach passion_fruit 0.44100413 0.214688098 0.21570296 0.29147272
#> 419          pear passion_fruit 0.37846667 0.173572387 0.20637620 0.21630807
#> 420     pineapple passion_fruit 0.42656441 0.056515135 0.43186006 0.45191586
#> 421          plum passion_fruit 0.40404595 0.263990868 0.26930433 0.31031287
#> 422     raspberry passion_fruit 0.48112166 0.392373644 0.40571994 0.40742242
#> 423    strawberry passion_fruit 0.51086622 0.448764564 0.48547884 0.48550423
#> 424     tangerine passion_fruit 0.26416778 0.075435277 0.09398100 0.09628223
#> 425   water_melon passion_fruit 0.57818570 0.332576181 0.45897854 0.49004286
#> 444          pear         peach 0.15344655 0.041115710 0.09962216 0.16476787
#> 445     pineapple         peach 0.53716700 0.271203232 0.48928632 0.49331247
#> 446          plum         peach 0.09756424 0.049302771 0.08902937 0.09837900
#> 447     raspberry         peach 0.30627914 0.177685546 0.19582292 0.30454404
#> 448    strawberry         peach 0.35218531 0.234076467 0.28598493 0.34393862
#> 449     tangerine         peach 0.30814949 0.139252821 0.14362147 0.26018714
#> 450   water_melon         peach 0.32000985 0.117888084 0.31806794 0.31899732
#> 470     pineapple          pear 0.42412449 0.230087522 0.39129995 0.39722327
#> 471          plum          pear 0.19563284 0.090418481 0.18803902 0.20820073
#> 472     raspberry          pear 0.26881660 0.218801257 0.21896377 0.24155536
#> 473    strawberry          pear 0.33896520 0.275192177 0.28485441 0.29106918
#> 474     tangerine          pear 0.15470294 0.098137111 0.11278713 0.14166275
#> 475   water_melon          pear 0.28254731 0.159003794 0.25917763 0.28036081
#> 496          plum     pineapple 0.60130633 0.320506003 0.57831567 0.57869821
#> 497     raspberry     pineapple 0.56566836 0.448888779 0.55415422 0.57974405
#> 498    strawberry     pineapple 0.50046343 0.505279699 0.56065058 0.57511423
#> 499     tangerine     pineapple 0.40987138 0.131950411 0.39479639 0.42379058
#> 500   water_melon     pineapple 0.38495463 0.389091316 0.40484385 0.40667652
#> 522     raspberry          plum 0.24374237 0.128382776 0.20237378 0.27853432
#> 523    strawberry          plum 0.31389097 0.184773696 0.30164954 0.33653414
#> 524     tangerine          plum 0.26311050 0.188555592 0.21793676 0.27956313
#> 525   water_melon          plum 0.34080642 0.068585313 0.37585632 0.37626475
#> 548    strawberry     raspberry 0.08136655 0.056390920 0.09951541 0.10808258
#> 549     tangerine     raspberry 0.34018620 0.316938367 0.32042687 0.32084029
#> 550   water_melon     raspberry 0.26556222 0.059797463 0.22133777 0.30435557
#> 574     tangerine    strawberry 0.41033480 0.373329288 0.39503722 0.39588472
#> 575   water_melon    strawberry 0.22459971 0.116188383 0.17518512 0.24184706
#> 600   water_melon     tangerine 0.41300783 0.257140905 0.36586481 0.41347444
#>        pcoa_4d    pcoa_5d    pcoa_6d    pcoa_7d    pcoa_8d    pcoa_9d
#> 2   0.18729888 0.27358298 0.27418323 0.28237134 0.29319702 0.29323770
#> 3   0.38163514 0.39407118 0.39458918 0.40290624 0.40741759 0.41102738
#> 4   0.29143356 0.34487369 0.34666097 0.37949628 0.37950727 0.37969988
#> 5   0.25763447 0.33370159 0.33374002 0.33988354 0.34090531 0.34458210
#> 6   0.35041846 0.36563735 0.38296881 0.38637890 0.38751549 0.38901230
#> 7   0.31656786 0.36158614 0.36682024 0.36850170 0.36904639 0.37146551
#> 8   0.26118318 0.26534821 0.35547486 0.35783171 0.35833001 0.36124247
#> 9   0.24559943 0.29866378 0.30400107 0.31749634 0.31822749 0.31865943
#> 10  0.15063391 0.15063864 0.15933745 0.16273683 0.17918477 0.18048570
#> 11  0.24611655 0.29985920 0.30583525 0.30586883 0.30608515 0.30628702
#> 12  0.35830055 0.43828439 0.43846987 0.43880970 0.44259792 0.44344325
#> 13  0.39883214 0.48631087 0.50836931 0.50908290 0.50909450 0.50911540
#> 14  0.41315044 0.46846568 0.47402671 0.47636985 0.48093649 0.48152759
#> 15  0.24455254 0.30694714 0.35432956 0.35434453 0.35617229 0.35631444
#> 16  0.10357576 0.16922492 0.17673299 0.20382614 0.20407627 0.20427416
#> 17  0.34482914 0.34482999 0.34580400 0.34886157 0.36004391 0.36147744
#> 18  0.15220028 0.21404961 0.21651435 0.21952394 0.21985064 0.21985078
#> 19  0.06447671 0.10125510 0.11118871 0.11245111 0.11728199 0.11732178
#> 20  0.46025311 0.49421187 0.50156563 0.50176598 0.50262242 0.50262397
#> 21  0.18025589 0.24592529 0.24629192 0.24714714 0.25049858 0.25050055
#> 22  0.23772226 0.30601245 0.30604612 0.30622330 0.30705671 0.31045880
#> 23  0.29858032 0.35813222 0.36033409 0.36449330 0.36684732 0.36765841
#> 24  0.18040030 0.22145972 0.22706912 0.23870776 0.23924550 0.24824005
#> 25  0.26487862 0.31964853 0.32244557 0.33678123 0.34155054 0.34767037
#> 28  0.51588714 0.52571918 0.52711571 0.52729989 0.54541453 0.54862214
#> 29  0.29358837 0.29397200 0.29876287 0.31114645 0.32030513 0.32038623
#> 30  0.28181015 0.28209490 0.28304704 0.28306486 0.28790354 0.29312644
#> 31  0.37906783 0.39079675 0.41249989 0.41282126 0.41575199 0.41677674
#> 32  0.17003748 0.17182062 0.18948065 0.19222134 0.21620819 0.22130223
#> 33  0.22914218 0.27529757 0.37503145 0.37596687 0.38073157 0.38408763
#> 34  0.38526514 0.38639078 0.39357340 0.39430908 0.39846173 0.39903975
#> 35  0.23929359 0.31225880 0.32002177 0.32186725 0.32189141 0.32298165
#> 36  0.37184896 0.37291073 0.37527191 0.38212384 0.39268093 0.39273042
#> 37  0.47140661 0.47437664 0.47440715 0.47705998 0.47752869 0.47861659
#> 38  0.40386296 0.41148800 0.44380636 0.44565606 0.45200592 0.45200600
#> 39  0.46917573 0.46966425 0.47830954 0.47874081 0.47891236 0.47977371
#> 40  0.39466554 0.39491078 0.44049470 0.44515574 0.44721805 0.44724804
#> 41  0.26063303 0.26875968 0.27075774 0.27288842 0.28143631 0.28146625
#> 42  0.42062348 0.46517323 0.46523868 0.46573150 0.46584098 0.46731163
#> 43  0.12253178 0.13193238 0.13272067 0.16844923 0.18126817 0.18132735
#> 44  0.21051821 0.24298587 0.25129114 0.25635549 0.27985390 0.27996650
#> 45  0.63470423 0.63500008 0.64341273 0.64561897 0.64752197 0.64755100
#> 46  0.04643004 0.05645901 0.05665423 0.07358877 0.08286553 0.08295682
#> 47  0.27573468 0.27581656 0.27615130 0.28199105 0.29971675 0.30397791
#> 48  0.35796695 0.35797079 0.36262414 0.38272861 0.38455592 0.38566993
#> 49  0.28770090 0.29632345 0.30409112 0.30415272 0.31058950 0.31658736
#> 50  0.37876310 0.37931692 0.38411323 0.38526027 0.38589119 0.39215839
#> 54  0.65039892 0.65608429 0.65625447 0.66029901 0.66333103 0.66665079
#> 55  0.60756132 0.61814008 0.61832566 0.61856238 0.62462788 0.62464167
#> 56  0.65768092 0.65770990 0.66434988 0.66503692 0.67111529 0.67691850
#> 57  0.59972933 0.60458982 0.60601478 0.60777932 0.60912180 0.60924058
#> 58  0.51441012 0.51697050 0.56040850 0.56186573 0.56744159 0.56750639
#> 59  0.29362652 0.30226015 0.30445635 0.30462516 0.31547422 0.31772700
#> 60  0.41415275 0.42591672 0.42709589 0.42982434 0.45066277 0.45184926
#> 61  0.35362493 0.36109690 0.36993407 0.37979317 0.38293549 0.38849236
#> 62  0.32997302 0.36422433 0.36571303 0.37130116 0.38967500 0.39060836
#> 63  0.36775129 0.40946050 0.42897517 0.43242297 0.43711873 0.44107815
#> 64  0.24095606 0.27035923 0.27534943 0.27747162 0.30498058 0.30650216
#> 65  0.45682539 0.46508902 0.49080924 0.49699767 0.50627929 0.51036109
#> 66  0.29751193 0.29963516 0.30797315 0.30862829 0.31659204 0.32286682
#> 67  0.31145004 0.32634037 0.32958715 0.33147803 0.36363032 0.36430698
#> 68  0.44754961 0.45059384 0.45367548 0.46868662 0.47425216 0.47738558
#> 69  0.37373259 0.37427513 0.37515826 0.38068622 0.38165301 0.38508524
#> 70  0.28729149 0.29871458 0.30578168 0.31309244 0.32571357 0.33001475
#> 71  0.49482921 0.49962726 0.50075885 0.50444838 0.51451925 0.51748766
#> 72  0.58246020 0.59007291 0.59059182 0.59484747 0.59605071 0.59611159
#> 73  0.62075998 0.62869019 0.62899446 0.64360310 0.65162839 0.65231614
#> 74  0.30995875 0.31142992 0.31286650 0.31296408 0.32217667 0.34399906
#> 75  0.53149822 0.53759109 0.53804791 0.53827900 0.55092022 0.55102201
#> 80  0.06722705 0.07270448 0.07868595 0.11960960 0.12189361 0.13688994
#> 81  0.20021290 0.21560845 0.22953958 0.25166894 0.25308867 0.25404307
#> 82  0.34762392 0.34775857 0.34877360 0.36859713 0.36930977 0.37329741
#> 83  0.33506510 0.36220844 0.41442930 0.42966718 0.42996507 0.43384224
#> 84  0.48627390 0.48648881 0.48696626 0.49100301 0.49135787 0.49219409
#> 85  0.34225591 0.38933916 0.38970019 0.40814952 0.41446962 0.41583933
#> 86  0.43689576 0.43709233 0.44736516 0.47476256 0.47498071 0.47498170
#> 87  0.52947800 0.53382842 0.53597390 0.55324265 0.55595947 0.55735810
#> 88  0.59670852 0.60404562 0.61451859 0.62759968 0.62759992 0.62764449
#> 89  0.69664157 0.69759342 0.69858671 0.70676748 0.70959026 0.71050006
#> 90  0.37867657 0.37867816 0.40437594 0.43170402 0.43297481 0.43297957
#> 91  0.37090576 0.37433833 0.38411735 0.38773918 0.38780626 0.38781867
#> 92  0.51469671 0.54647532 0.54987959 0.56044828 0.56703043 0.56875467
#> 93  0.35556612 0.35717810 0.36354358 0.41049641 0.41059709 0.41076779
#> 94  0.31999873 0.33720208 0.33737442 0.36436019 0.36615427 0.36646747
#> 95  0.65780472 0.65781923 0.65974841 0.67448793 0.67500598 0.67513778
#> 96  0.30138403 0.30186928 0.30575486 0.33377761 0.33592818 0.33611149
#> 97  0.08252358 0.08293926 0.09194950 0.17085159 0.17274351 0.18219749
#> 98  0.11822510 0.11897713 0.11906670 0.24081483 0.24388936 0.24660559
#> 99  0.39516344 0.39910465 0.39938655 0.40747519 0.40768711 0.41126421
#> 100 0.33394078 0.33398569 0.33406377 0.33892547 0.34319880 0.35173894
#> 106 0.24334825 0.26611439 0.28751223 0.28781097 0.28782971 0.29991965
#> 107 0.35726591 0.35921414 0.36365866 0.36482652 0.36776758 0.36785197
#> 108 0.28660039 0.33083232 0.40377294 0.40444631 0.40451554 0.40453973
#> 109 0.43261705 0.43466494 0.43772300 0.43857048 0.43859692 0.43988338
#> 110 0.29140891 0.36112046 0.36414848 0.36548622 0.36870573 0.36981117
#> 111 0.37579269 0.37800020 0.38358677 0.38971904 0.39155564 0.39632841
#> 112 0.46739714 0.46913374 0.46947193 0.47182511 0.47286941 0.47342048
#> 113 0.55643805 0.56035942 0.57833309 0.57954019 0.57999391 0.58257855
#> 114 0.65800463 0.65806273 0.66149843 0.66172045 0.66291294 0.66343649
#> 115 0.34534073 0.34636264 0.38669470 0.39148795 0.39160717 0.39621753
#> 116 0.32508997 0.33437785 0.33903961 0.34107546 0.34146347 0.34655606
#> 117 0.47512030 0.51999651 0.52091975 0.52123912 0.52499196 0.52530191
#> 118 0.33958689 0.34512529 0.34717184 0.36144133 0.36172782 0.36522974
#> 119 0.27772352 0.30836731 0.31106446 0.31467446 0.32028535 0.32373704
#> 120 0.62580195 0.62662229 0.63177301 0.63376052 0.63376741 0.63565546
#> 121 0.29165194 0.29507185 0.29565108 0.29887613 0.29922594 0.30357422
#> 122 0.02817742 0.03420539 0.03552812 0.06456861 0.08104546 0.08116311
#> 123 0.13532590 0.13608242 0.14043930 0.18422986 0.18484908 0.18664027
#> 124 0.34184612 0.35192836 0.35480580 0.35492759 0.35507825 0.37367678
#> 125 0.33063330 0.33229217 0.33438072 0.33599368 0.33737554 0.33769747
#> 132 0.38050019 0.38694625 0.39044592 0.39077599 0.39392398 0.40126767
#> 133 0.47953694 0.48297933 0.49830620 0.49841097 0.49852726 0.50488507
#> 134 0.57012183 0.57387807 0.57672071 0.57813133 0.57818805 0.58040483
#> 135 0.46162237 0.47354340 0.47758210 0.47792618 0.48007151 0.48329473
#> 136 0.54465902 0.54875278 0.57569926 0.57839263 0.57985557 0.58031051
#> 137 0.61486342 0.63243050 0.64498878 0.64588189 0.64649398 0.64940886
#> 138 0.65142508 0.67423031 0.67509815 0.67553446 0.67604314 0.67668559
#> 139 0.73394166 0.74312097 0.74427937 0.74429028 0.74518279 0.74743224
#> 140 0.43685054 0.44431649 0.44877551 0.45133136 0.45137643 0.45201557
#> 141 0.43274533 0.43374515 0.46402288 0.46674337 0.46715293 0.46782693
#> 142 0.52845738 0.53852160 0.55638220 0.55640589 0.55956511 0.56347374
#> 143 0.41513519 0.41768871 0.44262942 0.45118471 0.45153087 0.45279723
#> 144 0.40247110 0.40333071 0.40901733 0.41046288 0.41526528 0.41692352
#> 145 0.60941009 0.61408676 0.61474001 0.61585503 0.61585511 0.61686824
#> 146 0.36400042 0.36939539 0.39072773 0.39193021 0.39208966 0.39348425
#> 147 0.23559910 0.25160411 0.27808976 0.28106759 0.28588584 0.29684871
#> 148 0.14136844 0.16941802 0.18492683 0.21321551 0.21354366 0.22141191
#> 149 0.47476472 0.47537405 0.47962905 0.48015216 0.48034548 0.48141848
#> 150 0.33269350 0.34094040 0.34836246 0.35138548 0.35243823 0.36608916
#> 158 0.32566188 0.34987602 0.39111040 0.39115399 0.39308689 0.39310207
#> 159 0.53157537 0.53159688 0.53162069 0.53460612 0.53622242 0.53683980
#> 160 0.39777143 0.43493662 0.43504744 0.43505238 0.44531160 0.44579216
#> 161 0.51275112 0.51276260 0.52705541 0.52854812 0.52861707 0.53131177
#> 162 0.56926085 0.57453774 0.57934756 0.57962391 0.58482556 0.58501679
#> 163 0.42608868 0.43848799 0.44691511 0.44699065 0.44760697 0.45006126
#> 164 0.55061256 0.55253952 0.55264215 0.55277302 0.55944785 0.55975276
#> 165 0.53995082 0.54005831 0.55222189 0.55314285 0.55597842 0.55844120
#> 166 0.38616445 0.38832448 0.40435102 0.40976312 0.41087040 0.41406235
#> 167 0.44120166 0.47425865 0.48229596 0.48241967 0.49459762 0.49470204
#> 168 0.26212608 0.26324261 0.27963081 0.28859951 0.29037085 0.29347574
#> 169 0.35783970 0.37066317 0.37100009 0.37145412 0.37169109 0.37375977
#> 170 0.71423216 0.71425188 0.71464919 0.71495721 0.71665994 0.71783623
#> 171 0.14815358 0.14833993 0.16630299 0.16694438 0.17769944 0.18290318
#> 172 0.35240140 0.35285925 0.35903128 0.35988316 0.35989230 0.35990942
#> 173 0.40957654 0.41022363 0.41081163 0.42056940 0.42504731 0.42542476
#> 174 0.39422671 0.39693317 0.39710224 0.39896082 0.40058850 0.41503314
#> 175 0.49197319 0.49199111 0.49237216 0.49626670 0.50219281 0.50270201
#> 184 0.34337084 0.36477483 0.40669258 0.40982583 0.40983451 0.41087293
#> 185 0.20872757 0.21418078 0.28276586 0.28287656 0.28838616 0.28939541
#> 186 0.27290729 0.29995025 0.42190652 0.42435656 0.42544380 0.42923159
#> 187 0.32830345 0.38736153 0.46064873 0.46125983 0.46289664 0.46326264
#> 188 0.34248424 0.41334921 0.42270002 0.42293368 0.42321590 0.42620540
#> 189 0.49153680 0.52142617 0.54665297 0.54668799 0.54872419 0.54916213
#> 190 0.43326899 0.45491948 0.45879824 0.46034827 0.46066747 0.46403981
#> 191 0.26536542 0.27926152 0.40080946 0.40535563 0.40545093 0.40913328
#> 192 0.36958754 0.37244711 0.45564487 0.45567324 0.46104046 0.46124155
#> 193 0.29735569 0.31491166 0.41424003 0.42137738 0.42143403 0.42394044
#> 194 0.24561114 0.24759030 0.31245988 0.31339558 0.31771388 0.32057361
#> 195 0.65951420 0.67283172 0.68956115 0.69008267 0.69016162 0.69159682
#> 196 0.24721803 0.27500590 0.37164047 0.37220320 0.37284976 0.37577204
#> 197 0.28190592 0.31740814 0.39858375 0.39975561 0.40190382 0.40190383
#> 198 0.40355929 0.43085920 0.47366421 0.48327651 0.48380471 0.48427649
#> 199 0.21309426 0.22819182 0.29463606 0.29643690 0.29645069 0.31690089
#> 200 0.48934128 0.50685729 0.54277448 0.54567745 0.54699812 0.54733381
#> 210 0.17586646 0.24539153 0.24543824 0.25231215 0.25790802 0.25795736
#> 211 0.10882253 0.10883093 0.15970220 0.18639422 0.18930422 0.19132075
#> 212 0.25438860 0.26742329 0.27629864 0.28611845 0.28840246 0.28860392
#> 213 0.41803483 0.43184079 0.44140944 0.44611684 0.44648475 0.44698773
#> 214 0.34235977 0.34612108 0.34647613 0.34930838 0.35213900 0.35221395
#> 215 0.25283591 0.25331406 0.28042669 0.29400768 0.29436392 0.29556781
#> 216 0.16423691 0.16816166 0.19968332 0.19993150 0.20025959 0.20188612
#> 217 0.39531899 0.42999819 0.43786945 0.44022662 0.44536775 0.44564000
#> 218 0.32325079 0.32383465 0.33591943 0.35941158 0.35953919 0.35993326
#> 219 0.19234117 0.21315597 0.21342822 0.22614949 0.23271154 0.23310458
#> 220 0.42105464 0.42117553 0.42216300 0.42920038 0.42927070 0.42954467
#> 221 0.39050130 0.39051030 0.39676078 0.40307036 0.40353120 0.40391394
#> 222 0.41521036 0.41583325 0.42032066 0.42808587 0.43035828 0.43135108
#> 223 0.48496249 0.48575948 0.48605480 0.50764962 0.50804080 0.50810101
#> 224 0.19762884 0.20193730 0.20204368 0.20283988 0.20291518 0.21915963
#> 225 0.38669295 0.38679725 0.38706284 0.38710380 0.38871172 0.39170836
#> 236 0.15052570 0.22893321 0.25490143 0.25766274 0.27179407 0.27376003
#> 237 0.31206519 0.40212244 0.40729092 0.40759813 0.40796125 0.40800158
#> 238 0.41082690 0.49686564 0.50609244 0.50612973 0.51116282 0.51183608
#> 239 0.41758550 0.47294072 0.47338347 0.47359346 0.47367659 0.47368180
#> 240 0.24369028 0.30698486 0.33149010 0.33282973 0.33510175 0.33659731
#> 241 0.14112712 0.19531210 0.22075708 0.23112638 0.24006325 0.24200822
#> 242 0.40171555 0.40172033 0.40919665 0.40940304 0.40964379 0.40977910
#> 243 0.21443783 0.26267189 0.27593003 0.28450280 0.29139721 0.29221773
#> 244 0.09165498 0.12117648 0.12132437 0.12241258 0.16344772 0.16449990
#> 245 0.52763843 0.55789413 0.55890777 0.55922759 0.56108881 0.56145896
#> 246 0.25749734 0.30772168 0.31458647 0.31483654 0.31668362 0.31749078
#> 247 0.27985851 0.34046124 0.34511212 0.34585612 0.35936376 0.36017793
#> 248 0.37430856 0.42389713 0.42407131 0.43310459 0.43439830 0.43440718
#> 249 0.17412002 0.21708505 0.21709224 0.22084485 0.22858005 0.24487988
#> 250 0.33389719 0.37938145 0.37950159 0.38488109 0.38530712 0.38773407
#> 262 0.19191863 0.20835860 0.21368466 0.21479367 0.22569298 0.22895090
#> 263 0.40217625 0.41615875 0.46537336 0.46643671 0.46667601 0.46672135
#> 264 0.40785514 0.41085238 0.43170802 0.43479541 0.44167019 0.44305245
#> 265 0.28042681 0.28078638 0.36754971 0.36763223 0.37069358 0.37069508
#> 266 0.17129275 0.17534323 0.17558433 0.20513755 0.20627226 0.20628325
#> 267 0.37043014 0.40779770 0.40923158 0.41235034 0.42443041 0.42663131
#> 268 0.34422232 0.34484950 0.34595119 0.34739989 0.34819290 0.34836242
#> 269 0.19350758 0.21479221 0.23956904 0.24051716 0.24150386 0.24191942
#> 270 0.47130053 0.47138154 0.49339606 0.49375071 0.49543640 0.49559069
#> 271 0.38067030 0.38069133 0.38354868 0.38436747 0.38791533 0.38804744
#> 272 0.36010400 0.36073934 0.36500256 0.36530850 0.36547718 0.36988818
#> 273 0.44463205 0.44541885 0.45649258 0.45926354 0.46231129 0.46367431
#> 274 0.15780518 0.16351750 0.19725406 0.21217481 0.21395431 0.22093508
#> 275 0.41263716 0.41270771 0.42525593 0.43725832 0.44257340 0.44906229
#> 288 0.33207403 0.33307853 0.36989582 0.37002243 0.37399212 0.37535698
#> 289 0.42430504 0.42547893 0.43391545 0.43494604 0.43502589 0.43504010
#> 290 0.44197165 0.44700848 0.48561993 0.48582198 0.48630846 0.48774690
#> 291 0.26970606 0.29462725 0.29709503 0.30881654 0.31247611 0.31458411
#> 292 0.27300013 0.37129110 0.37152509 0.37264126 0.37394939 0.37398010
#> 293 0.46638230 0.47738660 0.47779803 0.48078306 0.48295989 0.48374901
#> 294 0.32986300 0.37310219 0.37769014 0.37769043 0.38852136 0.38928132
#> 295 0.47235502 0.47786857 0.48787930 0.48788908 0.48871813 0.48941568
#> 296 0.46754908 0.47523471 0.47523520 0.47524650 0.47554848 0.47639345
#> 297 0.44823302 0.45219392 0.45226849 0.45232038 0.45940889 0.45977991
#> 298 0.52633276 0.52916313 0.53176370 0.53663911 0.53688636 0.53689452
#> 299 0.19018570 0.22701872 0.23557462 0.24222300 0.24579506 0.26300827
#> 300 0.55695801 0.56178558 0.56448426 0.57011755 0.57011827 0.57135534
#> 314 0.31170945 0.31695627 0.32587751 0.32650577 0.33246794 0.33368423
#> 315 0.60356970 0.61065660 0.61133940 0.61179817 0.61266607 0.61269032
#> 316 0.34747786 0.37630261 0.42572380 0.43220978 0.43226114 0.43228329
#> 317 0.22025075 0.35428250 0.39473212 0.39519596 0.40435968 0.40602839
#> 318 0.41772751 0.43682756 0.47272768 0.47692976 0.47700627 0.47702620
#> 319 0.39898267 0.44638992 0.45793681 0.45804912 0.45952113 0.45958510
#> 320 0.60378383 0.61172221 0.61491208 0.61504460 0.61558934 0.61561723
#> 321 0.38510436 0.40077424 0.43211253 0.43215996 0.43377540 0.43379050
#> 322 0.53881432 0.54556655 0.56652407 0.56676516 0.56736330 0.56960170
#> 323 0.62233178 0.62751737 0.63680299 0.64204183 0.64316860 0.64382300
#> 324 0.27729132 0.31517302 0.33004746 0.33333253 0.33357036 0.33921127
#> 325 0.65770406 0.66516373 0.67351697 0.67717261 0.67927806 0.68283017
#> 340 0.53463106 0.53579722 0.54591473 0.54767970 0.54850521 0.54955275
#> 341 0.34874333 0.35943307 0.38001348 0.38388050 0.38794619 0.38933377
#> 342 0.40028012 0.45678598 0.46725013 0.46725144 0.46781306 0.46788685
#> 343 0.39755882 0.40373092 0.41715459 0.42541293 0.42884386 0.42952072
#> 344 0.39665014 0.42155850 0.42238796 0.42347964 0.43499636 0.43549320
#> 345 0.50938036 0.51101153 0.51118138 0.51224623 0.51356564 0.51406287
#> 346 0.45562489 0.45875967 0.46671962 0.46747972 0.46816274 0.46882155
#> 347 0.63849348 0.63911311 0.64372723 0.64477682 0.65085293 0.65122407
#> 348 0.70457076 0.70494846 0.70570288 0.71304998 0.71347512 0.71347533
#> 349 0.35598731 0.36777830 0.36844993 0.36939719 0.37277738 0.38350409
#> 350 0.60081344 0.60227312 0.60302070 0.60509197 0.60516251 0.60655615
#> 366 0.25736369 0.26250042 0.34767760 0.36130252 0.36223223 0.36223383
#> 367 0.52766153 0.55906514 0.59476159 0.59630175 0.59865230 0.60013956
#> 368 0.33020146 0.33205059 0.39266337 0.39464335 0.39537585 0.39549758
#> 369 0.22223978 0.24684275 0.27948413 0.27981220 0.28827834 0.28857669
#> 370 0.43528170 0.43531603 0.44481567 0.44494973 0.44500028 0.44514401
#> 371 0.40293326 0.40334443 0.44604467 0.44637950 0.44640533 0.44649747
#> 372 0.33205768 0.33213562 0.37851658 0.37858426 0.38309953 0.38715620
#> 373 0.35029211 0.35050651 0.37641574 0.38088285 0.38092191 0.38247851
#> 374 0.34231163 0.34703275 0.36949024 0.37613175 0.37666338 0.38082607
#> 375 0.16811020 0.16823873 0.21547643 0.23507015 0.23599227 0.24762621
#> 392 0.29345023 0.32220777 0.32317845 0.32790161 0.33726854 0.33976964
#> 393 0.21740951 0.21804857 0.21882260 0.25857871 0.25858554 0.25873312
#> 394 0.08461823 0.10133393 0.14021023 0.16382766 0.16948255 0.16991000
#> 395 0.41037284 0.41296633 0.43495274 0.44363963 0.44405614 0.44417411
#> 396 0.25320866 0.25541095 0.25815382 0.27056257 0.27230304 0.27242041
#> 397 0.30346961 0.30912715 0.31259386 0.32560528 0.32724402 0.33180471
#> 398 0.36949149 0.37498196 0.38580450 0.41631821 0.41750001 0.41883368
#> 399 0.10598414 0.10612020 0.14658699 0.14922098 0.14933880 0.15992840
#> 400 0.32645302 0.32955374 0.34251834 0.34254568 0.34572490 0.35354217
#> 418 0.43107835 0.45634525 0.45639356 0.46375685 0.47011461 0.47123066
#> 419 0.36266612 0.37081411 0.37771556 0.37884906 0.39811537 0.39917788
#> 420 0.45390393 0.48802268 0.50060002 0.50161592 0.50515672 0.50610174
#> 421 0.39183891 0.42575723 0.42594076 0.42670603 0.42941924 0.43069722
#> 422 0.45123105 0.49035369 0.49082027 0.49211482 0.50461976 0.50480499
#> 423 0.50189414 0.53916971 0.54315909 0.55246848 0.55451037 0.55456456
#> 424 0.23262220 0.26536226 0.27605736 0.27742775 0.28687213 0.30327191
#> 425 0.56944568 0.59666531 0.60056342 0.60273573 0.60359290 0.60448261
#> 444 0.16688636 0.18192790 0.19814729 0.20512029 0.21006177 0.21008780
#> 445 0.56486673 0.56563797 0.57784564 0.58003977 0.58029914 0.58030108
#> 446 0.12576461 0.12687996 0.12831593 0.14031258 0.14324948 0.14325140
#> 447 0.32869051 0.33138705 0.33257094 0.33582555 0.33760210 0.34073331
#> 448 0.39311844 0.39594803 0.40250396 0.40293723 0.40401648 0.40476836
#> 449 0.28089145 0.28175589 0.29365223 0.31352727 0.31355343 0.32041717
#> 450 0.32018550 0.32144424 0.33006772 0.35602064 0.35884049 0.36471523
#> 470 0.46863481 0.47959858 0.48123351 0.48124067 0.48530066 0.48530403
#> 471 0.21455768 0.23236892 0.23983399 0.23986315 0.25106371 0.25109635
#> 472 0.26035867 0.28447343 0.28891727 0.28898786 0.28918611 0.29233322
#> 473 0.33403441 0.35482913 0.35488267 0.36205507 0.36970230 0.37031842
#> 474 0.16236538 0.17000199 0.17005439 0.17929745 0.18596390 0.19844706
#> 475 0.28036278 0.29795077 0.29797204 0.30863126 0.32153994 0.32744211
#> 496 0.61125268 0.61138546 0.61935002 0.61938272 0.61948965 0.61949370
#> 497 0.59921572 0.59934947 0.60608588 0.60609756 0.60831875 0.60995064
#> 498 0.58132919 0.58159932 0.58339978 0.58747594 0.58760184 0.58805812
#> 499 0.45636546 0.45927127 0.46063457 0.46445498 0.46464529 0.46951748
#> 500 0.47612619 0.47612748 0.47808274 0.48524010 0.48602153 0.49017794
#> 522 0.28219511 0.28333616 0.28347570 0.28365655 0.29066569 0.29441375
#> 523 0.35469532 0.35600127 0.35995548 0.36777779 0.36777837 0.36865454
#> 524 0.28091197 0.28358441 0.29062724 0.29543559 0.29647528 0.30356499
#> 525 0.37967708 0.37985524 0.38393361 0.39151428 0.39184269 0.39735205
#> 548 0.12702558 0.12712643 0.13462856 0.14963879 0.16278990 0.16419399
#> 549 0.32133433 0.32769334 0.33222726 0.33818743 0.34038766 0.35835658
#> 550 0.31916569 0.31946260 0.32288964 0.33435257 0.34367096 0.34420189
#> 574 0.40481226 0.41070215 0.41083344 0.43047076 0.43122302 0.44064296
#> 575 0.29155468 0.29216229 0.29217390 0.32940254 0.32976170 0.33224398
#> 600 0.42082007 0.42383613 0.42390754 0.42456319 0.42652357 0.44623536
#>       pcoa_10d tree_average
#> 2   0.29423834   0.32109423
#> 3   0.41188119   0.46338132
#> 4   0.37970684   0.36489558
#> 5   0.34458337   0.36489558
#> 6   0.38930266   0.36489558
#> 7   0.37165906   0.32109423
#> 8   0.36130480   0.33713555
#> 9   0.31884085   0.20080550
#> 10  0.18079936   0.11842324
#> 11  0.30647948   0.20080550
#> 12  0.44459287   0.28433580
#> 13  0.50911558   0.46338132
#> 14  0.48164955   0.46338132
#> 15  0.35666205   0.36489558
#> 16  0.20472493   0.20080550
#> 17  0.36148728   0.46338132
#> 18  0.22094147   0.32109423
#> 19  0.12703348   0.06793207
#> 20  0.50266415   0.46338132
#> 21  0.25869067   0.32109423
#> 22  0.31319125   0.36489558
#> 23  0.36774951   0.36489558
#> 24  0.24948754   0.20080550
#> 25  0.34802882   0.36489558
#> 28  0.54862680   0.46338132
#> 29  0.32148398   0.36489558
#> 30  0.29405160   0.36489558
#> 31  0.41687855   0.36489558
#> 32  0.22424968   0.14028171
#> 33  0.38448766   0.33713555
#> 34  0.40057174   0.32109423
#> 35  0.32326787   0.32109423
#> 36  0.39429629   0.32109423
#> 37  0.47867862   0.32109423
#> 38  0.45267899   0.46338132
#> 39  0.47996101   0.46338132
#> 40  0.44732886   0.36489558
#> 41  0.28399629   0.32109423
#> 42  0.46780957   0.46338132
#> 43  0.18711351   0.10427836
#> 44  0.28103374   0.32109423
#> 45  0.64779810   0.46338132
#> 46  0.09224177   0.02899600
#> 47  0.30445486   0.36489558
#> 48  0.38703108   0.36489558
#> 49  0.31658808   0.32109423
#> 50  0.39224950   0.36489558
#> 54  0.66727283   0.46338132
#> 55  0.62516474   0.46338132
#> 56  0.67701573   0.46338132
#> 57  0.61045582   0.46338132
#> 58  0.56785153   0.46338132
#> 59  0.31990439   0.46338132
#> 60  0.45212757   0.46338132
#> 61  0.39028522   0.46338132
#> 62  0.39064631   0.46338132
#> 63  0.44189984   0.37606375
#> 64  0.30690239   0.37606375
#> 65  0.51047459   0.46338132
#> 66  0.32534562   0.46338132
#> 67  0.36508616   0.37606375
#> 68  0.47983616   0.46338132
#> 69  0.38572507   0.46338132
#> 70  0.33062944   0.15956682
#> 71  0.51888623   0.46338132
#> 72  0.59629464   0.46338132
#> 73  0.65323799   0.46338132
#> 74  0.34400273   0.46338132
#> 75  0.55112618   0.46338132
#> 80  0.13692821   0.17152720
#> 81  0.25463361   0.14929723
#> 82  0.37342326   0.36489558
#> 83  0.43393579   0.36489558
#> 84  0.49226671   0.36489558
#> 85  0.41604075   0.36489558
#> 86  0.47505884   0.36489558
#> 87  0.55840948   0.36489558
#> 88  0.62764727   0.46338132
#> 89  0.71062150   0.46338132
#> 90  0.43335529   0.29820443
#> 91  0.38798266   0.36489558
#> 92  0.56877635   0.46338132
#> 93  0.41123644   0.36489558
#> 94  0.37000105   0.36489558
#> 95  0.67519324   0.46338132
#> 96  0.34270027   0.36489558
#> 97  0.18733686   0.17152720
#> 98  0.24667583   0.17152720
#> 99  0.41216382   0.36489558
#> 100 0.35220385   0.29820443
#> 106 0.30025066   0.17152720
#> 107 0.36807919   0.36489558
#> 108 0.40458092   0.36489558
#> 109 0.44003875   0.36489558
#> 110 0.36993854   0.36489558
#> 111 0.39650397   0.36489558
#> 112 0.47443522   0.36489558
#> 113 0.58258015   0.46338132
#> 114 0.66351035   0.46338132
#> 115 0.39649401   0.29820443
#> 116 0.34685995   0.36489558
#> 117 0.52530475   0.46338132
#> 118 0.36594475   0.36489558
#> 119 0.32724339   0.36489558
#> 120 0.63567854   0.46338132
#> 121 0.31017309   0.36489558
#> 122 0.09063639   0.04387834
#> 123 0.18686313   0.08973041
#> 124 0.37444541   0.36489558
#> 125 0.33802397   0.29820443
#> 132 0.40217677   0.36489558
#> 133 0.50495365   0.36489558
#> 134 0.58097741   0.36489558
#> 135 0.48331465   0.36489558
#> 136 0.58088787   0.36489558
#> 137 0.64962920   0.36489558
#> 138 0.67686222   0.46338132
#> 139 0.74744401   0.46338132
#> 140 0.45201612   0.29820443
#> 141 0.46870102   0.36489558
#> 142 0.56360941   0.46338132
#> 143 0.45430312   0.36489558
#> 144 0.41828180   0.36489558
#> 145 0.61692927   0.46338132
#> 146 0.39659120   0.36489558
#> 147 0.29800686   0.17152720
#> 148 0.22262597   0.17152720
#> 149 0.48151996   0.36489558
#> 150 0.36608994   0.29820443
#> 158 0.39354680   0.33713555
#> 159 0.53684123   0.32109423
#> 160 0.44636660   0.32109423
#> 161 0.53131298   0.32109423
#> 162 0.58666496   0.32109423
#> 163 0.45020980   0.46338132
#> 164 0.56021819   0.46338132
#> 165 0.55912954   0.36489558
#> 166 0.41406539   0.32109423
#> 167 0.49491921   0.46338132
#> 168 0.29364380   0.14028171
#> 169 0.37865786   0.32109423
#> 170 0.71807068   0.46338132
#> 171 0.19828444   0.14028171
#> 172 0.36383082   0.36489558
#> 173 0.42544180   0.36489558
#> 174 0.41667126   0.32109423
#> 175 0.50346927   0.36489558
#> 184 0.41124397   0.33713555
#> 185 0.28942214   0.33713555
#> 186 0.42959108   0.33713555
#> 187 0.46394973   0.33713555
#> 188 0.42626521   0.46338132
#> 189 0.54917763   0.46338132
#> 190 0.46412769   0.36489558
#> 191 0.40963606   0.33713555
#> 192 0.46125928   0.46338132
#> 193 0.42490659   0.33713555
#> 194 0.32331366   0.33713555
#> 195 0.69159691   0.46338132
#> 196 0.38020155   0.33713555
#> 197 0.40338782   0.36489558
#> 198 0.48450555   0.36489558
#> 199 0.31742348   0.33713555
#> 200 0.54740913   0.36489558
#> 210 0.25884346   0.20080550
#> 211 0.19132078   0.11692821
#> 212 0.29174649   0.28433580
#> 213 0.44710697   0.46338132
#> 214 0.35287519   0.46338132
#> 215 0.29675314   0.36489558
#> 216 0.20190587   0.13424457
#> 217 0.44584205   0.46338132
#> 218 0.36010660   0.32109423
#> 219 0.24057074   0.20080550
#> 220 0.42988528   0.46338132
#> 221 0.41087887   0.32109423
#> 222 0.43447829   0.36489558
#> 223 0.50810750   0.36489558
#> 224 0.22204383   0.13424457
#> 225 0.39260683   0.36489558
#> 236 0.27460337   0.20080550
#> 237 0.40855752   0.28433580
#> 238 0.51195587   0.46338132
#> 239 0.47368184   0.46338132
#> 240 0.33663591   0.36489558
#> 241 0.24321750   0.20080550
#> 242 0.40985676   0.46338132
#> 243 0.29402740   0.32109423
#> 244 0.16884753   0.11842324
#> 245 0.56147535   0.46338132
#> 246 0.32203906   0.32109423
#> 247 0.36147845   0.36489558
#> 248 0.43481511   0.36489558
#> 249 0.24529548   0.20080550
#> 250 0.38776823   0.36489558
#> 262 0.23291911   0.28433580
#> 263 0.46683789   0.46338132
#> 264 0.44358344   0.46338132
#> 265 0.37164840   0.36489558
#> 266 0.20630116   0.13424457
#> 267 0.42684569   0.46338132
#> 268 0.34853816   0.32109423
#> 269 0.24914671   0.20080550
#> 270 0.49588958   0.46338132
#> 271 0.39531212   0.32109423
#> 272 0.37354500   0.36489558
#> 273 0.46368202   0.36489558
#> 274 0.22381322   0.13424457
#> 275 0.44985246   0.36489558
#> 288 0.37675107   0.46338132
#> 289 0.43555214   0.46338132
#> 290 0.48801616   0.36489558
#> 291 0.31786180   0.28433580
#> 292 0.37512492   0.46338132
#> 293 0.48674012   0.32109423
#> 294 0.38964206   0.28433580
#> 295 0.49008457   0.46338132
#> 296 0.47750939   0.32109423
#> 297 0.45987455   0.36489558
#> 298 0.53839269   0.36489558
#> 299 0.26310230   0.28433580
#> 300 0.57158383   0.36489558
#> 314 0.33387440   0.16014333
#> 315 0.61290371   0.46338132
#> 316 0.43248323   0.46338132
#> 317 0.40604020   0.36833202
#> 318 0.47751035   0.46338132
#> 319 0.46220510   0.46338132
#> 320 0.61565461   0.37606375
#> 321 0.43863440   0.46338132
#> 322 0.57112678   0.46338132
#> 323 0.64386972   0.46338132
#> 324 0.34015697   0.46338132
#> 325 0.68302280   0.46338132
#> 340 0.54957464   0.46338132
#> 341 0.39009863   0.46338132
#> 342 0.46795819   0.36833202
#> 343 0.43076857   0.46338132
#> 344 0.43713721   0.46338132
#> 345 0.51408242   0.37606375
#> 346 0.47189154   0.46338132
#> 347 0.65193522   0.46338132
#> 348 0.71372890   0.46338132
#> 349 0.38376247   0.46338132
#> 350 0.60657637   0.46338132
#> 366 0.36341858   0.36489558
#> 367 0.60028199   0.46338132
#> 368 0.39728743   0.36489558
#> 369 0.29045420   0.36489558
#> 370 0.44524299   0.46338132
#> 371 0.44916020   0.36489558
#> 372 0.38799756   0.29820443
#> 373 0.38322625   0.29820443
#> 374 0.38093659   0.36489558
#> 375 0.24762622   0.10600649
#> 392 0.34015778   0.46338132
#> 393 0.25886776   0.32109423
#> 394 0.18096902   0.20080550
#> 395 0.44462115   0.46338132
#> 396 0.28341012   0.32109423
#> 397 0.33630916   0.36489558
#> 398 0.41886840   0.36489558
#> 399 0.16449638   0.10966533
#> 400 0.35476000   0.36489558
#> 418 0.47187198   0.46338132
#> 419 0.40182504   0.46338132
#> 420 0.50611518   0.37606375
#> 421 0.43512423   0.46338132
#> 422 0.50627966   0.46338132
#> 423 0.55467074   0.46338132
#> 424 0.30408710   0.46338132
#> 425 0.60462507   0.46338132
#> 444 0.22164616   0.32109423
#> 445 0.58098983   0.46338132
#> 446 0.16734396   0.10427836
#> 447 0.34654626   0.36489558
#> 448 0.40500154   0.36489558
#> 449 0.32382336   0.32109423
#> 450 0.36666040   0.36489558
#> 470 0.48714919   0.46338132
#> 471 0.25159697   0.32109423
#> 472 0.29242773   0.36489558
#> 473 0.37466427   0.36489558
#> 474 0.19986872   0.20080550
#> 475 0.32909307   0.36489558
#> 496 0.62222374   0.46338132
#> 497 0.61094974   0.46338132
#> 498 0.58823785   0.46338132
#> 499 0.46988430   0.46338132
#> 500 0.49026877   0.46338132
#> 522 0.29533410   0.36489558
#> 523 0.37576680   0.36489558
#> 524 0.30614478   0.32109423
#> 525 0.40033583   0.36489558
#> 548 0.17148334   0.08973041
#> 549 0.35872994   0.36489558
#> 550 0.34514436   0.29820443
#> 574 0.44188466   0.36489558
#> 575 0.33310802   0.29820443
#> 600 0.44632867   0.36489558
#> 
#> 
#> $details_deviation
#> $details_deviation$dev_distsp
#>              sp.x          sp.y       pcoa_1d       pcoa_2d       pcoa_3d
#> 2         apricot         apple -0.1634793627 -0.0613131726 -0.0567510971
#> 3          banana         apple -0.0131320500 -0.0047392768 -0.0037427035
#> 4         currant         apple -0.1010576998 -0.1006330170 -0.0813228267
#> 5      blackberry         apple -0.1078357918 -0.1075952587 -0.0657691608
#> 6       blueberry         apple -0.0423032775 -0.0375635362 -0.0370317817
#> 7          cherry         apple -0.1616933564  0.0044655132  0.0185309251
#> 8           grape         apple -0.3014574226 -0.1180817603 -0.0492048338
#> 9      grapefruit         apple -0.0782134735 -0.0620468549 -0.0406330304
#> 10      kiwifruit         apple -0.1182215595 -0.1151277801 -0.0521533878
#> 11          lemon         apple -0.1163257141 -0.1107436173 -0.0223964456
#> 12           lime         apple -0.1717551260 -0.1715743017 -0.0544854936
#> 13         litchi         apple -0.2208832159 -0.1118805638 -0.1066939412
#> 14          mango         apple -0.0487058101 -0.0364587186 -0.0208028301
#> 15          melon         apple -0.2616967886 -0.0665607673 -0.0638025397
#> 16         orange         apple -0.0875671626 -0.0866583269 -0.0779119903
#> 17  passion_fruit         apple -0.1096418149 -0.0838538643 -0.0646337510
#> 18          peach         apple -0.1348560405 -0.0634209972 -0.0167038995
#> 19           pear         apple -0.0406116659 -0.0398174530 -0.0264415399
#> 20      pineapple         apple -0.1919552128 -0.0362306401 -0.0345022788
#> 21           plum         apple -0.1480900385 -0.0408337829 -0.0309493853
#> 22      raspberry         apple -0.0967829650 -0.0967745107 -0.0553955783
#> 23     strawberry         apple -0.0862982218 -0.0736475401 -0.0584288511
#> 24      tangerine         apple -0.0825254764 -0.0733131341 -0.0300916371
#> 25    water_melon         apple -0.1703111415 -0.0530135968 -0.0415572094
#> 28         banana       apricot -0.0857023219 -0.0242932667 -0.0241109329
#> 29        currant       apricot -0.1042450037 -0.0545054018 -0.0146009012
#> 30     blackberry       apricot -0.1110230958 -0.0436697614  0.0248648307
#> 31      blueberry       apricot -0.1356171216 -0.0461783026 -0.0453007506
#> 32         cherry       apricot -0.0730266810 -0.0090297993  0.0019372342
#> 33          grape       apricot -0.2956743988 -0.2906455598 -0.1377553518
#> 34     grapefruit       apricot -0.1183953004 -0.0243571666  0.0006936907
#> 35      kiwifruit       apricot -0.1584033864 -0.0853719907 -0.0419028383
#> 36          lemon       apricot -0.1565075411 -0.0769666856 -0.0021222401
#> 37           lime       apricot -0.1286036196 -0.0935940794  0.0201976926
#> 38         litchi       apricot -0.0373928817 -0.0221055863 -0.0075622973
#> 39          mango       apricot  0.0696330090  0.0740771360  0.0795749880
#> 40          melon       apricot -0.3439820414 -0.0187968599 -0.0114241535
#> 41         orange       apricot -0.0732035350  0.0003625967  0.0138623513
#> 42  passion_fruit       apricot -0.1554444099 -0.1497362756 -0.1188010163
#> 43          peach       apricot -0.0480830291 -0.0111865203  0.0054971105
#> 44           pear       apricot -0.1167117412 -0.0229469924 -0.0105328150
#> 45      pineapple       apricot -0.2945870897 -0.0399771690 -0.0399703870
#> 46           plum       apricot -0.0153893242 -0.0149825960 -0.0060163710
#> 47      raspberry       apricot -0.1423945114 -0.0604310211  0.0049237280
#> 48     strawberry       apricot -0.1561521924 -0.0310780438 -0.0034079774
#> 49      tangerine       apricot -0.0657018088 -0.0364410525  0.0117241279
#> 50    water_melon       apricot -0.2568317788  0.0651181273  0.0668086545
#> 54        currant        banana -0.0232806589 -0.0163820775 -0.0030818734
#> 55     blackberry        banana -0.0300587510 -0.0232900603  0.0002809691
#> 56      blueberry        banana -0.0130110851 -0.0124973002 -0.0124278354
#> 57         cherry        banana -0.0660341977  0.0420972235  0.0455975807
#> 58          grape        banana -0.1674144255 -0.0804373324 -0.0372737598
#> 59     grapefruit        banana -0.1024010940 -0.1024008480 -0.0646444426
#> 60      kiwifruit        banana -0.0946385403 -0.0816225538 -0.0610000317
#> 61          lemon        banana -0.1340038606 -0.1306800188 -0.0487136663
#> 62           lime        banana -0.2208557785 -0.1930523006 -0.0171950080
#> 63         litchi        banana -0.2912331831 -0.0275701602 -0.0163747357
#> 64          mango        banana -0.2022404267 -0.0462949627 -0.0269137086
#> 65          melon        banana -0.0864450002 -0.0626896105 -0.0580126769
#> 66         orange        banana -0.0147881642 -0.0069464042  0.0015288165
#> 67  passion_fruit        banana -0.1545003362 -0.0712841793 -0.0422261591
#> 68          peach        banana -0.0813214239 -0.0476578445 -0.0410890920
#> 69           pear        banana -0.0033340150  0.0042738603  0.0091095536
#> 70      pineapple        banana -0.0557616968  0.1063538199  0.1065677751
#> 71           plum        banana -0.0703129977 -0.0088275768 -0.0077939136
#> 72      raspberry        banana -0.0190059241 -0.0132293892  0.0091482188
#> 73     strawberry        banana -0.0327636052 -0.0327607052 -0.0215307530
#> 74      tangerine        banana -0.0705055635 -0.0383462584 -0.0023877290
#> 75    water_melon        banana -0.0092007673  0.0084254471  0.0107939978
#> 80     blackberry       currant -0.0502274025 -0.0500022593 -0.0389149021
#> 81      blueberry       currant -0.1287289277 -0.0786732777 -0.0090084589
#> 82         cherry       currant -0.1616628782 -0.0164758571  0.0438968236
#> 83          grape       currant -0.1499194025 -0.0975652203 -0.0894623777
#> 84     grapefruit       currant -0.0131622822 -0.0035252359 -0.0034422918
#> 85      kiwifruit       currant -0.0902582803 -0.0902582026 -0.0900754121
#> 86          lemon       currant -0.0512745228 -0.0475976058 -0.0401120522
#> 87           lime       currant -0.0233706014 -0.0233370624  0.0013152483
#> 88         litchi       currant -0.0583045521 -0.0027321360 -0.0013323560
#> 89          mango       currant -0.0346119947 -0.0294009509  0.0052273453
#> 90          melon       currant -0.1388558607 -0.0480657892 -0.0410986853
#> 91         orange       currant -0.0084924947 -0.0073959172 -0.0021501549
#> 92  passion_fruit       currant -0.0768611308 -0.0683091410 -0.0682377642
#> 93          peach       currant -0.1328683259 -0.1236499732 -0.0536370397
#> 94           pear       currant -0.0542900783 -0.0534915442 -0.0446038482
#> 95      pineapple       currant -0.1210543157 -0.0215406693 -0.0057361827
#> 96           plum       currant -0.1196343279 -0.0745078580 -0.0263173850
#> 97      raspberry       currant -0.0722692403 -0.0711190425 -0.0658724616
#> 98     strawberry       currant -0.1826471459 -0.1050775552 -0.1041838161
#> 99      tangerine       currant -0.0032801459 -0.0018235003 -0.0016200373
#> 100   water_melon       currant -0.1494334270 -0.0212950301  0.0341925798
#> 106     blueberry    blackberry -0.1425926938 -0.1198579992 -0.0319868888
#> 107        cherry    blackberry -0.1846025864 -0.0029153831  0.0834226480
#> 108         grape    blackberry -0.1566974945 -0.0904380138 -0.0881675515
#> 109    grapefruit    blackberry -0.0185667479 -0.0087531920 -0.0062425888
#> 110     kiwifruit    blackberry -0.0970363723 -0.0969911261 -0.0925690637
#> 111         lemon    blackberry -0.0010471203  0.0024633279  0.0049577139
#> 112          lime    blackberry  0.0268568011  0.0268576408  0.0433965702
#> 113        litchi    blackberry -0.0650826441 -0.0016304175  0.0038078583
#> 114         mango    blackberry -0.0413900867 -0.0350098427  0.0157181544
#> 115         melon    blackberry -0.1456339528 -0.0417186480 -0.0235214622
#> 116        orange    blackberry -0.0152705868 -0.0143936233  0.0009361898
#> 117 passion_fruit    blackberry -0.0836392229 -0.0729129995 -0.0705637247
#> 118         peach    blackberry -0.1396464180 -0.1263748437 -0.0178358171
#> 119          pear    blackberry -0.0610681703 -0.0604846474 -0.0367793261
#> 120     pineapple    blackberry -0.1278324077 -0.0217966063  0.0050104300
#> 121          plum    blackberry -0.1264124200 -0.0656793526  0.0133150833
#> 122     raspberry    blackberry -0.0220418378 -0.0205056631 -0.0193202051
#> 123    strawberry    blackberry -0.0635398526 -0.0013808365  0.0106468099
#> 124     tangerine    blackberry -0.0100582379 -0.0078413278 -0.0070103500
#> 125   water_melon    blackberry -0.1562115190 -0.0018367241  0.0823861359
#> 132        cherry     blueberry -0.1271508801  0.0667314715  0.0743952001
#> 133         grape     blueberry -0.1638922529 -0.0749047500 -0.0277515215
#> 134    grapefruit     blueberry -0.0873457554 -0.0866216698 -0.0734037705
#> 135     kiwifruit     blueberry -0.1273538414 -0.1203937378 -0.1022813840
#> 136         lemon     blueberry -0.1254579960 -0.1252654078 -0.0807582478
#> 137          lime     blueberry -0.0975540745 -0.0939381975 -0.0239826308
#> 138        litchi     blueberry -0.0722774025  0.0120965610  0.0171258743
#> 139         mango     blueberry -0.0485848452 -0.0320055949 -0.0255240022
#> 140         melon     blueberry -0.0907102741 -0.0450716630 -0.0407286918
#> 141        orange     blueberry -0.0421539899 -0.0401296328 -0.0355770832
#> 142 passion_fruit     blueberry -0.0181067085  0.0066244862  0.0188959426
#> 143         peach     blueberry -0.1352153022 -0.1030979994 -0.0920402991
#> 144          pear     blueberry -0.0497439477 -0.0464304175 -0.0428117339
#> 145     pineapple     blueberry  0.0211347275  0.0845822197  0.0849161857
#> 146          plum     blueberry -0.1111962560 -0.0279512489 -0.0252093743
#> 147     raspberry     blueberry -0.1222102892 -0.1080164379 -0.0322289736
#> 148    strawberry     blueberry -0.0972346594 -0.0886234114 -0.0189774446
#> 149     tangerine     blueberry -0.0172529963 -0.0049974611  0.0155522376
#> 150   water_melon     blueberry -0.1135172775 -0.0470372950 -0.0394305372
#> 158         grape        cherry -0.2034080395 -0.1819563757 -0.0216705651
#> 159    grapefruit        cherry -0.1914219815 -0.0392690809 -0.0029767167
#> 160     kiwifruit        cherry -0.2314300675 -0.0957104915 -0.0427964444
#> 161         lemon        cherry -0.2295342221 -0.0914128522 -0.0070527005
#> 162          lime        cherry -0.2016303006 -0.1213014649  0.0031759239
#> 163        litchi        cherry  0.0548734776  0.0548761653  0.0840397080
#> 164         mango        cherry  0.0785660350  0.1050981332  0.1054867477
#> 165         melon        cherry -0.4170087224 -0.0319605478 -0.0159188158
#> 166        orange        cherry -0.1462302160 -0.0126209943  0.0117021164
#> 167 passion_fruit        cherry -0.1189528314 -0.0823886704 -0.0371119584
#> 168         peach        cherry -0.0726248617  0.0268893826  0.0270746275
#> 169          pear        cherry -0.1538201738  0.0035825380  0.0267583808
#> 170     pineapple        cherry -0.3051637208  0.0123794817  0.0144561003
#> 171          plum        cherry -0.0486058154  0.0112353493  0.0155654633
#> 172     raspberry        cherry -0.2078931939 -0.0025180762  0.0785077694
#> 173    strawberry        cherry -0.1812468345  0.0625514695  0.1055026953
#> 174     tangerine        cherry -0.1080508342 -0.0287779733  0.0311242144
#> 175   water_melon        cherry -0.4299062189  0.0405236865  0.0407749698
#> 184    grapefruit         grape -0.1776506941 -0.0389963195 -0.0292206691
#> 185     kiwifruit         grape -0.1469517094 -0.0082199132  0.0098665282
#> 186         lemon         grape -0.2561669751 -0.1299319274 -0.1298461535
#> 187          lime         grape -0.1878590133 -0.1286605997 -0.1161421041
#> 188        litchi         grape -0.0428479090 -0.0323309533 -0.0108948220
#> 189         mango         grape -0.0894420095 -0.0788520328  0.0141223634
#> 190         melon         grape -0.3721746529  0.0165680287  0.0403708538
#> 191        orange         grape -0.2758932721 -0.1502588494 -0.1135193000
#> 192 passion_fruit         grape -0.1981108370 -0.1817451705 -0.1677844751
#> 193         peach         grape -0.3988692584 -0.3021089636 -0.1117148453
#> 194          pear         grape -0.2208569673 -0.0555072487 -0.0081946104
#> 195     pineapple         grape -0.4126045546 -0.1000071709 -0.0617766611
#> 196          plum         grape -0.3045274988 -0.2969118822 -0.1189354807
#> 197     raspberry         grape -0.1698870919 -0.0921123859 -0.0886479777
#> 198    strawberry         grape -0.1836447730 -0.0639164533 -0.0540593268
#> 199     tangerine         grape -0.1346835872 -0.0740454082 -0.0651176166
#> 200   water_melon         grape -0.3615970866 -0.0682179928  0.0015414264
#> 210     kiwifruit    grapefruit -0.0711030251 -0.0441897349 -0.0441808613
#> 211         lemon    grapefruit -0.0720068070 -0.0574729403 -0.0112342765
#> 212          lime    grapefruit -0.1768749858 -0.1277915443 -0.0334733482
#> 213        litchi    grapefruit -0.2965158962 -0.0276659703 -0.0262255156
#> 214         mango    grapefruit -0.0571667732  0.0095175572  0.0806111285
#> 215         melon    grapefruit  0.0033645187  0.0448281429  0.0520934779
#> 216        orange    grapefruit -0.0451917655 -0.0237047061 -0.0115933960
#> 217 passion_fruit    grapefruit -0.3872946972 -0.2206505970 -0.2206499567
#> 218         peach    grapefruit -0.1187971198 -0.0602496655  0.0060147092
#> 219          pear    grapefruit -0.0376018076 -0.0217725252 -0.0101165707
#> 220     pineapple    grapefruit -0.2931040437 -0.1118771239 -0.0795045717
#> 221          plum    grapefruit -0.1428161660 -0.0473697432 -0.0143895225
#> 222     raspberry    grapefruit -0.0317563453 -0.0232212604 -0.0213712782
#> 223    strawberry    grapefruit -0.0455140263 -0.0455111169 -0.0454878345
#> 224     tangerine    grapefruit -0.0833711473 -0.0010010376  0.0006271003
#> 225   water_melon    grapefruit  0.0079897040  0.0348185587  0.0753277173
#> 236         lemon     kiwifruit -0.1092152657 -0.0957822071 -0.0657676637
#> 237          lime     kiwifruit -0.1646446776 -0.1645673882 -0.1026883259
#> 238        litchi     kiwifruit -0.2458240496 -0.1371476849 -0.1358194774
#> 239         mango     kiwifruit -0.0736466437 -0.0639642335 -0.0092587302
#> 240         melon     kiwifruit -0.1862509044 -0.0049093272  0.0021948199
#> 241        orange     kiwifruit -0.1323756286 -0.1261557728 -0.1084963855
#> 242 passion_fruit     kiwifruit -0.1850876990 -0.1622803795 -0.1622669962
#> 243         peach     kiwifruit -0.1725414695 -0.1362741422 -0.0092101013
#> 244          pear     kiwifruit -0.0776098936 -0.0618923103 -0.0179709910
#> 245     pineapple     kiwifruit -0.2853414899 -0.1035086810 -0.0828889510
#> 246          plum     kiwifruit -0.1828242521 -0.1068380682 -0.0496642265
#> 247     raspberry     kiwifruit -0.1102259697 -0.1098479850 -0.1064294583
#> 248    strawberry     kiwifruit -0.1239836508 -0.1083939369 -0.1083815700
#> 249     tangerine     kiwifruit -0.1074663101 -0.1013098765 -0.0985915281
#> 250   water_melon     kiwifruit -0.1595117220 -0.0448971383  0.0023158473
#> 262          lime         lemon -0.0554294119 -0.0440033211 -0.0174707022
#> 263        litchi         lemon -0.2851893699 -0.0839240827 -0.0631119153
#> 264         mango         lemon -0.1130119641 -0.0771103121  0.0491707097
#> 265         melon         lemon -0.1584850956 -0.0807337454 -0.0385406786
#> 266        orange         lemon -0.0833040061 -0.0759364321  0.0221946807
#> 267 passion_fruit         lemon -0.3355641305 -0.2406276248 -0.2182694085
#> 268         peach         lemon -0.1569093604 -0.1128743083  0.0353315784
#> 269          pear         lemon -0.0757140482 -0.0708962607  0.0042115984
#> 270     pineapple         lemon -0.3247068103 -0.1364359942 -0.0658292112
#> 271          plum         lemon -0.1809284067 -0.0999191547 -0.0109778229
#> 272     raspberry         lemon -0.0032477068 -0.0006499275  0.0030460928
#> 273    strawberry         lemon -0.0282233365 -0.0261011387 -0.0153774821
#> 274     tangerine         lemon -0.1214833879 -0.0446382132 -0.0229531023
#> 275   water_melon         lemon -0.1376982942 -0.0877297884  0.0037470198
#> 288        litchi          lime -0.2483039141 -0.0217930487  0.0468610036
#> 289         mango          lime -0.1756214577 -0.1508903151  0.0913165327
#> 290         melon          lime -0.1735104671 -0.0896162194 -0.0202748993
#> 291        orange          lime -0.1387334180 -0.1369989612 -0.0171901107
#> 292 passion_fruit          lime -0.2665158530 -0.2004120061 -0.1061164965
#> 293         peach          lime -0.2123387723 -0.2014167640 -0.0149157567
#> 294          pear          lime -0.1311434601 -0.1305420922 -0.0272277478
#> 295     pineapple          lime -0.4115587282 -0.0987718372  0.0171440540
#> 296          plum          lime -0.1530244852 -0.1181001709  0.0115767364
#> 297     raspberry          lime  0.0246562147  0.0247211904  0.0442083322
#> 298    strawberry          lime -0.0003194151  0.0078614002  0.0372874012
#> 299     tangerine          lime -0.0935794665 -0.0867387481 -0.0164402005
#> 300   water_melon          lime -0.1688852818 -0.1079670551  0.0100307363
#> 314         mango        litchi -0.0465941004  0.0425935670  0.1077304908
#> 315         melon        litchi -0.2805598024 -0.0091424317 -0.0084356276
#> 316        orange        litchi -0.2164787240 -0.0624470945 -0.0617668873
#> 317 passion_fruit        litchi -0.2970999869 -0.1907680478 -0.1877670736
#> 318         peach        litchi -0.1405877413 -0.0873849630 -0.0496809370
#> 319          pear        litchi -0.1949235647 -0.0725905586 -0.0711589856
#> 320     pineapple        litchi -0.5135128627  0.0637928780  0.0723343382
#> 321          plum        litchi -0.0462459817 -0.0293441938 -0.0087295223
#> 322     raspberry        litchi -0.0782722415 -0.0080886709 -0.0035309912
#> 323    strawberry        litchi -0.0920299226  0.0129486694  0.0135573477
#> 324     tangerine        litchi -0.1383577395 -0.0220291761 -0.0165950275
#> 325   water_melon        litchi -0.2699822362 -0.0477992554 -0.0324275118
#> 340         melon         mango -0.0412106794  0.0657895239  0.0865077974
#> 341        orange         mango -0.0261195000 -0.0056275802  0.0305854672
#> 342 passion_fruit         mango -0.1900740962 -0.1899212168 -0.0957518626
#> 343         peach         mango  0.0497714827  0.0499996401  0.0501674011
#> 344          pear         mango -0.0227461588 -0.0077762284  0.0188106524
#> 345     pineapple         mango -0.2580021235  0.0558671409  0.0621341082
#> 346          plum         mango  0.0607799090  0.0648937934  0.0680529849
#> 347     raspberry         mango -0.0545796842 -0.0466460735  0.0027453467
#> 348    strawberry         mango -0.0683373652 -0.0437523447 -0.0131284688
#> 349     tangerine         mango -0.0495136670 -0.0449576967  0.0408009408
#> 350   water_melon         mango  0.0037103212  0.0893564581  0.0904279946
#> 366        orange         melon -0.2075154219 -0.0900011804 -0.0898954652
#> 367 passion_fruit         melon -0.3713386034 -0.2037118201 -0.1990172341
#> 368         peach         melon -0.3298710510 -0.0377726845 -0.0085874499
#> 369          pear         melon -0.1664443186 -0.0007544235 -0.0007145660
#> 370     pineapple         melon -0.0549257276 -0.0347226196 -0.0260552069
#> 371          plum         melon -0.3684029070 -0.0346912098 -0.0232303063
#> 372     raspberry         melon -0.1588235502 -0.0533585067 -0.0365142488
#> 373    strawberry         melon -0.0614701201 -0.0243536751 -0.0185412522
#> 374     tangerine         melon -0.1826061033 -0.0270876441 -0.0162623537
#> 375   water_melon         melon  0.0046251853  0.0046922425  0.0507853782
#> 392 passion_fruit        orange -0.2001868179 -0.1478148749 -0.1390260345
#> 393         peach        orange -0.1089575021 -0.0715811087 -0.0133868438
#> 394          pear        orange -0.0792787290 -0.0789659721 -0.0780444327
#> 395     pineapple        orange -0.2236729320 -0.0366574875 -0.0277580578
#> 396          plum        orange -0.0976244006 -0.0222391098 -0.0014300479
#> 397     raspberry        orange -0.0284601842 -0.0280716885 -0.0139766380
#> 398    strawberry        orange -0.0422178652 -0.0356935385 -0.0323076676
#> 399     tangerine        orange -0.0781209845 -0.0400577149 -0.0076529300
#> 400   water_melon        orange -0.0777459364 -0.0030403701  0.0193293975
#> 418         peach passion_fruit -0.2263160372 -0.2253011779 -0.1495314151
#> 419          pear passion_fruit -0.2048942848 -0.1720904684 -0.1621585973
#> 420     pineapple passion_fruit -0.3700492732  0.0052956504  0.0253514552
#> 421          plum passion_fruit -0.1400550857 -0.1347416264 -0.0937330859
#> 422     raspberry passion_fruit -0.0887480122 -0.0754017176 -0.0736992375
#> 423    strawberry passion_fruit -0.0621016528 -0.0253873805 -0.0253619862
#> 424     tangerine passion_fruit -0.1887325001 -0.1701867806 -0.1678855457
#> 425   water_melon passion_fruit -0.2456095220 -0.1192071653 -0.0881428395
#> 444          pear         peach -0.1123308433 -0.0538243938  0.0113213139
#> 445     pineapple         peach -0.2659637675 -0.0478806800 -0.0438545332
#> 446          plum         peach -0.0482614705 -0.0085348753  0.0008147577
#> 447     raspberry         peach -0.1285935911 -0.1104562209 -0.0017350935
#> 448    strawberry         peach -0.1181088479 -0.0662003872 -0.0082466993
#> 449     tangerine         peach -0.1688966685 -0.1645280209 -0.0479623498
#> 450   water_melon         peach -0.2021217676 -0.0019419155 -0.0010125348
#> 470     pineapple          pear -0.1940369646 -0.0328245368 -0.0269012183
#> 471          plum          pear -0.1052143584 -0.0075938212  0.0125678879
#> 472     raspberry          pear -0.0500153435 -0.0498528276 -0.0272612427
#> 473    strawberry          pear -0.0637730246 -0.0541107915 -0.0478960196
#> 474     tangerine          pear -0.0565658252 -0.0419158097 -0.0130401894
#> 475   water_melon          pear -0.1235435200 -0.0233696848 -0.0021865076
#> 496          plum     pineapple -0.2808003296 -0.0229906606 -0.0226081209
#> 497     raspberry     pineapple -0.1167795809 -0.0115141420  0.0140756941
#> 498    strawberry     pineapple  0.0048162734  0.0601871544  0.0746508002
#> 499     tangerine     pineapple -0.2779209674 -0.0150749859  0.0139192041
#> 500   water_melon     pineapple  0.0041366871  0.0198892254  0.0217218905
#> 522     raspberry          plum -0.1153595931 -0.0413685884  0.0347919512
#> 523    strawberry          plum -0.1291172742 -0.0122414327  0.0226431675
#> 524     tangerine          plum -0.0745549089 -0.0451737424  0.0164526336
#> 525   water_melon          plum -0.2722211030  0.0350499062  0.0354583298
#> 548    strawberry     raspberry -0.0249756298  0.0181488565  0.0267160315
#> 549     tangerine     raspberry -0.0232478353 -0.0197593321 -0.0193459093
#> 550   water_melon     raspberry -0.2057647528 -0.0442244406  0.0387933564
#> 574     tangerine    strawberry -0.0370055164 -0.0152975875 -0.0144500881
#> 575   water_melon    strawberry -0.1084113227 -0.0494145856  0.0172473566
#> 600   water_melon     tangerine -0.1558669209 -0.0471430131  0.0004666108
#>           pcoa_4d       pcoa_5d       pcoa_6d       pcoa_7d      pcoa_8d
#> 2   -0.0528852420  0.0333988631  0.0339991097  0.0421872222 0.0530128968
#> 3    0.0072900396  0.0197260797  0.0202440801  0.0285611434 0.0330724861
#> 4   -0.0772269489 -0.0237868122 -0.0219995315  0.0108357707 0.0108467623
#> 5   -0.0635186776  0.0125484403  0.0125868642  0.0187303827 0.0197521546
#> 6    0.0199440721  0.0351629611  0.0524944227  0.0559045095 0.0570411048
#> 7    0.0281458699  0.0731641406  0.0783982461  0.0800797013 0.0806243911
#> 8   -0.0457945661 -0.0416295350  0.0484971172  0.0508539651 0.0513522627
#> 9   -0.0156670849  0.0373972686  0.0427345586  0.0562298325 0.0569609754
#> 10  -0.0017553655 -0.0017506389  0.0069481718  0.0103475530 0.0267954936
#> 11  -0.0083407947  0.0454018507  0.0513779010  0.0514114852 0.0516278013
#> 12  -0.0506307980  0.0293530414  0.0295385242  0.0298783547 0.0336665740
#> 13  -0.0909142259 -0.0034354962  0.0186229480  0.0193365308 0.0193481370
#> 14  -0.0179677430  0.0373474919  0.0429085187  0.0452516629 0.0498183062
#> 15  -0.0381959608  0.0241986355  0.0715810607  0.0715960261 0.0734237918
#> 16  -0.0779045654 -0.0122554013 -0.0047473389  0.0223458167 0.0225959398
#> 17   0.0342945375  0.0342953875  0.0352693917  0.0383269678 0.0495093071
#> 18   0.0035489276  0.0653982569  0.0678630003  0.0708725951 0.0711992901
#> 19  -0.0034553576  0.0333230276  0.0432566380  0.0445190412 0.0493499262
#> 20   0.0108899739  0.0448487287  0.0522024909  0.0524028442 0.0532592880
#> 21  -0.0309322278  0.0347371715  0.0351037985  0.0359590227 0.0393104590
#> 22  -0.0505415553  0.0177486307  0.0177823051  0.0179594778 0.0187928866
#> 23  -0.0355896728  0.0239622226  0.0261640935  0.0303233066 0.0326773201
#> 24  -0.0275826903  0.0134767349  0.0190861287  0.0307247672 0.0312625144
#> 25  -0.0371159138  0.0176539919  0.0204510351  0.0347866971 0.0395560055
#> 28  -0.0077329903  0.0020990516  0.0034955847  0.0036797634 0.0217944008
#> 29  -0.0015546801 -0.0011710552  0.0036198210  0.0160033955 0.0251620784
#> 30   0.0341744537  0.0344592035  0.0354113439  0.0354291618 0.0402678424
#> 31   0.0319843570  0.0437132823  0.0654164211  0.0657377916 0.0686685199
#> 32   0.0469869234  0.0487700617  0.0664300893  0.0691707809 0.0931576272
#> 33  -0.1377166522 -0.0915612653  0.0081726135  0.0091080295 0.0138727363
#> 34   0.0071120459  0.0082376801  0.0154202998  0.0161559835 0.0203086324
#> 35  -0.0299822755  0.0429829372  0.0507459073  0.0525913916 0.0526155478
#> 36   0.0005050269  0.0015667955  0.0039279806  0.0107799101 0.0213369999
#> 37   0.0289220112  0.0318920444  0.0319225550  0.0345753813 0.0350440886
#> 38   0.0209021675  0.0285272105  0.0608455748  0.0626952680 0.0690451296
#> 39   0.0796916045  0.0801801264  0.0888254092  0.0892566850 0.0894282379
#> 40  -0.0049695465 -0.0047243118  0.0408596167  0.0455206549 0.0475829584
#> 41   0.0168115778  0.0249382267  0.0269362847  0.0290669650 0.0376148564
#> 42  -0.0124184745  0.0321312713  0.0321967185  0.0326895430 0.0327990199
#> 43   0.0115392987  0.0209398961  0.0217281882  0.0574567505 0.0702756952
#> 44  -0.0102186893  0.0222489718  0.0305542383  0.0356185873 0.0591169970
#> 45   0.0060044575  0.0063003090  0.0147129587  0.0169191983 0.0188221997
#> 46   0.0174340360  0.0274630074  0.0276582284  0.0445927639 0.0538695247
#> 47   0.0185640763  0.0186459541  0.0189806914  0.0248204455 0.0425461408
#> 48   0.0306477419  0.0306515791  0.0353049349  0.0554094056 0.0572367130
#> 49   0.0198368169  0.0284593657  0.0362270388  0.0362886362 0.0427254147
#> 50   0.0669526850  0.0675065130  0.0723028194  0.0734498571 0.0740807788
#> 54  -0.0016975985  0.0039877798  0.0041579581  0.0082024968 0.0112345104
#> 55   0.0029721605  0.0135509190  0.0137364947  0.0139732229 0.0200387142
#> 56  -0.0047143213 -0.0046853450  0.0019546314  0.0026416738 0.0087200462
#> 57   0.0457534418  0.0506139349  0.0520388937  0.0538034362 0.0551459167
#> 58  -0.0197376766 -0.0171772960  0.0262607057  0.0277179365 0.0332937948
#> 59   0.0130654136  0.0216990493  0.0238952483  0.0240640578 0.0349131135
#> 60  -0.0075311177  0.0042328486  0.0054120166  0.0081404697 0.0289788972
#> 61  -0.0034603457  0.0040116245  0.0128487945  0.0227078962 0.0258502187
#> 62  -0.0149195857  0.0193317261  0.0208204180  0.0264085495 0.0447823908
#> 63  -0.0158317942  0.0258774139  0.0453920861  0.0488398878 0.0535356433
#> 64   0.0175163094  0.0469194705  0.0519096760  0.0540318701 0.0815408262
#> 65  -0.0118843763 -0.0036207394  0.0220994821  0.0282879112 0.0375695258
#> 66   0.0154238798  0.0175471122  0.0258850979  0.0265402365 0.0345039902
#> 67  -0.0033705550  0.0115197749  0.0147665496  0.0166574363 0.0488097254
#> 68  -0.0087801702 -0.0057359402 -0.0026543013  0.0123568406 0.0179223822
#> 69   0.0365059324  0.0370484673  0.0379315936  0.0434595598 0.0444263454
#> 70   0.1277246690  0.1391477584  0.1462148572  0.1535256143 0.1661467438
#> 71   0.0002050878  0.0050031313  0.0061347194  0.0098242517 0.0198951254
#> 72   0.0107603745  0.0183730798  0.0188919895  0.0231476415 0.0243508819
#> 73  -0.0210884457 -0.0131582383 -0.0128539729  0.0017546657 0.0097799562
#> 74   0.0036976533  0.0051688202  0.0066053980  0.0067029777 0.0159155696
#> 75   0.0294010125  0.0354938837  0.0359506974  0.0361817928 0.0488230095
#> 80  -0.0372857996 -0.0318083668 -0.0258268946  0.0150967487 0.0173807569
#> 81   0.0509156727  0.0663112189  0.0802423538  0.1023717047 0.1037914359
#> 82   0.0450868738  0.0452215191  0.0462365559  0.0660600827 0.0667727236
#> 83  -0.0769367830 -0.0497934514  0.0024274126  0.0176652914 0.0179631830
#> 84   0.0224557785  0.0226706842  0.0231481372  0.0271848826 0.0275397470
#> 85  -0.0497728941 -0.0026896449 -0.0023286117  0.0161207144 0.0224408147
#> 86  -0.0201132026 -0.0199166375 -0.0096437987  0.0177535967 0.0179717422
#> 87   0.0013283676  0.0056787930  0.0078242754  0.0250930190 0.0278098422
#> 88   0.0019380086  0.0092751146  0.0197480841  0.0328291695 0.0328294098
#> 89   0.0120143890  0.0129662377  0.0139595315  0.0221403002 0.0249630821
#> 90  -0.0067303866 -0.0067287940  0.0189689848  0.0462970697 0.0475678596
#> 91   0.0008972994  0.0043298644  0.0141088837  0.0177307189 0.0177977928
#> 92  -0.0306600172  0.0011185949  0.0045228588  0.0150915564 0.0216737050
#> 93  -0.0311097019 -0.0294977247 -0.0231322479  0.0238205834 0.0239212708
#> 94  -0.0292145606 -0.0120112019 -0.0118388676  0.0151469004 0.0169409819
#> 95   0.0117396763  0.0117541792  0.0136833598  0.0284228795 0.0289409387
#> 96  -0.0227550292 -0.0222697781 -0.0183841943  0.0096385522 0.0117891294
#> 97  -0.0658676118 -0.0654519349 -0.0564416947  0.0224603967 0.0243523185
#> 98  -0.0841530730 -0.0834010497 -0.0833114742  0.0384366527 0.0415111874
#> 99  -0.0011770271  0.0027641861  0.0030460870  0.0111347234 0.0113466435
#> 100  0.0485879409  0.0486328477  0.0487109327  0.0535726310 0.0578459627
#> 106  0.0259018033  0.0486679525  0.0700657917  0.0703645325 0.0703832640
#> 107  0.0860745990  0.0880228350  0.0924673499  0.0936352073 0.0965762745
#> 108 -0.0778941383 -0.0336622093  0.0392784105  0.0399517748 0.0400210056
#> 109  0.0176799055  0.0197277931  0.0227858526  0.0236333385 0.0236597704
#> 110 -0.0531125407  0.0165990056  0.0196270334  0.0209647727 0.0241842799
#> 111  0.0232965771  0.0255040840  0.0310906593  0.0372229283 0.0390595244
#> 112  0.0437603586  0.0454969571  0.0458351490  0.0481883271 0.0492326297
#> 113  0.0091749011  0.0130962656  0.0310699343  0.0322770388 0.0327307538
#> 114  0.0208848022  0.0209429014  0.0243786025  0.0246006299 0.0257931131
#> 115  0.0074411263  0.0084630389  0.0487951011  0.0535883509 0.0537075742
#> 116  0.0025888558  0.0118767398  0.0165385043  0.0185743484 0.0189623582
#> 117 -0.0227290686  0.0221471368  0.0230703767  0.0233897467 0.0271425882
#> 118  0.0004184159  0.0059568204  0.0080033709  0.0222728633 0.0225593446
#> 119 -0.0239824088  0.0066613758  0.0093585266  0.0129685241 0.0185794215
#> 120  0.0272442579  0.0280645947  0.0332153187  0.0352028278 0.0352097158
#> 121  0.0150202368  0.0184401469  0.0190193781  0.0222444255 0.0225942404
#> 122 -0.0157009220 -0.0096729502 -0.0083502281  0.0206902689 0.0371671123
#> 123  0.0372316357  0.0379881517  0.0423450316  0.0861355966 0.0867548153
#> 124 -0.0069869904  0.0030952464  0.0059726918  0.0060944747 0.0062451361
#> 125  0.0927878080  0.0944466783  0.0965352349  0.0981481915 0.0995300531
#> 132  0.0919068375  0.0983529004  0.1018525667  0.1021826341 0.1053306331
#> 133  0.0329938991  0.0364362887  0.0517631644  0.0518679289 0.0519842208
#> 134  0.0115519269  0.0153081702  0.0181508078  0.0195614259 0.0196181494
#> 135  0.0119296989  0.0238507325  0.0278894313  0.0282335079 0.0303788407
#> 136 -0.0071017172 -0.0030079573  0.0239385187  0.0266318922 0.0280948263
#> 137 -0.0080379866  0.0095290984  0.0220873783  0.0229804864 0.0235925779
#> 138  0.0221134219  0.0449186494  0.0457864887  0.0462228014 0.0467314778
#> 139  0.0147733238  0.0239526360  0.0251110358  0.0251219480 0.0260144581
#> 140  0.0790208662  0.0864868180  0.0909458405  0.0935016884 0.0935467625
#> 141  0.0085070647  0.0095068859  0.0397846223  0.0425051083 0.0429146718
#> 142  0.0212867704  0.0313509943  0.0492115946  0.0492352785 0.0523945059
#> 143  0.0055440875  0.0080976043  0.0330383215  0.0415936080 0.0419397693
#> 144  0.0372356407  0.0380952513  0.0437818733  0.0452274254 0.0500298211
#> 145  0.0849657882  0.0896424496  0.0902957017  0.0914107196 0.0914108074
#> 146  0.0277311323  0.0331261024  0.0544584430  0.0556609241 0.0558203753
#> 147  0.0166985562  0.0327035671  0.0591892200  0.0621670412 0.0669852999
#> 148  0.0038344430  0.0318840292  0.0473928350  0.0756815182 0.0760096666
#> 149  0.0438830976  0.0444924323  0.0487474338  0.0492705451 0.0494638588
#> 150  0.0626885014  0.0709354036  0.0783574676  0.0813804817 0.0824332355
#> 158  0.0010455224  0.0252596679  0.0664940424  0.0665376298 0.0684705368
#> 159  0.0303717118  0.0303932240  0.0304170285  0.0334024665 0.0350187675
#> 160  0.0054450094  0.0426101993  0.0427210162  0.0427259601 0.0529851796
#> 161  0.0183566251  0.0183681027  0.0326609202  0.0341536247 0.0342225782
#> 162  0.0037256874  0.0090025759  0.0138124053  0.0140887517 0.0192904015
#> 163  0.0853703702  0.0977696812  0.1061967999  0.1062723431 0.1068886586
#> 164  0.1200375756  0.1219645369  0.1220671727  0.1221980399 0.1288728714
#> 165  0.0172651686  0.0173726589  0.0295362445  0.0304572018 0.0332927724
#> 166  0.0192924321  0.0214524680  0.0374790008  0.0428911038 0.0439983792
#> 167 -0.0053726016  0.0276843908  0.0357217049  0.0358454069 0.0480233584
#> 168  0.0765678909  0.0776844168  0.0940726171  0.1030413202 0.1048126539
#> 169  0.0499704846  0.0627939599  0.0631308755  0.0635849034 0.0638218793
#> 170  0.0249318723  0.0249515970  0.0253489115  0.0256569315 0.0273596540
#> 171  0.0359172087  0.0361035554  0.0540666115  0.0547080077 0.0654630632
#> 172  0.0797559859  0.0802138388  0.0863858651  0.0872377510 0.0872468880
#> 173  0.1071865730  0.1078336612  0.1084216564  0.1181794317 0.1226573418
#> 174  0.0339897297  0.0366961869  0.0368652533  0.0387238319 0.0403515155
#> 175  0.0571122177  0.0571301412  0.0575111829  0.0614057315 0.0673318407
#> 184 -0.0228532165 -0.0014492265  0.0404685237  0.0436017731 0.0436104610
#> 185  0.0220878220  0.0275410345  0.0961261067  0.0962368139 0.1017464086
#> 186 -0.1269116362 -0.0998686785  0.0220875897  0.0245376330 0.0256248662
#> 187 -0.1022521058 -0.0431940242  0.0300931793  0.0307042775 0.0323410840
#> 188  0.0252528557  0.0961178296  0.1054686386  0.1057023030 0.1059845231
#> 189  0.0141620864  0.0440514606  0.0692782578  0.0693132838 0.0713494781
#> 190  0.0455629456  0.0672134373  0.0710921931  0.0726422239 0.0729614251
#> 191 -0.1099613357 -0.0960652392  0.0254827027  0.0300288777 0.0301241706
#> 192 -0.0349364036 -0.0320768337  0.0511209238  0.0511492962 0.0565165158
#> 193 -0.1097885594 -0.0922325808  0.0070957846  0.0142331359 0.0142897877
#> 194 -0.0080865504 -0.0061073865  0.0587621870  0.0596978840 0.0640161937
#> 195 -0.0160185957 -0.0027010801  0.0140283543  0.0145498740 0.0146288170
#> 196 -0.1148872257 -0.0870993535  0.0095352150  0.0100979408 0.0107445075
#> 197 -0.0739417070 -0.0384394885  0.0427361263  0.0439079871 0.0460561972
#> 198 -0.0224369369  0.0048629694  0.0476679857  0.0572802827 0.0578084853
#> 199 -0.0525671659 -0.0374696013  0.0289746382  0.0307754823 0.0307892689
#> 200  0.0015811278  0.0190971370  0.0550143234  0.0579172886 0.0592379671
#> 210 -0.0441218884  0.0254031884  0.0254498941  0.0323238042 0.0379196733
#> 211 -0.0081056779 -0.0080972851  0.0427739899  0.0694660133 0.0723760074
#> 212  0.0233904348  0.0364251178  0.0453004678  0.0551202803 0.0574042885
#> 213  0.0357088259  0.0495147803  0.0590834341  0.0637908376 0.0641587474
#> 214  0.0858336578  0.0895949709  0.0899500139  0.0927822674 0.0956128904
#> 215  0.0520956772  0.0525738307  0.0796864596  0.0932674489 0.0936236862
#> 216  0.0299052651  0.0338300204  0.0653516764  0.0655998588 0.0659279481
#> 217 -0.0098154635  0.0248637399  0.0327349981  0.0350921673 0.0402333000
#> 218  0.0076053265  0.0081891821  0.0202739648  0.0437661137 0.0438937258
#> 219 -0.0009932691  0.0198215267  0.0200937750  0.0328150472 0.0393770973
#> 220  0.0535957122  0.0537166039  0.0547040744  0.0617414502 0.0618117739
#> 221  0.0015340210  0.0015430173  0.0077935014  0.0141030821 0.0145639135
#> 222  0.0089201189  0.0095430121  0.0140304266  0.0217956363 0.0240680456
#> 223  0.0085236510  0.0093206370  0.0096159605  0.0312107849 0.0316019582
#> 224  0.0566621687  0.0609706325  0.0610770120  0.0618732090 0.0619485102
#> 225  0.0799462276  0.0800505273  0.0803161146  0.0803570775 0.0819649976
#> 236 -0.0626534787  0.0157540250  0.0417222512  0.0444835542 0.0586148899
#> 237 -0.0555879871  0.0344692595  0.0396377350  0.0399449461 0.0403080676
#> 238 -0.0696925823  0.0163461624  0.0255729585  0.0256102455 0.0306433421
#> 239 -0.0043058042  0.0510494142  0.0514921680  0.0517021588 0.0517852912
#> 240  0.0022199478  0.0655145202  0.0900197615  0.0913593948 0.0936314166
#> 241 -0.0509939498  0.0031910302  0.0286360028  0.0390053078 0.0479421805
#> 242  0.0499027844  0.0499075606  0.0573838805  0.0575902744 0.0578310231
#> 243 -0.0060666672  0.0421673933  0.0554255360  0.0639983069 0.0708927161
#> 244  0.0071977713  0.0367192694  0.0368671650  0.0379553723 0.0789905063
#> 245  0.0190567357  0.0493124320  0.0503260729  0.0506458935 0.0525071136
#> 246 -0.0225927134  0.0276316313  0.0344964190  0.0347464888 0.0365935753
#> 247 -0.0560160275  0.0045866931  0.0092375765  0.0099815750 0.0234892210
#> 248 -0.0317145872  0.0178739835  0.0180481623  0.0270814451 0.0283751558
#> 249 -0.0246360838  0.0183289426  0.0183361388  0.0220887440 0.0298239424
#> 250  0.0085343606  0.0540186217  0.0541387581  0.0595182614 0.0599442890
#> 262  0.0374446302  0.0538846037  0.0592106603  0.0603196765 0.0712189807
#> 263 -0.0137446333  0.0002378642  0.0494524773  0.0505158241 0.0507551255
#> 264  0.0505624333  0.0535596725  0.0744153085  0.0775027048 0.0843774782
#> 265 -0.0372416304 -0.0368820667  0.0498812630  0.0499637875 0.0530251398
#> 266  0.0437702701  0.0478207551  0.0480618501  0.0776150769 0.0787497836
#> 267 -0.0278951472  0.0094724169  0.0109062992  0.0140250552 0.0261051275
#> 268  0.0353860212  0.0360132016  0.0371148909  0.0385635888 0.0393565929
#> 269  0.0069823040  0.0282669252  0.0530437555  0.0539918810 0.0549785749
#> 270  0.0273174279  0.0273984409  0.0494129591  0.0497676072 0.0514533038
#> 271 -0.0014878236 -0.0014667850  0.0013905628  0.0022093553 0.0057572097
#> 272  0.0272438087  0.0278791435  0.0321423637  0.0324483092 0.0326169800
#> 273  0.0304053031  0.0311921034  0.0422658341  0.0450367993 0.0480845491
#> 274  0.0236476721  0.0293599940  0.0630965462  0.0780172973 0.0797968034
#> 275  0.0051238435  0.0051943865  0.0177426067  0.0297450018 0.0350600826
#> 288  0.0520831908  0.0530876842  0.0899049824  0.0900315879 0.0940012766
#> 289  0.1034474266  0.1046213194  0.1130578383  0.1140884265 0.1141682799
#> 290  0.0102332483  0.0152700845  0.0538815344  0.0540835764 0.0545700613
#> 291 -0.0122904137  0.0126307792  0.0150985534  0.0268200604 0.0304796373
#> 292 -0.0297991527  0.0684918128  0.0687258012  0.0698419750 0.0711501019
#> 293  0.0030720014  0.0140762956  0.0144877312  0.0174727617 0.0196495928
#> 294 -0.0111362740  0.0321029093  0.0366908663  0.0366911511 0.0475220767
#> 295  0.0405645927  0.0460781344  0.0560888699  0.0560986435 0.0569277003
#> 296  0.0142502914  0.0219359301  0.0219364183  0.0219477118 0.0222496949
#> 297  0.0442321623  0.0481930606  0.0482676293  0.0483195174 0.0554080333
#> 298  0.0409653473  0.0437957194  0.0463962928  0.0512717032 0.0515189468
#> 299 -0.0151124736  0.0217205434  0.0302764439  0.0369248253 0.0404968845
#> 300  0.0192131161  0.0240406849  0.0267393709  0.0323726592 0.0323733778
#> 314  0.1515661191  0.1568129381  0.1657341815  0.1663624411 0.1723246112
#> 315  0.0330950356  0.0401819316  0.0408647367  0.0413235094 0.0421914104
#> 316 -0.0439508492 -0.0151261002  0.0342950940  0.0407810670 0.0408324254
#> 317 -0.1448195923 -0.0107878510  0.0296617698  0.0301256171 0.0392893366
#> 318 -0.0055186885  0.0135813607  0.0494814772  0.0536835613 0.0537600732
#> 319 -0.0374836449  0.0099236071  0.0214704995  0.0215828067 0.0230548157
#> 320  0.0788157459  0.0867541190  0.0899439918  0.0900765113 0.0906212506
#> 321  0.0068971503  0.0225670256  0.0539053219  0.0539527499 0.0555681909
#> 322  0.0001980778  0.0069503053  0.0279078207  0.0281489195 0.0287470562
#> 323  0.0135669296  0.0187525237  0.0280381462  0.0332769837 0.0344037534
#> 324 -0.0044720511  0.0334096406  0.0482840821  0.0515691524 0.0518069878
#> 325 -0.0128247199 -0.0053650431  0.0029881897  0.0066438370 0.0087492859
#> 340  0.0899562891  0.0911224471  0.1012399591  0.1030049283 0.1038304386
#> 341  0.0341246109  0.0448143573  0.0653947692  0.0692617896 0.0733274717
#> 342  0.0286864406  0.0851922922  0.0956564496  0.0956577525 0.0962193783
#> 343  0.0511226212  0.0572947177  0.0707183857  0.0789767246 0.0824076587
#> 344  0.0188120053  0.0437203643  0.0445498274  0.0456415039 0.0571582243
#> 345  0.1263737807  0.1280049549  0.1281748064  0.1292396538 0.1305590611
#> 346  0.0708943447  0.0740291258  0.0819890700  0.0827491774 0.0834321932
#> 347  0.0100205683  0.0106401949  0.0152543156  0.0163039032 0.0223800105
#> 348  0.0059492418  0.0063269476  0.0070813620  0.0144284600 0.0148535997
#> 349  0.0495187748  0.0613097679  0.0619814018  0.0629286579 0.0663088530
#> 350  0.0904279948  0.0918876758  0.0926352522  0.0947065185 0.0947770647
#> 366 -0.0651166053 -0.0599798771  0.0251973062  0.0388222197 0.0397519362
#> 367 -0.0656215764 -0.0342179615  0.0014784885  0.0030186449 0.0053691895
#> 368 -0.0069259922 -0.0050768668  0.0555359136  0.0575158925 0.0582483964
#> 369  0.0074233455  0.0320263189  0.0646676950  0.0649957689 0.0734619018
#> 370  0.1018963376  0.1019306604  0.1114303035  0.1115643618 0.1116149188
#> 371 -0.0075160167 -0.0071048472  0.0355953994  0.0359302259 0.0359560559
#> 372  0.0028049856  0.0028829278  0.0492638854  0.0493315725 0.0538468344
#> 373  0.0620019241  0.0622163237  0.0881255560  0.0925926658 0.0926317283
#> 374  0.0131962970  0.0179174232  0.0403749079  0.0470164239 0.0475480490
#> 375  0.0621037098  0.0622322385  0.1094699388  0.1290636583 0.1299857776
#> 392 -0.0137162184  0.0150413218  0.0160120078  0.0207351657 0.0301020954
#> 393  0.0007435419  0.0013826002  0.0021566303  0.0419127351 0.0419195700
#> 394 -0.0612532640 -0.0445375642 -0.0056612607  0.0179561669 0.0236110568
#> 395  0.0232051460  0.0257986350  0.0477850436  0.0564719378 0.0568884456
#> 396 -0.0014269781  0.0007753127  0.0035181825  0.0159269314 0.0176674008
#> 397 -0.0103845939 -0.0047270536 -0.0012603406  0.0117510808 0.0133898234
#> 398 -0.0145113097 -0.0090208435  0.0018016930  0.0323154086 0.0334972107
#> 399 -0.0036811989 -0.0035451340  0.0369216513  0.0395556443 0.0396734658
#> 400  0.0231105298  0.0262112507  0.0391758531  0.0392031887 0.0423824074
#> 418 -0.0099257844  0.0153411107  0.0153894217  0.0227527188 0.0291104771
#> 419 -0.0158005501 -0.0076525609 -0.0007511146  0.0003823916 0.0196486967
#> 420  0.0273395218  0.0614582718  0.0740356080  0.0750515166 0.0785923115
#> 421 -0.0122070414  0.0217112786  0.0218948102  0.0226600776 0.0253732839
#> 422 -0.0298906025  0.0092320302  0.0096986107  0.0109931613 0.0234980996
#> 423 -0.0089720781  0.0283034969  0.0322928702  0.0416022596 0.0436441524
#> 424 -0.0315455749  0.0011944828  0.0118895832  0.0132599729 0.0227043520
#> 425 -0.0087400196  0.0184796109  0.0223777190  0.0245500246 0.0254071970
#> 444  0.0134398038  0.0284813433  0.0447007346  0.0516737345 0.0566152168
#> 445  0.0276997324  0.0284709685  0.0406786390  0.0428727658 0.0431321420
#> 446  0.0282003716  0.0293157193  0.0307516889  0.0427483407 0.0456852423
#> 447  0.0224113737  0.0251079118  0.0262917999  0.0295464171 0.0313229662
#> 448  0.0409331297  0.0437627160  0.0503186467  0.0507519118 0.0518311632
#> 449 -0.0272580440 -0.0263935995 -0.0144972637  0.0053777769 0.0054039454
#> 450  0.0001756438  0.0014343851  0.0100578705  0.0360107864 0.0388306402
#> 470  0.0445103217  0.0554740975  0.0571090252  0.0571161807 0.0611761782
#> 471  0.0189248359  0.0367360847  0.0442011474  0.0442303093 0.0554308699
#> 472 -0.0084579310  0.0156568307  0.0201006698  0.0201712641 0.0203695087
#> 473 -0.0049307896  0.0158639237  0.0159174660  0.0230898647 0.0307370963
#> 474  0.0076624429  0.0152990563  0.0153514531  0.0245945182 0.0312609606
#> 475 -0.0021845352  0.0154034607  0.0154247225  0.0260839471 0.0389926290
#> 496  0.0099463505  0.0100791279  0.0180436858  0.0180763858 0.0181833190
#> 497  0.0335473556  0.0336811087  0.0404175182  0.0404291993 0.0426503925
#> 498  0.0808657651  0.0811358959  0.0829363563  0.0870125172 0.0871384097
#> 499  0.0464940775  0.0493998914  0.0507631885  0.0545835986 0.0547739133
#> 500  0.0911715635  0.0911728541  0.0931281075  0.1002854734 0.1010668986
#> 522  0.0384527449  0.0395937863  0.0397333333  0.0399141809 0.0469233185
#> 523  0.0408043494  0.0421102990  0.0460645064  0.0538868223 0.0538874019
#> 524  0.0178014675  0.0204739057  0.0275167351  0.0323250929 0.0333647773
#> 525  0.0388706604  0.0390488196  0.0431271990  0.0507078602 0.0510362758
#> 548  0.0456590319  0.0457598784  0.0532620114  0.0682722358 0.0814233480
#> 549 -0.0188518765 -0.0124928599 -0.0079589449 -0.0019987706 0.0002014544
#> 550  0.0536034734  0.0539003815  0.0573274218  0.0687903578 0.0781087460
#> 574 -0.0055225446  0.0003673425  0.0004986317  0.0201359528 0.0208882182
#> 575  0.0669549784  0.0675625842  0.0675741960  0.1048028304 0.1051619930
#> 600  0.0078122451  0.0108283024  0.0108997140  0.0115553601 0.0135157493
#>         pcoa_9d    pcoa_10d  tree_average
#> 2   0.053053576 0.054054217  8.091011e-02
#> 3   0.036682285 0.037536087  8.903622e-02
#> 4   0.011039372 0.011046331 -3.764924e-03
#> 5   0.023428945 0.023430222  4.374243e-02
#> 6   0.058537914 0.058828275  3.442120e-02
#> 7   0.083043516 0.083237066  3.267223e-02
#> 8   0.054264722 0.054327057  3.015781e-02
#> 9   0.057392916 0.057574335 -6.046101e-02
#> 10  0.028096420 0.028410081 -3.396603e-02
#> 11  0.051829673 0.052022130 -5.365185e-02
#> 12  0.034511902 0.035661523 -1.245955e-01
#> 13  0.019369034 0.019369215 -2.636504e-02
#> 14  0.050409401 0.050531358  3.226313e-02
#> 15  0.073565942 0.073913549  8.214708e-02
#> 16  0.022793831 0.023244607  1.932518e-02
#> 17  0.050942834 0.050952672  1.528467e-01
#> 18  0.071199434 0.072290120  1.724429e-01
#> 19  0.049389715 0.059101410  0.000000e+00
#> 20  0.053260835 0.053301009  1.401818e-02
#> 21  0.039312429 0.047502556  1.099061e-01
#> 22  0.022194983 0.024927428  7.663176e-02
#> 23  0.033488412 0.033579516  3.072559e-02
#> 24  0.040257060 0.041504555 -7.177487e-03
#> 25  0.045675834 0.046034289  6.290105e-02
#> 28  0.025002008 0.025006667 -6.023881e-02
#> 29  0.025243178 0.026340931  6.975253e-02
#> 30  0.045490744 0.046415906  1.172599e-01
#> 31  0.069693271 0.069795079  1.781211e-02
#> 32  0.098251669 0.101199118  1.723115e-02
#> 33  0.017228790 0.017628827 -2.972329e-02
#> 34  0.020886656 0.022418645 -5.705887e-02
#> 35  0.053705784 0.053992008  5.181836e-02
#> 36  0.021386485 0.022952352 -5.024971e-02
#> 37  0.036131991 0.036194022 -1.213904e-01
#> 38  0.069045211 0.069718200  8.042053e-02
#> 39  0.090289582 0.090476885  7.389719e-02
#> 40  0.047612951 0.047693771 -3.473951e-02
#> 41  0.037644794 0.040174835  7.727277e-02
#> 42  0.034269668 0.034767615  3.033936e-02
#> 43  0.070334875 0.076121032 -6.714119e-03
#> 44  0.059229597 0.060296841  1.003573e-01
#> 45  0.018851227 0.019098323 -1.653185e-01
#> 46  0.053960811 0.063245766  0.000000e+00
#> 47  0.046807305 0.047284250  1.077250e-01
#> 48  0.058350721 0.059711872  3.757637e-02
#> 49  0.048723284 0.048723999  5.323015e-02
#> 50  0.080347975 0.080439088  5.308517e-02
#> 54  0.014554275 0.015176319 -1.887152e-01
#> 55  0.020052507 0.020575581 -1.412078e-01
#> 56  0.014523260 0.014620490 -1.990139e-01
#> 57  0.055264695 0.056479933 -9.059457e-02
#> 58  0.033358590 0.033703732 -7.076648e-02
#> 59  0.037165891 0.039343286  1.828202e-01
#> 60  0.030165393 0.030443701  4.169745e-02
#> 61  0.031407089 0.033199943  1.062960e-01
#> 62  0.045715755 0.045753704  1.184887e-01
#> 63  0.057495067 0.058316754 -7.519333e-03
#> 64  0.083062410 0.083462635  1.526240e-01
#> 65  0.041651332 0.041764828 -5.328442e-03
#> 66  0.040778765 0.043257568  1.812933e-01
#> 67  0.049486380 0.050265561  6.124315e-02
#> 68  0.021055795 0.023506381  7.051539e-03
#> 69  0.047858577 0.048498410  1.261547e-01
#> 70  0.170447933 0.171062617  0.000000e+00
#> 71  0.022863535 0.024262107 -3.124281e-02
#> 72  0.024411766 0.024594817 -1.083185e-01
#> 73  0.010467713 0.011389563 -1.784671e-01
#> 74  0.037737958 0.037741626  1.571202e-01
#> 75  0.048924799 0.049028969 -3.871589e-02
#> 80  0.032377088 0.032415359  6.701435e-02
#> 81  0.104745836 0.105336382  0.000000e+00
#> 82  0.070760368 0.070886219  6.235854e-02
#> 83  0.021840356 0.021933907 -4.710630e-02
#> 84  0.028375966 0.028448587 -9.892254e-02
#> 85  0.023810530 0.024011942 -2.713322e-02
#> 86  0.017972741 0.018049880 -9.211338e-02
#> 87  0.029208472 0.030259854 -1.632540e-01
#> 88  0.032873980 0.032876766 -1.313892e-01
#> 89  0.025872881 0.025994323 -2.212459e-01
#> 90  0.047572613 0.047948337 -8.720252e-02
#> 91  0.017810208 0.017974196 -5.112882e-03
#> 92  0.023397945 0.023419623 -8.197541e-02
#> 93  0.024091963 0.024560619 -2.178024e-02
#> 94  0.017254181 0.020787761  1.568230e-02
#> 95  0.029072736 0.029128197 -1.826837e-01
#> 96  0.011972434 0.018561220  4.075653e-02
#> 97  0.033806294 0.038945670  2.313601e-02
#> 98  0.044227411 0.044297652 -3.085098e-02
#> 99  0.014923742 0.015823356 -3.144488e-02
#> 100 0.066386095 0.066851004  1.285159e-02
#> 106 0.082473208 0.082804219 -4.591924e-02
#> 107 0.096660658 0.096887877  9.370427e-02
#> 108 0.040045201 0.040086390  4.010489e-04
#> 109 0.024946234 0.025101607 -5.004156e-02
#> 110 0.025289715 0.025417090  2.037413e-02
#> 111 0.043832294 0.044007854  1.239947e-02
#> 112 0.049783699 0.050798439 -5.874120e-02
#> 113 0.035315392 0.035316997 -8.388183e-02
#> 114 0.026316666 0.026390526 -1.737385e-01
#> 115 0.058317927 0.058594414 -3.969517e-02
#> 116 0.024054947 0.024358844  4.239447e-02
#> 117 0.027452538 0.027455382 -3.446805e-02
#> 118 0.026061274 0.026776276  2.572711e-02
#> 119 0.022031103 0.025537459  6.318965e-02
#> 120 0.037097768 0.037120847 -1.351764e-01
#> 121 0.026942518 0.033541386  8.826388e-02
#> 122 0.037284763 0.046758046  0.000000e+00
#> 123 0.088546004 0.088768862 -8.363858e-03
#> 124 0.024843667 0.025612303  1.606247e-02
#> 125 0.099851981 0.100178482  6.035895e-02
#> 132 0.112674324 0.113583419  7.630223e-02
#> 133 0.058342033 0.058410612 -8.164746e-02
#> 134 0.021834924 0.022407506 -1.936743e-01
#> 135 0.033602065 0.033621983 -8.479709e-02
#> 136 0.028549774 0.029127134 -1.868652e-01
#> 137 0.026507451 0.026727797 -2.580058e-01
#> 138 0.047373932 0.047550561 -1.659303e-01
#> 139 0.028263910 0.028275681 -2.557870e-01
#> 140 0.094185897 0.094186454 -5.962524e-02
#> 141 0.043588673 0.044462755 -5.934268e-02
#> 142 0.056303130 0.056438807 -4.378929e-02
#> 143 0.043206122 0.044712018 -4.469552e-02
#> 144 0.051688065 0.053046343 -3.398769e-04
#> 145 0.092423930 0.092484967 -6.106299e-02
#> 146 0.057214966 0.060321909  2.862630e-02
#> 147 0.077948166 0.079106321 -4.737334e-02
#> 148 0.083877914 0.085091974  3.399321e-02
#> 149 0.050536865 0.050638341 -6.598604e-02
#> 150 0.096084161 0.096084947  2.819944e-02
#> 158 0.068485714 0.068930444  1.251919e-02
#> 159 0.035636144 0.035637573 -1.801094e-01
#> 160 0.053465741 0.054040178 -7.123220e-02
#> 161 0.036917276 0.036918485 -1.733003e-01
#> 162 0.019481635 0.021129803 -2.444409e-01
#> 163 0.109342952 0.109491490  1.226630e-01
#> 164 0.129177780 0.129643205  3.280634e-02
#> 165 0.035755550 0.036443895 -1.577901e-01
#> 166 0.047190335 0.047193369 -4.577779e-02
#> 167 0.048127782 0.048344951  1.680706e-02
#> 168 0.107917547 0.108085607 -4.527648e-02
#> 169 0.065890560 0.070788641  1.322501e-02
#> 170 0.028535944 0.028770395 -2.259190e-01
#> 171 0.070666807 0.086048065  2.804533e-02
#> 172 0.087264006 0.091185411  9.225017e-02
#> 173 0.123034793 0.123051833  6.250561e-02
#> 174 0.054796155 0.056434274 -3.914276e-02
#> 175 0.067841042 0.068608295 -6.996539e-02
#> 184 0.044648874 0.045019919 -2.908850e-02
#> 185 0.102755657 0.102782396  1.504958e-01
#> 186 0.029412661 0.029772144 -6.268338e-02
#> 187 0.032707087 0.033394177 -9.342001e-02
#> 188 0.108974021 0.109033828  1.461499e-01
#> 189 0.071787417 0.071802924 -1.399339e-02
#> 190 0.076333767 0.076421650 -2.281046e-02
#> 191 0.033806524 0.034309301 -3.819121e-02
#> 192 0.056717600 0.056735330  5.885737e-02
#> 193 0.016796199 0.017762344 -7.000870e-02
#> 194 0.066875920 0.069615964  8.343786e-02
#> 195 0.016064020 0.016064111 -2.121515e-01
#> 196 0.013666787 0.018096292 -2.496971e-02
#> 197 0.046056201 0.047540193  9.047958e-03
#> 198 0.058280266 0.058509320 -6.110064e-02
#> 199 0.051239463 0.051762055  7.147413e-02
#> 200 0.059573651 0.059648971 -1.228646e-01
#> 210 0.037969013 0.038855116 -1.918284e-02
#> 211 0.074392539 0.074392568  0.000000e+00
#> 212 0.057605751 0.060748323  5.333763e-02
#> 213 0.064661720 0.064780958  8.105531e-02
#> 214 0.095687839 0.096349082  2.068552e-01
#> 215 0.094827580 0.096012909  1.641554e-01
#> 216 0.067554483 0.067574231 -8.706571e-05
#> 217 0.040505547 0.040707603  5.824687e-02
#> 218 0.044287790 0.044461135  5.448761e-03
#> 219 0.039770135 0.047236299  7.471059e-03
#> 220 0.062085739 0.062426346  9.592239e-02
#> 221 0.014946656 0.021911591 -6.787306e-02
#> 222 0.025060840 0.028188051 -4.139466e-02
#> 223 0.031662169 0.031668664 -1.115433e-01
#> 224 0.078192959 0.081077162 -6.722097e-03
#> 225 0.084961630 0.085860105  5.814886e-02
#> 236 0.060580847 0.061424189 -1.237368e-02
#> 237 0.040348400 0.040904344 -8.331738e-02
#> 238 0.031316601 0.031436391 -1.713816e-02
#> 239 0.051790497 0.051790536  4.149002e-02
#> 240 0.095126978 0.095165575  1.234252e-01
#> 241 0.049887151 0.051096430  8.684429e-03
#> 242 0.057966332 0.058043992  1.115685e-01
#> 243 0.071713236 0.073522904  1.005897e-01
#> 244 0.080042695 0.084390322  3.396603e-02
#> 245 0.052877263 0.052893653 -4.520038e-02
#> 246 0.037400729 0.041949012  4.100418e-02
#> 247 0.024303390 0.025603908  2.902104e-02
#> 248 0.028384032 0.028791967 -4.112756e-02
#> 249 0.046123773 0.046539370  2.049397e-03
#> 250 0.062371241 0.062405394  3.953275e-02
#> 262 0.074476899 0.078445113  1.298618e-01
#> 263 0.050800469 0.050917008  4.746044e-02
#> 264 0.085759741 0.086290730  1.060886e-01
#> 265 0.053026640 0.053979955  4.722714e-02
#> 266 0.078760771 0.078778684  6.722097e-03
#> 267 0.028306027 0.028520407  6.505603e-02
#> 268 0.039526115 0.039701856  1.225792e-02
#> 269 0.055394142 0.062621427  1.428022e-02
#> 270 0.051607595 0.051906484  1.939822e-02
#> 271 0.005889322 0.013154003 -6.106389e-02
#> 272 0.037027984 0.040684805  3.203539e-02
#> 273 0.049447565 0.049455278 -4.933116e-02
#> 274 0.086777573 0.089655707  8.706571e-05
#> 275 0.041548972 0.042339137 -4.261774e-02
#> 288 0.095366137 0.096760229  1.833905e-01
#> 289 0.114182488 0.114694526  1.425237e-01
#> 290 0.056008504 0.056277755 -6.684282e-02
#> 291 0.032587638 0.035865327  2.339327e-03
#> 292 0.071180818 0.072325638  1.605820e-01
#> 293 0.020438712 0.023429818 -1.422161e-01
#> 294 0.048282037 0.048642785 -5.666348e-02
#> 295 0.057625245 0.058294143  3.159089e-02
#> 296 0.023094669 0.024210603 -1.322046e-01
#> 297 0.055779054 0.055873688 -3.910528e-02
#> 298 0.051527105 0.053025277 -1.204718e-01
#> 299 0.057710094 0.057804123  7.903763e-02
#> 300 0.033610448 0.033838932 -1.728493e-01
#> 314 0.173540898 0.173731070  0.000000e+00
#> 315 0.042215659 0.042429049 -1.070933e-01
#> 316 0.040854581 0.041054518  7.195261e-02
#> 317 0.040958045 0.040969849  3.261669e-03
#> 318 0.053780007 0.054264151  4.013512e-02
#> 319 0.023118790 0.025738791  2.691501e-02
#> 320 0.090649140 0.090686519 -1.489043e-01
#> 321 0.055583289 0.060427192  8.517411e-02
#> 322 0.030985455 0.032510531 -7.523492e-02
#> 323 0.035058154 0.035104872 -1.453835e-01
#> 324 0.057447891 0.058393593  1.816179e-01
#> 325 0.012301392 0.012494028 -2.071475e-01
#> 340 0.104877985 0.104899870  1.870655e-02
#> 341 0.074715057 0.075479920  1.487626e-01
#> 342 0.096293163 0.096364510 -3.261669e-03
#> 343 0.083084515 0.084332369  1.169451e-01
#> 344 0.057655066 0.059299079  8.554319e-02
#> 345 0.131056297 0.131075846 -6.942826e-03
#> 346 0.084091000 0.087160996  7.865077e-02
#> 347 0.022751158 0.023462307 -1.650916e-01
#> 348 0.014853816 0.015107383 -2.352402e-01
#> 349 0.077035562 0.077293940  1.569128e-01
#> 350 0.096170698 0.096190924 -4.700413e-02
#> 366 0.039753530 0.040938283  4.241528e-02
#> 367 0.006856455 0.006998886 -1.299018e-01
#> 368 0.058370126 0.060159975  2.776813e-02
#> 369 0.073760260 0.075637766  1.500791e-01
#> 370 0.111758643 0.111857623  1.299960e-01
#> 371 0.036048192 0.038710922 -4.555369e-02
#> 372 0.057903507 0.058744864 -3.104826e-02
#> 373 0.094188325 0.094936069  9.914252e-03
#> 374 0.051710746 0.051821261  3.578025e-02
#> 375 0.141619721 0.141619726  0.000000e+00
#> 392 0.032603199 0.032991335  1.562149e-01
#> 393 0.042067144 0.042201789  1.044283e-01
#> 394 0.024038507 0.035097530  5.493401e-02
#> 395 0.057006416 0.057453461  7.621363e-02
#> 396 0.017784769 0.028774476  6.645858e-02
#> 397 0.017950508 0.022454961  5.104138e-02
#> 398 0.034830874 0.034865593 -1.910722e-02
#> 399 0.050263068 0.054831045  0.000000e+00
#> 400 0.050199683 0.051417508  6.155309e-02
#> 418 0.030226528 0.030867849  2.237719e-02
#> 419 0.020711212 0.023358364  8.491465e-02
#> 420 0.079537335 0.079550773 -5.050066e-02
#> 421 0.026651269 0.031078276  5.933537e-02
#> 422 0.023683336 0.025158004 -1.774034e-02
#> 423 0.043698347 0.043804520 -4.748490e-02
#> 424 0.039104131 0.039919322  1.992135e-01
#> 425 0.026296902 0.026439369 -1.148044e-01
#> 444 0.056641249 0.068199604  1.676477e-01
#> 445 0.043134078 0.043822826 -7.378568e-02
#> 446 0.045687162 0.069779717  6.714119e-03
#> 447 0.034454169 0.040267121  5.861644e-02
#> 448 0.052583049 0.052816225  1.271027e-02
#> 449 0.012267685 0.015673875  1.294474e-02
#> 450 0.044705381 0.046650552  4.488573e-02
#> 470 0.061179547 0.063024702  3.925683e-02
#> 471 0.055463515 0.055964130  1.254614e-01
#> 472 0.023516622 0.023611127  9.607898e-02
#> 473 0.031353216 0.035699064  2.593038e-02
#> 474 0.043744125 0.045165781  4.610257e-02
#> 475 0.044894795 0.046545755  8.234827e-02
#> 496 0.018187371 0.020917409 -1.379250e-01
#> 497 0.044282285 0.045281382 -1.022870e-01
#> 498 0.087594690 0.087774423 -3.708211e-02
#> 499 0.059646103 0.060012920  5.350994e-02
#> 500 0.105223312 0.105314139  7.842669e-02
#> 522 0.050671382 0.051591735  1.211532e-01
#> 523 0.054763573 0.061875827  5.100461e-02
#> 524 0.040454493 0.043034278  5.798373e-02
#> 525 0.056545635 0.059529412  2.408917e-02
#> 548 0.082827439 0.090116791  8.363858e-03
#> 549 0.018170381 0.018543741  2.470938e-02
#> 550 0.078639673 0.079582140  3.264222e-02
#> 574 0.030308156 0.031549853 -4.543922e-02
#> 575 0.107644274 0.108508316  7.360473e-02
#> 600 0.033227530 0.033320848 -4.811224e-02
#> 
#> $details_deviation$abs_dev_distsp
#>              sp.x          sp.y      pcoa_1d      pcoa_2d      pcoa_3d
#> 2         apricot         apple 0.1634793627 0.0613131726 0.0567510971
#> 3          banana         apple 0.0131320500 0.0047392768 0.0037427035
#> 4         currant         apple 0.1010576998 0.1006330170 0.0813228267
#> 5      blackberry         apple 0.1078357918 0.1075952587 0.0657691608
#> 6       blueberry         apple 0.0423032775 0.0375635362 0.0370317817
#> 7          cherry         apple 0.1616933564 0.0044655132 0.0185309251
#> 8           grape         apple 0.3014574226 0.1180817603 0.0492048338
#> 9      grapefruit         apple 0.0782134735 0.0620468549 0.0406330304
#> 10      kiwifruit         apple 0.1182215595 0.1151277801 0.0521533878
#> 11          lemon         apple 0.1163257141 0.1107436173 0.0223964456
#> 12           lime         apple 0.1717551260 0.1715743017 0.0544854936
#> 13         litchi         apple 0.2208832159 0.1118805638 0.1066939412
#> 14          mango         apple 0.0487058101 0.0364587186 0.0208028301
#> 15          melon         apple 0.2616967886 0.0665607673 0.0638025397
#> 16         orange         apple 0.0875671626 0.0866583269 0.0779119903
#> 17  passion_fruit         apple 0.1096418149 0.0838538643 0.0646337510
#> 18          peach         apple 0.1348560405 0.0634209972 0.0167038995
#> 19           pear         apple 0.0406116659 0.0398174530 0.0264415399
#> 20      pineapple         apple 0.1919552128 0.0362306401 0.0345022788
#> 21           plum         apple 0.1480900385 0.0408337829 0.0309493853
#> 22      raspberry         apple 0.0967829650 0.0967745107 0.0553955783
#> 23     strawberry         apple 0.0862982218 0.0736475401 0.0584288511
#> 24      tangerine         apple 0.0825254764 0.0733131341 0.0300916371
#> 25    water_melon         apple 0.1703111415 0.0530135968 0.0415572094
#> 28         banana       apricot 0.0857023219 0.0242932667 0.0241109329
#> 29        currant       apricot 0.1042450037 0.0545054018 0.0146009012
#> 30     blackberry       apricot 0.1110230958 0.0436697614 0.0248648307
#> 31      blueberry       apricot 0.1356171216 0.0461783026 0.0453007506
#> 32         cherry       apricot 0.0730266810 0.0090297993 0.0019372342
#> 33          grape       apricot 0.2956743988 0.2906455598 0.1377553518
#> 34     grapefruit       apricot 0.1183953004 0.0243571666 0.0006936907
#> 35      kiwifruit       apricot 0.1584033864 0.0853719907 0.0419028383
#> 36          lemon       apricot 0.1565075411 0.0769666856 0.0021222401
#> 37           lime       apricot 0.1286036196 0.0935940794 0.0201976926
#> 38         litchi       apricot 0.0373928817 0.0221055863 0.0075622973
#> 39          mango       apricot 0.0696330090 0.0740771360 0.0795749880
#> 40          melon       apricot 0.3439820414 0.0187968599 0.0114241535
#> 41         orange       apricot 0.0732035350 0.0003625967 0.0138623513
#> 42  passion_fruit       apricot 0.1554444099 0.1497362756 0.1188010163
#> 43          peach       apricot 0.0480830291 0.0111865203 0.0054971105
#> 44           pear       apricot 0.1167117412 0.0229469924 0.0105328150
#> 45      pineapple       apricot 0.2945870897 0.0399771690 0.0399703870
#> 46           plum       apricot 0.0153893242 0.0149825960 0.0060163710
#> 47      raspberry       apricot 0.1423945114 0.0604310211 0.0049237280
#> 48     strawberry       apricot 0.1561521924 0.0310780438 0.0034079774
#> 49      tangerine       apricot 0.0657018088 0.0364410525 0.0117241279
#> 50    water_melon       apricot 0.2568317788 0.0651181273 0.0668086545
#> 54        currant        banana 0.0232806589 0.0163820775 0.0030818734
#> 55     blackberry        banana 0.0300587510 0.0232900603 0.0002809691
#> 56      blueberry        banana 0.0130110851 0.0124973002 0.0124278354
#> 57         cherry        banana 0.0660341977 0.0420972235 0.0455975807
#> 58          grape        banana 0.1674144255 0.0804373324 0.0372737598
#> 59     grapefruit        banana 0.1024010940 0.1024008480 0.0646444426
#> 60      kiwifruit        banana 0.0946385403 0.0816225538 0.0610000317
#> 61          lemon        banana 0.1340038606 0.1306800188 0.0487136663
#> 62           lime        banana 0.2208557785 0.1930523006 0.0171950080
#> 63         litchi        banana 0.2912331831 0.0275701602 0.0163747357
#> 64          mango        banana 0.2022404267 0.0462949627 0.0269137086
#> 65          melon        banana 0.0864450002 0.0626896105 0.0580126769
#> 66         orange        banana 0.0147881642 0.0069464042 0.0015288165
#> 67  passion_fruit        banana 0.1545003362 0.0712841793 0.0422261591
#> 68          peach        banana 0.0813214239 0.0476578445 0.0410890920
#> 69           pear        banana 0.0033340150 0.0042738603 0.0091095536
#> 70      pineapple        banana 0.0557616968 0.1063538199 0.1065677751
#> 71           plum        banana 0.0703129977 0.0088275768 0.0077939136
#> 72      raspberry        banana 0.0190059241 0.0132293892 0.0091482188
#> 73     strawberry        banana 0.0327636052 0.0327607052 0.0215307530
#> 74      tangerine        banana 0.0705055635 0.0383462584 0.0023877290
#> 75    water_melon        banana 0.0092007673 0.0084254471 0.0107939978
#> 80     blackberry       currant 0.0502274025 0.0500022593 0.0389149021
#> 81      blueberry       currant 0.1287289277 0.0786732777 0.0090084589
#> 82         cherry       currant 0.1616628782 0.0164758571 0.0438968236
#> 83          grape       currant 0.1499194025 0.0975652203 0.0894623777
#> 84     grapefruit       currant 0.0131622822 0.0035252359 0.0034422918
#> 85      kiwifruit       currant 0.0902582803 0.0902582026 0.0900754121
#> 86          lemon       currant 0.0512745228 0.0475976058 0.0401120522
#> 87           lime       currant 0.0233706014 0.0233370624 0.0013152483
#> 88         litchi       currant 0.0583045521 0.0027321360 0.0013323560
#> 89          mango       currant 0.0346119947 0.0294009509 0.0052273453
#> 90          melon       currant 0.1388558607 0.0480657892 0.0410986853
#> 91         orange       currant 0.0084924947 0.0073959172 0.0021501549
#> 92  passion_fruit       currant 0.0768611308 0.0683091410 0.0682377642
#> 93          peach       currant 0.1328683259 0.1236499732 0.0536370397
#> 94           pear       currant 0.0542900783 0.0534915442 0.0446038482
#> 95      pineapple       currant 0.1210543157 0.0215406693 0.0057361827
#> 96           plum       currant 0.1196343279 0.0745078580 0.0263173850
#> 97      raspberry       currant 0.0722692403 0.0711190425 0.0658724616
#> 98     strawberry       currant 0.1826471459 0.1050775552 0.1041838161
#> 99      tangerine       currant 0.0032801459 0.0018235003 0.0016200373
#> 100   water_melon       currant 0.1494334270 0.0212950301 0.0341925798
#> 106     blueberry    blackberry 0.1425926938 0.1198579992 0.0319868888
#> 107        cherry    blackberry 0.1846025864 0.0029153831 0.0834226480
#> 108         grape    blackberry 0.1566974945 0.0904380138 0.0881675515
#> 109    grapefruit    blackberry 0.0185667479 0.0087531920 0.0062425888
#> 110     kiwifruit    blackberry 0.0970363723 0.0969911261 0.0925690637
#> 111         lemon    blackberry 0.0010471203 0.0024633279 0.0049577139
#> 112          lime    blackberry 0.0268568011 0.0268576408 0.0433965702
#> 113        litchi    blackberry 0.0650826441 0.0016304175 0.0038078583
#> 114         mango    blackberry 0.0413900867 0.0350098427 0.0157181544
#> 115         melon    blackberry 0.1456339528 0.0417186480 0.0235214622
#> 116        orange    blackberry 0.0152705868 0.0143936233 0.0009361898
#> 117 passion_fruit    blackberry 0.0836392229 0.0729129995 0.0705637247
#> 118         peach    blackberry 0.1396464180 0.1263748437 0.0178358171
#> 119          pear    blackberry 0.0610681703 0.0604846474 0.0367793261
#> 120     pineapple    blackberry 0.1278324077 0.0217966063 0.0050104300
#> 121          plum    blackberry 0.1264124200 0.0656793526 0.0133150833
#> 122     raspberry    blackberry 0.0220418378 0.0205056631 0.0193202051
#> 123    strawberry    blackberry 0.0635398526 0.0013808365 0.0106468099
#> 124     tangerine    blackberry 0.0100582379 0.0078413278 0.0070103500
#> 125   water_melon    blackberry 0.1562115190 0.0018367241 0.0823861359
#> 132        cherry     blueberry 0.1271508801 0.0667314715 0.0743952001
#> 133         grape     blueberry 0.1638922529 0.0749047500 0.0277515215
#> 134    grapefruit     blueberry 0.0873457554 0.0866216698 0.0734037705
#> 135     kiwifruit     blueberry 0.1273538414 0.1203937378 0.1022813840
#> 136         lemon     blueberry 0.1254579960 0.1252654078 0.0807582478
#> 137          lime     blueberry 0.0975540745 0.0939381975 0.0239826308
#> 138        litchi     blueberry 0.0722774025 0.0120965610 0.0171258743
#> 139         mango     blueberry 0.0485848452 0.0320055949 0.0255240022
#> 140         melon     blueberry 0.0907102741 0.0450716630 0.0407286918
#> 141        orange     blueberry 0.0421539899 0.0401296328 0.0355770832
#> 142 passion_fruit     blueberry 0.0181067085 0.0066244862 0.0188959426
#> 143         peach     blueberry 0.1352153022 0.1030979994 0.0920402991
#> 144          pear     blueberry 0.0497439477 0.0464304175 0.0428117339
#> 145     pineapple     blueberry 0.0211347275 0.0845822197 0.0849161857
#> 146          plum     blueberry 0.1111962560 0.0279512489 0.0252093743
#> 147     raspberry     blueberry 0.1222102892 0.1080164379 0.0322289736
#> 148    strawberry     blueberry 0.0972346594 0.0886234114 0.0189774446
#> 149     tangerine     blueberry 0.0172529963 0.0049974611 0.0155522376
#> 150   water_melon     blueberry 0.1135172775 0.0470372950 0.0394305372
#> 158         grape        cherry 0.2034080395 0.1819563757 0.0216705651
#> 159    grapefruit        cherry 0.1914219815 0.0392690809 0.0029767167
#> 160     kiwifruit        cherry 0.2314300675 0.0957104915 0.0427964444
#> 161         lemon        cherry 0.2295342221 0.0914128522 0.0070527005
#> 162          lime        cherry 0.2016303006 0.1213014649 0.0031759239
#> 163        litchi        cherry 0.0548734776 0.0548761653 0.0840397080
#> 164         mango        cherry 0.0785660350 0.1050981332 0.1054867477
#> 165         melon        cherry 0.4170087224 0.0319605478 0.0159188158
#> 166        orange        cherry 0.1462302160 0.0126209943 0.0117021164
#> 167 passion_fruit        cherry 0.1189528314 0.0823886704 0.0371119584
#> 168         peach        cherry 0.0726248617 0.0268893826 0.0270746275
#> 169          pear        cherry 0.1538201738 0.0035825380 0.0267583808
#> 170     pineapple        cherry 0.3051637208 0.0123794817 0.0144561003
#> 171          plum        cherry 0.0486058154 0.0112353493 0.0155654633
#> 172     raspberry        cherry 0.2078931939 0.0025180762 0.0785077694
#> 173    strawberry        cherry 0.1812468345 0.0625514695 0.1055026953
#> 174     tangerine        cherry 0.1080508342 0.0287779733 0.0311242144
#> 175   water_melon        cherry 0.4299062189 0.0405236865 0.0407749698
#> 184    grapefruit         grape 0.1776506941 0.0389963195 0.0292206691
#> 185     kiwifruit         grape 0.1469517094 0.0082199132 0.0098665282
#> 186         lemon         grape 0.2561669751 0.1299319274 0.1298461535
#> 187          lime         grape 0.1878590133 0.1286605997 0.1161421041
#> 188        litchi         grape 0.0428479090 0.0323309533 0.0108948220
#> 189         mango         grape 0.0894420095 0.0788520328 0.0141223634
#> 190         melon         grape 0.3721746529 0.0165680287 0.0403708538
#> 191        orange         grape 0.2758932721 0.1502588494 0.1135193000
#> 192 passion_fruit         grape 0.1981108370 0.1817451705 0.1677844751
#> 193         peach         grape 0.3988692584 0.3021089636 0.1117148453
#> 194          pear         grape 0.2208569673 0.0555072487 0.0081946104
#> 195     pineapple         grape 0.4126045546 0.1000071709 0.0617766611
#> 196          plum         grape 0.3045274988 0.2969118822 0.1189354807
#> 197     raspberry         grape 0.1698870919 0.0921123859 0.0886479777
#> 198    strawberry         grape 0.1836447730 0.0639164533 0.0540593268
#> 199     tangerine         grape 0.1346835872 0.0740454082 0.0651176166
#> 200   water_melon         grape 0.3615970866 0.0682179928 0.0015414264
#> 210     kiwifruit    grapefruit 0.0711030251 0.0441897349 0.0441808613
#> 211         lemon    grapefruit 0.0720068070 0.0574729403 0.0112342765
#> 212          lime    grapefruit 0.1768749858 0.1277915443 0.0334733482
#> 213        litchi    grapefruit 0.2965158962 0.0276659703 0.0262255156
#> 214         mango    grapefruit 0.0571667732 0.0095175572 0.0806111285
#> 215         melon    grapefruit 0.0033645187 0.0448281429 0.0520934779
#> 216        orange    grapefruit 0.0451917655 0.0237047061 0.0115933960
#> 217 passion_fruit    grapefruit 0.3872946972 0.2206505970 0.2206499567
#> 218         peach    grapefruit 0.1187971198 0.0602496655 0.0060147092
#> 219          pear    grapefruit 0.0376018076 0.0217725252 0.0101165707
#> 220     pineapple    grapefruit 0.2931040437 0.1118771239 0.0795045717
#> 221          plum    grapefruit 0.1428161660 0.0473697432 0.0143895225
#> 222     raspberry    grapefruit 0.0317563453 0.0232212604 0.0213712782
#> 223    strawberry    grapefruit 0.0455140263 0.0455111169 0.0454878345
#> 224     tangerine    grapefruit 0.0833711473 0.0010010376 0.0006271003
#> 225   water_melon    grapefruit 0.0079897040 0.0348185587 0.0753277173
#> 236         lemon     kiwifruit 0.1092152657 0.0957822071 0.0657676637
#> 237          lime     kiwifruit 0.1646446776 0.1645673882 0.1026883259
#> 238        litchi     kiwifruit 0.2458240496 0.1371476849 0.1358194774
#> 239         mango     kiwifruit 0.0736466437 0.0639642335 0.0092587302
#> 240         melon     kiwifruit 0.1862509044 0.0049093272 0.0021948199
#> 241        orange     kiwifruit 0.1323756286 0.1261557728 0.1084963855
#> 242 passion_fruit     kiwifruit 0.1850876990 0.1622803795 0.1622669962
#> 243         peach     kiwifruit 0.1725414695 0.1362741422 0.0092101013
#> 244          pear     kiwifruit 0.0776098936 0.0618923103 0.0179709910
#> 245     pineapple     kiwifruit 0.2853414899 0.1035086810 0.0828889510
#> 246          plum     kiwifruit 0.1828242521 0.1068380682 0.0496642265
#> 247     raspberry     kiwifruit 0.1102259697 0.1098479850 0.1064294583
#> 248    strawberry     kiwifruit 0.1239836508 0.1083939369 0.1083815700
#> 249     tangerine     kiwifruit 0.1074663101 0.1013098765 0.0985915281
#> 250   water_melon     kiwifruit 0.1595117220 0.0448971383 0.0023158473
#> 262          lime         lemon 0.0554294119 0.0440033211 0.0174707022
#> 263        litchi         lemon 0.2851893699 0.0839240827 0.0631119153
#> 264         mango         lemon 0.1130119641 0.0771103121 0.0491707097
#> 265         melon         lemon 0.1584850956 0.0807337454 0.0385406786
#> 266        orange         lemon 0.0833040061 0.0759364321 0.0221946807
#> 267 passion_fruit         lemon 0.3355641305 0.2406276248 0.2182694085
#> 268         peach         lemon 0.1569093604 0.1128743083 0.0353315784
#> 269          pear         lemon 0.0757140482 0.0708962607 0.0042115984
#> 270     pineapple         lemon 0.3247068103 0.1364359942 0.0658292112
#> 271          plum         lemon 0.1809284067 0.0999191547 0.0109778229
#> 272     raspberry         lemon 0.0032477068 0.0006499275 0.0030460928
#> 273    strawberry         lemon 0.0282233365 0.0261011387 0.0153774821
#> 274     tangerine         lemon 0.1214833879 0.0446382132 0.0229531023
#> 275   water_melon         lemon 0.1376982942 0.0877297884 0.0037470198
#> 288        litchi          lime 0.2483039141 0.0217930487 0.0468610036
#> 289         mango          lime 0.1756214577 0.1508903151 0.0913165327
#> 290         melon          lime 0.1735104671 0.0896162194 0.0202748993
#> 291        orange          lime 0.1387334180 0.1369989612 0.0171901107
#> 292 passion_fruit          lime 0.2665158530 0.2004120061 0.1061164965
#> 293         peach          lime 0.2123387723 0.2014167640 0.0149157567
#> 294          pear          lime 0.1311434601 0.1305420922 0.0272277478
#> 295     pineapple          lime 0.4115587282 0.0987718372 0.0171440540
#> 296          plum          lime 0.1530244852 0.1181001709 0.0115767364
#> 297     raspberry          lime 0.0246562147 0.0247211904 0.0442083322
#> 298    strawberry          lime 0.0003194151 0.0078614002 0.0372874012
#> 299     tangerine          lime 0.0935794665 0.0867387481 0.0164402005
#> 300   water_melon          lime 0.1688852818 0.1079670551 0.0100307363
#> 314         mango        litchi 0.0465941004 0.0425935670 0.1077304908
#> 315         melon        litchi 0.2805598024 0.0091424317 0.0084356276
#> 316        orange        litchi 0.2164787240 0.0624470945 0.0617668873
#> 317 passion_fruit        litchi 0.2970999869 0.1907680478 0.1877670736
#> 318         peach        litchi 0.1405877413 0.0873849630 0.0496809370
#> 319          pear        litchi 0.1949235647 0.0725905586 0.0711589856
#> 320     pineapple        litchi 0.5135128627 0.0637928780 0.0723343382
#> 321          plum        litchi 0.0462459817 0.0293441938 0.0087295223
#> 322     raspberry        litchi 0.0782722415 0.0080886709 0.0035309912
#> 323    strawberry        litchi 0.0920299226 0.0129486694 0.0135573477
#> 324     tangerine        litchi 0.1383577395 0.0220291761 0.0165950275
#> 325   water_melon        litchi 0.2699822362 0.0477992554 0.0324275118
#> 340         melon         mango 0.0412106794 0.0657895239 0.0865077974
#> 341        orange         mango 0.0261195000 0.0056275802 0.0305854672
#> 342 passion_fruit         mango 0.1900740962 0.1899212168 0.0957518626
#> 343         peach         mango 0.0497714827 0.0499996401 0.0501674011
#> 344          pear         mango 0.0227461588 0.0077762284 0.0188106524
#> 345     pineapple         mango 0.2580021235 0.0558671409 0.0621341082
#> 346          plum         mango 0.0607799090 0.0648937934 0.0680529849
#> 347     raspberry         mango 0.0545796842 0.0466460735 0.0027453467
#> 348    strawberry         mango 0.0683373652 0.0437523447 0.0131284688
#> 349     tangerine         mango 0.0495136670 0.0449576967 0.0408009408
#> 350   water_melon         mango 0.0037103212 0.0893564581 0.0904279946
#> 366        orange         melon 0.2075154219 0.0900011804 0.0898954652
#> 367 passion_fruit         melon 0.3713386034 0.2037118201 0.1990172341
#> 368         peach         melon 0.3298710510 0.0377726845 0.0085874499
#> 369          pear         melon 0.1664443186 0.0007544235 0.0007145660
#> 370     pineapple         melon 0.0549257276 0.0347226196 0.0260552069
#> 371          plum         melon 0.3684029070 0.0346912098 0.0232303063
#> 372     raspberry         melon 0.1588235502 0.0533585067 0.0365142488
#> 373    strawberry         melon 0.0614701201 0.0243536751 0.0185412522
#> 374     tangerine         melon 0.1826061033 0.0270876441 0.0162623537
#> 375   water_melon         melon 0.0046251853 0.0046922425 0.0507853782
#> 392 passion_fruit        orange 0.2001868179 0.1478148749 0.1390260345
#> 393         peach        orange 0.1089575021 0.0715811087 0.0133868438
#> 394          pear        orange 0.0792787290 0.0789659721 0.0780444327
#> 395     pineapple        orange 0.2236729320 0.0366574875 0.0277580578
#> 396          plum        orange 0.0976244006 0.0222391098 0.0014300479
#> 397     raspberry        orange 0.0284601842 0.0280716885 0.0139766380
#> 398    strawberry        orange 0.0422178652 0.0356935385 0.0323076676
#> 399     tangerine        orange 0.0781209845 0.0400577149 0.0076529300
#> 400   water_melon        orange 0.0777459364 0.0030403701 0.0193293975
#> 418         peach passion_fruit 0.2263160372 0.2253011779 0.1495314151
#> 419          pear passion_fruit 0.2048942848 0.1720904684 0.1621585973
#> 420     pineapple passion_fruit 0.3700492732 0.0052956504 0.0253514552
#> 421          plum passion_fruit 0.1400550857 0.1347416264 0.0937330859
#> 422     raspberry passion_fruit 0.0887480122 0.0754017176 0.0736992375
#> 423    strawberry passion_fruit 0.0621016528 0.0253873805 0.0253619862
#> 424     tangerine passion_fruit 0.1887325001 0.1701867806 0.1678855457
#> 425   water_melon passion_fruit 0.2456095220 0.1192071653 0.0881428395
#> 444          pear         peach 0.1123308433 0.0538243938 0.0113213139
#> 445     pineapple         peach 0.2659637675 0.0478806800 0.0438545332
#> 446          plum         peach 0.0482614705 0.0085348753 0.0008147577
#> 447     raspberry         peach 0.1285935911 0.1104562209 0.0017350935
#> 448    strawberry         peach 0.1181088479 0.0662003872 0.0082466993
#> 449     tangerine         peach 0.1688966685 0.1645280209 0.0479623498
#> 450   water_melon         peach 0.2021217676 0.0019419155 0.0010125348
#> 470     pineapple          pear 0.1940369646 0.0328245368 0.0269012183
#> 471          plum          pear 0.1052143584 0.0075938212 0.0125678879
#> 472     raspberry          pear 0.0500153435 0.0498528276 0.0272612427
#> 473    strawberry          pear 0.0637730246 0.0541107915 0.0478960196
#> 474     tangerine          pear 0.0565658252 0.0419158097 0.0130401894
#> 475   water_melon          pear 0.1235435200 0.0233696848 0.0021865076
#> 496          plum     pineapple 0.2808003296 0.0229906606 0.0226081209
#> 497     raspberry     pineapple 0.1167795809 0.0115141420 0.0140756941
#> 498    strawberry     pineapple 0.0048162734 0.0601871544 0.0746508002
#> 499     tangerine     pineapple 0.2779209674 0.0150749859 0.0139192041
#> 500   water_melon     pineapple 0.0041366871 0.0198892254 0.0217218905
#> 522     raspberry          plum 0.1153595931 0.0413685884 0.0347919512
#> 523    strawberry          plum 0.1291172742 0.0122414327 0.0226431675
#> 524     tangerine          plum 0.0745549089 0.0451737424 0.0164526336
#> 525   water_melon          plum 0.2722211030 0.0350499062 0.0354583298
#> 548    strawberry     raspberry 0.0249756298 0.0181488565 0.0267160315
#> 549     tangerine     raspberry 0.0232478353 0.0197593321 0.0193459093
#> 550   water_melon     raspberry 0.2057647528 0.0442244406 0.0387933564
#> 574     tangerine    strawberry 0.0370055164 0.0152975875 0.0144500881
#> 575   water_melon    strawberry 0.1084113227 0.0494145856 0.0172473566
#> 600   water_melon     tangerine 0.1558669209 0.0471430131 0.0004666108
#>          pcoa_4d      pcoa_5d      pcoa_6d      pcoa_7d      pcoa_8d
#> 2   0.0528852420 0.0333988631 0.0339991097 0.0421872222 0.0530128968
#> 3   0.0072900396 0.0197260797 0.0202440801 0.0285611434 0.0330724861
#> 4   0.0772269489 0.0237868122 0.0219995315 0.0108357707 0.0108467623
#> 5   0.0635186776 0.0125484403 0.0125868642 0.0187303827 0.0197521546
#> 6   0.0199440721 0.0351629611 0.0524944227 0.0559045095 0.0570411048
#> 7   0.0281458699 0.0731641406 0.0783982461 0.0800797013 0.0806243911
#> 8   0.0457945661 0.0416295350 0.0484971172 0.0508539651 0.0513522627
#> 9   0.0156670849 0.0373972686 0.0427345586 0.0562298325 0.0569609754
#> 10  0.0017553655 0.0017506389 0.0069481718 0.0103475530 0.0267954936
#> 11  0.0083407947 0.0454018507 0.0513779010 0.0514114852 0.0516278013
#> 12  0.0506307980 0.0293530414 0.0295385242 0.0298783547 0.0336665740
#> 13  0.0909142259 0.0034354962 0.0186229480 0.0193365308 0.0193481370
#> 14  0.0179677430 0.0373474919 0.0429085187 0.0452516629 0.0498183062
#> 15  0.0381959608 0.0241986355 0.0715810607 0.0715960261 0.0734237918
#> 16  0.0779045654 0.0122554013 0.0047473389 0.0223458167 0.0225959398
#> 17  0.0342945375 0.0342953875 0.0352693917 0.0383269678 0.0495093071
#> 18  0.0035489276 0.0653982569 0.0678630003 0.0708725951 0.0711992901
#> 19  0.0034553576 0.0333230276 0.0432566380 0.0445190412 0.0493499262
#> 20  0.0108899739 0.0448487287 0.0522024909 0.0524028442 0.0532592880
#> 21  0.0309322278 0.0347371715 0.0351037985 0.0359590227 0.0393104590
#> 22  0.0505415553 0.0177486307 0.0177823051 0.0179594778 0.0187928866
#> 23  0.0355896728 0.0239622226 0.0261640935 0.0303233066 0.0326773201
#> 24  0.0275826903 0.0134767349 0.0190861287 0.0307247672 0.0312625144
#> 25  0.0371159138 0.0176539919 0.0204510351 0.0347866971 0.0395560055
#> 28  0.0077329903 0.0020990516 0.0034955847 0.0036797634 0.0217944008
#> 29  0.0015546801 0.0011710552 0.0036198210 0.0160033955 0.0251620784
#> 30  0.0341744537 0.0344592035 0.0354113439 0.0354291618 0.0402678424
#> 31  0.0319843570 0.0437132823 0.0654164211 0.0657377916 0.0686685199
#> 32  0.0469869234 0.0487700617 0.0664300893 0.0691707809 0.0931576272
#> 33  0.1377166522 0.0915612653 0.0081726135 0.0091080295 0.0138727363
#> 34  0.0071120459 0.0082376801 0.0154202998 0.0161559835 0.0203086324
#> 35  0.0299822755 0.0429829372 0.0507459073 0.0525913916 0.0526155478
#> 36  0.0005050269 0.0015667955 0.0039279806 0.0107799101 0.0213369999
#> 37  0.0289220112 0.0318920444 0.0319225550 0.0345753813 0.0350440886
#> 38  0.0209021675 0.0285272105 0.0608455748 0.0626952680 0.0690451296
#> 39  0.0796916045 0.0801801264 0.0888254092 0.0892566850 0.0894282379
#> 40  0.0049695465 0.0047243118 0.0408596167 0.0455206549 0.0475829584
#> 41  0.0168115778 0.0249382267 0.0269362847 0.0290669650 0.0376148564
#> 42  0.0124184745 0.0321312713 0.0321967185 0.0326895430 0.0327990199
#> 43  0.0115392987 0.0209398961 0.0217281882 0.0574567505 0.0702756952
#> 44  0.0102186893 0.0222489718 0.0305542383 0.0356185873 0.0591169970
#> 45  0.0060044575 0.0063003090 0.0147129587 0.0169191983 0.0188221997
#> 46  0.0174340360 0.0274630074 0.0276582284 0.0445927639 0.0538695247
#> 47  0.0185640763 0.0186459541 0.0189806914 0.0248204455 0.0425461408
#> 48  0.0306477419 0.0306515791 0.0353049349 0.0554094056 0.0572367130
#> 49  0.0198368169 0.0284593657 0.0362270388 0.0362886362 0.0427254147
#> 50  0.0669526850 0.0675065130 0.0723028194 0.0734498571 0.0740807788
#> 54  0.0016975985 0.0039877798 0.0041579581 0.0082024968 0.0112345104
#> 55  0.0029721605 0.0135509190 0.0137364947 0.0139732229 0.0200387142
#> 56  0.0047143213 0.0046853450 0.0019546314 0.0026416738 0.0087200462
#> 57  0.0457534418 0.0506139349 0.0520388937 0.0538034362 0.0551459167
#> 58  0.0197376766 0.0171772960 0.0262607057 0.0277179365 0.0332937948
#> 59  0.0130654136 0.0216990493 0.0238952483 0.0240640578 0.0349131135
#> 60  0.0075311177 0.0042328486 0.0054120166 0.0081404697 0.0289788972
#> 61  0.0034603457 0.0040116245 0.0128487945 0.0227078962 0.0258502187
#> 62  0.0149195857 0.0193317261 0.0208204180 0.0264085495 0.0447823908
#> 63  0.0158317942 0.0258774139 0.0453920861 0.0488398878 0.0535356433
#> 64  0.0175163094 0.0469194705 0.0519096760 0.0540318701 0.0815408262
#> 65  0.0118843763 0.0036207394 0.0220994821 0.0282879112 0.0375695258
#> 66  0.0154238798 0.0175471122 0.0258850979 0.0265402365 0.0345039902
#> 67  0.0033705550 0.0115197749 0.0147665496 0.0166574363 0.0488097254
#> 68  0.0087801702 0.0057359402 0.0026543013 0.0123568406 0.0179223822
#> 69  0.0365059324 0.0370484673 0.0379315936 0.0434595598 0.0444263454
#> 70  0.1277246690 0.1391477584 0.1462148572 0.1535256143 0.1661467438
#> 71  0.0002050878 0.0050031313 0.0061347194 0.0098242517 0.0198951254
#> 72  0.0107603745 0.0183730798 0.0188919895 0.0231476415 0.0243508819
#> 73  0.0210884457 0.0131582383 0.0128539729 0.0017546657 0.0097799562
#> 74  0.0036976533 0.0051688202 0.0066053980 0.0067029777 0.0159155696
#> 75  0.0294010125 0.0354938837 0.0359506974 0.0361817928 0.0488230095
#> 80  0.0372857996 0.0318083668 0.0258268946 0.0150967487 0.0173807569
#> 81  0.0509156727 0.0663112189 0.0802423538 0.1023717047 0.1037914359
#> 82  0.0450868738 0.0452215191 0.0462365559 0.0660600827 0.0667727236
#> 83  0.0769367830 0.0497934514 0.0024274126 0.0176652914 0.0179631830
#> 84  0.0224557785 0.0226706842 0.0231481372 0.0271848826 0.0275397470
#> 85  0.0497728941 0.0026896449 0.0023286117 0.0161207144 0.0224408147
#> 86  0.0201132026 0.0199166375 0.0096437987 0.0177535967 0.0179717422
#> 87  0.0013283676 0.0056787930 0.0078242754 0.0250930190 0.0278098422
#> 88  0.0019380086 0.0092751146 0.0197480841 0.0328291695 0.0328294098
#> 89  0.0120143890 0.0129662377 0.0139595315 0.0221403002 0.0249630821
#> 90  0.0067303866 0.0067287940 0.0189689848 0.0462970697 0.0475678596
#> 91  0.0008972994 0.0043298644 0.0141088837 0.0177307189 0.0177977928
#> 92  0.0306600172 0.0011185949 0.0045228588 0.0150915564 0.0216737050
#> 93  0.0311097019 0.0294977247 0.0231322479 0.0238205834 0.0239212708
#> 94  0.0292145606 0.0120112019 0.0118388676 0.0151469004 0.0169409819
#> 95  0.0117396763 0.0117541792 0.0136833598 0.0284228795 0.0289409387
#> 96  0.0227550292 0.0222697781 0.0183841943 0.0096385522 0.0117891294
#> 97  0.0658676118 0.0654519349 0.0564416947 0.0224603967 0.0243523185
#> 98  0.0841530730 0.0834010497 0.0833114742 0.0384366527 0.0415111874
#> 99  0.0011770271 0.0027641861 0.0030460870 0.0111347234 0.0113466435
#> 100 0.0485879409 0.0486328477 0.0487109327 0.0535726310 0.0578459627
#> 106 0.0259018033 0.0486679525 0.0700657917 0.0703645325 0.0703832640
#> 107 0.0860745990 0.0880228350 0.0924673499 0.0936352073 0.0965762745
#> 108 0.0778941383 0.0336622093 0.0392784105 0.0399517748 0.0400210056
#> 109 0.0176799055 0.0197277931 0.0227858526 0.0236333385 0.0236597704
#> 110 0.0531125407 0.0165990056 0.0196270334 0.0209647727 0.0241842799
#> 111 0.0232965771 0.0255040840 0.0310906593 0.0372229283 0.0390595244
#> 112 0.0437603586 0.0454969571 0.0458351490 0.0481883271 0.0492326297
#> 113 0.0091749011 0.0130962656 0.0310699343 0.0322770388 0.0327307538
#> 114 0.0208848022 0.0209429014 0.0243786025 0.0246006299 0.0257931131
#> 115 0.0074411263 0.0084630389 0.0487951011 0.0535883509 0.0537075742
#> 116 0.0025888558 0.0118767398 0.0165385043 0.0185743484 0.0189623582
#> 117 0.0227290686 0.0221471368 0.0230703767 0.0233897467 0.0271425882
#> 118 0.0004184159 0.0059568204 0.0080033709 0.0222728633 0.0225593446
#> 119 0.0239824088 0.0066613758 0.0093585266 0.0129685241 0.0185794215
#> 120 0.0272442579 0.0280645947 0.0332153187 0.0352028278 0.0352097158
#> 121 0.0150202368 0.0184401469 0.0190193781 0.0222444255 0.0225942404
#> 122 0.0157009220 0.0096729502 0.0083502281 0.0206902689 0.0371671123
#> 123 0.0372316357 0.0379881517 0.0423450316 0.0861355966 0.0867548153
#> 124 0.0069869904 0.0030952464 0.0059726918 0.0060944747 0.0062451361
#> 125 0.0927878080 0.0944466783 0.0965352349 0.0981481915 0.0995300531
#> 132 0.0919068375 0.0983529004 0.1018525667 0.1021826341 0.1053306331
#> 133 0.0329938991 0.0364362887 0.0517631644 0.0518679289 0.0519842208
#> 134 0.0115519269 0.0153081702 0.0181508078 0.0195614259 0.0196181494
#> 135 0.0119296989 0.0238507325 0.0278894313 0.0282335079 0.0303788407
#> 136 0.0071017172 0.0030079573 0.0239385187 0.0266318922 0.0280948263
#> 137 0.0080379866 0.0095290984 0.0220873783 0.0229804864 0.0235925779
#> 138 0.0221134219 0.0449186494 0.0457864887 0.0462228014 0.0467314778
#> 139 0.0147733238 0.0239526360 0.0251110358 0.0251219480 0.0260144581
#> 140 0.0790208662 0.0864868180 0.0909458405 0.0935016884 0.0935467625
#> 141 0.0085070647 0.0095068859 0.0397846223 0.0425051083 0.0429146718
#> 142 0.0212867704 0.0313509943 0.0492115946 0.0492352785 0.0523945059
#> 143 0.0055440875 0.0080976043 0.0330383215 0.0415936080 0.0419397693
#> 144 0.0372356407 0.0380952513 0.0437818733 0.0452274254 0.0500298211
#> 145 0.0849657882 0.0896424496 0.0902957017 0.0914107196 0.0914108074
#> 146 0.0277311323 0.0331261024 0.0544584430 0.0556609241 0.0558203753
#> 147 0.0166985562 0.0327035671 0.0591892200 0.0621670412 0.0669852999
#> 148 0.0038344430 0.0318840292 0.0473928350 0.0756815182 0.0760096666
#> 149 0.0438830976 0.0444924323 0.0487474338 0.0492705451 0.0494638588
#> 150 0.0626885014 0.0709354036 0.0783574676 0.0813804817 0.0824332355
#> 158 0.0010455224 0.0252596679 0.0664940424 0.0665376298 0.0684705368
#> 159 0.0303717118 0.0303932240 0.0304170285 0.0334024665 0.0350187675
#> 160 0.0054450094 0.0426101993 0.0427210162 0.0427259601 0.0529851796
#> 161 0.0183566251 0.0183681027 0.0326609202 0.0341536247 0.0342225782
#> 162 0.0037256874 0.0090025759 0.0138124053 0.0140887517 0.0192904015
#> 163 0.0853703702 0.0977696812 0.1061967999 0.1062723431 0.1068886586
#> 164 0.1200375756 0.1219645369 0.1220671727 0.1221980399 0.1288728714
#> 165 0.0172651686 0.0173726589 0.0295362445 0.0304572018 0.0332927724
#> 166 0.0192924321 0.0214524680 0.0374790008 0.0428911038 0.0439983792
#> 167 0.0053726016 0.0276843908 0.0357217049 0.0358454069 0.0480233584
#> 168 0.0765678909 0.0776844168 0.0940726171 0.1030413202 0.1048126539
#> 169 0.0499704846 0.0627939599 0.0631308755 0.0635849034 0.0638218793
#> 170 0.0249318723 0.0249515970 0.0253489115 0.0256569315 0.0273596540
#> 171 0.0359172087 0.0361035554 0.0540666115 0.0547080077 0.0654630632
#> 172 0.0797559859 0.0802138388 0.0863858651 0.0872377510 0.0872468880
#> 173 0.1071865730 0.1078336612 0.1084216564 0.1181794317 0.1226573418
#> 174 0.0339897297 0.0366961869 0.0368652533 0.0387238319 0.0403515155
#> 175 0.0571122177 0.0571301412 0.0575111829 0.0614057315 0.0673318407
#> 184 0.0228532165 0.0014492265 0.0404685237 0.0436017731 0.0436104610
#> 185 0.0220878220 0.0275410345 0.0961261067 0.0962368139 0.1017464086
#> 186 0.1269116362 0.0998686785 0.0220875897 0.0245376330 0.0256248662
#> 187 0.1022521058 0.0431940242 0.0300931793 0.0307042775 0.0323410840
#> 188 0.0252528557 0.0961178296 0.1054686386 0.1057023030 0.1059845231
#> 189 0.0141620864 0.0440514606 0.0692782578 0.0693132838 0.0713494781
#> 190 0.0455629456 0.0672134373 0.0710921931 0.0726422239 0.0729614251
#> 191 0.1099613357 0.0960652392 0.0254827027 0.0300288777 0.0301241706
#> 192 0.0349364036 0.0320768337 0.0511209238 0.0511492962 0.0565165158
#> 193 0.1097885594 0.0922325808 0.0070957846 0.0142331359 0.0142897877
#> 194 0.0080865504 0.0061073865 0.0587621870 0.0596978840 0.0640161937
#> 195 0.0160185957 0.0027010801 0.0140283543 0.0145498740 0.0146288170
#> 196 0.1148872257 0.0870993535 0.0095352150 0.0100979408 0.0107445075
#> 197 0.0739417070 0.0384394885 0.0427361263 0.0439079871 0.0460561972
#> 198 0.0224369369 0.0048629694 0.0476679857 0.0572802827 0.0578084853
#> 199 0.0525671659 0.0374696013 0.0289746382 0.0307754823 0.0307892689
#> 200 0.0015811278 0.0190971370 0.0550143234 0.0579172886 0.0592379671
#> 210 0.0441218884 0.0254031884 0.0254498941 0.0323238042 0.0379196733
#> 211 0.0081056779 0.0080972851 0.0427739899 0.0694660133 0.0723760074
#> 212 0.0233904348 0.0364251178 0.0453004678 0.0551202803 0.0574042885
#> 213 0.0357088259 0.0495147803 0.0590834341 0.0637908376 0.0641587474
#> 214 0.0858336578 0.0895949709 0.0899500139 0.0927822674 0.0956128904
#> 215 0.0520956772 0.0525738307 0.0796864596 0.0932674489 0.0936236862
#> 216 0.0299052651 0.0338300204 0.0653516764 0.0655998588 0.0659279481
#> 217 0.0098154635 0.0248637399 0.0327349981 0.0350921673 0.0402333000
#> 218 0.0076053265 0.0081891821 0.0202739648 0.0437661137 0.0438937258
#> 219 0.0009932691 0.0198215267 0.0200937750 0.0328150472 0.0393770973
#> 220 0.0535957122 0.0537166039 0.0547040744 0.0617414502 0.0618117739
#> 221 0.0015340210 0.0015430173 0.0077935014 0.0141030821 0.0145639135
#> 222 0.0089201189 0.0095430121 0.0140304266 0.0217956363 0.0240680456
#> 223 0.0085236510 0.0093206370 0.0096159605 0.0312107849 0.0316019582
#> 224 0.0566621687 0.0609706325 0.0610770120 0.0618732090 0.0619485102
#> 225 0.0799462276 0.0800505273 0.0803161146 0.0803570775 0.0819649976
#> 236 0.0626534787 0.0157540250 0.0417222512 0.0444835542 0.0586148899
#> 237 0.0555879871 0.0344692595 0.0396377350 0.0399449461 0.0403080676
#> 238 0.0696925823 0.0163461624 0.0255729585 0.0256102455 0.0306433421
#> 239 0.0043058042 0.0510494142 0.0514921680 0.0517021588 0.0517852912
#> 240 0.0022199478 0.0655145202 0.0900197615 0.0913593948 0.0936314166
#> 241 0.0509939498 0.0031910302 0.0286360028 0.0390053078 0.0479421805
#> 242 0.0499027844 0.0499075606 0.0573838805 0.0575902744 0.0578310231
#> 243 0.0060666672 0.0421673933 0.0554255360 0.0639983069 0.0708927161
#> 244 0.0071977713 0.0367192694 0.0368671650 0.0379553723 0.0789905063
#> 245 0.0190567357 0.0493124320 0.0503260729 0.0506458935 0.0525071136
#> 246 0.0225927134 0.0276316313 0.0344964190 0.0347464888 0.0365935753
#> 247 0.0560160275 0.0045866931 0.0092375765 0.0099815750 0.0234892210
#> 248 0.0317145872 0.0178739835 0.0180481623 0.0270814451 0.0283751558
#> 249 0.0246360838 0.0183289426 0.0183361388 0.0220887440 0.0298239424
#> 250 0.0085343606 0.0540186217 0.0541387581 0.0595182614 0.0599442890
#> 262 0.0374446302 0.0538846037 0.0592106603 0.0603196765 0.0712189807
#> 263 0.0137446333 0.0002378642 0.0494524773 0.0505158241 0.0507551255
#> 264 0.0505624333 0.0535596725 0.0744153085 0.0775027048 0.0843774782
#> 265 0.0372416304 0.0368820667 0.0498812630 0.0499637875 0.0530251398
#> 266 0.0437702701 0.0478207551 0.0480618501 0.0776150769 0.0787497836
#> 267 0.0278951472 0.0094724169 0.0109062992 0.0140250552 0.0261051275
#> 268 0.0353860212 0.0360132016 0.0371148909 0.0385635888 0.0393565929
#> 269 0.0069823040 0.0282669252 0.0530437555 0.0539918810 0.0549785749
#> 270 0.0273174279 0.0273984409 0.0494129591 0.0497676072 0.0514533038
#> 271 0.0014878236 0.0014667850 0.0013905628 0.0022093553 0.0057572097
#> 272 0.0272438087 0.0278791435 0.0321423637 0.0324483092 0.0326169800
#> 273 0.0304053031 0.0311921034 0.0422658341 0.0450367993 0.0480845491
#> 274 0.0236476721 0.0293599940 0.0630965462 0.0780172973 0.0797968034
#> 275 0.0051238435 0.0051943865 0.0177426067 0.0297450018 0.0350600826
#> 288 0.0520831908 0.0530876842 0.0899049824 0.0900315879 0.0940012766
#> 289 0.1034474266 0.1046213194 0.1130578383 0.1140884265 0.1141682799
#> 290 0.0102332483 0.0152700845 0.0538815344 0.0540835764 0.0545700613
#> 291 0.0122904137 0.0126307792 0.0150985534 0.0268200604 0.0304796373
#> 292 0.0297991527 0.0684918128 0.0687258012 0.0698419750 0.0711501019
#> 293 0.0030720014 0.0140762956 0.0144877312 0.0174727617 0.0196495928
#> 294 0.0111362740 0.0321029093 0.0366908663 0.0366911511 0.0475220767
#> 295 0.0405645927 0.0460781344 0.0560888699 0.0560986435 0.0569277003
#> 296 0.0142502914 0.0219359301 0.0219364183 0.0219477118 0.0222496949
#> 297 0.0442321623 0.0481930606 0.0482676293 0.0483195174 0.0554080333
#> 298 0.0409653473 0.0437957194 0.0463962928 0.0512717032 0.0515189468
#> 299 0.0151124736 0.0217205434 0.0302764439 0.0369248253 0.0404968845
#> 300 0.0192131161 0.0240406849 0.0267393709 0.0323726592 0.0323733778
#> 314 0.1515661191 0.1568129381 0.1657341815 0.1663624411 0.1723246112
#> 315 0.0330950356 0.0401819316 0.0408647367 0.0413235094 0.0421914104
#> 316 0.0439508492 0.0151261002 0.0342950940 0.0407810670 0.0408324254
#> 317 0.1448195923 0.0107878510 0.0296617698 0.0301256171 0.0392893366
#> 318 0.0055186885 0.0135813607 0.0494814772 0.0536835613 0.0537600732
#> 319 0.0374836449 0.0099236071 0.0214704995 0.0215828067 0.0230548157
#> 320 0.0788157459 0.0867541190 0.0899439918 0.0900765113 0.0906212506
#> 321 0.0068971503 0.0225670256 0.0539053219 0.0539527499 0.0555681909
#> 322 0.0001980778 0.0069503053 0.0279078207 0.0281489195 0.0287470562
#> 323 0.0135669296 0.0187525237 0.0280381462 0.0332769837 0.0344037534
#> 324 0.0044720511 0.0334096406 0.0482840821 0.0515691524 0.0518069878
#> 325 0.0128247199 0.0053650431 0.0029881897 0.0066438370 0.0087492859
#> 340 0.0899562891 0.0911224471 0.1012399591 0.1030049283 0.1038304386
#> 341 0.0341246109 0.0448143573 0.0653947692 0.0692617896 0.0733274717
#> 342 0.0286864406 0.0851922922 0.0956564496 0.0956577525 0.0962193783
#> 343 0.0511226212 0.0572947177 0.0707183857 0.0789767246 0.0824076587
#> 344 0.0188120053 0.0437203643 0.0445498274 0.0456415039 0.0571582243
#> 345 0.1263737807 0.1280049549 0.1281748064 0.1292396538 0.1305590611
#> 346 0.0708943447 0.0740291258 0.0819890700 0.0827491774 0.0834321932
#> 347 0.0100205683 0.0106401949 0.0152543156 0.0163039032 0.0223800105
#> 348 0.0059492418 0.0063269476 0.0070813620 0.0144284600 0.0148535997
#> 349 0.0495187748 0.0613097679 0.0619814018 0.0629286579 0.0663088530
#> 350 0.0904279948 0.0918876758 0.0926352522 0.0947065185 0.0947770647
#> 366 0.0651166053 0.0599798771 0.0251973062 0.0388222197 0.0397519362
#> 367 0.0656215764 0.0342179615 0.0014784885 0.0030186449 0.0053691895
#> 368 0.0069259922 0.0050768668 0.0555359136 0.0575158925 0.0582483964
#> 369 0.0074233455 0.0320263189 0.0646676950 0.0649957689 0.0734619018
#> 370 0.1018963376 0.1019306604 0.1114303035 0.1115643618 0.1116149188
#> 371 0.0075160167 0.0071048472 0.0355953994 0.0359302259 0.0359560559
#> 372 0.0028049856 0.0028829278 0.0492638854 0.0493315725 0.0538468344
#> 373 0.0620019241 0.0622163237 0.0881255560 0.0925926658 0.0926317283
#> 374 0.0131962970 0.0179174232 0.0403749079 0.0470164239 0.0475480490
#> 375 0.0621037098 0.0622322385 0.1094699388 0.1290636583 0.1299857776
#> 392 0.0137162184 0.0150413218 0.0160120078 0.0207351657 0.0301020954
#> 393 0.0007435419 0.0013826002 0.0021566303 0.0419127351 0.0419195700
#> 394 0.0612532640 0.0445375642 0.0056612607 0.0179561669 0.0236110568
#> 395 0.0232051460 0.0257986350 0.0477850436 0.0564719378 0.0568884456
#> 396 0.0014269781 0.0007753127 0.0035181825 0.0159269314 0.0176674008
#> 397 0.0103845939 0.0047270536 0.0012603406 0.0117510808 0.0133898234
#> 398 0.0145113097 0.0090208435 0.0018016930 0.0323154086 0.0334972107
#> 399 0.0036811989 0.0035451340 0.0369216513 0.0395556443 0.0396734658
#> 400 0.0231105298 0.0262112507 0.0391758531 0.0392031887 0.0423824074
#> 418 0.0099257844 0.0153411107 0.0153894217 0.0227527188 0.0291104771
#> 419 0.0158005501 0.0076525609 0.0007511146 0.0003823916 0.0196486967
#> 420 0.0273395218 0.0614582718 0.0740356080 0.0750515166 0.0785923115
#> 421 0.0122070414 0.0217112786 0.0218948102 0.0226600776 0.0253732839
#> 422 0.0298906025 0.0092320302 0.0096986107 0.0109931613 0.0234980996
#> 423 0.0089720781 0.0283034969 0.0322928702 0.0416022596 0.0436441524
#> 424 0.0315455749 0.0011944828 0.0118895832 0.0132599729 0.0227043520
#> 425 0.0087400196 0.0184796109 0.0223777190 0.0245500246 0.0254071970
#> 444 0.0134398038 0.0284813433 0.0447007346 0.0516737345 0.0566152168
#> 445 0.0276997324 0.0284709685 0.0406786390 0.0428727658 0.0431321420
#> 446 0.0282003716 0.0293157193 0.0307516889 0.0427483407 0.0456852423
#> 447 0.0224113737 0.0251079118 0.0262917999 0.0295464171 0.0313229662
#> 448 0.0409331297 0.0437627160 0.0503186467 0.0507519118 0.0518311632
#> 449 0.0272580440 0.0263935995 0.0144972637 0.0053777769 0.0054039454
#> 450 0.0001756438 0.0014343851 0.0100578705 0.0360107864 0.0388306402
#> 470 0.0445103217 0.0554740975 0.0571090252 0.0571161807 0.0611761782
#> 471 0.0189248359 0.0367360847 0.0442011474 0.0442303093 0.0554308699
#> 472 0.0084579310 0.0156568307 0.0201006698 0.0201712641 0.0203695087
#> 473 0.0049307896 0.0158639237 0.0159174660 0.0230898647 0.0307370963
#> 474 0.0076624429 0.0152990563 0.0153514531 0.0245945182 0.0312609606
#> 475 0.0021845352 0.0154034607 0.0154247225 0.0260839471 0.0389926290
#> 496 0.0099463505 0.0100791279 0.0180436858 0.0180763858 0.0181833190
#> 497 0.0335473556 0.0336811087 0.0404175182 0.0404291993 0.0426503925
#> 498 0.0808657651 0.0811358959 0.0829363563 0.0870125172 0.0871384097
#> 499 0.0464940775 0.0493998914 0.0507631885 0.0545835986 0.0547739133
#> 500 0.0911715635 0.0911728541 0.0931281075 0.1002854734 0.1010668986
#> 522 0.0384527449 0.0395937863 0.0397333333 0.0399141809 0.0469233185
#> 523 0.0408043494 0.0421102990 0.0460645064 0.0538868223 0.0538874019
#> 524 0.0178014675 0.0204739057 0.0275167351 0.0323250929 0.0333647773
#> 525 0.0388706604 0.0390488196 0.0431271990 0.0507078602 0.0510362758
#> 548 0.0456590319 0.0457598784 0.0532620114 0.0682722358 0.0814233480
#> 549 0.0188518765 0.0124928599 0.0079589449 0.0019987706 0.0002014544
#> 550 0.0536034734 0.0539003815 0.0573274218 0.0687903578 0.0781087460
#> 574 0.0055225446 0.0003673425 0.0004986317 0.0201359528 0.0208882182
#> 575 0.0669549784 0.0675625842 0.0675741960 0.1048028304 0.1051619930
#> 600 0.0078122451 0.0108283024 0.0108997140 0.0115553601 0.0135157493
#>         pcoa_9d    pcoa_10d tree_average
#> 2   0.053053576 0.054054217 8.091011e-02
#> 3   0.036682285 0.037536087 8.903622e-02
#> 4   0.011039372 0.011046331 3.764924e-03
#> 5   0.023428945 0.023430222 4.374243e-02
#> 6   0.058537914 0.058828275 3.442120e-02
#> 7   0.083043516 0.083237066 3.267223e-02
#> 8   0.054264722 0.054327057 3.015781e-02
#> 9   0.057392916 0.057574335 6.046101e-02
#> 10  0.028096420 0.028410081 3.396603e-02
#> 11  0.051829673 0.052022130 5.365185e-02
#> 12  0.034511902 0.035661523 1.245955e-01
#> 13  0.019369034 0.019369215 2.636504e-02
#> 14  0.050409401 0.050531358 3.226313e-02
#> 15  0.073565942 0.073913549 8.214708e-02
#> 16  0.022793831 0.023244607 1.932518e-02
#> 17  0.050942834 0.050952672 1.528467e-01
#> 18  0.071199434 0.072290120 1.724429e-01
#> 19  0.049389715 0.059101410 0.000000e+00
#> 20  0.053260835 0.053301009 1.401818e-02
#> 21  0.039312429 0.047502556 1.099061e-01
#> 22  0.022194983 0.024927428 7.663176e-02
#> 23  0.033488412 0.033579516 3.072559e-02
#> 24  0.040257060 0.041504555 7.177487e-03
#> 25  0.045675834 0.046034289 6.290105e-02
#> 28  0.025002008 0.025006667 6.023881e-02
#> 29  0.025243178 0.026340931 6.975253e-02
#> 30  0.045490744 0.046415906 1.172599e-01
#> 31  0.069693271 0.069795079 1.781211e-02
#> 32  0.098251669 0.101199118 1.723115e-02
#> 33  0.017228790 0.017628827 2.972329e-02
#> 34  0.020886656 0.022418645 5.705887e-02
#> 35  0.053705784 0.053992008 5.181836e-02
#> 36  0.021386485 0.022952352 5.024971e-02
#> 37  0.036131991 0.036194022 1.213904e-01
#> 38  0.069045211 0.069718200 8.042053e-02
#> 39  0.090289582 0.090476885 7.389719e-02
#> 40  0.047612951 0.047693771 3.473951e-02
#> 41  0.037644794 0.040174835 7.727277e-02
#> 42  0.034269668 0.034767615 3.033936e-02
#> 43  0.070334875 0.076121032 6.714119e-03
#> 44  0.059229597 0.060296841 1.003573e-01
#> 45  0.018851227 0.019098323 1.653185e-01
#> 46  0.053960811 0.063245766 0.000000e+00
#> 47  0.046807305 0.047284250 1.077250e-01
#> 48  0.058350721 0.059711872 3.757637e-02
#> 49  0.048723284 0.048723999 5.323015e-02
#> 50  0.080347975 0.080439088 5.308517e-02
#> 54  0.014554275 0.015176319 1.887152e-01
#> 55  0.020052507 0.020575581 1.412078e-01
#> 56  0.014523260 0.014620490 1.990139e-01
#> 57  0.055264695 0.056479933 9.059457e-02
#> 58  0.033358590 0.033703732 7.076648e-02
#> 59  0.037165891 0.039343286 1.828202e-01
#> 60  0.030165393 0.030443701 4.169745e-02
#> 61  0.031407089 0.033199943 1.062960e-01
#> 62  0.045715755 0.045753704 1.184887e-01
#> 63  0.057495067 0.058316754 7.519333e-03
#> 64  0.083062410 0.083462635 1.526240e-01
#> 65  0.041651332 0.041764828 5.328442e-03
#> 66  0.040778765 0.043257568 1.812933e-01
#> 67  0.049486380 0.050265561 6.124315e-02
#> 68  0.021055795 0.023506381 7.051539e-03
#> 69  0.047858577 0.048498410 1.261547e-01
#> 70  0.170447933 0.171062617 0.000000e+00
#> 71  0.022863535 0.024262107 3.124281e-02
#> 72  0.024411766 0.024594817 1.083185e-01
#> 73  0.010467713 0.011389563 1.784671e-01
#> 74  0.037737958 0.037741626 1.571202e-01
#> 75  0.048924799 0.049028969 3.871589e-02
#> 80  0.032377088 0.032415359 6.701435e-02
#> 81  0.104745836 0.105336382 0.000000e+00
#> 82  0.070760368 0.070886219 6.235854e-02
#> 83  0.021840356 0.021933907 4.710630e-02
#> 84  0.028375966 0.028448587 9.892254e-02
#> 85  0.023810530 0.024011942 2.713322e-02
#> 86  0.017972741 0.018049880 9.211338e-02
#> 87  0.029208472 0.030259854 1.632540e-01
#> 88  0.032873980 0.032876766 1.313892e-01
#> 89  0.025872881 0.025994323 2.212459e-01
#> 90  0.047572613 0.047948337 8.720252e-02
#> 91  0.017810208 0.017974196 5.112882e-03
#> 92  0.023397945 0.023419623 8.197541e-02
#> 93  0.024091963 0.024560619 2.178024e-02
#> 94  0.017254181 0.020787761 1.568230e-02
#> 95  0.029072736 0.029128197 1.826837e-01
#> 96  0.011972434 0.018561220 4.075653e-02
#> 97  0.033806294 0.038945670 2.313601e-02
#> 98  0.044227411 0.044297652 3.085098e-02
#> 99  0.014923742 0.015823356 3.144488e-02
#> 100 0.066386095 0.066851004 1.285159e-02
#> 106 0.082473208 0.082804219 4.591924e-02
#> 107 0.096660658 0.096887877 9.370427e-02
#> 108 0.040045201 0.040086390 4.010489e-04
#> 109 0.024946234 0.025101607 5.004156e-02
#> 110 0.025289715 0.025417090 2.037413e-02
#> 111 0.043832294 0.044007854 1.239947e-02
#> 112 0.049783699 0.050798439 5.874120e-02
#> 113 0.035315392 0.035316997 8.388183e-02
#> 114 0.026316666 0.026390526 1.737385e-01
#> 115 0.058317927 0.058594414 3.969517e-02
#> 116 0.024054947 0.024358844 4.239447e-02
#> 117 0.027452538 0.027455382 3.446805e-02
#> 118 0.026061274 0.026776276 2.572711e-02
#> 119 0.022031103 0.025537459 6.318965e-02
#> 120 0.037097768 0.037120847 1.351764e-01
#> 121 0.026942518 0.033541386 8.826388e-02
#> 122 0.037284763 0.046758046 0.000000e+00
#> 123 0.088546004 0.088768862 8.363858e-03
#> 124 0.024843667 0.025612303 1.606247e-02
#> 125 0.099851981 0.100178482 6.035895e-02
#> 132 0.112674324 0.113583419 7.630223e-02
#> 133 0.058342033 0.058410612 8.164746e-02
#> 134 0.021834924 0.022407506 1.936743e-01
#> 135 0.033602065 0.033621983 8.479709e-02
#> 136 0.028549774 0.029127134 1.868652e-01
#> 137 0.026507451 0.026727797 2.580058e-01
#> 138 0.047373932 0.047550561 1.659303e-01
#> 139 0.028263910 0.028275681 2.557870e-01
#> 140 0.094185897 0.094186454 5.962524e-02
#> 141 0.043588673 0.044462755 5.934268e-02
#> 142 0.056303130 0.056438807 4.378929e-02
#> 143 0.043206122 0.044712018 4.469552e-02
#> 144 0.051688065 0.053046343 3.398769e-04
#> 145 0.092423930 0.092484967 6.106299e-02
#> 146 0.057214966 0.060321909 2.862630e-02
#> 147 0.077948166 0.079106321 4.737334e-02
#> 148 0.083877914 0.085091974 3.399321e-02
#> 149 0.050536865 0.050638341 6.598604e-02
#> 150 0.096084161 0.096084947 2.819944e-02
#> 158 0.068485714 0.068930444 1.251919e-02
#> 159 0.035636144 0.035637573 1.801094e-01
#> 160 0.053465741 0.054040178 7.123220e-02
#> 161 0.036917276 0.036918485 1.733003e-01
#> 162 0.019481635 0.021129803 2.444409e-01
#> 163 0.109342952 0.109491490 1.226630e-01
#> 164 0.129177780 0.129643205 3.280634e-02
#> 165 0.035755550 0.036443895 1.577901e-01
#> 166 0.047190335 0.047193369 4.577779e-02
#> 167 0.048127782 0.048344951 1.680706e-02
#> 168 0.107917547 0.108085607 4.527648e-02
#> 169 0.065890560 0.070788641 1.322501e-02
#> 170 0.028535944 0.028770395 2.259190e-01
#> 171 0.070666807 0.086048065 2.804533e-02
#> 172 0.087264006 0.091185411 9.225017e-02
#> 173 0.123034793 0.123051833 6.250561e-02
#> 174 0.054796155 0.056434274 3.914276e-02
#> 175 0.067841042 0.068608295 6.996539e-02
#> 184 0.044648874 0.045019919 2.908850e-02
#> 185 0.102755657 0.102782396 1.504958e-01
#> 186 0.029412661 0.029772144 6.268338e-02
#> 187 0.032707087 0.033394177 9.342001e-02
#> 188 0.108974021 0.109033828 1.461499e-01
#> 189 0.071787417 0.071802924 1.399339e-02
#> 190 0.076333767 0.076421650 2.281046e-02
#> 191 0.033806524 0.034309301 3.819121e-02
#> 192 0.056717600 0.056735330 5.885737e-02
#> 193 0.016796199 0.017762344 7.000870e-02
#> 194 0.066875920 0.069615964 8.343786e-02
#> 195 0.016064020 0.016064111 2.121515e-01
#> 196 0.013666787 0.018096292 2.496971e-02
#> 197 0.046056201 0.047540193 9.047958e-03
#> 198 0.058280266 0.058509320 6.110064e-02
#> 199 0.051239463 0.051762055 7.147413e-02
#> 200 0.059573651 0.059648971 1.228646e-01
#> 210 0.037969013 0.038855116 1.918284e-02
#> 211 0.074392539 0.074392568 0.000000e+00
#> 212 0.057605751 0.060748323 5.333763e-02
#> 213 0.064661720 0.064780958 8.105531e-02
#> 214 0.095687839 0.096349082 2.068552e-01
#> 215 0.094827580 0.096012909 1.641554e-01
#> 216 0.067554483 0.067574231 8.706571e-05
#> 217 0.040505547 0.040707603 5.824687e-02
#> 218 0.044287790 0.044461135 5.448761e-03
#> 219 0.039770135 0.047236299 7.471059e-03
#> 220 0.062085739 0.062426346 9.592239e-02
#> 221 0.014946656 0.021911591 6.787306e-02
#> 222 0.025060840 0.028188051 4.139466e-02
#> 223 0.031662169 0.031668664 1.115433e-01
#> 224 0.078192959 0.081077162 6.722097e-03
#> 225 0.084961630 0.085860105 5.814886e-02
#> 236 0.060580847 0.061424189 1.237368e-02
#> 237 0.040348400 0.040904344 8.331738e-02
#> 238 0.031316601 0.031436391 1.713816e-02
#> 239 0.051790497 0.051790536 4.149002e-02
#> 240 0.095126978 0.095165575 1.234252e-01
#> 241 0.049887151 0.051096430 8.684429e-03
#> 242 0.057966332 0.058043992 1.115685e-01
#> 243 0.071713236 0.073522904 1.005897e-01
#> 244 0.080042695 0.084390322 3.396603e-02
#> 245 0.052877263 0.052893653 4.520038e-02
#> 246 0.037400729 0.041949012 4.100418e-02
#> 247 0.024303390 0.025603908 2.902104e-02
#> 248 0.028384032 0.028791967 4.112756e-02
#> 249 0.046123773 0.046539370 2.049397e-03
#> 250 0.062371241 0.062405394 3.953275e-02
#> 262 0.074476899 0.078445113 1.298618e-01
#> 263 0.050800469 0.050917008 4.746044e-02
#> 264 0.085759741 0.086290730 1.060886e-01
#> 265 0.053026640 0.053979955 4.722714e-02
#> 266 0.078760771 0.078778684 6.722097e-03
#> 267 0.028306027 0.028520407 6.505603e-02
#> 268 0.039526115 0.039701856 1.225792e-02
#> 269 0.055394142 0.062621427 1.428022e-02
#> 270 0.051607595 0.051906484 1.939822e-02
#> 271 0.005889322 0.013154003 6.106389e-02
#> 272 0.037027984 0.040684805 3.203539e-02
#> 273 0.049447565 0.049455278 4.933116e-02
#> 274 0.086777573 0.089655707 8.706571e-05
#> 275 0.041548972 0.042339137 4.261774e-02
#> 288 0.095366137 0.096760229 1.833905e-01
#> 289 0.114182488 0.114694526 1.425237e-01
#> 290 0.056008504 0.056277755 6.684282e-02
#> 291 0.032587638 0.035865327 2.339327e-03
#> 292 0.071180818 0.072325638 1.605820e-01
#> 293 0.020438712 0.023429818 1.422161e-01
#> 294 0.048282037 0.048642785 5.666348e-02
#> 295 0.057625245 0.058294143 3.159089e-02
#> 296 0.023094669 0.024210603 1.322046e-01
#> 297 0.055779054 0.055873688 3.910528e-02
#> 298 0.051527105 0.053025277 1.204718e-01
#> 299 0.057710094 0.057804123 7.903763e-02
#> 300 0.033610448 0.033838932 1.728493e-01
#> 314 0.173540898 0.173731070 0.000000e+00
#> 315 0.042215659 0.042429049 1.070933e-01
#> 316 0.040854581 0.041054518 7.195261e-02
#> 317 0.040958045 0.040969849 3.261669e-03
#> 318 0.053780007 0.054264151 4.013512e-02
#> 319 0.023118790 0.025738791 2.691501e-02
#> 320 0.090649140 0.090686519 1.489043e-01
#> 321 0.055583289 0.060427192 8.517411e-02
#> 322 0.030985455 0.032510531 7.523492e-02
#> 323 0.035058154 0.035104872 1.453835e-01
#> 324 0.057447891 0.058393593 1.816179e-01
#> 325 0.012301392 0.012494028 2.071475e-01
#> 340 0.104877985 0.104899870 1.870655e-02
#> 341 0.074715057 0.075479920 1.487626e-01
#> 342 0.096293163 0.096364510 3.261669e-03
#> 343 0.083084515 0.084332369 1.169451e-01
#> 344 0.057655066 0.059299079 8.554319e-02
#> 345 0.131056297 0.131075846 6.942826e-03
#> 346 0.084091000 0.087160996 7.865077e-02
#> 347 0.022751158 0.023462307 1.650916e-01
#> 348 0.014853816 0.015107383 2.352402e-01
#> 349 0.077035562 0.077293940 1.569128e-01
#> 350 0.096170698 0.096190924 4.700413e-02
#> 366 0.039753530 0.040938283 4.241528e-02
#> 367 0.006856455 0.006998886 1.299018e-01
#> 368 0.058370126 0.060159975 2.776813e-02
#> 369 0.073760260 0.075637766 1.500791e-01
#> 370 0.111758643 0.111857623 1.299960e-01
#> 371 0.036048192 0.038710922 4.555369e-02
#> 372 0.057903507 0.058744864 3.104826e-02
#> 373 0.094188325 0.094936069 9.914252e-03
#> 374 0.051710746 0.051821261 3.578025e-02
#> 375 0.141619721 0.141619726 0.000000e+00
#> 392 0.032603199 0.032991335 1.562149e-01
#> 393 0.042067144 0.042201789 1.044283e-01
#> 394 0.024038507 0.035097530 5.493401e-02
#> 395 0.057006416 0.057453461 7.621363e-02
#> 396 0.017784769 0.028774476 6.645858e-02
#> 397 0.017950508 0.022454961 5.104138e-02
#> 398 0.034830874 0.034865593 1.910722e-02
#> 399 0.050263068 0.054831045 0.000000e+00
#> 400 0.050199683 0.051417508 6.155309e-02
#> 418 0.030226528 0.030867849 2.237719e-02
#> 419 0.020711212 0.023358364 8.491465e-02
#> 420 0.079537335 0.079550773 5.050066e-02
#> 421 0.026651269 0.031078276 5.933537e-02
#> 422 0.023683336 0.025158004 1.774034e-02
#> 423 0.043698347 0.043804520 4.748490e-02
#> 424 0.039104131 0.039919322 1.992135e-01
#> 425 0.026296902 0.026439369 1.148044e-01
#> 444 0.056641249 0.068199604 1.676477e-01
#> 445 0.043134078 0.043822826 7.378568e-02
#> 446 0.045687162 0.069779717 6.714119e-03
#> 447 0.034454169 0.040267121 5.861644e-02
#> 448 0.052583049 0.052816225 1.271027e-02
#> 449 0.012267685 0.015673875 1.294474e-02
#> 450 0.044705381 0.046650552 4.488573e-02
#> 470 0.061179547 0.063024702 3.925683e-02
#> 471 0.055463515 0.055964130 1.254614e-01
#> 472 0.023516622 0.023611127 9.607898e-02
#> 473 0.031353216 0.035699064 2.593038e-02
#> 474 0.043744125 0.045165781 4.610257e-02
#> 475 0.044894795 0.046545755 8.234827e-02
#> 496 0.018187371 0.020917409 1.379250e-01
#> 497 0.044282285 0.045281382 1.022870e-01
#> 498 0.087594690 0.087774423 3.708211e-02
#> 499 0.059646103 0.060012920 5.350994e-02
#> 500 0.105223312 0.105314139 7.842669e-02
#> 522 0.050671382 0.051591735 1.211532e-01
#> 523 0.054763573 0.061875827 5.100461e-02
#> 524 0.040454493 0.043034278 5.798373e-02
#> 525 0.056545635 0.059529412 2.408917e-02
#> 548 0.082827439 0.090116791 8.363858e-03
#> 549 0.018170381 0.018543741 2.470938e-02
#> 550 0.078639673 0.079582140 3.264222e-02
#> 574 0.030308156 0.031549853 4.543922e-02
#> 575 0.107644274 0.108508316 7.360473e-02
#> 600 0.033227530 0.033320848 4.811224e-02
#> 
#> 
  
# Retrieve Species Coordinates
sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord
sp_faxes_coord_fruits
#>                         PC1          PC2          PC3          PC4          PC5
#> apple          0.0265193144  0.013054128 -0.043717435 -0.021211462 -0.148218537
#> apricot        0.1032240731 -0.148535464 -0.084372792 -0.059069034  0.051197559
#> banana        -0.3346937349  0.091371413 -0.070877547  0.069888722 -0.050000845
#> currant        0.2941221208 -0.002028080  0.059840313  0.027477102  0.036184138
#> blackberry     0.2398366750  0.002921130  0.096332589  0.012767053  0.063869252
#> blueberry      0.3146904236  0.065534392 -0.061375184  0.170320608 -0.043827091
#> cherry         0.1532479527 -0.250996898 -0.135570593  0.056216685  0.026507879
#> grape          0.0320396363 -0.175761176  0.131682624 -0.063280205 -0.101388849
#> grapefruit    -0.1565337233  0.091667465  0.051101664 -0.129100040  0.021725499
#> kiwifruit     -0.0076484035 -0.001811478  0.049335316 -0.133654078 -0.149411856
#> lemon         -0.1116123197  0.052718824  0.138487467 -0.103193749  0.023077059
#> lime          -0.2106569060  0.003790913  0.219519143  0.031204379  0.104197631
#> litchi        -0.2423438344 -0.252455139  0.019104502  0.089864709  0.130046060
#> mango         -0.3558930629 -0.084500319 -0.155978719 -0.069529132  0.072613533
#> melon          0.0475710274  0.228214444 -0.009073564 -0.130154617  0.037282406
#> orange        -0.0673938482  0.026151045 -0.002061839 -0.019971288 -0.014393513
#> passion_fruit -0.1743734750 -0.091951801  0.051587716  0.220532013 -0.147452934
#> peach          0.0403146225 -0.071052366 -0.144444375 -0.097071825  0.002288305
#> pear          -0.0008010876  0.019689428 -0.013204612 -0.070565107 -0.070145826
#> pineapple     -0.2308886096  0.336194377 -0.081546943  0.178095824  0.031816021
#> plum           0.0896173933 -0.145183768 -0.102585094 -0.018724455  0.019074804
#> raspberry      0.2180001689  0.011254750  0.088794687  0.026582440  0.044477454
#> strawberry     0.2743910893  0.093250971  0.046622086  0.093318143  0.049540096
#> tangerine     -0.0989381984 -0.035898731  0.072512352  0.008770777 -0.019765660
#> water_melon    0.1582027062  0.224361939 -0.120111760 -0.069513464  0.030707416
#>                        PC6          PC7          PC8          PC9          PC10
#> apple          0.040965006  0.039965979 -0.022577267 -0.012017348 -1.067572e-02
#> apricot        0.059097707 -0.027540550  0.056358883 -0.016901596  1.356995e-02
#> banana         0.020752971 -0.041476147 -0.083039226  0.042337245  1.583090e-02
#> currant        0.005808690 -0.114447728 -0.019688896 -0.024109923 -1.297464e-02
#> blackberry     0.035900843 -0.024364564  0.003797170  0.038186173 -9.737454e-03
#> blueberry     -0.072940389 -0.011254549  0.007080864 -0.046110072  4.357324e-03
#> cherry        -0.020780891  0.004803343 -0.042620591  0.030307476 -2.266872e-02
#> grape         -0.195578235 -0.001035913 -0.003686485  0.033761747 -3.964539e-03
#> grapefruit    -0.015750045 -0.051616104 -0.001017932  0.004568718 -2.143000e-02
#> kiwifruit     -0.010962090  0.006877412  0.052415441  0.009613759 -3.049371e-05
#> lemon          0.101128620  0.044498483 -0.034082720 -0.023135838 -2.153532e-02
#> lime           0.053717374  0.022699643  0.035206514  0.015350425  2.127588e-02
#> litchi        -0.107160270  0.013020916 -0.019139655 -0.016630151 -1.110479e-02
#> mango         -0.031431302 -0.007224234  0.043541389  0.011834427  1.624666e-04
#> melon         -0.136045993  0.036709357  0.013459634 -0.022081148  5.067072e-03
#> orange         0.091930421 -0.061574888 -0.012476483 -0.021006738 -2.425390e-02
#> passion_fruit  0.066901090 -0.006120818  0.066457679  0.020143504 -8.008755e-03
#> peach          0.073541526  0.076191641 -0.010596386 -0.012269290 -3.260209e-02
#> pear          -0.004973365  0.023163490 -0.055891141 -0.008962102  3.803879e-02
#> pineapple     -0.044607834  0.025787798  0.006751938 -0.010770300 -4.320707e-03
#> plum           0.054398550  0.019423316  0.018261782 -0.013010921  5.390245e-02
#> raspberry      0.045504907  0.029550728 -0.045185074  0.033817649  3.060490e-02
#> strawberry     0.001191018  0.094872368  0.018914728  0.012390681 -1.886095e-02
#> tangerine     -0.009194489 -0.033661520 -0.006545485 -0.078234499  1.424240e-02
#> water_melon   -0.001413820 -0.057247459  0.034301318  0.052928120  5.115953e-03
```
