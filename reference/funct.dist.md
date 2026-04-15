# Compute functional distance between species

For a given combination of traits, this function returns the functional
distance matrix between species.

## Usage

``` r
funct.dist(
  sp_tr,
  tr_cat,
  metric,
  scale_euclid = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA = TRUE
)
```

## Arguments

- sp_tr:

  a data frame of traits values (columns) for each species (rows).

- tr_cat:

  a data frame containing three columns for each trait (rows):

  - **trait_name**: the name of all traits as in `sp_tr` data frame;

  - **trait_type**: the category code for each trait as followed: `N`
    for Nominal traits (factor variable), `O` for Ordinal traits
    (ordered variable), `C` for Circular traits (integer values), `Q`
    for quantitative traits (numeric values) that is allowed **only** if
    there are at least 2 species with the same value, and `F` for fuzzy
    traits (i.e. described with several values defined with several
    column);

  - **fuzzy_name**: name of fuzzy-coded trait to which 'sub-trait'
    belongs (if trait is not fuzzy, ignored so could be trait name or
    NA).

  - **trait_weight**: Optional, a numeric vector of length n (traits
    number) to specify a weight for each trait.

- metric:

  the distance to be computed: `euclidean`, the Euclidean distance,
  `gower`, the Classical Gower distance as defined by Gower (1971),
  extent by de Bello *et al.* (2021) and based on the
  [`gawdis`](https://rdrr.io/pkg/gawdis/man/gawdis.html) function.

- scale_euclid:

  only when computing euclidean distance a string value to compute (or
  not) scaling of quantitative traits using the
  [`tr.cont.scale`](https://cmlmagneville.github.io/mFD/reference/tr.cont.scale.md)
  function. Possible options are: `range` (standardize by the range:
  \\({x' = x - min(x) )} / (max(x) - min (x))\\) `center` (use the
  center transformation: \\x' = x - mean(x)\\), `scale` (use the scale
  transformation: \\x' = \frac{x}{sd(x)}\\), `scale_center` (use the
  scale-center transformation: \\x' = \frac{x - mean(x)}{sd(x)}\\), or
  `noscale` traits are not scaled Default is `scale_center`.

- ordinal_var:

  a character string specifying the method to be used for ordinal
  variables (i.e. ordered). `classic` simply treats ordinal variables as
  continuous variables; `metric` refers to Eq. 3 of Podani (1999);
  `podani` refers to Eqs. 2a-b of Podani (1999), Both options convert
  ordinal variables to ranks. Default is `classic`.

- weight_type:

  the type of used method to weight traits. `user` user defined weights
  in tr_cat, `equal` all traits having the same weight. More methods are
  available using [`gawdis`](https://rdrr.io/pkg/gawdis/man/gawdis.html)
  from `gawdis` package. To compute gower distance with fuzzy trait and
  weight please refer to
  [`gawdis`](https://rdrr.io/pkg/gawdis/man/gawdis.html). Default is
  `equal`.

- stop_if_NA:

  a logical value to stop or not the process if the `sp_tr` data frame
  contains NA. Functional measures are sensitive to missing traits. For
  further explanations, see the Note section. Default is `TRUE`.

## Value

a `dist` object containing distance between each pair of species.

## Note

If the `sp_tr` data frame contains `NA` you can either chose to compute
anyway functional distances (but keep in mind that **Functional measures
are sensitive to missing traits!**) or you can delete species with
missing or extrapolate missing traits (see Johnson *et al.* (2020)).

## References

de Bello *et al.* (2021) Towards a more balanced combination of multiple
traits when computing functional differences between species. *Method in
Ecology and Evolution*, **12**, 443-448.  
Gower (1971 ) A general coefficient of similarity and some of its
properties. *Biometrics*, **27**, 857-871.  
Johnson *et al.* (2020) Handling missing values in trait data. *Global
Ecology and Biogeography*, **30**, 51-62.  
Podani (1999) Extending Gower's general coefficient of similarity to
ordinal characters, *Taxon*, **48**, 331-340.

## Author

Nicolas Loiseau and Sebastien Villeger

## Examples

``` r
# Load Species x Traits data
data("fruits_traits", package = "mFD")

# Load Traits x Categories data
data("fruits_traits_cat", package = "mFD")

# Remove fuzzy traits for this example and thus remove lat column:
fruits_traits     <- fruits_traits[ , -c(6:8)]
fruits_traits_cat <- fruits_traits_cat[-c(6:8), ]
fruits_traits_cat <- fruits_traits_cat[ , -3]

# Compute Functional Distance
sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
                                  tr_cat        = fruits_traits_cat,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)
sp_dist_fruits
#>                     apple     apricot      banana     currant  blackberry
#> apricot       0.165796703                                                
#> banana        0.375274725 0.541071429                                    
#> currant       0.391483516 0.425686813 0.766758242                        
#> blackberry    0.375686813 0.409890110 0.750961538 0.084203297            
#> blueberry     0.355357143 0.410439560 0.730631868 0.236126374 0.320329670
#> cherry        0.233379121 0.099175824 0.558104396 0.424862637 0.409065934
#> grape         0.380494505 0.446291209 0.705219780 0.371978022 0.356181319
#> grapefruit    0.192307692 0.326510989 0.267582418 0.500824176 0.483379121
#> kiwifruit     0.219230769 0.353434066 0.594505495 0.372252747 0.356456044
#> lemon         0.208379121 0.342582418 0.383653846 0.516895604 0.432692308
#> lime          0.369505495 0.403708791 0.344780220 0.578021978 0.493818681
#> litchi        0.466483516 0.332280220 0.391208791 0.657967033 0.642170330
#> mango         0.394917582 0.360714286 0.219642857 0.786401099 0.770604396
#> melon         0.284752747 0.418956044 0.560027473 0.406730769 0.390934066
#> orange        0.117170330 0.251373626 0.292445055 0.474313187 0.458516484
#> passion_fruit 0.461126374 0.526923077 0.414148352 0.552609890 0.536813187
#> peach         0.127472527 0.061675824 0.502747253 0.464010989 0.448214286
#> pear          0.008791209 0.157005495 0.384065934 0.382692308 0.366895604
#> pineapple     0.557417582 0.708379121 0.232692308 0.734065934 0.718269231
#> plum          0.156456044 0.009340659 0.531730769 0.435027473 0.419230769
#> raspberry     0.382280220 0.416483516 0.757554945 0.090796703 0.006593407
#> strawberry    0.375549451 0.409752747 0.750824176 0.284065934 0.200137363
#> tangerine     0.152609890 0.218406593 0.322664835 0.444093407 0.428296703
#> water_melon   0.281181319 0.415384615 0.556456044 0.410302198 0.394505495
#>                 blueberry      cherry       grape  grapefruit   kiwifruit
#> apricot                                                                  
#> banana                                                                   
#> currant                                                                  
#> blackberry                                                               
#> blueberry                                                                
#> cherry        0.388736264                                                
#> grape         0.335851648 0.347115385                                    
#> grapefruit    0.536950549 0.425686813 0.572802198                        
#> kiwifruit     0.363873626 0.452609890 0.199725275 0.373076923            
#> lemon         0.553021978 0.441758242 0.588873626 0.116071429 0.389148352
#> lime          0.614148352 0.502884615 0.650000000 0.277197802 0.550274725
#> litchi        0.621840659 0.233104396 0.514010989 0.458791209 0.685714286
#> mango         0.750274725 0.361538462 0.685576923 0.287225275 0.614148352
#> melon         0.229395604 0.518131868 0.465247253 0.307554945 0.265521978
#> orange        0.461813187 0.350549451 0.497664835 0.075137363 0.302060440
#> passion_fruit 0.516483516 0.572252747 0.319368132 0.453434066 0.280357143
#> peach         0.472115385 0.160851648 0.507967033 0.264835165 0.308241758
#> pear          0.353434066 0.242170330 0.389285714 0.183516484 0.210439560
#> pineapple     0.502060440 0.790796703 0.737912088 0.434890110 0.561813187
#> plum          0.401098901 0.089835165 0.436950549 0.335851648 0.362774725
#> raspberry     0.326923077 0.415659341 0.362774725 0.489972527 0.363049451
#> strawberry    0.120192308 0.408928571 0.356043956 0.483241758 0.356318681
#> tangerine     0.407967033 0.280769231 0.427884615 0.144917582 0.371840659
#> water_melon   0.225824176 0.514560440 0.461675824 0.311126374 0.261950549
#>                     lemon        lime      litchi       mango       melon
#> apricot                                                                  
#> banana                                                                   
#> currant                                                                  
#> blackberry                                                               
#> blueberry                                                                
#> cherry                                                                   
#> grape                                                                    
#> grapefruit                                                               
#> kiwifruit                                                                
#> lemon                                                                    
#> lime          0.161126374                                                
#> litchi        0.474862637 0.335989011                                    
#> mango         0.403296703 0.364423077 0.171565934                        
#> melon         0.423626374 0.584752747 0.751236264 0.579670330            
#> orange        0.091208791 0.252335165 0.383653846 0.312087912 0.367582418
#> passion_fruit 0.469505495 0.330631868 0.405357143 0.433791209 0.545879121
#> peach         0.280906593 0.442032967 0.393956044 0.322390110 0.357280220
#> pear          0.199587912 0.360714286 0.475274725 0.403708791 0.275961538
#> pineapple     0.550961538 0.512087912 0.623901099 0.452335165 0.327335165
#> plum          0.351923077 0.413049451 0.322939560 0.351373626 0.428296703
#> raspberry     0.426098901 0.487225275 0.648763736 0.777197802 0.397527473
#> strawberry    0.432829670 0.493956044 0.642032967 0.770467033 0.190796703
#> tangerine     0.160989011 0.222115385 0.313873626 0.342307692 0.437362637
#> water_melon   0.427197802 0.588324176 0.747664835 0.576098901 0.003571429
#>                    orange passion_fruit       peach        pear   pineapple
#> apricot                                                                    
#> banana                                                                     
#> currant                                                                    
#> blackberry                                                                 
#> blueberry                                                                  
#> cherry                                                                     
#> grape                                                                      
#> grapefruit                                                                 
#> kiwifruit                                                                  
#> lemon                                                                      
#> lime                                                                       
#> litchi                                                                     
#> mango                                                                      
#> melon                                                                      
#> orange                                                                     
#> passion_fruit 0.378296703                                                  
#> peach         0.210302198   0.588598901                                    
#> pear          0.108379121   0.469917582 0.118681319                        
#> pineapple     0.459752747   0.418543956 0.670054945 0.551373626            
#> plum          0.260714286   0.517582418 0.071016484 0.152335165 0.700961538
#> raspberry     0.465109890   0.543406593 0.454807692 0.373489011 0.724862637
#> strawberry    0.458379121   0.536675824 0.448076923 0.366758242 0.518131868
#> tangerine     0.069780220   0.308516484 0.280082418 0.161401099 0.510027473
#> water_melon   0.364010989   0.542307692 0.353708791 0.272390110 0.323763736
#>                      plum   raspberry  strawberry   tangerine
#> apricot                                                      
#> banana                                                       
#> currant                                                      
#> blackberry                                                   
#> blueberry                                                    
#> cherry                                                       
#> grape                                                        
#> grapefruit                                                   
#> kiwifruit                                                    
#> lemon                                                        
#> lime                                                         
#> litchi                                                       
#> mango                                                        
#> melon                                                        
#> orange                                                       
#> passion_fruit                                                
#> peach                                                        
#> pear                                                         
#> pineapple                                                    
#> plum                                                         
#> raspberry     0.425824176                                    
#> strawberry    0.419093407 0.206730769                        
#> tangerine     0.209065934 0.434890110 0.428159341            
#> water_melon   0.424725275 0.401098901 0.194368132 0.433791209
```
