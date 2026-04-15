# Compute the Minimum Spanning Tree (MST) linking species of a given assemblage

This function computes the MST linking species of a given assemblage and
is used to compute FEve index.

## Usage

``` r
mst.computation(sp_faxes_coord_k)
```

## Arguments

- sp_faxes_coord_k:

  a matrix relating species coordinates for species present in a given
  assemblage.

## Value

A dist object summarizing the MST for all species of a given assemblage
`mst_asb_k`.

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
# Load Species*Traits dataframe:
 data("fruits_traits", package = "mFD")

# Load Assemblages*Species dataframe:      
 data("baskets_fruits_weights", package = "mFD") 

# Load Traits categories dataframe:
 data("fruits_traits_cat", package = "mFD") 
 
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
  
# Compute functional spaces quality to retrieve species coordinates matrix:
 fspaces_quality_fruits <- mFD::quality.fspaces(
                                  sp_dist             = sp_dist_fruits, 
                                  maxdim_pcoa         = 10,
                                  deviation_weighting = "absolute",
                                  fdist_scaling       = FALSE,
                                  fdendro             = "average")
 
# Retrieve species coordinates matrix:
 sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord

# Compute the distance of "pear" to its nearest neighbor(s):
 mst_fruits <- mst.computation(sp_faxes_coord_fruits)
 mst_fruits
#>               apple apricot banana currant blackberry blueberry cherry grape
#> apricot           0                                                         
#> banana            0       0                                                 
#> currant           0       0      0                                          
#> blackberry        0       0      0       1                                  
#> blueberry         0       0      0       0          0                       
#> cherry            0       0      0       0          0         0             
#> grape             0       0      0       0          0         0      0      
#> grapefruit        0       0      1       0          0         0      0     0
#> kiwifruit         0       0      0       0          0         0      0     1
#> lemon             0       0      0       0          0         0      0     0
#> lime              0       0      0       0          0         0      0     0
#> litchi            0       0      0       0          0         0      0     0
#> mango             0       0      1       0          0         0      0     0
#> melon             0       0      0       0          0         0      0     0
#> orange            0       0      0       0          0         0      0     0
#> passion_fruit     0       0      0       0          0         0      0     0
#> peach             1       0      0       0          0         0      0     0
#> pear              1       0      0       0          0         0      0     0
#> pineapple         0       0      1       0          0         0      0     0
#> plum              0       1      0       0          0         0      1     0
#> raspberry         0       0      0       0          1         0      0     0
#> strawberry        0       0      0       0          0         1      0     0
#> tangerine         0       0      0       0          0         0      0     0
#> water_melon       0       0      0       0          0         0      0     0
#>               grapefruit kiwifruit lemon lime litchi mango melon orange
#> apricot                                                                
#> banana                                                                 
#> currant                                                                
#> blackberry                                                             
#> blueberry                                                              
#> cherry                                                                 
#> grape                                                                  
#> grapefruit                                                             
#> kiwifruit              0                                               
#> lemon                  1         0                                     
#> lime                   0         0     1                               
#> litchi                 0         0     0    0                          
#> mango                  0         0     0    0      1                   
#> melon                  0         0     0    0      0     0             
#> orange                 1         0     0    0      0     0     0       
#> passion_fruit          0         0     0    0      0     0     0      0
#> peach                  0         0     0    0      0     0     0      0
#> pear                   0         1     0    0      0     0     1      1
#> pineapple              0         0     0    0      0     0     0      0
#> plum                   0         0     0    0      0     0     0      0
#> raspberry              0         0     0    0      0     0     0      0
#> strawberry             0         0     0    0      0     0     0      0
#> tangerine              0         0     0    0      0     0     0      1
#> water_melon            0         0     0    0      0     0     1      0
#>               passion_fruit peach pear pineapple plum raspberry strawberry
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
#> peach                     0                                               
#> pear                      0     0                                         
#> pineapple                 0     0    0                                    
#> plum                      0     1    0         0                          
#> raspberry                 0     0    1         0    0                     
#> strawberry                0     0    0         0    0         1           
#> tangerine                 1     0    0         0    0         0          0
#> water_melon               0     0    0         0    0         0          0
#>               tangerine
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
#> raspberry              
#> strawberry             
#> tangerine              
#> water_melon           0
```
