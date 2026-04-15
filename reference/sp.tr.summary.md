# Summarize Species x Traits data frame

This function computes a summary data helping to choose the type of
analysis you can do with your data. For this function to work, there
must be no NA in your `sp_tr` data frame.

## Usage

``` r
sp.tr.summary(tr_cat, sp_tr, stop_if_NA = TRUE)
```

## Arguments

- tr_cat:

  a data frame containing three columns for each trait (rows):

  - **trait_name**: the name of all traits as in `sp_tr` data frame;

  - **trait_type**: the category code for each trait as followed: `N`
    for Nominal traits (factor variable), `O` for Ordinal traits
    (ordered variable), `C` for Circular traits (integer
    values)(circular traits can not be used in mFD function used to
    compute functional distance but ok for summary function and function
    to group species into Functional Entities), `Q` for quantitative
    traits (numeric values), and `F` for fuzzy traits (i.e. described
    with several values defined with several column);

  - **fuzzy_name**: name of fuzzy-coded trait to which 'sub-trait'
    belongs (if trait is not fuzzy, ignored so could be trait name or
    NA).

  - **trait_weight**: weights of each traits if the user wants to
    specify a weight for each trait.

- sp_tr:

  a data frame of traits values (columns) for each species (rows). Note
  that species names **must be** specified in the row names.

- stop_if_NA:

  a logical value indicating whether the process should stop if there is
  some NA in the `sp_tr` dataframe. If you continue with functional
  analysis, we remind you that functional measures, are sensitive to
  missing traits

## Value

If there is no fuzzy-coded trait, a three-elements list with:

- tr_summary_list:

  a table summarizing for each trait the number of species per modality
  for non-continuous trait and min, max, mean, median, and quartiles for
  continuous traits.

- tr_types:

  a list containing traits type.

- mod_list:

  a list containing modalities for all traits.

If there is fuzzy-coded trait, a four-elements list with:

- tr_summary_non_fuzzy_list:

  a table summarizing for each trait the number of species per modality
  for non-continuous trait and min, max, mean, median, and quartiles for
  continuous traits.

- tr_summary_fuzzy_list:

  a table summarizing for each subtrait min, max, mean, median and
  quartiles

- tr_types:

  a list containing traits type.

- mod_list:

  a list containing modalities for non-continuous trait.

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
# Load Species x Traits data
data('fruits_traits', package = 'mFD')

# Load Traits x Categories data
data('fruits_traits_cat', package = 'mFD')

# Summarize Species x Traits data
mFD::sp.tr.summary(tr_cat = fruits_traits_cat, sp_tr = fruits_traits)
#> $tr_summary_non_fuzzy_list
#>       Size     Plant           Climate     Seed        Sugar       
#>  0-1cm  :2   forb : 5   temperate  :15   none: 2   Min.   : 16.90  
#>  1-3cm  :6   shrub: 3   subtropical: 4   pip :17   1st Qu.: 73.70  
#>  3-5cm  :5   tree :14   tropical   : 6   pit : 6   Median : 92.40  
#>  5-10cm :6   vine : 3                              Mean   : 90.66  
#>  10-20cm:6                                         3rd Qu.:105.80  
#>                                                    Max.   :162.50  
#> 
#> $tr_summary_fuzzy_list
#>     Use.raw        Use.pastry      Use.jam  
#>  Min.   : 10.0   Min.   : 0.0   Min.   : 0  
#>  1st Qu.: 40.0   1st Qu.: 0.0   1st Qu.: 0  
#>  Median : 60.0   Median :10.0   Median :10  
#>  Mean   : 61.2   Mean   :16.8   Mean   :22  
#>  3rd Qu.: 90.0   3rd Qu.:30.0   3rd Qu.:40  
#>  Max.   :100.0   Max.   :60.0   Max.   :80  
#> 
#> $tr_types
#> $tr_types$Size
#> [1] "ordered" "factor" 
#> 
#> $tr_types$Plant
#> [1] "factor"
#> 
#> $tr_types$Climate
#> [1] "ordered" "factor" 
#> 
#> $tr_types$Seed
#> [1] "ordered" "factor" 
#> 
#> $tr_types$Sugar
#> [1] "numeric"
#> 
#> 
#> $mod_list
#> $mod_list$Size
#> [1] 5-10cm  3-5cm   10-20cm 0-1cm   1-3cm  
#> Levels: 0-1cm < 1-3cm < 3-5cm < 5-10cm < 10-20cm
#> 
#> $mod_list$Plant
#> [1] tree  shrub forb  vine 
#> Levels: forb shrub tree vine
#> 
#> $mod_list$Climate
#> [1] temperate   tropical    subtropical
#> Levels: temperate < subtropical < tropical
#> 
#> $mod_list$Seed
#> [1] pip  pit  none
#> Levels: none < pip < pit
#> 
#> 
```
