# Compute Functional Entities composition based on a Species x Traits matrix

Compute Functional Entities composition based on a Species x Traits
matrix

## Usage

``` r
sp.to.fe(sp_tr, tr_cat, fe_nm_type = "fe_rank", check_input = TRUE)
```

## Arguments

- sp_tr:

  a data frame containing species as rows and traits as columns.

- tr_cat:

  a data frame containing three columns for each trait (rows):

  - **trait_name**: names of all traits as in `sp_tr` data frame;

  - **trait_type**: category codes for each trait as followed: *N* for
    Nominal traits (factor variable), *O* for Ordinal traits (ordered
    variable), *C* for Circular traits (integer values), *Q* for
    Quantitative traits (numeric values) that is allowed **only** if
    there are at least 2 species with the same value, and *F* for
    fuzzy-coded traits (i.e. described with several 'sub-traits');

  - **fuzzy_name** name of fuzzy-coded trait to which 'sub-trait'
    belongs (if trait is not fuzzy, ignored so could be trait name or
    NA).

- fe_nm_type:

  a character string referring to the type of naming functional
  entities. Two possible values: *"fe_rank"* (FE are named after their
  decreasing rank in term of number of species *i.e.* fe_1 is the one
  gathering most species) and *"tr_val"* (FE are named after names of
  traits and of trait values for each FE, if possible, *see details
  below*). Default: `fe_nm_type = "fe_rank"`.

- check_input:

  a logical value indicating whether key features the inputs are checked
  (e.g. class and/or mode of objects, names of rows and/or columns,
  missing values). If an error is detected, a detailed message is
  returned. Default: `check.input = TRUE`.

## Value

A list of objects containing:

- **fe_nm**: a vector with names of all FE (following fe_nm_type). FE
  are ordered according to the decreasing number of species they gather.

- **sp_fe**: a vector containing for each species the name of the FE it
  belongs to. FE order is done according to decreasing number of
  species.

- **fe_tr**: a data frame containing traits values (variables in
  columns) for each FE (rows). FE order is done according to decreasing
  number of species.

- **fe_nb_sp**: a vector with species number per FE. If all FE have only
  one species, a warning message is returned. FE are ordered according
  to the decreasing number of species they gather.

- **details_fe**: a list containing: *fe_codes* a vector containing
  character referring to traits values (like a barcode) with names as in
  `fe_nm_type` and sorted according to `fe_nb_sp`; *tr_uval* a list
  containing for each trait a vector of its unique values or a data
  frame for fuzzy-coded traits; *fuzzy_E* a list with for each
  fuzzy-coded trait a data frame with names of entities (E) and names of
  species (sp); *tr_nb_uval* a vector with number of unique values per
  trait (or combinations for fuzzy-coded traits); *max_nb_fe* the
  maximum number of FE possible given number of unique values per trait.

## Details

`fe_nm_type = "tr_val"` is allowed **only** if:

- there are less than 7 traits;

- none of them is fuzzy-coded (so that names are not too long)

- all trait names and all trait values have different 2 first letters

If these 3 conditions are met, names of Functional Entities are made as
a character string of up to 2 letters for trait name in upper case font
then up to 2 letters for trait value in lower case font, separated by
"\_" between traits. Trait names are abbreviated to a single letter
whenever possible. *Examples:* ("TAc2_TBxx_TCyy", "TAc3_TBff_TCyy") or
("A2_Bx_Cy", "A3_Bf_Cy")

## Author

Sebastien Villeger, Nicolas Loiseau, and Camille Magneville

## Examples

``` r
# Load species traits data:
 data("fruits_traits", package = "mFD")

# Transform species traits data:
# Only keep the first 4 traits to illustrate FEs:
 fruits_traits <- fruits_traits[ , c(1:4)]   

# Load trait types data:
 data("fruits_traits_cat", package = "mFD")

# Transform the trait types data to only keep traits 1 - 4:
 fruits_traits_cat <- fruits_traits_cat[c(1:4), ]

# Gather species into FEs:
## gathering species into FEs (FEs named according to the decreasing...
## ...  number of species they gather):
 sp_FEs <- mFD::sp.to.fe(
      sp_tr      = fruits_traits, 
      tr_cat     = fruits_traits_cat, 
      fe_nm_type = "fe_rank")

## display FEs names:
sp_FEs$fe_nm
#>  [1] "fe_1"  "fe_2"  "fe_3"  "fe_4"  "fe_5"  "fe_6"  "fe_7"  "fe_8"  "fe_9" 
#> [10] "fe_10" "fe_11" "fe_12" "fe_13" "fe_14" "fe_15" "fe_16" "fe_17" "fe_18"
#> [19] "fe_19" "fe_20"

## display for each species the name of the FE it belongs to:
sp_FEs$sp_fe
#>         apple       apricot        banana       currant    blackberry 
#>        "fe_1"        "fe_2"        "fe_6"        "fe_7"        "fe_3" 
#>     blueberry        cherry         grape    grapefruit     kiwifruit 
#>        "fe_8"        "fe_9"       "fe_10"       "fe_11"       "fe_12" 
#>         lemon          lime        litchi         mango         melon 
#>        "fe_4"       "fe_13"       "fe_14"       "fe_15"        "fe_5" 
#>        orange passion_fruit         peach          pear     pineapple 
#>        "fe_4"       "fe_16"       "fe_17"        "fe_1"       "fe_18" 
#>          plum     raspberry    strawberry     tangerine   water_melon 
#>        "fe_2"        "fe_3"       "fe_19"       "fe_20"        "fe_5" 

## display trait values for each FE:
sp_FEs$fe_tr
#>          Size Plant     Climate Seed
#> fe_1   5-10cm  tree   temperate  pip
#> fe_2    3-5cm  tree   temperate  pit
#> fe_3    1-3cm shrub   temperate  pip
#> fe_4   5-10cm  tree subtropical  pip
#> fe_5  10-20cm  forb   temperate  pip
#> fe_6  10-20cm  tree    tropical none
#> fe_7    0-1cm shrub   temperate  pip
#> fe_8    0-1cm  forb   temperate  pip
#> fe_9    1-3cm  tree   temperate  pit
#> fe_10   1-3cm  vine   temperate  pip
#> fe_11 10-20cm  tree subtropical  pip
#> fe_12  5-10cm  vine   temperate  pip
#> fe_13   3-5cm  tree    tropical  pip
#> fe_14   1-3cm  tree    tropical  pit
#> fe_15 10-20cm  tree    tropical  pit
#> fe_16   3-5cm  vine    tropical  pip
#> fe_17  5-10cm  tree   temperate  pit
#> fe_18 10-20cm  forb    tropical none
#> fe_19   1-3cm  forb   temperate  pip
#> fe_20   3-5cm  tree subtropical  pip
 
## display the number of species per FEs:
sp_FEs$fe_nb_sp
#>  fe_1  fe_2  fe_3  fe_4  fe_5  fe_6  fe_7  fe_8  fe_9 fe_10 fe_11 fe_12 fe_13 
#>     2     2     2     2     2     1     1     1     1     1     1     1     1 
#> fe_14 fe_15 fe_16 fe_17 fe_18 fe_19 fe_20 
#>     1     1     1     1     1     1     1 
```
