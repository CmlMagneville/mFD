# Dataset: Baskets Composition in Fruits Species

This dataset represents the abundance of 25 fruits species in 10
baskets.

## Usage

``` r
baskets_fruits_weights
```

## Format

A matrix of integers with 10 rows (baskets) and 25 columns (species).

## See also

`fruits_traits`, `fruits_traits_cat`

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Assemblages x Species Matrix
data("baskets_fruits_weights", package = "mFD")
baskets_fruits_weights[1:5, 1:5]

# Summarize Assemblages Data
mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights)
} # }
```
