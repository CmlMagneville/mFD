# Dataset: Traits Values of Fruits Species

This dataset represents the value of 6 traits for 15 fruits species.
**Important:** Species names must be specified as the data frame row
names (not in an additional column).

## Usage

``` r
fruits_traits
```

## Format

A data frame with 25 rows (species) and the following columns (traits):

- Size:

  an ordered factor describing the size of fruits species

- Plant:

  an unordered factor describing the type of plant

- Climate:

  an ordered factor describing the climate regions

- Seed:

  an ordered factor describing the type of seed

- Sugar:

  a numeric describing the quantity of sugar

- Use.raw:

  an integer (percentage) describing the proportion of a raw use (fuzzy
  trait) of the fruit

- Use.pastry:

  an integer (percentage) describing the proportion of a pastry use
  (fuzzy trait) of the fruit

- Use.jam:

  an integer (percentage) describing the proportion of a jam use (fuzzy
  trait) of the fruit

## See also

`fruits_traits_cat`, `baskets_fruits_weights`

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Species x Traits Data Frame
data("fruits_traits", package = "mFD")
fruits_traits

# Load Traits Information
data("fruits_traits_cat", package = "mFD")
fruits_traits_cat

# Summarize Species x Traits Data
mFD::sp.tr.summary(tr_cat = fruits_traits_cat, sp_tr = fruits_traits)
} # }
```
