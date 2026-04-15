# Merge distance object(s) into a single data frame

This function merges distance object(s) into a single data frame which
rows are pairs of elements and column(s) distance metric(s). It stands
on the
[`dist_long`](https://talgalili.github.io/dendextend/reference/dist_long.html)
function.

## Usage

``` r
dist.to.df(list_dist)
```

## Arguments

- list_dist:

  a list of dist object(s). All dist objects should have a name (e.g.
  name of distance metric) and the same labels (i.e. names of sets
  between which distance was computed).

## Value

A data frame which first and second columns (names `x1` and `x2`)
contain names of the 2 sets involved in each pair, and with one column
for each dist object (named after its name in `list_dist`.

## Author

Sebastien Villeger

## Examples

``` r
# Create dist objects: 
dist_A <- round(dist(matrix(runif(10, 0, 1), 5, 2, 
                      dimnames = list(letters[1:5], NULL))), 2)
dist_B <- round(dist(matrix(runif(10, 0, 1), 5, 2, 
                      dimnames = list(letters[1:5], NULL))), 2)
dist_C <- round(dist(matrix(runif(10, 0, 1), 5, 2, 
                      dimnames = list(letters[1:5], NULL))), 2)

# First example with only 1 distance:
dist.to.df(list(dA = dist_A))
#>    x1 x2   dA
#> 1   a  b 0.36
#> 2   a  c 0.70
#> 3   a  d 0.63
#> 4   a  e 0.58
#> 5   b  c 0.34
#> 6   b  d 0.78
#> 7   b  e 0.22
#> 8   c  d 1.03
#> 9   c  e 0.14
#> 10  d  e 0.95

# Second example with 3 distances:
dist.to.df(list(d1 = dist_A, d2 = dist_B, d3 = dist_C))
#>    x1 x2   d1   d2   d3
#> 1   a  b 0.36 0.78 0.15
#> 2   a  c 0.70 0.62 0.53
#> 3   a  d 0.63 0.60 0.47
#> 4   a  e 0.58 0.09 0.74
#> 5   b  c 0.34 0.43 0.39
#> 6   b  d 0.78 0.24 0.33
#> 7   b  e 0.22 0.81 0.59
#> 8   c  d 1.03 0.20 0.07
#> 9   c  e 0.14 0.69 0.34
#> 10  d  e 0.95 0.66 0.33
```
