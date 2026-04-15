# Compute Functional beta-Diversity indices for pairs of assemblages in a multidimensional space

Computes a set of indices of pairwise functional beta-diversity
(dissimilarity and its turnover and nestedness-resultant components)
based on overlap between convex hulls in a multidimensional space. For
details about indices formulas see Villeger *et al.* (2013). This
functions stands on
[`functional.betapart.core.pairwise`](https://rdrr.io/pkg/betapart/man/functional.betapart.core.pairwise.html).

## Usage

``` r
beta.fd.multidim(
  sp_faxes_coord,
  asb_sp_occ,
  check_input = TRUE,
  beta_family = "Jaccard",
  details_returned = TRUE,
  betapart_step = TRUE,
  betapart_chullopt = list(conv1 = "Qt", conv2 = "QJ"),
  betapart_para = FALSE,
  betapart_para_opt = betapart::beta.para.control()
)
```

## Arguments

- sp_faxes_coord:

  a matrix with coordinates of species (rows) on functional axes
  (columns). Species coordinates have been retrieved thanks to
  [`tr.cont.fspace`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  or
  [`quality.fspaces`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md).

- asb_sp_occ:

  a matrix with presence/absence (coded as 0/1) of species (columns) in
  a set of assemblages (rows).

- check_input:

  a logical value defining whether inputs are checked before computation
  of indices. Possible error messages will thus may be more
  understandable for the user than R error messages. Default:
  `check_input = TRUE`.

- beta_family:

  a character string for the type of beta-diversity index to use,
  'Jaccard' (default) and/or 'Sorensen'.

- details_returned:

  a logical value indicating whether the user wants to details_returned.
  Details are used in the graphical function `beta.multidim.plot` and
  thus must be kept if the user want to have graphical outputs for the
  computed indices.

- betapart_step:

  a logical value indicating whether the computation progress should be
  displayed in the R console. Default: `betapart_step = TRUE`.

- betapart_chullopt:

  a A list of two named vectors of character conv1 and conv2 defining
  the options that will be used to compute the convex hulls (through the
  options of geometry::convhulln function). For further details see help
  of
  [`functional.betapart.core.pairwise`](https://rdrr.io/pkg/betapart/man/functional.betapart.core.pairwise.html).
  Default: `betapart_chullopt = c(conv1 = 'Qt', conv2 = 'QJ')`.

- betapart_para:

  a logical value indicating whether internal parallelization should be
  used to compute pairwise dissimilarities. Default:
  `betapart_para = FALSE`.

- betapart_para_opt:

  a list with details about parallelization. Default value means those
  parameters are set according to computer specifications. `nc` is the
  number of cores (default = 4), `type` is a character string with code
  of method used (default PSOCK), `LB` is a boolean specifying whether
  load-balancing is applied (default is TRUE) and `size` is a numeric
  value for number of tasks performed at each time (default is 1). See
  help of
  [`functional.betapart.core.pairwise`](https://rdrr.io/pkg/betapart/man/functional.betapart.core.pairwise.html)
  for more details.

## Value

A list with:

- *pairasb_fbd_indices* a list with for each index a *dist* object with
  values for all pairs of assemblages. Indices names start with the
  abbreviation of the type of dissimilarity index ('jac' for
  Jaccard-like and 'sor' for Sorensen-like dissimilarity) and end with
  abbreviation of component ('diss': dissimilarity, 'turn' its turnover
  component, and 'nest' its nestedness-resultant component).

- if *store_details* is TRUE,

- *details_beta* list: **inputs** a list with *sp_faxes_coord* and
  *asb_sp_occ* on which indices were computed (required for drawing
  graphics), **pool_vertices** a list of vectors (1 per assemblage) with
  names of species being vertices of the convex hull shaping all
  species; **asb_FRic** a vector with volume of the convex hull shaping
  each assemblage (relative to volume filled by all species) ;
  **asb_vertices** a list of vectors (1 per assemblage) with names of
  species being vertices of the convex hull

## Note

All assemblages should have a number of species strictly higher than the
number of functional axes. Computing intersection of convex hulls in
space of \>5 dimensions is yet impossible with most computers. This
function uses R libraries 'betapart' (\> =1.5.4) for indices
computation. Indices values are stored as *dist* objects to optimize
memory. See below example of how merging distance values in a
*dataframe* with `dist.to.df`.

## References

Villeger *et al.* (2013) Decomposing functional beta-diversity reveals
that low functional beta-diversity is driven by low functional turnover
in European fish assemblages. *Global Ecology and Biogeography*, **22**,
671-681.

## Author

Sebastien Villeger and Camille Magneville

## Examples

``` r
if (FALSE) { # \dontrun{
 # Load Species*Traits dataframe:
data('fruits_traits', package = 'mFD')

# Load Assemblages*Species dataframe:      
data('baskets_fruits_weights', package = 'mFD') 

# Load Traits categories dataframe:
data('fruits_traits_cat', package = 'mFD') 
 
# Compute functional distance 
sp_dist_fruits <- mFD::funct.dist(sp_tr         = fruits_traits,
                                  tr_cat        = fruits_traits_cat,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fruits <- mFD::quality.fspaces(
 sp_dist             = sp_dist_fruits, 
 maxdim_pcoa         = 10,
 deviation_weighting = 'absolute',
 fdist_scaling       = FALSE,
 fdendro             = 'average')
 
# Retrieve species coordinates matrix:
sp_faxes_coord_fruits <- fspaces_quality_fruits$details_fspaces$sp_pc_coord

# Get the occurrence dataframe:
asb_sp_fruits_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_fruits_weights) 
asb_sp_fruits_occ <- asb_sp_fruits_summ$'asb_sp_occ'

# Compute beta diversity indices:
beta_fd_fruits <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fruits[, c('PC1', 'PC2', 'PC3', 'PC4')], 
  asb_sp_occ       = asb_sp_fruits_occ,
  check_input      = TRUE,
  beta_family      = c('Jaccard'),
  details_returned = TRUE)

# merging pairwise beta-diversity indices in a data.frame
dist.to.df(beta_fd_fruits$pairasb_fbd_indices)
} # }
```
