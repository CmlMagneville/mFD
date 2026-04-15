# Package index

## Summarize data

These two functions are used to summarize up the species-traits and
assemblage data.

- [`sp.tr.summary()`](https://cmlmagneville.github.io/mFD/reference/sp.tr.summary.md)
  : Summarize Species x Traits data frame
- [`asb.sp.summary()`](https://cmlmagneville.github.io/mFD/reference/asb.sp.summary.md)
  : Summarize Assemblage x Species data frame

## Functional distances and Functional spaces

These functions compute the functional distances between species, build
functional spaces (if the data gather only continuous traits or not),
assess/plot the quality of functional spaces (and dendrogram if asked),
plot the chosen functional space and caracterize its functional axes
based on traits correlation with each functional axes.

- [`funct.dist()`](https://cmlmagneville.github.io/mFD/reference/funct.dist.md)
  : Compute functional distance between species
- [`tr.cont.scale()`](https://cmlmagneville.github.io/mFD/reference/tr.cont.scale.md)
  : Scale continuous traits
- [`tr.cont.fspace()`](https://cmlmagneville.github.io/mFD/reference/tr.cont.fspace.md)
  : Build a functional space based on continuous traits only
- [`quality.fspaces()`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.md)
  : Compute functional spaces and their quality
- [`quality.fspaces.plot()`](https://cmlmagneville.github.io/mFD/reference/quality.fspaces.plot.md)
  : Plot functional space quality with a chosen quality metric
- [`funct.space.plot()`](https://cmlmagneville.github.io/mFD/reference/funct.space.plot.md)
  : Plot species position in a functional space
- [`traits.faxes.cor()`](https://cmlmagneville.github.io/mFD/reference/traits.faxes.cor.md)
  : Correlation between Traits and Axes

## Compute Functional Entities

This function gathers species into Functional Entities (FEs)

- [`sp.to.fe()`](https://cmlmagneville.github.io/mFD/reference/sp.to.fe.md)
  : Compute Functional Entities composition based on a Species x Traits
  matrix

## Compute and Plot functional indices

### Based on FEs

These functions compute and plot functional indices based on FEs as in
`Mouillot et al. (2014)`. They compute/plot Functional Redundancy,
Functional Overredundancy and Functional Vulnerability. They also help
to work with FEs.

- [`alpha.fd.fe()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.fe.md)
  : Compute the set of indices based on number of species in Functional
  Entities
- [`alpha.fd.fe.plot()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.fe.plot.md)
  : Illustrate Functional Diversity indices based on Functional Entities
- [`fe.sp.df.computation()`](https://cmlmagneville.github.io/mFD/reference/fe.sp.df.computation.md)
  : Get a data frame linking Functional Entities names and species names
- [`from.fecoord.to.spcoord()`](https://cmlmagneville.github.io/mFD/reference/from.fecoord.to.spcoord.md)
  : Convert the data frame of FEs coordinates to a species coordinates
  one
- [`from.spfe.to.feasb()`](https://cmlmagneville.github.io/mFD/reference/from.spfe.to.feasb.md)
  : Build the assemblage-FEs dataframe from the assemblages-species one
- [`search.sp.nm()`](https://cmlmagneville.github.io/mFD/reference/search.sp.nm.md)
  : Get the names of species belonging to a specific Functional Entity
  (FE)

### Based on Hill numbers

These functions compute alpha and beta indices based on Hill numbers
according to the Chao et al. (2019) framework.

- [`alpha.fd.hill()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.hill.md)
  : Compute Functional alpha-Diversity indices based on Hill Numbers
- [`beta.fd.hill()`](https://cmlmagneville.github.io/mFD/reference/beta.fd.hill.md)
  : Compute Functional beta-Diversity indices based on Hill Numbers

### FUSE index

This function computes FUSE (Functionally Unique, Specialized, and
Endangered) index that combines functional uniqueness, specialisation
and global endangerment to identify threatened species of particular
importance for functional diversity based on `Pimiento et al. (2020)`.

- [`fuse()`](https://cmlmagneville.github.io/mFD/reference/fuse.md) :
  Compute FUSE (Functionally Unique, Specialized and Endangered)

### Based on multidimensional space

These functions compute/plot alpha and beta indices based on a given
multidimensional functional space.

- [`alpha.fd.multidim()`](https://cmlmagneville.github.io/mFD/reference/alpha.fd.multidim.md)
  : Compute a set of alpha functional indices for a set of assemblages
- [`beta.fd.multidim()`](https://cmlmagneville.github.io/mFD/reference/beta.fd.multidim.md)
  : Compute Functional beta-Diversity indices for pairs of assemblages
  in a multidimensional space
- [`alpha.multidim.plot()`](https://cmlmagneville.github.io/mFD/reference/alpha.multidim.plot.md)
  : Plot functional space and chosen functional indices
- [`beta.multidim.plot()`](https://cmlmagneville.github.io/mFD/reference/beta.multidim.plot.md)
  : Illustrate Functional beta-Diversity indices for pairs of
  assemblages in a multidimensional space

### Based on multidimensinal space for more complex graphs

These functions return `ggplot` layers for each index allowing users to
draw more complex graphs.

- [`background.plot()`](https://cmlmagneville.github.io/mFD/reference/background.plot.md)
  : Plot background of multidimensional plots
- [`fdiv.plot()`](https://cmlmagneville.github.io/mFD/reference/fdiv.plot.md)
  : Plot FDiv indice
- [`fdis.plot()`](https://cmlmagneville.github.io/mFD/reference/fdis.plot.md)
  : Plot FDis index
- [`feve.plot()`](https://cmlmagneville.github.io/mFD/reference/feve.plot.md)
  : Plot FEve index
- [`fide.plot()`](https://cmlmagneville.github.io/mFD/reference/fide.plot.md)
  : Plot FIde index
- [`fnnd.plot()`](https://cmlmagneville.github.io/mFD/reference/fnnd.plot.md)
  : Plot FNND index
- [`fori.plot()`](https://cmlmagneville.github.io/mFD/reference/fori.plot.md)
  : Plot FOri
- [`fric.plot()`](https://cmlmagneville.github.io/mFD/reference/fric.plot.md)
  : Plot FRic index
- [`fspe.plot()`](https://cmlmagneville.github.io/mFD/reference/fspe.plot.md)
  : Plot FSpe
- [`panels.to.patchwork()`](https://cmlmagneville.github.io/mFD/reference/panels.to.patchwork.md)
  : Plot individual plots along a pair of functional axes into a unique
  graph
- [`pool.plot()`](https://cmlmagneville.github.io/mFD/reference/pool.plot.md)
  : Plot species from the pool

## Other functions

Various functions that can be used by the user for diverse usage

- [`dist.nearneighb()`](https://cmlmagneville.github.io/mFD/reference/dist.nearneighb.md)
  : Compute distance of a given point to its nearest neighbor in the
  functional space and the identity of the nearest neighbor
- [`dist.point()`](https://cmlmagneville.github.io/mFD/reference/dist.point.md)
  : Compute distances of all points to a given point in the functional
  space
- [`dist.to.df()`](https://cmlmagneville.github.io/mFD/reference/dist.to.df.md)
  : Merge distance object(s) into a single data frame
- [`mst.computation()`](https://cmlmagneville.github.io/mFD/reference/mst.computation.md)
  : Compute the Minimum Spanning Tree (MST) linking species of a given
  assemblage
- [`sp.filter()`](https://cmlmagneville.github.io/mFD/reference/sp.filter.md)
  : Retrieve information about species in a given assemblage
- [`vertices()`](https://cmlmagneville.github.io/mFD/reference/vertices.md)
  : Compute vertices of the Minimal Convex Hull shaping species from a
  single assemblage in a multidimensional functional space

## Data sets

The three data sets used for examples and tutorials in the `mFD`
package.

- [`baskets_fruits_weights`](https://cmlmagneville.github.io/mFD/reference/baskets_fruits_weights.md)
  : Dataset: Baskets Composition in Fruits Species
- [`fruits_traits`](https://cmlmagneville.github.io/mFD/reference/fruits_traits.md)
  : Dataset: Traits Values of Fruits Species
- [`fruits_traits_cat`](https://cmlmagneville.github.io/mFD/reference/fruits_traits_cat.md)
  : Dataset: Fruits Traits Informations
