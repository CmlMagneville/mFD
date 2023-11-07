# mFD (development version)

# mFD 1.0.6
* Fix bug:
  * Correct the sp.plot() fct when plotting one asb (pb with vertices aes)
  * Correct the sp.plot() fct when several sp same place but only one as vert

# mFD 1.0.5
* Fix bug:
  * Correct funct.dist() function error when only quantitative traits (ISSUE 34)
* Remove FSGE and FUGE from fuse() outputs

# mFD 1.0.4
* Add warning messages in `funct.dist` function: when using only continuous 
traits, no weighting is realised for now.
* Add informations in the General Tutorial and in the Continuous Traits,
no weighting is realised for now.
* Fix bug:
  * when computing FRic and FDiv but not enough species at different
coordinates in the functional space to compute the convex-hull. Before it 
stopped with the error "Error in fdiv.computation(): Names of the vertices are 
not all present in species coordinates matrix. Please check." as FRic vert
were NA. Now FRic vert are NULL and pas the check.
  * when computing FEs with fuzzy traits, attribution of species to wrong FEs,
so when computing FEs, creating as many FEs as species.

# mFD 1.0.3

* Fix bug: 
  * scale center traits (parenthesis missing) in `tr.cont.scale` function
  * correct FDiv legend in plot
  * correct typo in Functional entities vignette

# mFD 1.0.2

* New vignette: create customised plots
* Review mFD general workflow: fixed typos and missing chunks
* Fix minor bugs in some functions:
  * correct beta functional axes names, `beta.multidim.plot()` function
  * correct check fuzzy names, `check.fuzzy()` function
* Review function documentation
* Add paper citation to the READ.ME file

# mFD 1.0.1

* Fix bug with fdist_scaling = TRUE in quality.fspaces function (Issue #19)


# mFD 1.0.0

* First release of the package.
* Submission to CRAN





