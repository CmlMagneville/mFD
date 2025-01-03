# mFD (development version)
* Include new functions to deal with FEs:
* from.spfe.to.feasb() which computes the assemblages × FEs dataframe with outputs of the mFD::sp.to.fe() function
* fe.sp.df.computation() which creates a dataframe linking FEs names to species names with outputs of the mFD::sp.to.fe() function
* search.sp.nm() which finds a species name given an FE name
* from.fecoord.to.spcoord() function which converts the dataframe of FEs coordinates to one with species coordinates.
* Update the website with a FAQ based on user's questions 


# mFD 1.0.7
* Traits names in the tr_cat dataframe must be the same that in the sp_tr
dataframe: otherwise can lead to bugs if using weights: check.sp.tr()
* Add a sentence explaining that in the General Workflow vignette

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





