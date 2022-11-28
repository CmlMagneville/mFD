## Resubmit comments

* Correct follow moved content ulr
* Fixed the M1Mac issue when building vignettes
* Add new vignettes:
  * `vignettes/Customised_plots.Rmd`
* Review old vignettes:
  * fixed typos and missing chunks in the general workflow vignette
* Fix minor bugs in some functions:
  * correct beta functional axes names
  * correct check fuzzy names
* Review function documentation
* For a complete overview of the new version, please see `NEWS.md`
* Add paper citation to the READ.ME file

## Test environments

* Local
  * Windows 11, R 4.0.5
* Github Actions
  * macOS 11.6.5, R-release (R 4.2.2)
  * Windows Server 2022, R-release (R 4.2.2)
  * Ubuntu 18.04.6 LTS, R-devel, R-release (R 4.2.2), R-oldrel
* WinBuilder
  * r-devel
  * r-release
  * r-oldrel
* R-hub
  * Windows Server 2022, R-devel 64 bit
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran


## R CMD check results

0 error | 0 warning | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Camille Magneville <camille.magneville@gmail.com>'
  
  Version contains large components (1.0.2.9000)
  
  Found the following (possibly) invalid URLs:
    URL: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343
      From: inst/doc/Compute_functional_hill_indices.html
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/10.1111/ecog.05904
      From: README.md
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05904
      From: inst/CITATION
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12021
      From: inst/doc/mFD_general_workflow.html
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12299
      From: inst/doc/Compute_and_interpret_quality_of_functional_spaces.html
            inst/doc/Continuous_traits_framework.html
            inst/doc/mFD_general_workflow.html
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13185
      From: inst/doc/mFD_general_workflow.html
      Status: 503
      Message: Service Unavailable
    URL: https://www.pnas.org/doi/abs/10.1073/pnas.1317625111
      From: inst/doc/How_to_deal_with_Functional_Entities.html
            inst/doc/mFD_general_workflow.html
      Status: 503
      Message: Service Unavailable
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1002/ecm.1343
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503
    DOI: 10.1073/pnas.1317625111
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503
    DOI: 10.1111/geb.12299
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503



## Downstream dependencies

There are currently no downstream dependencies for this package.

