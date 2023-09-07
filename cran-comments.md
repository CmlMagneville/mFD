## Resubmit comments

* Fix bug in some function:
  * correct funct.dist() function when only quantitative traits
* Remove FSGE and FUGE outputs from the fuse() function
* For a complete overview of the new version, please see `NEWS.md`

## Test environments

* Local
  * Windows 11, R 4.3.1
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

0 error | 0 warning | 2 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Camille Magneville <camille.magneville@gmail.com>'
  
  Version contains large components (1.0.4.9000)
  
  Found the following (possibly) invalid URLs:
    URL: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343
      From: inst/doc/Compute_functional_hill_indices.html
      Status: 403
      Message: Forbidden
    URL: https://onlinelibrary.wiley.com/doi/10.1111/ecog.05904
      From: README.md
      Status: 403
      Message: Forbidden
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05904
      From: inst/CITATION
      Status: 403
      Message: Forbidden
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12021
      From: inst/doc/mFD_general_workflow.html
      Status: 403
      Message: Forbidden
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12299
      From: inst/doc/Compute_and_interpret_quality_of_functional_spaces.html
            inst/doc/Continuous_traits_framework.html
            inst/doc/mFD_general_workflow.html
      Status: 403
      Message: Forbidden
    URL: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13185
      From: inst/doc/mFD_general_workflow.html
      Status: 403
      Message: Forbidden
    URL: https://www.pnas.org/doi/abs/10.1073/pnas.1317625111
      From: inst/doc/How_to_deal_with_Functional_Entities.html
            inst/doc/mFD_general_workflow.html
      Status: 403
      Message: Forbidden
  
* checking for future file timestamps ... NOTE
  unable to verify current time


## Downstream dependencies

There are currently no downstream dependencies for this package.

