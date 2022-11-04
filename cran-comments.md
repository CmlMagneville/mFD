## Test environments

* local Windows 11, R 4.0.5
* Mac OS X 10.15 (on GitHub Actions), r-release (R 4.1.1)
* Windows Server 2019 (on GitHub Actions), r-release (R 4.1.1)
* Ubuntu 20.04 (on GitHub Actions), r-release (R 4.1.1)


## R CMD check results

0 error | 0 warning | 1 note

* checking CRAN incoming feasibility ... NOTE

  Maintainer: 'Camille Magneville <camille.magneville@gmail.com>'
  
  Found the following (possibly) invalid URLs:
    URL: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1343
      From: inst/doc/Compute_functional_hill_indices.html
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/10.1111/ecog.05904
      From: README.md
      Status: 503
      Message: Service Unavailable
    URL: https://onlinelibrary.wiley.com/doi/10.1111/ecog.05904?af=R
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
    URL: https://www.pnas.org/content/111/38/13757.short (moved to https://www.pnas.org/doi/abs/10.1073/pnas.1317625111)
      From: inst/doc/How_to_deal_with_Functional_Entities.html
      Status: 503
      Message: Service Unavailable
    URL: https://www.pnas.org/doi/abs/10.1073/pnas.1317625111
      From: inst/doc/mFD_general_workflow.html
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


## Resubmit comments

* References have been added to the DESCRIPTION file
* Packages have been written in single quotes in title and description
* Examples are now executed
* User messages with print() function have been changed using stop() function
