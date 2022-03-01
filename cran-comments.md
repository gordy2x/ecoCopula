## Test environments
* local R installation, R 4.1.2
* ubuntu 16.04 (on travis-ci), R 4.1.2
* win-builder (devel)
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1111/2041-210X.13247
    From: inst/doc/the_basics.html
    Status: 503
    Message: Service Unavailable

Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.13247
    From: inst/CITATION
    Status: Service Unavailable
    Message: 503

The URL and DOI are valid and work in browser.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
