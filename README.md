## ecoCopula 
<!-- badges: start -->
  [![R-CMD-check](https://github.com/gordy2x/ecoCopula/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gordy2x/ecoCopula/actions/workflows/R-CMD-check.yaml)
[![License](http://img.shields.io/badge/license-LGPL%20%28%3E=%202.1%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![DOI](https://zenodo.org/badge/139233335.svg)](https://zenodo.org/badge/latestdoi/139233335)
[![Codecov test coverage](https://codecov.io/gh/gordy2x/ecoCopula/branch/master/graph/badge.svg)](https://app.codecov.io/gh/gordy2x/ecoCopula?branch=master)
<!-- badges: end -->


R package to find direct and indirect species associations from co-occurrence data

### Installation

To install `ecoCopula` from [CRAN](https://CRAN.R-project.org/package=ecoCopula):
```r
install.packages("ecoCopula")
```

For the development version with the latest bells and whistles:
```r
# install.packages("devtools")
devtools::install_github("gordy2x/ecoCopula")
```

For the development version with zero-inflated functionality:
```r
# install.packages("devtools")
# devtools::install_github("r-forge/countreg/pkg")
devtools::install_github("gordy2x/ecoCopula", ref="e401671")
```

If you have trouble installing, please email me at [g.popovic@unsw.edu.au](mailto:g.popovic@unsw.edu.au)
