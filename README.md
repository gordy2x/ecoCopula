## ecoCopula [![Build Status](https://travis-ci.com/gordy2x/ecoCopula.svg)](https://travis-ci.com/gordy2x/ecoCopula) [![License](http://img.shields.io/badge/license-LGPL%20%28%3E=%202.1%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![DOI](https://zenodo.org/badge/139233335.svg)](https://zenodo.org/badge/latestdoi/139233335)

=======

R package to find direct and indirect species associations from co-occurrence data

### Installation

#install.packages("devtools")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
devtools::install_github("gordy2x/ecoCopula",upgrade = "always")

#### Test it is working

library(ecoCopula)
example(cord)
example(cgr)

If you have a problem installing, please email me (g.popovic@unsw.edu.au). 

### Author

Gordana Popovic 

### Contributor

Michelle Lim

### License

LGPL (>= 2.1)

