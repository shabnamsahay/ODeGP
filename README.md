
# ODeGP

<!-- badges: start -->
<!-- badges: end -->

**O**scillation **De**tection with **G**aussian **P**rocesses

## Installation

You can install the development version of ODeGP like so:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("shabnamsahay/ODeGP")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("ODeGP")

fname_file <- './dataset.csv' 
threshold <- 10
nonstat <- TRUE

d <- extractData(fname_file, errorbars = FALSE)
bayesF <- oscOrNot(d, nonstat, threshold, plotting=TRUE)

```

The sample <code>dataset.csv</code> is provided in the top-level of this repository.

