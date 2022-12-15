
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

This is a basic example of how to run ODeGP. The sample <code>dataset.csv</code> is provided in the top-level of this repository.

``` r
library("ODeGP")

fname_file <- './dataset.csv' 
threshold <- 10
nonstat <- TRUE

d <- extractData(fname_file, errorbars = FALSE)
bayesF <- oscOrNot(d, nonstat, threshold, plotting=TRUE)

```

On completion, the calculated Bayes factor for the data will be printed to the console, along with its comparison to the given threshold to determine whether the data is oscillatory or not. The non-stationary kernel's posterior will be saved to the file nonstationaryPosteriorData.csv in the working directory.


