
# ODeGP

<!-- badges: start -->
<!-- badges: end -->

**O**scillation **De**tection with **G**aussian **P**rocesses

## Installation

If not installed, R can be downloaded from the [CRAN](https://cran.r-project.org/) homepage. RStudio is available for download [here](https://posit.co/products/open-source/rstudio/). A detailed guide to installation of both R and RStudio can be found on [this](https://rstudio-education.github.io/hopr/starting.html) webpage.

You can install the development version of ODeGP like so:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("shabnamsahay/ODeGP")
```

## Example

This is a basic example of how to run ODeGP. The sample file <code>dataset.csv</code> is provided in the top-level of this repository.

``` r
library("ODeGP")

filepath <- './dataset.csv'
threshold <- 14

d <- extractData(filepath, errorbars = FALSE)
bayesF <- oscOrNot(d, threshold, detrend = TRUE, plotting = TRUE)

```

On completion, the calculated Bayes factor for the data will be printed to the console, along with its comparison to the given `threshold` to determine whether the data is oscillatory or not. The non-stationary kernel's posterior will be saved to the file <code>nonstationaryPosteriorData.csv</code> in the working directory.

## Parameters

A short description of the remaining parameters:

- `errorbars`: set this to TRUE if replicates have been merged in the input data (resulting in errorbars already being present), else FALSE (i.e. distinct replicates are provided in the input data, so errorbars are not present)
- `detrend`: if you want to detrend your data via linear regression, set this to TRUE. If your data is already detrended OR you do not want to detrend the data at all, set this parameter to FALSE
- `plotting`: set this to TRUE if you want an output plot of the optimised non-stationary kernel posterior on the input data, else set to FALSE

An example dataset file that would require `errorbars = TRUE` to be set when `extractData()` is called is provided in the top-level of this repository as <code>dataset_errorbars.csv</code>. Column names used in any input file to ODeGP must be identical to those shown in the example files in either case. 

More comprehensive running examples will be added shortly.
