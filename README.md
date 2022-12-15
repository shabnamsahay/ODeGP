
# ODeGP

<!-- badges: start -->
<!-- badges: end -->

**O**scillation **De**tection with **G**aussian **P**rocesses

## Installation

If not installed, R can be downloaded from the [CRAN](https://cran.r-project.org/) homepage. RStudio is available for download [here](https://posit.co/products/open-source/rstudio/). A detailed guide to installation of both R and RStudio can be found at [this](https://rstudio-education.github.io/hopr/starting.html) webpage.

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

fname_file <- './dataset.csv' 
threshold <- 10
nonstat <- TRUE

d <- extractData(fname_file, errorbars = FALSE)
bayesF <- oscOrNot(d, nonstat, threshold, plotting=TRUE)

```

On completion, the calculated Bayes factor for the data will be printed to the console, along with its comparison to the given threshold to determine whether the data is oscillatory or not. The non-stationary kernel's posterior will be saved to the file <code>nonstationaryPosteriorData.csv</code> in the working directory.

More comprehensive running examples will be added shortly.


