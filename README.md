## DeepPheWAS

## Overview

DeepPheWAS is an R package for running phenome wide association studies (PheWAS). 
It allows user control of all the stages of PheWAS, from data wrangling, through
phenotype generation and association testing. 

## Installation

``` r
# The development version from GitHub:
# install.packages("devtools")
devtools::install_github("Richard-Packer/DeepPheWAS")
```
## Using DeepPheWAS
DeepPheWAS is designed to be used via the inbuilt R scripts. These R-scripts  accessed via 
using a bash interface and utilise an argument parser (docopt) to process command line 
inputs (arguments).
The inbuilt R-scripts have been split into two folders to represent the two main stages of PheWAS, phenotype generation and association testing. The R scripts are located within the DeepPheWAS package at:

/extdata/scripts/phenotype_generation/

and:

/extdata/scripts/association_testing/

For further details of the functionality of the scripts and how to use them, see [user guide](https://www.medrxiv.org/content/10.1101/2022.05.05.22274419v1.supplementary-material)
