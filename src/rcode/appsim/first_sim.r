#' # UK school closure scenarios for COVID-19 
#' 
#' ## Load functions and population data
#' 
#' First clear memory. 
rm(list=ls(all=TRUE))

#' Load functions to run the hybrid model from local packages, from github or directly from local copies
#' if that is how the code has been distributed. 
library(devtools)

## install_github("c97sr/idd")
## uninstall("~/git/idd")
install("~/gdrive/git/idd",dependencies=FALSE)
library("idd")

## Next load the network into r
fnNetwork <- "~/"
