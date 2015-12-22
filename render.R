args <- commandArgs(trailingOnly = TRUE)
library(knitr)
library(rmarkdown)
render(args[1], "all")
