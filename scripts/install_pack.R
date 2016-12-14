source("http://bioconductor.org/biocLite.R")
ipak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    biocLite(new.pkg)
sapply(pkg, require, character.only = TRUE)
}

packages <- c("S4Vectors", "IRanges", "BiocGenerics", "Biobase", "BiocParallel", 
"genefilter", "geneplotter", "vsn", "phyloseq", "Heatplus", "mmnet", "RCy3", "edgeR")

ipak(packages)

ipak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, repos='http://cran.us.r-project.org')
sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "vegan", "dplyr", "reshape2", "RColorBrewer", "rmarkdown", 
"knitr", "grid", "reshape", "data.table", "biom", "gplots", "gridExtra", "mvtnorm",
"pls", "QuantPsyc", "reshape2", "R.utils", "pheatmap", "optparse", "Hmisc", "locfit", 
"Rcpp", "httr", "RJSONIO", "RcppArmadillo", "ggthemes", "sitools")

ipak(packages)

devtools::install_github("gavinsimpson/ggvegan")