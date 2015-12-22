source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")

ipak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "vegan", "reshape2", "RColorBrewer", "rmarkdown", 
"knitr", "grid", "reshape", "data.table", "devtools", "ggvegan", "biom", "RCurl", 
"gplots", "gridExtra", "mvtnorm")

ipak(packages)
