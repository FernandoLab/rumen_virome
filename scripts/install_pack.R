source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")

ipak <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, repos='http://cran.us.r-project.org')
sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "vegan", "reshape2", "RColorBrewer", "rmarkdown", 
"knitr", "grid", "reshape", "data.table", "biom", "gplots", "gridExtra", "mvtnorm")

ipak(packages)

devtools::install_github("gavinsimpson/ggvegan")
