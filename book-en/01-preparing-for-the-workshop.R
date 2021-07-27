# Install the required packages
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("ggplot2")

# install mvpart from package archive file
install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Load the required packages
library(vegan)
library(labdsv)
library(MASS)
library(mvpart)
library(ggplot2)
