# Installez les paquets requis
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("ggplot2")

# installez mvpart de l'archive
install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Chargez les pacquets requis
library(vegan)
library(labdsv)
library(MASS)
library(mvpart)
library(ggplot2)
