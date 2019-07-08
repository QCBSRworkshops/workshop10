# Name: Amanda Winegardner, Xavier Giroux-Bougard, Bérenger Bourgeois, Emmanuelle Chrétien and Monica Granados 
# Date: March 2016
# Description: Multivariate analyses II
# Dataset: File names are "DoubsSpe.csv", "DoubsEnv.csv"
# Notes: Material in R script obtained from 
#        Borcard D., Gillet F. et Legendre P., 2011. Numerical Ecology with R. Springer.
#***************************************************************************************#


# 0. Loading the data ####

rm(list=ls())
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("plyr")
install.packages("rJava")
install.packages("venneuler")

# For the two following packages, upload the file provided on the wiki page.
# To do so, go to Packages tab on the bottom right panel of R Studio
# Click on Install Packages
# Choose to install from Package Archive file and upload these files
install.packages("mvpart")
install.packages("MVPARTwrap") 
install.packages("rdaTest") 

# If it is not working, you might want to try to load them from github.
# To do so, you need to install the package "devtools" and to perform these few steps
install.packages("devtools")
library(devtools)

# Install mvpart, MVPARTwrap and rdaTest from github
devtools::install_github("cran/mvpart")
devtools::install_github("cran/MVPARTwrap")
devtools::install_github("philippec/fonctions_R_git/rdaTest")

# Load packages
library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(MASS)
library(plyr)

# Species community data frame (fish abundance): "DoubsSpe.csv"
spe <- read.csv(file.choose(), row.names=1)
spe <- spe[-8,] # site number 8 contains no species and must be removed

# Environmental data frame: "DoubsEnv.csv"
env <- read.csv(file.choose(), row.names=1)
env <- env[-8,]

#***************************************************************************************#

# 1.Explore the data ####

## 1.1 Species data ####
names(spe)
dim(spe)
str(spe)
head(spe)
summary(spe) 

### Species distibution, all species confounded
(ab <- table(unlist(spe)))
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

### Number of absences
sum(spe==0)

### Proportion of zeros in the community data set
sum(spe==0)/(nrow(spe)*ncol(spe))

## Apply an appropriate transformation for community composition data

?decostand # this function provide some standardiztion methods for community domposition data (see method options)

# Hellinger transformation
spe.hel<-decostand(spe, method="hellinger")

## 1.2 Environmental data ####
names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Bivariate Plots of the Environmental Data" )

# Standardization of the 11 environmental data
env.z <- decostand(env, method="standardize")
apply(env.z, 2, mean) # the data are now centered (means~0)
apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)

# Note: explanatory variables expressed in different units
#       should be standardized before computing distances measures and ordination analysis. 
#       It is however not necessary for multivariate tree analysis.


#***************************************************************************************#
# 2. Canonical Analysis
## 2.1. Redundancy Analysis (RDA) ####

### Prepare the data
env.z <- subset(env.z, select = -das) # remove the "distance from the source" variable

### Run the RDA of all exlanatory variables of env on species abundances
?rda
spe.rda <- rda(spe.hel~., data=env.z) # the dot "." means that all variables of env.z will be included.

### Extract the results
summary(spe.rda, display=NULL)

### Select the significant explanatory variables by forward selection
?ordiR2step
ordiR2step(rda(spe.hel~1, data=env.z), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
env.signif <- subset(env.z, select = c("alt", "oxy", "dbo"))

### Re-run the RDA of significant explanatory variables on species abundances
spe.rda.signif <- rda(spe.hel~., data=env.signif)

### Extract the results
summary(spe.rda.signif, display=NULL)

### Calculate the adjusted R^2 of the RDA
(R2adj <- RsquareAdj(spe.rda.signif)$adj.r.squared)

### Test the significance of the model and the significance of axis
?anova.cca
anova(spe.rda.signif, step=1000)
anova(spe.rda.signif, step=1000, by="axis")

### Plot the results
#### quick plots scaling 1 and 2
windows()
plot(spe.rda.signif, scaling=1, main="Triplot RDA (scaling 1)")
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA (scaling 2)")

#### advanced plots scaling 1
windows()
plot(spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spe.rda.signif, display="species", choices=c(1), scaling=1),
       scores(spe.rda.signif, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(spe.rda.signif, display="species", choices=c(1), scaling=1),
     scores(spe.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spe.rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spe.rda.signif, display="bp", choices=c(1), scaling=1),
       scores(spe.rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(spe.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(spe.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(spe.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#### advanced plots scaling 2
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2,
       scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2.1,
     scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(spe.rda.signif, display="species", scaling=2)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spe.rda.signif, display="bp", choices=c(1), scaling=2),
       scores(spe.rda.signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(spe.rda.signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# Challenge 1 ####
# Run the RDA of the mite environmental variables of the mite species abundances
data(mite) #data from vegan package
mite.spe <- mite
data(mite.env) 
# Which are the significant explanatory variables ?  
# How much of the variation explain the significant explanatory variables ?
# Which are the significant axis ?
# Which group of sites can you identify ? 
# Which species are related to each group of sites ?
# (answer below)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
mite.spe.hel<-decostand(mite.spe,method="hell")
mite.spe.rda <- rda(mite.spe.hel~., data=mite.env)
ordiR2step(rda(mite.spe.hel~1, data=mite.env), 
           scope= formula(mite.spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
which(colnames(mite.env)=="WatrCont")
which(colnames(mite.env)=="Shrub")
which(colnames(mite.env)=="Substrate")
which(colnames(mite.env)=="Topo")
which(colnames(mite.env)=="SubsDens")
(mite.env.signif=mite.env[,c(2,4,3,5, 1)])
mite.env.signif <- subset(mite.env, 
                          select = c("WatrCont", "Shrub", "Substrate", "Topo", "SubsDens"))
mite.spe.rda.signif <- rda(mite.spe~., data=mite.env.signif) 
summary(mite.spe.rda.signif, display=NULL)
(R2adj <- RsquareAdj(mite.spe.rda.signif)$adj.r.squared)
anova(mite.spe.rda.signif, step=1000)
anova(mite.spe.rda.signif, step=1000, by="axis")          
windows()
plot(mite.spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.rda.signif, display="species", choices=c(1), scaling=1),
     scores(mite.spe.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1),
       scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

## 2.2  Partial RDA ####
## effect of water chemistry on species abundances partialling out the physiography 
### Divide the env data frame in two data frames
envtopo <- env.z[, c(1:3)] # Physiography : explanatory dataset 1
names(envtopo)
envchem <- env.z[, c(4:10)] # Water quality : explanatory dataset 2
names(envchem)
### Run the partial RDA
spechem.physio <- rda(spe.hel, envchem, envtopo)
# or : 
spechem.physio2 <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
                         Condition(alt + pen + deb), data=env.z)
### Extract the results
summary(spechem.physio, display=NULL)
summary(spechem.physio2, display=NULL)
### Calculate adjusted the R^2 of the partial RDA
(R2adj <- RsquareAdj(spechem.physio)$adj.r.squared)
### Test of the partial RDA and the significance of axis
anova(spechem.physio, step=1000)
anova(spechem.physio2, step=1000, by="axis")
### Construct the triplots
#### Scaling 1
windows(title="Partial RDA scaling 1")
plot(spechem.physio, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spechem.physio, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spechem.physio, display="species", choices=c(1), scaling=1),
       scores(spechem.physio, display="species", choices=c(2), scaling=1),
       col="grey",length=0)
text(scores(spechem.physio, display="species", choices=c(1), scaling=1),
     scores(spechem.physio, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spechem.physio, display="species", scaling=1)),
     col="grey", cex=0.8)  
arrows(0,0,
       scores(spechem.physio, display="bp", choices=c(1), scaling=1),
       scores(spechem.physio, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(spechem.physio, display="bp", choices=c(1), scaling=1)+0.05,
     scores(spechem.physio, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(spechem.physio, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 
#### Scaling 2
windows(title="Partial RDA scaling 2")
plot(spechem.physio, scaling=2, main="Triplot partial RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spechem.physio, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spechem.physio, display="species", choices=c(1), scaling=2),
       scores(spechem.physio, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(spechem.physio, display="species", choices=c(1), scaling=2),
     scores(spechem.physio, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(spechem.physio, display="species", scaling=2)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spechem.physio, display="bp", choices=c(1), scaling=2),
       scores(spechem.physio, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(spechem.physio, display="bp", choices=c(1), scaling=2)+0.05,
     scores(spechem.physio, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(spechem.physio, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# Challenge 2 ####
# Run the partial RDA of the mite environmental variables of the mite species abundances 
# partialling out for the substrate variables (SubsDens, WaterCont and Substrate)
# Is the model significant ?
# Which are the significant axis ?
# Interpret the obtained triplot.
# (answer below)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
mite.spe.subs<-rda(mite.spe.hel ~ Shrub + Topo
                  + Condition(SubsDens + WatrCont + Substrate), data=mite.env)
summary(mite.spe.subs, display=NULL)
(R2adj <- RsquareAdj(mite.spe.subs)$adj.r.squared)
anova(mite.spe.subs, step=1000)
anova(mite.spe.subs, step=1000, by="axis")
#### Triplot scaling 1
windows(title="Partial RDA scaling 1")
plot(mite.spe.subs, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(mite.spe.subs, display="species", choices=c(1), scaling=1),
       scores(mite.spe.subs, display="species", choices=c(2), scaling=1),
       col="grey",length=0)
text(scores(mite.spe.subs, display="species", choices=c(1), scaling=1),
     scores(mite.spe.subs, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.subs, display="species", scaling=1)),
     col="grey", cex=0.8)  
arrows(0,0,
       scores(mite.spe.subs, display="bp", choices=c(1), scaling=1),
       scores(mite.spe.subs, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(mite.spe.subs, display="bp", choices=c(1), scaling=1)+0.05,
     scores(mite.spe.subs, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(mite.spe.subs, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 
#### Triplot scaling 2
windows(title="Partial RDA scaling 2")
plot(mite.spe.subs, scaling=2, main="Triplot partial RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
       col="grey",length=0)
text(scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
     scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(mite.spe.subs, display="species", scaling=2)),
     col="grey", cex=0.8)  
arrows(0,0,
       scores(mite.spe.subs, display="bp", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(mite.spe.subs, display="bp", choices=c(1), scaling=2)+0.05,
     scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# 2.3. Variation partitioning by partial RDA ####
?varpart
vegandocs("partitioning")
# you can access the document using this line, as learned in first workshop
vignette("partitioning", package="vegan")

## Variation partitioning with all explanatory variables
spe.part.all <- varpart(spe.hel, envchem, envtopo)
spe.part.all
windows(title="Variation partitioning - all variables")
plot(spe.part.all, digits=2)  

## Variation partitioning with only significant variables
### Separate forward selection in each subset of environmental  variables
spe.chem <- rda(spe.hel~., data=envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envchem), 
           scope= formula(spe.chem), direction= "forward", R2scope=TRUE, pstep=1000)
names(envchem)
(envchem.pars <- envchem[, c( 4, 6, 7 )])
spe.topo <- rda(spe.hel~., data=envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envtopo), 
           scope= formula(spe.topo), direction= "forward", R2scope=TRUE, pstep=1000)
names(envtopo)
envtopo.pars <- envtopo[, c(1,2)]

### Variation partitioning
(spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars))
windows(title="Variation partitioning - parsimonious subsets")
plot(spe.part, digits=2)

# For a cool Venn diagram, use the package Venneuler
library(venneuler)
v <- venneuler(c(A=0.241,"A&B"= 0.233, B=0.112))
v$labels <- c("Chem\n24%","Topo \n11%") 
v$colors <- c(.78,.6) 
plot(v)

# Tests of all testable fractions
anova(rda(spe.hel, envchem.pars), step=1000) # Test of fractions [a+b]
anova(rda(spe.hel, envtopo.pars), step=1000) # Test of fractions [b+c]
env.pars <- cbind(envchem.pars, envtopo.pars) 
anova(rda(spe.hel, env.pars), step=1000) # Test of fractions [a+b+c]
anova(rda(spe.hel, envchem.pars, envtopo.pars), step=1000) # Test of fraction [a]
anova(rda(spe.hel, envtopo.pars, envchem.pars), step=1000) # Test of fraction [c]

# Challenge 3 ####
# Run the variation partitioning of the mite species abundances 
# with a first dataset for the significant substrate variables (SubsDens, WaterCont and Substrate)
# and a second dataset for the significant other variables (Shrud and Topo)
# Which proportion of the variation explain each dataset ?  
# Which are the significant fractions ?
# (answer below)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
str(mite.env)
(mite.subs<-mite.env[,c(1,2,3)])
(mite.other<-mite.env[,c(4,5)])
rda.mite.subs <- rda(mite.spe.hel~., data=mite.subs)
R2a.all.chem <- RsquareAdj(rda.mite.subs)$adj.r.squared
ordiR2step(rda(mite.spe.hel~1, data=mite.subs), 
           scope= formula(rda.mite.subs), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.subs)
(mite.subs.pars <- mite.subs[, c(2, 3)])
rda.mite.other <- rda(mite.spe.hel~., data=mite.other)
R2a.all.chem <- RsquareAdj(rda.mite.other)$adj.r.squared
ordiR2step(rda(mite.spe.hel~1, data=mite.other), 
           scope= formula(rda.mite.other), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.other)
(mite.other.pars <- mite.other[, c(1,2)])
### Variation partitioning
(mite.spe.part <- varpart(mite.spe.hel, ~WatrCont+Substrate, ~Shrub+Topo,
                          data=mite.env))
windows(title="Variation partitioning - parsimonious subsets")
plot(mite.spe.part, digits=2)
# Tests of all testable fractions
anova(rda(mite.spe.hel~ WatrCont+Substrate, data=mite.env), step=1000) # Test of fractions [a+b]
anova(rda(mite.spe.hel~Shrub+Topo, data=mite.env), step=1000) # Test of fractions [b+c]
(env.pars <- cbind(mite.env[,c(2,3,4,5)])) 
anova(rda(mite.spe.hel~ WatrCont+Substrate+Shrub+Topo, data=env.pars), step=1000) # Test of fractions [a+b+c]
anova(rda(mite.spe.hel~WatrCont+Substrate + Condition(Shrub+Topo), data=env.pars), step=1000) # Test of fraction [a]
anova(rda(mite.spe.hel~Shrub+Topo+ Condition(WatrCont+Substrate ), data=env.pars), step=1000) # Test of fraction [c]



# 3. MRT: Multivariate regression tree ####

# Prepare the data
env <- subset(env, select = -das) # remove the "distance from the source" variable

# Create the regression tree
# Using the 1 standard error from cross validated residual error criterion
doubs.mrt<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4)
# Click on the figure to select the number of partitions you want.
# When using the 1se from CVRE criterion, click on the orange dot
# The tree will then appear
summary(doubs.mrt)


# Using the CVRE criterion
doubs.mrt.cvre<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4)

# Choosing yourself the best number of partitions
doubs.mrt.4<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4)
summary(doubs.mrt.4)


# Find indicative species with MRT results
doubs.mrt.wrap<-MRT(doubs.mrt,percent=10,species=colnames(spe.hel))
summary(doubs.mrt.wrap)

doubs.mrt.indval<-indval(spe.hel,doubs.mrt$where)
doubs.mrt.indval$pval

doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)]
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]


# Challenge 5 ####
# Run the multivariate regression tree of the mite species abundances 
# What is the best sized tree according to Breiman et al. 2004 ?
# Which proportion of the variation is explained by the tree ?  
# What are the discriminant species ?
# (answer below)
#
#
#
#
#
#
#
#
#
mite.mrt<-mvpart(data.matrix(mite.spe.hel)~.,mite.env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(mite.spe.hel),xvmult=100,which=4)
summary(mite.mrt)

mite.mrt.wrap<-MRT(mite.mrt,percent=10,species=colnames(mite.spe.hel))
summary(mite.mrt.wrap)

mite.mrt.indval<-indval(mite.spe.hel,mite.mrt$where)
mite.mrt.indval$pval

mite.mrt.indval$maxcls[which(mite.mrt.indval$pval<=0.05)]
mite.mrt.indval$indcls[which(mite.mrt.indval$pval<=0.05)]

# 4. LDA: Linear Discriminat Analysis ####

#load spatial data to determine groups 
spa <- read.csv ('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpa.csv', row.names = 1)
spa <- spa[-8,]

#View spatial data 
View (spa)

#add site numbers
numbers<-(1:30)
numbers<-numbers[!numbers%in%8] 
spa$site<-numbers

#make groups based on lattitude y<82=group1, 82<y<156=group2, y>156=group3
spa.group<-ddply(.data=spa, .variables=.(x, y, site), .fun= summarise, group = if(y <= 82) 1 else if (y <= 156) 2 else 3)

#order by site
spa.group<-spa.group[with(spa.group, order(site)), ]

#run LDA
LDA<-lda(env,spa.group[,4])

#classification of the objects based on LDA
spe.class <- predict(LDA)$class

#posterior probabilities of the objects to belong to the groups
spe.post <- predict(LDA)$posterior

#table of prior versus predicted classifications
spe.table <- table(spa.group[,4], spe.class)

#proportion of correct classification
diag(prop.table(spe.table, 1))

#predicting classification of new data 
#read in new sites 
classify.me<-read.csv("classifyme.csv", header = T)

#predict grouping of new data
predict.group<-predict(LDA, newdata=classify.me)

#give classification for each new site
group.new<-predict.group$class

# Challenge 6 ####
# Run an LDA for the mite env data (only first two vars) based on four latitudinal groups you create from the mite.xy data set. 
data(mite.xy)
# What group was group 2 most incorrectly grouped into?
# What proportion of sites was correctly classified in group 1? group 2?
# (answer below)
#
#
#
#
#
#
#
#
#
mite.xy$site<-seq(1:70)
(max(mite.xy[,2])-min(mite.xy[,2]))/4
mite.xy.group<-ddply(.data=mite.xy, .variables=.(x, y, site), .fun= summarise, group = if(y <= 2.5) 1 else if (y <= 4.9) 2 else if (y <= 7.3) 3 else 4)
mite.xy.group<-mite.xy.group[with(mite.xy.group, order(site)), ]
LDA.mite<-lda(mite.env[,1:2],mite.xy.group[,4])
mite.class <- predict(LDA.mite)$class
mite.post <- predict(LDA.mite)$posterior
mite.table <- table(mite.xy.group[,4], mite.class)
diag(prop.table(mite.table, 1))

# 5. Other useful ordination methods####

?cca # Constrained Correspondence Analysis (CCA)

?metaMDS # Nonmetric Multidimesional Scaling

?CCorA # Canonical Correlation Analysis

help(coinertia, package=ade4) # Coinertia Analysis 

help(mfa, package=ade4) # Multiple Factorial Analysis

https://r-forge.r-project.org/R/?group_id=195  # packages AEM and PCNM contain functions to perform spatial analysis 

