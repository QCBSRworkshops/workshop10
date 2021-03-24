##Section: 01-preparing-for-the-workshop.R 

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


##Section: 02-introduction.R 




##Section: 03-data-exploration.R 

# Make sure the files are in your working directory! 
# If R cannot find the dataset, set your working directory with setwd()
# to the folder in which your data is stored (e.g. setwd("~/Desktop/workshop10"))

# Species community data frame (fish abundance)
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # Site number 8 contains no species, so we remove row 8 (site 8) 
# Be careful to only run this command line once as you are overwriting "spe" each time! 

# Environmental data frame: “DoubsEnv.csv”
env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # Remove corresponding abiotic data for site 8 (because removed from fish data). 
# Again, be careful to only run the last line once.

names(spe) # names of objects (species)
dim(spe) # dataset dimensions
head(spe) # look at first 5 rows

str(spe) # structure of objects in dataset
summary(spe) # summary statistics for all objects (min, mean, max, etc.)

# Count number of species frequencies in each abundance class
ab <- table(unlist(spe))
# Plot distribution of species frequencies
barplot(ab, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars

# Count the number of zeros in the dataset
sum(spe == 0) 

# Calculate proportion of zeros in the dataset
sum(spe == 0)/(nrow(spe)*ncol(spe))

# Apply Hellinger transformation to correct for the double zero problem
spe.hel <- decostand(spe, method = "hellinger")

names(env)
dim(env)
head(env)

str(env)
summary(env)

# We can visually look for correlations between variables:
heatmap(abs(cor(env)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0
round(apply(env.z, 2, mean), 1)

# and scaled to have a standard deviation of 1
apply(env.z, 2, sd)


##Section: 04-canonical-analysis.R 




##Section: 05-redundancy-analysis.R 

# sometimes cache needs to be set to true in the knitr setup chunk for this to take effect
# in xaringan::infinite_moon_reader()
library(knitr)
hook_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(more, x[lines], more)
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::include_graphics("images/RDA.png")

knitr::include_graphics("images/constrained_ord_diagram.png")

# We'll use our standardized environmental data, but we will remove 'das', which was correlated with many other variables:
env.z <- subset(env.z, select = -das)

# Model the effect of all environmental variables on fish community composition
spe.rda <- rda(spe.hel ~ ., data = env.z)

summary(spe.rda)

summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # lower model limit (simple!)
               scope = formula(spe.rda), # upper model limit (the "full" model)
               direction = "forward",
               R2scope = TRUE, # can't surpass the "full" model's R2
               pstep = 1000,
               trace = FALSE) # change to TRUE to see the selection process!

# Check the new model with forward-selected variables
fwd.sel$call

# Write our new model
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)

anova.cca(spe.rda.signif, step = 1000, by = "term")

anova.cca(spe.rda.signif, step = 1000, by = "axis")

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")

extract % explained by the first 2 axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

extract scores - these are coordinates in the RDA space
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(spe.rda.signif,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

# Challenge 1: Run an RDA to model the effects of environmental variables on mite species abundances.

# Load mite species abundance data
data("mite")

# Load environmental data
data("mite.env")

decostand()
rda()
ordiR2step()
anova.cca()
ordiplot()

# Challenge 1: Solution! Spoilers ahead!!

# Hellinger transform the community data
mite.spe.hel <- decostand(mite, method = "hellinger")

# Standardize quantitative environmental data
mite.env$SubsDens <- decostand(mite.env$SubsDens, method = "standardize")
mite.env$WatrCont <- decostand(mite.env$WatrCont, method = "standardize")

# Initial RDA with ALL of the environmental data
mite.spe.rda <- rda(mite.spe.hel ~ ., data = mite.env)

# Forward selection of environmental variables
fwd.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.env),
                      scope = formula(mite.spe.rda),
                      direction = "forward",
                      R2scope = TRUE, pstep = 1000, trace = FALSE)
fwd.sel$call

# Re-run the RDA with the significant variables
mite.spe.rda.signif <- rda(mite.spe.hel ~ WatrCont + Shrub +
                           Substrate + Topo + SubsDens,
                           data = mite.env)

# Find the adjusted R2 of the model with the retained env variables
RsquareAdj(mite.spe.rda.signif)$adj.r.squared


anova.cca(mite.spe.rda.signif, step = 1000)

# Scaling 1
ordiplot(mite.spe.rda.signif,
         scaling = 1,
         main = "Mite RDA - Scaling 1")
# Scaling 2
ordiplot(mite.spe.rda.signif,
         scaling = 2,
         main = "Mite RDA - Scaling 2")


##Section: 06-partial-redundancy-analysis.R 

knitr::include_graphics("images/PartialRDA.png")

# Subset environmental data into topography variables and chemistry variables
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Run a partial RDA
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

# Alternative syntax for the partial RDA:
spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo + # these are the effects we are interested in
                       Condition(alt + pen + deb), # these are the covariates
                       data = env.z)

summary(spe.partial.rda)

# Extract the model's adjusted R2
RsquareAdj(spe.partial.rda)$adj.r.squared

# Test whether the model is statistically significant
anova.cca(spe.partial.rda, step = 1000)

ordiplot(spe.partial.rda, 
         scaling = 2,
         main = "Doubs River partial RDA - Scaling 2")

# Challenge 2:
# Run a partial RDA to model the effects of environmental variables on mite species abundances (`mite.spe.hel`), while controlling for substrate variables (`SubsDens`, `WatrCont`, and `Substrate`).

rda()
summary()
RsquareAdj()
anova.cca() # hint: see the 'by' argument in ?anova.cca

# Challenge 2: Solution! Spoilers ahead!!

# Compute partial RDA
mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Check summary
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")


##Section: 07-variation-partitioning.R 

?varpart
vegandocs("partitioning.pdf")

#Variation partitioning with all explanatory variables
spe.part.all <- varpart(spe.hel, envchem, envtopo)
spe.part.all
windows(title="Variation partitioning - all variables")
plot(spe.part.all, digits=2)

#RDA of chemistry variables
spe.chem <- rda(spe.hel~., data=envchem)

#Select significant chemistry variables
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envchem),
           scope= formula(spe.chem), direction= "forward", R2scope=TRUE, pstep=1000)
names(envchem)
(envchem.pars <- envchem[, c( 4, 6, 7 )])

#RDA with other environmental variables
spe.topo <- rda(spe.hel~., data=envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envtopo),
           scope= formula(spe.topo), direction= "forward", R2scope=TRUE, pstep=1000)
names(envtopo)
envtopo.pars <- envtopo[, c(1,2)]

#Varpart
spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars)
windows(title="Variation partitioning - parsimonious subsets")
plot(spe.part, digits=2)

#Tests of significance
anova.cca(rda(spe.hel, envchem.pars), step=1000) # Test of fractions [a+b]
anova.cca(rda(spe.hel, envtopo.pars), step=1000) # Test of fractions [b+c]
env.pars <- cbind(envchem.pars, envtopo.pars)
anova.cca(rda(spe.hel, env.pars), step=1000) # Test of fractions [a+b+c]
anova.cca(rda(spe.hel, envchem.pars, envtopo.pars), step=1000) # Test of fraction [a]
anova.cca(rda(spe.hel, envtopo.pars, envchem.pars), step=1000) # Test of fraction [c]

str(mite.env)
(mite.subs=mite.env[,c(1,2,3)]) #First set of variables outlined in challenge
(mite.other=mite.env[,c(4,5)]) #Second set of variables outlined in challenge

#RDA for mite.subs
rda.mite.subs <- rda(mite.spe.hel~., data=mite.subs)
R2a.all.subs <- RsquareAdj(rda.mite.subs)$adj.r.squared

#Forward selection for mite.subs
ordiR2step(rda(mite.spe.hel~1, data=mite.subs),
           scope= formula(rda.mite.subs), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.subs)
(mite.subs.pars <- mite.subs[, c(2, 3)])

#RDA for mite.other
rda.mite.other <- rda(mite.spe.hel~., data=mite.other)
R2a.all.other <- RsquareAdj(rda.mite.other)$adj.r.squared

#Forward selection for mite.other
ordiR2step(rda(mite.spe.hel~1, data=mite.other),
           scope= formula(rda.mite.other), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.other)
(mite.other.pars <- mite.other[, c(1,2)])

#Variation partitioning
(mite.spe.part <- varpart(mite.spe.hel, ~WatrCont+Substrate, ~Shrub+Topo,
                          data=mite.env))
windows(title="Variation partitioning - parsimonious subsets")
plot(mite.spe.part, digits=2)

# Tests of all testable fractions
anova.cca(rda(mite.spe.hel~ WatrCont+Substrate, data=mite.env), step=1000) # Test of fractions [a+b]
anova.cca(rda(mite.spe.hel~Shrub+Topo, data=mite.env), step=1000) # Test of fractions [b+c]
(env.pars <- cbind(mite.env[,c(2,3,4,5)]))
anova.cca(rda(mite.spe.hel~ WatrCont+Substrate+Shrub+Topo, data=env.pars), step=1000) # Test of fractions [a+b+c]
anova.cca(rda(mite.spe.hel~WatrCont+Substrate + Condition(Shrub+Topo), data=env.pars), step=1000) # Test of fraction [a]
anova.cca(rda(mite.spe.hel~Shrub+Topo+ Condition(WatrCont+Substrate ), data=env.pars), step=1000) # Test of fraction [c]


##Section: 08-multivariate-regression-tree.R 

?mvpart

#Prepare the data: remove “distance from source”
env <- subset(env, select = -das)

# Create the regression tree
doubs.mrt <- mvpart(as.matrix(spe.hel) ~. ,env,
            legend=FALSE, margin=0.01, cp=0, xv="pick",
            xval=nrow(spe.hel), xvmult=100, which=4)

# Using the CVRE criterion
doubs.mrt.cvre <- mvpart(as.matrix(spe.hel)~., env,
                 legend=FALSE, margin=0.01, cp=0,xv="pick",
                 xval=nrow(spe.hel), xvmult=100,which=4)

# Choosing ourself the best number of partitions
doubs.mrt.4 <- mvpart(as.matrix(spe.hel)~., env,
              legend=FALSE, margin=0.01, cp=0, xv="pick",
              xval=nrow(spe.hel), xvmult=100,which=4)

summary(doubs.mrt)

# Find discriminant species with MRT results
doubs.mrt.wrap<-MRT(doubs.mrt,percent=10,species=colnames(spe.hel))
summary(doubs.mrt.wrap)

# Extract indval p-values
doubs.mrt.indval<-indval(spe.hel,doubs.mrt$where)
doubs.mrt.indval$pval

# Extract indicator species of each node, with its indval
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)]
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]

mite.mrt<-mvpart(data.matrix(mite.spe.hel)~.,mite.env,
legend=FALSE,margin=0.01,cp=0,xv="pick",
xval=nrow(mite.spe.hel),xvmult=100,which=4)
summary(mite.mrt)

mite.mrt.wrap<-MRT(mite.mrt,percent=10,species=colnames(mite.spe.hel))
summary(mite.mrt.wrap)

mite.mrt.indval<-indval(mite.spe.hel,mite.mrt$where)
mite.mrt.indval$pval

mite.mrt.indval$maxcls[which(mite.mrt.indval$pval<=0.05)]
mite.mrt.indval$indcls[which(mite.mrt.indval$pval<=0.05)]


##Section: 09-linear-discriminant-analysis.R 

#load spatial data to determine groups
spa <- read.csv ('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpa.csv', row.names = 1)
spa <- spa[,-8]

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

mite.xy$site<-seq(1:70)
(max(mite.xy[,2])-min(mite.xy[,2]))/4

mite.xy.group<-ddply(.data=mite.xy, .variables=.(x, y, site), .fun= summarise, group = if(y <= 2.5) 1 else if (y <= 4.9) 2 else if (y <= 7.3) 3 else 4)
mite.xy.group<-mite.xy.group[with(mite.xy.group, order(site)), ]

LDA.mite<-lda(mite.env[,1:2],mite.xy.group[,4])
mite.class <- predict(LDA.mite)$class
mite.post <- predict(LDA.mite)$posterior
mite.table <- table(mite.xy.group[,4], mite.class)
diag(prop.table(mite.table, 1))


##Section: 10-final-considerations.R 

?cca #(constrained correspondence analysis)
# Constrained Correspondence Analysis (CCA) is a canonical ordination method similar to RDA that preserve
# Chi-square distances among object (instead of Euclidean distances in RDA). This method is well suited for the
# analysis of large ecological gradients.



?CCorA # Canonical Correlation Analysis

# Canonical Correlation Analysis (CCorA) differs from RDA given that the two matrices are considered symmetric
# while in RDA the Y matrix is dependent on the X matrix. The main use of this technique is to test the
# significance of the correlation between two multidimensional data sets, then explore the structure of the data by
# computing the correlations (which are the square roots of the CCorA eigenvalues) that can be found between
# linear functions of two groups of descriptors.


help(coinertia, package=ade4) # Coinertia Analysis

#Coinertia Analysis (CoIA) is a symmetric canonical ordination method that is appropriate to compare pairs
# of data sets that play equivalent roles in the analysis. The method finds a common space onto which the objects
# and variables of these data sets can be projected and compared. Compared to CCorA, co-inertia analysis
# imposes no constraint regarding the number of variables in the two sets, so that it can be used to compare
# ecological communities even when they are species-rich. Co-inertia analysis is not well-suited, however, to
# analyse pairs of data sets that contain the same variables, because the analysis does not establish one-to-one
# correspondences between variables in the two data sets; the method does not ‘know’ that the first variable is the
# same in the first and the second data sets, and likewise for the other variables.


help(mfa, package=ade4) # Multiple Factorial Analysis

# Multiple factor analysis (MFA) can be used to compare several data sets describing the same objects. MFA
# consists in projecting objects and variables of two or more data sets on a global PCA, computed from all data
# sets, in which the sets receive equal weights.


# Spatial analysis can be performed using packages AEM and PCNM : http://r-forge.r-project.org/R/?group_id=195


##Section: 11-references.R 




