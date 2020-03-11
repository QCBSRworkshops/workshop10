## ----setup, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  comment = "#",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = TRUE,
  fig.width=6, fig.height=6,
  fig.retina = 3,
  fig.align = 'center'
)
options(repos=structure(c(CRAN="http://cran.r-project.org")))


## ----install_pkgs, echo = FALSE, results = "asis"-----------------------------
cat(
  qcbsRworkshops::first_slides(10, c('Hmisc', 'labdsv', 'MASS', 'vegan'))
)


## ---- echo = TRUE-------------------------------------------------------------
# Make sure the files are in your working directory!
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # remove site with no data


## ---- echo = TRUE-------------------------------------------------------------
env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # remove site with no data


## ---- echo = TRUE-------------------------------------------------------------
names(spe) # names of objects (species)
dim(spe) # dataset dimensions


## ---- echo = TRUE, results = 'hide'-------------------------------------------
head(spe) # look at first 5 rows
str(spe) # structure of objects in dataset
summary(spe) # summary statistics for all objects (min, mean, max, etc.)


## ---- echo = TRUE, results = 'hide', fig.width = 6, fig.height = 3.5----------
# Count number of species frequencies in each abundance class
ab <- table(unlist(spe))
# Plot distribution of species frequencies
barplot(ab, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars


## -----------------------------------------------------------------------------
sum(spe == 0)


## -----------------------------------------------------------------------------
sum(spe==0)/(nrow(spe)*ncol(spe))


## -----------------------------------------------------------------------------
# The decostand() function in the vegan package makes this easy for us:
library(vegan)
spe.hel <- decostand(spe, method = "hellinger")


## ---- echo = TRUE-------------------------------------------------------------
names(env) # names of objects (environmental variables)
dim(env) # dataset dimensions
head(env) # look at first 5 rows


## ---- echo = TRUE, results = 'hide'-------------------------------------------
str(env) # structure of objects in dataset
summary(env) # summary statistics for all objects (min, mean, max, etc.)


## ---- fig.height = 5, fig.width = 9-------------------------------------------
# We can visually look for correlations between variables:
pairs(env)


## -----------------------------------------------------------------------------
# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0:
round(apply(env.z, 2, mean), 1)

# and scaled to have a standard deviation of 1
apply(env.z, 2, sd)


## -----------------------------------------------------------------------------
# We'll use our standardized environmental data
# But we will remove 'das', which was correlated with many other variables:
env.z <- subset(env.z, select = -das)


## -----------------------------------------------------------------------------
# Model the effect of all environmental variables on fish community composition
spe.rda <- rda(spe.hel ~ ., data = env.z)


## ---- eval = FALSE, results = 'hide'------------------------------------------
## summary(spe.rda, display = NULL)


## -----------------------------------------------------------------------------
# Forward selection of variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # lower model limit (simple!)
               scope = formula(spe.rda), # upper model limit (the "full" model)
               direction = "forward",
               R2scope = TRUE, # can't surpass the "full" model's R2
               pstep = 1000,
               trace = FALSE) # change to TRUE to see the selection process!


## -----------------------------------------------------------------------------
# Check the new model with forward-selected variables
fwd.sel$call


## -----------------------------------------------------------------------------
# Write our new model
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
# check the adjusted R2
RsquareAdj(spe.rda.signif)


## -----------------------------------------------------------------------------
anova.cca(spe.rda.signif, step = 1000)


## ---- results = 'hide'--------------------------------------------------------
anova.cca(spe.rda.signif, step = 1000, by = "term")


## ---- fig.height = 6.5, fig.width = 6, strip.white = TRUE---------------------
ordiplot(spe.rda.signif,
         scaling = 1,
         type = "text")


## ---- fig.height = 6.5, fig.width = 6, strip.white = TRUE---------------------
ordiplot(spe.rda.signif,
         scaling = 2,
         type = "text")


## ---- echo = FALSE------------------------------------------------------------
## extract % explained
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)
## scores
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)
## plot
plot(spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)"), xlim=c(-1,1), ylim=c(-1,1))
points(sc_si, pch=21, col="black", bg="steelblue", cex=1.2)
points(sc_sp, pch=22, col="black", bg = "#f2bd33", cex=1.2)
text(sc_sp, labels = rownames(sc_sp), col="black", cex=0.6)
arrows(0,0, sc_bp[,1], sc_bp[,2], col="red", lwd = 3)
text(x = sc_bp[,1] -0.1, y = sc_bp[,2] - 0.03, labels=rownames(sc_bp), col="red", cex=1, font = 2)


## -----------------------------------------------------------------------------
# Load mite species abundance data
data("mite")

# Load environmental data
data("mite.env")


## ---- eval = FALSE------------------------------------------------------------
## decostand()
## rda()
## ordiR2step()
## anova.cca()
## ordiplot()


## -----------------------------------------------------------------------------
# Hellinger transform the community data
mite.spe.hel <- decostand(mite, method = "hellinger")

# Standardize quantiative environmental data
mite.env$SubsDens <- decostand(mite.env$SubsDens, method = "standardize")
mite.env$WatrCont <- decostand(mite.env$WatrCont, method = "standardize")


## -----------------------------------------------------------------------------
# Initial RDA with ALL of the environmental data
mite.spe.rda <- rda(mite.spe.hel ~ ., data = mite.env)

# Forward selection of environmental variables
fwd.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.env),
                      scope = formula(mite.spe.rda),
                      direction = "forward",
                      R2scope = TRUE, pstep = 1000, trace = FALSE)
fwd.sel$call


## -----------------------------------------------------------------------------
# Re-run the RDA with the significant variables
mite.spe.rda.signif <- rda(mite.spe.hel ~ WatrCont + Shrub +
                           Substrate + Topo + SubsDens,
                           data = mite.env)

# Find the adjusted R2 of the model with the retained env variables
RsquareAdj(mite.spe.rda.signif)$adj.r.squared



## -----------------------------------------------------------------------------
anova.cca(mite.spe.rda.signif, step = 1000)


## ---- fig.height = 6.5--------------------------------------------------------
ordiplot(mite.spe.rda.signif,
         scaling = 1,
         main = "Mite RDA - Scaling 1")


## ---- fig.height = 6.5--------------------------------------------------------
ordiplot(mite.spe.rda.signif,
         scaling = 2,
         main = "Mite RDA - Scaling 2")


## -----------------------------------------------------------------------------
# Subset environmental data into topography variables and chemistry variables
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Partial RDA
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)


## ---- eval = FALSE------------------------------------------------------------
## spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
##                        Condition(alt + pen + deb),
##                        data = env.z)


## ---- eval = FALSE, results = 'hide'------------------------------------------
## summary(spe.partial.rda, display = NULL)


## -----------------------------------------------------------------------------
RsquareAdj(spe.partial.rda)$adj.r.squared


## -----------------------------------------------------------------------------
anova.cca(spe.partial.rda, step = 1000)


## -----------------------------------------------------------------------------
ordiplot(spe.partial.rda, scaling = 2,
         main = "Doubs River partial RDA - Scaling 2")


## ---- eval = FALSE------------------------------------------------------------
## rda()
## summary()
## RsquareAdj()
## anova.cca() # hint: see the 'by' argument in ?anova.cca


## ---- results = 'hide'--------------------------------------------------------
# Compute partial RDA
mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Check summary
summary(mite.spe.subs, display = NULL)


## -----------------------------------------------------------------------------
RsquareAdj(mite.spe.subs)$adj.r.squared


## -----------------------------------------------------------------------------
anova.cca(mite.spe.subs, step = 1000)


## -----------------------------------------------------------------------------
anova.cca(mite.spe.subs, step = 1000, by = "axis")


## -----------------------------------------------------------------------------
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!


## ---- strip.white = TRUE, fig.width = 6, fig.height = 6-----------------------
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)


## -----------------------------------------------------------------------------
# [a+b] Chemistry without controlling for topography
anova.cca(rda(spe.hel, env.chem))


## -----------------------------------------------------------------------------
# [b+c] Topography without controlling for chemistry
anova.cca(rda(spe.hel, env.topo))


## -----------------------------------------------------------------------------
# [a] Chemistry alone
anova.cca(rda(spe.hel, env.chem, env.topo))


## -----------------------------------------------------------------------------
# [c] Topography alone
anova.cca(rda(spe.hel, env.topo, env.chem))


## -----------------------------------------------------------------------------
data("mite.pcnm")


## ---- eval = FALSE------------------------------------------------------------
## ordiR2step()
## varpart()
## anova.cca(rda())
## plot()


## -----------------------------------------------------------------------------
# Write full RDA model with all variables
full.spat <- rda(mite.spe.hel ~ ., data = mite.pcnm)

# Forward selection of spatial variables
spat.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.pcnm),
               scope = formula(full.spat),
               R2scope = RsquareAdj(full.spat)$adj.r.squared,
               direction = "forward",
               trace = FALSE)
spat.sel$call


## -----------------------------------------------------------------------------
# Subset environmental data to retain only substrate variables
mite.subs <- subset(mite.env, select = c(SubsDens, WatrCont))

# Subset to keep only selected spatial variables
mite.spat <- subset(mite.pcnm,
                    select = names(spat.sel$terminfo$ordered))
                    # a faster way to access the selected variables!


## -----------------------------------------------------------------------------
mite.part <- varpart(mite.spe.hel, mite.subs, mite.spat)
mite.part$part$indfract # access results!


## ---- results = 'hide'--------------------------------------------------------
# [a]: Substrate only
anova.cca(rda(mite.spe.hel, mite.subs, mite.spat))
# p = 0.001 ***

# [c]: Space only
anova.cca(rda(mite.spe.hel, mite.spat, mite.subs))
# p = 0.001 ***


## ---- fig.height=6, fig.width=6-----------------------------------------------
plot(mite.part, digits = 2, cex = 1.5,
     bg = c("pink", "skyblue"), alpha = 90) # add colour!


## ----mvpart_install-----------------------------------------------------------
remotes::install_github("cran/mvpart")
library(mvpart)


## ---- results = 'hide', fig.show = 'hide', eval = F---------------------------
## # First, remove the "distance from source" variable
## env <- subset(env, select = -das)
## 
## # Create multivariate regression tree
## # library(mvpart)
## doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
##                     xv = "pick", # interactively select best tree
##                     xval = nrow(spe.hel), # number of cross-validations
##                     xvmult = 100, # number of multiple cross-validations
##                     which = 4, # plot both node labels
##                     legend = FALSE, margin = 0.01, cp = 0)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5.5, fig.align = 'center'----
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # interactively select tree size
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0, plot.add = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5.5, fig.align = 'center', eval =F----
## doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
##                     xv = "pick", # interactively select tree size
##                     xval = nrow(spe.hel), # number of cross-validations
##                     xvmult = 100, # number of multiple cross-validations
##                     which = 4, # plot both node labels
##                     legend = FALSE, margin = 0.01, cp = 0,
##                     plot.add = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5.5, fig.width = 5.5-------
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "1se", # interactively select tree size
                    xval = nrow(spe.hel), # number of cross-validations
                    xvse = 1,
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 4.5, fig.width = 12--------
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 10, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 4.5, fig.width =8----------
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 4, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## -----------------------------------------------------------------------------
doubs.mrt$cptable


## -----------------------------------------------------------------------------
summary(doubs.mrt)


## ----MVPARTwrap_install-------------------------------------------------------
remotes::install_github("cran/MVPARTwrap")
library(MVPARTwrap)


## ---- results = 'hide'--------------------------------------------------------
# Generate a nicer and more informative output
doubs.mrt.wrap <- MRT(doubs.mrt, percent = 10, species = colnames(spe.hel))

# Access the full output:
summary(doubs.mrt.wrap)


## -----------------------------------------------------------------------------
summary(doubs.mrt.wrap)


## ----labdsv-------------------------------------------------------------------
library(labdsv)

# Calculate indicator values (indval) for each species
doubs.mrt.indval <- indval(spe.hel, doubs.mrt$where)

# Extract the significant indicator species (and which node they represent)
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval <= 0.05)]

# Extract their indicator values
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval <= 0.05)]


## -----------------------------------------------------------------------------
data("mite")
data("mite.env")


## ---- eval = FALSE------------------------------------------------------------
## ?mvpart() # hint: pay attention to the 'xv' argument!
## ?MRT()
## summary()


## ---- results = 'hide', fig.height = 4.5, fig.width = 4.5---------------------
mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01, cp = 0,
                   prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5-----------
mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01,
                   cp = 0, prn = FALSE)


## ---- results = 'hide'--------------------------------------------------------
# Generate nicer MRT output
mite.mrt.wrap <- MRT(mite.mrt,
                     percent = 10,
                     species = colnames(mite.spe.hel))

# Look at discriminant species table from MRT output
summary(mite.mrt.wrap)


## -----------------------------------------------------------------------------
# load spatial data for Doubs sites
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # add site numbers
spa <- spa[-8,] # remove site 8


## -----------------------------------------------------------------------------
# group sites based on latitude
spa$group <- NA # create "group" column
spa$group[which(spa$y < 82)] <- 1
spa$group[which(spa$y > 82 & spa$y < 156)] <- 2
spa$group[which(spa$y > 156)] <- 3


## ---- fig.width = 5.5, fig.height = 5.5---------------------------------------
plot(spa$x, spa$y, col = spa$group, pch = 16, cex = 1.5)


## -----------------------------------------------------------------------------
# load required library
library(MASS)

# run the LDA grouping sites into latitude groups based on env data
LDA <- lda(env, spa$group)


## -----------------------------------------------------------------------------
# Classification of the objects based on the LDA
spe.class <- predict(LDA)$class

# Posterior probabilities that the objects belong to those groups
spe.post <- predict(LDA)$posterior

# Table of prior vs. predicted classifications
(spe.table <- table(spa$group, spe.class))

# Proportion of corrected classification
diag(prop.table(spe.table, 1))


## -----------------------------------------------------------------------------
# Load the new site data
classify.me <- read.csv("data/classifyme.csv", header = TRUE)
# classify.me <- classify.me[,-1] # remove das variable

# Predict grouping of new sites
predict.group <- predict(LDA, newdata = classify.me)

# View site classification
predict.group$class


## -----------------------------------------------------------------------------
data(mite.xy)


## ---- eval = FALSE------------------------------------------------------------
## lda()
## predict()
## table()
## diag()


## -----------------------------------------------------------------------------
# assign numbers to sites
mite.xy$site <- 1:nrow(mite.xy)

# find latitudinal range for each group
(max(mite.xy[,2])-min(mite.xy[,2]))/4

# group sites into 4 latitude groups
# group sites based on latitude
mite.xy$group <- NA # create "group" column
mite.xy$group[which(mite.xy$y < 2.5)] <- 1
mite.xy$group[which(mite.xy$y >= 2.5 & mite.xy$y < 4.9)] <- 2
mite.xy$group[which(mite.xy$y >= 4.9 & mite.xy$y < 7.3)] <- 3
mite.xy$group[which(mite.xy$y >= 7.3)] <- 4


## -----------------------------------------------------------------------------
LDA.mite <- lda(mite.env[,1:2], mite.xy$group)


## -----------------------------------------------------------------------------
# classification of the objects based on LDA
mite.class <- predict(LDA.mite)$class
# table of prior versus predicted classifications
(mite.table <- table(mite.xy$group, mite.class))
# proportion of correct classification
diag(prop.table(mite.table, 1))

