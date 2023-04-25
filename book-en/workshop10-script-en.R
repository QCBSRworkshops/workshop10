##Section: 01-preparing-for-the-workshop.R 

###Notice ###
#                                                                            #
#This is an automatically generated script based on the code chunks from the #
#book for this workshop.                                                     #
#                                                                            #
#It is minimally annotated to allow participants to provide their comments:  #
#a practice that we highly encourage.                                        #
#                                                                            #
#Note that the solutions to the challenges are also included in this script. #
#When solving the challenges by yourself, attempt to not scroll and peek at  #
#the solutions.                                                              #
#                                                                            #
#Happy coding!                                                               #


# Install the required packages
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("ggplot2")

# install mvpart from package archive file
install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Load the required packages
library(labdsv)
library(vegan)
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

# Custom triplot code!

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

# Partition the variation in fish community composition
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!

# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

# Significance testing

# [a+b] Chemistry without controlling for topography
anova.cca(rda(spe.hel, env.chem))

# [b+c] Topography without controlling for chemistry
anova.cca(rda(spe.hel, env.topo))

# [a] Chemistry alone
anova.cca(rda(spe.hel, env.chem, env.topo))

# [c] Topography alone
anova.cca(rda(spe.hel, env.topo, env.chem))

# Challenge 3

# Partition the variation in the mite species data according to substrate variables (`SubsDens`, `WatrCont`) and significant spatial variables.
        # What proportion of the variation is explained by substrate variables? By space?
        # Which individual fractions are significant?
        # Plot your results!

data("mite.pcnm")

ordiR2step()
varpart()
anova.cca(rda())
plot()

# Challenge 3 - Solution! Spoilers ahead!

# Step 1: Forward selection!

# Write full RDA model with all variables
full.spat <- rda(mite.spe.hel ~ ., data = mite.pcnm)

# Forward selection of spatial variables
spat.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.pcnm),
               scope = formula(full.spat),
               R2scope = RsquareAdj(full.spat)$adj.r.squared,
               direction = "forward",
               trace = FALSE)
spat.sel$call

# Step 2: Group variables of interest.

# Subset environmental data to retain only substrate variables
mite.subs <- subset(mite.env, select = c(SubsDens, WatrCont))

# Subset to keep only selected spatial variables
mite.spat <- subset(mite.pcnm,
                    select = names(spat.sel$terminfo$ordered))
                    # a faster way to access the selected variables!

# Step 3: Partition the variation in species abundances.
mite.part <- varpart(mite.spe.hel, mite.subs, mite.spat)
mite.part$part$indfract # access results!

# Step 4: Significance testing
# [a]: Substrate only
anova.cca(rda(mite.spe.hel, mite.subs, mite.spat))

# [c]: Space only
anova.cca(rda(mite.spe.hel, mite.spat, mite.subs))

# Step 5: Plot
plot(mite.part, 
     digits = 2, 
     Xnames = c("Subs", "Space"), # label the fractions
     cex = 1.5,
     bg = c("seagreen3", "mediumpurple"), # add colour!
     alpha = 80) # adjust transparency


##Section: 08-multivariate-regression-tree.R 

# First, remove the "distance from source" variable
env <- subset(env, select = -das)

# Create multivariate regression tree
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # interactively select best tree
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0)

# Select the solution we want
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "1se", # select smallest tree within 1 se
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0)

# Trying 10 groups
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 10, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)

# Trying fewer groups
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 4, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)

# Checking the complexity parameter
doubs.mrt$cptable

# Checking the tree result summary
summary(doubs.mrt)

# Calculate indicator values (indval) for each species
doubs.mrt.indval <- indval(spe.hel, doubs.mrt$where)

# Extract the significant indicator species (and which node they represent)
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval <= 0.05)]

# Extract their indicator values
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval <= 0.05)]

# Challenge 4
#
# Create a multivariate regression tree for the mite data.
# * Select the smallest tree within 1 SE of the CVRE.
# * What is the proportion of variance (R2) explained by this tree?
# * How many leaves does it have?
# * What are the top 3 discriminant species?

data("mite")
data("mite.env")

?mvpart() # hint: pay attention to the 'xv' argument!
summary()

# Challenge 4: Solution! Spoilers ahead!

mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01, cp = 0,
                   prn = FALSE)

# Calculate indicator values (indval) for each species
mite.mrt.indval <- indval(mite.spe.hel, mite.mrt$where)

# Extract the significant indicator species (and which node they represent)
mite.mrt.indval$maxcls[which(mite.mrt.indval$pval <= 0.05)]

# Extract their indicator values
mite.mrt.indval$indcls[which(mite.mrt.indval$pval <= 0.05)]


##Section: 09-linear-discriminant-analysis.R 

# load spatial data for Doubs sites
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # add site numbers
spa <- spa[-8,] # remove site 8

# group sites based on latitude
spa$group <- NA # create "group" column
spa$group[which(spa$y < 82)] <- 1
spa$group[which(spa$y > 82 & spa$y < 156)] <- 2
spa$group[which(spa$y > 156)] <- 3

ggplot(data = spa) +
  geom_point(aes(x = x, 
                 y = y, 
                 col = as.factor(group)), 
             size = 4) +
  labs(color = "Groups", 
       x = "Longitude", 
       y = "Latitude") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # formatting the plot to make it pretty
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# run the LDA grouping sites into latitude groups based on env data
LDA <- lda(env, spa$group)

# predict the groupings
lda.plotdf <- data.frame(group = spa$group, lda = predict(LDA)$x)

# Plot the newly reorganised sites according to the LDA
ggplot(lda.plotdf) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 col = factor(group)), 
             size = 4) +
  labs(color = "Groups") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # formatting the plot to make it pretty
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# Classification of the objects based on the LDA
spe.class <- predict(LDA)$class

# Posterior probabilities that the objects belong to those groups
spe.post <- predict(LDA)$posterior

# Table of prior vs. predicted classifications
(spe.table <- table(spa$group, spe.class))

# Proportion of corrected classification
diag(prop.table(spe.table, 1))

# Load the new site data
classify.me <- read.csv("data/classifyme.csv", header = TRUE)
# remove das
classify.me <- subset(classify.me, select = -das)
  
# Predict grouping of new sites
predict.group <- predict(LDA, newdata = classify.me)

# View site classification
predict.group$class

data(mite.xy)

lda()
predict()
table()
diag()

# assign numbers to sites
mite.xy$site <- 1:nrow(mite.xy)

# figure out the spacing to make 4 equally spaced latitude groups 
(max(mite.xy[,2])-min(mite.xy[,2]))/4

# use this to group sites into 4 latitude groups
mite.xy$group <- NA # create "group" column
mite.xy$group[which(mite.xy$y < 2.5)] <- 1
mite.xy$group[which(mite.xy$y >= 2.5 & mite.xy$y < 4.9)] <- 2
mite.xy$group[which(mite.xy$y >= 4.9 & mite.xy$y < 7.3)] <- 3
mite.xy$group[which(mite.xy$y >= 7.3)] <- 4

LDA.mite <- lda(mite.env[,1:2], mite.xy$group)

# group sites based on the LDA
mite.class <- predict(LDA.mite)$class

# get the table of prior versus predicted classifications
(mite.table <- table(mite.xy$group, mite.class))

# proportion of correct classification
diag(prop.table(mite.table, 1))


##Section: 10-final-considerations.R 




##Section: 11-references.R 




