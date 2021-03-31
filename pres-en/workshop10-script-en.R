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

# Standard procedure to check and install packages and their dependencies, if needed.

list.of.packages <- c("remotes", "Hmisc", "labdsv", "MASS", "vegan")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE)
  print(paste0("The following package was installed:", new.packages))
} else if(length(new.packages) == 0) {
    print("All packages were already installed previously")
  }

# Load all required libraries at once
lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE)

# Make sure the files are in your working directory!
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # remove site with no data

env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # remove site with no data

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

sum(spe == 0)

sum(spe==0)/(nrow(spe)*ncol(spe))

# The decostand() function in the vegan package makes this easy for us:
library(vegan)
spe.hel <- decostand(spe, method = "hellinger")

names(env) # names of objects (environmental variables)
dim(env) # dataset dimensions
head(env) # look at first 5 rows

str(env) # structure of objects in dataset
summary(env) # summary statistics for all objects (min, mean, max, etc.)

# We can visually look for correlations between variables:
heatmap(abs(cor(env)), Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0:
round(apply(env.z, 2, mean), 1)

# and scaled to have a standard deviation of 1
apply(env.z, 2, sd)

# We'll use our standardized environmental data
# But we will remove 'das', which was correlated with many other variables:
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
# check the adjusted R2
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)

anova.cca(spe.rda.signif, step = 1000, by = "term")

ordiplot(spe.rda.signif,
         scaling = 1,
         type = "text")

ordiplot(spe.rda.signif,
         scaling = 2,
         type = "text")

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

# Load mite species abundance data
data("mite")

# Load environmental data
data("mite.env")

decostand()
rda()
ordiR2step()
anova.cca()
ordiplot()

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

ordiplot(mite.spe.rda.signif,
         scaling = 1,
         frame = F,
         main = "Mite RDA - Scaling 1")

ordiplot(mite.spe.rda.signif,
         scaling = 2,
         frame = F,  
         main = "Mite RDA - Scaling 2")

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
         frame = FALSE, 
         main = "Doubs River partial RDA - Scaling 2")

rda()
summary()
RsquareAdj()
anova.cca() # hint: see the 'by' argument in ?anova.cca

# Compute partial RDA
mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Check summary
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")

spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!

plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

# [a+b] Chemistry without controlling for topography
anova.cca(rda(spe.hel, env.chem))

# [b+c] Topography without controlling for chemistry
anova.cca(rda(spe.hel, env.topo))

# [a] Chemistry alone
anova.cca(rda(spe.hel, env.chem, env.topo))

# [c] Topography alone
anova.cca(rda(spe.hel, env.topo, env.chem))

data("mite.pcnm")

ordiR2step()
varpart()
anova.cca(rda())
plot()

# Write full RDA model with all variables
full.spat <- rda(mite.spe.hel ~ ., data = mite.pcnm)

# Forward selection of spatial variables
spat.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.pcnm),
               scope = formula(full.spat),
               R2scope = RsquareAdj(full.spat)$adj.r.squared,
               direction = "forward",
               trace = FALSE)
spat.sel$call

# Subset environmental data to retain only substrate variables
mite.subs <- subset(mite.env, select = c(SubsDens, WatrCont))

# Subset to keep only selected spatial variables
mite.spat <- subset(mite.pcnm,
                    select = names(spat.sel$terminfo$ordered))
                    # a faster way to access the selected variables!

mite.part <- varpart(mite.spe.hel, mite.subs, mite.spat)
mite.part$part$indfract # access results!

anova.cca(rda(mite.spe.hel, mite.subs, mite.spat))

anova.cca(rda(mite.spe.hel, mite.spat, mite.subs))

plot(mite.part, 
     digits = 2, 
     Xnames = c("Subs", "Space"), # label the fractions
     cex = 1.5,
     bg = c("seagreen3", "mediumpurple"), # add colour!
     alpha = 80) # adjust transparency

remotes::install_github("cran/mvpart")
library(mvpart)

# remove this chunk when the MRT section is back in. It is needed for the LDA section.
env <- subset(env, select = -das)

# First, remove the "distance from source" variable
env <- subset(env, select = -das)

# Create multivariate regression tree
# library(mvpart)
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # interactively select best tree
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # interactively select tree size
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0, plot.add = FALSE)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # interactively select tree size
                    xval = nrow(spe.hel), # number of cross-validations
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0,
                    plot.add = FALSE)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "1se", # interactively select tree size
                    xval = nrow(spe.hel), # number of cross-validations
                    xvse = 1,
                    xvmult = 100, # number of multiple cross-validations
                    which = 4, # plot both node labels
                    legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)

mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 10, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)

mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 4, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)

doubs.mrt$cptable

summary(doubs.mrt)

library(labdsv)

# Calculate indicator values (indval) for each species
doubs.mrt.indval <- indval(spe.hel, doubs.mrt$where)

# Extract the significant indicator species (and which node they represent)
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval <= 0.05)]

# Extract their indicator values
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval <= 0.05)]

data("mite")
data("mite.env")

?mvpart() # hint: pay attention to the 'xv' argument!
summary()

mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01, cp = 0,
                   prn = FALSE)

mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01,
                   cp = 0, prn = FALSE)

# Calculate indicator values (indval) for each species
mite.mrt.indval <- indval(mite.spe.hel, mite.mrt$where)

# Extract the significant indicator species (and which node they represent)
mite.mrt.indval$maxcls[which(mite.mrt.indval$pval <= 0.05)]

# Extract their indicator values
mite.mrt.indval$indcls[which(mite.mrt.indval$pval <= 0.05)]

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

# load required library
library(MASS)

# run the LDA grouping sites into latitude groups based on env data
LDA <- lda(env, spa$group)

# predict the groupings
lda.plotdf <- data.frame(group = spa$group, lda = predict(LDA)$x)

library(ggplot2)
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

# proportion of correct classification
diag(prop.table(mite.table, 1))
