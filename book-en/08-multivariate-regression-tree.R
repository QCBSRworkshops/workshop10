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

## # Challenge 4
## #
## # Create a multivariate regression tree for the mite data.
## # * Select the smallest tree within 1 SE of the CVRE.
## # * What is the proportion of variance (R2) explained by this tree?
## # * How many leaves does it have?
## # * What are the top 3 discriminant species?

data("mite")
data("mite.env")

## ?mvpart() # hint: pay attention to the 'xv' argument!
## summary()

## # Challenge 4: Solution! Spoilers ahead!

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
