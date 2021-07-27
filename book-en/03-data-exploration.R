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

## str(spe) # structure of objects in dataset
## summary(spe) # summary statistics for all objects (min, mean, max, etc.)

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

## str(env)
## summary(env)

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
