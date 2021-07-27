# Partition the variation in fish community composition
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!

# plot the variation partitioning Venn diagram
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

## # Significance testing

# [a+b] Chemistry without controlling for topography
anova.cca(rda(spe.hel, env.chem))

# [b+c] Topography without controlling for chemistry
anova.cca(rda(spe.hel, env.topo))

# [a] Chemistry alone
anova.cca(rda(spe.hel, env.chem, env.topo))

# [c] Topography alone
anova.cca(rda(spe.hel, env.topo, env.chem))

## # Challenge 3
## 
## # Partition the variation in the mite species data according to substrate variables (`SubsDens`, `WatrCont`) and significant spatial variables.
##         # What proportion of the variation is explained by substrate variables? By space?
##         # Which individual fractions are significant?
##         # Plot your results!

data("mite.pcnm")

## ordiR2step()
## varpart()
## anova.cca(rda())
## plot()

## # Challenge 3 - Solution! Spoilers ahead!

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
