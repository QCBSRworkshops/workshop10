knitr::include_graphics("images/PartialRDA.png")

# Subset environmental data into topography variables and chemistry variables
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Run a partial RDA
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

## # Alternative syntax for the partial RDA:
## spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo + # these are the effects we are interested in
##                        Condition(alt + pen + deb), # these are the covariates
##                        data = env.z)

summary(spe.partial.rda)

# Extract the model's adjusted R2
RsquareAdj(spe.partial.rda)$adj.r.squared

# Test whether the model is statistically significant
anova.cca(spe.partial.rda, step = 1000)

ordiplot(spe.partial.rda, 
         scaling = 2,
         main = "Doubs River partial RDA - Scaling 2")

## # Challenge 2:
## # Run a partial RDA to model the effects of environmental variables on mite species abundances (`mite.spe.hel`), while controlling for substrate variables (`SubsDens`, `WatrCont`, and `Substrate`).

## rda()
## summary()
## RsquareAdj()
## anova.cca() # hint: see the 'by' argument in ?anova.cca

## # Challenge 2: Solution! Spoilers ahead!!

# Compute partial RDA
mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Check summary
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")
