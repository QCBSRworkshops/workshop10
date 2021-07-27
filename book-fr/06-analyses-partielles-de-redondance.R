knitr::include_graphics("images/PartialRDA.png")

# Divisez le tableau de données environnementales en deux:
# variables topographiques et chimiques
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Faire la RDA partielle
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

## # Syntaxe alternative
## spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
##                        Condition(alt + pen + deb), # covariables ici
##                        data = env.z)

summary(spe.partial.rda)

# Extraire le R2 ajusté du modèle
RsquareAdj(spe.partial.rda)$adj.r.squared

# Évaluer la significativité statistique du modèle
anova.cca(spe.partial.rda, step = 1000)

ordiplot(spe.partial.rda, 
         scaling = 2,
         main = "Rivière Doubs - Cadrage 2")

## # Défi 2:
## # Effectuez une RDA partielle de l’abondance des espèces de mites (`mite.spe.hel`) en fonction des variables environnementales, tenant compte de l’effet du substrat (`SubsDens`, `WaterCont` and `Substrate`).
## # * Quel pourcentage de variance est expliqué par les variables environnementales?
## # * Le modèle est-il significatif?
## # * Quels sont les axes significatifs?

## rda()
## summary()
## RsquareAdj()
## anova.cca() # voir l'argument 'by' dans ?anova.cca

## # Défi 2: Solution! Spoilers ci-dessous!!

mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Extraire les résultats
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")
