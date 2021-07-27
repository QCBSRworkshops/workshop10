# Partitionner la variation de la composition des espèces de poissons 
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!

# Visualiser les résultats avec un diagramme Venn
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # noms des matrices explicatives
     bg = c("seagreen3", "mediumpurple"), alpha = 80,
     digits = 2,
     cex = 1.5)

## # Tester la significativité

# [a+b] Chimie sans tenir compte de topographie
anova.cca(rda(spe.hel, env.chem))

# [b+c] Topographie sans tenir compte de chimie
anova.cca(rda(spe.hel, env.topo))

# [a] Chimie seulement
anova.cca(rda(spe.hel, env.chem, env.topo))

# [c] Topographie seulement
anova.cca(rda(spe.hel, env.chem, env.topo))

## # Défi 3
## 
## # Partitionnez la variation de l’abondance des espèces de mites entres des variables de substrat (`SubsDens`, `WatrCont`) et des variables spatiales significatives.
## # * Quelle est la proportion de variance expliquée par le substrat? par l'espace?
## # * Quelles sont les fractions significatives?
## # * Diagramme Venn des résultats!

data("mite.pcnm")

## ordiR2step()
## varpart()
## anova.cca(rda())
## plot()

## # Défi 3 - Solution! Spoilers ci-dessous!

# Modèle RDA avec tous les variables spatiales
full.spat <- rda(mite.spe.hel ~ ., data = mite.pcnm)

# Sélection progressive des variables spatiales
spat.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.pcnm),
               scope = formula(full.spat),
               R2scope = RsquareAdj(full.spat)$adj.r.squared,
               direction = "forward",
               trace = FALSE)
spat.sel$call

# Variables de substrat
mite.subs <- subset(mite.env, select = c(SubsDens, WatrCont))

# Variables spatiales significatives
mite.spat <- subset(mite.pcnm,
                    select = names(spat.sel$terminfo$ordered))
                    # pour rapidement accèder aux variables sélectionnées

mite.part <- varpart(mite.spe.hel, mite.subs, mite.spat)
mite.part$part$indfract # extraire résultats

anova.cca(rda(mite.spe.hel, mite.subs, mite.spat))

anova.cca(rda(mite.spe.hel, mite.spat, mite.subs))

plot(mite.part,
     digits = 2,
     Xnames = c("Subs", "Space"), # titre des fractions
     cex = 1.5,
     bg = c("seagreen3", "mediumpurple"), # ajoutez des couleurs!
     alpha = 80)
