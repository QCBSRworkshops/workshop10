# On utilisera nos données explicatives standardisées
# Enlever la variable "distance from the source" (colinéarité avec autres variables)
env.z <- subset(env.z, select = -das)

# Modèlise l'effect de tous les variables environnementales sur la composition en espèces des communautés
spe.rda <- rda(spe.hel ~ ., data = env.z)

summary(spe.rda)

summary(spe.rda)

# Sélection progressive de variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # modèle le plus simple
               scope = formula(spe.rda), # modèle "complet"
               direction = "forward",
               R2scope = TRUE, # limité par le R2 du modèle "complet"
               pstep = 1000,
               trace = FALSE) # mettre TRUE pour voir le processus du sélection!

# Vérifier le nouveau modèle avec les variables sélectionnée
fwd.sel$call

# Écrire notre nouveau modèle
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
# vérifier son R2 ajusté
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, permutations = 1000)

anova.cca(spe.rda.signif, permutations = 1000, by = "term")

anova.cca(spe.rda.signif, step = 1000, by = "axis")

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")

# Configuration des triplots RDA!

## extrait le % expliqué par les 2 premiers axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

## scores - ceux-ci sont des coordonnées dans l'espace RDA
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

## Configuration du graphique

# Commencer avec un graphique vide avec le cadrage, les axes, et des titres  
plot(spe.rda.signif,
     scaling = 1, # type de cadrage
     type = "none", # garder le graphique vide pour l'instant
     frame = FALSE,
     # fixer les limites des axes
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # ajouter des titres au graphique et aux axes
     main = "Triplot RDA - cadrage 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# ajouter des points pour les scores des sites
points(sc_si, 
       pch = 21, # fixer le symbole (ici, un cercle rempli avec une couleur)
       col = "black", # couleur de la bordure du cercle
       bg = "steelblue", # couleur pour remplir le cercle
       cex = 1.2) # taille du cercle
# ajouter des points pour les scores des espèces
points(sc_sp, 
       pch = 22, # fixer le symbole (ici, un carré rempli avec une couleur)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# ajouter du texte pour identifier les espèces 
text(sc_sp + c(0.03, 0.09), # ajuster les coordonnées pour éviter des chevauchements 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # gras
     cex = 0.6)
# ajouter des flèches pour les effets des variables explicatives
arrows(0,0, # chaque flèche commence à (0,0)
       sc_bp[,1], sc_bp[,2], # et finit au score de la variable
       col = "red", 
       lwd = 3)
# ajouter du texte pour identifier les variables explicatives
text(x = sc_bp[,1] -0.1, # ajuster les coordonnées pour éviter des chevauchements
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

## # Défi 1: Effectuer une RDA pour modèliser les effect des variables environnementales sur l’abondance des espèces d'acariens.

# Charger les données d'abondance des espèces d'acariens
data("mite")

# Charger les données environnementales
data("mite.env")

## decostand()
## rda()
## ordiR2step()
## anova.cca()
## ordiplot()

## # Défi 1: Solution! Spoilers ci-dessous!!

# Transformer les données d'abondances
mite.spe.hel <- decostand(mite, method = "hellinger")

# Standardiser les données environmentales quantiatives
mite.env$SubsDens <- decostand(mite.env$SubsDens, method = "standardize")
mite.env$WatrCont <- decostand(mite.env$WatrCont, method = "standardize")

# RDA avec tous les variables environnementales
mite.spe.rda <- rda(mite.spe.hel ~ ., data = mite.env)

# Sélection progressive des variables environnementales significatives
fwd.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.env),
                      scope = formula(mite.spe.rda),
                      direction = "forward",
                      R2scope = TRUE, pstep = 1000, trace = FALSE)
fwd.sel$call

# Refaire la RDA avec seulement les variables significatives
mite.spe.rda.signif <- rda(mite.spe.hel ~ WatrCont + Shrub +
                           Substrate + Topo + SubsDens,
                           data = mite.env)

# Calculer le R2 ajusté
RsquareAdj(mite.spe.rda.signif)$adj.r.squared


anova.cca(mite.spe.rda.signif, step = 1000)

# Cadrage 1
ordiplot(mite.spe.rda.signif,
         scaling = 1,
         main = "Cadrage 1")
# Cadrage 2
ordiplot(mite.spe.rda.signif,
         scaling = 2,
         main = "Cadrage 2")
