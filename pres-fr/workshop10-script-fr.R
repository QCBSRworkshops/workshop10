# Standard procedure to check and install packages and their dependencies, if needed.

# remotes::install_github("cran/mvpart")

list.of.packages <- c("remotes", "Hmisc", "labdsv", "MASS", "vegan"
                      #"mvpart"
                      )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE)
  print(paste0("The following package was installed:", new.packages))
} else if(length(new.packages) == 0) {
    print("All packages were already installed previously")
  }

# Load all required libraries at once
lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE)

# Assurez vous que les fichiers se trouvent dans votre répertoire de travail!
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # Supprimer site 8 (pas d'espèces).

env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # Supprimer site 8 (pas d'espèces).

names(spe) # noms d'objets (espèces)
dim(spe) # dimensions de la matrice

head(spe) # 5 premières lignes
str(spe) # structure d'objets de la matrice
summary(spe) # statistiques descriptives des objets (min, moyenne, max, etc.)

# Compter la fréquence d'espèces dans chaque classe d'abondance
ab <- table(unlist(spe))
# Visualiser cette distribution
barplot(ab, las = 1,
        xlab = "Abundance class", ylab = "Frequency",
        col = grey(5:0/5))

sum(spe == 0)

sum(spe == 0)/(nrow(spe)*ncol(spe))

# la fonction decostand() dans la libraire vegan nous facilite la tâche:
library(vegan)
spe.hel <- decostand(spe, method = "hellinger")

names(env) # noms des objets (variables environnementales)
dim(env) # dimensions de la matrice
head(env) # 5 premières lignes

str(env) # structure des objets
summary(env) # statistiques descriptives (min, moyenne, max, etc.)

# On peut également détecter (visuellement) les colinéarités entres variables:
heatmap(abs(cor(env)), # corrélation de Pearson (note: ce sont des valeurs absolues!)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topright",
       title = "R de Pearson",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

# standardiser les données
env.z <- decostand(env, method = "standardize")

# centrer les données (moyenne ~ 0)
round(apply(env.z, 2, mean), 1)

# réduire les données (écart type = 1)
apply(env.z, 2, sd)

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

fwd.sel$call

# Écrire notre nouveau modèle
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
# vérifier son R2 ajusté
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, permutations = 1000)

anova.cca(spe.rda.signif, permutations = 1000, by = "term")

ordiplot(spe.rda.signif,
         scaling = 1,
         type = "text")

ordiplot(spe.rda.signif,
         scaling = 2,
         type = "text")

extrait le % expliqué par les 2 premiers axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

scores - ceux-ci sont des coordonnées dans l'espace RDA
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

Configuration du graphique

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

# Charger les données d'abondance des espèces d'acariens
data("mite")

# Charger les données environnementales
data("mite.env")

decostand()
rda()
ordiR2step()
anova.cca()
ordiplot()

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

ordiplot(mite.spe.rda.signif,
         scaling = 1,
         frame = F,
         main = "Cadrage 1")

ordiplot(mite.spe.rda.signif,
         scaling = 2,
         frame = F,
         main = "Cadrage 2")

# Divisez le tableau de données environnementales en deux:
# variables topographiques et chimiques
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Faire la RDA partielle
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
                       Condition(alt + pen + deb), # covariables ici
                       data = env.z)

summary(spe.partial.rda)

# Extraire le R2 ajusté du modèle
RsquareAdj(spe.partial.rda)$adj.r.squared

# Évaluer la significativité statistique du modèle
anova.cca(spe.partial.rda, step = 1000)

ordiplot(spe.partial.rda,
         scaling = 2,
         main = "Rivière Doubs - Cadrage 2")

rda()
summary()
RsquareAdj()
anova.cca() # voir l'argument 'by' dans ?anova.cca

mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Extraire les résultats
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")

spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # extraire résultats

plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # noms des matrices explicatives
     bg = c("seagreen3", "mediumpurple"), alpha = 80,
     digits = 2,
     cex = 1.5)

anova.cca(rda(spe.hel, env.chem))

anova.cca(rda(spe.hel, env.topo))

anova.cca(rda(spe.hel, env.chem, env.topo))

anova.cca(rda(spe.hel, env.topo, env.chem))

data("mite.pcnm")

ordiR2step()
varpart()
anova.cca(rda())
plot()

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

remotes::install_github("cran/mvpart")

library(mvpart)

# Enlever la variable “distance from source”
# Remove this chunk when LDA is back in. It is needed for the LDA (just switch next chunk to true).
env <- subset(env, select = -das)

# Enlever la variable “distance from source”
env <- subset(env, select = -das)

# Créer l'arbre de regression multivarié
# library(mvpart)
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick", # selection graphique intéractive
                    xval = nrow(spe.hel), # nombre de validations
                    xvmult = 100, # nombre de validations multiples
                    which = 4, # identifier les noeuds
                    legend = FALSE, margin = 0.01, cp = 0)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick",
                    xval = nrow(spe.hel),
                    xvmult = 100,
                    which = 4,
                    legend = FALSE, margin = 0.01, cp = 0, plot.add = FALSE)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick",
                    xval = nrow(spe.hel),
                    xvmult = 100,  cross-validations
                    which = 4,
                    legend = FALSE, margin = 0.01, cp = 0,
                    plot.add = FALSE)

doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "1se",
                    xval = nrow(spe.hel),
                    xvse = 1,
                    xvmult = 100,
                    which = 4,
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

# Calcul d'une valeur indval pour chaque espèce
doubs.mrt.indval <- indval(spe.hel, doubs.mrt$where)

# Extraire les espèces indicatrices à chaque noeud
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval <= 0.05)]

# Extraire leur valeur indval
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval <= 0.05)]

data("mite")
data("mite.env")

?mvpart() # argument 'xv'!
summary()

mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se",
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

# Calcul d'une valeur indicatrice pour chaque espèce
mite.mrt.indval <- indval(mite.spe.hel, mite.mrt$where)

# Extraire les espèces indicatrices à chaque noeud
mite.mrt.indval$maxcls[which(mite.mrt.indval$pval <= 0.05)]

# Extraire leur valeur indval
mite.mrt.indval$indcls[which(mite.mrt.indval$pval <= 0.05)]

# charger les données spatiales des sites Doubs:
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # assigner un chiffre par site
spa <- spa[-8,] # enlever le site #8

spa$group <- NA # créer colonne "group"
spa$group[which(spa$y < 82)] <- 1
spa$group[which(spa$y > 82 & spa$y < 156)] <- 2
spa$group[which(spa$y > 156)] <- 3

ggplot(data = spa) +
  geom_point(aes(x = x,
                 y = y,
                 col = as.factor(group)),
             size = 4) +
  labs(color = "Groupes",
       x = "Longitude",
       y = "Latitude") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # configuration
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# charger la libraire requise
library(MASS)

# faire la LDA
LDA <- lda(env, spa$group)

# prédire les groupes à partir de la LDA
lda.plotdf <- data.frame(group = spa$group, lda = predict(LDA)$x)

library(ggplot2)
# Visualiser les sites réorganisés à partir de la LDA
ggplot(lda.plotdf) +
  geom_point(aes(x = lda.LD1,
                 y = lda.LD2,
                 col = factor(group)),
             size = 4) +
  labs(color = "Groupes") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # configuration de la figure pour la rendre plus belle
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# classification des objets en fonction de la LDA
spe.class <- predict(LDA)$class

# probabilités que les objets appartiennent à chaque groupe a posteriori
spe.post <- predict(LDA)$posterior

# tableau des classifications a priori et prédites
(spe.table <- table(spa$group, spe.class))

# proportion de classification correcte
diag(prop.table(spe.table, 1))

# charger les nouvelles données
classify.me <- read.csv("data/classifyme.csv", header = TRUE)
# enlever das
classify.me <- subset(classify.me, select = -das)

# prédire le groupement des nouvelles données
predict.group <- predict(LDA, newdata = classify.me)

# prédire la classification pour chaque site
predict.group$class

data(mite.xy)

lda()
predict()
table()
diag()

# numéroter les sites
mite.xy$site <- 1:nrow(mite.xy)

# trouver une étendue égale de latitudes par groupe
(max(mite.xy[,2])-min(mite.xy[,2]))/4

# classifier les sites dans 4 groupes de latitude
mite.xy$group <- NA # nouvelle colonne "group"
mite.xy$group[which(mite.xy$y < 2.5)] <- 1
mite.xy$group[which(mite.xy$y >= 2.5 & mite.xy$y < 4.9)] <- 2
mite.xy$group[which(mite.xy$y >= 4.9 & mite.xy$y < 7.3)] <- 3
mite.xy$group[which(mite.xy$y >= 7.3)] <- 4

LDA.mite <- lda(mite.env[,1:2], mite.xy$group)

# classification des objects en fonction de la LDA
mite.class <- predict(LDA.mite)$class
# tableeau de classifications  (prior versus predicted)
(mite.table <- table(mite.xy$group, mite.class))
# proportion de classifications exactes
diag(prop.table(mite.table, 1))

# proportion of correct classification
diag(prop.table(mite.table, 1))
