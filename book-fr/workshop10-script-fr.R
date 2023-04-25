##Section: 01-preparation-pour-l-atelier.R 

###Avis ###
#                                                                            #
#Ceci est un script généré automatiquement basé sur les morceaux de code du  #
#livre pour cet atelier.                                                     #
#                                                                            #
#Il est minimalement annoté pour permettre aux participants de fournir leurs #
#commentaires : une pratique que nous encourageons vivement.                 #
#                                                                            #
#Notez que les solutions aux défis sont également incluses dans ce script.   #
#Lorsque vous résolvez les défis par vous-méme, essayez de ne pas parcourir  #
#le code et de regarder les solutions.                                       #
#                                                                            #
#Bon codage !                                                               #


# Installez les paquets requis
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("ggplot2")

# installez mvpart de l'archive
install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Chargez les pacquets requis
library(labdsv)
library(vegan)
library(MASS)
library(mvpart)
library(ggplot2)


##Section: 02-introduction-fr.R 




##Section: 03-exploration-des-donnees.R 

# Assurez vous que les fichiers se trouvent dans votre répertoire de travail!
# Si R ne trouve pas le jeu de données, définissez votre répertoire de travail avec setwd()
# au dossier dans lequel vos données sont sauvegardées (par exemple setwd("~/Desktop/workshop10"))
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # Supprimer site 8 (pas d'espèces).

# Matrice d'abondances d'espèces de poissons: “DoubsSpe.csv”
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # Supprimer site 8 (pas d'espèces).
# Attention! Exécuter cette ligne une seule fois. 

# Matrice de données environnementales: “DoubsEnv.csv”
env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # Supprimer le site 8 puisqu'on l'a supprimé de la matrice d'abondance. # N'exécuter qu'une seule fois.

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

# Appliquer la transformation de Hellinger pour corriger le problème de double zéro
spe.hel <- decostand(spe, method = "hellinger")

names(env)
dim(env)
head(env)

str(env)
summary(env)

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


##Section: 04-analyses-canoniques.R 




##Section: 05-analyses-de-redondance.R 

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

# Défi 1: Effectuer une RDA pour modèliser les effect des variables environnementales sur l’abondance des espèces d'acariens.

# Charger les données d'abondance des espèces d'acariens
data("mite")

# Charger les données environnementales
data("mite.env")

decostand()
rda()
ordiR2step()
anova.cca()
ordiplot()

# Défi 1: Solution! Spoilers ci-dessous!!

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


##Section: 06-analyses-partielles-de-redondance.R 

knitr::include_graphics("images/PartialRDA.png")

# Divisez le tableau de données environnementales en deux:
# variables topographiques et chimiques
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Faire la RDA partielle
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)

# Syntaxe alternative
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

# Défi 2:
# Effectuez une RDA partielle de l’abondance des espèces de mites (`mite.spe.hel`) en fonction des variables environnementales, tenant compte de l’effet du substrat (`SubsDens`, `WaterCont` and `Substrate`).
# * Quel pourcentage de variance est expliqué par les variables environnementales?
# * Le modèle est-il significatif?
# * Quels sont les axes significatifs?

rda()
summary()
RsquareAdj()
anova.cca() # voir l'argument 'by' dans ?anova.cca

# Défi 2: Solution! Spoilers ci-dessous!!

mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Extraire les résultats
summary(mite.spe.subs)

RsquareAdj(mite.spe.subs)$adj.r.squared

anova.cca(mite.spe.subs, step = 1000)

anova.cca(mite.spe.subs, step = 1000, by = "axis")


##Section: 07-partition-de-la-variance.R 

# Partitionner la variation de la composition des espèces de poissons 
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # access results!

# Visualiser les résultats avec un diagramme Venn
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # noms des matrices explicatives
     bg = c("seagreen3", "mediumpurple"), alpha = 80,
     digits = 2,
     cex = 1.5)

# Tester la significativité

# [a+b] Chimie sans tenir compte de topographie
anova.cca(rda(spe.hel, env.chem))

# [b+c] Topographie sans tenir compte de chimie
anova.cca(rda(spe.hel, env.topo))

# [a] Chimie seulement
anova.cca(rda(spe.hel, env.chem, env.topo))

# [c] Topographie seulement
anova.cca(rda(spe.hel, env.chem, env.topo))

# Défi 3

# Partitionnez la variation de l’abondance des espèces de mites entres des variables de substrat (`SubsDens`, `WatrCont`) et des variables spatiales significatives.
# * Quelle est la proportion de variance expliquée par le substrat? par l'espace?
# * Quelles sont les fractions significatives?
# * Diagramme Venn des résultats!

data("mite.pcnm")

ordiR2step()
varpart()
anova.cca(rda())
plot()

# Défi 3 - Solution! Spoilers ci-dessous!

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


##Section: 08-arbre-de-regression-multivarie.R 

# Enlever la variable “distance from source”
env <- subset(env, select = -das)


##Section: 09-analyse-de-discrimination.R 

# charger les données spatiales des sites Doubs:
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # assigner un chiffre par site
spa <- spa[-8,] # enlever le site #8

# classification a priori
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

# faire la LDA
LDA <- lda(env, spa$group)

# prédire les groupes à partir de la LDA
lda.plotdf <- data.frame(group = spa$group, lda = predict(LDA)$x)

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

# Défi 5

# Créez quatre groupes de latitude avec des étendues égales à partir des données `mite.xy`. Ensuite, faites une LDA sur les données environnementales `mite.env` des acariens (`SubsDens` et `WatrCont`). **Quelle proportion de sites ont été classifiés correctement au groupe 1? Au groupe 2?**

data(mite.xy)

lda()
predict()
table()
diag()

# Défi 5: Solution! Spoilers ci-dessous!

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


##Section: 10-considerations-finales.R 




##Section: 11-references-fr.R 




