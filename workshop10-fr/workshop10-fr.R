## ----setup, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  comment = "#",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = TRUE,
  fig.width = 6, fig.height = 6,
  fig.retina = 3,
  fig.retina = 3,
  fig.align = 'center'
)
options(repos=structure(c(CRAN="http://cran.r-project.org")))


## ----install_pkgs, echo = FALSE, results = "asis"-----------------------------
cat(
  qcbsRworkshops::first_slides(10, c('Hmisc', 'labdsv', 'MASS', 'vegan'),
    lang = "fr")
)


## ---- echo = TRUE-------------------------------------------------------------
# Assurez vous que les fichiers se trouvent dans votre répertoire de travail!
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <- spe[-8,] # Supprimer site 8 (pas d'espèces).


## ---- echo = TRUE-------------------------------------------------------------
env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # Supprimer site 8 (pas d'espèces).


## ---- echo = TRUE, results = 'hide'-------------------------------------------
names(spe) # voir les noms des colonnes (espèces)
dim(spe) # dimensions de la matrice
head(spe) # 5 premières lignes
str(spe) # structure interne de la matrice
summary(spe) # statistiques descriptives des objets (min, moyenne, max, etc.)


## ---- echo = FALSE------------------------------------------------------------
dim(spe)


## ---- echo = TRUE, results = 'hide', fig.width = 6, fig.height = 3.5----------
# Compter la fréquence d'espèces dans chaque classe d'abondance
ab <- table(unlist(spe))
# Visualiser cette distribution
barplot(ab, las = 1,
        xlab = "Abundance class", ylab = "Frequency",
        col = grey(5:0/5))


## -----------------------------------------------------------------------------
sum(spe == 0)


## -----------------------------------------------------------------------------
sum(spe == 0)/(nrow(spe)*ncol(spe))


## -----------------------------------------------------------------------------
# a fonction decostand() dans la libraire vegan nous facilite la tâche:
library(vegan)
spe.hel <- decostand(spe, method = "hellinger")


## ---- echo = TRUE, results = 'hide'-------------------------------------------
names(env) # noms des objets (variables environnementales)
dim(env) # dimensions de la matrice
head(env) # 5 premières lignes
str(env) # structure des objets
summary(env) # statistiques descriptives (min, moyenne, max, etc.)


## ---- echo = FALSE------------------------------------------------------------
names(env) # noms des objets (variables environnementales)
dim(env) # structure des objets


## ---- fig.height = 5, fig.width = 9-------------------------------------------
# On peut également détecter (visuellement) les colinéarités entres variables:
pairs(env)


## -----------------------------------------------------------------------------
# standardiser les données
env.z <- decostand(env, method = "standardize")

# centrer les données (moyenne ~ 0)
round(apply(env.z, 2, mean), 1)

# réduire les données (écart type = 1)
apply(env.z, 2, sd)


## -----------------------------------------------------------------------------
# On utilisera nos données explicatives standardisées
# Enlever la variable "distance from the source" (colinéarity avec plusieurs variables)
env.z <- subset(env.z, select = -das)


## -----------------------------------------------------------------------------
# Modèlise l'effect de tous les variables environnementales sur la composition en espèces des communautés
spe.rda <- rda(spe.hel ~ ., data = env.z)


## ---- eval = FALSE, results = 'hide'------------------------------------------
## summary(spe.rda, display = NULL)


## -----------------------------------------------------------------------------
# Sélection progressive de variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # modèle le plus simple
               scope = formula(spe.rda), # modèle "complet"
               direction = "forward",
               R2scope = TRUE, # limité par le R2 du modèle "complet"
               pstep = 1000,
               trace = FALSE) # mettre TRUE pour voir le processus du sélection!


## -----------------------------------------------------------------------------
fwd.sel$call


## -----------------------------------------------------------------------------
spe.rda.signif <- rda(spe.hel ~ alt + oxy + dbo, data = env.z)
RsquareAdj(spe.rda.signif)


## -----------------------------------------------------------------------------
anova.cca(spe.rda.signif, permutations = 1000)


## ---- results = 'hide'--------------------------------------------------------
anova.cca(spe.rda.signif, permutations = 1000, by = "axis")


## ---- fig.height = 6.5, fig.width = 6, strip.white = TRUE---------------------
ordiplot(spe.rda.signif,
         scaling = 1,
         type = "text")


## ---- fig.height = 6.5, fig.width = 6, strip.white = TRUE---------------------
ordiplot(spe.rda.signif,
         scaling = 2,
         type = "text")


## ---- echo = FALSE------------------------------------------------------------
## extrait le % expliqué par les 2 axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)
## scores => coordonnées
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)
## plot
plot(spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=paste0("RDA1 (", perc[1], "%)"), ylab = paste0("RDA2 (", perc[2], "%)"), xlim=c(-1,1), ylim=c(-1,1))
points(sc_si, pch=21, col="black", bg="steelblue", cex=1.2)
points(sc_sp, pch=22, col="black", bg = "#f2bd33", cex=1.2)
text(sc_sp, labels = rownames(sc_sp), col="black", cex=0.6)
arrows(0,0, sc_bp[,1], sc_bp[,2], col="#f04f6c", lwd = 3)
text(x = sc_bp[,1] -0.1, y = sc_bp[,2] - 0.03, labels=rownames(sc_bp), col="#f04f6c", cex=1, font = 2)


## -----------------------------------------------------------------------------
# Charger les données d'abondance des espèces d'acariens
data("mite")

# Charger les données environnementales
data("mite.env")


## ---- eval = FALSE------------------------------------------------------------
## decostand()
## rda()
## ordiR2step()
## anova.cca()
## ordiplot()


## -----------------------------------------------------------------------------
# Transformer les données d'abondances
mite.spe.hel <- decostand(mite, method = "hellinger")

# Standardiser les données environmentales quantiatives
mite.env$SubsDens <- decostand(mite.env$SubsDens, method = "standardize")
mite.env$WatrCont <- decostand(mite.env$WatrCont, method = "standardize")


## -----------------------------------------------------------------------------
# RDA avec tous les variables environnementales
mite.spe.rda <- rda(mite.spe.hel ~ ., data = mite.env)

# Sélection progressive des variables environnementales significatives
fwd.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.env),
                      scope = formula(mite.spe.rda),
                      direction = "forward",
                      R2scope = TRUE, pstep = 1000, trace = FALSE)
fwd.sel$call


## -----------------------------------------------------------------------------
# Refaire la RDA avec seulement les variables significatives
mite.spe.rda.signif <- rda(mite.spe.hel ~ WatrCont + Shrub +
                           Substrate + Topo + SubsDens,
                           data = mite.env)

# Calculer le R2 ajusté
RsquareAdj(mite.spe.rda.signif)$adj.r.squared



## -----------------------------------------------------------------------------
anova.cca(mite.spe.rda.signif, step = 1000)


## ---- fig.height = 6.5--------------------------------------------------------
ordiplot(mite.spe.rda.signif,
         scaling = 1,
         main = "Cadrage 1")


## ---- fig.height = 6.5--------------------------------------------------------
ordiplot(mite.spe.rda.signif,
         scaling = 2,
         main = "Cadrage 2")


## -----------------------------------------------------------------------------
# Divisez le tableau de données environnementales en deux:
# variables topographiques et chimiques
env.topo <- subset(env.z, select = c(alt, pen, deb))
env.chem <- subset(env.z, select = c(pH, dur, pho, nit, amm, oxy, dbo))

# Faire la RDA partielle
spe.partial.rda <- rda(spe.hel, env.chem, env.topo)


## ---- eval = FALSE------------------------------------------------------------
## spe.partial.rda <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
##                        Condition(alt + pen + deb),
##                        data = env.z)


## ---- eval = FALSE, results = 'hide'------------------------------------------
## # Extraire les résultats
## summary(spe.partial.rda, display = NULL)


## -----------------------------------------------------------------------------
RsquareAdj(spe.partial.rda)$adj.r.squared


## -----------------------------------------------------------------------------
anova.cca(spe.partial.rda, step = 1000)


## -----------------------------------------------------------------------------
ordiplot(spe.partial.rda, scaling = 2,
         main = "Doubs River partial RDA - Scaling 2")


## ---- eval = FALSE------------------------------------------------------------
## mite.spe.hel
## mite.env
## rda()
## summary()
## RsquareAdj()
## anova.cca()


## ---- results = 'hide'--------------------------------------------------------
mite.spe.subs <- rda(mite.spe.hel ~ Shrub + Topo
                     + Condition(SubsDens + WatrCont + Substrate),
                     data = mite.env)

# Extraire les résultats
summary(mite.spe.subs, display = NULL)


## -----------------------------------------------------------------------------
anova.cca(mite.spe.subs, step = 1000)


## -----------------------------------------------------------------------------
anova.cca(mite.spe.subs, step = 1000, by = "axis")


## -----------------------------------------------------------------------------
spe.part.all <- varpart(spe.hel, env.chem, env.topo)
spe.part.all$part # extraire résultats


## ---- strip.white = TRUE, fig.width = 6, fig.height = 6-----------------------
plot(spe.part.all,
     Xnames = c("Chem", "Topo"), # noms des matrices explicatives
     bg = c("seagreen3", "mediumpurple"), alpha = 80,
     digits = 2,
     cex = 1.5)


## -----------------------------------------------------------------------------
anova.cca(rda(spe.hel, env.chem))


## -----------------------------------------------------------------------------
anova.cca(rda(spe.hel, env.topo))


## -----------------------------------------------------------------------------
anova.cca(rda(spe.hel, env.chem, env.topo))


## -----------------------------------------------------------------------------
anova.cca(rda(spe.hel, env.topo, env.chem))


## -----------------------------------------------------------------------------
data("mite.pcnm")


## ---- eval = FALSE------------------------------------------------------------
## ordiR2step()
## varpart()
## anova.cca(rda())
## plot()


## -----------------------------------------------------------------------------
# Modèle RDA avec tous les variables spatiales
full.spat <- rda(mite.spe.hel ~ ., data = mite.pcnm)

# Sélection progressive des variables spatiales
spat.sel <- ordiR2step(rda(mite.spe.hel ~ 1, data = mite.pcnm),
               scope = formula(full.spat),
               R2scope = RsquareAdj(full.spat)$adj.r.squared,
               direction = "forward",
               trace = FALSE)
spat.sel$call


## -----------------------------------------------------------------------------
# Variables de substrat
mite.subs <- subset(mite.env, select = c(SubsDens, WatrCont))

# Variables spatiales significatives
mite.spat <- subset(mite.pcnm,
                    select = names(spat.sel$terminfo$ordered))
                    # pour rapidement accèder aux variables sélectionnées


## -----------------------------------------------------------------------------
mite.part <- varpart(mite.spe.hel, mite.subs, mite.spat)
mite.part$part$indfract # extraire résultats


## ---- results = 'hide'--------------------------------------------------------
# [a]: Substrat seulement
anova.cca(rda(mite.spe.hel, mite.subs, mite.spat))
# p = 0.001 ***

# [c]: Espace seulement
anova.cca(rda(mite.spe.hel, mite.spat, mite.subs))
# p = 0.001 ***


## ---- fig.height=6, fig.width=6-----------------------------------------------
plot(mite.part, digits = 2, cex = 1.5,
     bg = c("pink", "skyblue"), alpha = 90)


## ----mvpart_install-----------------------------------------------------------
remotes::install_github("cran/mvpart")
library(mvpart)


## ---- results = 'hide', fig.show = 'hide', eval = F---------------------------
## # Enlever la variable “distance from source”
## env <- subset(env, select = -das)
## 
## # Créer l'arbre de regression multivarié
## # library(mvpart)
## doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
##                     xv = "pick", # selection graphique intéractive
##                     xval = nrow(spe.hel), # nombre de validations
##                     xvmult = 100, # nombre de validations multiples
##                     which = 4, # identifier les noeuds
##                     legend = FALSE, margin = 0.01, cp = 0)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5.5, fig.align = 'center'----
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "pick",
                    xval = nrow(spe.hel),
                    xvmult = 100,
                    which = 4,
                    legend = FALSE, margin = 0.01, cp = 0, plot.add = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5.5, fig.align = 'center', eval =FALSE----
## doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
##                     xv = "pick",
##                     xval = nrow(spe.hel),
##                     xvmult = 100,  cross-validations
##                     which = 4,
##                     legend = FALSE, margin = 0.01, cp = 0,
##                     plot.add = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5.5, fig.width = 5.5-------
doubs.mrt <- mvpart(as.matrix(spe.hel) ~ ., data = env,
                    xv = "1se",
                    xval = nrow(spe.hel),
                    xvse = 1,
                    xvmult = 100,
                    which = 4,
                    legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 4.5, fig.width = 12--------
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 10, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 4.5, fig.width =8----------
mvpart(as.matrix(spe.hel) ~ ., data = env,
        xv = "none", # no cross-validation
        size = 4, # set tree size
        which = 4,
        legend = FALSE, margin = 0.01, cp = 0, prn = FALSE)


## -----------------------------------------------------------------------------
doubs.mrt$cptable


## -----------------------------------------------------------------------------
summary(doubs.mrt)


## ----MVPARTwrap_install-------------------------------------------------------
remotes::install_github("cran/MVPARTwrap")
library(MVPARTwrap)


## ---- results = 'hide'--------------------------------------------------------
# Créer un sommaire plus informatif et moins dense
doubs.mrt.wrap <- MRT(doubs.mrt, percent = 10, species = colnames(spe.hel))

# Voir le sommaire
summary(doubs.mrt.wrap)


## -----------------------------------------------------------------------------
summary(doubs.mrt.wrap)


## -----------------------------------------------------------------------------
library(labdsv)

# Calcul d'une valeur indval pour chaque espèce
doubs.mrt.indval <- indval(spe.hel, doubs.mrt$where)

# Extraire les espèces indicatrices à chaque noeud
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval <= 0.05)]

# Extraire leur valeur indval
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval <= 0.05)]


## -----------------------------------------------------------------------------
data("mite")
data("mite.env")


## ---- eval = FALSE------------------------------------------------------------
## ?mvpart() # argument 'xv'!
## ?MRT()
## summary()


## ---- results = 'hide', fig.height = 4.5, fig.width = 4.5---------------------
mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se",
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01, cp = 0,
                   prn = FALSE)


## ---- echo = FALSE, results = 'hide', fig.height = 5, fig.width = 5-----------
mite.mrt <- mvpart(as.matrix(mite.spe.hel) ~ ., data = mite.env,
                   xv = "1se", # choose smallest tree within 1 SE
                   xval = nrow(mite.spe.hel),
                   xvmult = 100,
                   which = 4, legend = FALSE, margin = 0.01,
                   cp = 0, prn = FALSE)


## ---- results = 'hide'--------------------------------------------------------
# Créer sommaire plus informatif
mite.mrt.wrap <- MRT(mite.mrt,
                     percent = 10,
                     species = colnames(mite.spe.hel))

# Voir sommaire (pour voir les espèces discriminantes)
summary(mite.mrt.wrap)


## -----------------------------------------------------------------------------
# charger les données spatiales des sites Doubs:
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # assigner un chiffre par site
spa <- spa[-8,] # enlever le site #8


## -----------------------------------------------------------------------------
spa$group <- NA # créer colonne "group"
spa$group[which(spa$y < 82)] <- 1
spa$group[which(spa$y > 82 & spa$y < 156)] <- 2
spa$group[which(spa$y > 156)] <- 3


## ---- fig.width = 5.5, fig.height = 5.5---------------------------------------
plot(spa$x, spa$y, col = spa$group, pch = 16, cex = 1.5)


## -----------------------------------------------------------------------------
# charger la libraire requise
library(MASS)

# faire la LDA
LDA <- lda(env, spa$group)


## -----------------------------------------------------------------------------
# classification des objets en fonction de la LDA
spe.class <- predict(LDA)$class

# probabilités que les objets appartiennent à chaque groupe a posteriori
spe.post <- predict(LDA)$posterior

# tableau des classifications a priori et prédites
(spe.table <- table(spa$group, spe.class))

# proportion de classification correcte
diag(prop.table(spe.table, 1))


## -----------------------------------------------------------------------------
# charger les nouvelles données
classify.me <- read.csv("data/classifyme.csv", header = TRUE)
# classify.me <- classify.me[,-1] # remove das variable

# prédire le groupement des nouvelles données
predict.group <- predict(LDA, newdata = classify.me)

# donner la classification pour chaque site
predict.group$class


## -----------------------------------------------------------------------------
data(mite.xy)


## ---- eval = FALSE------------------------------------------------------------
## lda()
## predict()
## table()
## diag()


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
LDA.mite <- lda(mite.env[,1:2], mite.xy$group)


## -----------------------------------------------------------------------------
# classification des objects en fonction de la LDA
mite.class <- predict(LDA.mite)$class
# tableeau de classifications  (prior versus predicted)
(mite.table <- table(mite.xy$group, mite.class))
# proportion de classifications exactes
diag(prop.table(mite.table, 1))

