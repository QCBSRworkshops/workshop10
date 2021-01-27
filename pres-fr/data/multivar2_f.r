# Nom: Amanda Winegardner, Xavier Giroux-Bougard, Bérenger Bourgeois, Emmanuelle Chrétien and Monica Granados 
# Date: Mars 2016
# Description: Analyses multivariées avancées
# Données:"DoubsSpe.csv", "DoubsEnv.csv"
# Notes: Matériel pour scripts R obtenu (en partie) de 
#        Borcard D., Gillet F. et Legendre P., 2011. Numerical Ecology avec R. Springer.
#***************************************************************************************#


# 0. Charger les données ####

rm(list=ls())
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("plyr")
install.packages("rJava")
install.packages("venneuler")

# Pour les librairies suivantes, chargez le fichier fourni sur la page wiki de l'atelier.
# Pour ce faire, cliquez sur l'onglet Packages" de la section en bas à droite de R Studio
# Cliquez sur Install Package
# Choisissez d'installer depuis Package Archive file et chargez ces fichiers
install.packages("mvpart") 
install.packages("MVPARTwrap") 
install.packages("rdaTest") 

# Si ça ne fonctionne pas, vous pouvez les installer depuis github.
# Pour ce faire, vous devez installer la librairie "devtools" et rouler ces quelques lignes de code
install.packages("devtools")
library(devtools)

# Installer les librairies depuis github
devtools::install_github("cran/mvpart")
devtools::install_github("cran/MVPARTwrap")
devtools::install_github("philippec/fonctions_R_git/rdaTest")

# Charger les librairies
library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(MASS)
library(plyr)

# Fichier d'abondance des espèces de poissons: "DoubsSpe.csv"
spe <- read.csv(file.choose(), row.names=1)
spe <- spe[-8,] # le site 8 ne contient aucune espèce et doit être enlevé

# Fichier de variables environnementales: "DoubsEnv.csv"
env <- read.csv(file.choose(), row.names=1)
env <- env[-8,]

#***************************************************************************************#

# 1.Exploration des données ####

## 1.1 Abondances d'espèces ####
names(spe)
dim(spe)
str(spe)
head(spe)
summary(spe) 

### Distribution des espèces, toutes espèces confondues
(ab <- table(unlist(spe)))
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

### Nombre d'absences
sum(spe==0)

### Proportion de zéros dans le fichier
sum(spe==0)/(nrow(spe)*ncol(spe))

## Appliquer une transformation appropriées aux données d'abondances d'espèces

?decostand # cette fonction rassemble des méthodes standards pour des données d'abondances d'espèces (voir l'argument "method")

# Transformation de Hellinger
spe.hel<-decostand(spe, method="hellinger")

## 1.2 Variables environnementales  ####
names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Bivariate Plots of the Environmental Data" )

# Centrer-réduire les 11 variables environnementales
env.z <- decostand(env, method="standardize")
apply(env.z, 2, mean) # centrage des données (moyennes~0)
apply(env.z, 2, sd)   # réduction des données (écarts-types=1)

# Note: les variables explicatives sont exprimées en unités différentes
#       et doivent donc être centrées-réduites avant la création de matrice de distances et avant leur 
#       utilisation dans des ordinations
#       Ce n'est cependant pas une étape nécessaire pour l'arbre de régression.


#***************************************************************************************#
# 2. Analyses canoniques
## 2.1. Analyse canonique de redondances (ACR ou RDA) ####

### Préparer les données
env.z <- subset(env.z, select = -das) # enlever la variable "distance from the source"

### Faire la RDA en utilisant toutes les variables environnementales
?rda
spe.rda <- rda(spe.hel~., data=env.z) # le point "." signifie que toutes les variables de env.z seront incluses

### Extraire les résultats
summary(spe.rda, display=NULL)

### Sélectionner les variables explicatives significatrices par sélection progressive
?ordiR2step
ordiR2step(rda(spe.hel~1, data=env.z), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
env.signif <- subset(env.z, select = c("alt", "oxy", "dbo"))

### Refaire la RDA en utilisant seulement les variables explicatives significatrices
spe.rda.signif <- rda(spe.hel~., data=env.signif)

### Extraire les résultats
summary(spe.rda.signif, display=NULL)

### Calculer le R2 ajusté de la RDA
(R2adj <- RsquareAdj(spe.rda.signif)$adj.r.squared)

### Tester la significativité du modèle et des axes
?anova.cca
anova(spe.rda.signif, step=1000)
anova(spe.rda.signif, step=1000, by="axis")

### Représentation graphique des résultats
#### graphiques rapides 
# cadrage 1
windows()
plot(spe.rda.signif, scaling=1, main="Triplot RDA (scaling 1)")
# cadrage 2
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA (scaling 2)")

#### graphiques avancés - cadrage 1
windows()
plot(spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spe.rda.signif, display="species", choices=c(1), scaling=1),
       scores(spe.rda.signif, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(spe.rda.signif, display="species", choices=c(1), scaling=1),
     scores(spe.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spe.rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spe.rda.signif, display="bp", choices=c(1), scaling=1),
       scores(spe.rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(spe.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(spe.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(spe.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#### graphiques avancés - cadrage 2
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2,
       scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(spe.rda.signif, display="species", choices=c(1), scaling=2)*2.1,
     scores(spe.rda.signif, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(spe.rda.signif, display="species", scaling=2)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spe.rda.signif, display="bp", choices=c(1), scaling=2),
       scores(spe.rda.signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(spe.rda.signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(spe.rda.signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# Défi 1 ####
# Faire une RDA sur les données d'acariens
data(mite) #charger les données depuis la librairie vegan
mite.spe <- mite
data(mite.env) 
# Quelles sont les variables explicatives significatrices?  
# Quelle proportion de la variance est expliquée par les variables explicatives?
# Quels sont les axes significatifs?
# Quels groupes de sites pouvez-vous identifier? 
# Quelles espèces sont liées à quels groupes de sites?
# (réponse ci-bas)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
mite.spe.hel<-decostand(mite.spe,method="hell")
mite.spe.rda <- rda(mite.spe.hel~., data=mite.env)
ordiR2step(rda(mite.spe.hel~1, data=mite.env), 
           scope= formula(mite.spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
which(colnames(mite.env)=="WatrCont")
which(colnames(mite.env)=="Shrub")
which(colnames(mite.env)=="Substrate")
which(colnames(mite.env)=="Topo")
which(colnames(mite.env)=="SubsDens")
(mite.env.signif=mite.env[,c(2,4,3,5, 1)])
mite.env.signif <- subset(mite.env, 
                          select = c("WatrCont", "Shrub", "Substrate", "Topo", "SubsDens"))
mite.spe.rda.signif <- rda(mite.spe~., data=mite.env.signif) 
summary(mite.spe.rda.signif, display=NULL)
(R2adj <- RsquareAdj(mite.spe.rda.signif)$adj.r.squared)
anova(mite.spe.rda.signif, step=1000)
anova(mite.spe.rda.signif, step=1000, by="axis")          
windows()
plot(mite.spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.rda.signif, display="species", choices=c(1), scaling=1),
     scores(mite.spe.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1),
       scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

## 2.2  RDA partielle ####
## Effet de la chimie de l'eau sur les abondances d'espèces en tenant compte de la physiographie
### Diviser les données environnementales en deux fichiers
envtopo <- env.z[, c(1:3)] # Physiographie : tableau 1
names(envtopo)
envchem <- env.z[, c(4:10)] # Physico-chimie de l'eau : tableau 2
names(envchem)

### Faire la RDA partielle
spechem.physio <- rda(spe.hel, envchem, envtopo)
# ou : 
spechem.physio2 <- rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo +
                         Condition(alt + pen + deb), data=env.z)

### Extraire les résultats
summary(spechem.physio, display=NULL)
summary(spechem.physio2, display=NULL)

### Calculer le R2 ajusté de la RDA partielle
(R2adj <- RsquareAdj(spechem.physio)$adj.r.squared)

### Tester la significativité de la RDA partielle et des axes
anova(spechem.physio, step=1000)
anova(spechem.physio2, step=1000, by="axis")

### Construire les triplots
#### Cadrage 1
windows(title="Partial RDA scaling 1")
plot(spechem.physio, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spechem.physio, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spechem.physio, display="species", choices=c(1), scaling=1),
       scores(spechem.physio, display="species", choices=c(2), scaling=1),
       col="grey",length=0)
text(scores(spechem.physio, display="species", choices=c(1), scaling=1),
     scores(spechem.physio, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spechem.physio, display="species", scaling=1)),
     col="grey", cex=0.8)  
arrows(0,0,
       scores(spechem.physio, display="bp", choices=c(1), scaling=1),
       scores(spechem.physio, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(spechem.physio, display="bp", choices=c(1), scaling=1)+0.05,
     scores(spechem.physio, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(spechem.physio, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#### Cadrage 2
windows(title="Partial RDA scaling 2")
plot(spechem.physio, scaling=2, main="Triplot partial RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spechem.physio, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spechem.physio, display="species", choices=c(1), scaling=2),
       scores(spechem.physio, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(spechem.physio, display="species", choices=c(1), scaling=2),
     scores(spechem.physio, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(spechem.physio, display="species", scaling=2)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(spechem.physio, display="bp", choices=c(1), scaling=2),
       scores(spechem.physio, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(spechem.physio, display="bp", choices=c(1), scaling=2)+0.05,
     scores(spechem.physio, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(spechem.physio, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# Défi 2 ####
# Faire la RDA partielle sur les données d'acariens 
# en contrôlant pour les variables de substrat (SubsDens, WaterCont and Substrate)
# Est-ce que le modèle est significatif ?
# Quels sont les axes significatifs ?
# Interpréter le triplot
# (réponse ci-bas)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
mite.spe.subs<-rda(mite.spe.hel ~ Shrub + Topo
                  + Condition(SubsDens + WatrCont + Substrate), data=mite.env)
summary(mite.spe.subs, display=NULL)
(R2adj <- RsquareAdj(mite.spe.subs)$adj.r.squared)
anova(mite.spe.subs, step=1000)
anova(mite.spe.subs, step=1000, by="axis")

#### Triplot cadrage 1
windows(title="Partial RDA scaling 1")
plot(mite.spe.subs, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1)
arrows(0,0,
       scores(mite.spe.subs, display="species", choices=c(1), scaling=1),
       scores(mite.spe.subs, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(mite.spe.subs, display="species", choices=c(1), scaling=1),
     scores(mite.spe.subs, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.subs, display="species", scaling=1)),
     col="black", cex=0.8)  
arrows(0,0,
       scores(mite.spe.subs, display="bp", choices=c(1), scaling=1),
       scores(mite.spe.subs, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(mite.spe.subs, display="bp", choices=c(1), scaling=1)+0.05,
     scores(mite.spe.subs, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(mite.spe.subs, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#### Triplot cadrage 2
windows(title="Partial RDA scaling 2")
plot(mite.spe.subs, scaling=2, main="Triplot partial RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
       col="grey",length=0)
text(scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
     scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(mite.spe.subs, display="species", scaling=2)),
     col="grey", cex=0.8)  
arrows(0,0,
       scores(mite.spe.subs, display="bp", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(mite.spe.subs, display="bp", choices=c(1), scaling=2)+0.05,
     scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  

# 2.3. Partitionnement de la variation par RDA partielle ####
?varpart
vegandocs("partitioning")
# vous pouvez aussi accéder à la documentation de cette manière, tel que vu dans le premier atelier
vignette("partitioning", package="vegan")

## Partitionnement de la variation avec toutes les variables explicatives
spe.part.all <- varpart(spe.hel, envchem, envtopo)
spe.part.all
windows(title="Variation partitioning - all variables")
plot(spe.part.all, digits=2)  

## Partitionnement de la variation avec les variables explicatives significatrices seulement
### Sélection progressive des variables dans chaque tableau
spe.chem <- rda(spe.hel~., data=envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envchem), 
           scope= formula(spe.chem), direction= "forward", R2scope=TRUE, pstep=1000)
names(envchem)
(envchem.pars <- envchem[, c( 4, 6, 7 )])
spe.topo <- rda(spe.hel~., data=envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envtopo), 
           scope= formula(spe.topo), direction= "forward", R2scope=TRUE, pstep=1000)
names(envtopo)
envtopo.pars <- envtopo[, c(1,2)]

### Partitionnement de la variation
(spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars))
windows(title="Variation partitioning - parsimonious subsets")
plot(spe.part, digits=2)

# Tester les fractions
anova(rda(spe.hel, envchem.pars), step=1000) # Test of fractions [a+b]
anova(rda(spe.hel, envtopo.pars), step=1000) # Test of fractions [b+c]
env.pars <- cbind(envchem.pars, envtopo.pars) 
anova(rda(spe.hel, env.pars), step=1000) # Test of fractions [a+b+c]
anova(rda(spe.hel, envchem.pars, envtopo.pars), step=1000) # Test of fraction [a]
anova(rda(spe.hel, envtopo.pars, envchem.pars), step=1000) # Test of fraction [c]

# Pour faire un diagramme de Venn plus cool
install.packages("venneuler")
library(venneuler)
v <- venneuler(c(A=0.241,"A&B"= 0.233, B=0.112))
v$labels <- c("Chem\n24%","Topo \n11%") 
v$colors <- c(.78,.6) 
plot(v)


# Défi 3 ####
# Faire le partitionnement de variation sur les données d'acariens
# avec un premier tableau avec les variables significatives de substrat (SubsDens, WaterCont and Substrate)
# et un second tableau avec les autres variables significatives (Shrud and Topo)
# Quelle proportion de variance est expliqué par chaque fichier de données?  
# Quelles sont les fractions significatives?
# (réponse ci-bas)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
str(mite.env)
(mite.subs<-mite.env[,c(1,2,3)])
(mite.other<-mite.env[,c(4,5)])
rda.mite.subs <- rda(mite.spe.hel~., data=mite.subs)
R2a.all.chem <- RsquareAdj(rda.mite.subs)$adj.r.squared
ordiR2step(rda(mite.spe.hel~1, data=mite.subs), 
           scope= formula(rda.mite.subs), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.subs)
(mite.subs.pars <- mite.subs[, c(2, 3)])
rda.mite.other <- rda(mite.spe.hel~., data=mite.other)
R2a.all.chem <- RsquareAdj(rda.mite.other)$adj.r.squared
ordiR2step(rda(mite.spe.hel~1, data=mite.other), 
           scope= formula(rda.mite.other), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.other)
(mite.other.pars <- mite.other[, c(1,2)])

### Partitionnement de la variation
(mite.spe.part <- varpart(mite.spe.hel, ~WatrCont+Substrate, ~Shrub+Topo,
                          data=mite.env))
windows(title="Variation partitioning - parsimonious subsets")
plot(mite.spe.part, digits=2)

# Tester les fractions
anova(rda(mite.spe.hel~ WatrCont+Substrate, data=mite.env), step=1000) # Test of fractions [a+b]
anova(rda(mite.spe.hel~Shrub+Topo, data=mite.env), step=1000) # Test of fractions [b+c]
(env.pars <- cbind(mite.env[,c(2,3,4,5)])) 
anova(rda(mite.spe.hel~ WatrCont+Substrate+Shrub+Topo, data=env.pars), step=1000) # Test of fractions [a+b+c]
anova(rda(mite.spe.hel~WatrCont+Substrate + Condition(Shrub+Topo), data=env.pars), step=1000) # Test of fraction [a]
anova(rda(mite.spe.hel~Shrub+Topo+ Condition(WatrCont+Substrate ), data=env.pars), step=1000) # Test of fraction [c]



# 3. MRT: Arbre de régression multivarié ####

# Préparer les données
env <- subset(env, select = -das) # Enlever "distance from the source" 

# Créer l'arbre de régression
# Utiliser le critère de l'arbre de taille minimale à 1 erreur-type de la CVRE minimale
doubs.mrt<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4) # cliquer sur le point orange sur le graphique
summary(doubs.mrt)


# Utiliser le critère de la CVRE minimale
doubs.mrt.cvre<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4) # cliquer sur le point rouge sur le graphique

# En choisissant un arbre de 4 branches
doubs.mrt.4<-mvpart(data.matrix(spe.hel)~.,env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(spe.hel),xvmult=100,which=4) # cliquer sur le point bleu sous le chiffre 4 sur le graphique
summary(doubs.mrt.4)


# Trouver des espèces indicatrices avec les résultats
doubs.mrt.wrap<-MRT(doubs.mrt,percent=10,species=colnames(spe.hel))
summary(doubs.mrt.wrap)

doubs.mrt.indval<-indval(spe.hel,doubs.mrt$where)
doubs.mrt.indval$pval

doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)]
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]


# Défi 5 ####
# Faire le MRT avec les données d'acariens
# Quelle est la taille idéale selon le critère de Breiman et al. 1984?
# Quelle proportion de la variation est expliquée par l'arbre?  
# Quelles sont les espèces discriminantes?
# (réponse ci-bas)
#
#
#
#
#
#
#
#
#
mite.mrt<-mvpart(data.matrix(mite.spe.hel)~.,mite.env,legend=FALSE,margin=0.01,cp=0,xv="pick",
                  xval=nrow(mite.spe.hel),xvmult=100,which=4)
summary(mite.mrt)

mite.mrt.wrap<-MRT(mite.mrt,percent=10,species=colnames(mite.spe.hel))
summary(mite.mrt.wrap)

mite.mrt.indval<-indval(mite.spe.hel,mite.mrt$where)
mite.mrt.indval$pval

mite.mrt.indval$maxcls[which(mite.mrt.indval$pval<=0.05)]
mite.mrt.indval$indcls[which(mite.mrt.indval$pval<=0.05)]

# 4. LDA: Analyse discriminante linéaire ####

#Chargez les données spatiales pour déterminer les groupes 
spa <- read.csv ('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpa.csv', row.names = 1)
spa <- spa[-8,] # enlever site 8 car on l'a enlevé des deux autres jeux de données Doubs

#Visualiser les données spatiales
View (spa)

#Ajouter les numéros de sites
numbers<-(1:30)
numbers<-numbers[!numbers%in%8] 
spa$site<-numbers

#Faire des groupes en fonction de la latitude y<82=group1, 82<y<156=group2, y>156=group3
spa.group<-ddply(.data=spa, .variables=.(x, y, site), .fun= summarise, group = if(y <= 82) 1 else if (y <= 156) 2 else 3)

#Ordonner les données par site
spa.group<-spa.group[with(spa.group, order(site)), ]

#faire la LDA
LDA<-lda(env,spa.group[,4])

#classification des objets en fonction de la LDA
spe.class <- predict(LDA)$class

#probabilités que les objets appartiennent à chaque groupe a posteriori
spe.post <- predict(LDA)$posterior

#tableau des classifications a priori et prédites
spe.table <- table(spa.group[,4], spe.class)

#proportion de classification correcte
diag(prop.table(spe.table, 1))

#prédire la classification des nouvelles données 
#harger les nouvelles données 
classify.me<-read.csv("classifyme.csv", header = T)

#prédire le groupement des nouvelles données
predict.group<-predict(LDA, newdata=classify.me)

#donner la classification pour chaque site
group.new<-predict.group$class

# Défi 6 ####
# Faire une LDA sur les données environnementales des acariens (deux premières variables) 
# en se basant sur 4 groupes de latitude à créer à partir des données mite.xy. 
data(mite.xy)
# Quelle proportion de sites ont été classifiés correctement au groupe 1? Au groupe 2?
# (réponse ci-bas)
#
#
#
#
#
#
#
#
#
mite.xy$site<-seq(1:70)
(max(mite.xy[,2])-min(mite.xy[,2]))/4
mite.xy.group<-ddply(.data=mite.xy, .variables=.(x, y, site), .fun= summarise, group = if(y <= 2.5) 1 else if (y <= 4.9) 2 else if (y <= 7.3) 3 else 4)
mite.xy.group<-mite.xy.group[with(mite.xy.group, order(site)), ]
LDA.mite<-lda(mite.env[,1:2],mite.xy.group[,4])
mite.class <- predict(LDA.mite)$class
mite.post <- predict(LDA.mite)$posterior
mite.table <- table(mite.xy.group[,4], mite.class)
diag(prop.table(mite.table, 1))

# 5. autres méthodes d'ordination utiles####

?cca # analyse canonique des correspondances (CCA)

?CCorA # analyse canonique des corrélations

help(coinertia, package=ade4) # analyse de coinertie

help(mfa, package=ade4) # analyse factorielle multiple

https://r-forge.r-project.org/R/?group_id=195  # Les packages AEM et PCNM permettent d’effectuer diverses formes d’analyses spatiales

