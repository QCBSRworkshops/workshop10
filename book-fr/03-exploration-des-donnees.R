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

## str(env)
## summary(env)

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
