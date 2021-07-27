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

## # Défi 5
## 
## # Créez quatre groupes de latitude avec des étendues égales à partir des données `mite.xy`. Ensuite, faites une LDA sur les données environnementales `mite.env` des acariens (`SubsDens` et `WatrCont`). **Quelle proportion de sites ont été classifiés correctement au groupe 1? Au groupe 2?**

data(mite.xy)

## lda()
## predict()
## table()
## diag()

## # Défi 5: Solution! Spoilers ci-dessous!

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
