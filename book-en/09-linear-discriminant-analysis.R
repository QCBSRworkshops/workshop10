# load spatial data for Doubs sites
spa <- read.csv("data/doubsspa.csv", row.names = 1)
spa$site <- 1:nrow(spa) # add site numbers
spa <- spa[-8,] # remove site 8

# group sites based on latitude
spa$group <- NA # create "group" column
spa$group[which(spa$y < 82)] <- 1
spa$group[which(spa$y > 82 & spa$y < 156)] <- 2
spa$group[which(spa$y > 156)] <- 3

ggplot(data = spa) +
  geom_point(aes(x = x, 
                 y = y, 
                 col = as.factor(group)), 
             size = 4) +
  labs(color = "Groups", 
       x = "Longitude", 
       y = "Latitude") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # formatting the plot to make it pretty
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# run the LDA grouping sites into latitude groups based on env data
LDA <- lda(env, spa$group)

# predict the groupings
lda.plotdf <- data.frame(group = spa$group, lda = predict(LDA)$x)

# Plot the newly reorganised sites according to the LDA
ggplot(lda.plotdf) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 col = factor(group)), 
             size = 4) +
  labs(color = "Groups") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + # formatting the plot to make it pretty
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# Classification of the objects based on the LDA
spe.class <- predict(LDA)$class

# Posterior probabilities that the objects belong to those groups
spe.post <- predict(LDA)$posterior

# Table of prior vs. predicted classifications
(spe.table <- table(spa$group, spe.class))

# Proportion of corrected classification
diag(prop.table(spe.table, 1))

# Load the new site data
classify.me <- read.csv("data/classifyme.csv", header = TRUE)
# remove das
classify.me <- subset(classify.me, select = -das)
  
# Predict grouping of new sites
predict.group <- predict(LDA, newdata = classify.me)

# View site classification
predict.group$class

data(mite.xy)

## lda()
## predict()
## table()
## diag()

# assign numbers to sites
mite.xy$site <- 1:nrow(mite.xy)

# figure out the spacing to make 4 equally spaced latitude groups 
(max(mite.xy[,2])-min(mite.xy[,2]))/4

# use this to group sites into 4 latitude groups
mite.xy$group <- NA # create "group" column
mite.xy$group[which(mite.xy$y < 2.5)] <- 1
mite.xy$group[which(mite.xy$y >= 2.5 & mite.xy$y < 4.9)] <- 2
mite.xy$group[which(mite.xy$y >= 4.9 & mite.xy$y < 7.3)] <- 3
mite.xy$group[which(mite.xy$y >= 7.3)] <- 4

LDA.mite <- lda(mite.env[,1:2], mite.xy$group)

# group sites based on the LDA
mite.class <- predict(LDA.mite)$class

# get the table of prior versus predicted classifications
(mite.table <- table(mite.xy$group, mite.class))

# proportion of correct classification
diag(prop.table(mite.table, 1))
