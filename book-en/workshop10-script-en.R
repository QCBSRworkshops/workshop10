##Section: 01-preparing-for-the-workshop.R 

install.packages("vegan")
install.packages("mvpart")
install.packages("labdsv")
install.packages("plyr")
install.packages("MASS")

# For the two following packages, upload the file provided on the wiki page.
# To do so, go to Packages tab on the bottom right panel of R Studio
# Click on Install Packages
# Choose to install from Package Archive file and upload these two files
install.packages("MVPARTwrap")
install.packages("rdaTest")

library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(plyr)
library(MASS)


##Section: 02-introduction.R 




##Section: 03-data-exploration.R 

#Species community data frame (fish abundance): “DoubsSpe.csv”
spe<- read.csv(file.choose(), row.names=1)
spe<- spe[-8,] #Site number 8 contains no species and so row 8 (site 8) is removed. Be careful to
#only run this command line once as you are overwriting "spe" each time.

#Environmental data frame: “DoubsEnv.csv”
env<- read.csv(file.choose(), row.names=1)
env<- env[-8,] #Remove corresponding abiotic data for site 8 (since removed from fish data).
#Again, be careful to only run the last line once.

names(spe) #see names of columns in spe
dim(spe) #dimensions of spe; number of columns and rows
str(spe) #displays internal structure of objects
head(spe) #first few rows of the data frame
summary(spe) #summary statistics for each column; min value, median value, max value, mean value etc.

#Species distribution
(ab <- table(unlist(spe))) #note that when you put an entire line of code in brackets like this, the output for that operation is displayed right away in the R console

barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

sum(spe==0)

sum(spe==0)/(nrow(spe)*ncol(spe))

spe.hel <- decostand(spe, method="hellinger") # you can also use method="hell"

names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Bivariate Plots of the Environmental Data" )

env.z <- decostand(env, method="standardize")
apply(env.z, 2, mean) # the data are now centered (means~0)
apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)


##Section: 04-canonical-analysis.R 




##Section: 05-redundancy-analysis.R 

#Preparing the data prior to RDA
env.z <- subset(env.z, select = -das) # remove the "distance from the source" variable

#Running the RDA
?rda
spe.rda <- rda(spe.hel~., data=env.z)

#Extract the results
summary(spe.rda, display=NULL)

#The results are called using summary:
summary(spe.rda, display=NULL) #display = NULL optional

#Select the significant explanatory variables by forward selection
?ordiR2step
ordiR2step(rda(spe.hel~1, data=env.z), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
env.signif <- subset(env.z, select = c("alt", "oxy", "dbo"))

spe.rda.signif <- rda(spe.hel~., data=env.signif)
summary(spe.rda.signif, display=NULL)

(R2adj <- RsquareAdj(spe.rda.signif)$adj.r.squared)
#Here the strength of the relationship between X and Y corrected for the number of X variables is 0.54.

?anova.cca
anova.cca(spe.rda.signif, step=1000)
anova.cca(spe.rda.signif, step=1000, by="axis")
#In this case, the RDA model is highly significant (p=0.001) as well as all three canonical axes.

#Quick plots scaling 1 and 2
windows()
plot(spe.rda.signif, scaling=1, main="Triplot RDA (scaling 1)")
windows()
plot(spe.rda.signif, scaling=2, main="Triplot RDA (scaling 2)")

#Advanced plots scaling 1
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

#Advanced plots scaling 2
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

#Load the mite species and environmental data from vegan package
data(mite)
mite.spe<-mite
mite.spe.hel <- decostand(mite.spe, method="hellinger")

data(mite.env)

#Initial RDA with ALL of the environmental data
mite.spe.rda<-rda(mite.spe.hel~., data=mite.env)

#Select significant environmental variables
ordiR2step(rda(mite.spe.hel~1, data=mite.env),
           scope= formula(mite.spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)

#Create a new dataframe with only the significant variables that you identified above
mite.env.signif <- subset(mite.env,
                          select = c("WatrCont", "Shrub", "Substrate", "Topo", "SubsDens"))

#Re-run the RDA with the significant variables and look at the summary
mite.spe.rda.signif=rda(mite.spe~., data=mite.env.signif)
summary(mite.spe.rda.signif, display=NULL)

#Find the R2 adjusted of the model with the retained environmental variables
(R2adj <- RsquareAdj(mite.spe.rda.signif)$adj.r.squared)

#Determine the significant canonical (constrained) axes)
anova.cca(mite.spe.rda.signif, step=1000)
anova.cca(mite.spe.rda.signif, step=1000, by="axis")

#Plot the RDA
windows()
plot(mite.spe.rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.rda.signif, display="species", choices=c(1), scaling=1),
     scores(mite.spe.rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.rda.signif, display="species", scaling=1)),
     col="grey", cex=0.8)
arrows(0,0,
      scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1),
      scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1),
      col="red")
text(scores(mite.spe.rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(mite.spe.rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1)


##Section: 06-partial-redundancy-analysis.R 

#Divide the env2 dataframe into two dataframes:
envtopo <- env[, c(1:3)] # Physiography : explanatory dataset 1
names(envtopo)
envchem <- env[, c(4:10)] # Water quality : explanatory dataset 2
names(envchem)

#Run the partial RDA
spechem.physio=rda(spe.hel, envchem, envtopo)
summary(spechem.physio, display=NULL)
#or
spechem.physio2=rda(spe.hel ~ pH + dur + pho + nit + amm + oxy + dbo
                    + Condition(alt + pen + deb), data=env)

#Extract the results
summary(spechem.physio, display=NULL)

#Calculate the adjusted R2 of the partial RDA
(R2adj <- RsquareAdj(spechem.physio)$adj.r.squared)

#Test the significance of the axes in partial RDA
anova.cca(spechem.physio, step=1000)
anova.cca(spechem.physio2, step=1000, by="axis")

#Construct the triplots
#Scaling 1
windows(title="Partial RDA scaling 1")
plot(spechem.physio, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(spechem.physio, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(spechem.physio, display="species", choices=c(1), scaling=1),
       scores(spechem.physio, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(spechem.physio, display="species", choices=c(1), scaling=1),
     scores(spechem.physio, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spechem.physio, display="species", scaling=1)),
     col="black", cex=0.8)
arrows(0,0,
      scores(spechem.physio, display="bp", choices=c(1), scaling=1),
      scores(spechem.physio, display="bp", choices=c(2), scaling=1),
      col="red")
text(scores(spechem.physio, display="bp", choices=c(1), scaling=1)+0.05,
     scores(spechem.physio, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(spechem.physio, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1)

#Scaling 2
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

#Partial RDA
mite.spe.subs=rda(mite.spe.hel ~ Shrub + Topo
                  + Condition(SubsDens + WatrCont + Substrate), data=mite.env)

#Summary
summary(mite.spe.subs, display=NULL)
(R2adj <- RsquareAdj(mite.spe.subs)$adj.r.squared)

#Significant axes
anova.cca(mite.spe.subs, step=1000)
anova.cca(mite.spe.subs, step=1000, by="axis")

#Triplot scaling 1
windows(title="Partial RDA scaling 1")
plot(mite.spe.subs, scaling=1, main="Triplot partial RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
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

#Triplot scaling 2
windows(title="Partial RDA scaling 2")
plot(mite.spe.subs, scaling=2, main="Triplot partial RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(mite.spe.subs, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(mite.spe.subs, display="species", choices=c(1), scaling=2),
     scores(mite.spe.subs, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(mite.spe.subs, display="species", scaling=2)),
     col="black", cex=0.8)
arrows(0,0,
       scores(mite.spe.subs, display="bp", choices=c(1), scaling=2),
       scores(mite.spe.subs, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(mite.spe.subs, display="bp", choices=c(1), scaling=2)+0.05,
     scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(mite.spe.subs, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)


##Section: 07-variation-partitioning.R 

?varpart
vegandocs("partitioning.pdf")

#Variation partitioning with all explanatory variables
spe.part.all <- varpart(spe.hel, envchem, envtopo)
spe.part.all
windows(title="Variation partitioning - all variables")
plot(spe.part.all, digits=2)

#RDA of chemistry variables
spe.chem <- rda(spe.hel~., data=envchem)

#Select significant chemistry variables
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envchem),
           scope= formula(spe.chem), direction= "forward", R2scope=TRUE, pstep=1000)
names(envchem)
(envchem.pars <- envchem[, c( 4, 6, 7 )])

#RDA with other environmental variables
spe.topo <- rda(spe.hel~., data=envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
ordiR2step(rda(spe.hel~1, data=envtopo),
           scope= formula(spe.topo), direction= "forward", R2scope=TRUE, pstep=1000)
names(envtopo)
envtopo.pars <- envtopo[, c(1,2)]

#Varpart
spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars)
windows(title="Variation partitioning - parsimonious subsets")
plot(spe.part, digits=2)

#Tests of significance
anova.cca(rda(spe.hel, envchem.pars), step=1000) # Test of fractions [a+b]
anova.cca(rda(spe.hel, envtopo.pars), step=1000) # Test of fractions [b+c]
env.pars <- cbind(envchem.pars, envtopo.pars)
anova.cca(rda(spe.hel, env.pars), step=1000) # Test of fractions [a+b+c]
anova.cca(rda(spe.hel, envchem.pars, envtopo.pars), step=1000) # Test of fraction [a]
anova.cca(rda(spe.hel, envtopo.pars, envchem.pars), step=1000) # Test of fraction [c]

str(mite.env)
(mite.subs=mite.env[,c(1,2,3)]) #First set of variables outlined in challenge
(mite.other=mite.env[,c(4,5)]) #Second set of variables outlined in challenge

#RDA for mite.subs
rda.mite.subs <- rda(mite.spe.hel~., data=mite.subs)
R2a.all.subs <- RsquareAdj(rda.mite.subs)$adj.r.squared

#Forward selection for mite.subs
ordiR2step(rda(mite.spe.hel~1, data=mite.subs),
           scope= formula(rda.mite.subs), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.subs)
(mite.subs.pars <- mite.subs[, c(2, 3)])

#RDA for mite.other
rda.mite.other <- rda(mite.spe.hel~., data=mite.other)
R2a.all.other <- RsquareAdj(rda.mite.other)$adj.r.squared

#Forward selection for mite.other
ordiR2step(rda(mite.spe.hel~1, data=mite.other),
           scope= formula(rda.mite.other), direction= "forward", R2scope=TRUE, pstep=1000)
names(mite.other)
(mite.other.pars <- mite.other[, c(1,2)])

#Variation partitioning
(mite.spe.part <- varpart(mite.spe.hel, ~WatrCont+Substrate, ~Shrub+Topo,
                          data=mite.env))
windows(title="Variation partitioning - parsimonious subsets")
plot(mite.spe.part, digits=2)

# Tests of all testable fractions
anova.cca(rda(mite.spe.hel~ WatrCont+Substrate, data=mite.env), step=1000) # Test of fractions [a+b]
anova.cca(rda(mite.spe.hel~Shrub+Topo, data=mite.env), step=1000) # Test of fractions [b+c]
(env.pars <- cbind(mite.env[,c(2,3,4,5)]))
anova.cca(rda(mite.spe.hel~ WatrCont+Substrate+Shrub+Topo, data=env.pars), step=1000) # Test of fractions [a+b+c]
anova.cca(rda(mite.spe.hel~WatrCont+Substrate + Condition(Shrub+Topo), data=env.pars), step=1000) # Test of fraction [a]
anova.cca(rda(mite.spe.hel~Shrub+Topo+ Condition(WatrCont+Substrate ), data=env.pars), step=1000) # Test of fraction [c]


##Section: 08-multivariate-regression-tree.R 

?mvpart

#Prepare the data: remove “distance from source”
env <- subset(env, select = -das)

# Create the regression tree
doubs.mrt <- mvpart(as.matrix(spe.hel) ~. ,env,
            legend=FALSE, margin=0.01, cp=0, xv="pick",
            xval=nrow(spe.hel), xvmult=100, which=4)

# Using the CVRE criterion
doubs.mrt.cvre <- mvpart(as.matrix(spe.hel)~., env,
                 legend=FALSE, margin=0.01, cp=0,xv="pick",
                 xval=nrow(spe.hel), xvmult=100,which=4)

# Choosing ourself the best number of partitions
doubs.mrt.4 <- mvpart(as.matrix(spe.hel)~., env,
              legend=FALSE, margin=0.01, cp=0, xv="pick",
              xval=nrow(spe.hel), xvmult=100,which=4)

summary(doubs.mrt)

# Find discriminant species with MRT results
doubs.mrt.wrap<-MRT(doubs.mrt,percent=10,species=colnames(spe.hel))
summary(doubs.mrt.wrap)

# Extract indval p-values
doubs.mrt.indval<-indval(spe.hel,doubs.mrt$where)
doubs.mrt.indval$pval

# Extract indicator species of each node, with its indval
doubs.mrt.indval$maxcls[which(doubs.mrt.indval$pval<=0.05)]
doubs.mrt.indval$indcls[which(doubs.mrt.indval$pval<=0.05)]

mite.mrt<-mvpart(data.matrix(mite.spe.hel)~.,mite.env,
legend=FALSE,margin=0.01,cp=0,xv="pick",
xval=nrow(mite.spe.hel),xvmult=100,which=4)
summary(mite.mrt)

mite.mrt.wrap<-MRT(mite.mrt,percent=10,species=colnames(mite.spe.hel))
summary(mite.mrt.wrap)

mite.mrt.indval<-indval(mite.spe.hel,mite.mrt$where)
mite.mrt.indval$pval

mite.mrt.indval$maxcls[which(mite.mrt.indval$pval<=0.05)]
mite.mrt.indval$indcls[which(mite.mrt.indval$pval<=0.05)]


##Section: 09-linear-discriminant-analysis.R 

#load spatial data to determine groups
spa <- read.csv ('http://www.davidzeleny.net/anadat-r/data-download/DoubsSpa.csv', row.names = 1)
spa <- spa[,-8]

#View spatial data
View (spa)

#add site numbers
numbers<-(1:30)
numbers<-numbers[!numbers%in%8]
spa$site<-numbers

#make groups based on lattitude y<82=group1, 82<y<156=group2, y>156=group3
spa.group<-ddply(.data=spa, .variables=.(x, y, site), .fun= summarise, group = if(y <= 82) 1 else if (y <= 156) 2 else 3)

#order by site
spa.group<-spa.group[with(spa.group, order(site)), ]

#run LDA
LDA<-lda(env,spa.group[,4])

#classification of the objects based on LDA
spe.class <- predict(LDA)$class

#posterior probabilities of the objects to belong to the groups
spe.post <- predict(LDA)$posterior

#table of prior versus predicted classifications
spe.table <- table(spa.group[,4], spe.class)

#proportion of correct classification
diag(prop.table(spe.table, 1))

#predicting classification of new data
#read in new sites
classify.me<-read.csv("classifyme.csv", header = T)

#predict grouping of new data
predict.group<-predict(LDA, newdata=classify.me)

#give classification for each new site
group.new<-predict.group$class

mite.xy$site<-seq(1:70)
(max(mite.xy[,2])-min(mite.xy[,2]))/4

mite.xy.group<-ddply(.data=mite.xy, .variables=.(x, y, site), .fun= summarise, group = if(y <= 2.5) 1 else if (y <= 4.9) 2 else if (y <= 7.3) 3 else 4)
mite.xy.group<-mite.xy.group[with(mite.xy.group, order(site)), ]

LDA.mite<-lda(mite.env[,1:2],mite.xy.group[,4])
mite.class <- predict(LDA.mite)$class
mite.post <- predict(LDA.mite)$posterior
mite.table <- table(mite.xy.group[,4], mite.class)
diag(prop.table(mite.table, 1))


##Section: 10-final-considerations.R 

?cca #(constrained correspondence analysis)
# Constrained Correspondence Analysis (CCA) is a canonical ordination method similar to RDA that preserve
# Chi-square distances among object (instead of Euclidean distances in RDA). This method is well suited for the
# analysis of large ecological gradients.



?CCorA # Canonical Correlation Analysis

# Canonical Correlation Analysis (CCorA) differs from RDA given that the two matrices are considered symmetric
# while in RDA the Y matrix is dependent on the X matrix. The main use of this technique is to test the
# significance of the correlation between two multidimensional data sets, then explore the structure of the data by
# computing the correlations (which are the square roots of the CCorA eigenvalues) that can be found between
# linear functions of two groups of descriptors.


help(coinertia, package=ade4) # Coinertia Analysis

#Coinertia Analysis (CoIA) is a symmetric canonical ordination method that is appropriate to compare pairs
# of data sets that play equivalent roles in the analysis. The method finds a common space onto which the objects
# and variables of these data sets can be projected and compared. Compared to CCorA, co-inertia analysis
# imposes no constraint regarding the number of variables in the two sets, so that it can be used to compare
# ecological communities even when they are species-rich. Co-inertia analysis is not well-suited, however, to
# analyse pairs of data sets that contain the same variables, because the analysis does not establish one-to-one
# correspondences between variables in the two data sets; the method does not ‘know’ that the first variable is the
# same in the first and the second data sets, and likewise for the other variables.


help(mfa, package=ade4) # Multiple Factorial Analysis

# Multiple factor analysis (MFA) can be used to compare several data sets describing the same objects. MFA
# consists in projecting objects and variables of two or more data sets on a global PCA, computed from all data
# sets, in which the sets receive equal weights.


# Spatial analysis can be performed using packages AEM and PCNM : http://r-forge.r-project.org/R/?group_id=195


##Section: 11-references.R 




