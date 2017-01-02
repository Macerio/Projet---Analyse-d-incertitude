##################################
## Projet ANALYSE D'INCERTITUDE ##
##################################

### Import package used
install.packages('DiceKriging')
install.packages('sensitivity')
library(DiceKriging) #for km
library(sensitivity) #for morris
library(ggplot2)

### Import data
#chemin = '/Users/HUGO/Documents/Ecole/INSA/5GM/Incertitude/Projet'
chemin = "~/Google Drive/INSA/Analyse d'incertitude/Projet"
jdd.1 = read.csv(paste(chemin,'/jdd_1.csv',sep=""),header= T )

### Look up and clean the data 
head(jdd.1)
# Pour chaque espece : 
# - taux de dissolution
# - taux de précipitation
# - taux de surface de réaction
summary(jdd.1)
X = jdd.1[,1:(length(jdd.1)-2)]
for (i in 1:dim(X)[2]){
  X[,i]=as.numeric(X[,i])
}
boxplot(X) # pensez à scale

Y = jdd.1[,(length(jdd.1)-1):(length(jdd.1))]
Y1 = data.frame(Y$X.wgt_calcite)
Y2 = data.frame(Y$X.wgt_chlorite)
boxplot(Y1)
boxplot(Y2)
X = scale(X)
X = data.frame(X)
boxplot(X)
data = data.frame(y= Y1,X)
### Create models
#### régression linéaire
modeleRL <- lm(data$Y.X.wgt_calcite~.,data = data)
summary(modeleRL)

### Selection des bonnes variables
model.empty <- lm(data$Y.X.wgt_calcite~1,data = data)
model.both = step(model.empty,scope=list(lower= model.empty,upper= modeleRL) ,direction = "both")
summary(model.both)

#### métamodélisation par processus gaussien
# Attention très très long...
# Normal, il y a bcp de variable. Ce n'est pas un modèle appropié
#modelePG <- km(Y1~X,design=X,response=Y1)

#### Morris
binf = apply(X,2,min)
bsup = apply(X,2,max)
scale.unif <- function(x) {
  return(x*(bsup-binf) + binf)
}
RL =function(X){
  X.scaled <- data.frame(t(apply(X,1,scale.unif)))
  return(predict(modeleRL,X.scaled))
}


model.morris = morris(RL,factors =colnames(X),design = list(type = "oat", levels = 10, grid.jump = 3), r = 4)
#model.morris = morris(modeleRL,factors =colnames(X),design = list(type = "oat", levels = 5, grid.jump = 3), r = 4)
#quartz()
plot(model.morris)

# Pour connaitre les effets d'interactions : 
annova = aov(data$Y.X.wgt_calcite~(.)^2,data=data)
interaction = annova$coefficients[-1] # On enlève l'intercept
interaction = interaction[which(interaction>10^-3)]
barplot(sort(interaction, decreasing = T), las = 2)
plot(annova)


install.package("ggplot2")
library("ggplot2")
qplot(interaction, geom="bar", stat="identity")
names(interaction)
