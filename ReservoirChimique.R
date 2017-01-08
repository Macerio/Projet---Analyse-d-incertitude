##################################
## Projet ANALYSE D'INCERTITUDE ##
##################################



### Import package used
library(DiceKriging) #for km
library(sensitivity) #for morris
library(ggplot2) # plot 
library(randomForest) # RandomForest model
library(rpart) # Tree model
library(rpart.plot) # plot tree
library(corrplot) # plot correlation
library(randtoolbox) # for sobol 
library(DiceDesign) # For mindist
library(plotly) # plot


### Import data
chemin = '/Users/HUGO/Documents/Ecole/INSA/5GM/Incertitude/Projet'
#chemin = "~/Google Drive/INSA/Analyse d'incertitude/Projet"
jdd.1 = read.csv(paste(chemin,'/jdd_1.csv',sep=""),header= T )
jdd.2 = read.csv(paste(chemin,'/jdd_2.csv',sep=""),header= T )


############
## MORRIS ##
############


### Look up and clean the data 
head(jdd.1)
# Pour chaque espece : 
# - taux de dissolution
# - taux de précipitation
# - taux de surface de réaction
summary(jdd.1)
X.1 = jdd.1[,1:(length(jdd.1)-2)]
for (i in 1:dim(X.1)[2]){
  X.1[,i]=as.numeric(X.1[,i])
}
boxplot(X.1) # pensez à scale


Y.1 = jdd.1[,(length(jdd.1)-1):(length(jdd.1))]
Y.1.1 = data.frame(Y.1$X.wgt_calcite)
Y.1.2 = data.frame(Y.1$X.wgt_chlorite)
boxplot(Y.1.1)
boxplot(Y.1.2)
X.1 = scale(X.1)
X.1 = data.frame(X.1)
boxplot(X.1)

#### Correlation between variables
data_full.1 = data.frame(y= Y.1,X.1)
cor <- cor(data_full.1)
quartz() # Pour l'afficher dans une autre fenetre
corrplot(cor, type="upper", order = "AOE", tl.col="black", tl.srt=55, tl.cex = 0.7)

#### Distribution of each variables
par(mfrow=c(2,2))
for (i in 3:dim(data_full.1)[2]){
  hist(data_full.1[,i],breaks = 30)
}
par(mfrow=c(1,1))

# There is no particulary form on each variable. 
# For morris, we can use a uniform distribution for each variable. 

####################
## Calcite weight ##
####################

data = data.frame(y= Y.1.1,X.1)

### Create models
#### régression linéaire
model.RL <- lm(data$Y.1.X.wgt_calcite~.,data = data)
summary(model.RL)

#### Tree
model.tree = rpart(data$Y.1.X.wgt_calcite~.,data = data,cp = 0.00001)
printcp(model.tree)
param=model.tree$cptable
cp.optim=model.tree$cptable[which.min(model.tree$cptable[,"xerror"]),"CP"]
arbre=prune(model.tree,cp=cp.optim)
rpart.plot(arbre) # long execution (20 secondes)

#### RandomForest
model.RF = randomForest(data$Y.1.X.wgt_calcite~.,data = data)
Imp = data.frame(name = rownames(model.RF$importance), IncNodePurity = model.RF$importance)
VarImp <-ggplot(data=Imp, aes(x=reorder(name,IncNodePurity), y=IncNodePurity)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="steelblue")
VarImp + coord_flip()

### Selection des bonnes variables
model.empty <- lm(data$Y.1.X.wgt_calcite~1,data = data)
model.both = step(model.empty,scope=list(lower= model.empty,upper= model.RL) ,direction = "both")
summary(model.both)
# We keep 21 variables. 
names(model.both$coefficients)[which(names(model.both$coefficients) %in% names(X.2))]
# Most of theme are in the file jdd.2 

#### Morris

# Parameters : 
# Dans ‘type’ on a Level qui représente le nombre de niveau de la grille dans lequel peut se déplacer la variable
# grid.jump c’est le nombre de niveau dans lequel la variable va de déplacer en diminuant ou augmentant (l/2)
# r c’est le nombre de variations appliquées à chaque variable pour voir les effets

binf = apply(X.1,2,min)
bsup = apply(X.1,2,max)
scale.unif <- function(x) {
  return(x*(bsup-binf) + binf)
}
RL =function(X){
  X.scaled <- data.frame(t(apply(X,1,scale.unif)))
  return(predict(model.RL,X.scaled))
}
RF =function(X){
  X.scaled <- data.frame(t(apply(X,1,scale.unif)))
  return(predict(model.RF,X.scaled))
}

model.morris = morris(RF,factors =colnames(X.1),design = list(type = "oat", levels = 10, grid.jump = 3), r = 20)
# Better plot
mu <- apply(model.morris$ee, 2, mean)
mu.star <- apply(model.morris$ee, 2, function(x) mean(abs(x)))
sigma <- apply(model.morris$ee, 2, sd)

marice= data.frame(name = model.morris$factors,mu.star,sigma)
d=data.frame(x1=c(0,0.003,0.002), x2=c(0.002,0.005,max(mu.star)+0.001), y1=c(0,0.001,0.008), y2=c(0.004,0.006,max(sigma)+0.001),
             t =c('Non influents','Effets linéaires','Effets non linéaires ou interactions') , 
             r=c(1,2,3))

g <- ggplot(marice) +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
  geom_point(aes(x=mu.star, y=sigma,label = model.morris$factors),size=1) +
  geom_text(aes(x=mu.star, y=sigma,label = model.morris$factors),size=3,hjust=1.1)
g
ggplotly(g) # Use for zoom
# Remark : 
# With the last graph, we can see the influent variable with a non linear effect : 
# - Kine, Kppt and K2diss Calcite
# - Kine and K2diss Dolomite. 
# Sometimes, there is 2 variables with an linear effect : 
# - K2diss and Kine illite. 

#### ANNOVA
# To identify the different interaction between variables and the important variable for the output Y1 : 
annova_Morris = aov(data$Y.1.X.wgt_calcite~(.)^2,data=data)
interaction = annova_Morris$coefficients[-1] # On enlève l'intercept
#plot(annova)
interaction = data.frame(Mixing = names(interaction),value = interaction)
interaction = interaction[order(interaction[,2], decreasing = T),]
interaction_Morris = interaction [1:25,]
anova_Morris<-ggplot(data=interaction_Morris, aes(x=reorder(Mixing,value), y=value)) +
  geom_bar(stat="identity") 
anova_Morris <- anova_Morris + geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(interaction_Morris$value,5)), hjust=1.6,color="white", size=3.5)
anova_Morris + coord_flip() + ggtitle("Anova with Morris")

# Remark : The important variables are the following : 
#                                   value
# Kine_Dolomite                    0.002696630
# Kppt_Calcite                     0.002254608
# Kine_Illite                      0.002154100
# K2diss_Dolomite                  0.002005415
# K2diss_Illite:Kine_Dolomite      0.001909999
# K2diss_Illite                    0.001612849
# Kine_Illite:Kine_Kaolinite       0.001427935
# K2diss_Dolomite:Kine_Gibbsite    0.001376154
# K2diss_Smectite:Kppt_Quartz      0.001266825
# Kppt_Smectite:Kine_Illite        0.001150967
# K2diss_Kaolinite:Kine_Microcline 0.001132863


#### discrepance ####
mindist(X.1)



########
# LHS ##
########

### Look up and clean the data 
head(jdd.2)
summary(jdd.2)
X.2 = jdd.2[,1:(length(jdd.2)-2)]
Y.2 = jdd.2[,(length(jdd.2)-1):(length(jdd.2))]
Y.2.1 = data.frame(Y.2$X.wgt_calcite)
Y.2.2 = data.frame(Y.2$X.wgt_chlorite)
X.2 = scale(X.2)
X.2 = data.frame(X.2)

#### Correlation between variables
data_full.2 = data.frame(y= Y.2,X.2)
cor <- cor(data_full.2)
quartz() # Pour l'afficher dans une autre fenetre
corrplot(cor, type="upper", order = "AOE", tl.col="black", tl.srt=55, tl.cex = 0.7)

data = data.frame(y= Y.2.1,X.2)

annova_LHS = aov(data$Y.2.X.wgt_calcite~(.)^2,data=data)
interaction = annova_LHS$coefficients[-1] # On enlève l'intercept
interaction = data.frame(Mixing = names(interaction),value = interaction)
interaction = interaction[order(interaction[,2], decreasing = T),]
interaction_LHS = interaction [1:25,]
anova_LHS<-ggplot(data=interaction_LHS, aes(x=reorder(Mixing,value), y=value)) +
  geom_bar(stat="identity") 
anova_LHS <- anova_LHS + geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(interaction_LHS$value,5)), hjust=1.6,color="white", size=3.5)
anova_LHS + coord_flip() + ggtitle("Anova with LHS data")

# Top interaction result of anova : 
#K2diss_Illite 0.0018723245
#Kine_Illite 0.0017312695
#K2diss_Dolomite 0.0017251225
#Kine_Dolomite 0.0017105872
#Kppt_Calcite 0.0009016229
#K2diss_Illite:Kine_Illite 0.0008265193
#Kine_Chlorite 0.0004496297

# To compare with the LHS anova, we have some similitude as : 
# We find some variable again in the top importante variables : 
interaction_LHS[which(rownames(interaction_LHS) %in% rownames(interaction_Morris)),]
# - K2diss_Illite
# - Kine_Illite
# - K2diss_Dolomite
# - Kine_Dolomite
# - Kppt_Calcite


#### kriging model as processus gaussien
# Time to comput : 10min 
model.PG <- km(Y.2.1~.,design=X.2,response=Y.2.1)
#getwd(chemin)
#save(model.PG, file = "model_PG.RData")
#load(file = "model_PG.RData")
# validation croisée
LOO<-leaveOneOut.km(model.PG, type="UK", trend.reestim=FALSE)
err.VC1=mean((LOO$mean-zi_1)^2/LOO$sd^2)

#### RandomForest
model.RF_LHS =randomForest(data$Y.2.X.wgt_calcite~.,data = data)
Imp_LHS = data.frame(name = rownames(model.RF_LHS$importance), IncNodePurity = model.RF_LHS$importance)
VarImp_LHS <-ggplot(data=Imp_LHS, aes(x=reorder(name,IncNodePurity), y=IncNodePurity)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="steelblue")
VarImp_LHS + coord_flip()

#### Sobol indice with RF model 
n <- 1000
set.seed(1234)
X1 <- matrix(runif(19 * n,min = -1.7,max= 1.7), nrow=n)
colnames(X1) = names(X.2)
X2 <- matrix(runif(19 * n,min = -1.7,max= 1.7), nrow=n)
colnames(X2) = names(X.2)
wwww.sobol <- sobol2002(model.RF_LHS,X1=X1, X2=X2)
# First order indice sobol
indice_sobol= data.frame(name = rownames(wwww.sobol$S), original =  c(wwww.sobol$S))
plot_sobol = ggplot(data=indice_sobol,aes(x=reorder(name,-original), y=original)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="steelblue")
plot_sobol + coord_flip() + ggtitle("sobol indice with LSH")

# indice_sobol[order(indice_sobol[,2], decreasing = F),]
#       Kine_Illite 0.3093727
#     Kine_Chlorite 0.3199128
#      Kppt_Gibbsite 0.4224477
#     Kppt_Kaolinite 0.4309436
#     Kine_Gibbsite 0.4376032
#     K2diss_Calcite 0.4389109
#      Kine_Calcite 0.4401204
#    Kine_Kaolinite 0.4431233
#   K2diss_Anhydrite 0.4433615
#  K2diss_Microcline 0.4434484
#    Kine_Anhydrite 0.4451381
#     Kppt_Smectite 0.4455828
#   Kine_Microcline 0.4466226
#     Kine_Smectite 0.4518645
#    K2diss_Chlorite 0.4579467
#      K2diss_Illite 0.4853016
#       Kppt_Calcite 0.5499403







                #####################
                ## Chlorite weight ##
                #####################

data = data.frame(y= Y.1.2,X.1)
# We have to erase the both variables K2diss and Kine chlorite because there 
#are too correlate and this create a bad model 
data = data[,-which(names(data) %in% c("K2diss_Chlorite","Kine_Chlorite"))]


#### Tree
model.tree = rpart(data$Y.1.X.wgt_chlorite~.,data = data,cp = 0.00001)
printcp(model.tree)
param=model.tree$cptable
cp.optim=model.tree$cptable[which.min(model.tree$cptable[,"xerror"]),"CP"]
arbre=prune(model.tree,cp=cp.optim)
#rpart.plot(arbre) # long execution (20 secondes)

#### RandomForest
model.RF = randomForest(data$Y.1.X.wgt_chlorite~.,data = data)
Imp = data.frame(name = rownames(model.RF$importance), IncNodePurity = model.RF$importance)
VarImp <-ggplot(data=Imp, aes(x=reorder(name,IncNodePurity), y=IncNodePurity)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="steelblue")
VarImp + coord_flip()

#### Morris

model.morris = morris(RF,factors =colnames(X.1),design = list(type = "oat", levels = 10, grid.jump = 3), r = 20)
mu <- apply(model.morris$ee, 2, mean)
mu.star <- apply(model.morris$ee, 2, function(x) mean(abs(x)))
sigma <- apply(model.morris$ee, 2, sd)

marice= data.frame(name = model.morris$factors,mu.star,sigma)
d=data.frame(x1=c(0,0.00075,0.00055), x2=c(0.0005,max(mu.star)+0.0001,max(mu.star)+0.0001), y1=c(0,0,0.001), y2=c(0.0005,0.0005,max(sigma)+0.0001),
             t =c('Non influents','Effets linéaires','Effets non linéaires ou interactions') , 
             r=c(1,2,3))

g <- ggplot(marice) +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
  geom_point(aes(x=mu.star, y=sigma,label = model.morris$factors),size=1) +
  geom_text(aes(x=mu.star, y=sigma,label = model.morris$factors),size=3,hjust=1.1)
g
ggplotly(g) # to zoom
# Variable with a non linear effect : 
# - Kppt calcite
# - K2diss dolomite 
# - Kppt Kaolinite
# - K2diss Microcline 
# - Kine Quartz

# We don't have other effects. 

#### ANOVA
# To identify the different interaction between variables and the important variable for the output Y1 : 
annova = aov(data$Y.1.X.wgt_chlorite~(.)^2,data=data)
interaction = annova$coefficients[-1] # On enlève l'intercept
interaction = data.frame(Mixing = names(interaction),value = interaction)
interaction = interaction[order(interaction[,2], decreasing = T),]
interaction = interaction [1:25,]
p<-ggplot(data=interaction, aes(x=reorder(Mixing,value), y=value)) +
  geom_bar(stat="identity") 
p <- p + geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(interaction$value,5)), hjust=1.6,color="white", size=3.5)
p + coord_flip() + ggtitle("Anova of interaction and influent variables")

# We can see that the variable Kine Quartz is present in a lot of interaction 












#### Shity code
result =matrix(nrow = 4,ncol=30)
k=1
good = c("Kppt_Calcite","K2diss_Dolomite","Kppt_Kaolinite","K2diss_Anhydrite")
ind = which(names(X)==good[4])
summary(X[,ind])
change = seq(-1.5,1.5,length.out = 30)
for (i in change){
  vector = data.frame(X[1,1:(ind-1)],K2diss_Anhydrite = change[k],X[1,(ind+1):dim(X)[2]])
  result[4,k]=predict(model.RF,vector)
  k=k+1
}
par(mfrow=c(2,2))
for (i in 1:4) {
  plot(result[i,])
}
par(mfrow=c(1,1))
