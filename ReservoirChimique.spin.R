##################################
## Projet ANALYSE D'INCERTITUDE ##
##################################


## Questions 
# - Comment étudier la réponse global des deux sorties ? 
# - Est-il utile de maximiser le plan d'experience LHS avec les observations que l'on a ? 
# - Variabilité du Morris 

### Import package used
library(DiceKriging) #for km
library(sensitivity) #for morris
library(ggplot2)
library(randomForest)
library(rpart)
library(rpart.plot)

### Import data
chemin = '/Users/HUGO/Documents/Ecole/INSA/5GM/Incertitude/Projet'
#chemin = "~/Google Drive/INSA/Analyse d'incertitude/Projet"
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

#### Tree
model.tree = rpart(data$Y.X.wgt_calcite~.,data = data,cp = 0.00001)
printcp(model.tree)
param=model.tree$cptable
cp.optim=model.tree$cptable[which.min(model.tree$cptable[,"xerror"]),"CP"]
arbre=prune(model.tree,cp=cp.optim)
rpart.plot(arbre) # long execution (20 secondes)

#### RandomForest
model.RF = randomForest(data$Y.X.wgt_calcite~.,data = data)
varImpPlot(model.RF)
Imp = data.frame(name = rownames(model.RF$importance), IncNodePurity = model.RF$importance)
VarImp <-ggplot(data=Imp, aes(x=reorder(name,IncNodePurity), y=IncNodePurity)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="steelblue")
VarImp + coord_flip()

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
RF =function(X){
  X.scaled <- data.frame(t(apply(X,1,scale.unif)))
  return(predict(model.RF,X.scaled))
}

model.morris = morris(RF,factors =colnames(X),design = list(type = "oat", levels = 10, grid.jump = 3), r = 20)

#model.morris = morris(modeleRL,factors =colnames(X),design = list(type = "oat", levels = 5, grid.jump = 3), r = 4)

# Better plot
mu <- apply(model.morris$ee, 2, mean)
mu.star <- apply(model.morris$ee, 2, function(x) mean(abs(x)))
sigma <- apply(model.morris$ee, 2, sd)

#plot(mu.star, sigma)
marice= data.frame(name = model.morris$factors,mu.star,sigma)
d=data.frame(x1=c(0,0.003,0.002), x2=c(0.002,0.005,max(mu.star)+0.001), y1=c(0,0.001,0.008), y2=c(0.004,0.006,max(sigma)+0.001),
             t =c('Non influents','Effets linéaires','Effets non linéaires ou interactions') , 
             r=c(1,2,3))

g <- ggplot(marice) +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
  geom_point(aes(x=mu.star, y=sigma,label = model.morris$factors),size=1) +
  geom_text(aes(x=mu.star, y=sigma,label = model.morris$factors),size=3,hjust=1.1)
g
# Remark : 
# With the last graph, we can see the influent variable with a non linear effect : 
# - Kine, Kppt and K2diss Calcite
# - Kine and K2diss Dolomite. 
# Sometimes, there is 2 variables with an linear effect : 
# - K2diss and Kine illite. 

#### ANNOVA
# Pour connaitre les effets d'interactions : 
annova = aov(data$Y.X.wgt_calcite~(.)^2,data=data)
interaction = annova$coefficients[-1] # On enlève l'intercept
interaction = interaction[which(interaction>10^-3)]
#plot(annova)
interaction = data.frame(Mixing = names(interaction),value = interaction)
interaction = interaction[order(interaction[,2], decreasing = T),]
p<-ggplot(data=interaction, aes(x=reorder(Mixing,value), y=value)) +
  geom_bar(stat="identity") 
p <- p + geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(interaction$value,5)), hjust=1.6,color="white", size=3.5)
p + coord_flip()

#### Correlation between variables
data_full = data.frame(y= Y,X)
install.packages("corrplot")
library(corrplot)
cor <- cor(data_full)
quartz() # Pour l'afficher dans une autre fenetre
corrplot(cor, type="upper", order = "AOE", tl.col="black", tl.srt=55, tl.cex = 0.7)



#### R shiny
server=function(input, output) {
  
  output$trendPlot <- renderPlot({
    
    model.morris = morris(RF,factors =colnames(X),design = list(type = "oat", levels = 10, grid.jump = 3), r = 20)
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
  })
}
ui=fluidPage(
  titlePanel("MORRIS"),
  mainPanel(
    plotOutput("trendPlot")
  )
)

shinyApp(ui,server)
    
    