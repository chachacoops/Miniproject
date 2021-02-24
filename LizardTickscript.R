rm(list=ls())	
install.packages("lmtest")
install.packages("usdm")
install.packages("psych")
install.packages("lmerTest")
install.packages("dplyr")
install.packages("lme4")
install.packages("car")
install.packages("ggplot2")
install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(lmtest)
library(usdm)
library(psych)
library(lmerTest)
library(dplyr)
library(lme4)
library(car)
library(ggplot2)
require(glmmADMB)
setwd("~/Desktop/Masters/Masters Scripts/Mini Project")

##### Data preparation

Lizard<-read.csv("Lizards.csv", header=T) #load data
View(Lizard)


Liz1<- dplyr::select(Lizard,Intensity,HABITAT,Season,
                     Month,Hab_RD,Elevation,Floodness,ID,Year,Canopy) #select variables (maybe not necessary I just like a clean set!)
Liz1 %>% group_by(ID) %>% summarise(count=length(ID)) %>% count(count) #test for repeated individuals
#although lizards can be repeated, they usually collected more ticks by next recapture
Liz1[["Intensity"]][is.na(Liz1[["Intensity"]])] <- 0 #convert NAs to 0s

# check for normality in response variable

hist1 <- ggplot(Liz1, aes(Intensity))+
  geom_histogram(fill="grey", binwidth=5)+  
  theme(axis.text = element_text(size = 11))+
  theme(axis.title = element_text(size = 15))+
  labs(x="Tick intensity")
hist1

Liz1<- filter(Liz1, Intensity<=200) #remove outlier of 200+ ticks and rerun histogram

#recheck histogram
hist2 <- ggplot(Liz1, aes(Intensity))+
  geom_histogram(fill="grey", binwidth=5)+  
  theme(axis.text = element_text(size = 11))+
  theme(axis.title = element_text(size = 15))+
  labs(x="Tick intensity")
hist2


###### boxplot

boxplot(Liz1$Intensity ~ Liz1$HABITAT) #doesn't look great
Liz1$sqrt.I <- sqrt(Liz1$Intensity) #square root intensity for better boxplot
#make fancy boxplot
bplot1 <- ggplot(Liz1, aes(x=HABITAT, y=sqrt.I, fill =HABITAT))+
  geom_boxplot()+
  labs(x="Habitat type", y="Tick intensity (sqrt)", fill="Habitat type")+
  scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))+
  theme(axis.text = element_text(size = 11))+
  theme(axis.title = element_text(size = 13))
bplot1
#looks better, easier to intepret!

##### comparison between oil palm and forest
wilcox.test(Liz1$Intensity ~ Liz1$HABITAT,conf.int=TRUE) #test habitat influence hypothesis
#W = 41271, p-value = 7.919e-12
test<-wilcox.test(Liz1$Intensity ~ Liz1$HABITAT) #calculate effect size
zstat<-qnorm(test$p.value/2)
abs(zstat)/sqrt(20)
#1.529472
    
Liz1$Year <- as.factor(Liz1$Year) #change to factors to check as random effect

summary(Liz1)
sum(Liz1$Intensity)
range(Liz1$Intensity)
#SE function
sE <- function(x){SD2 <- sd(x)
N2 <- sum(complete.cases(x))
SE2 <- SD2/sqrt(N2); SE2
}
sE(Liz1$Intensity)

###### figuring out which model to use
mean(Liz1$Intensity)
var(Liz1$Intensity)
#variance is not equal to the mean. Therefore poisson model will not be suitable


##### oil palm plantation (OPP) analysis

Liz2<- Liz1 #create new dataframe for OPP
Liz2<- filter(Liz2, HABITAT=="OPP") #filter all which habitat is equal to OPP
View(Liz2)
summary(Liz2)

#descriptive statistics for oil palm
var(Liz2$Intensity)
mean(Liz2$Intensity)
hist(Liz2$Intensity)
mean(Liz2$Intensity)
range(Liz2$Intensity)

#calculate means and SE for transect location
aggregate(Liz2[, 1,5], list(Liz2$Hab_RD), mean) 
aggregate(Liz2[, 1,5], list(Liz2$Hab_RD), sE) 
#calculate means and SE for canopy heights
aggregate(Liz2[, 1,10], list(Liz2$Canopy), mean) 
aggregate(Liz2[, 1,10], list(Liz2$Canopy), sE) 
#calculate means and SE for season
aggregate(Liz2[, 1,3], list(Liz2$Season), mean) 
aggregate(Liz2[, 1,3], list(Liz2$Season), sE) 


#test poisson model to see which other models to try
oppP1 <- glm(Intensity~Season+Hab_RD+Floodness+Elevation+Canopy, data=Liz2, family="poisson") #remove floodness
overdispersion <- oppP1$deviance/oppP1$df.residual
overdispersion #high at 9.18


#quassipoisson model test
QPopp1 <- glm(Intensity~Season+Hab_RD+Floodness+Elevation+Canopy, data=Liz2, family=quasipoisson) #remove floodness
summary(QPopp1)

#checking residual distributions
devresid <- resid(QPopp1, type = "deviance")
plot(devresid ~ QPopp1$fitted.values)
#overdispersion statistic
od <- QPopp1$deviance/QPopp1$df.residual #over dispersion stat
od
# overdispersion stat is same as poisson model 9.18


#negative binomial models
opp1 <- glm.nb(Intensity~Season+Hab_RD+Floodness+Elevation+Canopy, data=Liz2) #remove floodness
summary(opp1)

opp2 <- glm.nb(Intensity~Season+Hab_RD+Elevation+Canopy, data=Liz2) #remove floodness
summary(opp2)
lrtest(opp1,opp2)

opp3 <- glm.nb(Intensity~Season+Hab_RD+Canopy, data=Liz2) #remove elevation
summary(opp3)
lrtest(opp2,opp3)

#checking residual distributions for opp3 model
devresid <- resid(opp3, type = "deviance")
plot(devresid ~ opp3$fitted.values)
#overdispersion statistic
od <- opp3$deviance/opp3$df.residual #over dispersion stat
od #1.07 looks good!

#pseudo R^2 for opp3
pR2 <- (opp3$null.deviance-opp3$deviance)/opp3$null.deviance
pR2 #0.12 - not a huge amount but still worthwhile

##### forest analysis

Liz3<- Liz1
Liz3<- filter(Liz3, HABITAT=="FOREST")
View(Liz3)
summary(Liz3)
mean(Liz3$Intensity)
range(Liz3$Intensity)

#quassipoisson model test

QPfor1<- glm(Intensity~Season+Hab_RD+Elevation+Canopy, data=Liz3, family=quasipoisson) #remove floodness
summary(QPfor1)

#checking residual distributions
devresid <- resid(QPfor1, type = "deviance")
plot(devresid ~ QPfor1$fitted.values)
#overdispersion statistic
od <- QPfor1$deviance/QPfor1$df.residual #over dispersion stat
od
# overdispersion stat for quassipoisson model - 5.56267 bit too high


#negative binomial
for1 <- glm.nb(Intensity~Season+Hab_RD+Elevation+Canopy, data=Liz3) #remove floodness
summary(for1)
for2 <- glm.nb(Intensity~Season+Hab_RD+Elevation, data=Liz3) #remove canopy
summary(for2)
lrtest(for1,for2)

for3 <- glm.nb(Intensity~Season+Hab_RD, data=Liz3) #remove elevation
summary(for3)

for4 <- glm.nb(Intensity~Hab_RD, data=Liz3) #remove elevation
summary(for4)

#non significance for all forest models! 


