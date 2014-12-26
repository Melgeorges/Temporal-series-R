#Chemin des donn�es
Chemin="c:/FormationsSuivies/Inge+Stat/ANIMATION_ORGANISATION/Orga_journees_2014/SupportsDavidInge+Stat/"
setwd(Chemin)

#Si absents, installer les packages suivants : astsa
#la commande est :
#install.packages("astsa")

####################################################################################################
#Exercice 1 : Analyse d'une série n'incluant pas de saisonalité (rendement du blé dans le Vaucluse)#
####################################################################################################

##Lecture du fichier
TAB<-read.table("Rendement_Ble_Vaucluse.txt", header=T, sep="\t")
plot(TAB$Time, TAB$Yield, xlab="Annee", ylab="Rendement (t ha-1)", type="l")
Time2<-TAB$Time^2
Time3<-TAB$Time^3

##Methode 1 : Regression

#Ajustement des modeles lineaire, quadratique et cubique
Mod.1<-lm(Yield~Time,data=TAB)
summary(Mod.1)
AIC(Mod.1)
Mod.2<-lm(Yield~Time+Time2,data=TAB)
summary(Mod.2)
AIC(Mod.2)
Mod.3<-lm(Yield~Time+Time2+Time3,data=TAB)
summary(Mod.3)
AIC(Mod.3)

#Analyse des residus
dev.new()
par(mfrow=c(2,2))
plot(TAB$Time,Mod.1$residuals)
abline(0,0)
acf(Mod.1$residuals) # autocorrelation
plot(TAB$Time,Mod.3$residuals)
abline(0,0)
acf(Mod.3$residuals)

#Prediction du rendements jusqu'en 2020
Time<-1950:2020
Pred<-predict(Mod.3, newdata=data.frame(Time=Time, Time2=Time^2, Time3=Time^3), se.fit=T)
dev.new()
plot(TAB$Time, TAB$Yield, xlab="Annee", ylab="Rendement (t ha-1)", type="l", xlim=c(1950,2020), ylim=c(1,6))
lines(Time,Pred$fit,lwd=2, col="red")
lines(Time,Pred$fit-1.96*Pred$se.fit,lty=2,col="red")
lines(Time,Pred$fit+1.96*Pred$se.fit,lty=2,col="red")

##Methode 2 : HoltWinters
Yield<-TAB$Yield
Yield.ts<-ts(Yield, start=1950, frequency=1)

dev.new()
par(mfrow=c(1,2))

#HW sans tendance
HW_Yield<-HoltWinters(Yield.ts, beta=F, gamma=F)
print(HW_Yield)
plot(HW_Yield, predict(HW_Yield, prediction.interval=T, n.ahead=7))

#HW avec tendance linéaire
HW_Yield<-HoltWinters(Yield.ts, gamma=F)
print(HW_Yield)
plot(HW_Yield, predict(HW_Yield, prediction.interval=T, n.ahead=7))

##Methode 3 : ARIMA

#ACF sur données initiales et données différenciées
dev.new()
par(mfrow=c(1,2))
acf(Yield.ts)    # reflette la tendance
acf(diff(Yield.ts))      # MA1 pic caract�ristique, parametre theta negatif

Nobs<-length(TAB$Yield)

Mod.1<-arima(Yield.ts,order=c(0,1,1),xreg=1:Nobs) # MA 1 donnees diff�renciees, xreg dedifferencie (ajuste une tendance lineaire)
print(Mod.1)
tsdiag(Mod.1)

Mod.2<-arima(Yield.ts,order=c(1,1,0),xreg=1:Nobs)   # AR 1 donnees differenciees
print(Mod.2)
tsdiag(Mod.2)

Mod.3<-arima(Yield.ts,order=c(1,1,1),xreg=1:Nobs)    # ARMA donnees differenciees
print(Mod.3)
tsdiag(Mod.3)

#Prediction  avec Mod.1
Pred<-predict(Mod.1, n.ahead=7,newxreg=(Nobs+1):(Nobs+7)) # indice des observations futures

dev.new()
plot(Yield.ts,xlim=c(1950, 2020),ylim=c(1,6))
lines(Pred$pred,col="red",lwd=2)
lines(Pred$pred-1.96*Pred$se,col="red",lwd=1,lty=2)
lines(Pred$pred+1.96*Pred$se,col="red",lwd=1,lty=2)

#package alternatif pour arima : library astsa
library(astsa)

Mod.1<-sarima(Yield.ts, 0,1,1,0,0,0)
print(Mod.1)

Mod.2<-sarima(Yield.ts, 1,1,0,0,0,0)
print(Mod.2)

Mod.3<-sarima(Yield.ts, 0,1,9,0,0,0)
print(Mod.3)

sarima.for(Yield.ts,n.ahead=7, 0,1,1,0,0,0)

################################################################################################
#Exercice 2 : Analyse d'une série incluant une saisonalité (bilan hydrique mensuel du Vaucluse)# 
################################################################################################

##Lecture du fichier 
EVT=read.csv("ETPpluvioVC.csv",sep=";",skip=9)

head(EVT)

#Evapotranspiration potentielle
ETP = EVT[,5]
#Precipitation
Precip = EVT[,4]

##Création d'objets "time series"
ETP.ts = ts(ETP,start=c(1970,10),frequency = 12)
Precip.ts = ts(Precip,start=c(1970,10),frequency = 12)
bilan = Precip - ETP
bilan.ts = ts(bilan,start=c(1970,10),frequency = 12)

##Graphiques des séries chronologiques
dev.new()
par(mfrow=c(3,1))
plot(ETP.ts)
plot(Precip.ts)
plot(bilan.ts)

par(mfrow=c(1,1))

##Décomposition des series (tendance, saison, résidus)
dev.new()
plot(decompose(ETP.ts))

dev.new()
plot(decompose(Precip.ts))

dev.new()
plot(decompose(bilan.ts))

##Extraction des composantes tendance et saison du bilan (Precipitation - ETP)
decompose(ETP.ts)$trend
decompose(bilan.ts)$seasonal

##Modèles sarima sur le bilan : on suppose d'abord qu'il n'y a pas de tendance

#ACF sur série initiale et sur la série désaisonnalisée
dev.new()
par(mfrow=c(1,2))
acf(bilan.ts)
acf(diff(bilan.ts,lag=12))

#MA(1) sur série désaisonalisée : Y(t)-Y(t-12) = Theta * e(t-12) + e(t)
Mod.1<-arima(bilan.ts, order=c(0,0,0), seasonal=list(order = c(0,1,1), period = 12))
Mod.1
dev.new()
tsdiag(Mod.1)

#AR(1) sur série désaisonalisée : Y(t)-Y(t-12) = Phi * (Y(t-1)-Y(t-13)) + e(t)
Mod.2<-arima(bilan.ts, order=c(0,0,0), seasonal=list(order = c(1,1,0), period = 12))
Mod.2
dev.new()
tsdiag(Mod.2)

#ARMA(1,1) sur série désaisonalisée : Y(t)-Y(t-12) = Phi * (Y(t-1)-Y(t-13)) + Theta * e(t-12) + e(t)
Mod.3<-arima(bilan.ts, order=c(0,0,0), seasonal=list(order = c(1,1,1), period = 12))
Mod.3
dev.new()
tsdiag(Mod.3)

#Prediction sur 3 ans avec Mod.1
Pred<-predict(Mod.1, n.ahead=400)

dev.new()
plot(bilan.ts,xlim=c(1970, 2050),ylim=c(-300,250))
lines(Pred$pred,col="red",lwd=2)
lines(Pred$pred-1.96*Pred$se,col="red",lwd=1,lty=2)
lines(Pred$pred+1.96*Pred$se,col="red",lwd=1,lty=2)

##Modèles sarima : on suppose maintenant qu'il y a une tendance linéaire

#ACF sur série initiale et sur la série différenciée
dev.new()
par(mfrow=c(2,2))
acf(bilan.ts)
acf(diff(bilan.ts))
acf(diff(bilan.ts,lag=12))
acf(diff(diff(bilan.ts,lag=12)))

#ARIMA(0,1,1)(0,1,1)12 : Y(t)-Y(t-12) - (Y(t-1)-Y(t-13)) = Theta*e(t-1)+ THETA*e(t-12)+ Theta*THETA*e(t-13) + e(t)
Mod.1<-arima(bilan.ts, order=c(0,1,1), seasonal=list(order = c(0,1,1), period = 12))
Mod.1
dev.new()
tsdiag(Mod.1)

#ARIMA(0,1,0)(0,1,1)12 : Y(t)-Y(t-12) - (Y(t-1)-Y(t-13)) = THETA*e(t-12)+ e(t)
Mod.2<-arima(bilan.ts, order=c(0,1,0), seasonal=list(order = c(0,1,1), period = 12))
Mod.2
dev.new()
tsdiag(Mod.2)

#ARIMA(0,1,1)(0,1,0)12 : Y(t)-Y(t-12) - (Y(t-1)-Y(t-13)) = Theta*e(t-1)+ e(t)
Mod.3<-arima(bilan.ts, order=c(0,1,1), seasonal=list(order = c(0,1,0), period = 12))
Mod.3
dev.new()
tsdiag(Mod.3)

#Prediction sur 3 ans avec Mod.1
Pred<-predict(Mod.1, n.ahead=36)

#dev.new()
#plot(bilan.ts,xlim=c(1970, 2016),ylim=c(-300,250))
#lines(Pred$pred,col="red",lwd=2)
#lines(Pred$pred-1.96*Pred$se,col="red",lwd=1,lty=2)
#lines(Pred$pred+1.96*Pred$se,col="red",lwd=1,lty=2)

#package alternatif pour arima : library astsa
library(astsa)

Mod.1<-sarima(bilan.ts, 0,1,1,0,1,1, 12)
print(Mod.1)

Mod.2<-sarima(bilan.ts, 0,1,0,0,1,1, 12)
print(Mod.2)

Mod.3<-sarima(bilan.ts, 0,1,1,0,1,0, 12)
print(Mod.3)

sarima.for(bilan.ts,n.ahead=36, 0,1,1,0,1,1, 12)

