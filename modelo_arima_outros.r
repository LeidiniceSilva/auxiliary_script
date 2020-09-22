##################################################################
## script modelo para ARIMA e outros                           ###
##                                                             ###
## desenvolvido por: Amanda Sampaio                            ###
## contato: amandaafs34@gmail.com                              ###
##################################################################

####################################################
## Carregando os pacotes                         ###
##                                               ###
## certifique-se de que todos os pacotes         ###
## usados abaixo estão instalados                ###
####################################################

library(tseries)
library(timeSeries)
library(forecast)
library(fBasics)
library(lmtest)#Para accuracy
library(timetools)
library(fBasics)
library(hydroGOF)
library(xts)

################################# Modelo de previsão ARIMA #################################
dir()
setwd('C:\\data//P3')#Definindo diretório de trabalho
getwd()
dir()

#####################################
#### Lendo e manipulando os dados####
#####################################

#Fortaleza cidade utilizada no estudo
fortal192meses <- read.table("F_mensal_2003_2018.csv",header=T,sep=",",dec=",") #dados sazonal
attach(fortal192meses)
names(fortal192meses)
fix(fortal192meses)#Facilita a visualização e pode editar também
length(Fortaleza)#Aqui essa função indica o número de elementos na coluna que quero no caso Fortaleza.

#####################################
###       Modelo ARIMA            ###
#####################################
#Série Temporal#
(velfort.ts=ts(fortal192meses$Fortaleza, start=c(2003,1), end=c(2018, 12),freq=12))
#ARIMA
fortal.arima=auto.arima(velfort.ts)
fortal.arima
Series: velfort.ts 
ARIMA(2,0,0)(2,1,2)[12] 
#TIPO:SARIMA(2,0,0)(2,1,2)
names(fortal.arima)
attach(fortal.arima)
fortaleza.fit=arima(velfort.ts,order=c(2,0,0),seasonal=list(order=c(2,1,2),period=12))
names(fortaleza.fit)
attach(fortaleza.fit)
###The following objects are masked from fortal.arima:

#aic, arma, call, code, coef, loglik, mask, model, n.cond, nobs, residuals, series, sigma2, var.coef

fitted(fortaleza.fit)
#######################################################################################################################
fortal.arima.fit=c(fitted(fortaleza.fit))
fortal.arima.fit.ts = ts(fortal.arima.fit, start=c(2003,1), end=c(2018, 12),freq=12)
length(fortal.arima.fit.ts)#192 meses

length(velfort.ts)#192 meses
length(fortal.arima.fit.ts)#192 meses
#Erros ARIMA
#Total hora
#01/2003 - 12/2018 de jan até dez
accuracy(velfort.ts,fortal.arima.fit.ts)
#RESULTADO DA ACURACIA
#               ME     RMSE       MAE      MPE     MAPE        ACF1 Theil's U
Test set 0.03591845 0.296938 0.2220829 1.291337 8.658903 -0.02015318 0.7222914
#Coeficiente eficiencia#precisa dessa-> library(hydroGOF)
NSE(fortal.arima.fit.ts,velfort.ts)
0.8377462
#####################################################################################################################