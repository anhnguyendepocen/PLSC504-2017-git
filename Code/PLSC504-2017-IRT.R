#########################################################
# PLSC 504 -- Fall 2017
#
# Item Response Models
#
########################################################
# Load packages (install as needed), set options:

library(RCurl)
library(lme4)
library(plm)
library(gtools)
library(plyr)
library(texreg)
library(statmod)
library(ltm)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

#######################
# SCOTUS voting data

url <- getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/SCOTUS-IRT.csv")
SCOTUS <- read.csv(text = url) 
rm(url)

head(SCOTUS,10)
summary(SCOTUS)

# 1PLM:

OnePLM<-rasch(SCOTUS[c(2:10)])
summary(OnePLM)

coef(OnePLM, prob=TRUE, order=TRUE)

# Alternative model constraining alpha = 1.0:

IRTData <- SCOTUS[c(2:10)]

AltOnePLM<-rasch(IRTData, constraint=cbind(length(IRTData)+1,1))
summary(AltOnePLM)

# 2PLM:

TwoPLM<-ltm(IRTData ~ z1)
summary(TwoPLM)

# 2PLM Probabilities and testing:

coef(TwoPLM, prob=TRUE, order=TRUE)
anova(OnePLM, TwoPLM)

# 3PLM:

ThreePLM<-tpm(IRTData)
summary(ThreePLM)

anova(TwoPLM, ThreePLM)

# Plots:

pdf("1PLMIRFsR.pdf",6,5)
par(mar=c(4,4,2,2))
plot(OnePLM,lty=c(1,2,3,4,5,6,7,8,9), lwd=3, 
     zrange=c(-2.5,2.5),xlab="Liberalism",
     legend=TRUE,main="1PLM ICCs")
dev.off()

pdf("2PLMIRFsR.pdf",6,5)
par(mar=c(4,4,2,2))
plot(TwoPLM,lty=c(1,2,3,4,5,6,7,8,9), lwd=3, 
     zrange=c(-2.5,2.5),xlab="Liberalism",
     legend=TRUE,main="2PLM ICCs")
dev.off()

pdf("3PLMIRFsR.pdf",6,5)
par(mar=c(4,4,2,2))
plot(ThreePLM,lty=c(1,2,3,4,5,6,7,8,9), lwd=3, 
     zrange=c(-2.5,2.5),xlab="Liberalism",
     legend=TRUE,main="3PLM ICCs")
dev.off()

