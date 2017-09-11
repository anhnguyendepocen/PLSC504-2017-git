##########################################
# Code for PLSC 504 - Fall 2017
#
# Nominal Response Models
#
##########################################
# Options:

options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

# Packages, etc.:

require(RCurl)

###################################################
# MNL, etc.

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/Election1992small.csv")
nes92<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(nes92)

nes92.mlogit<-vglm(presvote~partyid, multinomial, nes92)
summary(nes92.mlogit)

Bush.nes92.mlogit<-vglm(formula=presvote~partyid, 
                        family=multinomial(refLevel=1),data=nes92) 
summary(Bush.nes92.mlogit)

Clinton.nes92.mlogit<-vglm(formula=presvote~partyid,
                           family=multinomial(refLevel=2),data=nes92)
summary(Clinton.nes92.mlogit)

# CL

library(mlogit)
colnames(nes92)<-c("caseid","presvote","partyid","FT.Bush","FT.Clinton","FT.Perot")
nes92$PVote<-factor(nes92$presvote,labels=c("Bush","Clinton","Perot"))
head(nes92)

nes92CL<-mlogit.data(nes92,shape="wide",choice="PVote",varying=4:6)
head(nes92CL,6)

# CL regression:

nes92.clogit<-mlogit(PVote~FT|partyid,data=nes92CL)
summary(nes92.clogit)

# Interpretation:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/Election1992.csv")
BigNES92<-read.csv(text=temp, header=TRUE)
rm(temp)


# library(VGAM)
NES.MNL<-vglm(presvote~partyid+age+white+female,data=BigNES92,
              multinomial(refLevel=1))
summaryvglm(NES.MNL)

library(aod)
wald.test(b=c(t(coef(NES.MNL))),Sigma=vcov(NES.MNL),Terms=c(5,6))
wald.test(b=c(t(coef(NES.MNL))),Sigma=vcov(NES.MNL),Terms=c(1,3,5,7,9))

# Hats

PickBush<-ifelse(fitted.values(NES.MNL)[,1]>fitted.values(NES.MNL)[,2] 
                 & fitted.values(NES.MNL)[,1]>fitted.values(NES.MNL)[,3], 1,0)
PickWJC<-ifelse(fitted.values(NES.MNL)[,2]>fitted.values(NES.MNL)[,1] 
                & fitted.values(NES.MNL)[,2]>fitted.values(NES.MNL)[,3], 2, 0)
PickHRP<-ifelse(fitted.values(NES.MNL)[,3]>fitted.values(NES.MNL)[,1] 
                & fitted.values(NES.MNL)[,3]>fitted.values(NES.MNL)[,2], 3, 0)
OutHat<-PickBush+PickWJC+PickHRP
table(BigNES92$presvote,OutHat)

# Odds ratios:

mnl.or <- function(model) { 
  coeffs <- c(t(coef(NES.MNL))) 
  lci <- exp(coeffs - 1.96 * diag(vcov(NES.MNL))^0.5) 
  or <- exp(coeffs) 
  uci <- exp(coeffs + 1.96* diag(vcov(NES.MNL))^0.5) 
  lreg.or <- cbind(lci, or, uci) 
  lreg.or 
} 

mnl.or(NES.MNL)

# In-sample predictions:

hats<-as.data.frame(fitted.values(NES.MNL))
names(hats)<-c("Bush","Clinton","Perot")
attach(hats)

library(car)

pdf("InSampleRScatterplotMatrix.pdf",8,7)
spm(~Bush+Clinton+Perot,pch=20,plot.points=TRUE,
    diagonal="histogram",col=c("black","grey"))
dev.off()

pdf("InSampleMNLPredProbsR.pdf",8,6)
par(mfrow=c(1,3))
plot(BigNES92$partyid,Bush,xlab="Party ID")
plot(BigNES92$partyid,Clinton,xlab="Party ID")
plot(BigNES92$partyid,Perot,xlab="Party ID")
par(mfrow=c(1,1))
dev.off()

# CL example:

nes92.clogit<-mlogit(PVote~FT|partyid,data=nes92CL)
summary(nes92.clogit)

# In-sample predictions:

CLhats<-predict(nes92.clogit,nes92CL)

pdf("InSampleCLHatsR.pdf",7,6)
plot(nes92$FT.Bush,CLhats[,1],pch=19,
     col=rgb(100,0,0,100,maxColorValue=255),
     xlab="Feeling Thermometer",
     ylab="Predicted Probability")
points(nes92$FT.Clinton+runif(nrow(CLhats),-1,1),
       CLhats[,2],pch=4,col=rgb(0,0,100,100,maxColorValue=255))
points(nes92$FT.Perot+runif(nrow(CLhats),-1,1),
       CLhats[,3],pch=17,col=rgb(0,100,0,50,maxColorValue=255))
lines(lowess(nes92$FT.Bush,CLhats[,1]),lwd=2,col="red")
lines(lowess(nes92$FT.Clinton,CLhats[,2]),lwd=2,col="blue")
lines(lowess(nes92$FT.Perot,CLhats[,3]),lwd=2,col="darkgreen")
legend("topleft",bty="n",c("Bush","Clinton","Perot"),
       col=c("red","blue","darkgreen"),pch=c(19,4,17))
dev.off()

# Multinomial Probit

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/Sweden2002.csv")
Sweden<-read.csv(text=temp, header=TRUE)
rm(temp)

library(mlogit)
Sweden.Long<-mlogit.data(Sweden,choice="partychoice",
                         shape="wide")
Sweden.MNL<-mlogit(partychoice~1|female+union+leftright+
                     age,data=Sweden.Long)
summary(Sweden.MNL)

# Restricted model (omitting Social Democrats)

Sweden.MNL.Restr<-mlogit(partychoice~1|female+union+leftright+age,
                         Sweden.Long,alt.subset=c("Conservatives","Liberals","Left Party"))

hmftest(Sweden.MNL,Sweden.MNL.Restr)

# Heteroscedastic logit:

Sweden.Het<-mlogit(partychoice~1|female+union+leftright+
                     age,data=Sweden.Long,heterosc=TRUE)
summary(Sweden.Het)

# MNP:

library(MNP)
Sweden.MNP<-mnp(partychoice~female+union+leftright+age,
                data=Sweden)
summary(Sweden.MNP)

# Alternative, using -mlogit-:
#
# Warning: This takes a long time...

Sweden.MNP2<-mlogit(partychoice~1|female+union+leftright+
                      age,data=Sweden.Long,probit=TRUE)
summary(Sweden.MNP2)

# Predictions:

MNL.hats<-data.frame(MNL.P = fitted(Sweden.MNL,outcome=FALSE))
foo<-predict(Sweden.MNP,type="prob")
MNP.hats<-data.frame(MNP.P=foo$p)
hats<-cbind(MNL.hats,MNP.hats)

pdf("MNL-MNP-ProbhatsR.pdf",7,7)
par(mfrow=c(2,2))
par(mar=c(4,4,4,2))
with(hats, plot(MNL.P.Conservatives,MNP.P.Conservatives,pch=19,
                xlab="MNL Predictions",ylab="MNP Predictions",
                main="Conservatives"))
with(hats, plot(MNL.P.Liberals,MNP.P.Liberals,pch=19,
                xlab="MNL Predictions",ylab="MNP Predictions",
                main="Liberals"))
with(hats, plot(MNL.P.Social.Democrats,MNP.P.Social.Democrats,pch=19,
                xlab="MNL Predictions",ylab="MNP Predictions",
                main="Social Democrats"))
with(hats, plot(MNL.P.Left.Party,MNP.P.Left.Party,pch=19,
                xlab="MNL Predictions",ylab="MNP Predictions",
                main="Left Party"))
dev.off()
