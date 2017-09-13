#####################################
# PLSC 504 -- Fall 2017
#
# Event Count materials - Day One
#####################################
# Options:

options(scipen=7)
options(digits=4)
par(mar=c(4,4,2,2))

# Packages:
library(RCurl)

# Various Poisson histograms

set.seed(7222009)
N<-1000
LP05<-rpois(N,0.5)
LP1<-rpois(N,1)
LP5<-rpois(N,5)
LP10<-rpois(N,10)

pdf("PoissonHistogramsR.pdf",7,6)
par(mfrow=c(2,2))
hist(LP05,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 0.5")
hist(LP1,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 1.0")
hist(LP5,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 5")
hist(LP10,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 10")
dev.off()

# Get SCOTUS nullifications data:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/nulls.csv")
Nulls<-read.csv(text=temp, header=TRUE)
rm(temp)

# Histogram:

pdf("NullsHist.pdf",6,5)
par(mar=c(4,4,2,2))
with(Nulls, 
     hist(nulls,main="",xlab="Number of Nullifications",
     col="grey"))
dev.off()

# Poisson regression:

nulls.poisson<-glm(nulls~tenure+unified,family="poisson",
                data=Nulls)
summary(nulls.poisson)

# IRRs:

library(mfx)
nulls.poisson.IRR<-poissonirr(nulls~tenure+unified,
                   data=Nulls)
nulls.poisson.IRR

# Predictions:

tenure<-seq(0,20,1)
unified<-1
simdata<-as.data.frame(cbind(tenure,unified))
nullhats<-predict(nulls.poisson,newdata=simdata,se.fit=TRUE)

# NOTE: These are XBs, not predicted counts.
# Transforming:

nullhats$Yhat<-exp(nullhats$fit)
nullhats$UB<-exp(nullhats$fit + 1.96*(nullhats$se.fit))
nullhats$LB<-exp(nullhats$fit - 1.96*(nullhats$se.fit))

# Plot...

pdf("NullsOutOfSampleHatsR.pdf",6,5)
plot(simdata$tenure,nullhats$Yhat,t="l",lwd=3,ylim=c(0,5),ylab=
           "Predicted Count", xlab="Mean Tenure")
lines(simdata$tenure,nullhats$UB,lwd=2,lty=2)
lines(simdata$tenure,nullhats$LB,lwd=2,lty=2)
dev.off()

# Offsets with dyadic data...Aggregated counts
# of conflicts between the countries in each
# dyad, 1950-1985...

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/offsetIR.csv")
IR<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(IR)

cor(IR,use="complete.obs")

IR.fit1<-glm(disputes~allies+openness,data=IR,family="poisson")
summary(IR.fit1)

IR.fit2<-glm(disputes~allies+openness,data=IR,family="poisson",
               offset=log(Ndyads))
summary(IR.fit2)

IR.fit3<-glm(disputes~allies+openness+log(Ndyads),data=IR,
             family="poisson")
summary(IR.fit3)

# z-test:
2*pnorm((0.811-1)/.071)

# Wald test:
wald.test(b=coef(IR.fit3),Sigma=vcov(IR.fit3),Terms=4,H0=1)

# CPB Max(Y) figure:

L.CPB <- seq(0.1,10,by=0.1)
MaxY8 <- (-L.CPB) / (0.8-1)
MaxY5 <- (-L.CPB) / (0.5-1)
MaxY2 <- (-L.CPB) / (0.2-1)

pdf("CPBMaxR.pdf",6,5)
par(mfrow=c(1,1))
par(mar=c(4,4,2,2))
plot(L.CPB,MaxY8,t="l",lwd=3,col="black",xlab="Lambda",
     ylab="Maximum Value of Y")
lines(L.CPB,MaxY5,lwd=3,col="red",lty=2)
lines(L.CPB,MaxY2,lwd=3,col="darkgreen",lty=3)
legend("topleft",bty="n",lty=c(1,2,3),lwd=3,
       col=c("black","red","darkgreen"),
       legend=c("Alpha = 0.8","Alpha = 0.5","Alpha = 0.2"))
dev.off()


# SCOTUS amici data:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/Amici.csv")
amici<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(amici)

amici.poisson<-glm(namici~term+civlibs,data=amici,
                   family="poisson")
summary(amici.poisson)

Phats<-fitted.values(amici.poisson)
Uhats<-((amici$namici-Phats)^2 - amici$namici) / (Phats * sqrt(2))
summary(lm(Uhats~Phats))

library(MASS)
amici.NB<-glm.nb(namici~term+civlibs,data=amici)
summary(amici.NB)
1 / amici.NB$theta

# Plot:

pdf("PoissonNBYHatsR.pdf",6,5)
par(mar=c(4,4,2,2))
plot(amici.poisson$fitted.values,amici.NB$fitted.values,pch=20,
     xlab="Poisson",ylab="Negative Binomial",main="")
abline(a=0,b=1,lwd=2)
dev.off()

