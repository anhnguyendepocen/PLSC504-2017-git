#########################################################
# PLSC 504 -- Fall 2017
#
# Duration Dependence
#
########################################################
# Load packages (install as needed), set options:

library(RCurl)
library(foreign)
library(gtools)
library(plyr)
library(survival)
library(flexsurv)
library(nnet)
library(mstate)
library(texreg)

###########################################################
# Duration dependence...

# Unobserved heterogeneity plot:

pdf("UnobsHet.pdf",6,5)
par(mar=c(4,4,2,2))
plot(seq(0:20),rep(0.05,times=21),pch=NA,ylim=c(0.015,0.055),
     xlim=c(0,20),xlab="Time",ylab="Hazard")
abline(h=0.02,lwd=3,col="red")
abline(h=0.05,lwd=3,col="black")
abline(0.045,-0.001,lwd=3,lty=2,col="blue")
text(18,0.017,label=c("h(t) | Z=1"),col="red")
text(18,0.053,label=c("h(t) | Z=0"),col="black")
text(15.5,0.035,label=c("Estimated Hazard"),col="blue")
dev.off()

# Heterogeneity sim:

set.seed(7222009)
W<-rnorm(500)
X<-rnorm(500)
Z<-rnorm(500)
T<-rexp(500,rate=(exp(0+0.5*W+0.5*X-0.6*Z)))
C<-rep(1,times=500)
S<-Surv(T,C)

summary(survreg(S~W,dist="weibull"))
summary(survreg(S~W+X,dist="weibull"))
summary(survreg(S~W+X+Z,dist="weibull"))

M1<-survreg(S~W,dist="weibull")
M2<-survreg(S~W+X,dist="weibull")
M3<-survreg(S~W+X+Z,dist="weibull")

# Plot:

t<-cbind(1:60,1:60,1:60)
P<-c(1/(M1$scale),1/(M2$scale),1/(M3$scale))
DurDepHs<-t(apply(t,1,function(t) 0.02*P*((0.02*t)^(P-1))))

pdf("DurDepHs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t[,1],DurDepHs[,1],t="l",lwd=3,lty=1,col="green",
     xlab="Time",ylab="Weibull Hazard",ylim=c(0.01,0.04))
lines(t[,2],DurDepHs[,2],t="l",lwd=3,lty=2,col="blue")
lines(t[,3],DurDepHs[,3],t="l",lwd=3,lty=3,col="red")
abline(h=0.02,lty=4,lwd=2)
legend("topright",inset=.02,
       c("One Covariate","Two Covariates","Correct Specification","True Hazard"),
       lty=c(1,2,3,4),lwd=c(3,3,3,2),col=c("green","blue","red","black"),
       cex=1.2,bty="n")
dev.off()

# SCOTUS data...

scotusURL<-"https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/scotus2.csv"
temp<-getURL(scotusURL)
scotus<-read.csv(textConnection(temp))
rm(scotusURL,temp)

scotus.S<-Surv(scotus$svcstart,scotus$service,scotus$retire)

# Parameterized duration dependence:

ct.weib<-flexsurvreg(scotus.S~age+pension+pagree,
                     data=scotus,dist="weibull")
ct.weib.DD<-flexsurvreg(scotus.S~age+pension+pagree+shape(age),
                        data=scotus,dist="weibull")

# Plots (this code is a hotter mess than usual, and could use
# a solid smack with a few Hadleyverse tools):

t<-1:max(scotus$service)

Age<-seq(min(scotus$age),max(scotus$age),by=1)
P.vary<-exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*Age))

# Plot:

pdf("PbyAge.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Age,P.vary,t="l",lwd=3,lty=1,ylab="Estimate of p")
abline(h=1,lty=2,lwd=2)
dev.off()            

age55<-55
age75<-75

P<-c(exp(scotus.weib$coefficients[1]),
     exp(ct.weib$coefficients[1]),
     exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age55)),
     exp(ct.weib.DD$coefficients[1]+(ct.weib.DD$coefficients[6]*age75)))
XB55<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age55)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB75<-ct.weib$coefficients[2]+(ct.weib$coefficients[3]*age75)+
  (ct.weib$coefficients[4]*median(scotus$pension))+
  (ct.weib$coefficients[5]*median(scotus$pagree)) 
XB55DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age55)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB75DD<-ct.weib.DD$coefficients[2]+(ct.weib.DD$coefficients[3]*age75)+
  (ct.weib.DD$coefficients[4]*median(scotus$pension))+
  (ct.weib.DD$coefficients[5]*median(scotus$pagree)) 
XB<-c(XB55,XB75,XB55DD,XB75DD)

h55<-dweibull(t,P[1],scale=XB[1]) / (1 - pweibull(t,P[1],scale=XB[1]))
h75<-dweibull(t,P[2],scale=XB[2]) / (1 - pweibull(t,P[2],scale=XB[2]))
h55DD<-dweibull(t,P[3],scale=XB[3]) / (1 - pweibull(t,P[3],scale=XB[3]))
h75DD<-dweibull(t,P[4],scale=XB[4]) / (1 - pweibull(t,P[4],scale=XB[4]))

# Plot:

pdf("DDHazards.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t,h75DD,t="l",lwd=3,lty=1,col="red",ylim=c(0.2,0.8),
     xlab="Time (in years)",ylab="Hazard")
lines(t,h55DD,lwd=3,lty=1,col="blue")
lines(t,h75,lwd=3,lty=2,col="red")
lines(t,h55,lwd=3,lty=2,col="blue")
legend("topleft",inset=0.02,
       c("Age 55 (p varying)","Age 75 (p varying)","Age 55 (p fixed)",
         "Age 75 (p fixed)"),lty=c(1,1,2,2),lwd=c(3,3,3,3),
       col=c("blue","red","blue","red"),bty="n")
dev.off()

