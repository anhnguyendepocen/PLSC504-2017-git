##########################################
# Code for PLSC 504 - Fall 2017
#
# Ordinal Response Models
#
##########################################
# Options:

options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

# Packages, etc.:

require(RCurl)

# Ordered sim:

set.seed(7222009)
X<-runif(1000,0,10)
Ystar<-0 + 1*X + rnorm(1000)
Y1<-Ystar
Y1[Ystar<2.5]<-1
Y1[Ystar>=2.5 & Ystar<5]<-2
Y1[Ystar>=5 & Ystar<7.5]<-3
Y1[Ystar>=7.5]<-4
table(Y1)

summary(lm(Ystar~X))
summary(lm(Y1~X))

pdf("OrdinalOneR.pdf",7,5)
par(mar=c(4,4,2,2))
par(mfrow=c(1,2))
plot(X,Ystar,pch=20,xlab="X",ylab="Y*")
abline(lm(Ystar~X),lwd=3,col="red")
abline(h=c(2.5,5,7.5),lty=2)
plot(X,Y1,pch=20,xlab="X",ylab="Y1")
abline(lm(Y1~X),lwd=3,col="red")
dev.off()

Y2<-Ystar
Y2[Ystar<2]<-1
Y2[Ystar>=2 & Ystar<8]<-2
Y2[Ystar>=8 & Ystar<9]<-3
Y2[Ystar>9]<-4
table(Y2)

summary(lm(Y2~X))

pdf("OrdinalTwoR.pdf",7,5)
par(mar=c(4,4,2,2))
par(mfrow=c(1,2))
plot(X,Ystar,pch=20,xlab="X",ylab="Y*")
abline(lm(Ystar~X),lwd=3,col="red")
abline(h=c(2,8,9),lty=2)
plot(X,Y2,pch=20,xlab="X",ylab="Y2")
abline(lm(Y2~X),lwd=3,col="red")
dev.off()

# Best Example Ever...

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/Beer.csv")
beer<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(beer)

library(MASS)
beer.logit<-polr(as.factor(quality)~price+calories+craftbeer
                 +bitter+malty,data=beer)
summary(beer.logit)

beer.probit<-polr(as.factor(quality)~price+calories+craftbeer+
                    bitter+malty,data=beer,method="probit")
summary(beer.probit)

# Odds Ratios

olreg.or <- function(model) { 
   coeffs <- coef(summary(beer.logit)) 
   lci <- exp(coeffs[ ,1] - 1.96 * coeffs[ ,2]) 
   or <- exp(coeffs[ ,1]) 
   uci <- exp(coeffs[ ,1] + 1.96 * coeffs[ ,2]) 
   lreg.or <- cbind(lci, or, uci) 
   lreg.or 
   } 

olreg.or(beer.logit)

# Predicted probs

calories<-seq(60,200,1)
price<-mean(beer$price)
craftbeer<-median(beer$craftbeer)
bitter<-mean(beer$bitter)
malty<-mean(beer$malty)
beersim<-cbind(calories,price,craftbeer,bitter,malty)
beer.hat<-predict(beer.logit,beersim,type='probs')

pdf("ROrdinalProbs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(c(60,200), c(0,1), type='n', xlab="Calories", ylab='Fitted 
  Probability')
lines(60:200, beer.hat[1:141, 1], lty=1, lwd=3)
lines(60:200, beer.hat[1:141, 2], lty=2, lwd=3)
lines(60:200, beer.hat[1:141, 3], lty=3, lwd=3)
lines(60:200, beer.hat[1:141, 4], lty=4, lwd=3)
dev.off()

# Cumulative probs:

xaxis<-c(60,60:200,200)
yaxis1<-c(0,beer.hat[,1],0)
yaxis2<-c(0,beer.hat[,2]+beer.hat[,1],0)
yaxis3<-c(0,beer.hat[,3]+beer.hat[,2]+beer.hat[,1],0)
yaxis4<-c(0,beer.hat[,4]+beer.hat[,3]+beer.hat[,2]+beer.hat[,1],0)

pdf("ROrdinalCumProbs.pdf",6,5)
par(mar=c(4,4,2,2))
plot(c(60,200), c(0,1), type='n', xlab="Calories", 
 ylab="Cumulative Probability")
polygon(xaxis,yaxis4,col="white")
polygon(xaxis,yaxis3,col="grey80")
polygon(xaxis,yaxis2,col="grey50")
polygon(xaxis,yaxis1,col="grey10")
dev.off()

# Variants: Series of Binary Regressions:

beer$goodplus<-as.factor(quality>1)
beer.good<-glm(goodplus~price+calories+craftbeer+bitter+
                 malty,family="binomial",data=beer)
summary(beer.good)

beer$VGplus<-as.factor(quality>2)
beer.VG<-glm(VGplus~price+calories+craftbeer+bitter+malty,
             family="binomial",data=beer)
summary(beer.VG)

