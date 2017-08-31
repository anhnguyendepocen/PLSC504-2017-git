##########################################
# Code for PLSC 504 - Fall 2017
#
# Binary Response Models, Part II
#
##########################################
# Packages, etc.:

require(RCurl)

# Options:

options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places


####################################
# Separation:

# Table

Yeas<-t(c(rep(0,times=212),rep(1,times=219)))
Dems<-t(c(rep(0,times=178),rep(1,times=253)))
table(Yeas,Dems)

# Simulated Logits:

set.seed(7222009)
X<-runif(100,min=-5,max=5)
X<-X[order(X)]
Z<-runif(100,min=-5,max=5)
Y<-ifelse(plogis(X+Z)>0.5,1,0)
Y2<-ifelse(plogis(X+0.5*Z)>0.5,1,0)
Y3<-ifelse(plogis(X+0.1*Z)>0.5,1,0)
Ysep<-ifelse(plogis(X)>0.5,1,0)
Yfit<-glm(Y~X,family="binomial")
Y2fit<-glm(Y2~X,family="binomial")
Y3fit<-glm(Y3~X,family="binomial")
Ysepfit<-glm(Ysep~X,family="binomial")

# Plots:

pdf("Separation.pdf",8,7)
par(mar=c(4,4,2,2))
par(mfrow=c(2,2))
plot(X,Y,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Yfit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Yfit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Yfit))[4],digits=2))))
plot(X,Y2,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Y2fit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Y2fit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Y2fit))[4],digits=2))))
plot(X,Y3,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Y3fit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Y3fit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Y3fit))[4],digits=2))))
plot(X,Ysep,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Ysepfit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Ysepfit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Ysepfit))[4],digits=2))))
dev.off()

# Toy data:

rm(X,Y,Z)
set.seed(7222009)
Z<-rnorm(500)
W<-rnorm(500)
Y<-rbinom(500,size=1,prob=plogis((0.2+0.5*W-0.5*Z)))
X<-rbinom(500,1,0.5)
X<-ifelse(Y==0,0,X)

summary(glm(Y~W+Z+X,family="binomial"))

data<-as.data.frame(cbind(W,X,Y,Z))
write.dta(data,"SepSim.dta") # for the Stata illustration

# Pets data:

PetsURL<-"https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2016-git/master/Data/Pets.csv"
temp<-getURL(PetsURL)
Pets<-read.csv(textConnection(temp))
rm(temp)

Pets.1<-glm(petfamily~female+as.factor(married)+as.factor(partyid)
            +as.factor(education),data=Pets,family=binomial)

Pets.2<-glm(petfamily~female+as.factor(married)*female+as.factor(partyid)+
              as.factor(education),data=Pets,family=binomial)

with(Pets, xtabs(~petfamily+as.factor(married)+female))

Pets.Firth<-logistf(petfamily~female+
                      as.factor(married)*female+as.factor(partyid)+
                      as.factor(education),data=Pets)


