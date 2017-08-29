##########################################
# Code for PLSC 504 - Fall 2017
#
# Binary Response Models, Part I
#
##########################################
# Packages, etc.:

require(RCurl)

# Options:

options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

###########################
# "Toy" example:

set.seed(7222009)
ystar<-rnorm(100)
y<-ifelse(ystar>0,1,0)
x<-ystar+(0.5*rnorm(100))
data<-data.frame(ystar,y,x)
head(data)

pdf("YstarYX-R.pdf",6,5)
par(mar=c(4,4,2,2))
plot(x,ystar,pch=19,ylab="Y* / Y",xlab="X")
points(x,y,pch=4,col="red")
abline(h=0)
legend("topleft",bty="n",pch=c(19,4),col=c("black","red"),
       legend=c("Y*","Y"))
dev.off()

# probits and logits...

myprobit<-glm(y~x,family=binomial(link="probit"),
              data=data)
summary(myprobit)

mylogit<-glm(y~x,family=binomial(link="logit"),
             data=data)
summary(mylogit)

pdf("LogitProbitHats.pdf",6,5)
plot(mylogit$fitted.values,myprobit$fitted.values,
     pch=20,xlab="Logit Predictions",
     ylab="Probit Predictions")
dev.off()

#################################
# NAFTA example...

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2016-git/master/Data/NAFTA.csv")
NAFTA<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(NAFTA)

# Logit:

NAFTA.GLM.fit<-glm(vote~democrat+pcthispc+cope93+DemXCOPE,
                   NAFTA,family=binomial)
summary(NAFTA.GLM.fit)

# Interactions...

NAFTA.GLM.fit$coeff[4]+NAFTA.GLM.fit$coeff[5]
(NAFTA.GLM.fit$coeff[4]+NAFTA.GLM.fit$coeff[5]) / 
  (sqrt(vcov(NAFTA.GLM.fit)[4,4] + 
  (1)^2*vcov(NAFTA.GLM.fit)[5,5] + 
  2*1*vcov(NAFTA.GLM.fit)[4,5]))

# Same thing, using -car-:

library(car)
linear.hypothesis(NAFTA.GLM.fit,"cope93+DemXCOPE=0")

# Predicted values:

preds<-NAFTA.GLM.fit$fitted.values
hats<-predict(NAFTA.GLM.fit,se.fit=TRUE)

# Plotting in-sample predictions:

XBUB<-hats$fit + (1.96*hats$se.fit) 
XBLB<-hats$fit - (1.96*hats$se.fit)
plotdata<-cbind(as.data.frame(hats),XBUB,XBLB)
plotdata<-data.frame(lapply(plotdata,binomial(link="logit")$linkinv))
par(mfrow=c(1,2))
library(plotrix)
plotCI(cope93[democrat==1],plotdata$fit[democrat==1],ui=plotdata$XBUB[democrat==1],
         li=plotdata$XBLB[democrat==1],pch=20,xlab="COPE Score",ylab="Predicted 
         Pr(Pro-NAFTA Vote)")
plotCI(cope93[democrat==0],plotdata$fit[democrat==0],ui=plotdata$XBUB[democrat==0],
         li=plotdata$XBLB[democrat==0],pch=20,xlab="COPE Score",ylab="Predicted 
         Pr(Pro-NAFTA Vote)")

# Plotting Out-of-sample Predictions:

sim.data<-data.frame(pcthispc=mean(nafta$pcthispc),democrat=rep(0:1,101),
                       cope93=seq(from=0,to=100,length.out=101))
sim.data$DemXCOPE<-sim.data$democrat*sim.data$cope93

OutHats<-predict(NAFTA.GLM.fit,se.fit=TRUE,newdata=sim.data)
OutHatsUB<-OutHats$fit+(1.96*OutHats$se.fit)
OutHatsLB<-OutHats$fit-(1.96*OutHats$se.fit)
OutHats<-cbind(as.data.frame(OutHats),OutHatsUB,OutHatsLB)
OutHats<-data.frame(lapply(OutHats,binomial(link="logit")$linkinv))

par(mfrow=c(1,2))
both<-cbind(sim.data,OutHats)
both<-both[order(both$cope93,both$democrat),]

plot(both$cope93[democrat==1],both$fit[democrat==1],t="l",lwd=2,ylim=c(0,1),
       xlab="COPE Score",ylab="Predicted Pr(Pro-NAFTA Vote)")
lines(both$cope93[democrat==1],both$OutHatsUB[democrat==1],lty=2)
lines(both$cope93[democrat==1],both$OutHatsLB[democrat==1],lty=2)
text(locator(1),label="Democrats")

plot(both$cope93[democrat==0],both$fit[democrat==0],t="l",lwd=2,ylim=c(0,1),
       xlab="COPE Score",ylab="Predicted Pr(Pro-NAFTA Vote)")
lines(both$cope93[democrat==0],both$OutHatsUB[democrat==0],lty=2)
lines(both$cope93[democrat==0],both$OutHatsLB[democrat==0],lty=2)
text(locator(1),label="Republicans")

# Odds Ratios:

lreg.or <- function(model)
       {
        coeffs <- coef(summary(NAFTA.GLM.fit))
        lci <- exp(coeffs[ ,1] - 1.96 * coeffs[ ,2])
        or <- exp(coeffs[ ,1])
        uci <- exp(coeffs[ ,1] + 1.96 * coeffs[ ,2])
        lreg.or <- cbind(lci, or, uci)        
        lreg.or
        }

lreg.or(NAFTA.GLM.fit)

####################
# Goodness of fit:

table(NAFTA.GLM.fit$fitted.values>0.5,nafta$vote==1)
chisq.test(NAFTA.GLM.fit$fitted.values>0.5,nafta$vote==1)

# ROC curves, plots, etc.:

library(ROCR)
NAFTA.GLM.logithats<-predict(NAFTA.GLM.fit,
                       type="response")
preds<-prediction(NAFTA.GLM.logithats,NAFTA$vote)
plot(performance(preds,"tpr","fpr"),lwd=2,lty=2,
       col="red",xlab="1 - Specificity",ylab="Sensitivity")
abline(a=0,b=1,lwd=3)




