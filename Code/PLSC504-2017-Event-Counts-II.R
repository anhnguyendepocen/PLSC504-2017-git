#####################################
# PLSC 504 -- Fall 2017
#
# Event Count materials - Day Two
#####################################
# Options:

options(scipen=7)
options(digits=4)
par(mar=c(4,4,2,2))

# Packages:
library(RCurl)
library(mfx)
library(MASS)

###############################
# Hurdles, Zero-Inflation, etc.
#
# Conflict data:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/wars.csv")
wars<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(wars)

# Poisson:

wars.poisson<-glm(conflicts~polity+politysq+logPopulation+logGDP+
                    GDPGrowth+logOpenness+govshareGDP,family="poisson",
                  data=wars)
summary.glm(wars.poisson)

# Negative Binomial:

wars.nb<-glm.nb(conflicts~polity+politysq+logPopulation+logGDP+GDPGrowth+
                  logOpenness+govshareGDP,data=wars)
summary(wars.nb)

# Truncation E(Y) / Var(Y) plot:

lambda<-seq(0.1,6,by=0.1)
EY <- (lambda)/(1-exp(-lambda))
VY <- ((lambda)/(1-exp(-lambda))) * (1-((lambda)/(exp(lambda)-1)))

pdf("TruncatedPoissonEYVarYR.pdf",6,5)
plot(lambda,EY,t="l",lwd=2,lty=2,col="red",ylim=c(0,6),
     xlab="Lambda",ylab="E(Y|Y>0) and Var(Y|Y>0)")
lines(lambda,VY,lwd=2,lty=3,col="darkgreen")
abline(a=0,b=1,lwd=1)
legend("bottomright",bty="n",lwd=c(2,2,1),lty=c(2,3,1),
       legend=c("E(Y)","Var(Y)","Mean-Variance Equality"),
       col=c("red","darkgreen","black"))
dev.off()

# Truncation barplot:

outcomes<-seq(0,10,by=1)
outsPos<-seq(1,10,by=1)
upperT <- 4  # value at which dist. is upper-truncated...
outsUpper<-seq(0,upperT,by=1)
L <- 2.0
PoissonPr<- dpois(outcomes,L)
ZTPr <- (exp(-L)*(L^outsPos)) / ((factorial(outsPos)) * (1-exp(-L)))
ZTPr<-append(0,ZTPr)
T4Pr <- numeric(upperT+1)
PoisProb <- (exp(-L)*(L^outsUpper))/(factorial(outsUpper))
for (i in 1:(upperT+1)) {
  T4Pr[i] <- (exp(-L)*(L^outsUpper[i])) / 
    ((factorial(outsUpper[i])) * sum(PoisProb))
}
T4Pr<-append(T4Pr, rep(0,times=(10-upperT)))

df<-data.frame(Count = seq(0,10,by=1),
               Poisson = PoissonPr,
               ZeroTruncated = ZTPr,
               TruncAtFour = T4Pr)

pdf("PoissonTruncatedDensitiesR.pdf",7,5)
with(df, plot(Count,Poisson,pch=20,cex=1.2,col="black",lab=c(11,5,7),
              ylim=c(0,0.35),xlim=c(-0.3,10.3),ylab="Probability"))
with(df, points(Count-0.2,ZTPr,pch=4,col="red"))
with(df, points(Count+0.2,T4Pr,pch=17,col="darkgreen"))
with(df, segments(Count,0,Count,Poisson,lwd=1))
with(df, segments(Count-0.2,0,Count-0.2,ZTPr,lwd=1,col="red"))
with(df, segments(Count+0.2,0,Count+0.2,T4Pr,lwd=1,col="darkgreen"))
legend("topright",bty="n",col=c("black","red","darkgreen"),
       pch=c(20,4,17),lwd=1,lty=1,legend=c("Poisson","Truncated At Zero",
                                           "Upper Truncated at Four"))
dev.off()

# Truncation example. Standard Poisson (no zeros):

wars.poisNo0s<-glm(conflicts_no_zeros~polity+politysq+logPopulation+
                     logGDP+GDPGrowth+logOpenness+govshareGDP,
                   family="poisson",data=wars)

summary(wars.poisNo0s)

# Truncated at Zero:

library(VGAM)
wars.0tpois<-vglm(conflicts_no_zeros~polity+politysq+logPopulation+
                    logGDP+GDPGrowth+logOpenness+govshareGDP,
                  pospoisson,data=wars)

summary(wars.0tpois)

# Upper-censored Poisson:

wars$censoredconflicts<-wars$conflicts
wars$censoredconflicts<-ifelse(wars$conflicts>3,4,wars$censoredconflicts)
wars$censindicator<-ifelse(wars$censoredconflicts==4,1,0)

# Incorrect:

wars.poisCensored<-glm(censoredconflicts~polity+politysq+
                         logPopulation+logGDP+GDPGrowth+logOpenness+
                         govshareGDP,family="poisson",data=wars)

summary(wars.poisCensored)

# Censored:

wars.censpois<-vglm(SurvS4(censoredconflicts,censindicator)~polity+
                      politysq+logPopulation+logGDP+GDPGrowth+
                      logOpenness+govshareGDP,
                    cens.poisson,data=wars)
summary(wars.censpois)

# Zero-Inflated models, etc.

require(pscl)

wars.ZIP<-zeroinfl(conflicts~polity+politysq+logPopulation+
                     logGDP+GDPGrowth+logOpenness+govshareGDP,
                   data=wars,dist="poisson",link="logit")
summary(wars.ZIP)

wars.ZINB<-zeroinfl(conflicts~polity+politysq+logPopulation+
                      logGDP+GDPGrowth+logOpenness+govshareGDP,
                    data=wars,dist="negbin",link="logit")
summary(wars.ZINB)

wars.hurdle<-hurdle(conflicts~polity+politysq+logPopulation+
                      logGDP+GDPGrowth+logOpenness+govshareGDP,
                    data=wars,dist=c("poisson"),zero.dist=c("poisson"),
                    link=c("log"))
summary(wars.hurdle)

