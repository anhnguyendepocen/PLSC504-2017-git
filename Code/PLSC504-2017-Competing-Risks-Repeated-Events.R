#########################################################
# PLSC 504 -- Fall 2017
#
# Competing Risks + Repeated Events
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


#########################################
# Competing Risks...

# SCOTUS data...

scotusURL<-"https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/scotus2.csv"
temp<-getURL(scotusURL)
scotus<-read.csv(textConnection(temp))
rm(scotusURL,temp)

scotus.S<-Surv(scotus$svcstart,scotus$service,scotus$retire)

# Illustrating competing risks: 

scotus.SR<-Surv(scotus$svcstart,scotus$service,scotus$retire)
scotus.SD<-Surv(scotus$svcstart,scotus$service,scotus$death)

# Plot:

pdf("CR-KMs.pdf",6,5)
plot(survfit(scotus.SR~1),mark.time=F,lwd=c(3,1,1),
     xlab="Time (in years)",ylab="Survival")
par(new=TRUE)
plot(survfit(scotus.SD~1),mark.time=F,lwd=c(4,2,2),
     lty=c(3,3,3),col=c("red","red","red"),
     xlab="Time (in years)",ylab="Survival")
legend("topright",bty="n",inset=0.02,lty=c(1,2),lwd=c(3,4),
       col=c("black","red"),c("Retirement","Death"))
dev.off()

scotus$C<-scotus$retire+scotus$death
scotus.SC<-Surv(scotus$svcstart,scotus$service,scotus$C)

scotus.C<-coxph(scotus.SC~age+chief+south+pension+pagree,
                data=scotus,method="efron")
scotus.R<-coxph(scotus.SR~age+chief+south+pension+pagree,
                data=scotus,method="efron")
scotus.D<-coxph(scotus.SD~age+chief+south+pension+pagree,
                data=scotus,method="efron")

# # Pretty table:
# scotus.C.out<-extract(scotus.C,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# scotus.R.out<-extract(scotus.R,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# scotus.D.out<-extract(scotus.D,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# texreg(l=list(scotus.C.out,scotus.R.out,scotus.D.out),
#        file="scotusCRTable.tex",
#        custom.model.names=c("Combined","Retirement","Death"),
#        custom.coef.names=c("Age","Chief","South",
#                            "Pension Eligibility","Party Agreement"),
#        stars=numeric(0),caption="",label="",table=F)

# MNL:

scotus$lnT<-log(scotus$service)
scotus.MNL<-multinom(threecat~chief+south+age+pension+pagree+lnT,
                     data=scotus,na.action=na.omit)

# # Pretty table:
# scotus.MNL.out<-extract(scotus.MNL,include.deviance=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# texreg(l=scotus.MNL.out,
#        file="scotusMNLTable.tex",       
#        custom.model.names=c("Retirement","Death"),
#        custom.coef.names=c("Intercept","Age","Chief","South",
#                            "Pension Eligibility","Party Agreement",
#                            "log(Time)"),
#        stars=numeric(0),caption="",label="",table=F)

##########################################
# Repeated Events...
#
# Back to the Oneal-Russett data...

ORURL<-"https://raw.githubusercontent.com/PrisonRodeo/PLSC504-2017-git/master/Data/OR.csv"
temp<-getURL(ORURL)
OR<-read.csv(textConnection(temp))
rm(ORURL,temp)


OR<-OR[order(OR$dyadid,OR$year),]
OR$one<-rep(1,times=nrow(OR))
OR<-ddply(OR,"dyadid",mutate,eventno=cumsum(dispute)-dispute+1,
          altstart=cumsum(one)-1,altstop=cumsum(one))

listvars<-c("dyadid","year","start","stop",
            "altstart","altstop","dispute","eventno")
ORsm<-OR[listvars]
print(ORsm[which(ORsm$dyadid==2130 & ORsm$year<1966),])

# First events:

OR1st<-OR[OR$eventno==1,]
OR.1st<-Surv(OR1st$altstart,OR1st$altstop,OR1st$dispute)
OR.Cox.1st<-coxph(OR.1st~allies+contig+capratio+growth+democracy+
                    trade+cluster(dyadid),data=OR1st,method="efron")

# AG:

OR.AGS<-Surv(OR$altstart,OR$altstop,OR$dispute)
OR.Cox.AG<-coxph(OR.AGS~allies+contig+capratio+growth+democracy+
                   trade+cluster(dyadid),data=OR,method="efron")

# PWP Elapsed:

OR.PWPES<-Surv(OR$altstart,OR$altstop,OR$dispute)
OR.Cox.PWPE<-coxph(OR.PWPES~allies+contig+capratio+growth+democracy+
                     trade+strata(eventno)+cluster(dyadid),data=OR,
                   method="efron")

# PWP Gap:

OR.PWPGS<-Surv(OR$start,OR$stop,OR$dispute)
OR.Cox.PWPG<-coxph(OR.PWPGS~allies+contig+capratio+growth+democracy+
                     trade+strata(eventno)+cluster(dyadid),data=OR,
                   method="efron")

# WLW:

OR.expand<-OR[rep(1:nrow(OR),each=max(OR$eventno)),]
OR.expand<-ddply(OR.expand,c("dyadid","year"),mutate,
                 eventrisk=cumsum(one))
OR.expand$dispute<-ifelse(OR.expand$eventno==OR.expand$eventrisk 
                          & OR.expand$dispute==1,1,0)

OR.expand.S<-Surv(OR.expand$altstart,OR.expand$altstop,
                  OR.expand$dispute)
OR.Cox.WLW<-coxph(OR.expand.S~allies+contig+capratio+growth+
                    democracy+trade+strata(eventno)+
                    cluster(dyadid),data=OR.expand,
                  method="efron")

# # Pretty table:
# 
# OR.1st.out<-extract(OR.Cox.1st,include.rsquared=F,include.maxrs=F,
#                       include.zph=F,include.nobs=F,
#                       include.missings=F)
# OR.AG.out<-extract(OR.Cox.AG,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.PWPE.out<-extract(OR.Cox.PWPE,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.PWPG.out<-extract(OR.Cox.PWPG,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# OR.WLW.out<-extract(OR.Cox.WLW,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
# 
# texreg(l=list(OR.1st.out,OR.AG.out,OR.PWPE.out,OR.PWPG.out,OR.WLW.out),
#        file="RepeatedEventsTable.tex",
#        custom.model.names=c("First","AG","PWP-E","PWP-G","WLW"),
#        custom.coef.names=c("Allies","Contiguity","Capability Ratio",
#                            "Growth","Democracy","Trade"),
#        stars=numeric(0),caption="",label="",table=F)

# Parameter change:

OR$capXevent<-OR$capratio*OR$eventno
OR.Cox.BVary<-coxph(OR.PWPGS~allies+contig+growth+democracy+
                      trade+capratio+capXevent+strata(eventno)+
                      cluster(dyadid),data=OR,
                      method="efron")

