# Code for SPRINT clinical decision score replication
# Sanjay Basu, basus@stanford.edu

# INSTRUCTIONS: set all working directory calls (setwd) and load/save commands to appropriate directories on your local machine, to load/save data from/to your desired location. 
# This code does not contain the data itself, which requires IRB and NIH-BioLINCC approval, but is built to analyze the NIH-BioLINCC versions of the SPRINT-Pop and ACCORD-BP datasets.
# Note that we do not use the composite Framingham risk score included in the SPRINT-Pop dataset, as we found this to be erroneously calculated in the original release of the SPRINT-Pop dataset.

setwd("~/Data/sprint_pop/data")

#### load packages ####

install.packages('sas7bdat')
install.packages('gdata')
install.packages('glmnet') # https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
install.packages('Hmisc')
install.packages('survival')
install.packages('matrixStats')
install.packages('doMC')
install.packages('doBy')
install.packages('cvAUC')
install.packages('survivalROC')
install.packages('psych')
install.packages('metafor')
install.packages('survminer')
install.packages('mvmeta')
install.packages('survMisc')
install.packages('RTCGA.clinical')
install.packages('survIDINRI')
install.packages('corrgram')

library(sas7bdat)
library(gdata)
library(glmnet) # https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
library(Hmisc)
library(survival)
library(matrixStats)
library(doMC)
library(doBy)
library(cvAUC)
library(survivalROC)
library(psych)
library(metafor)
library(survminer)
library(mvmeta)
library(survMisc)
library(RTCGA.clinical)
library(survIDINRI)
library(corrgram)
registerDoMC(cores=4)

#### merge sprintpop data ####
rm(list=ls())
baseline = read.sas7bdat("baseline.sas7bdat")
bp = read.sas7bdat("bp.sas7bdat")
outcomes = read.sas7bdat("outcomes.sas7bdat")
retention = read.sas7bdat("retention.sas7bdat")
safety = read.sas7bdat("safety.sas7bdat")
bp_cut = bp[which(bp$VISITCODE=="RZ"),]
sprint_set = merge(baseline,bp_cut,by="MASKID")
sprint_set = merge(sprint_set,outcomes,by="MASKID")
sprint_set = merge(sprint_set,retention,by="MASKID")
sprint_set = merge(sprint_set,safety,by="MASKID")
save.image("sprint_cut.RData")

##### load derivation data ####
rm(list=ls())
load("sprint_cut.RData")
keep(sprint_set, sure=TRUE)
attach(sprint_set)

#### gen MI/stroke/CVD death model ####
kmdec=function(dec.num,dec.name, datain, adm.cens){
  stopped=0
  data.sub=datain[datain[,dec.name]==dec.num,]
  if (sum(data.sub$out)>1){
    avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
    avsurv.est=ifelse(min(avsurv$time)<=adm.cens,avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],1)
    avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    avsurv.stderr=avsurv.stderr*avsurv.est
    avsurv.num=ifelse(min(avsurv$time)<=adm.cens,avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
  } else {
    return(c(0,0,0,0,stopped=-1))
  }
  if (sum(data.sub$out)<5) stopped=1
  c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped) 
}
GND.calib = function(pred, tvar, out, cens.t, groups, adm.cens){
  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t=ifelse(tvar>adm.cens, 0, out)
  datause=data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec))
  groups=sort(unique(datause$dec))
  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\
                                        (see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
                                        Consider collapsing some groups to avoid this problem.")
  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1], 
                   kmvar=kmtab[,2]^2, 
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))
  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  print(hltab[c(1,2,3,4,10,5,6,9,7,11)], digits=4)
  plot(tapply(datause$pred,datause$dec,mean),1-kmtab[,1],xlab="Expected K-M rate",ylab="Observed K-M rate",xlim=c(0,1),ylim=c(0,1))
  abline(a=0,b=1, col = "gray60")
  calline = lm(hltab$kmperc~hltab$expectedperc)
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1),slope=calline$coefficients[2],intercept = calline$coefficients[1])
}
hisp = (RACE4=="HISPANIC")
currentsmoker = (SMOKE_3CAT==3)
formersmoker = (SMOKE_3CAT==2)
cvd = (EVENT_MI==1)|(EVENT_STROKE==1)|(EVENT_CVDDEATH==1)|(EVENT_HF==1)
t_censor = rowMaxs(cbind(T_MI,T_STROKE,T_CVDDEATH,T_HF))
t_cvds = rowMaxs(cbind(T_MI*EVENT_MI,T_STROKE*EVENT_STROKE,T_CVDDEATH*EVENT_CVDDEATH,T_HF*EVENT_HF))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)
sae = (HYP_SAE_EVNT==1)|(SYN_SAE_EVNT==1)|(ELE_SAE_EVNT==1)|(AKI_SAE_EVNT==1)|(BRA_SAE_EVNT==1)|(INJ_SAE_EVNT==1)
t_censor = rowMaxs(cbind(HYP_SAE_DAYS,SYN_SAE_DAYS,ELE_SAE_DAYS,AKI_SAE_DAYS,BRA_SAE_DAYS,INJ_SAE_DAYS)) 
t_saes = rowMaxs(cbind(HYP_SAE_DAYS*HYP_SAE_EVNT,SYN_SAE_DAYS*SYN_SAE_EVNT,ELE_SAE_DAYS*ELE_SAE_EVNT,AKI_SAE_DAYS*AKI_SAE_EVNT,BRA_SAE_DAYS*BRA_SAE_EVNT,INJ_SAE_DAYS*INJ_SAE_EVNT)) 
t_saes[t_saes==0] = t_censor[t_saes==0]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
testsubset = data.frame(cOutcome,
                        INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        SCREAT,CHR,HDL,TRR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]
cvdmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox", parallel = TRUE)
plot(cvdmodel.cv)
coef.min = coef(cvdmodel.cv, s = 'lambda.min')
coef.min
c<-data.frame(cvd,t_cvds,sae,t_saes,
              INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
              SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
              ASPIRIN,STATIN,
              SCREAT,CHR,HDL,TRR,BMI,
              INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
              INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
              INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
              INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
c=c[complete.cases(c),]
adm.cens=5*365.25
c$fu.time <- pmin(c$t_cvds, adm.cens)
c$status <- ifelse(as.numeric(adm.cens < c$t_cvds), 0, c$cvd)



calt = c
survcox_c<-coxph(data=c, Surv(fu.time, status)~AGE+FEMALE+RACE_BLACK+hisp+SBP.y+DBP.y+
                   N_AGENTS+currentsmoker+formersmoker+ASPIRIN+STATIN+SCREAT+CHR+
                   HDL+TRR+BMI+INTENSIVE*AGE+INTENSIVE*RACE_BLACK+
                   INTENSIVE*DBP.y+INTENSIVE*currentsmoker+INTENSIVE*HDL+
                   INTENSIVE*TRR)
summary(survcox_c)
survfit_c=survfit(survcox_c, newdata=c, se.fit=FALSE)
estinc_c=1-survfit_c$surv[dim(survfit_c$surv)[1],]
c$dec=as.numeric(cut2(estinc_c, g=10))
GND.result=GND.calib(pred=estinc_c, tvar=c$fu.time, out=c$status, 
                     cens.t=adm.cens, groups=c$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_c,c$cvd)



#### gen adverse events outcome model ####
testsubset = data.frame(dOutcome,
                        INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        SCREAT,CHR,HDL,TRR,BMI,UMALCR,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI,INTENSIVE*UMALCR)
testsubset=testsubset[complete.cases(testsubset),]
saemodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox", parallel = TRUE)
plot(saemodel.cv)
coef.min = coef(saemodel.cv, s = 'lambda.min')
coef.min
d<-c
adm.cens=5*365.25
d$fu.time <- pmin(d$t_saes, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_saes), 0, d$sae)
dalt = d

survcox_d<-coxph(data=d, Surv(fu.time, status)~AGE+FEMALE+hisp+SBP.y+DBP.y+
                   N_AGENTS+currentsmoker+formersmoker+ASPIRIN+SCREAT+CHR+HDL+
                   INTENSIVE*FEMALE+INTENSIVE*currentsmoker+INTENSIVE*STATIN+INTENSIVE*SCREAT+
                   INTENSIVE*CHR+INTENSIVE*TRR)
summary(survcox_d)
survfit_d=survfit(survcox_d, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$sae)



#### C/GDN ####
INTENSIVE = sprint_set$INTENSIVE
INTENSIVE[INTENSIVE==0]=1
cint<-data.frame(cvd,t_cvds,sae,t_saes,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 SCREAT,CHR,HDL,TRR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
cint=cint[complete.cases(cint),]
adm.cens=5*365.25
cint$fu.time <- pmin(cint$t_cvds, adm.cens)
cint$status <- ifelse(as.numeric(adm.cens < cint$t_cvds), 0, cint$cvd)
survfit_cint=survfit(survcox_c, newdata=cint, se.fit=FALSE)
estinc_cint=1-survfit_cint$surv[dim(survfit_cint$surv)[1],]
INTENSIVE = sprint_set$INTENSIVE
INTENSIVE[INTENSIVE==1]=0
cstd<-data.frame(cvd,t_cvds,sae,t_saes,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 SCREAT,CHR,HDL,TRR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
cstd=cstd[complete.cases(cstd),]
adm.cens=5*365.25
cstd$fu.time <- pmin(cstd$t_cvds, adm.cens)
cstd$status <- ifelse(as.numeric(adm.cens < cstd$t_cvds), 0, cstd$cvd)
survfit_cstd=survfit(survcox_c, newdata=cstd, se.fit=FALSE)
estinc_cstd=1-survfit_cstd$surv[dim(survfit_cstd$surv)[1],]
benefit = estinc_cstd-estinc_cint
hist(benefit,xlab="Predicted reduction in CVD events/CVD death \n from intensive treatment (probability)")
survfit_dint=survfit(survcox_d, newdata=cint, se.fit=FALSE)
estinc_dint=1-survfit_dint$surv[dim(survfit_dint$surv)[1],]
survfit_dstd=survfit(survcox_d, newdata=cstd, se.fit=FALSE)
estinc_dstd=1-survfit_dstd$surv[dim(survfit_dstd$surv)[1],]
risk = estinc_dint-estinc_dstd
hist(risk,xlab="Predicted increase in serious adverse events \n from intensive treatment (probability)")
plot(benefit,risk, xlim=c(-.1,.5),ylim=c(-.1,.5),xlab="Predicted absolute risk reduction in CVD events/death", ylab="Predicted absolute risk increase in serious adverse events", col=rgb(0,0,100,50,maxColorValue=255), pch=19)
abline(a=0,b=0, col = "gray60")
abline(v=0, col = "gray60")

sprintbenefit = benefit
sprintrisk = risk

INTENSIVE = sprint_set$INTENSIVE
pv = c
pv$fu.time <- pmin(c$t_cvds, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < c$t_cvds), 0, c$cvd)

betax=(survcox_c$coefficients[1]*pv$AGE+
         survcox_c$coefficients[2]*pv$FEMALE+
         survcox_c$coefficients[3]*pv$RACE_BLACK+
         survcox_c$coefficients[4]*pv$hisp+
         survcox_c$coefficients[5]*pv$SBP.y+
         survcox_c$coefficients[6]*pv$DBP.y+
         survcox_c$coefficients[7]*pv$N_AGENTS+
         survcox_c$coefficients[8]*pv$currentsmoker+
         survcox_c$coefficients[9]*pv$formersmoker+
         survcox_c$coefficients[10]*pv$ASPIRIN+
         survcox_c$coefficients[11]*pv$STATIN+
         survcox_c$coefficients[12]*pv$SCREAT+
         survcox_c$coefficients[13]*pv$CHR+
         survcox_c$coefficients[14]*pv$HDL+
         survcox_c$coefficients[15]*pv$TRR+
         survcox_c$coefficients[16]*pv$BMI+
         survcox_c$coefficients[17]*pv$INTENSIVE+
         survcox_c$coefficients[18]*pv$INTENSIVE*pv$AGE+
         survcox_c$coefficients[19]*pv$INTENSIVE*pv$RACE_BLACK+
         survcox_c$coefficients[20]*pv$INTENSIVE*pv$DBP.y+
         survcox_c$coefficients[21]*pv$INTENSIVE*pv$currentsmoker+
         survcox_c$coefficients[22]*pv$INTENSIVE*pv$HDL+
         survcox_c$coefficients[23]*pv$INTENSIVE*pv$TRR)
cvdpred = 1 - .943^exp(betax-mean(na.omit(betax)))
estinc_e=cvdpred
pv$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=pv$fu.time, out=pv$status, 
                     cens.t=adm.cens, groups=pv$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,pv$cvd)

benpredbetaint = (survcox_c$coefficients[1]*pv$AGE+
                    survcox_c$coefficients[2]*pv$FEMALE+
                    survcox_c$coefficients[3]*pv$RACE_BLACK+
                    survcox_c$coefficients[4]*pv$hisp+
                    survcox_c$coefficients[5]*pv$SBP.y+
                    survcox_c$coefficients[6]*pv$DBP.y+
                    survcox_c$coefficients[7]*pv$N_AGENTS+
                    survcox_c$coefficients[8]*pv$currentsmoker+
                    survcox_c$coefficients[9]*pv$formersmoker+
                    survcox_c$coefficients[10]*pv$ASPIRIN+
                    survcox_c$coefficients[11]*pv$STATIN+
                    survcox_c$coefficients[12]*pv$SCREAT+
                    survcox_c$coefficients[13]*pv$CHR+
                    survcox_c$coefficients[14]*pv$HDL+
                    survcox_c$coefficients[15]*pv$TRR+
                    survcox_c$coefficients[16]*pv$BMI+
                    survcox_c$coefficients[17]*1+
                    survcox_c$coefficients[18]*1*pv$AGE+
                    survcox_c$coefficients[19]*1*pv$RACE_BLACK+
                    survcox_c$coefficients[20]*1*pv$DBP.y+
                    survcox_c$coefficients[21]*1*pv$currentsmoker+
                    survcox_c$coefficients[22]*1*pv$HDL+
                    survcox_c$coefficients[23]*1*pv$TRR)
benintpred = 1 - .943^exp(benpredbetaint-mean(na.omit(betax)))

benpredbetastd = (survcox_c$coefficients[1]*pv$AGE+
                    survcox_c$coefficients[2]*pv$FEMALE+
                    survcox_c$coefficients[3]*pv$RACE_BLACK+
                    survcox_c$coefficients[4]*pv$hisp+
                    survcox_c$coefficients[5]*pv$SBP.y+
                    survcox_c$coefficients[6]*pv$DBP.y+
                    survcox_c$coefficients[7]*pv$N_AGENTS+
                    survcox_c$coefficients[8]*pv$currentsmoker+
                    survcox_c$coefficients[9]*pv$formersmoker+
                    survcox_c$coefficients[10]*pv$ASPIRIN+
                    survcox_c$coefficients[11]*pv$STATIN+
                    survcox_c$coefficients[12]*pv$SCREAT+
                    survcox_c$coefficients[13]*pv$CHR+
                    survcox_c$coefficients[14]*pv$HDL+
                    survcox_c$coefficients[15]*pv$TRR+
                    survcox_c$coefficients[16]*pv$BMI+
                    survcox_c$coefficients[17]*0+
                    survcox_c$coefficients[18]*0*pv$AGE+
                    survcox_c$coefficients[19]*0*pv$RACE_BLACK+
                    survcox_c$coefficients[20]*0*pv$DBP.y+
                    survcox_c$coefficients[21]*0*pv$currentsmoker+
                    survcox_c$coefficients[22]*0*pv$HDL+
                    survcox_c$coefficients[23]*0*pv$TRR)
benstdpred = 1 - .943^exp(benpredbetastd-mean(na.omit(betax)))

netben = benstdpred-benintpred 
hist(netben)






pv$fu.time <- pmin(c$t_saes, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < c$t_saes), 0, c$sae)
betax=(survcox_d$coefficients[1]*pv$AGE+
         survcox_d$coefficients[2]*pv$FEMALE+
         survcox_d$coefficients[3]*pv$hisp+
         survcox_d$coefficients[4]*pv$SBP.y+
         survcox_d$coefficients[5]*pv$DBP.y+
         survcox_d$coefficients[6]*pv$N_AGENTS+
         survcox_d$coefficients[7]*pv$currentsmoker+
         survcox_d$coefficients[8]*pv$formersmoker+
         survcox_d$coefficients[9]*pv$ASPIRIN+
         survcox_d$coefficients[10]*pv$SCREAT+
         survcox_d$coefficients[11]*pv$CHR+
         survcox_d$coefficients[12]*pv$HDL+
         survcox_d$coefficients[13]*pv$INTENSIVE+
         survcox_d$coefficients[14]*pv$STATIN+
         survcox_d$coefficients[15]*pv$TRR+
         survcox_d$coefficients[16]*pv$FEMALE*pv$INTENSIVE+
         survcox_d$coefficients[17]*pv$INTENSIVE*pv$currentsmoker+
         survcox_d$coefficients[18]*pv$STATIN*pv$INTENSIVE+
         survcox_d$coefficients[19]*pv$INTENSIVE*pv$SCREAT+
         survcox_d$coefficients[20]*pv$INTENSIVE*pv$CHR+
         survcox_d$coefficients[21]*pv$INTENSIVE*pv$TRR)
saepred = 1 - .897^exp(betax-mean(na.omit(betax)))
estinc_e=saepred
pv$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=pv$fu.time, out=pv$status, 
                     cens.t=adm.cens, groups=pv$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,pv$sae)

harmpredbetaint = (survcox_d$coefficients[1]*pv$AGE+
                     survcox_d$coefficients[2]*pv$FEMALE+
                     survcox_d$coefficients[3]*pv$hisp+
                     survcox_d$coefficients[4]*pv$SBP.y+
                     survcox_d$coefficients[5]*pv$DBP.y+
                     survcox_d$coefficients[6]*pv$N_AGENTS+
                     survcox_d$coefficients[7]*pv$currentsmoker+
                     survcox_d$coefficients[8]*pv$formersmoker+
                     survcox_d$coefficients[9]*pv$ASPIRIN+
                     survcox_d$coefficients[10]*pv$SCREAT+
                     survcox_d$coefficients[11]*pv$CHR+
                     survcox_d$coefficients[12]*pv$HDL+
                     survcox_d$coefficients[13]*1+
                     survcox_d$coefficients[14]*pv$STATIN+
                     survcox_d$coefficients[15]*pv$TRR+
                     survcox_d$coefficients[16]*pv$FEMALE*1+
                     survcox_d$coefficients[17]*1*pv$currentsmoker+
                     survcox_d$coefficients[18]*pv$STATIN*1+
                     survcox_d$coefficients[19]*1*pv$SCREAT+
                     survcox_d$coefficients[20]*1*pv$CHR+
                     survcox_d$coefficients[21]*1*pv$TRR)
harmintpred = 1 - .897^exp(harmpredbetaint-mean(na.omit(betax)))

harmpredbetastd = (survcox_d$coefficients[1]*pv$AGE+
                     survcox_d$coefficients[2]*pv$FEMALE+
                     survcox_d$coefficients[3]*pv$hisp+
                     survcox_d$coefficients[4]*pv$SBP.y+
                     survcox_d$coefficients[5]*pv$DBP.y+
                     survcox_d$coefficients[6]*pv$N_AGENTS+
                     survcox_d$coefficients[7]*pv$currentsmoker+
                     survcox_d$coefficients[8]*pv$formersmoker+
                     survcox_d$coefficients[9]*pv$ASPIRIN+
                     survcox_d$coefficients[10]*pv$SCREAT+
                     survcox_d$coefficients[11]*pv$CHR+
                     survcox_d$coefficients[12]*pv$HDL+
                     survcox_d$coefficients[13]*0+
                     survcox_d$coefficients[14]*pv$STATIN+
                     survcox_d$coefficients[15]*pv$TRR+
                     survcox_d$coefficients[16]*pv$FEMALE*0+
                     survcox_d$coefficients[17]*0*pv$currentsmoker+
                     survcox_d$coefficients[18]*pv$STATIN*0+
                     survcox_d$coefficients[19]*0*pv$SCREAT+
                     survcox_d$coefficients[20]*0*pv$CHR+
                     survcox_d$coefficients[21]*0*pv$TRR)
harmstdpred = 1 - .897^exp(harmpredbetastd-mean(na.omit(betax)))

netharm = harmintpred-harmstdpred 
hist(netharm)

netbensprint = netben
netharmsprint = netharm


plot(netben,netharm, xlim=c(-.2,.4),ylim=c(-.2,.4),xlab="Predicted reduction in CVD events/death", ylab="Predicted increase in serious adverse events")
abline(a=0,b=0, col = "gray60")
abline(v=0, col = "gray60")




#### score vs SPRINT outcomes ####
bencats = c(.01,.03)
bencat = 1*(netben<bencats[1])+
  2*((netben>=bencats[1])&(netben<bencats[2]))+
  3*((netben>=bencats[2]))
cvdtable = describeBy(pv$cvd,list(bencat,pv$INTENSIVE),mat=TRUE)
cvdtable
prop.test(x=c(cvdtable[1,5]*cvdtable[1,6],cvdtable[4,5]*cvdtable[4,6]), n=c(cvdtable[1,5],cvdtable[4,5]), correct=FALSE)
prop.test(x=c(cvdtable[2,5]*cvdtable[2,6],cvdtable[5,5]*cvdtable[5,6]), n=c(cvdtable[2,5],cvdtable[5,5]), correct=FALSE)
prop.test(x=c(cvdtable[3,5]*cvdtable[3,6],cvdtable[6,5]*cvdtable[6,6]), n=c(cvdtable[3,5],cvdtable[6,5]), correct=FALSE)
cvdtabletimes = describeBy(pv$t_cvds,list(bencat,pv$INTENSIVE),mat=TRUE)
x1i = cvdtable[1:3,5]*cvdtable[1:3,6]
x2i = cvdtable[4:6,5]*cvdtable[4:6,6]
t1i = cvdtabletimes[1:3,5]*cvdtabletimes[1:3,6]
t2i = cvdtabletimes[4:6,5]*cvdtabletimes[4:6,6]
rma(measure="IR", x1i, x2i, t1i, t2i, method="REML")

harmcats = c(.005,.04)
harmcat = 1*(netharm<harmcats[1])+
  2*((netharm>=harmcats[1])&(netharm<harmcats[2]))+
  3*((netharm>=harmcats[2]))
saetable = describeBy(pv$sae,list(harmcat,pv$INTENSIVE),mat=TRUE)
saetable
prop.test(x=c(saetable[1,5]*saetable[1,6],saetable[4,5]*saetable[4,6]), n=c(saetable[1,5],saetable[4,5]), correct=FALSE)
prop.test(x=c(saetable[2,5]*saetable[2,6],saetable[5,5]*saetable[5,6]), n=c(saetable[2,5],saetable[5,5]), correct=FALSE)
prop.test(x=c(saetable[3,5]*saetable[3,6],saetable[6,5]*saetable[6,6]), n=c(saetable[3,5],saetable[6,5]), correct=FALSE)
saetabletimes = describeBy(pv$t_saes,list(harmcat,pv$INTENSIVE),mat=TRUE)
x1i = saetable[1:3,5]*saetable[1:3,6]
x2i = saetable[4:6,5]*saetable[4:6,6]
t1i = saetabletimes[1:3,5]*saetabletimes[1:3,6]
t2i = saetabletimes[4:6,5]*saetabletimes[4:6,6]
rma(measure="IR", x1i, x2i, t1i, t2i, method="REML")




save.image("~/Data/sprint_pop/data/sprint_run.RData")

#### merge accord data ####
setwd("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets")
accord_key = read.sas7bdat("accord_key.sas7bdat")
accord_key = accord_key[(accord_key$arm>=1)&(accord_key$arm<=4),]
bloodpressure = read.sas7bdat("bloodpressure.sas7bdat")
concomitantmeds = read.sas7bdat("concomitantmeds.sas7bdat")
cvdoutcomes = read.sas7bdat("cvdoutcomes.sas7bdat")
lipids = read.sas7bdat("lipids.sas7bdat")
otherlabs = read.sas7bdat("otherlabs.sas7bdat")
setwd("~/Data/accord/4-Data_Sets-CRFs/4a-CRF_Data_Sets")
f13_intensivebpmanagement = read.sas7bdat("f13_intensivebpmanagement.sas7bdat")
f14_standardbpmanagement = read.sas7bdat("f14_standardbpmanagement.sas7bdat")
f07_baselinehistoryphysicalexam = read.sas7bdat("f07_baselinehistoryphysicalexam.sas7bdat")
accord_key_cut = accord_key[which((accord_key$arm==1)|(accord_key$arm==2)|(accord_key$arm==3)|(accord_key$arm==4)),]
accord_key_cut$INTENSIVE = (accord_key_cut$arm==1)|(accord_key_cut$arm==3)
bloodpressure_cut = bloodpressure[which(bloodpressure$Visit=="BLR"),]
concomitantmeds_cut = concomitantmeds[which(concomitantmeds$Visit=="BLR"),]
f13_intensivebpmanagement_cut = f13_intensivebpmanagement[which(f13_intensivebpmanagement$Visit=="BLR"|
                                                                  f13_intensivebpmanagement$Visit=="F01"|
                                                                  f13_intensivebpmanagement$Visit=="F02"|
                                                                  f13_intensivebpmanagement$Visit=="F03"|
                                                                  f13_intensivebpmanagement$Visit=="F04"|
                                                                  f13_intensivebpmanagement$Visit=="F05"|
                                                                  f13_intensivebpmanagement$Visit=="F06"|
                                                                  f13_intensivebpmanagement$Visit=="F07"|
                                                                  f13_intensivebpmanagement$Visit=="F08"|
                                                                  f13_intensivebpmanagement$Visit=="F09"|
                                                                  f13_intensivebpmanagement$Visit=="F10"|
                                                                  f13_intensivebpmanagement$Visit=="F11"|
                                                                  f13_intensivebpmanagement$Visit=="F12"|
                                                                  f13_intensivebpmanagement$Visit=="F13"|
                                                                  f13_intensivebpmanagement$Visit=="F14"|
                                                                  f13_intensivebpmanagement$Visit=="F15"|
                                                                  f13_intensivebpmanagement$Visit=="F16"|
                                                                  f13_intensivebpmanagement$Visit=="F17"|
                                                                  f13_intensivebpmanagement$Visit=="F18"|
                                                                  f13_intensivebpmanagement$Visit=="F19"|
                                                                  f13_intensivebpmanagement$Visit=="F20"|
                                                                  f13_intensivebpmanagement$Visit=="F21"|
                                                                  f13_intensivebpmanagement$Visit=="F22"|
                                                                  f13_intensivebpmanagement$Visit=="F23"|
                                                                  f13_intensivebpmanagement$Visit=="F24"|
                                                                  f13_intensivebpmanagement$Visit=="F25"|
                                                                  f13_intensivebpmanagement$Visit=="F26"|
                                                                  f13_intensivebpmanagement$Visit=="F27"|
                                                                  f13_intensivebpmanagement$Visit=="F28"|
                                                                  f13_intensivebpmanagement$Visit=="F29"|
                                                                  f13_intensivebpmanagement$Visit=="F30"|
                                                                  f13_intensivebpmanagement$Visit=="F31"|
                                                                  f13_intensivebpmanagement$Visit=="F32"|
                                                                  f13_intensivebpmanagement$Visit=="F33"|
                                                                  f13_intensivebpmanagement$Visit=="F34"|
                                                                  f13_intensivebpmanagement$Visit=="F35"|
                                                                  f13_intensivebpmanagement$Visit=="F36"|
                                                                  f13_intensivebpmanagement$Visit=="F37"|
                                                                  f13_intensivebpmanagement$Visit=="F38"|
                                                                  f13_intensivebpmanagement$Visit=="F39"|
                                                                  f13_intensivebpmanagement$Visit=="F40"|
                                                                  f13_intensivebpmanagement$Visit=="F41"|
                                                                  f13_intensivebpmanagement$Visit=="F42"|
                                                                  f13_intensivebpmanagement$Visit=="F43"|
                                                                  f13_intensivebpmanagement$Visit=="F44"|
                                                                  f13_intensivebpmanagement$Visit=="F45"|
                                                                  f13_intensivebpmanagement$Visit=="F46"|
                                                                  f13_intensivebpmanagement$Visit=="F47"|
                                                                  f13_intensivebpmanagement$Visit=="F48"|
                                                                  f13_intensivebpmanagement$Visit=="F49"|
                                                                  f13_intensivebpmanagement$Visit=="F50"|
                                                                  f13_intensivebpmanagement$Visit=="F51"|
                                                                  f13_intensivebpmanagement$Visit=="F52"|
                                                                  f13_intensivebpmanagement$Visit=="F53"|
                                                                  f13_intensivebpmanagement$Visit=="F54"|
                                                                  f13_intensivebpmanagement$Visit=="F55"|
                                                                  f13_intensivebpmanagement$Visit=="F56"|
                                                                  f13_intensivebpmanagement$Visit=="F57"|
                                                                  f13_intensivebpmanagement$Visit=="F58"|
                                                                  f13_intensivebpmanagement$Visit=="F59"|
                                                                  f13_intensivebpmanagement$Visit=="F60"),]
f13_intensivebpmanagement_cut$visvec = as.numeric(f13_intensivebpmanagement_cut$Visit)
f13_intensivebpmanagement_cut$visnum = 0*(f13_intensivebpmanagement_cut$visvec==1)+
  1*(f13_intensivebpmanagement_cut$visvec==3)+
  2*(f13_intensivebpmanagement_cut$visvec==4)+
  3*(f13_intensivebpmanagement_cut$visvec==5)+
  4*(f13_intensivebpmanagement_cut$visvec==6)+
  6*(f13_intensivebpmanagement_cut$visvec==7)+
  8*(f13_intensivebpmanagement_cut$visvec==8)+
  12*(f13_intensivebpmanagement_cut$visvec==10)+
  16*(f13_intensivebpmanagement_cut$visvec==12)+
  20*(f13_intensivebpmanagement_cut$visvec==14)+
  24*(f13_intensivebpmanagement_cut$visvec==16)+
  28*(f13_intensivebpmanagement_cut$visvec==18)+
  30*(f13_intensivebpmanagement_cut$visvec==19)+
  32*(f13_intensivebpmanagement_cut$visvec==20)+
  34*(f13_intensivebpmanagement_cut$visvec==21)+
  36*(f13_intensivebpmanagement_cut$visvec==22)+
  38*(f13_intensivebpmanagement_cut$visvec==23)+
  40*(f13_intensivebpmanagement_cut$visvec==24)+
  42*(f13_intensivebpmanagement_cut$visvec==25)+
  44*(f13_intensivebpmanagement_cut$visvec==26)+
  46*(f13_intensivebpmanagement_cut$visvec==27)+
  48*(f13_intensivebpmanagement_cut$visvec==28)+
  50*(f13_intensivebpmanagement_cut$visvec==29)+
  52*(f13_intensivebpmanagement_cut$visvec==30)+
  54*(f13_intensivebpmanagement_cut$visvec==31)+
  56*(f13_intensivebpmanagement_cut$visvec==32)+
  58*(f13_intensivebpmanagement_cut$visvec==33)+
  60*(f13_intensivebpmanagement_cut$visvec==34)+
  10*(f13_intensivebpmanagement_cut$visvec==9)+
  14*(f13_intensivebpmanagement_cut$visvec==11)+
  18*(f13_intensivebpmanagement_cut$visvec==13)+
  22*(f13_intensivebpmanagement_cut$visvec==15)+
  26*(f13_intensivebpmanagement_cut$visvec==17)
f13_intensivebpmanagement_cut$visadv = f13_intensivebpmanagement_cut$visnum*f13_intensivebpmanagement_cut$advexp
f13_intensivebpmanagement_cut1 = summaryBy(advexp~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=max)
f13_intensivebpmanagement_cut1$advexp.max[f13_intensivebpmanagement_cut1$advexp.max=="-Inf"]=""
f13_intensivebpmanagement_cut2 = summaryBy(visadv~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=min)
f13_intensivebpmanagement_cut2$visadv.min[f13_intensivebpmanagement_cut2$visadv.min=="Inf"]=""
f13_intensivebpmanagement_cut3 = summaryBy(visnum~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=max)
f13_intensivebpmanagement_cut = merge(f13_intensivebpmanagement_cut1,f13_intensivebpmanagement_cut2,by="MaskID")
f13_intensivebpmanagement_cut = merge(f13_intensivebpmanagement_cut,f13_intensivebpmanagement_cut3,by="MaskID")
f14_standardbpmanagement_cut = f14_standardbpmanagement[which(f14_standardbpmanagement$Visit=="BLR"|
                                                                f14_standardbpmanagement$Visit=="F01"|
                                                                f14_standardbpmanagement$Visit=="F02"|
                                                                f14_standardbpmanagement$Visit=="F03"|
                                                                f14_standardbpmanagement$Visit=="F04"|
                                                                f14_standardbpmanagement$Visit=="F05"|
                                                                f14_standardbpmanagement$Visit=="F06"|
                                                                f14_standardbpmanagement$Visit=="F07"|
                                                                f14_standardbpmanagement$Visit=="F08"|
                                                                f14_standardbpmanagement$Visit=="F09"|
                                                                f14_standardbpmanagement$Visit=="F10"|
                                                                f14_standardbpmanagement$Visit=="F11"|
                                                                f14_standardbpmanagement$Visit=="F12"|
                                                                f14_standardbpmanagement$Visit=="F13"|
                                                                f14_standardbpmanagement$Visit=="F14"|
                                                                f14_standardbpmanagement$Visit=="F15"|
                                                                f14_standardbpmanagement$Visit=="F16"|
                                                                f14_standardbpmanagement$Visit=="F17"|
                                                                f14_standardbpmanagement$Visit=="F18"|
                                                                f14_standardbpmanagement$Visit=="F19"|
                                                                f14_standardbpmanagement$Visit=="F20"|
                                                                f14_standardbpmanagement$Visit=="F21"|
                                                                f14_standardbpmanagement$Visit=="F22"|
                                                                f14_standardbpmanagement$Visit=="F23"|
                                                                f14_standardbpmanagement$Visit=="F24"|
                                                                f14_standardbpmanagement$Visit=="F25"|
                                                                f14_standardbpmanagement$Visit=="F26"|
                                                                f14_standardbpmanagement$Visit=="F27"|
                                                                f14_standardbpmanagement$Visit=="F28"|
                                                                f14_standardbpmanagement$Visit=="F29"|
                                                                f14_standardbpmanagement$Visit=="F30"|
                                                                f14_standardbpmanagement$Visit=="F31"|
                                                                f14_standardbpmanagement$Visit=="F32"|
                                                                f14_standardbpmanagement$Visit=="F33"|
                                                                f14_standardbpmanagement$Visit=="F34"|
                                                                f14_standardbpmanagement$Visit=="F35"|
                                                                f14_standardbpmanagement$Visit=="F36"|
                                                                f14_standardbpmanagement$Visit=="F31"|
                                                                f14_standardbpmanagement$Visit=="F32"|
                                                                f14_standardbpmanagement$Visit=="F33"|
                                                                f14_standardbpmanagement$Visit=="F34"|
                                                                f14_standardbpmanagement$Visit=="F35"|
                                                                f14_standardbpmanagement$Visit=="F36"|
                                                                f14_standardbpmanagement$Visit=="F37"|
                                                                f14_standardbpmanagement$Visit=="F38"|
                                                                f14_standardbpmanagement$Visit=="F39"|
                                                                f14_standardbpmanagement$Visit=="F40"|
                                                                f14_standardbpmanagement$Visit=="F41"|
                                                                f14_standardbpmanagement$Visit=="F42"|
                                                                f14_standardbpmanagement$Visit=="F43"|
                                                                f14_standardbpmanagement$Visit=="F44"|
                                                                f14_standardbpmanagement$Visit=="F45"|
                                                                f14_standardbpmanagement$Visit=="F46"|
                                                                f14_standardbpmanagement$Visit=="F47"|
                                                                f14_standardbpmanagement$Visit=="F48"|
                                                                f14_standardbpmanagement$Visit=="F49"|
                                                                f14_standardbpmanagement$Visit=="F50"|
                                                                f14_standardbpmanagement$Visit=="F51"|
                                                                f14_standardbpmanagement$Visit=="F52"|
                                                                f14_standardbpmanagement$Visit=="F53"|
                                                                f14_standardbpmanagement$Visit=="F54"|
                                                                f14_standardbpmanagement$Visit=="F55"|
                                                                f14_standardbpmanagement$Visit=="F56"|
                                                                f14_standardbpmanagement$Visit=="F57"|
                                                                f14_standardbpmanagement$Visit=="F58"|
                                                                f14_standardbpmanagement$Visit=="F59"|
                                                                f14_standardbpmanagement$Visit=="F60"),]
f14_standardbpmanagement_cut$visvec = as.numeric(f14_standardbpmanagement_cut$Visit)
f14_standardbpmanagement_cut$visnum = 0*(f14_standardbpmanagement_cut$visvec==1)+
  1*(f14_standardbpmanagement_cut$visvec==3)+
  2*(f14_standardbpmanagement_cut$visvec==4)+
  3*(f14_standardbpmanagement_cut$visvec==5)+
  4*(f14_standardbpmanagement_cut$visvec==6)+
  6*(f14_standardbpmanagement_cut$visvec==7)+
  8*(f14_standardbpmanagement_cut$visvec==8)+
  12*(f14_standardbpmanagement_cut$visvec==9)+
  16*(f14_standardbpmanagement_cut$visvec==10)+
  20*(f14_standardbpmanagement_cut$visvec==11)+
  24*(f14_standardbpmanagement_cut$visvec==12)+
  28*(f14_standardbpmanagement_cut$visvec==13)+
  32*(f14_standardbpmanagement_cut$visvec==14)+
  36*(f14_standardbpmanagement_cut$visvec==15)+
  40*(f14_standardbpmanagement_cut$visvec==16)+
  44*(f14_standardbpmanagement_cut$visvec==17)+
  48*(f14_standardbpmanagement_cut$visvec==18)+
  52*(f14_standardbpmanagement_cut$visvec==19)+
  56*(f14_standardbpmanagement_cut$visvec==20)+
  60*(f14_standardbpmanagement_cut$visvec==21)
f14_standardbpmanagement_cut$visadv = f14_standardbpmanagement_cut$visnum*f14_standardbpmanagement_cut$advexp
f14_standardbpmanagement_cut1 = summaryBy(advexp~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=max)
f14_standardbpmanagement_cut1$advexp.max[f14_standardbpmanagement_cut1$advexp.max=="-Inf"]=""
f14_standardbpmanagement_cut2 = summaryBy(visadv~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=min)
f14_standardbpmanagement_cut2$visadv.min[f14_standardbpmanagement_cut2$visadv.min=="Inf"]=""
f14_standardbpmanagement_cut3 = summaryBy(visnum~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=max)
f14_standardbpmanagement_cut = merge(f14_standardbpmanagement_cut1,f14_standardbpmanagement_cut2,by="MaskID")
f14_standardbpmanagement_cut = merge(f14_standardbpmanagement_cut,f14_standardbpmanagement_cut3,by="MaskID")
bpmg = rbind(f13_intensivebpmanagement_cut,f14_standardbpmanagement_cut)
lipids_cut = lipids[which(lipids$Visit=="BLR"),]
otherlabs_cut = otherlabs[which(otherlabs$Visit=="BLR"),]
accord_set = merge(accord_key_cut,bloodpressure_cut,by="MaskID")
accord_set = merge(accord_set,concomitantmeds_cut,by=c("MaskID","Visit"))
accord_set = merge(accord_set,bpmg,by="MaskID")
accord_set = merge(accord_set,f07_baselinehistoryphysicalexam,by="MaskID")
accord_set = merge(accord_set,cvdoutcomes,by="MaskID")
accord_set = merge(accord_set,lipids_cut,by="MaskID")
accord_set = merge(accord_set,otherlabs_cut,by="MaskID")
save.image("~/Data/sprint_pop/data/accord_cut.RData")

load("~/Data/sprint_pop/data/accord_cut.RData")
cvd = (accord_set$censor_nmi==0)|(accord_set$censor_nst==0)|(accord_set$censor_cm==0)|(accord_set$censor_chf==0)|(accord_set$censor_maj==0)
t_censor = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25,accord_set$fuyrs_nst*365.25,accord_set$fuyrs_cm*365.25,accord_set$fuyrs_chf*365.25,accord_set$fuyrs_maj*365.25))
t_cvds = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25*(1-accord_set$censor_nmi),accord_set$fuyrs_nst*365.25*(1-accord_set$censor_nst),accord_set$fuyrs_cm*365.25*(1-accord_set$censor_cm),accord_set$fuyrs_chf*365.25*(1-accord_set$censor_chf),accord_set$fuyrs_maj*365.25*(1-accord_set$censor_maj)))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)
sae = (accord_set$advexp.max==1)
accord_set$visadv.min=as.numeric(accord_set$visadv.min)
t_censor = rowMaxs(cbind(accord_set$visnum.max*30.42))
t_saes = rowMaxs(cbind(accord_set$visadv.min*30.42))
t_saes[is.na(t_saes)] = t_censor[is.na(t_saes)]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
INTENSIVE = as.numeric(accord_set$INTENSIVE)
AGE = accord_set$baseline_age
FEMALE = accord_set$female
RACE_BLACK = as.numeric(accord_set$raceclass=="Black")
hisp = (accord_set$raceclass=="Hispanic")
SBP.y = accord_set$sbp
DBP.y = accord_set$dbp
N_AGENTS = (accord_set$loop+accord_set$thiazide+accord_set$ksparing+accord_set$a2rb+accord_set$acei+accord_set$dhp_ccb+accord_set$nondhp_ccb+accord_set$alpha_blocker+accord_set$central_agent+accord_set$beta_blocker+accord_set$vasodilator+accord_set$reserpine+accord_set$other_bpmed)
currentsmoker = (accord_set$cigarett==1)
formersmoker = (accord_set$smokelif==1)
ASPIRIN= accord_set$aspirin
STATIN = accord_set$statin
SUB_SENIOR = as.numeric(AGE>=75)
SUB_CKD = as.numeric(accord_set$gfr<60)
CHR = accord_set$chol
GLUR = accord_set$fpg
HDL = accord_set$hdl
TRR = accord_set$trig
UMALCR = accord_set$uacr
EGFR = accord_set$gfr
SCREAT = accord_set$screat
BMI = accord_set$wt_kg/((accord_set$ht_cm/1000)^2)/100
c2<-data.frame(cvd,t_cvds,sae,t_saes,
              INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
              SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
              ASPIRIN,STATIN,
              SCREAT,CHR,HDL,TRR,BMI,
              INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
              INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
              INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
              INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
c2=c2[complete.cases(c2),]
adm.cens=5*365.25
c2$fu.time <- pmin(c2$t_cvds, adm.cens)
c2$status <- ifelse(as.numeric(adm.cens < c2$t_cvds), 0, c2$cvd)

#### test CVD model on accord ####
survcox_c2<-coxph(data=c2, Surv(fu.time, status)~AGE+FEMALE+RACE_BLACK+hisp+SBP.y+DBP.y+
                    N_AGENTS+currentsmoker+formersmoker+ASPIRIN+STATIN+SCREAT+CHR+
                    HDL+TRR+BMI+INTENSIVE*AGE+INTENSIVE*RACE_BLACK+
                    INTENSIVE*DBP.y+INTENSIVE*currentsmoker+INTENSIVE*HDL+
                    INTENSIVE*TRR)
summary(survcox_c2)
survfit_c2=survfit(survcox_c2, newdata=c2, se.fit=FALSE)
estinc_c2=1-survfit_c2$surv[dim(survfit_c2$surv)[1],]
c2$dec=as.numeric(cut2(estinc_c2, g=10))
GND.result=GND.calib(pred=estinc_c2, tvar=c2$fu.time, out=c2$status, 
                     cens.t=adm.cens, groups=c2$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_c2,c2$cvd)


d2<-c2
d2=d2[complete.cases(d2),]
adm.cens=5*365.25
d2$fu.time <- pmin(d2$t_saes, adm.cens)
d2$status <- ifelse(as.numeric(adm.cens < d2$t_saes), 0, d2$sae)

#### test SAE model on accord  ####
survcox_d2<-coxph(data=d2, Surv(fu.time, status)~AGE+FEMALE+hisp+SBP.y+DBP.y+
                    N_AGENTS+currentsmoker+formersmoker+ASPIRIN+SCREAT+CHR+HDL+
                    INTENSIVE*FEMALE+INTENSIVE*currentsmoker+INTENSIVE*STATIN+INTENSIVE*SCREAT+INTENSIVE*CHR+INTENSIVE*TRR)
summary(survcox_d2)
survfit_d2=survfit(survcox_d2, newdata=d2, se.fit=FALSE)
estinc_d2=1-survfit_d2$surv[dim(survfit_d2$surv)[1],]
d2$dec=as.numeric(cut2(estinc_d2, g=10))
GND.result=GND.calib(pred=estinc_d2, tvar=d2$fu.time, out=d2$status, 
                     cens.t=adm.cens, groups=d2$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d2,d2$sae)


#### C/GDN ####
INTENSIVE = accord_set$INTENSIVE
INTENSIVE[INTENSIVE==0]=1
cint<-data.frame(cvd,t_cvds,sae,t_saes,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 SCREAT,CHR,HDL,TRR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
cint=cint[complete.cases(cint),]
adm.cens=5*365.25
cint$fu.time <- pmin(cint$t_cvds, adm.cens)
cint$status <- ifelse(as.numeric(adm.cens < cint$t_cvds), 0, cint$cvd)
survfit_cint=survfit(survcox_c, newdata=cint, se.fit=FALSE)
estinc_cint=1-survfit_cint$surv[dim(survfit_cint$surv)[1],]
INTENSIVE = accord_set$INTENSIVE
INTENSIVE[INTENSIVE==1]=0
cstd<-data.frame(cvd,t_cvds,sae,t_saes,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 SCREAT,CHR,HDL,TRR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
cstd=cstd[complete.cases(cstd),]
adm.cens=5*365.25
cstd$fu.time <- pmin(cstd$t_cvds, adm.cens)
cstd$status <- ifelse(as.numeric(adm.cens < cstd$t_cvds), 0, cstd$cvd)
survfit_cstd=survfit(survcox_c, newdata=cstd, se.fit=FALSE)
estinc_cstd=1-survfit_cstd$surv[dim(survfit_cstd$surv)[1],]
benefit = estinc_cstd-estinc_cint
hist(benefit,xlab="Predicted reduction in CVD events/CVD death \n from intensive treatment (probability)")
survfit_dint=survfit(survcox_d, newdata=cint, se.fit=FALSE)
estinc_dint=1-survfit_dint$surv[dim(survfit_dint$surv)[1],]
survfit_dstd=survfit(survcox_d, newdata=cstd, se.fit=FALSE)
estinc_dstd=1-survfit_dstd$surv[dim(survfit_dstd$surv)[1],]
risk = estinc_dint-estinc_dstd
hist(risk,xlab="Predicted increase in serious adverse events \n from intensive treatment (probability)")
plot(sprintbenefit,sprintrisk, xlim=c(-.1,.5),ylim=c(-.1,.5),xlab="Predicted absolute risk reduction in CVD events/death", ylab="Predicted absolute risk increase in serious adverse events", col=rgb(0,0,100,50,maxColorValue=255), pch=19)
points(benefit,risk, xlim=c(-.1,.5),ylim=c(-.1,.5),xlab="Predicted absolute risk reduction in CVD events/death", ylab="Predicted absolute risk increase in serious adverse events", col=rgb(255,140,0,50,maxColorValue=255), pch=19)
abline(a=0,b=0, col = "gray60")
abline(v=0, col = "gray60")

INTENSIVE = accord_set$INTENSIVE
pv = c2
pv$fu.time <- pmin(c2$t_cvds, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < c2$t_cvds), 0, c2$cvd)

betax=(survcox_c2$coefficients[1]*pv$AGE+
         survcox_c2$coefficients[2]*pv$FEMALE+
         survcox_c2$coefficients[3]*pv$RACE_BLACK+
         survcox_c2$coefficients[4]*pv$hisp+
         survcox_c2$coefficients[5]*pv$SBP.y+
         survcox_c2$coefficients[6]*pv$DBP.y+
         survcox_c2$coefficients[7]*pv$N_AGENTS+
         survcox_c2$coefficients[8]*pv$currentsmoker+
         survcox_c2$coefficients[9]*pv$formersmoker+
         survcox_c2$coefficients[10]*pv$ASPIRIN+
         survcox_c2$coefficients[11]*pv$STATIN+
         survcox_c2$coefficients[12]*pv$SCREAT+
         survcox_c2$coefficients[13]*pv$CHR+
         survcox_c2$coefficients[14]*pv$HDL+
         survcox_c2$coefficients[15]*pv$TRR+
         survcox_c2$coefficients[16]*pv$BMI+
         survcox_c2$coefficients[17]*pv$INTENSIVE+
         survcox_c2$coefficients[18]*pv$INTENSIVE*pv$AGE+
         survcox_c2$coefficients[19]*pv$INTENSIVE*pv$RACE_BLACK+
         survcox_c2$coefficients[20]*pv$INTENSIVE*pv$DBP.y+
         survcox_c2$coefficients[21]*pv$INTENSIVE*pv$currentsmoker+
         survcox_c2$coefficients[22]*pv$INTENSIVE*pv$HDL+
         survcox_c2$coefficients[23]*pv$INTENSIVE*pv$TRR)
cvdpred = 1 - .881^exp(betax-mean(na.omit(betax)))
estinc_e=cvdpred
pv$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=pv$fu.time, out=pv$status, 
                     cens.t=adm.cens, groups=pv$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,pv$cvd)

benpredbetaint = (survcox_c2$coefficients[1]*pv$AGE+
                    survcox_c2$coefficients[2]*pv$FEMALE+
                    survcox_c2$coefficients[3]*pv$RACE_BLACK+
                    survcox_c2$coefficients[4]*pv$hisp+
                    survcox_c2$coefficients[5]*pv$SBP.y+
                    survcox_c2$coefficients[6]*pv$DBP.y+
                    survcox_c2$coefficients[7]*pv$N_AGENTS+
                    survcox_c2$coefficients[8]*pv$currentsmoker+
                    survcox_c2$coefficients[9]*pv$formersmoker+
                    survcox_c2$coefficients[10]*pv$ASPIRIN+
                    survcox_c2$coefficients[11]*pv$STATIN+
                    survcox_c2$coefficients[12]*pv$SCREAT+
                    survcox_c2$coefficients[13]*pv$CHR+
                    survcox_c2$coefficients[14]*pv$HDL+
                    survcox_c2$coefficients[15]*pv$TRR+
                    survcox_c2$coefficients[16]*pv$BMI+
                    survcox_c2$coefficients[17]*1+
                    survcox_c2$coefficients[18]*1*pv$AGE+
                    survcox_c2$coefficients[19]*1*pv$RACE_BLACK+
                    survcox_c2$coefficients[20]*1*pv$DBP.y+
                    survcox_c2$coefficients[21]*1*pv$currentsmoker+
                    survcox_c2$coefficients[22]*1*pv$HDL+
                    survcox_c2$coefficients[23]*1*pv$TRR)
benintpred = 1 - .881^exp(benpredbetaint-mean(na.omit(betax)))

benpredbetastd = (survcox_c2$coefficients[1]*pv$AGE+
                    survcox_c2$coefficients[2]*pv$FEMALE+
                    survcox_c2$coefficients[3]*pv$RACE_BLACK+
                    survcox_c2$coefficients[4]*pv$hisp+
                    survcox_c2$coefficients[5]*pv$SBP.y+
                    survcox_c2$coefficients[6]*pv$DBP.y+
                    survcox_c2$coefficients[7]*pv$N_AGENTS+
                    survcox_c2$coefficients[8]*pv$currentsmoker+
                    survcox_c2$coefficients[9]*pv$formersmoker+
                    survcox_c2$coefficients[10]*pv$ASPIRIN+
                    survcox_c2$coefficients[11]*pv$STATIN+
                    survcox_c2$coefficients[12]*pv$SCREAT+
                    survcox_c2$coefficients[13]*pv$CHR+
                    survcox_c2$coefficients[14]*pv$HDL+
                    survcox_c2$coefficients[15]*pv$TRR+
                    survcox_c2$coefficients[16]*pv$BMI+
                    survcox_c2$coefficients[17]*0+
                    survcox_c2$coefficients[18]*0*pv$AGE+
                    survcox_c2$coefficients[19]*0*pv$RACE_BLACK+
                    survcox_c2$coefficients[20]*0*pv$DBP.y+
                    survcox_c2$coefficients[21]*0*pv$currentsmoker+
                    survcox_c2$coefficients[22]*0*pv$HDL+
                    survcox_c2$coefficients[23]*0*pv$TRR)
benstdpred = 1 - .881^exp(benpredbetastd-mean(na.omit(betax)))

netben = benstdpred-benintpred 
hist(netben)


pv$fu.time <- pmin(c2$t_saes, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < c2$t_saes), 0, c2$sae)
betax=(survcox_d$coefficients[1]*pv$AGE+
         survcox_d$coefficients[2]*pv$FEMALE+
         survcox_d$coefficients[3]*pv$hisp+
         survcox_d$coefficients[4]*pv$SBP.y+
         survcox_d$coefficients[5]*pv$DBP.y+
         survcox_d$coefficients[6]*pv$N_AGENTS+
         survcox_d$coefficients[7]*pv$currentsmoker+
         survcox_d$coefficients[8]*pv$formersmoker+
         survcox_d$coefficients[9]*pv$ASPIRIN+
         survcox_d$coefficients[10]*pv$SCREAT+
         survcox_d$coefficients[11]*pv$CHR+
         survcox_d$coefficients[12]*pv$HDL+
         survcox_d$coefficients[13]*pv$INTENSIVE+
         survcox_d$coefficients[14]*pv$STATIN+
         survcox_d$coefficients[15]*pv$TRR+
         survcox_d$coefficients[16]*pv$FEMALE*pv$INTENSIVE+
         survcox_d$coefficients[17]*pv$INTENSIVE*pv$currentsmoker+
         survcox_d$coefficients[18]*pv$STATIN*pv$INTENSIVE+
         survcox_d$coefficients[19]*pv$INTENSIVE*pv$SCREAT+
         survcox_d$coefficients[20]*pv$INTENSIVE*pv$CHR+
         survcox_d$coefficients[21]*pv$INTENSIVE*pv$TRR)
saepred = 1 - .887^exp(betax-mean(na.omit(betax)))
estinc_e=saepred
pv$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=pv$fu.time, out=pv$status, 
                     cens.t=adm.cens, groups=pv$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,pv$sae)

harmpredbetaint = (survcox_d$coefficients[1]*pv$AGE+
                     survcox_d$coefficients[2]*pv$FEMALE+
                     survcox_d$coefficients[3]*pv$hisp+
                     survcox_d$coefficients[4]*pv$SBP.y+
                     survcox_d$coefficients[5]*pv$DBP.y+
                     survcox_d$coefficients[6]*pv$N_AGENTS+
                     survcox_d$coefficients[7]*pv$currentsmoker+
                     survcox_d$coefficients[8]*pv$formersmoker+
                     survcox_d$coefficients[9]*pv$ASPIRIN+
                     survcox_d$coefficients[10]*pv$SCREAT+
                     survcox_d$coefficients[11]*pv$CHR+
                     survcox_d$coefficients[12]*pv$HDL+
                     survcox_d$coefficients[13]*1+
                     survcox_d$coefficients[14]*pv$STATIN+
                     survcox_d$coefficients[15]*pv$TRR+
                     survcox_d$coefficients[16]*pv$FEMALE*1+
                     survcox_d$coefficients[17]*1*pv$currentsmoker+
                     survcox_d$coefficients[18]*pv$STATIN*1+
                     survcox_d$coefficients[19]*1*pv$SCREAT+
                     survcox_d$coefficients[20]*1*pv$CHR+
                     survcox_d$coefficients[21]*1*pv$TRR)
harmintpred = 1 - .887^exp(harmpredbetaint-mean(na.omit(betax)))

harmpredbetastd = (survcox_d$coefficients[1]*pv$AGE+
                     survcox_d$coefficients[2]*pv$FEMALE+
                     survcox_d$coefficients[3]*pv$hisp+
                     survcox_d$coefficients[4]*pv$SBP.y+
                     survcox_d$coefficients[5]*pv$DBP.y+
                     survcox_d$coefficients[6]*pv$N_AGENTS+
                     survcox_d$coefficients[7]*pv$currentsmoker+
                     survcox_d$coefficients[8]*pv$formersmoker+
                     survcox_d$coefficients[9]*pv$ASPIRIN+
                     survcox_d$coefficients[10]*pv$SCREAT+
                     survcox_d$coefficients[11]*pv$CHR+
                     survcox_d$coefficients[12]*pv$HDL+
                     survcox_d$coefficients[13]*0+
                     survcox_d$coefficients[14]*pv$STATIN+
                     survcox_d$coefficients[15]*pv$TRR+
                     survcox_d$coefficients[16]*pv$FEMALE*0+
                     survcox_d$coefficients[17]*0*pv$currentsmoker+
                     survcox_d$coefficients[18]*pv$STATIN*0+
                     survcox_d$coefficients[19]*0*pv$SCREAT+
                     survcox_d$coefficients[20]*0*pv$CHR+
                     survcox_d$coefficients[21]*0*pv$TRR)
harmstdpred = 1 - .887^exp(harmpredbetastd-mean(na.omit(betax)))

netharm = harmintpred-harmstdpred 
hist(netharm)

plot(netben,netharm, xlim=c(-.2,.4),ylim=c(-.2,.4),xlab="Predicted reduction in CVD events/death", ylab="Predicted increase in serious adverse events")
abline(a=0,b=0, col = "gray60")
abline(v=0, col = "gray60")

netbenaccord = netben
netharmaccord = netharm

#### scatterplots ARD####
plot(-netbenaccord,netharmaccord, 
     xlim=c(min(-netbenaccord,-netbensprint),max(-netbenaccord,-netbensprint)), 
     ylim=c(min(netharmaccord,netharmsprint),max(netharmaccord,netharmsprint)),
     col=rgb(255,140,0,70,maxColorValue=255), pch=16,
     xlab = c("Change in predicted probability of CVD events/death w/ int Rx"),
     ylab = c("\n Change in predicted probability of serious adverse events w/ int Rx"))
points(-netbensprint,netharmsprint, col=rgb(30,144,255,50,maxColorValue=255), pch=16)
abline(h=0,v=0,col="gray60")
legend(-.4, .5, c("SPRINT", "ACCORD-BP"), col = c(rgb(30,144,255,255,maxColorValue=255),rgb(255,140,0,255,maxColorValue=255)),
       text.col = "black", bty = "n", pch = c(16,16),
       bg = "white")


#### score vs ACCORD outcomes ####
bencats = c(.01,.03)
bencat = 1*(netben<bencats[1])+
  2*((netben>=bencats[1])&(netben<bencats[2]))+
  3*((netben>=bencats[2]))
cvdtable = describeBy(pv$cvd,list(bencat,pv$INTENSIVE),mat=TRUE)
cvdtable
prop.test(x=c(cvdtable[1,5]*cvdtable[1,6],cvdtable[4,5]*cvdtable[4,6]), n=c(cvdtable[1,5],cvdtable[4,5]), correct=FALSE)
prop.test(x=c(cvdtable[2,5]*cvdtable[2,6],cvdtable[5,5]*cvdtable[5,6]), n=c(cvdtable[2,5],cvdtable[5,5]), correct=FALSE)
prop.test(x=c(cvdtable[3,5]*cvdtable[3,6],cvdtable[6,5]*cvdtable[6,6]), n=c(cvdtable[3,5],cvdtable[6,5]), correct=FALSE)
cvdtabletimes = describeBy(pv$t_cvds,list(bencat,pv$INTENSIVE),mat=TRUE)
x1i = cvdtable[1:3,5]*cvdtable[1:3,6]
x2i = cvdtable[4:6,5]*cvdtable[4:6,6]
t1i = cvdtabletimes[1:3,5]*cvdtabletimes[1:3,6]
t2i = cvdtabletimes[4:6,5]*cvdtabletimes[4:6,6]
rma(measure="IR", x1i, x2i, t1i, t2i, method="REML")

harmcats = c(.005,.04)
harmcat = 1*(netharm<harmcats[1])+
  2*((netharm>=harmcats[1])&(netharm<harmcats[2]))+
  3*((netharm>=harmcats[2]))
saetable = describeBy(pv$sae,list(harmcat,pv$INTENSIVE),mat=TRUE)
saetable
prop.test(x=c(saetable[1,5]*saetable[1,6],saetable[4,5]*saetable[4,6]), n=c(saetable[1,5],saetable[4,5]), correct=FALSE)
prop.test(x=c(saetable[2,5]*saetable[2,6],saetable[5,5]*saetable[5,6]), n=c(saetable[2,5],saetable[5,5]), correct=FALSE)
prop.test(x=c(saetable[3,5]*saetable[3,6],saetable[6,5]*saetable[6,6]), n=c(saetable[3,5],saetable[6,5]), correct=FALSE)
saetabletimes = describeBy(pv$t_saes,list(harmcat,pv$INTENSIVE),mat=TRUE)
x1i = saetable[1:3,5]*saetable[1:3,6]
x2i = saetable[4:6,5]*saetable[4:6,6]
t1i = saetabletimes[1:3,5]*saetabletimes[1:3,6]
t2i = saetabletimes[4:6,5]*saetabletimes[4:6,6]
rma(measure="IR", x1i, x2i, t1i, t2i, method="REML")



