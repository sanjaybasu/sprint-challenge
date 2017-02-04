# SPRINT Challenge entry, Sanjay Basu et al., basus@stanford.edu

# INSTRUCTIONS: 
# (1) install all packages listed under the "load packages" header
# (2) adjust command "registerDoMC(cores=4)" under "load packages" heading to specify the number of processor cores available on your machine, if you wish to implement parallel processing for speed (which is recommended)
# (3) change calls to appropriate directory on your local machine, to load/save data from your desired location. this code does not contain the data itself, which requires IRB and NIH-BioLINCC approval, but is built to analyze the NIH-BioLINCC versions of the SPRINT-Pop and ACCORD-BP datasets 
# (4) note the Framingham risk score is recalculated here for SPRINT participants, rather than using the score listed in the SPRINT-Pop data, as the latter was found to be incorrectly calculated


#### load packages ####

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
registerDoMC(cores=4)

#### merge sprintpop data ####
rm(list=ls())
setwd("~/Data/sprint_pop/data")
baseline = read.sas7bdat("baseline.sas7bdat")
bp = read.sas7bdat("bp.sas7bdat")
outcomes = read.sas7bdat("outcomes.sas7bdat")
retention = read.sas7bdat("retention.sas7bdat")
safety = read.sas7bdat("safety.sas7bdat")
save.image("~/Data/sprint_pop/data/sprint.RData")

rm(list=ls())
load("~/Data/sprint_pop/data/sprint.RData")
bp_cut = bp[which(bp$VISITCODE=="RZ"),]
sprint_set = merge(baseline,bp_cut,by="MASKID")
sprint_set = merge(sprint_set,outcomes,by="MASKID")
sprint_set = merge(sprint_set,retention,by="MASKID")
sprint_set = merge(sprint_set,safety,by="MASKID")
save.image("~/Data/sprint_pop/data/sprint_cut.RData")

##### load derivation data ####
rm(list=ls())
load("~/Data/sprint_pop/data/sprint_cut.RData")
keep(sprint_set, sure=TRUE)
attach(sprint_set)
hisp = (RACE4=="HISPANIC")
currentsmoker = (SMOKE_3CAT==3)
formersmoker = (SMOKE_3CAT==2)
frammen = 3.06117*log(AGE)+1.12370*log(CHR)-0.93263*log(HDL)+1.93303*log(SBP.y)*(N_AGENTS==0)+1.99881*log(SBP.y)*(N_AGENTS>0)+0.65451*currentsmoker+0.57367*1
framfem = 2.32888*log(AGE)+1.20904*log(CHR)-0.70833*log(HDL)+2.76157*log(SBP.y)*(N_AGENTS==0) + 2.82263*log(SBP.y)*(N_AGENTS>0)+0.52873*currentsmoker+0.69154*1
RISK10YRS =100*((1-FEMALE)*(1-0.88936^exp(frammen - 23.9802))+FEMALE*(1-0.95012^exp(framfem - 26.1931)))


#### gen MI/stroke/CVD death model ####
cvd = (EVENT_MI==1)|(EVENT_STROKE==1)|(EVENT_CVDDEATH==1)
t_censor = rowMaxs(cbind(T_MI,T_STROKE,T_CVDDEATH))
t_cvds = rowMaxs(cbind(T_MI*EVENT_MI,T_STROKE*EVENT_STROKE,T_CVDDEATH*EVENT_CVDDEATH))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)

testsubset = data.frame(cOutcome,
                        INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]
cvdmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox", parallel = TRUE)
plot(cvdmodel.cv)
cvdmodel.cv
coef.min = coef(cvdmodel.cv, s = 'lambda.min')
coef.min
c<-data.frame(cvd,t_cvds,
              INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
              SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
              ASPIRIN,STATIN,
              EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
              INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
              INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
              INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
              INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
c=c[complete.cases(c),]
adm.cens=3*365.25
c$fu.time <- pmin(c$t_cvds, adm.cens)
c$status <- ifelse(as.numeric(adm.cens < c$t_cvds), 0, c$cvd)
mean(c$status)
survcox_c<-coxph(data=c, Surv(fu.time, status)~AGE+FEMALE+RACE_BLACK+hisp+SBP.y+
                   N_AGENTS+currentsmoker+formersmoker+ASPIRIN+STATIN+EGFR+SCREAT+
                   CHR+GLUR+HDL+TRR+UMALCR+BMI+
                   INTENSIVE*FEMALE+INTENSIVE*RACE_BLACK+INTENSIVE*DBP.y+
                   INTENSIVE*currentsmoker+INTENSIVE*ASPIRIN+INTENSIVE*STATIN+
                   INTENSIVE*EGFR+INTENSIVE*CHR+INTENSIVE*HDL+INTENSIVE*TRR+
                   INTENSIVE*UMALCR)
summary(survcox_c)
survfit_c=survfit(survcox_c, newdata=c, se.fit=FALSE)
estinc_c=1-survfit_c$surv[dim(survfit_c$surv)[1],]
c$dec=as.numeric(cut2(estinc_c, g=9))
table(c$dec, c$status)
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
}#kmdec

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
  plot(tapply(datause$pred,datause$dec,mean),1-kmtab[,1],xlab="Expected K-M rate",ylab="Observed K-M rate",xlim=c(0,.15),ylim=c(0,.15),main="CVD events or CVD death")
  abline(a=0,b=1, col = "gray60")
  
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
}#GND.calib


GND.result=GND.calib(pred=estinc_c, tvar=c$fu.time, out=c$status, 
                     cens.t=adm.cens, groups=c$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_c,c$cvd)


#### gen adverse events outcome model ####
sae = (HYP_SAE_EVNT==1)|(SYN_SAE_EVNT==1)|(ELE_SAE_EVNT==1)|(AKI_SAE_EVNT==1) #|(BRA_SAE_EVNT==1)
t_censor = rowMaxs(cbind(HYP_SAE_DAYS,SYN_SAE_DAYS,ELE_SAE_DAYS,AKI_SAE_DAYS)) #,BRA_SAE_DAYS
t_saes = rowMaxs(cbind(HYP_SAE_DAYS*HYP_SAE_EVNT,SYN_SAE_DAYS*SYN_SAE_EVNT,ELE_SAE_DAYS*ELE_SAE_EVNT,AKI_SAE_DAYS*AKI_SAE_EVNT)) #,BRA_SAE_DAYS*BRA_SAE_EVNT
t_saes[t_saes==0] = t_censor[t_saes==0]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
hisp = (RACE4=="HISPANIC")
currentsmoker = (SMOKE_3CAT==3)
formersmoker = (SMOKE_3CAT==2)
testsubset = data.frame(dOutcome,
                        INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]
saemodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox", parallel = TRUE)
plot(saemodel.cv)
saemodel.cv
coef.min = coef(saemodel.cv, s = 'lambda.min')
coef.min
d<-data.frame(sae,t_saes,
              INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
              SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
              ASPIRIN,STATIN,
              EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
              INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
              INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
              INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
              INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
d=d[complete.cases(d),]
adm.cens=3*365.25
d$fu.time <- pmin(d$t_saes, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_saes), 0, d$sae)
mean(d$status)
chartable =describeBy(d,d$INTENSIVE,mat=TRUE)
chartable
chartable[is.na(chartable)]=0

survcox_d<-coxph(data=d, Surv(fu.time, status)~AGE+FEMALE+hisp+SBP.y+DBP.y+
                   N_AGENTS+currentsmoker+formersmoker+EGFR+SCREAT+CHR+HDL+UMALCR+BMI+
                   INTENSIVE*FEMALE+INTENSIVE*STATIN+INTENSIVE*SCREAT+INTENSIVE*TRR+INTENSIVE*UMALCR)
summary(survcox_d)
survfit_d=survfit(survcox_d, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
table(d$dec, d$status)
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
}#kmdec

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
  plot(tapply(datause$pred,datause$dec,mean),1-kmtab[,1],xlab="Expected K-M rate",ylab="Observed K-M rate",xlim=c(0,.3),ylim=c(0,.3), main="Serious adverse events")
  abline(a=0,b=1, col = "gray60")
  
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
}#GND.calib


GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$sae)


#### clin pred score development ####
# CVD model predict benefit as CVD risk under standard minus under intensive Rx
INTENSIVE = sprint_set$INTENSIVE
INTENSIVE[INTENSIVE==0]=1
cint<-data.frame(cvd,t_cvds,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
cint=cint[complete.cases(cint),]
adm.cens=3*365.25
cint$fu.time <- pmin(cint$t_cvds, adm.cens)
cint$status <- ifelse(as.numeric(adm.cens < cint$t_cvds), 0, cint$cvd)
mean(cint$status)
survfit_cint=survfit(survcox_c, newdata=cint, se.fit=FALSE)
estinc_cint=1-survfit_cint$surv[dim(survfit_cint$surv)[1],]

INTENSIVE = sprint_set$INTENSIVE
INTENSIVE[INTENSIVE==1]=0
cstd<-data.frame(cvd,t_cvds,
                 INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                 SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                 ASPIRIN,STATIN,
                 EGFR,SCREAT,CHR,GLUR,HDL,TRR,UMALCR,BMI,
                 INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                 INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                 INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                 INTENSIVE*EGFR,INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*GLUR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*UMALCR,INTENSIVE*BMI)
cstd=cstd[complete.cases(cstd),]
adm.cens=3*365.25
cstd$fu.time <- pmin(cstd$t_cvds, adm.cens)
cstd$status <- ifelse(as.numeric(adm.cens < cstd$t_cvds), 0, cstd$cvd)
mean(cstd$status)
survfit_cstd=survfit(survcox_c, newdata=cstd, se.fit=FALSE)
estinc_cstd=1-survfit_cstd$surv[dim(survfit_cstd$surv)[1],]

benefit = estinc_cstd-estinc_cint
summary(benefit)
hist(benefit,xlab="Predicted reduction in CVD events/CVD death \n from intensive treatment (probability)")

# SAE model predict risk as SAE risk under intensive minus under standard Rx
survfit_dint=survfit(survcox_d, newdata=cint, se.fit=FALSE)
estinc_dint=1-survfit_dint$surv[dim(survfit_dint$surv)[1],]

survfit_dstd=survfit(survcox_d, newdata=cstd, se.fit=FALSE)
estinc_dstd=1-survfit_dstd$surv[dim(survfit_dstd$surv)[1],]

risk = estinc_dint-estinc_dstd
summary(risk)
hist(risk,xlab="Predicted increase in serious adverse events \n from intensive treatment (probability)")

benriskdiff = benefit - risk
summary(benriskdiff)
hist(benriskdiff)
plot(benefit,risk, xlim=c(-.2,.4),ylim=c(-.2,.4),xlab="Predicted reduction in CVD events/death", ylab="Predicted increase in serious adverse events")
abline(a=0,b=1, col = "gray60")
abline(a=0,b=0, col = "gray60")
abline(v=0, col = "gray60")

benriskmodel = lm(benriskdiff~AGE + FEMALE + RACE_BLACK + 
                    hisp + SBP.y + N_AGENTS + currentsmoker + formersmoker + 
                    ASPIRIN + STATIN + EGFR + SCREAT + CHR + GLUR + HDL + TRR + 
                    UMALCR + BMI + INTENSIVE * FEMALE + INTENSIVE * RACE_BLACK + 
                    INTENSIVE * DBP.y + INTENSIVE * currentsmoker + INTENSIVE * 
                    ASPIRIN + INTENSIVE * STATIN + INTENSIVE * EGFR + INTENSIVE * 
                    CHR + INTENSIVE * HDL + INTENSIVE * TRR + INTENSIVE * UMALCR+ 
                    AGE + FEMALE + hisp + 
                    SBP.y + DBP.y + N_AGENTS + currentsmoker + formersmoker + 
                    EGFR + SCREAT + CHR + HDL + UMALCR + BMI + INTENSIVE * FEMALE + 
                    INTENSIVE * STATIN + INTENSIVE * SCREAT + INTENSIVE * TRR + 
                    INTENSIVE * UMALCR, data=cint)
summary(benriskmodel)

benriskmodel = lm(benriskdiff~AGE+FEMALE+RACE_BLACK+hisp+
                    SBP.y+N_AGENTS+currentsmoker+formersmoker+
                    ASPIRIN+STATIN+EGFR+SCREAT+CHR+TRR+UMALCR+DBP.y, data=cint)
summary(benriskmodel)
benriskmodelfinal = benriskmodel

#### score vs SPRINT outcomes ####
benriskmodel_pred = predict(benriskmodel,sprint_set)
summary(benriskmodel_pred)
hist(benriskmodel_pred)
quantile(benriskmodel_pred,na.rm=TRUE)
sprintrisk = benriskmodel_pred
riskcatquants = c(quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.25)),quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.5)),quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.75)))
sprintcats = riskcatquants
riskcat = 1*(benriskmodel_pred<riskcatquants[1])+
  2*((benriskmodel_pred>=riskcatquants[1])&(benriskmodel_pred<riskcatquants[2]))+
  3*((benriskmodel_pred>=riskcatquants[2])&(benriskmodel_pred<riskcatquants[3]))+
  4*(benriskmodel_pred>=riskcatquants[3])
intensive = sprint_set$INTENSIVE

benriskquants = c(quantile(benriskdiff,na.rm=TRUE,probs=c(.25)),quantile(benriskdiff,na.rm=TRUE,probs=c(.5)),quantile(benriskdiff,na.rm=TRUE,probs=c(.75)))
benriskcat = 1*(benriskdiff<benriskquants[1])+
  2*((benriskdiff>=benriskquants[1])&(benriskdiff<benriskquants[2]))+
  3*((benriskdiff>=benriskquants[2])&(benriskdiff<benriskquants[3]))+
  4*(benriskdiff>=benriskquants[3])
benriskmodel_pred2 = predict(benriskmodel,data=cint)
riskcat2 = 1*(benriskmodel_pred2<riskcatquants[1])+
  2*((benriskmodel_pred2>=riskcatquants[1])&(benriskmodel_pred2<riskcatquants[2]))+
  3*((benriskmodel_pred2>=riskcatquants[2])&(benriskmodel_pred2<riskcatquants[3]))+
  4*(benriskmodel_pred2>=riskcatquants[3])
table(riskcat2,benriskcat)

cvdtable = describeBy(cvd,list(riskcat,intensive),mat=TRUE)
cvdtable
prop.test(x=c(cvdtable[1,5]*cvdtable[1,6],cvdtable[5,5]*cvdtable[5,6]), n=c(cvdtable[1,5],cvdtable[5,5]), correct=FALSE)
prop.test(x=c(cvdtable[2,5]*cvdtable[2,6],cvdtable[6,5]*cvdtable[6,6]), n=c(cvdtable[2,5],cvdtable[6,5]), correct=FALSE)
prop.test(x=c(cvdtable[3,5]*cvdtable[3,6],cvdtable[7,5]*cvdtable[7,6]), n=c(cvdtable[3,5],cvdtable[7,5]), correct=FALSE)
prop.test(x=c(cvdtable[4,5]*cvdtable[4,6],cvdtable[8,5]*cvdtable[8,6]), n=c(cvdtable[4,5],cvdtable[8,5]), correct=FALSE)
survdiff(Surv(t_cvds, cvd)~ riskcat+strata(sprint_set$INTENSIVE))


saetable = describeBy(sae,list(riskcat,intensive),mat=TRUE)
saetable
prop.test(x=c(saetable[1,5]*saetable[1,6],saetable[5,5]*saetable[5,6]), n=c(saetable[1,5],saetable[5,5]), correct=FALSE)
prop.test(x=c(saetable[2,5]*saetable[2,6],saetable[6,5]*saetable[6,6]), n=c(saetable[2,5],saetable[6,5]), correct=FALSE)
prop.test(x=c(saetable[3,5]*saetable[3,6],saetable[7,5]*saetable[7,6]), n=c(saetable[3,5],saetable[7,5]), correct=FALSE)
prop.test(x=c(saetable[4,5]*saetable[4,6],saetable[8,5]*saetable[8,6]), n=c(saetable[4,5],saetable[8,5]), correct=FALSE)
survdiff(Surv(t_saes, sae)~ riskcat+strata(sprint_set$INTENSIVE))


fit0 <- survfit(Surv(t_cvds[riskcat==1], cvd[riskcat==1]) ~ sprint_set$INTENSIVE[riskcat==1])
ggsurvplot(
  fit0,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 0, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,5*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


fit1 <- survfit(Surv(t_cvds[riskcat==2], cvd[riskcat==2]) ~ sprint_set$INTENSIVE[riskcat==2])
ggsurvplot(
  fit1,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 1, without diabetes",
  legend = c(0.2,.8),
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)



fit2 <- survfit(Surv(t_cvds[riskcat==3], cvd[riskcat==3]) ~ sprint_set$INTENSIVE[riskcat==3])
ggsurvplot(
  fit2,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 2, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)



fit3 <- survfit(Surv(t_cvds[riskcat==4], cvd[riskcat==4]) ~ sprint_set$INTENSIVE[riskcat==4])
ggsurvplot(
  fit3,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 3, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


sfit0 <- survfit(Surv(t_saes[riskcat==1], sae[riskcat==1]) ~ sprint_set$INTENSIVE[riskcat==1])
ggsurvplot(
  sfit0,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 0, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit1 <- survfit(Surv(t_saes[riskcat==2], sae[riskcat==2]) ~ sprint_set$INTENSIVE[riskcat==2])
ggsurvplot(
  sfit1,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 1, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit2 <- survfit(Surv(t_saes[riskcat==3], sae[riskcat==3]) ~ sprint_set$INTENSIVE[riskcat==3])
ggsurvplot(
  sfit2,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 2, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit3 <- survfit(Surv(t_saes[riskcat==4], sae[riskcat==4]) ~ sprint_set$INTENSIVE[riskcat==4])
ggsurvplot(
  sfit3,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend = c(0.2,.8),
  legend.labs=c("Standard","Intensive"),
  main = "Score = 3, without diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


#### merge accord data ####
setwd("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets")
accord_key = read.sas7bdat("accord_key.sas7bdat")
bloodpressure = read.sas7bdat("bloodpressure.sas7bdat")
concomitantmeds = read.sas7bdat("concomitantmeds.sas7bdat")
cvdoutcomes = read.sas7bdat("cvdoutcomes.sas7bdat")
lipids = read.sas7bdat("lipids.sas7bdat")
otherlabs = read.sas7bdat("otherlabs.sas7bdat")
setwd("~/Data/accord/4-Data_Sets-CRFs/4a-CRF_Data_Sets")
f13_intensivebpmanagement = read.sas7bdat("f13_intensivebpmanagement.sas7bdat")
f14_standardbpmanagement = read.sas7bdat("f14_standardbpmanagement.sas7bdat")
f07_baselinehistoryphysicalexam = read.sas7bdat("f07_baselinehistoryphysicalexam.sas7bdat")
save.image("~/Data/sprint_pop/data/accord.RData")

load("~/Data/sprint_pop/data/accord.RData")
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
                                                                  f13_intensivebpmanagement$Visit=="F36"),]
f13_intensivebpmanagement_cut$visvec = as.numeric(f13_intensivebpmanagement_cut$Visit)
f13_intensivebpmanagement_cut$visnum = 0*(f13_intensivebpmanagement_cut$visvec==1)+1*(f13_intensivebpmanagement_cut$visvec==3)+2*(f13_intensivebpmanagement_cut$visvec==4)+3*(f13_intensivebpmanagement_cut$visvec==5)+4*(f13_intensivebpmanagement_cut$visvec==6)+6*(f13_intensivebpmanagement_cut$visvec==7)+8*(f13_intensivebpmanagement_cut$visvec==8)+12*(f13_intensivebpmanagement_cut$visvec==10)+16*(f13_intensivebpmanagement_cut$visvec==12)+20*(f13_intensivebpmanagement_cut$visvec==14)+24*(f13_intensivebpmanagement_cut$visvec==16)+28*(f13_intensivebpmanagement_cut$visvec==18)+30*(f13_intensivebpmanagement_cut$visvec==19)+32*(f13_intensivebpmanagement_cut$visvec==20)+34*(f13_intensivebpmanagement_cut$visvec==21)+36*(f13_intensivebpmanagement_cut$visvec==22)+10*(f13_intensivebpmanagement_cut$visvec==9)+14*(f13_intensivebpmanagement_cut$visvec==11)+18*(f13_intensivebpmanagement_cut$visvec==13)+22*(f13_intensivebpmanagement_cut$visvec==15)+26*(f13_intensivebpmanagement_cut$visvec==17)
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
                                                                f14_standardbpmanagement$Visit=="F36"),]
f14_standardbpmanagement_cut$visvec = as.numeric(f14_standardbpmanagement_cut$Visit)
f14_standardbpmanagement_cut$visnum = 0*(f14_standardbpmanagement_cut$visvec==1)+1*(f14_standardbpmanagement_cut$visvec==3)+2*(f14_standardbpmanagement_cut$visvec==4)+3*(f14_standardbpmanagement_cut$visvec==5)+4*(f14_standardbpmanagement_cut$visvec==6)+6*(f14_standardbpmanagement_cut$visvec==7)+8*(f14_standardbpmanagement_cut$visvec==8)+12*(f14_standardbpmanagement_cut$visvec==9)+16*(f14_standardbpmanagement_cut$visvec==10)+20*(f14_standardbpmanagement_cut$visvec==11)+24*(f14_standardbpmanagement_cut$visvec==12)+28*(f14_standardbpmanagement_cut$visvec==13)+32*(f14_standardbpmanagement_cut$visvec==14)+36*(f14_standardbpmanagement_cut$visvec==15)
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


#### score vs ACCORD outcomes ####
# CVD model predict benefit as CVD risk under standard minus under intensive Rx
load("~/Data/sprint_pop/data/accord_cut.RData")
cvd = (accord_set$censor_nmi==0)|(accord_set$censor_nst==0)|(accord_set$censor_cm==0)
t_censor = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25,accord_set$fuyrs_nst*365.25,accord_set$fuyrs_cm*365.25))
t_cvds = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25*(1-accord_set$censor_nmi),accord_set$fuyrs_nst*365.25*(1-accord_set$censor_nst),accord_set$fuyrs_cm*365.25*(1-accord_set$censor_cm)))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)
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
UMALCR = accord_set$ualb
EGFR = accord_set$gfr
SCREAT = accord_set$screat
BMI = accord_set$wt_kg/((accord_set$ht_cm/1000)^2)/100
frammen = 3.06117*log(AGE)+1.12370*log(CHR)-0.93263*log(HDL)+1.93303*log(SBP.y)*(N_AGENTS==0)+1.99881*log(SBP.y)*(N_AGENTS>0)+0.65451*currentsmoker+0.57367*1
framfem = 2.32888*log(AGE)+1.20904*log(CHR)-0.70833*log(HDL)+2.76157*log(SBP.y)*(N_AGENTS==0) + 2.82263*log(SBP.y)*(N_AGENTS>0)+0.52873*currentsmoker+0.69154*1
RISK10YRS =100*((1-FEMALE)*(1-0.88936^exp(frammen - 23.9802))+FEMALE*(1-0.95012^exp(framfem - 26.1931)))
sae = (accord_set$advexp.max==1)
accord_set$visadv.min=as.numeric(accord_set$visadv.min)
t_censor = rowMaxs(cbind(accord_set$visnum.max*30.42))
t_saes = rowMaxs(cbind(accord_set$visadv.min*30.42))
t_saes[is.na(t_saes)] = t_censor[is.na(t_saes)]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)

benriskmodel_pred = predict(benriskmodelfinal,accord_set)
summary(benriskmodel_pred)
accordrisk = benriskmodel_pred
hist(benriskmodel_pred)
quantile(benriskmodel_pred,na.rm=TRUE)
riskcatquants = c(quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.25)),quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.5)),quantile(benriskmodel_pred,na.rm=TRUE,probs=c(.75)))
accordcats = riskcatquants
riskcat = 1*(benriskmodel_pred<riskcatquants[1])+
  2*((benriskmodel_pred>=riskcatquants[1])&(benriskmodel_pred<riskcatquants[2]))+
  3*((benriskmodel_pred>=riskcatquants[2])&(benriskmodel_pred<riskcatquants[3]))+
  4*(benriskmodel_pred>=riskcatquants[3])

intensive = accord_set$INTENSIVE
cvdtable = describeBy(cvd,list(riskcat,intensive),mat=TRUE)
saetable = describeBy(sae,list(riskcat,intensive),mat=TRUE)

fit0 <- survfit(Surv(t_cvds[riskcat==1], cvd[riskcat==1]) ~ accord_set$INTENSIVE[riskcat==1])
ggsurvplot(
  fit0,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 0, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


fit1 <- survfit(Surv(t_cvds[riskcat==2], cvd[riskcat==2]) ~ accord_set$INTENSIVE[riskcat==2])
ggsurvplot(
  fit1,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 1, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)



fit2 <- survfit(Surv(t_cvds[riskcat==3], cvd[riskcat==3]) ~ accord_set$INTENSIVE[riskcat==3])
ggsurvplot(
  fit2,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 2, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)



fit3 <- survfit(Surv(t_cvds[riskcat==4], cvd[riskcat==4]) ~ accord_set$INTENSIVE[riskcat==4])
ggsurvplot(
  fit3,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 3, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of CVD events/deaths",
  risk.table = TRUE,       # show risk table.
  pval = T,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


sfit0 <- survfit(Surv(t_saes[riskcat==1], sae[riskcat==1]) ~ accord_set$INTENSIVE[riskcat==1])
ggsurvplot(
  sfit0,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 0, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit1 <- survfit(Surv(t_saes[riskcat==2], sae[riskcat==2]) ~ accord_set$INTENSIVE[riskcat==2])
ggsurvplot(
  sfit1,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 1, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit2 <- survfit(Surv(t_saes[riskcat==3], sae[riskcat==3]) ~ accord_set$INTENSIVE[riskcat==3])
ggsurvplot(
  sfit2,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 2, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)

sfit3 <- survfit(Surv(t_saes[riskcat==4], sae[riskcat==4]) ~ accord_set$INTENSIVE[riskcat==4])
ggsurvplot(
  sfit3,                     # survfit object with calculated statistics.
  fun = "cumhaz",
  legend.labs=c("Standard","Intensive"),
  main = "Score = 3, with diabetes",
  font.legend=c(14),font.x=c(14),font.y=c(14),font.tickslab=c(14),font.main = c(14),
  xlab = "Days",
  ylab = "Cumulative hazard \n of severe adverse events",
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  ylim = c(0,.2),
  xlim = c(0,3*365.25),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 365,     # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T # show bars instead of names in text annotations
  # in legend of risk table
)


save.image("~/Data/sprint_pop/data/challenge.RData")



