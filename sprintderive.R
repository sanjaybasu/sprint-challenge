library(gdata)
library(glmnet) 
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
attach(sprint_set)
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
sae = (HYP_SAE_EVNT==1)|(SYN_SAE_EVNT==1)|(ELE_SAE_EVNT==1)|(AKI_SAE_EVNT==1)|(BRA_SAE_EVNT==1)
t_censor = rowMaxs(cbind(HYP_SAE_DAYS,SYN_SAE_DAYS,ELE_SAE_DAYS,AKI_SAE_DAYS,BRA_SAE_DAYS)) 
t_saes = rowMaxs(cbind(HYP_SAE_DAYS*HYP_SAE_EVNT,SYN_SAE_DAYS*SYN_SAE_EVNT,ELE_SAE_DAYS*ELE_SAE_EVNT,AKI_SAE_DAYS*AKI_SAE_EVNT,BRA_SAE_DAYS*BRA_SAE_EVNT)) 
t_saes[t_saes==0] = t_censor[t_saes==0]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
testsubset = data.frame(cOutcome,
                        AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        SCREAT,CHR,HDL,TRR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]
fit.lasso <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=1)
fit.ridge <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=0)
fit.elnet <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=.5)
fit.lasso.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=1, 
                          family="cox")
fit.ridge.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=0,
                          family="cox")
fit.elnet.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=.5,
                          family="cox")
for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), 
                                            alpha=i/10,family="cox"))
}
par(mfrow=c(3,2))
plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")
plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")
plot(fit.elnet, xvar="lambda")
plot(fit5, main="Elastic Net")

testsubset = data.frame(dOutcome,
                        AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        SCREAT,CHR,HDL,TRR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]
fit.lasso <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=1)
fit.ridge <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=0)
fit.elnet <- glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), family="cox", alpha=.5)
fit.lasso.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=1, 
                          family="cox")
fit.ridge.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=0,
                          family="cox")
fit.elnet.cv <- cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), alpha=.5,
                          family="cox")
for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(as.matrix(testsubset[,-c(1)]), as.matrix(testsubset[,1]), 
                                            alpha=i/10,family="cox"))
}
par(mfrow=c(3,2))
plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")
plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")
plot(fit.elnet, xvar="lambda")
plot(fit5, main="Elastic Net")
