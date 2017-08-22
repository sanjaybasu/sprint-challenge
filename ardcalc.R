
pv = data.frame(cvd,t_cvds,sae,t_saes,
                INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                ASPIRIN,STATIN,
                SCREAT,CHR,HDL,TRR,BMI,
                INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
pv=pv[complete.cases(pv),]
adm.cens=5*365.25
pv$fu.time <- pmin(pv$t_cvds, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < pv$t_cvds), 0, pv$cvd)

betax=(survcox_c[1]*pv$AGE+
                 survcox_c[2]*pv$FEMALE+
                 survcox_c[3]*pv$RACE_BLACK+
                 survcox_c[4]*pv$hisp+
                 survcox_c[5]*pv$SBP.y+
                 survcox_c[6]*pv$DBP.y+
                 survcox_c[7]*pv$N_AGENTS+
                 survcox_c[8]*pv$currentsmoker+
                 survcox_c[9]*pv$formersmoker+
                 survcox_c[10]*pv$ASPIRIN+
                 survcox_c[11]*pv$STATIN+
                 survcox_c[12]*pv$SCREAT+
                 survcox_c[13]*pv$CHR+
                 survcox_c[14]*pv$HDL+
                 survcox_c[15]*pv$TRR+
                 survcox_c[16]*pv$BMI+
                 survcox_c[17]*pv$INTENSIVE+
                 survcox_c[18]*pv$INTENSIVE*pv$AGE+
                 survcox_c[19]*pv$INTENSIVE*pv$RACE_BLACK+
                 survcox_c[20]*pv$INTENSIVE*pv$DBP.y+
                 survcox_c[21]*pv$INTENSIVE*pv$currentsmoker+
                 survcox_c[22]*pv$INTENSIVE*pv$HDL+
                 survcox_c[23]*pv$INTENSIVE*pv$TRR)

benpredbetaint = (survcox_c[1]*pv$AGE+
                    survcox_c[2]*pv$FEMALE+
                    survcox_c[3]*pv$RACE_BLACK+
                    survcox_c[4]*pv$hisp+
                    survcox_c[5]*pv$SBP.y+
                    survcox_c[6]*pv$DBP.y+
                    survcox_c[7]*pv$N_AGENTS+
                    survcox_c[8]*pv$currentsmoker+
                    survcox_c[9]*pv$formersmoker+
                    survcox_c[10]*pv$ASPIRIN+
                    survcox_c[11]*pv$STATIN+
                    survcox_c[12]*pv$SCREAT+
                    survcox_c[13]*pv$CHR+
                    survcox_c[14]*pv$HDL+
                    survcox_c[15]*pv$TRR+
                    survcox_c[16]*pv$BMI+
                    survcox_c[17]*1+
                    survcox_c[18]*1*pv$AGE+
                    survcox_c[19]*1*pv$RACE_BLACK+
                    survcox_c[20]*1*pv$DBP.y+
                    survcox_c[21]*1*pv$currentsmoker+
                    survcox_c[22]*1*pv$HDL+
                    survcox_c[23]*1*pv$TRR)

benintpred = 1 - .943^exp(benpredbetaint-mean(na.omit(betax)))

benpredbetastd = (survcox_c[1]*pv$AGE+
                    survcox_c[2]*pv$FEMALE+
                    survcox_c[3]*pv$RACE_BLACK+
                    survcox_c[4]*pv$hisp+
                    survcox_c[5]*pv$SBP.y+
                    survcox_c[6]*pv$DBP.y+
                    survcox_c[7]*pv$N_AGENTS+
                    survcox_c[8]*pv$currentsmoker+
                    survcox_c[9]*pv$formersmoker+
                    survcox_c[10]*pv$ASPIRIN+
                    survcox_c[11]*pv$STATIN+
                    survcox_c[12]*pv$SCREAT+
                    survcox_c[13]*pv$CHR+
                    survcox_c[14]*pv$HDL+
                    survcox_c[15]*pv$TRR+
                    survcox_c[16]*pv$BMI+
                    survcox_c[17]*0+
                    survcox_c[18]*0*pv$AGE+
                    survcox_c[19]*0*pv$RACE_BLACK+
                    survcox_c[20]*0*pv$DBP.y+
                    survcox_c[21]*0*pv$currentsmoker+
                    survcox_c[22]*0*pv$HDL+
                    survcox_c[23]*0*pv$TRR)
benstdpred = 1 - .943^exp(benpredbetastd-mean(na.omit(betax)))

netben = benstdpred-benintpred 

pv$fu.time <- pmin(pv$t_saes, adm.cens)
pv$status <- ifelse(as.numeric(adm.cens < pv$t_saes), 0, pv$sae)

betax=(survcox_d[1]*pv$AGE+
         survcox_d[2]*pv$FEMALE+
         survcox_d[3]*pv$hisp+
         survcox_d[4]*pv$SBP.y+
         survcox_d[5]*pv$DBP.y+
         survcox_d[6]*pv$N_AGENTS+
         survcox_d[7]*pv$currentsmoker+
         survcox_d[8]*pv$formersmoker+
         survcox_d[9]*pv$ASPIRIN+
         survcox_d[10]*pv$SCREAT+
         survcox_d[11]*pv$CHR+
         survcox_d[12]*pv$HDL+
         survcox_d[13]*pv$INTENSIVE+
         survcox_d[14]*pv$STATIN+
         survcox_d[15]*pv$TRR+
         survcox_d[16]*pv$FEMALE*pv$INTENSIVE+
         survcox_d[17]*pv$INTENSIVE*pv$currentsmoker+
         survcox_d[18]*pv$STATIN*pv$INTENSIVE+
         survcox_d[19]*pv$INTENSIVE*pv$SCREAT+
         survcox_d[20]*pv$INTENSIVE*pv$CHR+
         survcox_d[21]*pv$INTENSIVE*pv$TRR)

harmpredbetaint = (survcox_d[1]*pv$AGE+
                     survcox_d[2]*pv$FEMALE+
                     survcox_d[3]*pv$hisp+
                     survcox_d[4]*pv$SBP.y+
                     survcox_d[5]*pv$DBP.y+
                     survcox_d[6]*pv$N_AGENTS+
                     survcox_d[7]*pv$currentsmoker+
                     survcox_d[8]*pv$formersmoker+
                     survcox_d[9]*pv$ASPIRIN+
                     survcox_d[10]*pv$SCREAT+
                     survcox_d[11]*pv$CHR+
                     survcox_d[12]*pv$HDL+
                     survcox_d[13]*1+
                     survcox_d[14]*pv$STATIN+
                     survcox_d[15]*pv$TRR+
                     survcox_d[16]*pv$FEMALE*1+
                     survcox_d[17]*1*pv$currentsmoker+
                     survcox_d[18]*pv$STATIN*1+
                     survcox_d[19]*1*pv$SCREAT+
                     survcox_d[20]*1*pv$CHR+
                     survcox_d[21]*1*pv$TRR)
harmintpred = 1 - .897^exp(harmpredbetaint-mean(na.omit(betax)))

harmpredbetastd = (survcox_d[1]*pv$AGE+
                     survcox_d[2]*pv$FEMALE+
                     survcox_d[3]*pv$hisp+
                     survcox_d[4]*pv$SBP.y+
                     survcox_d[5]*pv$DBP.y+
                     survcox_d[6]*pv$N_AGENTS+
                     survcox_d[7]*pv$currentsmoker+
                     survcox_d[8]*pv$formersmoker+
                     survcox_d[9]*pv$ASPIRIN+
                     survcox_d[10]*pv$SCREAT+
                     survcox_d[11]*pv$CHR+
                     survcox_d[12]*pv$HDL+
                     survcox_d[13]*0+
                     survcox_d[14]*pv$STATIN+
                     survcox_d[15]*pv$TRR+
                     survcox_d[16]*pv$FEMALE*0+
                     survcox_d[17]*0*pv$currentsmoker+
                     survcox_d[18]*pv$STATIN*0+
                     survcox_d[19]*0*pv$SCREAT+
                     survcox_d[20]*0*pv$CHR+
                     survcox_d[21]*0*pv$TRR)
harmstdpred = 1 - .897^exp(harmpredbetastd-mean(na.omit(betax)))

netharm = harmintpred-harmstdpred 

bencats = c(.01,.03)
bencat = 1*(netben<bencats[1])+
  2*((netben>=bencats[1])&(netben<bencats[2]))+
  3*((netben>=bencats[2]))
cvdtable = describeBy(pv$cvd,list(bencat,pv$INTENSIVE),mat=TRUE)
cvdtable
prop.test(x=c(cvdtable[1,5]*cvdtable[1,6],cvdtable[4,5]*cvdtable[4,6]), n=c(cvdtable[1,5],cvdtable[4,5]), correct=FALSE)
prop.test(x=c(cvdtable[2,5]*cvdtable[2,6],cvdtable[5,5]*cvdtable[5,6]), n=c(cvdtable[2,5],cvdtable[5,5]), correct=FALSE)
prop.test(x=c(cvdtable[3,5]*cvdtable[3,6],cvdtable[6,5]*cvdtable[6,6]), n=c(cvdtable[3,5],cvdtable[6,5]), correct=FALSE)

harmcats = c(.005,.04)
harmcat = 1*(netharm<harmcats[1])+
  2*((netharm>=harmcats[1])&(netharm<harmcats[2]))+
  3*((netharm>=harmcats[2]))
saetable = describeBy(pv$sae,list(harmcat,pv$INTENSIVE),mat=TRUE)
saetable
prop.test(x=c(saetable[1,5]*saetable[1,6],saetable[4,5]*saetable[4,6]), n=c(saetable[1,5],saetable[4,5]), correct=FALSE)
prop.test(x=c(saetable[2,5]*saetable[2,6],saetable[5,5]*saetable[5,6]), n=c(saetable[2,5],saetable[5,5]), correct=FALSE)
prop.test(x=c(saetable[3,5]*saetable[3,6],saetable[6,5]*saetable[6,6]), n=c(saetable[3,5],saetable[6,5]), correct=FALSE)
