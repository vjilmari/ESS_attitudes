cor.dat<-read.csv2("associations_within_voting_groups.csv",
                   stringsAsFactors = F)

library(metafor)
library(numform)

cor.dat$observed_r.z<-transf.rtoz(cor.dat$observed_r)
cor.dat$r.SE<-1/sqrt(cor.dat$n.y-3)
cor.dat$LL<-transf.ztor(cor.dat$observed_r.z+qnorm(.025)*cor.dat$r.SE)
cor.dat$UL<-transf.ztor(cor.dat$observed_r.z+qnorm(.975)*cor.dat$r.SE)
cor.dat$z<-cor.dat$observed_r.z/cor.dat$r.SE
cor.dat$p<-(1-pnorm(abs(cor.dat$z)))*2




anti.imm.r<-cbind.data.frame(group=cor.dat[cor.dat$anti.imm==1,c("voting.group")],
                            round(cor.dat[cor.dat$anti.imm==1,c("observed_r","LL","UL")],2),
                            p=round(cor.dat[cor.dat$anti.imm==1,c("p")],3))

anti.imm.r$r<-f_num(x=anti.imm.r[,"observed_r"],digits=2)


anti.imm.r$CI<-paste0("[",
                     f_num(x=anti.imm.r[,"LL"],digits=2),", ",
                     f_num(x=anti.imm.r[,"UL"],digits=2),"]")

anti.imm.r$pno0<-f_num(x=anti.imm.r[,"p"],digits=3)
anti.imm.r

write.csv2(anti.imm.r,"anti.imm.r.csv")

pro.env.r<-cbind.data.frame(group=cor.dat[cor.dat$pro.env==1,c("voting.group")],
                             round(cor.dat[cor.dat$pro.env==1,c("observed_r","LL","UL")],2),
                             p=round(cor.dat[cor.dat$pro.env==1,c("p")],3))

pro.env.r$r<-f_num(x=pro.env.r[,"observed_r"],digits=2)


pro.env.r$CI<-paste0("[",
                     f_num(x=pro.env.r[,"LL"],digits=2),", ",
                     f_num(x=pro.env.r[,"UL"],digits=2),"]")

pro.env.r$pno0<-f_num(x=pro.env.r[,"p"],digits=3)
pro.env.r

write.csv2(pro.env.r,"pro.env.r.csv")
      




# For table 1


tab1.mod2<-read.csv2("FE.H1.mod1.csv",
                   stringsAsFactors = F)       

tab1.mod2<-data.frame(Eff=tab1.mod2$Eff,
                      Est=tab1.mod2$Est,
                      CI=paste0("[",tab1.mod2[,"LL"],", ",tab1.mod2[,"UL"],"]"),
                      p=tab1.mod2$p)

tab1.mod3<-read.csv2("FE.H1.mod2.csv",
                     stringsAsFactors = F)

tab1.mod3<-data.frame(Eff=tab1.mod3$Eff,
                      Est=tab1.mod3$Est,
                      CI=paste0("[",tab1.mod3[,"LL"],", ",tab1.mod3[,"UL"],"]"),
                      p=tab1.mod3$p)

tab1.mod4<-read.csv2("FE.H1.mod3.csv",
                     stringsAsFactors = F)

tab1.mod4<-data.frame(Eff=tab1.mod4$Eff,
                      Est=tab1.mod4$Est,
                      CI=paste0("[",tab1.mod4[,"LL"],", ",tab1.mod4[,"UL"],"]"),
                      p=tab1.mod4$p)


tab1<-full_join(x=tab1.mod2,
          y=tab1.mod3,
          by="Eff")
tab1<-full_join(x=tab1,
                y=tab1.mod4,
                by="Eff")
tab1

write.csv2(tab1,"tab1.csv")

confint.merMod(H1.mod1,parm="theta_",method="profile",level=.95,oldNames=F)
confint.merMod(H1.mod2,parm="theta_",method="profile",level=.95,oldNames=F)
#don't run the random slopes -model again, CIs below
#confint.merMod(H1.mod3,parm="theta_",method="profile",level=.95,oldNames=F)

#d_(Intercept)|voting.group                0.2793160 0.34313941
#cor_environ.lvl1.(Intercept)|voting.group -0.1898851 0.44957436
#sd_environ.lvl1|voting.group               0.0242823 0.06094131
#sd_(Intercept)|cntry                       0.3720908 0.71144406
#cor_environ.lvl1.(Intercept)|cntry        -0.5288256 0.45582915
#sd_environ.lvl1|cntry                      0.0242841 0.06260749
#sigma                                      1.0215694 1.03677305
