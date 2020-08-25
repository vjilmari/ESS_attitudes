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
      






       
