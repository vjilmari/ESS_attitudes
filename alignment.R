#testing alignment with sirt package

library(sirt)

dat<-read.csv2("C:/Users/vjilmari/R/ESS_attitudes/dat.no.miss.csv",stringsAsFactors = F)

imm.vars<-c("gvrfgap.R","rfgfrpc","rfgbfml.R")
env.vars<-c("inctxff.R","sbsrnen.R","banhhap.R")

#standardization function

na.standardize<-function(var){
  (var-mean(var,na.rm=T))/sd(var,na.rm=T)
}

dat$gvrfgap.R.z<-na.standardize(dat$gvrfgap.R)
dat$rfgfrpc.z<-na.standardize(dat$rfgfrpc)
dat$rfgbfml.R.z<-na.standardize(dat$rfgbfml.R)

imm.vars.z<-paste0(imm.vars,".z")
imm.vars.z

describeBy(dat.imm[,imm.vars],group = dat.imm$cntry)

table(dat.imm[,imm.vars[1]],dat.imm$cntry)
table(dat.imm[,imm.vars[2]],dat.imm$cntry)
table(dat.imm[,imm.vars[3]],dat.imm$cntry)

dat.imm<-dat %>% 
  filter(cntry!="HU") %>% 
  dplyr::select(cntry,gvrfgap.R,rfgfrpc,rfgbfml.R) %>%
  na.omit()

par <- invariance_alignment_cfa_config(dat = dat.imm[,imm.vars],
                                       group = dat.imm$cntry)
par

#Czech republic, lithuania, and poland have negative variances on the first item


mod1 <- invariance.alignment(lambda = par$lambda,
                             nu =par$nu,
                             align.scale = c(0.2, 0.4),
                             align.pow = c(0.25, 0.25))

mod1$es.invariance["R2",]

#The loadings is impossible, it should be between 0 and 1. This is likely because of weird results for certain
#countries

cmod1 <- invariance_alignment_constraints(mod1,
                                          lambda_parm_tol = 0.4,
                                          nu_parm_tol = 0.2)
summary(cmod1)

#the same problems are found with different item order

#look if excluding Lithuania would solve some problems

dat.imm<-dat %>% 
  filter(cntry!="HU" & cntry!="LT") %>% 
  dplyr::select(cntry,gvrfgap.R,rfgfrpc,rfgbfml.R) %>%
  na.omit()

par <- invariance_alignment_cfa_config(dat = dat.imm[,imm.vars],
                                       group = dat.imm$cntry)
par

config.lambdas<-data.frame(par$lambda)

#I think these are the standardized factor loadings from configural model
#values beyond 1 are Heywood cases, but it may not matter if we are able to align them somehow
round(config.lambdas,2)


mod1 <- invariance.alignment(lambda = par$lambda,
                             nu =par$nu,
                             align.scale = c(0.2, 0.4),
                             align.pow = c(0.25, 0.25))

mod1$es.invariance["R2",]
mod1

#Now most of the invariance can be absorbed by group varying estimates

aligned.lambdas<-data.frame(mod1$lambda.aligned)
round(aligned.lambdas,2)

cbind(round(config.lambdas,2),
      round(aligned.lambdas,2))

cmod1 <- invariance_alignment_constraints(mod1,
                                          lambda_parm_tol = 0.4,
                                          nu_parm_tol = 0.2)
summary(cmod1)

#I guess that these tell about the item parameters that would be obtained with the given alignment


## see what the weighted sum scores would look like

config.lambdas$cntry<-row.names(config.lambdas)
aligned.lambdas$cntry<-row.names(aligned.lambdas)

## weighing based on configural
dat<-left_join(x=dat,
               y=config.lambdas,
               by="cntry",
               suffix=c("",".conf.weight"))

dat<-left_join(x=dat,
               y=aligned.lambdas,
               by="cntry",
               suffix=c("",".alig.weight"))

na.standardize<-function(var){
  (var-mean(var,na.rm=T))/sd(var,na.rm=T)
}

dat$gvrfgap.R.z<-na.standardize(dat$gvrfgap.R)
dat$rfgfrpc.z<-na.standardize(dat$rfgfrpc)
dat$rfgbfml.R.z<-na.standardize(dat$rfgbfml.R)

z.imm.vars<-paste0(imm.vars,".z")


dat$imm.unweighted.mean<-rowMeans(dat[,imm.vars])
dat$imm.unweighted.sum<-rowSums(dat[,imm.vars])
dat$imm.conf.mean<-rowMeans(data.frame(dat$gvrfgap.R*dat$gvrfgap.R.conf.weight,
                              dat$rfgfrpc*dat$rfgfrpc.conf.weight,
                              dat$rfgbfml.R*dat$rfgbfml.R.conf.weight))
dat$imm.conf.sum<-rowSums(data.frame(dat$gvrfgap.R*dat$gvrfgap.R.conf.weight,
                                         dat$rfgfrpc*dat$rfgfrpc.conf.weight,
                                         dat$rfgbfml.R*dat$rfgbfml.R.conf.weight))

dat$imm.alig.mean<-rowMeans(data.frame(dat$gvrfgap.R*dat$gvrfgap.R.alig.weight,
                                       dat$rfgfrpc*dat$rfgfrpc.alig.weight,
                                       dat$rfgbfml.R*dat$rfgbfml.R.alig.weight))
dat$imm.alig.sum<-rowSums(data.frame(dat$gvrfgap.R*dat$gvrfgap.R.alig.weight,
                                     dat$rfgfrpc*dat$rfgfrpc.alig.weight,
                                     dat$rfgbfml.R*dat$rfgbfml.R.alig.weight))

dat$imm.unweighted.mean.z<-rowMeans(dat[,z.imm.vars])
dat$imm.unweighted.sum.z<-rowSums(dat[,z.imm.vars])
dat$imm.conf.mean.z<-rowMeans(data.frame(dat$gvrfgap.R.z*dat$gvrfgap.R.conf.weight,
                                       dat$rfgfrpc.z*dat$rfgfrpc.conf.weight,
                                       dat$rfgbfml.R.z*dat$rfgbfml.R.conf.weight))
dat$imm.conf.sum.z<-rowSums(data.frame(dat$gvrfgap.R.z*dat$gvrfgap.R.conf.weight,
                                     dat$rfgfrpc.z*dat$rfgfrpc.conf.weight,
                                     dat$rfgbfml.R.z*dat$rfgbfml.R.conf.weight))

dat$imm.alig.mean.z<-rowMeans(data.frame(dat$gvrfgap.R.z*dat$gvrfgap.R.alig.weight,
                                       dat$rfgfrpc.z*dat$rfgfrpc.alig.weight,
                                       dat$rfgbfml.R.z*dat$rfgbfml.R.alig.weight))
dat$imm.alig.sum.z<-rowSums(data.frame(dat$gvrfgap.R.z*dat$gvrfgap.R.alig.weight,
                                     dat$rfgfrpc.z*dat$rfgfrpc.alig.weight,
                                     dat$rfgbfml.R.z*dat$rfgbfml.R.alig.weight))

imm.sum.vars<-c("imm.unweighted.mean",
                "imm.unweighted.sum",
                "imm.conf.mean",
                "imm.conf.sum",
                "imm.alig.mean",
                "imm.alig.sum")

imm.sum.vars<-c(imm.sum.vars,paste0(imm.sum.vars,".z"))

describe(dat[,imm.sum.vars])

round(cor(dat[,imm.sum.vars],
      use="pairwise.complete.obs"),2)
