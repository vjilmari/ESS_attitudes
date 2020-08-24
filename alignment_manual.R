#Try the manual alignment approach found in https://cran.r-project.org/web/packages/sirt/sirt.pdf page 153

library(semTools)
library(sirt)

dat<-read.csv2("C:/Users/vjilmari/R/ESS_attitudes/dat.no.miss.csv",stringsAsFactors = F)

imm.vars<-c("gvrfgap.R","rfgfrpc","rfgbfml.R")
env.vars<-c("inctxff.R","sbsrnen.R","banhhap.R")

dat.imm<-dat %>% 
  filter(cntry!="HU" & cntry!="LT") %>% 
  dplyr::select(cntry,gvrfgap.R,rfgfrpc,rfgbfml.R) %>%
  na.omit()

describeBy(dat.imm[,imm.vars],group = dat.imm$cntry)

table(dat.imm[,imm.vars[1]],dat.imm$cntry)
table(dat.imm[,imm.vars[2]],dat.imm$cntry)
table(dat.imm[,imm.vars[3]],dat.imm$cntry)




model1 <- "
F.imm=~ gvrfgap.R + rfgfrpc + rfgbfml.R
F.imm ~~ 1*F.imm
"


res <- cfa(model=model1, std.lv=TRUE, data=dat.imm, group=c("cntry"))
summary(res,standardized=T,fit=T)


# extract item parameters separate group analyses
ipars <- lavaan::parameterEstimates(res)
ipars
# extract lambda's: groups are in rows, items in columns
lambda <- matrix( ipars[ ipars$op=="=~", "est"], nrow=18, byrow=TRUE)
colnames(lambda) <- colnames(dat.imm)[-1]

# extract nu's
nu <- matrix( ipars[ ipars$op=="~1" & ipars$se !=0, "est" ], nrow=18, byrow=TRUE)
colnames(nu) <- colnames(dat.imm)[-1]
nu


# Model 1: least squares optimization (I guess this means alignment scale values at 1, and power at .50)
mod1 <- invariance.alignment( lambda=lambda, nu=nu )
summary(mod1)

mod2 <- invariance.alignment( lambda=lambda, nu=nu, align.pow=c(.5,.5) )
summary(mod2)

mod3 <- invariance.alignment( lambda=lambda, nu=nu, align.scale = c(0.2, 0.4),
                              align.pow = c(0.25, 0.25))

cmod3 <- invariance_alignment_constraints(mod3,
                                              lambda_parm_tol = 0.2,
                                              nu_parm_tol = 0.4)
summary(cmod3)
str(cmod3)


mod.constraints<-data.frame(cmod3$lambda_list$parm_est)
rownames(mod.constraints)<-unique(dat.imm$cntry)
mod.constraints

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



mod.free<-mod.constraints
mod.free[,1]<-ifelse(mod.free[,1]==getmode(mod.constraints[,1]),0,1)
mod.free[,2]<-ifelse(mod.free[,2]==getmode(mod.constraints[,2]),0,1)
mod.free[,3]<-ifelse(mod.free[,3]==getmode(mod.constraints[,3]),0,1)
mod.free

mod.pars<-mod.free

#partially invariant model
paste0("F.imm=~","(",paste0(rep("l1",18),collapse=","),")*gvrfgap.R",
       " + ",
       )


model2 <- "
F.imm=~ 
c(l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1,l1)*gvrfgap.R + 
c(l2,l2,l2,l2a,l2,l2,l2,l2,l2,l2,l2,l2b,l2,l2,l2c,l2,l2,l2d)*rfgfrpc + 
c(l3,l3,l3,l3a,l3,l3,l3b,l3c,l3d,l3,l3,l3,l3,l3,l3e,l3,l3,l3)*rfgbfml.R

F.imm~c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18)*1

F.imm~~c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18)*F.imm
m1==0
v1==1
"


metric <- cfa(model=model2,
           data=dat.imm,
           std.lv=F,
           #meanstructure=F,
           group=c("cntry"))
summary(metric)
fitMeasures(metric)

scores<-lavPredict(metric,type="lv",method = "regression")
str(scores)


## multigroup models return a list of factor scores (one per group)


idx <- lavInspect(metric, "case.idx") # list: 1 vector per group
fscores<-lavPredict(metric,type="lv",method = "regression")        # list: 1 matrix per group
## loop over groups and factors
for (g in seq_along(fscores)) {
  for (fs in colnames(fscores[[g]])) {
    dat.imm[ idx[[g]], fs] <- fscores[[g]][ , fs]
  }
}
head(dat.imm)

describeBy(dat.imm[,5],group=dat.imm[,"cntry"])
