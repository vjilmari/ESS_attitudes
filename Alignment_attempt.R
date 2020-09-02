#' ---
#' title: "Alignment attempt"
#' output: 
#'   html_document: 
#'     keep_md: yes
#'     number_sections: yes
#'     toc: yes
#'     toc_depth: 5
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Preparations
#' 
#' ### Load packages
#' 
## --------------------------------------------------------------------------------------------------------------
library(sirt)
library(psych)
library(dplyr)
library(lavaan)

#' 
#' ### Read data
#' 
## --------------------------------------------------------------------------------------------------------------
dat<-read.csv2("dat.no.miss.csv",stringsAsFactors = F)

# data for refugee attitudes
imm.vars<-c("gvrfgap.R","rfgfrpc","rfgbfml.R")

# standardization of the indicators (as suggested by Asparouhov & Muthen, 2014)
na.standardize<-function(var){
  (var-mean(var,na.rm=T))/sd(var,na.rm=T)
}

dat$imm.1.z<-na.standardize(dat$gvrfgap.R)
dat$imm.2.z<-na.standardize(dat$rfgfrpc)
dat$imm.3.z<-na.standardize(dat$rfgbfml.R)

imm.z.vars<-c("imm.1.z","imm.2.z","imm.3.z")

describe(dat[,imm.z.vars])
describeBy(dat[,imm.z.vars],group=dat[,"cntry"])

#exclude Hungary
dat.imm<-dat %>%
  filter(cntry!="HU") %>%
  dplyr::select(cntry,all_of(imm.z.vars),all_of(imm.vars))
dat.imm<-na.omit(dat.imm)


env.vars<-c("inctxff.R","sbsrnen.R","banhhap.R")

dat$env.1.z<-na.standardize(dat$inctxff.R)
dat$env.2.z<-na.standardize(dat$sbsrnen.R)
dat$env.3.z<-na.standardize(dat$banhhap.R)

env.z.vars<-c("env.1.z","env.2.z","env.3.z")

describe(dat[,env.z.vars])
describeBy(dat[,env.z.vars],group=dat[,"cntry"])


dat.env<-dat %>%
  dplyr::select(cntry,all_of(env.z.vars),all_of(env.vars))
dat.env<-na.omit(dat.env)




#' 
#' \newpage
#' 
#' # Alignment for refugee attitudes
#' 
#' Approach by Fischer and Karl (2019)
#' 
#' ### Configural models
#' 
## --------------------------------------------------------------------------------------------------------------
par.imm <- invariance_alignment_cfa_config(dat = dat.imm[,imm.z.vars],
                                       group = dat.imm[,"cntry"])
par.imm

#' 
#' For Lithuania, solution was not found. For Poland and Czech, there are negative variances. Also for France, the model seems to produce Heywood loadings and error variances (>1)
#' 
#' \newpage
#' 
#' ### Construct separate configural models for the countries with problems
#' 
## --------------------------------------------------------------------------------------------------------------
conf.mod<-"
F.imm=~imm.1.z+imm.2.z+imm.3.z
"

#' 
#' \newpage
#' 
#' #### Configural model for Lithuania
#' 
## --------------------------------------------------------------------------------------------------------------
dat.imm.LT<-dat.imm %>%
  filter(cntry=="LT") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.LT<-cfa(data=dat.imm.LT,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(imm.conf.LT)

#' 
#' \newpage
#' 
#' #### Configural model for Poland
#' 
## --------------------------------------------------------------------------------------------------------------
dat.imm.PL<-dat.imm %>%
  filter(cntry=="PL") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.PL<-cfa(data=dat.imm.PL,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(imm.conf.PL)

#' 
#' \newpage
#' 
#' #### Configural model for Czech Republic
#' 
## --------------------------------------------------------------------------------------------------------------
dat.imm.CZ<-dat.imm %>%
  filter(cntry=="CZ") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.CZ<-cfa(data=dat.imm.CZ,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(imm.conf.CZ)

#' 
#' \newpage
#' 
#' #### Configural model for France
#' 
## --------------------------------------------------------------------------------------------------------------
dat.imm.FR<-dat.imm %>%
  filter(cntry=="FR") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.FR<-cfa(data=dat.imm.FR,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(imm.conf.FR)

#' 
#' 
#' \newpage
#' 
#' 
#' # Alignment for environment attitudes
#' 
#' 
#' ### Configural models
#' 
## --------------------------------------------------------------------------------------------------------------
par.env <- invariance_alignment_cfa_config(dat = dat.env[,env.z.vars],
                                       group = dat.env[,"cntry"])
par.env

#' 
#' Configural model does not converge for Hungary. For Czech Republic, one error variance seems to be > 1.
#' 
#' 
#' \newpage
#' 
#' ### Construct separate configural models for the countries with problems
#' 
## --------------------------------------------------------------------------------------------------------------
conf.mod<-"
F.env=~env.1.z+env.2.z+env.3.z
"

#' 
#' \newpage
#' 
#' #### Configural model for Hungary
#' 
## --------------------------------------------------------------------------------------------------------------
dat.env.HU<-dat.env %>%
  filter(cntry=="HU") %>%
  mutate(env.1.z=na.standardize(inctxff.R),
         env.2.z=na.standardize(sbsrnen.R),
         env.3.z=na.standardize(banhhap.R))


env.conf.HU<-cfa(data=dat.env.HU,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(env.conf.HU)

#' 
#' \newpage
#' 
#' #### Configural model for Czech Republic
#' 
## --------------------------------------------------------------------------------------------------------------
dat.env.CZ<-dat.env %>%
  filter(cntry=="CZ") %>%
  mutate(env.1.z=na.standardize(inctxff.R),
         env.2.z=na.standardize(sbsrnen.R),
         env.3.z=na.standardize(banhhap.R))


env.conf.CZ<-cfa(data=dat.env.CZ,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(env.conf.CZ)

#' 
