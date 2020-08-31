---
title: "Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_depth: 5
---

# Preparations





## Load packages


```r
library(lme4)
library(lmerTest)
library(dplyr)
library(psych)
library(emmeans)
library(ggplot2)
library(metafor)
library(merTools)
```

## Session information about the packages


```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Finnish_Finland.1252  LC_CTYPE=Finnish_Finland.1252   
## [3] LC_MONETARY=Finnish_Finland.1252 LC_NUMERIC=C                    
## [5] LC_TIME=Finnish_Finland.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] merTools_0.5.0  arm_1.10-1      MASS_7.3-51.5   metafor_2.4-0  
##  [5] ggplot2_3.3.2   emmeans_1.4.6   psych_1.9.12.31 dplyr_0.8.5    
##  [9] lmerTest_3.1-2  lme4_1.1-23     Matrix_1.2-18  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6        mvtnorm_1.1-0       lattice_0.20-38    
##  [4] tidyr_1.1.0         zoo_1.8-7           assertthat_0.2.1   
##  [7] digest_0.6.25       foreach_1.5.0       mime_0.9           
## [10] R6_2.4.1            backports_1.1.6     evaluate_0.14      
## [13] coda_0.19-3         pillar_1.4.3        rlang_0.4.6        
## [16] multcomp_1.4-13     minqa_1.2.4         nloptr_1.2.2.1     
## [19] rmarkdown_2.1       splines_3.6.3       statmod_1.4.34     
## [22] stringr_1.4.0       munsell_0.5.0       shiny_1.4.0.2      
## [25] broom_0.5.6         httpuv_1.5.2        compiler_3.6.3     
## [28] numDeriv_2016.8-1.1 xfun_0.13           pkgconfig_2.0.3    
## [31] mnormt_1.5-6        htmltools_0.4.0     tidyselect_1.1.0   
## [34] tibble_3.0.1        codetools_0.2-16    later_1.0.0        
## [37] crayon_1.3.4        withr_2.2.0         grid_3.6.3         
## [40] nlme_3.1-144        xtable_1.8-4        gtable_0.3.0       
## [43] lifecycle_0.2.0     magrittr_1.5        scales_1.1.1       
## [46] estimability_1.3    stringi_1.4.6       promises_1.1.0     
## [49] ellipsis_0.3.1      vctrs_0.3.0         generics_0.0.2     
## [52] boot_1.3-24         sandwich_2.5-1      blme_1.0-4         
## [55] TH.data_1.0-10      iterators_1.0.12    tools_3.6.3        
## [58] glue_1.4.1          purrr_0.3.4         fastmap_1.0.1      
## [61] abind_1.4-5         parallel_3.6.3      survival_3.1-8     
## [64] yaml_2.2.1          colorspace_1.4-1    knitr_1.28
```

\newpage

## Custom functions


```r
#to extract fixed effects
getFE<-function(model){
  coefs<-data.frame(summary(model)$coefficients)
  coefs$lower<-coefs[,1]-qt(p=.975,df=coefs[,"df"])*coefs[,2]
  coefs$upper<-coefs[,1]+qt(p=.975,df=coefs[,"df"])*coefs[,2]
  coefs<-cbind.data.frame(round(coefs[,1:4],2),
                          p=round(coefs[,5],3),
                          LL=round(coefs$lower,2),
                          UL=round(coefs$upper,2))
  #row.names(coefs)<-substr(row.names(coefs),1,25)
  return(coefs)
}



#to extract random effects
getVC<-function(model){
  VC<-as.data.frame(VarCorr(model))
  VC<-cbind(VC[,c(1:3)],est_SD=VC[,5],est_SD2=VC[,4])
  return(VC)
}


#to extract model deviance
getDEV<-function(model){
  DEV<-unname(summary(model)$devcomp$cmp["dev"])
  return(DEV)
}


#partial correlation test
pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){


	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(use == "mat"){
		p.use <- "Var-Cov matrix"
		pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
	}else if(use == "rec"){
		p.use <- "Recursive formula"
		pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
	}else{
		stop("\'use\' should be either \"rec\" or \"mat\"!\n")
	}

	# print the method
	if(gregexpr("p",method)[[1]][1] == 1){
		p.method <- "Pearson"
	}else if(gregexpr("s",method)[[1]][1] == 1){
		p.method <- "Spearman"
	}else if(gregexpr("k",method)[[1]][1] == 1){
		p.method <- "Kendall"
	}else{
		stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
	}

	# sample number
	n <- dim(na.omit(data.frame(x,y,z)))[1]
	
	# given variables' number
	gn <- dim(z)[2]

	# p-value
	if(p.method == "Kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  		p.value <- 2*pnorm(-abs(statistic))
	}

	data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			


# By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	xdata <- na.omit(data.frame(data[,c(1,2)]))
	Sxx <- cov(xdata,xdata,m=method)

	xzdata <- na.omit(data)
	xdata <- data.frame(xzdata[,c(1,2)])
	zdata <- data.frame(xzdata[,-c(1,2)])
	Sxz <- cov(xdata,zdata,m=method)

	zdata <- na.omit(data.frame(data[,-c(1,2)]))
	Szz <- cov(zdata,zdata,m=method)

	# is Szz positive definite?
	zz.ev <- eigen(Szz)$values
	if(min(zz.ev)[1]<0){
		stop("\'Szz\' is not positive definite!\n")
	}

	# partial correlation
	Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
	
	rxx.z <- cov2cor(Sxx.z)[1,2]

	rxx.z
}

# By using recursive formula
pcor.rec <- function(x,y,z,method="p",na.rm=T){
	# 

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	# recursive formula
	if(dim(z)[2] == 1){
		tdata <- na.omit(data.frame(data[,1],data[,2]))
		rxy <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
		rxz <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
		ryz <- cor(tdata[,1],tdata[,2],m=method)

		rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
		
		return(rxy.z)
	}else{
		x <- c(data[,1])
		y <- c(data[,2])
		z0 <- c(data[,3])
		zc <- as.data.frame(data[,-c(1,2,3)])

		rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
		rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
		ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
		
		rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
		return(rxy.z)
	}			
}	
```

\newpage

## Load data


```r
dat<-read.csv2("alig.dat.csv",stringsAsFactors = F)
```

### Variable transformations

#### Country


```r
table(dat$cntry)
```

```
## 
##   AT   BE   CH   DE   EE   ES   FI   FR   GB   IE   IT   NL   NO   PT   SE   SI 
## 1973 1753 1503 2819 1974 1817 1862 2015 1876 2676 2317 1661 1538 1228 1525 1276
```

#### Voting group


```r
#make voting group variable names unique to each country
dat$voting.group<-paste0(dat$cntry,": ",dat$vote.group.combined)
```

#### Centering Attitudes towards the Environment


```r
#rename the variable
dat$environ<-dat$F.env
                                 
describe(dat$environ,fast=T)
```

```
##    vars     n mean   sd   min  max range se
## X1    1 28884    0 0.76 -2.95 1.91  4.85  0
```

```r
#grand mean center
dat$environ.gmc<-dat$environ-mean(dat$environ,na.rm=T)

#obtain dataframe with country means and add to data

environ.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(environ.cntry=mean(environ.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=environ.cntry,
               by=c("cntry"))

#center individuals around country means

dat$environ.cntrymc<-dat$environ.gmc-dat$environ.cntry

#obtain dataframe with voting group means and add to data

environ.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(environ.voting.group=mean(environ.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=environ.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$environ.vgmc<-dat$environ.cntrymc-dat$environ.voting.group

#describe the variable

describe(dat$environ.vgmc,fast=T)
```

```
##    vars     n mean   sd   min  max range se
## X1    1 28884    0 0.74 -3.04 2.23  5.27  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$environ.lvl1<-dat$environ.vgmc
dat$environ.lvl2<-dat$environ.voting.group
dat$environ.lvl3<-dat$environ.cntry
```

\newpage

#### Centering Political Engagement


```r
#correlation between the variables

corr.test(dat$nwspol.4,dat$polintr.R,adjust="none")
```

```
## Call:corr.test(x = dat$nwspol.4, y = dat$polintr.R, adjust = "none")
## Correlation matrix 
## [1] 0.31
## Sample Size 
## [1] 29637
## [1] 0
## 
##  To see confidence intervals of the correlations, print with the short=FALSE option
```

```r
#rename the variable

dat$engagement<-dat$polint.agg

#descriptive statistics
psych::describe(dat$engagement,fast=T)
```

```
##    vars     n mean   sd min max range se
## X1    1 29813 2.54 0.79   1   4     3  0
```

```r
#grand mean center
dat$engagement.gmc<-dat$engagement-mean(dat$engagement,na.rm=T)

#obtain dataframe with country means and add to data

engagement.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(engagement.cntry=mean(engagement.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=engagement.cntry,
               by=c("cntry"))

#center individuals around country means

dat$engagement.cntrymc<-dat$engagement.gmc-dat$engagement.cntry

#obtain dataframe with voting group means and add to data

engagement.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(engagement.voting.group=mean(engagement.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=engagement.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$engagement.vgmc<-dat$engagement.cntrymc-dat$engagement.voting.group

#describe the centered variable

describe(dat$engagement.vgmc,fast=T)
```

```
##    vars     n mean   sd   min max range se
## X1    1 29813    0 0.74 -2.06 2.1  4.17  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$engagement.lvl1<-dat$engagement.vgmc
dat$engagement.lvl2<-dat$engagement.voting.group
dat$engagement.lvl3<-dat$engagement.cntry
```

\newpage

#### Centering Political Interest Item (for exploratory analysis)


```r
#rename the variable (this will replace the previous item with same name)
dat$polintr<-dat$polintr.R

#descriptive statistics
psych::describe(dat$polintr,fast=T)
```

```
##    vars     n mean   sd min max range   se
## X1    1 29793  2.5 0.91   1   4     3 0.01
```

```r
#grand mean center
dat$polintr.gmc<-dat$polintr-mean(dat$polintr,na.rm=T)

#obtain dataframe with country means and add to data

polintr.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(polintr.cntry=mean(polintr.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=polintr.cntry,
               by=c("cntry"))

#center individuals around country means

dat$polintr.cntrymc<-dat$polintr.gmc-dat$polintr.cntry

#obtain dataframe with voting group means and add to data

polintr.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(polintr.voting.group=mean(polintr.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=polintr.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$polintr.vgmc<-dat$polintr.cntrymc-dat$polintr.voting.group

#describe the centered variable

describe(dat$polintr.vgmc,fast=T)
```

```
##    vars     n mean   sd   min  max range se
## X1    1 29793    0 0.83 -2.19 2.36  4.55  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$polintr.lvl1<-dat$polintr.vgmc
dat$polintr.lvl2<-dat$polintr.voting.group
dat$polintr.lvl3<-dat$polintr.cntry
```


\newpage

#### Centering "Time used for consuming political media" (for exploratory analysis)


```r
#rename the variable
dat$polnews<-dat$nwspol.4

#descriptive statistics
psych::describe(dat$polnews,fast=T)
```

```
##    vars     n mean   sd min max range   se
## X1    1 29657 2.59 1.04   1   4     3 0.01
```

```r
#grand mean center
dat$polnews.gmc<-dat$polnews-mean(dat$polnews,na.rm=T)

#obtain dataframe with country means and add to data

polnews.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(polnews.cntry=mean(polnews.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=polnews.cntry,
               by=c("cntry"))

#center individuals around country means

dat$polnews.cntrymc<-dat$polnews.gmc-dat$polnews.cntry

#obtain dataframe with voting group means and add to data

polnews.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(polnews.voting.group=mean(polnews.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=polnews.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$polnews.vgmc<-dat$polnews.cntrymc-dat$polnews.voting.group

#describe the centered variable

describe(dat$polnews.vgmc,fast=T)
```

```
##    vars     n mean sd   min  max range   se
## X1    1 29657    0  1 -2.19 2.31   4.5 0.01
```

```r
#rename as lvl1, lvl2, and lvl3

dat$polnews.lvl1<-dat$polnews.vgmc
dat$polnews.lvl2<-dat$polnews.voting.group
dat$polnews.lvl3<-dat$polnews.cntry
```

\newpage

#### Rename and grand mean center the Attitudes towards refugees (pro-refugee attitudes indicate high scores)


```r
#calculate the sum score
dat$refugees<-dat$F.imm
                                 
describe(dat$refugees,fast=T)
```

```
##    vars     n mean   sd   min  max range   se
## X1    1 27492    0 1.27 -3.96 3.96  7.92 0.01
```

```r
#grand mean center
dat$refugees<-dat$refugees-mean(dat$refugees,na.rm=T)

#rename
dat$refugees.gmc<-dat$refugees

#obtain dataframe with country means and add to data

refugees.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(refugees.cntry=mean(refugees.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=refugees.cntry,
               by=c("cntry"))

#center individuals around country means

dat$refugees.cntrymc<-dat$refugees.gmc-dat$refugees.cntry

#obtain dataframe with voting group means and add to data

refugees.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(refugees.voting.group=mean(refugees.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=refugees.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$refugees.vgmc<-dat$refugees.cntrymc-dat$refugees.voting.group

#describe the variable

describe(dat$refugees.vgmc,fast=T)
```

```
##    vars     n mean  sd   min  max range   se
## X1    1 27492    0 1.2 -4.13 4.35  8.48 0.01
```

```r
#rename as lvl1, lvl2, and lvl3

dat$refugees.lvl1<-dat$refugees.vgmc
dat$refugees.lvl2<-dat$refugees.voting.group
dat$refugees.lvl3<-dat$refugees.cntry
```

\newpage

#### Rename and Center the covariates aroung grand mean or logical middle points if applicable


```r
#grand-mean center age
dat$age<-dat$agea-mean(dat$agea,na.rm=T)
#sex around zero
dat$gender<-dat$gndr-1.5 #-0.5 males, 0.5 females
#rename occupation variable
dat$occup<-dat$isco.13
#grand-mean center education years
dat$educ<-dat$eduyrs-mean(dat$eduyrs,na.rm=T)
#residence around zero
dat$resid<-dat$rural-0.5 #-0.5 urban, 0.5 rural
```

\newpage

#### Voting group dummy-coded variables


```r
#recode if the party voted is =1, or not =0 anti-immigration
dat$anti.imm.party.dummy<-ifelse(is.na(dat$anti.imm.party.rule2),0,1)
#recode if the party voted is =1, or not =0 pro-environment
dat$pro.env.party.dummy<-ifelse(is.na(dat$pro.env.party.manual),0,1)

#dat$other.party.dummy<-ifelse(grepl("Other",dat$vote.group.combined),1,0)

#dummy-code not voting
dat$did.not.vote.dummy<-ifelse(grepl("did not vote",dat$vote.group.combined),1,0)
table(dat$did.not.vote.dummy)
```

```
## 
##     0     1 
## 24330  5483
```

```r
#dummy-code "don't know"
dat$dont.know.dummy<-ifelse(grepl("Don't know",dat$vote.group.combined),1,0)
table(dat$dont.know.dummy)
```

```
## 
##     0     1 
## 28722  1091
```

```r
#dummy-code invalid vote
dat$invalid.vote.dummy<-ifelse(grepl("Invalid vote",dat$vote.group.combined),1,0)
table(dat$invalid.vote.dummy)
```

```
## 
##     0     1 
## 29798    15
```

```r
#dummy-code "no answer"
dat$no.answer.dummy<-ifelse(grepl("No answer",dat$vote.group.combined),1,0)
table(dat$no.answer.dummy)
```

```
## 
##     0     1 
## 29801    12
```

```r
#dummy-code not-eligible: age
dat$not.eligible.age.dummy<-ifelse(grepl("not eligible: age",dat$vote.group.combined),1,0)
table(dat$not.eligible.age.dummy)
```

```
## 
##     0     1 
## 28390  1423
```

```r
#dummy code not-eligible: citizenship
dat$not.eligible.citizenship.dummy<-ifelse(grepl("not eligible: citizenship",dat$vote.group.combined),1,0)
table(dat$not.eligible.citizenship.dummy)
```

```
## 
##     0     1 
## 28611  1202
```

```r
#dummy-code not-eligible: other reasons
dat$not.eligible.other.dummy<-ifelse(grepl("not eligible: other",dat$vote.group.combined),1,0)
table(dat$not.eligible.other.dummy)
```

```
## 
##     0     1 
## 29594   219
```

```r
#add dummy-variable for other_party voting

dat<- dat %>%
  mutate(other.party.dummy:=case_when(
    anti.imm.party.dummy==1 |
      pro.env.party.dummy==1 |
      did.not.vote.dummy==1 |
      dont.know.dummy==1 |
      invalid.vote.dummy==1 |
      no.answer.dummy==1 |
      not.eligible.age.dummy==1 |
      not.eligible.citizenship.dummy==1 |
      not.eligible.other.dummy==1 ~0,
    TRUE~1
  ))

table(dat$other.party.dummy)
```

```
## 
##     0     1 
## 14441 15372
```

```r
#recode the names for a new multi-category variable: all.parties.lvl2

dat<-dat %>%
  mutate(all.parties.lvl2:=case_when(
    did.not.vote.dummy==1~"Did not vote",
    dont.know.dummy==1~"Don't know",
    no.answer.dummy==1~"No answer",
    invalid.vote.dummy==1~"Invalid vote",
    not.eligible.age.dummy==1~"NE age",
    not.eligible.citizenship.dummy==1~"NE citizen",
    not.eligible.other.dummy==1~"NE other",
    other.party.dummy==1~"Other party",
    anti.imm.party.dummy==1~"Anti-immigration party",
    pro.env.party.dummy==1~"Pro-environment party",
  ),
  party:=case_when(
    other.party.dummy==1~"Other party",
    anti.imm.party.dummy==1~"Anti-immigration party",
    pro.env.party.dummy==1~"Pro-environment party",
    TRUE~NA_character_
  ))
```


#### Omit missing variables


```r
#missing values per each row
dat$analysis.miss<-
  is.na(dat$cntry)+
  is.na(dat$voting.group)+
  is.na(dat$refugees)+
  is.na(dat$environ)+
  is.na(dat$vote.group.combined)+
  is.na(dat$age)+
  is.na(dat$gender)+
  is.na(dat$occup)+
  is.na(dat$educ)+
  is.na(dat$resid)+
  is.na(dat$engagement)
table(dat$analysis.miss)
```

```
## 
##     0     1     2 
## 26886  2604   323
```

```r
#include only those without any missing values
dat<-dat %>%
  filter(analysis.miss ==0)
```


# Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes

### Model 0: Intercepts only


```r
H1.mod0<-lmer(refugees~(1|voting.group),#+(1|cntry),
              data=dat,REML=F)

(FE.H1.mod0<-getFE(H1.mod0))
```

```
##             Estimate Std..Error     df t.value     p   LL   UL
## (Intercept)     0.08       0.03 251.81    2.44 0.015 0.01 0.14
```

```r
(VC.H1.mod0<-getVC(H1.mod0))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.4568751 0.2087348
## 2     Residual        <NA> <NA> 1.2016464 1.4439540
```

```r
getDEV(H1.mod0)
```

```
## [1] 86754.07
```

```r
#ICC

##voting group

VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.1263001
```

```r
##country

VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## numeric(0)
```


\newpage

### Model 1: Covariates


```r
H1.mod1<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

#model comparison
anova(H1.mod0,
      H1.mod1)
```

```
## Data: dat
## Models:
## H1.mod0: refugees ~ (1 | voting.group)
## H1.mod1: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.mod1:     occup
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod0    3 86760 86785 -43377    86754                         
## H1.mod1   19 85874 86030 -42918    85836 918.15 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod1<-getFE(H1.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.00       0.11
## occupUnemployed                                            -0.02       0.14
##                                                               df t.value     p
## (Intercept)                                             19271.69    0.28 0.778
## age                                                     26347.70   -2.70 0.007
## gender                                                  26757.71    6.49 0.000
## educ                                                    26811.38   14.59 0.000
## resid                                                   26849.56   -7.87 0.000
## occupClerical support workers                           26727.20   -0.17 0.863
## occupCraft and related trades workers                   26731.95   -1.28 0.202
## occupElementary occupations                             26736.38   -0.04 0.969
## occupManagers                                           26728.08    0.60 0.549
## occupOther: Not in paid work                            26834.31    1.58 0.115
## occupPlant and machine operators, and assemblers        26734.02   -0.65 0.517
## occupProfessionals                                      26724.72    1.75 0.080
## occupRetired                                            26732.28    0.22 0.830
## occupService and sales workers                          26727.59   -0.63 0.527
## occupSkilled agricultural, forestry and fishery workers 26730.64   -0.18 0.859
## occupTechnicians and associate professionals            26723.36    0.04 0.970
## occupUnemployed                                         26752.26   -0.13 0.898
##                                                            LL    UL
## (Intercept)                                             -0.19  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.22  0.27
## occupService and sales workers                          -0.28  0.14
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.28  0.25
```

```r
(VC.H1.mod1<-getVC(H1.mod1))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.4133723 0.1708767
## 2     Residual        <NA> <NA> 1.1820745 1.3973001
```

```r
getDEV(H1.mod1)
```

```
## [1] 85835.92
```

```r
write.csv2(FE.H1.mod1,"FE.H1.mod1.csv")

#variance explained

##lvl 1: individuals

(VC.H1.mod0[VC.H1.mod0$grp=="Residual","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="Residual","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="Residual","est_SD2"]
```

```
## [1] 0.03230982
```

```r
##lvl 2: voting group

(VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.1813696
```

```r
##lvl 3: country

(VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="cntry","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]
```

```
## numeric(0)
```

```r
##total

(sum(VC.H1.mod0$est_SD2)-sum(VC.H1.mod1$est_SD2))/
  sum(VC.H1.mod0$est_SD2)
```

```
## [1] 0.05113609
```

```r
#individual contributions of covariates
anova(H1.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age     10.178  10.178     1 26348   7.2843  0.006961 ** 
## gender  58.764  58.764     1 26758  42.0557 9.026e-11 ***
## educ   297.411 297.411     1 26811 212.8470 < 2.2e-16 ***
## resid   86.588  86.588     1 26850  61.9680 3.622e-15 ***
## occup  243.412  20.284    12 26690  14.5168 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage

### Model 2: Fixed effects for Attitudes towards the Environment


```r
H1.mod2<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1,data=dat,REML=F)

#model comparison
anova(H1.mod1,
      H1.mod2)
```

```
## Data: dat
## Models:
## H1.mod1: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.mod1:     occup
## H1.mod2: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.mod2:     occup + environ.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod1   19 85874 86030 -42918    85836                         
## H1.mod2   20 85196 85360 -42578    85156 679.54  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod2<-getFE(H1.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.02       0.13
## environ.lvl1                                                0.26       0.01
##                                                               df t.value     p
## (Intercept)                                             18787.02    0.35 0.725
## age                                                     26421.71   -1.65 0.099
## gender                                                  26754.01    5.82 0.000
## educ                                                    26825.81   12.34 0.000
## resid                                                   26844.75   -7.18 0.000
## occupClerical support workers                           26723.85   -0.25 0.803
## occupCraft and related trades workers                   26728.58   -1.20 0.231
## occupElementary occupations                             26732.79    0.09 0.930
## occupManagers                                           26724.75    0.45 0.654
## occupOther: Not in paid work                            26829.27    1.61 0.107
## occupPlant and machine operators, and assemblers        26730.45   -0.67 0.502
## occupProfessionals                                      26721.52    1.57 0.116
## occupRetired                                            26729.05    0.43 0.670
## occupService and sales workers                          26724.33   -0.60 0.549
## occupSkilled agricultural, forestry and fishery workers 26727.42   -0.13 0.894
## occupTechnicians and associate professionals            26720.11   -0.11 0.912
## occupUnemployed                                         26749.23    0.12 0.906
## environ.lvl1                                            26661.62   26.24 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.16  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.28
## environ.lvl1                                             0.24  0.28
```

```r
(VC.H1.mod2<-getVC(H1.mod2))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.4189837 0.1755473
## 2     Residual        <NA> <NA> 1.1669922 1.3618708
```

```r
getDEV(H1.mod2)
```

```
## [1] 85156.38
```

```r
write.csv2(FE.H1.mod2,"FE.H1.mod2.csv")

#variance explained

##lvl 1: individuals

(VC.H1.mod1[VC.H1.mod1$grp=="Residual","est_SD2"]-
     VC.H1.mod2[VC.H1.mod2$grp=="Residual","est_SD2"])/
  VC.H1.mod1[VC.H1.mod1$grp=="Residual","est_SD2"]
```

```
## [1] 0.02535559
```

```r
##total

(sum(VC.H1.mod1$est_SD2)-sum(VC.H1.mod2$est_SD2))/
  sum(VC.H1.mod1$est_SD2)
```

```
## [1] 0.01961433
```

\newpage

### Model 3: Random effects for Attitudes towards the Environment


```r
H1.mod3<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


#model comparison
anova(H1.mod2,
      H1.mod3)
```

```
## Data: dat
## Models:
## H1.mod2: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.mod2:     occup + environ.lvl1
## H1.mod3: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H1.mod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod2   20 85196 85360 -42578    85156                         
## H1.mod3   23 85078 85267 -42516    85032 124.29  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod3<-getFE(H1.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.27       0.03
##                                                               df t.value     p
## (Intercept)                                             18762.80    0.39 0.695
## age                                                     25948.59   -1.57 0.117
## gender                                                  26742.17    5.64 0.000
## educ                                                    26806.67   12.26 0.000
## resid                                                   26830.44   -7.22 0.000
## occupClerical support workers                           26695.53   -0.26 0.792
## occupCraft and related trades workers                   26702.49   -1.23 0.219
## occupElementary occupations                             26705.26    0.07 0.943
## occupManagers                                           26699.38    0.40 0.692
## occupOther: Not in paid work                            26800.07    1.58 0.113
## occupPlant and machine operators, and assemblers        26704.59   -0.72 0.471
## occupProfessionals                                      26694.69    1.54 0.124
## occupRetired                                            26700.61    0.37 0.709
## occupService and sales workers                          26696.12   -0.60 0.552
## occupSkilled agricultural, forestry and fishery workers 26704.73   -0.20 0.843
## occupTechnicians and associate professionals            26692.75   -0.14 0.891
## occupUnemployed                                         26720.62    0.08 0.940
## environ.lvl1                                               15.95    8.22 0.000
##                                                            LL    UL
## (Intercept)                                             -0.17  0.26
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.21  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.34
```

```r
(VC.H1.mod3<-getVC(H1.mod3))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41857780 0.175207375
## 2 voting.group environ.lvl1         <NA> 0.09120255 0.008317904
## 3 voting.group  (Intercept) environ.lvl1 0.63782258 0.024349107
## 4        cntry environ.lvl1         <NA> 0.11995018 0.014388045
## 5     Residual         <NA>         <NA> 1.16262160 1.351688977
```

```r
getDEV(H1.mod3)
```

```
## [1] 85032.09
```

```r
write.csv2(FE.H1.mod3,"FE.H1.mod3.csv")
```

\newpage


#### Describe the correlation between refugee and environment attitudes by country


```r
#model implied associations (posterior modes)

#fit also a model without covariates


slope.mod.no.cov<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

#correlations from data

cntry.cor.dat<-dat %>%
  dplyr::select(cntry,environ.gmc,refugees,age,gender,educ,resid) %>%
  group_by(cntry) %>%
  summarize(observed_r=cor(environ.gmc,refugees,use="pairwise.complete.obs"),
            partial_r=pcor.test(x=environ.gmc,
                         y=refugees,
                         z=data.frame(age,
                                      gender,
                                      educ,
                                      resid))$estimate,
            n=n())

#posterior modes from multilevel models

country.effs.dat<-
  cbind.data.frame(
    cntry=rownames(
      coefficients(H1.mod3)$cntry),
    slope_with_covs=coefficients(
      H1.mod3)$cntry[,"environ.lvl1"],
    slope_without_covs=coefficients(
      slope.mod.no.cov)$cntry[,"environ.lvl1"],
    cntry.cor.dat[,2:4])


country.effs.dat$partial_r.se<-1/sqrt(country.effs.dat$n-3)
country.effs.dat$country<-
  c("Austria",
    "Belgium",
    "Switzerland",
    #"Czech Republic",
    "Germany",
    "Estonia",
    "Spain",
    "Finland",
    "France",
    "Great Britain",
    #"Hungary",
    "Ireland",
    "Italy",
    #"Lithuania",
    "Netherlands",
    "Norway",
    #"Poland",
    "Portugal",
    "Sweden",
    "Slovenia")


country.effs.dat <- country.effs.dat[order(country.effs.dat$country),]

write.csv2(country.effs.dat,"associations_within_countries.csv")
```

\newpage

#### Describe the correlation between refugee and environment attitudes by voting group

* It was not possible to calculate the partial coefficients for all voting groups because of small group size and lack of variance in either attitudes or covariates, so only zero-order correlations were calculated as well as posterior mode slopes


```r
#get size of each voting group

voting.group.n<-dat %>%
  group_by(voting.group) %>%
  summarize(voting.group.n=n())

#add voting group size to the main data

dat<-left_join(x=dat,
               y=voting.group.n,
               by=c("voting.group"))

#calculate observed correlations


all.voting.group.cor.dat<-dat %>%
  #filter(voting.group.n>12) %>%
  dplyr::select(voting.group,environ.gmc,age,
                                      gender,
                                      educ,
                                      resid,
         refugees,anti.imm.party.dummy,
         pro.env.party.dummy) %>%
  group_by(voting.group) %>%
  summarize(observed_r=cor(environ.gmc,
                           refugees,use="pairwise.complete.obs"),

            n=n(),
            anti.imm=sum(anti.imm.party.dummy)/n(),
            pro.env=sum(pro.env.party.dummy)/n())

#posterior modes for voting group specific associations

voting.group.effs.dat<-
  cbind.data.frame(
    voting.group=rownames(
      coefficients(H1.mod3)$voting.group),
    slope_with_covs=coefficients(
      H1.mod3)$voting.group[,"environ.lvl1"],
    slope_without_covs=coefficients(
      slope.mod.no.cov)$voting.group[,"environ.lvl1"])

#add country level variability to slope-estimates

voting.group.effs.dat$cntry<-substr(voting.group.effs.dat$voting.group,1,2)
voting.group.effs.dat.w.country<-left_join(x=voting.group.effs.dat,
                                           y=country.effs.dat,
                                           by="cntry",suffix=c("",".cntry"))
```

```
## Warning: Column `cntry` joining character vector and factor, coercing into
## character vector
```

```r
voting.group.effs.dat.w.country$slope_with_covs_add_cntry<-
  voting.group.effs.dat.w.country$slope_with_covs+
  voting.group.effs.dat.w.country$slope_with_covs.cntry
  
voting.group.effs.dat.w.country$slope_without_covs_add_cntry<-
  voting.group.effs.dat.w.country$slope_without_covs+
  voting.group.effs.dat.w.country$slope_without_covs.cntry

voting.group.effs.dat.w.country<-
  voting.group.effs.dat.w.country %>%
  dplyr::select(voting.group,
                slope_with_covs,
                slope_without_covs,
                slope_with_covs_add_cntry,
                slope_without_covs_add_cntry,
                n)


voting.group.slopes.dat<-left_join(x=voting.group.effs.dat.w.country,
                                 y=all.voting.group.cor.dat,
                                 by=c("voting.group"))
```

```
## Warning: Column `voting.group` joining factor and character vector, coercing
## into character vector
```

```r
write.csv2(voting.group.slopes.dat,"associations_within_voting_groups.csv")
```


\newpage

##### Print a forest plot


```r
#Create forest plot of country estimates
library(dmetar)
```

```
## Extensive documentation for the dmetar package can be found at: 
##  www.bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/
```

```r
library(meta)
```

```
## Loading 'meta' package (version 4.13-0).
## Type 'help(meta)' for a brief overview.
```

```r
# country-level
#nrow(country.effs.dat)
#col.vect.cntry<-rep(c("#4CC542","#05B3D8"),each=10)

m <- metacor(cor=partial_r,
             n=n,
             data=country.effs.dat,
             studlab=country.effs.dat$country,
             comb.fixed = FALSE,
             comb.random = TRUE,
             prediction=TRUE,
             sm="ZCOR",
             backtransf=TRUE,
             level=.95)
m
```

```
##                   COR            95%-CI %W(random)
## Austria        0.2195 [ 0.1747; 0.2633]        6.3
## Belgium        0.1586 [ 0.1117; 0.2047]        6.3
## Estonia       -0.0150 [-0.0601; 0.0302]        6.3
## Finland        0.2875 [ 0.2446; 0.3292]        6.3
## France         0.2320 [ 0.1888; 0.2743]        6.3
## Germany        0.2636 [ 0.2279; 0.2986]        6.5
## Great Britain  0.2597 [ 0.2159; 0.3025]        6.3
## Ireland        0.1877 [ 0.1477; 0.2271]        6.4
## Italy          0.1335 [ 0.0901; 0.1764]        6.4
## Netherlands    0.1706 [ 0.1193; 0.2210]        6.2
## Norway         0.2698 [ 0.2222; 0.3161]        6.2
## Portugal       0.1488 [ 0.0912; 0.2054]        6.0
## Slovenia       0.0196 [-0.0383; 0.0774]        6.1
## Spain          0.1622 [ 0.1081; 0.2153]        6.1
## Sweden         0.2160 [ 0.1645; 0.2664]        6.2
## Switzerland    0.2167 [ 0.1651; 0.2671]        6.2
## 
## Number of studies combined: k = 16
## 
##                         COR           95%-CI    z  p-value
## Random effects model 0.1847 [0.1431; 0.2257] 8.56 < 0.0001
## Prediction interval         [0.0016; 0.3559]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0070 [0.0035; 0.0176]; tau = 0.0836 [0.0593; 0.1326];
##  I^2 = 92.1% [88.8%; 94.5%]; H = 3.56 [2.98; 4.25]
## 
## Test of heterogeneity:
##       Q d.f.  p-value
##  190.13   15 < 0.0001
## 
## Details on meta-analytical method:
## - Inverse variance method
## - DerSimonian-Laird estimator for tau^2
## - Jackson method for confidence interval of tau^2 and tau
## - Fisher's z transformation of correlations
```

```r
grDevices::pdf(file = "forestplot.pdf",family="sans") 

forest.m<-meta::forest(m,overall=T,
                       #layout = "JAMA",
                       prediction=F,
                       leftlabs=c("Country","n"),
                       print.I2=F,
                       print.tau2=F,
                       het.stat=F,
                       overall.hetstat=F,
                       text.random="Overall",
                       #weights=F,
                       #label.right="Partial Correlation Coefficient",
                       #bottom.lr=F
                       rightcols=c("effect", "ci"),
                       rightlabs = c("","95% CI"),
                       smlab = "Partial Correlation Coefficient",
                       weight.study="random"#,
                       #col.study=col.vect.cntry,
                       #col.square=col.vect.cntry
                       
)

graphics.off()


# print separate plots for voting groups within each country


cntry.vect<-as.character(unique(country.effs.dat$cntry))

#drop groups smaller than 6
voting.group.plot.dat<-voting.group.slopes.dat %>%
  filter(n.y>4)

voting.group.plot.dat$color<-ifelse(voting.group.plot.dat$anti.imm==1,
                                    "#05B3D8",
                                    ifelse(voting.group.plot.dat$pro.env==1,
                                           "#4CC542","darkgray"))

table(voting.group.plot.dat$color)
```

```
## 
##  #05B3D8  #4CC542 darkgray 
##       25       25      199
```

```r
pdf("forest_for_each_country.pdf",family = "sans",width = 12,height=10)

for (i in cntry.vect){
  
  dat.temp<-voting.group.plot.dat[grepl(paste0(i,": "),voting.group.plot.dat$voting.group),]
  
  m.temp <- metacor(cor=observed_r,
                    n=n.y,
                    data=dat.temp,
                    studlab=dat.temp[,"voting.group"],
                    comb.fixed = FALSE,
                    comb.random = TRUE,
                    prediction=TRUE,
                    sm="ZCOR",
                    backtransf=TRUE,
                    level=.95)
  
  forest.m.temp<-meta::forest(m.temp,overall=T,
                              prediction=F,
                              leftlabs=c("Voting group","n"),
                              print.I2=F,
                              print.tau2=F,
                              het.stat=F,
                              overall.hetstat=F,
                              text.random="Overall",
                              rightcols=c("effect", "ci"),
                              rightlabs = c("","95% CI"),
                              smlab = "Correlation Coefficient",
                              weight.study="random",
                              xlim=c(-0.8,0.8),
                              col.study=dat.temp[,"color"],
                              col.square=dat.temp[,"color"])
  
}


dev.off()
```

```
## null device 
##           1
```

```r
# print separate plots for anti-immigration parties


#drop groups smaller than 6 and filter only anti-refugee parties
anti.imm.group.plot.dat<-voting.group.slopes.dat %>%
  filter(n.y>4,anti.imm==1)

anti.imm.group.plot.dat$color<-"#05B3D8"

pdf("forest_for_anti_refugee.pdf",family = "sans",width = 12,height=10)

m.anti.imm <- metacor(cor=observed_r,
                    n=n.y,
                    data=anti.imm.group.plot.dat,
                    studlab=anti.imm.group.plot.dat[,"voting.group"],
                    comb.fixed = FALSE,
                    comb.random = TRUE,
                    prediction=TRUE,
                    sm="ZCOR",
                    backtransf=TRUE,
                    level=.95)

m.anti.imm
```

```
##                                                 COR             95%-CI
## AT: BZÖ                                      0.2746 [-0.6032;  0.8516]
## AT: FPÖ                                      0.2317 [ 0.1035;  0.3523]
## BE: N-VA                                    -0.0224 [-0.1451;  0.1010]
## BE: Parti Populaire                         -0.2968 [-0.9344;  0.7932]
## BE: Vlaams Belang                            0.2186 [-0.1538;  0.5366]
## CH: Swiss People's Party                     0.1508 [-0.0187;  0.3118]
## DE: AfD                                      0.0508 [-0.2126;  0.3073]
## DE: NPD                                      0.5766 [-0.1418;  0.8972]
## EE: Eesti Konservatiivne Rahvaerakond        0.1494 [-0.0887;  0.3713]
## ES: Partido Popular - PP                     0.0966 [-0.0312;  0.2213]
## FI: True Finns                               0.1808 [ 0.0390;  0.3153]
## FR: FN (Front National)                      0.1966 [ 0.0140;  0.3665]
## FR: MPF (Mouvement pour la France)           0.4571 [ 0.0037;  0.7546]
## FR: UMP (Union pour un Mouvement Populaire)  0.1545 [ 0.0375;  0.2673]
## GB: Conservative                             0.2588 [ 0.1729;  0.3407]
## GB: Democratic Unionist Party (nir)          0.3725 [-0.5290;  0.8790]
## GB: UK Independence Party                    0.0081 [-0.1857;  0.2013]
## IT: Fratelli d'Italia                        0.0918 [-0.3065;  0.4627]
## IT: Lega Nord                               -0.2438 [-0.4514; -0.0111]
## IT: Popolo delle Libertà (PdL)               0.2902 [ 0.0973;  0.4620]
## NL: Party for Freedom                        0.0408 [-0.1689;  0.2470]
## NL: Reformed Political Party                 0.0484 [-0.4149;  0.4918]
## NO: Progress Party (FRP)                     0.1817 [ 0.0026;  0.3496]
## SE: Sverigedemokraterna                      0.0396 [-0.1737;  0.2494]
## SI: SDS - Slovenska demokratska stranka      0.0140 [-0.1642;  0.1913]
##                                             %W(random)
## AT: BZÖ                                            0.3
## AT: FPÖ                                            7.0
## BE: N-VA                                           7.4
## BE: Parti Populaire                                0.2
## BE: Vlaams Belang                                  1.8
## CH: Swiss People's Party                           5.5
## DE: AfD                                            3.1
## DE: NPD                                            0.4
## EE: Eesti Konservatiivne Rahvaerakond              3.6
## ES: Partido Popular - PP                           7.2
## FI: True Finns                                     6.5
## FR: FN (Front National)                            5.0
## FR: MPF (Mouvement pour la France)                 1.1
## FR: UMP (Union pour un Mouvement Populaire)        7.6
## GB: Conservative                                   8.9
## GB: Democratic Unionist Party (nir)                0.3
## GB: UK Independence Party                          4.7
## IT: Fratelli d'Italia                              1.5
## IT: Lega Nord                                      3.6
## IT: Popolo delle Libertà (PdL)                     4.5
## NL: Party for Freedom                              4.3
## NL: Reformed Political Party                       1.1
## NO: Progress Party (FRP)                           5.1
## SE: Sverigedemokraterna                            4.2
## SI: SDS - Slovenska demokratska stranka            5.2
## 
## Number of studies combined: k = 25
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1280 [ 0.0741; 0.1811] 4.63 < 0.0001
## Prediction interval         [-0.0478; 0.2961]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0065 [0.0000; 0.0227]; tau = 0.0807 [0.0000; 0.1506];
##  I^2 = 40.9% [4.7%; 63.4%]; H = 1.30 [1.02; 1.65]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  40.62   24  0.0183
## 
## Details on meta-analytical method:
## - Inverse variance method
## - DerSimonian-Laird estimator for tau^2
## - Jackson method for confidence interval of tau^2 and tau
## - Fisher's z transformation of correlations
```

```r
forest.m.anti.imm<-meta::forest(m.anti.imm,overall=T,
                              prediction=F,
                              leftlabs=c("Voting group","n"),
                              print.I2=F,
                              print.tau2=F,
                              het.stat=F,
                              overall.hetstat=F,
                              text.random="Overall",
                              rightcols=c("effect", "ci"),
                              rightlabs = c("","95% CI"),
                              smlab = "Correlation Coefficient",
                              weight.study="random",
                              xlim=c(-0.8,0.8),
                              col.study=anti.imm.group.plot.dat[,"color"],
                              col.square=anti.imm.group.plot.dat[,"color"])
  

dev.off()
```

```
## null device 
##           1
```

```r
#drop groups smaller than 6 and filter only pro-environment parties
pro.env.group.plot.dat<-voting.group.slopes.dat %>%
  filter(n.y>4,pro.env==1)

pro.env.group.plot.dat$color<-"#4CC542"

pdf("forest_for_pro_environment.pdf",family = "sans",width = 12,height=10)

m.pro.env <- metacor(cor=observed_r,
                      n=n.y,
                      data=pro.env.group.plot.dat,
                      studlab=pro.env.group.plot.dat[,"voting.group"],
                      comb.fixed = FALSE,
                      comb.random = TRUE,
                      prediction=TRUE,
                      sm="ZCOR",
                      backtransf=TRUE,
                      level=.95)

m.pro.env
```

```
##                                                        COR            95%-CI
## AT: Grüne                                           0.4411 [ 0.3117; 0.5545]
## BE: Ecolo                                           0.1651 [-0.0973; 0.4061]
## BE: Groen!                                          0.0050 [-0.2192; 0.2287]
## CH: Green Party                                     0.1653 [-0.0948; 0.4042]
## CH: Social Democratic Party                         0.2520 [ 0.0792; 0.4101]
## DE: Bündnis 90/ Die Grünen                          0.1388 [ 0.0087; 0.2642]
## EE: Erakond Eestimaa Rohelised                      0.5466 [-0.0398; 0.8529]
## FI: Green League                                    0.3189 [ 0.1736; 0.4505]
## FR: Autres mouvements écologistes                  -0.2693 [-0.6360; 0.1967]
## FR: EELV (Europe Ecologie Les Verts)                0.2375 [-0.0246; 0.4691]
## GB: Green Party                                     0.2902 [-0.0476; 0.5685]
## IE: Green Party                                     0.0425 [-0.3763; 0.4469]
## IT: Movimento 5 Stelle                              0.1271 [-0.0085; 0.2580]
## IT: Sinistra Ecologia e Libertà (SEL)               0.4806 [ 0.0958; 0.7404]
## NL: Green Left                                      0.1772 [-0.0850; 0.4164]
## NL: Party for the Animals                           0.2571 [-0.1447; 0.5861]
## NO: Green Party (MDG)                               0.1821 [-0.1663; 0.4901]
## NO: Liberal Party (V)                              -0.0212 [-0.2977; 0.2587]
## NO: Socialist Left Party (SV)                       0.2948 [ 0.0442; 0.5105]
## PT: B.E. - Bloco de Esquerda                        0.2588 [ 0.0029; 0.4829]
## PT: PAN - Pessoas-Animais-Natureza                  0.2388 [-0.4215; 0.7336]
## SE: FI (Feministiskt initiativ)                     0.1285 [-0.2644; 0.4848]
## SE: Miljöpartiet de gröna                           0.3086 [ 0.0982; 0.4926]
## SE: Vänsterpartiet                                  0.3888 [ 0.1775; 0.5658]
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS) -0.2527 [-0.5830; 0.1493]
##                                                    %W(random)
## AT: Grüne                                                 7.2
## BE: Ecolo                                                 4.1
## BE: Groen!                                                5.0
## CH: Green Party                                           4.2
## CH: Social Democratic Party                               6.3
## DE: Bündnis 90/ Die Grünen                                8.0
## EE: Erakond Eestimaa Rohelised                            1.0
## FI: Green League                                          7.1
## FR: Autres mouvements écologistes                         1.7
## FR: EELV (Europe Ecologie Les Verts)                      4.1
## GB: Green Party                                           2.8
## IE: Green Party                                           2.0
## IT: Movimento 5 Stelle                                    7.8
## IT: Sinistra Ecologia e Libertà (SEL)                     2.0
## NL: Green Left                                            4.1
## NL: Party for the Animals                                 2.2
## NO: Green Party (MDG)                                     2.8
## NO: Liberal Party (V)                                     3.7
## NO: Socialist Left Party (SV)                             4.2
## PT: B.E. - Bloco de Esquerda                              4.2
## PT: PAN - Pessoas-Animais-Natureza                        0.9
## SE: FI (Feministiskt initiativ)                           2.3
## SE: Miljöpartiet de gröna                                 5.1
## SE: Vänsterpartiet                                        4.9
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)        2.2
## 
## Number of studies combined: k = 25
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.2144 [ 0.1491; 0.2778] 6.32 < 0.0001
## Prediction interval         [-0.0055; 0.4145]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0105 [0.0000; 0.0410]; tau = 0.1023 [0.0000; 0.2024];
##  I^2 = 40.7% [4.3%; 63.2%]; H = 1.30 [1.02; 1.65]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  40.45   24  0.0191
## 
## Details on meta-analytical method:
## - Inverse variance method
## - DerSimonian-Laird estimator for tau^2
## - Jackson method for confidence interval of tau^2 and tau
## - Fisher's z transformation of correlations
```

```r
forest.m.pro.env<-meta::forest(m.pro.env,overall=T,
                                prediction=F,
                                leftlabs=c("Voting group","n"),
                                print.I2=F,
                                print.tau2=F,
                                het.stat=F,
                                overall.hetstat=F,
                                text.random="Overall",
                                rightcols=c("effect", "ci"),
                                rightlabs = c("","95% CI"),
                                smlab = "Correlation Coefficient",
                                weight.study="random",
                                xlim=c(-0.8,0.8),
                                col.study=pro.env.group.plot.dat[,"color"],
                                col.square=pro.env.group.plot.dat[,"color"])


dev.off()
```

```
## null device 
##           1
```

\newpage

## Alternative (exploratory) approach for Hypothesis 1 with Environment attitudes as dependent variable, and immigrant attitudes as independent

### Model 0: Intercepts only


```r
H1.env.mod0<-lmer(environ.gmc~(1|voting.group),#+(1|cntry),
                  data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##             Estimate Std..Error     df t.value    p   LL   UL
## (Intercept)     0.03       0.01 219.95    2.61 0.01 0.01 0.06
```

```r
(VC.H1.env.mod0<-getVC(H1.env.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1750141 0.03062994
## 2     Residual        <NA> <NA> 0.7409473 0.54900283
```

```r
getDEV(H1.env.mod0)
```

```
## [1] 60567.31
```

```r
#ICC

##voting group

VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.0528437
```

```r
##country

VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## numeric(0)
```


\newpage

### Model 1: Covariates


```r
H1.env.mod1<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

#model comparison
anova(H1.env.mod0,
      H1.env.mod1)
```

```
## Data: dat
## Models:
## H1.env.mod0: environ.gmc ~ (1 | voting.group)
## H1.env.mod1: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.env.mod1:     occup
##             npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H1.env.mod0    3 60573 60598 -30284    60567                        
## H1.env.mod1   19 59726 59882 -29844    59688 879.3 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod1<-getFE(H1.env.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers            0.01       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.09       0.08
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.07       0.07
## occupUnemployed                                            -0.12       0.08
##                                                               df t.value     p
## (Intercept)                                             24922.59   -0.13 0.894
## age                                                     22321.00   -6.87 0.000
## gender                                                  26838.20    4.84 0.000
## educ                                                    25203.66   15.41 0.000
## resid                                                   26831.48   -5.11 0.000
## occupClerical support workers                           26800.75    0.61 0.539
## occupCraft and related trades workers                   26807.69   -0.47 0.640
## occupElementary occupations                             26817.30   -0.68 0.497
## occupManagers                                           26800.42    1.14 0.256
## occupOther: Not in paid work                            26881.47    0.06 0.949
## occupPlant and machine operators, and assemblers        26813.57    0.19 0.849
## occupProfessionals                                      26794.04    1.45 0.146
## occupRetired                                            26803.99   -1.16 0.247
## occupService and sales workers                          26799.57   -0.13 0.895
## occupSkilled agricultural, forestry and fishery workers 26801.96   -0.15 0.882
## occupTechnicians and associate professionals            26794.67    1.08 0.279
## occupUnemployed                                         26820.32   -1.40 0.160
##                                                            LL    UL
## (Intercept)                                             -0.14  0.12
## age                                                      0.00  0.00
## gender                                                   0.03  0.06
## educ                                                     0.02  0.03
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.16  0.10
## occupElementary occupations                             -0.18  0.09
## occupManagers                                           -0.06  0.21
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.12  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.24  0.06
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.15  0.13
## occupTechnicians and associate professionals            -0.06  0.20
## occupUnemployed                                         -0.28  0.05
```

```r
(VC.H1.env.mod1<-getVC(H1.env.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1473832 0.02172182
## 2     Residual        <NA> <NA> 0.7296591 0.53240246
```

```r
getDEV(H1.env.mod1)
```

```
## [1] 59688.01
```

```r
write.csv2(FE.H1.env.mod1,"FE.H1.env.mod1.csv")

#variance explained

##lvl 1: individuals

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="Residual","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="Residual","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="Residual","est_SD2"]
```

```
## [1] 0.03023731
```

```r
##lvl 2: voting group

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.2908304
```

```r
##lvl 3: country

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]
```

```
## numeric(0)
```

```r
##total

(sum(VC.H1.env.mod0$est_SD2)-sum(VC.H1.env.mod1$est_SD2))/
  sum(VC.H1.env.mod0$est_SD2)
```

```
## [1] 0.04400802
```

```r
#individual contributions of covariates
anova(H1.env.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age     25.098  25.098     1 22321  47.1418 6.777e-12 ***
## gender  12.470  12.470     1 26838  23.4228 1.307e-06 ***
## educ   126.411 126.411     1 25204 237.4347 < 2.2e-16 ***
## resid   13.913  13.913     1 26832  26.1324 3.210e-07 ***
## occup   55.892   4.658    12 26249   8.7484 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage

### Model 2: Fixed effects for Attitudes towards Immigrants


```r
H1.env.mod2<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                refugees.lvl1,data=dat,REML=F)

#model comparison
anova(H1.env.mod1,
      H1.env.mod2)
```

```
## Data: dat
## Models:
## H1.env.mod1: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.env.mod1:     occup
## H1.env.mod2: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.env.mod2:     occup + refugees.lvl1
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod1   19 59726 59882 -29844    59688                         
## H1.env.mod2   20 59054 59218 -29507    59014 674.01  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod2<-getFE(H1.env.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.00       0.07
## age                                                         0.00       0.00
## gender                                                      0.04       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.04       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.07       0.07
## occupOther: Not in paid work                               -0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.08       0.07
## occupRetired                                               -0.09       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.07       0.07
## occupUnemployed                                            -0.11       0.08
## refugees.lvl1                                               0.10       0.00
##                                                               df t.value     p
## (Intercept)                                             24676.53   -0.07 0.948
## age                                                     22838.64   -6.58 0.000
## gender                                                  26831.22    3.86 0.000
## educ                                                    25427.43   13.26 0.000
## resid                                                   26850.28   -3.95 0.000
## occupClerical support workers                           26792.99    0.64 0.521
## occupCraft and related trades workers                   26799.90   -0.28 0.778
## occupElementary occupations                             26809.20   -0.68 0.497
## occupManagers                                           26792.72    1.05 0.295
## occupOther: Not in paid work                            26884.56   -0.16 0.870
## occupPlant and machine operators, and assemblers        26805.54    0.29 0.773
## occupProfessionals                                      26786.50    1.19 0.234
## occupRetired                                            26796.49   -1.21 0.226
## occupService and sales workers                          26791.98   -0.04 0.972
## occupSkilled agricultural, forestry and fishery workers 26794.35   -0.13 0.897
## occupTechnicians and associate professionals            26786.97    1.08 0.278
## occupUnemployed                                         26813.28   -1.39 0.165
## refugees.lvl1                                           26679.13   26.13 0.000
##                                                            LL    UL
## (Intercept)                                             -0.13  0.13
## age                                                      0.00  0.00
## gender                                                   0.02  0.05
## educ                                                     0.02  0.02
## resid                                                   -0.06 -0.02
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.15  0.11
## occupElementary occupations                             -0.18  0.09
## occupManagers                                           -0.06  0.20
## occupOther: Not in paid work                            -0.15  0.12
## occupPlant and machine operators, and assemblers        -0.11  0.15
## occupProfessionals                                      -0.05  0.21
## occupRetired                                            -0.24  0.06
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.13
## occupTechnicians and associate professionals            -0.06  0.20
## occupUnemployed                                         -0.28  0.05
## refugees.lvl1                                            0.09  0.10
```

```r
(VC.H1.env.mod2<-getVC(H1.env.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1522856 0.02319091
## 2     Residual        <NA> <NA> 0.7203663 0.51892756
```

```r
getDEV(H1.env.mod2)
```

```
## [1] 59014
```

```r
write.csv2(FE.H1.env.mod2,"FE.H1.env.mod2.csv")

#variance explained

##lvl 1: individuals

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="Residual","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="Residual","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="Residual","est_SD2"]
```

```
## [1] 0.02530962
```

```r
##lvl 2: voting group

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]
```

```
## [1] -0.06763193
```

```r
##lvl 3: country

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="cntry","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]
```

```
## numeric(0)
```

```r
##total

(sum(VC.H1.env.mod1$est_SD2)-sum(VC.H1.env.mod2$est_SD2))/
  sum(VC.H1.env.mod1$est_SD2)
```

```
## [1] 0.02166629
```

\newpage

### Model 3: Random effects for Attitudes towards the Immigrant


```r
H1.env.mod3<-lmer(environ.gmc~(refugees.lvl1|voting.group)+
                (0+refugees.lvl1|cntry)+
                age+gender+educ+resid+occup+
                refugees.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


#model comparison
anova(H1.env.mod2,
      H1.env.mod3)
```

```
## Data: dat
## Models:
## H1.env.mod2: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H1.env.mod2:     occup + refugees.lvl1
## H1.env.mod3: environ.gmc ~ (refugees.lvl1 | voting.group) + (0 + refugees.lvl1 | 
## H1.env.mod3:     cntry) + age + gender + educ + resid + occup + refugees.lvl1
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod2   20 59054 59218 -29507    59014                         
## H1.env.mod3   23 58931 59120 -29443    58885 128.98  3  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod3<-getFE(H1.env.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.04       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.04       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.08       0.07
## occupRetired                                               -0.09       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.10       0.08
## refugees.lvl1                                               0.10       0.01
##                                                               df t.value     p
## (Intercept)                                             24635.10   -0.12 0.901
## age                                                     22778.47   -6.51 0.000
## gender                                                  26816.36    3.84 0.000
## educ                                                    25432.75   13.01 0.000
## resid                                                   26835.34   -4.23 0.000
## occupClerical support workers                           26774.32    0.70 0.485
## occupCraft and related trades workers                   26777.52   -0.23 0.816
## occupElementary occupations                             26788.60   -0.64 0.522
## occupManagers                                           26773.91    1.13 0.259
## occupOther: Not in paid work                            26864.19   -0.05 0.959
## occupPlant and machine operators, and assemblers        26786.79    0.36 0.717
## occupProfessionals                                      26766.41    1.22 0.221
## occupRetired                                            26766.99   -1.15 0.250
## occupService and sales workers                          26772.01    0.03 0.978
## occupSkilled agricultural, forestry and fishery workers 26773.71   -0.07 0.941
## occupTechnicians and associate professionals            26767.24    1.16 0.247
## occupUnemployed                                         26795.44   -1.25 0.212
## refugees.lvl1                                              15.93    8.02 0.000
##                                                            LL    UL
## (Intercept)                                             -0.14  0.12
## age                                                      0.00  0.00
## gender                                                   0.02  0.05
## educ                                                     0.02  0.02
## resid                                                   -0.06 -0.02
## occupClerical support workers                           -0.08  0.18
## occupCraft and related trades workers                   -0.15  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.06  0.21
## occupOther: Not in paid work                            -0.14  0.13
## occupPlant and machine operators, and assemblers        -0.11  0.16
## occupProfessionals                                      -0.05  0.21
## occupRetired                                            -0.24  0.06
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.14  0.13
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.26  0.06
## refugees.lvl1                                            0.07  0.13
```

```r
(VC.H1.env.mod3<-getVC(H1.env.mod3))
```

```
##            grp          var1          var2     est_SD      est_SD2
## 1 voting.group   (Intercept)          <NA> 0.15262567 0.0232945943
## 2 voting.group refugees.lvl1          <NA> 0.02877705 0.0008281183
## 3 voting.group   (Intercept) refugees.lvl1 0.29060987 0.0012763922
## 4        cntry refugees.lvl1          <NA> 0.04703414 0.0022122108
## 5     Residual          <NA>          <NA> 0.71755521 0.5148854828
```

```r
getDEV(H1.env.mod3)
```

```
## [1] 58885.02
```

```r
write.csv2(FE.H1.env.mod3,"FE.H1.env.mod3.csv")
```

