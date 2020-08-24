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
##   AT   BE   CH   DE   EE   ES   FI   GB   IE   IT   NL   NO   PT   SE   SI 
## 1973 1753 1503 2819 1974 1817 1862 1876 2676 2317 1661 1538 1228 1525 1276
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
##    vars     n mean   sd  min  max range   se
## X1    1 26887 -0.3 0.86 -3.6 1.67  5.27 0.01
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
## X1    1 26887    0 0.75 -3.01 2.47  5.48  0
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
## [1] 0.3
## Sample Size 
## [1] 27625
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
## X1    1 27798 2.54 0.79   1   4     3  0
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
## X1    1 27798    0 0.74 -2.06 2.1  4.17  0
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
## X1    1 27778 2.51 0.91   1   4     3 0.01
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
## X1    1 27778    0 0.82 -2.19 2.36  4.55  0
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
## X1    1 27645 2.59 1.04   1   4     3 0.01
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
## X1    1 27645    0  1 -2.19 2.31   4.5 0.01
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
## X1    1 25599 0.45 1.02 -2.23 2.93  5.15 0.01
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
##    vars     n mean   sd   min  max range   se
## X1    1 25599    0 0.81 -3.02 3.63  6.66 0.01
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
## 22913  4885
```

```r
#dummy-code "don't know"
dat$dont.know.dummy<-ifelse(grepl("Don't know",dat$vote.group.combined),1,0)
table(dat$dont.know.dummy)
```

```
## 
##     0     1 
## 26778  1020
```

```r
#dummy-code invalid vote
dat$invalid.vote.dummy<-ifelse(grepl("Invalid vote",dat$vote.group.combined),1,0)
table(dat$invalid.vote.dummy)
```

```
## 
##     0     1 
## 27783    15
```

```r
#dummy-code "no answer"
dat$no.answer.dummy<-ifelse(grepl("No answer",dat$vote.group.combined),1,0)
table(dat$no.answer.dummy)
```

```
## 
##     0     1 
## 27786    12
```

```r
#dummy-code not-eligible: age
dat$not.eligible.age.dummy<-ifelse(grepl("not eligible: age",dat$vote.group.combined),1,0)
table(dat$not.eligible.age.dummy)
```

```
## 
##     0     1 
## 26501  1297
```

```r
#dummy code not-eligible: citizenship
dat$not.eligible.citizenship.dummy<-ifelse(grepl("not eligible: citizenship",dat$vote.group.combined),1,0)
table(dat$not.eligible.citizenship.dummy)
```

```
## 
##     0     1 
## 26659  1139
```

```r
#dummy-code not-eligible: other reasons
dat$not.eligible.other.dummy<-ifelse(grepl("not eligible: other",dat$vote.group.combined),1,0)
table(dat$not.eligible.other.dummy)
```

```
## 
##     0     1 
## 27636   162
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
## 13004 14794
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
## 25004  2478   316
```

```r
#include only those without any missing values
dat<-dat %>%
  filter(analysis.miss ==0)
```


# Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes

### Model 0: Intercepts only


```r
H1.mod0<-lmer(refugees~(1|voting.group)+(1|cntry),
              data=dat,REML=F)

(FE.H1.mod0<-getFE(H1.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL  UL
## (Intercept)     0.09       0.15 14.99    0.59 0.563 -0.23 0.4
```

```r
(VC.H1.mod0<-getVC(H1.mod0))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.3295584 0.1086087
## 2        cntry (Intercept) <NA> 0.5592081 0.3127137
## 3     Residual        <NA> <NA> 0.8146453 0.6636470
```

```r
getDEV(H1.mod0)
```

```
## [1] 61322.27
```

```r
#ICC

##voting group

VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.100103
```

```r
##country

VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.2882235
```


\newpage

### Model 1: Covariates


```r
H1.mod1<-lmer(refugees~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

#model comparison
anova(H1.mod0,
      H1.mod1)
```

```
## Data: dat
## Models:
## H1.mod0: refugees ~ (1 | voting.group) + (1 | cntry)
## H1.mod1: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.mod1:     educ + resid + occup
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod0    4 61330 61363 -30661    61322                         
## H1.mod1   20 60562 60725 -30261    60522 800.22 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod1<-getFE(H1.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.06       0.16
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.08       0.01
## occupClerical support workers                              -0.02       0.08
## occupCraft and related trades workers                      -0.10       0.08
## occupElementary occupations                                 0.00       0.08
## occupManagers                                               0.05       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.13       0.08
## occupRetired                                                0.01       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals                0.00       0.08
## occupUnemployed                                            -0.01       0.09
##                                                               df t.value     p
## (Intercept)                                                24.50    0.37 0.711
## age                                                     24677.79   -2.69 0.007
## gender                                                  24867.92    5.87 0.000
## educ                                                    24989.31   12.52 0.000
## resid                                                   24901.20   -7.18 0.000
## occupClerical support workers                           24838.29   -0.25 0.802
## occupCraft and related trades workers                   24842.39   -1.32 0.188
## occupElementary occupations                             24842.22    0.02 0.986
## occupManagers                                           24839.36    0.58 0.561
## occupOther: Not in paid work                            24924.54    1.44 0.151
## occupPlant and machine operators, and assemblers        24842.98   -0.69 0.492
## occupProfessionals                                      24839.51    1.66 0.096
## occupRetired                                            24840.54    0.07 0.946
## occupService and sales workers                          24839.72   -0.81 0.419
## occupSkilled agricultural, forestry and fishery workers 24844.05   -0.13 0.900
## occupTechnicians and associate professionals            24835.97    0.01 0.995
## occupUnemployed                                         24853.70   -0.14 0.886
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.04  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.10 -0.06
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.26  0.05
## occupElementary occupations                             -0.15  0.16
## occupManagers                                           -0.11  0.20
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.02  0.28
## occupRetired                                            -0.17  0.18
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.15  0.15
## occupUnemployed                                         -0.20  0.17
```

```r
(VC.H1.mod1<-getVC(H1.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2998019 0.08988117
## 2        cntry (Intercept) <NA> 0.5521168 0.30483291
## 3     Residual        <NA> <NA> 0.8021827 0.64349707
```

```r
getDEV(H1.mod1)
```

```
## [1] 60522.05
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
## [1] 0.03036237
```

```r
##lvl 2: voting group

(VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.1724315
```

```r
##lvl 3: country

(VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="cntry","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]
```

```
## [1] 0.02520132
```

```r
##total

(sum(VC.H1.mod0$est_SD2)-sum(VC.H1.mod1$est_SD2))/
  sum(VC.H1.mod0$est_SD2)
```

```
## [1] 0.04309638
```

```r
#individual contributions of covariates
anova(H1.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age      4.665   4.665     1 24678   7.2491  0.007098 ** 
## gender  22.204  22.204     1 24868  34.5050 4.306e-09 ***
## educ   100.905 100.905     1 24989 156.8074 < 2.2e-16 ***
## resid   33.157  33.157     1 24901  51.5259 7.264e-13 ***
## occup  112.886   9.407    12 24824  14.6188 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage

### Model 2: Fixed effects for Attitudes towards the Environment


```r
H1.mod2<-lmer(refugees~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1,data=dat,REML=F)

#model comparison
anova(H1.mod1,
      H1.mod2)
```

```
## Data: dat
## Models:
## H1.mod1: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.mod1:     educ + resid + occup
## H1.mod2: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.mod2:     educ + resid + occup + environ.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod1   20 60562 60725 -30261    60522                         
## H1.mod2   21 59992 60162 -29975    59950 572.36  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod2<-getFE(H1.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.06       0.16
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.02       0.08
## occupCraft and related trades workers                      -0.09       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.04       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.12       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers     0.00       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.01       0.09
## environ.lvl1                                                0.16       0.01
##                                                               df t.value     p
## (Intercept)                                                24.27    0.38 0.708
## age                                                     24721.71   -1.79 0.074
## gender                                                  24865.24    5.33 0.000
## educ                                                    24987.26   10.46 0.000
## resid                                                   24897.56   -6.62 0.000
## occupClerical support workers                           24836.16   -0.26 0.794
## occupCraft and related trades workers                   24840.21   -1.18 0.239
## occupElementary occupations                             24840.00    0.17 0.861
## occupManagers                                           24837.23    0.49 0.623
## occupOther: Not in paid work                            24920.56    1.50 0.133
## occupPlant and machine operators, and assemblers        24840.67   -0.65 0.517
## occupProfessionals                                      24837.36    1.54 0.123
## occupRetired                                            24838.41    0.31 0.753
## occupService and sales workers                          24837.60   -0.72 0.469
## occupSkilled agricultural, forestry and fishery workers 24841.91   -0.02 0.986
## occupTechnicians and associate professionals            24833.87   -0.09 0.932
## occupUnemployed                                         24851.78    0.08 0.933
## environ.lvl1                                            24793.19   24.07 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.04  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.17
## occupManagers                                           -0.11  0.19
## occupOther: Not in paid work                            -0.04  0.28
## occupPlant and machine operators, and assemblers        -0.20  0.10
## occupProfessionals                                      -0.03  0.27
## occupRetired                                            -0.14  0.20
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.16  0.16
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.19
## environ.lvl1                                             0.15  0.18
```

```r
(VC.H1.mod2<-getVC(H1.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.3036170 0.09218326
## 2        cntry (Intercept) <NA> 0.5514434 0.30408984
## 3     Residual        <NA> <NA> 0.7929074 0.62870215
```

```r
getDEV(H1.mod2)
```

```
## [1] 59949.69
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
## [1] 0.02299143
```

```r
##total

(sum(VC.H1.mod1$est_SD2)-sum(VC.H1.mod2$est_SD2))/
  sum(VC.H1.mod1$est_SD2)
```

```
## [1] 0.01274876
```

\newpage

### Model 3: Random effects for Attitudes towards the Environment


```r
H1.mod3<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## H1.mod2: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.mod2:     educ + resid + occup + environ.lvl1
## H1.mod3: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H1.mod3:     age + gender + educ + resid + occup + environ.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod2   21 59992 60162 -29975    59950                         
## H1.mod3   25 59850 60053 -29900    59800 149.46  4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod3<-getFE(H1.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.07       0.16
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.02       0.08
## occupCraft and related trades workers                      -0.09       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.02       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.17       0.03
##                                                               df t.value     p
## (Intercept)                                                24.16    0.40 0.689
## age                                                     24472.93   -1.74 0.082
## gender                                                  24848.33    5.24 0.000
## educ                                                    24981.72   10.46 0.000
## resid                                                   24880.41   -6.71 0.000
## occupClerical support workers                           24806.23   -0.29 0.771
## occupCraft and related trades workers                   24814.17   -1.22 0.221
## occupElementary occupations                             24811.50    0.14 0.889
## occupManagers                                           24812.38    0.42 0.673
## occupOther: Not in paid work                            24892.44    1.46 0.145
## occupPlant and machine operators, and assemblers        24815.43   -0.70 0.482
## occupProfessionals                                      24811.25    1.48 0.139
## occupRetired                                            24809.13    0.27 0.788
## occupService and sales workers                          24810.57   -0.74 0.460
## occupSkilled agricultural, forestry and fishery workers 24815.14   -0.10 0.918
## occupTechnicians and associate professionals            24806.13   -0.12 0.904
## occupUnemployed                                         24824.02    0.00 0.999
## environ.lvl1                                               14.55    6.78 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.04  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.18
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.04  0.26
## occupRetired                                            -0.15  0.20
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.12  0.23
```

```r
(VC.H1.mod3<-getVC(H1.mod3))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.30340904 0.092057048
## 2 voting.group environ.lvl1         <NA> 0.06796968 0.004619877
## 3 voting.group  (Intercept) environ.lvl1 0.51590618 0.010639335
## 4        cntry  (Intercept)         <NA> 0.55222959 0.304957519
## 5        cntry environ.lvl1         <NA> 0.09172306 0.008413120
## 6        cntry  (Intercept) environ.lvl1 0.24468664 0.012393914
## 7     Residual         <NA>         <NA> 0.78911031 0.622695077
```

```r
getDEV(H1.mod3)
```

```
## [1] 59800.23
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
                (environ.lvl1|cntry)+
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
    #"France",
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
## Austria        0.1952 [ 0.1500; 0.2395]        6.7
## Belgium        0.1645 [ 0.1177; 0.2105]        6.7
## Estonia       -0.0189 [-0.0640; 0.0263]        6.7
## Finland        0.2866 [ 0.2437; 0.3284]        6.7
## Germany        0.2579 [ 0.2221; 0.2930]        6.9
## Great Britain  0.2510 [ 0.2070; 0.2939]        6.7
## Ireland        0.1844 [ 0.1444; 0.2239]        6.8
## Italy          0.1111 [ 0.0675; 0.1543]        6.8
## Netherlands    0.1698 [ 0.1185; 0.2202]        6.6
## Norway         0.3019 [ 0.2552; 0.3472]        6.6
## Portugal       0.1487 [ 0.0911; 0.2053]        6.5
## Slovenia       0.0182 [-0.0397; 0.0760]        6.5
## Spain          0.1465 [ 0.0921; 0.1999]        6.5
## Sweden         0.2140 [ 0.1624; 0.2644]        6.6
## Switzerland    0.2185 [ 0.1669; 0.2688]        6.6
## 
## Number of studies combined: k = 15
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1782 [ 0.1323; 0.2233] 7.50 < 0.0001
## Prediction interval         [-0.0201; 0.3630]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0080 [0.0040; 0.0208]; tau = 0.0895 [0.0631; 0.1442];
##  I^2 = 93.0% [90.0%; 95.1%]; H = 3.78 [3.17; 4.51]
## 
## Test of heterogeneity:
##       Q d.f.  p-value
##  199.79   14 < 0.0001
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
##       22       23      186
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
##                                             COR            95%-CI %W(random)
## AT: BZÖ                                  0.0359 [-0.7371; 0.7682]        0.3
## AT: FPÖ                                  0.2447 [ 0.1171; 0.3643]        8.2
## BE: N-VA                                -0.0186 [-0.1414; 0.1047]        8.6
## BE: Parti Populaire                     -0.3576 [-0.9425; 0.7665]        0.2
## BE: Vlaams Belang                        0.2188 [-0.1536; 0.5368]        2.0
## CH: Swiss People's Party                 0.1408 [-0.0289; 0.3025]        6.3
## DE: AfD                                  0.0242 [-0.2379; 0.2830]        3.5
## DE: NPD                                  0.6411 [-0.0401; 0.9154]        0.5
## EE: Eesti Konservatiivne Rahvaerakond    0.1247 [-0.1136; 0.3495]        4.1
## ES: Partido Popular - PP                 0.0872 [-0.0406; 0.2123]        8.4
## FI: True Finns                           0.1694 [ 0.0273; 0.3048]        7.6
## GB: Conservative                         0.2321 [ 0.1452; 0.3154]       10.6
## GB: Democratic Unionist Party (nir)      0.4390 [-0.4691; 0.8959]        0.3
## GB: UK Independence Party               -0.0145 [-0.2074; 0.1796]        5.4
## IT: Fratelli d'Italia                    0.1834 [-0.2195; 0.5329]        1.7
## IT: Lega Nord                           -0.1356 [-0.3576; 0.1009]        4.2
## IT: Popolo delle Libertà (PdL)           0.2532 [ 0.0577; 0.4300]        5.2
## NL: Party for Freedom                    0.0085 [-0.2001; 0.2164]        4.9
## NL: Reformed Political Party             0.1026 [-0.3688; 0.5320]        1.3
## NO: Progress Party (FRP)                 0.2586 [ 0.0833; 0.4185]        5.9
## SE: Sverigedemokraterna                  0.0426 [-0.1709; 0.2522]        4.8
## SI: SDS - Slovenska demokratska stranka -0.0048 [-0.1824; 0.1731]        6.0
## 
## Number of studies combined: k = 22
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1174 [ 0.0604; 0.1736] 4.02 < 0.0001
## Prediction interval         [-0.0542; 0.2822]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0060 [0.0000; 0.0244]; tau = 0.0772 [0.0000; 0.1562];
##  I^2 = 37.9% [0.0%; 62.8%]; H = 1.27 [1.00; 1.64]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  33.79   21  0.0381
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
## AT: Grüne                                           0.3902 [ 0.2550; 0.5104]
## BE: Ecolo                                           0.1782 [-0.0840; 0.4173]
## BE: Groen!                                          0.0014 [-0.2227; 0.2253]
## CH: Green Party                                     0.1534 [-0.1069; 0.3940]
## CH: Social Democratic Party                         0.2456 [ 0.0725; 0.4044]
## DE: Bündnis 90/ Die Grünen                          0.1254 [-0.0049; 0.2515]
## EE: Erakond Eestimaa Rohelised                      0.5566 [-0.0255; 0.8568]
## FI: Green League                                    0.3160 [ 0.1706; 0.4480]
## GB: Green Party                                     0.3826 [ 0.0566; 0.6349]
## IE: Green Party                                     0.0320 [-0.3852; 0.4385]
## IT: Movimento 5 Stelle                              0.1105 [-0.0253; 0.2422]
## IT: Sinistra Ecologia e Libertà (SEL)               0.3879 [-0.0184; 0.6842]
## NL: Green Left                                      0.1912 [-0.0706; 0.4283]
## NL: Party for the Animals                           0.2307 [-0.1720; 0.5673]
## NO: Green Party (MDG)                               0.1174 [-0.2298; 0.4382]
## NO: Liberal Party (V)                              -0.0039 [-0.2820; 0.2747]
## NO: Socialist Left Party (SV)                       0.2476 [-0.0067; 0.4719]
## PT: B.E. - Bloco de Esquerda                        0.2333 [-0.0242; 0.4618]
## PT: PAN - Pessoas-Animais-Natureza                  0.2806 [-0.3839; 0.7536]
## SE: FI (Feministiskt initiativ)                     0.1298 [-0.2632; 0.4859]
## SE: Miljöpartiet de gröna                           0.3021 [ 0.0911; 0.4872]
## SE: Vänsterpartiet                                  0.3868 [ 0.1752; 0.5642]
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS) -0.3630 [-0.6579; 0.0283]
##                                                    %W(random)
## AT: Grüne                                                 7.9
## BE: Ecolo                                                 4.3
## BE: Groen!                                                5.2
## CH: Green Party                                           4.4
## CH: Social Democratic Party                               6.9
## DE: Bündnis 90/ Die Grünen                                8.9
## EE: Erakond Eestimaa Rohelised                            1.0
## FI: Green League                                          7.8
## GB: Green Party                                           2.9
## IE: Green Party                                           2.0
## IT: Movimento 5 Stelle                                    8.6
## IT: Sinistra Ecologia e Libertà (SEL)                     2.1
## NL: Green Left                                            4.3
## NL: Party for the Animals                                 2.2
## NO: Green Party (MDG)                                     2.8
## NO: Liberal Party (V)                                     3.9
## NO: Socialist Left Party (SV)                             4.4
## PT: B.E. - Bloco de Esquerda                              4.4
## PT: PAN - Pessoas-Animais-Natureza                        0.9
## SE: FI (Feministiskt initiativ)                           2.3
## SE: Miljöpartiet de gröna                                 5.4
## SE: Vänsterpartiet                                        5.1
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)        2.2
## 
## Number of studies combined: k = 23
## 
##                         COR           95%-CI    z  p-value
## Random effects model 0.2081 [0.1434; 0.2710] 6.20 < 0.0001
## Prediction interval         [0.0052; 0.3945]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0086 [0.0000; 0.0387]; tau = 0.0930 [0.0000; 0.1967];
##  I^2 = 37.1% [0.0%; 61.9%]; H = 1.26 [1.00; 1.62]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  34.96   22  0.0391
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
H1.env.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.07       0.09 14.98    0.73 0.476 -0.13 0.26
```

```r
(VC.H1.env.mod0<-getVC(H1.env.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1899440 0.03607874
## 2        cntry (Intercept) <NA> 0.3472832 0.12060562
## 3     Residual        <NA> <NA> 0.7506191 0.56342898
```

```r
getDEV(H1.env.mod0)
```

```
## [1] 57054.11
```

```r
#ICC

##voting group

VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.05010147
```

```r
##country

VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.1674814
```


\newpage

### Model 1: Covariates


```r
H1.env.mod1<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

#model comparison
anova(H1.env.mod0,
      H1.env.mod1)
```

```
## Data: dat
## Models:
## H1.env.mod0: environ.gmc ~ (1 | voting.group) + (1 | cntry)
## H1.env.mod1: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.env.mod1:     educ + resid + occup
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod0    4 57062 57095 -28527    57054                         
## H1.env.mod1   20 56269 56432 -28115    56229 824.82 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod1<-getFE(H1.env.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## age                                                         0.00       0.00
## gender                                                      0.04       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.04       0.01
## occupClerical support workers                               0.01       0.07
## occupCraft and related trades workers                      -0.06       0.07
## occupElementary occupations                                -0.07       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                               -0.02       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.08       0.07
## occupRetired                                               -0.12       0.08
## occupService and sales workers                             -0.04       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.05       0.08
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.12       0.09
##                                                               df t.value     p
## (Intercept)                                                38.05    0.42 0.676
## age                                                     21667.53   -6.18 0.000
## gender                                                  24931.54    4.19 0.000
## educ                                                    24833.90   14.70 0.000
## resid                                                   24985.71   -4.33 0.000
## occupClerical support workers                           24884.73    0.18 0.856
## occupCraft and related trades workers                   24891.04   -0.90 0.368
## occupElementary occupations                             24892.83   -0.92 0.359
## occupManagers                                           24884.87    0.78 0.436
## occupOther: Not in paid work                            24989.36   -0.22 0.829
## occupPlant and machine operators, and assemblers        24895.01   -0.21 0.831
## occupProfessionals                                      24886.54    1.11 0.269
## occupRetired                                            24886.72   -1.51 0.131
## occupService and sales workers                          24886.18   -0.49 0.622
## occupSkilled agricultural, forestry and fishery workers 24891.60   -0.59 0.553
## occupTechnicians and associate professionals            24881.47    0.75 0.455
## occupUnemployed                                         24894.49   -1.38 0.169
##                                                            LL    UL
## (Intercept)                                             -0.18  0.28
## age                                                      0.00  0.00
## gender                                                   0.02  0.06
## educ                                                     0.02  0.03
## resid                                                   -0.06 -0.02
## occupClerical support workers                           -0.13  0.15
## occupCraft and related trades workers                   -0.21  0.08
## occupElementary occupations                             -0.21  0.08
## occupManagers                                           -0.09  0.20
## occupOther: Not in paid work                            -0.16  0.13
## occupPlant and machine operators, and assemblers        -0.16  0.13
## occupProfessionals                                      -0.06  0.22
## occupRetired                                            -0.28  0.04
## occupService and sales workers                          -0.17  0.10
## occupSkilled agricultural, forestry and fishery workers -0.20  0.11
## occupTechnicians and associate professionals            -0.09  0.19
## occupUnemployed                                         -0.29  0.05
```

```r
(VC.H1.env.mod1<-getVC(H1.env.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1590318 0.02529112
## 2        cntry (Intercept) <NA> 0.3484034 0.12138491
## 3     Residual        <NA> <NA> 0.7390657 0.54621805
```

```r
getDEV(H1.env.mod1)
```

```
## [1] 56229.29
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
## [1] 0.03054675
```

```r
##lvl 2: voting group

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.2990019
```

```r
##lvl 3: country

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]
```

```
## [1] -0.006461495
```

```r
##total

(sum(VC.H1.env.mod0$est_SD2)-sum(VC.H1.env.mod1$est_SD2))/
  sum(VC.H1.env.mod0$est_SD2)
```

```
## [1] 0.03779856
```

```r
#individual contributions of covariates
anova(H1.env.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age     20.864  20.864     1 21668  38.1968 6.510e-10 ***
## gender   9.574   9.574     1 24932  17.5277 2.841e-05 ***
## educ   118.039 118.039     1 24834 216.1026 < 2.2e-16 ***
## resid   10.240  10.240     1 24986  18.7473 1.498e-05 ***
## occup   55.679   4.640    12 24517   8.4946 2.533e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage

### Model 2: Fixed effects for Attitudes towards Immigrants


```r
H1.env.mod2<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup+
                refugees.lvl1,data=dat,REML=F)

#model comparison
anova(H1.env.mod1,
      H1.env.mod2)
```

```
## Data: dat
## Models:
## H1.env.mod1: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.env.mod1:     educ + resid + occup
## H1.env.mod2: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.env.mod2:     educ + resid + occup + refugees.lvl1
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod1   20 56269 56432 -28115    56229                         
## H1.env.mod2   21 55704 55875 -27831    55662 567.41  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod2<-getFE(H1.env.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## age                                                         0.00       0.00
## gender                                                      0.03       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.03       0.01
## occupClerical support workers                               0.02       0.07
## occupCraft and related trades workers                      -0.05       0.07
## occupElementary occupations                                -0.07       0.07
## occupManagers                                               0.05       0.07
## occupOther: Not in paid work                               -0.03       0.07
## occupPlant and machine operators, and assemblers           -0.01       0.07
## occupProfessionals                                          0.06       0.07
## occupRetired                                               -0.12       0.08
## occupService and sales workers                             -0.03       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.04       0.08
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.12       0.09
## refugees.lvl1                                               0.14       0.01
##                                                               df t.value     p
## (Intercept)                                                37.52    0.45 0.653
## age                                                     22070.85   -5.90 0.000
## gender                                                  24926.24    3.33 0.001
## educ                                                    24865.82   12.94 0.000
## resid                                                   24981.31   -3.32 0.001
## occupClerical support workers                           24879.93    0.22 0.829
## occupCraft and related trades workers                   24886.13   -0.72 0.471
## occupElementary occupations                             24887.70   -0.93 0.353
## occupManagers                                           24880.16    0.69 0.488
## occupOther: Not in paid work                            24987.46   -0.42 0.678
## occupPlant and machine operators, and assemblers        24889.81   -0.12 0.905
## occupProfessionals                                      24881.72    0.86 0.388
## occupRetired                                            24882.04   -1.54 0.125
## occupService and sales workers                          24881.43   -0.38 0.705
## occupSkilled agricultural, forestry and fishery workers 24886.71   -0.59 0.557
## occupTechnicians and associate professionals            24876.68    0.75 0.453
## occupUnemployed                                         24890.24   -1.36 0.173
## refugees.lvl1                                           24803.57   23.96 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.28
## age                                                      0.00  0.00
## gender                                                   0.01  0.05
## educ                                                     0.02  0.02
## resid                                                   -0.05 -0.01
## occupClerical support workers                           -0.12  0.16
## occupCraft and related trades workers                   -0.19  0.09
## occupElementary occupations                             -0.21  0.07
## occupManagers                                           -0.09  0.19
## occupOther: Not in paid work                            -0.17  0.11
## occupPlant and machine operators, and assemblers        -0.15  0.13
## occupProfessionals                                      -0.08  0.20
## occupRetired                                            -0.28  0.03
## occupService and sales workers                          -0.16  0.11
## occupSkilled agricultural, forestry and fishery workers -0.19  0.10
## occupTechnicians and associate professionals            -0.09  0.19
## occupUnemployed                                         -0.29  0.05
## refugees.lvl1                                            0.13  0.15
```

```r
(VC.H1.env.mod2<-getVC(H1.env.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1642361 0.02697348
## 2        cntry (Intercept) <NA> 0.3474618 0.12072973
## 3     Residual        <NA> <NA> 0.7305329 0.53367837
```

```r
getDEV(H1.env.mod2)
```

```
## [1] 55661.88
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
## [1] 0.02295729
```

```r
##lvl 2: voting group

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]
```

```
## [1] -0.06651973
```

```r
##lvl 3: country

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="cntry","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]
```

```
## [1] 0.00539753
```

```r
##total

(sum(VC.H1.env.mod1$est_SD2)-sum(VC.H1.env.mod2$est_SD2))/
  sum(VC.H1.env.mod1$est_SD2)
```

```
## [1] 0.0166151
```

\newpage

### Model 3: Random effects for Attitudes towards the Immigrant


```r
H1.env.mod3<-lmer(environ.gmc~(refugees.lvl1|voting.group)+
                (refugees.lvl1|cntry)+
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
## H1.env.mod2: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H1.env.mod2:     educ + resid + occup + refugees.lvl1
## H1.env.mod3: environ.gmc ~ (refugees.lvl1 | voting.group) + (refugees.lvl1 | 
## H1.env.mod3:     cntry) + age + gender + educ + resid + occup + refugees.lvl1
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod2   21 55704 55875 -27831    55662                         
## H1.env.mod3   25 55596 55799 -27773    55546 115.78  4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod3<-getFE(H1.env.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## age                                                         0.00       0.00
## gender                                                      0.03       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.03       0.01
## occupClerical support workers                               0.02       0.07
## occupCraft and related trades workers                      -0.05       0.07
## occupElementary occupations                                -0.06       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                               -0.02       0.07
## occupPlant and machine operators, and assemblers            0.00       0.07
## occupProfessionals                                          0.06       0.07
## occupRetired                                               -0.12       0.08
## occupService and sales workers                             -0.02       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.04       0.08
## occupTechnicians and associate professionals                0.06       0.07
## occupUnemployed                                            -0.10       0.09
## refugees.lvl1                                               0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                37.34    0.40 0.688
## age                                                     22080.01   -5.93 0.000
## gender                                                  24911.45    3.23 0.001
## educ                                                    24825.57   12.68 0.000
## resid                                                   24967.99   -3.49 0.000
## occupClerical support workers                           24855.68    0.30 0.762
## occupCraft and related trades workers                   24858.75   -0.66 0.507
## occupElementary occupations                             24862.11   -0.88 0.377
## occupManagers                                           24857.28    0.78 0.434
## occupOther: Not in paid work                            24961.15   -0.28 0.778
## occupPlant and machine operators, and assemblers        24864.76   -0.03 0.977
## occupProfessionals                                      24855.77    0.92 0.358
## occupRetired                                            24847.24   -1.46 0.144
## occupService and sales workers                          24856.06   -0.30 0.763
## occupSkilled agricultural, forestry and fishery workers 24860.15   -0.53 0.599
## occupTechnicians and associate professionals            24851.84    0.86 0.392
## occupUnemployed                                         24870.01   -1.21 0.225
## refugees.lvl1                                              14.12    7.05 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.28
## age                                                      0.00  0.00
## gender                                                   0.01  0.05
## educ                                                     0.02  0.02
## resid                                                   -0.05 -0.02
## occupClerical support workers                           -0.12  0.16
## occupCraft and related trades workers                   -0.19  0.09
## occupElementary occupations                             -0.20  0.08
## occupManagers                                           -0.08  0.20
## occupOther: Not in paid work                            -0.16  0.12
## occupPlant and machine operators, and assemblers        -0.14  0.14
## occupProfessionals                                      -0.07  0.20
## occupRetired                                            -0.28  0.04
## occupService and sales workers                          -0.16  0.12
## occupSkilled agricultural, forestry and fishery workers -0.19  0.11
## occupTechnicians and associate professionals            -0.08  0.20
## occupUnemployed                                         -0.27  0.06
## refugees.lvl1                                            0.10  0.19
```

```r
(VC.H1.env.mod3<-getVC(H1.env.mod3))
```

```
##            grp          var1          var2     est_SD      est_SD2
## 1 voting.group   (Intercept)          <NA> 0.16474058 0.0271394599
## 2 voting.group refugees.lvl1          <NA> 0.04773257 0.0022783983
## 3 voting.group   (Intercept) refugees.lvl1 0.11732372 0.0009225741
## 4        cntry   (Intercept)          <NA> 0.34759066 0.1208192693
## 5        cntry refugees.lvl1          <NA> 0.07595944 0.0057698370
## 6        cntry   (Intercept) refugees.lvl1 0.05380809 0.0014206839
## 7     Residual          <NA>          <NA> 0.72759257 0.5293909471
```

```r
getDEV(H1.env.mod3)
```

```
## [1] 55546.1
```

```r
write.csv2(FE.H1.env.mod3,"FE.H1.env.mod3.csv")
```

