---
title: "Exploratory non-linearity"
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


# Exploratory analysis: Is the association linear?

## Data re-preparations


```r
##Remove two of the smallest voting groups

Rdat<-dat %>%
  filter(all.parties.lvl2!="Invalid vote" & all.parties.lvl2!="No answer")
```

## Model 1 (H1 selected model without random effect correlations)


```r
EX9RR.NCmod1<-lmer(refugees~
                     (environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1,data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod1)
```

```
## [1] FALSE
```

```r
(VC.EX9RR.NCmod1<-getVC(EX9RR.NCmod1))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.30461158 0.092788214
## 2 voting.group.1 environ.lvl1 <NA> 0.07051058 0.004971742
## 3          cntry  (Intercept) <NA> 0.55259526 0.305361518
## 4        cntry.1 environ.lvl1 <NA> 0.08975493 0.008055948
## 5       Residual         <NA> <NA> 0.78921911 0.622866807
```

```r
getFE(EX9RR.NCmod1)
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
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.05       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.17       0.03
##                                                               df t.value     p
## (Intercept)                                                24.15    0.39 0.702
## age                                                     24702.48   -1.72 0.085
## gender                                                  24822.96    5.21 0.000
## educ                                                    24962.23   10.52 0.000
## resid                                                   24856.71   -6.71 0.000
## occupClerical support workers                           24791.11   -0.27 0.785
## occupCraft and related trades workers                   24799.23   -1.21 0.228
## occupElementary occupations                             24796.41    0.16 0.876
## occupManagers                                           24797.03    0.45 0.653
## occupOther: Not in paid work                            24876.65    1.48 0.140
## occupPlant and machine operators, and assemblers        24799.62   -0.68 0.499
## occupProfessionals                                      24797.21    1.51 0.131
## occupRetired                                            24793.55    0.31 0.753
## occupService and sales workers                          24796.18   -0.72 0.471
## occupSkilled agricultural, forestry and fishery workers 24799.28   -0.07 0.945
## occupTechnicians and associate professionals            24791.24   -0.10 0.923
## occupUnemployed                                         24806.90    0.02 0.988
## environ.lvl1                                               14.54    6.73 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.04  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.19
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.20  0.10
## occupProfessionals                                      -0.03  0.26
## occupRetired                                            -0.14  0.20
## occupService and sales workers                          -0.20  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.16
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.12  0.22
```

\newpage

## Model 2 (add squared fixed effect)


```r
EX9RR.NCmod2<-lmer(refugees~(environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2),data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))


anova(EX9RR.NCmod1,EX9RR.NCmod2)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod1: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1
## EX9RR.NCmod2: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod2:     I(environ.lvl1^2)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX9RR.NCmod1   23 59832 60019 -29893    59786                     
## EX9RR.NCmod2   24 59834 60029 -29893    59786 0.0789  1     0.7789
```

```r
isSingular(EX9RR.NCmod2)
```

```
## [1] FALSE
```

```r
(VC.EX9RR.NCmod2<-getVC(EX9RR.NCmod2))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.30468277 0.092831593
## 2 voting.group.1 environ.lvl1 <NA> 0.07054960 0.004977246
## 3          cntry  (Intercept) <NA> 0.55249814 0.305254193
## 4        cntry.1 environ.lvl1 <NA> 0.08973973 0.008053220
## 5       Residual         <NA> <NA> 0.78921571 0.622861432
```

```r
getFE(EX9RR.NCmod2)
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
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.17       0.03
## I(environ.lvl1^2)                                           0.00       0.01
##                                                               df t.value     p
## (Intercept)                                                24.17    0.38 0.706
## age                                                     24703.98   -1.73 0.084
## gender                                                  24822.93    5.21 0.000
## educ                                                    24962.51   10.51 0.000
## resid                                                   24855.91   -6.70 0.000
## occupClerical support workers                           24791.07   -0.27 0.785
## occupCraft and related trades workers                   24799.16   -1.21 0.228
## occupElementary occupations                             24796.40    0.15 0.877
## occupManagers                                           24797.02    0.45 0.653
## occupOther: Not in paid work                            24876.77    1.47 0.140
## occupPlant and machine operators, and assemblers        24799.55   -0.68 0.499
## occupProfessionals                                      24797.28    1.51 0.131
## occupRetired                                            24793.65    0.31 0.755
## occupService and sales workers                          24796.13   -0.72 0.471
## occupSkilled agricultural, forestry and fishery workers 24799.28   -0.07 0.944
## occupTechnicians and associate professionals            24791.18   -0.10 0.923
## occupUnemployed                                         24806.87    0.01 0.989
## environ.lvl1                                               14.83    6.72 0.000
## I(environ.lvl1^2)                                       24618.66    0.28 0.779
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.04  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.19
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.20  0.10
## occupProfessionals                                      -0.03  0.26
## occupRetired                                            -0.14  0.20
## occupService and sales workers                          -0.20  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.16
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.12  0.22
## I(environ.lvl1^2)                                       -0.01  0.01
```

\newpage

## Model 3 (add fixed cubic effect)


```r
EX9RR.NCmod3<-lmer(refugees~(environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
            I(environ.lvl1^2)+
              I(environ.lvl1^3),
            data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod3)
```

```
## [1] FALSE
```

```r
#test for cubic term only
anova(EX9RR.NCmod2,EX9RR.NCmod3)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod2: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod2:     I(environ.lvl1^2)
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX9RR.NCmod2   24 59834 60029 -29893    59786                     
## EX9RR.NCmod3   25 59834 60037 -29892    59784 2.1119  1     0.1462
```

```r
#test for both non-linear terms
anova(EX9RR.NCmod1,EX9RR.NCmod3)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod1: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX9RR.NCmod1   23 59832 60019 -29893    59786                     
## EX9RR.NCmod3   25 59834 60037 -29892    59784 2.1908  2     0.3344
```

```r
(VC.EX9RR.NCmod3<-getVC(EX9RR.NCmod3))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.30486956 0.092945448
## 2 voting.group.1 environ.lvl1 <NA> 0.07102279 0.005044237
## 3          cntry  (Intercept) <NA> 0.55220002 0.304924859
## 4        cntry.1 environ.lvl1 <NA> 0.09059896 0.008208172
## 5       Residual         <NA> <NA> 0.78916354 0.622779089
```

```r
getFE(EX9RR.NCmod3)
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
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.16       0.03
## I(environ.lvl1^2)                                           0.01       0.01
## I(environ.lvl1^3)                                           0.01       0.01
##                                                               df t.value     p
## (Intercept)                                                24.19    0.37 0.713
## age                                                     24704.92   -1.73 0.084
## gender                                                  24823.52    5.19 0.000
## educ                                                    24962.15   10.51 0.000
## resid                                                   24855.83   -6.69 0.000
## occupClerical support workers                           24790.68   -0.28 0.780
## occupCraft and related trades workers                   24798.73   -1.21 0.225
## occupElementary occupations                             24796.04    0.15 0.881
## occupManagers                                           24796.61    0.45 0.656
## occupOther: Not in paid work                            24876.38    1.47 0.142
## occupPlant and machine operators, and assemblers        24799.15   -0.68 0.495
## occupProfessionals                                      24797.00    1.51 0.131
## occupRetired                                            24793.27    0.31 0.756
## occupService and sales workers                          24795.73   -0.73 0.466
## occupSkilled agricultural, forestry and fishery workers 24798.96   -0.08 0.940
## occupTechnicians and associate professionals            24790.85   -0.10 0.921
## occupUnemployed                                         24806.67    0.02 0.984
## environ.lvl1                                               17.53    6.01 0.000
## I(environ.lvl1^2)                                       24108.72    1.11 0.267
## I(environ.lvl1^3)                                       19441.88    1.45 0.146
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.19
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.03  0.26
## occupRetired                                            -0.14  0.20
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.16
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.10  0.22
## I(environ.lvl1^2)                                       -0.01  0.02
## I(environ.lvl1^3)                                        0.00  0.02
```

## Model 4 (add random effects for non-linear terms by country)


```r
EX9RR.NCmod4<-lmer(refugees~(environ.lvl1||voting.group)+
                (environ.lvl1+
                   I(environ.lvl1^2)+
                   I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                  I(environ.lvl1^2)+
                  I(environ.lvl1^3),
                data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(EX9RR.NCmod4)
```

```
## [1] TRUE
```

```r
anova(EX9RR.NCmod3,EX9RR.NCmod4)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod4: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 + I(environ.lvl1^2) + 
## EX9RR.NCmod4:     I(environ.lvl1^3) || cntry) + age + gender + educ + resid + 
## EX9RR.NCmod4:     occup + environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX9RR.NCmod3   25 59834 60037 -29892    59784                        
## EX9RR.NCmod4   27 59827 60047 -29887    59773 10.447  2   0.005389 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(VC.EX9RR.NCmod4<-getVC(EX9RR.NCmod4))
```

```
##              grp              var1 var2     est_SD     est_SD2
## 1   voting.group       (Intercept) <NA> 0.30476354 0.092880818
## 2 voting.group.1      environ.lvl1 <NA> 0.07161810 0.005129152
## 3          cntry       (Intercept) <NA> 0.55254157 0.305302191
## 4        cntry.1      environ.lvl1 <NA> 0.08589362 0.007377714
## 5        cntry.2 I(environ.lvl1^2) <NA> 0.00000000 0.000000000
## 6        cntry.3 I(environ.lvl1^3) <NA> 0.01829238 0.000334611
## 7       Residual              <NA> <NA> 0.78881101 0.622222813
```

```r
getFE(EX9RR.NCmod4)
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
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.12       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.05       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.01       0.09
## environ.lvl1                                                0.17       0.03
## I(environ.lvl1^2)                                           0.01       0.01
## I(environ.lvl1^3)                                           0.00       0.01
##                                                               df t.value     p
## (Intercept)                                                24.17    0.37 0.714
## age                                                     24704.51   -1.70 0.090
## gender                                                  24818.74    5.16 0.000
## educ                                                    24962.75   10.50 0.000
## resid                                                   24857.01   -6.65 0.000
## occupClerical support workers                           24784.35   -0.28 0.779
## occupCraft and related trades workers                   24791.65   -1.21 0.227
## occupElementary occupations                             24790.26    0.17 0.868
## occupManagers                                           24790.24    0.45 0.656
## occupOther: Not in paid work                            24870.79    1.49 0.135
## occupPlant and machine operators, and assemblers        24792.97   -0.68 0.496
## occupProfessionals                                      24790.54    1.51 0.130
## occupRetired                                            24790.22    0.30 0.764
## occupService and sales workers                          24789.18   -0.71 0.476
## occupSkilled agricultural, forestry and fishery workers 24794.67   -0.08 0.939
## occupTechnicians and associate professionals            24785.28   -0.10 0.923
## occupUnemployed                                         24803.24    0.10 0.922
## environ.lvl1                                               16.54    6.47 0.000
## I(environ.lvl1^2)                                       21318.67    0.67 0.503
## I(environ.lvl1^3)                                          17.86    0.41 0.689
##                                                            LL    UL
## (Intercept)                                             -0.28  0.40
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.19
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.03  0.26
## occupRetired                                            -0.15  0.20
## occupService and sales workers                          -0.20  0.10
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.17  0.19
## environ.lvl1                                             0.11  0.22
## I(environ.lvl1^2)                                       -0.01  0.02
## I(environ.lvl1^3)                                       -0.01  0.02
```
## Model 5 (add random effects for non-linear terms by voting group)


```r
EX9RR.NCmod6<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3),data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(EX9RR.NCmod6)
```

```
## [1] TRUE
```

```r
#test against fixed effects only
anova(EX9RR.NCmod3,EX9RR.NCmod6)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX9RR.NCmod3   25 59834 60037 -29892    59784                        
## EX9RR.NCmod6   29 59826 60062 -29884    59768 15.646  4   0.003533 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#test against random effects by country
anova(EX9RR.NCmod4,EX9RR.NCmod6)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod4: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 + I(environ.lvl1^2) + 
## EX9RR.NCmod4:     I(environ.lvl1^3) || cntry) + age + gender + educ + resid + 
## EX9RR.NCmod4:     occup + environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod4   27 59827 60047 -29887    59773                       
## EX9RR.NCmod6   29 59826 60062 -29884    59768 5.1991  2    0.07431 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(VC.EX9RR.NCmod6<-getVC(EX9RR.NCmod6))
```

```
##              grp              var1 var2     est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.30371822 0.0922447574
## 2 voting.group.1      environ.lvl1 <NA> 0.07271645 0.0052876828
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.04248299 0.0018048046
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.00000000 0.0000000000
## 5          cntry       (Intercept) <NA> 0.55145445 0.3041020070
## 6        cntry.1      environ.lvl1 <NA> 0.08549299 0.0073090510
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.00000000 0.0000000000
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.01893597 0.0003585711
## 9       Residual              <NA> <NA> 0.78809129 0.6210878869
```

```r
getFE(EX9RR.NCmod6)
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
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.12       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.03       0.09
## occupService and sales workers                             -0.05       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.01       0.09
## environ.lvl1                                                0.17       0.03
## I(environ.lvl1^2)                                           0.01       0.01
## I(environ.lvl1^3)                                           0.00       0.01
##                                                               df t.value     p
## (Intercept)                                                24.20    0.37 0.715
## age                                                     24690.96   -1.69 0.090
## gender                                                  24810.17    5.15 0.000
## educ                                                    24962.44   10.53 0.000
## resid                                                   24841.91   -6.65 0.000
## occupClerical support workers                           24772.03   -0.27 0.785
## occupCraft and related trades workers                   24779.44   -1.21 0.226
## occupElementary occupations                             24779.12    0.16 0.871
## occupManagers                                           24774.58    0.44 0.657
## occupOther: Not in paid work                            24861.15    1.49 0.137
## occupPlant and machine operators, and assemblers        24777.75   -0.69 0.492
## occupProfessionals                                      24776.89    1.51 0.132
## occupRetired                                            24780.87    0.31 0.756
## occupService and sales workers                          24774.89   -0.71 0.476
## occupSkilled agricultural, forestry and fishery workers 24787.33   -0.10 0.924
## occupTechnicians and associate professionals            24773.57   -0.09 0.924
## occupUnemployed                                         24789.60    0.09 0.929
## environ.lvl1                                               16.58    6.47 0.000
## I(environ.lvl1^2)                                         130.87    0.68 0.498
## I(environ.lvl1^3)                                          16.82    0.45 0.656
##                                                            LL    UL
## (Intercept)                                             -0.27  0.40
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.19
## occupOther: Not in paid work                            -0.04  0.27
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.03  0.26
## occupRetired                                            -0.14  0.20
## occupService and sales workers                          -0.20  0.10
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.17  0.19
## environ.lvl1                                             0.11  0.22
## I(environ.lvl1^2)                                       -0.01  0.03
## I(environ.lvl1^3)                                       -0.01  0.02
```

\newpage

### Marginal effects for non-linearity at level-1 


```r
EX9RR.NCmod6.trends<-
  emtrends(EX9RR.NCmod6,
           specs = c("environ.lvl1"),
           var=c("environ.lvl1"),
           at=list(environ.lvl1=
                     c(mean(Rdat$environ.lvl1)-1*sd(Rdat$environ.lvl1),
                       mean(Rdat$environ.lvl1)-0*sd(Rdat$environ.lvl1),
                       mean(Rdat$environ.lvl1)+1*sd(Rdat$environ.lvl1))))

(EX9RR.NCmod6.trends.tab<-data.frame(EX9RR.NCmod6.trends))
```

```
##   environ.lvl1 environ.lvl1.trend         SE  df asymp.LCL asymp.UCL
## 1 -0.741717423          0.1621068 0.02792410 Inf 0.1073766 0.2168370
## 2  0.005748891          0.1658870 0.02562385 Inf 0.1156652 0.2161088
## 3  0.753215205          0.1815229 0.03108767 Inf 0.1205922 0.2424537
```

```r
EX9RR.NCmod6.trends.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod6.trends.tab$environ.lvl1.trend/
                   EX9RR.NCmod6.trends.tab$SE)))
EX9RR.NCmod6.trends.tab$adj.p<-
  p.adjust(EX9RR.NCmod6.trends.tab$p,method="holm")

EX9RR.NCmod6.trends.tab<-
  cbind(group=EX9RR.NCmod6.trends.tab[,1],
        round(EX9RR.NCmod6.trends.tab[,c(2,3)],2),
        round(EX9RR.NCmod6.trends.tab[,c(7,8)],4),
        round(EX9RR.NCmod6.trends.tab[,c(5,6)],2))
EX9RR.NCmod6.trends.tab
```

```
##          group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -0.741717423               0.16 0.03 0     0      0.11      0.22
## 2  0.005748891               0.17 0.03 0     0      0.12      0.22
## 3  0.753215205               0.18 0.03 0     0      0.12      0.24
```

```r
write.csv2(EX9RR.NCmod6.trends.tab,"EX9RR.NCmod6.trends.tab.csv")
```

\newpage

## Model 6 (include voting group main effects at level-2)


```r
EX9RR.NCmod7<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3)+
                     all.parties.lvl2,data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(EX9RR.NCmod7)
```

```
## [1] TRUE
```

```r
anova(EX9RR.NCmod6,EX9RR.NCmod7)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod7: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod7:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX9RR.NCmod6   29 59826 60062 -29884    59768                         
## EX9RR.NCmod7   36 59676 59968 -29802    59604 164.08  7  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod7)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.51       0.17
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.02       0.08
## occupCraft and related trades workers                      -0.09       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.11       0.08
## occupPlant and machine operators, and assemblers           -0.05       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.02       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.17       0.03
## I(environ.lvl1^2)                                           0.01       0.01
## I(environ.lvl1^3)                                           0.00       0.01
## all.parties.lvl2Did not vote                                0.46       0.07
## all.parties.lvl2Don't know                                  0.43       0.07
## all.parties.lvl2NE age                                      0.73       0.08
## all.parties.lvl2NE citizen                                  0.83       0.08
## all.parties.lvl2NE other                                    0.82       0.10
## all.parties.lvl2Other party                                 0.57       0.05
## all.parties.lvl2Pro-environment party                       0.95       0.07
##                                                               df t.value     p
## (Intercept)                                                28.96   -3.09 0.004
## age                                                     24929.16   -1.36 0.173
## gender                                                  24833.76    5.11 0.000
## educ                                                    24935.94   10.47 0.000
## resid                                                   24889.34   -6.56 0.000
## occupClerical support workers                           24816.86   -0.29 0.768
## occupCraft and related trades workers                   24824.10   -1.22 0.223
## occupElementary occupations                             24821.00    0.13 0.897
## occupManagers                                           24818.62    0.42 0.677
## occupOther: Not in paid work                            24848.03    1.36 0.173
## occupPlant and machine operators, and assemblers        24824.06   -0.71 0.481
## occupProfessionals                                      24820.30    1.48 0.138
## occupRetired                                            24821.11    0.28 0.782
## occupService and sales workers                          24818.04   -0.72 0.470
## occupSkilled agricultural, forestry and fishery workers 24832.89   -0.12 0.906
## occupTechnicians and associate professionals            24817.68   -0.11 0.910
## occupUnemployed                                         24817.48    0.01 0.994
## environ.lvl1                                               16.58    6.46 0.000
## I(environ.lvl1^2)                                         130.01    0.68 0.496
## I(environ.lvl1^3)                                          17.16    0.46 0.648
## all.parties.lvl2Did not vote                              162.46    6.65 0.000
## all.parties.lvl2Don't know                                209.08    5.73 0.000
## all.parties.lvl2NE age                                    215.17    9.75 0.000
## all.parties.lvl2NE citizen                                216.68   10.73 0.000
## all.parties.lvl2NE other                                  599.83    7.93 0.000
## all.parties.lvl2Other party                               195.02   10.79 0.000
## all.parties.lvl2Pro-environment party                     206.16   14.12 0.000
##                                                            LL    UL
## (Intercept)                                             -0.85 -0.17
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.24  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.18
## occupOther: Not in paid work                            -0.05  0.26
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.04  0.26
## occupRetired                                            -0.15  0.20
## occupService and sales workers                          -0.20  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.11  0.22
## I(environ.lvl1^2)                                       -0.01  0.02
## I(environ.lvl1^3)                                       -0.01  0.02
## all.parties.lvl2Did not vote                             0.33  0.60
## all.parties.lvl2Don't know                               0.28  0.58
## all.parties.lvl2NE age                                   0.59  0.88
## all.parties.lvl2NE citizen                               0.67  0.98
## all.parties.lvl2NE other                                 0.61  1.02
## all.parties.lvl2Other party                              0.46  0.67
## all.parties.lvl2Pro-environment party                    0.81  1.08
```

```r
getVC(EX9RR.NCmod7)
```

```
##              grp              var1 var2     est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.18818707 0.0354143733
## 2 voting.group.1      environ.lvl1 <NA> 0.07285455 0.0053077857
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.04153995 0.0017255678
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.00000000 0.0000000000
## 5          cntry       (Intercept) <NA> 0.54006403 0.2916691529
## 6        cntry.1      environ.lvl1 <NA> 0.08568104 0.0073412399
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.00000000 0.0000000000
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.01869883 0.0003496462
## 9       Residual              <NA> <NA> 0.78807924 0.6210688897
```

## Model 7 (voting group interactions with linear term)


```r
EX9RR.NCmod8<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3)+
                     all.parties.lvl2+
                     all.parties.lvl2:environ.lvl1,
                   data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(EX9RR.NCmod8)
```

```
## [1] TRUE
```

```r
anova(EX9RR.NCmod7,EX9RR.NCmod8)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod7: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod7:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
##              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX9RR.NCmod7   36 59676 59968 -29802    59604                        
## EX9RR.NCmod8   43 59666 60016 -29790    59580 23.726  7   0.001273 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod8)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.51       0.17
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.02       0.08
## occupCraft and related trades workers                      -0.10       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.11       0.08
## occupPlant and machine operators, and assemblers           -0.06       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.02       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.08       0.04
## I(environ.lvl1^2)                                           0.01       0.01
## I(environ.lvl1^3)                                           0.00       0.01
## all.parties.lvl2Did not vote                                0.46       0.07
## all.parties.lvl2Don't know                                  0.43       0.07
## all.parties.lvl2NE age                                      0.73       0.08
## all.parties.lvl2NE citizen                                  0.83       0.08
## all.parties.lvl2NE other                                    0.82       0.10
## all.parties.lvl2Other party                                 0.56       0.05
## all.parties.lvl2Pro-environment party                       0.94       0.07
## environ.lvl1:all.parties.lvl2Did not vote                   0.03       0.04
## environ.lvl1:all.parties.lvl2Don't know                     0.12       0.05
## environ.lvl1:all.parties.lvl2NE age                         0.13       0.05
## environ.lvl1:all.parties.lvl2NE citizen                     0.02       0.05
## environ.lvl1:all.parties.lvl2NE other                       0.01       0.10
## environ.lvl1:all.parties.lvl2Other party                    0.11       0.03
## environ.lvl1:all.parties.lvl2Pro-environment party          0.09       0.04
##                                                               df t.value     p
## (Intercept)                                                28.97   -3.07 0.005
## age                                                     24924.77   -1.31 0.189
## gender                                                  24836.30    5.09 0.000
## educ                                                    24936.95   10.47 0.000
## resid                                                   24896.16   -6.54 0.000
## occupClerical support workers                           24816.80   -0.31 0.757
## occupCraft and related trades workers                   24822.99   -1.24 0.214
## occupElementary occupations                             24818.70    0.11 0.909
## occupManagers                                           24813.34    0.40 0.691
## occupOther: Not in paid work                            24847.16    1.36 0.175
## occupPlant and machine operators, and assemblers        24821.01   -0.73 0.468
## occupProfessionals                                      24818.18    1.46 0.144
## occupRetired                                            24820.93    0.24 0.809
## occupService and sales workers                          24816.07   -0.74 0.459
## occupSkilled agricultural, forestry and fishery workers 24834.90   -0.12 0.901
## occupTechnicians and associate professionals            24816.23   -0.13 0.896
## occupUnemployed                                         24820.15   -0.02 0.986
## environ.lvl1                                               57.89    2.11 0.039
## I(environ.lvl1^2)                                         134.79    0.82 0.414
## I(environ.lvl1^3)                                          16.99    0.62 0.541
## all.parties.lvl2Did not vote                              162.44    6.65 0.000
## all.parties.lvl2Don't know                                209.15    5.70 0.000
## all.parties.lvl2NE age                                    215.09    9.72 0.000
## all.parties.lvl2NE citizen                                216.88   10.72 0.000
## all.parties.lvl2NE other                                  603.68    7.94 0.000
## all.parties.lvl2Other party                               195.08   10.73 0.000
## all.parties.lvl2Pro-environment party                     206.17   14.10 0.000
## environ.lvl1:all.parties.lvl2Did not vote                  87.41    0.85 0.397
## environ.lvl1:all.parties.lvl2Don't know                   277.09    2.31 0.022
## environ.lvl1:all.parties.lvl2NE age                       214.97    2.62 0.009
## environ.lvl1:all.parties.lvl2NE citizen                   162.39    0.43 0.669
## environ.lvl1:all.parties.lvl2NE other                     698.74    0.07 0.948
## environ.lvl1:all.parties.lvl2Other party                  120.24    3.71 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        204.02    2.15 0.033
##                                                            LL    UL
## (Intercept)                                             -0.85 -0.17
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.17  0.13
## occupCraft and related trades workers                   -0.25  0.06
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.18
## occupOther: Not in paid work                            -0.05  0.26
## occupPlant and machine operators, and assemblers        -0.21  0.10
## occupProfessionals                                      -0.04  0.26
## occupRetired                                            -0.15  0.19
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                             0.00  0.15
## I(environ.lvl1^2)                                       -0.01  0.03
## I(environ.lvl1^3)                                       -0.01  0.02
## all.parties.lvl2Did not vote                             0.33  0.60
## all.parties.lvl2Don't know                               0.28  0.57
## all.parties.lvl2NE age                                   0.58  0.88
## all.parties.lvl2NE citizen                               0.67  0.98
## all.parties.lvl2NE other                                 0.62  1.02
## all.parties.lvl2Other party                              0.46  0.67
## all.parties.lvl2Pro-environment party                    0.81  1.08
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.10
## environ.lvl1:all.parties.lvl2Don't know                  0.02  0.22
## environ.lvl1:all.parties.lvl2NE age                      0.03  0.22
## environ.lvl1:all.parties.lvl2NE citizen                 -0.07  0.11
## environ.lvl1:all.parties.lvl2NE other                   -0.19  0.20
## environ.lvl1:all.parties.lvl2Other party                 0.05  0.17
## environ.lvl1:all.parties.lvl2Pro-environment party       0.01  0.18
```

```r
getVC(EX9RR.NCmod8)
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.188126526 3.539159e-02
## 2 voting.group.1      environ.lvl1 <NA> 0.056955407 3.243918e-03
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.045190637 2.042194e-03
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.007411094 5.492431e-05
## 5          cntry       (Intercept) <NA> 0.539783385 2.913661e-01
## 6        cntry.1      environ.lvl1 <NA> 0.089152995 7.948257e-03
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.000000000 0.000000e+00
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.018554349 3.442639e-04
## 9       Residual              <NA> <NA> 0.787883565 6.207605e-01
```

\newpage

## Model 8 (interactions with non-linear terms)


```r
EX9RR.NCmod10<-lmer(refugees~
                      (environ.lvl1+
                         I(environ.lvl1^2)+
                         I(environ.lvl1^3)||voting.group)+
                      (environ.lvl1+
                         I(environ.lvl1^2)+
                         I(environ.lvl1^3)||cntry)+
                      age+gender+educ+resid+occup+
                      environ.lvl1+
                      I(environ.lvl1^2)+
                      I(environ.lvl1^3)+
                      all.parties.lvl2+
                      all.parties.lvl2:environ.lvl1+
                      all.parties.lvl2:I(environ.lvl1^2)+
                      all.parties.lvl2:I(environ.lvl1^3),
                    data=Rdat,REML=F,
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
#test against linear interaction by voting groups
anova(EX9RR.NCmod8,EX9RR.NCmod10)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
## EX9RR.NCmod10: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod10:     all.parties.lvl2:environ.lvl1 + all.parties.lvl2:I(environ.lvl1^2) + 
## EX9RR.NCmod10:     all.parties.lvl2:I(environ.lvl1^3)
##               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod8    43 59666 60016 -29790    59580                       
## EX9RR.NCmod10   57 59670 60133 -29778    59556 24.571 14    0.03904 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod10)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.49       0.17
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.03       0.08
## occupCraft and related trades workers                      -0.10       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.11       0.08
## occupPlant and machine operators, and assemblers           -0.06       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.02       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.08       0.05
## I(environ.lvl1^2)                                          -0.04       0.03
## I(environ.lvl1^3)                                          -0.01       0.02
## all.parties.lvl2Did not vote                                0.46       0.07
## all.parties.lvl2Don't know                                  0.41       0.08
## all.parties.lvl2NE age                                      0.69       0.08
## all.parties.lvl2NE citizen                                  0.82       0.08
## all.parties.lvl2NE other                                    0.77       0.11
## all.parties.lvl2Other party                                 0.53       0.05
## all.parties.lvl2Pro-environment party                       0.94       0.07
## environ.lvl1:all.parties.lvl2Did not vote                   0.01       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.15       0.07
## environ.lvl1:all.parties.lvl2NE age                         0.03       0.07
## environ.lvl1:all.parties.lvl2NE citizen                     0.02       0.07
## environ.lvl1:all.parties.lvl2NE other                      -0.16       0.15
## environ.lvl1:all.parties.lvl2Other party                    0.11       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.12       0.06
## I(environ.lvl1^2):all.parties.lvl2Did not vote              0.01       0.03
## I(environ.lvl1^2):all.parties.lvl2Don't know                0.03       0.05
## I(environ.lvl1^2):all.parties.lvl2NE age                    0.11       0.05
## I(environ.lvl1^2):all.parties.lvl2NE citizen                0.01       0.05
## I(environ.lvl1^2):all.parties.lvl2NE other                  0.14       0.10
## I(environ.lvl1^2):all.parties.lvl2Other party               0.07       0.03
## I(environ.lvl1^2):all.parties.lvl2Pro-environment party    -0.02       0.06
## I(environ.lvl1^3):all.parties.lvl2Did not vote              0.02       0.02
## I(environ.lvl1^3):all.parties.lvl2Don't know               -0.01       0.03
## I(environ.lvl1^3):all.parties.lvl2NE age                    0.08       0.03
## I(environ.lvl1^3):all.parties.lvl2NE citizen                0.01       0.03
## I(environ.lvl1^3):all.parties.lvl2NE other                  0.11       0.07
## I(environ.lvl1^3):all.parties.lvl2Other party               0.01       0.02
## I(environ.lvl1^3):all.parties.lvl2Pro-environment party    -0.03       0.03
##                                                               df t.value     p
## (Intercept)                                                29.24   -2.94 0.006
## age                                                     24929.36   -1.35 0.177
## gender                                                  24839.79    5.10 0.000
## educ                                                    24931.47   10.38 0.000
## resid                                                   24901.00   -6.56 0.000
## occupClerical support workers                           24820.82   -0.32 0.745
## occupCraft and related trades workers                   24827.47   -1.26 0.208
## occupElementary occupations                             24822.32    0.09 0.930
## occupManagers                                           24816.61    0.38 0.703
## occupOther: Not in paid work                            24849.83    1.34 0.181
## occupPlant and machine operators, and assemblers        24825.65   -0.75 0.451
## occupProfessionals                                      24822.28    1.46 0.146
## occupRetired                                            24825.38    0.25 0.799
## occupService and sales workers                          24820.79   -0.76 0.449
## occupSkilled agricultural, forestry and fishery workers 24838.75   -0.14 0.892
## occupTechnicians and associate professionals            24820.17   -0.14 0.887
## occupUnemployed                                         24815.30    0.02 0.987
## environ.lvl1                                              117.81    1.85 0.067
## I(environ.lvl1^2)                                         173.47   -1.42 0.158
## I(environ.lvl1^3)                                         211.25   -0.33 0.745
## all.parties.lvl2Did not vote                              175.63    6.43 0.000
## all.parties.lvl2Don't know                                245.64    5.26 0.000
## all.parties.lvl2NE age                                    247.43    8.88 0.000
## all.parties.lvl2NE citizen                                252.85   10.30 0.000
## all.parties.lvl2NE other                                  929.83    6.86 0.000
## all.parties.lvl2Other party                               218.08    9.89 0.000
## all.parties.lvl2Pro-environment party                     248.07   13.46 0.000
## environ.lvl1:all.parties.lvl2Did not vote                 250.54    0.25 0.802
## environ.lvl1:all.parties.lvl2Don't know                   955.46    2.01 0.045
## environ.lvl1:all.parties.lvl2NE age                       725.82    0.52 0.603
## environ.lvl1:all.parties.lvl2NE citizen                   620.06    0.26 0.795
## environ.lvl1:all.parties.lvl2NE other                    8884.13   -1.07 0.286
## environ.lvl1:all.parties.lvl2Other party                  344.92    2.70 0.007
## environ.lvl1:all.parties.lvl2Pro-environment party        518.83    2.13 0.034
## I(environ.lvl1^2):all.parties.lvl2Did not vote             98.85    0.36 0.716
## I(environ.lvl1^2):all.parties.lvl2Don't know              458.85    0.64 0.522
## I(environ.lvl1^2):all.parties.lvl2NE age                  386.85    2.05 0.041
## I(environ.lvl1^2):all.parties.lvl2NE citizen              193.29    0.31 0.757
## I(environ.lvl1^2):all.parties.lvl2NE other                445.66    1.39 0.165
## I(environ.lvl1^2):all.parties.lvl2Other party             152.78    2.31 0.022
## I(environ.lvl1^2):all.parties.lvl2Pro-environment party   950.74   -0.37 0.712
## I(environ.lvl1^3):all.parties.lvl2Did not vote            100.19    0.73 0.468
## I(environ.lvl1^3):all.parties.lvl2Don't know              339.82   -0.29 0.770
## I(environ.lvl1^3):all.parties.lvl2NE age                  248.08    2.37 0.019
## I(environ.lvl1^3):all.parties.lvl2NE citizen              362.29    0.17 0.863
## I(environ.lvl1^3):all.parties.lvl2NE other                464.03    1.69 0.092
## I(environ.lvl1^3):all.parties.lvl2Other party             215.53    0.62 0.535
## I(environ.lvl1^3):all.parties.lvl2Pro-environment party   393.97   -1.02 0.310
##                                                            LL    UL
## (Intercept)                                             -0.83 -0.15
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.18  0.13
## occupCraft and related trades workers                   -0.25  0.05
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.18
## occupOther: Not in paid work                            -0.05  0.26
## occupPlant and machine operators, and assemblers        -0.21  0.09
## occupProfessionals                                      -0.04  0.26
## occupRetired                                            -0.15  0.19
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                            -0.01  0.17
## I(environ.lvl1^2)                                       -0.10  0.02
## I(environ.lvl1^3)                                       -0.04  0.03
## all.parties.lvl2Did not vote                             0.32  0.60
## all.parties.lvl2Don't know                               0.26  0.56
## all.parties.lvl2NE age                                   0.54  0.85
## all.parties.lvl2NE citizen                               0.67  0.98
## all.parties.lvl2NE other                                 0.55  0.99
## all.parties.lvl2Other party                              0.43  0.64
## all.parties.lvl2Pro-environment party                    0.81  1.08
## environ.lvl1:all.parties.lvl2Did not vote               -0.09  0.11
## environ.lvl1:all.parties.lvl2Don't know                  0.00  0.29
## environ.lvl1:all.parties.lvl2NE age                     -0.10  0.17
## environ.lvl1:all.parties.lvl2NE citizen                 -0.12  0.15
## environ.lvl1:all.parties.lvl2NE other                   -0.45  0.13
## environ.lvl1:all.parties.lvl2Other party                 0.03  0.20
## environ.lvl1:all.parties.lvl2Pro-environment party       0.01  0.24
## I(environ.lvl1^2):all.parties.lvl2Did not vote          -0.06  0.08
## I(environ.lvl1^2):all.parties.lvl2Don't know            -0.07  0.14
## I(environ.lvl1^2):all.parties.lvl2NE age                 0.00  0.22
## I(environ.lvl1^2):all.parties.lvl2NE citizen            -0.08  0.11
## I(environ.lvl1^2):all.parties.lvl2NE other              -0.06  0.34
## I(environ.lvl1^2):all.parties.lvl2Other party            0.01  0.13
## I(environ.lvl1^2):all.parties.lvl2Pro-environment party -0.14  0.09
## I(environ.lvl1^3):all.parties.lvl2Did not vote          -0.03  0.06
## I(environ.lvl1^3):all.parties.lvl2Don't know            -0.07  0.05
## I(environ.lvl1^3):all.parties.lvl2NE age                 0.01  0.14
## I(environ.lvl1^3):all.parties.lvl2NE citizen            -0.05  0.06
## I(environ.lvl1^3):all.parties.lvl2NE other              -0.02  0.25
## I(environ.lvl1^3):all.parties.lvl2Other party           -0.03  0.05
## I(environ.lvl1^3):all.parties.lvl2Pro-environment party -0.09  0.03
```

```r
getVC(EX9RR.NCmod10)
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.187898363 3.530579e-02
## 2 voting.group.1      environ.lvl1 <NA> 0.058762365 3.453016e-03
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.035635879 1.269916e-03
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.005707735 3.257824e-05
## 5          cntry       (Intercept) <NA> 0.539884711 2.914755e-01
## 6        cntry.1      environ.lvl1 <NA> 0.089733237 8.052054e-03
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.000000000 0.000000e+00
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.016330226 2.666763e-04
## 9       Residual              <NA> <NA> 0.787748002 6.205469e-01
```

\newpage

### Marginal trends for each voting group

#### Linear coefficients


```r
EX9RR.NCmod10.linear<-
  emtrends(EX9RR.NCmod10,specs = c("all.parties.lvl2"),var=c("environ.lvl1"))
(EX9RR.NCmod10.linear.tab<-data.frame(EX9RR.NCmod10.linear))
```

```
##         all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL
## 1 Anti-immigration party         0.08306674 0.04524963 Inf -0.005620906
## 2           Did not vote         0.09570158 0.03843323 Inf  0.020373837
## 3             Don't know         0.22924118 0.06529439 Inf  0.101266533
## 4                 NE age         0.11985396 0.05942261 Inf  0.003387789
## 5             NE citizen         0.10137912 0.06179487 Inf -0.019736597
## 6               NE other        -0.07207832 0.14404533 Inf -0.354401982
## 7            Other party         0.19868954 0.02869801 Inf  0.142442481
## 8  Pro-environment party         0.20734184 0.04945930 Inf  0.110403402
##   asymp.UCL
## 1 0.1717544
## 2 0.1710293
## 3 0.3572158
## 4 0.2363201
## 5 0.2224948
## 6 0.2102453
## 7 0.2549366
## 8 0.3042803
```

```r
EX9RR.NCmod10.linear.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.linear.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.linear.tab$SE)))
EX9RR.NCmod10.linear.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.linear.tab$p,method="holm")

EX9RR.NCmod10.linear.tab<-
  cbind(group=EX9RR.NCmod10.linear.tab[,1],
        round(EX9RR.NCmod10.linear.tab[,c(2,3)],2),
        round(EX9RR.NCmod10.linear.tab[,c(7,8)],4),
        round(EX9RR.NCmod10.linear.tab[,c(5,6)],2))
EX9RR.NCmod10.linear.tab
```

```
##                    group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 1 Anti-immigration party               0.08 0.05 0.0664 0.1992     -0.01
## 2           Did not vote               0.10 0.04 0.0128 0.0639      0.02
## 3             Don't know               0.23 0.07 0.0004 0.0027      0.10
## 4                 NE age               0.12 0.06 0.0437 0.1748      0.00
## 5             NE citizen               0.10 0.06 0.1009 0.2018     -0.02
## 6               NE other              -0.07 0.14 0.6168 0.6168     -0.35
## 7            Other party               0.20 0.03 0.0000 0.0000      0.14
## 8  Pro-environment party               0.21 0.05 0.0000 0.0002      0.11
##   asymp.UCL
## 1      0.17
## 2      0.17
## 3      0.36
## 4      0.24
## 5      0.22
## 6      0.21
## 7      0.25
## 8      0.30
```

```r
write.csv2(EX9RR.NCmod10.linear.tab,"EX9RR.NCmod10.linear.tab.csv")

#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.linear, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1           Other party - Anti-immigration party 0.115622804 0.04233000 Inf
## 2 Pro-environment party - Anti-immigration party 0.124275104 0.05831638 Inf
## 3            Pro-environment party - Other party 0.008652301 0.04677753 Inf
##    z.ratio    p.value
## 1 2.731463 0.00630539
## 2 2.131050 0.03308504
## 3 0.184967 0.85325491
```

### Quadratic and cubic coefficients

* Non-linear variables must be manually recoded to new variables for extracting marginal effects


```r
#recoding
Rdat$environ.lvl1.sq<-Rdat$environ.lvl1^2
Rdat$environ.lvl1.cu<-Rdat$environ.lvl1^3

#refitting the model
EX9RR.NCmod10.ma<-lmer(refugees~
                         (environ.lvl1+
                            environ.lvl1.sq+
                            environ.lvl1.cu||voting.group)+
                         (environ.lvl1+
                            environ.lvl1.sq+
                            environ.lvl1.cu||cntry)+
                         age+gender+educ+resid+occup+
                         environ.lvl1+
                         environ.lvl1.sq+
                         environ.lvl1.cu+
                         all.parties.lvl2+
                         all.parties.lvl2:environ.lvl1+
                         all.parties.lvl2:environ.lvl1.sq+
                       all.parties.lvl2:environ.lvl1.cu,
                       data=Rdat,REML=F,
                       control=lmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(EX9RR.NCmod10.ma)
```

```
## [1] TRUE
```

```r
#check if identical
anova(EX9RR.NCmod10,EX9RR.NCmod10.ma)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod10: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod10:     all.parties.lvl2:environ.lvl1 + all.parties.lvl2:I(environ.lvl1^2) + 
## EX9RR.NCmod10:     all.parties.lvl2:I(environ.lvl1^3)
## EX9RR.NCmod10.ma: refugees ~ (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     voting.group) + (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10.ma:     environ.lvl1.sq + environ.lvl1.cu + all.parties.lvl2 + all.parties.lvl2:environ.lvl1 + 
## EX9RR.NCmod10.ma:     all.parties.lvl2:environ.lvl1.sq + all.parties.lvl2:environ.lvl1.cu
##                  npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)
## EX9RR.NCmod10      57 59670 60133 -29778    59556                    
## EX9RR.NCmod10.ma   57 59670 60133 -29778    59556     0  0          1
```

```r
#
anova(EX9RR.NCmod8,EX9RR.NCmod10.ma)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
## EX9RR.NCmod10.ma: refugees ~ (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     voting.group) + (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10.ma:     environ.lvl1.sq + environ.lvl1.cu + all.parties.lvl2 + all.parties.lvl2:environ.lvl1 + 
## EX9RR.NCmod10.ma:     all.parties.lvl2:environ.lvl1.sq + all.parties.lvl2:environ.lvl1.cu
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod8       43 59666 60016 -29790    59580                       
## EX9RR.NCmod10.ma   57 59670 60133 -29778    59556 24.571 14    0.03904 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod10.ma)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.49       0.17
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.07       0.01
## occupClerical support workers                              -0.03       0.08
## occupCraft and related trades workers                      -0.10       0.08
## occupElementary occupations                                 0.01       0.08
## occupManagers                                               0.03       0.08
## occupOther: Not in paid work                                0.11       0.08
## occupPlant and machine operators, and assemblers           -0.06       0.08
## occupProfessionals                                          0.11       0.08
## occupRetired                                                0.02       0.09
## occupService and sales workers                             -0.06       0.08
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.08
## occupTechnicians and associate professionals               -0.01       0.08
## occupUnemployed                                             0.00       0.09
## environ.lvl1                                                0.08       0.05
## environ.lvl1.sq                                            -0.04       0.03
## environ.lvl1.cu                                            -0.01       0.02
## all.parties.lvl2Did not vote                                0.46       0.07
## all.parties.lvl2Don't know                                  0.41       0.08
## all.parties.lvl2NE age                                      0.69       0.08
## all.parties.lvl2NE citizen                                  0.82       0.08
## all.parties.lvl2NE other                                    0.77       0.11
## all.parties.lvl2Other party                                 0.53       0.05
## all.parties.lvl2Pro-environment party                       0.94       0.07
## environ.lvl1:all.parties.lvl2Did not vote                   0.01       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.15       0.07
## environ.lvl1:all.parties.lvl2NE age                         0.03       0.07
## environ.lvl1:all.parties.lvl2NE citizen                     0.02       0.07
## environ.lvl1:all.parties.lvl2NE other                      -0.16       0.15
## environ.lvl1:all.parties.lvl2Other party                    0.11       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.12       0.06
## environ.lvl1.sq:all.parties.lvl2Did not vote                0.01       0.03
## environ.lvl1.sq:all.parties.lvl2Don't know                  0.03       0.05
## environ.lvl1.sq:all.parties.lvl2NE age                      0.11       0.05
## environ.lvl1.sq:all.parties.lvl2NE citizen                  0.01       0.05
## environ.lvl1.sq:all.parties.lvl2NE other                    0.14       0.10
## environ.lvl1.sq:all.parties.lvl2Other party                 0.07       0.03
## environ.lvl1.sq:all.parties.lvl2Pro-environment party      -0.02       0.06
## environ.lvl1.cu:all.parties.lvl2Did not vote                0.02       0.02
## environ.lvl1.cu:all.parties.lvl2Don't know                 -0.01       0.03
## environ.lvl1.cu:all.parties.lvl2NE age                      0.08       0.03
## environ.lvl1.cu:all.parties.lvl2NE citizen                  0.01       0.03
## environ.lvl1.cu:all.parties.lvl2NE other                    0.11       0.07
## environ.lvl1.cu:all.parties.lvl2Other party                 0.01       0.02
## environ.lvl1.cu:all.parties.lvl2Pro-environment party      -0.03       0.03
##                                                               df t.value     p
## (Intercept)                                                29.24   -2.94 0.006
## age                                                     24929.36   -1.35 0.177
## gender                                                  24839.79    5.10 0.000
## educ                                                    24931.47   10.38 0.000
## resid                                                   24901.00   -6.56 0.000
## occupClerical support workers                           24820.82   -0.32 0.745
## occupCraft and related trades workers                   24827.47   -1.26 0.208
## occupElementary occupations                             24822.32    0.09 0.930
## occupManagers                                           24816.61    0.38 0.703
## occupOther: Not in paid work                            24849.83    1.34 0.181
## occupPlant and machine operators, and assemblers        24825.65   -0.75 0.451
## occupProfessionals                                      24822.28    1.46 0.146
## occupRetired                                            24825.38    0.25 0.799
## occupService and sales workers                          24820.79   -0.76 0.449
## occupSkilled agricultural, forestry and fishery workers 24838.75   -0.14 0.892
## occupTechnicians and associate professionals            24820.17   -0.14 0.887
## occupUnemployed                                         24815.30    0.02 0.987
## environ.lvl1                                              117.81    1.85 0.067
## environ.lvl1.sq                                           173.47   -1.42 0.158
## environ.lvl1.cu                                           211.25   -0.33 0.745
## all.parties.lvl2Did not vote                              175.63    6.43 0.000
## all.parties.lvl2Don't know                                245.64    5.26 0.000
## all.parties.lvl2NE age                                    247.43    8.88 0.000
## all.parties.lvl2NE citizen                                252.85   10.30 0.000
## all.parties.lvl2NE other                                  929.83    6.86 0.000
## all.parties.lvl2Other party                               218.08    9.89 0.000
## all.parties.lvl2Pro-environment party                     248.07   13.46 0.000
## environ.lvl1:all.parties.lvl2Did not vote                 250.54    0.25 0.802
## environ.lvl1:all.parties.lvl2Don't know                   955.46    2.01 0.045
## environ.lvl1:all.parties.lvl2NE age                       725.82    0.52 0.603
## environ.lvl1:all.parties.lvl2NE citizen                   620.06    0.26 0.795
## environ.lvl1:all.parties.lvl2NE other                    8884.13   -1.07 0.286
## environ.lvl1:all.parties.lvl2Other party                  344.92    2.70 0.007
## environ.lvl1:all.parties.lvl2Pro-environment party        518.83    2.13 0.034
## environ.lvl1.sq:all.parties.lvl2Did not vote               98.85    0.36 0.716
## environ.lvl1.sq:all.parties.lvl2Don't know                458.85    0.64 0.522
## environ.lvl1.sq:all.parties.lvl2NE age                    386.85    2.05 0.041
## environ.lvl1.sq:all.parties.lvl2NE citizen                193.29    0.31 0.757
## environ.lvl1.sq:all.parties.lvl2NE other                  445.66    1.39 0.165
## environ.lvl1.sq:all.parties.lvl2Other party               152.78    2.31 0.022
## environ.lvl1.sq:all.parties.lvl2Pro-environment party     950.74   -0.37 0.712
## environ.lvl1.cu:all.parties.lvl2Did not vote              100.19    0.73 0.468
## environ.lvl1.cu:all.parties.lvl2Don't know                339.82   -0.29 0.770
## environ.lvl1.cu:all.parties.lvl2NE age                    248.08    2.37 0.019
## environ.lvl1.cu:all.parties.lvl2NE citizen                362.29    0.17 0.863
## environ.lvl1.cu:all.parties.lvl2NE other                  464.03    1.69 0.092
## environ.lvl1.cu:all.parties.lvl2Other party               215.53    0.62 0.535
## environ.lvl1.cu:all.parties.lvl2Pro-environment party     393.97   -1.02 0.310
##                                                            LL    UL
## (Intercept)                                             -0.83 -0.15
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.09 -0.05
## occupClerical support workers                           -0.18  0.13
## occupCraft and related trades workers                   -0.25  0.05
## occupElementary occupations                             -0.14  0.16
## occupManagers                                           -0.12  0.18
## occupOther: Not in paid work                            -0.05  0.26
## occupPlant and machine operators, and assemblers        -0.21  0.09
## occupProfessionals                                      -0.04  0.26
## occupRetired                                            -0.15  0.19
## occupService and sales workers                          -0.21  0.09
## occupSkilled agricultural, forestry and fishery workers -0.17  0.15
## occupTechnicians and associate professionals            -0.16  0.14
## occupUnemployed                                         -0.18  0.18
## environ.lvl1                                            -0.01  0.17
## environ.lvl1.sq                                         -0.10  0.02
## environ.lvl1.cu                                         -0.04  0.03
## all.parties.lvl2Did not vote                             0.32  0.60
## all.parties.lvl2Don't know                               0.26  0.56
## all.parties.lvl2NE age                                   0.54  0.85
## all.parties.lvl2NE citizen                               0.67  0.98
## all.parties.lvl2NE other                                 0.55  0.99
## all.parties.lvl2Other party                              0.43  0.64
## all.parties.lvl2Pro-environment party                    0.81  1.08
## environ.lvl1:all.parties.lvl2Did not vote               -0.09  0.11
## environ.lvl1:all.parties.lvl2Don't know                  0.00  0.29
## environ.lvl1:all.parties.lvl2NE age                     -0.10  0.17
## environ.lvl1:all.parties.lvl2NE citizen                 -0.12  0.15
## environ.lvl1:all.parties.lvl2NE other                   -0.45  0.13
## environ.lvl1:all.parties.lvl2Other party                 0.03  0.20
## environ.lvl1:all.parties.lvl2Pro-environment party       0.01  0.24
## environ.lvl1.sq:all.parties.lvl2Did not vote            -0.06  0.08
## environ.lvl1.sq:all.parties.lvl2Don't know              -0.07  0.14
## environ.lvl1.sq:all.parties.lvl2NE age                   0.00  0.22
## environ.lvl1.sq:all.parties.lvl2NE citizen              -0.08  0.11
## environ.lvl1.sq:all.parties.lvl2NE other                -0.06  0.34
## environ.lvl1.sq:all.parties.lvl2Other party              0.01  0.13
## environ.lvl1.sq:all.parties.lvl2Pro-environment party   -0.14  0.09
## environ.lvl1.cu:all.parties.lvl2Did not vote            -0.03  0.06
## environ.lvl1.cu:all.parties.lvl2Don't know              -0.07  0.05
## environ.lvl1.cu:all.parties.lvl2NE age                   0.01  0.14
## environ.lvl1.cu:all.parties.lvl2NE citizen              -0.05  0.06
## environ.lvl1.cu:all.parties.lvl2NE other                -0.02  0.25
## environ.lvl1.cu:all.parties.lvl2Other party             -0.03  0.05
## environ.lvl1.cu:all.parties.lvl2Pro-environment party   -0.09  0.03
```

```r
getVC(EX9RR.NCmod10.ma)
```

```
##              grp            var1 var2      est_SD      est_SD2
## 1   voting.group     (Intercept) <NA> 0.187898363 3.530579e-02
## 2 voting.group.1    environ.lvl1 <NA> 0.058762365 3.453016e-03
## 3 voting.group.2 environ.lvl1.sq <NA> 0.035635879 1.269916e-03
## 4 voting.group.3 environ.lvl1.cu <NA> 0.005707735 3.257824e-05
## 5          cntry     (Intercept) <NA> 0.539884711 2.914755e-01
## 6        cntry.1    environ.lvl1 <NA> 0.089733237 8.052054e-03
## 7        cntry.2 environ.lvl1.sq <NA> 0.000000000 0.000000e+00
## 8        cntry.3 environ.lvl1.cu <NA> 0.016330226 2.666763e-04
## 9       Residual            <NA> <NA> 0.787748002 6.205469e-01
```

#### Marginal quadratic effects


```r
EX9RR.NCmod10.quadratic<-
  emtrends(EX9RR.NCmod10.ma,specs = c("all.parties.lvl2"),var=c("environ.lvl1.sq"))
(EX9RR.NCmod10.quadratic.tab<-data.frame(EX9RR.NCmod10.quadratic))
```

```
##         all.parties.lvl2 environ.lvl1.sq.trend         SE  df    asymp.LCL
## 1 Anti-immigration party           -0.04004512 0.02826913 Inf -0.095451592
## 2           Did not vote           -0.02732427 0.02059382 Inf -0.067687404
## 3             Don't know           -0.00551286 0.04603177 Inf -0.095733467
## 4                 NE age            0.07015228 0.04568164 Inf -0.019382101
## 5             NE citizen           -0.02565880 0.03690659 Inf -0.097994389
## 6               NE other            0.09932585 0.09613015 Inf -0.089085785
## 7            Other party            0.03078653 0.01212030 Inf  0.007031181
## 8  Pro-environment party           -0.06150590 0.05074020 Inf -0.160954870
##    asymp.UCL
## 1 0.01536135
## 2 0.01303887
## 3 0.08470775
## 4 0.15968665
## 5 0.04667679
## 6 0.28773749
## 7 0.05454189
## 8 0.03794307
```

```r
EX9RR.NCmod10.quadratic.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.quadratic.tab$environ.lvl1.sq.trend/
                   EX9RR.NCmod10.quadratic.tab$SE)))
EX9RR.NCmod10.quadratic.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.quadratic.tab$p,method="holm")

EX9RR.NCmod10.quadratic.tab<-
  cbind(group=EX9RR.NCmod10.quadratic.tab[,1],
        round(EX9RR.NCmod10.quadratic.tab[,c(2,3)],2),
        round(EX9RR.NCmod10.quadratic.tab[,c(7,8)],4),
        round(EX9RR.NCmod10.quadratic.tab[,c(5,6)],2))
EX9RR.NCmod10.quadratic.tab
```

```
##                    group environ.lvl1.sq.trend   SE      p  adj.p asymp.LCL
## 1 Anti-immigration party                 -0.04 0.03 0.1566 0.9397     -0.10
## 2           Did not vote                 -0.03 0.02 0.1846 0.9397     -0.07
## 3             Don't know                 -0.01 0.05 0.9047 0.9738     -0.10
## 4                 NE age                  0.07 0.05 0.1246 0.8723     -0.02
## 5             NE citizen                 -0.03 0.04 0.4869 0.9738     -0.10
## 6               NE other                  0.10 0.10 0.3015 0.9397     -0.09
## 7            Other party                  0.03 0.01 0.0111 0.0887      0.01
## 8  Pro-environment party                 -0.06 0.05 0.2254 0.9397     -0.16
##   asymp.UCL
## 1      0.02
## 2      0.01
## 3      0.08
## 4      0.16
## 5      0.05
## 6      0.29
## 7      0.05
## 8      0.04
```

```r
write.csv2(EX9RR.NCmod10.quadratic.tab,"EX9RR.NCmod10.quadratic.tab.csv")

#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.quadratic, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1           Other party - Anti-immigration party  0.07083165 0.03069187 Inf
## 2 Pro-environment party - Anti-immigration party -0.02146078 0.05803643 Inf
## 3            Pro-environment party - Other party -0.09229244 0.05205859 Inf
##      z.ratio    p.value
## 1  2.3078310 0.02100854
## 2 -0.3697812 0.71154550
## 3 -1.7728568 0.07625243
```

```r
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.quadratic, "del.eff",
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##        all.parties.lvl2_del.eff    estimate         SE  df    z.ratio
## 1 Anti-immigration party effect -0.02468544 0.03843936 Inf -0.6421916
## 2            Other party effect  0.08156204 0.03136857 Inf  2.6001196
## 3  Pro-environment party effect -0.05687661 0.05294964 Inf -1.0741642
##       p.value
## 1 0.520748791
## 2 0.009319127
## 3 0.282749089
```

\newpage

#### Marginal cubic effects


```r
EX9RR.NCmod10.cubic<-
  emtrends(EX9RR.NCmod10.ma,specs = c("all.parties.lvl2"),var=c("environ.lvl1.cu"))
(EX9RR.NCmod10.cubic.tab<-data.frame(EX9RR.NCmod10.cubic))
```

```
##         all.parties.lvl2 environ.lvl1.cu.trend          SE  df   asymp.LCL
## 1 Anti-immigration party         -0.0058793709 0.018017522 Inf -0.04119307
## 2           Did not vote          0.0096346446 0.013120965 Inf -0.01608197
## 3             Don't know         -0.0152933277 0.027720604 Inf -0.06962471
## 4                 NE age          0.0702186332 0.027633108 Inf  0.01605874
## 5             NE citizen         -0.0007571017 0.024378256 Inf -0.04853761
## 6               NE other          0.1080511334 0.065318131 Inf -0.01997005
## 7            Other party          0.0058757166 0.008870916 Inf -0.01151096
## 8  Pro-environment party         -0.0376678540 0.026482679 Inf -0.08957295
##    asymp.UCL
## 1 0.02943432
## 2 0.03535126
## 3 0.03903806
## 4 0.12437853
## 5 0.04702340
## 6 0.23607232
## 7 0.02326239
## 8 0.01423724
```

```r
EX9RR.NCmod10.cubic.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.cubic.tab$environ.lvl1.cu.trend/
                   EX9RR.NCmod10.cubic.tab$SE)))
EX9RR.NCmod10.cubic.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.cubic.tab$p,method="holm")

EX9RR.NCmod10.cubic.tab<-
  cbind(group=EX9RR.NCmod10.cubic.tab[,1],
        round(EX9RR.NCmod10.cubic.tab[,c(2,3)],2),
        round(EX9RR.NCmod10.cubic.tab[,c(7,8)],4),
        round(EX9RR.NCmod10.cubic.tab[,c(5,6)],2))
EX9RR.NCmod10.cubic.tab
```

```
##                    group environ.lvl1.cu.trend   SE      p  adj.p asymp.LCL
## 1 Anti-immigration party                 -0.01 0.02 0.7442 1.0000     -0.04
## 2           Did not vote                  0.01 0.01 0.4628 1.0000     -0.02
## 3             Don't know                 -0.02 0.03 0.5812 1.0000     -0.07
## 4                 NE age                  0.07 0.03 0.0111 0.0884      0.02
## 5             NE citizen                  0.00 0.02 0.9752 1.0000     -0.05
## 6               NE other                  0.11 0.07 0.0981 0.6866     -0.02
## 7            Other party                  0.01 0.01 0.5077 1.0000     -0.01
## 8  Pro-environment party                 -0.04 0.03 0.1549 0.9295     -0.09
##   asymp.UCL
## 1      0.03
## 2      0.04
## 3      0.04
## 4      0.12
## 5      0.05
## 6      0.24
## 7      0.02
## 8      0.01
```

```r
write.csv2(EX9RR.NCmod10.cubic.tab,"EX9RR.NCmod10.cubic.tab.csv")

#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.cubic, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1           Other party - Anti-immigration party  0.01175509 0.01892451 Inf
## 2 Pro-environment party - Anti-immigration party -0.03178848 0.03128832 Inf
## 3            Pro-environment party - Other party -0.04354357 0.02720831 Inf
##      z.ratio   p.value
## 1  0.6211566 0.5344966
## 2 -1.0159858 0.3096362
## 3 -1.6003779 0.1095148
```

```r
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.cubic, "del.eff",
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##        all.parties.lvl2_del.eff    estimate         SE  df    z.ratio   p.value
## 1 Anti-immigration party effect  0.01001670 0.02198806 Inf  0.4555516 0.6487125
## 2            Other party effect  0.02764933 0.01744921 Inf  1.5845606 0.1130662
## 3  Pro-environment party effect -0.03766603 0.02775051 Inf -1.3573092 0.1746830
```

\newpage

### Simple slopes for each voting group at -1SD, mean, and +1SD (grand mean used to enhance comparisons between voting groups)

* Because the location of points on grand-mean vary for group-means, the locations of grand mean points for each voting group are calculated separately

#### This script would produce the slopes for group means (not interpreted in the main text)


```r
EX9RR.NCmod10.slopes<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                               at=list(environ.lvl1=c(mean(dat$environ.lvl1)-sd(dat$environ.lvl1),
                                                      mean(dat$environ.lvl1)-0*sd(dat$environ.lvl1),
                                                      mean(dat$environ.lvl1)+sd(dat$environ.lvl1))))
(EX9RR.NCmod10.slopes.tab<-data.frame(EX9RR.NCmod10.slopes))
```

```
##    environ.lvl1       all.parties.lvl2 environ.lvl1.trend         SE  df
## 1  -0.741754824 Anti-immigration party         0.13330284 0.05025507 Inf
## 2   0.005700588 Anti-immigration party         0.08307062 0.04525014 Inf
## 3   0.753156000 Anti-immigration party         0.01312994 0.06357303 Inf
## 4  -0.741754824           Did not vote         0.15233514 0.04300194 Inf
## 5   0.005700588           Did not vote         0.09570419 0.03843351 Inf
## 6   0.753156000           Did not vote         0.07136990 0.04896252 Inf
## 7  -0.741754824             Don't know         0.21242924 0.07410781 Inf
## 8   0.005700588             Don't know         0.22924175 0.06529539 Inf
## 9   0.753156000             Don't know         0.19478892 0.09724522 Inf
## 10 -0.741754824                 NE age         0.13000874 0.07178727 Inf
## 11  0.005700588                 NE age         0.11984701 0.05942368 Inf
## 12  0.753156000                 NE age         0.34506779 0.09580586 Inf
## 13 -0.741754824             NE citizen         0.13849892 0.06660817 Inf
## 14  0.005700588             NE citizen         0.10138160 0.06179577 Inf
## 15  0.753156000             NE citizen         0.06172637 0.07775873 Inf
## 16 -0.741754824               NE other        -0.04356062 0.16346230 Inf
## 17  0.005700588               NE other        -0.07208818 0.14404670 Inf
## 18  0.753156000               NE other         0.26158652 0.18645956 Inf
## 19 -0.741754824            Other party         0.16228911 0.03096371 Inf
## 20  0.005700588            Other party         0.19868655 0.02869814 Inf
## 21  0.753156000            Other party         0.25478021 0.03536920 Inf
## 22 -0.741754824  Pro-environment party         0.23758572 0.06505555 Inf
## 23  0.005700588  Pro-environment party         0.20734788 0.04946009 Inf
## 24  0.753156000  Pro-environment party         0.05084221 0.11120028 Inf
##       asymp.LCL asymp.UCL
## 1   0.034804722 0.2318010
## 2  -0.005618022 0.1717593
## 3  -0.111470903 0.1377308
## 4   0.068052879 0.2366174
## 5   0.020375890 0.1710325
## 6  -0.024594884 0.1673347
## 7   0.067180597 0.3576779
## 8   0.101265146 0.3572184
## 9   0.004191783 0.3853861
## 10 -0.010691712 0.2707092
## 11  0.003378736 0.2363153
## 12  0.157291749 0.5328438
## 13  0.007949316 0.2690485
## 14 -0.019735887 0.2224991
## 15 -0.090677940 0.2141307
## 16 -0.363940841 0.2768196
## 17 -0.354414528 0.2102382
## 18 -0.103867514 0.6270405
## 19  0.101601360 0.2229769
## 20  0.142439228 0.2549339
## 21  0.185457849 0.3241026
## 22  0.110079185 0.3650923
## 23  0.110407886 0.3042879
## 24 -0.167106339 0.2687908
```

```r
EX9RR.NCmod10.slopes.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.tab$SE)))
EX9RR.NCmod10.slopes.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.tab$p,method="holm")

EX9RR.NCmod10.slopes.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.tab[,1],group=EX9RR.NCmod10.slopes.tab[,2],
        round(EX9RR.NCmod10.slopes.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.tab
```

```
##       env_point                  group environ.lvl1.trend   SE      p  adj.p
## 1  -0.741754824 Anti-immigration party               0.13 0.05 0.0080 0.1198
## 2   0.005700588 Anti-immigration party               0.08 0.05 0.0664 0.6639
## 3   0.753156000 Anti-immigration party               0.01 0.06 0.8364 1.0000
## 4  -0.741754824           Did not vote               0.15 0.04 0.0004 0.0071
## 5   0.005700588           Did not vote               0.10 0.04 0.0128 0.1788
## 6   0.753156000           Did not vote               0.07 0.05 0.1449 1.0000
## 7  -0.741754824             Don't know               0.21 0.07 0.0042 0.0664
## 8   0.005700588             Don't know               0.23 0.07 0.0004 0.0076
## 9   0.753156000             Don't know               0.19 0.10 0.0452 0.5246
## 10 -0.741754824                 NE age               0.13 0.07 0.0701 0.6639
## 11  0.005700588                 NE age               0.12 0.06 0.0437 0.5246
## 12  0.753156000                 NE age               0.35 0.10 0.0003 0.0060
## 13 -0.741754824             NE citizen               0.14 0.07 0.0376 0.4887
## 14  0.005700588             NE citizen               0.10 0.06 0.1009 0.8071
## 15  0.753156000             NE citizen               0.06 0.08 0.4273 1.0000
## 16 -0.741754824               NE other              -0.04 0.16 0.7899 1.0000
## 17  0.005700588               NE other              -0.07 0.14 0.6168 1.0000
## 18  0.753156000               NE other               0.26 0.19 0.1606 1.0000
## 19 -0.741754824            Other party               0.16 0.03 0.0000 0.0000
## 20  0.005700588            Other party               0.20 0.03 0.0000 0.0000
## 21  0.753156000            Other party               0.25 0.04 0.0000 0.0000
## 22 -0.741754824  Pro-environment party               0.24 0.07 0.0003 0.0052
## 23  0.005700588  Pro-environment party               0.21 0.05 0.0000 0.0006
## 24  0.753156000  Pro-environment party               0.05 0.11 0.6475 1.0000
##    asymp.LCL asymp.UCL
## 1       0.03      0.23
## 2      -0.01      0.17
## 3      -0.11      0.14
## 4       0.07      0.24
## 5       0.02      0.17
## 6      -0.02      0.17
## 7       0.07      0.36
## 8       0.10      0.36
## 9       0.00      0.39
## 10     -0.01      0.27
## 11      0.00      0.24
## 12      0.16      0.53
## 13      0.01      0.27
## 14     -0.02      0.22
## 15     -0.09      0.21
## 16     -0.36      0.28
## 17     -0.35      0.21
## 18     -0.10      0.63
## 19      0.10      0.22
## 20      0.14      0.25
## 21      0.19      0.32
## 22      0.11      0.37
## 23      0.11      0.30
## 24     -0.17      0.27
```

```r
write.csv2(EX9RR.NCmod10.slopes.tab,"EX9RR.NCmod10.slopes.tab.csv")
```

\newpage

#### Locating the grand-mean points for group-means


```r
all.parties.lvl2.means<-Rdat %>%
  group_by(all.parties.lvl2) %>%
  summarize(env.lvl2.mean=mean(environ.cntrymc),
            env.lvl1.sd=sd(environ.lvl1),
            ref.lvl2.mean=mean(refugees.cntrymc),
            ref.lvl1.sd=sd(refugees.lvl1),
            n=n())

all.parties.lvl2.means<-data.frame(all.parties.lvl2.means)

all.parties.lvl2.means$env.num<-
  (all.parties.lvl2.means$n-1)*
  (all.parties.lvl2.means$env.lvl1.sd^2)

all.parties.lvl2.means$env.pooled.sd<-
  sqrt(sum(all.parties.lvl2.means$env.num)/
         (sum(all.parties.lvl2.means$n)-
            nrow(all.parties.lvl2.means)))

all.parties.lvl2.means$ref.num<-
  (all.parties.lvl2.means$n-1)*
  (all.parties.lvl2.means$ref.lvl1.sd^2)

all.parties.lvl2.means$ref.pooled.sd<-
  sqrt(sum(all.parties.lvl2.means$ref.num)/
         (sum(all.parties.lvl2.means$n)-
            nrow(all.parties.lvl2.means)))
all.parties.lvl2.means
```

```
##         all.parties.lvl2 env.lvl2.mean env.lvl1.sd ref.lvl2.mean ref.lvl1.sd
## 1 Anti-immigration party   -0.13763144   0.7466511   -0.49793401   0.8030285
## 2           Did not vote   -0.11260263   0.7622766   -0.09031957   0.8381866
## 3             Don't know   -0.06464144   0.7025777   -0.09262905   0.7908895
## 4                 NE age    0.01546465   0.6598148    0.21977398   0.7528518
## 5             NE citizen    0.02895292   0.8084471    0.30236880   0.8707147
## 6               NE other    0.19384372   0.7753446    0.30818352   0.7603034
## 7            Other party    0.02509211   0.7579057    0.02283530   0.8118594
## 8  Pro-environment party    0.37454011   0.6664397    0.48916011   0.7665512
##       n    env.num env.pooled.sd    ref.num ref.pooled.sd
## 1  2437 1358.04036     0.7475447 1570.86628     0.8112792
## 2  4295 2495.09569     0.7475447 3016.77871     0.8112792
## 3   884  435.86246     0.7475447  552.32196     0.8112792
## 4  1187  516.33163     0.7475447  672.20807     0.8112792
## 5   969  632.67197     0.7475447  733.88346     0.8112792
## 6   134   79.95418     0.7475447   76.88216     0.8112792
## 7 13386 7688.62522     0.7475447 8822.26325     0.8112792
## 8  1696  752.82044     0.7475447  995.98315     0.8112792
```

##### group.1: Anti-immigration party


```r
x.points.group.1<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]))

x.points.group.1
```

```
## [1] -0.5985121  0.1490326  0.8965774
```

```r
EX9RR.NCmod10.slopes.group.1<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.1)))
EX9RR.NCmod10.slopes.group.1.tab<-data.frame(EX9RR.NCmod10.slopes.group.1)
(EX9RR.NCmod10.slopes.group.1.tab<-
    EX9RR.NCmod10.slopes.group.1.tab[EX9RR.NCmod10.slopes.group.1.tab$all.parties.lvl2=="Anti-immigration party",])
```

```
##   environ.lvl1       all.parties.lvl2 environ.lvl1.trend         SE  df
## 1   -0.5985121 Anti-immigration party        0.125202877 0.04948049 Inf
## 2    0.1490326 Anti-immigration party        0.071186088 0.04397956 Inf
## 3    0.8965774 Anti-immigration party       -0.002543871 0.07739172 Inf
##     asymp.LCL asymp.UCL
## 1  0.02822289 0.2221829
## 2 -0.01501226 0.1573844
## 3 -0.15422886 0.1491411
```

```r
EX9RR.NCmod10.slopes.group.1.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.1.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.1.tab$SE)))
EX9RR.NCmod10.slopes.group.1.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.1.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.1.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.1.tab[,1],group=EX9RR.NCmod10.slopes.group.1.tab[,2],
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.1.tab
```

```
##    env_point                  group environ.lvl1.trend   SE      p  adj.p
## 1 -0.5985121 Anti-immigration party               0.13 0.05 0.0114 0.0342
## 2  0.1490326 Anti-immigration party               0.07 0.04 0.1055 0.2111
## 3  0.8965774 Anti-immigration party               0.00 0.08 0.9738 0.9738
##   asymp.LCL asymp.UCL
## 1      0.03      0.22
## 2     -0.02      0.16
## 3     -0.15      0.15
```

\newpage

##### group.2: Did not vote -group


```r
x.points.group.2<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]))

x.points.group.2
```

```
## [1] -0.6235409  0.1240038  0.8715485
```

```r
EX9RR.NCmod10.slopes.group.2<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.2)))
EX9RR.NCmod10.slopes.group.2.tab<-data.frame(EX9RR.NCmod10.slopes.group.2)
(EX9RR.NCmod10.slopes.group.2.tab<-
    EX9RR.NCmod10.slopes.group.2.tab[EX9RR.NCmod10.slopes.group.2.tab$all.parties.lvl2=="Did not vote",])
```

```
##   environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL
## 4   -0.6235409     Did not vote         0.14122864 0.04189849 Inf  0.05910911
## 5    0.1240038     Did not vote         0.08970135 0.03781345 Inf  0.01558834
## 6    0.8715485     Did not vote         0.07047843 0.05592312 Inf -0.03912888
##   asymp.UCL
## 4 0.2233482
## 5 0.1638144
## 6 0.1800857
```

```r
EX9RR.NCmod10.slopes.group.2.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.2.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.2.tab$SE)))
EX9RR.NCmod10.slopes.group.2.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.2.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.2.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.2.tab[,1],group=EX9RR.NCmod10.slopes.group.2.tab[,2],
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.2.tab
```

```
##    env_point        group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 4 -0.6235409 Did not vote               0.14 0.04 0.0007 0.0022      0.06
## 5  0.1240038 Did not vote               0.09 0.04 0.0177 0.0354      0.02
## 6  0.8715485 Did not vote               0.07 0.06 0.2076 0.2076     -0.04
##   asymp.UCL
## 4      0.22
## 5      0.16
## 6      0.18
```

\newpage

##### group.3: Don't know -group


```r
x.points.group.3<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]))

x.points.group.3
```

```
## [1] -0.67150212  0.07604262  0.82358736
```

```r
EX9RR.NCmod10.slopes.group.3<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.3)))
EX9RR.NCmod10.slopes.group.3.tab<-data.frame(EX9RR.NCmod10.slopes.group.3)
(EX9RR.NCmod10.slopes.group.3.tab<-
    EX9RR.NCmod10.slopes.group.3.tab[EX9RR.NCmod10.slopes.group.3.tab$all.parties.lvl2=="Don't know",])
```

```
##   environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL
## 7  -0.67150212       Don't know          0.2161922 0.07372206 Inf  0.07169961
## 8   0.07604262       Don't know          0.2281847 0.06392488 Inf  0.10289421
## 9   0.82358736       Don't know          0.1888996 0.10793927 Inf -0.02265750
##   asymp.UCL
## 7 0.3606848
## 8 0.3534752
## 9 0.4004567
```

```r
EX9RR.NCmod10.slopes.group.3.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.3.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.3.tab$SE)))
EX9RR.NCmod10.slopes.group.3.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.3.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.3.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.3.tab[,1],group=EX9RR.NCmod10.slopes.group.3.tab[,2],
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.3.tab
```

```
##     env_point      group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 7 -0.67150212 Don't know               0.22 0.07 0.0034 0.0067      0.07
## 8  0.07604262 Don't know               0.23 0.06 0.0004 0.0011      0.10
## 9  0.82358736 Don't know               0.19 0.11 0.0801 0.0801     -0.02
##   asymp.UCL
## 7      0.36
## 8      0.35
## 9      0.40
```

\newpage

##### group.4: Not eligible Age -group


```r
x.points.group.4<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]))

x.points.group.4
```

```
## [1] -0.751608216 -0.004063477  0.743481262
```

```r
EX9RR.NCmod10.slopes.group.4<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.4)))
EX9RR.NCmod10.slopes.group.4.tab<-data.frame(EX9RR.NCmod10.slopes.group.4)
(EX9RR.NCmod10.slopes.group.4.tab<-
    EX9RR.NCmod10.slopes.group.4.tab[EX9RR.NCmod10.slopes.group.4.tab$all.parties.lvl2=="NE age",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL
## 10 -0.751608216           NE age          0.1317146 0.07193465 Inf -0.009274686
## 11 -0.004063477           NE age          0.1184624 0.05964229 Inf  0.001565680
## 12  0.743481262           NE age          0.3406490 0.09436243 Inf  0.155702035
##    asymp.UCL
## 10 0.2727039
## 11 0.2353592
## 12 0.5255960
```

```r
EX9RR.NCmod10.slopes.group.4.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.4.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.4.tab$SE)))
EX9RR.NCmod10.slopes.group.4.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.4.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.4.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.4.tab[,1],group=EX9RR.NCmod10.slopes.group.4.tab[,2],
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.4.tab
```

```
##       env_point  group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 10 -0.751608216 NE age               0.13 0.07 0.0671 0.0940     -0.01
## 11 -0.004063477 NE age               0.12 0.06 0.0470 0.0940      0.00
## 12  0.743481262 NE age               0.34 0.09 0.0003 0.0009      0.16
##    asymp.UCL
## 10      0.27
## 11      0.24
## 12      0.53
```

\newpage

##### group.5: Not eligible citizenship -group


```r
x.points.group.5<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]))

x.points.group.5
```

```
## [1] -0.76509648 -0.01755174  0.72999300
```

```r
EX9RR.NCmod10.slopes.group.5<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.5)))
EX9RR.NCmod10.slopes.group.5.tab<-data.frame(EX9RR.NCmod10.slopes.group.5)
(EX9RR.NCmod10.slopes.group.5.tab<-
    EX9RR.NCmod10.slopes.group.5.tab[EX9RR.NCmod10.slopes.group.5.tab$all.parties.lvl2=="NE citizen",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL
## 13  -0.76509648       NE citizen         0.13961717 0.06669871 Inf  0.008890088
## 14  -0.01755174       NE citizen         0.10257452 0.06222661 Inf -0.019387389
## 15   0.72999300       NE citizen         0.06299336 0.07542719 Inf -0.084841223
##    asymp.UCL
## 13 0.2703442
## 14 0.2245364
## 15 0.2108279
```

```r
EX9RR.NCmod10.slopes.group.5.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.5.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.5.tab$SE)))
EX9RR.NCmod10.slopes.group.5.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.5.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.5.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.5.tab[,1],group=EX9RR.NCmod10.slopes.group.5.tab[,2],
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.5.tab
```

```
##      env_point      group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 13 -0.76509648 NE citizen               0.14 0.07 0.0363 0.1090      0.01
## 14 -0.01755174 NE citizen               0.10 0.06 0.0993 0.1985     -0.02
## 15  0.72999300 NE citizen               0.06 0.08 0.4036 0.4036     -0.08
##    asymp.UCL
## 13      0.27
## 14      0.22
## 15      0.21
```

\newpage

##### group.6: Not eligible Other -group


```r
x.points.group.6<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]))

x.points.group.6
```

```
## [1] -0.9299873 -0.1824425  0.5651022
```

```r
EX9RR.NCmod10.slopes.group.6<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.6)))
EX9RR.NCmod10.slopes.group.6.tab<-data.frame(EX9RR.NCmod10.slopes.group.6)
(EX9RR.NCmod10.slopes.group.6.tab<-
    EX9RR.NCmod10.slopes.group.6.tab[EX9RR.NCmod10.slopes.group.6.tab$all.parties.lvl2=="NE other",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend        SE  df  asymp.LCL
## 16   -0.9299873         NE other         0.02071562 0.1866336 Inf -0.3450795
## 17   -0.1824425         NE other        -0.09901837 0.1484839 Inf -0.3900415
## 18    0.5651022         NE other         0.14353647 0.1518875 Inf -0.1541576
##    asymp.UCL
## 16 0.3865107
## 17 0.1920048
## 18 0.4412305
```

```r
EX9RR.NCmod10.slopes.group.6.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.6.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.6.tab$SE)))
EX9RR.NCmod10.slopes.group.6.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.6.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.6.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.6.tab[,1],group=EX9RR.NCmod10.slopes.group.6.tab[,2],
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.6.tab
```

```
##     env_point    group environ.lvl1.trend   SE      p adj.p asymp.LCL asymp.UCL
## 16 -0.9299873 NE other               0.02 0.19 0.9116     1     -0.35      0.39
## 17 -0.1824425 NE other              -0.10 0.15 0.5049     1     -0.39      0.19
## 18  0.5651022 NE other               0.14 0.15 0.3446     1     -0.15      0.44
```

##### group.7: Other party -group


```r
x.points.group.7<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]))

x.points.group.7
```

```
## [1] -0.76123567 -0.01369093  0.73385381
```

```r
EX9RR.NCmod10.slopes.group.7<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.7)))
EX9RR.NCmod10.slopes.group.7.tab<-data.frame(EX9RR.NCmod10.slopes.group.7)
(EX9RR.NCmod10.slopes.group.7.tab<-
    EX9RR.NCmod10.slopes.group.7.tab[EX9RR.NCmod10.slopes.group.7.tab$all.parties.lvl2=="Other party",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df asymp.LCL
## 19  -0.76123567      Other party          0.1616038 0.03110747 Inf 0.1006343
## 20  -0.01369093      Other party          0.1974934 0.02875290 Inf 0.1411388
## 21   0.73385381      Other party          0.2530839 0.03477012 Inf 0.1849357
##    asymp.UCL
## 19 0.2225734
## 20 0.2538481
## 21 0.3212321
```

```r
EX9RR.NCmod10.slopes.group.7.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.7.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.7.tab$SE)))
EX9RR.NCmod10.slopes.group.7.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.7.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.7.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.7.tab[,1],group=EX9RR.NCmod10.slopes.group.7.tab[,2],
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.7.tab
```

```
##      env_point       group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 19 -0.76123567 Other party               0.16 0.03 0     0      0.10      0.22
## 20 -0.01369093 Other party               0.20 0.03 0     0      0.14      0.25
## 21  0.73385381 Other party               0.25 0.03 0     0      0.18      0.32
```

\newpage

##### group.8: Pro-environment party -group


```r
x.points.group.8<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]))

x.points.group.8
```

```
## [1] -1.1106837 -0.3631389  0.3844058
```

```r
EX9RR.NCmod10.slopes.group.8<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.8)))
EX9RR.NCmod10.slopes.group.8.tab<-data.frame(EX9RR.NCmod10.slopes.group.8)
(EX9RR.NCmod10.slopes.group.8.tab<-
    EX9RR.NCmod10.slopes.group.8.tab[EX9RR.NCmod10.slopes.group.8.tab$all.parties.lvl2=="Pro-environment party",])
```

```
##    environ.lvl1      all.parties.lvl2 environ.lvl1.trend         SE  df
## 22   -1.1106837 Pro-environment party          0.2059681 0.07165753 Inf
## 23   -0.3631389 Pro-environment party          0.2380498 0.05892478 Inf
## 24    0.3844058 Pro-environment party          0.1438335 0.06037701 Inf
##     asymp.LCL asymp.UCL
## 22 0.06552197 0.3464143
## 23 0.12255937 0.3535403
## 24 0.02549673 0.2621703
```

```r
EX9RR.NCmod10.slopes.group.8.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.8.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.8.tab$SE)))
EX9RR.NCmod10.slopes.group.8.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.8.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.8.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.8.tab[,1],group=EX9RR.NCmod10.slopes.group.8.tab[,2],
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.8.tab
```

```
##     env_point                 group environ.lvl1.trend   SE      p  adj.p
## 22 -1.1106837 Pro-environment party               0.21 0.07 0.0040 0.0081
## 23 -0.3631389 Pro-environment party               0.24 0.06 0.0001 0.0002
## 24  0.3844058 Pro-environment party               0.14 0.06 0.0172 0.0172
##    asymp.LCL asymp.UCL
## 22      0.07      0.35
## 23      0.12      0.35
## 24      0.03      0.26
```

### Plotting the effects

#### Make a separate long format data file for the focal groups + young non-voters


```r
EX9.plot.dat<-Rdat %>%
  filter(anti.imm.party.dummy==1 |
           pro.env.party.dummy==1 |
           other.party.dummy==1 |
           not.eligible.age.dummy==1)

#rename the groups
EX9.plot.dat$party<-ifelse(EX9.plot.dat$anti.imm.party.dummy==1,"Anti-immigration voters",
                           ifelse(EX9.plot.dat$pro.env.party.dummy==1,"Pro-environment voters",
                                  ifelse(EX9.plot.dat$other.party.dummy==1,"Other party voters",
                                         ifelse(EX9.plot.dat$not.eligible.age.dummy==1,"Not eligible to vote: Age",NA))))
```

\newpage

#### Plot the curves


```r
non.lin.plot<-
  ggplot(data=EX9.plot.dat,aes(x=environ.gmc,y=refugees))+
  geom_smooth(method="lm",
              se=T,
              formula=y ~ poly(x, 3),
              size=2,color="black")+
  geom_smooth(method="lm",
              se=F,
              formula=y ~ x,
              size=1,
              color="red",linetype="dashed")+
  xlab("Attitudes towards the Environment")+
  ylab("Attitudes towards Refugees")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.25,
                                        linetype = 'dotted',
                                        colour = "black"))+
  coord_cartesian(xlim=c(-2,2))+
  facet_wrap(~party,ncol=2)+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.background = element_blank())

non.lin.plot
```

![](Exploratory_non_linearity_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```r
ggsave(plot = non.lin.plot,
       filename="non.lin.plot.png",device = "png",
       units = "cm",width=12,height=18,dpi = 600)
```

