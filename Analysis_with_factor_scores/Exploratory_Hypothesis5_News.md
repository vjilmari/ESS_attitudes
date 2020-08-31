---
title: "Exlopratory Hypothesis 5 with Political News Consumption"
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


# Hypothesis 5: The strength of the association between environment and refugee attitudes is stronger among more politically engaged individuals. 

Omit missing values on political interest item


```r
dat.H5.news<-dat %>%
  filter(!is.na(polnews))
```

### Model 1: without interactions (only main effects)


```r
H5.exp.news.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


isSingular(H5.exp.news.mod1)
```

```
## [1] FALSE
```

```r
(FE.H5.exp.news.mod1<-getFE(H5.exp.news.mod1))
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
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.27       0.03
## polnews.lvl1                                                0.02       0.01
##                                                               df t.value     p
## (Intercept)                                             18709.32    0.36 0.719
## age                                                     25665.01   -2.07 0.038
## gender                                                  26614.61    5.66 0.000
## educ                                                    26681.09   12.26 0.000
## resid                                                   26707.37   -7.22 0.000
## occupClerical support workers                           26572.28   -0.28 0.780
## occupCraft and related trades workers                   26579.82   -1.20 0.231
## occupElementary occupations                             26582.36    0.09 0.924
## occupManagers                                           26576.50    0.42 0.676
## occupOther: Not in paid work                            26678.03    1.60 0.109
## occupPlant and machine operators, and assemblers        26581.67   -0.68 0.498
## occupProfessionals                                      26571.71    1.57 0.117
## occupRetired                                            26575.11    0.39 0.695
## occupService and sales workers                          26573.08   -0.57 0.571
## occupSkilled agricultural, forestry and fishery workers 26581.33   -0.17 0.867
## occupTechnicians and associate professionals            26569.73   -0.12 0.906
## occupUnemployed                                         26597.91    0.09 0.925
## environ.lvl1                                               15.90    8.23 0.000
## polnews.lvl1                                            26600.02    2.35 0.019
##                                                            LL    UL
## (Intercept)                                             -0.17  0.25
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.34
## polnews.lvl1                                             0.00  0.03
```

```r
(VC.H5.exp.news.mod1<-getVC(H5.exp.news.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41807641 0.174787881
## 2 voting.group environ.lvl1         <NA> 0.08888826 0.007901123
## 3 voting.group  (Intercept) environ.lvl1 0.62444649 0.023205733
## 4        cntry environ.lvl1         <NA> 0.11977394 0.014345797
## 5     Residual         <NA>         <NA> 1.16240969 1.351196289
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political polnews


```r
H5.exp.news.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.news.mod2)
```

```
## [1] FALSE
```

```r
anova(H5.exp.news.mod1,H5.exp.news.mod2)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.news.mod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod1:     polnews.lvl1
## H5.exp.news.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.news.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod2:     polnews.lvl1 + environ.lvl1:polnews.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.news.mod1   24 84682 84879 -42317    84634                     
## H5.exp.news.mod2   25 84684 84889 -42317    84634 0.1967  1     0.6574
```

```r
(FE.H5.exp.news.mod2<-getFE(H5.exp.news.mod2))
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
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.27       0.03
## polnews.lvl1                                                0.02       0.01
## environ.lvl1:polnews.lvl1                                   0.00       0.01
##                                                               df t.value     p
## (Intercept)                                             18710.81    0.36 0.717
## age                                                     25664.52   -2.08 0.038
## gender                                                  26614.58    5.66 0.000
## educ                                                    26681.05   12.26 0.000
## resid                                                   26707.36   -7.22 0.000
## occupClerical support workers                           26572.09   -0.28 0.777
## occupCraft and related trades workers                   26579.66   -1.20 0.229
## occupElementary occupations                             26582.21    0.09 0.928
## occupManagers                                           26576.37    0.41 0.679
## occupOther: Not in paid work                            26677.87    1.60 0.109
## occupPlant and machine operators, and assemblers        26581.53   -0.68 0.496
## occupProfessionals                                      26571.58    1.56 0.118
## occupRetired                                            26575.07    0.39 0.696
## occupService and sales workers                          26572.91   -0.57 0.569
## occupSkilled agricultural, forestry and fishery workers 26581.27   -0.17 0.865
## occupTechnicians and associate professionals            26569.65   -0.12 0.904
## occupUnemployed                                         26597.57    0.09 0.928
## environ.lvl1                                               15.90    8.24 0.000
## polnews.lvl1                                            26599.74    2.35 0.019
## environ.lvl1:polnews.lvl1                               26428.11    0.44 0.657
##                                                            LL    UL
## (Intercept)                                             -0.17  0.25
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.34
## polnews.lvl1                                             0.00  0.03
## environ.lvl1:polnews.lvl1                               -0.01  0.02
```

```r
(VC.H5.exp.news.mod2<-getVC(H5.exp.news.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41804977 0.174765612
## 2 voting.group environ.lvl1         <NA> 0.08896189 0.007914218
## 3 voting.group  (Intercept) environ.lvl1 0.62446234 0.023224066
## 4        cntry environ.lvl1         <NA> 0.11964974 0.014316060
## 5     Residual         <NA>         <NA> 1.16240540 1.351186313
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political polnews, allow polnews effect to vary between voting groups and countries


```r
H5.exp.news.mod3<-lmer(refugees~(environ.lvl1+polnews.lvl1|voting.group)+
                (0+environ.lvl1+polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.news.mod3)
```

```
## [1] FALSE
```

```r
anova(H5.exp.news.mod2,H5.exp.news.mod3)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.news.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod2:     polnews.lvl1 + environ.lvl1:polnews.lvl1
## H5.exp.news.mod3: refugees ~ (environ.lvl1 + polnews.lvl1 | voting.group) + (0 + 
## H5.exp.news.mod3:     environ.lvl1 + polnews.lvl1 | cntry) + age + gender + educ + 
## H5.exp.news.mod3:     resid + occup + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.exp.news.mod2   25 84684 84889 -42317    84634                       
## H5.exp.news.mod3   30 84682 84928 -42311    84622 11.961  5    0.03532 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod3<-getFE(H5.exp.news.mod3))
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
## occupRetired                                                0.04       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.27       0.03
## polnews.lvl1                                                0.02       0.01
## environ.lvl1:polnews.lvl1                                   0.00       0.01
##                                                               df t.value     p
## (Intercept)                                             18679.42    0.35 0.723
## age                                                     25071.96   -1.88 0.059
## gender                                                  26608.90    5.72 0.000
## educ                                                    26661.06   12.25 0.000
## resid                                                   26693.48   -7.13 0.000
## occupClerical support workers                           26549.63   -0.27 0.791
## occupCraft and related trades workers                   26560.87   -1.20 0.231
## occupElementary occupations                             26564.69    0.11 0.912
## occupManagers                                           26555.13    0.42 0.672
## occupOther: Not in paid work                            26656.45    1.60 0.110
## occupPlant and machine operators, and assemblers        26561.30   -0.65 0.517
## occupProfessionals                                      26543.16    1.57 0.117
## occupRetired                                            26550.31    0.36 0.720
## occupService and sales workers                          26551.67   -0.56 0.574
## occupSkilled agricultural, forestry and fishery workers 26565.39   -0.16 0.877
## occupTechnicians and associate professionals            26545.30   -0.11 0.913
## occupUnemployed                                         26572.61    0.07 0.945
## environ.lvl1                                               15.90    8.26 0.000
## polnews.lvl1                                               17.55    1.43 0.170
## environ.lvl1:polnews.lvl1                               26425.03    0.46 0.643
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.20  0.28
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.34
## polnews.lvl1                                            -0.01  0.04
## environ.lvl1:polnews.lvl1                               -0.01  0.02
```

```r
(VC.H5.exp.news.mod3<-getVC(H5.exp.news.mod3))
```

```
##             grp         var1         var2       est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.418398210  1.750571e-01
## 2  voting.group environ.lvl1         <NA>  0.088998932  7.920810e-03
## 3  voting.group polnews.lvl1         <NA>  0.043504207  1.892616e-03
## 4  voting.group  (Intercept) environ.lvl1  0.620903228  2.312057e-02
## 5  voting.group  (Intercept) polnews.lvl1  0.022324974  4.063610e-04
## 6  voting.group environ.lvl1 polnews.lvl1 -0.007011017 -2.714545e-05
## 7         cntry environ.lvl1         <NA>  0.119258034  1.422248e-02
## 8         cntry polnews.lvl1         <NA>  0.029899045  8.939529e-04
## 9         cntry environ.lvl1 polnews.lvl1  0.169321146  6.037486e-04
## 10     Residual         <NA>         <NA>  1.161222565  1.348438e+00
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.exp.news.mod4<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+
                   environ.lvl1:polnews.lvl1|voting.group)+
                (0+environ.lvl1+polnews.lvl1+
                   environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(H5.exp.news.mod4)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod3,H5.exp.news.mod4)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod3: refugees ~ (environ.lvl1 + polnews.lvl1 | voting.group) + (0 + 
## H5.exp.news.mod3:     environ.lvl1 + polnews.lvl1 | cntry) + age + gender + educ + 
## H5.exp.news.mod3:     resid + occup + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1
## H5.exp.news.mod4: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod4:     polnews.lvl1 + environ.lvl1:polnews.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.exp.news.mod3   30 84682 84928 -42311    84622                       
## H5.exp.news.mod4   37 84684 84987 -42305    84610 12.036  7    0.09939 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod4<-getFE(H5.exp.news.mod4))
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
## occupRetired                                                0.04       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.27       0.03
## polnews.lvl1                                                0.02       0.01
## environ.lvl1:polnews.lvl1                                   0.00       0.01
##                                                               df t.value     p
## (Intercept)                                             18660.83    0.35 0.728
## age                                                     25072.62   -1.91 0.056
## gender                                                  26600.28    5.72 0.000
## educ                                                    26633.87   12.22 0.000
## resid                                                   26690.63   -7.09 0.000
## occupClerical support workers                           26506.42   -0.26 0.795
## occupCraft and related trades workers                   26510.45   -1.19 0.232
## occupElementary occupations                             26524.79    0.11 0.909
## occupManagers                                           26516.49    0.44 0.664
## occupOther: Not in paid work                            26620.97    1.62 0.105
## occupPlant and machine operators, and assemblers        26520.60   -0.63 0.527
## occupProfessionals                                      26500.41    1.57 0.116
## occupRetired                                            26508.17    0.36 0.715
## occupService and sales workers                          26513.57   -0.55 0.581
## occupSkilled agricultural, forestry and fishery workers 26539.12   -0.15 0.883
## occupTechnicians and associate professionals            26501.74   -0.10 0.918
## occupUnemployed                                         26557.95    0.06 0.955
## environ.lvl1                                               15.79    8.14 0.000
## polnews.lvl1                                               18.11    1.36 0.192
## environ.lvl1:polnews.lvl1                                  16.04    0.16 0.877
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.20  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.34
## polnews.lvl1                                            -0.01  0.04
## environ.lvl1:polnews.lvl1                               -0.03  0.03
```

```r
(VC.H5.exp.news.mod4<-getVC(H5.exp.news.mod4))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.41838878
## 2  voting.group              environ.lvl1                      <NA>  0.08857233
## 3  voting.group              polnews.lvl1                      <NA>  0.04147059
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.03367763
## 5  voting.group               (Intercept)              environ.lvl1  0.60050944
## 6  voting.group               (Intercept)              polnews.lvl1  0.01492126
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1  0.21659263
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.05546997
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.54398363
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.41962970
## 11        cntry              environ.lvl1                      <NA>  0.12021439
## 12        cntry              polnews.lvl1                      <NA>  0.03136057
## 13        cntry environ.lvl1:polnews.lvl1                      <NA>  0.03741365
## 14        cntry              environ.lvl1              polnews.lvl1  0.19935410
## 15        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.84152801
## 16        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1  0.69713212
## 17     Residual                      <NA>                      <NA>  1.16070177
##          est_SD2
## 1   0.1750491706
## 2   0.0078450571
## 3   0.0017198099
## 4   0.0011341825
## 5   0.0222534793
## 6   0.0002588962
## 7   0.0030518641
## 8  -0.0002037493
## 9  -0.0016226519
## 10 -0.0005860679
## 11  0.0144515005
## 12  0.0009834855
## 13  0.0013997809
## 14  0.0007515634
## 15  0.0037849059
## 16  0.0008179544
## 17  1.3472285947
```

```r
theta <- getME(H5.exp.news.mod4,"theta")

## diagonal elements are identifiable because they are fitted
##  with a lower bound of zero ...
diag.element <- getME(H5.exp.news.mod4,"lower")==0
any(theta[diag.element]<1e-5)
```

```
## [1] TRUE
```

```r
round(theta,5)
```

```
##                            voting.group.(Intercept) 
##                                             0.36046 
##               voting.group.environ.lvl1.(Intercept) 
##                                             0.04582 
##               voting.group.polnews.lvl1.(Intercept) 
##                                             0.00053 
##  voting.group.environ.lvl1:polnews.lvl1.(Intercept) 
##                                             0.00628 
##                           voting.group.environ.lvl1 
##                                             0.06102 
##              voting.group.polnews.lvl1.environ.lvl1 
##                                            -0.00288 
## voting.group.environ.lvl1:polnews.lvl1.environ.lvl1 
##                                            -0.02446 
##                           voting.group.polnews.lvl1 
##                                             0.03561 
## voting.group.environ.lvl1:polnews.lvl1.polnews.lvl1 
##                                            -0.01429 
##              voting.group.environ.lvl1:polnews.lvl1 
##                                             0.00000 
##                                  cntry.environ.lvl1 
##                                             0.10357 
##                     cntry.polnews.lvl1.environ.lvl1 
##                                             0.00539 
##        cntry.environ.lvl1:polnews.lvl1.environ.lvl1 
##                                             0.02713 
##                                  cntry.polnews.lvl1 
##                                             0.02648 
##        cntry.environ.lvl1:polnews.lvl1.polnews.lvl1 
##                                             0.01741 
##                     cntry.environ.lvl1:polnews.lvl1 
##                                             0.00000
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and polnews


```r
H5.exp.news.mod5<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|voting.group)+
                (0+environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
isSingular(H5.exp.news.mod5)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod4,H5.exp.news.mod5)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod4: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod4:     polnews.lvl1 + environ.lvl1:polnews.lvl1
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.news.mod4   37 84684 84987 -42305    84610                         
## H5.exp.news.mod5   64 84528 85052 -42200    84400 209.82 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod5<-getFE(H5.exp.news.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.76       0.12
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.06       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals                0.00       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.15       0.05
## polnews.lvl1                                                0.01       0.03
## all.parties.lvl2Did not vote                                0.63       0.09
## all.parties.lvl2Don't know                                  0.59       0.10
## all.parties.lvl2Invalid vote                                0.83       0.50
## all.parties.lvl2NE age                                      1.03       0.10
## all.parties.lvl2NE citizen                                  1.19       0.10
## all.parties.lvl2NE other                                    1.01       0.13
## all.parties.lvl2No answer                                   1.04       0.50
## all.parties.lvl2Other party                                 0.79       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:polnews.lvl1                                   0.00       0.01
## environ.lvl1:all.parties.lvl2Did not vote                   0.06       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   0.24       0.78
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.06
## environ.lvl1:all.parties.lvl2NE citizen                     0.04       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.04       0.13
## environ.lvl1:all.parties.lvl2No answer                      0.24       0.58
## environ.lvl1:all.parties.lvl2Other party                    0.15       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.13       0.06
## polnews.lvl1:all.parties.lvl2Did not vote                  -0.01       0.03
## polnews.lvl1:all.parties.lvl2Don't know                    -0.05       0.05
## polnews.lvl1:all.parties.lvl2Invalid vote                   0.34       0.62
## polnews.lvl1:all.parties.lvl2NE age                        -0.01       0.04
## polnews.lvl1:all.parties.lvl2NE citizen                    -0.02       0.04
## polnews.lvl1:all.parties.lvl2NE other                      -0.07       0.09
## polnews.lvl1:all.parties.lvl2No answer                     -0.44       0.59
## polnews.lvl1:all.parties.lvl2Other party                    0.02       0.03
## polnews.lvl1:all.parties.lvl2Pro-environment party          0.01       0.04
##                                                               df t.value     p
## (Intercept)                                              3085.66   -6.28 0.000
## age                                                     26371.73   -1.53 0.127
## gender                                                  26637.86    5.71 0.000
## educ                                                    26169.01   12.14 0.000
## resid                                                   26705.78   -6.91 0.000
## occupClerical support workers                           26597.87   -0.21 0.836
## occupCraft and related trades workers                   26606.35   -1.13 0.259
## occupElementary occupations                             26617.12    0.16 0.870
## occupManagers                                           26606.44    0.49 0.624
## occupOther: Not in paid work                            26648.86    1.59 0.112
## occupPlant and machine operators, and assemblers        26615.57   -0.58 0.560
## occupProfessionals                                      26587.71    1.63 0.104
## occupRetired                                            26594.91    0.42 0.672
## occupService and sales workers                          26603.05   -0.50 0.620
## occupSkilled agricultural, forestry and fishery workers 26624.64   -0.11 0.915
## occupTechnicians and associate professionals            26592.27   -0.04 0.967
## occupUnemployed                                         26628.66    0.06 0.950
## environ.lvl1                                               54.40    3.20 0.002
## polnews.lvl1                                              106.30    0.42 0.674
## all.parties.lvl2Did not vote                              168.56    7.01 0.000
## all.parties.lvl2Don't know                                234.29    6.00 0.000
## all.parties.lvl2Invalid vote                             1617.35    1.67 0.096
## all.parties.lvl2NE age                                    242.69   10.41 0.000
## all.parties.lvl2NE citizen                                238.56   11.81 0.000
## all.parties.lvl2NE other                                  579.81    7.62 0.000
## all.parties.lvl2No answer                                1650.01    2.06 0.040
## all.parties.lvl2Other party                               211.78   11.92 0.000
## all.parties.lvl2Pro-environment party                     230.73   15.20 0.000
## environ.lvl1:polnews.lvl1                                  15.29    0.03 0.977
## environ.lvl1:all.parties.lvl2Did not vote                  78.55    1.20 0.233
## environ.lvl1:all.parties.lvl2Don't know                   306.92    2.07 0.039
## environ.lvl1:all.parties.lvl2Invalid vote               21590.07    0.30 0.762
## environ.lvl1:all.parties.lvl2NE age                       268.87    3.26 0.001
## environ.lvl1:all.parties.lvl2NE citizen                   206.09    0.67 0.501
## environ.lvl1:all.parties.lvl2NE other                    1078.45    0.29 0.769
## environ.lvl1:all.parties.lvl2No answer                  11418.78    0.41 0.680
## environ.lvl1:all.parties.lvl2Other party                  121.12    3.70 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        225.52    2.31 0.022
## polnews.lvl1:all.parties.lvl2Did not vote                  62.30   -0.34 0.737
## polnews.lvl1:all.parties.lvl2Don't know                   287.96   -1.08 0.280
## polnews.lvl1:all.parties.lvl2Invalid vote               23529.93    0.56 0.576
## polnews.lvl1:all.parties.lvl2NE age                       190.75   -0.32 0.752
## polnews.lvl1:all.parties.lvl2NE citizen                   191.83   -0.49 0.625
## polnews.lvl1:all.parties.lvl2NE other                    1873.69   -0.79 0.429
## polnews.lvl1:all.parties.lvl2No answer                  23363.21   -0.75 0.453
## polnews.lvl1:all.parties.lvl2Other party                  102.95    0.73 0.466
## polnews.lvl1:all.parties.lvl2Pro-environment party        178.87    0.28 0.778
##                                                            LL    UL
## (Intercept)                                             -1.00 -0.52
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.08
## occupClerical support workers                           -0.23  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.20  0.23
## occupManagers                                           -0.16  0.27
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.21  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.06  0.24
## polnews.lvl1                                            -0.04  0.06
## all.parties.lvl2Did not vote                             0.45  0.81
## all.parties.lvl2Don't know                               0.39  0.78
## all.parties.lvl2Invalid vote                            -0.15  1.81
## all.parties.lvl2NE age                                   0.83  1.22
## all.parties.lvl2NE citizen                               0.99  1.38
## all.parties.lvl2NE other                                 0.75  1.28
## all.parties.lvl2No answer                                0.05  2.03
## all.parties.lvl2Other party                              0.66  0.92
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:polnews.lvl1                               -0.03  0.03
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.15
## environ.lvl1:all.parties.lvl2Don't know                  0.01  0.28
## environ.lvl1:all.parties.lvl2Invalid vote               -1.29  1.76
## environ.lvl1:all.parties.lvl2NE age                      0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.08  0.17
## environ.lvl1:all.parties.lvl2NE other                   -0.21  0.28
## environ.lvl1:all.parties.lvl2No answer                  -0.90  1.38
## environ.lvl1:all.parties.lvl2Other party                 0.07  0.22
## environ.lvl1:all.parties.lvl2Pro-environment party       0.02  0.25
## polnews.lvl1:all.parties.lvl2Did not vote               -0.07  0.05
## polnews.lvl1:all.parties.lvl2Don't know                 -0.14  0.04
## polnews.lvl1:all.parties.lvl2Invalid vote               -0.86  1.55
## polnews.lvl1:all.parties.lvl2NE age                     -0.10  0.07
## polnews.lvl1:all.parties.lvl2NE citizen                 -0.11  0.07
## polnews.lvl1:all.parties.lvl2NE other                   -0.23  0.10
## polnews.lvl1:all.parties.lvl2No answer                  -1.60  0.71
## polnews.lvl1:all.parties.lvl2Other party                -0.03  0.07
## polnews.lvl1:all.parties.lvl2Pro-environment party      -0.06  0.09
```

```r
(VC.H5.exp.news.mod5<-getVC(H5.exp.news.mod5))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.25618171
## 2  voting.group              environ.lvl1                      <NA>  0.07369466
## 3  voting.group              polnews.lvl1                      <NA>  0.03526321
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.02305376
## 5  voting.group               (Intercept)              environ.lvl1  0.56674590
## 6  voting.group               (Intercept)              polnews.lvl1 -0.04390335
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1 -0.52370758
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.24651432
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.83988323
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.31892647
## 11        cntry              environ.lvl1                      <NA>  0.12094848
## 12        cntry              polnews.lvl1                      <NA>  0.03331485
## 13        cntry environ.lvl1:polnews.lvl1                      <NA>  0.03907540
## 14        cntry              environ.lvl1              polnews.lvl1  0.17838799
## 15        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.81526381
## 16        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1  0.71523454
## 17     Residual                      <NA>                      <NA>  1.16066661
##          est_SD2
## 1   0.0656290666
## 2   0.0054309035
## 3   0.0012434941
## 4   0.0005314760
## 5   0.0106997231
## 6  -0.0003966136
## 7  -0.0030929920
## 8  -0.0006406194
## 9  -0.0014269106
## 10 -0.0002592712
## 11  0.0146285354
## 12  0.0011098795
## 13  0.0015268865
## 14  0.0007187932
## 15  0.0038530263
## 16  0.0009310860
## 17  1.3471469907
```

\newpage

#### Look among which voting group there is strongest association between polnews and refugee attitudes


```r
H5.exp.news.mod5.trends<-emtrends(H5.exp.news.mod5,specs = c("all.parties.lvl2"),var=c("polnews.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 26763' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 26763)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 26763' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 26763)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.exp.news.mod5.trends.tab<-data.frame(H5.exp.news.mod5.trends))
```

```
##          all.parties.lvl2 polnews.lvl1.trend         SE  df    asymp.LCL
## 1  Anti-immigration party       0.0108471603 0.02572134 Inf -0.039565744
## 2            Did not vote       0.0006122752 0.02035014 Inf -0.039273260
## 3              Don't know      -0.0383947242 0.03951906 Inf -0.115850668
## 4            Invalid vote       0.3551649038 0.61520059 Inf -0.850606087
## 5                  NE age      -0.0023231165 0.03486805 Inf -0.070663231
## 6              NE citizen      -0.0107573307 0.03775427 Inf -0.084754345
## 7                NE other      -0.0565142136 0.08219992 Inf -0.217623090
## 8               No answer      -0.4323222178 0.58983729 Inf -1.588382063
## 9             Other party       0.0304959285 0.01436349 Inf  0.002344002
## 10  Pro-environment party       0.0216755797 0.03085533 Inf -0.038799762
##     asymp.UCL
## 1  0.06126006
## 2  0.04049781
## 3  0.03906122
## 4  1.56093589
## 5  0.06601700
## 6  0.06323968
## 7  0.10459466
## 8  0.72373763
## 9  0.05864785
## 10 0.08215092
```

```r
H5.exp.news.mod5.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.news.mod5.trends.tab$polnews.lvl1.trend/
                   H5.exp.news.mod5.trends.tab$SE)))
H5.exp.news.mod5.trends.tab$adj.p<-
  p.adjust(H5.exp.news.mod5.trends.tab$p,method="holm")

H5.exp.news.mod5.trends.tab<-
  cbind(group=H5.exp.news.mod5.trends.tab[,1],
      round(H5.exp.news.mod5.trends.tab[,c(2,3)],2),
      round(H5.exp.news.mod5.trends.tab[,c(7,8)],4),
      round(H5.exp.news.mod5.trends.tab[,c(5,6)],2))
H5.exp.news.mod5.trends.tab
```

```
##                     group polnews.lvl1.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party               0.01 0.03 0.6732 1.0000     -0.04
## 2            Did not vote               0.00 0.02 0.9760 1.0000     -0.04
## 3              Don't know              -0.04 0.04 0.3313 1.0000     -0.12
## 4            Invalid vote               0.36 0.62 0.5637 1.0000     -0.85
## 5                  NE age               0.00 0.03 0.9469 1.0000     -0.07
## 6              NE citizen              -0.01 0.04 0.7757 1.0000     -0.08
## 7                NE other              -0.06 0.08 0.4918 1.0000     -0.22
## 8               No answer              -0.43 0.59 0.4636 1.0000     -1.59
## 9             Other party               0.03 0.01 0.0337 0.3374      0.00
## 10  Pro-environment party               0.02 0.03 0.4824 1.0000     -0.04
##    asymp.UCL
## 1       0.06
## 2       0.04
## 3       0.04
## 4       1.56
## 5       0.07
## 6       0.06
## 7       0.10
## 8       0.72
## 9       0.06
## 10      0.08
```

```r
write.csv2(H5.exp.news.mod5.trends.tab,"H5.exp.news.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.exp.news.contrast<-data.frame(pairs(H5.exp.news.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.01082842 0.03826971 Inf
##     z.ratio   p.value
## 1 0.2829501 0.7772151
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.news.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  0.02555 0.0985 Inf  0.259  0.7953 
##  Did not vote effect            0.01418 0.0972 Inf  0.146  0.8841 
##  Don't know effect             -0.02916 0.1029 Inf -0.283  0.7769 
##  Invalid vote effect            0.40813 0.6187 Inf  0.660  0.5095 
##  NE age effect                  0.01092 0.1013 Inf  0.108  0.9141 
##  NE citizen effect              0.00155 0.1023 Inf  0.015  0.9879 
##  NE other effect               -0.04929 0.1254 Inf -0.393  0.6942 
##  No answer effect              -0.46686 0.5939 Inf -0.786  0.4318 
##  Other party effect             0.04739 0.0962 Inf  0.493  0.6222 
##  Pro-environment party effect   0.03759 0.1000 Inf  0.376  0.7070 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic
```

```r
#contrast for three voting groups
(H4.more.contrasts<-data.frame(pairs(H5.exp.news.mod5.trends, 
         exclude=c(2:8), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast     estimate         SE  df
## 1           Other party - Anti-immigration party  0.019648768 0.02685432 Inf
## 2 Pro-environment party - Anti-immigration party  0.010828419 0.03826971 Inf
## 3            Pro-environment party - Other party -0.008820349 0.03174610 Inf
##      z.ratio   p.value
## 1  0.7316798 0.4643640
## 2  0.2829501 0.7772151
## 3 -0.2778404 0.7811349
```

\newpage

### Model 6: Enter three-way interaction voting group x polnews x environment attitudes


```r
H5.exp.news.mod6<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|voting.group)+
                (0+environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2+
                environ.lvl1:polnews.lvl1:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
isSingular(H5.exp.news.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod5,H5.exp.news.mod6)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
## H5.exp.news.mod6: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod6:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod6:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod6:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2 + 
## H5.exp.news.mod6:     environ.lvl1:polnews.lvl1:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)  
## H5.exp.news.mod5   64 84528 85052 -42200    84400                      
## H5.exp.news.mod6   73 84530 85128 -42192    84384 16.16  9    0.06362 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod6<-getFE(H5.exp.news.mod6))
```

```
##                                                                 Estimate
## (Intercept)                                                        -0.76
## age                                                                 0.00
## gender                                                              0.09
## educ                                                                0.03
## resid                                                              -0.10
## occupClerical support workers                                      -0.02
## occupCraft and related trades workers                              -0.12
## occupElementary occupations                                         0.02
## occupManagers                                                       0.05
## occupOther: Not in paid work                                        0.18
## occupPlant and machine operators, and assemblers                   -0.06
## occupProfessionals                                                  0.17
## occupRetired                                                        0.05
## occupService and sales workers                                     -0.05
## occupSkilled agricultural, forestry and fishery workers            -0.01
## occupTechnicians and associate professionals                       -0.01
## occupUnemployed                                                     0.01
## environ.lvl1                                                        0.15
## polnews.lvl1                                                        0.01
## all.parties.lvl2Did not vote                                        0.63
## all.parties.lvl2Don't know                                          0.58
## all.parties.lvl2Invalid vote                                        0.69
## all.parties.lvl2NE age                                              1.02
## all.parties.lvl2NE citizen                                          1.16
## all.parties.lvl2NE other                                            1.02
## all.parties.lvl2No answer                                           0.97
## all.parties.lvl2Other party                                         0.79
## all.parties.lvl2Pro-environment party                               1.31
## environ.lvl1:polnews.lvl1                                          -0.01
## environ.lvl1:all.parties.lvl2Did not vote                           0.05
## environ.lvl1:all.parties.lvl2Don't know                             0.14
## environ.lvl1:all.parties.lvl2Invalid vote                           0.55
## environ.lvl1:all.parties.lvl2NE age                                 0.21
## environ.lvl1:all.parties.lvl2NE citizen                             0.03
## environ.lvl1:all.parties.lvl2NE other                               0.03
## environ.lvl1:all.parties.lvl2No answer                              0.02
## environ.lvl1:all.parties.lvl2Other party                            0.14
## environ.lvl1:all.parties.lvl2Pro-environment party                  0.13
## polnews.lvl1:all.parties.lvl2Did not vote                          -0.01
## polnews.lvl1:all.parties.lvl2Don't know                            -0.05
## polnews.lvl1:all.parties.lvl2Invalid vote                           0.79
## polnews.lvl1:all.parties.lvl2NE age                                -0.01
## polnews.lvl1:all.parties.lvl2NE citizen                            -0.03
## polnews.lvl1:all.parties.lvl2NE other                              -0.06
## polnews.lvl1:all.parties.lvl2No answer                             -0.51
## polnews.lvl1:all.parties.lvl2Other party                            0.02
## polnews.lvl1:all.parties.lvl2Pro-environment party                  0.01
## environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote              0.01
## environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know                0.02
## environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote             -1.55
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE age                    0.01
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen                0.18
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE other                 -0.03
## environ.lvl1:polnews.lvl1:all.parties.lvl2No answer                -0.34
## environ.lvl1:polnews.lvl1:all.parties.lvl2Other party               0.00
## environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party     0.03
##                                                                 Std..Error
## (Intercept)                                                           0.12
## age                                                                   0.00
## gender                                                                0.02
## educ                                                                  0.00
## resid                                                                 0.02
## occupClerical support workers                                         0.11
## occupCraft and related trades workers                                 0.11
## occupElementary occupations                                           0.11
## occupManagers                                                         0.11
## occupOther: Not in paid work                                          0.11
## occupPlant and machine operators, and assemblers                      0.11
## occupProfessionals                                                    0.11
## occupRetired                                                          0.12
## occupService and sales workers                                        0.11
## occupSkilled agricultural, forestry and fishery workers               0.12
## occupTechnicians and associate professionals                          0.11
## occupUnemployed                                                       0.13
## environ.lvl1                                                          0.05
## polnews.lvl1                                                          0.03
## all.parties.lvl2Did not vote                                          0.09
## all.parties.lvl2Don't know                                            0.10
## all.parties.lvl2Invalid vote                                          0.53
## all.parties.lvl2NE age                                                0.10
## all.parties.lvl2NE citizen                                            0.10
## all.parties.lvl2NE other                                              0.13
## all.parties.lvl2No answer                                             0.54
## all.parties.lvl2Other party                                           0.07
## all.parties.lvl2Pro-environment party                                 0.09
## environ.lvl1:polnews.lvl1                                             0.03
## environ.lvl1:all.parties.lvl2Did not vote                             0.05
## environ.lvl1:all.parties.lvl2Don't know                               0.07
## environ.lvl1:all.parties.lvl2Invalid vote                             0.89
## environ.lvl1:all.parties.lvl2NE age                                   0.06
## environ.lvl1:all.parties.lvl2NE citizen                               0.06
## environ.lvl1:all.parties.lvl2NE other                                 0.13
## environ.lvl1:all.parties.lvl2No answer                                0.86
## environ.lvl1:all.parties.lvl2Other party                              0.04
## environ.lvl1:all.parties.lvl2Pro-environment party                    0.06
## polnews.lvl1:all.parties.lvl2Did not vote                             0.03
## polnews.lvl1:all.parties.lvl2Don't know                               0.05
## polnews.lvl1:all.parties.lvl2Invalid vote                             0.86
## polnews.lvl1:all.parties.lvl2NE age                                   0.04
## polnews.lvl1:all.parties.lvl2NE citizen                               0.04
## polnews.lvl1:all.parties.lvl2NE other                                 0.09
## polnews.lvl1:all.parties.lvl2No answer                                0.62
## polnews.lvl1:all.parties.lvl2Other party                              0.03
## polnews.lvl1:all.parties.lvl2Pro-environment party                    0.04
## environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote                0.04
## environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know                  0.06
## environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote                2.11
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE age                      0.06
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen                  0.05
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE other                    0.12
## environ.lvl1:polnews.lvl1:all.parties.lvl2No answer                   1.05
## environ.lvl1:polnews.lvl1:all.parties.lvl2Other party                 0.03
## environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party       0.05
##                                                                       df
## (Intercept)                                                      3096.26
## age                                                             26403.06
## gender                                                          26639.16
## educ                                                            26191.83
## resid                                                           26710.56
## occupClerical support workers                                   26600.06
## occupCraft and related trades workers                           26607.49
## occupElementary occupations                                     26618.45
## occupManagers                                                   26607.49
## occupOther: Not in paid work                                    26649.50
## occupPlant and machine operators, and assemblers                26615.57
## occupProfessionals                                              26588.72
## occupRetired                                                    26596.19
## occupService and sales workers                                  26604.35
## occupSkilled agricultural, forestry and fishery workers         26624.61
## occupTechnicians and associate professionals                    26594.51
## occupUnemployed                                                 26627.71
## environ.lvl1                                                       55.14
## polnews.lvl1                                                      109.12
## all.parties.lvl2Did not vote                                      175.71
## all.parties.lvl2Don't know                                        234.62
## all.parties.lvl2Invalid vote                                     2108.63
## all.parties.lvl2NE age                                            245.39
## all.parties.lvl2NE citizen                                        239.18
## all.parties.lvl2NE other                                          581.08
## all.parties.lvl2No answer                                        2272.00
## all.parties.lvl2Other party                                       214.23
## all.parties.lvl2Pro-environment party                             232.20
## environ.lvl1:polnews.lvl1                                         204.09
## environ.lvl1:all.parties.lvl2Did not vote                          81.09
## environ.lvl1:all.parties.lvl2Don't know                           311.31
## environ.lvl1:all.parties.lvl2Invalid vote                       23627.84
## environ.lvl1:all.parties.lvl2NE age                               271.89
## environ.lvl1:all.parties.lvl2NE citizen                           208.87
## environ.lvl1:all.parties.lvl2NE other                            1104.66
## environ.lvl1:all.parties.lvl2No answer                          23229.74
## environ.lvl1:all.parties.lvl2Other party                          124.88
## environ.lvl1:all.parties.lvl2Pro-environment party                228.89
## polnews.lvl1:all.parties.lvl2Did not vote                          65.51
## polnews.lvl1:all.parties.lvl2Don't know                           300.73
## polnews.lvl1:all.parties.lvl2Invalid vote                       25968.78
## polnews.lvl1:all.parties.lvl2NE age                               199.03
## polnews.lvl1:all.parties.lvl2NE citizen                           199.66
## polnews.lvl1:all.parties.lvl2NE other                            2240.15
## polnews.lvl1:all.parties.lvl2No answer                          24163.24
## polnews.lvl1:all.parties.lvl2Other party                          108.47
## polnews.lvl1:all.parties.lvl2Pro-environment party                187.74
## environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote            237.43
## environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know             1246.02
## environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote          26232.04
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE age                 1297.08
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen              851.13
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE other               6115.77
## environ.lvl1:polnews.lvl1:all.parties.lvl2No answer             26130.28
## environ.lvl1:polnews.lvl1:all.parties.lvl2Other party             406.03
## environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party   750.68
##                                                                 t.value     p
## (Intercept)                                                       -6.25 0.000
## age                                                               -1.49 0.136
## gender                                                             5.70 0.000
## educ                                                              12.18 0.000
## resid                                                             -6.85 0.000
## occupClerical support workers                                     -0.20 0.839
## occupCraft and related trades workers                             -1.13 0.256
## occupElementary occupations                                        0.16 0.872
## occupManagers                                                      0.48 0.630
## occupOther: Not in paid work                                       1.60 0.109
## occupPlant and machine operators, and assemblers                  -0.58 0.559
## occupProfessionals                                                 1.63 0.104
## occupRetired                                                       0.42 0.675
## occupService and sales workers                                    -0.49 0.625
## occupSkilled agricultural, forestry and fishery workers           -0.10 0.920
## occupTechnicians and associate professionals                      -0.05 0.963
## occupUnemployed                                                    0.05 0.957
## environ.lvl1                                                       3.25 0.002
## polnews.lvl1                                                       0.40 0.687
## all.parties.lvl2Did not vote                                       6.94 0.000
## all.parties.lvl2Don't know                                         5.96 0.000
## all.parties.lvl2Invalid vote                                       1.30 0.195
## all.parties.lvl2NE age                                            10.34 0.000
## all.parties.lvl2NE citizen                                        11.58 0.000
## all.parties.lvl2NE other                                           7.61 0.000
## all.parties.lvl2No answer                                          1.78 0.075
## all.parties.lvl2Other party                                       11.88 0.000
## all.parties.lvl2Pro-environment party                             15.11 0.000
## environ.lvl1:polnews.lvl1                                         -0.43 0.670
## environ.lvl1:all.parties.lvl2Did not vote                          1.14 0.258
## environ.lvl1:all.parties.lvl2Don't know                            2.03 0.043
## environ.lvl1:all.parties.lvl2Invalid vote                          0.62 0.536
## environ.lvl1:all.parties.lvl2NE age                                3.22 0.001
## environ.lvl1:all.parties.lvl2NE citizen                            0.51 0.610
## environ.lvl1:all.parties.lvl2NE other                              0.28 0.783
## environ.lvl1:all.parties.lvl2No answer                             0.02 0.981
## environ.lvl1:all.parties.lvl2Other party                           3.65 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party                 2.24 0.026
## polnews.lvl1:all.parties.lvl2Did not vote                         -0.33 0.746
## polnews.lvl1:all.parties.lvl2Don't know                           -1.07 0.287
## polnews.lvl1:all.parties.lvl2Invalid vote                          0.92 0.360
## polnews.lvl1:all.parties.lvl2NE age                               -0.30 0.763
## polnews.lvl1:all.parties.lvl2NE citizen                           -0.61 0.540
## polnews.lvl1:all.parties.lvl2NE other                             -0.69 0.489
## polnews.lvl1:all.parties.lvl2No answer                            -0.82 0.413
## polnews.lvl1:all.parties.lvl2Other party                           0.76 0.450
## polnews.lvl1:all.parties.lvl2Pro-environment party                 0.30 0.767
## environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote             0.37 0.712
## environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know               0.36 0.715
## environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote            -0.73 0.464
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE age                   0.24 0.808
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen               3.38 0.001
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE other                -0.24 0.812
## environ.lvl1:polnews.lvl1:all.parties.lvl2No answer               -0.33 0.742
## environ.lvl1:polnews.lvl1:all.parties.lvl2Other party             -0.02 0.985
## environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party    0.57 0.572
##                                                                    LL    UL
## (Intercept)                                                     -0.99 -0.52
## age                                                              0.00  0.00
## gender                                                           0.06  0.12
## educ                                                             0.02  0.03
## resid                                                           -0.13 -0.07
## occupClerical support workers                                   -0.23  0.19
## occupCraft and related trades workers                           -0.33  0.09
## occupElementary occupations                                     -0.20  0.23
## occupManagers                                                   -0.16  0.26
## occupOther: Not in paid work                                    -0.04  0.40
## occupPlant and machine operators, and assemblers                -0.28  0.15
## occupProfessionals                                              -0.04  0.38
## occupRetired                                                    -0.19  0.29
## occupService and sales workers                                  -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers         -0.24  0.21
## occupTechnicians and associate professionals                    -0.21  0.20
## occupUnemployed                                                 -0.25  0.27
## environ.lvl1                                                     0.06  0.25
## polnews.lvl1                                                    -0.04  0.06
## all.parties.lvl2Did not vote                                     0.45  0.80
## all.parties.lvl2Don't know                                       0.39  0.78
## all.parties.lvl2Invalid vote                                    -0.35  1.73
## all.parties.lvl2NE age                                           0.83  1.22
## all.parties.lvl2NE citizen                                       0.97  1.36
## all.parties.lvl2NE other                                         0.75  1.28
## all.parties.lvl2No answer                                       -0.10  2.03
## all.parties.lvl2Other party                                      0.66  0.93
## all.parties.lvl2Pro-environment party                            1.14  1.48
## environ.lvl1:polnews.lvl1                                       -0.08  0.05
## environ.lvl1:all.parties.lvl2Did not vote                       -0.04  0.15
## environ.lvl1:all.parties.lvl2Don't know                          0.00  0.27
## environ.lvl1:all.parties.lvl2Invalid vote                       -1.19  2.29
## environ.lvl1:all.parties.lvl2NE age                              0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                         -0.09  0.16
## environ.lvl1:all.parties.lvl2NE other                           -0.21  0.28
## environ.lvl1:all.parties.lvl2No answer                          -1.66  1.70
## environ.lvl1:all.parties.lvl2Other party                         0.07  0.22
## environ.lvl1:all.parties.lvl2Pro-environment party               0.02  0.25
## polnews.lvl1:all.parties.lvl2Did not vote                       -0.07  0.05
## polnews.lvl1:all.parties.lvl2Don't know                         -0.14  0.04
## polnews.lvl1:all.parties.lvl2Invalid vote                       -0.90  2.47
## polnews.lvl1:all.parties.lvl2NE age                             -0.09  0.07
## polnews.lvl1:all.parties.lvl2NE citizen                         -0.11  0.06
## polnews.lvl1:all.parties.lvl2NE other                           -0.23  0.11
## polnews.lvl1:all.parties.lvl2No answer                          -1.72  0.71
## polnews.lvl1:all.parties.lvl2Other party                        -0.03  0.07
## polnews.lvl1:all.parties.lvl2Pro-environment party              -0.06  0.09
## environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote          -0.06  0.09
## environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know            -0.10  0.14
## environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote          -5.68  2.59
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE age                -0.10  0.13
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen             0.08  0.29
## environ.lvl1:polnews.lvl1:all.parties.lvl2NE other              -0.26  0.20
## environ.lvl1:polnews.lvl1:all.parties.lvl2No answer             -2.40  1.71
## environ.lvl1:polnews.lvl1:all.parties.lvl2Other party           -0.07  0.06
## environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party -0.07  0.13
```

```r
(VC.H5.exp.news.mod6<-getVC(H5.exp.news.mod6))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.25615535
## 2  voting.group              environ.lvl1                      <NA>  0.07366816
## 3  voting.group              polnews.lvl1                      <NA>  0.03575911
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.02260760
## 5  voting.group               (Intercept)              environ.lvl1  0.56875475
## 6  voting.group               (Intercept)              polnews.lvl1 -0.03833647
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1 -0.50322445
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.21499434
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.83697265
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.35398019
## 11        cntry              environ.lvl1                      <NA>  0.12088630
## 12        cntry              polnews.lvl1                      <NA>  0.03398824
## 13        cntry environ.lvl1:polnews.lvl1                      <NA>  0.03964676
## 14        cntry              environ.lvl1              polnews.lvl1  0.18228064
## 15        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.81937872
## 16        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1  0.71300549
## 17     Residual                      <NA>                      <NA>  1.16029572
##          est_SD2
## 1   0.0656155631
## 2   0.0054269985
## 3   0.0012787137
## 4   0.0005111036
## 5   0.0107326835
## 6  -0.0003511577
## 7  -0.0029142020
## 8  -0.0005663613
## 9  -0.0013939449
## 10 -0.0002861674
## 11  0.0146134972
## 12  0.0011552003
## 13  0.0015718654
## 14  0.0007489387
## 15  0.0039270772
## 16  0.0009607916
## 17  1.3462861596
```

#### Refit with manually coded level-1 interaction


```r
dat.H5.news$env.news.int<-dat.H5.news$environ.lvl1*dat.H5.news$polnews.lvl1

H5.exp.news.mod6<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+env.news.int|voting.group)+
                (0+environ.lvl1+polnews.lvl1+env.news.int|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                env.news.int+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2+
                env.news.int:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
isSingular(H5.exp.news.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod5,H5.exp.news.mod6)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (0 + environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
## H5.exp.news.mod6: refugees ~ (environ.lvl1 + polnews.lvl1 + env.news.int | voting.group) + 
## H5.exp.news.mod6:     (0 + environ.lvl1 + polnews.lvl1 + env.news.int | cntry) + 
## H5.exp.news.mod6:     age + gender + educ + resid + occup + environ.lvl1 + polnews.lvl1 + 
## H5.exp.news.mod6:     env.news.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.exp.news.mod6:     polnews.lvl1:all.parties.lvl2 + env.news.int:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)  
## H5.exp.news.mod5   64 84528 85052 -42200    84400                      
## H5.exp.news.mod6   73 84530 85128 -42192    84384 16.16  9    0.06362 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod6<-getFE(H5.exp.news.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.76       0.12
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.06       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.15       0.05
## polnews.lvl1                                                0.01       0.03
## env.news.int                                               -0.01       0.03
## all.parties.lvl2Did not vote                                0.63       0.09
## all.parties.lvl2Don't know                                  0.58       0.10
## all.parties.lvl2Invalid vote                                0.69       0.53
## all.parties.lvl2NE age                                      1.02       0.10
## all.parties.lvl2NE citizen                                  1.16       0.10
## all.parties.lvl2NE other                                    1.02       0.13
## all.parties.lvl2No answer                                   0.97       0.54
## all.parties.lvl2Other party                                 0.79       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:all.parties.lvl2Did not vote                   0.05       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   0.55       0.89
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.06
## environ.lvl1:all.parties.lvl2NE citizen                     0.03       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.03       0.13
## environ.lvl1:all.parties.lvl2No answer                      0.02       0.86
## environ.lvl1:all.parties.lvl2Other party                    0.14       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.13       0.06
## polnews.lvl1:all.parties.lvl2Did not vote                  -0.01       0.03
## polnews.lvl1:all.parties.lvl2Don't know                    -0.05       0.05
## polnews.lvl1:all.parties.lvl2Invalid vote                   0.79       0.86
## polnews.lvl1:all.parties.lvl2NE age                        -0.01       0.04
## polnews.lvl1:all.parties.lvl2NE citizen                    -0.03       0.04
## polnews.lvl1:all.parties.lvl2NE other                      -0.06       0.09
## polnews.lvl1:all.parties.lvl2No answer                     -0.51       0.62
## polnews.lvl1:all.parties.lvl2Other party                    0.02       0.03
## polnews.lvl1:all.parties.lvl2Pro-environment party          0.01       0.04
## env.news.int:all.parties.lvl2Did not vote                   0.01       0.04
## env.news.int:all.parties.lvl2Don't know                     0.02       0.06
## env.news.int:all.parties.lvl2Invalid vote                  -1.55       2.11
## env.news.int:all.parties.lvl2NE age                         0.01       0.06
## env.news.int:all.parties.lvl2NE citizen                     0.18       0.05
## env.news.int:all.parties.lvl2NE other                      -0.03       0.12
## env.news.int:all.parties.lvl2No answer                     -0.34       1.05
## env.news.int:all.parties.lvl2Other party                    0.00       0.03
## env.news.int:all.parties.lvl2Pro-environment party          0.03       0.05
##                                                               df t.value     p
## (Intercept)                                              3096.25   -6.25 0.000
## age                                                     26403.06   -1.49 0.136
## gender                                                  26639.16    5.70 0.000
## educ                                                    26191.83   12.18 0.000
## resid                                                   26710.56   -6.85 0.000
## occupClerical support workers                           26600.06   -0.20 0.839
## occupCraft and related trades workers                   26607.49   -1.13 0.256
## occupElementary occupations                             26618.45    0.16 0.872
## occupManagers                                           26607.49    0.48 0.630
## occupOther: Not in paid work                            26649.50    1.60 0.109
## occupPlant and machine operators, and assemblers        26615.57   -0.58 0.559
## occupProfessionals                                      26588.72    1.63 0.104
## occupRetired                                            26596.19    0.42 0.675
## occupService and sales workers                          26604.36   -0.49 0.625
## occupSkilled agricultural, forestry and fishery workers 26624.61   -0.10 0.920
## occupTechnicians and associate professionals            26594.51   -0.05 0.963
## occupUnemployed                                         26627.71    0.05 0.957
## environ.lvl1                                               55.14    3.25 0.002
## polnews.lvl1                                              109.12    0.40 0.687
## env.news.int                                              204.09   -0.43 0.670
## all.parties.lvl2Did not vote                              175.71    6.94 0.000
## all.parties.lvl2Don't know                                234.61    5.96 0.000
## all.parties.lvl2Invalid vote                             2108.62    1.30 0.195
## all.parties.lvl2NE age                                    245.39   10.34 0.000
## all.parties.lvl2NE citizen                                239.18   11.58 0.000
## all.parties.lvl2NE other                                  581.08    7.61 0.000
## all.parties.lvl2No answer                                2271.99    1.78 0.075
## all.parties.lvl2Other party                               214.23   11.88 0.000
## all.parties.lvl2Pro-environment party                     232.20   15.11 0.000
## environ.lvl1:all.parties.lvl2Did not vote                  81.09    1.14 0.258
## environ.lvl1:all.parties.lvl2Don't know                   311.30    2.03 0.043
## environ.lvl1:all.parties.lvl2Invalid vote               23627.81    0.62 0.536
## environ.lvl1:all.parties.lvl2NE age                       271.89    3.22 0.001
## environ.lvl1:all.parties.lvl2NE citizen                   208.87    0.51 0.610
## environ.lvl1:all.parties.lvl2NE other                    1104.65    0.28 0.783
## environ.lvl1:all.parties.lvl2No answer                  23229.70    0.02 0.981
## environ.lvl1:all.parties.lvl2Other party                  124.88    3.65 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        228.89    2.24 0.026
## polnews.lvl1:all.parties.lvl2Did not vote                  65.51   -0.33 0.746
## polnews.lvl1:all.parties.lvl2Don't know                   300.73   -1.07 0.287
## polnews.lvl1:all.parties.lvl2Invalid vote               25968.77    0.92 0.360
## polnews.lvl1:all.parties.lvl2NE age                       199.03   -0.30 0.763
## polnews.lvl1:all.parties.lvl2NE citizen                   199.66   -0.61 0.540
## polnews.lvl1:all.parties.lvl2NE other                    2240.13   -0.69 0.489
## polnews.lvl1:all.parties.lvl2No answer                  24163.22   -0.82 0.413
## polnews.lvl1:all.parties.lvl2Other party                  108.47    0.76 0.450
## polnews.lvl1:all.parties.lvl2Pro-environment party        187.74    0.30 0.767
## env.news.int:all.parties.lvl2Did not vote                 237.43    0.37 0.712
## env.news.int:all.parties.lvl2Don't know                  1246.00    0.36 0.715
## env.news.int:all.parties.lvl2Invalid vote               26232.04   -0.73 0.464
## env.news.int:all.parties.lvl2NE age                      1297.07    0.24 0.808
## env.news.int:all.parties.lvl2NE citizen                   851.12    3.38 0.001
## env.news.int:all.parties.lvl2NE other                    6115.71   -0.24 0.812
## env.news.int:all.parties.lvl2No answer                  26130.28   -0.33 0.742
## env.news.int:all.parties.lvl2Other party                  406.02   -0.02 0.985
## env.news.int:all.parties.lvl2Pro-environment party        750.67    0.57 0.572
##                                                            LL    UL
## (Intercept)                                             -0.99 -0.52
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.23  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.20  0.23
## occupManagers                                           -0.16  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.21  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.06  0.25
## polnews.lvl1                                            -0.04  0.06
## env.news.int                                            -0.08  0.05
## all.parties.lvl2Did not vote                             0.45  0.80
## all.parties.lvl2Don't know                               0.39  0.78
## all.parties.lvl2Invalid vote                            -0.35  1.73
## all.parties.lvl2NE age                                   0.83  1.22
## all.parties.lvl2NE citizen                               0.97  1.36
## all.parties.lvl2NE other                                 0.75  1.28
## all.parties.lvl2No answer                               -0.10  2.03
## all.parties.lvl2Other party                              0.66  0.93
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.15
## environ.lvl1:all.parties.lvl2Don't know                  0.00  0.27
## environ.lvl1:all.parties.lvl2Invalid vote               -1.19  2.29
## environ.lvl1:all.parties.lvl2NE age                      0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.09  0.16
## environ.lvl1:all.parties.lvl2NE other                   -0.21  0.28
## environ.lvl1:all.parties.lvl2No answer                  -1.66  1.70
## environ.lvl1:all.parties.lvl2Other party                 0.07  0.22
## environ.lvl1:all.parties.lvl2Pro-environment party       0.02  0.25
## polnews.lvl1:all.parties.lvl2Did not vote               -0.07  0.05
## polnews.lvl1:all.parties.lvl2Don't know                 -0.14  0.04
## polnews.lvl1:all.parties.lvl2Invalid vote               -0.90  2.47
## polnews.lvl1:all.parties.lvl2NE age                     -0.09  0.07
## polnews.lvl1:all.parties.lvl2NE citizen                 -0.11  0.06
## polnews.lvl1:all.parties.lvl2NE other                   -0.23  0.11
## polnews.lvl1:all.parties.lvl2No answer                  -1.72  0.71
## polnews.lvl1:all.parties.lvl2Other party                -0.03  0.07
## polnews.lvl1:all.parties.lvl2Pro-environment party      -0.06  0.09
## env.news.int:all.parties.lvl2Did not vote               -0.06  0.09
## env.news.int:all.parties.lvl2Don't know                 -0.10  0.14
## env.news.int:all.parties.lvl2Invalid vote               -5.68  2.59
## env.news.int:all.parties.lvl2NE age                     -0.10  0.13
## env.news.int:all.parties.lvl2NE citizen                  0.08  0.29
## env.news.int:all.parties.lvl2NE other                   -0.26  0.20
## env.news.int:all.parties.lvl2No answer                  -2.40  1.71
## env.news.int:all.parties.lvl2Other party                -0.07  0.06
## env.news.int:all.parties.lvl2Pro-environment party      -0.07  0.13
```

```r
(VC.H5.exp.news.mod6<-getVC(H5.exp.news.mod6))
```

```
##             grp         var1         var2      est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.25615555  0.0656156659
## 2  voting.group environ.lvl1         <NA>  0.07366860  0.0054270627
## 3  voting.group polnews.lvl1         <NA>  0.03575940  0.0012787347
## 4  voting.group env.news.int         <NA>  0.02260782  0.0005111136
## 5  voting.group  (Intercept) environ.lvl1  0.56874975  0.0107326610
## 6  voting.group  (Intercept) polnews.lvl1 -0.03833964 -0.0003511900
## 7  voting.group  (Intercept) env.news.int -0.50321899 -0.0029142010
## 8  voting.group environ.lvl1 polnews.lvl1 -0.21498186 -0.0005663364
## 9  voting.group environ.lvl1 env.news.int -0.83697554 -0.0013939715
## 10 voting.group polnews.lvl1 env.news.int -0.35398741 -0.0002861783
## 11        cntry environ.lvl1         <NA>  0.12088644  0.0146135321
## 12        cntry polnews.lvl1         <NA>  0.03398813  0.0011551931
## 13        cntry env.news.int         <NA>  0.03964677  0.0015718662
## 14        cntry environ.lvl1 polnews.lvl1  0.18228101  0.0007489388
## 15        cntry environ.lvl1 env.news.int  0.81937631  0.0039270714
## 16        cntry polnews.lvl1 env.news.int  0.71300870  0.0009607932
## 17     Residual         <NA>         <NA>  1.16029569  1.3462860981
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.exp.news.mod6.trends<-emtrends(H5.exp.news.mod6,specs = c("all.parties.lvl2"),var=c("env.news.int"))
(H5.exp.news.mod6.trends.tab<-data.frame(H5.exp.news.mod6.trends))
```

```
##          all.parties.lvl2 env.news.int.trend         SE  df   asymp.LCL
## 1  Anti-immigration party      -1.343582e-02 0.03144273 Inf -0.07506244
## 2            Did not vote       5.254399e-05 0.02341364 Inf -0.04583734
## 3              Don't know       8.421066e-03 0.05289599 Inf -0.09525317
## 4            Invalid vote      -1.558455e+00 2.11149875 Inf -5.69691647
## 5                  NE age       7.887817e-04 0.05116279 Inf -0.09948844
## 6              NE citizen       1.705658e-01 0.04665679 Inf  0.07912019
## 7                NE other      -4.172579e-02 0.11555513 Inf -0.26820969
## 8               No answer      -3.575803e-01 1.04709048 Inf -2.40983988
## 9             Other party      -1.407327e-02 0.01728033 Inf -0.04794209
## 10  Pro-environment party       1.521679e-02 0.04215225 Inf -0.06740011
##     asymp.UCL
## 1  0.04819081
## 2  0.04594243
## 3  0.11209531
## 4  2.58000653
## 5  0.10106600
## 6  0.26201145
## 7  0.18475811
## 8  1.69467936
## 9  0.01979555
## 10 0.09783369
```

```r
H5.exp.news.mod6.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.news.mod6.trends.tab$env.news.int.trend/
                   H5.exp.news.mod6.trends.tab$SE)))
H5.exp.news.mod6.trends.tab$adj.p<-
  p.adjust(H5.exp.news.mod6.trends.tab$p,method="holm")

H5.exp.news.mod6.trends.tab<-
  cbind(group=H5.exp.news.mod6.trends.tab[,1],
      round(H5.exp.news.mod6.trends.tab[,c(2,3)],2),
      round(H5.exp.news.mod6.trends.tab[,c(7,8)],4),
      round(H5.exp.news.mod6.trends.tab[,c(5,6)],2))
H5.exp.news.mod6.trends.tab
```

```
##                     group env.news.int.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party              -0.01 0.03 0.6692 1.0000     -0.08
## 2            Did not vote               0.00 0.02 0.9982 1.0000     -0.05
## 3              Don't know               0.01 0.05 0.8735 1.0000     -0.10
## 4            Invalid vote              -1.56 2.11 0.4605 1.0000     -5.70
## 5                  NE age               0.00 0.05 0.9877 1.0000     -0.10
## 6              NE citizen               0.17 0.05 0.0003 0.0026      0.08
## 7                NE other              -0.04 0.12 0.7180 1.0000     -0.27
## 8               No answer              -0.36 1.05 0.7327 1.0000     -2.41
## 9             Other party              -0.01 0.02 0.4154 1.0000     -0.05
## 10  Pro-environment party               0.02 0.04 0.7181 1.0000     -0.07
##    asymp.UCL
## 1       0.05
## 2       0.05
## 3       0.11
## 4       2.58
## 5       0.10
## 6       0.26
## 7       0.18
## 8       1.69
## 9       0.02
## 10      0.10
```

```r
write.csv2(H5.exp.news.mod6.trends.tab,"H5.exp.news.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.exp.news.contrast<-data.frame(pairs(H5.exp.news.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.02865261 0.05066405 Inf
##     z.ratio   p.value
## 1 0.5655412 0.5717057
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.news.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect    0.184 0.264 Inf  0.697  1.0000 
##  Did not vote effect              0.199 0.263 Inf  0.756  1.0000 
##  Don't know effect                0.208 0.267 Inf  0.779  1.0000 
##  Invalid vote effect             -1.533 2.115 Inf -0.725  1.0000 
##  NE age effect                    0.200 0.267 Inf  0.748  1.0000 
##  NE citizen effect                0.388 0.266 Inf  1.459  1.0000 
##  NE other effect                  0.153 0.286 Inf  0.533  1.0000 
##  No answer effect                -0.198 1.073 Inf -0.185  1.0000 
##  Other party effect               0.183 0.263 Inf  0.697  1.0000 
##  Pro-environment party effect     0.216 0.266 Inf  0.813  1.0000 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H5.exp.news.more.contrasts<-data.frame(pairs(H5.exp.news.mod6.trends, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast      estimate         SE  df
## 1           Other party - Anti-immigration party -0.0006374564 0.03301471 Inf
## 2 Pro-environment party - Anti-immigration party  0.0286526070 0.05066405 Inf
## 3            Pro-environment party - Other party  0.0292900633 0.04328038 Inf
##       z.ratio p.value
## 1 -0.01930825       1
## 2  0.56554122       1
## 3  0.67675157       1
```

