---
title: "Hypothesis 5"
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

### Model 1: without interactions (only main effects)


```r
H5.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


isSingular(H5.mod1)
```

```
## [1] FALSE
```

```r
(FE.H5.mod1<-getFE(H5.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.02       0.13
## environ.lvl1                                                0.26       0.03
## engagement.lvl1                                             0.07       0.01
##                                                               df t.value     p
## (Intercept)                                             18811.69    0.33 0.744
## age                                                     25790.71   -3.21 0.001
## gender                                                  26736.11    6.43 0.000
## educ                                                    26788.06   11.50 0.000
## resid                                                   26831.31   -6.97 0.000
## occupClerical support workers                           26696.39   -0.22 0.827
## occupCraft and related trades workers                   26703.54   -1.15 0.252
## occupElementary occupations                             26706.78    0.19 0.848
## occupManagers                                           26700.12    0.40 0.691
## occupOther: Not in paid work                            26801.21    1.60 0.110
## occupPlant and machine operators, and assemblers        26705.69   -0.64 0.522
## occupProfessionals                                      26695.51    1.56 0.120
## occupRetired                                            26701.77    0.42 0.674
## occupService and sales workers                          26697.30   -0.52 0.601
## occupSkilled agricultural, forestry and fishery workers 26705.67   -0.12 0.904
## occupTechnicians and associate professionals            26693.62   -0.11 0.915
## occupUnemployed                                         26721.95    0.13 0.900
## environ.lvl1                                               15.93    8.07 0.000
## engagement.lvl1                                         26753.06    6.85 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.24  0.28
## environ.lvl1                                             0.20  0.33
## engagement.lvl1                                          0.05  0.09
```

```r
(VC.H5.mod1<-getVC(H5.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41719664 0.174053034
## 2 voting.group environ.lvl1         <NA> 0.08901292 0.007923299
## 3 voting.group  (Intercept) environ.lvl1 0.62813789 0.023326459
## 4        cntry environ.lvl1         <NA> 0.11989870 0.014375698
## 5     Residual         <NA>         <NA> 1.16165188 1.349435086
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political engagement


```r
H5.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                environ.lvl1:engagement.lvl1
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.mod2)
```

```
## [1] FALSE
```

```r
anova(H5.mod1,H5.mod2)
```

```
## Data: dat
## Models:
## H5.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.mod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod1:     engagement.lvl1
## H5.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod2:     engagement.lvl1 + environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.mod1   24 85033 85230 -42493    84985                        
## H5.mod2   25 85026 85231 -42488    84976 9.2338  1   0.002376 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod2<-getFE(H5.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.26       0.03
## engagement.lvl1                                             0.07       0.01
## environ.lvl1:engagement.lvl1                                0.04       0.01
##                                                               df t.value     p
## (Intercept)                                             18840.05    0.33 0.743
## age                                                     25772.74   -3.24 0.001
## gender                                                  26736.22    6.42 0.000
## educ                                                    26786.42   11.52 0.000
## resid                                                   26831.46   -6.95 0.000
## occupClerical support workers                           26696.41   -0.24 0.814
## occupCraft and related trades workers                   26703.54   -1.16 0.244
## occupElementary occupations                             26706.80    0.17 0.865
## occupManagers                                           26700.16    0.38 0.704
## occupOther: Not in paid work                            26801.27    1.59 0.113
## occupPlant and machine operators, and assemblers        26705.72   -0.66 0.509
## occupProfessionals                                      26695.52    1.54 0.124
## occupRetired                                            26701.85    0.41 0.683
## occupService and sales workers                          26697.31   -0.54 0.590
## occupSkilled agricultural, forestry and fishery workers 26705.69   -0.13 0.894
## occupTechnicians and associate professionals            26693.65   -0.12 0.907
## occupUnemployed                                         26721.66    0.10 0.920
## environ.lvl1                                               15.95    8.10 0.000
## engagement.lvl1                                         26753.14    6.89 0.000
## environ.lvl1:engagement.lvl1                            26479.49    3.04 0.002
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.17  0.25
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.19  0.33
## engagement.lvl1                                          0.05  0.09
## environ.lvl1:engagement.lvl1                             0.01  0.06
```

```r
(VC.H5.mod2<-getVC(H5.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41646428 0.173442494
## 2 voting.group environ.lvl1         <NA> 0.09008067 0.008114527
## 3 voting.group  (Intercept) environ.lvl1 0.63107327 0.023674953
## 4        cntry environ.lvl1         <NA> 0.11910701 0.014186480
## 5     Residual         <NA>         <NA> 1.16145695 1.348982243
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political engagement, allow engagement effect to vary between voting groups and countries


```r
H5.mod3<-lmer(refugees~(environ.lvl1+engagement.lvl1|voting.group)+
                (0+environ.lvl1+engagement.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                environ.lvl1:engagement.lvl1
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.mod3)
```

```
## [1] FALSE
```

```r
anova(H5.mod2,H5.mod3)
```

```
## Data: dat
## Models:
## H5.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod2:     engagement.lvl1 + environ.lvl1:engagement.lvl1
## H5.mod3: refugees ~ (environ.lvl1 + engagement.lvl1 | voting.group) + 
## H5.mod3:     (0 + environ.lvl1 + engagement.lvl1 | cntry) + age + gender + 
## H5.mod3:     educ + resid + occup + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.mod2   25 85026 85231 -42488    84976                         
## H5.mod3   30 85009 85255 -42475    84949 26.846  5  6.111e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod3<-getFE(H5.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.26       0.03
## engagement.lvl1                                             0.07       0.02
## environ.lvl1:engagement.lvl1                                0.04       0.01
##                                                               df t.value     p
## (Intercept)                                             18791.40    0.31 0.758
## age                                                     25456.65   -3.11 0.002
## gender                                                  26725.36    6.46 0.000
## educ                                                    26783.99   11.49 0.000
## resid                                                   26806.70   -6.91 0.000
## occupClerical support workers                           26671.79   -0.21 0.834
## occupCraft and related trades workers                   26681.63   -1.15 0.252
## occupElementary occupations                             26685.11    0.20 0.845
## occupManagers                                           26676.68    0.41 0.678
## occupOther: Not in paid work                            26778.74    1.61 0.108
## occupPlant and machine operators, and assemblers        26685.48   -0.62 0.536
## occupProfessionals                                      26666.75    1.54 0.123
## occupRetired                                            26674.30    0.38 0.704
## occupService and sales workers                          26672.48   -0.51 0.613
## occupSkilled agricultural, forestry and fishery workers 26688.35   -0.10 0.918
## occupTechnicians and associate professionals            26668.34   -0.09 0.927
## occupUnemployed                                         26694.04    0.10 0.920
## environ.lvl1                                               15.91    8.12 0.000
## engagement.lvl1                                            17.23    4.45 0.000
## environ.lvl1:engagement.lvl1                            26473.11    3.07 0.002
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.08
## occupClerical support workers                           -0.23  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.19  0.33
## engagement.lvl1                                          0.04  0.10
## environ.lvl1:engagement.lvl1                             0.01  0.06
```

```r
(VC.H5.mod3<-getVC(H5.mod3))
```

```
##             grp            var1            var2     est_SD      est_SD2
## 1  voting.group     (Intercept)            <NA> 0.41691343 0.1738168060
## 2  voting.group    environ.lvl1            <NA> 0.08914639 0.0079470787
## 3  voting.group engagement.lvl1            <NA> 0.07659566 0.0058668956
## 4  voting.group     (Intercept)    environ.lvl1 0.61768280 0.0229570006
## 5  voting.group     (Intercept) engagement.lvl1 0.29126274 0.0093011144
## 6  voting.group    environ.lvl1 engagement.lvl1 0.05783774 0.0003949292
## 7         cntry    environ.lvl1            <NA> 0.11781675 0.0138807855
## 8         cntry engagement.lvl1            <NA> 0.04143857 0.0017171548
## 9         cntry    environ.lvl1 engagement.lvl1 0.19180347 0.0009364147
## 10     Residual            <NA>            <NA> 1.15964026 1.3447655431
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.mod4<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+
                   environ.lvl1:engagement.lvl1|voting.group)+
                (0+environ.lvl1+engagement.lvl1+
                   environ.lvl1:engagement.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                environ.lvl1:engagement.lvl1
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(H5.mod4)
```

```
## [1] TRUE
```

```r
anova(H5.mod3,H5.mod4)
```

```
## Data: dat
## Models:
## H5.mod3: refugees ~ (environ.lvl1 + engagement.lvl1 | voting.group) + 
## H5.mod3:     (0 + environ.lvl1 + engagement.lvl1 | cntry) + age + gender + 
## H5.mod3:     educ + resid + occup + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1
## H5.mod4: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod4:     engagement.lvl1 + environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod3   30 85009 85255 -42475    84949                       
## H5.mod4   37 85006 85309 -42466    84932 17.041  7    0.01714 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod4<-getFE(H5.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.26       0.03
## engagement.lvl1                                             0.07       0.02
## environ.lvl1:engagement.lvl1                                0.03       0.02
##                                                               df t.value     p
## (Intercept)                                             18774.52    0.32 0.753
## age                                                     25449.11   -3.17 0.002
## gender                                                  26708.36    6.47 0.000
## educ                                                    26708.83   11.47 0.000
## resid                                                   26801.12   -6.93 0.000
## occupClerical support workers                           26637.09   -0.22 0.826
## occupCraft and related trades workers                   26646.05   -1.16 0.248
## occupElementary occupations                             26659.29    0.19 0.851
## occupManagers                                           26642.17    0.42 0.671
## occupOther: Not in paid work                            26752.45    1.60 0.110
## occupPlant and machine operators, and assemblers        26656.13   -0.62 0.534
## occupProfessionals                                      26629.25    1.53 0.125
## occupRetired                                            26643.03    0.39 0.700
## occupService and sales workers                          26644.51   -0.52 0.607
## occupSkilled agricultural, forestry and fishery workers 26661.85   -0.11 0.914
## occupTechnicians and associate professionals            26631.73   -0.10 0.923
## occupUnemployed                                         26687.13    0.10 0.918
## environ.lvl1                                               15.91    8.13 0.000
## engagement.lvl1                                            17.83    4.35 0.000
## environ.lvl1:engagement.lvl1                               17.53    1.70 0.107
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.19  0.33
## engagement.lvl1                                          0.04  0.10
## environ.lvl1:engagement.lvl1                            -0.01  0.07
```

```r
(VC.H5.mod4<-getVC(H5.mod4))
```

```
##             grp                         var1                         var2
## 1  voting.group                  (Intercept)                         <NA>
## 2  voting.group                 environ.lvl1                         <NA>
## 3  voting.group              engagement.lvl1                         <NA>
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>
## 5  voting.group                  (Intercept)                 environ.lvl1
## 6  voting.group                  (Intercept)              engagement.lvl1
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1
## 8  voting.group                 environ.lvl1              engagement.lvl1
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1
## 11        cntry                 environ.lvl1                         <NA>
## 12        cntry              engagement.lvl1                         <NA>
## 13        cntry environ.lvl1:engagement.lvl1                         <NA>
## 14        cntry                 environ.lvl1              engagement.lvl1
## 15        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 16        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 17     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.41677632  0.1737025050
## 2   0.08897592  0.0079167142
## 3   0.07541843  0.0056879391
## 4   0.06406261  0.0041040174
## 5   0.60079275  0.0222792316
## 6   0.28638770  0.0090019142
## 7  -0.02160911 -0.0005769585
## 8   0.05916292  0.0003970082
## 9  -0.77219661 -0.0044015432
## 10 -0.16805275 -0.0008119470
## 11  0.11620809  0.0135043210
## 12  0.04286872  0.0018377268
## 13  0.05024416  0.0025244758
## 14  0.19635384  0.0009781743
## 15  0.92099401  0.0053774798
## 16  0.56283359  0.0012122888
## 17  1.15882977  1.3428864343
```

```r
theta <- getME(H5.mod4,"theta")

## diagonal elements are identifiable because they are fitted
##  with a lower bound of zero ...
diag.element <- getME(H5.mod4,"lower")==0
any(theta[diag.element]<1e-5)
```

```
## [1] TRUE
```

```r
round(theta,5)
```

```
##                                  voting.group.(Intercept) 
##                                                   0.35965 
##                     voting.group.environ.lvl1.(Intercept) 
##                                                   0.04613 
##                  voting.group.engagement.lvl1.(Intercept) 
##                                                   0.01864 
##     voting.group.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                  -0.00119 
##                                 voting.group.environ.lvl1 
##                                                   0.06138 
##                 voting.group.engagement.lvl1.environ.lvl1 
##                                                  -0.00919 
##    voting.group.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                  -0.05250 
##                              voting.group.engagement.lvl1 
##                                                   0.06167 
## voting.group.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                  -0.01727 
##                 voting.group.environ.lvl1:engagement.lvl1 
##                                                   0.00000 
##                                        cntry.environ.lvl1 
##                                                   0.10028 
##                        cntry.engagement.lvl1.environ.lvl1 
##                                                   0.00726 
##           cntry.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                   0.03993 
##                                     cntry.engagement.lvl1 
##                                                   0.03627 
##        cntry.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                   0.01689 
##                        cntry.environ.lvl1:engagement.lvl1 
##                                                   0.00000
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and engagement


```r
H5.mod5<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|voting.group)+
                (0+environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                environ.lvl1:engagement.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                engagement.lvl1:all.parties.lvl2
                ,data=dat,REML=F,
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
isSingular(H5.mod5)
```

```
## [1] TRUE
```

```r
anova(H5.mod4,H5.mod5)
```

```
## Data: dat
## Models:
## H5.mod4: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod4:     engagement.lvl1 + environ.lvl1:engagement.lvl1
## H5.mod5: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H5.mod4   37 85006 85309 -42466    84932                        
## H5.mod5   64 84841 85366 -42356    84713 219.3 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod5<-getFE(H5.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.76       0.12
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
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
## occupTechnicians and associate professionals                0.00       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.15       0.05
## engagement.lvl1                                             0.02       0.04
## all.parties.lvl2Did not vote                                0.62       0.09
## all.parties.lvl2Don't know                                  0.59       0.10
## all.parties.lvl2Invalid vote                                0.61       0.52
## all.parties.lvl2NE age                                      1.00       0.10
## all.parties.lvl2NE citizen                                  1.20       0.10
## all.parties.lvl2NE other                                    1.00       0.13
## all.parties.lvl2No answer                                   1.05       0.51
## all.parties.lvl2Other party                                 0.79       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:engagement.lvl1                                0.03       0.02
## environ.lvl1:all.parties.lvl2Did not vote                   0.05       0.04
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   0.45       0.86
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.06
## environ.lvl1:all.parties.lvl2NE citizen                     0.06       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.03       0.13
## environ.lvl1:all.parties.lvl2No answer                      0.22       0.61
## environ.lvl1:all.parties.lvl2Other party                    0.14       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.12       0.06
## engagement.lvl1:all.parties.lvl2Did not vote                0.02       0.04
## engagement.lvl1:all.parties.lvl2Don't know                  0.00       0.06
## engagement.lvl1:all.parties.lvl2Invalid vote                0.67       0.86
## engagement.lvl1:all.parties.lvl2NE age                     -0.01       0.06
## engagement.lvl1:all.parties.lvl2NE citizen                 -0.05       0.06
## engagement.lvl1:all.parties.lvl2NE other                    0.00       0.12
## engagement.lvl1:all.parties.lvl2No answer                  -0.72       1.32
## engagement.lvl1:all.parties.lvl2Other party                 0.08       0.04
## engagement.lvl1:all.parties.lvl2Pro-environment party       0.12       0.05
##                                                               df t.value     p
## (Intercept)                                              3137.48   -6.30 0.000
## age                                                     26443.90   -2.79 0.005
## gender                                                  26761.32    6.45 0.000
## educ                                                    25701.39   11.43 0.000
## resid                                                   26800.08   -6.71 0.000
## occupClerical support workers                           26713.52   -0.17 0.865
## occupCraft and related trades workers                   26723.10   -1.10 0.272
## occupElementary occupations                             26735.86    0.22 0.823
## occupManagers                                           26717.48    0.47 0.639
## occupOther: Not in paid work                            26768.05    1.57 0.116
## occupPlant and machine operators, and assemblers        26738.06   -0.59 0.554
## occupProfessionals                                      26702.11    1.57 0.117
## occupRetired                                            26713.14    0.43 0.665
## occupService and sales workers                          26718.64   -0.47 0.639
## occupSkilled agricultural, forestry and fishery workers 26739.97   -0.08 0.934
## occupTechnicians and associate professionals            26708.61   -0.04 0.969
## occupUnemployed                                         26742.87    0.10 0.924
## environ.lvl1                                               55.65    3.17 0.003
## engagement.lvl1                                           127.10    0.53 0.598
## all.parties.lvl2Did not vote                              159.46    7.04 0.000
## all.parties.lvl2Don't know                                236.21    6.04 0.000
## all.parties.lvl2Invalid vote                             1725.56    1.17 0.243
## all.parties.lvl2NE age                                    243.53   10.26 0.000
## all.parties.lvl2NE citizen                                235.50   12.08 0.000
## all.parties.lvl2NE other                                  584.14    7.55 0.000
## all.parties.lvl2No answer                                1659.21    2.06 0.039
## all.parties.lvl2Other party                               210.19   12.02 0.000
## all.parties.lvl2Pro-environment party                     230.70   15.28 0.000
## environ.lvl1:engagement.lvl1                               17.56    1.75 0.098
## environ.lvl1:all.parties.lvl2Did not vote                  88.95    1.14 0.255
## environ.lvl1:all.parties.lvl2Don't know                   342.85    2.10 0.037
## environ.lvl1:all.parties.lvl2Invalid vote               23801.15    0.52 0.603
## environ.lvl1:all.parties.lvl2NE age                       305.65    3.29 0.001
## environ.lvl1:all.parties.lvl2NE citizen                   230.16    1.02 0.311
## environ.lvl1:all.parties.lvl2NE other                    1181.84    0.26 0.797
## environ.lvl1:all.parties.lvl2No answer                  14312.84    0.36 0.717
## environ.lvl1:all.parties.lvl2Other party                  134.23    3.60 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        255.42    2.13 0.034
## engagement.lvl1:all.parties.lvl2Did not vote               75.65    0.43 0.666
## engagement.lvl1:all.parties.lvl2Don't know                344.49   -0.06 0.956
## engagement.lvl1:all.parties.lvl2Invalid vote            23307.11    0.78 0.437
## engagement.lvl1:all.parties.lvl2NE age                    245.39   -0.19 0.849
## engagement.lvl1:all.parties.lvl2NE citizen                189.62   -0.90 0.368
## engagement.lvl1:all.parties.lvl2NE other                 1952.82    0.01 0.993
## engagement.lvl1:all.parties.lvl2No answer               26219.60   -0.55 0.586
## engagement.lvl1:all.parties.lvl2Other party               120.44    2.15 0.034
## engagement.lvl1:all.parties.lvl2Pro-environment party     209.65    2.25 0.025
##                                                            LL    UL
## (Intercept)                                             -1.00 -0.52
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.23  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.24
## occupManagers                                           -0.16  0.26
## occupOther: Not in paid work                            -0.04  0.39
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.22
## occupTechnicians and associate professionals            -0.21  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.05  0.24
## engagement.lvl1                                         -0.05  0.09
## all.parties.lvl2Did not vote                             0.45  0.80
## all.parties.lvl2Don't know                               0.39  0.78
## all.parties.lvl2Invalid vote                            -0.42  1.64
## all.parties.lvl2NE age                                   0.81  1.19
## all.parties.lvl2NE citizen                               1.00  1.39
## all.parties.lvl2NE other                                 0.74  1.26
## all.parties.lvl2No answer                                0.05  2.05
## all.parties.lvl2Other party                              0.66  0.92
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:engagement.lvl1                            -0.01  0.07
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.14
## environ.lvl1:all.parties.lvl2Don't know                  0.01  0.28
## environ.lvl1:all.parties.lvl2Invalid vote               -1.24  2.13
## environ.lvl1:all.parties.lvl2NE age                      0.09  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.06  0.19
## environ.lvl1:all.parties.lvl2NE other                   -0.22  0.28
## environ.lvl1:all.parties.lvl2No answer                  -0.97  1.42
## environ.lvl1:all.parties.lvl2Other party                 0.06  0.22
## environ.lvl1:all.parties.lvl2Pro-environment party       0.01  0.24
## engagement.lvl1:all.parties.lvl2Did not vote            -0.07  0.10
## engagement.lvl1:all.parties.lvl2Don't know              -0.13  0.12
## engagement.lvl1:all.parties.lvl2Invalid vote            -1.02  2.36
## engagement.lvl1:all.parties.lvl2NE age                  -0.13  0.11
## engagement.lvl1:all.parties.lvl2NE citizen              -0.17  0.06
## engagement.lvl1:all.parties.lvl2NE other                -0.23  0.24
## engagement.lvl1:all.parties.lvl2No answer               -3.31  1.87
## engagement.lvl1:all.parties.lvl2Other party              0.01  0.15
## engagement.lvl1:all.parties.lvl2Pro-environment party    0.02  0.22
```

```r
(VC.H5.mod5<-getVC(H5.mod5))
```

```
##             grp                         var1                         var2
## 1  voting.group                  (Intercept)                         <NA>
## 2  voting.group                 environ.lvl1                         <NA>
## 3  voting.group              engagement.lvl1                         <NA>
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>
## 5  voting.group                  (Intercept)                 environ.lvl1
## 6  voting.group                  (Intercept)              engagement.lvl1
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1
## 8  voting.group                 environ.lvl1              engagement.lvl1
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1
## 11        cntry                 environ.lvl1                         <NA>
## 12        cntry              engagement.lvl1                         <NA>
## 13        cntry environ.lvl1:engagement.lvl1                         <NA>
## 14        cntry                 environ.lvl1              engagement.lvl1
## 15        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 16        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 17     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.25638896  0.0657352997
## 2   0.07500612  0.0056259178
## 3   0.06061981  0.0036747619
## 4   0.05245696  0.0027517322
## 5   0.56928938  0.0109478566
## 6   0.23734322  0.0036888479
## 7  -0.58225751 -0.0078310051
## 8  -0.28635412 -0.0013020112
## 9  -0.99726128 -0.0039238169
## 10  0.21770582  0.0006922894
## 11  0.11762307  0.0138351872
## 12  0.04851450  0.0023536567
## 13  0.05047429  0.0025476537
## 14  0.18990027  0.0010836516
## 15  0.88263005  0.0052401223
## 16  0.62912632  0.0015405635
## 17  1.15880374  1.3428261063
```

\newpage

#### Look among which voting group there is strongest association between engagement and refugee attitudes


```r
H5.mod5.trends<-emtrends(H5.mod5,specs = c("all.parties.lvl2"),var=c("engagement.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 26886' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 26886)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 26886' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 26886)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.mod5.trends.tab<-data.frame(H5.mod5.trends))
```

```
##          all.parties.lvl2 engagement.lvl1.trend         SE  df   asymp.LCL
## 1  Anti-immigration party            0.01910620 0.03583922 Inf -0.05113738
## 2            Did not vote            0.03774589 0.02973676 Inf -0.02053708
## 3              Don't know            0.01558112 0.05577807 Inf -0.09374189
## 4            Invalid vote            0.68913879 0.86158073 Inf -0.99952842
## 5                  NE age            0.00782722 0.05046247 Inf -0.09107740
## 6              NE citizen           -0.03448363 0.05054434 Inf -0.13354871
## 7                NE other            0.02021220 0.11612225 Inf -0.20738322
## 8               No answer           -0.70011699 1.31900325 Inf -3.28531586
## 9             Other party            0.09908194 0.02037265 Inf  0.05915229
## 10  Pro-environment party            0.13878791 0.04305092 Inf  0.05440965
##     asymp.UCL
## 1  0.08934979
## 2  0.09602887
## 3  0.12490413
## 4  2.37780600
## 5  0.10673184
## 6  0.06458145
## 7  0.24780762
## 8  1.88508188
## 9  0.13901160
## 10 0.22316616
```

```r
H5.mod5.trends.tab$p<-
  2*(1-pnorm(abs(H5.mod5.trends.tab$engagement.lvl1.trend/
                   H5.mod5.trends.tab$SE)))
H5.mod5.trends.tab$adj.p<-
  p.adjust(H5.mod5.trends.tab$p,method="holm")

H5.mod5.trends.tab<-
  cbind(group=H5.mod5.trends.tab[,1],
      round(H5.mod5.trends.tab[,c(2,3)],2),
      round(H5.mod5.trends.tab[,c(7,8)],4),
      round(H5.mod5.trends.tab[,c(5,6)],2))
H5.mod5.trends.tab
```

```
##                     group engagement.lvl1.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party                  0.02 0.04 0.5940 1.0000     -0.05
## 2            Did not vote                  0.04 0.03 0.2043 1.0000     -0.02
## 3              Don't know                  0.02 0.06 0.7800 1.0000     -0.09
## 4            Invalid vote                  0.69 0.86 0.4238 1.0000     -1.00
## 5                  NE age                  0.01 0.05 0.8767 1.0000     -0.09
## 6              NE citizen                 -0.03 0.05 0.4951 1.0000     -0.13
## 7                NE other                  0.02 0.12 0.8618 1.0000     -0.21
## 8               No answer                 -0.70 1.32 0.5956 1.0000     -3.29
## 9             Other party                  0.10 0.02 0.0000 0.0000      0.06
## 10  Pro-environment party                  0.14 0.04 0.0013 0.0114      0.05
##    asymp.UCL
## 1       0.09
## 2       0.10
## 3       0.12
## 4       2.38
## 5       0.11
## 6       0.06
## 7       0.25
## 8       1.89
## 9       0.14
## 10      0.22
```

```r
write.csv2(H5.mod5.trends.tab,"H5.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.1196817 0.05308428 Inf
##   z.ratio    p.value
## 1 2.25456 0.02416096
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect  -0.0113 0.179 Inf -0.063  0.9496 
##  Did not vote effect             0.0094 0.178 Inf  0.053  0.9579 
##  Don't know effect              -0.0152 0.184 Inf -0.083  0.9341 
##  Invalid vote effect             0.7332 0.874 Inf  0.839  0.4016 
##  NE age effect                  -0.0238 0.183 Inf -0.131  0.8960 
##  NE citizen effect              -0.0709 0.183 Inf -0.388  0.6979 
##  NE other effect                -0.0101 0.210 Inf -0.048  0.9617 
##  No answer effect               -0.8105 1.323 Inf -0.613  0.5400 
##  Other party effect              0.0775 0.177 Inf  0.439  0.6606 
##  Pro-environment party effect    0.1217 0.181 Inf  0.673  0.5006 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic
```

```r
#contrast for three voting groups
(H4.more.contrasts<-data.frame(pairs(H5.mod5.trends, 
         exclude=c(2:8), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1           Other party - Anti-immigration party 0.07997574 0.03722896 Inf
## 2 Pro-environment party - Anti-immigration party 0.11968170 0.05308428 Inf
## 3            Pro-environment party - Other party 0.03970596 0.04417894 Inf
##     z.ratio    p.value
## 1 2.1482132 0.03169682
## 2 2.2545599 0.02416096
## 3 0.8987532 0.36878413
```

\newpage

### Model 6: Enter three-way interaction voting group x engagement x environment attitudes


```r
H5.mod6<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|voting.group)+
                (0+environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                environ.lvl1:engagement.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                engagement.lvl1:all.parties.lvl2+
                environ.lvl1:engagement.lvl1:all.parties.lvl2
                ,data=dat,REML=F,
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
isSingular(H5.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.mod5,H5.mod6)
```

```
## Data: dat
## Models:
## H5.mod5: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
## H5.mod6: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod6:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod6:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod6:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2 + 
## H5.mod6:     environ.lvl1:engagement.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   64 84841 85366 -42356    84713                       
## H5.mod6   73 84839 85438 -42347    84693 19.567  9    0.02078 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                                    Estimate
## (Intercept)                                                           -0.76
## age                                                                    0.00
## gender                                                                 0.10
## educ                                                                   0.03
## resid                                                                 -0.10
## occupClerical support workers                                         -0.02
## occupCraft and related trades workers                                 -0.12
## occupElementary occupations                                            0.02
## occupManagers                                                          0.05
## occupOther: Not in paid work                                           0.18
## occupPlant and machine operators, and assemblers                      -0.07
## occupProfessionals                                                     0.17
## occupRetired                                                           0.05
## occupService and sales workers                                        -0.05
## occupSkilled agricultural, forestry and fishery workers               -0.01
## occupTechnicians and associate professionals                          -0.01
## occupUnemployed                                                        0.01
## environ.lvl1                                                           0.15
## engagement.lvl1                                                        0.02
## all.parties.lvl2Did not vote                                           0.62
## all.parties.lvl2Don't know                                             0.58
## all.parties.lvl2Invalid vote                                           0.29
## all.parties.lvl2NE age                                                 0.99
## all.parties.lvl2NE citizen                                             1.14
## all.parties.lvl2NE other                                               1.01
## all.parties.lvl2No answer                                              0.90
## all.parties.lvl2Other party                                            0.80
## all.parties.lvl2Pro-environment party                                  1.31
## environ.lvl1:engagement.lvl1                                           0.02
## environ.lvl1:all.parties.lvl2Did not vote                              0.05
## environ.lvl1:all.parties.lvl2Don't know                                0.14
## environ.lvl1:all.parties.lvl2Invalid vote                              1.15
## environ.lvl1:all.parties.lvl2NE age                                    0.21
## environ.lvl1:all.parties.lvl2NE citizen                                0.05
## environ.lvl1:all.parties.lvl2NE other                                  0.03
## environ.lvl1:all.parties.lvl2No answer                                -0.04
## environ.lvl1:all.parties.lvl2Other party                               0.14
## environ.lvl1:all.parties.lvl2Pro-environment party                     0.12
## engagement.lvl1:all.parties.lvl2Did not vote                           0.02
## engagement.lvl1:all.parties.lvl2Don't know                             0.00
## engagement.lvl1:all.parties.lvl2Invalid vote                           1.07
## engagement.lvl1:all.parties.lvl2NE age                                -0.01
## engagement.lvl1:all.parties.lvl2NE citizen                            -0.06
## engagement.lvl1:all.parties.lvl2NE other                               0.02
## engagement.lvl1:all.parties.lvl2No answer                             -1.08
## engagement.lvl1:all.parties.lvl2Other party                            0.08
## engagement.lvl1:all.parties.lvl2Pro-environment party                  0.12
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote              0.01
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know                0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote             -1.43
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                    0.05
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen                0.26
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                 -0.09
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer                -1.15
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party              -0.01
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party    -0.01
##                                                                    Std..Error
## (Intercept)                                                              0.12
## age                                                                      0.00
## gender                                                                   0.02
## educ                                                                     0.00
## resid                                                                    0.02
## occupClerical support workers                                            0.11
## occupCraft and related trades workers                                    0.11
## occupElementary occupations                                              0.11
## occupManagers                                                            0.11
## occupOther: Not in paid work                                             0.11
## occupPlant and machine operators, and assemblers                         0.11
## occupProfessionals                                                       0.11
## occupRetired                                                             0.12
## occupService and sales workers                                           0.11
## occupSkilled agricultural, forestry and fishery workers                  0.12
## occupTechnicians and associate professionals                             0.11
## occupUnemployed                                                          0.13
## environ.lvl1                                                             0.05
## engagement.lvl1                                                          0.04
## all.parties.lvl2Did not vote                                             0.09
## all.parties.lvl2Don't know                                               0.10
## all.parties.lvl2Invalid vote                                             0.67
## all.parties.lvl2NE age                                                   0.10
## all.parties.lvl2NE citizen                                               0.10
## all.parties.lvl2NE other                                                 0.13
## all.parties.lvl2No answer                                                0.59
## all.parties.lvl2Other party                                              0.07
## all.parties.lvl2Pro-environment party                                    0.09
## environ.lvl1:engagement.lvl1                                             0.04
## environ.lvl1:all.parties.lvl2Did not vote                                0.05
## environ.lvl1:all.parties.lvl2Don't know                                  0.07
## environ.lvl1:all.parties.lvl2Invalid vote                                1.27
## environ.lvl1:all.parties.lvl2NE age                                      0.07
## environ.lvl1:all.parties.lvl2NE citizen                                  0.06
## environ.lvl1:all.parties.lvl2NE other                                    0.13
## environ.lvl1:all.parties.lvl2No answer                                   0.81
## environ.lvl1:all.parties.lvl2Other party                                 0.04
## environ.lvl1:all.parties.lvl2Pro-environment party                       0.06
## engagement.lvl1:all.parties.lvl2Did not vote                             0.04
## engagement.lvl1:all.parties.lvl2Don't know                               0.06
## engagement.lvl1:all.parties.lvl2Invalid vote                             1.01
## engagement.lvl1:all.parties.lvl2NE age                                   0.06
## engagement.lvl1:all.parties.lvl2NE citizen                               0.06
## engagement.lvl1:all.parties.lvl2NE other                                 0.12
## engagement.lvl1:all.parties.lvl2No answer                                1.51
## engagement.lvl1:all.parties.lvl2Other party                              0.04
## engagement.lvl1:all.parties.lvl2Pro-environment party                    0.05
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote                0.05
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know                  0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote                1.91
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                      0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen                  0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                    0.16
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer                   2.38
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party                 0.04
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party       0.07
##                                                                          df
## (Intercept)                                                         3077.29
## age                                                                26481.26
## gender                                                             26761.57
## educ                                                               25806.59
## resid                                                              26805.66
## occupClerical support workers                                      26723.71
## occupCraft and related trades workers                              26732.63
## occupElementary occupations                                        26743.73
## occupManagers                                                      26726.53
## occupOther: Not in paid work                                       26774.68
## occupPlant and machine operators, and assemblers                   26744.32
## occupProfessionals                                                 26712.08
## occupRetired                                                       26721.81
## occupService and sales workers                                     26727.25
## occupSkilled agricultural, forestry and fishery workers            26746.50
## occupTechnicians and associate professionals                       26718.84
## occupUnemployed                                                    26748.35
## environ.lvl1                                                          58.06
## engagement.lvl1                                                      126.56
## all.parties.lvl2Did not vote                                         175.50
## all.parties.lvl2Don't know                                           233.96
## all.parties.lvl2Invalid vote                                        4947.03
## all.parties.lvl2NE age                                               245.54
## all.parties.lvl2NE citizen                                           238.15
## all.parties.lvl2NE other                                             583.35
## all.parties.lvl2No answer                                           3120.04
## all.parties.lvl2Other party                                          213.92
## all.parties.lvl2Pro-environment party                                232.51
## environ.lvl1:engagement.lvl1                                         170.31
## environ.lvl1:all.parties.lvl2Did not vote                             86.68
## environ.lvl1:all.parties.lvl2Don't know                              331.76
## environ.lvl1:all.parties.lvl2Invalid vote                          26025.00
## environ.lvl1:all.parties.lvl2NE age                                  291.19
## environ.lvl1:all.parties.lvl2NE citizen                              223.73
## environ.lvl1:all.parties.lvl2NE other                               1197.41
## environ.lvl1:all.parties.lvl2No answer                             22679.08
## environ.lvl1:all.parties.lvl2Other party                             134.18
## environ.lvl1:all.parties.lvl2Pro-environment party                   246.85
## engagement.lvl1:all.parties.lvl2Did not vote                          75.57
## engagement.lvl1:all.parties.lvl2Don't know                           344.77
## engagement.lvl1:all.parties.lvl2Invalid vote                       25500.84
## engagement.lvl1:all.parties.lvl2NE age                               245.48
## engagement.lvl1:all.parties.lvl2NE citizen                           190.97
## engagement.lvl1:all.parties.lvl2NE other                            2345.18
## engagement.lvl1:all.parties.lvl2No answer                          26430.44
## engagement.lvl1:all.parties.lvl2Other party                          120.60
## engagement.lvl1:all.parties.lvl2Pro-environment party                209.36
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote            154.04
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know              959.38
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote          26355.15
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                  825.11
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen              475.43
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other               3966.26
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer             26379.11
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party             249.64
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party   668.25
##                                                                    t.value
## (Intercept)                                                          -6.25
## age                                                                  -2.76
## gender                                                                6.47
## educ                                                                 11.43
## resid                                                                -6.67
## occupClerical support workers                                        -0.18
## occupCraft and related trades workers                                -1.11
## occupElementary occupations                                           0.21
## occupManagers                                                         0.45
## occupOther: Not in paid work                                          1.58
## occupPlant and machine operators, and assemblers                     -0.60
## occupProfessionals                                                    1.56
## occupRetired                                                          0.42
## occupService and sales workers                                       -0.47
## occupSkilled agricultural, forestry and fishery workers              -0.08
## occupTechnicians and associate professionals                         -0.05
## occupUnemployed                                                       0.08
## environ.lvl1                                                          3.16
## engagement.lvl1                                                       0.51
## all.parties.lvl2Did not vote                                          6.83
## all.parties.lvl2Don't know                                            5.95
## all.parties.lvl2Invalid vote                                          0.43
## all.parties.lvl2NE age                                               10.02
## all.parties.lvl2NE citizen                                           11.36
## all.parties.lvl2NE other                                              7.54
## all.parties.lvl2No answer                                             1.52
## all.parties.lvl2Other party                                          11.91
## all.parties.lvl2Pro-environment party                                15.11
## environ.lvl1:engagement.lvl1                                          0.48
## environ.lvl1:all.parties.lvl2Did not vote                             1.04
## environ.lvl1:all.parties.lvl2Don't know                               2.04
## environ.lvl1:all.parties.lvl2Invalid vote                             0.91
## environ.lvl1:all.parties.lvl2NE age                                   3.20
## environ.lvl1:all.parties.lvl2NE citizen                               0.80
## environ.lvl1:all.parties.lvl2NE other                                 0.23
## environ.lvl1:all.parties.lvl2No answer                               -0.06
## environ.lvl1:all.parties.lvl2Other party                              3.58
## environ.lvl1:all.parties.lvl2Pro-environment party                    2.12
## engagement.lvl1:all.parties.lvl2Did not vote                          0.46
## engagement.lvl1:all.parties.lvl2Don't know                           -0.04
## engagement.lvl1:all.parties.lvl2Invalid vote                          1.06
## engagement.lvl1:all.parties.lvl2NE age                               -0.17
## engagement.lvl1:all.parties.lvl2NE citizen                           -0.98
## engagement.lvl1:all.parties.lvl2NE other                              0.12
## engagement.lvl1:all.parties.lvl2No answer                            -0.71
## engagement.lvl1:all.parties.lvl2Other party                           2.13
## engagement.lvl1:all.parties.lvl2Pro-environment party                 2.24
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote             0.29
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know               0.19
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote            -0.75
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                   0.57
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen               3.59
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                -0.54
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer               -0.48
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party             -0.17
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party   -0.19
##                                                                        p    LL
## (Intercept)                                                        0.000 -0.99
## age                                                                0.006  0.00
## gender                                                             0.000  0.07
## educ                                                               0.000  0.02
## resid                                                              0.000 -0.13
## occupClerical support workers                                      0.861 -0.23
## occupCraft and related trades workers                              0.267 -0.33
## occupElementary occupations                                        0.831 -0.19
## occupManagers                                                      0.651 -0.16
## occupOther: Not in paid work                                       0.114 -0.04
## occupPlant and machine operators, and assemblers                   0.549 -0.28
## occupProfessionals                                                 0.119 -0.04
## occupRetired                                                       0.673 -0.19
## occupService and sales workers                                     0.637 -0.26
## occupSkilled agricultural, forestry and fishery workers            0.938 -0.23
## occupTechnicians and associate professionals                       0.957 -0.21
## occupUnemployed                                                    0.933 -0.25
## environ.lvl1                                                       0.003  0.05
## engagement.lvl1                                                    0.613 -0.05
## all.parties.lvl2Did not vote                                       0.000  0.44
## all.parties.lvl2Don't know                                         0.000  0.39
## all.parties.lvl2Invalid vote                                       0.670 -1.03
## all.parties.lvl2NE age                                             0.000  0.80
## all.parties.lvl2NE citizen                                         0.000  0.94
## all.parties.lvl2NE other                                           0.000  0.75
## all.parties.lvl2No answer                                          0.129 -0.26
## all.parties.lvl2Other party                                        0.000  0.66
## all.parties.lvl2Pro-environment party                              0.000  1.14
## environ.lvl1:engagement.lvl1                                       0.632 -0.06
## environ.lvl1:all.parties.lvl2Did not vote                          0.299 -0.04
## environ.lvl1:all.parties.lvl2Don't know                            0.042  0.01
## environ.lvl1:all.parties.lvl2Invalid vote                          0.365 -1.34
## environ.lvl1:all.parties.lvl2NE age                                0.002  0.08
## environ.lvl1:all.parties.lvl2NE citizen                            0.424 -0.07
## environ.lvl1:all.parties.lvl2NE other                              0.819 -0.22
## environ.lvl1:all.parties.lvl2No answer                             0.956 -1.64
## environ.lvl1:all.parties.lvl2Other party                           0.000  0.06
## environ.lvl1:all.parties.lvl2Pro-environment party                 0.035  0.01
## engagement.lvl1:all.parties.lvl2Did not vote                       0.649 -0.07
## engagement.lvl1:all.parties.lvl2Don't know                         0.964 -0.13
## engagement.lvl1:all.parties.lvl2Invalid vote                       0.290 -0.91
## engagement.lvl1:all.parties.lvl2NE age                             0.862 -0.13
## engagement.lvl1:all.parties.lvl2NE citizen                         0.326 -0.18
## engagement.lvl1:all.parties.lvl2NE other                           0.902 -0.22
## engagement.lvl1:all.parties.lvl2No answer                          0.475 -4.05
## engagement.lvl1:all.parties.lvl2Other party                        0.035  0.01
## engagement.lvl1:all.parties.lvl2Pro-environment party              0.026  0.01
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote          0.772 -0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know            0.851 -0.14
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote          0.454 -5.16
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                0.566 -0.11
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen            0.000  0.12
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other              0.591 -0.41
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer             0.629 -5.81
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party           0.861 -0.09
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party 0.848 -0.15
##                                                                       UL
## (Intercept)                                                        -0.52
## age                                                                 0.00
## gender                                                              0.13
## educ                                                                0.03
## resid                                                              -0.07
## occupClerical support workers                                       0.19
## occupCraft and related trades workers                               0.09
## occupElementary occupations                                         0.24
## occupManagers                                                       0.26
## occupOther: Not in paid work                                        0.40
## occupPlant and machine operators, and assemblers                    0.15
## occupProfessionals                                                  0.37
## occupRetired                                                        0.29
## occupService and sales workers                                      0.16
## occupSkilled agricultural, forestry and fishery workers             0.22
## occupTechnicians and associate professionals                        0.20
## occupUnemployed                                                     0.27
## environ.lvl1                                                        0.24
## engagement.lvl1                                                     0.09
## all.parties.lvl2Did not vote                                        0.80
## all.parties.lvl2Don't know                                          0.78
## all.parties.lvl2Invalid vote                                        1.61
## all.parties.lvl2NE age                                              1.19
## all.parties.lvl2NE citizen                                          1.34
## all.parties.lvl2NE other                                            1.27
## all.parties.lvl2No answer                                           2.06
## all.parties.lvl2Other party                                         0.93
## all.parties.lvl2Pro-environment party                               1.48
## environ.lvl1:engagement.lvl1                                        0.10
## environ.lvl1:all.parties.lvl2Did not vote                           0.14
## environ.lvl1:all.parties.lvl2Don't know                             0.28
## environ.lvl1:all.parties.lvl2Invalid vote                           3.64
## environ.lvl1:all.parties.lvl2NE age                                 0.34
## environ.lvl1:all.parties.lvl2NE citizen                             0.18
## environ.lvl1:all.parties.lvl2NE other                               0.28
## environ.lvl1:all.parties.lvl2No answer                              1.55
## environ.lvl1:all.parties.lvl2Other party                            0.22
## environ.lvl1:all.parties.lvl2Pro-environment party                  0.24
## engagement.lvl1:all.parties.lvl2Did not vote                        0.11
## engagement.lvl1:all.parties.lvl2Don't know                          0.12
## engagement.lvl1:all.parties.lvl2Invalid vote                        3.05
## engagement.lvl1:all.parties.lvl2NE age                              0.11
## engagement.lvl1:all.parties.lvl2NE citizen                          0.06
## engagement.lvl1:all.parties.lvl2NE other                            0.25
## engagement.lvl1:all.parties.lvl2No answer                           1.88
## engagement.lvl1:all.parties.lvl2Other party                         0.15
## engagement.lvl1:all.parties.lvl2Pro-environment party               0.22
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote           0.11
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know             0.17
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote           2.31
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                 0.20
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen             0.40
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other               0.23
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer              3.51
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party            0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party  0.12
```

```r
(VC.H5.mod6<-getVC(H5.mod6))
```

```
##             grp                         var1                         var2
## 1  voting.group                  (Intercept)                         <NA>
## 2  voting.group                 environ.lvl1                         <NA>
## 3  voting.group              engagement.lvl1                         <NA>
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>
## 5  voting.group                  (Intercept)                 environ.lvl1
## 6  voting.group                  (Intercept)              engagement.lvl1
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1
## 8  voting.group                 environ.lvl1              engagement.lvl1
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1
## 11        cntry                 environ.lvl1                         <NA>
## 12        cntry              engagement.lvl1                         <NA>
## 13        cntry environ.lvl1:engagement.lvl1                         <NA>
## 14        cntry                 environ.lvl1              engagement.lvl1
## 15        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 16        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 17     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.25657750  0.0658320148
## 2   0.07513531  0.0056453141
## 3   0.05998835  0.0035986023
## 4   0.04505391  0.0020298552
## 5   0.55882005  0.0107729491
## 6   0.24743023  0.0038083623
## 7  -0.62603675 -0.0072368725
## 8  -0.26471943 -0.0011931548
## 9  -0.98622016 -0.0033384929
## 10  0.10182750  0.0002752102
## 11  0.11757623  0.0138241699
## 12  0.04889703  0.0023909194
## 13  0.05202864  0.0027069798
## 14  0.21446415  0.0012329819
## 15  0.90015092  0.0055065219
## 16  0.61849316  0.0015734751
## 17  1.15841886  1.3419342610
```

#### Refit with manually coded level-1 interaction


```r
dat$env.eng.int<-dat$environ.lvl1*dat$engagement.lvl1

H5.mod6<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+env.eng.int|voting.group)+
                (0+environ.lvl1+engagement.lvl1+env.eng.int|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                engagement.lvl1+
                env.eng.int+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                engagement.lvl1:all.parties.lvl2+
                env.eng.int:all.parties.lvl2
                ,data=dat,REML=F,
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
isSingular(H5.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.mod5,H5.mod6)
```

```
## Data: dat
## Models:
## H5.mod5: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     voting.group) + (0 + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
## H5.mod6: refugees ~ (environ.lvl1 + engagement.lvl1 + env.eng.int | voting.group) + 
## H5.mod6:     (0 + environ.lvl1 + engagement.lvl1 + env.eng.int | cntry) + 
## H5.mod6:     age + gender + educ + resid + occup + environ.lvl1 + engagement.lvl1 + 
## H5.mod6:     env.eng.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.mod6:     engagement.lvl1:all.parties.lvl2 + env.eng.int:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   64 84841 85366 -42356    84713                       
## H5.mod6   73 84839 85438 -42347    84693 19.567  9    0.02078 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.76       0.12
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.15       0.05
## engagement.lvl1                                             0.02       0.04
## env.eng.int                                                 0.02       0.04
## all.parties.lvl2Did not vote                                0.62       0.09
## all.parties.lvl2Don't know                                  0.58       0.10
## all.parties.lvl2Invalid vote                                0.29       0.67
## all.parties.lvl2NE age                                      0.99       0.10
## all.parties.lvl2NE citizen                                  1.14       0.10
## all.parties.lvl2NE other                                    1.01       0.13
## all.parties.lvl2No answer                                   0.90       0.59
## all.parties.lvl2Other party                                 0.80       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:all.parties.lvl2Did not vote                   0.05       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   1.15       1.27
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.07
## environ.lvl1:all.parties.lvl2NE citizen                     0.05       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.03       0.13
## environ.lvl1:all.parties.lvl2No answer                     -0.04       0.81
## environ.lvl1:all.parties.lvl2Other party                    0.14       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.12       0.06
## engagement.lvl1:all.parties.lvl2Did not vote                0.02       0.04
## engagement.lvl1:all.parties.lvl2Don't know                  0.00       0.06
## engagement.lvl1:all.parties.lvl2Invalid vote                1.07       1.01
## engagement.lvl1:all.parties.lvl2NE age                     -0.01       0.06
## engagement.lvl1:all.parties.lvl2NE citizen                 -0.06       0.06
## engagement.lvl1:all.parties.lvl2NE other                    0.02       0.12
## engagement.lvl1:all.parties.lvl2No answer                  -1.08       1.51
## engagement.lvl1:all.parties.lvl2Other party                 0.08       0.04
## engagement.lvl1:all.parties.lvl2Pro-environment party       0.12       0.05
## env.eng.int:all.parties.lvl2Did not vote                    0.01       0.05
## env.eng.int:all.parties.lvl2Don't know                      0.02       0.08
## env.eng.int:all.parties.lvl2Invalid vote                   -1.43       1.91
## env.eng.int:all.parties.lvl2NE age                          0.05       0.08
## env.eng.int:all.parties.lvl2NE citizen                      0.26       0.07
## env.eng.int:all.parties.lvl2NE other                       -0.09       0.16
## env.eng.int:all.parties.lvl2No answer                      -1.15       2.38
## env.eng.int:all.parties.lvl2Other party                    -0.01       0.04
## env.eng.int:all.parties.lvl2Pro-environment party          -0.01       0.07
##                                                               df t.value     p
## (Intercept)                                              3077.30   -6.25 0.000
## age                                                     26481.26   -2.76 0.006
## gender                                                  26761.57    6.47 0.000
## educ                                                    25806.59   11.43 0.000
## resid                                                   26805.66   -6.67 0.000
## occupClerical support workers                           26723.71   -0.18 0.861
## occupCraft and related trades workers                   26732.63   -1.11 0.267
## occupElementary occupations                             26743.73    0.21 0.831
## occupManagers                                           26726.53    0.45 0.651
## occupOther: Not in paid work                            26774.68    1.58 0.114
## occupPlant and machine operators, and assemblers        26744.32   -0.60 0.549
## occupProfessionals                                      26712.08    1.56 0.119
## occupRetired                                            26721.81    0.42 0.673
## occupService and sales workers                          26727.25   -0.47 0.637
## occupSkilled agricultural, forestry and fishery workers 26746.50   -0.08 0.938
## occupTechnicians and associate professionals            26718.84   -0.05 0.957
## occupUnemployed                                         26748.35    0.08 0.933
## environ.lvl1                                               58.06    3.16 0.003
## engagement.lvl1                                           126.56    0.51 0.613
## env.eng.int                                               170.30    0.48 0.632
## all.parties.lvl2Did not vote                              175.50    6.83 0.000
## all.parties.lvl2Don't know                                233.96    5.95 0.000
## all.parties.lvl2Invalid vote                             4947.05    0.43 0.670
## all.parties.lvl2NE age                                    245.54   10.02 0.000
## all.parties.lvl2NE citizen                                238.15   11.36 0.000
## all.parties.lvl2NE other                                  583.35    7.54 0.000
## all.parties.lvl2No answer                                3120.05    1.52 0.129
## all.parties.lvl2Other party                               213.92   11.91 0.000
## all.parties.lvl2Pro-environment party                     232.51   15.11 0.000
## environ.lvl1:all.parties.lvl2Did not vote                  86.68    1.04 0.299
## environ.lvl1:all.parties.lvl2Don't know                   331.76    2.04 0.042
## environ.lvl1:all.parties.lvl2Invalid vote               26025.00    0.91 0.365
## environ.lvl1:all.parties.lvl2NE age                       291.19    3.20 0.002
## environ.lvl1:all.parties.lvl2NE citizen                   223.73    0.80 0.424
## environ.lvl1:all.parties.lvl2NE other                    1197.41    0.23 0.819
## environ.lvl1:all.parties.lvl2No answer                  22679.09   -0.06 0.956
## environ.lvl1:all.parties.lvl2Other party                  134.18    3.58 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        246.85    2.12 0.035
## engagement.lvl1:all.parties.lvl2Did not vote               75.57    0.46 0.649
## engagement.lvl1:all.parties.lvl2Don't know                344.77   -0.04 0.964
## engagement.lvl1:all.parties.lvl2Invalid vote            25500.84    1.06 0.290
## engagement.lvl1:all.parties.lvl2NE age                    245.48   -0.17 0.862
## engagement.lvl1:all.parties.lvl2NE citizen                190.97   -0.98 0.326
## engagement.lvl1:all.parties.lvl2NE other                 2345.18    0.12 0.902
## engagement.lvl1:all.parties.lvl2No answer               26430.44   -0.71 0.475
## engagement.lvl1:all.parties.lvl2Other party               120.60    2.13 0.035
## engagement.lvl1:all.parties.lvl2Pro-environment party     209.36    2.24 0.026
## env.eng.int:all.parties.lvl2Did not vote                  154.04    0.29 0.772
## env.eng.int:all.parties.lvl2Don't know                    959.36    0.19 0.851
## env.eng.int:all.parties.lvl2Invalid vote                26355.15   -0.75 0.454
## env.eng.int:all.parties.lvl2NE age                        825.09    0.57 0.566
## env.eng.int:all.parties.lvl2NE citizen                    475.43    3.59 0.000
## env.eng.int:all.parties.lvl2NE other                     3966.21   -0.54 0.591
## env.eng.int:all.parties.lvl2No answer                   26379.11   -0.48 0.629
## env.eng.int:all.parties.lvl2Other party                   249.64   -0.17 0.861
## env.eng.int:all.parties.lvl2Pro-environment party         668.23   -0.19 0.848
##                                                            LL    UL
## (Intercept)                                             -0.99 -0.52
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.23  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.24
## occupManagers                                           -0.16  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.15
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.23  0.22
## occupTechnicians and associate professionals            -0.21  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.05  0.24
## engagement.lvl1                                         -0.05  0.09
## env.eng.int                                             -0.06  0.10
## all.parties.lvl2Did not vote                             0.44  0.80
## all.parties.lvl2Don't know                               0.39  0.78
## all.parties.lvl2Invalid vote                            -1.03  1.61
## all.parties.lvl2NE age                                   0.80  1.19
## all.parties.lvl2NE citizen                               0.94  1.34
## all.parties.lvl2NE other                                 0.75  1.27
## all.parties.lvl2No answer                               -0.26  2.06
## all.parties.lvl2Other party                              0.66  0.93
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.14
## environ.lvl1:all.parties.lvl2Don't know                  0.01  0.28
## environ.lvl1:all.parties.lvl2Invalid vote               -1.34  3.64
## environ.lvl1:all.parties.lvl2NE age                      0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.07  0.18
## environ.lvl1:all.parties.lvl2NE other                   -0.22  0.28
## environ.lvl1:all.parties.lvl2No answer                  -1.64  1.55
## environ.lvl1:all.parties.lvl2Other party                 0.06  0.22
## environ.lvl1:all.parties.lvl2Pro-environment party       0.01  0.24
## engagement.lvl1:all.parties.lvl2Did not vote            -0.07  0.11
## engagement.lvl1:all.parties.lvl2Don't know              -0.13  0.12
## engagement.lvl1:all.parties.lvl2Invalid vote            -0.91  3.05
## engagement.lvl1:all.parties.lvl2NE age                  -0.13  0.11
## engagement.lvl1:all.parties.lvl2NE citizen              -0.18  0.06
## engagement.lvl1:all.parties.lvl2NE other                -0.22  0.25
## engagement.lvl1:all.parties.lvl2No answer               -4.05  1.88
## engagement.lvl1:all.parties.lvl2Other party              0.01  0.15
## engagement.lvl1:all.parties.lvl2Pro-environment party    0.01  0.22
## env.eng.int:all.parties.lvl2Did not vote                -0.08  0.11
## env.eng.int:all.parties.lvl2Don't know                  -0.14  0.17
## env.eng.int:all.parties.lvl2Invalid vote                -5.16  2.31
## env.eng.int:all.parties.lvl2NE age                      -0.11  0.20
## env.eng.int:all.parties.lvl2NE citizen                   0.12  0.40
## env.eng.int:all.parties.lvl2NE other                    -0.41  0.23
## env.eng.int:all.parties.lvl2No answer                   -5.81  3.51
## env.eng.int:all.parties.lvl2Other party                 -0.09  0.08
## env.eng.int:all.parties.lvl2Pro-environment party       -0.15  0.12
```

```r
(VC.H5.mod6<-getVC(H5.mod6))
```

```
##             grp            var1            var2      est_SD       est_SD2
## 1  voting.group     (Intercept)            <NA>  0.25657729  0.0658319058
## 2  voting.group    environ.lvl1            <NA>  0.07513514  0.0056452891
## 3  voting.group engagement.lvl1            <NA>  0.05998821  0.0035985859
## 4  voting.group     env.eng.int            <NA>  0.04505408  0.0020298706
## 5  voting.group     (Intercept)    environ.lvl1  0.55882314  0.0107729758
## 6  voting.group     (Intercept) engagement.lvl1  0.24742652  0.0038082934
## 7  voting.group     (Intercept)     env.eng.int -0.62603501 -0.0072368739
## 8  voting.group    environ.lvl1 engagement.lvl1 -0.26473164 -0.0011932045
## 9  voting.group    environ.lvl1     env.eng.int -0.98621665 -0.0033384863
## 10 voting.group engagement.lvl1     env.eng.int  0.10181793  0.0002751848
## 11        cntry    environ.lvl1            <NA>  0.11757609  0.0138241373
## 12        cntry engagement.lvl1            <NA>  0.04889704  0.0023909205
## 13        cntry     env.eng.int            <NA>  0.05202879  0.0027069952
## 14        cntry    environ.lvl1 engagement.lvl1  0.21446783  0.0012330019
## 15        cntry    environ.lvl1     env.eng.int  0.90015244  0.0055065403
## 16        cntry engagement.lvl1     env.eng.int  0.61849337  0.0015734805
## 17     Residual            <NA>            <NA>  1.15841888  1.3419342926
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.mod6.trends<-emtrends(H5.mod6,specs = c("all.parties.lvl2"),var=c("env.eng.int"))
(H5.mod6.trends.tab<-data.frame(H5.mod6.trends))
```

```
##          all.parties.lvl2 env.eng.int.trend         SE  df   asymp.LCL
## 1  Anti-immigration party       0.020123261 0.04198007 Inf -0.06215615
## 2            Did not vote       0.034574347 0.03274632 Inf -0.02960725
## 3              Don't know       0.035305365 0.07130445 Inf -0.10444878
## 4            Invalid vote      -1.406159530 1.90559325 Inf -5.14105368
## 5                  NE age       0.065382387 0.06907325 Inf -0.06999869
## 6              NE citizen       0.278011175 0.06122316 Inf  0.15801599
## 7                NE other      -0.067677722 0.15906837 Inf -0.37944600
## 8               No answer      -1.127339718 2.37764052 Inf -5.78742950
## 9             Other party       0.012412460 0.02297133 Inf -0.03261052
## 10  Pro-environment party       0.006731737 0.05866624 Inf -0.10825197
##     asymp.UCL
## 1  0.10240268
## 2  0.09875594
## 3  0.17505951
## 4  2.32873462
## 5  0.20076347
## 6  0.39800636
## 7  0.24409055
## 8  3.53275006
## 9  0.05743544
## 10 0.12171545
```

```r
H5.mod6.trends.tab$p<-
  2*(1-pnorm(abs(H5.mod6.trends.tab$env.eng.int.trend/
                   H5.mod6.trends.tab$SE)))
H5.mod6.trends.tab$adj.p<-
  p.adjust(H5.mod6.trends.tab$p,method="holm")

H5.mod6.trends.tab<-
  cbind(group=H5.mod6.trends.tab[,1],
      round(H5.mod6.trends.tab[,c(2,3)],2),
      round(H5.mod6.trends.tab[,c(7,8)],4),
      round(H5.mod6.trends.tab[,c(5,6)],2))
H5.mod6.trends.tab
```

```
##                     group env.eng.int.trend   SE      p adj.p asymp.LCL
## 1  Anti-immigration party              0.02 0.04 0.6317 1e+00     -0.06
## 2            Did not vote              0.03 0.03 0.2910 1e+00     -0.03
## 3              Don't know              0.04 0.07 0.6205 1e+00     -0.10
## 4            Invalid vote             -1.41 1.91 0.4606 1e+00     -5.14
## 5                  NE age              0.07 0.07 0.3439 1e+00     -0.07
## 6              NE citizen              0.28 0.06 0.0000 1e-04      0.16
## 7                NE other             -0.07 0.16 0.6705 1e+00     -0.38
## 8               No answer             -1.13 2.38 0.6354 1e+00     -5.79
## 9             Other party              0.01 0.02 0.5890 1e+00     -0.03
## 10  Pro-environment party              0.01 0.06 0.9086 1e+00     -0.11
##    asymp.UCL
## 1       0.10
## 2       0.10
## 3       0.18
## 4       2.33
## 5       0.20
## 6       0.40
## 7       0.24
## 8       3.53
## 9       0.06
## 10      0.12
```

```r
write.csv2(H5.mod6.trends.tab,"H5.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1 Pro-environment party - Anti-immigration party -0.01339152 0.06972524 Inf
##      z.ratio   p.value
## 1 -0.1920614 0.8476941
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect    0.261 0.342 Inf  0.764  1.0000 
##  Did not vote effect              0.277 0.341 Inf  0.814  1.0000 
##  Don't know effect                0.278 0.346 Inf  0.802  1.0000 
##  Invalid vote effect             -1.324 1.924 Inf -0.688  1.0000 
##  NE age effect                    0.311 0.346 Inf  0.900  1.0000 
##  NE citizen effect                0.548 0.345 Inf  1.590  1.0000 
##  NE other effect                  0.164 0.374 Inf  0.437  1.0000 
##  No answer effect                -1.014 2.387 Inf -0.425  1.0000 
##  Other party effect               0.253 0.340 Inf  0.743  1.0000 
##  Pro-environment party effect     0.246 0.344 Inf  0.716  1.0000 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H5.more.contrasts<-data.frame(pairs(H5.mod6.trends, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast     estimate         SE  df
## 1           Other party - Anti-immigration party -0.007710801 0.04412518 Inf
## 2 Pro-environment party - Anti-immigration party -0.013391524 0.06972524 Inf
## 3            Pro-environment party - Other party -0.005680723 0.06019848 Inf
##       z.ratio p.value
## 1 -0.17474833       1
## 2 -0.19206136       1
## 3 -0.09436655       1
```

