---
title: 'Exploratory analysis: covariate moderation'
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


# Exploratory analyses for moderators

## Does the association vary by age?

### Center the age variable


```r
describe(dat$age)
```

```
##    vars     n  mean    sd median trimmed   mad   min  max range skew kurtosis
## X1    1 26886 -0.34 18.44    0.4   -0.45 22.24 -34.6 50.4    85 0.03     -0.9
##      se
## X1 0.11
```

```r
#already grand-mean centered

#obtain dataframe with country means and add to data

age.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(age.cntry=mean(age,na.rm=T))

dat<-left_join(x=dat,
               y=age.cntry,
               by=c("cntry"))

#center individuals around country means

dat$age.cntrymc<-dat$age-dat$age.cntry

#obtain dataframe with voting group means and add to data

age.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(age.voting.group=mean(age.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=age.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$age.vgmc<-dat$age.cntrymc-dat$age.voting.group

#describe the variable

describe(dat$age.vgmc)
```

```
##    vars     n mean    sd median trimmed   mad    min   max range skew kurtosis
## X1    1 26886    0 16.19  -0.05    -0.2 17.65 -41.74 55.49 97.23  0.1    -0.56
##     se
## X1 0.1
```

```r
#rename as lvl1, lvl2, and lvl3

dat$age.lvl1<-dat$age.vgmc
dat$age.lvl2<-dat$age.voting.group
dat$age.lvl3<-dat$age.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered age)

* Divide age variable by 10 to give interpretation by a decade


```r
dat$age.lvl1.10<-dat$age.lvl1/10

EX3.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+educ+resid+
                 age.lvl1.10+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX3.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX3.mod1<-getVC(EX3.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42043360 0.176764413
## 2 voting.group environ.lvl1         <NA> 0.09112496 0.008303758
## 3 voting.group  (Intercept) environ.lvl1 0.63814039 0.024448431
## 4        cntry environ.lvl1         <NA> 0.12006884 0.014416525
## 5     Residual         <NA>         <NA> 1.16262084 1.351687219
```

```r
getFE(EX3.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## gender                                                      0.09       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.04       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.02       0.13
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## age.lvl1.10                                                 0.00       0.00
## environ.lvl1                                                0.27       0.03
##                                                               df t.value     p
## (Intercept)                                             18678.17    0.39 0.696
## gender                                                  26746.65    5.61 0.000
## occupClerical support workers                           26695.69   -0.25 0.802
## occupCraft and related trades workers                   26703.25   -1.21 0.225
## occupElementary occupations                             26707.87    0.09 0.928
## occupManagers                                           26698.90    0.39 0.695
## occupOther: Not in paid work                            26823.37    1.64 0.101
## occupPlant and machine operators, and assemblers        26704.76   -0.71 0.476
## occupProfessionals                                      26694.18    1.54 0.125
## occupRetired                                            26697.55    0.34 0.737
## occupService and sales workers                          26699.16   -0.57 0.568
## occupSkilled agricultural, forestry and fishery workers 26704.38   -0.20 0.844
## occupTechnicians and associate professionals            26692.73   -0.13 0.898
## occupUnemployed                                         26733.57    0.13 0.900
## educ                                                    26862.33   12.48 0.000
## resid                                                   26830.34   -7.24 0.000
## age.lvl1.10                                             26721.90   -0.84 0.402
## environ.lvl1                                               15.94    8.23 0.000
##                                                            LL    UL
## (Intercept)                                             -0.17  0.26
## gender                                                   0.06  0.12
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.20  0.28
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.24  0.28
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## age.lvl1.10                                             -0.01  0.01
## environ.lvl1                                             0.20  0.34
```


\newpage

### Model 2 (interaction between age and environmental attitudes)


```r
EX3.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+educ+resid+
                 age.lvl1.10+
                environ.lvl1+
                 age.lvl1.10:environ.lvl1,
                 data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX3.mod1,EX3.mod2)
```

```
## Data: dat
## Models:
## EX3.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX3.mod1:     cntry) + gender + occup + educ + resid + age.lvl1.10 + environ.lvl1
## EX3.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX3.mod2:     cntry) + gender + occup + educ + resid + age.lvl1.10 + environ.lvl1 + 
## EX3.mod2:     age.lvl1.10:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX3.mod1   23 85080 85268 -42517    85034                         
## EX3.mod2   24 85068 85265 -42510    85020 13.883  1  0.0001945 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX3.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX3.mod2<-getVC(EX3.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42097276 0.177218068
## 2 voting.group environ.lvl1         <NA> 0.09081816 0.008247938
## 3 voting.group  (Intercept) environ.lvl1 0.63693050 0.024351109
## 4        cntry environ.lvl1         <NA> 0.12118443 0.014685667
## 5     Residual         <NA>         <NA> 1.16230398 1.350950543
```

```r
getFE(EX3.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## gender                                                      0.09       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.04       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.02       0.13
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## age.lvl1.10                                                 0.00       0.00
## environ.lvl1                                                0.27       0.03
## age.lvl1.10:environ.lvl1                                   -0.02       0.01
##                                                               df t.value     p
## (Intercept)                                             18647.92    0.37 0.714
## gender                                                  26746.16    5.57 0.000
## occupClerical support workers                           26695.37   -0.24 0.807
## occupCraft and related trades workers                   26702.84   -1.19 0.233
## occupElementary occupations                             26707.49    0.11 0.911
## occupManagers                                           26698.59    0.40 0.688
## occupOther: Not in paid work                            26823.12    1.63 0.104
## occupPlant and machine operators, and assemblers        26704.39   -0.70 0.484
## occupProfessionals                                      26693.89    1.54 0.124
## occupRetired                                            26697.06    0.29 0.775
## occupService and sales workers                          26698.79   -0.55 0.579
## occupSkilled agricultural, forestry and fishery workers 26704.08   -0.21 0.837
## occupTechnicians and associate professionals            26692.42   -0.12 0.907
## occupUnemployed                                         26732.98    0.15 0.880
## educ                                                    26862.45   12.55 0.000
## resid                                                   26830.06   -7.20 0.000
## age.lvl1.10                                             26721.61   -0.87 0.385
## environ.lvl1                                               15.95    8.21 0.000
## age.lvl1.10:environ.lvl1                                26569.04   -3.73 0.000
##                                                            LL    UL
## (Intercept)                                             -0.17  0.25
## gender                                                   0.06  0.12
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.23
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.21  0.28
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.24  0.28
## educ                                                     0.03  0.03
## resid                                                   -0.14 -0.08
## age.lvl1.10                                             -0.01  0.01
## environ.lvl1                                             0.20  0.34
## age.lvl1.10:environ.lvl1                                -0.03 -0.01
```


\newpage

#### Marginal effects for ages at -1SD and +1SD


```r
EX3.mod2.trends<-
  emtrends(EX3.mod2,specs = c("age.lvl1.10"),var=c("environ.lvl1"),
           at=list(age.lvl1.10=c(
             
             mean(dat$age.lvl1.10)-sd(dat$age.lvl1.10),
             mean(dat$age.lvl1.10),
             mean(dat$age.lvl1.10)+sd(dat$age.lvl1.10)
             )))

(EX3.mod2.trends.tab<-data.frame(EX3.mod2.trends))
```

```
##     age.lvl1.10 environ.lvl1.trend         SE  df asymp.LCL asymp.UCL
## 1 -1.619436e+00          0.3081306 0.03464194 Inf 0.2402337 0.3760276
## 2  1.108176e-18          0.2718712 0.03310899 Inf 0.2069787 0.3367636
## 3  1.619436e+00          0.2356117 0.03437513 Inf 0.1682377 0.3029857
```

```r
EX3.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX3.mod2.trends.tab$environ.lvl1.trend/
                   EX3.mod2.trends.tab$SE)))
EX3.mod2.trends.tab$adj.p<-
  p.adjust(EX3.mod2.trends.tab$p,method="holm")

EX3.mod2.trends.tab<-
  cbind(group=round(EX3.mod2.trends.tab[,1],2),
      round(EX3.mod2.trends.tab[,c(2,3)],2),
      round(EX3.mod2.trends.tab[,c(7,8)],4),
      round(EX3.mod2.trends.tab[,c(5,6)],2))
EX3.mod2.trends.tab
```

```
##   group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -1.62               0.31 0.03 0     0      0.24      0.38
## 2  0.00               0.27 0.03 0     0      0.21      0.34
## 3  1.62               0.24 0.03 0     0      0.17      0.30
```

```r
pairs(EX3.mod2.trends,adjust="none")
```

```
##  contrast                                 estimate      SE  df z.ratio p.value
##  -1.61943561763092 - 1.10817592059019e-18   0.0363 0.00973 Inf 3.727   0.0002 
##  -1.61943561763092 - 1.61943561763092       0.0725 0.01946 Inf 3.727   0.0002 
##  1.10817592059019e-18 - 1.61943561763092    0.0363 0.00973 Inf 3.727   0.0002 
## 
## Results are averaged over the levels of: gender, occup, resid 
## Degrees-of-freedom method: asymptotic
```

\newpage

## Does the association vary by sex?

### Center the sex variable


```r
describe(dat$gender)
```

```
##    vars     n mean  sd median trimmed mad  min max range  skew kurtosis se
## X1    1 26886 0.01 0.5    0.5    0.01   0 -0.5 0.5     1 -0.03       -2  0
```

```r
#grand mean center
dat$gender.gmc<-dat$gender-mean(dat$gender,na.rm=T)

#obtain dataframe with country means and add to data

gender.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(gender.cntry=mean(gender.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=gender.cntry,
               by=c("cntry"))

#center individuals around country means

dat$gender.cntrymc<-dat$gender.gmc-dat$gender.cntry

#obtain dataframe with voting group means and add to data

gender.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(gender.voting.group=mean(gender.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=gender.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$gender.vgmc<-dat$gender.cntrymc-dat$gender.voting.group

#describe the variable

describe(dat$gender.vgmc)
```

```
##    vars     n mean   sd median trimmed  mad  min  max range  skew kurtosis se
## X1    1 26886    0 0.49    0.3       0 0.52 -0.8 0.86  1.66 -0.03    -1.91  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$gender.lvl1<-dat$gender.vgmc
dat$gender.lvl2<-dat$gender.voting.group
dat$gender.lvl3<-dat$gender.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered sex)


```r
EX1.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+occup+educ+resid+
                 gender.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX1.mod1)
```

```
## [1] FALSE
```

```r
getVC(EX1.mod1)
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42080481 0.177076691
## 2 voting.group environ.lvl1         <NA> 0.09126277 0.008328893
## 3 voting.group  (Intercept) environ.lvl1 0.63874478 0.024530235
## 4        cntry environ.lvl1         <NA> 0.12012462 0.014429925
## 5     Residual         <NA>         <NA> 1.16261292 1.351668799
```

```r
getFE(EX1.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## gender.lvl1                                                 0.08       0.02
## environ.lvl1                                                0.27       0.03
##                                                               df t.value     p
## (Intercept)                                             18672.18    0.39 0.695
## age                                                     25964.98   -1.56 0.120
## occupClerical support workers                           26695.84   -0.25 0.802
## occupCraft and related trades workers                   26702.05   -1.23 0.219
## occupElementary occupations                             26705.20    0.08 0.936
## occupManagers                                           26699.00    0.40 0.689
## occupOther: Not in paid work                            26798.69    1.59 0.111
## occupPlant and machine operators, and assemblers        26704.17   -0.72 0.472
## occupProfessionals                                      26694.91    1.55 0.122
## occupRetired                                            26700.82    0.38 0.701
## occupService and sales workers                          26696.33   -0.58 0.560
## occupSkilled agricultural, forestry and fishery workers 26704.28   -0.20 0.845
## occupTechnicians and associate professionals            26692.57   -0.13 0.898
## occupUnemployed                                         26720.25    0.08 0.936
## educ                                                    26810.22   12.26 0.000
## resid                                                   26829.46   -7.22 0.000
## gender.lvl1                                             26649.01    5.47 0.000
## environ.lvl1                                               15.95    8.21 0.000
##                                                            LL    UL
## (Intercept)                                             -0.17  0.26
## age                                                      0.00  0.00
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## gender.lvl1                                              0.05  0.11
## environ.lvl1                                             0.20  0.34
```

\newpage

### Model 2 (interaction between sex and environmental attitudes)


```r
EX1.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+occup+educ+resid+
                 gender.lvl1+
                environ.lvl1+
                 gender.lvl1:environ.lvl1,
                 data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX1.mod1,EX1.mod2)
```

```
## Data: dat
## Models:
## EX1.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX1.mod1:     cntry) + age + occup + educ + resid + gender.lvl1 + environ.lvl1
## EX1.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX1.mod2:     cntry) + age + occup + educ + resid + gender.lvl1 + environ.lvl1 + 
## EX1.mod2:     gender.lvl1:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX1.mod1   23 85080 85269 -42517    85034                     
## EX1.mod2   24 85082 85279 -42517    85034 0.1721  1     0.6783
```

```r
isSingular(EX1.mod2)
```

```
## [1] FALSE
```

```r
getVC(EX1.mod2)
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42074353 0.177025114
## 2 voting.group environ.lvl1         <NA> 0.09139018 0.008352165
## 3 voting.group  (Intercept) environ.lvl1 0.63760448 0.024517057
## 4        cntry environ.lvl1         <NA> 0.12001477 0.014403545
## 5     Residual         <NA>         <NA> 1.16260739 1.351655932
```

```r
getFE(EX1.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## gender.lvl1                                                 0.08       0.02
## environ.lvl1                                                0.27       0.03
## gender.lvl1:environ.lvl1                                   -0.01       0.02
##                                                               df t.value     p
## (Intercept)                                             18674.20    0.39 0.695
## age                                                     25964.07   -1.56 0.119
## occupClerical support workers                           26695.81   -0.25 0.801
## occupCraft and related trades workers                   26702.20   -1.23 0.220
## occupElementary occupations                             26705.14    0.08 0.936
## occupManagers                                           26698.98    0.40 0.688
## occupOther: Not in paid work                            26798.62    1.59 0.112
## occupPlant and machine operators, and assemblers        26704.20   -0.72 0.473
## occupProfessionals                                      26694.95    1.55 0.121
## occupRetired                                            26700.75    0.38 0.703
## occupService and sales workers                          26696.29   -0.58 0.560
## occupSkilled agricultural, forestry and fishery workers 26704.33   -0.19 0.848
## occupTechnicians and associate professionals            26692.55   -0.13 0.898
## occupUnemployed                                         26720.25    0.08 0.935
## educ                                                    26810.79   12.26 0.000
## resid                                                   26829.45   -7.22 0.000
## gender.lvl1                                             26650.89    5.48 0.000
## environ.lvl1                                               15.95    8.21 0.000
## gender.lvl1:environ.lvl1                                26566.47   -0.41 0.678
##                                                            LL    UL
## (Intercept)                                             -0.17  0.26
## age                                                      0.00  0.00
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## gender.lvl1                                              0.05  0.11
## environ.lvl1                                             0.20  0.34
## gender.lvl1:environ.lvl1                                -0.05  0.03
```

\newpage


## Does the association vary by education (years)?

### Center the education variable


```r
describe(dat$educ)
```

```
##    vars     n mean   sd median trimmed  mad    min   max range skew kurtosis
## X1    1 26886 0.09 3.94  -0.15    0.04 4.45 -13.15 40.85    54 0.27     1.34
##      se
## X1 0.02
```

```r
#already grand-mean centered

#obtain dataframe with country means and add to data

educ.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(educ.cntry=mean(educ,na.rm=T))

dat<-left_join(x=dat,
               y=educ.cntry,
               by=c("cntry"))

#center individuals around country means

dat$educ.cntrymc<-dat$educ-dat$educ.cntry

#obtain dataframe with voting group means and add to data

educ.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(educ.voting.group=mean(educ.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=educ.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$educ.vgmc<-dat$educ.cntrymc-dat$educ.voting.group

#describe the variable

describe(dat$educ.vgmc)
```

```
##    vars     n mean   sd median trimmed  mad   min   max range skew kurtosis
## X1    1 26886    0 3.63  -0.15   -0.08 3.38 -15.9 39.72 55.62 0.37     1.79
##      se
## X1 0.02
```

```r
#rename as lvl1, lvl2, and lvl3

dat$educ.lvl1<-dat$educ.vgmc
dat$educ.lvl2<-dat$educ.voting.group
dat$educ.lvl3<-dat$educ.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered education)


```r
EX4.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+age+resid+
                 educ.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX4.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX4.mod1<-getVC(EX4.mod1))
```

```
##            grp         var1         var2     est_SD    est_SD2
## 1 voting.group  (Intercept)         <NA> 0.43257180 0.18711836
## 2 voting.group environ.lvl1         <NA> 0.09113935 0.00830638
## 3 voting.group  (Intercept) environ.lvl1 0.65216270 0.02571107
## 4        cntry environ.lvl1         <NA> 0.11826157 0.01398580
## 5     Residual         <NA>         <NA> 1.16263514 1.35172046
```

```r
getFE(EX4.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## gender                                                      0.09       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.04       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.12
## occupTechnicians and associate professionals               -0.02       0.11
## occupUnemployed                                             0.00       0.13
## age                                                         0.00       0.00
## resid                                                      -0.11       0.02
## educ.lvl1                                                   0.03       0.00
## environ.lvl1                                                0.27       0.03
##                                                               df t.value     p
## (Intercept)                                             18121.47    0.46 0.644
## gender                                                  26737.47    5.63 0.000
## occupClerical support workers                           26692.34   -0.28 0.778
## occupCraft and related trades workers                   26701.49   -1.27 0.205
## occupElementary occupations                             26704.92    0.02 0.981
## occupManagers                                           26695.33    0.40 0.690
## occupOther: Not in paid work                            26801.97    1.51 0.131
## occupPlant and machine operators, and assemblers        26703.49   -0.76 0.446
## occupProfessionals                                      26691.52    1.56 0.119
## occupRetired                                            26699.92    0.34 0.736
## occupService and sales workers                          26694.08   -0.63 0.531
## occupSkilled agricultural, forestry and fishery workers 26702.56   -0.23 0.817
## occupTechnicians and associate professionals            26689.06   -0.15 0.884
## occupUnemployed                                         26719.92    0.04 0.970
## age                                                     26258.30   -1.66 0.096
## resid                                                   26830.28   -7.27 0.000
## educ.lvl1                                               26757.52   11.72 0.000
## environ.lvl1                                               15.98    8.36 0.000
##                                                            LL    UL
## (Intercept)                                             -0.16  0.26
## gender                                                   0.06  0.12
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.35  0.07
## occupElementary occupations                             -0.21  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.30  0.13
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.20  0.28
## occupService and sales workers                          -0.28  0.14
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.23  0.19
## occupUnemployed                                         -0.26  0.27
## age                                                      0.00  0.00
## resid                                                   -0.14 -0.08
## educ.lvl1                                                0.02  0.03
## environ.lvl1                                             0.20  0.34
```

\newpage

### Model 2 (interaction between education and environment attitudes)


```r
EX4.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+age+resid+
                 educ.lvl1+
                environ.lvl1+
                 environ.lvl1:educ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

anova(EX4.mod1,EX4.mod2)
```

```
## Data: dat
## Models:
## EX4.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX4.mod1:     cntry) + gender + occup + age + resid + educ.lvl1 + environ.lvl1
## EX4.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX4.mod2:     cntry) + gender + occup + age + resid + educ.lvl1 + environ.lvl1 + 
## EX4.mod2:     environ.lvl1:educ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX4.mod1   23 85091 85280 -42523    85045                         
## EX4.mod2   24 85050 85247 -42501    85002 43.238  1  4.847e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX4.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX4.mod2<-getVC(EX4.mod2))
```

```
##            grp         var1         var2    est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.4311365 0.185878646
## 2 voting.group environ.lvl1         <NA> 0.0923120 0.008521506
## 3 voting.group  (Intercept) environ.lvl1 0.6620128 0.026347492
## 4        cntry environ.lvl1         <NA> 0.1177848 0.013873263
## 5     Residual         <NA>         <NA> 1.1617204 1.349594343
```

```r
getFE(EX4.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## gender                                                      0.08       0.02
## occupClerical support workers                              -0.04       0.11
## occupCraft and related trades workers                      -0.15       0.11
## occupElementary occupations                                -0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.16       0.11
## occupPlant and machine operators, and assemblers           -0.09       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.04       0.12
## occupTechnicians and associate professionals               -0.02       0.11
## occupUnemployed                                             0.00       0.13
## age                                                         0.00       0.00
## resid                                                      -0.11       0.02
## educ.lvl1                                                   0.03       0.00
## environ.lvl1                                                0.27       0.03
## educ.lvl1:environ.lvl1                                      0.02       0.00
##                                                               df t.value     p
## (Intercept)                                             18167.84    0.47 0.636
## gender                                                  26736.04    5.40 0.000
## occupClerical support workers                           26692.49   -0.33 0.744
## occupCraft and related trades workers                   26701.76   -1.35 0.176
## occupElementary occupations                             26705.08   -0.06 0.954
## occupManagers                                           26695.54    0.34 0.738
## occupOther: Not in paid work                            26802.47    1.41 0.159
## occupPlant and machine operators, and assemblers        26703.80   -0.83 0.408
## occupProfessionals                                      26691.55    1.47 0.142
## occupRetired                                            26700.38    0.22 0.829
## occupService and sales workers                          26694.22   -0.67 0.501
## occupSkilled agricultural, forestry and fishery workers 26703.00   -0.36 0.719
## occupTechnicians and associate professionals            26689.22   -0.20 0.839
## occupUnemployed                                         26720.14   -0.02 0.985
## age                                                     26235.59   -1.88 0.061
## resid                                                   26830.38   -7.15 0.000
## educ.lvl1                                               26754.36   11.49 0.000
## environ.lvl1                                               15.99    8.46 0.000
## educ.lvl1:environ.lvl1                                  26350.53    6.58 0.000
##                                                            LL    UL
## (Intercept)                                             -0.16  0.27
## gender                                                   0.05  0.11
## occupClerical support workers                           -0.25  0.18
## occupCraft and related trades workers                   -0.36  0.07
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.18  0.25
## occupOther: Not in paid work                            -0.06  0.38
## occupPlant and machine operators, and assemblers        -0.30  0.12
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.21  0.27
## occupService and sales workers                          -0.28  0.14
## occupSkilled agricultural, forestry and fishery workers -0.27  0.18
## occupTechnicians and associate professionals            -0.23  0.19
## occupUnemployed                                         -0.26  0.26
## age                                                      0.00  0.00
## resid                                                   -0.14 -0.08
## educ.lvl1                                                0.02  0.03
## environ.lvl1                                             0.20  0.34
## educ.lvl1:environ.lvl1                                   0.01  0.02
```


\newpage

#### Marginal effects for different levels of education


```r
EX4.mod2.trends<-
  emtrends(EX4.mod2,specs = c("educ.lvl1"),var=c("environ.lvl1"),
           at=
             list(educ.lvl1=
                    c(mean(dat$educ.lvl1)-sd(dat$educ.lvl1),
                      mean(dat$educ.lvl1),                              mean(dat$educ.lvl1)+sd(dat$educ.lvl1))))
(EX4.mod2.trends.tab<-data.frame(EX4.mod2.trends))
```

```
##       educ.lvl1 environ.lvl1.trend         SE  df asymp.LCL asymp.UCL
## 1 -3.631213e+00          0.2113345 0.03359684 Inf 0.1454859 0.2771831
## 2  1.581291e-18          0.2734103 0.03233354 Inf 0.2100378 0.3367829
## 3  3.631213e+00          0.3354862 0.03376648 Inf 0.2693051 0.4016673
```

```r
EX4.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX4.mod2.trends.tab$environ.lvl1.trend/
                   EX4.mod2.trends.tab$SE)))
EX4.mod2.trends.tab$adj.p<-
  p.adjust(EX4.mod2.trends.tab$p,method="holm")

EX4.mod2.trends.tab<-
  cbind(group=round(EX4.mod2.trends.tab[,1],2),
        round(EX4.mod2.trends.tab[,c(2,3)],2),
        round(EX4.mod2.trends.tab[,c(7,8)],4),
        round(EX4.mod2.trends.tab[,c(5,6)],2))
EX4.mod2.trends.tab
```

```
##   group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -3.63               0.21 0.03 0     0      0.15      0.28
## 2  0.00               0.27 0.03 0     0      0.21      0.34
## 3  3.63               0.34 0.03 0     0      0.27      0.40
```

\newpage

## Does the association vary by place of residence (urban/rural)?

### Center the residence variable


```r
describe(dat$resid)
```

```
##    vars     n mean   sd median trimmed mad  min max range skew kurtosis se
## X1    1 26886 -0.1 0.49   -0.5   -0.12   0 -0.5 0.5     1  0.4    -1.84  0
```

```r
dat$resid.gmc<-dat$resid-mean(dat$resid,na.rm=T)


#obtain dataframe with country means and add to data

resid.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(resid.cntry=mean(resid.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=resid.cntry,
               by=c("cntry"))

#center individuals around country means

dat$resid.cntrymc<-dat$resid.gmc-dat$resid.cntry

#obtain dataframe with voting group means and add to data

resid.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(resid.voting.group=mean(resid.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=resid.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$resid.vgmc<-dat$resid.cntrymc-dat$resid.voting.group

#describe the variable

describe(dat$resid.vgmc)
```

```
##    vars     n mean   sd median trimmed  mad   min  max range skew kurtosis se
## X1    1 26886    0 0.47  -0.25   -0.02 0.34 -0.84 0.98  1.82 0.36    -1.58  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$resid.lvl1<-dat$resid.vgmc
dat$resid.lvl2<-dat$resid.voting.group
dat$resid.lvl3<-dat$resid.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered residence)


```r
EX5.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+age+educ+
                 resid.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX5.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX5.mod1<-getVC(EX5.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42273069 0.178701240
## 2 voting.group environ.lvl1         <NA> 0.09095172 0.008272216
## 3 voting.group  (Intercept) environ.lvl1 0.63640585 0.024468587
## 4        cntry environ.lvl1         <NA> 0.11998503 0.014396408
## 5     Residual         <NA>         <NA> 1.16262125 1.351688179
```

```r
getFE(EX5.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## gender                                                      0.09       0.02
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
## age                                                         0.00       0.00
## educ                                                        0.03       0.00
## resid.lvl1                                                 -0.11       0.02
## environ.lvl1                                                0.27       0.03
##                                                               df t.value     p
## (Intercept)                                             18580.43    0.49 0.626
## gender                                                  26740.96    5.63 0.000
## occupClerical support workers                           26694.56   -0.26 0.792
## occupCraft and related trades workers                   26701.43   -1.23 0.219
## occupElementary occupations                             26704.12    0.07 0.943
## occupManagers                                           26698.42    0.40 0.692
## occupOther: Not in paid work                            26798.28    1.58 0.114
## occupPlant and machine operators, and assemblers        26703.44   -0.72 0.470
## occupProfessionals                                      26693.80    1.54 0.124
## occupRetired                                            26699.63    0.37 0.711
## occupService and sales workers                          26695.16   -0.60 0.552
## occupSkilled agricultural, forestry and fishery workers 26704.41   -0.21 0.831
## occupTechnicians and associate professionals            26691.83   -0.14 0.890
## occupUnemployed                                         26719.59    0.07 0.942
## age                                                     25981.86   -1.57 0.116
## educ                                                    26807.92   12.29 0.000
## resid.lvl1                                              26638.58   -6.94 0.000
## environ.lvl1                                               15.96    8.22 0.000
##                                                            LL    UL
## (Intercept)                                             -0.16  0.27
## gender                                                   0.06  0.12
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.21  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.20  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.19
## occupUnemployed                                         -0.25  0.27
## age                                                      0.00  0.00
## educ                                                     0.02  0.03
## resid.lvl1                                              -0.14 -0.08
## environ.lvl1                                             0.20  0.34
```

\newpage

### Model 2 (interaction between residence and environment attitudes)


```r
EX5.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                gender+occup+age+educ+
                 resid.lvl1+
                environ.lvl1+
                 environ.lvl1:resid.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

anova(EX5.mod1,EX5.mod2)
```

```
## Data: dat
## Models:
## EX5.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX5.mod1:     cntry) + gender + occup + age + educ + resid.lvl1 + environ.lvl1
## EX5.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX5.mod2:     cntry) + gender + occup + age + educ + resid.lvl1 + environ.lvl1 + 
## EX5.mod2:     environ.lvl1:resid.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX5.mod1   23 85082 85271 -42518    85036                     
## EX5.mod2   24 85084 85281 -42518    85036 0.0317  1     0.8588
```

```r
isSingular(EX5.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX5.mod2<-getVC(EX5.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.42271829 0.178690754
## 2 voting.group environ.lvl1         <NA> 0.09090207 0.008263186
## 3 voting.group  (Intercept) environ.lvl1 0.63681574 0.024470260
## 4        cntry environ.lvl1         <NA> 0.11997201 0.014393283
## 5     Residual         <NA>         <NA> 1.16262226 1.351690519
```

```r
getFE(EX5.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.05       0.11
## gender                                                      0.09       0.02
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
## age                                                         0.00       0.00
## educ                                                        0.03       0.00
## resid.lvl1                                                 -0.11       0.02
## environ.lvl1                                                0.27       0.03
## resid.lvl1:environ.lvl1                                     0.00       0.02
##                                                               df t.value     p
## (Intercept)                                             18582.66    0.49 0.627
## gender                                                  26740.98    5.63 0.000
## occupClerical support workers                           26694.58   -0.26 0.793
## occupCraft and related trades workers                   26701.45   -1.23 0.219
## occupElementary occupations                             26704.16    0.07 0.942
## occupManagers                                           26698.46    0.40 0.691
## occupOther: Not in paid work                            26798.41    1.58 0.113
## occupPlant and machine operators, and assemblers        26703.45   -0.72 0.470
## occupProfessionals                                      26693.87    1.54 0.124
## occupRetired                                            26699.59    0.37 0.711
## occupService and sales workers                          26695.21   -0.59 0.553
## occupSkilled agricultural, forestry and fishery workers 26704.42   -0.21 0.831
## occupTechnicians and associate professionals            26691.88   -0.14 0.891
## occupUnemployed                                         26719.71    0.07 0.941
## age                                                     25980.98   -1.57 0.116
## educ                                                    26807.27   12.28 0.000
## resid.lvl1                                              26637.77   -6.94 0.000
## environ.lvl1                                               15.95    8.22 0.000
## resid.lvl1:environ.lvl1                                 26520.53   -0.18 0.859
##                                                            LL    UL
## (Intercept)                                             -0.16  0.27
## gender                                                   0.06  0.12
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.21  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.20  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## age                                                      0.00  0.00
## educ                                                     0.02  0.03
## resid.lvl1                                              -0.14 -0.08
## environ.lvl1                                             0.20  0.34
## resid.lvl1:environ.lvl1                                 -0.04  0.04
```


\newpage

## Does the association vary by occupational groups

### Model 1 (Same as H1 selected model)


```r
EX2.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX2.mod1)
```

```
## [1] FALSE
```

```r
getVC(EX2.mod1)
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41857780 0.175207375
## 2 voting.group environ.lvl1         <NA> 0.09120255 0.008317904
## 3 voting.group  (Intercept) environ.lvl1 0.63782258 0.024349107
## 4        cntry environ.lvl1         <NA> 0.11995018 0.014388045
## 5     Residual         <NA>         <NA> 1.16262160 1.351688977
```

\newpage


### Model 2 (Interaction between environment and occupational groups)


```r
EX2.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                 occup:environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX2.mod1,EX2.mod2)
```

```
## Data: dat
## Models:
## EX2.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX2.mod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1
## EX2.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## EX2.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX2.mod2:     occup:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX2.mod1   23 85078 85267 -42516    85032                       
## EX2.mod2   35 85080 85367 -42505    85010 21.587 12    0.04242 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX2.mod2)
```

```
## [1] FALSE
```

```r
getVC(EX2.mod2)
```

```
##            grp         var1         var2     est_SD    est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41813605 0.17483775
## 2 voting.group environ.lvl1         <NA> 0.08973199 0.00805183
## 3 voting.group  (Intercept) environ.lvl1 0.61298129 0.02299917
## 4        cntry environ.lvl1         <NA> 0.11802202 0.01392920
## 5     Residual         <NA>         <NA> 1.16216050 1.35061703
```

\newpage

#### Marginal effects for each occupation group


```r
EX2.mod2.trends<-emtrends(EX2.mod2,specs = c("occup"),var=c("environ.lvl1"))
(EX2.mod2.trends.tab<-data.frame(EX2.mod2.trends))
```

```
##                                                 occup environ.lvl1.trend
## 1                                        Armed forces          0.1488196
## 2                            Clerical support workers          0.3132550
## 3                    Craft and related trades workers          0.2862463
## 4                              Elementary occupations          0.2354232
## 5                                            Managers          0.2838370
## 6                             Other: Not in paid work          0.3073434
## 7         Plant and machine operators, and assemblers          0.1685377
## 8                                       Professionals          0.2955657
## 9                                             Retired          0.2765820
## 10                          Service and sales workers          0.2167257
## 11 Skilled agricultural, forestry and fishery workers          0.2307368
## 12            Technicians and associate professionals          0.3210390
## 13                                         Unemployed          0.2112246
##            SE  df    asymp.LCL asymp.UCL
## 1  0.14957075 Inf -0.144333696 0.4419729
## 2  0.04616785 Inf  0.222767650 0.4037423
## 3  0.04365816 Inf  0.200677898 0.3718147
## 4  0.04537882 Inf  0.146482335 0.3243640
## 5  0.04718220 Inf  0.191361604 0.3763124
## 6  0.05317431 Inf  0.203123628 0.4115631
## 7  0.04933858 Inf  0.071835883 0.2652396
## 8  0.03839569 Inf  0.220311552 0.3708199
## 9  0.08250632 Inf  0.114872622 0.4382915
## 10 0.03923800 Inf  0.139820656 0.2936308
## 11 0.06863929 Inf  0.096206299 0.3652674
## 12 0.04076030 Inf  0.241150256 0.4009277
## 13 0.10557905 Inf  0.004293424 0.4181557
```

```r
EX2.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX2.mod2.trends.tab$environ.lvl1.trend/
                   EX2.mod2.trends.tab$SE)))
EX2.mod2.trends.tab$adj.p<-
  p.adjust(EX2.mod2.trends.tab$p,method="holm")

EX2.mod2.trends.tab<-
  cbind(group=EX2.mod2.trends.tab[,1],
        round(EX2.mod2.trends.tab[,c(2,3)],2),
        round(EX2.mod2.trends.tab[,c(7,8)],4),
        round(EX2.mod2.trends.tab[,c(5,6)],2))
EX2.mod2.trends.tab
```

```
##                                                 group environ.lvl1.trend   SE
## 1                                        Armed forces               0.15 0.15
## 2                            Clerical support workers               0.31 0.05
## 3                    Craft and related trades workers               0.29 0.04
## 4                              Elementary occupations               0.24 0.05
## 5                                            Managers               0.28 0.05
## 6                             Other: Not in paid work               0.31 0.05
## 7         Plant and machine operators, and assemblers               0.17 0.05
## 8                                       Professionals               0.30 0.04
## 9                                             Retired               0.28 0.08
## 10                          Service and sales workers               0.22 0.04
## 11 Skilled agricultural, forestry and fishery workers               0.23 0.07
## 12            Technicians and associate professionals               0.32 0.04
## 13                                         Unemployed               0.21 0.11
##         p  adj.p asymp.LCL asymp.UCL
## 1  0.3197 0.3197     -0.14      0.44
## 2  0.0000 0.0000      0.22      0.40
## 3  0.0000 0.0000      0.20      0.37
## 4  0.0000 0.0000      0.15      0.32
## 5  0.0000 0.0000      0.19      0.38
## 6  0.0000 0.0000      0.20      0.41
## 7  0.0006 0.0032      0.07      0.27
## 8  0.0000 0.0000      0.22      0.37
## 9  0.0008 0.0032      0.11      0.44
## 10 0.0000 0.0000      0.14      0.29
## 11 0.0008 0.0032      0.10      0.37
## 12 0.0000 0.0000      0.24      0.40
## 13 0.0454 0.0909      0.00      0.42
```

```r
#contrast for all groups against mean of other groups
contrast(EX2.mod2.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                                                  estimate     SE  df
##  Armed forces effect                                        -0.1134 0.1472 Inf
##  Clerical support workers effect                             0.0647 0.0390 Inf
##  Craft and related trades workers effect                     0.0355 0.0361 Inf
##  Elementary occupations effect                              -0.0196 0.0380 Inf
##  Managers effect                                             0.0329 0.0402 Inf
##  Other: Not in paid work effect                              0.0583 0.0469 Inf
##  Plant and machine operators, and assemblers effect         -0.0920 0.0429 Inf
##  Professionals effect                                        0.0456 0.0297 Inf
##  Retired effect                                              0.0250 0.0784 Inf
##  Service and sales workers effect                           -0.0398 0.0306 Inf
##  Skilled agricultural, forestry and fishery workers effect  -0.0246 0.0638 Inf
##  Technicians and associate professionals effect              0.0732 0.0326 Inf
##  Unemployed effect                                          -0.0458 0.1022 Inf
##  z.ratio p.value
##  -0.770  1.0000 
##   1.658  1.0000 
##   0.983  1.0000 
##  -0.515  1.0000 
##   0.817  1.0000 
##   1.244  1.0000 
##  -2.147  0.3816 
##   1.535  1.0000 
##   0.319  1.0000 
##  -1.302  1.0000 
##  -0.386  1.0000 
##   2.245  0.3221 
##  -0.448  1.0000 
## 
## Results are averaged over the levels of: gender, resid 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 13 tests
```




