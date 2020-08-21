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
## Running under: Windows 10 x64 (build 18362)
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
##  [5] ggplot2_3.3.0   emmeans_1.4.5   psych_1.9.12.31 dplyr_0.8.5    
##  [9] lmerTest_3.1-2  lme4_1.1-23     Matrix_1.2-18  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6        mvtnorm_1.1-0       lattice_0.20-38    
##  [4] tidyr_1.0.2         zoo_1.8-7           foreach_1.5.0      
##  [7] assertthat_0.2.1    digest_0.6.25       mime_0.9           
## [10] R6_2.4.1            backports_1.1.6     evaluate_0.14      
## [13] coda_0.19-3         pillar_1.4.3        rlang_0.4.5        
## [16] multcomp_1.4-13     minqa_1.2.4         nloptr_1.2.2.1     
## [19] rmarkdown_2.1       splines_3.6.3       statmod_1.4.34     
## [22] stringr_1.4.0       munsell_0.5.0       shiny_1.4.0.2      
## [25] broom_0.5.5         httpuv_1.5.2        compiler_3.6.3     
## [28] numDeriv_2016.8-1.1 xfun_0.13           pkgconfig_2.0.3    
## [31] mnormt_1.5-6        htmltools_0.4.0     tidyselect_1.0.0   
## [34] tibble_3.0.0        codetools_0.2-16    fansi_0.4.1        
## [37] later_1.0.0         crayon_1.3.4        withr_2.1.2        
## [40] grid_3.6.3          nlme_3.1-144        xtable_1.8-4       
## [43] gtable_0.3.0        lifecycle_0.2.0     magrittr_1.5       
## [46] scales_1.1.0        estimability_1.3    cli_2.0.2          
## [49] stringi_1.4.6       promises_1.1.0      generics_0.0.2     
## [52] ellipsis_0.3.0      vctrs_0.2.4         boot_1.3-24        
## [55] sandwich_2.5-1      blme_1.0-4          TH.data_1.0-10     
## [58] iterators_1.0.12    tools_3.6.3         glue_1.4.0         
## [61] purrr_0.3.3         fastmap_1.0.1       abind_1.4-5        
## [64] parallel_3.6.3      survival_3.1-8      yaml_2.2.0         
## [67] colorspace_1.4-1    knitr_1.28
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
dat<-read.csv2("dat.no.miss.csv",stringsAsFactors = F)
```

### Variable transformations

#### Country


```r
table(dat$cntry)
```

```
## 
##   AT   BE   CH   CZ   DE   EE   ES   FI   FR   GB   HU   IE   IT   LT   NL   NO 
## 1973 1753 1503 2156 2819 1974 1817 1862 2015 1876 1391 2676 2317 1927 1661 1538 
##   PL   PT   SE   SI 
## 1589 1228 1525 1276
```

#### Voting group


```r
#make voting group variable names unique to each country
dat$voting.group<-paste0(dat$cntry,": ",dat$vote.group.combined)
```

#### Centering Attitudes towards the Environment


```r
#calculate the sum score
dat$environ<-rowMeans(data.frame(dat$inctxff.R,
                                 dat$sbsrnen.R,
                                 dat$banhhap.R),na.rm=T)
                                 
describe(dat$environ,fast=T)
```

```
##    vars     n mean   sd min max range se
## X1    1 36876 3.44 0.83   1   5     4  0
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
## X1    1 36876    0 0.79 -3.06 2.11  5.17  0
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
## [1] 0.32
## Sample Size 
## [1] 36545
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
##    vars     n mean  sd min max range se
## X1    1 36876 2.49 0.8   1   4     3  0
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
##    vars     n mean   sd   min  max range se
## X1    1 36876    0 0.74 -2.06 2.63  4.69  0
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
##    vars     n mean   sd min max range se
## X1    1 36831 2.42 0.91   1   4     3  0
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
## X1    1 36831    0 0.81 -2.19 2.65  4.85  0
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
## X1    1 36590 2.57 1.05   1   4     3 0.01
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
##    vars     n mean sd   min max range   se
## X1    1 36590    0  1 -2.19 2.6  4.79 0.01
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
dat$refugees<-rowMeans(data.frame(dat$gvrfgap.R,
                                 dat$rfgfrpc,
                                 dat$rfgbfml.R),na.rm=T)
                                 
describe(dat$refugees,fast=T)
```

```
##    vars     n mean   sd min max range se
## X1    1 36876 3.02 0.88   1   5     4  0
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
##    vars     n mean   sd   min  max range se
## X1    1 36876    0 0.77 -2.92 3.15  6.07  0
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
## 29063  7813
```

```r
#dummy-code "don't know"
dat$dont.know.dummy<-ifelse(grepl("Don't know",dat$vote.group.combined),1,0)
table(dat$dont.know.dummy)
```

```
## 
##     0     1 
## 35670  1206
```

```r
#dummy-code invalid vote
dat$invalid.vote.dummy<-ifelse(grepl("Invalid vote",dat$vote.group.combined),1,0)
table(dat$invalid.vote.dummy)
```

```
## 
##     0     1 
## 36861    15
```

```r
#dummy-code "no answer"
dat$no.answer.dummy<-ifelse(grepl("No answer",dat$vote.group.combined),1,0)
table(dat$no.answer.dummy)
```

```
## 
##     0     1 
## 36864    12
```

```r
#dummy-code not-eligible: age
dat$not.eligible.age.dummy<-ifelse(grepl("not eligible: age",dat$vote.group.combined),1,0)
table(dat$not.eligible.age.dummy)
```

```
## 
##     0     1 
## 35063  1813
```

```r
#dummy code not-eligible: citizenship
dat$not.eligible.citizenship.dummy<-ifelse(grepl("not eligible: citizenship",dat$vote.group.combined),1,0)
table(dat$not.eligible.citizenship.dummy)
```

```
## 
##     0     1 
## 35663  1213
```

```r
#dummy-code not-eligible: other reasons
dat$not.eligible.other.dummy<-ifelse(grepl("not eligible: other",dat$vote.group.combined),1,0)
table(dat$not.eligible.other.dummy)
```

```
## 
##     0     1 
## 36611   265
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
## 18517 18359
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
##     0 
## 36876
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
##    vars     n mean    sd median trimmed   mad    min   max range skew kurtosis
## X1    1 36876    0 18.49   0.64   -0.08 22.24 -34.36 50.64    85 0.02    -0.91
##     se
## X1 0.1
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
## X1    1 36876    0 16.18  -0.01   -0.18 17.71 -43.21 55.39  98.6 0.09    -0.56
##      se
## X1 0.08
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
                (environ.lvl1|cntry)+
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
## 1 voting.group  (Intercept)         <NA> 0.25793316 0.066529513
## 2 voting.group environ.lvl1         <NA> 0.06017025 0.003620459
## 3 voting.group  (Intercept) environ.lvl1 0.32665879 0.005069713
## 4        cntry  (Intercept)         <NA> 0.35903356 0.128905101
## 5        cntry environ.lvl1         <NA> 0.06041360 0.003649803
## 6        cntry  (Intercept) environ.lvl1 0.28937886 0.006276775
## 7     Residual         <NA>         <NA> 0.75432268 0.569002699
```

```r
getFE(EX3.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## age.lvl1.10                                                -0.01       0.00
## environ.lvl1                                                0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                49.49    0.19 0.847
## gender                                                  36694.83    7.33 0.000
## occupClerical support workers                           36604.49   -0.21 0.833
## occupCraft and related trades workers                   36617.70   -0.77 0.439
## occupElementary occupations                             36617.85    0.35 0.728
## occupManagers                                           36609.60    0.67 0.504
## occupOther: Not in paid work                            36783.74    2.04 0.042
## occupPlant and machine operators, and assemblers        36616.22   -0.33 0.743
## occupProfessionals                                      36612.40    1.67 0.096
## occupRetired                                            36610.89    0.40 0.686
## occupService and sales workers                          36613.54   -0.22 0.825
## occupSkilled agricultural, forestry and fishery workers 36619.41   -0.41 0.682
## occupTechnicians and associate professionals            36605.32    0.27 0.787
## occupUnemployed                                         36636.85   -1.02 0.307
## educ                                                    36835.85   13.35 0.000
## resid                                                   36735.51   -5.82 0.000
## age.lvl1.10                                             36645.83   -2.56 0.011
## environ.lvl1                                               19.92   10.04 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.23
## gender                                                   0.05  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## age.lvl1.10                                             -0.01  0.00
## environ.lvl1                                             0.12  0.18
```


\newpage

### Model 2 (interaction between age and environmental attitudes)


```r
EX3.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## EX3.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX3.mod1:     gender + occup + educ + resid + age.lvl1.10 + environ.lvl1
## EX3.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX3.mod2:     gender + occup + educ + resid + age.lvl1.10 + environ.lvl1 + 
## EX3.mod2:     age.lvl1.10:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX3.mod1   25 84768 84981 -42359    84718                         
## EX3.mod2   26 84751 84973 -42350    84699 18.638  1   1.58e-05 ***
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
## 1 voting.group  (Intercept)         <NA> 0.25820911 0.066671944
## 2 voting.group environ.lvl1         <NA> 0.05974302 0.003569229
## 3 voting.group  (Intercept) environ.lvl1 0.33075491 0.005102289
## 4        cntry  (Intercept)         <NA> 0.35889446 0.128805233
## 5        cntry environ.lvl1         <NA> 0.06068514 0.003682686
## 6        cntry  (Intercept) environ.lvl1 0.29974543 0.006528324
## 7     Residual         <NA>         <NA> 0.75413634 0.568721625
```

```r
getFE(EX3.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## age.lvl1.10                                                -0.01       0.00
## environ.lvl1                                                0.15       0.02
## age.lvl1.10:environ.lvl1                                   -0.01       0.00
##                                                               df t.value     p
## (Intercept)                                                49.49    0.18 0.857
## gender                                                  36694.79    7.31 0.000
## occupClerical support workers                           36604.53   -0.21 0.835
## occupCraft and related trades workers                   36617.69   -0.76 0.449
## occupElementary occupations                             36617.89    0.36 0.720
## occupManagers                                           36609.69    0.68 0.498
## occupOther: Not in paid work                            36783.80    2.02 0.043
## occupPlant and machine operators, and assemblers        36616.27   -0.32 0.748
## occupProfessionals                                      36612.41    1.66 0.097
## occupRetired                                            36610.87    0.34 0.732
## occupService and sales workers                          36613.55   -0.21 0.831
## occupSkilled agricultural, forestry and fishery workers 36619.48   -0.41 0.680
## occupTechnicians and associate professionals            36605.36    0.28 0.782
## occupUnemployed                                         36636.84   -0.98 0.329
## educ                                                    36835.94   13.43 0.000
## resid                                                   36735.91   -5.78 0.000
## age.lvl1.10                                             36645.70   -2.59 0.010
## environ.lvl1                                               19.92   10.07 0.000
## age.lvl1.10:environ.lvl1                                36627.74   -4.32 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.23
## gender                                                   0.05  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.08
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## age.lvl1.10                                             -0.01  0.00
## environ.lvl1                                             0.12  0.19
## age.lvl1.10:environ.lvl1                                -0.02 -0.01
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
## 1 -1.618158e+00          0.1752827 0.01612289 Inf 0.1436824 0.2068830
## 2 -3.759247e-18          0.1537212 0.01526370 Inf 0.1238049 0.1836375
## 3  1.618158e+00          0.1321598 0.01599597 Inf 0.1008082 0.1635113
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
## 1 -1.62               0.18 0.02 0     0      0.14      0.21
## 2  0.00               0.15 0.02 0     0      0.12      0.18
## 3  1.62               0.13 0.02 0     0      0.10      0.16
```

```r
pairs(EX3.mod2.trends,adjust="none")
```

```
##  contrast                                 estimate      SE  df z.ratio p.value
##  -1.61815846642672 - -3.7592467012089e-18   0.0216 0.00499 Inf 4.318   <.0001 
##  -1.61815846642672 - 1.61815846642672       0.0431 0.00999 Inf 4.318   <.0001 
##  -3.7592467012089e-18 - 1.61815846642672    0.0216 0.00499 Inf 4.318   <.0001 
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
## X1    1 36876 0.02 0.5    0.5    0.03   0 -0.5 0.5     1 -0.09    -1.99  0
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
##    vars     n mean   sd median trimmed  mad   min  max range  skew kurtosis se
## X1    1 36876    0 0.49   0.34       0 0.38 -0.79 0.86  1.64 -0.09    -1.91  0
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
                (environ.lvl1|cntry)+
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
## 1 voting.group  (Intercept)         <NA> 0.25762748 0.066371921
## 2 voting.group environ.lvl1         <NA> 0.06014494 0.003617414
## 3 voting.group  (Intercept) environ.lvl1 0.31969901 0.004953733
## 4        cntry  (Intercept)         <NA> 0.35890996 0.128816361
## 5        cntry environ.lvl1         <NA> 0.06041543 0.003650024
## 6        cntry  (Intercept) environ.lvl1 0.28612262 0.006204197
## 7     Residual         <NA>         <NA> 0.75432208 0.569001804
```

```r
getFE(EX1.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## gender.lvl1                                                 0.06       0.01
## environ.lvl1                                                0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                49.52    0.20 0.846
## age                                                     35789.18   -3.17 0.002
## occupClerical support workers                           36604.12   -0.20 0.838
## occupCraft and related trades workers                   36616.61   -0.79 0.431
## occupElementary occupations                             36615.24    0.35 0.730
## occupManagers                                           36610.01    0.68 0.498
## occupOther: Not in paid work                            36750.49    1.99 0.047
## occupPlant and machine operators, and assemblers        36615.75   -0.33 0.739
## occupProfessionals                                      36613.46    1.68 0.093
## occupRetired                                            36616.20    0.45 0.656
## occupService and sales workers                          36610.10   -0.23 0.821
## occupSkilled agricultural, forestry and fishery workers 36619.62   -0.41 0.684
## occupTechnicians and associate professionals            36605.15    0.27 0.785
## occupUnemployed                                         36618.71   -1.06 0.289
## educ                                                    36863.14   13.17 0.000
## resid                                                   36734.66   -5.80 0.000
## gender.lvl1                                             36579.83    7.14 0.000
## environ.lvl1                                               19.92   10.03 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.23
## age                                                      0.00  0.00
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## gender.lvl1                                              0.04  0.08
## environ.lvl1                                             0.12  0.18
```

\newpage

### Model 2 (interaction between sex and environmental attitudes)


```r
EX1.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## EX1.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX1.mod1:     age + occup + educ + resid + gender.lvl1 + environ.lvl1
## EX1.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX1.mod2:     age + occup + educ + resid + gender.lvl1 + environ.lvl1 + 
## EX1.mod2:     gender.lvl1:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX1.mod1   25 84768 84981 -42359    84718                     
## EX1.mod2   26 84769 84991 -42359    84717 0.4486  1      0.503
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
## 1 voting.group  (Intercept)         <NA> 0.25757040 0.066342510
## 2 voting.group environ.lvl1         <NA> 0.06016994 0.003620421
## 3 voting.group  (Intercept) environ.lvl1 0.31917201 0.004946526
## 4        cntry  (Intercept)         <NA> 0.35889735 0.128807307
## 5        cntry environ.lvl1         <NA> 0.06037731 0.003645420
## 6        cntry  (Intercept) environ.lvl1 0.28658871 0.006210165
## 7     Residual         <NA>         <NA> 0.75431800 0.568995649
```

```r
getFE(EX1.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## gender.lvl1                                                 0.06       0.01
## environ.lvl1                                                0.15       0.02
## gender.lvl1:environ.lvl1                                   -0.01       0.01
##                                                               df t.value     p
## (Intercept)                                                49.52    0.20 0.846
## age                                                     35788.22   -3.18 0.001
## occupClerical support workers                           36604.12   -0.20 0.838
## occupCraft and related trades workers                   36616.67   -0.78 0.434
## occupElementary occupations                             36615.24    0.35 0.730
## occupManagers                                           36610.01    0.68 0.498
## occupOther: Not in paid work                            36750.56    1.98 0.047
## occupPlant and machine operators, and assemblers        36615.76   -0.33 0.741
## occupProfessionals                                      36613.47    1.68 0.092
## occupRetired                                            36616.27    0.44 0.658
## occupService and sales workers                          36610.12   -0.23 0.821
## occupSkilled agricultural, forestry and fishery workers 36619.61   -0.40 0.686
## occupTechnicians and associate professionals            36605.13    0.27 0.783
## occupUnemployed                                         36618.68   -1.06 0.290
## educ                                                    36863.20   13.18 0.000
## resid                                                   36734.62   -5.80 0.000
## gender.lvl1                                             36581.11    7.15 0.000
## environ.lvl1                                               19.92   10.02 0.000
## gender.lvl1:environ.lvl1                                36615.61   -0.67 0.503
##                                                            LL    UL
## (Intercept)                                             -0.19  0.23
## age                                                      0.00  0.00
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## gender.lvl1                                              0.04  0.08
## environ.lvl1                                             0.12  0.18
## gender.lvl1:environ.lvl1                                -0.03  0.01
```

\newpage


## Does the association vary by education (years)?

### Center the education variable


```r
describe(dat$educ)
```

```
##    vars     n mean   sd median trimmed  mad    min   max range skew kurtosis
## X1    1 36876    0 3.87  -0.02   -0.06 2.97 -13.02 40.98    54 0.32     1.85
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
##    vars     n mean   sd median trimmed  mad    min   max range skew kurtosis
## X1    1 36876    0 3.57  -0.15   -0.08 3.24 -15.82 39.72 55.54 0.42     2.48
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
                (environ.lvl1|cntry)+
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
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.26536114 0.070416533
## 2 voting.group environ.lvl1         <NA> 0.06015460 0.003618576
## 3 voting.group  (Intercept) environ.lvl1 0.32648038 0.005211506
## 4        cntry  (Intercept)         <NA> 0.36099355 0.130316343
## 5        cntry environ.lvl1         <NA> 0.06036018 0.003643351
## 6        cntry  (Intercept) environ.lvl1 0.30441525 0.006633097
## 7     Residual         <NA>         <NA> 0.75432316 0.569003427
```

```r
getFE(EX4.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.02       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.09       0.08
## age                                                         0.00       0.00
## resid                                                      -0.05       0.01
## educ.lvl1                                                   0.02       0.00
## environ.lvl1                                                0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                49.01    0.23 0.821
## gender                                                  36680.04    7.35 0.000
## occupClerical support workers                           36600.03   -0.24 0.812
## occupCraft and related trades workers                   36614.70   -0.83 0.409
## occupElementary occupations                             36615.55    0.28 0.782
## occupManagers                                           36606.47    0.68 0.498
## occupOther: Not in paid work                            36754.86    1.89 0.059
## occupPlant and machine operators, and assemblers        36614.44   -0.38 0.704
## occupProfessionals                                      36611.50    1.69 0.090
## occupRetired                                            36614.27    0.40 0.693
## occupService and sales workers                          36607.35   -0.28 0.783
## occupSkilled agricultural, forestry and fishery workers 36617.94   -0.45 0.649
## occupTechnicians and associate professionals            36601.11    0.26 0.798
## occupUnemployed                                         36618.85   -1.11 0.266
## age                                                     36207.24   -3.32 0.001
## resid                                                   36736.47   -5.86 0.000
## educ.lvl1                                               36726.39   12.51 0.000
## environ.lvl1                                               19.91   10.07 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.23
## gender                                                   0.05  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.25
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.07
## age                                                      0.00  0.00
## resid                                                   -0.07 -0.03
## educ.lvl1                                                0.01  0.02
## environ.lvl1                                             0.12  0.18
```

\newpage

### Model 2 (interaction between education and environment attitudes)


```r
EX4.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## EX4.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX4.mod1:     gender + occup + age + resid + educ.lvl1 + environ.lvl1
## EX4.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX4.mod2:     gender + occup + age + resid + educ.lvl1 + environ.lvl1 + 
## EX4.mod2:     environ.lvl1:educ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX4.mod1   25 84782 84994 -42366    84732                         
## EX4.mod2   26 84751 84973 -42350    84699 32.337  1  1.296e-08 ***
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
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.26480998 0.070124327
## 2 voting.group environ.lvl1         <NA> 0.06009609 0.003611540
## 3 voting.group  (Intercept) environ.lvl1 0.33592326 0.005345898
## 4        cntry  (Intercept)         <NA> 0.35971231 0.129392946
## 5        cntry environ.lvl1         <NA> 0.06020536 0.003624685
## 6        cntry  (Intercept) environ.lvl1 0.32697954 0.007081268
## 7     Residual         <NA>         <NA> 0.75401009 0.568531210
```

```r
getFE(EX4.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.02       0.06
## occupCraft and related trades workers                      -0.06       0.06
## occupElementary occupations                                 0.01       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.06
## occupProfessionals                                          0.10       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.01       0.06
## occupUnemployed                                            -0.09       0.08
## age                                                         0.00       0.00
## resid                                                      -0.05       0.01
## educ.lvl1                                                   0.02       0.00
## environ.lvl1                                                0.15       0.02
## educ.lvl1:environ.lvl1                                      0.01       0.00
##                                                               df t.value     p
## (Intercept)                                                49.23    0.22 0.823
## gender                                                  36678.65    7.18 0.000
## occupClerical support workers                           36600.38   -0.26 0.792
## occupCraft and related trades workers                   36615.08   -0.87 0.383
## occupElementary occupations                             36615.88    0.22 0.826
## occupManagers                                           36606.82    0.64 0.519
## occupOther: Not in paid work                            36755.51    1.83 0.067
## occupPlant and machine operators, and assemblers        36614.87   -0.43 0.670
## occupProfessionals                                      36611.73    1.63 0.103
## occupRetired                                            36614.77    0.28 0.776
## occupService and sales workers                          36607.70   -0.30 0.765
## occupSkilled agricultural, forestry and fishery workers 36618.45   -0.52 0.605
## occupTechnicians and associate professionals            36601.42    0.23 0.821
## occupUnemployed                                         36619.42   -1.15 0.251
## age                                                     36197.48   -3.50 0.000
## resid                                                   36737.33   -5.79 0.000
## educ.lvl1                                               36726.15   12.29 0.000
## environ.lvl1                                               19.91   10.15 0.000
## educ.lvl1:environ.lvl1                                  36609.49    5.69 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.23
## gender                                                   0.04  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                            -0.01  0.25
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.17  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.06
## age                                                      0.00  0.00
## resid                                                   -0.07 -0.03
## educ.lvl1                                                0.01  0.02
## environ.lvl1                                             0.12  0.19
## educ.lvl1:environ.lvl1                                   0.01  0.01
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
## 1 -3.572954e+00          0.1264271 0.01588559 Inf 0.0952919 0.1575623
## 2  2.759112e-18          0.1540018 0.01517330 Inf 0.1242627 0.1837410
## 3  3.572954e+00          0.1815766 0.01597138 Inf 0.1502733 0.2128799
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
## 1 -3.57               0.13 0.02 0     0      0.10      0.16
## 2  0.00               0.15 0.02 0     0      0.12      0.18
## 3  3.57               0.18 0.02 0     0      0.15      0.21
```

\newpage

## Does the association vary by place of residence (urban/rural)?

### Center the residence variable


```r
describe(dat$resid)
```

```
##    vars     n  mean   sd median trimmed mad  min max range skew kurtosis se
## X1    1 36876 -0.11 0.49   -0.5   -0.14   0 -0.5 0.5     1 0.45    -1.79  0
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
## X1    1 36876    0 0.47  -0.25   -0.02 0.35 -0.86 0.97  1.83 0.41    -1.54  0
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
                (environ.lvl1|cntry)+
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
## 1 voting.group  (Intercept)         <NA> 0.25826547 0.066701051
## 2 voting.group environ.lvl1         <NA> 0.06012703 0.003615259
## 3 voting.group  (Intercept) environ.lvl1 0.31807394 0.004939286
## 4        cntry  (Intercept)         <NA> 0.35806598 0.128211246
## 5        cntry environ.lvl1         <NA> 0.06039734 0.003647839
## 6        cntry  (Intercept) environ.lvl1 0.28716657 0.006210331
## 7     Residual         <NA>         <NA> 0.75432146 0.569000870
```

```r
getFE(EX5.mod1)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## age                                                         0.00       0.00
## educ                                                        0.02       0.00
## resid.lvl1                                                 -0.05       0.01
## environ.lvl1                                                0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                49.65    0.25 0.807
## gender                                                  36685.27    7.36 0.000
## occupClerical support workers                           36603.29   -0.22 0.824
## occupCraft and related trades workers                   36616.32   -0.79 0.429
## occupElementary occupations                             36614.42    0.33 0.742
## occupManagers                                           36609.66    0.67 0.503
## occupOther: Not in paid work                            36750.45    1.97 0.049
## occupPlant and machine operators, and assemblers        36615.42   -0.34 0.734
## occupProfessionals                                      36612.48    1.67 0.096
## occupRetired                                            36615.10    0.43 0.667
## occupService and sales workers                          36609.18   -0.24 0.807
## occupSkilled agricultural, forestry and fishery workers 36621.25   -0.43 0.665
## occupTechnicians and associate professionals            36604.55    0.26 0.795
## occupUnemployed                                         36618.15   -1.07 0.284
## age                                                     35803.71   -3.19 0.001
## educ                                                    36862.99   13.21 0.000
## resid.lvl1                                              36563.00   -5.43 0.000
## environ.lvl1                                               19.92   10.04 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.23
## gender                                                   0.05  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## age                                                      0.00  0.00
## educ                                                     0.02  0.02
## resid.lvl1                                              -0.06 -0.03
## environ.lvl1                                             0.12  0.18
```

\newpage

### Model 2 (interaction between residence and environment attitudes)


```r
EX5.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## EX5.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX5.mod1:     gender + occup + age + educ + resid.lvl1 + environ.lvl1
## EX5.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX5.mod2:     gender + occup + age + educ + resid.lvl1 + environ.lvl1 + 
## EX5.mod2:     environ.lvl1:resid.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX5.mod1   25 84769 84982 -42359    84719                     
## EX5.mod2   26 84770 84992 -42359    84718 0.4625  1     0.4965
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
## 1 voting.group  (Intercept)         <NA> 0.25824349 0.066689702
## 2 voting.group environ.lvl1         <NA> 0.06002198 0.003602638
## 3 voting.group  (Intercept) environ.lvl1 0.31923956 0.004948304
## 4        cntry  (Intercept)         <NA> 0.35805193 0.128201183
## 5        cntry environ.lvl1         <NA> 0.06038924 0.003646860
## 6        cntry  (Intercept) environ.lvl1 0.28781783 0.006223336
## 7     Residual         <NA>         <NA> 0.75432011 0.568998827
```

```r
getFE(EX5.mod2)
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## gender                                                      0.06       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## age                                                         0.00       0.00
## educ                                                        0.02       0.00
## resid.lvl1                                                 -0.05       0.01
## environ.lvl1                                                0.15       0.02
## resid.lvl1:environ.lvl1                                    -0.01       0.01
##                                                               df t.value     p
## (Intercept)                                                49.66    0.24 0.809
## gender                                                  36685.36    7.36 0.000
## occupClerical support workers                           36603.41   -0.22 0.826
## occupCraft and related trades workers                   36616.44   -0.79 0.430
## occupElementary occupations                             36614.56    0.33 0.740
## occupManagers                                           36609.80    0.67 0.502
## occupOther: Not in paid work                            36750.56    1.97 0.049
## occupPlant and machine operators, and assemblers        36615.52   -0.34 0.735
## occupProfessionals                                      36612.58    1.67 0.096
## occupRetired                                            36615.14    0.43 0.669
## occupService and sales workers                          36609.29   -0.24 0.808
## occupSkilled agricultural, forestry and fishery workers 36621.32   -0.43 0.665
## occupTechnicians and associate professionals            36604.69    0.26 0.793
## occupUnemployed                                         36618.25   -1.07 0.284
## age                                                     35802.46   -3.19 0.001
## educ                                                    36862.95   13.20 0.000
## resid.lvl1                                              36560.76   -5.44 0.000
## environ.lvl1                                               19.92   10.03 0.000
## resid.lvl1:environ.lvl1                                 36603.51   -0.68 0.496
##                                                            LL    UL
## (Intercept)                                             -0.18  0.23
## gender                                                   0.05  0.08
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## age                                                      0.00  0.00
## educ                                                     0.02  0.02
## resid.lvl1                                              -0.06 -0.03
## environ.lvl1                                             0.12  0.18
## resid.lvl1:environ.lvl1                                 -0.03  0.01
```


\newpage

## Does the association vary by occupational groups

### Model 1 (Same as H1 selected model)


```r
EX2.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## 1 voting.group  (Intercept)         <NA> 0.25587191 0.065470436
## 2 voting.group environ.lvl1         <NA> 0.06016336 0.003619630
## 3 voting.group  (Intercept) environ.lvl1 0.31876624 0.004907124
## 4        cntry  (Intercept)         <NA> 0.35892016 0.128823682
## 5        cntry environ.lvl1         <NA> 0.06040006 0.003648167
## 6        cntry  (Intercept) environ.lvl1 0.28974115 0.006281240
## 7     Residual         <NA>         <NA> 0.75432125 0.569000555
```

\newpage


### Model 2 (Interaction between environment and occupational groups)


```r
EX2.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## EX2.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX2.mod1:     age + gender + educ + resid + occup + environ.lvl1
## EX2.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX2.mod2:     age + gender + educ + resid + occup + environ.lvl1 + occup:environ.lvl1
##          npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX2.mod1   25 84764 84977 -42357    84714                        
## EX2.mod2   37 84760 85075 -42343    84686 28.362 12   0.004896 **
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
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.25592354 0.065496859
## 2 voting.group environ.lvl1         <NA> 0.05693231 0.003241288
## 3 voting.group  (Intercept) environ.lvl1 0.28996773 0.004224922
## 4        cntry  (Intercept)         <NA> 0.35839191 0.128444761
## 5        cntry environ.lvl1         <NA> 0.05953179 0.003544033
## 6        cntry  (Intercept) environ.lvl1 0.29663639 0.006328948
## 7     Residual         <NA>         <NA> 0.75409475 0.568658890
```

\newpage

#### Marginal effects for each occupation group


```r
EX2.mod2.trends<-emtrends(EX2.mod2,specs = c("occup"),var=c("environ.lvl1"))
(EX2.mod2.trends.tab<-data.frame(EX2.mod2.trends))
```

```
##                                                 occup environ.lvl1.trend
## 1                                        Armed forces         0.08226301
## 2                            Clerical support workers         0.18341909
## 3                    Craft and related trades workers         0.15353627
## 4                              Elementary occupations         0.12067892
## 5                                            Managers         0.15412258
## 6                             Other: Not in paid work         0.20920355
## 7         Plant and machine operators, and assemblers         0.11719694
## 8                                       Professionals         0.16183601
## 9                                             Retired         0.11913605
## 10                          Service and sales workers         0.12411022
## 11 Skilled agricultural, forestry and fishery workers         0.12130618
## 12            Technicians and associate professionals         0.17656513
## 13                                         Unemployed         0.16922309
##            SE  df   asymp.LCL asymp.UCL
## 1  0.08601659 Inf -0.08632642 0.2508524
## 2  0.02292783 Inf  0.13848138 0.2283568
## 3  0.02066707 Inf  0.11302956 0.1940430
## 4  0.02189452 Inf  0.07776645 0.1635914
## 5  0.02345812 Inf  0.10814551 0.2000996
## 6  0.02520195 Inf  0.15980865 0.2585985
## 7  0.02316605 Inf  0.07179233 0.1626016
## 8  0.01883895 Inf  0.12491235 0.1987597
## 9  0.03908144 Inf  0.04253784 0.1957343
## 10 0.01892106 Inf  0.08702563 0.1611948
## 11 0.03510658 Inf  0.05249855 0.1901138
## 12 0.01995901 Inf  0.13744620 0.2156841
## 13 0.05648591 Inf  0.05851274 0.2799334
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
## 1                                        Armed forces               0.08 0.09
## 2                            Clerical support workers               0.18 0.02
## 3                    Craft and related trades workers               0.15 0.02
## 4                              Elementary occupations               0.12 0.02
## 5                                            Managers               0.15 0.02
## 6                             Other: Not in paid work               0.21 0.03
## 7         Plant and machine operators, and assemblers               0.12 0.02
## 8                                       Professionals               0.16 0.02
## 9                                             Retired               0.12 0.04
## 10                          Service and sales workers               0.12 0.02
## 11 Skilled agricultural, forestry and fishery workers               0.12 0.04
## 12            Technicians and associate professionals               0.18 0.02
## 13                                         Unemployed               0.17 0.06
##         p  adj.p asymp.LCL asymp.UCL
## 1  0.3389 0.3389     -0.09      0.25
## 2  0.0000 0.0000      0.14      0.23
## 3  0.0000 0.0000      0.11      0.19
## 4  0.0000 0.0000      0.08      0.16
## 5  0.0000 0.0000      0.11      0.20
## 6  0.0000 0.0000      0.16      0.26
## 7  0.0000 0.0000      0.07      0.16
## 8  0.0000 0.0000      0.12      0.20
## 9  0.0023 0.0069      0.04      0.20
## 10 0.0000 0.0000      0.09      0.16
## 11 0.0005 0.0022      0.05      0.19
## 12 0.0000 0.0000      0.14      0.22
## 13 0.0027 0.0069      0.06      0.28
```

```r
#contrast for all groups against mean of other groups
contrast(EX2.mod2.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                                                  estimate     SE  df
##  Armed forces effect                                       -0.06860 0.0852 Inf
##  Clerical support workers effect                            0.04099 0.0207 Inf
##  Craft and related trades workers effect                    0.00861 0.0182 Inf
##  Elementary occupations effect                             -0.02698 0.0195 Inf
##  Managers effect                                            0.00925 0.0213 Inf
##  Other: Not in paid work effect                             0.06892 0.0232 Inf
##  Plant and machine operators, and assemblers effect        -0.03075 0.0210 Inf
##  Professionals effect                                       0.01761 0.0162 Inf
##  Retired effect                                            -0.02865 0.0377 Inf
##  Service and sales workers effect                          -0.02326 0.0161 Inf
##  Skilled agricultural, forestry and fishery workers effect -0.02630 0.0336 Inf
##  Technicians and associate professionals effect             0.03356 0.0174 Inf
##  Unemployed effect                                          0.02561 0.0554 Inf
##  z.ratio p.value
##  -0.805  1.0000 
##   1.980  0.5720 
##   0.474  1.0000 
##  -1.385  1.0000 
##   0.434  1.0000 
##   2.964  0.0394 
##  -1.468  1.0000 
##   1.086  1.0000 
##  -0.761  1.0000 
##  -1.441  1.0000 
##  -0.783  1.0000 
##   1.926  0.5957 
##   0.462  1.0000 
## 
## Results are averaged over the levels of: gender, resid 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 13 tests
```




