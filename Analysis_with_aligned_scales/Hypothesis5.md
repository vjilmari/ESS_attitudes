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

# Hypothesis 5: The strength of the association between environment and refugee attitudes is stronger among more politically engaged individuals. 

### Model 1: without interactions (only main effects)


```r
H5.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.03       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.04       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## environ.lvl1                                                0.15       0.02
## engagement.lvl1                                             0.04       0.01
##                                                               df t.value     p
## (Intercept)                                                49.47    0.15 0.882
## age                                                     35483.71   -4.84 0.000
## gender                                                  36678.99    8.20 0.000
## educ                                                    36858.06   12.39 0.000
## resid                                                   36736.44   -5.58 0.000
## occupClerical support workers                           36605.28   -0.18 0.857
## occupCraft and related trades workers                   36618.42   -0.71 0.478
## occupElementary occupations                             36617.15    0.45 0.652
## occupManagers                                           36611.54    0.68 0.498
## occupOther: Not in paid work                            36754.25    2.01 0.045
## occupPlant and machine operators, and assemblers        36617.49   -0.26 0.795
## occupProfessionals                                      36614.45    1.69 0.090
## occupRetired                                            36617.59    0.50 0.618
## occupService and sales workers                          36611.39   -0.17 0.863
## occupSkilled agricultural, forestry and fishery workers 36621.42   -0.33 0.739
## occupTechnicians and associate professionals            36606.50    0.29 0.770
## occupUnemployed                                         36620.54   -1.01 0.313
## environ.lvl1                                               19.95    9.93 0.000
## engagement.lvl1                                         36710.18    7.07 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.10  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.15  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## environ.lvl1                                             0.12  0.18
## engagement.lvl1                                          0.03  0.05
```

```r
(VC.H5.mod1<-getVC(H5.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.25529255 0.065174288
## 2 voting.group environ.lvl1         <NA> 0.05941664 0.003530338
## 3 voting.group  (Intercept) environ.lvl1 0.30027063 0.004554693
## 4        cntry  (Intercept)         <NA> 0.35902650 0.128900029
## 5        cntry environ.lvl1         <NA> 0.06001899 0.003602279
## 6        cntry  (Intercept) environ.lvl1 0.28040364 0.006042252
## 7     Residual         <NA>         <NA> 0.75382726 0.568255545
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political engagement


```r
H5.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## H5.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.mod1:     age + gender + educ + resid + occup + environ.lvl1 + engagement.lvl1
## H5.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.mod2:     age + gender + educ + resid + occup + environ.lvl1 + engagement.lvl1 + 
## H5.mod2:     environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.mod1   26 84717 84938 -42332    84665                        
## H5.mod2   27 84712 84942 -42329    84658 6.8511  1   0.008859 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod2<-getFE(H5.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.03       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## environ.lvl1                                                0.15       0.02
## engagement.lvl1                                             0.04       0.01
## environ.lvl1:engagement.lvl1                                0.02       0.01
##                                                               df t.value     p
## (Intercept)                                                49.52    0.15 0.884
## age                                                     35473.19   -4.86 0.000
## gender                                                  36679.13    8.20 0.000
## educ                                                    36858.02   12.42 0.000
## resid                                                   36736.61   -5.58 0.000
## occupClerical support workers                           36605.35   -0.19 0.852
## occupCraft and related trades workers                   36618.46   -0.72 0.474
## occupElementary occupations                             36617.21    0.44 0.659
## occupManagers                                           36611.60    0.67 0.504
## occupOther: Not in paid work                            36754.44    2.00 0.045
## occupPlant and machine operators, and assemblers        36617.52   -0.27 0.788
## occupProfessionals                                      36614.52    1.69 0.092
## occupRetired                                            36617.62    0.48 0.629
## occupService and sales workers                          36611.47   -0.18 0.860
## occupSkilled agricultural, forestry and fishery workers 36621.51   -0.33 0.738
## occupTechnicians and associate professionals            36606.57    0.29 0.772
## occupUnemployed                                         36620.24   -1.03 0.302
## environ.lvl1                                               19.95    9.96 0.000
## engagement.lvl1                                         36709.23    7.12 0.000
## environ.lvl1:engagement.lvl1                            36578.40    2.62 0.009
##                                                            LL    UL
## (Intercept)                                             -0.19  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.10  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.15  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.23  0.07
## environ.lvl1                                             0.12  0.18
## engagement.lvl1                                          0.03  0.05
## environ.lvl1:engagement.lvl1                             0.00  0.03
```

```r
(VC.H5.mod2<-getVC(H5.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.25495174 0.065000389
## 2 voting.group environ.lvl1         <NA> 0.05959337 0.003551369
## 3 voting.group  (Intercept) environ.lvl1 0.30154898 0.004581564
## 4        cntry  (Intercept)         <NA> 0.35873329 0.128689574
## 5        cntry environ.lvl1         <NA> 0.05970281 0.003564426
## 6        cntry  (Intercept) environ.lvl1 0.28330454 0.006067643
## 7     Residual         <NA>         <NA> 0.75376191 0.568157017
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political engagement, allow engagement effect to vary between voting groups and countries


```r
H5.mod3<-lmer(refugees~(environ.lvl1+engagement.lvl1|voting.group)+
                (environ.lvl1+engagement.lvl1|cntry)+
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
## H5.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.mod2:     age + gender + educ + resid + occup + environ.lvl1 + engagement.lvl1 + 
## H5.mod2:     environ.lvl1:engagement.lvl1
## H5.mod3: refugees ~ (environ.lvl1 + engagement.lvl1 | voting.group) + 
## H5.mod3:     (environ.lvl1 + engagement.lvl1 | cntry) + age + gender + 
## H5.mod3:     educ + resid + occup + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.mod2   27 84712 84942 -42329    84658                         
## H5.mod3   33 84674 84955 -42304    84608 50.177  6  4.333e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod3<-getFE(H5.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.03       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## environ.lvl1                                                0.15       0.01
## engagement.lvl1                                             0.04       0.01
## environ.lvl1:engagement.lvl1                                0.02       0.01
##                                                               df t.value     p
## (Intercept)                                                49.43    0.15 0.879
## age                                                     35128.34   -4.75 0.000
## gender                                                  36665.93    8.25 0.000
## educ                                                    36846.43   12.31 0.000
## resid                                                   36706.96   -5.56 0.000
## occupClerical support workers                           36582.60   -0.19 0.849
## occupCraft and related trades workers                   36598.42   -0.71 0.475
## occupElementary occupations                             36597.31    0.43 0.665
## occupManagers                                           36589.95    0.68 0.498
## occupOther: Not in paid work                            36734.10    1.99 0.047
## occupPlant and machine operators, and assemblers        36598.84   -0.26 0.792
## occupProfessionals                                      36591.68    1.66 0.097
## occupRetired                                            36597.25    0.47 0.642
## occupService and sales workers                          36590.77   -0.17 0.865
## occupSkilled agricultural, forestry and fishery workers 36603.31   -0.32 0.748
## occupTechnicians and associate professionals            36584.66    0.29 0.771
## occupUnemployed                                         36591.92   -1.09 0.278
## environ.lvl1                                               19.89    9.97 0.000
## engagement.lvl1                                            22.72    4.81 0.000
## environ.lvl1:engagement.lvl1                            36562.78    2.67 0.008
##                                                            LL    UL
## (Intercept)                                             -0.19  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.15  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.07
## environ.lvl1                                             0.12  0.18
## engagement.lvl1                                          0.02  0.06
## environ.lvl1:engagement.lvl1                             0.00  0.03
```

```r
(VC.H5.mod3<-getVC(H5.mod3))
```

```
##             grp            var1            var2     est_SD      est_SD2
## 1  voting.group     (Intercept)            <NA> 0.25533021 0.0651935174
## 2  voting.group    environ.lvl1            <NA> 0.05784013 0.0033454804
## 3  voting.group engagement.lvl1            <NA> 0.05914156 0.0034977242
## 4  voting.group     (Intercept)    environ.lvl1 0.29024318 0.0042864077
## 5  voting.group     (Intercept) engagement.lvl1 0.30986622 0.0046791742
## 6  voting.group    environ.lvl1 engagement.lvl1 0.05491270 0.0001878429
## 7         cntry     (Intercept)            <NA> 0.35894957 0.1288447928
## 8         cntry    environ.lvl1            <NA> 0.05909483 0.0034921993
## 9         cntry engagement.lvl1            <NA> 0.02363258 0.0005584989
## 10        cntry     (Intercept)    environ.lvl1 0.28122690 0.0059654033
## 11        cntry     (Intercept) engagement.lvl1 0.35326397 0.0029967047
## 12        cntry    environ.lvl1 engagement.lvl1 0.11896721 0.0001661453
## 13     Residual            <NA>            <NA> 0.75233876 0.5660136074
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.mod4<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+
                   environ.lvl1:engagement.lvl1|voting.group)+
                (environ.lvl1+engagement.lvl1+
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

```
## Warning: Model failed to converge with 1 negative eigenvalue: -4.9e+02
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
## H5.mod3:     (environ.lvl1 + engagement.lvl1 | cntry) + age + gender + 
## H5.mod3:     educ + resid + occup + environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1
## H5.mod4: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod4:     engagement.lvl1 + environ.lvl1:engagement.lvl1
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.mod3   33 84674 84955 -42304    84608                        
## H5.mod4   41 84669 85018 -42294    84587 20.209  8   0.009573 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod4<-getFE(H5.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.03       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.10       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.09       0.08
## environ.lvl1                                                0.15       0.01
## engagement.lvl1                                             0.04       0.01
## environ.lvl1:engagement.lvl1                                0.01       0.01
##                                                               df t.value     p
## (Intercept)                                                49.42    0.16 0.872
## age                                                     35128.59   -4.80 0.000
## gender                                                  36655.00    8.25 0.000
## educ                                                    36800.74   12.31 0.000
## resid                                                   36702.20   -5.59 0.000
## occupClerical support workers                           36570.94   -0.21 0.837
## occupCraft and related trades workers                   36586.89   -0.73 0.466
## occupElementary occupations                             36589.34    0.41 0.681
## occupManagers                                           36575.31    0.67 0.505
## occupOther: Not in paid work                            36729.08    1.96 0.050
## occupPlant and machine operators, and assemblers        36589.17   -0.27 0.790
## occupProfessionals                                      36577.84    1.64 0.100
## occupRetired                                            36587.02    0.47 0.635
## occupService and sales workers                          36582.19   -0.19 0.850
## occupSkilled agricultural, forestry and fishery workers 36595.06   -0.33 0.742
## occupTechnicians and associate professionals            36571.47    0.28 0.779
## occupUnemployed                                         36589.65   -1.12 0.263
## environ.lvl1                                               20.04    9.98 0.000
## engagement.lvl1                                            22.78    4.81 0.000
## environ.lvl1:engagement.lvl1                               22.52    1.51 0.144
##                                                            LL    UL
## (Intercept)                                             -0.19  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.10  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.15  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.07
## environ.lvl1                                             0.12  0.18
## engagement.lvl1                                          0.02  0.06
## environ.lvl1:engagement.lvl1                            -0.01  0.03
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
## 11        cntry                  (Intercept)                         <NA>
## 12        cntry                 environ.lvl1                         <NA>
## 13        cntry              engagement.lvl1                         <NA>
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>
## 15        cntry                  (Intercept)                 environ.lvl1
## 16        cntry                  (Intercept)              engagement.lvl1
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1
## 18        cntry                 environ.lvl1              engagement.lvl1
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 21     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.25514592  6.509944e-02
## 2   0.05784974  3.346593e-03
## 3   0.05891045  3.470441e-03
## 4   0.02028314  4.114056e-04
## 5   0.28713069  4.238085e-03
## 6   0.30779857  4.626447e-03
## 7   0.09933381  5.140683e-04
## 8   0.04869774  1.659597e-04
## 9  -0.91303438 -1.071331e-03
## 10 -0.07579976 -9.057227e-05
## 11  0.35884243  1.287679e-01
## 12  0.05861115  3.435267e-03
## 13  0.02355699  5.549320e-04
## 14  0.03094585  9.576457e-04
## 15  0.29721321  6.251038e-03
## 16  0.34607903  2.925492e-03
## 17  0.31759929  3.526840e-03
## 18  0.11999810  1.656817e-04
## 19  0.99746393  1.809172e-03
## 20  0.06344070  4.624772e-05
## 21  0.75202320  5.655389e-01
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
##                                                   0.33928 
##                     voting.group.environ.lvl1.(Intercept) 
##                                                   0.02209 
##                  voting.group.engagement.lvl1.(Intercept) 
##                                                   0.02411 
##     voting.group.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                   0.00268 
##                                 voting.group.environ.lvl1 
##                                                   0.07369 
##                 voting.group.engagement.lvl1.environ.lvl1 
##                                                  -0.00325 
##    voting.group.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                  -0.02651 
##                              voting.group.engagement.lvl1 
##                                                   0.07446 
## voting.group.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                  -0.00417 
##                 voting.group.environ.lvl1:engagement.lvl1 
##                                                   0.00004 
##                                         cntry.(Intercept) 
##                                                   0.47717 
##                            cntry.environ.lvl1.(Intercept) 
##                                                   0.02316 
##                         cntry.engagement.lvl1.(Intercept) 
##                                                   0.01084 
##            cntry.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                   0.01307 
##                                        cntry.environ.lvl1 
##                                                   0.07442 
##                        cntry.engagement.lvl1.environ.lvl1 
##                                                   0.00056 
##           cntry.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                   0.03892 
##                                     cntry.engagement.lvl1 
##                                                   0.02938 
##        cntry.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                  -0.00278 
##                        cntry.environ.lvl1:engagement.lvl1 
##                                                   0.00000
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and engagement


```r
H5.mod5<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|voting.group)+
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|cntry)+
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
## H5.mod4:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod4:     engagement.lvl1 + environ.lvl1:engagement.lvl1
## H5.mod5: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.mod4   41 84669 85018 -42294    84587                         
## H5.mod5   68 84497 85076 -42180    84361 226.78 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod5<-getFE(H5.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.44       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.10       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.09       0.08
## environ.lvl1                                                0.12       0.02
## engagement.lvl1                                             0.01       0.02
## all.parties.lvl2Did not vote                                0.36       0.05
## all.parties.lvl2Don't know                                  0.34       0.06
## all.parties.lvl2Invalid vote                                0.29       0.26
## all.parties.lvl2NE age                                      0.55       0.06
## all.parties.lvl2NE citizen                                  0.70       0.06
## all.parties.lvl2NE other                                    0.55       0.07
## all.parties.lvl2No answer                                   0.29       0.28
## all.parties.lvl2Other party                                 0.47       0.04
## all.parties.lvl2Pro-environment party                       0.78       0.05
## environ.lvl1:engagement.lvl1                                0.01       0.01
## environ.lvl1:all.parties.lvl2Did not vote                   0.01       0.03
## environ.lvl1:all.parties.lvl2Don't know                     0.02       0.04
## environ.lvl1:all.parties.lvl2Invalid vote                  -0.20       0.29
## environ.lvl1:all.parties.lvl2NE age                         0.10       0.03
## environ.lvl1:all.parties.lvl2NE citizen                     0.00       0.04
## environ.lvl1:all.parties.lvl2NE other                       0.00       0.06
## environ.lvl1:all.parties.lvl2No answer                     -0.10       0.26
## environ.lvl1:all.parties.lvl2Other party                    0.04       0.02
## environ.lvl1:all.parties.lvl2Pro-environment party          0.04       0.03
## engagement.lvl1:all.parties.lvl2Did not vote                0.00       0.03
## engagement.lvl1:all.parties.lvl2Don't know                  0.00       0.04
## engagement.lvl1:all.parties.lvl2Invalid vote               -0.11       0.33
## engagement.lvl1:all.parties.lvl2NE age                      0.00       0.03
## engagement.lvl1:all.parties.lvl2NE citizen                 -0.04       0.04
## engagement.lvl1:all.parties.lvl2NE other                    0.06       0.06
## engagement.lvl1:all.parties.lvl2No answer                  -0.75       0.49
## engagement.lvl1:all.parties.lvl2Other party                 0.05       0.02
## engagement.lvl1:all.parties.lvl2Pro-environment party       0.10       0.03
##                                                               df t.value     p
## (Intercept)                                                66.41   -4.18 0.000
## age                                                     36680.85   -4.56 0.000
## gender                                                  36684.44    8.20 0.000
## educ                                                    36778.09   12.28 0.000
## resid                                                   36765.30   -5.48 0.000
## occupClerical support workers                           36632.78   -0.21 0.832
## occupCraft and related trades workers                   36650.28   -0.72 0.469
## occupElementary occupations                             36646.11    0.39 0.698
## occupManagers                                           36640.88    0.67 0.500
## occupOther: Not in paid work                            36687.19    1.86 0.063
## occupPlant and machine operators, and assemblers        36653.44   -0.28 0.783
## occupProfessionals                                      36640.66    1.64 0.102
## occupRetired                                            36639.02    0.43 0.668
## occupService and sales workers                          36642.92   -0.19 0.851
## occupSkilled agricultural, forestry and fishery workers 36663.59   -0.36 0.723
## occupTechnicians and associate professionals            36635.02    0.29 0.772
## occupUnemployed                                         36630.36   -1.19 0.236
## environ.lvl1                                               90.90    4.91 0.000
## engagement.lvl1                                           164.27    0.54 0.589
## all.parties.lvl2Did not vote                              194.57    7.11 0.000
## all.parties.lvl2Don't know                                276.47    5.99 0.000
## all.parties.lvl2Invalid vote                              799.52    1.11 0.268
## all.parties.lvl2NE age                                    278.14    9.89 0.000
## all.parties.lvl2NE citizen                                280.96   11.65 0.000
## all.parties.lvl2NE other                                  617.63    7.36 0.000
## all.parties.lvl2No answer                                1043.04    1.04 0.297
## all.parties.lvl2Other party                               248.02   12.19 0.000
## all.parties.lvl2Pro-environment party                     270.24   15.20 0.000
## environ.lvl1:engagement.lvl1                               23.73    1.55 0.134
## environ.lvl1:all.parties.lvl2Did not vote                 102.99    0.28 0.778
## environ.lvl1:all.parties.lvl2Don't know                   391.71    0.53 0.593
## environ.lvl1:all.parties.lvl2Invalid vote               10065.22   -0.70 0.483
## environ.lvl1:all.parties.lvl2NE age                       267.92    3.03 0.003
## environ.lvl1:all.parties.lvl2NE citizen                   248.16   -0.02 0.987
## environ.lvl1:all.parties.lvl2NE other                    1575.65    0.04 0.969
## environ.lvl1:all.parties.lvl2No answer                   8286.55   -0.40 0.692
## environ.lvl1:all.parties.lvl2Other party                  165.30    1.79 0.076
## environ.lvl1:all.parties.lvl2Pro-environment party        256.43    1.25 0.212
## engagement.lvl1:all.parties.lvl2Did not vote              103.64    0.18 0.855
## engagement.lvl1:all.parties.lvl2Don't know                425.45    0.00 0.998
## engagement.lvl1:all.parties.lvl2Invalid vote            13148.96   -0.35 0.730
## engagement.lvl1:all.parties.lvl2NE age                    297.94    0.12 0.901
## engagement.lvl1:all.parties.lvl2NE citizen                229.51   -1.09 0.276
## engagement.lvl1:all.parties.lvl2NE other                 1823.75    0.89 0.372
## engagement.lvl1:all.parties.lvl2No answer               29067.55   -1.53 0.127
## engagement.lvl1:all.parties.lvl2Other party               164.70    2.11 0.036
## engagement.lvl1:all.parties.lvl2Pro-environment party     264.95    2.98 0.003
##                                                            LL    UL
## (Intercept)                                             -0.65 -0.23
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                            -0.01  0.25
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.06
## environ.lvl1                                             0.07  0.16
## engagement.lvl1                                         -0.03  0.05
## all.parties.lvl2Did not vote                             0.26  0.46
## all.parties.lvl2Don't know                               0.23  0.45
## all.parties.lvl2Invalid vote                            -0.22  0.80
## all.parties.lvl2NE age                                   0.44  0.66
## all.parties.lvl2NE citizen                               0.58  0.82
## all.parties.lvl2NE other                                 0.40  0.69
## all.parties.lvl2No answer                               -0.25  0.83
## all.parties.lvl2Other party                              0.39  0.55
## all.parties.lvl2Pro-environment party                    0.68  0.88
## environ.lvl1:engagement.lvl1                             0.00  0.03
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.06
## environ.lvl1:all.parties.lvl2Don't know                 -0.05  0.09
## environ.lvl1:all.parties.lvl2Invalid vote               -0.77  0.36
## environ.lvl1:all.parties.lvl2NE age                      0.04  0.17
## environ.lvl1:all.parties.lvl2NE citizen                 -0.07  0.07
## environ.lvl1:all.parties.lvl2NE other                   -0.12  0.13
## environ.lvl1:all.parties.lvl2No answer                  -0.62  0.41
## environ.lvl1:all.parties.lvl2Other party                 0.00  0.08
## environ.lvl1:all.parties.lvl2Pro-environment party      -0.02  0.11
## engagement.lvl1:all.parties.lvl2Did not vote            -0.05  0.06
## engagement.lvl1:all.parties.lvl2Don't know              -0.08  0.08
## engagement.lvl1:all.parties.lvl2Invalid vote            -0.76  0.53
## engagement.lvl1:all.parties.lvl2NE age                  -0.06  0.07
## engagement.lvl1:all.parties.lvl2NE citizen              -0.11  0.03
## engagement.lvl1:all.parties.lvl2NE other                -0.07  0.19
## engagement.lvl1:all.parties.lvl2No answer               -1.70  0.21
## engagement.lvl1:all.parties.lvl2Other party              0.00  0.09
## engagement.lvl1:all.parties.lvl2Pro-environment party    0.03  0.16
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
## 11        cntry                  (Intercept)                         <NA>
## 12        cntry                 environ.lvl1                         <NA>
## 13        cntry              engagement.lvl1                         <NA>
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>
## 15        cntry                  (Intercept)                 environ.lvl1
## 16        cntry                  (Intercept)              engagement.lvl1
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1
## 18        cntry                 environ.lvl1              engagement.lvl1
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 21     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.16355980  2.675181e-02
## 2   0.05176991  2.680123e-03
## 3   0.05004058  2.504060e-03
## 4   0.01728683  2.988346e-04
## 5   0.21275829  1.801526e-03
## 6   0.17964390  1.470318e-03
## 7  -0.54201440 -1.532508e-03
## 8  -0.10612025 -2.749147e-04
## 9  -0.88435933 -7.914466e-04
## 10  0.30536185  2.641512e-04
## 11  0.34137410  1.165363e-01
## 12  0.05954648  3.545783e-03
## 13  0.02740233  7.508876e-04
## 14  0.02945741  8.677392e-04
## 15  0.30243163  6.147717e-03
## 16  0.26596623  2.487966e-03
## 17  0.35174864  3.537184e-03
## 18  0.19446665  3.173136e-04
## 19  0.98378475  1.725642e-03
## 20  0.04097690  3.307662e-05
## 21  0.75199104  5.654905e-01
```

\newpage

#### Look among which voting group there is strongest association between engagement and refugee attitudes


```r
H5.mod5.trends<-emtrends(H5.mod5,specs = c("all.parties.lvl2"),var=c("engagement.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 36876' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 36876)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 36876' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 36876)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.mod5.trends.tab<-data.frame(H5.mod5.trends))
```

```
##          all.parties.lvl2 engagement.lvl1.trend         SE  df   asymp.LCL
## 1  Anti-immigration party            0.01142262 0.02108929 Inf -0.02991163
## 2            Did not vote            0.01615834 0.01753044 Inf -0.01820070
## 3              Don't know            0.01133892 0.03360036 Inf -0.05451657
## 4            Invalid vote           -0.10251447 0.32929852 Inf -0.74792770
## 5                  NE age            0.01573999 0.02883423 Inf -0.04077406
## 6              NE citizen           -0.02869090 0.03135627 Inf -0.09014806
## 7                NE other            0.06935223 0.06209503 Inf -0.05235179
## 8               No answer           -0.73442322 0.48783159 Inf -1.69055557
## 9             Other party            0.05863312 0.01162927 Inf  0.03584018
## 10  Pro-environment party            0.10906114 0.02656268 Inf  0.05699925
##     asymp.UCL
## 1  0.05275688
## 2  0.05051737
## 3  0.07719442
## 4  0.54289877
## 5  0.07225404
## 6  0.03276626
## 7  0.19105626
## 8  0.22170913
## 9  0.08142606
## 10 0.16112303
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
##                     group engagement.lvl1.trend   SE      p adj.p asymp.LCL
## 1  Anti-immigration party                  0.01 0.02 0.5881 1e+00     -0.03
## 2            Did not vote                  0.02 0.02 0.3567 1e+00     -0.02
## 3              Don't know                  0.01 0.03 0.7358 1e+00     -0.05
## 4            Invalid vote                 -0.10 0.33 0.7556 1e+00     -0.75
## 5                  NE age                  0.02 0.03 0.5851 1e+00     -0.04
## 6              NE citizen                 -0.03 0.03 0.3602 1e+00     -0.09
## 7                NE other                  0.07 0.06 0.2640 1e+00     -0.05
## 8               No answer                 -0.73 0.49 0.1322 1e+00     -1.69
## 9             Other party                  0.06 0.01 0.0000 0e+00      0.04
## 10  Pro-environment party                  0.11 0.03 0.0000 4e-04      0.06
##    asymp.UCL
## 1       0.05
## 2       0.05
## 3       0.08
## 4       0.54
## 5       0.07
## 6       0.03
## 7       0.19
## 8       0.22
## 9       0.08
## 10      0.16
```

```r
write.csv2(H5.mod5.trends.tab,"H5.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.09763852 0.03274783 Inf
##    z.ratio     p.value
## 1 2.981526 0.002868157
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect   0.0765 0.0691 Inf  1.107  0.2683 
##  Did not vote effect             0.0817 0.0681 Inf  1.200  0.2303 
##  Don't know effect               0.0764 0.0738 Inf  1.034  0.3009 
##  Invalid vote effect            -0.0501 0.3338 Inf -0.150  0.8806 
##  NE age effect                   0.0813 0.0718 Inf  1.131  0.2580 
##  NE citizen effect               0.0319 0.0728 Inf  0.438  0.6615 
##  NE other effect                 0.1408 0.0903 Inf  1.560  0.1187 
##  No answer effect               -0.7523 0.4893 Inf -1.538  0.1242 
##  Other party effect              0.1289 0.0668 Inf  1.929  0.0537 
##  Pro-environment party effect    0.1849 0.0710 Inf  2.606  0.0092 
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
## 1           Other party - Anti-immigration party 0.04721050 0.02237756 Inf
## 2 Pro-environment party - Anti-immigration party 0.09763852 0.03274783 Inf
## 3            Pro-environment party - Other party 0.05042802 0.02756586 Inf
##    z.ratio     p.value
## 1 2.109725 0.034882076
## 2 2.981526 0.002868157
## 3 1.829365 0.067344979
```

\newpage

### Model 6: Enter three-way interaction voting group x engagement x environment attitudes


```r
H5.mod6<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|voting.group)+
                (environ.lvl1+engagement.lvl1+environ.lvl1:engagement.lvl1|cntry)+
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
## H5.mod5:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
## H5.mod6: refugees ~ (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod6:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod6:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod6:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2 + 
## H5.mod6:     environ.lvl1:engagement.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   68 84497 85076 -42180    84361                       
## H5.mod6   77 84497 85153 -42171    84343 17.602  9    0.04008 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                                    Estimate
## (Intercept)                                                           -0.43
## age                                                                    0.00
## gender                                                                 0.07
## educ                                                                   0.02
## resid                                                                 -0.05
## occupClerical support workers                                         -0.01
## occupCraft and related trades workers                                 -0.05
## occupElementary occupations                                            0.02
## occupManagers                                                          0.04
## occupOther: Not in paid work                                           0.12
## occupPlant and machine operators, and assemblers                      -0.02
## occupProfessionals                                                     0.10
## occupRetired                                                           0.03
## occupService and sales workers                                        -0.01
## occupSkilled agricultural, forestry and fishery workers               -0.02
## occupTechnicians and associate professionals                           0.02
## occupUnemployed                                                       -0.09
## environ.lvl1                                                           0.12
## engagement.lvl1                                                        0.01
## all.parties.lvl2Did not vote                                           0.36
## all.parties.lvl2Don't know                                             0.33
## all.parties.lvl2Invalid vote                                           0.28
## all.parties.lvl2NE age                                                 0.55
## all.parties.lvl2NE citizen                                             0.68
## all.parties.lvl2NE other                                               0.55
## all.parties.lvl2No answer                                              0.30
## all.parties.lvl2Other party                                            0.47
## all.parties.lvl2Pro-environment party                                  0.78
## environ.lvl1:engagement.lvl1                                          -0.01
## environ.lvl1:all.parties.lvl2Did not vote                              0.00
## environ.lvl1:all.parties.lvl2Don't know                                0.02
## environ.lvl1:all.parties.lvl2Invalid vote                             -0.24
## environ.lvl1:all.parties.lvl2NE age                                    0.10
## environ.lvl1:all.parties.lvl2NE citizen                               -0.01
## environ.lvl1:all.parties.lvl2NE other                                  0.01
## environ.lvl1:all.parties.lvl2No answer                                -0.01
## environ.lvl1:all.parties.lvl2Other party                               0.04
## environ.lvl1:all.parties.lvl2Pro-environment party                     0.04
## engagement.lvl1:all.parties.lvl2Did not vote                           0.01
## engagement.lvl1:all.parties.lvl2Don't know                             0.00
## engagement.lvl1:all.parties.lvl2Invalid vote                          -0.21
## engagement.lvl1:all.parties.lvl2NE age                                 0.01
## engagement.lvl1:all.parties.lvl2NE citizen                            -0.04
## engagement.lvl1:all.parties.lvl2NE other                               0.06
## engagement.lvl1:all.parties.lvl2No answer                             -0.64
## engagement.lvl1:all.parties.lvl2Other party                            0.05
## engagement.lvl1:all.parties.lvl2Pro-environment party                  0.10
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote              0.03
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know                0.03
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote             -0.15
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                    0.00
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen                0.14
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                 -0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer                -0.53
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party               0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party     0.00
##                                                                    Std..Error
## (Intercept)                                                              0.10
## age                                                                      0.00
## gender                                                                   0.01
## educ                                                                     0.00
## resid                                                                    0.01
## occupClerical support workers                                            0.06
## occupCraft and related trades workers                                    0.06
## occupElementary occupations                                              0.06
## occupManagers                                                            0.06
## occupOther: Not in paid work                                             0.07
## occupPlant and machine operators, and assemblers                         0.06
## occupProfessionals                                                       0.06
## occupRetired                                                             0.07
## occupService and sales workers                                           0.06
## occupSkilled agricultural, forestry and fishery workers                  0.07
## occupTechnicians and associate professionals                             0.06
## occupUnemployed                                                          0.08
## environ.lvl1                                                             0.02
## engagement.lvl1                                                          0.02
## all.parties.lvl2Did not vote                                             0.05
## all.parties.lvl2Don't know                                               0.06
## all.parties.lvl2Invalid vote                                             0.26
## all.parties.lvl2NE age                                                   0.06
## all.parties.lvl2NE citizen                                               0.06
## all.parties.lvl2NE other                                                 0.07
## all.parties.lvl2No answer                                                0.28
## all.parties.lvl2Other party                                              0.04
## all.parties.lvl2Pro-environment party                                    0.05
## environ.lvl1:engagement.lvl1                                             0.02
## environ.lvl1:all.parties.lvl2Did not vote                                0.03
## environ.lvl1:all.parties.lvl2Don't know                                  0.04
## environ.lvl1:all.parties.lvl2Invalid vote                                0.31
## environ.lvl1:all.parties.lvl2NE age                                      0.03
## environ.lvl1:all.parties.lvl2NE citizen                                  0.04
## environ.lvl1:all.parties.lvl2NE other                                    0.06
## environ.lvl1:all.parties.lvl2No answer                                   0.28
## environ.lvl1:all.parties.lvl2Other party                                 0.02
## environ.lvl1:all.parties.lvl2Pro-environment party                       0.03
## engagement.lvl1:all.parties.lvl2Did not vote                             0.03
## engagement.lvl1:all.parties.lvl2Don't know                               0.04
## engagement.lvl1:all.parties.lvl2Invalid vote                             0.40
## engagement.lvl1:all.parties.lvl2NE age                                   0.03
## engagement.lvl1:all.parties.lvl2NE citizen                               0.04
## engagement.lvl1:all.parties.lvl2NE other                                 0.06
## engagement.lvl1:all.parties.lvl2No answer                                0.50
## engagement.lvl1:all.parties.lvl2Other party                              0.02
## engagement.lvl1:all.parties.lvl2Pro-environment party                    0.03
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote                0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know                  0.04
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote                0.43
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                      0.04
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen                  0.04
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                    0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer                   0.62
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party                 0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party       0.04
##                                                                          df
## (Intercept)                                                           66.36
## age                                                                36684.64
## gender                                                             36683.58
## educ                                                               36778.42
## resid                                                              36766.59
## occupClerical support workers                                      36634.73
## occupCraft and related trades workers                              36651.04
## occupElementary occupations                                        36646.91
## occupManagers                                                      36641.88
## occupOther: Not in paid work                                       36685.96
## occupPlant and machine operators, and assemblers                   36654.28
## occupProfessionals                                                 36641.77
## occupRetired                                                       36639.19
## occupService and sales workers                                     36643.91
## occupSkilled agricultural, forestry and fishery workers            36664.39
## occupTechnicians and associate professionals                       36636.69
## occupUnemployed                                                    36632.25
## environ.lvl1                                                          91.91
## engagement.lvl1                                                      164.72
## all.parties.lvl2Did not vote                                         200.34
## all.parties.lvl2Don't know                                           275.65
## all.parties.lvl2Invalid vote                                         806.19
## all.parties.lvl2NE age                                               278.73
## all.parties.lvl2NE citizen                                           281.49
## all.parties.lvl2NE other                                             614.15
## all.parties.lvl2No answer                                           1043.15
## all.parties.lvl2Other party                                          248.22
## all.parties.lvl2Pro-environment party                                271.95
## environ.lvl1:engagement.lvl1                                         414.95
## environ.lvl1:all.parties.lvl2Did not vote                            104.24
## environ.lvl1:all.parties.lvl2Don't know                              390.60
## environ.lvl1:all.parties.lvl2Invalid vote                          13166.97
## environ.lvl1:all.parties.lvl2NE age                                  265.81
## environ.lvl1:all.parties.lvl2NE citizen                              246.05
## environ.lvl1:all.parties.lvl2NE other                               1540.19
## environ.lvl1:all.parties.lvl2No answer                             10865.61
## environ.lvl1:all.parties.lvl2Other party                             165.74
## environ.lvl1:all.parties.lvl2Pro-environment party                   255.52
## engagement.lvl1:all.parties.lvl2Did not vote                         103.98
## engagement.lvl1:all.parties.lvl2Don't know                           424.30
## engagement.lvl1:all.parties.lvl2Invalid vote                       23542.13
## engagement.lvl1:all.parties.lvl2NE age                               298.07
## engagement.lvl1:all.parties.lvl2NE citizen                           229.16
## engagement.lvl1:all.parties.lvl2NE other                            1859.27
## engagement.lvl1:all.parties.lvl2No answer                          30398.91
## engagement.lvl1:all.parties.lvl2Other party                          165.55
## engagement.lvl1:all.parties.lvl2Pro-environment party                264.56
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote            777.50
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know             5133.68
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote          35036.90
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                 3547.26
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen             2314.10
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other              11828.20
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer             35971.42
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party            1444.52
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party  3102.31
##                                                                    t.value
## (Intercept)                                                          -4.15
## age                                                                  -4.54
## gender                                                                8.22
## educ                                                                 12.28
## resid                                                                -5.47
## occupClerical support workers                                        -0.22
## occupCraft and related trades workers                                -0.73
## occupElementary occupations                                           0.39
## occupManagers                                                         0.67
## occupOther: Not in paid work                                          1.86
## occupPlant and machine operators, and assemblers                     -0.28
## occupProfessionals                                                    1.62
## occupRetired                                                          0.42
## occupService and sales workers                                       -0.19
## occupSkilled agricultural, forestry and fishery workers              -0.35
## occupTechnicians and associate professionals                          0.28
## occupUnemployed                                                      -1.19
## environ.lvl1                                                          5.00
## engagement.lvl1                                                       0.46
## all.parties.lvl2Did not vote                                          6.94
## all.parties.lvl2Don't know                                            5.90
## all.parties.lvl2Invalid vote                                          1.08
## all.parties.lvl2NE age                                                9.84
## all.parties.lvl2NE citizen                                           11.26
## all.parties.lvl2NE other                                              7.33
## all.parties.lvl2No answer                                             1.08
## all.parties.lvl2Other party                                          12.08
## all.parties.lvl2Pro-environment party                                15.12
## environ.lvl1:engagement.lvl1                                         -0.38
## environ.lvl1:all.parties.lvl2Did not vote                             0.14
## environ.lvl1:all.parties.lvl2Don't know                               0.48
## environ.lvl1:all.parties.lvl2Invalid vote                            -0.80
## environ.lvl1:all.parties.lvl2NE age                                   3.00
## environ.lvl1:all.parties.lvl2NE citizen                              -0.23
## environ.lvl1:all.parties.lvl2NE other                                 0.15
## environ.lvl1:all.parties.lvl2No answer                               -0.03
## environ.lvl1:all.parties.lvl2Other party                              1.71
## environ.lvl1:all.parties.lvl2Pro-environment party                    1.23
## engagement.lvl1:all.parties.lvl2Did not vote                          0.27
## engagement.lvl1:all.parties.lvl2Don't know                            0.04
## engagement.lvl1:all.parties.lvl2Invalid vote                         -0.51
## engagement.lvl1:all.parties.lvl2NE age                                0.20
## engagement.lvl1:all.parties.lvl2NE citizen                           -1.10
## engagement.lvl1:all.parties.lvl2NE other                              0.99
## engagement.lvl1:all.parties.lvl2No answer                            -1.27
## engagement.lvl1:all.parties.lvl2Other party                           2.17
## engagement.lvl1:all.parties.lvl2Pro-environment party                 3.01
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote             1.22
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know               0.75
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote            -0.35
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                  -0.12
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen               3.54
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other                -0.98
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer               -0.85
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party              0.95
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party   -0.02
##                                                                        p    LL
## (Intercept)                                                        0.000 -0.64
## age                                                                0.000  0.00
## gender                                                             0.000  0.05
## educ                                                               0.000  0.01
## resid                                                              0.000 -0.06
## occupClerical support workers                                      0.828 -0.14
## occupCraft and related trades workers                              0.463 -0.17
## occupElementary occupations                                        0.700 -0.10
## occupManagers                                                      0.504 -0.08
## occupOther: Not in paid work                                       0.063 -0.01
## occupPlant and machine operators, and assemblers                   0.776 -0.14
## occupProfessionals                                                 0.104 -0.02
## occupRetired                                                       0.672 -0.11
## occupService and sales workers                                     0.847 -0.14
## occupSkilled agricultural, forestry and fishery workers            0.725 -0.16
## occupTechnicians and associate professionals                       0.778 -0.11
## occupUnemployed                                                    0.234 -0.24
## environ.lvl1                                                       0.000  0.07
## engagement.lvl1                                                    0.646 -0.03
## all.parties.lvl2Did not vote                                       0.000  0.26
## all.parties.lvl2Don't know                                         0.000  0.22
## all.parties.lvl2Invalid vote                                       0.282 -0.23
## all.parties.lvl2NE age                                             0.000  0.44
## all.parties.lvl2NE citizen                                         0.000  0.56
## all.parties.lvl2NE other                                           0.000  0.40
## all.parties.lvl2No answer                                          0.278 -0.24
## all.parties.lvl2Other party                                        0.000  0.39
## all.parties.lvl2Pro-environment party                              0.000  0.68
## environ.lvl1:engagement.lvl1                                       0.701 -0.05
## environ.lvl1:all.parties.lvl2Did not vote                          0.888 -0.05
## environ.lvl1:all.parties.lvl2Don't know                            0.633 -0.06
## environ.lvl1:all.parties.lvl2Invalid vote                          0.424 -0.84
## environ.lvl1:all.parties.lvl2NE age                                0.003  0.04
## environ.lvl1:all.parties.lvl2NE citizen                            0.817 -0.08
## environ.lvl1:all.parties.lvl2NE other                              0.884 -0.11
## environ.lvl1:all.parties.lvl2No answer                             0.979 -0.56
## environ.lvl1:all.parties.lvl2Other party                           0.090 -0.01
## environ.lvl1:all.parties.lvl2Pro-environment party                 0.219 -0.02
## engagement.lvl1:all.parties.lvl2Did not vote                       0.788 -0.04
## engagement.lvl1:all.parties.lvl2Don't know                         0.966 -0.07
## engagement.lvl1:all.parties.lvl2Invalid vote                       0.607 -1.00
## engagement.lvl1:all.parties.lvl2NE age                             0.843 -0.06
## engagement.lvl1:all.parties.lvl2NE citizen                         0.272 -0.11
## engagement.lvl1:all.parties.lvl2NE other                           0.322 -0.06
## engagement.lvl1:all.parties.lvl2No answer                          0.205 -1.62
## engagement.lvl1:all.parties.lvl2Other party                        0.031  0.00
## engagement.lvl1:all.parties.lvl2Pro-environment party              0.003  0.03
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote          0.221 -0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know            0.456 -0.05
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote          0.723 -1.00
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                0.901 -0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen            0.000  0.06
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other              0.325 -0.22
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer             0.395 -1.74
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party           0.341 -0.02
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party 0.987 -0.07
##                                                                       UL
## (Intercept)                                                        -0.22
## age                                                                 0.00
## gender                                                              0.09
## educ                                                                0.02
## resid                                                              -0.03
## occupClerical support workers                                       0.11
## occupCraft and related trades workers                               0.08
## occupElementary occupations                                         0.15
## occupManagers                                                       0.17
## occupOther: Not in paid work                                        0.25
## occupPlant and machine operators, and assemblers                    0.11
## occupProfessionals                                                  0.23
## occupRetired                                                        0.17
## occupService and sales workers                                      0.11
## occupSkilled agricultural, forestry and fishery workers             0.11
## occupTechnicians and associate professionals                        0.14
## occupUnemployed                                                     0.06
## environ.lvl1                                                        0.16
## engagement.lvl1                                                     0.05
## all.parties.lvl2Did not vote                                        0.46
## all.parties.lvl2Don't know                                          0.44
## all.parties.lvl2Invalid vote                                        0.79
## all.parties.lvl2NE age                                              0.66
## all.parties.lvl2NE citizen                                          0.80
## all.parties.lvl2NE other                                            0.69
## all.parties.lvl2No answer                                           0.85
## all.parties.lvl2Other party                                         0.54
## all.parties.lvl2Pro-environment party                               0.89
## environ.lvl1:engagement.lvl1                                        0.03
## environ.lvl1:all.parties.lvl2Did not vote                           0.05
## environ.lvl1:all.parties.lvl2Don't know                             0.09
## environ.lvl1:all.parties.lvl2Invalid vote                           0.35
## environ.lvl1:all.parties.lvl2NE age                                 0.17
## environ.lvl1:all.parties.lvl2NE citizen                             0.06
## environ.lvl1:all.parties.lvl2NE other                               0.13
## environ.lvl1:all.parties.lvl2No answer                              0.55
## environ.lvl1:all.parties.lvl2Other party                            0.08
## environ.lvl1:all.parties.lvl2Pro-environment party                  0.11
## engagement.lvl1:all.parties.lvl2Did not vote                        0.06
## engagement.lvl1:all.parties.lvl2Don't know                          0.08
## engagement.lvl1:all.parties.lvl2Invalid vote                        0.58
## engagement.lvl1:all.parties.lvl2NE age                              0.07
## engagement.lvl1:all.parties.lvl2NE citizen                          0.03
## engagement.lvl1:all.parties.lvl2NE other                            0.19
## engagement.lvl1:all.parties.lvl2No answer                           0.35
## engagement.lvl1:all.parties.lvl2Other party                         0.09
## engagement.lvl1:all.parties.lvl2Pro-environment party               0.16
## environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote           0.08
## environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know             0.12
## environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote           0.69
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE age                 0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen             0.21
## environ.lvl1:engagement.lvl1:all.parties.lvl2NE other               0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2No answer              0.69
## environ.lvl1:engagement.lvl1:all.parties.lvl2Other party            0.07
## environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party  0.07
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
## 11        cntry                  (Intercept)                         <NA>
## 12        cntry                 environ.lvl1                         <NA>
## 13        cntry              engagement.lvl1                         <NA>
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>
## 15        cntry                  (Intercept)                 environ.lvl1
## 16        cntry                  (Intercept)              engagement.lvl1
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1
## 18        cntry                 environ.lvl1              engagement.lvl1
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1
## 21     Residual                         <NA>                         <NA>
##         est_SD       est_SD2
## 1   0.16364275  2.677895e-02
## 2   0.05138266  2.640178e-03
## 3   0.04983214  2.483242e-03
## 4   0.01350661  1.824286e-04
## 5   0.21368560  1.796754e-03
## 6   0.18425552  1.502543e-03
## 7  -0.67819859 -1.498995e-03
## 8  -0.09495344 -2.431290e-04
## 9  -0.84444019 -5.860464e-04
## 10  0.13432092  9.040649e-05
## 11  0.34173977  1.167861e-01
## 12  0.05931879  3.518719e-03
## 13  0.02767161  7.657182e-04
## 14  0.03011049  9.066416e-04
## 15  0.30069673  6.095601e-03
## 16  0.26567119  2.512317e-03
## 17  0.33983873  3.496924e-03
## 18  0.19889010  3.264675e-04
## 19  0.99013518  1.768498e-03
## 20  0.07933824  6.610508e-05
## 21  0.75183742  5.652595e-01
```

#### Refit with manually coded level-1 interaction


```r
dat$env.eng.int<-dat$environ.lvl1*dat$engagement.lvl1

H5.mod6<-lmer(refugees~
                (environ.lvl1+engagement.lvl1+env.eng.int|voting.group)+
                (environ.lvl1+engagement.lvl1+env.eng.int|cntry)+
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
## H5.mod5:     voting.group) + (environ.lvl1 + engagement.lvl1 + environ.lvl1:engagement.lvl1 | 
## H5.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.mod5:     engagement.lvl1 + environ.lvl1:engagement.lvl1 + all.parties.lvl2 + 
## H5.mod5:     environ.lvl1:all.parties.lvl2 + engagement.lvl1:all.parties.lvl2
## H5.mod6: refugees ~ (environ.lvl1 + engagement.lvl1 + env.eng.int | voting.group) + 
## H5.mod6:     (environ.lvl1 + engagement.lvl1 + env.eng.int | cntry) + 
## H5.mod6:     age + gender + educ + resid + occup + environ.lvl1 + engagement.lvl1 + 
## H5.mod6:     env.eng.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.mod6:     engagement.lvl1:all.parties.lvl2 + env.eng.int:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   68 84497 85076 -42180    84361                       
## H5.mod6   77 84497 85153 -42171    84343 17.602  9    0.04008 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.43       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.10       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.09       0.08
## environ.lvl1                                                0.12       0.02
## engagement.lvl1                                             0.01       0.02
## env.eng.int                                                -0.01       0.02
## all.parties.lvl2Did not vote                                0.36       0.05
## all.parties.lvl2Don't know                                  0.33       0.06
## all.parties.lvl2Invalid vote                                0.28       0.26
## all.parties.lvl2NE age                                      0.55       0.06
## all.parties.lvl2NE citizen                                  0.68       0.06
## all.parties.lvl2NE other                                    0.55       0.07
## all.parties.lvl2No answer                                   0.30       0.28
## all.parties.lvl2Other party                                 0.47       0.04
## all.parties.lvl2Pro-environment party                       0.78       0.05
## environ.lvl1:all.parties.lvl2Did not vote                   0.00       0.03
## environ.lvl1:all.parties.lvl2Don't know                     0.02       0.04
## environ.lvl1:all.parties.lvl2Invalid vote                  -0.24       0.31
## environ.lvl1:all.parties.lvl2NE age                         0.10       0.03
## environ.lvl1:all.parties.lvl2NE citizen                    -0.01       0.04
## environ.lvl1:all.parties.lvl2NE other                       0.01       0.06
## environ.lvl1:all.parties.lvl2No answer                     -0.01       0.28
## environ.lvl1:all.parties.lvl2Other party                    0.04       0.02
## environ.lvl1:all.parties.lvl2Pro-environment party          0.04       0.03
## engagement.lvl1:all.parties.lvl2Did not vote                0.01       0.03
## engagement.lvl1:all.parties.lvl2Don't know                  0.00       0.04
## engagement.lvl1:all.parties.lvl2Invalid vote               -0.21       0.40
## engagement.lvl1:all.parties.lvl2NE age                      0.01       0.03
## engagement.lvl1:all.parties.lvl2NE citizen                 -0.04       0.04
## engagement.lvl1:all.parties.lvl2NE other                    0.06       0.06
## engagement.lvl1:all.parties.lvl2No answer                  -0.64       0.50
## engagement.lvl1:all.parties.lvl2Other party                 0.05       0.02
## engagement.lvl1:all.parties.lvl2Pro-environment party       0.10       0.03
## env.eng.int:all.parties.lvl2Did not vote                    0.03       0.02
## env.eng.int:all.parties.lvl2Don't know                      0.03       0.04
## env.eng.int:all.parties.lvl2Invalid vote                   -0.15       0.43
## env.eng.int:all.parties.lvl2NE age                          0.00       0.04
## env.eng.int:all.parties.lvl2NE citizen                      0.14       0.04
## env.eng.int:all.parties.lvl2NE other                       -0.07       0.07
## env.eng.int:all.parties.lvl2No answer                      -0.53       0.62
## env.eng.int:all.parties.lvl2Other party                     0.02       0.02
## env.eng.int:all.parties.lvl2Pro-environment party           0.00       0.04
##                                                               df t.value     p
## (Intercept)                                                66.35   -4.15 0.000
## age                                                     36684.64   -4.54 0.000
## gender                                                  36683.58    8.22 0.000
## educ                                                    36778.42   12.28 0.000
## resid                                                   36766.59   -5.47 0.000
## occupClerical support workers                           36634.73   -0.22 0.828
## occupCraft and related trades workers                   36651.04   -0.73 0.463
## occupElementary occupations                             36646.91    0.39 0.700
## occupManagers                                           36641.88    0.67 0.504
## occupOther: Not in paid work                            36685.96    1.86 0.063
## occupPlant and machine operators, and assemblers        36654.28   -0.28 0.776
## occupProfessionals                                      36641.77    1.62 0.104
## occupRetired                                            36639.19    0.42 0.672
## occupService and sales workers                          36643.92   -0.19 0.847
## occupSkilled agricultural, forestry and fishery workers 36664.39   -0.35 0.725
## occupTechnicians and associate professionals            36636.70    0.28 0.778
## occupUnemployed                                         36632.25   -1.19 0.234
## environ.lvl1                                               91.91    5.00 0.000
## engagement.lvl1                                           164.72    0.46 0.646
## env.eng.int                                               414.93   -0.38 0.701
## all.parties.lvl2Did not vote                              200.34    6.94 0.000
## all.parties.lvl2Don't know                                275.65    5.90 0.000
## all.parties.lvl2Invalid vote                              806.19    1.08 0.282
## all.parties.lvl2NE age                                    278.73    9.84 0.000
## all.parties.lvl2NE citizen                                281.49   11.26 0.000
## all.parties.lvl2NE other                                  614.15    7.33 0.000
## all.parties.lvl2No answer                                1043.15    1.08 0.278
## all.parties.lvl2Other party                               248.22   12.08 0.000
## all.parties.lvl2Pro-environment party                     271.95   15.12 0.000
## environ.lvl1:all.parties.lvl2Did not vote                 104.24    0.14 0.888
## environ.lvl1:all.parties.lvl2Don't know                   390.60    0.48 0.633
## environ.lvl1:all.parties.lvl2Invalid vote               13166.93   -0.80 0.424
## environ.lvl1:all.parties.lvl2NE age                       265.81    3.00 0.003
## environ.lvl1:all.parties.lvl2NE citizen                   246.05   -0.23 0.817
## environ.lvl1:all.parties.lvl2NE other                    1540.19    0.15 0.884
## environ.lvl1:all.parties.lvl2No answer                  10865.57   -0.03 0.979
## environ.lvl1:all.parties.lvl2Other party                  165.74    1.71 0.090
## environ.lvl1:all.parties.lvl2Pro-environment party        255.51    1.23 0.219
## engagement.lvl1:all.parties.lvl2Did not vote              103.98    0.27 0.788
## engagement.lvl1:all.parties.lvl2Don't know                424.30    0.04 0.966
## engagement.lvl1:all.parties.lvl2Invalid vote            23542.09   -0.51 0.607
## engagement.lvl1:all.parties.lvl2NE age                    298.07    0.20 0.843
## engagement.lvl1:all.parties.lvl2NE citizen                229.16   -1.10 0.272
## engagement.lvl1:all.parties.lvl2NE other                 1859.26    0.99 0.322
## engagement.lvl1:all.parties.lvl2No answer               30398.88   -1.27 0.205
## engagement.lvl1:all.parties.lvl2Other party               165.55    2.17 0.031
## engagement.lvl1:all.parties.lvl2Pro-environment party     264.56    3.01 0.003
## env.eng.int:all.parties.lvl2Did not vote                  777.51    1.22 0.221
## env.eng.int:all.parties.lvl2Don't know                   5133.72    0.75 0.456
## env.eng.int:all.parties.lvl2Invalid vote                35036.92   -0.35 0.723
## env.eng.int:all.parties.lvl2NE age                       3547.31   -0.12 0.901
## env.eng.int:all.parties.lvl2NE citizen                   2314.13    3.54 0.000
## env.eng.int:all.parties.lvl2NE other                    11828.18   -0.98 0.325
## env.eng.int:all.parties.lvl2No answer                   35971.43   -0.85 0.395
## env.eng.int:all.parties.lvl2Other party                  1444.55    0.95 0.341
## env.eng.int:all.parties.lvl2Pro-environment party        3102.37   -0.02 0.987
##                                                            LL    UL
## (Intercept)                                             -0.64 -0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.01  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                            -0.01  0.25
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.06
## environ.lvl1                                             0.07  0.16
## engagement.lvl1                                         -0.03  0.05
## env.eng.int                                             -0.05  0.03
## all.parties.lvl2Did not vote                             0.26  0.46
## all.parties.lvl2Don't know                               0.22  0.44
## all.parties.lvl2Invalid vote                            -0.23  0.79
## all.parties.lvl2NE age                                   0.44  0.66
## all.parties.lvl2NE citizen                               0.56  0.80
## all.parties.lvl2NE other                                 0.40  0.69
## all.parties.lvl2No answer                               -0.24  0.85
## all.parties.lvl2Other party                              0.39  0.54
## all.parties.lvl2Pro-environment party                    0.68  0.89
## environ.lvl1:all.parties.lvl2Did not vote               -0.05  0.05
## environ.lvl1:all.parties.lvl2Don't know                 -0.06  0.09
## environ.lvl1:all.parties.lvl2Invalid vote               -0.84  0.35
## environ.lvl1:all.parties.lvl2NE age                      0.04  0.17
## environ.lvl1:all.parties.lvl2NE citizen                 -0.08  0.06
## environ.lvl1:all.parties.lvl2NE other                   -0.11  0.13
## environ.lvl1:all.parties.lvl2No answer                  -0.56  0.55
## environ.lvl1:all.parties.lvl2Other party                -0.01  0.08
## environ.lvl1:all.parties.lvl2Pro-environment party      -0.02  0.11
## engagement.lvl1:all.parties.lvl2Did not vote            -0.04  0.06
## engagement.lvl1:all.parties.lvl2Don't know              -0.07  0.08
## engagement.lvl1:all.parties.lvl2Invalid vote            -1.00  0.58
## engagement.lvl1:all.parties.lvl2NE age                  -0.06  0.07
## engagement.lvl1:all.parties.lvl2NE citizen              -0.11  0.03
## engagement.lvl1:all.parties.lvl2NE other                -0.06  0.19
## engagement.lvl1:all.parties.lvl2No answer               -1.62  0.35
## engagement.lvl1:all.parties.lvl2Other party              0.00  0.09
## engagement.lvl1:all.parties.lvl2Pro-environment party    0.03  0.16
## env.eng.int:all.parties.lvl2Did not vote                -0.02  0.08
## env.eng.int:all.parties.lvl2Don't know                  -0.05  0.12
## env.eng.int:all.parties.lvl2Invalid vote                -1.00  0.69
## env.eng.int:all.parties.lvl2NE age                      -0.08  0.07
## env.eng.int:all.parties.lvl2NE citizen                   0.06  0.21
## env.eng.int:all.parties.lvl2NE other                    -0.22  0.07
## env.eng.int:all.parties.lvl2No answer                   -1.74  0.69
## env.eng.int:all.parties.lvl2Other party                 -0.02  0.07
## env.eng.int:all.parties.lvl2Pro-environment party       -0.07  0.07
```

```r
(VC.H5.mod6<-getVC(H5.mod6))
```

```
##             grp            var1            var2      est_SD       est_SD2
## 1  voting.group     (Intercept)            <NA>  0.16364269  0.0267789316
## 2  voting.group    environ.lvl1            <NA>  0.05138278  0.0026401899
## 3  voting.group engagement.lvl1            <NA>  0.04983220  0.0024832478
## 4  voting.group     env.eng.int            <NA>  0.01350645  0.0001824243
## 5  voting.group     (Intercept)    environ.lvl1  0.21368602  0.0017967610
## 6  voting.group     (Intercept) engagement.lvl1  0.18425932  0.0015025749
## 7  voting.group     (Intercept)     env.eng.int -0.67820928 -0.0014990001
## 8  voting.group    environ.lvl1 engagement.lvl1 -0.09495130 -0.0002431244
## 9  voting.group    environ.lvl1     env.eng.int -0.84443281 -0.0005860356
## 10 voting.group engagement.lvl1     env.eng.int  0.13431390  0.0000904008
## 11        cntry     (Intercept)            <NA>  0.34175223  0.1167945841
## 12        cntry    environ.lvl1            <NA>  0.05931960  0.0035188149
## 13        cntry engagement.lvl1            <NA>  0.02767161  0.0007657178
## 14        cntry     env.eng.int            <NA>  0.03011101  0.0009066729
## 15        cntry     (Intercept)    environ.lvl1  0.30072003  0.0060963784
## 16        cntry     (Intercept) engagement.lvl1  0.26568452  0.0025125342
## 17        cntry     (Intercept)     env.eng.int  0.33986379  0.0034973698
## 18        cntry    environ.lvl1 engagement.lvl1  0.19890239  0.0003264920
## 19        cntry    environ.lvl1     env.eng.int  0.99013511  0.0017685526
## 20        cntry engagement.lvl1     env.eng.int  0.07935203  0.0000661177
## 21     Residual            <NA>            <NA>  0.75183740  0.5652594822
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.mod6.trends<-emtrends(H5.mod6,specs = c("all.parties.lvl2"),var=c("env.eng.int"))
(H5.mod6.trends.tab<-data.frame(H5.mod6.trends))
```

```
##          all.parties.lvl2 env.eng.int.trend         SE  df    asymp.LCL
## 1  Anti-immigration party      -0.008198955 0.02131248 Inf -0.049970655
## 2            Did not vote       0.021848757 0.01546376 Inf -0.008459653
## 3              Don't know       0.023780915 0.03838105 Inf -0.051444552
## 4            Invalid vote      -0.160948280 0.43045084 Inf -1.004616426
## 5                  NE age      -0.013019676 0.03364955 Inf -0.078971588
## 6              NE citizen       0.127233400 0.03311840 Inf  0.062322532
## 7                NE other      -0.080170391 0.07053638 Inf -0.218419160
## 8               No answer      -0.534388253 0.61813635 Inf -1.745913234
## 9             Other party       0.013148610 0.01174939 Inf -0.009879769
## 10  Pro-environment party      -0.008773747 0.03094462 Inf -0.069424091
##     asymp.UCL
## 1  0.03357274
## 2  0.05215717
## 3  0.09900638
## 4  0.68271987
## 5  0.05293224
## 6  0.19214427
## 7  0.05807838
## 8  0.67713673
## 9  0.03617699
## 10 0.05187660
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
##                     group env.eng.int.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party             -0.01 0.02 0.7005 1.0000     -0.05
## 2            Did not vote              0.02 0.02 0.1577 1.0000     -0.01
## 3              Don't know              0.02 0.04 0.5355 1.0000     -0.05
## 4            Invalid vote             -0.16 0.43 0.7085 1.0000     -1.00
## 5                  NE age             -0.01 0.03 0.6988 1.0000     -0.08
## 6              NE citizen              0.13 0.03 0.0001 0.0012      0.06
## 7                NE other             -0.08 0.07 0.2557 1.0000     -0.22
## 8               No answer             -0.53 0.62 0.3873 1.0000     -1.75
## 9             Other party              0.01 0.01 0.2631 1.0000     -0.01
## 10  Pro-environment party             -0.01 0.03 0.7768 1.0000     -0.07
##    asymp.UCL
## 1       0.03
## 2       0.05
## 3       0.10
## 4       0.68
## 5       0.05
## 6       0.19
## 7       0.06
## 8       0.68
## 9       0.04
## 10      0.05
```

```r
write.csv2(H5.mod6.trends.tab,"H5.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast      estimate         SE  df
## 1 Pro-environment party - Anti-immigration party -0.0005747923 0.03641241 Inf
##       z.ratio   p.value
## 1 -0.01578561 0.9874054
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect   0.0597 0.0868 Inf  0.688  1.0000 
##  Did not vote effect             0.0931 0.0856 Inf  1.088  1.0000 
##  Don't know effect               0.0953 0.0924 Inf  1.031  1.0000 
##  Invalid vote effect            -0.1100 0.4360 Inf -0.252  1.0000 
##  NE age effect                   0.0544 0.0906 Inf  0.600  1.0000 
##  NE citizen effect               0.2102 0.0904 Inf  2.326  0.2000 
##  NE other effect                -0.0202 0.1095 Inf -0.185  1.0000 
##  No answer effect               -0.5249 0.6201 Inf -0.847  1.0000 
##  Other party effect              0.0834 0.0850 Inf  0.982  1.0000 
##  Pro-environment party effect    0.0591 0.0896 Inf  0.659  1.0000 
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
##                                         contrast      estimate         SE  df
## 1           Other party - Anti-immigration party  0.0213475652 0.02242258 Inf
## 2 Pro-environment party - Anti-immigration party -0.0005747923 0.03641241 Inf
## 3            Pro-environment party - Other party -0.0219223575 0.03169840 Inf
##       z.ratio p.value
## 1  0.95205677       1
## 2 -0.01578561       1
## 3 -0.69159194       1
```

