---
title: "Hypothesis 2: Those who voted for pro-environment parties will report higher pro-refugee attitudes than those who voted for anti-immigration parties"
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


\newpage

# Hypothesis 2: Those who voted for pro-environment parties will report higher pro-refugee attitudes than those who voted for anti-immigration parties

### Model 1: random intercepts + covariates (same as in H1)


```r
H2.mod1<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)


(FE.H2.mod1<-getFE(H2.mod1))
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
(VC.H2.mod1<-getVC(H2.mod1))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.4133723 0.1708767
## 2     Residual        <NA> <NA> 1.1820745 1.3973001
```

\newpage


### Model 2: Categorical predictor at level-2


```r
H2.mod2<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                all.parties.lvl2
                ,data=dat,REML=F)


(FE.H2.mod2<-getFE(H2.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.75       0.12
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.01       0.11
## occupUnemployed                                            -0.02       0.14
## all.parties.lvl2Did not vote                                0.63       0.09
## all.parties.lvl2Don't know                                  0.58       0.10
## all.parties.lvl2Invalid vote                                0.77       0.49
## all.parties.lvl2NE age                                      1.02       0.10
## all.parties.lvl2NE citizen                                  1.15       0.10
## all.parties.lvl2NE other                                    1.02       0.13
## all.parties.lvl2No answer                                   1.00       0.49
## all.parties.lvl2Other party                                 0.79       0.07
## all.parties.lvl2Pro-environment party                       1.29       0.09
##                                                               df t.value     p
## (Intercept)                                              3229.56   -6.15 0.000
## age                                                     26885.99   -2.28 0.022
## gender                                                  26796.75    6.48 0.000
## educ                                                    26553.32   14.52 0.000
## resid                                                   26880.87   -7.70 0.000
## occupClerical support workers                           26795.71   -0.15 0.882
## occupCraft and related trades workers                   26800.78   -1.24 0.214
## occupElementary occupations                             26805.92   -0.02 0.984
## occupManagers                                           26795.33    0.62 0.538
## occupOther: Not in paid work                            26839.52    1.50 0.132
## occupPlant and machine operators, and assemblers        26805.29   -0.62 0.532
## occupProfessionals                                      26788.45    1.77 0.077
## occupRetired                                            26795.24    0.25 0.804
## occupService and sales workers                          26794.02   -0.60 0.548
## occupSkilled agricultural, forestry and fishery workers 26796.88   -0.18 0.860
## occupTechnicians and associate professionals            26790.15    0.06 0.951
## occupUnemployed                                         26807.62   -0.14 0.886
## all.parties.lvl2Did not vote                              174.91    6.97 0.000
## all.parties.lvl2Don't know                                235.82    5.93 0.000
## all.parties.lvl2Invalid vote                             1594.66    1.56 0.118
## all.parties.lvl2NE age                                    244.96   10.33 0.000
## all.parties.lvl2NE citizen                                239.02   11.41 0.000
## all.parties.lvl2NE other                                  588.95    7.60 0.000
## all.parties.lvl2No answer                                1595.33    2.03 0.042
## all.parties.lvl2Other party                               214.12   11.77 0.000
## all.parties.lvl2Pro-environment party                     232.78   14.86 0.000
##                                                            LL    UL
## (Intercept)                                             -0.99 -0.51
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.21  0.28
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.28  0.25
## all.parties.lvl2Did not vote                             0.45  0.81
## all.parties.lvl2Don't know                               0.39  0.77
## all.parties.lvl2Invalid vote                            -0.20  1.74
## all.parties.lvl2NE age                                   0.83  1.22
## all.parties.lvl2NE citizen                               0.95  1.35
## all.parties.lvl2NE other                                 0.76  1.29
## all.parties.lvl2No answer                                0.03  1.97
## all.parties.lvl2Other party                              0.65  0.92
## all.parties.lvl2Pro-environment party                    1.12  1.46
```

```r
(VC.H2.mod2<-getVC(H2.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2549113 0.06497977
## 2     Residual        <NA> <NA> 1.1820631 1.39727308
```

```r
anova(H2.mod1,H2.mod2)
```

```
## Data: dat
## Models:
## H2.mod1: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H2.mod1:     occup
## H2.mod2: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H2.mod2:     occup + all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H2.mod1   19 85874 86030 -42918    85836                         
## H2.mod2   28 85708 85937 -42826    85652 184.21  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H2.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                  Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age                7.29   7.293     1 26886   5.2196   0.02234 *  
## gender            58.66  58.664     1 26797  41.9844 9.360e-11 ***
## educ             294.40 294.403     1 26553 210.6985 < 2.2e-16 ***
## resid             82.85  82.850     1 26881  59.2944 1.404e-14 ***
## occup            236.63  19.719    12 26808  14.1124 < 2.2e-16 ***
## all.parties.lvl2 389.10  43.233     9   314  30.9409 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H2.total.eff<-(VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"]-
     VC.H2.mod2[VC.H2.mod2$grp=="voting.group","est_SD2"])/
  VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.6197271
```

\newpage

#### Marginal means for and contrasts between Pro-environment and Anti-immigration parties


```r
H2.mod2.mmeans<-emmeans(H2.mod2,specs="all.parties.lvl2")

H2.mod2.mmeans.tab<-cbind(group=data.frame(H2.mod2.mmeans)[,1],
      data.frame(H2.mod2.mmeans)[,2:6])
H2.mod2.mmeans.tab$p<-
  2*(1-pnorm(abs(H2.mod2.mmeans.tab$emmean/
                   H2.mod2.mmeans.tab$SE)))
H2.mod2.mmeans.tab$adj.p<-
  p.adjust(H2.mod2.mmeans.tab$p,method="holm")

H2.mod2.mmeans.tab<-
  cbind(group=H2.mod2.mmeans.tab[,1],
      round(H2.mod2.mmeans.tab[,c(2,3)],2),
      round(H2.mod2.mmeans.tab[,c(7,8)],4),
      round(H2.mod2.mmeans.tab[,c(5,6)],2))
H2.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.74 0.06 0.0000 0.0000     -0.86     -0.62
## 2            Did not vote  -0.11 0.07 0.0996 0.3985     -0.24      0.02
## 3              Don't know  -0.16 0.08 0.0405 0.2023     -0.31     -0.01
## 4            Invalid vote   0.03 0.49 0.9500 1.0000     -0.93      0.99
## 5                  NE age   0.28 0.08 0.0003 0.0021      0.13      0.44
## 6              NE citizen   0.41 0.08 0.0000 0.0000      0.25      0.57
## 7                NE other   0.28 0.12 0.0187 0.1125      0.05      0.52
## 8               No answer   0.26 0.49 0.5920 1.0000     -0.70      1.22
## 9             Other party   0.05 0.03 0.1231 0.3985     -0.01      0.11
## 10  Pro-environment party   0.55 0.06 0.0000 0.0000      0.42      0.67
```

```r
write.csv2(H2.mod2.mmeans.tab,"H2.mod2.mmeans.tab.csv")


#contrast between anti-immigration and pro-environment
(H2.contrast<-data.frame(pairs(H2.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 1.286806 0.08657952 Inf
##    z.ratio      p.value
## 1 14.86271 5.754788e-50
```

```r
#contrast for all groups against mean of other groups
contrast(H2.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.9165 0.1007 Inf -9.106  <.0001 
##  Did not vote effect            -0.2180 0.1042 Inf -2.092  0.2185 
##  Don't know effect              -0.2713 0.1110 Inf -2.445  0.1015 
##  Invalid vote effect            -0.0603 0.4932 Inf -0.122  1.0000 
##  NE age effect                   0.2191 0.1115 Inf  1.965  0.2469 
##  NE citizen effect               0.3595 0.1132 Inf  3.176  0.0120 
##  NE other effect                 0.2196 0.1437 Inf  1.528  0.5056 
##  No answer effect                0.1972 0.4933 Inf  0.400  1.0000 
##  Other party effect             -0.0425 0.0853 Inf -0.499  1.0000 
##  Pro-environment party effect    0.5133 0.1012 Inf  5.072  <.0001 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H2.more.contrasts<-data.frame(pairs(H2.mod2.mmeans, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1           Other party - Anti-immigration party 0.7865943 0.06683493 Inf
## 2 Pro-environment party - Anti-immigration party 1.2868064 0.08657952 Inf
## 3            Pro-environment party - Other party 0.5002121 0.06764793 Inf
##     z.ratio      p.value
## 1 11.769210 1.124977e-31
## 2 14.862712 1.726437e-49
## 3  7.394345 1.421070e-13
```

\newpage

#### Effect size for the difference between Anti-immigration and Pro-environment party voters

Pool the standard deviations first within both groups and then across


```r
H2.anti.imm.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Anti-immigration party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(refugees),
            pro.ref.sd=sd(refugees),
            n=n())
H2.anti.imm.sd.dat$numerator<-
  (H2.anti.imm.sd.dat$n-1)*H2.anti.imm.sd.dat$pro.ref.sd^2

H2.anti.imm.sd<-sqrt(sum(H2.anti.imm.sd.dat$numerator)/
  ((sum(H2.anti.imm.sd.dat$n)-nrow(H2.anti.imm.sd.dat))))

H2.anti.imm.sd
```

```
## [1] 1.165751
```

```r
H2.pro.env.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Pro-environment party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(refugees),
            pro.ref.sd=sd(refugees),
            n=n())
H2.pro.env.sd.dat$numerator<-
  (H2.pro.env.sd.dat$n-1)*H2.pro.env.sd.dat$pro.ref.sd^2

H2.pro.env.sd<-sqrt(sum(H2.pro.env.sd.dat$numerator)/
  ((sum(H2.pro.env.sd.dat$n)-nrow(H2.pro.env.sd.dat))))
H2.pro.env.sd
```

```
## [1] 1.152943
```

```r
H2.pooled.sd<-sqrt(
  ((nrow(H2.anti.imm.sd.dat)-1)*H2.anti.imm.sd^2+
  (nrow(H2.pro.env.sd.dat)-1)*H2.pro.env.sd^2)/
  (nrow(H2.anti.imm.sd.dat)+
     nrow(H2.pro.env.sd.dat)-2))
H2.pooled.sd  
```

```
## [1] 1.159128
```

```r
(H2.effect.size<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd)
```

```
## [1] 1.112906
```

```r
H2.other.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Other party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(refugees),
            pro.ref.sd=sd(refugees),
            n=n())
H2.other.sd.dat$numerator<-
  (H2.other.sd.dat$n-1)*H2.other.sd.dat$pro.ref.sd^2

H2.other.sd<-sqrt(sum(H2.other.sd.dat$numerator)/
  ((sum(H2.other.sd.dat$n)-nrow(H2.other.sd.dat))))
H2.other.sd
```

```
## [1] 1.224836
```

```r
H2.pooled.sd.other.env<-sqrt(
  ((nrow(H2.other.sd.dat)-1)*H2.other.sd^2+
  (nrow(H2.pro.env.sd.dat)-1)*H2.pro.env.sd^2)/
  (nrow(H2.other.sd.dat)+
     nrow(H2.pro.env.sd.dat)-2))
H2.pooled.sd.other.env 
```

```
## [1] 1.190671
```

```r
(H2.effect.size.env.other<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[9,2])/H2.pooled.sd.other.env)
```

```
## [1] 0.4199312
```

```r
H2.pooled.sd.other.imm<-sqrt(
  ((nrow(H2.other.sd.dat)-1)*H2.other.sd^2+
  (nrow(H2.anti.imm.sd.dat)-1)*H2.anti.imm.sd^2)/
  (nrow(H2.other.sd.dat)+
     nrow(H2.anti.imm.sd.dat)-2))
H2.pooled.sd.other.imm 
```

```
## [1] 1.197766
```

```r
(H2.effect.size.imm.other<-(H2.mod2.mmeans.tab[9,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd.other.imm)
```

```
## [1] 0.659561
```



\newpage


### Model 3: Dummy-predictors at level-2 


```r
#did not vote left as reference

H2.mod3<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                other.party.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy+
                anti.imm.party.dummy+
                 pro.env.party.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


(FE.H2.mod3<-getFE(H2.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.13       0.13
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.01       0.11
## occupUnemployed                                            -0.02       0.14
## other.party.dummy                                           0.16       0.07
## dont.know.dummy                                            -0.05       0.10
## invalid.vote.dummy                                          0.14       0.49
## no.answer.dummy                                             0.37       0.49
## not.eligible.age.dummy                                      0.39       0.10
## not.eligible.citizenship.dummy                              0.52       0.10
## not.eligible.other.dummy                                    0.39       0.14
## anti.imm.party.dummy                                       -0.63       0.09
## pro.env.party.dummy                                         0.66       0.09
##                                                               df t.value     p
## (Intercept)                                              1841.24   -1.00 0.319
## age                                                     26885.99   -2.28 0.022
## gender                                                  26796.75    6.48 0.000
## educ                                                    26553.32   14.52 0.000
## resid                                                   26880.87   -7.70 0.000
## occupClerical support workers                           26795.71   -0.15 0.882
## occupCraft and related trades workers                   26800.78   -1.24 0.214
## occupElementary occupations                             26805.92   -0.02 0.984
## occupManagers                                           26795.33    0.62 0.538
## occupOther: Not in paid work                            26839.52    1.50 0.132
## occupPlant and machine operators, and assemblers        26805.29   -0.62 0.532
## occupProfessionals                                      26788.45    1.77 0.077
## occupRetired                                            26795.24    0.25 0.804
## occupService and sales workers                          26794.02   -0.60 0.548
## occupSkilled agricultural, forestry and fishery workers 26796.88   -0.18 0.860
## occupTechnicians and associate professionals            26790.15    0.06 0.951
## occupUnemployed                                         26807.62   -0.14 0.886
## other.party.dummy                                         160.59    2.18 0.031
## dont.know.dummy                                           199.17   -0.47 0.638
## invalid.vote.dummy                                       1549.87    0.29 0.774
## no.answer.dummy                                          1551.48    0.76 0.450
## not.eligible.age.dummy                                    200.14    3.86 0.000
## not.eligible.citizenship.dummy                            202.05    4.99 0.000
## not.eligible.other.dummy                                  484.85    2.87 0.004
## anti.imm.party.dummy                                      174.91   -6.97 0.000
## pro.env.party.dummy                                       189.61    7.24 0.000
##                                                            LL    UL
## (Intercept)                                             -0.37  0.12
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.21  0.28
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.28  0.25
## other.party.dummy                                        0.01  0.30
## dont.know.dummy                                         -0.25  0.15
## invalid.vote.dummy                                      -0.83  1.11
## no.answer.dummy                                         -0.60  1.34
## not.eligible.age.dummy                                   0.19  0.59
## not.eligible.citizenship.dummy                           0.31  0.73
## not.eligible.other.dummy                                 0.12  0.66
## anti.imm.party.dummy                                    -0.81 -0.45
## pro.env.party.dummy                                      0.48  0.84
```

```r
(VC.H2.mod3<-getVC(H2.mod3))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2549113 0.06497977
## 2     Residual        <NA> <NA> 1.1820631 1.39727308
```

```r
#this just confirms that the dummy and categorical
#models are identical
anova(H2.mod2,H2.mod3)
```

```
## Data: dat
## Models:
## H2.mod2: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H2.mod2:     occup + all.parties.lvl2
## H2.mod3: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H2.mod3:     occup + other.party.dummy + dont.know.dummy + invalid.vote.dummy + 
## H2.mod3:     no.answer.dummy + not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H2.mod3:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H2.mod2   28 85708 85937 -42826    85652                        
## H2.mod3   28 85708 85937 -42826    85652     0  0  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


\newpage


### Model 4: Dummy-predictors (anti-immigration and pro-environment) at level-2 allowed to vary between countries


```r
#did not vote left as reference

H2.mod4<-lmer(refugees~(1|voting.group)+
                (0+anti.imm.party.dummy+pro.env.party.dummy||cntry)+
                age+gender+educ+resid+occup+
                other.party.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy+
                anti.imm.party.dummy+
                 pro.env.party.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

isSingular(H2.mod4)
```

```
## [1] FALSE
```

```r
(FE.H2.mod4<-getFE(H2.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.13       0.12
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.01       0.11
## occupUnemployed                                            -0.02       0.14
## other.party.dummy                                           0.16       0.07
## dont.know.dummy                                            -0.05       0.10
## invalid.vote.dummy                                          0.14       0.49
## no.answer.dummy                                             0.37       0.49
## not.eligible.age.dummy                                      0.39       0.10
## not.eligible.citizenship.dummy                              0.52       0.10
## not.eligible.other.dummy                                    0.39       0.13
## anti.imm.party.dummy                                       -0.65       0.10
## pro.env.party.dummy                                         0.67       0.10
##                                                               df t.value     p
## (Intercept)                                              1726.13   -1.02 0.310
## age                                                     26881.61   -2.31 0.021
## gender                                                  26798.21    6.47 0.000
## educ                                                    26453.94   14.52 0.000
## resid                                                   26870.38   -7.71 0.000
## occupClerical support workers                           26797.24   -0.14 0.885
## occupCraft and related trades workers                   26801.73   -1.24 0.216
## occupElementary occupations                             26807.67   -0.01 0.988
## occupManagers                                           26796.45    0.61 0.541
## occupOther: Not in paid work                            26840.76    1.51 0.131
## occupPlant and machine operators, and assemblers        26806.91   -0.62 0.536
## occupProfessionals                                      26790.24    1.77 0.077
## occupRetired                                            26796.82    0.25 0.803
## occupService and sales workers                          26795.63   -0.60 0.550
## occupSkilled agricultural, forestry and fishery workers 26800.14   -0.18 0.858
## occupTechnicians and associate professionals            26790.88    0.06 0.950
## occupUnemployed                                         26809.25   -0.13 0.893
## other.party.dummy                                         126.05    2.30 0.023
## dont.know.dummy                                           160.19   -0.49 0.623
## invalid.vote.dummy                                       1445.33    0.29 0.770
## no.answer.dummy                                          1446.88    0.77 0.441
## not.eligible.age.dummy                                    161.04    4.07 0.000
## not.eligible.citizenship.dummy                            161.38    5.26 0.000
## not.eligible.other.dummy                                  409.62    2.93 0.004
## anti.imm.party.dummy                                       32.94   -6.50 0.000
## pro.env.party.dummy                                        30.54    6.46 0.000
##                                                            LL    UL
## (Intercept)                                             -0.37  0.12
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.22
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.21  0.28
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.28  0.25
## other.party.dummy                                        0.02  0.29
## dont.know.dummy                                         -0.24  0.14
## invalid.vote.dummy                                      -0.81  1.09
## no.answer.dummy                                         -0.58  1.33
## not.eligible.age.dummy                                   0.20  0.58
## not.eligible.citizenship.dummy                           0.32  0.71
## not.eligible.other.dummy                                 0.13  0.65
## anti.imm.party.dummy                                    -0.85 -0.44
## pro.env.party.dummy                                      0.46  0.88
```

```r
(VC.H2.mod4<-getVC(H2.mod4))
```

```
##            grp                 var1 var2    est_SD    est_SD2
## 1 voting.group          (Intercept) <NA> 0.2381506 0.05671572
## 2        cntry anti.imm.party.dummy <NA> 0.1827981 0.03341514
## 3      cntry.1  pro.env.party.dummy <NA> 0.2089535 0.04366157
## 4     Residual                 <NA> <NA> 1.1821198 1.39740727
```

```r
anova(H2.mod3,H2.mod4)
```

```
## Data: dat
## Models:
## H2.mod3: refugees ~ (1 | voting.group) + age + gender + educ + resid + 
## H2.mod3:     occup + other.party.dummy + dont.know.dummy + invalid.vote.dummy + 
## H2.mod3:     no.answer.dummy + not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H2.mod3:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
## H2.mod4: refugees ~ (1 | voting.group) + (0 + anti.imm.party.dummy + pro.env.party.dummy || 
## H2.mod4:     cntry) + age + gender + educ + resid + occup + other.party.dummy + 
## H2.mod4:     dont.know.dummy + invalid.vote.dummy + no.answer.dummy + 
## H2.mod4:     not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H2.mod4:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## H2.mod3   28 85708 85937 -42826    85652                     
## H2.mod4   30 85708 85954 -42824    85648 3.9917  2     0.1359
```

\newpage

### Model 5: explained variance by the focus groups


```r
#leave the focus group dummies out

H2.mod5<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                did.not.vote.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

(FE.H2.mod5<-getFE(H2.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.00       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.04       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                -0.01       0.11
## occupManagers                                               0.06       0.11
## occupOther: Not in paid work                                0.16       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.02       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.00       0.11
## occupUnemployed                                            -0.02       0.14
## did.not.vote.dummy                                         -0.12       0.11
## dont.know.dummy                                            -0.17       0.11
## invalid.vote.dummy                                          0.02       0.57
## no.answer.dummy                                             0.25       0.57
## not.eligible.age.dummy                                      0.28       0.11
## not.eligible.citizenship.dummy                              0.40       0.12
## not.eligible.other.dummy                                    0.30       0.15
##                                                               df t.value     p
## (Intercept)                                             15816.36   -0.02 0.987
## age                                                     26847.00   -2.33 0.020
## gender                                                  26747.27    6.55 0.000
## educ                                                    26881.23   14.69 0.000
## resid                                                   26846.56   -7.79 0.000
## occupClerical support workers                           26735.24   -0.17 0.865
## occupCraft and related trades workers                   26739.49   -1.29 0.198
## occupElementary occupations                             26741.40   -0.07 0.943
## occupManagers                                           26735.94    0.59 0.557
## occupOther: Not in paid work                            26769.28    1.44 0.149
## occupPlant and machine operators, and assemblers        26741.09   -0.66 0.508
## occupProfessionals                                      26732.62    1.75 0.081
## occupRetired                                            26737.54    0.20 0.844
## occupService and sales workers                          26734.81   -0.65 0.517
## occupSkilled agricultural, forestry and fishery workers 26738.23   -0.19 0.852
## occupTechnicians and associate professionals            26731.40    0.04 0.972
## occupUnemployed                                         26751.82   -0.18 0.854
## did.not.vote.dummy                                        192.89   -1.14 0.257
## dont.know.dummy                                           246.54   -1.51 0.131
## invalid.vote.dummy                                        787.82    0.04 0.969
## no.answer.dummy                                           788.04    0.44 0.659
## not.eligible.age.dummy                                    254.47    2.44 0.015
## not.eligible.citizenship.dummy                            262.33    3.43 0.001
## not.eligible.other.dummy                                  594.99    2.03 0.042
##                                                            LL    UL
## (Intercept)                                             -0.22  0.22
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.36  0.07
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.06  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.22  0.27
## occupService and sales workers                          -0.28  0.14
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.29  0.24
## did.not.vote.dummy                                      -0.33  0.09
## dont.know.dummy                                         -0.39  0.05
## invalid.vote.dummy                                      -1.10  1.15
## no.answer.dummy                                         -0.87  1.38
## not.eligible.age.dummy                                   0.05  0.50
## not.eligible.citizenship.dummy                           0.17  0.62
## not.eligible.other.dummy                                 0.01  0.60
```

```r
(VC.H2.mod5<-getVC(H2.mod5))
```

```
##            grp        var1 var2    est_SD  est_SD2
## 1 voting.group (Intercept) <NA> 0.3908401 0.152756
## 2     Residual        <NA> <NA> 1.1820293 1.397193
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H2.total.eff<-(VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"]-
     VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])/
  VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1060454
```

```r
#see how much residual variance was explained at level-2 by anti-immigrants

H2.mod6<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                did.not.vote.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy+
                anti.imm.party.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

(FE.H2.mod6<-getFE(H2.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.12       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                 0.00       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.01       0.11
## occupUnemployed                                            -0.02       0.14
## did.not.vote.dummy                                         -0.25       0.08
## dont.know.dummy                                            -0.29       0.09
## invalid.vote.dummy                                         -0.10       0.51
## no.answer.dummy                                             0.13       0.51
## not.eligible.age.dummy                                      0.15       0.09
## not.eligible.citizenship.dummy                              0.27       0.09
## not.eligible.other.dummy                                    0.16       0.13
## anti.imm.party.dummy                                       -0.87       0.07
##                                                               df t.value     p
## (Intercept)                                             19152.33    1.09 0.274
## age                                                     26883.16   -2.42 0.015
## gender                                                  26781.81    6.53 0.000
## educ                                                    26741.47   14.61 0.000
## resid                                                   26885.04   -7.77 0.000
## occupClerical support workers                           26775.09   -0.16 0.876
## occupCraft and related trades workers                   26779.88   -1.25 0.212
## occupElementary occupations                             26783.94   -0.04 0.967
## occupManagers                                           26775.09    0.60 0.551
## occupOther: Not in paid work                            26817.04    1.49 0.137
## occupPlant and machine operators, and assemblers        26783.42   -0.63 0.527
## occupProfessionals                                      26770.07    1.78 0.076
## occupRetired                                            26775.62    0.23 0.815
## occupService and sales workers                          26773.82   -0.62 0.538
## occupSkilled agricultural, forestry and fishery workers 26776.78   -0.18 0.855
## occupTechnicians and associate professionals            26770.15    0.06 0.952
## occupUnemployed                                         26788.39   -0.15 0.880
## did.not.vote.dummy                                        171.67   -3.03 0.003
## dont.know.dummy                                           248.90   -3.28 0.001
## invalid.vote.dummy                                       1283.23   -0.20 0.839
## no.answer.dummy                                          1283.62    0.25 0.802
## not.eligible.age.dummy                                    261.52    1.62 0.106
## not.eligible.citizenship.dummy                            258.03    2.93 0.004
## not.eligible.other.dummy                                  722.47    1.22 0.223
## anti.imm.party.dummy                                      225.47  -11.88 0.000
##                                                            LL    UL
## (Intercept)                                             -0.10  0.34
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.41
## occupRetired                                            -0.22  0.27
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.29  0.24
## did.not.vote.dummy                                      -0.41 -0.09
## dont.know.dummy                                         -0.47 -0.12
## invalid.vote.dummy                                      -1.11  0.90
## no.answer.dummy                                         -0.87  1.13
## not.eligible.age.dummy                                  -0.03  0.33
## not.eligible.citizenship.dummy                           0.09  0.46
## not.eligible.other.dummy                                -0.10  0.41
## anti.imm.party.dummy                                    -1.02 -0.73
```

```r
(VC.H2.mod6<-getVC(H2.mod6))
```

```
##            grp        var1 var2   est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.293254 0.08599792
## 2     Residual        <NA> <NA> 1.182050 1.39724248
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod6[VC.H2.mod6$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.4370242
```

```r
#see how much residual variance was explained at level-2 by pro-environments

H2.mod7<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                did.not.vote.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy+
                pro.env.party.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

(FE.H2.mod7<-getFE(H2.mod7))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.10       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.14       0.11
## occupElementary occupations                                -0.01       0.11
## occupManagers                                               0.07       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.00       0.11
## occupUnemployed                                            -0.02       0.14
## did.not.vote.dummy                                         -0.02       0.09
## dont.know.dummy                                            -0.07       0.10
## invalid.vote.dummy                                          0.12       0.54
## no.answer.dummy                                             0.35       0.54
## not.eligible.age.dummy                                      0.37       0.10
## not.eligible.citizenship.dummy                              0.49       0.11
## not.eligible.other.dummy                                    0.39       0.14
## pro.env.party.dummy                                         0.63       0.08
##                                                               df t.value     p
## (Intercept)                                             16171.21   -0.87 0.382
## age                                                     26860.08   -2.21 0.027
## gender                                                  26758.93    6.51 0.000
## educ                                                    26851.37   14.62 0.000
## resid                                                   26867.13   -7.73 0.000
## occupClerical support workers                           26751.64   -0.17 0.868
## occupCraft and related trades workers                   26756.15   -1.29 0.199
## occupElementary occupations                             26759.00   -0.06 0.954
## occupManagers                                           26752.03    0.60 0.547
## occupOther: Not in paid work                            26789.76    1.46 0.145
## occupPlant and machine operators, and assemblers        26758.62   -0.66 0.510
## occupProfessionals                                      26747.26    1.74 0.083
## occupRetired                                            26753.26    0.21 0.837
## occupService and sales workers                          26750.86   -0.64 0.522
## occupSkilled agricultural, forestry and fishery workers 26754.11   -0.18 0.855
## occupTechnicians and associate professionals            26747.14    0.03 0.973
## occupUnemployed                                         26766.89   -0.18 0.855
## did.not.vote.dummy                                        183.07   -0.26 0.798
## dont.know.dummy                                           245.76   -0.73 0.468
## invalid.vote.dummy                                        963.21    0.22 0.828
## no.answer.dummy                                           963.50    0.64 0.521
## not.eligible.age.dummy                                    255.60    3.62 0.000
## not.eligible.citizenship.dummy                            259.35    4.69 0.000
## not.eligible.other.dummy                                  642.96    2.79 0.005
## pro.env.party.dummy                                       252.77    7.55 0.000
##                                                            LL    UL
## (Intercept)                                             -0.32  0.12
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.36  0.07
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.06  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.22  0.27
## occupService and sales workers                          -0.28  0.14
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.29  0.24
## did.not.vote.dummy                                      -0.21  0.16
## dont.know.dummy                                         -0.27  0.13
## invalid.vote.dummy                                      -0.95  1.18
## no.answer.dummy                                         -0.72  1.41
## not.eligible.age.dummy                                   0.17  0.57
## not.eligible.citizenship.dummy                           0.29  0.70
## not.eligible.other.dummy                                 0.12  0.67
## pro.env.party.dummy                                      0.47  0.80
```

```r
(VC.H2.mod7<-getVC(H2.mod7))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.3440338 0.1183592
## 2     Residual        <NA> <NA> 1.1820235 1.3971794
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod7[VC.H2.mod7$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.2251745
```

```r
#see how much residual variance was explained at level-2 by both focus groups

H2.mod8<-lmer(refugees~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                did.not.vote.dummy+
                dont.know.dummy+
                 invalid.vote.dummy+
                no.answer.dummy+
                not.eligible.age.dummy+
                 not.eligible.citizenship.dummy+
                not.eligible.other.dummy+
                anti.imm.party.dummy+
                pro.env.party.dummy
                ,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

(FE.H2.mod8<-getFE(H2.mod8))
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
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.19       0.11
## occupRetired                                                0.03       0.12
## occupService and sales workers                             -0.07       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals                0.01       0.11
## occupUnemployed                                            -0.02       0.14
## did.not.vote.dummy                                         -0.16       0.07
## dont.know.dummy                                            -0.21       0.08
## invalid.vote.dummy                                         -0.02       0.49
## no.answer.dummy                                             0.22       0.49
## not.eligible.age.dummy                                      0.24       0.08
## not.eligible.citizenship.dummy                              0.36       0.08
## not.eligible.other.dummy                                    0.24       0.12
## anti.imm.party.dummy                                       -0.79       0.07
## pro.env.party.dummy                                         0.50       0.07
##                                                               df t.value     p
## (Intercept)                                             19007.55    0.29 0.768
## age                                                     26885.99   -2.28 0.022
## gender                                                  26796.75    6.48 0.000
## educ                                                    26553.32   14.52 0.000
## resid                                                   26880.87   -7.70 0.000
## occupClerical support workers                           26795.71   -0.15 0.882
## occupCraft and related trades workers                   26800.78   -1.24 0.214
## occupElementary occupations                             26805.92   -0.02 0.984
## occupManagers                                           26795.33    0.62 0.538
## occupOther: Not in paid work                            26839.52    1.50 0.132
## occupPlant and machine operators, and assemblers        26805.29   -0.62 0.532
## occupProfessionals                                      26788.45    1.77 0.077
## occupRetired                                            26795.24    0.25 0.804
## occupService and sales workers                          26794.02   -0.60 0.548
## occupSkilled agricultural, forestry and fishery workers 26796.88   -0.18 0.860
## occupTechnicians and associate professionals            26790.15    0.06 0.951
## occupUnemployed                                         26807.62   -0.14 0.886
## did.not.vote.dummy                                        160.59   -2.18 0.031
## dont.know.dummy                                           249.60   -2.52 0.012
## invalid.vote.dummy                                       1663.54   -0.03 0.974
## no.answer.dummy                                          1664.01    0.44 0.660
## not.eligible.age.dummy                                    265.18    2.83 0.005
## not.eligible.citizenship.dummy                            253.30    4.26 0.000
## not.eligible.other.dummy                                  777.60    1.91 0.056
## anti.imm.party.dummy                                      214.12  -11.77 0.000
## pro.env.party.dummy                                       248.54    7.39 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.03  0.04
## resid                                                   -0.15 -0.09
## occupClerical support workers                           -0.23  0.20
## occupCraft and related trades workers                   -0.35  0.08
## occupElementary occupations                             -0.22  0.21
## occupManagers                                           -0.15  0.28
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.15
## occupProfessionals                                      -0.02  0.40
## occupRetired                                            -0.21  0.28
## occupService and sales workers                          -0.28  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.21
## occupTechnicians and associate professionals            -0.21  0.22
## occupUnemployed                                         -0.28  0.25
## did.not.vote.dummy                                      -0.30 -0.01
## dont.know.dummy                                         -0.37 -0.04
## invalid.vote.dummy                                      -0.98  0.95
## no.answer.dummy                                         -0.75  1.18
## not.eligible.age.dummy                                   0.07  0.40
## not.eligible.citizenship.dummy                           0.19  0.53
## not.eligible.other.dummy                                -0.01  0.48
## anti.imm.party.dummy                                    -0.92 -0.65
## pro.env.party.dummy                                      0.37  0.63
```

```r
(VC.H2.mod8<-getVC(H2.mod8))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2549113 0.06497977
## 2     Residual        <NA> <NA> 1.1820631 1.39727308
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod8[VC.H2.mod8$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5746172
```

```r
#how much variance was left at level-2
1-H2.total.eff
```

```
## [1] 0.4253828
```

