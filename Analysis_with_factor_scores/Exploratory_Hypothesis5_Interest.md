---
title: "Exlopratory Hypothesis 5 with Political Interest"
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
dat.H5.intr<-dat %>%
  filter(!is.na(polintr))
```

### Model 1: without interactions (only main effects)


```r
H5.exp.intr.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1
                ,data=dat.H5.intr,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


isSingular(H5.exp.intr.mod1)
```

```
## [1] FALSE
```

```r
(FE.H5.exp.intr.mod1<-getFE(H5.exp.intr.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.11       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.03       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.03       0.13
## environ.lvl1                                                0.26       0.03
## polintr.lvl1                                                0.09       0.01
##                                                               df t.value     p
## (Intercept)                                             18697.43    0.34 0.731
## age                                                     25950.69   -2.84 0.005
## gender                                                  26723.18    6.83 0.000
## educ                                                    26777.13   10.79 0.000
## resid                                                   26818.15   -6.87 0.000
## occupClerical support workers                           26682.70   -0.24 0.811
## occupCraft and related trades workers                   26689.94   -1.12 0.262
## occupElementary occupations                             26692.92    0.22 0.823
## occupManagers                                           26686.32    0.32 0.751
## occupOther: Not in paid work                            26787.83    1.64 0.101
## occupPlant and machine operators, and assemblers        26691.87   -0.65 0.518
## occupProfessionals                                      26681.76    1.49 0.135
## occupRetired                                            26688.02    0.43 0.669
## occupService and sales workers                          26683.56   -0.51 0.611
## occupSkilled agricultural, forestry and fishery workers 26692.12   -0.14 0.891
## occupTechnicians and associate professionals            26679.94   -0.13 0.895
## occupUnemployed                                         26707.38    0.24 0.809
## environ.lvl1                                               15.74    7.99 0.000
## polintr.lvl1                                            26705.97    9.47 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.08  0.14
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.24
## occupManagers                                           -0.18  0.25
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.23  0.29
## environ.lvl1                                             0.19  0.33
## polintr.lvl1                                             0.07  0.10
```

```r
(VC.H5.exp.intr.mod1<-getVC(H5.exp.intr.mod1))
```

```
##            grp         var1         var2     est_SD    est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41911659 0.17565872
## 2 voting.group environ.lvl1         <NA> 0.08684613 0.00754225
## 3 voting.group  (Intercept) environ.lvl1 0.60539205 0.02203545
## 4        cntry environ.lvl1         <NA> 0.11887070 0.01413024
## 5     Residual         <NA>         <NA> 1.16044527 1.34663322
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political polintr


```r
H5.exp.intr.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1
                ,data=dat.H5.intr,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.intr.mod2)
```

```
## [1] FALSE
```

```r
anova(H5.exp.intr.mod1,H5.exp.intr.mod2)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.intr.mod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod1:     polintr.lvl1
## H5.exp.intr.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.intr.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod2:     polintr.lvl1 + environ.lvl1:polintr.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.intr.mod1   24 84938 85135 -42445    84890                         
## H5.exp.intr.mod2   25 84917 85122 -42434    84867 23.326  1  1.367e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod2<-getFE(H5.exp.intr.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.03       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.03       0.13
## environ.lvl1                                                0.26       0.03
## polintr.lvl1                                                0.09       0.01
## environ.lvl1:polintr.lvl1                                   0.05       0.01
##                                                               df t.value     p
## (Intercept)                                             18764.55    0.32 0.751
## age                                                     25915.37   -2.88 0.004
## gender                                                  26723.31    6.81 0.000
## educ                                                    26773.15   10.84 0.000
## resid                                                   26818.36   -6.79 0.000
## occupClerical support workers                           26682.96   -0.24 0.808
## occupCraft and related trades workers                   26690.20   -1.13 0.257
## occupElementary occupations                             26693.21    0.22 0.829
## occupManagers                                           26686.55    0.31 0.757
## occupOther: Not in paid work                            26788.30    1.62 0.105
## occupPlant and machine operators, and assemblers        26692.15   -0.65 0.513
## occupProfessionals                                      26681.95    1.48 0.138
## occupRetired                                            26688.31    0.41 0.679
## occupService and sales workers                          26683.81   -0.51 0.609
## occupSkilled agricultural, forestry and fishery workers 26692.29   -0.14 0.886
## occupTechnicians and associate professionals            26680.10   -0.13 0.898
## occupUnemployed                                         26707.46    0.22 0.822
## environ.lvl1                                               15.78    7.98 0.000
## polintr.lvl1                                            26708.56    9.57 0.000
## environ.lvl1:polintr.lvl1                               26159.72    4.83 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.24
## occupManagers                                           -0.18  0.25
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.28  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.23  0.29
## environ.lvl1                                             0.19  0.33
## polintr.lvl1                                             0.07  0.10
## environ.lvl1:polintr.lvl1                                0.03  0.08
```

```r
(VC.H5.exp.intr.mod2<-getVC(H5.exp.intr.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.41742694 0.174245253
## 2 voting.group environ.lvl1         <NA> 0.08857118 0.007844854
## 3 voting.group  (Intercept) environ.lvl1 0.61506265 0.022740095
## 4        cntry environ.lvl1         <NA> 0.11822937 0.013978183
## 5     Residual         <NA>         <NA> 1.15996033 1.345507968
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political polintr, allow polintr effect to vary between voting groups and countries


```r
H5.exp.intr.mod3<-lmer(refugees~(environ.lvl1+polintr.lvl1|voting.group)+
                (0+environ.lvl1+polintr.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1
                ,data=dat.H5.intr,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Model failed to converge with 1 negative eigenvalue: -6.0e+03
```

```r
isSingular(H5.exp.intr.mod3)
```

```
## [1] TRUE
```

```r
anova(H5.exp.intr.mod2,H5.exp.intr.mod3)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H5.exp.intr.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod2:     polintr.lvl1 + environ.lvl1:polintr.lvl1
## H5.exp.intr.mod3: refugees ~ (environ.lvl1 + polintr.lvl1 | voting.group) + (0 + 
## H5.exp.intr.mod3:     environ.lvl1 + polintr.lvl1 | cntry) + age + gender + educ + 
## H5.exp.intr.mod3:     resid + occup + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.intr.mod2   25 84917 85122 -42434    84867                         
## H5.exp.intr.mod3   30 84878 85124 -42409    84818 48.839  5  2.394e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod3<-getFE(H5.exp.intr.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.03       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.03       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.03       0.13
## environ.lvl1                                                0.25       0.03
## polintr.lvl1                                                0.09       0.02
## environ.lvl1:polintr.lvl1                                   0.05       0.01
##                                                               df t.value     p
## (Intercept)                                             18715.41    0.32 0.751
## age                                                     25811.97   -2.97 0.003
## gender                                                  26717.39    6.79 0.000
## educ                                                    26726.48   10.79 0.000
## resid                                                   26789.73   -6.78 0.000
## occupClerical support workers                           26648.23   -0.24 0.811
## occupCraft and related trades workers                   26655.47   -1.13 0.260
## occupElementary occupations                             26658.00    0.23 0.815
## occupManagers                                           26656.35    0.34 0.735
## occupOther: Not in paid work                            26758.91    1.62 0.104
## occupPlant and machine operators, and assemblers        26664.81   -0.66 0.506
## occupProfessionals                                      26651.06    1.48 0.140
## occupRetired                                            26655.15    0.42 0.672
## occupService and sales workers                          26646.94   -0.48 0.629
## occupSkilled agricultural, forestry and fishery workers 26663.12   -0.11 0.910
## occupTechnicians and associate professionals            26646.84   -0.13 0.900
## occupUnemployed                                         26681.24    0.23 0.815
## environ.lvl1                                               15.58    7.98 0.000
## polintr.lvl1                                               17.33    5.57 0.000
## environ.lvl1:polintr.lvl1                               26636.59    4.85 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.24
## occupManagers                                           -0.18  0.25
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.23  0.29
## environ.lvl1                                             0.18  0.32
## polintr.lvl1                                             0.06  0.13
## environ.lvl1:polintr.lvl1                                0.03  0.08
```

```r
(VC.H5.exp.intr.mod3<-getVC(H5.exp.intr.mod3))
```

```
##             grp         var1         var2     est_SD      est_SD2
## 1  voting.group  (Intercept)         <NA> 0.41763282 0.1744171699
## 2  voting.group environ.lvl1         <NA> 0.05248939 0.0027551363
## 3  voting.group polintr.lvl1         <NA> 0.08693883 0.0075583610
## 4  voting.group  (Intercept) environ.lvl1 1.00000000 0.0219212927
## 5  voting.group  (Intercept) polintr.lvl1 0.43921958 0.0159474088
## 6  voting.group environ.lvl1 polintr.lvl1 0.43921958 0.0020043200
## 7         cntry environ.lvl1         <NA> 0.11670510 0.0136200794
## 8         cntry polintr.lvl1         <NA> 0.04883098 0.0023844645
## 9         cntry environ.lvl1 polintr.lvl1 0.10832747 0.0006173392
## 10     Residual         <NA>         <NA> 1.15804094 1.3410588282
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.exp.intr.mod4<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+
                   environ.lvl1:polintr.lvl1|voting.group)+
                (0+environ.lvl1+polintr.lvl1+
                   environ.lvl1:polintr.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1
                ,data=dat.H5.intr,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```r
isSingular(H5.exp.intr.mod4)
```

```
## [1] TRUE
```

```r
anova(H5.exp.intr.mod3,H5.exp.intr.mod4)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod3: refugees ~ (environ.lvl1 + polintr.lvl1 | voting.group) + (0 + 
## H5.exp.intr.mod3:     environ.lvl1 + polintr.lvl1 | cntry) + age + gender + educ + 
## H5.exp.intr.mod3:     resid + occup + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1
## H5.exp.intr.mod4: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod4:     polintr.lvl1 + environ.lvl1:polintr.lvl1
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.exp.intr.mod3   30 84878 85124 -42409    84818                        
## H5.exp.intr.mod4   37 84870 85174 -42398    84796 21.943  7   0.002599 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod4<-getFE(H5.exp.intr.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.11
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.18       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.15       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.02       0.11
## occupUnemployed                                             0.04       0.13
## environ.lvl1                                                0.25       0.03
## polintr.lvl1                                                0.09       0.02
## environ.lvl1:polintr.lvl1                                   0.05       0.02
##                                                               df t.value     p
## (Intercept)                                             18662.41    0.35 0.728
## age                                                     25830.07   -3.05 0.002
## gender                                                  26687.81    6.80 0.000
## educ                                                    26636.61   10.85 0.000
## resid                                                   26764.67   -6.88 0.000
## occupClerical support workers                           26629.94   -0.27 0.786
## occupCraft and related trades workers                   26639.72   -1.16 0.248
## occupElementary occupations                             26643.54    0.20 0.845
## occupManagers                                           26638.86    0.33 0.743
## occupOther: Not in paid work                            26746.03    1.59 0.112
## occupPlant and machine operators, and assemblers        26650.75   -0.68 0.494
## occupProfessionals                                      26629.51    1.44 0.149
## occupRetired                                            26635.34    0.39 0.694
## occupService and sales workers                          26632.51   -0.52 0.605
## occupSkilled agricultural, forestry and fishery workers 26647.22   -0.15 0.880
## occupTechnicians and associate professionals            26629.67   -0.14 0.885
## occupUnemployed                                         26664.27    0.27 0.784
## environ.lvl1                                               15.77    8.21 0.000
## polintr.lvl1                                               17.16    5.37 0.000
## environ.lvl1:polintr.lvl1                                  18.89    2.97 0.008
##                                                            LL    UL
## (Intercept)                                             -0.18  0.25
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.24  0.18
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.18  0.25
## occupOther: Not in paid work                            -0.04  0.40
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.06  0.36
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.15
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.19
## occupUnemployed                                         -0.22  0.30
## environ.lvl1                                             0.19  0.32
## polintr.lvl1                                             0.06  0.13
## environ.lvl1:polintr.lvl1                                0.01  0.08
```

```r
(VC.H5.exp.intr.mod4<-getVC(H5.exp.intr.mod4))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.41817800
## 2  voting.group              environ.lvl1                      <NA>  0.08621675
## 3  voting.group              polintr.lvl1                      <NA>  0.08813827
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.05692954
## 5  voting.group               (Intercept)              environ.lvl1  0.59062673
## 6  voting.group               (Intercept)              polintr.lvl1  0.42642296
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1 -0.16488856
## 8  voting.group              environ.lvl1              polintr.lvl1  0.09152983
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.86334554
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1 -0.11820118
## 11        cntry              environ.lvl1                      <NA>  0.11163389
## 12        cntry              polintr.lvl1                      <NA>  0.05022285
## 13        cntry environ.lvl1:polintr.lvl1                      <NA>  0.03973149
## 14        cntry              environ.lvl1              polintr.lvl1  0.14846748
## 15        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.94302988
## 16        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.46903000
## 17     Residual                      <NA>                      <NA>  1.15634476
##          est_SD2
## 1   0.1748728373
## 2   0.0074333277
## 3   0.0077683552
## 4   0.0032409728
## 5   0.0212944250
## 6   0.0157168786
## 7  -0.0039254495
## 8   0.0006955347
## 9  -0.0042375417
## 10 -0.0005930947
## 11  0.0124621246
## 12  0.0025223342
## 13  0.0015785910
## 14  0.0008323935
## 15  0.0041826960
## 16  0.0009359157
## 17  1.3371332136
```

```r
theta <- getME(H5.exp.intr.mod4,"theta")

## diagonal elements are identifiable because they are fitted
##  with a lower bound of zero ...
diag.element <- getME(H5.exp.intr.mod4,"lower")==0
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
##                                             0.36164 
##               voting.group.environ.lvl1.(Intercept) 
##                                             0.04404 
##               voting.group.polintr.lvl1.(Intercept) 
##                                             0.03250 
##  voting.group.environ.lvl1:polintr.lvl1.(Intercept) 
##                                            -0.00812 
##                           voting.group.environ.lvl1 
##                                             0.06017 
##              voting.group.polintr.lvl1.environ.lvl1 
##                                            -0.01514 
## voting.group.environ.lvl1:polintr.lvl1.environ.lvl1 
##                                            -0.04673 
##                           voting.group.polintr.lvl1 
##                                             0.06726 
## voting.group.environ.lvl1:polintr.lvl1.polintr.lvl1 
##                                            -0.01319 
##              voting.group.environ.lvl1:polintr.lvl1 
##                                             0.00000 
##                                  cntry.environ.lvl1 
##                                             0.09654 
##                     cntry.polintr.lvl1.environ.lvl1 
##                                             0.00645 
##        cntry.environ.lvl1:polintr.lvl1.environ.lvl1 
##                                             0.03240 
##                                  cntry.polintr.lvl1 
##                                             0.04295 
##        cntry.environ.lvl1:polintr.lvl1.polintr.lvl1 
##                                             0.01143 
##                     cntry.environ.lvl1:polintr.lvl1 
##                                             0.00000
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and polintr


```r
H5.exp.intr.mod5<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|voting.group)+
                (0+environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polintr.lvl1:all.parties.lvl2
                ,data=dat.H5.intr,REML=F,
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
isSingular(H5.exp.intr.mod5)
```

```
## [1] TRUE
```

```r
anova(H5.exp.intr.mod4,H5.exp.intr.mod5)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod4: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod4:     polintr.lvl1 + environ.lvl1:polintr.lvl1
## H5.exp.intr.mod5: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H5.exp.intr.mod4   37 84870 85174 -42398    84796                        
## H5.exp.intr.mod5   64 84697 85222 -42285    84569 227.1 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod5<-getFE(H5.exp.intr.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.75       0.12
## age                                                         0.00       0.00
## gender                                                      0.10       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.15       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.11
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.03       0.13
## environ.lvl1                                                0.15       0.05
## polintr.lvl1                                                0.01       0.03
## all.parties.lvl2Did not vote                                0.62       0.09
## all.parties.lvl2Don't know                                  0.58       0.10
## all.parties.lvl2Invalid vote                                0.51       0.63
## all.parties.lvl2NE age                                      1.00       0.10
## all.parties.lvl2NE citizen                                  1.19       0.10
## all.parties.lvl2NE other                                    1.00       0.13
## all.parties.lvl2No answer                                   1.12       0.50
## all.parties.lvl2Other party                                 0.80       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:polintr.lvl1                                   0.05       0.02
## environ.lvl1:all.parties.lvl2Did not vote                   0.05       0.04
## environ.lvl1:all.parties.lvl2Don't know                     0.13       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   0.32       0.84
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.07
## environ.lvl1:all.parties.lvl2NE citizen                     0.06       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.02       0.13
## environ.lvl1:all.parties.lvl2No answer                      0.40       0.59
## environ.lvl1:all.parties.lvl2Other party                    0.13       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.10       0.06
## polintr.lvl1:all.parties.lvl2Did not vote                   0.05       0.04
## polintr.lvl1:all.parties.lvl2Don't know                     0.09       0.06
## polintr.lvl1:all.parties.lvl2Invalid vote                   0.40       0.67
## polintr.lvl1:all.parties.lvl2NE age                         0.01       0.06
## polintr.lvl1:all.parties.lvl2NE citizen                    -0.04       0.05
## polintr.lvl1:all.parties.lvl2NE other                       0.12       0.11
## polintr.lvl1:all.parties.lvl2No answer                      0.45       0.98
## polintr.lvl1:all.parties.lvl2Other party                    0.11       0.03
## polintr.lvl1:all.parties.lvl2Pro-environment party          0.20       0.05
##                                                               df t.value     p
## (Intercept)                                              3126.34   -6.26 0.000
## age                                                     26713.82   -2.57 0.010
## gender                                                  26720.75    6.80 0.000
## educ                                                    25600.05   10.83 0.000
## resid                                                   26775.85   -6.62 0.000
## occupClerical support workers                           26683.89   -0.23 0.816
## occupCraft and related trades workers                   26693.72   -1.11 0.265
## occupElementary occupations                             26696.96    0.21 0.837
## occupManagers                                           26701.81    0.34 0.731
## occupOther: Not in paid work                            26733.82    1.54 0.123
## occupPlant and machine operators, and assemblers        26708.26   -0.67 0.505
## occupProfessionals                                      26680.68    1.45 0.146
## occupRetired                                            26681.40    0.44 0.661
## occupService and sales workers                          26684.38   -0.49 0.624
## occupSkilled agricultural, forestry and fishery workers 26706.51   -0.16 0.875
## occupTechnicians and associate professionals            26684.83   -0.10 0.917
## occupUnemployed                                         26696.48    0.25 0.805
## environ.lvl1                                               55.18    3.22 0.002
## polintr.lvl1                                              110.65    0.18 0.855
## all.parties.lvl2Did not vote                              164.74    6.98 0.000
## all.parties.lvl2Don't know                                237.02    5.97 0.000
## all.parties.lvl2Invalid vote                             2861.86    0.82 0.412
## all.parties.lvl2NE age                                    245.21   10.23 0.000
## all.parties.lvl2NE citizen                                235.37   11.96 0.000
## all.parties.lvl2NE other                                  581.27    7.53 0.000
## all.parties.lvl2No answer                                1606.84    2.24 0.025
## all.parties.lvl2Other party                               213.04   12.02 0.000
## all.parties.lvl2Pro-environment party                     232.22   15.28 0.000
## environ.lvl1:polintr.lvl1                                  18.33    3.17 0.005
## environ.lvl1:all.parties.lvl2Did not vote                  67.81    1.15 0.256
## environ.lvl1:all.parties.lvl2Don't know                   283.73    1.99 0.047
## environ.lvl1:all.parties.lvl2Invalid vote               20290.47    0.38 0.703
## environ.lvl1:all.parties.lvl2NE age                       271.67    3.18 0.002
## environ.lvl1:all.parties.lvl2NE citizen                   204.94    0.98 0.330
## environ.lvl1:all.parties.lvl2NE other                    1080.01    0.12 0.905
## environ.lvl1:all.parties.lvl2No answer                  15409.84    0.67 0.501
## environ.lvl1:all.parties.lvl2Other party                  107.31    3.43 0.001
## environ.lvl1:all.parties.lvl2Pro-environment party        224.99    1.67 0.096
## polintr.lvl1:all.parties.lvl2Did not vote                  74.45    1.34 0.184
## polintr.lvl1:all.parties.lvl2Don't know                   323.38    1.50 0.135
## polintr.lvl1:all.parties.lvl2Invalid vote               18044.46    0.59 0.555
## polintr.lvl1:all.parties.lvl2NE age                       229.24    0.14 0.885
## polintr.lvl1:all.parties.lvl2NE citizen                   173.85   -0.77 0.442
## polintr.lvl1:all.parties.lvl2NE other                     899.16    1.11 0.269
## polintr.lvl1:all.parties.lvl2No answer                  25300.60    0.46 0.645
## polintr.lvl1:all.parties.lvl2Other party                  112.93    3.12 0.002
## polintr.lvl1:all.parties.lvl2Pro-environment party        197.24    3.95 0.000
##                                                            LL    UL
## (Intercept)                                             -0.99 -0.52
## age                                                      0.00  0.00
## gender                                                   0.07  0.13
## educ                                                     0.02  0.03
## resid                                                   -0.13 -0.07
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.33  0.09
## occupElementary occupations                             -0.19  0.23
## occupManagers                                           -0.17  0.25
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.05  0.36
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers -0.24  0.21
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.23  0.29
## environ.lvl1                                             0.06  0.24
## polintr.lvl1                                            -0.06  0.07
## all.parties.lvl2Did not vote                             0.44  0.80
## all.parties.lvl2Don't know                               0.39  0.77
## all.parties.lvl2Invalid vote                            -0.71  1.74
## all.parties.lvl2NE age                                   0.81  1.20
## all.parties.lvl2NE citizen                               0.99  1.38
## all.parties.lvl2NE other                                 0.74  1.26
## all.parties.lvl2No answer                                0.14  2.10
## all.parties.lvl2Other party                              0.67  0.93
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:polintr.lvl1                                0.02  0.08
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.14
## environ.lvl1:all.parties.lvl2Don't know                  0.00  0.27
## environ.lvl1:all.parties.lvl2Invalid vote               -1.33  1.97
## environ.lvl1:all.parties.lvl2NE age                      0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.06  0.19
## environ.lvl1:all.parties.lvl2NE other                   -0.23  0.26
## environ.lvl1:all.parties.lvl2No answer                  -0.76  1.55
## environ.lvl1:all.parties.lvl2Other party                 0.06  0.21
## environ.lvl1:all.parties.lvl2Pro-environment party      -0.02  0.21
## polintr.lvl1:all.parties.lvl2Did not vote               -0.03  0.13
## polintr.lvl1:all.parties.lvl2Don't know                 -0.03  0.21
## polintr.lvl1:all.parties.lvl2Invalid vote               -0.92  1.71
## polintr.lvl1:all.parties.lvl2NE age                     -0.10  0.12
## polintr.lvl1:all.parties.lvl2NE citizen                 -0.15  0.07
## polintr.lvl1:all.parties.lvl2NE other                   -0.09  0.33
## polintr.lvl1:all.parties.lvl2No answer                  -1.47  2.37
## polintr.lvl1:all.parties.lvl2Other party                 0.04  0.18
## polintr.lvl1:all.parties.lvl2Pro-environment party       0.10  0.29
```

```r
(VC.H5.exp.intr.mod5<-getVC(H5.exp.intr.mod5))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.25729842
## 2  voting.group              environ.lvl1                      <NA>  0.07546274
## 3  voting.group              polintr.lvl1                      <NA>  0.06461187
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.05027388
## 5  voting.group               (Intercept)              environ.lvl1  0.51416752
## 6  voting.group               (Intercept)              polintr.lvl1  0.40320882
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1 -0.50194524
## 8  voting.group              environ.lvl1              polintr.lvl1 -0.24109258
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.95614761
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1  0.29144639
## 11        cntry              environ.lvl1                      <NA>  0.11369716
## 12        cntry              polintr.lvl1                      <NA>  0.05767173
## 13        cntry environ.lvl1:polintr.lvl1                      <NA>  0.03870668
## 14        cntry              environ.lvl1              polintr.lvl1  0.16058866
## 15        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.91613508
## 16        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.54278789
## 17     Residual                      <NA>                      <NA>  1.15626171
##          est_SD2
## 1   0.0662024744
## 2   0.0056946250
## 3   0.0041746933
## 4   0.0025274631
## 5   0.0099833045
## 6   0.0067031574
## 7  -0.0064928573
## 8  -0.0011755164
## 9  -0.0036274373
## 10  0.0009467022
## 11  0.0129270444
## 12  0.0033260286
## 13  0.0014982069
## 14  0.0010529978
## 15  0.0040317634
## 16  0.0012116552
## 17  1.3369411441
```

\newpage

#### Look among which voting group there is strongest association between polintr and refugee attitudes


```r
H5.exp.intr.mod5.trends<-emtrends(H5.exp.intr.mod5,specs = c("all.parties.lvl2"),var=c("polintr.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 26873' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 26873)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 26873' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 26873)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.exp.intr.mod5.trends.tab<-data.frame(H5.exp.intr.mod5.trends))
```

```
##          all.parties.lvl2 polintr.lvl1.trend         SE  df    asymp.LCL
## 1  Anti-immigration party        0.006569538 0.03456168 Inf -0.061170117
## 2            Did not vote        0.060870688 0.02970721 Inf  0.002645631
## 3              Don't know        0.097057832 0.05368755 Inf -0.008167823
## 4            Invalid vote        0.403101111 0.67190809 Inf -0.913814543
## 5                  NE age        0.014574802 0.04803716 Inf -0.079576299
## 6              NE citizen       -0.035713015 0.04732294 Inf -0.128464273
## 7                NE other        0.125085943 0.10354806 Inf -0.077864518
## 8               No answer        0.457441113 0.97924648 Inf -1.461846713
## 9             Other party        0.115032214 0.02095351 Inf  0.073964086
## 10  Pro-environment party        0.202674496 0.04129123 Inf  0.121745162
##     asymp.UCL
## 1  0.07430919
## 2  0.11909575
## 3  0.20228349
## 4  1.72001676
## 5  0.10872590
## 6  0.05703824
## 7  0.32803640
## 8  2.37672894
## 9  0.15610034
## 10 0.28360383
```

```r
H5.exp.intr.mod5.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.intr.mod5.trends.tab$polintr.lvl1.trend/
                   H5.exp.intr.mod5.trends.tab$SE)))
H5.exp.intr.mod5.trends.tab$adj.p<-
  p.adjust(H5.exp.intr.mod5.trends.tab$p,method="holm")

H5.exp.intr.mod5.trends.tab<-
  cbind(group=H5.exp.intr.mod5.trends.tab[,1],
      round(H5.exp.intr.mod5.trends.tab[,c(2,3)],2),
      round(H5.exp.intr.mod5.trends.tab[,c(7,8)],4),
      round(H5.exp.intr.mod5.trends.tab[,c(5,6)],2))
H5.exp.intr.mod5.trends.tab
```

```
##                     group polintr.lvl1.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party               0.01 0.03 0.8492 1.0000     -0.06
## 2            Did not vote               0.06 0.03 0.0405 0.3237      0.00
## 3              Don't know               0.10 0.05 0.0706 0.4944     -0.01
## 4            Invalid vote               0.40 0.67 0.5485 1.0000     -0.91
## 5                  NE age               0.01 0.05 0.7616 1.0000     -0.08
## 6              NE citizen              -0.04 0.05 0.4504 1.0000     -0.13
## 7                NE other               0.13 0.10 0.2270 1.0000     -0.08
## 8               No answer               0.46 0.98 0.6404 1.0000     -1.46
## 9             Other party               0.12 0.02 0.0000 0.0000      0.07
## 10  Pro-environment party               0.20 0.04 0.0000 0.0000      0.12
##    asymp.UCL
## 1       0.07
## 2       0.12
## 3       0.20
## 4       1.72
## 5       0.11
## 6       0.06
## 7       0.33
## 8       2.38
## 9       0.16
## 10      0.28
```

```r
write.csv2(H5.exp.intr.mod5.trends.tab,"H5.exp.intr.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.196105 0.04966265 Inf
##    z.ratio      p.value
## 1 3.948741 7.856322e-05
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.intr.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect  -0.1534 0.136 Inf -1.125  0.2607 
##  Did not vote effect            -0.0931 0.135 Inf -0.688  0.4915 
##  Don't know effect              -0.0529 0.143 Inf -0.371  0.7105 
##  Invalid vote effect             0.2871 0.681 Inf  0.422  0.6731 
##  NE age effect                  -0.1445 0.140 Inf -1.029  0.3034 
##  NE citizen effect              -0.2004 0.140 Inf -1.430  0.1528 
##  NE other effect                -0.0218 0.167 Inf -0.130  0.8966 
##  No answer effect                0.3475 0.982 Inf  0.354  0.7235 
##  Other party effect             -0.0329 0.134 Inf -0.246  0.8054 
##  Pro-environment party effect    0.0644 0.138 Inf  0.466  0.6414 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic
```

```r
#contrast for three voting groups
(H4.more.contrasts<-data.frame(pairs(H5.exp.intr.mod5.trends, 
         exclude=c(2:8), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1           Other party - Anti-immigration party 0.10846268 0.03470886 Inf
## 2 Pro-environment party - Anti-immigration party 0.19610496 0.04966265 Inf
## 3            Pro-environment party - Other party 0.08764228 0.04140901 Inf
##    z.ratio      p.value
## 1 3.124928 1.778486e-03
## 2 3.948741 7.856322e-05
## 3 2.116503 3.430208e-02
```

\newpage

### Model 6: Enter three-way interaction voting group x polintr x environment attitudes


```r
H5.exp.intr.mod6<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|voting.group)+
                (0+environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polintr.lvl1:all.parties.lvl2+
                environ.lvl1:polintr.lvl1:all.parties.lvl2
                ,data=dat.H5.intr,REML=F,
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
isSingular(H5.exp.intr.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.intr.mod5,H5.exp.intr.mod6)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod5: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
## H5.exp.intr.mod6: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod6:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod6:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod6:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2 + 
## H5.exp.intr.mod6:     environ.lvl1:polintr.lvl1:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.intr.mod5   64 84697 85222 -42285    84569                     
## H5.exp.intr.mod6   73 84704 85302 -42279    84558 11.614  9     0.2359
```

```r
(FE.H5.exp.intr.mod6<-getFE(H5.exp.intr.mod6))
```

```
##                                                                 Estimate
## (Intercept)                                                        -0.75
## age                                                                 0.00
## gender                                                              0.11
## educ                                                                0.03
## resid                                                              -0.10
## occupClerical support workers                                      -0.03
## occupCraft and related trades workers                              -0.12
## occupElementary occupations                                         0.02
## occupManagers                                                       0.04
## occupOther: Not in paid work                                        0.17
## occupPlant and machine operators, and assemblers                   -0.07
## occupProfessionals                                                  0.15
## occupRetired                                                        0.05
## occupService and sales workers                                     -0.05
## occupSkilled agricultural, forestry and fishery workers            -0.02
## occupTechnicians and associate professionals                       -0.01
## occupUnemployed                                                     0.03
## environ.lvl1                                                        0.15
## polintr.lvl1                                                        0.01
## all.parties.lvl2Did not vote                                        0.62
## all.parties.lvl2Don't know                                          0.58
## all.parties.lvl2Invalid vote                                       -0.24
## all.parties.lvl2NE age                                              0.99
## all.parties.lvl2NE citizen                                          1.15
## all.parties.lvl2NE other                                            1.01
## all.parties.lvl2No answer                                           1.04
## all.parties.lvl2Other party                                         0.80
## all.parties.lvl2Pro-environment party                               1.32
## environ.lvl1:polintr.lvl1                                           0.05
## environ.lvl1:all.parties.lvl2Did not vote                           0.05
## environ.lvl1:all.parties.lvl2Don't know                             0.14
## environ.lvl1:all.parties.lvl2Invalid vote                           1.81
## environ.lvl1:all.parties.lvl2NE age                                 0.20
## environ.lvl1:all.parties.lvl2NE citizen                             0.06
## environ.lvl1:all.parties.lvl2NE other                               0.01
## environ.lvl1:all.parties.lvl2No answer                              0.92
## environ.lvl1:all.parties.lvl2Other party                            0.14
## environ.lvl1:all.parties.lvl2Pro-environment party                  0.10
## polintr.lvl1:all.parties.lvl2Did not vote                           0.06
## polintr.lvl1:all.parties.lvl2Don't know                             0.09
## polintr.lvl1:all.parties.lvl2Invalid vote                           1.13
## polintr.lvl1:all.parties.lvl2NE age                                 0.01
## polintr.lvl1:all.parties.lvl2NE citizen                            -0.04
## polintr.lvl1:all.parties.lvl2NE other                               0.13
## polintr.lvl1:all.parties.lvl2No answer                              0.73
## polintr.lvl1:all.parties.lvl2Other party                            0.11
## polintr.lvl1:all.parties.lvl2Pro-environment party                  0.19
## environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote              0.01
## environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know               -0.01
## environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote             -1.70
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE age                    0.07
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen                0.14
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE other                 -0.10
## environ.lvl1:polintr.lvl1:all.parties.lvl2No answer                -2.22
## environ.lvl1:polintr.lvl1:all.parties.lvl2Other party              -0.01
## environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party    -0.06
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
## occupSkilled agricultural, forestry and fishery workers               0.11
## occupTechnicians and associate professionals                          0.11
## occupUnemployed                                                       0.13
## environ.lvl1                                                          0.05
## polintr.lvl1                                                          0.03
## all.parties.lvl2Did not vote                                          0.09
## all.parties.lvl2Don't know                                            0.10
## all.parties.lvl2Invalid vote                                          0.98
## all.parties.lvl2NE age                                                0.10
## all.parties.lvl2NE citizen                                            0.10
## all.parties.lvl2NE other                                              0.13
## all.parties.lvl2No answer                                             0.53
## all.parties.lvl2Other party                                           0.07
## all.parties.lvl2Pro-environment party                                 0.09
## environ.lvl1:polintr.lvl1                                             0.04
## environ.lvl1:all.parties.lvl2Did not vote                             0.05
## environ.lvl1:all.parties.lvl2Don't know                               0.07
## environ.lvl1:all.parties.lvl2Invalid vote                             1.71
## environ.lvl1:all.parties.lvl2NE age                                   0.07
## environ.lvl1:all.parties.lvl2NE citizen                               0.06
## environ.lvl1:all.parties.lvl2NE other                                 0.13
## environ.lvl1:all.parties.lvl2No answer                                1.31
## environ.lvl1:all.parties.lvl2Other party                              0.04
## environ.lvl1:all.parties.lvl2Pro-environment party                    0.06
## polintr.lvl1:all.parties.lvl2Did not vote                             0.04
## polintr.lvl1:all.parties.lvl2Don't know                               0.06
## polintr.lvl1:all.parties.lvl2Invalid vote                             0.99
## polintr.lvl1:all.parties.lvl2NE age                                   0.06
## polintr.lvl1:all.parties.lvl2NE citizen                               0.05
## polintr.lvl1:all.parties.lvl2NE other                                 0.11
## polintr.lvl1:all.parties.lvl2No answer                                1.17
## polintr.lvl1:all.parties.lvl2Other party                              0.03
## polintr.lvl1:all.parties.lvl2Pro-environment party                    0.05
## environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote                0.04
## environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know                  0.07
## environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote                1.70
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE age                      0.07
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen                  0.06
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE other                    0.14
## environ.lvl1:polintr.lvl1:all.parties.lvl2No answer                   4.95
## environ.lvl1:polintr.lvl1:all.parties.lvl2Other party                 0.04
## environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party       0.06
##                                                                       df
## (Intercept)                                                      3022.65
## age                                                             26726.70
## gender                                                          26733.28
## educ                                                            25692.65
## resid                                                           26780.29
## occupClerical support workers                                   26695.83
## occupCraft and related trades workers                           26705.67
## occupElementary occupations                                     26708.68
## occupManagers                                                   26706.28
## occupOther: Not in paid work                                    26749.66
## occupPlant and machine operators, and assemblers                26720.37
## occupProfessionals                                              26691.05
## occupRetired                                                    26692.92
## occupService and sales workers                                  26697.35
## occupSkilled agricultural, forestry and fishery workers         26713.20
## occupTechnicians and associate professionals                    26695.98
## occupUnemployed                                                 26711.21
## environ.lvl1                                                       59.80
## polintr.lvl1                                                      109.83
## all.parties.lvl2Did not vote                                      174.89
## all.parties.lvl2Don't know                                        232.60
## all.parties.lvl2Invalid vote                                    14046.40
## all.parties.lvl2NE age                                            243.52
## all.parties.lvl2NE citizen                                        238.26
## all.parties.lvl2NE other                                          579.17
## all.parties.lvl2No answer                                        2037.94
## all.parties.lvl2Other party                                       212.99
## all.parties.lvl2Pro-environment party                             232.26
## environ.lvl1:polintr.lvl1                                         166.38
## environ.lvl1:all.parties.lvl2Did not vote                          83.41
## environ.lvl1:all.parties.lvl2Don't know                           318.90
## environ.lvl1:all.parties.lvl2Invalid vote                       26379.56
## environ.lvl1:all.parties.lvl2NE age                               276.51
## environ.lvl1:all.parties.lvl2NE citizen                           215.81
## environ.lvl1:all.parties.lvl2NE other                            1131.23
## environ.lvl1:all.parties.lvl2No answer                          26035.61
## environ.lvl1:all.parties.lvl2Other party                          129.85
## environ.lvl1:all.parties.lvl2Pro-environment party                237.65
## polintr.lvl1:all.parties.lvl2Did not vote                          73.83
## polintr.lvl1:all.parties.lvl2Don't know                           314.82
## polintr.lvl1:all.parties.lvl2Invalid vote                       25520.36
## polintr.lvl1:all.parties.lvl2NE age                               224.08
## polintr.lvl1:all.parties.lvl2NE citizen                           171.16
## polintr.lvl1:all.parties.lvl2NE other                            1061.75
## polintr.lvl1:all.parties.lvl2No answer                          26132.10
## polintr.lvl1:all.parties.lvl2Other party                          112.31
## polintr.lvl1:all.parties.lvl2Pro-environment party                193.05
## environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote            124.05
## environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know              835.75
## environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote          26315.53
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE age                  786.26
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen              332.89
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE other               1766.17
## environ.lvl1:polintr.lvl1:all.parties.lvl2No answer             26334.04
## environ.lvl1:polintr.lvl1:all.parties.lvl2Other party             193.56
## environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party   581.17
##                                                                 t.value     p
## (Intercept)                                                       -6.23 0.000
## age                                                               -2.58 0.010
## gender                                                             6.84 0.000
## educ                                                              10.81 0.000
## resid                                                             -6.62 0.000
## occupClerical support workers                                     -0.24 0.812
## occupCraft and related trades workers                             -1.12 0.264
## occupElementary occupations                                        0.20 0.840
## occupManagers                                                      0.34 0.735
## occupOther: Not in paid work                                       1.54 0.123
## occupPlant and machine operators, and assemblers                  -0.67 0.502
## occupProfessionals                                                 1.45 0.148
## occupRetired                                                       0.44 0.663
## occupService and sales workers                                    -0.49 0.621
## occupSkilled agricultural, forestry and fishery workers           -0.16 0.875
## occupTechnicians and associate professionals                      -0.11 0.910
## occupUnemployed                                                    0.25 0.800
## environ.lvl1                                                       3.18 0.002
## polintr.lvl1                                                       0.17 0.865
## all.parties.lvl2Did not vote                                       6.79 0.000
## all.parties.lvl2Don't know                                         5.93 0.000
## all.parties.lvl2Invalid vote                                      -0.25 0.804
## all.parties.lvl2NE age                                             9.96 0.000
## all.parties.lvl2NE citizen                                        11.34 0.000
## all.parties.lvl2NE other                                           7.54 0.000
## all.parties.lvl2No answer                                          1.96 0.050
## all.parties.lvl2Other party                                       11.91 0.000
## all.parties.lvl2Pro-environment party                             15.19 0.000
## environ.lvl1:polintr.lvl1                                          1.27 0.205
## environ.lvl1:all.parties.lvl2Did not vote                          1.05 0.296
## environ.lvl1:all.parties.lvl2Don't know                            1.98 0.048
## environ.lvl1:all.parties.lvl2Invalid vote                          1.06 0.289
## environ.lvl1:all.parties.lvl2NE age                                3.08 0.002
## environ.lvl1:all.parties.lvl2NE citizen                            0.87 0.387
## environ.lvl1:all.parties.lvl2NE other                              0.09 0.928
## environ.lvl1:all.parties.lvl2No answer                             0.70 0.482
## environ.lvl1:all.parties.lvl2Other party                           3.39 0.001
## environ.lvl1:all.parties.lvl2Pro-environment party                 1.70 0.090
## polintr.lvl1:all.parties.lvl2Did not vote                          1.36 0.179
## polintr.lvl1:all.parties.lvl2Don't know                            1.50 0.135
## polintr.lvl1:all.parties.lvl2Invalid vote                          1.14 0.255
## polintr.lvl1:all.parties.lvl2NE age                                0.16 0.872
## polintr.lvl1:all.parties.lvl2NE citizen                           -0.78 0.434
## polintr.lvl1:all.parties.lvl2NE other                              1.21 0.226
## polintr.lvl1:all.parties.lvl2No answer                             0.63 0.529
## polintr.lvl1:all.parties.lvl2Other party                           3.09 0.003
## polintr.lvl1:all.parties.lvl2Pro-environment party                 3.88 0.000
## environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote             0.14 0.886
## environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know              -0.15 0.879
## environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote            -1.00 0.317
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE age                   0.89 0.374
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen               2.29 0.023
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE other                -0.69 0.491
## environ.lvl1:polintr.lvl1:all.parties.lvl2No answer               -0.45 0.655
## environ.lvl1:polintr.lvl1:all.parties.lvl2Other party             -0.20 0.841
## environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party   -0.86 0.389
##                                                                     LL    UL
## (Intercept)                                                      -0.99 -0.52
## age                                                               0.00  0.00
## gender                                                            0.07  0.14
## educ                                                              0.02  0.03
## resid                                                            -0.13 -0.07
## occupClerical support workers                                    -0.24  0.19
## occupCraft and related trades workers                            -0.33  0.09
## occupElementary occupations                                      -0.19  0.23
## occupManagers                                                    -0.18  0.25
## occupOther: Not in paid work                                     -0.05  0.39
## occupPlant and machine operators, and assemblers                 -0.29  0.14
## occupProfessionals                                               -0.05  0.36
## occupRetired                                                     -0.19  0.29
## occupService and sales workers                                   -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers          -0.24  0.21
## occupTechnicians and associate professionals                     -0.22  0.20
## occupUnemployed                                                  -0.23  0.29
## environ.lvl1                                                      0.05  0.24
## polintr.lvl1                                                     -0.06  0.07
## all.parties.lvl2Did not vote                                      0.44  0.80
## all.parties.lvl2Don't know                                        0.39  0.78
## all.parties.lvl2Invalid vote                                     -2.17  1.68
## all.parties.lvl2NE age                                            0.79  1.19
## all.parties.lvl2NE citizen                                        0.95  1.34
## all.parties.lvl2NE other                                          0.75  1.27
## all.parties.lvl2No answer                                         0.00  2.08
## all.parties.lvl2Other party                                       0.67  0.93
## all.parties.lvl2Pro-environment party                             1.15  1.49
## environ.lvl1:polintr.lvl1                                        -0.03  0.12
## environ.lvl1:all.parties.lvl2Did not vote                        -0.04  0.14
## environ.lvl1:all.parties.lvl2Don't know                           0.00  0.27
## environ.lvl1:all.parties.lvl2Invalid vote                        -1.54  5.17
## environ.lvl1:all.parties.lvl2NE age                               0.07  0.33
## environ.lvl1:all.parties.lvl2NE citizen                          -0.07  0.18
## environ.lvl1:all.parties.lvl2NE other                            -0.24  0.26
## environ.lvl1:all.parties.lvl2No answer                           -1.65  3.49
## environ.lvl1:all.parties.lvl2Other party                          0.06  0.21
## environ.lvl1:all.parties.lvl2Pro-environment party               -0.02  0.22
## polintr.lvl1:all.parties.lvl2Did not vote                        -0.03  0.14
## polintr.lvl1:all.parties.lvl2Don't know                          -0.03  0.21
## polintr.lvl1:all.parties.lvl2Invalid vote                        -0.81  3.07
## polintr.lvl1:all.parties.lvl2NE age                              -0.10  0.12
## polintr.lvl1:all.parties.lvl2NE citizen                          -0.15  0.07
## polintr.lvl1:all.parties.lvl2NE other                            -0.08  0.34
## polintr.lvl1:all.parties.lvl2No answer                           -1.55  3.02
## polintr.lvl1:all.parties.lvl2Other party                          0.04  0.18
## polintr.lvl1:all.parties.lvl2Pro-environment party                0.10  0.29
## environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote           -0.08  0.09
## environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know             -0.16  0.13
## environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote           -5.04  1.63
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE age                 -0.08  0.21
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen              0.02  0.27
## environ.lvl1:polintr.lvl1:all.parties.lvl2NE other               -0.37  0.18
## environ.lvl1:polintr.lvl1:all.parties.lvl2No answer             -11.92  7.49
## environ.lvl1:polintr.lvl1:all.parties.lvl2Other party            -0.08  0.07
## environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party  -0.18  0.07
```

```r
(VC.H5.exp.intr.mod6<-getVC(H5.exp.intr.mod6))
```

```
##             grp                      var1                      var2      est_SD
## 1  voting.group               (Intercept)                      <NA>  0.25767494
## 2  voting.group              environ.lvl1                      <NA>  0.07618001
## 3  voting.group              polintr.lvl1                      <NA>  0.06498842
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.04515735
## 5  voting.group               (Intercept)              environ.lvl1  0.49810096
## 6  voting.group               (Intercept)              polintr.lvl1  0.39897315
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1 -0.52604503
## 8  voting.group              environ.lvl1              polintr.lvl1 -0.25425445
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.99619799
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1  0.29274361
## 11        cntry              environ.lvl1                      <NA>  0.11379091
## 12        cntry              polintr.lvl1                      <NA>  0.05756919
## 13        cntry environ.lvl1:polintr.lvl1                      <NA>  0.03844036
## 14        cntry              environ.lvl1              polintr.lvl1  0.17453607
## 15        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.94516567
## 16        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.48654385
## 17     Residual                      <NA>                      <NA>  1.15602079
##         est_SD2
## 1   0.066396375
## 2   0.005803394
## 3   0.004223494
## 4   0.002039186
## 5   0.009777562
## 6   0.006681159
## 7  -0.006121016
## 8  -0.001258768
## 9  -0.003427008
## 10  0.000859116
## 11  0.012948371
## 12  0.003314212
## 13  0.001477661
## 14  0.001143360
## 15  0.004134309
## 16  0.001076712
## 17  1.336384057
```

#### Refit with manually coded level-1 interaction


```r
dat.H5.intr$env.intr.int<-dat.H5.intr$environ.lvl1*dat.H5.intr$polintr.lvl1

H5.exp.intr.mod6<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+env.intr.int|voting.group)+
                (0+environ.lvl1+polintr.lvl1+env.intr.int|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                env.intr.int+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polintr.lvl1:all.parties.lvl2+
                env.intr.int:all.parties.lvl2
                ,data=dat.H5.intr,REML=F,
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

```
## Warning: Model failed to converge with 2 negative eigenvalues: -2.3e+00 -2.0e+03
```

```r
isSingular(H5.exp.intr.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.intr.mod5,H5.exp.intr.mod6)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod5: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     voting.group) + (0 + environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
## H5.exp.intr.mod6: refugees ~ (environ.lvl1 + polintr.lvl1 + env.intr.int | voting.group) + 
## H5.exp.intr.mod6:     (0 + environ.lvl1 + polintr.lvl1 + env.intr.int | cntry) + 
## H5.exp.intr.mod6:     age + gender + educ + resid + occup + environ.lvl1 + polintr.lvl1 + 
## H5.exp.intr.mod6:     env.intr.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.exp.intr.mod6:     polintr.lvl1:all.parties.lvl2 + env.intr.int:all.parties.lvl2
##                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.intr.mod5   64 84697 85222 -42285    84569                     
## H5.exp.intr.mod6   73 84713 85311 -42283    84567 2.6279  9     0.9772
```

```r
(FE.H5.exp.intr.mod6<-getFE(H5.exp.intr.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.76       0.12
## age                                                         0.00       0.00
## gender                                                      0.11       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.10       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.12       0.11
## occupElementary occupations                                 0.02       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.16       0.11
## occupRetired                                                0.06       0.12
## occupService and sales workers                             -0.05       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.11
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.03       0.13
## environ.lvl1                                                0.14       0.04
## polintr.lvl1                                                0.01       0.03
## env.intr.int                                                0.05       0.04
## all.parties.lvl2Did not vote                                0.62       0.09
## all.parties.lvl2Don't know                                  0.58       0.10
## all.parties.lvl2Invalid vote                               -0.24       0.98
## all.parties.lvl2NE age                                      0.99       0.10
## all.parties.lvl2NE citizen                                  1.15       0.10
## all.parties.lvl2NE other                                    1.01       0.13
## all.parties.lvl2No answer                                   1.04       0.53
## all.parties.lvl2Other party                                 0.80       0.07
## all.parties.lvl2Pro-environment party                       1.32       0.09
## environ.lvl1:all.parties.lvl2Did not vote                   0.06       0.04
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.06
## environ.lvl1:all.parties.lvl2Invalid vote                   1.81       1.71
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.06
## environ.lvl1:all.parties.lvl2NE citizen                     0.05       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.02       0.12
## environ.lvl1:all.parties.lvl2No answer                      0.91       1.31
## environ.lvl1:all.parties.lvl2Other party                    0.13       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.11       0.05
## polintr.lvl1:all.parties.lvl2Did not vote                   0.05       0.04
## polintr.lvl1:all.parties.lvl2Don't know                     0.09       0.06
## polintr.lvl1:all.parties.lvl2Invalid vote                   1.13       0.99
## polintr.lvl1:all.parties.lvl2NE age                         0.01       0.06
## polintr.lvl1:all.parties.lvl2NE citizen                    -0.04       0.05
## polintr.lvl1:all.parties.lvl2NE other                       0.13       0.11
## polintr.lvl1:all.parties.lvl2No answer                      0.74       1.17
## polintr.lvl1:all.parties.lvl2Other party                    0.11       0.03
## polintr.lvl1:all.parties.lvl2Pro-environment party          0.19       0.05
## env.intr.int:all.parties.lvl2Did not vote                   0.00       0.04
## env.intr.int:all.parties.lvl2Don't know                    -0.01       0.07
## env.intr.int:all.parties.lvl2Invalid vote                  -1.70       1.70
## env.intr.int:all.parties.lvl2NE age                         0.06       0.07
## env.intr.int:all.parties.lvl2NE citizen                     0.14       0.06
## env.intr.int:all.parties.lvl2NE other                      -0.10       0.14
## env.intr.int:all.parties.lvl2No answer                     -2.21       4.96
## env.intr.int:all.parties.lvl2Other party                    0.00       0.04
## env.intr.int:all.parties.lvl2Pro-environment party         -0.06       0.06
##                                                               df t.value     p
## (Intercept)                                              3029.00   -6.25 0.000
## age                                                     26678.08   -2.49 0.013
## gender                                                  26709.67    6.85 0.000
## educ                                                    25352.10   10.78 0.000
## resid                                                   26733.51   -6.59 0.000
## occupClerical support workers                           26676.61   -0.22 0.828
## occupCraft and related trades workers                   26692.22   -1.09 0.275
## occupElementary occupations                             26693.94    0.23 0.818
## occupManagers                                           26700.11    0.35 0.724
## occupOther: Not in paid work                            26733.93    1.57 0.117
## occupPlant and machine operators, and assemblers        26706.70   -0.66 0.510
## occupProfessionals                                      26680.58    1.47 0.141
## occupRetired                                            26678.91    0.46 0.648
## occupService and sales workers                          26683.16   -0.48 0.634
## occupSkilled agricultural, forestry and fishery workers 26690.37   -0.13 0.898
## occupTechnicians and associate professionals            26684.14   -0.10 0.921
## occupUnemployed                                         26705.82    0.25 0.801
## environ.lvl1                                               55.12    3.36 0.001
## polintr.lvl1                                              108.39    0.16 0.873
## env.intr.int                                               89.58    1.37 0.174
## all.parties.lvl2Did not vote                              174.74    6.80 0.000
## all.parties.lvl2Don't know                                232.54    5.93 0.000
## all.parties.lvl2Invalid vote                            14073.81   -0.25 0.804
## all.parties.lvl2NE age                                    243.50    9.97 0.000
## all.parties.lvl2NE citizen                                238.27   11.35 0.000
## all.parties.lvl2NE other                                  579.68    7.53 0.000
## all.parties.lvl2No answer                                2044.86    1.96 0.051
## all.parties.lvl2Other party                               212.80   11.92 0.000
## all.parties.lvl2Pro-environment party                     232.19   15.20 0.000
## environ.lvl1:all.parties.lvl2Did not vote                 381.58    1.45 0.149
## environ.lvl1:all.parties.lvl2Don't know                  2172.15    2.13 0.033
## environ.lvl1:all.parties.lvl2Invalid vote               26345.43    1.06 0.289
## environ.lvl1:all.parties.lvl2NE age                      1737.69    3.45 0.001
## environ.lvl1:all.parties.lvl2NE citizen                  1353.49    0.91 0.361
## environ.lvl1:all.parties.lvl2NE other                    8965.88    0.18 0.856
## environ.lvl1:all.parties.lvl2No answer                  26293.17    0.70 0.485
## environ.lvl1:all.parties.lvl2Other party                  801.76    3.63 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party       1682.91    1.99 0.046
## polintr.lvl1:all.parties.lvl2Did not vote                  73.04    1.36 0.178
## polintr.lvl1:all.parties.lvl2Don't know                   315.72    1.52 0.130
## polintr.lvl1:all.parties.lvl2Invalid vote               25601.37    1.14 0.255
## polintr.lvl1:all.parties.lvl2NE age                       223.99    0.15 0.882
## polintr.lvl1:all.parties.lvl2NE citizen                   168.89   -0.68 0.496
## polintr.lvl1:all.parties.lvl2NE other                    1070.00    1.18 0.237
## polintr.lvl1:all.parties.lvl2No answer                  26187.29    0.63 0.529
## polintr.lvl1:all.parties.lvl2Other party                  111.47    3.20 0.002
## polintr.lvl1:all.parties.lvl2Pro-environment party        192.29    3.89 0.000
## env.intr.int:all.parties.lvl2Did not vote                  60.83   -0.05 0.962
## env.intr.int:all.parties.lvl2Don't know                   397.00   -0.18 0.859
## env.intr.int:all.parties.lvl2Invalid vote               26442.94   -1.00 0.317
## env.intr.int:all.parties.lvl2NE age                       445.65    0.80 0.426
## env.intr.int:all.parties.lvl2NE citizen                   168.11    2.18 0.030
## env.intr.int:all.parties.lvl2NE other                    1075.94   -0.74 0.459
## env.intr.int:all.parties.lvl2No answer                  26370.75   -0.45 0.655
## env.intr.int:all.parties.lvl2Other party                   88.89   -0.13 0.898
## env.intr.int:all.parties.lvl2Pro-environment party        296.52   -0.98 0.330
##                                                             LL    UL
## (Intercept)                                              -0.99 -0.52
## age                                                       0.00  0.00
## gender                                                    0.08  0.14
## educ                                                      0.02  0.03
## resid                                                    -0.13 -0.07
## occupClerical support workers                            -0.23  0.19
## occupCraft and related trades workers                    -0.33  0.09
## occupElementary occupations                              -0.19  0.24
## occupManagers                                            -0.17  0.25
## occupOther: Not in paid work                             -0.04  0.39
## occupPlant and machine operators, and assemblers         -0.28  0.14
## occupProfessionals                                       -0.05  0.36
## occupRetired                                             -0.18  0.30
## occupService and sales workers                           -0.26  0.16
## occupSkilled agricultural, forestry and fishery workers  -0.24  0.21
## occupTechnicians and associate professionals             -0.22  0.20
## occupUnemployed                                          -0.23  0.29
## environ.lvl1                                              0.06  0.23
## polintr.lvl1                                             -0.06  0.07
## env.intr.int                                             -0.02  0.12
## all.parties.lvl2Did not vote                              0.44  0.80
## all.parties.lvl2Don't know                                0.39  0.77
## all.parties.lvl2Invalid vote                             -2.17  1.68
## all.parties.lvl2NE age                                    0.79  1.19
## all.parties.lvl2NE citizen                                0.95  1.34
## all.parties.lvl2NE other                                  0.75  1.27
## all.parties.lvl2No answer                                 0.00  2.08
## all.parties.lvl2Other party                               0.67  0.93
## all.parties.lvl2Pro-environment party                     1.15  1.49
## environ.lvl1:all.parties.lvl2Did not vote                -0.02  0.14
## environ.lvl1:all.parties.lvl2Don't know                   0.01  0.26
## environ.lvl1:all.parties.lvl2Invalid vote                -1.54  5.17
## environ.lvl1:all.parties.lvl2NE age                       0.09  0.33
## environ.lvl1:all.parties.lvl2NE citizen                  -0.06  0.17
## environ.lvl1:all.parties.lvl2NE other                    -0.22  0.26
## environ.lvl1:all.parties.lvl2No answer                   -1.65  3.48
## environ.lvl1:all.parties.lvl2Other party                  0.06  0.20
## environ.lvl1:all.parties.lvl2Pro-environment party        0.00  0.21
## polintr.lvl1:all.parties.lvl2Did not vote                -0.03  0.14
## polintr.lvl1:all.parties.lvl2Don't know                  -0.03  0.21
## polintr.lvl1:all.parties.lvl2Invalid vote                -0.81  3.07
## polintr.lvl1:all.parties.lvl2NE age                      -0.10  0.12
## polintr.lvl1:all.parties.lvl2NE citizen                  -0.14  0.07
## polintr.lvl1:all.parties.lvl2NE other                    -0.08  0.34
## polintr.lvl1:all.parties.lvl2No answer                   -1.55  3.03
## polintr.lvl1:all.parties.lvl2Other party                  0.04  0.18
## polintr.lvl1:all.parties.lvl2Pro-environment party        0.10  0.29
## env.intr.int:all.parties.lvl2Did not vote                -0.09  0.09
## env.intr.int:all.parties.lvl2Don't know                  -0.16  0.13
## env.intr.int:all.parties.lvl2Invalid vote                -5.04  1.63
## env.intr.int:all.parties.lvl2NE age                      -0.09  0.20
## env.intr.int:all.parties.lvl2NE citizen                   0.01  0.26
## env.intr.int:all.parties.lvl2NE other                    -0.38  0.17
## env.intr.int:all.parties.lvl2No answer                  -11.93  7.50
## env.intr.int:all.parties.lvl2Other party                 -0.08  0.07
## env.intr.int:all.parties.lvl2Pro-environment party       -0.19  0.06
```

```r
(VC.H5.exp.intr.mod6<-getVC(H5.exp.intr.mod6))
```

```
##             grp         var1         var2      est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.25740792  0.0662588350
## 2  voting.group environ.lvl1         <NA>  0.04132206  0.0017075126
## 3  voting.group polintr.lvl1         <NA>  0.06294352  0.0039618861
## 4  voting.group env.intr.int         <NA>  0.04426321  0.0019592322
## 5  voting.group  (Intercept) environ.lvl1  1.00000000  0.0106366253
## 6  voting.group  (Intercept) polintr.lvl1  0.39985226  0.0064784700
## 7  voting.group  (Intercept) env.intr.int -0.60053605 -0.0068423288
## 8  voting.group environ.lvl1 polintr.lvl1  0.39985226  0.0010399980
## 9  voting.group environ.lvl1 env.intr.int -0.60053605 -0.0010984088
## 10 voting.group polintr.lvl1 env.intr.int  0.21688476  0.0006042588
## 11        cntry environ.lvl1         <NA>  0.11527209  0.0132876554
## 12        cntry polintr.lvl1         <NA>  0.05782657  0.0033439119
## 13        cntry env.intr.int         <NA>  0.03881572  0.0015066604
## 14        cntry environ.lvl1 polintr.lvl1  0.14940596  0.0009959086
## 15        cntry environ.lvl1 env.intr.int  0.92913159  0.0041572782
## 16        cntry polintr.lvl1 env.intr.int  0.50441693  0.0011322042
## 17     Residual         <NA>         <NA>  1.15688300  1.3383782805
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.exp.intr.mod6.trends<-emtrends(H5.exp.intr.mod6,specs = c("all.parties.lvl2"),var=c("env.intr.int"))
(H5.exp.intr.mod6.trends.tab<-data.frame(H5.exp.intr.mod6.trends))
```

```
##          all.parties.lvl2 env.intr.int.trend         SE  df     asymp.LCL
## 1  Anti-immigration party         0.04986683 0.03639405 Inf  -0.021464195
## 2            Did not vote         0.04777783 0.02871579 Inf  -0.008504081
## 3              Don't know         0.03665599 0.06611392 Inf  -0.092924924
## 4            Invalid vote        -1.65458669 1.70202498 Inf  -4.990494343
## 5                  NE age         0.10895220 0.06614026 Inf  -0.020680324
## 6              NE citizen         0.18652934 0.05267541 Inf   0.083287430
## 7                NE other        -0.05351129 0.13537753 Inf  -0.318846377
## 8               No answer        -2.16477552 4.95574468 Inf -11.877856599
## 9             Other party         0.04488365 0.01929926 Inf   0.007057793
## 10  Pro-environment party        -0.01262596 0.05441715 Inf  -0.119281623
##     asymp.UCL
## 1  0.12119786
## 2  0.10405974
## 3  0.16623690
## 4  1.68132097
## 5  0.23858472
## 6  0.28977126
## 7  0.21182379
## 8  7.54830557
## 9  0.08270951
## 10 0.09402970
```

```r
H5.exp.intr.mod6.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.intr.mod6.trends.tab$env.intr.int.trend/
                   H5.exp.intr.mod6.trends.tab$SE)))
H5.exp.intr.mod6.trends.tab$adj.p<-
  p.adjust(H5.exp.intr.mod6.trends.tab$p,method="holm")

H5.exp.intr.mod6.trends.tab<-
  cbind(group=H5.exp.intr.mod6.trends.tab[,1],
      round(H5.exp.intr.mod6.trends.tab[,c(2,3)],2),
      round(H5.exp.intr.mod6.trends.tab[,c(7,8)],4),
      round(H5.exp.intr.mod6.trends.tab[,c(5,6)],2))
H5.exp.intr.mod6.trends.tab
```

```
##                     group env.intr.int.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party               0.05 0.04 0.1706 1.0000     -0.02
## 2            Did not vote               0.05 0.03 0.0961 0.7692     -0.01
## 3              Don't know               0.04 0.07 0.5793 1.0000     -0.09
## 4            Invalid vote              -1.65 1.70 0.3310 1.0000     -4.99
## 5                  NE age               0.11 0.07 0.0995 0.7692     -0.02
## 6              NE citizen               0.19 0.05 0.0004 0.0040      0.08
## 7                NE other              -0.05 0.14 0.6926 1.0000     -0.32
## 8               No answer              -2.16 4.96 0.6622 1.0000    -11.88
## 9             Other party               0.04 0.02 0.0200 0.1803      0.01
## 10  Pro-environment party              -0.01 0.05 0.8165 1.0000     -0.12
##    asymp.UCL
## 1       0.12
## 2       0.10
## 3       0.17
## 4       1.68
## 5       0.24
## 6       0.29
## 7       0.21
## 8       7.55
## 9       0.08
## 10      0.09
```

```r
write.csv2(H5.exp.intr.mod6.trends.tab,"H5.exp.intr.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1 Pro-environment party - Anti-immigration party -0.06249279 0.06400678 Inf
##      z.ratio   p.value
## 1 -0.9763464 0.3288928
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.intr.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect    0.434 0.584 Inf  0.744  1.0000 
##  Did not vote effect              0.432 0.583 Inf  0.741  1.0000 
##  Don't know effect                0.420 0.586 Inf  0.716  1.0000 
##  Invalid vote effect             -1.459 1.789 Inf -0.816  1.0000 
##  NE age effect                    0.500 0.586 Inf  0.853  1.0000 
##  NE citizen effect                0.586 0.585 Inf  1.002  1.0000 
##  NE other effect                  0.320 0.598 Inf  0.534  1.0000 
##  No answer effect                -2.026 4.959 Inf -0.409  1.0000 
##  Other party effect               0.429 0.583 Inf  0.736  1.0000 
##  Pro-environment party effect     0.365 0.585 Inf  0.624  1.0000 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H5.exp.intr.more.contrasts<-data.frame(pairs(H5.exp.intr.mod6.trends, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast     estimate         SE  df
## 1           Other party - Anti-immigration party -0.004983177 0.03882330 Inf
## 2 Pro-environment party - Anti-immigration party -0.062492794 0.06400678 Inf
## 3            Pro-environment party - Other party -0.057509617 0.05604240 Inf
##      z.ratio   p.value
## 1 -0.1283553 0.9144194
## 2 -0.9763464 0.9144194
## 3 -1.0261806 0.9144194
```

