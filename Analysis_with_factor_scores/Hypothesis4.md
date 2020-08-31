---
title: "Hypothesis 4"
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




# Hypothesis 4: The strength of the association between environment and refugee attitudes is stronger among those who voted for pro-environment or anti-immigration parties in the previous national elections

### Model 1: without interactions (only main effects, combine H1 and H2 final models in terms of predictors)


```r
H4.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                all.parties.lvl2,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


isSingular(H4.mod1)
```

```
## [1] FALSE
```

```r
(FE.H4.mod1<-getFE(H4.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.70       0.12
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.02       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.05       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.07       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.26       0.03
## all.parties.lvl2Did not vote                                0.63       0.09
## all.parties.lvl2Don't know                                  0.53       0.10
## all.parties.lvl2Invalid vote                                0.72       0.49
## all.parties.lvl2NE age                                      0.95       0.10
## all.parties.lvl2NE citizen                                  1.14       0.10
## all.parties.lvl2NE other                                    0.97       0.13
## all.parties.lvl2No answer                                   0.99       0.48
## all.parties.lvl2Other party                                 0.73       0.06
## all.parties.lvl2Pro-environment party                       1.26       0.08
##                                                               df t.value     p
## (Intercept)                                              3281.55   -5.84 0.000
## age                                                     26833.17   -1.13 0.258
## gender                                                  26776.47    5.59 0.000
## educ                                                    26557.78   12.24 0.000
## resid                                                   26854.89   -7.05 0.000
## occupClerical support workers                           26759.43   -0.22 0.823
## occupCraft and related trades workers                   26767.26   -1.19 0.236
## occupElementary occupations                             26770.77    0.10 0.922
## occupManagers                                           26763.15    0.43 0.664
## occupOther: Not in paid work                            26804.87    1.54 0.124
## occupPlant and machine operators, and assemblers        26771.37   -0.68 0.494
## occupProfessionals                                      26754.96    1.58 0.114
## occupRetired                                            26758.98    0.45 0.655
## occupService and sales workers                          26758.29   -0.55 0.579
## occupSkilled agricultural, forestry and fishery workers 26765.86   -0.19 0.852
## occupTechnicians and associate professionals            26755.61   -0.09 0.926
## occupUnemployed                                         26771.92    0.09 0.930
## environ.lvl1                                               15.81    8.11 0.000
## all.parties.lvl2Did not vote                              146.35    7.32 0.000
## all.parties.lvl2Don't know                                236.69    5.58 0.000
## all.parties.lvl2Invalid vote                             1503.56    1.47 0.141
## all.parties.lvl2NE age                                    243.88    9.88 0.000
## all.parties.lvl2NE citizen                                233.04   11.70 0.000
## all.parties.lvl2NE other                                  600.86    7.32 0.000
## all.parties.lvl2No answer                                1601.97    2.06 0.040
## all.parties.lvl2Other party                               201.20   11.31 0.000
## all.parties.lvl2Pro-environment party                     229.20   14.89 0.000
##                                                            LL    UL
## (Intercept)                                             -0.94 -0.47
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.38
## occupRetired                                            -0.19  0.30
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.20  0.33
## all.parties.lvl2Did not vote                             0.46  0.80
## all.parties.lvl2Don't know                               0.35  0.72
## all.parties.lvl2Invalid vote                            -0.24  1.68
## all.parties.lvl2NE age                                   0.76  1.14
## all.parties.lvl2NE citizen                               0.95  1.34
## all.parties.lvl2NE other                                 0.71  1.23
## all.parties.lvl2No answer                                0.05  1.94
## all.parties.lvl2Other party                              0.60  0.86
## all.parties.lvl2Pro-environment party                    1.09  1.42
```

```r
(VC.H4.mod1<-getVC(H4.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.25782789 0.066475222
## 2 voting.group environ.lvl1         <NA> 0.09412447 0.008859416
## 3 voting.group  (Intercept) environ.lvl1 0.55797450 0.013540877
## 4        cntry environ.lvl1         <NA> 0.11868378 0.014085839
## 5     Residual         <NA>         <NA> 1.16255916 1.351543809
```

\newpage

### Model 2: Cross-level interaction between environmental attitudes and voting group categories


```r
H4.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (0+environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
isSingular(H4.mod2)
```

```
## [1] FALSE
```

```r
anova(H4.mod1,H4.mod2)
```

```
## Data: dat
## Models:
## H4.mod1: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H4.mod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H4.mod1:     all.parties.lvl2
## H4.mod2: refugees ~ (environ.lvl1 | voting.group) + (0 + environ.lvl1 | 
## H4.mod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H4.mod2:     all.parties.lvl2 + environ.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H4.mod1   32 84915 85178 -42426    84851                        
## H4.mod2   41 84910 85246 -42414    84828 23.641  9   0.004906 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H4.mod2<-getFE(H4.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.75       0.12
## age                                                         0.00       0.00
## gender                                                      0.09       0.02
## educ                                                        0.03       0.00
## resid                                                      -0.11       0.02
## occupClerical support workers                              -0.03       0.11
## occupCraft and related trades workers                      -0.13       0.11
## occupElementary occupations                                 0.01       0.11
## occupManagers                                               0.04       0.11
## occupOther: Not in paid work                                0.17       0.11
## occupPlant and machine operators, and assemblers           -0.08       0.11
## occupProfessionals                                          0.17       0.11
## occupRetired                                                0.05       0.12
## occupService and sales workers                             -0.06       0.11
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.12
## occupTechnicians and associate professionals               -0.01       0.11
## occupUnemployed                                             0.01       0.13
## environ.lvl1                                                0.15       0.05
## all.parties.lvl2Did not vote                                0.63       0.09
## all.parties.lvl2Don't know                                  0.59       0.10
## all.parties.lvl2Invalid vote                                0.77       0.49
## all.parties.lvl2NE age                                      1.03       0.10
## all.parties.lvl2NE citizen                                  1.17       0.10
## all.parties.lvl2NE other                                    1.02       0.13
## all.parties.lvl2No answer                                   1.09       0.50
## all.parties.lvl2Other party                                 0.80       0.07
## all.parties.lvl2Pro-environment party                       1.31       0.09
## environ.lvl1:all.parties.lvl2Did not vote                   0.05       0.05
## environ.lvl1:all.parties.lvl2Don't know                     0.14       0.07
## environ.lvl1:all.parties.lvl2Invalid vote                   0.08       0.74
## environ.lvl1:all.parties.lvl2NE age                         0.21       0.07
## environ.lvl1:all.parties.lvl2NE citizen                     0.03       0.06
## environ.lvl1:all.parties.lvl2NE other                       0.04       0.13
## environ.lvl1:all.parties.lvl2No answer                      0.32       0.58
## environ.lvl1:all.parties.lvl2Other party                    0.15       0.04
## environ.lvl1:all.parties.lvl2Pro-environment party          0.14       0.06
##                                                               df t.value     p
## (Intercept)                                              3103.52   -6.19 0.000
## age                                                     26835.83   -1.10 0.272
## gender                                                  26781.09    5.60 0.000
## educ                                                    26554.04   12.24 0.000
## resid                                                   26859.61   -7.04 0.000
## occupClerical support workers                           26763.63   -0.25 0.805
## occupCraft and related trades workers                   26772.22   -1.21 0.228
## occupElementary occupations                             26775.77    0.08 0.938
## occupManagers                                           26767.39    0.40 0.686
## occupOther: Not in paid work                            26808.03    1.53 0.127
## occupPlant and machine operators, and assemblers        26775.44   -0.71 0.479
## occupProfessionals                                      26758.35    1.55 0.122
## occupRetired                                            26762.19    0.41 0.684
## occupService and sales workers                          26763.31   -0.57 0.566
## occupSkilled agricultural, forestry and fishery workers 26771.95   -0.19 0.849
## occupTechnicians and associate professionals            26759.68   -0.12 0.908
## occupUnemployed                                         26779.15    0.06 0.951
## environ.lvl1                                               56.52    3.19 0.002
## all.parties.lvl2Did not vote                              175.97    6.96 0.000
## all.parties.lvl2Don't know                                235.10    6.02 0.000
## all.parties.lvl2Invalid vote                             1526.33    1.58 0.115
## all.parties.lvl2NE age                                    244.11   10.42 0.000
## all.parties.lvl2NE citizen                                238.83   11.60 0.000
## all.parties.lvl2NE other                                  581.77    7.63 0.000
## all.parties.lvl2No answer                                1660.74    2.18 0.030
## all.parties.lvl2Other party                               214.64   11.91 0.000
## all.parties.lvl2Pro-environment party                     232.62   15.15 0.000
## environ.lvl1:all.parties.lvl2Did not vote                  80.20    1.15 0.253
## environ.lvl1:all.parties.lvl2Don't know                   303.49    2.02 0.045
## environ.lvl1:all.parties.lvl2Invalid vote               20555.74    0.11 0.909
## environ.lvl1:all.parties.lvl2NE age                       264.49    3.17 0.002
## environ.lvl1:all.parties.lvl2NE citizen                   204.41    0.50 0.618
## environ.lvl1:all.parties.lvl2NE other                    1058.69    0.28 0.778
## environ.lvl1:all.parties.lvl2No answer                  14169.85    0.55 0.580
## environ.lvl1:all.parties.lvl2Other party                  123.07    3.76 0.000
## environ.lvl1:all.parties.lvl2Pro-environment party        223.94    2.35 0.020
##                                                            LL    UL
## (Intercept)                                             -0.99 -0.51
## age                                                      0.00  0.00
## gender                                                   0.06  0.12
## educ                                                     0.02  0.03
## resid                                                   -0.14 -0.08
## occupClerical support workers                           -0.24  0.19
## occupCraft and related trades workers                   -0.34  0.08
## occupElementary occupations                             -0.20  0.22
## occupManagers                                           -0.17  0.26
## occupOther: Not in paid work                            -0.05  0.39
## occupPlant and machine operators, and assemblers        -0.29  0.14
## occupProfessionals                                      -0.04  0.37
## occupRetired                                            -0.19  0.29
## occupService and sales workers                          -0.27  0.15
## occupSkilled agricultural, forestry and fishery workers -0.25  0.20
## occupTechnicians and associate professionals            -0.22  0.20
## occupUnemployed                                         -0.25  0.27
## environ.lvl1                                             0.06  0.24
## all.parties.lvl2Did not vote                             0.45  0.81
## all.parties.lvl2Don't know                               0.40  0.78
## all.parties.lvl2Invalid vote                            -0.19  1.73
## all.parties.lvl2NE age                                   0.83  1.22
## all.parties.lvl2NE citizen                               0.97  1.36
## all.parties.lvl2NE other                                 0.76  1.28
## all.parties.lvl2No answer                                0.11  2.07
## all.parties.lvl2Other party                              0.66  0.93
## all.parties.lvl2Pro-environment party                    1.14  1.48
## environ.lvl1:all.parties.lvl2Did not vote               -0.04  0.15
## environ.lvl1:all.parties.lvl2Don't know                  0.00  0.27
## environ.lvl1:all.parties.lvl2Invalid vote               -1.36  1.53
## environ.lvl1:all.parties.lvl2NE age                      0.08  0.34
## environ.lvl1:all.parties.lvl2NE citizen                 -0.09  0.16
## environ.lvl1:all.parties.lvl2NE other                   -0.21  0.28
## environ.lvl1:all.parties.lvl2No answer                  -0.81  1.45
## environ.lvl1:all.parties.lvl2Other party                 0.07  0.23
## environ.lvl1:all.parties.lvl2Pro-environment party       0.02  0.25
```

```r
(VC.H4.mod2<-getVC(H4.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.25633692 0.065708617
## 2 voting.group environ.lvl1         <NA> 0.07592881 0.005765184
## 3 voting.group  (Intercept) environ.lvl1 0.59871614 0.011653026
## 4        cntry environ.lvl1         <NA> 0.12043334 0.014504188
## 5     Residual         <NA>         <NA> 1.16243385 1.351252446
```

\newpage

#### Marginal effect for pro-environment and anti-immigration voters


```r
H4.mod2.trends<-emtrends(H4.mod2,specs = c("all.parties.lvl2"),var=c("environ.lvl1"))
(H4.mod2.trends.tab<-data.frame(H4.mod2.trends))
```

```
##          all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL
## 1  Anti-immigration party          0.1501048 0.04705690 Inf  0.05787496
## 2            Did not vote          0.2040939 0.04289907 Inf  0.12001326
## 3              Don't know          0.2891729 0.06607949 Inf  0.15965950
## 4            Invalid vote          0.2342102 0.73837543 Inf -1.21297903
## 5                  NE age          0.3569415 0.06245242 Inf  0.23453705
## 6              NE citizen          0.1819405 0.06070804 Inf  0.06295487
## 7                NE other          0.1858491 0.12518458 Inf -0.05950812
## 8               No answer          0.4683672 0.57467837 Inf -0.65798171
## 9             Other party          0.3003528 0.03483645 Inf  0.23207460
## 10  Pro-environment party          0.2882629 0.05550016 Inf  0.17948463
##    asymp.UCL
## 1  0.2423346
## 2  0.2881745
## 3  0.4186864
## 4  1.6813995
## 5  0.4793460
## 6  0.3009260
## 7  0.4312064
## 8  1.5947161
## 9  0.3686310
## 10 0.3970413
```

```r
H4.mod2.trends.tab$p<-
  2*(1-pnorm(abs(H4.mod2.trends.tab$environ.lvl1.trend/
                   H4.mod2.trends.tab$SE)))
H4.mod2.trends.tab$adj.p<-
  p.adjust(H4.mod2.trends.tab$p,method="holm")

H4.mod2.trends.tab<-
  cbind(group=H4.mod2.trends.tab[,1],
      round(H4.mod2.trends.tab[,c(2,3)],2),
      round(H4.mod2.trends.tab[,c(7,8)],4),
      round(H4.mod2.trends.tab[,c(5,6)],2))
H4.mod2.trends.tab
```

```
##                     group environ.lvl1.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party               0.15 0.05 0.0014 0.0071      0.06
## 2            Did not vote               0.20 0.04 0.0000 0.0000      0.12
## 3              Don't know               0.29 0.07 0.0000 0.0001      0.16
## 4            Invalid vote               0.23 0.74 0.7511 0.8301     -1.21
## 5                  NE age               0.36 0.06 0.0000 0.0000      0.23
## 6              NE citizen               0.18 0.06 0.0027 0.0109      0.06
## 7                NE other               0.19 0.13 0.1376 0.4129     -0.06
## 8               No answer               0.47 0.57 0.4151 0.8301     -0.66
## 9             Other party               0.30 0.03 0.0000 0.0000      0.23
## 10  Pro-environment party               0.29 0.06 0.0000 0.0000      0.18
##    asymp.UCL
## 1       0.24
## 2       0.29
## 3       0.42
## 4       1.68
## 5       0.48
## 6       0.30
## 7       0.43
## 8       1.59
## 9       0.37
## 10      0.40
```

```r
write.csv2(H4.mod2.trends.tab,"H4.mod2.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H4.contrast<-data.frame(pairs(H4.mod2.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.1381582 0.05879975 Inf
##    z.ratio    p.value
## 1 2.349639 0.01879163
```

```r
#contrast for all groups against mean of other groups
contrast(H4.mod2.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate    SE  df z.ratio p.value
##  Anti-immigration party effect  -0.1287 0.111 Inf -1.156  0.2478 
##  Did not vote effect            -0.0687 0.110 Inf -0.626  0.5311 
##  Don't know effect               0.0258 0.121 Inf  0.214  0.8304 
##  Invalid vote effect            -0.0352 0.741 Inf -0.048  0.9621 
##  NE age effect                   0.1011 0.119 Inf  0.853  0.3938 
##  NE citizen effect              -0.0933 0.118 Inf -0.793  0.4280 
##  NE other effect                -0.0890 0.160 Inf -0.555  0.5789 
##  No answer effect                0.2249 0.580 Inf  0.388  0.6982 
##  Other party effect              0.0382 0.107 Inf  0.358  0.7202 
##  Pro-environment party effect    0.0248 0.115 Inf  0.215  0.8295 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic
```

```r
#contrast for three voting groups
(H4.more.contrasts<-data.frame(pairs(H4.mod2.trends, 
         exclude=c(2:8), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast    estimate         SE  df
## 1           Other party - Anti-immigration party  0.15024801 0.04000306 Inf
## 2 Pro-environment party - Anti-immigration party  0.13815817 0.05879975 Inf
## 3            Pro-environment party - Other party -0.01208984 0.04932082 Inf
##      z.ratio      p.value
## 1  3.7559133 0.0001727105
## 2  2.3496389 0.0187916323
## 3 -0.2451265 0.8063584866
```



\newpage
