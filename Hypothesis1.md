---
title: "Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_depth: 5
---


\newpage

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
#recode the first variable to represent this attitude
dat$environ<-dat$inctxff.R
describe(dat$environ,fast=T)
```

```
##    vars     n mean   sd min max range   se
## X1    1 36131 2.78 1.24   1   5     4 0.01
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
##    vars     n mean   sd   min max range   se
## X1    1 36131    0 1.18 -3.29   3  6.29 0.01
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
#grand mean center
dat$refugees<-dat$gvrfgap.R-mean(dat$gvrfgap.R,na.rm=T)
describe(dat$refugees,fast=T)
```

```
##    vars     n mean   sd   min  max range   se
## X1    1 36425    0 1.19 -1.93 2.07     4 0.01
```

```r
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
## X1    1 36425    0 1.04 -3.27 3.26  6.53 0.01
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
##     0     1     2 
## 35740  1076    60
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
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.08       0.11 19.92    0.69 0.499 -0.16 0.32
```

```r
(VC.H1.mod0<-getVC(H1.mod0))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.3247280 0.1054483
## 2        cntry (Intercept) <NA> 0.5023414 0.2523469
## 3     Residual        <NA> <NA> 1.0475008 1.0972579
```

```r
getDEV(H1.mod0)
```

```
## [1] 105426.7
```

```r
#ICC

##voting group

VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.07247039
```

```r
##country

VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.173428
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod0    4 105435 105469 -52713   105427                         
## H1.mod1   20 105044 105214 -52502   105004 422.87 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod1<-getFE(H1.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.06       0.14
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.08       0.01
## occupClerical support workers                              -0.04       0.09
## occupCraft and related trades workers                      -0.08       0.09
## occupElementary occupations                                 0.01       0.09
## occupManagers                                               0.00       0.09
## occupOther: Not in paid work                                0.16       0.09
## occupPlant and machine operators, and assemblers           -0.07       0.09
## occupProfessionals                                          0.10       0.09
## occupRetired                                               -0.03       0.10
## occupService and sales workers                             -0.04       0.09
## occupSkilled agricultural, forestry and fishery workers    -0.05       0.09
## occupTechnicians and associate professionals               -0.03       0.09
## occupUnemployed                                            -0.03       0.11
##                                                               df t.value     p
## (Intercept)                                                50.42    0.44 0.663
## age                                                     34095.35    3.30 0.001
## gender                                                  35598.09    4.77 0.000
## educ                                                    35706.37    9.28 0.000
## resid                                                   35657.25   -6.58 0.000
## occupClerical support workers                           35532.93   -0.40 0.687
## occupCraft and related trades workers                   35539.44   -0.85 0.398
## occupElementary occupations                             35540.85    0.14 0.892
## occupManagers                                           35533.47    0.00 0.998
## occupOther: Not in paid work                            35696.52    1.71 0.088
## occupPlant and machine operators, and assemblers        35539.06   -0.77 0.442
## occupProfessionals                                      35535.19    1.10 0.273
## occupRetired                                            35538.37   -0.28 0.782
## occupService and sales workers                          35534.90   -0.40 0.691
## occupSkilled agricultural, forestry and fishery workers 35544.87   -0.52 0.601
## occupTechnicians and associate professionals            35532.10   -0.29 0.774
## occupUnemployed                                         35546.28   -0.23 0.816
##                                                            LL    UL
## (Intercept)                                             -0.22  0.35
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.10 -0.05
## occupClerical support workers                           -0.21  0.14
## occupCraft and related trades workers                   -0.25  0.10
## occupElementary occupations                             -0.16  0.19
## occupManagers                                           -0.18  0.18
## occupOther: Not in paid work                            -0.02  0.34
## occupPlant and machine operators, and assemblers        -0.24  0.11
## occupProfessionals                                      -0.08  0.27
## occupRetired                                            -0.22  0.17
## occupService and sales workers                          -0.21  0.14
## occupSkilled agricultural, forestry and fishery workers -0.23  0.14
## occupTechnicians and associate professionals            -0.20  0.15
## occupUnemployed                                         -0.24  0.19
```

```r
(VC.H1.mod1<-getVC(H1.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.3034866 0.09210413
## 2        cntry (Intercept) <NA> 0.4986022 0.24860415
## 3     Residual        <NA> <NA> 1.0417310 1.08520339
```

```r
getDEV(H1.mod1)
```

```
## [1] 105003.8
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
## [1] 0.01098601
```

```r
##lvl 2: voting group

(VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.1265467
```

```r
##lvl 3: country

(VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="cntry","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]
```

```
## [1] 0.01483183
```

```r
##total

(sum(VC.H1.mod0$est_SD2)-sum(VC.H1.mod1$est_SD2))/
  sum(VC.H1.mod0$est_SD2)
```

```
## [1] 0.02002771
```

```r
#individual contributions of covariates
anova(H1.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age     11.807  11.807     1 34095 10.8797 0.0009732 ***
## gender  24.687  24.687     1 35598 22.7491 1.853e-06 ***
## educ    93.417  93.417     1 35706 86.0829 < 2.2e-16 ***
## resid   47.042  47.042     1 35657 43.3489 4.643e-11 ***
## occup  125.754  10.479    12 35223  9.6567 < 2.2e-16 ***
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod1   20 105044 105214 -52502   105004                         
## H1.mod2   21 104337 104515 -52147   104295 708.83  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod2<-getFE(H1.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.07       0.14
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.01       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.04       0.09
## occupCraft and related trades workers                      -0.07       0.09
## occupElementary occupations                                 0.01       0.09
## occupManagers                                              -0.02       0.09
## occupOther: Not in paid work                                0.14       0.09
## occupPlant and machine operators, and assemblers           -0.06       0.09
## occupProfessionals                                          0.07       0.09
## occupRetired                                               -0.03       0.10
## occupService and sales workers                             -0.04       0.09
## occupSkilled agricultural, forestry and fishery workers    -0.05       0.09
## occupTechnicians and associate professionals               -0.03       0.09
## occupUnemployed                                            -0.02       0.11
## environ.lvl1                                                0.13       0.00
##                                                               df t.value     p
## (Intercept)                                                49.59    0.51 0.611
## age                                                     34270.64    4.39 0.000
## gender                                                  35593.07    4.67 0.000
## educ                                                    35714.29    7.51 0.000
## resid                                                   35651.76   -5.23 0.000
## occupClerical support workers                           35528.94   -0.46 0.646
## occupCraft and related trades workers                   35535.31   -0.76 0.446
## occupElementary occupations                             35536.62    0.16 0.871
## occupManagers                                           35529.47   -0.18 0.859
## occupOther: Not in paid work                            35691.43    1.54 0.123
## occupPlant and machine operators, and assemblers        35534.89   -0.72 0.474
## occupProfessionals                                      35531.19    0.83 0.407
## occupRetired                                            35534.45   -0.34 0.734
## occupService and sales workers                          35530.91   -0.43 0.668
## occupSkilled agricultural, forestry and fishery workers 35540.59   -0.58 0.564
## occupTechnicians and associate professionals            35528.10   -0.36 0.718
## occupUnemployed                                         35542.75   -0.19 0.852
## environ.lvl1                                            35459.04   26.76 0.000
##                                                            LL    UL
## (Intercept)                                             -0.21  0.36
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.08 -0.04
## occupClerical support workers                           -0.21  0.13
## occupCraft and related trades workers                   -0.24  0.11
## occupElementary occupations                             -0.16  0.19
## occupManagers                                           -0.19  0.16
## occupOther: Not in paid work                            -0.04  0.32
## occupPlant and machine operators, and assemblers        -0.24  0.11
## occupProfessionals                                      -0.10  0.24
## occupRetired                                            -0.23  0.16
## occupService and sales workers                          -0.21  0.13
## occupSkilled agricultural, forestry and fishery workers -0.24  0.13
## occupTechnicians and associate professionals            -0.20  0.14
## occupUnemployed                                         -0.23  0.19
## environ.lvl1                                             0.12  0.13
```

```r
(VC.H1.mod2<-getVC(H1.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.3088935 0.09541521
## 2        cntry (Intercept) <NA> 0.4990552 0.24905606
## 3     Residual        <NA> <NA> 1.0312663 1.06351024
```

```r
getDEV(H1.mod2)
```

```
## [1] 104295
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
## [1] 0.01998994
```

```r
##total

(sum(VC.H1.mod1$est_SD2)-sum(VC.H1.mod2$est_SD2))/
  sum(VC.H1.mod1$est_SD2)
```

```
## [1] 0.01257453
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
##         npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)    
## H1.mod2   21 104337 104515 -52147   104295                        
## H1.mod3   25 104289 104502 -52120   104239 55.58  4  2.456e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod3<-getFE(H1.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.08       0.14
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.01       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.04       0.09
## occupCraft and related trades workers                      -0.07       0.09
## occupElementary occupations                                 0.01       0.09
## occupManagers                                              -0.02       0.09
## occupOther: Not in paid work                                0.14       0.09
## occupPlant and machine operators, and assemblers           -0.07       0.09
## occupProfessionals                                          0.07       0.09
## occupRetired                                               -0.04       0.10
## occupService and sales workers                             -0.04       0.09
## occupSkilled agricultural, forestry and fishery workers    -0.06       0.09
## occupTechnicians and associate professionals               -0.03       0.09
## occupUnemployed                                            -0.03       0.11
## environ.lvl1                                                0.12       0.01
##                                                               df t.value     p
## (Intercept)                                                49.47    0.53 0.598
## age                                                     34211.38    4.49 0.000
## gender                                                  35569.60    4.70 0.000
## educ                                                    35698.92    7.49 0.000
## resid                                                   35634.38   -5.24 0.000
## occupClerical support workers                           35508.21   -0.48 0.632
## occupCraft and related trades workers                   35518.38   -0.77 0.439
## occupElementary occupations                             35520.28    0.13 0.900
## occupManagers                                           35507.43   -0.22 0.827
## occupOther: Not in paid work                            35674.75    1.52 0.127
## occupPlant and machine operators, and assemblers        35516.07   -0.74 0.460
## occupProfessionals                                      35512.07    0.79 0.432
## occupRetired                                            35508.44   -0.36 0.716
## occupService and sales workers                          35513.92   -0.44 0.657
## occupSkilled agricultural, forestry and fishery workers 35523.04   -0.63 0.529
## occupTechnicians and associate professionals            35508.10   -0.38 0.705
## occupUnemployed                                         35525.32   -0.24 0.812
## environ.lvl1                                               19.40   11.53 0.000
##                                                            LL    UL
## (Intercept)                                             -0.21  0.36
## age                                                      0.00  0.00
## gender                                                   0.03  0.08
## educ                                                     0.01  0.02
## resid                                                   -0.08 -0.04
## occupClerical support workers                           -0.22  0.13
## occupCraft and related trades workers                   -0.24  0.10
## occupElementary occupations                             -0.16  0.18
## occupManagers                                           -0.19  0.15
## occupOther: Not in paid work                            -0.04  0.32
## occupPlant and machine operators, and assemblers        -0.24  0.11
## occupProfessionals                                      -0.10  0.24
## occupRetired                                            -0.23  0.16
## occupService and sales workers                          -0.21  0.13
## occupSkilled agricultural, forestry and fishery workers -0.24  0.12
## occupTechnicians and associate professionals            -0.20  0.14
## occupUnemployed                                         -0.24  0.18
## environ.lvl1                                             0.10  0.15
```

```r
(VC.H1.mod3<-getVC(H1.mod3))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30918363  0.0955945141
## 2 voting.group environ.lvl1         <NA>  0.04317550  0.0018641240
## 3 voting.group  (Intercept) environ.lvl1  0.12926379  0.0017255628
## 4        cntry  (Intercept)         <NA>  0.49932256  0.2493230223
## 5        cntry environ.lvl1         <NA>  0.04000565  0.0016004517
## 6        cntry  (Intercept) environ.lvl1 -0.04865068 -0.0009718325
## 7     Residual         <NA>         <NA>  1.02912455  1.0590973317
```

```r
getDEV(H1.mod3)
```

```
## [1] 104239.4
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
    "Czech Republic",
    "Germany",
    "Estonia",
    "Spain",
    "Finland",
    "France",
    "Great Britain",
    "Hungary",
    "Ireland",
    "Italy",
    "Lithuania",
    "Netherlands",
    "Norway",
    "Poland",
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
```

```
## Warning in cor(environ.gmc, refugees, use = "pairwise.complete.obs"): the
## standard deviation is zero

## Warning in cor(environ.gmc, refugees, use = "pairwise.complete.obs"): the
## standard deviation is zero
```

```r
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
## Austria        0.2229 [ 0.1801; 0.2648]        5.1
## Belgium        0.1463 [ 0.1001; 0.1919]        5.0
## Czech Republic 0.1659 [ 0.1237; 0.2075]        5.1
## Estonia        0.0472 [ 0.0025; 0.0917]        5.1
## Finland        0.2230 [ 0.1791; 0.2660]        5.1
## France         0.1974 [ 0.1549; 0.2392]        5.1
## Germany        0.2335 [ 0.1981; 0.2683]        5.3
## Great Britain  0.1816 [ 0.1372; 0.2252]        5.1
## Hungary        0.1737 [ 0.1213; 0.2252]        4.8
## Ireland        0.1739 [ 0.1363; 0.2110]        5.3
## Italy          0.1357 [ 0.0943; 0.1767]        5.2
## Lithuania      0.0951 [ 0.0486; 0.1412]        5.0
## Netherlands    0.1479 [ 0.0999; 0.1952]        5.0
## Norway         0.2278 [ 0.1797; 0.2748]        4.9
## Poland         0.0229 [-0.0279; 0.0735]        4.9
## Portugal       0.0864 [ 0.0303; 0.1420]        4.7
## Slovenia       0.1530 [ 0.0980; 0.2071]        4.7
## Spain          0.1390 [ 0.0916; 0.1858]        5.0
## Sweden         0.2231 [ 0.1742; 0.2708]        4.9
## Switzerland    0.1760 [ 0.1260; 0.2251]        4.9
## 
## Number of studies combined: k = 20
## 
##                         COR           95%-CI     z  p-value
## Random effects model 0.1596 [0.1338; 0.1853] 11.94 < 0.0001
## Prediction interval         [0.0415; 0.2733]               
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0031 [0.0016; 0.0072]; tau = 0.0552 [0.0395; 0.0851];
##  I^2 = 84.4% [77.2%; 89.4%]; H = 2.54 [2.09; 3.07]
## 
## Test of heterogeneity:
##       Q d.f.  p-value
##  122.14   19 < 0.0001
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
##       31       28      240
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
##                                                                       COR
## AT: BZÖ                                                            0.0000
## AT: FPÖ                                                            0.2587
## BE: N-VA                                                          -0.0320
## BE: Parti Populaire                                               -0.1021
## BE: Vlaams Belang                                                  0.1874
## CH: Swiss People's Party                                           0.0096
## CZ: ODS                                                            0.1914
## CZ: Usvit                                                          0.3619
## DE: AfD                                                           -0.0016
## DE: NPD                                                            0.1633
## EE: Eesti Konservatiivne Rahvaerakond                              0.0235
## ES: Partido Popular - PP                                           0.0959
## FI: True Finns                                                     0.1420
## FR: FN (Front National)                                            0.2846
## FR: MPF (Mouvement pour la France)                                -0.2187
## FR: UMP (Union pour un Mouvement Populaire)                        0.1587
## GB: Conservative                                                   0.2133
## GB: Democratic Unionist Party (nir)                                0.1931
## GB: UK Independence Party                                          0.0790
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)  0.2420
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                       0.2087
## IT: Fratelli d'Italia                                              0.0080
## IT: Lega Nord                                                     -0.1013
## IT: Popolo delle Libertà (PdL)                                     0.2035
## NL: Party for Freedom                                             -0.0299
## NL: Reformed Political Party                                       0.4091
## NO: Progress Party (FRP)                                           0.1942
## PL: KORWIN                                                        -0.1311
## PL: Kukiz'15                                                       0.0307
## SE: Sverigedemokraterna                                            0.0849
## SI: SDS - Slovenska demokratska stranka                            0.1381
##                                                                              95%-CI
## AT: BZÖ                                                           [-0.6296; 0.6296]
## AT: FPÖ                                                           [ 0.1373; 0.3724]
## BE: N-VA                                                          [-0.1533; 0.0902]
## BE: Parti Populaire                                               [-0.9030; 0.8574]
## BE: Vlaams Belang                                                 [-0.1854; 0.5130]
## CH: Swiss People's Party                                          [-0.1531; 0.1718]
## CZ: ODS                                                           [-0.0042; 0.3729]
## CZ: Usvit                                                         [ 0.0379; 0.6171]
## DE: AfD                                                           [-0.2598; 0.2568]
## DE: NPD                                                           [-0.5197; 0.7190]
## EE: Eesti Konservatiivne Rahvaerakond                             [-0.2127; 0.2570]
## ES: Partido Popular - PP                                          [-0.0145; 0.2041]
## FI: True Finns                                                    [-0.0007; 0.2791]
## FR: FN (Front National)                                           [ 0.1110; 0.4413]
## FR: MPF (Mouvement pour la France)                                [-0.6029; 0.2478]
## FR: UMP (Union pour un Mouvement Populaire)                       [ 0.0463; 0.2672]
## GB: Conservative                                                  [ 0.1269; 0.2964]
## GB: Democratic Unionist Party (nir)                               [-0.6552; 0.8260]
## GB: UK Independence Party                                         [-0.1134; 0.2658]
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd) [ 0.1583; 0.3223]
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                      [-0.0033; 0.4028]
## IT: Fratelli d'Italia                                             [-0.3661; 0.3800]
## IT: Lega Nord                                                     [-0.3209; 0.1286]
## IT: Popolo delle Libertà (PdL)                                    [ 0.0114; 0.3812]
## NL: Party for Freedom                                             [-0.2270; 0.1695]
## NL: Reformed Political Party                                      [-0.0151; 0.7085]
## NO: Progress Party (FRP)                                          [ 0.0163; 0.3602]
## PL: KORWIN                                                        [-0.4493; 0.2167]
## PL: Kukiz'15                                                      [-0.1677; 0.2268]
## SE: Sverigedemokraterna                                           [-0.1221; 0.2847]
## SI: SDS - Slovenska demokratska stranka                           [-0.0343; 0.3024]
##                                                                   %W(random)
## AT: BZÖ                                                                  0.3
## AT: FPÖ                                                                  5.8
## BE: N-VA                                                                 5.9
## BE: Parti Populaire                                                      0.1
## BE: Vlaams Belang                                                        1.2
## CH: Swiss People's Party                                                 4.4
## CZ: ODS                                                                  3.4
## CZ: Usvit                                                                1.5
## DE: AfD                                                                  2.2
## DE: NPD                                                                  0.3
## EE: Eesti Konservatiivne Rahvaerakond                                    2.6
## ES: Partido Popular - PP                                                 6.5
## FI: True Finns                                                           5.1
## FR: FN (Front National)                                                  3.9
## FR: MPF (Mouvement pour la France)                                       0.8
## FR: UMP (Union pour un Mouvement Populaire)                              6.3
## GB: Conservative                                                         7.6
## GB: Democratic Unionist Party (nir)                                      0.2
## GB: UK Independence Party                                                3.6
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)        7.6
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                             3.1
## IT: Fratelli d'Italia                                                    1.1
## IT: Lega Nord                                                            2.8
## IT: Popolo delle Libertà (PdL)                                           3.5
## NL: Party for Freedom                                                    3.4
## NL: Reformed Political Party                                             0.9
## NO: Progress Party (FRP)                                                 3.9
## PL: KORWIN                                                               1.4
## PL: Kukiz'15                                                             3.4
## SE: Sverigedemokraterna                                                  3.2
## SI: SDS - Slovenska demokratska stranka                                  4.1
## 
## Number of studies combined: k = 31
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1293 [ 0.0855; 0.1727] 5.74 < 0.0001
## Prediction interval         [-0.0180; 0.2711]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0047 [0.0000; 0.0136]; tau = 0.0687 [0.0000; 0.1164];
##  I^2 = 35.9% [0.8%; 58.6%]; H = 1.25 [1.00; 1.55]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  46.83   30  0.0259
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
## AT: Grüne                                           0.2685 [ 0.1276; 0.3987]
## BE: Ecolo                                           0.1091 [-0.1489; 0.3533]
## BE: Groen!                                         -0.0325 [-0.2532; 0.1914]
## CH: Green Party                                     0.0892 [-0.1601; 0.3278]
## CH: Social Democratic Party                         0.3476 [ 0.1866; 0.4904]
## DE: Bündnis 90/ Die Grünen                          0.1805 [ 0.0545; 0.3007]
## EE: Erakond Eestimaa Rohelised                      0.4930 [-0.1128; 0.8316]
## FI: Green League                                    0.1873 [ 0.0351; 0.3311]
## FR: Autres mouvements écologistes                  -0.0066 [-0.4478; 0.4372]
## FR: EELV (Europe Ecologie Les Verts)               -0.0099 [-0.2632; 0.2446]
## GB: Green Party                                     0.3931 [ 0.0840; 0.6332]
## HU: LMP (Lehet Más A Politika)                     -0.0781 [-0.5258; 0.4035]
## IE: Green Party                                     0.1208 [-0.2797; 0.4854]
## IT: Movimento 5 Stelle                              0.3121 [ 0.1884; 0.4260]
## IT: Sinistra Ecologia e Libertà (SEL)               0.3445 [-0.0586; 0.6510]
## LT: Lithuanian Greens Party (LZP)                   0.3536 [-0.4676; 0.8472]
## LT: Lithuanian Peasant and Greens Union (LVZS)      0.1456 [ 0.0199; 0.2668]
## NL: Green Left                                      0.0306 [-0.2229; 0.2803]
## NL: Party for the Animals                          -0.2584 [-0.5660; 0.1123]
## NO: Green Party (MDG)                              -0.0569 [-0.3876; 0.2868]
## NO: Liberal Party (V)                               0.0484 [-0.2303; 0.3197]
## NO: Socialist Left Party (SV)                       0.1245 [-0.1336; 0.3669]
## PT: B.E. - Bloco de Esquerda                        0.1221 [-0.1338; 0.3628]
## PT: PAN - Pessoas-Animais-Natureza                 -0.0356 [-0.6222; 0.5766]
## SE: FI (Feministiskt initiativ)                     0.3334 [-0.0453; 0.6283]
## SE: Miljöpartiet de gröna                           0.4540 [ 0.2779; 0.6006]
## SE: Vänsterpartiet                                  0.4020 [ 0.1985; 0.5722]
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)  0.0579 [-0.3294; 0.4285]
##                                                    %W(random)
## AT: Grüne                                                 6.6
## BE: Ecolo                                                 3.7
## BE: Groen!                                                4.4
## CH: Green Party                                           3.9
## CH: Social Democratic Party                               5.8
## DE: Bündnis 90/ Die Grünen                                7.2
## EE: Erakond Eestimaa Rohelised                            0.8
## FI: Green League                                          6.3
## FR: Autres mouvements écologistes                         1.5
## FR: EELV (Europe Ecologie Les Verts)                      3.7
## GB: Green Party                                           2.7
## HU: LMP (Lehet Más A Politika)                            1.3
## IE: Green Party                                           1.9
## IT: Movimento 5 Stelle                                    7.1
## IT: Sinistra Ecologia e Libertà (SEL)                     1.9
## LT: Lithuanian Greens Party (LZP)                         0.5
## LT: Lithuanian Peasant and Greens Union (LVZS)            7.2
## NL: Green Left                                            3.8
## NL: Party for the Animals                                 2.2
## NO: Green Party (MDG)                                     2.4
## NO: Liberal Party (V)                                     3.3
## NO: Socialist Left Party (SV)                             3.7
## PT: B.E. - Bloco de Esquerda                              3.8
## PT: PAN - Pessoas-Animais-Natureza                        0.8
## SE: FI (Feministiskt initiativ)                           2.1
## SE: Miljöpartiet de gröna                                 4.9
## SE: Vänsterpartiet                                        4.4
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)        2.0
## 
## Number of studies combined: k = 28
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1826 [ 0.1212; 0.2426] 5.76 < 0.0001
## Prediction interval         [-0.0311; 0.3804]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0100 [0.0000; 0.0360]; tau = 0.1000 [0.0000; 0.1898];
##  I^2 = 41.3% [7.8%; 62.6%]; H = 1.31 [1.04; 1.64]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  46.01   27  0.0127
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

##### Confidence intervals for supplementary tables


```r
cor.dat<-read.csv2("associations_within_voting_groups.csv",
                   stringsAsFactors = F)

library(metafor)
library(numform)
```

```
## 
## Attaching package: 'numform'
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```r
cor.dat$observed_r.z<-transf.rtoz(cor.dat$observed_r)
cor.dat$r.SE<-1/sqrt(cor.dat$n.y-3)
```

```
## Warning in sqrt(cor.dat$n.y - 3): NaNs produced
```

```r
cor.dat$LL<-transf.ztor(cor.dat$observed_r.z+qnorm(.025)*cor.dat$r.SE)
cor.dat$UL<-transf.ztor(cor.dat$observed_r.z+qnorm(.975)*cor.dat$r.SE)
cor.dat$z<-cor.dat$observed_r.z/cor.dat$r.SE
cor.dat$p<-(1-pnorm(abs(cor.dat$z)))*2




anti.imm.r<-cbind.data.frame(group=cor.dat[cor.dat$anti.imm==1,c("voting.group")],
                            round(cor.dat[cor.dat$anti.imm==1,c("observed_r","LL","UL")],2),
                            p=round(cor.dat[cor.dat$anti.imm==1,c("p")],3))

anti.imm.r$r<-f_num(x=anti.imm.r[,"observed_r"],digits=2)


anti.imm.r$CI<-paste0("[",
                     f_num(x=anti.imm.r[,"LL"],digits=2),", ",
                     f_num(x=anti.imm.r[,"UL"],digits=2),"]")

anti.imm.r$pno0<-f_num(x=anti.imm.r[,"p"],digits=3)
anti.imm.r
```

```
##                                                                 group
## 1                                                             AT: BZÖ
## 4                                                             AT: FPÖ
## 22                                                           BE: N-VA
## 28                                                BE: Parti Populaire
## 34                                                  BE: Vlaams Belang
## 42                                       CH: Federal Democratic Union
## 51                                           CH: Swiss People's Party
## 52                                                  CH: Ticino League
## 62                                                            CZ: ODS
## 66                                                          CZ: Usvit
## 67                                                            DE: AfD
## 77                                                            DE: NPD
## 85                              EE: Eesti Konservatiivne Rahvaerakond
## 111                                          ES: Partido Popular - PP
## 129                                                    FI: True Finns
## 136                                           FR: FN (Front National)
## 138                                FR: MPF (Mouvement pour la France)
## 148                       FR: UMP (Union pour un Mouvement Populaire)
## 149                                                  GB: Conservative
## 150                               GB: Democratic Unionist Party (nir)
## 163                                         GB: UK Independence Party
## 166 HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)
## 167                      HU: Jobbik (Jobbik Magyarországért Mozgalom)
## 190                                             IT: Fratelli d'Italia
## 192                                                     IT: Lega Nord
## 199                                    IT: Popolo delle Libertà (PdL)
## 234                                             NL: Party for Freedom
## 237                                      NL: Reformed Political Party
## 252                                          NO: Progress Party (FRP)
## 258                                                        PL: KORWIN
## 259                                                      PL: Kukiz'15
## 298                                           SE: Sverigedemokraterna
## 311                           SI: SDS - Slovenska demokratska stranka
##     observed_r    LL   UL     p    r            CI  pno0
## 1         0.00 -0.63 0.63 1.000  .00   [-.63, .63] 1.000
## 4         0.26  0.14 0.37 0.000  .26    [.14, .37]  .000
## 22       -0.03 -0.15 0.09 0.608 -.03   [-.15, .09]  .608
## 28       -0.10 -0.90 0.86 0.885 -.10   [-.90, .86]  .885
## 34        0.19 -0.19 0.51 0.325  .19   [-.19, .51]  .325
## 42        0.87 -1.00 1.00 1.000  .87 [-1.00, 1.00] 1.000
## 51        0.01 -0.15 0.17 0.909  .01   [-.15, .17]  .909
## 52          NA    NA   NA    NA <NA>      [NA, NA]  <NA>
## 62        0.19  0.00 0.37 0.055  .19    [.00, .37]  .055
## 66        0.36  0.04 0.62 0.029  .36    [.04, .62]  .029
## 67        0.00 -0.26 0.26 0.990  .00   [-.26, .26]  .990
## 77        0.16 -0.52 0.72 0.663  .16   [-.52, .72]  .663
## 85        0.02 -0.21 0.26 0.848  .02   [-.21, .26]  .848
## 111       0.10 -0.01 0.20 0.089  .10   [-.01, .20]  .089
## 129       0.14  0.00 0.28 0.051  .14    [.00, .28]  .051
## 136       0.28  0.11 0.44 0.002  .28    [.11, .44]  .002
## 138      -0.22 -0.60 0.25 0.359 -.22   [-.60, .25]  .359
## 148       0.16  0.05 0.27 0.006  .16    [.05, .27]  .006
## 149       0.21  0.13 0.30 0.000  .21    [.13, .30]  .000
## 150       0.19 -0.66 0.83 0.696  .19   [-.66, .83]  .696
## 163       0.08 -0.11 0.27 0.422  .08   [-.11, .27]  .422
## 166       0.24  0.16 0.32 0.000  .24    [.16, .32]  .000
## 167       0.21  0.00 0.40 0.054  .21    [.00, .40]  .054
## 190       0.01 -0.37 0.38 0.968  .01   [-.37, .38]  .968
## 192      -0.10 -0.32 0.13 0.388 -.10   [-.32, .13]  .388
## 199       0.20  0.01 0.38 0.038  .20    [.01, .38]  .038
## 234      -0.03 -0.23 0.17 0.771 -.03   [-.23, .17]  .771
## 237       0.41 -0.02 0.71 0.058  .41   [-.02, .71]  .058
## 252       0.19  0.02 0.36 0.033  .19    [.02, .36]  .033
## 258      -0.13 -0.45 0.22 0.463 -.13   [-.45, .22]  .463
## 259       0.03 -0.17 0.23 0.763  .03   [-.17, .23]  .763
## 298       0.08 -0.12 0.28 0.422  .08   [-.12, .28]  .422
## 311       0.14 -0.03 0.30 0.116  .14   [-.03, .30]  .116
```

```r
write.csv2(anti.imm.r,"anti.imm.r.csv")

pro.env.r<-cbind.data.frame(group=cor.dat[cor.dat$pro.env==1,c("voting.group")],
                             round(cor.dat[cor.dat$pro.env==1,c("observed_r","LL","UL")],2),
                             p=round(cor.dat[cor.dat$pro.env==1,c("p")],3))

pro.env.r$r<-f_num(x=pro.env.r[,"observed_r"],digits=2)


pro.env.r$CI<-paste0("[",
                     f_num(x=pro.env.r[,"LL"],digits=2),", ",
                     f_num(x=pro.env.r[,"UL"],digits=2),"]")

pro.env.r$pno0<-f_num(x=pro.env.r[,"p"],digits=3)
pro.env.r
```

```
##                                                  group observed_r    LL   UL
## 5                                            AT: Grüne       0.27  0.13 0.40
## 19                                           BE: Ecolo       0.11 -0.15 0.35
## 20                                          BE: Groen!      -0.03 -0.25 0.19
## 44                                     CH: Green Party       0.09 -0.16 0.33
## 50                         CH: Social Democratic Party       0.35  0.19 0.49
## 68                          DE: Bündnis 90/ Die Grünen       0.18  0.05 0.30
## 88                      EE: Erakond Eestimaa Rohelised       0.49 -0.11 0.83
## 118                                   FI: Green League       0.19  0.04 0.33
## 130                  FR: Autres mouvements écologistes      -0.01 -0.45 0.44
## 134               FR: EELV (Europe Ecologie Les Verts)      -0.01 -0.26 0.24
## 153                                    GB: Green Party       0.39  0.08 0.63
## 168                     HU: LMP (Lehet Más A Politika)      -0.08 -0.53 0.40
## 180                                    IE: Green Party       0.12 -0.28 0.49
## 193                             IT: Movimento 5 Stelle       0.31  0.19 0.43
## 201                   IT: Rivoluzione Civile (Ingroia)       0.76 -1.00 1.00
## 203              IT: Sinistra Ecologia e Libertà (SEL)       0.34 -0.06 0.65
## 211                  LT: Lithuanian Greens Party (LZP)       0.35 -0.47 0.85
## 212     LT: Lithuanian Peasant and Greens Union (LVZS)       0.15  0.02 0.27
## 227                                     NL: Green Left       0.03 -0.22 0.28
## 235                          NL: Party for the Animals      -0.26 -0.57 0.11
## 245                              NO: Green Party (MDG)      -0.06 -0.39 0.29
## 247                              NO: Liberal Party (V)       0.05 -0.23 0.32
## 254                      NO: Socialist Left Party (SV)       0.12 -0.13 0.37
## 269                       PT: B.E. - Bloco de Esquerda       0.12 -0.13 0.36
## 273                         PT: MPT - Partido da Terra         NA    NA   NA
## 278                 PT: PAN - Pessoas-Animais-Natureza      -0.04 -0.62 0.58
## 286                    SE: FI (Feministiskt initiativ)       0.33 -0.05 0.63
## 289                          SE: Miljöpartiet de gröna       0.45  0.28 0.60
## 299                                 SE: Vänsterpartiet       0.40  0.20 0.57
## 314 SI: ZL - Združena levica (DSD, IDS in Stranka TRS)       0.06 -0.33 0.43
##         p    r            CI  pno0
## 5   0.000  .27    [.13, .40]  .000
## 19  0.408  .11   [-.15, .35]  .408
## 20  0.778 -.03   [-.25, .19]  .778
## 44  0.485  .09   [-.16, .33]  .485
## 50  0.000  .35    [.19, .49]  .000
## 68  0.005  .18    [.05, .30]  .005
## 88  0.105  .49   [-.11, .83]  .105
## 118 0.016  .19    [.04, .33]  .016
## 130 0.978 -.01   [-.45, .44]  .978
## 134 0.940 -.01   [-.26, .24]  .940
## 153 0.014  .39    [.08, .63]  .014
## 168 0.762 -.08   [-.53, .40]  .762
## 180 0.561  .12   [-.28, .49]  .561
## 193 0.000  .31    [.19, .43]  .000
## 201 1.000  .76 [-1.00, 1.00] 1.000
## 203 0.092  .34   [-.06, .65]  .092
## 211 0.409  .35   [-.47, .85]  .409
## 212 0.023  .15    [.02, .27]  .023
## 227 0.815  .03   [-.22, .28]  .815
## 235 0.169 -.26   [-.57, .11]  .169
## 245 0.751 -.06   [-.39, .29]  .751
## 247 0.737  .05   [-.23, .32]  .737
## 254 0.345  .12   [-.13, .37]  .345
## 269 0.350  .12   [-.13, .36]  .350
## 273    NA <NA>      [NA, NA]  <NA>
## 278 0.920 -.04   [-.62, .58]  .920
## 286 0.083  .33   [-.05, .63]  .083
## 289 0.000  .45    [.28, .60]  .000
## 299 0.000  .40    [.20, .57]  .000
## 314 0.777  .06   [-.33, .43]  .777
```

```r
write.csv2(pro.env.r,"pro.env.r.csv")
```

\newpage

## Alternative (exploratory) approach for Hypothesis 1 with Environment attitudes as dependent variable, and immigrant attitudes as independent

### Model 0: Intercepts only


```r
H1.env.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL  UL
## (Intercept)     0.05       0.07 20.08    0.78 0.443 -0.09 0.2
```

```r
(VC.H1.env.mod0<-getVC(H1.env.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2851622 0.08131749
## 2        cntry (Intercept) <NA> 0.3004455 0.09026750
## 3     Residual        <NA> <NA> 1.1844251 1.40286286
```

```r
getDEV(H1.env.mod0)
```

```
## [1] 114073.9
```

```r
#ICC

##voting group

VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.05164826
```

```r
##country

VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.0573328
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
##             npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod0    4 114082 114116 -57037   114074                         
## H1.env.mod1   20 113129 113299 -56545   113089 984.48 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod1<-getFE(H1.env.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.04       0.12
## age                                                         0.00       0.00
## gender                                                      0.02       0.01
## educ                                                        0.03       0.00
## resid                                                      -0.14       0.01
## occupClerical support workers                               0.04       0.10
## occupCraft and related trades workers                      -0.06       0.10
## occupElementary occupations                                -0.01       0.10
## occupManagers                                               0.13       0.10
## occupOther: Not in paid work                                0.16       0.10
## occupPlant and machine operators, and assemblers           -0.04       0.10
## occupProfessionals                                          0.21       0.10
## occupRetired                                                0.06       0.11
## occupService and sales workers                              0.02       0.10
## occupSkilled agricultural, forestry and fishery workers     0.05       0.11
## occupTechnicians and associate professionals                0.06       0.10
## occupUnemployed                                            -0.03       0.12
##                                                               df t.value     p
## (Intercept)                                               188.50   -0.32 0.751
## age                                                     29829.87   -8.11 0.000
## gender                                                  35659.58    1.14 0.253
## educ                                                    35342.84   13.39 0.000
## resid                                                   35731.67  -10.19 0.000
## occupClerical support workers                           35578.76    0.42 0.672
## occupCraft and related trades workers                   35588.19   -0.60 0.548
## occupElementary occupations                             35591.13   -0.11 0.915
## occupManagers                                           35579.35    1.33 0.184
## occupOther: Not in paid work                            35719.64    1.55 0.121
## occupPlant and machine operators, and assemblers        35588.53   -0.41 0.684
## occupProfessionals                                      35581.91    2.09 0.037
## occupRetired                                            35582.94    0.52 0.604
## occupService and sales workers                          35580.41    0.25 0.806
## occupSkilled agricultural, forestry and fishery workers 35595.50    0.43 0.664
## occupTechnicians and associate professionals            35577.81    0.59 0.557
## occupUnemployed                                         35586.90   -0.25 0.803
##                                                            LL    UL
## (Intercept)                                             -0.27  0.20
## age                                                      0.00  0.00
## gender                                                  -0.01  0.04
## educ                                                     0.02  0.03
## resid                                                   -0.16 -0.11
## occupClerical support workers                           -0.15  0.24
## occupCraft and related trades workers                   -0.26  0.14
## occupElementary occupations                             -0.21  0.19
## occupManagers                                           -0.06  0.33
## occupOther: Not in paid work                            -0.04  0.36
## occupPlant and machine operators, and assemblers        -0.24  0.16
## occupProfessionals                                       0.01  0.40
## occupRetired                                            -0.16  0.28
## occupService and sales workers                          -0.17  0.22
## occupSkilled agricultural, forestry and fishery workers -0.16  0.25
## occupTechnicians and associate professionals            -0.14  0.25
## occupUnemployed                                         -0.27  0.21
```

```r
(VC.H1.env.mod1<-getVC(H1.env.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2357854 0.05559476
## 2        cntry (Intercept) <NA> 0.2942470 0.08658130
## 3     Residual        <NA> <NA> 1.1693790 1.36744715
```

```r
getDEV(H1.env.mod1)
```

```
## [1] 113089.5
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
## [1] 0.02524532
```

```r
##lvl 2: voting group

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.3163247
```

```r
##lvl 3: country

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]
```

```
## [1] 0.04083635
```

```r
##total

(sum(VC.H1.env.mod0$est_SD2)-sum(VC.H1.env.mod1$est_SD2))/
  sum(VC.H1.env.mod0$est_SD2)
```

```
## [1] 0.04117294
```

```r
#individual contributions of covariates
anova(H1.env.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age     89.975  89.975     1 29830  65.7977 5.187e-16 ***
## gender   1.791   1.791     1 35660   1.3094    0.2525    
## educ   245.226 245.226     1 35343 179.3314 < 2.2e-16 ***
## resid  142.111 142.111     1 35732 103.9243 < 2.2e-16 ***
## occup  205.203  17.100    12 34368  12.5052 < 2.2e-16 ***
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
##             npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod1   20 113129 113299 -56545   113089                         
## H1.env.mod2   21 112422 112600 -56190   112380 709.58  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod2<-getFE(H1.env.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.04       0.12
## age                                                         0.00       0.00
## gender                                                      0.01       0.01
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.01
## occupClerical support workers                               0.05       0.10
## occupCraft and related trades workers                      -0.05       0.10
## occupElementary occupations                                -0.01       0.10
## occupManagers                                               0.13       0.10
## occupOther: Not in paid work                                0.14       0.10
## occupPlant and machine operators, and assemblers           -0.03       0.10
## occupProfessionals                                          0.19       0.10
## occupRetired                                                0.06       0.11
## occupService and sales workers                              0.03       0.10
## occupSkilled agricultural, forestry and fishery workers     0.05       0.10
## occupTechnicians and associate professionals                0.06       0.10
## occupUnemployed                                            -0.02       0.12
## refugees.lvl1                                               0.16       0.01
##                                                               df t.value     p
## (Intercept)                                               182.49   -0.30 0.763
## age                                                     30300.28   -8.71 0.000
## gender                                                  35654.32    0.49 0.622
## educ                                                    35392.50   12.22 0.000
## resid                                                   35728.36   -9.40 0.000
## occupClerical support workers                           35574.16    0.48 0.629
## occupCraft and related trades workers                   35583.40   -0.49 0.625
## occupElementary occupations                             35586.20   -0.12 0.903
## occupManagers                                           35574.73    1.34 0.181
## occupOther: Not in paid work                            35722.40    1.35 0.178
## occupPlant and machine operators, and assemblers        35583.69   -0.30 0.762
## occupProfessionals                                      35577.21    1.96 0.050
## occupRetired                                            35578.52    0.56 0.573
## occupService and sales workers                          35575.86    0.31 0.760
## occupSkilled agricultural, forestry and fishery workers 35590.53    0.51 0.610
## occupTechnicians and associate professionals            35573.19    0.63 0.527
## occupUnemployed                                         35582.73   -0.20 0.838
## refugees.lvl1                                           35453.62   26.78 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.20
## age                                                      0.00  0.00
## gender                                                  -0.02  0.03
## educ                                                     0.02  0.03
## resid                                                   -0.15 -0.10
## occupClerical support workers                           -0.15  0.24
## occupCraft and related trades workers                   -0.24  0.15
## occupElementary occupations                             -0.21  0.18
## occupManagers                                           -0.06  0.33
## occupOther: Not in paid work                            -0.06  0.34
## occupPlant and machine operators, and assemblers        -0.23  0.17
## occupProfessionals                                       0.00  0.38
## occupRetired                                            -0.15  0.28
## occupService and sales workers                          -0.16  0.22
## occupSkilled agricultural, forestry and fishery workers -0.15  0.26
## occupTechnicians and associate professionals            -0.13  0.25
## occupUnemployed                                         -0.26  0.21
## refugees.lvl1                                            0.15  0.17
```

```r
(VC.H1.env.mod2<-getVC(H1.env.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2407043 0.05793855
## 2        cntry (Intercept) <NA> 0.2946595 0.08682425
## 3     Residual        <NA> <NA> 1.1576308 1.34010908
```

```r
getDEV(H1.env.mod2)
```

```
## [1] 112379.9
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
## [1] 0.01999204
```

```r
##lvl 2: voting group

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]
```

```
## [1] -0.04215841
```

```r
##lvl 3: country

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="cntry","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]
```

```
## [1] -0.00280595
```

```r
##total

(sum(VC.H1.env.mod1$est_SD2)-sum(VC.H1.env.mod2$est_SD2))/
  sum(VC.H1.env.mod1$est_SD2)
```

```
## [1] 0.0163957
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
##             npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod2   21 112422 112600 -56190   112380                         
## H1.env.mod3   25 112389 112602 -56170   112339 40.414  4  3.554e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod3<-getFE(H1.env.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.04       0.12
## age                                                         0.00       0.00
## gender                                                      0.01       0.01
## educ                                                        0.03       0.00
## resid                                                      -0.12       0.01
## occupClerical support workers                               0.05       0.10
## occupCraft and related trades workers                      -0.04       0.10
## occupElementary occupations                                -0.01       0.10
## occupManagers                                               0.14       0.10
## occupOther: Not in paid work                                0.14       0.10
## occupPlant and machine operators, and assemblers           -0.03       0.10
## occupProfessionals                                          0.19       0.10
## occupRetired                                                0.06       0.11
## occupService and sales workers                              0.03       0.10
## occupSkilled agricultural, forestry and fishery workers     0.05       0.10
## occupTechnicians and associate professionals                0.07       0.10
## occupUnemployed                                            -0.02       0.12
## refugees.lvl1                                               0.16       0.01
##                                                               df t.value     p
## (Intercept)                                               181.41   -0.33 0.740
## age                                                     30229.24   -8.64 0.000
## gender                                                  35647.27    0.51 0.613
## educ                                                    35243.03   12.16 0.000
## resid                                                   35688.34   -9.42 0.000
## occupClerical support workers                           35558.96    0.52 0.604
## occupCraft and related trades workers                   35566.00   -0.45 0.652
## occupElementary occupations                             35571.70   -0.10 0.921
## occupManagers                                           35560.44    1.37 0.170
## occupOther: Not in paid work                            35708.21    1.40 0.163
## occupPlant and machine operators, and assemblers        35566.12   -0.27 0.788
## occupProfessionals                                      35558.85    1.99 0.047
## occupRetired                                            35554.14    0.58 0.560
## occupService and sales workers                          35560.22    0.34 0.735
## occupSkilled agricultural, forestry and fishery workers 35572.65    0.50 0.615
## occupTechnicians and associate professionals            35559.38    0.69 0.489
## occupUnemployed                                         35571.95   -0.19 0.848
## refugees.lvl1                                              18.18   13.89 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.19
## age                                                      0.00  0.00
## gender                                                  -0.02  0.03
## educ                                                     0.02  0.03
## resid                                                   -0.15 -0.10
## occupClerical support workers                           -0.14  0.25
## occupCraft and related trades workers                   -0.24  0.15
## occupElementary occupations                             -0.20  0.19
## occupManagers                                           -0.06  0.33
## occupOther: Not in paid work                            -0.06  0.34
## occupPlant and machine operators, and assemblers        -0.22  0.17
## occupProfessionals                                       0.00  0.39
## occupRetired                                            -0.15  0.28
## occupService and sales workers                          -0.16  0.23
## occupSkilled agricultural, forestry and fishery workers -0.15  0.26
## occupTechnicians and associate professionals            -0.12  0.26
## occupUnemployed                                         -0.26  0.21
## refugees.lvl1                                            0.13  0.18
```

```r
(VC.H1.env.mod3<-getVC(H1.env.mod3))
```

```
##            grp          var1          var2      est_SD       est_SD2
## 1 voting.group   (Intercept)          <NA>  0.24114570  0.0581512497
## 2 voting.group refugees.lvl1          <NA>  0.05466987  0.0029887952
## 3 voting.group   (Intercept) refugees.lvl1 -0.05527857 -0.0007287598
## 4        cntry   (Intercept)          <NA>  0.29511810  0.0870946937
## 5        cntry refugees.lvl1          <NA>  0.03814282  0.0014548745
## 6        cntry   (Intercept) refugees.lvl1  0.80966258  0.0091140767
## 7     Residual          <NA>          <NA>  1.15572880  1.3357090675
```

```r
getDEV(H1.env.mod3)
```

```
## [1] 112339.5
```

```r
write.csv2(FE.H1.env.mod3,"FE.H1.env.mod3.csv")
```

