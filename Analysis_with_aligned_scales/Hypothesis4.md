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


# Hypothesis 4: The strength of the association between environment and refugee attitudes is stronger among those who voted for pro-environment or anti-immigration parties in the previous national elections

### Model 1: without interactions (only main effects, combine H1 and H2 final models in terms of predictors)


```r
H4.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## (Intercept)                                                -0.43       0.10
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
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
## environ.lvl1                                                0.15       0.02
## all.parties.lvl2Did not vote                                0.37       0.05
## all.parties.lvl2Don't know                                  0.33       0.06
## all.parties.lvl2Invalid vote                                0.30       0.26
## all.parties.lvl2NE age                                      0.56       0.06
## all.parties.lvl2NE citizen                                  0.69       0.06
## all.parties.lvl2NE other                                    0.55       0.07
## all.parties.lvl2No answer                                   0.28       0.28
## all.parties.lvl2Other party                                 0.46       0.04
## all.parties.lvl2Pro-environment party                       0.78       0.05
##                                                               df t.value     p
## (Intercept)                                                66.43   -4.10 0.000
## age                                                     36830.97   -2.90 0.004
## gender                                                  36711.50    7.30 0.000
## educ                                                    36837.41   13.16 0.000
## resid                                                   36790.03   -5.73 0.000
## occupClerical support workers                           36661.06   -0.23 0.820
## occupCraft and related trades workers                   36675.57   -0.78 0.435
## occupElementary occupations                             36668.99    0.32 0.749
## occupManagers                                           36669.47    0.68 0.498
## occupOther: Not in paid work                            36708.71    1.87 0.062
## occupPlant and machine operators, and assemblers        36674.52   -0.34 0.735
## occupProfessionals                                      36670.78    1.67 0.095
## occupRetired                                            36666.96    0.41 0.679
## occupService and sales workers                          36666.39   -0.24 0.813
## occupSkilled agricultural, forestry and fishery workers 36681.79   -0.43 0.668
## occupTechnicians and associate professionals            36663.05    0.27 0.789
## occupUnemployed                                         36659.01   -1.11 0.266
## environ.lvl1                                               19.91    9.92 0.000
## all.parties.lvl2Did not vote                              194.25    7.20 0.000
## all.parties.lvl2Don't know                                275.54    5.92 0.000
## all.parties.lvl2Invalid vote                              809.03    1.13 0.258
## all.parties.lvl2NE age                                    276.25    9.96 0.000
## all.parties.lvl2NE citizen                                279.75   11.52 0.000
## all.parties.lvl2NE other                                  619.57    7.39 0.000
## all.parties.lvl2No answer                                1052.44    1.00 0.317
## all.parties.lvl2Other party                               245.17   11.97 0.000
## all.parties.lvl2Pro-environment party                     269.47   15.06 0.000
##                                                            LL    UL
## (Intercept)                                             -0.64 -0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.08
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.11  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                            -0.01  0.25
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.07
## environ.lvl1                                             0.12  0.18
## all.parties.lvl2Did not vote                             0.27  0.47
## all.parties.lvl2Don't know                               0.22  0.44
## all.parties.lvl2Invalid vote                            -0.22  0.81
## all.parties.lvl2NE age                                   0.45  0.67
## all.parties.lvl2NE citizen                               0.57  0.81
## all.parties.lvl2NE other                                 0.40  0.70
## all.parties.lvl2No answer                               -0.27  0.82
## all.parties.lvl2Other party                              0.38  0.54
## all.parties.lvl2Pro-environment party                    0.67  0.88
```

```r
(VC.H4.mod1<-getVC(H4.mod1))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.16326230 0.026654578
## 2 voting.group environ.lvl1         <NA> 0.05990535 0.003588651
## 3 voting.group  (Intercept) environ.lvl1 0.20041429 0.001960109
## 4        cntry  (Intercept)         <NA> 0.34186247 0.116869946
## 5        cntry environ.lvl1         <NA> 0.06030115 0.003636229
## 6        cntry  (Intercept) environ.lvl1 0.32179050 0.006633615
## 7     Residual         <NA>         <NA> 0.75433188 0.569016585
```

\newpage

### Model 2: Cross-level interaction between environmental attitudes and voting group categories


```r
H4.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## H4.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H4.mod1:     age + gender + educ + resid + occup + environ.lvl1 + all.parties.lvl2
## H4.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H4.mod2:     age + gender + educ + resid + occup + environ.lvl1 + all.parties.lvl2 + 
## H4.mod2:     environ.lvl1:all.parties.lvl2
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)  
## H4.mod1   34 84589 84879 -42261    84521                      
## H4.mod2   43 84591 84957 -42252    84505 16.28  9    0.06126 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H4.mod2<-getFE(H4.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.43       0.10
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.02       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.06
## occupManagers                                               0.04       0.06
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.06
## occupProfessionals                                          0.10       0.06
## occupRetired                                                0.03       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.03       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.09       0.08
## environ.lvl1                                                0.12       0.02
## all.parties.lvl2Did not vote                                0.36       0.05
## all.parties.lvl2Don't know                                  0.34       0.06
## all.parties.lvl2Invalid vote                                0.30       0.26
## all.parties.lvl2NE age                                      0.57       0.06
## all.parties.lvl2NE citizen                                  0.69       0.06
## all.parties.lvl2NE other                                    0.55       0.07
## all.parties.lvl2No answer                                   0.28       0.28
## all.parties.lvl2Other party                                 0.47       0.04
## all.parties.lvl2Pro-environment party                       0.78       0.05
## environ.lvl1:all.parties.lvl2Did not vote                   0.00       0.03
## environ.lvl1:all.parties.lvl2Don't know                     0.02       0.04
## environ.lvl1:all.parties.lvl2Invalid vote                  -0.20       0.29
## environ.lvl1:all.parties.lvl2NE age                         0.10       0.03
## environ.lvl1:all.parties.lvl2NE citizen                    -0.01       0.04
## environ.lvl1:all.parties.lvl2NE other                       0.00       0.06
## environ.lvl1:all.parties.lvl2No answer                     -0.15       0.26
## environ.lvl1:all.parties.lvl2Other party                    0.04       0.02
## environ.lvl1:all.parties.lvl2Pro-environment party          0.05       0.03
##                                                               df t.value     p
## (Intercept)                                                66.57   -4.14 0.000
## age                                                     36829.26   -2.91 0.004
## gender                                                  36714.50    7.29 0.000
## educ                                                    36838.24   13.17 0.000
## resid                                                   36795.24   -5.72 0.000
## occupClerical support workers                           36662.17   -0.24 0.808
## occupCraft and related trades workers                   36676.56   -0.80 0.422
## occupElementary occupations                             36670.12    0.30 0.765
## occupManagers                                           36670.85    0.66 0.511
## occupOther: Not in paid work                            36709.35    1.85 0.065
## occupPlant and machine operators, and assemblers        36675.52   -0.36 0.720
## occupProfessionals                                      36670.76    1.65 0.100
## occupRetired                                            36665.97    0.38 0.704
## occupService and sales workers                          36667.84   -0.26 0.798
## occupSkilled agricultural, forestry and fishery workers 36683.83   -0.45 0.656
## occupTechnicians and associate professionals            36663.68    0.25 0.805
## occupUnemployed                                         36662.46   -1.14 0.255
## environ.lvl1                                               90.39    4.96 0.000
## all.parties.lvl2Did not vote                              200.19    7.08 0.000
## all.parties.lvl2Don't know                                276.05    5.97 0.000
## all.parties.lvl2Invalid vote                              807.05    1.13 0.258
## all.parties.lvl2NE age                                    277.19   10.18 0.000
## all.parties.lvl2NE citizen                                281.70   11.47 0.000
## all.parties.lvl2NE other                                  617.84    7.43 0.000
## all.parties.lvl2No answer                                1049.14    1.00 0.317
## all.parties.lvl2Other party                               248.17   12.10 0.000
## all.parties.lvl2Pro-environment party                     271.03   15.14 0.000
## environ.lvl1:all.parties.lvl2Did not vote                 101.65    0.14 0.887
## environ.lvl1:all.parties.lvl2Don't know                   369.13    0.40 0.688
## environ.lvl1:all.parties.lvl2Invalid vote               10196.39   -0.68 0.495
## environ.lvl1:all.parties.lvl2NE age                       252.59    2.98 0.003
## environ.lvl1:all.parties.lvl2NE citizen                   233.86   -0.33 0.744
## environ.lvl1:all.parties.lvl2NE other                    1432.78    0.05 0.956
## environ.lvl1:all.parties.lvl2No answer                   7634.48   -0.55 0.581
## environ.lvl1:all.parties.lvl2Other party                  160.34    1.89 0.061
## environ.lvl1:all.parties.lvl2Pro-environment party        242.53    1.50 0.136
##                                                            LL    UL
## (Intercept)                                             -0.64 -0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.08
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.14  0.11
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.15
## occupManagers                                           -0.08  0.17
## occupOther: Not in paid work                            -0.01  0.25
## occupPlant and machine operators, and assemblers        -0.15  0.10
## occupProfessionals                                      -0.02  0.23
## occupRetired                                            -0.11  0.17
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.10
## occupTechnicians and associate professionals            -0.11  0.14
## occupUnemployed                                         -0.24  0.06
## environ.lvl1                                             0.07  0.17
## all.parties.lvl2Did not vote                             0.26  0.47
## all.parties.lvl2Don't know                               0.22  0.45
## all.parties.lvl2Invalid vote                            -0.22  0.81
## all.parties.lvl2NE age                                   0.46  0.68
## all.parties.lvl2NE citizen                               0.57  0.81
## all.parties.lvl2NE other                                 0.41  0.70
## all.parties.lvl2No answer                               -0.27  0.82
## all.parties.lvl2Other party                              0.39  0.54
## all.parties.lvl2Pro-environment party                    0.68  0.88
## environ.lvl1:all.parties.lvl2Did not vote               -0.05  0.06
## environ.lvl1:all.parties.lvl2Don't know                 -0.06  0.09
## environ.lvl1:all.parties.lvl2Invalid vote               -0.76  0.37
## environ.lvl1:all.parties.lvl2NE age                      0.03  0.17
## environ.lvl1:all.parties.lvl2NE citizen                 -0.09  0.06
## environ.lvl1:all.parties.lvl2NE other                   -0.12  0.13
## environ.lvl1:all.parties.lvl2No answer                  -0.66  0.37
## environ.lvl1:all.parties.lvl2Other party                 0.00  0.09
## environ.lvl1:all.parties.lvl2Pro-environment party      -0.02  0.12
```

```r
(VC.H4.mod2<-getVC(H4.mod2))
```

```
##            grp         var1         var2     est_SD     est_SD2
## 1 voting.group  (Intercept)         <NA> 0.16317440 0.026625884
## 2 voting.group environ.lvl1         <NA> 0.05338228 0.002849668
## 3 voting.group  (Intercept) environ.lvl1 0.23117368 0.002013667
## 4        cntry  (Intercept)         <NA> 0.34169248 0.116753750
## 5        cntry environ.lvl1         <NA> 0.06133034 0.003761411
## 6        cntry  (Intercept) environ.lvl1 0.29900958 0.006266080
## 7     Residual         <NA>         <NA> 0.75431583 0.568992376
```

\newpage

#### Marginal effect for pro-environment and anti-immigration voters


```r
H4.mod2.trends<-emtrends(H4.mod2,specs = c("all.parties.lvl2"),var=c("environ.lvl1"))
(H4.mod2.trends.tab<-data.frame(H4.mod2.trends))
```

```
##          all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL
## 1  Anti-immigration party         0.11994501 0.02417544 Inf  0.072562016
## 2            Did not vote         0.12362704 0.02149149 Inf  0.081504503
## 3              Don't know         0.13529942 0.03526300 Inf  0.066185223
## 4            Invalid vote        -0.07732986 0.28908379 Inf -0.643923679
## 5                  NE age         0.22163595 0.03099405 Inf  0.160888730
## 6              NE citizen         0.10779812 0.03410530 Inf  0.040952961
## 7                NE other         0.12342125 0.06179884 Inf  0.002297747
## 8               No answer        -0.02514293 0.26228743 Inf -0.539216832
## 9             Other party         0.16186177 0.01669941 Inf  0.129131521
## 10  Pro-environment party         0.16977729 0.02976333 Inf  0.111442224
##    asymp.UCL
## 1  0.1673280
## 2  0.1657496
## 3  0.2044136
## 4  0.4892640
## 5  0.2823832
## 6  0.1746433
## 7  0.2445448
## 8  0.4889310
## 9  0.1945920
## 10 0.2281123
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
## 1  Anti-immigration party               0.12 0.02 0.0000 0.0000      0.07
## 2            Did not vote               0.12 0.02 0.0000 0.0000      0.08
## 3              Don't know               0.14 0.04 0.0001 0.0006      0.07
## 4            Invalid vote              -0.08 0.29 0.7891 1.0000     -0.64
## 5                  NE age               0.22 0.03 0.0000 0.0000      0.16
## 6              NE citizen               0.11 0.03 0.0016 0.0063      0.04
## 7                NE other               0.12 0.06 0.0458 0.1374      0.00
## 8               No answer              -0.03 0.26 0.9236 1.0000     -0.54
## 9             Other party               0.16 0.02 0.0000 0.0000      0.13
## 10  Pro-environment party               0.17 0.03 0.0000 0.0000      0.11
##    asymp.UCL
## 1       0.17
## 2       0.17
## 3       0.20
## 4       0.49
## 5       0.28
## 6       0.17
## 7       0.24
## 8       0.49
## 9       0.19
## 10      0.23
```

```r
write.csv2(H4.mod2.trends.tab,"H4.mod2.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H4.contrast<-data.frame(pairs(H4.mod2.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.04983227 0.03328707 Inf
##    z.ratio   p.value
## 1 1.497046 0.1343813
```

```r
#contrast for all groups against mean of other groups
contrast(H4.mod2.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect   0.0154 0.0486 Inf  0.317  0.7515 
##  Did not vote effect             0.0195 0.0474 Inf  0.412  0.6807 
##  Don't know effect               0.0325 0.0549 Inf  0.591  0.5545 
##  Invalid vote effect            -0.2038 0.2904 Inf -0.702  0.4828 
##  NE age effect                   0.1284 0.0523 Inf  2.455  0.0141 
##  NE citizen effect               0.0019 0.0541 Inf  0.035  0.9720 
##  NE other effect                 0.0193 0.0746 Inf  0.258  0.7962 
##  No answer effect               -0.1458 0.2641 Inf -0.552  0.5809 
##  Other party effect              0.0620 0.0453 Inf  1.367  0.1717 
##  Pro-environment party effect    0.0708 0.0517 Inf  1.370  0.1707 
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
## 1           Other party - Anti-immigration party 0.041916758 0.02223034 Inf
## 2 Pro-environment party - Anti-immigration party 0.049832275 0.03328707 Inf
## 3            Pro-environment party - Other party 0.007915517 0.02803816 Inf
##     z.ratio    p.value
## 1 1.8855654 0.05935355
## 2 1.4970460 0.13438128
## 3 0.2823123 0.77770405
```



\newpage
