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

\newpage

# Hypothesis 2: Those who voted for pro-environment parties will report higher pro-refugee attitudes than those who voted for anti-immigration parties

### Model 1: random intercepts + covariates (same as in H1)


```r
H2.mod1<-lmer(refugees~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)


(FE.H2.mod1<-getFE(H2.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.01       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.10       0.08
##                                                               df t.value     p
## (Intercept)                                                50.82    0.10 0.919
## age                                                     35883.24   -4.31 0.000
## gender                                                  36713.19    8.24 0.000
## educ                                                    36866.31   15.48 0.000
## resid                                                   36766.97   -6.71 0.000
## occupClerical support workers                           36653.24   -0.10 0.917
## occupCraft and related trades workers                   36658.87   -0.81 0.417
## occupElementary occupations                             36660.30    0.28 0.783
## occupManagers                                           36654.18    0.88 0.377
## occupOther: Not in paid work                            36800.93    2.00 0.046
## occupPlant and machine operators, and assemblers        36658.99   -0.31 0.756
## occupProfessionals                                      36654.96    1.93 0.054
## occupRetired                                            36659.71    0.28 0.781
## occupService and sales workers                          36654.90   -0.24 0.812
## occupSkilled agricultural, forestry and fishery workers 36663.80   -0.33 0.739
## occupTechnicians and associate professionals            36651.27    0.42 0.671
## occupUnemployed                                         36669.21   -1.29 0.197
##                                                            LL    UL
## (Intercept)                                             -0.20  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.13  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.15
## occupManagers                                           -0.07  0.19
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
```

```r
(VC.H2.mod1<-getVC(H2.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2519575 0.06348257
## 2        cntry (Intercept) <NA> 0.3579166 0.12810428
## 3     Residual        <NA> <NA> 0.7657696 0.58640315
```

\newpage


### Model 2: Categorical predictor at level-2


```r
H2.mod2<-lmer(refugees~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup+
                all.parties.lvl2
                ,data=dat,REML=F)


(FE.H2.mod2<-getFE(H2.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.44       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.11       0.08
## all.parties.lvl2Did not vote                                0.37       0.05
## all.parties.lvl2Don't know                                  0.33       0.06
## all.parties.lvl2Invalid vote                                0.30       0.26
## all.parties.lvl2NE age                                      0.57       0.06
## all.parties.lvl2NE citizen                                  0.69       0.06
## all.parties.lvl2NE other                                    0.55       0.07
## all.parties.lvl2No answer                                   0.27       0.28
## all.parties.lvl2Other party                                 0.46       0.04
## all.parties.lvl2Pro-environment party                       0.77       0.05
##                                                               df t.value     p
## (Intercept)                                                67.99   -4.17 0.000
## age                                                     36859.99   -3.98 0.000
## gender                                                  36739.30    8.19 0.000
## educ                                                    36846.07   15.43 0.000
## resid                                                   36824.13   -6.63 0.000
## occupClerical support workers                           36712.37   -0.12 0.904
## occupCraft and related trades workers                   36719.16   -0.81 0.416
## occupElementary occupations                             36715.81    0.25 0.801
## occupManagers                                           36714.98    0.88 0.381
## occupOther: Not in paid work                            36756.02    1.86 0.063
## occupPlant and machine operators, and assemblers        36719.38   -0.32 0.746
## occupProfessionals                                      36713.86    1.91 0.056
## occupRetired                                            36712.43    0.24 0.811
## occupService and sales workers                          36713.29   -0.24 0.811
## occupSkilled agricultural, forestry and fishery workers 36727.61   -0.36 0.720
## occupTechnicians and associate professionals            36710.81    0.42 0.676
## occupUnemployed                                         36711.56   -1.35 0.177
## all.parties.lvl2Did not vote                              198.51    7.17 0.000
## all.parties.lvl2Don't know                                277.33    5.97 0.000
## all.parties.lvl2Invalid vote                              844.97    1.15 0.250
## all.parties.lvl2NE age                                    278.66   10.23 0.000
## all.parties.lvl2NE citizen                                281.66   11.41 0.000
## all.parties.lvl2NE other                                  631.13    7.34 0.000
## all.parties.lvl2No answer                                1105.01    0.97 0.333
## all.parties.lvl2Other party                               247.31   12.07 0.000
## all.parties.lvl2Pro-environment party                     271.16   14.98 0.000
##                                                            LL    UL
## (Intercept)                                             -0.65 -0.23
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## all.parties.lvl2Did not vote                             0.26  0.47
## all.parties.lvl2Don't know                               0.22  0.44
## all.parties.lvl2Invalid vote                            -0.21  0.81
## all.parties.lvl2NE age                                   0.46  0.68
## all.parties.lvl2NE citizen                               0.57  0.80
## all.parties.lvl2NE other                                 0.40  0.69
## all.parties.lvl2No answer                               -0.28  0.82
## all.parties.lvl2Other party                              0.39  0.54
## all.parties.lvl2Pro-environment party                    0.67  0.87
```

```r
(VC.H2.mod2<-getVC(H2.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1608237 0.02586425
## 2        cntry (Intercept) <NA> 0.3409973 0.11627913
## 3     Residual        <NA> <NA> 0.7657731 0.58640845
```

```r
anova(H2.mod1,H2.mod2)
```

```
## Data: dat
## Models:
## H2.mod1: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H2.mod1:     educ + resid + occup
## H2.mod2: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H2.mod2:     educ + resid + occup + all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H2.mod1   20 85726 85896 -42843    85686                         
## H2.mod2   29 85549 85796 -42745    85491 194.83  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H2.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age                9.276   9.276     1 36860  15.818 6.988e-05 ***
## gender            39.371  39.371     1 36739  67.140 2.611e-16 ***
## educ             139.699 139.699     1 36846 238.227 < 2.2e-16 ***
## resid             25.766  25.766     1 36824  43.938 3.436e-11 ***
## occup            104.542   8.712    12 36744  14.856 < 2.2e-16 ***
## all.parties.lvl2 167.017  18.557     9   347  31.646 < 2.2e-16 ***
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
## [1] 0.5925772
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
## 1  Anti-immigration party  -0.43 0.08 0.0000 0.0000     -0.59     -0.26
## 2            Did not vote  -0.06 0.09 0.4683 1.0000     -0.23      0.11
## 3              Don't know  -0.09 0.09 0.2895 1.0000     -0.27      0.08
## 4            Invalid vote  -0.13 0.27 0.6402 1.0000     -0.66      0.40
## 5                  NE age   0.14 0.09 0.1033 0.7232     -0.03      0.32
## 6              NE citizen   0.26 0.09 0.0044 0.0356      0.08      0.44
## 7                NE other   0.12 0.10 0.2320 1.0000     -0.08      0.32
## 8               No answer  -0.16 0.29 0.5862 1.0000     -0.72      0.41
## 9             Other party   0.04 0.08 0.6500 1.0000     -0.12      0.19
## 10  Pro-environment party   0.34 0.09 0.0001 0.0005      0.18      0.51
```

```r
write.csv2(H2.mod2.mmeans.tab,"H2.mod2.mmeans.tab.csv")


#contrast between anti-immigration and pro-environment
(H2.contrast<-data.frame(pairs(H2.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.7691774 0.05133488 Inf
##    z.ratio      p.value
## 1 14.98352 9.409446e-51
```

```r
#contrast for all groups against mean of other groups
contrast(H2.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.4782 0.0561 Inf -8.532  <.0001 
##  Did not vote effect            -0.0726 0.0578 Inf -1.257  0.8351 
##  Don't know effect              -0.1078 0.0621 Inf -1.736  0.4958 
##  Invalid vote effect            -0.1442 0.2612 Inf -0.552  1.0000 
##  NE age effect                   0.1552 0.0619 Inf  2.508  0.0849 
##  NE citizen effect               0.2833 0.0657 Inf  4.313  0.0001 
##  NE other effect                 0.1303 0.0792 Inf  1.645  0.5001 
##  No answer effect               -0.1779 0.2791 Inf -0.638  1.0000 
##  Other party effect              0.0354 0.0470 Inf  0.755  1.0000 
##  Pro-environment party effect    0.3764 0.0582 Inf  6.468  <.0001 
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
## 1           Other party - Anti-immigration party 0.4622940 0.03829935 Inf
## 2 Pro-environment party - Anti-immigration party 0.7691774 0.05133488 Inf
## 3            Pro-environment party - Other party 0.3068834 0.04111839 Inf
##     z.ratio      p.value
## 1 12.070544 3.022645e-33
## 2 14.983523 2.822834e-50
## 3  7.463409 8.431203e-14
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
## [1] 0.8143721
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
## [1] 0.7342258
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
## [1] 0.7753352
```

```r
(H2.effect.size<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd)
```

```
## [1] 0.9931188
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
## [1] 0.7737921
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
## [1] 0.7559617
```

```r
(H2.effect.size.env.other<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[9,2])/H2.pooled.sd.other.env)
```

```
## [1] 0.3968455
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
## [1] 0.7926008
```

```r
(H2.effect.size.imm.other<-(H2.mod2.mmeans.tab[9,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd.other.imm)
```

```
## [1] 0.5929845
```



\newpage


### Model 3: Dummy-predictors at level-2 


```r
#did not vote left as reference

H2.mod3<-lmer(refugees~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                -0.07       0.11
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.11       0.08
## other.party.dummy                                           0.10       0.04
## dont.know.dummy                                            -0.03       0.06
## invalid.vote.dummy                                         -0.06       0.26
## no.answer.dummy                                            -0.09       0.28
## not.eligible.age.dummy                                      0.21       0.06
## not.eligible.citizenship.dummy                              0.32       0.06
## not.eligible.other.dummy                                    0.18       0.08
## anti.imm.party.dummy                                       -0.37       0.05
## pro.env.party.dummy                                         0.40       0.05
##                                                               df t.value     p
## (Intercept)                                                70.49   -0.68 0.497
## age                                                     36859.99   -3.98 0.000
## gender                                                  36739.30    8.19 0.000
## educ                                                    36846.07   15.43 0.000
## resid                                                   36824.13   -6.63 0.000
## occupClerical support workers                           36712.37   -0.12 0.904
## occupCraft and related trades workers                   36719.16   -0.81 0.416
## occupElementary occupations                             36715.81    0.25 0.801
## occupManagers                                           36714.98    0.88 0.381
## occupOther: Not in paid work                            36756.02    1.86 0.063
## occupPlant and machine operators, and assemblers        36719.38   -0.32 0.746
## occupProfessionals                                      36713.86    1.91 0.056
## occupRetired                                            36712.43    0.24 0.811
## occupService and sales workers                          36713.29   -0.24 0.811
## occupSkilled agricultural, forestry and fishery workers 36727.61   -0.36 0.720
## occupTechnicians and associate professionals            36710.81    0.42 0.676
## occupUnemployed                                         36711.56   -1.35 0.177
## other.party.dummy                                         180.75    2.39 0.018
## dont.know.dummy                                           231.81   -0.55 0.583
## invalid.vote.dummy                                        830.56   -0.25 0.805
## no.answer.dummy                                          1080.15   -0.34 0.735
## not.eligible.age.dummy                                    227.67    3.59 0.000
## not.eligible.citizenship.dummy                            239.48    5.20 0.000
## not.eligible.other.dummy                                  522.59    2.40 0.017
## anti.imm.party.dummy                                      198.51   -7.17 0.000
## pro.env.party.dummy                                       221.68    7.59 0.000
##                                                            LL    UL
## (Intercept)                                             -0.28  0.14
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## other.party.dummy                                        0.02  0.18
## dont.know.dummy                                         -0.15  0.08
## invalid.vote.dummy                                      -0.58  0.45
## no.answer.dummy                                         -0.64  0.45
## not.eligible.age.dummy                                   0.09  0.32
## not.eligible.citizenship.dummy                           0.20  0.44
## not.eligible.other.dummy                                 0.03  0.33
## anti.imm.party.dummy                                    -0.47 -0.26
## pro.env.party.dummy                                      0.30  0.51
```

```r
(VC.H2.mod3<-getVC(H2.mod3))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1608236 0.02586425
## 2        cntry (Intercept) <NA> 0.3409973 0.11627915
## 3     Residual        <NA> <NA> 0.7657731 0.58640845
```

```r
#this just confirms that the dummy and categorical
#models are identical
anova(H2.mod2,H2.mod3)
```

```
## Data: dat
## Models:
## H2.mod2: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H2.mod2:     educ + resid + occup + all.parties.lvl2
## H2.mod3: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H2.mod3:     educ + resid + occup + other.party.dummy + dont.know.dummy + 
## H2.mod3:     invalid.vote.dummy + no.answer.dummy + not.eligible.age.dummy + 
## H2.mod3:     not.eligible.citizenship.dummy + not.eligible.other.dummy + 
## H2.mod3:     anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H2.mod2   29 85549 85796 -42745    85491                        
## H2.mod3   29 85549 85796 -42745    85491     0  0  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


\newpage


### Model 4: Dummy-predictors (anti-immigration and pro-environment) at level-2 allowed to vary between countries


```r
#did not vote left as reference

H2.mod4<-lmer(refugees~(1|voting.group)+
                (anti.imm.party.dummy+pro.env.party.dummy||cntry)+
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
## (Intercept)                                                -0.07       0.11
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.11       0.08
## other.party.dummy                                           0.10       0.04
## dont.know.dummy                                            -0.03       0.05
## invalid.vote.dummy                                         -0.06       0.25
## no.answer.dummy                                            -0.11       0.27
## not.eligible.age.dummy                                      0.21       0.05
## not.eligible.citizenship.dummy                              0.32       0.06
## not.eligible.other.dummy                                    0.18       0.07
## anti.imm.party.dummy                                       -0.37       0.06
## pro.env.party.dummy                                         0.40       0.06
##                                                               df t.value     p
## (Intercept)                                                67.07   -0.69 0.493
## age                                                     36852.35   -3.96 0.000
## gender                                                  36740.17    8.20 0.000
## educ                                                    36846.17   15.43 0.000
## resid                                                   36829.47   -6.66 0.000
## occupClerical support workers                           36716.18   -0.11 0.909
## occupCraft and related trades workers                   36723.05   -0.81 0.420
## occupElementary occupations                             36719.52    0.26 0.796
## occupManagers                                           36718.77    0.88 0.381
## occupOther: Not in paid work                            36762.43    1.86 0.063
## occupPlant and machine operators, and assemblers        36723.30   -0.32 0.752
## occupProfessionals                                      36717.71    1.92 0.055
## occupRetired                                            36715.61    0.24 0.808
## occupService and sales workers                          36717.05   -0.23 0.814
## occupSkilled agricultural, forestry and fishery workers 36733.33   -0.35 0.724
## occupTechnicians and associate professionals            36714.77    0.42 0.675
## occupUnemployed                                         36715.02   -1.34 0.181
## other.party.dummy                                         140.51    2.55 0.012
## dont.know.dummy                                           186.82   -0.61 0.545
## invalid.vote.dummy                                        775.44   -0.25 0.806
## no.answer.dummy                                          1028.51   -0.41 0.679
## not.eligible.age.dummy                                    183.33    3.85 0.000
## not.eligible.citizenship.dummy                            191.13    5.55 0.000
## not.eligible.other.dummy                                  449.61    2.47 0.014
## anti.imm.party.dummy                                       34.82   -6.47 0.000
## pro.env.party.dummy                                        29.86    6.28 0.000
##                                                            LL    UL
## (Intercept)                                             -0.28  0.14
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## other.party.dummy                                        0.02  0.17
## dont.know.dummy                                         -0.14  0.07
## invalid.vote.dummy                                      -0.56  0.43
## no.answer.dummy                                         -0.65  0.42
## not.eligible.age.dummy                                   0.10  0.31
## not.eligible.citizenship.dummy                           0.21  0.43
## not.eligible.other.dummy                                 0.04  0.32
## anti.imm.party.dummy                                    -0.48 -0.25
## pro.env.party.dummy                                      0.27  0.53
```

```r
(VC.H2.mod4<-getVC(H2.mod4))
```

```
##            grp                 var1 var2    est_SD    est_SD2
## 1 voting.group          (Intercept) <NA> 0.1471037 0.02163951
## 2        cntry          (Intercept) <NA> 0.3436041 0.11806377
## 3      cntry.1 anti.imm.party.dummy <NA> 0.1261482 0.01591337
## 4      cntry.2  pro.env.party.dummy <NA> 0.1588515 0.02523380
## 5     Residual                 <NA> <NA> 0.7658168 0.58647533
```

```r
anova(H2.mod3,H2.mod4)
```

```
## Data: dat
## Models:
## H2.mod3: refugees ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H2.mod3:     educ + resid + occup + other.party.dummy + dont.know.dummy + 
## H2.mod3:     invalid.vote.dummy + no.answer.dummy + not.eligible.age.dummy + 
## H2.mod3:     not.eligible.citizenship.dummy + not.eligible.other.dummy + 
## H2.mod3:     anti.imm.party.dummy + pro.env.party.dummy
## H2.mod4: refugees ~ (1 | voting.group) + (anti.imm.party.dummy + pro.env.party.dummy || 
## H2.mod4:     cntry) + age + gender + educ + resid + occup + other.party.dummy + 
## H2.mod4:     dont.know.dummy + invalid.vote.dummy + no.answer.dummy + 
## H2.mod4:     not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H2.mod4:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H2.mod3   29 85549 85796 -42745    85491                       
## H2.mod4   31 85545 85809 -42742    85483 7.6774  2    0.02152 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage

### Model 5: explained variance by the focus groups


```r
#leave the focus group dummies out

H2.mod5<-lmer(refugees~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.00       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.10       0.08
## did.not.vote.dummy                                         -0.07       0.06
## dont.know.dummy                                            -0.10       0.06
## invalid.vote.dummy                                         -0.12       0.32
## no.answer.dummy                                            -0.17       0.33
## not.eligible.age.dummy                                      0.14       0.06
## not.eligible.citizenship.dummy                              0.25       0.07
## not.eligible.other.dummy                                    0.13       0.08
##                                                               df t.value     p
## (Intercept)                                                51.90   -0.03 0.973
## age                                                     36823.68   -3.95 0.000
## gender                                                  36694.85    8.29 0.000
## educ                                                    36793.65   15.56 0.000
## resid                                                   36748.84   -6.63 0.000
## occupClerical support workers                           36659.76   -0.10 0.919
## occupCraft and related trades workers                   36664.75   -0.82 0.413
## occupElementary occupations                             36662.28    0.25 0.800
## occupManagers                                           36660.95    0.88 0.380
## occupOther: Not in paid work                            36689.96    1.88 0.061
## occupPlant and machine operators, and assemblers        36664.19   -0.32 0.750
## occupProfessionals                                      36661.88    1.93 0.054
## occupRetired                                            36662.19    0.26 0.798
## occupService and sales workers                          36660.74   -0.25 0.803
## occupSkilled agricultural, forestry and fishery workers 36670.59   -0.34 0.737
## occupTechnicians and associate professionals            36658.11    0.43 0.670
## occupUnemployed                                         36666.33   -1.33 0.183
## did.not.vote.dummy                                        213.12   -1.19 0.235
## dont.know.dummy                                           283.39   -1.63 0.104
## invalid.vote.dummy                                        531.24   -0.37 0.710
## no.answer.dummy                                           638.43   -0.52 0.604
## not.eligible.age.dummy                                    284.60    2.18 0.030
## not.eligible.citizenship.dummy                            309.61    3.59 0.000
## not.eligible.other.dummy                                  639.75    1.55 0.122
##                                                            LL    UL
## (Intercept)                                             -0.21  0.20
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.13  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## did.not.vote.dummy                                      -0.18  0.04
## dont.know.dummy                                         -0.22  0.02
## invalid.vote.dummy                                      -0.74  0.50
## no.answer.dummy                                         -0.82  0.48
## not.eligible.age.dummy                                   0.01  0.26
## not.eligible.citizenship.dummy                           0.11  0.38
## not.eligible.other.dummy                                -0.03  0.29
```

```r
(VC.H2.mod5<-getVC(H2.mod5))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.2393027 0.0572658
## 2        cntry (Intercept) <NA> 0.3580601 0.1282070
## 3     Residual        <NA> <NA> 0.7657626 0.5863923
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H2.total.eff<-(VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"]-
     VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])/
  VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.09792889
```

```r
#see how much residual variance was explained at level-2 by anti-immigrants

H2.mod6<-lmer(refugees~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.07       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.10       0.08
## did.not.vote.dummy                                         -0.15       0.04
## dont.know.dummy                                            -0.18       0.05
## invalid.vote.dummy                                         -0.18       0.27
## no.answer.dummy                                            -0.25       0.29
## not.eligible.age.dummy                                      0.06       0.05
## not.eligible.citizenship.dummy                              0.17       0.06
## not.eligible.other.dummy                                    0.04       0.07
## anti.imm.party.dummy                                       -0.51       0.04
##                                                               df t.value     p
## (Intercept)                                                55.96    0.71 0.478
## age                                                     36858.12   -4.11 0.000
## gender                                                  36726.26    8.25 0.000
## educ                                                    36838.49   15.54 0.000
## resid                                                   36803.36   -6.66 0.000
## occupClerical support workers                           36693.49   -0.10 0.917
## occupCraft and related trades workers                   36700.21   -0.79 0.428
## occupElementary occupations                             36696.97    0.27 0.788
## occupManagers                                           36695.33    0.88 0.380
## occupOther: Not in paid work                            36734.43    1.89 0.059
## occupPlant and machine operators, and assemblers        36699.94   -0.31 0.760
## occupProfessionals                                      36696.20    1.94 0.052
## occupRetired                                            36694.50    0.26 0.794
## occupService and sales workers                          36694.43   -0.23 0.819
## occupSkilled agricultural, forestry and fishery workers 36707.52   -0.34 0.734
## occupTechnicians and associate professionals            36692.01    0.44 0.663
## occupUnemployed                                         36694.90   -1.33 0.184
## did.not.vote.dummy                                        187.40   -3.24 0.001
## dont.know.dummy                                           286.25   -3.53 0.000
## invalid.vote.dummy                                        723.82   -0.66 0.512
## no.answer.dummy                                           921.03   -0.85 0.396
## not.eligible.age.dummy                                    291.06    1.15 0.253
## not.eligible.citizenship.dummy                            298.12    3.09 0.002
## not.eligible.other.dummy                                  760.86    0.56 0.578
## anti.imm.party.dummy                                      250.98  -12.30 0.000
##                                                            LL    UL
## (Intercept)                                             -0.13  0.27
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.13  0.12
## occupCraft and related trades workers                   -0.18  0.08
## occupElementary occupations                             -0.11  0.15
## occupManagers                                           -0.07  0.19
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## did.not.vote.dummy                                      -0.23 -0.06
## dont.know.dummy                                         -0.28 -0.08
## invalid.vote.dummy                                      -0.72  0.36
## no.answer.dummy                                         -0.82  0.32
## not.eligible.age.dummy                                  -0.04  0.16
## not.eligible.citizenship.dummy                           0.06  0.28
## not.eligible.other.dummy                                -0.10  0.18
## anti.imm.party.dummy                                    -0.59 -0.43
```

```r
(VC.H2.mod6<-getVC(H2.mod6))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1812905 0.03286626
## 2        cntry (Intercept) <NA> 0.3439275 0.11828613
## 3     Residual        <NA> <NA> 0.7658013 0.58645168
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod6[VC.H2.mod6$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.4260752
```

```r
#see how much residual variance was explained at level-2 by pro-environments

H2.mod7<-lmer(refugees~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                -0.06       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.11       0.08
## did.not.vote.dummy                                         -0.02       0.05
## dont.know.dummy                                            -0.05       0.06
## invalid.vote.dummy                                         -0.10       0.30
## no.answer.dummy                                            -0.11       0.31
## not.eligible.age.dummy                                      0.19       0.06
## not.eligible.citizenship.dummy                              0.30       0.06
## not.eligible.other.dummy                                    0.18       0.08
## pro.env.party.dummy                                         0.39       0.05
##                                                               df t.value     p
## (Intercept)                                                53.46   -0.54 0.590
## age                                                     36838.47   -3.82 0.000
## gender                                                  36706.84    8.25 0.000
## educ                                                    36802.72   15.47 0.000
## resid                                                   36770.48   -6.60 0.000
## occupClerical support workers                           36676.45   -0.12 0.907
## occupCraft and related trades workers                   36681.82   -0.84 0.402
## occupElementary occupations                             36679.08    0.24 0.811
## occupManagers                                           36678.31    0.87 0.382
## occupOther: Not in paid work                            36709.87    1.86 0.064
## occupPlant and machine operators, and assemblers        36681.51   -0.33 0.738
## occupProfessionals                                      36677.68    1.90 0.057
## occupRetired                                            36678.18    0.23 0.815
## occupService and sales workers                          36677.55   -0.26 0.796
## occupSkilled agricultural, forestry and fishery workers 36688.51   -0.35 0.725
## occupTechnicians and associate professionals            36674.70    0.41 0.681
## occupUnemployed                                         36680.80   -1.35 0.176
## did.not.vote.dummy                                        206.61   -0.30 0.765
## dont.know.dummy                                           289.61   -0.84 0.401
## invalid.vote.dummy                                        610.85   -0.34 0.730
## no.answer.dummy                                           751.59   -0.36 0.722
## not.eligible.age.dummy                                    291.84    3.34 0.001
## not.eligible.citizenship.dummy                            311.14    4.81 0.000
## not.eligible.other.dummy                                  697.68    2.31 0.021
## pro.env.party.dummy                                       295.52    7.79 0.000
##                                                            LL    UL
## (Intercept)                                             -0.26  0.15
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## did.not.vote.dummy                                      -0.12  0.09
## dont.know.dummy                                         -0.16  0.06
## invalid.vote.dummy                                      -0.68  0.48
## no.answer.dummy                                         -0.72  0.50
## not.eligible.age.dummy                                   0.08  0.30
## not.eligible.citizenship.dummy                           0.18  0.43
## not.eligible.other.dummy                                 0.03  0.33
## pro.env.party.dummy                                      0.29  0.49
```

```r
(VC.H2.mod7<-getVC(H2.mod7))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2129598 0.04535187
## 2        cntry (Intercept) <NA> 0.3534071 0.12489661
## 3     Residual        <NA> <NA> 0.7657395 0.58635697
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod7[VC.H2.mod7$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.2080461
```

```r
#see how much residual variance was explained at level-2 by both focus groups

H2.mod8<-lmer(refugees~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.07       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                              -0.01       0.07
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.02       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.12       0.07
## occupPlant and machine operators, and assemblers           -0.02       0.07
## occupProfessionals                                          0.12       0.06
## occupRetired                                                0.02       0.07
## occupService and sales workers                             -0.02       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.03       0.06
## occupUnemployed                                            -0.11       0.08
## did.not.vote.dummy                                         -0.10       0.04
## dont.know.dummy                                            -0.13       0.05
## invalid.vote.dummy                                         -0.16       0.26
## no.answer.dummy                                            -0.19       0.28
## not.eligible.age.dummy                                      0.11       0.05
## not.eligible.citizenship.dummy                              0.22       0.05
## not.eligible.other.dummy                                    0.09       0.07
## anti.imm.party.dummy                                       -0.46       0.04
## pro.env.party.dummy                                         0.31       0.04
##                                                               df t.value     p
## (Intercept)                                                57.16    0.25 0.805
## age                                                     36859.99   -3.98 0.000
## gender                                                  36739.30    8.19 0.000
## educ                                                    36846.07   15.43 0.000
## resid                                                   36824.13   -6.63 0.000
## occupClerical support workers                           36712.37   -0.12 0.904
## occupCraft and related trades workers                   36719.16   -0.81 0.416
## occupElementary occupations                             36715.81    0.25 0.801
## occupManagers                                           36714.98    0.88 0.381
## occupOther: Not in paid work                            36756.02    1.86 0.063
## occupPlant and machine operators, and assemblers        36719.38   -0.32 0.746
## occupProfessionals                                      36713.86    1.91 0.056
## occupRetired                                            36712.43    0.24 0.811
## occupService and sales workers                          36713.29   -0.24 0.811
## occupSkilled agricultural, forestry and fishery workers 36727.61   -0.36 0.720
## occupTechnicians and associate professionals            36710.81    0.42 0.676
## occupUnemployed                                         36711.56   -1.35 0.177
## did.not.vote.dummy                                        180.75   -2.39 0.018
## dont.know.dummy                                           295.52   -2.76 0.006
## invalid.vote.dummy                                        871.64   -0.62 0.533
## no.answer.dummy                                          1134.09   -0.69 0.489
## not.eligible.age.dummy                                    302.90    2.31 0.022
## not.eligible.citizenship.dummy                            298.81    4.32 0.000
## not.eligible.other.dummy                                  828.86    1.25 0.211
## anti.imm.party.dummy                                      247.31  -12.07 0.000
## pro.env.party.dummy                                       287.90    7.46 0.000
##                                                            LL    UL
## (Intercept)                                             -0.18  0.23
## age                                                      0.00  0.00
## gender                                                   0.05  0.09
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.18  0.07
## occupElementary occupations                             -0.11  0.14
## occupManagers                                           -0.07  0.18
## occupOther: Not in paid work                            -0.01  0.26
## occupPlant and machine operators, and assemblers        -0.15  0.11
## occupProfessionals                                       0.00  0.25
## occupRetired                                            -0.12  0.16
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.16  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.26  0.05
## did.not.vote.dummy                                      -0.18 -0.02
## dont.know.dummy                                         -0.22 -0.04
## invalid.vote.dummy                                      -0.67  0.35
## no.answer.dummy                                         -0.74  0.35
## not.eligible.age.dummy                                   0.02  0.20
## not.eligible.citizenship.dummy                           0.12  0.32
## not.eligible.other.dummy                                -0.05  0.22
## anti.imm.party.dummy                                    -0.54 -0.39
## pro.env.party.dummy                                      0.23  0.39
```

```r
(VC.H2.mod8<-getVC(H2.mod8))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1608236 0.02586424
## 2        cntry (Intercept) <NA> 0.3409975 0.11627929
## 3     Residual        <NA> <NA> 0.7657731 0.58640845
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod8[VC.H2.mod8$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5483475
```

```r
#how much variance was left at level-2
1-H2.total.eff
```

```
## [1] 0.4516525
```

