---
title: "Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes"
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

# Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes

### Model 0: Intercepts only


```r
H1.mod0<-lmer(refugees~(1|voting.group)+(1|cntry),
              data=dat,REML=F)

(FE.H1.mod0<-getFE(H1.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.06       0.08 19.78    0.66 0.514 -0.12 0.23
```

```r
(VC.H1.mod0<-getVC(H1.mod0))
```

```
##            grp        var1 var2    est_SD   est_SD2
## 1 voting.group (Intercept) <NA> 0.2800427 0.0784239
## 2        cntry (Intercept) <NA> 0.3622383 0.1312166
## 3     Residual        <NA> <NA> 0.7761931 0.6024757
```

```r
getDEV(H1.mod0)
```

```
## [1] 86726.56
```

```r
#ICC

##voting group

VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.09656733
```

```r
##country

VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.mod0[,"est_SD2"])
```

```
## [1] 0.1615736
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
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod0    4 86735 86769 -43363    86727                         
## H1.mod1   20 85726 85896 -42843    85686 1040.8 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod1<-getFE(H1.mod1))
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
(VC.H1.mod1<-getVC(H1.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2519575 0.06348257
## 2        cntry (Intercept) <NA> 0.3579166 0.12810428
## 3     Residual        <NA> <NA> 0.7657696 0.58640315
```

```r
getDEV(H1.mod1)
```

```
## [1] 85685.75
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
## [1] 0.02667757
```

```r
##lvl 2: voting group

(VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.19052
```

```r
##lvl 3: country

(VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]-
     VC.H1.mod1[VC.H1.mod1$grp=="cntry","est_SD2"])/
  VC.H1.mod0[VC.H1.mod0$grp=="cntry","est_SD2"]
```

```
## [1] 0.0237186
```

```r
##total

(sum(VC.H1.mod0$est_SD2)-sum(VC.H1.mod1$est_SD2))/
  sum(VC.H1.mod0$est_SD2)
```

```
## [1] 0.04202131
```

```r
#individual contributions of covariates
anova(H1.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age     10.887  10.887     1 35883  18.566 1.646e-05 ***
## gender  39.781  39.781     1 36713  67.839 < 2.2e-16 ***
## educ   140.436 140.436     1 36866 239.487 < 2.2e-16 ***
## resid   26.420  26.420     1 36767  45.055 1.944e-11 ***
## occup  108.350   9.029    12 36504  15.398 < 2.2e-16 ***
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
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod1   20 85726 85896 -42843    85686                         
## H1.mod2   21 84891 85070 -42425    84849 836.28  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod2<-getFE(H1.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.01       0.10
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                              -0.01       0.06
## occupCraft and related trades workers                      -0.05       0.06
## occupElementary occupations                                 0.03       0.06
## occupManagers                                               0.05       0.06
## occupOther: Not in paid work                                0.13       0.07
## occupPlant and machine operators, and assemblers           -0.01       0.06
## occupProfessionals                                          0.11       0.06
## occupRetired                                                0.04       0.07
## occupService and sales workers                             -0.01       0.06
## occupSkilled agricultural, forestry and fishery workers    -0.02       0.07
## occupTechnicians and associate professionals                0.02       0.06
## occupUnemployed                                            -0.08       0.08
## environ.lvl1                                                0.15       0.01
##                                                               df t.value     p
## (Intercept)                                                49.75    0.14 0.889
## age                                                     36006.74   -3.24 0.001
## gender                                                  36708.56    7.50 0.000
## educ                                                    36867.54   13.29 0.000
## resid                                                   36761.15   -5.73 0.000
## occupClerical support workers                           36649.62   -0.17 0.865
## occupCraft and related trades workers                   36655.10   -0.71 0.476
## occupElementary occupations                             36656.44    0.41 0.683
## occupManagers                                           36650.52    0.77 0.444
## occupOther: Not in paid work                            36794.59    2.03 0.042
## occupPlant and machine operators, and assemblers        36655.17   -0.22 0.823
## occupProfessionals                                      36651.30    1.74 0.081
## occupRetired                                            36656.02    0.52 0.600
## occupService and sales workers                          36651.26   -0.20 0.843
## occupSkilled agricultural, forestry and fishery workers 36659.89   -0.33 0.741
## occupTechnicians and associate professionals            36647.66    0.32 0.747
## occupUnemployed                                         36665.95   -0.98 0.328
## environ.lvl1                                            36591.72   29.09 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.22
## age                                                      0.00  0.00
## gender                                                   0.05  0.08
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.14  0.12
## occupCraft and related trades workers                   -0.17  0.08
## occupElementary occupations                             -0.10  0.15
## occupManagers                                           -0.08  0.18
## occupOther: Not in paid work                             0.00  0.26
## occupPlant and machine operators, and assemblers        -0.14  0.11
## occupProfessionals                                      -0.01  0.24
## occupRetired                                            -0.10  0.18
## occupService and sales workers                          -0.14  0.11
## occupSkilled agricultural, forestry and fishery workers -0.15  0.11
## occupTechnicians and associate professionals            -0.10  0.15
## occupUnemployed                                         -0.23  0.08
## environ.lvl1                                             0.14  0.16
```

```r
(VC.H1.mod2<-getVC(H1.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2559323 0.06550132
## 2        cntry (Intercept) <NA> 0.3587319 0.12868857
## 3     Residual        <NA> <NA> 0.7569991 0.57304761
```

```r
getDEV(H1.mod2)
```

```
## [1] 84849.47
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
## [1] 0.02277536
```

```r
##total

(sum(VC.H1.mod1$est_SD2)-sum(VC.H1.mod2$est_SD2))/
  sum(VC.H1.mod1$est_SD2)
```

```
## [1] 0.01382089
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
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.mod2   21 84891 85070 -42425    84849                         
## H1.mod3   25 84764 84977 -42357    84714 135.01  4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.mod3<-getFE(H1.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.10
## age                                                         0.00       0.00
## gender                                                      0.06       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
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
## environ.lvl1                                                0.15       0.02
##                                                               df t.value     p
## (Intercept)                                                49.53    0.19 0.852
## age                                                     35762.51   -3.19 0.001
## gender                                                  36686.67    7.36 0.000
## educ                                                    36862.77   13.18 0.000
## resid                                                   36736.08   -5.81 0.000
## occupClerical support workers                           36604.39   -0.22 0.825
## occupCraft and related trades workers                   36617.42   -0.79 0.431
## occupElementary occupations                             36615.56    0.33 0.739
## occupManagers                                           36610.76    0.67 0.502
## occupOther: Not in paid work                            36752.58    1.97 0.048
## occupPlant and machine operators, and assemblers        36616.53   -0.34 0.737
## occupProfessionals                                      36613.60    1.67 0.096
## occupRetired                                            36616.24    0.43 0.666
## occupService and sales workers                          36610.26   -0.24 0.808
## occupSkilled agricultural, forestry and fishery workers 36620.36   -0.41 0.680
## occupTechnicians and associate professionals            36605.61    0.26 0.793
## occupUnemployed                                         36619.16   -1.07 0.285
## environ.lvl1                                               19.92   10.03 0.000
##                                                            LL    UL
## (Intercept)                                             -0.19  0.23
## age                                                      0.00  0.00
## gender                                                   0.05  0.08
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
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
## environ.lvl1                                             0.12  0.18
```

```r
(VC.H1.mod3<-getVC(H1.mod3))
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

```r
getDEV(H1.mod3)
```

```
## [1] 84714.46
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
## Austria        0.2338 [ 0.1917; 0.2751]        5.1
## Belgium        0.1680 [ 0.1221; 0.2131]        5.0
## Czech Republic 0.1683 [ 0.1269; 0.2090]        5.1
## Estonia        0.0075 [-0.0366; 0.0516]        5.1
## Finland        0.2746 [ 0.2321; 0.3161]        5.0
## France         0.2516 [ 0.2102; 0.2921]        5.1
## Germany        0.2491 [ 0.2141; 0.2834]        5.2
## Great Britain  0.2334 [ 0.1902; 0.2758]        5.0
## Hungary        0.1092 [ 0.0569; 0.1608]        4.9
## Ireland        0.1756 [ 0.1386; 0.2121]        5.2
## Italy          0.1620 [ 0.1221; 0.2014]        5.1
## Lithuania      0.0636 [ 0.0190; 0.1079]        5.0
## Netherlands    0.1491 [ 0.1017; 0.1958]        5.0
## Norway         0.2933 [ 0.2469; 0.3383]        4.9
## Poland         0.1586 [ 0.1103; 0.2062]        4.9
## Portugal       0.1279 [ 0.0725; 0.1825]        4.8
## Slovenia       0.0807 [ 0.0259; 0.1349]        4.8
## Spain          0.1547 [ 0.1095; 0.1993]        5.0
## Sweden         0.1989 [ 0.1502; 0.2466]        4.9
## Switzerland    0.2253 [ 0.1767; 0.2727]        4.9
## 
## Number of studies combined: k = 20
## 
##                         COR           95%-CI     z  p-value
## Random effects model 0.1755 [0.1432; 0.2074] 10.49 < 0.0001
## Prediction interval         [0.0224; 0.3205]               
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0051 [0.0027; 0.0117]; tau = 0.0717 [0.0522; 0.1081];
##  I^2 = 90.4% [86.7%; 93.1%]; H = 3.23 [2.74; 3.82]
## 
## Test of heterogeneity:
##       Q d.f.  p-value
##  198.57   19 < 0.0001
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
##       31       28      243
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
## AT: BZÖ                                                            0.2651
## AT: FPÖ                                                            0.2654
## BE: N-VA                                                           0.0009
## BE: Parti Populaire                                               -0.2924
## BE: Vlaams Belang                                                  0.2946
## CH: Swiss People's Party                                           0.1853
## CZ: ODS                                                            0.2117
## CZ: Usvit                                                          0.3232
## DE: AfD                                                            0.0888
## DE: NPD                                                            0.4192
## EE: Eesti Konservatiivne Rahvaerakond                              0.1618
## ES: Partido Popular - PP                                           0.0655
## FI: True Finns                                                     0.1688
## FR: FN (Front National)                                            0.2445
## FR: MPF (Mouvement pour la France)                                 0.1995
## FR: UMP (Union pour un Mouvement Populaire)                        0.1390
## GB: Conservative                                                   0.2580
## GB: Democratic Unionist Party (nir)                                0.2644
## GB: UK Independence Party                                         -0.0072
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)  0.1993
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                       0.0530
## IT: Fratelli d'Italia                                              0.1696
## IT: Lega Nord                                                     -0.1957
## IT: Popolo delle Libertà (PdL)                                     0.3133
## NL: Party for Freedom                                             -0.0066
## NL: Reformed Political Party                                       0.0929
## NO: Progress Party (FRP)                                           0.2785
## PL: KORWIN                                                         0.1909
## PL: Kukiz'15                                                       0.1178
## SE: Sverigedemokraterna                                            0.0212
## SI: SDS - Slovenska demokratska stranka                            0.0428
##                                                                              95%-CI
## AT: BZÖ                                                           [-0.4375; 0.7668]
## AT: FPÖ                                                           [ 0.1454; 0.3777]
## BE: N-VA                                                          [-0.1210; 0.1228]
## BE: Parti Populaire                                               [-0.9338; 0.7950]
## BE: Vlaams Belang                                                 [-0.0735; 0.5920]
## CH: Swiss People's Party                                          [ 0.0247; 0.3366]
## CZ: ODS                                                           [ 0.0179; 0.3901]
## CZ: Usvit                                                         [ 0.0039; 0.5827]
## DE: AfD                                                           [-0.1735; 0.3393]
## DE: NPD                                                           [-0.2859; 0.8298]
## EE: Eesti Konservatiivne Rahvaerakond                             [-0.0761; 0.3822]
## ES: Partido Popular - PP                                          [-0.0395; 0.1690]
## FI: True Finns                                                    [ 0.0267; 0.3042]
## FR: FN (Front National)                                           [ 0.0691; 0.4053]
## FR: MPF (Mouvement pour la France)                                [-0.2665; 0.5900]
## FR: UMP (Union pour un Mouvement Populaire)                       [ 0.0263; 0.2481]
## GB: Conservative                                                  [ 0.1734; 0.3387]
## GB: Democratic Unionist Party (nir)                               [-0.6101; 0.8485]
## GB: UK Independence Party                                         [-0.1976; 0.1838]
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd) [ 0.1152; 0.2807]
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                      [-0.1582; 0.2596]
## IT: Fratelli d'Italia                                             [-0.2100; 0.5047]
## IT: Lega Nord                                                     [-0.4033; 0.0312]
## IT: Popolo delle Libertà (PdL)                                    [ 0.1330; 0.4735]
## NL: Party for Freedom                                             [-0.2028; 0.1901]
## NL: Reformed Political Party                                      [-0.3421; 0.4951]
## NO: Progress Party (FRP)                                          [ 0.1053; 0.4354]
## PL: KORWIN                                                        [-0.1371; 0.4812]
## PL: Kukiz'15                                                      [-0.0737; 0.3009]
## SE: Sverigedemokraterna                                           [-0.1833; 0.2239]
## SI: SDS - Slovenska demokratska stranka                           [-0.1270; 0.2103]
##                                                                   %W(random)
## AT: BZÖ                                                                  0.3
## AT: FPÖ                                                                  5.7
## BE: N-VA                                                                 5.8
## BE: Parti Populaire                                                      0.1
## BE: Vlaams Belang                                                        1.2
## CH: Swiss People's Party                                                 4.4
## CZ: ODS                                                                  3.4
## CZ: Usvit                                                                1.5
## DE: AfD                                                                  2.2
## DE: NPD                                                                  0.3
## EE: Eesti Konservatiivne Rahvaerakond                                    2.6
## ES: Partido Popular - PP                                                 6.6
## FI: True Finns                                                           5.0
## FR: FN (Front National)                                                  3.8
## FR: MPF (Mouvement pour la France)                                       0.8
## FR: UMP (Union pour un Mouvement Populaire)                              6.2
## GB: Conservative                                                         7.4
## GB: Democratic Unionist Party (nir)                                      0.2
## GB: UK Independence Party                                                3.5
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)        7.5
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                             3.1
## IT: Fratelli d'Italia                                                    1.2
## IT: Lega Nord                                                            2.8
## IT: Popolo delle Libertà (PdL)                                           3.6
## NL: Party for Freedom                                                    3.4
## NL: Reformed Political Party                                             0.9
## NO: Progress Party (FRP)                                                 3.8
## PL: KORWIN                                                               1.5
## PL: Kukiz'15                                                             3.5
## SE: Sverigedemokraterna                                                  3.2
## SI: SDS - Slovenska demokratska stranka                                  4.1
## 
## Number of studies combined: k = 31
## 
##                         COR            95%-CI    z  p-value
## Random effects model 0.1462 [ 0.1025; 0.1894] 6.50 < 0.0001
## Prediction interval         [-0.0031; 0.2891]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0049 [0.0000; 0.0125]; tau = 0.0699 [0.0000; 0.1119];
##  I^2 = 37.3% [3.0%; 59.4%]; H = 1.26 [1.02; 1.57]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  47.82   30  0.0207
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
## AT: Grüne                                           0.4217 [ 0.2947; 0.5341]
## BE: Ecolo                                           0.1838 [-0.0735; 0.4182]
## BE: Groen!                                          0.0746 [-0.1505; 0.2922]
## CH: Green Party                                     0.0735 [-0.1735; 0.3119]
## CH: Social Democratic Party                         0.3328 [ 0.1704; 0.4776]
## DE: Bündnis 90/ Die Grünen                          0.1550 [ 0.0290; 0.2762]
## EE: Erakond Eestimaa Rohelised                      0.5230 [-0.0728; 0.8437]
## FI: Green League                                    0.3421 [ 0.1997; 0.4703]
## FR: Autres mouvements écologistes                  -0.3769 [-0.7023; 0.0788]
## FR: EELV (Europe Ecologie Les Verts)                0.2064 [-0.0502; 0.4374]
## GB: Green Party                                     0.2912 [-0.0314; 0.5588]
## HU: LMP (Lehet Más A Politika)                      0.0028 [-0.4646; 0.4691]
## IE: Green Party                                     0.1822 [-0.2207; 0.5320]
## IT: Movimento 5 Stelle                              0.1755 [ 0.0469; 0.2983]
## IT: Sinistra Ecologia e Libertà (SEL)               0.4277 [ 0.0569; 0.6948]
## LT: Lithuanian Greens Party (LZP)                   0.1665 [-0.5595; 0.7479]
## LT: Lithuanian Peasant and Greens Union (LVZS)      0.0921 [-0.0285; 0.2101]
## NL: Green Left                                      0.2437 [-0.0065; 0.4651]
## NL: Party for the Animals                           0.2755 [-0.0749; 0.5654]
## NO: Green Party (MDG)                               0.1243 [-0.2232; 0.4438]
## NO: Liberal Party (V)                               0.0136 [-0.2629; 0.2881]
## NO: Socialist Left Party (SV)                       0.2315 [-0.0215; 0.4567]
## PT: B.E. - Bloco de Esquerda                        0.2491 [-0.0007; 0.4697]
## PT: PAN - Pessoas-Animais-Natureza                  0.2181 [-0.4392; 0.7234]
## SE: FI (Feministiskt initiativ)                     0.0280 [-0.3487; 0.3970]
## SE: Miljöpartiet de gröna                           0.3519 [ 0.1618; 0.5168]
## SE: Vänsterpartiet                                  0.3476 [ 0.1370; 0.5281]
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS) -0.1817 [-0.5254; 0.2131]
##                                                    %W(random)
## AT: Grüne                                                 6.7
## BE: Ecolo                                                 3.7
## BE: Groen!                                                4.3
## CH: Green Party                                           3.9
## CH: Social Democratic Party                               5.8
## DE: Bündnis 90/ Die Grünen                                7.4
## EE: Erakond Eestimaa Rohelised                            0.8
## FI: Green League                                          6.4
## FR: Autres mouvements écologistes                         1.4
## FR: EELV (Europe Ecologie Les Verts)                      3.7
## GB: Green Party                                           2.6
## HU: LMP (Lehet Más A Politika)                            1.3
## IE: Green Party                                           1.8
## IT: Movimento 5 Stelle                                    7.3
## IT: Sinistra Ecologia e Libertà (SEL)                     1.9
## LT: Lithuanian Greens Party (LZP)                         0.5
## LT: Lithuanian Peasant and Greens Union (LVZS)            7.6
## NL: Green Left                                            3.7
## NL: Party for the Animals                                 2.3
## NO: Green Party (MDG)                                     2.3
## NO: Liberal Party (V)                                     3.2
## NO: Socialist Left Party (SV)                             3.7
## PT: B.E. - Bloco de Esquerda                              3.7
## PT: PAN - Pessoas-Animais-Natureza                        0.7
## SE: FI (Feministiskt initiativ)                           2.0
## SE: Miljöpartiet de gröna                                 4.9
## SE: Vänsterpartiet                                        4.4
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)        1.9
## 
## Number of studies combined: k = 28
## 
##                         COR           95%-CI    z  p-value
## Random effects model 0.2115 [0.1527; 0.2688] 6.92 < 0.0001
## Prediction interval         [0.0116; 0.3951]              
## 
## Quantifying heterogeneity:
##  tau^2 = 0.0088 [0.0000; 0.0325]; tau = 0.0938 [0.0000; 0.1802];
##  I^2 = 38.8% [3.4%; 61.2%]; H = 1.28 [1.02; 1.60]
## 
## Test of heterogeneity:
##      Q d.f. p-value
##  44.09   27  0.0202
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

\newpage

## Alternative (exploratory) approach for Hypothesis 1 with Environment attitudes as dependent variable, and immigrant attitudes as independent

### Model 0: Intercepts only


```r
H1.env.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.04       0.04 19.97     1.1 0.284 -0.04 0.12
```

```r
(VC.H1.env.mod0<-getVC(H1.env.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1796220 0.03226405
## 2        cntry (Intercept) <NA> 0.1618768 0.02620409
## 3     Residual        <NA> <NA> 0.7969943 0.63519992
```

```r
getDEV(H1.env.mod0)
```

```
## [1] 88436.47
```

```r
#ICC

##voting group

VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.04651223
```

```r
##country

VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H1.env.mod0[,"est_SD2"])
```

```
## [1] 0.03777612
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
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod0    4 88444 88479 -44218    88436                         
## H1.env.mod1   20 87387 87557 -43673    87347 1089.5 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod1<-getFE(H1.env.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.01       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers           -0.04       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.16       0.08
##                                                               df t.value     p
## (Intercept)                                               291.25    0.13 0.895
## age                                                     29940.74   -7.55 0.000
## gender                                                  36805.07    5.62 0.000
## educ                                                    36362.19   15.82 0.000
## resid                                                   36874.11   -7.18 0.000
## occupClerical support workers                           36722.23    0.53 0.599
## occupCraft and related trades workers                   36732.27   -0.63 0.526
## occupElementary occupations                             36736.18   -0.79 0.431
## occupManagers                                           36724.35    0.93 0.353
## occupOther: Not in paid work                            36855.64    0.13 0.896
## occupPlant and machine operators, and assemblers        36733.98   -0.54 0.592
## occupProfessionals                                      36725.51    1.52 0.127
## occupRetired                                            36728.84   -1.50 0.134
## occupService and sales workers                          36723.79   -0.20 0.839
## occupSkilled agricultural, forestry and fishery workers 36740.01    0.06 0.952
## occupTechnicians and associate professionals            36720.21    0.83 0.409
## occupUnemployed                                         36731.58   -2.02 0.043
##                                                            LL    UL
## (Intercept)                                             -0.14  0.16
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.05
## occupClerical support workers                           -0.10  0.17
## occupCraft and related trades workers                   -0.17  0.09
## occupElementary occupations                             -0.18  0.08
## occupManagers                                           -0.07  0.19
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.17  0.10
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.25  0.03
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.13  0.14
## occupTechnicians and associate professionals            -0.07  0.18
## occupUnemployed                                         -0.32  0.00
```

```r
(VC.H1.env.mod1<-getVC(H1.env.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1504843 0.02264552
## 2        cntry (Intercept) <NA> 0.1658533 0.02750733
## 3     Residual        <NA> <NA> 0.7859685 0.61774648
```

```r
getDEV(H1.env.mod1)
```

```
## [1] 87346.94
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
## [1] 0.02747708
```

```r
##lvl 2: voting group

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.2981192
```

```r
##lvl 3: country

(VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]-
     VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"])/
  VC.H1.env.mod0[VC.H1.env.mod0$grp=="cntry","est_SD2"]
```

```
## [1] -0.04973449
```

```r
##total

(sum(VC.H1.env.mod0$est_SD2)-sum(VC.H1.env.mod1$est_SD2))/
  sum(VC.H1.env.mod0$est_SD2)
```

```
## [1] 0.0371485
```

```r
#individual contributions of covariates
anova(H1.env.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age     35.237  35.237     1 29941  57.041 4.391e-14 ***
## gender  19.540  19.540     1 36805  31.631 1.878e-08 ***
## educ   154.651 154.651     1 36362 250.347 < 2.2e-16 ***
## resid   31.824  31.824     1 36874  51.516 7.235e-13 ***
## occup   84.050   7.004    12 35387  11.338 < 2.2e-16 ***
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
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod1   20 87387 87557 -43673    87347                         
## H1.env.mod2   21 86556 86735 -43257    86514 832.92  1  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod2<-getFE(H1.env.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.02       0.07
## age                                                         0.00       0.00
## gender                                                      0.04       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.05       0.07
## occupOther: Not in paid work                               -0.01       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.08       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.15       0.08
## refugees.lvl1                                               0.15       0.01
##                                                               df t.value     p
## (Intercept)                                               285.22    0.22 0.825
## age                                                     30671.87   -7.05 0.000
## gender                                                  36797.79    4.44 0.000
## educ                                                    36443.92   13.63 0.000
## resid                                                   36871.11   -6.28 0.000
## occupClerical support workers                           36715.79    0.55 0.585
## occupCraft and related trades workers                   36725.51   -0.52 0.602
## occupElementary occupations                             36729.24   -0.83 0.406
## occupManagers                                           36717.82    0.81 0.420
## occupOther: Not in paid work                            36860.02   -0.15 0.884
## occupPlant and machine operators, and assemblers        36727.08   -0.50 0.620
## occupProfessionals                                      36718.98    1.25 0.209
## occupRetired                                            36722.51   -1.55 0.121
## occupService and sales workers                          36717.37   -0.17 0.867
## occupSkilled agricultural, forestry and fishery workers 36733.00    0.11 0.915
## occupTechnicians and associate professionals            36713.70    0.77 0.440
## occupUnemployed                                         36725.64   -1.84 0.066
## refugees.lvl1                                           36610.66   29.03 0.000
##                                                            LL    UL
## (Intercept)                                             -0.13  0.16
## age                                                      0.00  0.00
## gender                                                   0.02  0.06
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.16  0.09
## occupElementary occupations                             -0.18  0.07
## occupManagers                                           -0.08  0.18
## occupOther: Not in paid work                            -0.14  0.12
## occupPlant and machine operators, and assemblers        -0.16  0.10
## occupProfessionals                                      -0.05  0.21
## occupRetired                                            -0.26  0.03
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.13  0.14
## occupTechnicians and associate professionals            -0.08  0.18
## occupUnemployed                                         -0.30  0.01
## refugees.lvl1                                            0.14  0.16
```

```r
(VC.H1.env.mod2<-getVC(H1.env.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1554673 0.02417008
## 2        cntry (Intercept) <NA> 0.1647756 0.02715100
## 3     Residual        <NA> <NA> 0.7769649 0.60367438
```

```r
getDEV(H1.env.mod2)
```

```
## [1] 86514.02
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
## [1] 0.02277973
```

```r
##lvl 2: voting group

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="voting.group","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="voting.group","est_SD2"]
```

```
## [1] -0.06732308
```

```r
##lvl 3: country

(VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]-
     VC.H1.env.mod2[VC.H1.env.mod2$grp=="cntry","est_SD2"])/
  VC.H1.env.mod1[VC.H1.env.mod1$grp=="cntry","est_SD2"]
```

```
## [1] 0.01295424
```

```r
##total

(sum(VC.H1.env.mod1$est_SD2)-sum(VC.H1.env.mod2$est_SD2))/
  sum(VC.H1.env.mod1$est_SD2)
```

```
## [1] 0.01932008
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
##             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H1.env.mod2   21 86556 86735 -43257    86514                         
## H1.env.mod3   25 86467 86680 -43208    86417 97.037  4  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H1.env.mod3<-getFE(H1.env.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.04       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.08       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.06       0.07
## occupUnemployed                                            -0.14       0.08
## refugees.lvl1                                               0.16       0.02
##                                                               df t.value     p
## (Intercept)                                               282.71    0.17 0.868
## age                                                     30473.94   -6.90 0.000
## gender                                                  36784.19    4.30 0.000
## educ                                                    36418.83   13.42 0.000
## resid                                                   36853.40   -6.44 0.000
## occupClerical support workers                           36684.93    0.63 0.531
## occupCraft and related trades workers                   36683.93   -0.44 0.662
## occupElementary occupations                             36695.28   -0.79 0.432
## occupManagers                                           36687.21    0.88 0.378
## occupOther: Not in paid work                            36823.47   -0.07 0.947
## occupPlant and machine operators, and assemblers        36693.16   -0.44 0.662
## occupProfessionals                                      36687.18    1.30 0.194
## occupRetired                                            36678.04   -1.53 0.127
## occupService and sales workers                          36683.91   -0.10 0.924
## occupSkilled agricultural, forestry and fishery workers 36699.67    0.15 0.878
## occupTechnicians and associate professionals            36681.44    0.86 0.391
## occupUnemployed                                         36695.97   -1.77 0.077
## refugees.lvl1                                              19.60   10.54 0.000
##                                                            LL    UL
## (Intercept)                                             -0.13  0.16
## age                                                      0.00  0.00
## gender                                                   0.02  0.05
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.04
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.16  0.10
## occupElementary occupations                             -0.18  0.08
## occupManagers                                           -0.07  0.19
## occupOther: Not in paid work                            -0.14  0.13
## occupPlant and machine operators, and assemblers        -0.16  0.10
## occupProfessionals                                      -0.04  0.21
## occupRetired                                            -0.25  0.03
## occupService and sales workers                          -0.13  0.12
## occupSkilled agricultural, forestry and fishery workers -0.13  0.15
## occupTechnicians and associate professionals            -0.07  0.18
## occupUnemployed                                         -0.30  0.02
## refugees.lvl1                                            0.13  0.19
```

```r
(VC.H1.env.mod3<-getVC(H1.env.mod3))
```

```
##            grp          var1          var2     est_SD     est_SD2
## 1 voting.group   (Intercept)          <NA> 0.15573973 0.024254865
## 2 voting.group refugees.lvl1          <NA> 0.04977941 0.002477989
## 3 voting.group   (Intercept) refugees.lvl1 0.27983922 0.002169490
## 4        cntry   (Intercept)          <NA> 0.16498211 0.027219097
## 5        cntry refugees.lvl1          <NA> 0.05981851 0.003578255
## 6        cntry   (Intercept) refugees.lvl1 0.12797080 0.001262942
## 7     Residual          <NA>          <NA> 0.77489226 0.600458009
```

```r
getDEV(H1.env.mod3)
```

```
## [1] 86416.98
```

```r
write.csv2(FE.H1.env.mod3,"FE.H1.env.mod3.csv")
```

