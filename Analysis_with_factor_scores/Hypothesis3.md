---
title: "Hypothesis 3"
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

# Hypothesis 3: Those who voted for anti-immigration parties will report lower pro-environment attitudes than those who voted for pro-environment parties.

### Model 0: random intercepts for environment attitudes



```r
H3.mod0<-lmer(environ.gmc~(1|voting.group),#+(1|cntry),
              data=dat,REML=F)

(FE.H3.mod0<-getFE(H3.mod0))
```

```
##             Estimate Std..Error     df t.value    p   LL   UL
## (Intercept)     0.03       0.01 219.95    2.61 0.01 0.01 0.06
```

```r
(VC.H3.mod0<-getVC(H3.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1750141 0.03062994
## 2     Residual        <NA> <NA> 0.7409473 0.54900283
```

```r
#ICC

##voting group

VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## [1] 0.0528437
```

```r
##country

VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## numeric(0)
```

\newpage

### Model 1: random intercepts + covariates


```r
H3.mod1<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

anova(H3.mod0,H3.mod1)
```

```
## Data: dat
## Models:
## H3.mod0: environ.gmc ~ (1 | voting.group)
## H3.mod1: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod1:     occup
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H3.mod0    3 60573 60598 -30284    60567                        
## H3.mod1   19 59726 59882 -29844    59688 879.3 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H3.mod1<-getFE(H3.mod1))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers            0.01       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.09       0.08
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.07       0.07
## occupUnemployed                                            -0.12       0.08
##                                                               df t.value     p
## (Intercept)                                             24922.59   -0.13 0.894
## age                                                     22321.00   -6.87 0.000
## gender                                                  26838.20    4.84 0.000
## educ                                                    25203.66   15.41 0.000
## resid                                                   26831.48   -5.11 0.000
## occupClerical support workers                           26800.75    0.61 0.539
## occupCraft and related trades workers                   26807.69   -0.47 0.640
## occupElementary occupations                             26817.30   -0.68 0.497
## occupManagers                                           26800.42    1.14 0.256
## occupOther: Not in paid work                            26881.47    0.06 0.949
## occupPlant and machine operators, and assemblers        26813.57    0.19 0.849
## occupProfessionals                                      26794.04    1.45 0.146
## occupRetired                                            26803.99   -1.16 0.247
## occupService and sales workers                          26799.57   -0.13 0.895
## occupSkilled agricultural, forestry and fishery workers 26801.96   -0.15 0.882
## occupTechnicians and associate professionals            26794.67    1.08 0.279
## occupUnemployed                                         26820.32   -1.40 0.160
##                                                            LL    UL
## (Intercept)                                             -0.14  0.12
## age                                                      0.00  0.00
## gender                                                   0.03  0.06
## educ                                                     0.02  0.03
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.16  0.10
## occupElementary occupations                             -0.18  0.09
## occupManagers                                           -0.06  0.21
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.12  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.24  0.06
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.15  0.13
## occupTechnicians and associate professionals            -0.06  0.20
## occupUnemployed                                         -0.28  0.05
```

```r
(VC.H3.mod1<-getVC(H3.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1473832 0.02172182
## 2     Residual        <NA> <NA> 0.7296591 0.53240246
```

```r
#variance explained

##lvl 1: individuals

(VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="Residual","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]
```

```
## [1] 0.03023731
```

```r
##lvl 2: voting group

(VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.2908304
```

```r
##lvl 3: country

(VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="cntry","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]
```

```
## numeric(0)
```

```r
##total

(sum(VC.H3.mod0$est_SD2)-sum(VC.H3.mod1$est_SD2))/
  sum(VC.H3.mod0$est_SD2)
```

```
## [1] 0.04400802
```

```r
#individual contributions of covariates
anova(H3.mod1)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age     25.098  25.098     1 22321  47.1418 6.777e-12 ***
## gender  12.470  12.470     1 26838  23.4228 1.307e-06 ***
## educ   126.411 126.411     1 25204 237.4347 < 2.2e-16 ***
## resid   13.913  13.913     1 26832  26.1324 3.210e-07 ***
## occup   55.892   4.658    12 26249   8.7484 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\newpage


### Model 2: Categorical predictor at level-2


```r
H3.mod2<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
                age+gender+educ+resid+occup+
                all.parties.lvl2
                ,data=dat,REML=F)


(FE.H3.mod2<-getFE(H3.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.17       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.11       0.08
## all.parties.lvl2Did not vote                                0.03       0.04
## all.parties.lvl2Don't know                                  0.05       0.04
## all.parties.lvl2Invalid vote                               -0.23       0.27
## all.parties.lvl2NE age                                      0.16       0.04
## all.parties.lvl2NE citizen                                  0.15       0.04
## all.parties.lvl2NE other                                    0.26       0.07
## all.parties.lvl2No answer                                   0.10       0.28
## all.parties.lvl2Other party                                 0.15       0.03
## all.parties.lvl2Pro-environment party                       0.45       0.04
##                                                               df t.value     p
## (Intercept)                                              7166.05   -2.35 0.019
## age                                                     26695.16   -6.70 0.000
## gender                                                  26867.90    4.90 0.000
## educ                                                    23779.10   14.96 0.000
## resid                                                   26366.56   -4.92 0.000
## occupClerical support workers                           26866.97    0.70 0.484
## occupCraft and related trades workers                   26871.80   -0.36 0.719
## occupElementary occupations                             26876.94   -0.58 0.565
## occupManagers                                           26866.05    1.24 0.213
## occupOther: Not in paid work                            26885.59    0.14 0.889
## occupPlant and machine operators, and assemblers        26876.47    0.30 0.764
## occupProfessionals                                      26857.68    1.52 0.129
## occupRetired                                            26864.78   -1.04 0.300
## occupService and sales workers                          26864.46   -0.03 0.975
## occupSkilled agricultural, forestry and fishery workers 26868.14   -0.06 0.951
## occupTechnicians and associate professionals            26861.43    1.16 0.245
## occupUnemployed                                         26875.78   -1.35 0.176
## all.parties.lvl2Did not vote                              126.75    0.95 0.342
## all.parties.lvl2Don't know                                224.18    1.12 0.263
## all.parties.lvl2Invalid vote                             4317.35   -0.84 0.399
## all.parties.lvl2NE age                                    239.63    3.61 0.000
## all.parties.lvl2NE citizen                                208.72    3.50 0.001
## all.parties.lvl2NE other                                  773.05    3.91 0.000
## all.parties.lvl2No answer                                4318.42    0.37 0.708
## all.parties.lvl2Other party                               165.07    5.28 0.000
## all.parties.lvl2Pro-environment party                     200.99   11.89 0.000
##                                                            LL    UL
## (Intercept)                                             -0.30 -0.03
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.05  0.22
## occupOther: Not in paid work                            -0.13  0.15
## occupPlant and machine operators, and assemblers        -0.11  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.23  0.07
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.14
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.28  0.05
## all.parties.lvl2Did not vote                            -0.04  0.11
## all.parties.lvl2Don't know                              -0.04  0.13
## all.parties.lvl2Invalid vote                            -0.77  0.31
## all.parties.lvl2NE age                                   0.07  0.24
## all.parties.lvl2NE citizen                               0.07  0.24
## all.parties.lvl2NE other                                 0.13  0.40
## all.parties.lvl2No answer                               -0.44  0.64
## all.parties.lvl2Other party                              0.09  0.20
## all.parties.lvl2Pro-environment party                    0.37  0.52
```

```r
(VC.H3.mod2<-getVC(H3.mod2))
```

```
##            grp        var1 var2     est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.09154595 0.00838066
## 2     Residual        <NA> <NA> 0.72965308 0.53239361
```

```r
anova(H3.mod1,H3.mod2)
```

```
## Data: dat
## Models:
## H3.mod1: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod1:     occup
## H3.mod2: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod2:     occup + all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H3.mod1   19 59726 59882 -29844    59688                         
## H3.mod2   28 59607 59837 -29776    59551 136.78  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H3.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                   Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
## age               23.896  23.896     1 26695.2  44.883 2.133e-11 ***
## gender            12.783  12.783     1 26867.9  24.010 9.640e-07 ***
## educ             119.224 119.224     1 23779.1 223.940 < 2.2e-16 ***
## resid             12.882  12.882     1 26366.6  24.197 8.751e-07 ***
## occup             54.349   4.529    12 26777.1   8.507 2.351e-16 ***
## all.parties.lvl2  98.307  10.923     9   289.8  20.517 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H3.total.eff<-(VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"]-
     VC.H3.mod2[VC.H3.mod2$grp=="voting.group","est_SD2"])/
  VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.6141824
```

\newpage

#### Marginal means for and contrasts between Pro-environment and Anti-immigration parties


```r
H3.mod2.mmeans<-emmeans(H3.mod2,specs="all.parties.lvl2")

H3.mod2.mmeans.tab<-cbind(group=data.frame(H3.mod2.mmeans)[,1],
      data.frame(H3.mod2.mmeans)[,2:6])
H3.mod2.mmeans.tab$p<-
  2*(1-pnorm(abs(H3.mod2.mmeans.tab$emmean/
                   H3.mod2.mmeans.tab$SE)))
H3.mod2.mmeans.tab$adj.p<-
  p.adjust(H3.mod2.mmeans.tab$p,method="holm")

H3.mod2.mmeans.tab<-
  cbind(group=H3.mod2.mmeans.tab[,1],
      round(H3.mod2.mmeans.tab[,c(2,3)],2),
      round(H3.mod2.mmeans.tab[,c(7,8)],4),
      round(H3.mod2.mmeans.tab[,c(5,6)],2))
H3.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.16 0.03 0.0000 0.0000     -0.21     -0.11
## 2            Did not vote  -0.12 0.03 0.0000 0.0000     -0.17     -0.07
## 3              Don't know  -0.11 0.03 0.0018 0.0127     -0.18     -0.04
## 4            Invalid vote  -0.39 0.27 0.1558 0.7789     -0.93      0.15
## 5                  NE age   0.00 0.04 0.9926 1.0000     -0.07      0.07
## 6              NE citizen   0.00 0.04 0.9089 1.0000     -0.07      0.07
## 7                NE other   0.11 0.06 0.0873 0.5241     -0.02      0.23
## 8               No answer  -0.05 0.27 0.8439 1.0000     -0.59      0.48
## 9             Other party  -0.01 0.01 0.5590 1.0000     -0.04      0.02
## 10  Pro-environment party   0.29 0.03 0.0000 0.0000      0.23      0.34
```

```r
write.csv2(H3.mod2.mmeans.tab,"H3.mod2.mmeans.tab.csv")

#contrast between anti-immigration and pro-environment
(H3.contrast<-data.frame(pairs(H3.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.4450365 0.03743501 Inf
##    z.ratio      p.value
## 1 11.88824 1.362418e-32
```

```r
#contrast for all groups against mean of other groups
contrast(H3.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.1248 0.0512 Inf -2.437  0.1333 
##  Did not vote effect            -0.0863 0.0512 Inf -1.684  0.6448 
##  Don't know effect              -0.0714 0.0560 Inf -1.276  1.0000 
##  Invalid vote effect            -0.3822 0.2757 Inf -1.386  0.9937 
##  NE age effect                   0.0501 0.0563 Inf  0.889  1.0000 
##  NE citizen effect               0.0451 0.0566 Inf  0.797  1.0000 
##  NE other effect                 0.1693 0.0764 Inf  2.215  0.2143 
##  No answer effect               -0.0102 0.2758 Inf -0.037  1.0000 
##  Other party effect              0.0407 0.0460 Inf  0.884  1.0000 
##  Pro-environment party effect    0.3697 0.0520 Inf  7.104  <.0001 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H3.more.contrasts<-data.frame(pairs(H3.mod2.mmeans, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1           Other party - Anti-immigration party 0.1488763 0.02819157 Inf
## 2 Pro-environment party - Anti-immigration party 0.4450365 0.03743501 Inf
## 3            Pro-environment party - Other party 0.2961601 0.02974614 Inf
##     z.ratio      p.value
## 1  5.280882 1.285633e-07
## 2 11.888243 4.087254e-32
## 3  9.956255 4.736302e-23
```

\newpage

#### Effect size for the difference between Anti-immigration and Pro-environment party voters


```r
H3.anti.imm.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Anti-immigration party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(environ.gmc),
            pro.ref.sd=sd(environ.gmc),
            n=n())
H3.anti.imm.sd.dat$numerator<-
  (H3.anti.imm.sd.dat$n-1)*H3.anti.imm.sd.dat$pro.ref.sd^2

H3.anti.imm.sd<-sqrt(sum(H3.anti.imm.sd.dat$numerator)/
  ((sum(H3.anti.imm.sd.dat$n)-nrow(H3.anti.imm.sd.dat))))
H3.anti.imm.sd
```

```
## [1] 0.7521715
```

```r
H3.pro.env.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Pro-environment party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(environ.gmc),
            pro.ref.sd=sd(environ.gmc),
            n=n())
H3.pro.env.sd.dat$numerator<-
  (H3.pro.env.sd.dat$n-1)*H3.pro.env.sd.dat$pro.ref.sd^2

H3.pro.env.sd<-sqrt(sum(H3.pro.env.sd.dat$numerator)/
  ((sum(H3.pro.env.sd.dat$n)-nrow(H3.pro.env.sd.dat))))
H3.pro.env.sd
```

```
## [1] 0.6766856
```

```r
H3.pooled.sd<-sqrt(
  ((nrow(H3.anti.imm.sd.dat)-1)*H3.anti.imm.sd^2+
  (nrow(H3.pro.env.sd.dat)-1)*H3.pro.env.sd^2)/
  (nrow(H3.anti.imm.sd.dat)+
     nrow(H3.pro.env.sd.dat)-2))
H3.pooled.sd  
```

```
## [1] 0.7140275
```

```r
(H3.effect.size<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd)
```

```
## [1] 0.6302278
```

```r
H3.other.sd.dat<-dat %>%
  filter(all.parties.lvl2=="Other party") %>%
  group_by(all.parties.lvl2,cntry) %>%
  summarize(pro.ref.mean=mean(environ.gmc),
            pro.ref.sd=sd(environ.gmc),
            n=n())
H3.other.sd.dat$numerator<-
  (H3.other.sd.dat$n-1)*H3.other.sd.dat$pro.ref.sd^2

H3.other.sd<-sqrt(sum(H3.other.sd.dat$numerator)/
  ((sum(H3.other.sd.dat$n)-nrow(H3.other.sd.dat))))
H3.other.sd
```

```
## [1] 0.750495
```

```r
H3.pooled.sd.other.env<-sqrt(
  ((nrow(H3.other.sd.dat)-1)*H3.other.sd^2+
  (nrow(H3.pro.env.sd.dat)-1)*H3.pro.env.sd^2)/
  (nrow(H3.other.sd.dat)+
     nrow(H3.pro.env.sd.dat)-2))
H3.pooled.sd.other.env 
```

```
## [1] 0.7158137
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.16 0.03 0.0000 0.0000     -0.21     -0.11
## 2            Did not vote  -0.12 0.03 0.0000 0.0000     -0.17     -0.07
## 3              Don't know  -0.11 0.03 0.0018 0.0127     -0.18     -0.04
## 4            Invalid vote  -0.39 0.27 0.1558 0.7789     -0.93      0.15
## 5                  NE age   0.00 0.04 0.9926 1.0000     -0.07      0.07
## 6              NE citizen   0.00 0.04 0.9089 1.0000     -0.07      0.07
## 7                NE other   0.11 0.06 0.0873 0.5241     -0.02      0.23
## 8               No answer  -0.05 0.27 0.8439 1.0000     -0.59      0.48
## 9             Other party  -0.01 0.01 0.5590 1.0000     -0.04      0.02
## 10  Pro-environment party   0.29 0.03 0.0000 0.0000      0.23      0.34
```

```r
(H3.effect.size.env.other<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[9,2])/H3.pooled.sd.other.env)
```

```
## [1] 0.4191035
```

```r
H3.pooled.sd.other.imm<-sqrt(
  ((nrow(H3.other.sd.dat)-1)*H3.other.sd^2+
  (nrow(H3.anti.imm.sd.dat)-1)*H3.anti.imm.sd^2)/
  (nrow(H3.other.sd.dat)+
     nrow(H3.anti.imm.sd.dat)-2))
H3.pooled.sd.other.imm 
```

```
## [1] 0.7512739
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.16 0.03 0.0000 0.0000     -0.21     -0.11
## 2            Did not vote  -0.12 0.03 0.0000 0.0000     -0.17     -0.07
## 3              Don't know  -0.11 0.03 0.0018 0.0127     -0.18     -0.04
## 4            Invalid vote  -0.39 0.27 0.1558 0.7789     -0.93      0.15
## 5                  NE age   0.00 0.04 0.9926 1.0000     -0.07      0.07
## 6              NE citizen   0.00 0.04 0.9089 1.0000     -0.07      0.07
## 7                NE other   0.11 0.06 0.0873 0.5241     -0.02      0.23
## 8               No answer  -0.05 0.27 0.8439 1.0000     -0.59      0.48
## 9             Other party  -0.01 0.01 0.5590 1.0000     -0.04      0.02
## 10  Pro-environment party   0.29 0.03 0.0000 0.0000      0.23      0.34
```

```r
(H3.effect.size.imm.other<-(H3.mod2.mmeans.tab[9,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd.other.imm)
```

```
## [1] 0.1996609
```


\newpage


### Model 3: Dummy-predictors at level-2 


```r
#did not vote left as reference

H3.mod3<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
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


(FE.H3.mod3<-getFE(H3.mod3))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.13       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.11       0.08
## other.party.dummy                                           0.11       0.03
## dont.know.dummy                                             0.01       0.04
## invalid.vote.dummy                                         -0.27       0.27
## no.answer.dummy                                             0.07       0.28
## not.eligible.age.dummy                                      0.12       0.04
## not.eligible.citizenship.dummy                              0.12       0.04
## not.eligible.other.dummy                                    0.23       0.07
## anti.imm.party.dummy                                       -0.03       0.04
## pro.env.party.dummy                                         0.41       0.04
##                                                               df t.value     p
## (Intercept)                                              4845.34   -1.86 0.064
## age                                                     26695.16   -6.70 0.000
## gender                                                  26867.90    4.90 0.000
## educ                                                    23779.10   14.96 0.000
## resid                                                   26366.56   -4.92 0.000
## occupClerical support workers                           26866.97    0.70 0.484
## occupCraft and related trades workers                   26871.80   -0.36 0.719
## occupElementary occupations                             26876.94   -0.58 0.565
## occupManagers                                           26866.05    1.24 0.213
## occupOther: Not in paid work                            26885.59    0.14 0.889
## occupPlant and machine operators, and assemblers        26876.47    0.30 0.764
## occupProfessionals                                      26857.68    1.52 0.129
## occupRetired                                            26864.78   -1.04 0.300
## occupService and sales workers                          26864.46   -0.03 0.975
## occupSkilled agricultural, forestry and fishery workers 26868.14   -0.06 0.951
## occupTechnicians and associate professionals            26861.43    1.16 0.245
## occupUnemployed                                         26875.78   -1.35 0.176
## other.party.dummy                                         112.02    4.01 0.000
## dont.know.dummy                                           179.75    0.31 0.756
## invalid.vote.dummy                                       4216.54   -0.97 0.333
## no.answer.dummy                                          4220.80    0.25 0.804
## not.eligible.age.dummy                                    180.20    2.86 0.005
## not.eligible.citizenship.dummy                            168.42    2.70 0.008
## not.eligible.other.dummy                                  648.92    3.40 0.001
## anti.imm.party.dummy                                      126.75   -0.95 0.341
## pro.env.party.dummy                                       154.49   10.91 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.01
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.05  0.22
## occupOther: Not in paid work                            -0.13  0.15
## occupPlant and machine operators, and assemblers        -0.11  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.23  0.07
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.14
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.28  0.05
## other.party.dummy                                        0.06  0.17
## dont.know.dummy                                         -0.07  0.10
## invalid.vote.dummy                                      -0.81  0.27
## no.answer.dummy                                         -0.47  0.61
## not.eligible.age.dummy                                   0.04  0.21
## not.eligible.citizenship.dummy                           0.03  0.20
## not.eligible.other.dummy                                 0.10  0.36
## anti.imm.party.dummy                                    -0.11  0.04
## pro.env.party.dummy                                      0.34  0.48
```

```r
(VC.H3.mod3<-getVC(H3.mod3))
```

```
##            grp        var1 var2     est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.09154594 0.00838066
## 2     Residual        <NA> <NA> 0.72965308 0.53239361
```

```r
anova(H3.mod2,H3.mod3)
```

```
## Data: dat
## Models:
## H3.mod2: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod2:     occup + all.parties.lvl2
## H3.mod3: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod3:     occup + other.party.dummy + dont.know.dummy + invalid.vote.dummy + 
## H3.mod3:     no.answer.dummy + not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H3.mod3:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H3.mod2   28 59607 59837 -29776    59551                        
## H3.mod3   28 59607 59837 -29776    59551     0  0  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


\newpage


### Model 4: Dummy-predictors at level-2 allowed to vary between countries


```r
#did not vote left as reference

H3.mod4<-lmer(environ.gmc~(1|voting.group)+(0+anti.imm.party.dummy+pro.env.party.dummy||cntry)+
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

isSingular(H3.mod4)
```

```
## [1] TRUE
```

```r
(FE.H3.mod4<-getFE(H3.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.13       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.11       0.08
## other.party.dummy                                           0.11       0.03
## dont.know.dummy                                             0.01       0.04
## invalid.vote.dummy                                         -0.27       0.27
## no.answer.dummy                                             0.07       0.27
## not.eligible.age.dummy                                      0.12       0.04
## not.eligible.citizenship.dummy                              0.12       0.04
## not.eligible.other.dummy                                    0.23       0.07
## anti.imm.party.dummy                                       -0.04       0.04
## pro.env.party.dummy                                         0.41       0.04
##                                                               df t.value     p
## (Intercept)                                              4769.49   -1.86 0.063
## age                                                     26702.66   -6.70 0.000
## gender                                                  26867.09    4.90 0.000
## educ                                                    23725.79   14.98 0.000
## resid                                                   26371.82   -4.92 0.000
## occupClerical support workers                           26865.46    0.71 0.481
## occupCraft and related trades workers                   26870.96   -0.35 0.723
## occupElementary occupations                             26875.98   -0.57 0.568
## occupManagers                                           26864.80    1.25 0.211
## occupOther: Not in paid work                            26878.74    0.14 0.888
## occupPlant and machine operators, and assemblers        26875.74    0.31 0.758
## occupProfessionals                                      26856.23    1.52 0.128
## occupRetired                                            26862.89   -1.04 0.301
## occupService and sales workers                          26863.13   -0.03 0.978
## occupSkilled agricultural, forestry and fishery workers 26866.72   -0.06 0.954
## occupTechnicians and associate professionals            26859.34    1.17 0.243
## occupUnemployed                                         26872.50   -1.35 0.176
## other.party.dummy                                         105.98    4.05 0.000
## dont.know.dummy                                           171.58    0.31 0.755
## invalid.vote.dummy                                       4148.49   -0.97 0.332
## no.answer.dummy                                          4152.70    0.25 0.803
## not.eligible.age.dummy                                    171.98    2.88 0.004
## not.eligible.citizenship.dummy                            160.31    2.72 0.007
## not.eligible.other.dummy                                  626.12    3.41 0.001
## anti.imm.party.dummy                                       32.91   -0.94 0.354
## pro.env.party.dummy                                       146.90   11.01 0.000
##                                                            LL    UL
## (Intercept)                                             -0.27  0.01
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.05  0.22
## occupOther: Not in paid work                            -0.13  0.15
## occupPlant and machine operators, and assemblers        -0.11  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.23  0.07
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.14
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.28  0.05
## other.party.dummy                                        0.06  0.17
## dont.know.dummy                                         -0.07  0.10
## invalid.vote.dummy                                      -0.80  0.27
## no.answer.dummy                                         -0.47  0.61
## not.eligible.age.dummy                                   0.04  0.21
## not.eligible.citizenship.dummy                           0.03  0.20
## not.eligible.other.dummy                                 0.10  0.36
## anti.imm.party.dummy                                    -0.12  0.04
## pro.env.party.dummy                                      0.34  0.48
```

```r
(VC.H3.mod4<-getVC(H3.mod4))
```

```
##            grp                 var1 var2     est_SD     est_SD2
## 1 voting.group          (Intercept) <NA> 0.09022949 0.008141361
## 2        cntry anti.imm.party.dummy <NA> 0.05163445 0.002666116
## 3      cntry.1  pro.env.party.dummy <NA> 0.00000000 0.000000000
## 4     Residual                 <NA> <NA> 0.72964489 0.532381665
```

```r
anova(H3.mod3,H3.mod4)
```

```
## Data: dat
## Models:
## H3.mod3: environ.gmc ~ (1 | voting.group) + age + gender + educ + resid + 
## H3.mod3:     occup + other.party.dummy + dont.know.dummy + invalid.vote.dummy + 
## H3.mod3:     no.answer.dummy + not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H3.mod3:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
## H3.mod4: environ.gmc ~ (1 | voting.group) + (0 + anti.imm.party.dummy + 
## H3.mod4:     pro.env.party.dummy || cntry) + age + gender + educ + resid + 
## H3.mod4:     occup + other.party.dummy + dont.know.dummy + invalid.vote.dummy + 
## H3.mod4:     no.answer.dummy + not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H3.mod4:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
## H3.mod3   28 59607 59837 -29776    59551                     
## H3.mod4   30 59611 59857 -29775    59551 0.4798  2     0.7867
```


\newpage

### Model 5: explained variance by the focus groups


```r
H3.mod5<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
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

(FE.H3.mod5<-getFE(H3.mod5))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.09       0.08
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.07       0.07
## occupUnemployed                                            -0.12       0.08
## did.not.vote.dummy                                         -0.14       0.04
## dont.know.dummy                                            -0.12       0.05
## invalid.vote.dummy                                         -0.40       0.29
## no.answer.dummy                                            -0.07       0.29
## not.eligible.age.dummy                                     -0.01       0.05
## not.eligible.citizenship.dummy                             -0.01       0.05
## not.eligible.other.dummy                                    0.10       0.07
##                                                               df t.value     p
## (Intercept)                                             23357.87    0.13 0.897
## age                                                     26867.83   -6.86 0.000
## gender                                                  26821.13    4.95 0.000
## educ                                                    26110.66   15.24 0.000
## resid                                                   26840.10   -5.08 0.000
## occupClerical support workers                           26810.09    0.64 0.525
## occupCraft and related trades workers                   26815.79   -0.44 0.661
## occupElementary occupations                             26821.87   -0.66 0.512
## occupManagers                                           26809.24    1.15 0.250
## occupOther: Not in paid work                            26856.42    0.06 0.951
## occupPlant and machine operators, and assemblers        26821.14    0.22 0.824
## occupProfessionals                                      26803.02    1.47 0.140
## occupRetired                                            26808.48   -1.12 0.261
## occupService and sales workers                          26807.65   -0.11 0.913
## occupSkilled agricultural, forestry and fishery workers 26810.90   -0.13 0.900
## occupTechnicians and associate professionals            26803.99    1.10 0.269
## occupUnemployed                                         26821.77   -1.40 0.162
## did.not.vote.dummy                                        137.36   -3.50 0.001
## dont.know.dummy                                           236.83   -2.69 0.008
## invalid.vote.dummy                                       2010.55   -1.37 0.171
## no.answer.dummy                                          2011.18   -0.23 0.818
## not.eligible.age.dummy                                    254.06   -0.18 0.854
## not.eligible.citizenship.dummy                            234.18   -0.29 0.775
## not.eligible.other.dummy                                  809.88    1.35 0.178
##                                                            LL    UL
## (Intercept)                                             -0.12  0.14
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.03
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.10
## occupElementary occupations                             -0.18  0.09
## occupManagers                                           -0.06  0.21
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.12  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.24  0.06
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.15  0.13
## occupTechnicians and associate professionals            -0.06  0.21
## occupUnemployed                                         -0.28  0.05
## did.not.vote.dummy                                      -0.21 -0.06
## dont.know.dummy                                         -0.21 -0.03
## invalid.vote.dummy                                      -0.98  0.17
## no.answer.dummy                                         -0.64  0.51
## not.eligible.age.dummy                                  -0.10  0.08
## not.eligible.citizenship.dummy                          -0.11  0.08
## not.eligible.other.dummy                                -0.04  0.24
```

```r
(VC.H3.mod5<-getVC(H3.mod5))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1382869 0.01912327
## 2     Residual        <NA> <NA> 0.7296465 0.53238402
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H3.total.eff<-(VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"]-
     VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])/
  VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1196288
```

```r
#see how much residual variance was explained at level-2 by anti-immigrants

H3.mod6<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
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

(FE.H3.mod6<-getFE(H3.mod6))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                 0.04       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers    -0.01       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.12       0.08
## did.not.vote.dummy                                         -0.17       0.04
## dont.know.dummy                                            -0.15       0.04
## invalid.vote.dummy                                         -0.43       0.29
## no.answer.dummy                                            -0.10       0.29
## not.eligible.age.dummy                                     -0.04       0.04
## not.eligible.citizenship.dummy                             -0.04       0.04
## not.eligible.other.dummy                                    0.07       0.07
## anti.imm.party.dummy                                       -0.20       0.03
##                                                               df t.value     p
## (Intercept)                                             23090.59    0.56 0.576
## age                                                     26837.41   -6.92 0.000
## gender                                                  26829.96    4.92 0.000
## educ                                                    25781.82   15.12 0.000
## resid                                                   26785.65   -5.04 0.000
## occupClerical support workers                           26825.92    0.65 0.517
## occupCraft and related trades workers                   26831.83   -0.40 0.686
## occupElementary occupations                             26838.23   -0.64 0.523
## occupManagers                                           26824.89    1.16 0.244
## occupOther: Not in paid work                            26869.45    0.08 0.937
## occupPlant and machine operators, and assemblers        26837.51    0.25 0.802
## occupProfessionals                                      26817.88    1.49 0.137
## occupRetired                                            26823.83   -1.10 0.273
## occupService and sales workers                          26823.34   -0.09 0.931
## occupSkilled agricultural, forestry and fishery workers 26826.55   -0.12 0.908
## occupTechnicians and associate professionals            26819.75    1.12 0.262
## occupUnemployed                                         26837.03   -1.39 0.166
## did.not.vote.dummy                                        131.21   -4.61 0.000
## dont.know.dummy                                           240.49   -3.56 0.000
## invalid.vote.dummy                                       2458.20   -1.50 0.133
## no.answer.dummy                                          2458.88   -0.34 0.735
## not.eligible.age.dummy                                    259.99   -0.93 0.354
## not.eligible.citizenship.dummy                            232.37   -1.01 0.315
## not.eligible.other.dummy                                  856.66    0.95 0.343
## anti.imm.party.dummy                                      183.92   -5.94 0.000
##                                                            LL    UL
## (Intercept)                                             -0.09  0.17
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.03
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.18  0.09
## occupManagers                                           -0.05  0.21
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.12  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.24  0.07
## occupService and sales workers                          -0.14  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.13
## occupTechnicians and associate professionals            -0.06  0.21
## occupUnemployed                                         -0.28  0.05
## did.not.vote.dummy                                      -0.24 -0.09
## dont.know.dummy                                         -0.24 -0.07
## invalid.vote.dummy                                      -0.99  0.13
## no.answer.dummy                                         -0.66  0.47
## not.eligible.age.dummy                                  -0.13  0.05
## not.eligible.citizenship.dummy                          -0.13  0.04
## not.eligible.other.dummy                                -0.07  0.20
## anti.imm.party.dummy                                    -0.27 -0.14
```

```r
(VC.H3.mod6<-getVC(H3.mod6))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1249494 0.01561235
## 2     Residual        <NA> <NA> 0.7296225 0.53234906
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod6[VC.H3.mod6$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1835938
```

```r
#see how much residual variance was explained at level-2 by pro-environments

H3.mod7<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
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

(FE.H3.mod7<-getFE(H3.mod7))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.04       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.03       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.11       0.08
## did.not.vote.dummy                                         -0.09       0.03
## dont.know.dummy                                            -0.07       0.04
## invalid.vote.dummy                                         -0.35       0.28
## no.answer.dummy                                            -0.02       0.28
## not.eligible.age.dummy                                      0.04       0.04
## not.eligible.citizenship.dummy                              0.03       0.04
## not.eligible.other.dummy                                    0.14       0.07
## pro.env.party.dummy                                         0.32       0.03
##                                                               df t.value     p
## (Intercept)                                             24262.43   -0.63 0.528
## age                                                     26766.13   -6.65 0.000
## gender                                                  26860.83    4.94 0.000
## educ                                                    24457.92   15.10 0.000
## resid                                                   26538.11   -4.98 0.000
## occupClerical support workers                           26855.74    0.69 0.491
## occupCraft and related trades workers                   26861.12   -0.39 0.694
## occupElementary occupations                             26867.23   -0.59 0.553
## occupManagers                                           26854.73    1.23 0.218
## occupOther: Not in paid work                            26885.16    0.12 0.902
## occupPlant and machine operators, and assemblers        26866.58    0.27 0.785
## occupProfessionals                                      26846.09    1.51 0.132
## occupRetired                                            26853.22   -1.06 0.287
## occupService and sales workers                          26852.98   -0.05 0.957
## occupSkilled agricultural, forestry and fishery workers 26856.86   -0.07 0.946
## occupTechnicians and associate professionals            26849.64    1.15 0.251
## occupUnemployed                                         26865.65   -1.36 0.173
## did.not.vote.dummy                                        116.01   -2.90 0.005
## dont.know.dummy                                           252.53   -1.97 0.050
## invalid.vote.dummy                                       3758.07   -1.28 0.202
## no.answer.dummy                                          3758.85   -0.07 0.943
## not.eligible.age.dummy                                    277.58    0.95 0.345
## not.eligible.citizenship.dummy                            231.20    0.80 0.425
## not.eligible.other.dummy                                  997.74    2.19 0.029
## pro.env.party.dummy                                       229.54   10.41 0.000
##                                                            LL    UL
## (Intercept)                                             -0.17  0.09
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.07 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.05  0.22
## occupOther: Not in paid work                            -0.13  0.15
## occupPlant and machine operators, and assemblers        -0.12  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.23  0.07
## occupService and sales workers                          -0.14  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.14
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.28  0.05
## did.not.vote.dummy                                      -0.15 -0.03
## dont.know.dummy                                         -0.15  0.00
## invalid.vote.dummy                                      -0.90  0.19
## no.answer.dummy                                         -0.56  0.52
## not.eligible.age.dummy                                  -0.04  0.11
## not.eligible.citizenship.dummy                          -0.05  0.11
## not.eligible.other.dummy                                 0.01  0.27
## pro.env.party.dummy                                      0.26  0.38
```

```r
(VC.H3.mod7<-getVC(H3.mod7))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1011237 0.01022599
## 2     Residual        <NA> <NA> 0.7296738 0.53242388
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod7[VC.H3.mod7$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.465259
```

```r
#see how much residual variance was explained at level-2 by both focus #parties

H3.mod8<-lmer(environ.gmc~(1|voting.group)+#(1|cntry)+
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

(FE.H3.mod8<-getFE(H3.mod8))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.02       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.05       0.01
## occupClerical support workers                               0.05       0.07
## occupCraft and related trades workers                      -0.02       0.07
## occupElementary occupations                                -0.04       0.07
## occupManagers                                               0.08       0.07
## occupOther: Not in paid work                                0.01       0.07
## occupPlant and machine operators, and assemblers            0.02       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.08       0.08
## occupService and sales workers                              0.00       0.07
## occupSkilled agricultural, forestry and fishery workers     0.00       0.07
## occupTechnicians and associate professionals                0.08       0.07
## occupUnemployed                                            -0.11       0.08
## did.not.vote.dummy                                         -0.11       0.03
## dont.know.dummy                                            -0.10       0.04
## invalid.vote.dummy                                         -0.38       0.27
## no.answer.dummy                                            -0.05       0.27
## not.eligible.age.dummy                                      0.01       0.04
## not.eligible.citizenship.dummy                              0.00       0.04
## not.eligible.other.dummy                                    0.12       0.06
## anti.imm.party.dummy                                       -0.15       0.03
## pro.env.party.dummy                                         0.30       0.03
##                                                               df t.value     p
## (Intercept)                                             23722.02   -0.25 0.800
## age                                                     26695.16   -6.70 0.000
## gender                                                  26867.90    4.90 0.000
## educ                                                    23779.10   14.96 0.000
## resid                                                   26366.56   -4.92 0.000
## occupClerical support workers                           26866.97    0.70 0.484
## occupCraft and related trades workers                   26871.80   -0.36 0.719
## occupElementary occupations                             26876.94   -0.58 0.565
## occupManagers                                           26866.05    1.24 0.213
## occupOther: Not in paid work                            26885.59    0.14 0.889
## occupPlant and machine operators, and assemblers        26876.47    0.30 0.764
## occupProfessionals                                      26857.68    1.52 0.129
## occupRetired                                            26864.78   -1.04 0.300
## occupService and sales workers                          26864.46   -0.03 0.975
## occupSkilled agricultural, forestry and fishery workers 26868.14   -0.06 0.951
## occupTechnicians and associate professionals            26861.43    1.16 0.245
## occupUnemployed                                         26875.78   -1.35 0.176
## did.not.vote.dummy                                        112.02   -4.01 0.000
## dont.know.dummy                                           260.03   -2.78 0.006
## invalid.vote.dummy                                       4555.01   -1.39 0.165
## no.answer.dummy                                          4555.65   -0.17 0.867
## not.eligible.age.dummy                                    287.48    0.23 0.821
## not.eligible.citizenship.dummy                            231.91    0.11 0.915
## not.eligible.other.dummy                                 1065.08    1.81 0.070
## anti.imm.party.dummy                                      165.07   -5.28 0.000
## pro.env.party.dummy                                       228.62    9.96 0.000
##                                                            LL    UL
## (Intercept)                                             -0.15  0.11
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.06 -0.03
## occupClerical support workers                           -0.09  0.18
## occupCraft and related trades workers                   -0.16  0.11
## occupElementary occupations                             -0.17  0.09
## occupManagers                                           -0.05  0.22
## occupOther: Not in paid work                            -0.13  0.15
## occupPlant and machine operators, and assemblers        -0.11  0.15
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.23  0.07
## occupService and sales workers                          -0.13  0.13
## occupSkilled agricultural, forestry and fishery workers -0.15  0.14
## occupTechnicians and associate professionals            -0.05  0.21
## occupUnemployed                                         -0.28  0.05
## did.not.vote.dummy                                      -0.17 -0.06
## dont.know.dummy                                         -0.17 -0.03
## invalid.vote.dummy                                      -0.92  0.16
## no.answer.dummy                                         -0.58  0.49
## not.eligible.age.dummy                                  -0.07  0.08
## not.eligible.citizenship.dummy                          -0.07  0.08
## not.eligible.other.dummy                                -0.01  0.24
## anti.imm.party.dummy                                    -0.20 -0.09
## pro.env.party.dummy                                      0.24  0.35
```

```r
(VC.H3.mod8<-getVC(H3.mod8))
```

```
##            grp        var1 var2     est_SD     est_SD2
## 1 voting.group (Intercept) <NA> 0.09154591 0.008380654
## 2     Residual        <NA> <NA> 0.72965308 0.532393616
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod8[VC.H3.mod8$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5617562
```

```r
#variance not accounted for at level-2
1-H3.total.eff
```

```
## [1] 0.4382438
```
