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

# Hypothesis 3: Those who voted for anti-immigration parties will report lower pro-environment attitudes than those who voted for pro-environment parties.

### Model 0: random intercepts for environment attitudes



```r
H3.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H3.mod0<-getFE(H3.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.04       0.04 19.97     1.1 0.284 -0.04 0.12
```

```r
(VC.H3.mod0<-getVC(H3.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1796220 0.03226405
## 2        cntry (Intercept) <NA> 0.1618768 0.02620409
## 3     Residual        <NA> <NA> 0.7969943 0.63519992
```

```r
#ICC

##voting group

VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## [1] 0.04651223
```

```r
##country

VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## [1] 0.03777612
```

\newpage

### Model 1: random intercepts + covariates


```r
H3.mod1<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)

anova(H3.mod0,H3.mod1)
```

```
## Data: dat
## Models:
## H3.mod0: environ.gmc ~ (1 | voting.group) + (1 | cntry)
## H3.mod1: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod1:     educ + resid + occup
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H3.mod0    4 88444 88479 -44218    88436                         
## H3.mod1   20 87387 87557 -43673    87347 1089.5 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H3.mod1<-getFE(H3.mod1))
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
(VC.H3.mod1<-getVC(H3.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1504843 0.02264552
## 2        cntry (Intercept) <NA> 0.1658533 0.02750733
## 3     Residual        <NA> <NA> 0.7859685 0.61774648
```

```r
#variance explained

##lvl 1: individuals

(VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="Residual","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]
```

```
## [1] 0.02747708
```

```r
##lvl 2: voting group

(VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.2981192
```

```r
##lvl 3: country

(VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="cntry","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]
```

```
## [1] -0.04973449
```

```r
##total

(sum(VC.H3.mod0$est_SD2)-sum(VC.H3.mod1$est_SD2))/
  sum(VC.H3.mod0$est_SD2)
```

```
## [1] 0.0371485
```

```r
#individual contributions of covariates
anova(H3.mod1)
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


### Model 2: Categorical predictor at level-2


```r
H3.mod2<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup+
                all.parties.lvl2
                ,data=dat,REML=F)


(FE.H3.mod2<-getFE(H3.mod2))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.15       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.03       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.16       0.08
## all.parties.lvl2Did not vote                                0.05       0.03
## all.parties.lvl2Don't know                                  0.07       0.04
## all.parties.lvl2Invalid vote                               -0.22       0.23
## all.parties.lvl2NE age                                      0.20       0.04
## all.parties.lvl2NE citizen                                  0.16       0.04
## all.parties.lvl2NE other                                    0.22       0.06
## all.parties.lvl2No answer                                   0.38       0.25
## all.parties.lvl2Other party                                 0.16       0.03
## all.parties.lvl2Pro-environment party                       0.47       0.04
##                                                               df t.value     p
## (Intercept)                                               406.79   -1.89 0.059
## age                                                     36665.93   -7.26 0.000
## gender                                                  36816.07    5.65 0.000
## educ                                                    36825.59   15.44 0.000
## resid                                                   36865.53   -7.04 0.000
## occupClerical support workers                           36782.33    0.52 0.607
## occupCraft and related trades workers                   36791.62   -0.61 0.541
## occupElementary occupations                             36787.49   -0.78 0.435
## occupManagers                                           36786.06    0.95 0.344
## occupOther: Not in paid work                            36836.85    0.05 0.957
## occupPlant and machine operators, and assemblers        36793.11   -0.52 0.605
## occupProfessionals                                      36783.46    1.51 0.132
## occupRetired                                            36779.23   -1.50 0.133
## occupService and sales workers                          36781.86   -0.19 0.850
## occupSkilled agricultural, forestry and fishery workers 36802.21    0.07 0.940
## occupTechnicians and associate professionals            36780.63    0.83 0.408
## occupUnemployed                                         36779.35   -2.05 0.041
## all.parties.lvl2Did not vote                              151.04    1.53 0.129
## all.parties.lvl2Don't know                                280.68    1.65 0.100
## all.parties.lvl2Invalid vote                             1913.82   -0.98 0.327
## all.parties.lvl2NE age                                    292.01    4.79 0.000
## all.parties.lvl2NE citizen                                255.46    3.69 0.000
## all.parties.lvl2NE other                                  856.14    3.51 0.000
## all.parties.lvl2No answer                                2718.74    1.51 0.132
## all.parties.lvl2Other party                               205.79    5.81 0.000
## all.parties.lvl2Pro-environment party                     246.29   12.49 0.000
##                                                            LL    UL
## (Intercept)                                             -0.30  0.01
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.04
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
## occupUnemployed                                         -0.32 -0.01
## all.parties.lvl2Did not vote                            -0.02  0.12
## all.parties.lvl2Don't know                              -0.01  0.15
## all.parties.lvl2Invalid vote                            -0.67  0.22
## all.parties.lvl2NE age                                   0.12  0.28
## all.parties.lvl2NE citizen                               0.08  0.25
## all.parties.lvl2NE other                                 0.10  0.34
## all.parties.lvl2No answer                               -0.11  0.87
## all.parties.lvl2Other party                              0.10  0.21
## all.parties.lvl2Pro-environment party                    0.39  0.54
```

```r
(VC.H3.mod2<-getVC(H3.mod2))
```

```
##            grp        var1 var2     est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.09806513 0.00961677
## 2        cntry (Intercept) <NA> 0.15860669 0.02515608
## 3     Residual        <NA> <NA> 0.78588249 0.61761129
```

```r
anova(H3.mod1,H3.mod2)
```

```
## Data: dat
## Models:
## H3.mod1: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod1:     educ + resid + occup
## H3.mod2: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod2:     educ + resid + occup + all.parties.lvl2
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H3.mod1   20 87387 87557 -43673    87347                         
## H3.mod2   29 87255 87502 -43599    87197 149.62  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H3.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age               32.573  32.573     1 36666  52.741 3.881e-13 ***
## gender            19.689  19.689     1 36816  31.879 1.653e-08 ***
## educ             147.306 147.306     1 36826 238.509 < 2.2e-16 ***
## resid             30.567  30.567     1 36866  49.492 2.026e-12 ***
## occup             82.405   6.867    12 36726  11.119 < 2.2e-16 ***
## all.parties.lvl2 119.611  13.290     9   338  21.518 < 2.2e-16 ***
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
## [1] 0.5753345
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
## 1  Anti-immigration party  -0.16 0.04 0.0003 0.0024     -0.24     -0.07
## 2            Did not vote  -0.11 0.04 0.0141 0.1126     -0.19     -0.02
## 3              Don't know  -0.09 0.05 0.0622 0.4357     -0.19      0.00
## 4            Invalid vote  -0.38 0.23 0.0962 0.5772     -0.84      0.07
## 5                  NE age   0.04 0.05 0.4326 1.0000     -0.06      0.13
## 6              NE citizen   0.00 0.05 0.9457 1.0000     -0.10      0.10
## 7                NE other   0.06 0.07 0.3784 1.0000     -0.07      0.19
## 8               No answer   0.22 0.25 0.3854 1.0000     -0.28      0.71
## 9             Other party   0.00 0.04 0.9807 1.0000     -0.08      0.07
## 10  Pro-environment party   0.31 0.05 0.0000 0.0000      0.22      0.40
```

```r
write.csv2(H3.mod2.mmeans.tab,"H3.mod2.mmeans.tab.csv")

#contrast between anti-immigration and pro-environment
(H3.contrast<-data.frame(pairs(H3.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.4671259 0.03739546 Inf
##    z.ratio      p.value
## 1 12.49152 8.305668e-36
```

```r
#contrast for all groups against mean of other groups
contrast(H3.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.1642 0.0460 Inf -3.572  0.0032 
##  Did not vote effect            -0.1058 0.0457 Inf -2.313  0.1658 
##  Don't know effect              -0.0887 0.0509 Inf -1.741  0.5021 
##  Invalid vote effect            -0.4137 0.2297 Inf -1.801  0.5021 
##  NE age effect                   0.0550 0.0508 Inf  1.083  1.0000 
##  NE citizen effect               0.0165 0.0531 Inf  0.310  1.0000 
##  NE other effect                 0.0789 0.0691 Inf  1.143  1.0000 
##  No answer effect                0.2556 0.2509 Inf  1.019  1.0000 
##  Other party effect              0.0116 0.0406 Inf  0.285  1.0000 
##  Pro-environment party effect    0.3548 0.0480 Inf  7.387  <.0001 
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
## 1           Other party - Anti-immigration party 0.1582143 0.02725281 Inf
## 2 Pro-environment party - Anti-immigration party 0.4671259 0.03739546 Inf
## 3            Pro-environment party - Other party 0.3089116 0.03034150 Inf
##     z.ratio      p.value
## 1  5.805431 6.420048e-09
## 2 12.491515 2.491700e-35
## 3 10.181159 4.813565e-24
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
## [1] 0.7977597
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
## [1] 0.7467473
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
## [1] 0.7726746
```

```r
(H3.effect.size<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd)
```

```
## [1] 0.6082768
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
## [1] 0.8005313
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
## [1] 0.7764067
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.16 0.04 0.0003 0.0024     -0.24     -0.07
## 2            Did not vote  -0.11 0.04 0.0141 0.1126     -0.19     -0.02
## 3              Don't know  -0.09 0.05 0.0622 0.4357     -0.19      0.00
## 4            Invalid vote  -0.38 0.23 0.0962 0.5772     -0.84      0.07
## 5                  NE age   0.04 0.05 0.4326 1.0000     -0.06      0.13
## 6              NE citizen   0.00 0.05 0.9457 1.0000     -0.10      0.10
## 7                NE other   0.06 0.07 0.3784 1.0000     -0.07      0.19
## 8               No answer   0.22 0.25 0.3854 1.0000     -0.28      0.71
## 9             Other party   0.00 0.04 0.9807 1.0000     -0.08      0.07
## 10  Pro-environment party   0.31 0.05 0.0000 0.0000      0.22      0.40
```

```r
(H3.effect.size.env.other<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[9,2])/H3.pooled.sd.other.env)
```

```
## [1] 0.3992753
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
## [1] 0.7992655
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party  -0.16 0.04 0.0003 0.0024     -0.24     -0.07
## 2            Did not vote  -0.11 0.04 0.0141 0.1126     -0.19     -0.02
## 3              Don't know  -0.09 0.05 0.0622 0.4357     -0.19      0.00
## 4            Invalid vote  -0.38 0.23 0.0962 0.5772     -0.84      0.07
## 5                  NE age   0.04 0.05 0.4326 1.0000     -0.06      0.13
## 6              NE citizen   0.00 0.05 0.9457 1.0000     -0.10      0.10
## 7                NE other   0.06 0.07 0.3784 1.0000     -0.07      0.19
## 8               No answer   0.22 0.25 0.3854 1.0000     -0.28      0.71
## 9             Other party   0.00 0.04 0.9807 1.0000     -0.08      0.07
## 10  Pro-environment party   0.31 0.05 0.0000 0.0000      0.22      0.40
```

```r
(H3.effect.size.imm.other<-(H3.mod2.mmeans.tab[9,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd.other.imm)
```

```
## [1] 0.2001838
```


\newpage


### Model 3: Dummy-predictors at level-2 


```r
#did not vote left as reference

H3.mod3<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                -0.09       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.03       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.16       0.08
## other.party.dummy                                           0.11       0.03
## dont.know.dummy                                             0.02       0.04
## invalid.vote.dummy                                         -0.28       0.23
## no.answer.dummy                                             0.33       0.25
## not.eligible.age.dummy                                      0.14       0.04
## not.eligible.citizenship.dummy                              0.11       0.04
## not.eligible.other.dummy                                    0.17       0.06
## anti.imm.party.dummy                                       -0.05       0.03
## pro.env.party.dummy                                         0.41       0.04
##                                                               df t.value     p
## (Intercept)                                               402.82   -1.22 0.225
## age                                                     36665.93   -7.26 0.000
## gender                                                  36816.07    5.65 0.000
## educ                                                    36825.59   15.44 0.000
## resid                                                   36865.53   -7.04 0.000
## occupClerical support workers                           36782.33    0.52 0.607
## occupCraft and related trades workers                   36791.62   -0.61 0.541
## occupElementary occupations                             36787.49   -0.78 0.435
## occupManagers                                           36786.06    0.95 0.344
## occupOther: Not in paid work                            36836.85    0.05 0.957
## occupPlant and machine operators, and assemblers        36793.12   -0.52 0.605
## occupProfessionals                                      36783.46    1.51 0.132
## occupRetired                                            36779.24   -1.50 0.133
## occupService and sales workers                          36781.86   -0.19 0.850
## occupSkilled agricultural, forestry and fishery workers 36802.21    0.07 0.940
## occupTechnicians and associate professionals            36780.63    0.83 0.408
## occupUnemployed                                         36779.35   -2.05 0.041
## other.party.dummy                                         130.15    3.93 0.000
## dont.know.dummy                                           219.53    0.38 0.708
## invalid.vote.dummy                                       1884.32   -1.21 0.226
## no.answer.dummy                                          2659.21    1.30 0.194
## not.eligible.age.dummy                                    216.16    3.57 0.000
## not.eligible.citizenship.dummy                            205.82    2.52 0.013
## not.eligible.other.dummy                                  712.99    2.67 0.008
## anti.imm.party.dummy                                      151.04   -1.53 0.129
## pro.env.party.dummy                                       187.13   11.15 0.000
##                                                            LL    UL
## (Intercept)                                             -0.25  0.06
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.04
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
## occupUnemployed                                         -0.32 -0.01
## other.party.dummy                                        0.05  0.16
## dont.know.dummy                                         -0.07  0.10
## invalid.vote.dummy                                      -0.73  0.17
## no.answer.dummy                                         -0.17  0.82
## not.eligible.age.dummy                                   0.06  0.22
## not.eligible.citizenship.dummy                           0.02  0.20
## not.eligible.other.dummy                                 0.04  0.29
## anti.imm.party.dummy                                    -0.12  0.02
## pro.env.party.dummy                                      0.34  0.49
```

```r
(VC.H3.mod3<-getVC(H3.mod3))
```

```
##            grp        var1 var2     est_SD     est_SD2
## 1 voting.group (Intercept) <NA> 0.09806513 0.009616769
## 2        cntry (Intercept) <NA> 0.15860702 0.025156186
## 3     Residual        <NA> <NA> 0.78588249 0.617611292
```

```r
anova(H3.mod2,H3.mod3)
```

```
## Data: dat
## Models:
## H3.mod2: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod2:     educ + resid + occup + all.parties.lvl2
## H3.mod3: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod3:     educ + resid + occup + other.party.dummy + dont.know.dummy + 
## H3.mod3:     invalid.vote.dummy + no.answer.dummy + not.eligible.age.dummy + 
## H3.mod3:     not.eligible.citizenship.dummy + not.eligible.other.dummy + 
## H3.mod3:     anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)    
## H3.mod2   29 87255 87502 -43599    87197                        
## H3.mod3   29 87255 87502 -43599    87197     0  0  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


\newpage


### Model 4: Dummy-predictors at level-2 allowed to vary between countries


```r
#did not vote left as reference

H3.mod4<-lmer(environ.gmc~(1|voting.group)+(anti.imm.party.dummy+pro.env.party.dummy||cntry)+
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
## [1] FALSE
```

```r
(FE.H3.mod4<-getFE(H3.mod4))
```

```
##                                                         Estimate Std..Error
## (Intercept)                                                -0.10       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.06       0.07
## occupUnemployed                                            -0.16       0.08
## other.party.dummy                                           0.11       0.02
## dont.know.dummy                                             0.02       0.04
## invalid.vote.dummy                                         -0.27       0.22
## no.answer.dummy                                             0.32       0.25
## not.eligible.age.dummy                                      0.14       0.04
## not.eligible.citizenship.dummy                              0.11       0.04
## not.eligible.other.dummy                                    0.17       0.06
## anti.imm.party.dummy                                       -0.05       0.04
## pro.env.party.dummy                                         0.41       0.05
##                                                               df t.value     p
## (Intercept)                                               416.85   -1.25 0.212
## age                                                     36575.55   -7.26 0.000
## gender                                                  36819.39    5.64 0.000
## educ                                                    36778.05   15.44 0.000
## resid                                                   36838.49   -7.03 0.000
## occupClerical support workers                           36785.96    0.54 0.592
## occupCraft and related trades workers                   36795.29   -0.59 0.553
## occupElementary occupations                             36790.99   -0.76 0.446
## occupManagers                                           36789.56    0.97 0.334
## occupOther: Not in paid work                            36838.05    0.06 0.949
## occupPlant and machine operators, and assemblers        36796.96   -0.50 0.616
## occupProfessionals                                      36786.78    1.52 0.128
## occupRetired                                            36783.04   -1.48 0.138
## occupService and sales workers                          36785.28   -0.18 0.860
## occupSkilled agricultural, forestry and fishery workers 36806.72    0.09 0.930
## occupTechnicians and associate professionals            36783.97    0.84 0.398
## occupUnemployed                                         36784.04   -2.03 0.042
## other.party.dummy                                         103.64    4.26 0.000
## dont.know.dummy                                           187.35    0.41 0.681
## invalid.vote.dummy                                       1918.12   -1.22 0.222
## no.answer.dummy                                          2725.96    1.32 0.188
## not.eligible.age.dummy                                    184.68    3.78 0.000
## not.eligible.citizenship.dummy                            172.03    2.70 0.008
## not.eligible.other.dummy                                  657.12    2.76 0.006
## anti.imm.party.dummy                                       30.62   -1.40 0.172
## pro.env.party.dummy                                        26.19    8.81 0.000
##                                                            LL    UL
## (Intercept)                                             -0.25  0.06
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.04
## occupClerical support workers                           -0.09  0.17
## occupCraft and related trades workers                   -0.17  0.09
## occupElementary occupations                             -0.18  0.08
## occupManagers                                           -0.07  0.20
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.16  0.10
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.25  0.04
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.13  0.14
## occupTechnicians and associate professionals            -0.07  0.18
## occupUnemployed                                         -0.32 -0.01
## other.party.dummy                                        0.06  0.16
## dont.know.dummy                                         -0.06  0.09
## invalid.vote.dummy                                      -0.72  0.17
## no.answer.dummy                                         -0.16  0.81
## not.eligible.age.dummy                                   0.07  0.22
## not.eligible.citizenship.dummy                           0.03  0.19
## not.eligible.other.dummy                                 0.05  0.29
## anti.imm.party.dummy                                    -0.12  0.02
## pro.env.party.dummy                                      0.31  0.50
```

```r
(VC.H3.mod4<-getVC(H3.mod4))
```

```
##            grp                 var1 var2     est_SD     est_SD2
## 1 voting.group          (Intercept) <NA> 0.08868514 0.007865054
## 2        cntry          (Intercept) <NA> 0.15504735 0.024039680
## 3      cntry.1 anti.imm.party.dummy <NA> 0.06347956 0.004029655
## 4      cntry.2  pro.env.party.dummy <NA> 0.11858148 0.014061568
## 5     Residual                 <NA> <NA> 0.78592300 0.617674956
```

```r
anova(H3.mod3,H3.mod4)
```

```
## Data: dat
## Models:
## H3.mod3: environ.gmc ~ (1 | voting.group) + (1 | cntry) + age + gender + 
## H3.mod3:     educ + resid + occup + other.party.dummy + dont.know.dummy + 
## H3.mod3:     invalid.vote.dummy + no.answer.dummy + not.eligible.age.dummy + 
## H3.mod3:     not.eligible.citizenship.dummy + not.eligible.other.dummy + 
## H3.mod3:     anti.imm.party.dummy + pro.env.party.dummy
## H3.mod4: environ.gmc ~ (1 | voting.group) + (anti.imm.party.dummy + pro.env.party.dummy || 
## H3.mod4:     cntry) + age + gender + educ + resid + occup + other.party.dummy + 
## H3.mod4:     dont.know.dummy + invalid.vote.dummy + no.answer.dummy + 
## H3.mod4:     not.eligible.age.dummy + not.eligible.citizenship.dummy + 
## H3.mod4:     not.eligible.other.dummy + anti.imm.party.dummy + pro.env.party.dummy
##         npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H3.mod3   29 87255 87502 -43599    87197                       
## H3.mod4   31 87253 87517 -43595    87191 6.5836  2    0.03719 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


\newpage

### Model 5: explained variance by the focus groups


```r
H3.mod5<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.02       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.06       0.07
## occupUnemployed                                            -0.16       0.08
## did.not.vote.dummy                                         -0.12       0.04
## dont.know.dummy                                            -0.11       0.04
## invalid.vote.dummy                                         -0.38       0.25
## no.answer.dummy                                             0.20       0.27
## not.eligible.age.dummy                                      0.03       0.04
## not.eligible.citizenship.dummy                             -0.01       0.05
## not.eligible.other.dummy                                    0.05       0.07
##                                                               df t.value     p
## (Intercept)                                               295.34    0.32 0.748
## age                                                     36839.41   -7.37 0.000
## gender                                                  36776.29    5.76 0.000
## educ                                                    36861.67   15.73 0.000
## resid                                                   36862.40   -7.15 0.000
## occupClerical support workers                           36728.22    0.53 0.593
## occupCraft and related trades workers                   36736.57   -0.61 0.541
## occupElementary occupations                             36733.20   -0.76 0.445
## occupManagers                                           36730.37    0.93 0.350
## occupOther: Not in paid work                            36782.22    0.06 0.951
## occupPlant and machine operators, and assemblers        36737.01   -0.51 0.610
## occupProfessionals                                      36731.66    1.53 0.125
## occupRetired                                            36726.73   -1.49 0.136
## occupService and sales workers                          36727.98   -0.19 0.848
## occupSkilled agricultural, forestry and fishery workers 36745.86    0.08 0.934
## occupTechnicians and associate professionals            36726.50    0.84 0.402
## occupUnemployed                                         36725.93   -2.02 0.043
## did.not.vote.dummy                                        152.41   -3.34 0.001
## dont.know.dummy                                           283.42   -2.48 0.014
## invalid.vote.dummy                                       1022.22   -1.51 0.131
## no.answer.dummy                                          1372.17    0.75 0.456
## not.eligible.age.dummy                                    293.12    0.64 0.525
## not.eligible.citizenship.dummy                            274.40   -0.14 0.885
## not.eligible.other.dummy                                  878.87    0.82 0.415
##                                                            LL    UL
## (Intercept)                                             -0.12  0.17
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
## did.not.vote.dummy                                      -0.19 -0.05
## dont.know.dummy                                         -0.19 -0.02
## invalid.vote.dummy                                      -0.87  0.11
## no.answer.dummy                                         -0.33  0.73
## not.eligible.age.dummy                                  -0.06  0.11
## not.eligible.citizenship.dummy                          -0.10  0.09
## not.eligible.other.dummy                                -0.08  0.18
```

```r
(VC.H3.mod5<-getVC(H3.mod5))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1425061 0.02030800
## 2        cntry (Intercept) <NA> 0.1663984 0.02768842
## 3     Residual        <NA> <NA> 0.7859618 0.61773597
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H3.total.eff<-(VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"]-
     VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])/
  VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.103222
```

```r
#see how much residual variance was explained at level-2 by anti-immigrants

H3.mod6<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.06       0.08
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.04       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.06       0.07
## occupUnemployed                                            -0.16       0.08
## did.not.vote.dummy                                         -0.15       0.03
## dont.know.dummy                                            -0.14       0.04
## invalid.vote.dummy                                         -0.41       0.24
## no.answer.dummy                                             0.17       0.26
## not.eligible.age.dummy                                     -0.01       0.04
## not.eligible.citizenship.dummy                             -0.04       0.04
## not.eligible.other.dummy                                    0.02       0.06
## anti.imm.party.dummy                                       -0.21       0.03
##                                                               df t.value     p
## (Intercept)                                               290.67    0.74 0.459
## age                                                     36800.81   -7.49 0.000
## gender                                                  36782.18    5.71 0.000
## educ                                                    36860.64   15.64 0.000
## resid                                                   36872.37   -7.12 0.000
## occupClerical support workers                           36742.91    0.53 0.595
## occupCraft and related trades workers                   36752.44   -0.59 0.557
## occupElementary occupations                             36748.37   -0.76 0.448
## occupManagers                                           36745.52    0.94 0.348
## occupOther: Not in paid work                            36799.88    0.07 0.945
## occupPlant and machine operators, and assemblers        36752.96   -0.50 0.618
## occupProfessionals                                      36746.14    1.54 0.125
## occupRetired                                            36740.77   -1.48 0.138
## occupService and sales workers                          36742.81   -0.18 0.857
## occupSkilled agricultural, forestry and fishery workers 36761.76    0.08 0.934
## occupTechnicians and associate professionals            36741.34    0.84 0.400
## occupUnemployed                                         36739.46   -2.02 0.043
## did.not.vote.dummy                                        142.67   -4.63 0.000
## dont.know.dummy                                           286.30   -3.45 0.001
## invalid.vote.dummy                                       1212.66   -1.67 0.095
## no.answer.dummy                                          1657.07    0.65 0.517
## not.eligible.age.dummy                                    300.27   -0.19 0.848
## not.eligible.citizenship.dummy                            267.20   -0.92 0.359
## not.eligible.other.dummy                                  930.25    0.27 0.789
## anti.imm.party.dummy                                      211.51   -6.47 0.000
##                                                            LL    UL
## (Intercept)                                             -0.09  0.21
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.05
## occupClerical support workers                           -0.10  0.17
## occupCraft and related trades workers                   -0.17  0.09
## occupElementary occupations                             -0.18  0.08
## occupManagers                                           -0.07  0.19
## occupOther: Not in paid work                            -0.13  0.14
## occupPlant and machine operators, and assemblers        -0.16  0.10
## occupProfessionals                                      -0.03  0.23
## occupRetired                                            -0.25  0.04
## occupService and sales workers                          -0.14  0.12
## occupSkilled agricultural, forestry and fishery workers -0.13  0.14
## occupTechnicians and associate professionals            -0.07  0.18
## occupUnemployed                                         -0.32  0.00
## did.not.vote.dummy                                      -0.22 -0.09
## dont.know.dummy                                         -0.22 -0.06
## invalid.vote.dummy                                      -0.88  0.07
## no.answer.dummy                                         -0.35  0.69
## not.eligible.age.dummy                                  -0.09  0.07
## not.eligible.citizenship.dummy                          -0.13  0.05
## not.eligible.other.dummy                                -0.11  0.14
## anti.imm.party.dummy                                    -0.27 -0.14
```

```r
(VC.H3.mod6<-getVC(H3.mod6))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1271714 0.01617256
## 2        cntry (Intercept) <NA> 0.1686866 0.02845518
## 3     Residual        <NA> <NA> 0.7859656 0.61774185
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod6[VC.H3.mod6$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.2036361
```

```r
#see how much residual variance was explained at level-2 by pro-environments

H3.mod7<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                -0.02       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.03       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.04       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.16       0.08
## did.not.vote.dummy                                         -0.08       0.03
## dont.know.dummy                                            -0.06       0.04
## invalid.vote.dummy                                         -0.36       0.23
## no.answer.dummy                                             0.25       0.25
## not.eligible.age.dummy                                      0.07       0.04
## not.eligible.citizenship.dummy                              0.04       0.04
## not.eligible.other.dummy                                    0.09       0.06
## pro.env.party.dummy                                         0.34       0.03
##                                                               df t.value     p
## (Intercept)                                               367.46   -0.25 0.806
## age                                                     36731.47   -7.14 0.000
## gender                                                  36810.63    5.70 0.000
## educ                                                    36832.25   15.54 0.000
## resid                                                   36874.70   -7.07 0.000
## occupClerical support workers                           36770.09    0.52 0.604
## occupCraft and related trades workers                   36778.70   -0.63 0.527
## occupElementary occupations                             36775.03   -0.78 0.433
## occupManagers                                           36773.46    0.95 0.344
## occupOther: Not in paid work                            36824.84    0.05 0.962
## occupPlant and machine operators, and assemblers        36780.02   -0.53 0.598
## occupProfessionals                                      36771.45    1.51 0.131
## occupRetired                                            36767.16   -1.51 0.131
## occupService and sales workers                          36769.50   -0.20 0.842
## occupSkilled agricultural, forestry and fishery workers 36789.39    0.08 0.938
## occupTechnicians and associate professionals            36768.30    0.83 0.409
## occupUnemployed                                         36766.73   -2.04 0.041
## did.not.vote.dummy                                        137.68   -2.64 0.009
## dont.know.dummy                                           318.38   -1.67 0.097
## invalid.vote.dummy                                       1671.82   -1.54 0.123
## no.answer.dummy                                          2329.74    0.98 0.330
## not.eligible.age.dummy                                    338.90    1.92 0.055
## not.eligible.citizenship.dummy                            281.63    0.87 0.386
## not.eligible.other.dummy                                 1107.40    1.56 0.120
## pro.env.party.dummy                                       276.51   10.67 0.000
##                                                            LL    UL
## (Intercept)                                             -0.16  0.13
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.04
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
## occupUnemployed                                         -0.32 -0.01
## did.not.vote.dummy                                      -0.13 -0.02
## dont.know.dummy                                         -0.13  0.01
## invalid.vote.dummy                                      -0.82  0.10
## no.answer.dummy                                         -0.25  0.75
## not.eligible.age.dummy                                   0.00  0.14
## not.eligible.citizenship.dummy                          -0.04  0.11
## not.eligible.other.dummy                                -0.02  0.21
## pro.env.party.dummy                                      0.28  0.40
```

```r
(VC.H3.mod7<-getVC(H3.mod7))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1094496 0.01197922
## 2        cntry (Intercept) <NA> 0.1556582 0.02422948
## 3     Residual        <NA> <NA> 0.7858840 0.61761374
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod7[VC.H3.mod7$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.4101233
```

```r
#see how much residual variance was explained at level-2 by both focus #parties

H3.mod8<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
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
## (Intercept)                                                 0.01       0.07
## age                                                         0.00       0.00
## gender                                                      0.05       0.01
## educ                                                        0.02       0.00
## resid                                                      -0.06       0.01
## occupClerical support workers                               0.03       0.07
## occupCraft and related trades workers                      -0.04       0.07
## occupElementary occupations                                -0.05       0.07
## occupManagers                                               0.06       0.07
## occupOther: Not in paid work                                0.00       0.07
## occupPlant and machine operators, and assemblers           -0.03       0.07
## occupProfessionals                                          0.10       0.07
## occupRetired                                               -0.11       0.07
## occupService and sales workers                             -0.01       0.07
## occupSkilled agricultural, forestry and fishery workers     0.01       0.07
## occupTechnicians and associate professionals                0.05       0.07
## occupUnemployed                                            -0.16       0.08
## did.not.vote.dummy                                         -0.11       0.03
## dont.know.dummy                                            -0.09       0.03
## invalid.vote.dummy                                         -0.38       0.23
## no.answer.dummy                                             0.22       0.25
## not.eligible.age.dummy                                      0.04       0.04
## not.eligible.citizenship.dummy                              0.00       0.04
## not.eligible.other.dummy                                    0.06       0.06
## anti.imm.party.dummy                                       -0.16       0.03
## pro.env.party.dummy                                         0.31       0.03
##                                                               df t.value     p
## (Intercept)                                               354.25    0.15 0.883
## age                                                     36665.93   -7.26 0.000
## gender                                                  36816.07    5.65 0.000
## educ                                                    36825.59   15.44 0.000
## resid                                                   36865.53   -7.04 0.000
## occupClerical support workers                           36782.33    0.52 0.607
## occupCraft and related trades workers                   36791.62   -0.61 0.541
## occupElementary occupations                             36787.49   -0.78 0.435
## occupManagers                                           36786.06    0.95 0.344
## occupOther: Not in paid work                            36836.85    0.05 0.957
## occupPlant and machine operators, and assemblers        36793.11   -0.52 0.605
## occupProfessionals                                      36783.46    1.51 0.132
## occupRetired                                            36779.23   -1.50 0.133
## occupService and sales workers                          36781.86   -0.19 0.850
## occupSkilled agricultural, forestry and fishery workers 36802.21    0.07 0.940
## occupTechnicians and associate professionals            36780.63    0.83 0.408
## occupUnemployed                                         36779.35   -2.05 0.041
## did.not.vote.dummy                                        130.16   -3.93 0.000
## dont.know.dummy                                           324.62   -2.58 0.010
## invalid.vote.dummy                                       2004.91   -1.68 0.093
## no.answer.dummy                                          2822.41    0.88 0.379
## not.eligible.age.dummy                                    351.27    1.11 0.269
## not.eligible.citizenship.dummy                            276.74    0.12 0.908
## not.eligible.other.dummy                                 1175.51    1.03 0.302
## anti.imm.party.dummy                                      205.79   -5.81 0.000
## pro.env.party.dummy                                       269.32   10.18 0.000
##                                                            LL    UL
## (Intercept)                                             -0.14  0.16
## age                                                      0.00  0.00
## gender                                                   0.03  0.07
## educ                                                     0.02  0.02
## resid                                                   -0.08 -0.04
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
## occupUnemployed                                         -0.32 -0.01
## did.not.vote.dummy                                      -0.16 -0.05
## dont.know.dummy                                         -0.16 -0.02
## invalid.vote.dummy                                      -0.83  0.06
## no.answer.dummy                                         -0.27  0.71
## not.eligible.age.dummy                                  -0.03  0.11
## not.eligible.citizenship.dummy                          -0.07  0.08
## not.eligible.other.dummy                                -0.05  0.18
## anti.imm.party.dummy                                    -0.21 -0.10
## pro.env.party.dummy                                      0.25  0.37
```

```r
(VC.H3.mod8<-getVC(H3.mod8))
```

```
##            grp        var1 var2     est_SD     est_SD2
## 1 voting.group (Intercept) <NA> 0.09806508 0.009616759
## 2        cntry (Intercept) <NA> 0.15860739 0.025156304
## 3     Residual        <NA> <NA> 0.78588249 0.617611293
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod8[VC.H3.mod8$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5264547
```

```r
#variance not accounted for at level-2
1-H3.total.eff
```

```
## [1] 0.4735453
```
