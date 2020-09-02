---
title: "Hypothesis 5"
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
library(finalfit)
```

## Session information about the packages


```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19041)
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
##  [1] finalfit_1.0.2  merTools_0.5.0  arm_1.10-1      MASS_7.3-51.5  
##  [5] metafor_2.4-0   ggplot2_3.3.0   emmeans_1.4.5   psych_1.9.12.31
##  [9] dplyr_0.8.5     lmerTest_3.1-2  lme4_1.1-23     Matrix_1.2-18  
## 
## loaded via a namespace (and not attached):
##  [1] tidyr_1.0.2         splines_3.6.3       foreach_1.5.0      
##  [4] shiny_1.4.0.2       assertthat_0.2.1    statmod_1.4.34     
##  [7] yaml_2.2.0          numDeriv_2016.8-1.1 pillar_1.4.3       
## [10] backports_1.1.6     lattice_0.20-38     glue_1.4.0         
## [13] digest_0.6.25       promises_1.1.0      minqa_1.2.4        
## [16] colorspace_1.4-1    sandwich_2.5-1      htmltools_0.4.0    
## [19] httpuv_1.5.2        pkgconfig_2.0.3     broom_0.5.5        
## [22] purrr_0.3.3         xtable_1.8-4        mvtnorm_1.1-0      
## [25] scales_1.1.0        later_1.0.0         tibble_3.0.0       
## [28] generics_0.0.2      ellipsis_0.3.0      TH.data_1.0-10     
## [31] withr_2.1.2         cli_2.0.2           mnormt_1.5-6       
## [34] survival_3.1-8      magrittr_1.5        crayon_1.3.4       
## [37] mime_0.9            estimability_1.3    evaluate_0.14      
## [40] mice_3.8.0          fansi_0.4.1         nlme_3.1-144       
## [43] forcats_0.5.0       tools_3.6.3         lifecycle_0.2.0    
## [46] multcomp_1.4-13     stringr_1.4.0       munsell_0.5.0      
## [49] compiler_3.6.3      rlang_0.4.5         blme_1.0-4         
## [52] grid_3.6.3          nloptr_1.2.2.1      iterators_1.0.12   
## [55] rmarkdown_2.1       boot_1.3-24         gtable_0.3.0       
## [58] codetools_0.2-16    abind_1.4-5         R6_2.4.1           
## [61] zoo_1.8-7           knitr_1.28          fastmap_1.0.1      
## [64] stringi_1.4.6       parallel_3.6.3      Rcpp_1.0.4.6       
## [67] vctrs_0.2.4         tidyselect_1.0.0    xfun_0.13          
## [70] coda_0.19-3
```

\newpage

## Custom functions


```r
#to extract fixed effects
getFE<-function(model){
  
  coefs<-data.frame(summary(model)$coefficients)
  coefs$lower<-coefs[,1]-qt(p=.975,df=coefs[,"df"])*coefs[,2]
  coefs$upper<-coefs[,1]+qt(p=.975,df=coefs[,"df"])*coefs[,2]
  
  
  coefs<-cbind.data.frame(Eff=rownames(coefs),
                          Est=round_tidy(coefs[,1],2),
                          SE=round_tidy(coefs[,2],2),
                          df=round_tidy(coefs[,3],2),
                          t=round_tidy(coefs[,4],2),
                          p=round_tidy(coefs[,5],3),
                          LL=round_tidy(coefs$lower,2),
                          UL=round_tidy(coefs$upper,2))
  
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
dat$voting.group<-
  paste0(dat$cntry,": ",dat$vote.group.combined)
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

\newpage


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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.43
## 2                                                      age  0.00 0.00 33758.85
## 3                                                   gender  0.07 0.01 35565.52
## 4                                                     educ  0.01 0.00 35675.86
## 5                                                    resid -0.06 0.01 35635.00
## 6                            occupClerical support workers -0.04 0.09 35509.19
## 7                    occupCraft and related trades workers -0.06 0.09 35519.51
## 8                              occupElementary occupations  0.02 0.09 35521.98
## 9                                            occupManagers -0.02 0.09 35508.37
## 10                            occupOther: Not in paid work  0.14 0.09 35676.02
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.21
## 12                                      occupProfessionals  0.07 0.09 35513.08
## 13                                            occupRetired -0.03 0.10 35510.01
## 14                          occupService and sales workers -0.03 0.09 35515.16
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35524.18
## 16            occupTechnicians and associate professionals -0.03 0.09 35509.06
## 17                                         occupUnemployed -0.02 0.11 35526.75
## 18                                            environ.lvl1  0.12 0.01    19.45
## 19                                         engagement.lvl1  0.06 0.01 35609.62
##        t     p    LL    UL
## 1   0.48 0.634 -0.22  0.35
## 2   2.39 0.017  0.00  0.00
## 3   5.66 0.000  0.04  0.09
## 4   6.63 0.000  0.01  0.02
## 5  -4.98 0.000 -0.08 -0.04
## 6  -0.42 0.673 -0.21  0.14
## 7  -0.67 0.502 -0.23  0.11
## 8   0.27 0.787 -0.15  0.20
## 9  -0.20 0.842 -0.19  0.16
## 10  1.56 0.118 -0.04  0.32
## 11 -0.64 0.522 -0.23  0.12
## 12  0.83 0.409 -0.10  0.24
## 13 -0.28 0.776 -0.22  0.16
## 14 -0.35 0.725 -0.20  0.14
## 15 -0.53 0.598 -0.23  0.13
## 16 -0.33 0.738 -0.20  0.14
## 17 -0.16 0.876 -0.23  0.19
## 18 11.43 0.000  0.10  0.14
## 19  7.86 0.000  0.05  0.08
```

```r
(VC.H5.mod1<-getVC(H5.mod1))
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30760759  0.094622430
## 2 voting.group environ.lvl1         <NA>  0.04320212  0.001866423
## 3 voting.group  (Intercept) environ.lvl1  0.10968739  0.001457669
## 4        cntry  (Intercept)         <NA>  0.49927636  0.249276886
## 5        cntry environ.lvl1         <NA>  0.03945867  0.001556987
## 6        cntry  (Intercept) environ.lvl1 -0.06954526 -0.001370096
## 7     Residual         <NA>         <NA>  1.02826413  1.057327131
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.mod1   26 104230 104450 -52089   104178                     
## H5.mod2   27 104231 104460 -52088   104177 1.0709  1     0.3007
```

```r
(FE.H5.mod2<-getFE(H5.mod2))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.44
## 2                                                      age  0.00 0.00 33750.64
## 3                                                   gender  0.07 0.01 35565.68
## 4                                                     educ  0.01 0.00 35675.75
## 5                                                    resid -0.06 0.01 35635.19
## 6                            occupClerical support workers -0.04 0.09 35509.24
## 7                    occupCraft and related trades workers -0.06 0.09 35519.58
## 8                              occupElementary occupations  0.02 0.09 35522.03
## 9                                            occupManagers -0.02 0.09 35508.41
## 10                            occupOther: Not in paid work  0.14 0.09 35676.21
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.29
## 12                                      occupProfessionals  0.07 0.09 35513.13
## 13                                            occupRetired -0.03 0.10 35510.04
## 14                          occupService and sales workers -0.03 0.09 35515.24
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35524.27
## 16            occupTechnicians and associate professionals -0.03 0.09 35509.12
## 17                                         occupUnemployed -0.02 0.11 35526.77
## 18                                            environ.lvl1  0.12 0.01    19.46
## 19                                         engagement.lvl1  0.06 0.01 35610.30
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01 35543.53
##        t     p    LL    UL
## 1   0.48 0.637 -0.22  0.35
## 2   2.38 0.017  0.00  0.00
## 3   5.67 0.000  0.04  0.09
## 4   6.64 0.000  0.01  0.02
## 5  -4.97 0.000 -0.08 -0.04
## 6  -0.42 0.674 -0.21  0.14
## 7  -0.67 0.504 -0.23  0.11
## 8   0.27 0.786 -0.15  0.20
## 9  -0.20 0.843 -0.19  0.16
## 10  1.57 0.117 -0.04  0.32
## 11 -0.64 0.525 -0.23  0.12
## 12  0.83 0.408 -0.10  0.24
## 13 -0.28 0.777 -0.22  0.16
## 14 -0.35 0.728 -0.20  0.14
## 15 -0.52 0.600 -0.23  0.13
## 16 -0.33 0.741 -0.20  0.14
## 17 -0.16 0.875 -0.23  0.19
## 18 11.42 0.000  0.10  0.14
## 19  7.87 0.000  0.05  0.08
## 20  1.03 0.301 -0.01  0.02
```

```r
(VC.H5.mod2<-getVC(H5.mod2))
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30746398  0.094534097
## 2 voting.group environ.lvl1         <NA>  0.04325593  0.001871075
## 3 voting.group  (Intercept) environ.lvl1  0.10913820  0.001451499
## 4        cntry  (Intercept)         <NA>  0.49918132  0.249181994
## 5        cntry environ.lvl1         <NA>  0.03942038  0.001553966
## 6        cntry  (Intercept) environ.lvl1 -0.07008820 -0.001379190
## 7     Residual         <NA>         <NA>  1.02825024  1.057298557
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
##         npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)    
## H5.mod2   27 104231 104460 -52088   104177                        
## H5.mod3   33 104203 104483 -52069   104137 39.18  6  6.599e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod3<-getFE(H5.mod3))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.06 0.14    49.36
## 2                                                      age  0.00 0.00 33343.34
## 3                                                   gender  0.07 0.01 35558.07
## 4                                                     educ  0.01 0.00 35664.79
## 5                                                    resid -0.06 0.01 35611.45
## 6                            occupClerical support workers -0.03 0.09 35484.87
## 7                    occupCraft and related trades workers -0.06 0.09 35493.52
## 8                              occupElementary occupations  0.03 0.09 35493.55
## 9                                            occupManagers -0.01 0.09 35483.48
## 10                            occupOther: Not in paid work  0.15 0.09 35654.04
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35496.23
## 12                                      occupProfessionals  0.07 0.09 35488.49
## 13                                            occupRetired -0.03 0.10 35480.11
## 14                          occupService and sales workers -0.03 0.09 35487.34
## 15 occupSkilled agricultural, forestry and fishery workers -0.04 0.09 35503.28
## 16            occupTechnicians and associate professionals -0.03 0.09 35485.04
## 17                                         occupUnemployed -0.02 0.11 35505.62
## 18                                            environ.lvl1  0.12 0.01    19.45
## 19                                         engagement.lvl1  0.06 0.01    21.71
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01 35519.96
##        t     p    LL    UL
## 1   0.46 0.651 -0.22  0.35
## 2   2.59 0.010  0.00  0.00
## 3   5.68 0.000  0.04  0.09
## 4   6.57 0.000  0.01  0.02
## 5  -4.94 0.000 -0.08 -0.03
## 6  -0.39 0.699 -0.21  0.14
## 7  -0.63 0.527 -0.23  0.12
## 8   0.31 0.755 -0.15  0.20
## 9  -0.15 0.877 -0.19  0.16
## 10  1.62 0.106 -0.03  0.32
## 11 -0.59 0.555 -0.23  0.12
## 12  0.84 0.402 -0.10  0.24
## 13 -0.27 0.785 -0.22  0.16
## 14 -0.31 0.758 -0.20  0.14
## 15 -0.48 0.634 -0.23  0.14
## 16 -0.30 0.767 -0.20  0.15
## 17 -0.16 0.870 -0.23  0.19
## 18 11.34 0.000  0.10  0.14
## 19  4.38 0.000  0.03  0.09
## 20  1.05 0.295 -0.01  0.02
```

```r
(VC.H5.mod3<-getVC(H5.mod3))
```

```
##             grp            var1            var2      est_SD       est_SD2
## 1  voting.group     (Intercept)            <NA>  0.30798224  9.485306e-02
## 2  voting.group    environ.lvl1            <NA>  0.04340675  1.884146e-03
## 3  voting.group engagement.lvl1            <NA>  0.06210001  3.856411e-03
## 4  voting.group     (Intercept)    environ.lvl1  0.08789352  1.175005e-03
## 5  voting.group     (Intercept) engagement.lvl1  0.31575322  6.039001e-03
## 6  voting.group    environ.lvl1 engagement.lvl1  0.04962916  1.337784e-04
## 7         cntry     (Intercept)            <NA>  0.49938528  2.493857e-01
## 8         cntry    environ.lvl1            <NA>  0.03943128  1.554826e-03
## 9         cntry engagement.lvl1            <NA>  0.04626914  2.140833e-03
## 10        cntry     (Intercept)    environ.lvl1 -0.08175091 -1.609790e-03
## 11        cntry     (Intercept) engagement.lvl1  0.54595401  1.261488e-02
## 12        cntry    environ.lvl1 engagement.lvl1  0.04572445  8.342203e-05
## 13     Residual            <NA>            <NA>  1.02670523  1.054124e+00
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.mod3   33 104203 104483 -52069   104137                        
## H5.mod4   41 104199 104547 -52059   104117 20.312  8   0.009218 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod4<-getFE(H5.mod4))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.31
## 2                                                      age  0.00 0.00 33320.82
## 3                                                   gender  0.07 0.01 35527.23
## 4                                                     educ  0.01 0.00 35628.81
## 5                                                    resid -0.06 0.01 35604.24
## 6                            occupClerical support workers -0.04 0.09 35454.99
## 7                    occupCraft and related trades workers -0.06 0.09 35464.39
## 8                              occupElementary occupations  0.03 0.09 35459.54
## 9                                            occupManagers -0.02 0.09 35443.30
## 10                            occupOther: Not in paid work  0.15 0.09 35629.00
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35458.67
## 12                                      occupProfessionals  0.07 0.09 35453.25
## 13                                            occupRetired -0.02 0.10 35456.55
## 14                          occupService and sales workers -0.03 0.09 35452.88
## 15 occupSkilled agricultural, forestry and fishery workers -0.04 0.09 35474.75
## 16            occupTechnicians and associate professionals -0.03 0.09 35446.84
## 17                                         occupUnemployed -0.02 0.11 35488.30
## 18                                            environ.lvl1  0.12 0.01    19.65
## 19                                         engagement.lvl1  0.06 0.01    21.72
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01    27.73
##        t     p    LL    UL
## 1   0.46 0.648 -0.22  0.35
## 2   2.60 0.009  0.00  0.00
## 3   5.66 0.000  0.04  0.09
## 4   6.60 0.000  0.01  0.02
## 5  -4.95 0.000 -0.08 -0.04
## 6  -0.40 0.688 -0.21  0.14
## 7  -0.65 0.518 -0.23  0.12
## 8   0.30 0.767 -0.15  0.20
## 9  -0.17 0.863 -0.19  0.16
## 10  1.62 0.105 -0.03  0.32
## 11 -0.60 0.545 -0.23  0.12
## 12  0.82 0.413 -0.10  0.24
## 13 -0.25 0.801 -0.22  0.17
## 14 -0.32 0.747 -0.20  0.14
## 15 -0.48 0.628 -0.23  0.14
## 16 -0.31 0.755 -0.20  0.14
## 17 -0.18 0.855 -0.23  0.19
## 18 11.47 0.000  0.10  0.14
## 19  4.37 0.000  0.03  0.09
## 20  0.86 0.396 -0.01  0.03
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
## 1   0.30666977  0.0940463483
## 2   0.04316499  0.0018632163
## 3   0.06352955  0.0040360037
## 4   0.04537679  0.0020590531
## 5   0.08618317  0.0011408409
## 6   0.29709444  0.0057881699
## 7   0.30736238  0.0042771596
## 8   0.04133690  0.0001133562
## 9  -0.22779132 -0.0004461723
## 10 -0.17800013 -0.0005131329
## 11  0.49968388  0.2496839810
## 12  0.03883750  0.0015083511
## 13  0.04584691  0.0021019391
## 14  0.02252269  0.0005072716
## 15 -0.05677543 -0.0011018106
## 16  0.56550304  0.0129550876
## 17 -0.29825728 -0.0033566547
## 18  0.05592332  0.0000995759
## 19  0.90661876  0.0007930420
## 20 -0.37051502 -0.0003825922
## 21  1.02577908  1.0522227306
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
##                                                   0.29896 
##                     voting.group.environ.lvl1.(Intercept) 
##                                                   0.00363 
##                  voting.group.engagement.lvl1.(Intercept) 
##                                                   0.01840 
##     voting.group.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                   0.01360 
##                                 voting.group.environ.lvl1 
##                                                   0.04192 
##                 voting.group.engagement.lvl1.environ.lvl1 
##                                                   0.00098 
##    voting.group.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                  -0.01129 
##                              voting.group.engagement.lvl1 
##                                                   0.05913 
## voting.group.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                  -0.01229 
##                 voting.group.environ.lvl1:engagement.lvl1 
##                                                   0.03864 
##                                         cntry.(Intercept) 
##                                                   0.48713 
##                            cntry.environ.lvl1.(Intercept) 
##                                                  -0.00215 
##                         cntry.engagement.lvl1.(Intercept) 
##                                                   0.02528 
##            cntry.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                  -0.00655 
##                                        cntry.environ.lvl1 
##                                                   0.03780 
##                        cntry.engagement.lvl1.environ.lvl1 
##                                                   0.00394 
##           cntry.environ.lvl1:engagement.lvl1.environ.lvl1 
##                                                   0.01957 
##                                     cntry.engagement.lvl1 
##                                                   0.03665 
##        cntry.environ.lvl1:engagement.lvl1.engagement.lvl1 
##                                                  -0.00751 
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.mod4   41 104199 104547 -52059   104117                         
## H5.mod5   68 104053 104630 -51959   103917 199.98 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod5<-getFE(H5.mod5))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.48 0.14    61.87
## 2                                                      age  0.00 0.00 35562.27
## 3                                                   gender  0.07 0.01 35550.62
## 4                                                     educ  0.01 0.00 35668.50
## 5                                                    resid -0.06 0.01 35653.24
## 6                            occupClerical support workers -0.04 0.09 35509.60
## 7                    occupCraft and related trades workers -0.06 0.09 35519.06
## 8                              occupElementary occupations  0.02 0.09 35510.18
## 9                                            occupManagers -0.02 0.09 35496.94
## 10                            occupOther: Not in paid work  0.13 0.09 35568.62
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35515.24
## 12                                      occupProfessionals  0.07 0.09 35509.52
## 13                                            occupRetired -0.03 0.10 35499.27
## 14                          occupService and sales workers -0.03 0.09 35507.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35535.86
## 16            occupTechnicians and associate professionals -0.03 0.09 35502.94
## 17                                         occupUnemployed -0.03 0.11 35531.23
## 18                                            environ.lvl1  0.11 0.02   107.98
## 19                                         engagement.lvl1  0.01 0.03   142.46
## 20                            all.parties.lvl2Did not vote  0.44 0.06   182.83
## 21                              all.parties.lvl2Don't know  0.42 0.07   268.91
## 22                            all.parties.lvl2Invalid vote  0.48 0.35  1090.05
## 23                                  all.parties.lvl2NE age  0.72 0.07   273.86
## 24                              all.parties.lvl2NE citizen  0.86 0.08   268.92
## 25                                all.parties.lvl2NE other  0.71 0.10   657.74
## 26                               all.parties.lvl2No answer  0.54 0.38  1390.63
## 27                             all.parties.lvl2Other party  0.55 0.05   231.13
## 28                   all.parties.lvl2Pro-environment party  0.92 0.06   256.15
## 29                            environ.lvl1:engagement.lvl1  0.01 0.01    26.79
## 30               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    86.20
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   339.50
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25331.83
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.26
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   236.37
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1243.41
## 36                  environ.lvl1:all.parties.lvl2No answer  0.24 0.23  7550.15
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   134.21
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.02 0.03   195.88
## 39            engagement.lvl1:all.parties.lvl2Did not vote  0.03 0.03   105.78
## 40              engagement.lvl1:all.parties.lvl2Don't know -0.00 0.05   513.61
## 41            engagement.lvl1:all.parties.lvl2Invalid vote  0.11 0.45 18420.88
## 42                  engagement.lvl1:all.parties.lvl2NE age -0.01 0.05   354.94
## 43              engagement.lvl1:all.parties.lvl2NE citizen  0.01 0.05   258.49
## 44                engagement.lvl1:all.parties.lvl2NE other -0.02 0.09  2570.22
## 45               engagement.lvl1:all.parties.lvl2No answer -0.68 0.72 31757.51
## 46             engagement.lvl1:all.parties.lvl2Other party  0.07 0.03   171.19
## 47   engagement.lvl1:all.parties.lvl2Pro-environment party  0.10 0.04   296.05
##        t     p    LL    UL
## 1  -3.30 0.002 -0.77 -0.19
## 2   2.89 0.004  0.00  0.00
## 3   5.68 0.000  0.04  0.09
## 4   6.59 0.000  0.01  0.02
## 5  -4.88 0.000 -0.08 -0.03
## 6  -0.42 0.673 -0.21  0.14
## 7  -0.64 0.521 -0.23  0.12
## 8   0.26 0.797 -0.15  0.20
## 9  -0.19 0.850 -0.19  0.16
## 10  1.43 0.153 -0.05  0.31
## 11 -0.61 0.543 -0.23  0.12
## 12  0.80 0.425 -0.10  0.24
## 13 -0.27 0.786 -0.22  0.16
## 14 -0.33 0.743 -0.20  0.14
## 15 -0.51 0.613 -0.23  0.14
## 16 -0.32 0.753 -0.20  0.14
## 17 -0.24 0.812 -0.24  0.18
## 18  5.49 0.000  0.07  0.15
## 19  0.41 0.686 -0.04  0.07
## 20  6.98 0.000  0.32  0.57
## 21  6.00 0.000  0.29  0.56
## 22  1.35 0.176 -0.22  1.17
## 23 10.15 0.000  0.58  0.86
## 24 11.30 0.000  0.71  1.01
## 25  7.29 0.000  0.52  0.90
## 26  1.44 0.150 -0.20  1.28
## 27 11.42 0.000  0.46  0.65
## 28 14.13 0.000  0.79  1.04
## 29  0.63 0.535 -0.01  0.02
## 30  0.04 0.967 -0.04  0.05
## 31  0.91 0.361 -0.04  0.10
## 32  0.15 0.882 -0.70  0.82
## 33  1.34 0.183 -0.02  0.10
## 34 -1.06 0.290 -0.10  0.03
## 35  0.61 0.545 -0.08  0.14
## 36  1.06 0.290 -0.21  0.69
## 37  0.92 0.361 -0.02  0.06
## 38  0.80 0.423 -0.03  0.08
## 39  0.76 0.450 -0.04  0.09
## 40 -0.03 0.979 -0.10  0.10
## 41  0.25 0.799 -0.77  1.00
## 42 -0.13 0.900 -0.10  0.08
## 43  0.21 0.837 -0.09  0.11
## 44 -0.25 0.806 -0.20  0.16
## 45 -0.94 0.346 -2.09  0.73
## 46  2.46 0.015  0.01  0.13
## 47  2.33 0.020  0.02  0.19
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
## 1   0.19753217  3.901896e-02
## 2   0.04094707  1.676662e-03
## 3   0.05449208  2.969387e-03
## 4   0.04316977  1.863629e-03
## 5   0.06573845  5.317164e-04
## 6   0.15485894  1.666892e-03
## 7  -0.07651028 -6.524352e-04
## 8  -0.01103897 -2.463116e-05
## 9  -0.23564127 -4.165373e-04
## 10 -0.36840660 -8.666436e-04
## 11  0.48255877  2.328630e-01
## 12  0.03970632  1.576592e-03
## 13  0.04760744  2.266468e-03
## 14  0.02373098  5.631596e-04
## 15 -0.06057173 -1.160593e-03
## 16  0.51021906  1.172146e-02
## 17 -0.30364533 -3.477223e-03
## 18  0.10900597  2.060558e-04
## 19  0.89680996  8.450372e-04
## 20 -0.34142628 -3.857336e-04
## 21  1.02573145  1.052125e+00
```

\newpage

#### Look among which voting group there is strongest association between engagement and refugee attitudes


```r
H5.mod5.trends<-emtrends(H5.mod5,specs = c("all.parties.lvl2"),var=c("engagement.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 35740' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 35740)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 35740' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 35740)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.mod5.trends.tab<-data.frame(H5.mod5.trends))
```

```
##          all.parties.lvl2 engagement.lvl1.trend         SE  df   asymp.LCL
## 1  Anti-immigration party           0.011511487 0.02833794 Inf -0.04402985
## 2            Did not vote           0.036598204 0.02300658 Inf -0.00849387
## 3              Don't know           0.010135998 0.04581401 Inf -0.07965782
## 4            Invalid vote           0.126363928 0.45130400 Inf -0.75817565
## 5                  NE age           0.005759236 0.03915849 Inf -0.07098999
## 6              NE citizen           0.021501952 0.04214949 Inf -0.06110953
## 7                NE other          -0.010700005 0.08728812 Inf -0.18178157
## 8               No answer          -0.668262981 0.72111011 Inf -2.08161282
## 9             Other party           0.083001023 0.01669505 Inf  0.05027933
## 10  Pro-environment party           0.111997005 0.03582636 Inf  0.04177864
##     asymp.UCL
## 1  0.06705282
## 2  0.08169028
## 3  0.09992981
## 4  1.01090351
## 5  0.08250846
## 6  0.10411343
## 7  0.16038156
## 8  0.74508686
## 9  0.11572272
## 10 0.18221537
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
##                     group engagement.lvl1.trend   SE      p  adj.p asymp.LCL
## 1  Anti-immigration party                  0.01 0.03 0.6846 1.0000     -0.04
## 2            Did not vote                  0.04 0.02 0.1117 0.8933     -0.01
## 3              Don't know                  0.01 0.05 0.8249 1.0000     -0.08
## 4            Invalid vote                  0.13 0.45 0.7795 1.0000     -0.76
## 5                  NE age                  0.01 0.04 0.8831 1.0000     -0.07
## 6              NE citizen                  0.02 0.04 0.6100 1.0000     -0.06
## 7                NE other                 -0.01 0.09 0.9024 1.0000     -0.18
## 8               No answer                 -0.67 0.72 0.3541 1.0000     -2.08
## 9             Other party                  0.08 0.02 0.0000 0.0000      0.05
## 10  Pro-environment party                  0.11 0.04 0.0018 0.0159      0.04
##    asymp.UCL
## 1       0.07
## 2       0.08
## 3       0.10
## 4       1.01
## 5       0.08
## 6       0.10
## 7       0.16
## 8       0.75
## 9       0.12
## 10      0.18
```

```r
write.csv2(H5.mod5.trends.tab,"H5.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.1004855 0.04306793 Inf
##    z.ratio    p.value
## 1 2.333187 0.01963836
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect   0.0430 0.0989 Inf  0.435  0.6636 
##  Did not vote effect             0.0709 0.0976 Inf  0.726  0.4676 
##  Don't know effect               0.0415 0.1052 Inf  0.394  0.6933 
##  Invalid vote effect             0.1706 0.4584 Inf  0.372  0.7097 
##  NE age effect                   0.0366 0.1026 Inf  0.357  0.7210 
##  NE citizen effect               0.0541 0.1037 Inf  0.522  0.6016 
##  NE other effect                 0.0183 0.1286 Inf  0.143  0.8865 
##  No answer effect               -0.7123 0.7229 Inf -0.985  0.3245 
##  Other party effect              0.1225 0.0963 Inf  1.272  0.2034 
##  Pro-environment party effect    0.1547 0.1014 Inf  1.526  0.1270 
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
## 1           Other party - Anti-immigration party 0.07148954 0.02911709 Inf
## 2 Pro-environment party - Anti-immigration party 0.10048552 0.04306793 Inf
## 3            Pro-environment party - Other party 0.02899598 0.03638382 Inf
##     z.ratio    p.value
## 1 2.4552429 0.01407894
## 2 2.3331865 0.01963836
## 3 0.7969471 0.42548176
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   68 104053 104630 -51959   103917                       
## H5.mod6   77 104056 104710 -51951   103902 14.763  9    0.09765 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                                   Eff   Est
## 1                                                         (Intercept) -0.47
## 2                                                                 age  0.00
## 3                                                              gender  0.07
## 4                                                                educ  0.01
## 5                                                               resid -0.06
## 6                                       occupClerical support workers -0.04
## 7                               occupCraft and related trades workers -0.06
## 8                                         occupElementary occupations  0.02
## 9                                                       occupManagers -0.02
## 10                                       occupOther: Not in paid work  0.13
## 11                   occupPlant and machine operators, and assemblers -0.06
## 12                                                 occupProfessionals  0.07
## 13                                                       occupRetired -0.03
## 14                                     occupService and sales workers -0.03
## 15            occupSkilled agricultural, forestry and fishery workers -0.05
## 16                       occupTechnicians and associate professionals -0.03
## 17                                                    occupUnemployed -0.03
## 18                                                       environ.lvl1  0.11
## 19                                                    engagement.lvl1  0.02
## 20                                       all.parties.lvl2Did not vote  0.44
## 21                                         all.parties.lvl2Don't know  0.42
## 22                                       all.parties.lvl2Invalid vote  0.47
## 23                                             all.parties.lvl2NE age  0.72
## 24                                         all.parties.lvl2NE citizen  0.85
## 25                                           all.parties.lvl2NE other  0.71
## 26                                          all.parties.lvl2No answer  0.54
## 27                                        all.parties.lvl2Other party  0.55
## 28                              all.parties.lvl2Pro-environment party  0.91
## 29                                       environ.lvl1:engagement.lvl1 -0.04
## 30                          environ.lvl1:all.parties.lvl2Did not vote -0.00
## 31                            environ.lvl1:all.parties.lvl2Don't know  0.03
## 32                          environ.lvl1:all.parties.lvl2Invalid vote  0.06
## 33                                environ.lvl1:all.parties.lvl2NE age  0.04
## 34                            environ.lvl1:all.parties.lvl2NE citizen -0.04
## 35                              environ.lvl1:all.parties.lvl2NE other  0.04
## 36                             environ.lvl1:all.parties.lvl2No answer  0.24
## 37                           environ.lvl1:all.parties.lvl2Other party  0.02
## 38                 environ.lvl1:all.parties.lvl2Pro-environment party  0.02
## 39                       engagement.lvl1:all.parties.lvl2Did not vote  0.02
## 40                         engagement.lvl1:all.parties.lvl2Don't know -0.00
## 41                       engagement.lvl1:all.parties.lvl2Invalid vote  0.12
## 42                             engagement.lvl1:all.parties.lvl2NE age -0.01
## 43                         engagement.lvl1:all.parties.lvl2NE citizen -0.00
## 44                           engagement.lvl1:all.parties.lvl2NE other -0.01
## 45                          engagement.lvl1:all.parties.lvl2No answer -0.68
## 46                        engagement.lvl1:all.parties.lvl2Other party  0.07
## 47              engagement.lvl1:all.parties.lvl2Pro-environment party  0.10
## 48          environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote  0.03
## 49            environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know  0.04
## 50          environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote  0.34
## 51                environ.lvl1:engagement.lvl1:all.parties.lvl2NE age  0.02
## 52            environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen  0.13
## 53              environ.lvl1:engagement.lvl1:all.parties.lvl2NE other -0.04
## 54             environ.lvl1:engagement.lvl1:all.parties.lvl2No answer -0.08
## 55           environ.lvl1:engagement.lvl1:all.parties.lvl2Other party  0.05
## 56 environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party  0.05
##      SE       df     t     p    LL    UL
## 1  0.14    61.81 -3.26 0.002 -0.76 -0.18
## 2  0.00 35561.36  2.89 0.004  0.00  0.00
## 3  0.01 35553.22  5.72 0.000  0.04  0.09
## 4  0.00 35669.76  6.60 0.000  0.01  0.02
## 5  0.01 35654.28 -4.87 0.000 -0.08 -0.03
## 6  0.09 35512.26 -0.45 0.655 -0.21  0.13
## 7  0.09 35521.65 -0.67 0.503 -0.23  0.11
## 8  0.09 35512.45  0.23 0.816 -0.15  0.19
## 9  0.09 35499.42 -0.22 0.830 -0.19  0.15
## 10 0.09 35570.57  1.40 0.160 -0.05  0.31
## 11 0.09 35517.32 -0.64 0.521 -0.23  0.12
## 12 0.09 35511.77  0.76 0.444 -0.10  0.24
## 13 0.10 35501.08 -0.31 0.757 -0.22  0.16
## 14 0.09 35509.95 -0.36 0.721 -0.20  0.14
## 15 0.09 35539.17 -0.53 0.598 -0.23  0.13
## 16 0.09 35505.94 -0.34 0.731 -0.20  0.14
## 17 0.11 35532.88 -0.24 0.808 -0.24  0.18
## 18 0.02   107.96  5.60 0.000  0.07  0.15
## 19 0.03   142.26  0.53 0.597 -0.04  0.07
## 20 0.06   182.69  6.95 0.000  0.32  0.57
## 21 0.07   268.25  5.95 0.000  0.28  0.56
## 22 0.35  1087.13  1.34 0.180 -0.22  1.17
## 23 0.07   273.89 10.12 0.000  0.58  0.86
## 24 0.08   268.13 11.12 0.000  0.70  1.00
## 25 0.10   656.69  7.25 0.000  0.52  0.90
## 26 0.38  1390.54  1.44 0.149 -0.20  1.28
## 27 0.05   230.67 11.33 0.000  0.45  0.64
## 28 0.06   256.51 14.03 0.000  0.78  1.04
## 29 0.02   144.61 -1.70 0.091 -0.08  0.01
## 30 0.02    86.02 -0.01 0.991 -0.04  0.04
## 31 0.03   338.63  0.86 0.390 -0.04  0.10
## 32 0.39 25427.50  0.16 0.870 -0.70  0.83
## 33 0.03   237.66  1.31 0.193 -0.02  0.10
## 34 0.03   233.93 -1.31 0.191 -0.11  0.02
## 35 0.06  1246.82  0.67 0.505 -0.07  0.15
## 36 0.23  7561.32  1.05 0.292 -0.21  0.69
## 37 0.02   134.16  0.79 0.432 -0.02  0.05
## 38 0.03   195.40  0.74 0.463 -0.04  0.08
## 39 0.03   104.77  0.70 0.487 -0.04  0.09
## 40 0.05   502.12 -0.09 0.926 -0.11  0.10
## 41 0.45 18649.67  0.28 0.783 -0.76  1.01
## 42 0.05   348.95 -0.15 0.881 -0.10  0.08
## 43 0.05   254.20 -0.10 0.919 -0.10  0.09
## 44 0.09  2622.74 -0.12 0.901 -0.19  0.17
## 45 0.72 31703.95 -0.94 0.345 -2.10  0.73
## 46 0.03   168.71  2.29 0.023  0.01  0.12
## 47 0.04   291.53  2.25 0.025  0.01  0.18
## 48 0.03    89.15  1.25 0.216 -0.02  0.08
## 49 0.04   450.31  0.93 0.354 -0.04  0.12
## 50 0.55 31115.62  0.62 0.537 -0.73  1.41
## 51 0.04   338.20  0.65 0.515 -0.05  0.10
## 52 0.04   228.36  3.47 0.001  0.06  0.21
## 53 0.07  1447.94 -0.55 0.579 -0.17  0.10
## 54 0.54 30988.47 -0.14 0.889 -1.14  0.99
## 55 0.02   146.84  2.17 0.032  0.00  0.10
## 56 0.04   281.08  1.33 0.185 -0.02  0.12
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
## 1   0.19774809  3.910431e-02
## 2   0.04074286  1.659981e-03
## 3   0.05485458  3.009024e-03
## 4   0.03587722  1.287175e-03
## 5   0.05926903  4.775201e-04
## 6   0.15360858  1.666252e-03
## 7  -0.07280733 -5.165427e-04
## 8  -0.01718070 -3.839769e-05
## 9  -0.17335628 -2.534019e-04
## 10 -0.44021252 -8.663514e-04
## 11  0.48293468  2.332259e-01
## 12  0.03971972  1.577656e-03
## 13  0.04737523  2.244412e-03
## 14  0.02368416  5.609392e-04
## 15 -0.06217027 -1.192552e-03
## 16  0.51179009  1.170932e-02
## 17 -0.30965167 -3.541765e-03
## 18  0.11015106  2.072746e-04
## 19  0.91162995  8.575958e-04
## 20 -0.30542525 -3.427000e-04
## 21  1.02567279  1.052005e+00
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.mod5   68 104053 104630 -51959   103917                       
## H5.mod6   77 104056 104710 -51951   103902 14.763  9    0.09765 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.mod6<-getFE(H5.mod6))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.14    61.81
## 2                                                      age  0.00 0.00 35561.36
## 3                                                   gender  0.07 0.01 35553.22
## 4                                                     educ  0.01 0.00 35669.76
## 5                                                    resid -0.06 0.01 35654.28
## 6                            occupClerical support workers -0.04 0.09 35512.26
## 7                    occupCraft and related trades workers -0.06 0.09 35521.65
## 8                              occupElementary occupations  0.02 0.09 35512.44
## 9                                            occupManagers -0.02 0.09 35499.41
## 10                            occupOther: Not in paid work  0.13 0.09 35570.56
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.31
## 12                                      occupProfessionals  0.07 0.09 35511.76
## 13                                            occupRetired -0.03 0.10 35501.08
## 14                          occupService and sales workers -0.03 0.09 35509.95
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35539.17
## 16            occupTechnicians and associate professionals -0.03 0.09 35505.93
## 17                                         occupUnemployed -0.03 0.11 35532.88
## 18                                            environ.lvl1  0.11 0.02   107.96
## 19                                         engagement.lvl1  0.02 0.03   142.26
## 20                                             env.eng.int -0.04 0.02   144.61
## 21                            all.parties.lvl2Did not vote  0.44 0.06   182.69
## 22                              all.parties.lvl2Don't know  0.42 0.07   268.25
## 23                            all.parties.lvl2Invalid vote  0.47 0.35  1087.13
## 24                                  all.parties.lvl2NE age  0.72 0.07   273.89
## 25                              all.parties.lvl2NE citizen  0.85 0.08   268.13
## 26                                all.parties.lvl2NE other  0.71 0.10   656.69
## 27                               all.parties.lvl2No answer  0.54 0.38  1390.55
## 28                             all.parties.lvl2Other party  0.55 0.05   230.67
## 29                   all.parties.lvl2Pro-environment party  0.91 0.06   256.51
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    86.02
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   338.63
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25427.47
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   237.66
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   233.93
## 35                   environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1246.82
## 36                  environ.lvl1:all.parties.lvl2No answer  0.24 0.23  7561.30
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   134.16
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.02 0.03   195.40
## 39            engagement.lvl1:all.parties.lvl2Did not vote  0.02 0.03   104.77
## 40              engagement.lvl1:all.parties.lvl2Don't know -0.00 0.05   502.12
## 41            engagement.lvl1:all.parties.lvl2Invalid vote  0.12 0.45 18649.66
## 42                  engagement.lvl1:all.parties.lvl2NE age -0.01 0.05   348.95
## 43              engagement.lvl1:all.parties.lvl2NE citizen -0.00 0.05   254.20
## 44                engagement.lvl1:all.parties.lvl2NE other -0.01 0.09  2622.73
## 45               engagement.lvl1:all.parties.lvl2No answer -0.68 0.72 31703.95
## 46             engagement.lvl1:all.parties.lvl2Other party  0.07 0.03   168.71
## 47   engagement.lvl1:all.parties.lvl2Pro-environment party  0.10 0.04   291.53
## 48                env.eng.int:all.parties.lvl2Did not vote  0.03 0.03    89.15
## 49                  env.eng.int:all.parties.lvl2Don't know  0.04 0.04   450.31
## 50                env.eng.int:all.parties.lvl2Invalid vote  0.34 0.55 31115.63
## 51                      env.eng.int:all.parties.lvl2NE age  0.02 0.04   338.21
## 52                  env.eng.int:all.parties.lvl2NE citizen  0.13 0.04   228.36
## 53                    env.eng.int:all.parties.lvl2NE other -0.04 0.07  1447.94
## 54                   env.eng.int:all.parties.lvl2No answer -0.08 0.54 30988.48
## 55                 env.eng.int:all.parties.lvl2Other party  0.05 0.02   146.84
## 56       env.eng.int:all.parties.lvl2Pro-environment party  0.05 0.04   281.08
##        t     p    LL    UL
## 1  -3.26 0.002 -0.76 -0.18
## 2   2.89 0.004  0.00  0.00
## 3   5.72 0.000  0.04  0.09
## 4   6.60 0.000  0.01  0.02
## 5  -4.87 0.000 -0.08 -0.03
## 6  -0.45 0.655 -0.21  0.13
## 7  -0.67 0.503 -0.23  0.11
## 8   0.23 0.816 -0.15  0.19
## 9  -0.22 0.830 -0.19  0.15
## 10  1.40 0.160 -0.05  0.31
## 11 -0.64 0.521 -0.23  0.12
## 12  0.76 0.444 -0.10  0.24
## 13 -0.31 0.757 -0.22  0.16
## 14 -0.36 0.721 -0.20  0.14
## 15 -0.53 0.598 -0.23  0.13
## 16 -0.34 0.731 -0.20  0.14
## 17 -0.24 0.808 -0.24  0.18
## 18  5.60 0.000  0.07  0.15
## 19  0.53 0.597 -0.04  0.07
## 20 -1.70 0.091 -0.08  0.01
## 21  6.95 0.000  0.32  0.57
## 22  5.95 0.000  0.28  0.56
## 23  1.34 0.180 -0.22  1.17
## 24 10.12 0.000  0.58  0.86
## 25 11.12 0.000  0.70  1.00
## 26  7.25 0.000  0.52  0.90
## 27  1.44 0.149 -0.20  1.28
## 28 11.33 0.000  0.45  0.64
## 29 14.03 0.000  0.78  1.04
## 30 -0.01 0.991 -0.04  0.04
## 31  0.86 0.390 -0.04  0.10
## 32  0.16 0.870 -0.70  0.83
## 33  1.31 0.193 -0.02  0.10
## 34 -1.31 0.191 -0.11  0.02
## 35  0.67 0.505 -0.07  0.15
## 36  1.05 0.292 -0.21  0.69
## 37  0.79 0.432 -0.02  0.05
## 38  0.74 0.463 -0.04  0.08
## 39  0.70 0.487 -0.04  0.09
## 40 -0.09 0.926 -0.11  0.10
## 41  0.28 0.783 -0.76  1.01
## 42 -0.15 0.881 -0.10  0.08
## 43 -0.10 0.919 -0.10  0.09
## 44 -0.12 0.901 -0.19  0.17
## 45 -0.94 0.345 -2.10  0.73
## 46  2.29 0.023  0.01  0.12
## 47  2.25 0.025  0.01  0.18
## 48  1.25 0.216 -0.02  0.08
## 49  0.93 0.354 -0.04  0.12
## 50  0.62 0.537 -0.73  1.41
## 51  0.65 0.515 -0.05  0.10
## 52  3.47 0.001  0.06  0.21
## 53 -0.55 0.579 -0.17  0.10
## 54 -0.14 0.889 -1.14  0.99
## 55  2.17 0.032  0.00  0.10
## 56  1.33 0.185 -0.02  0.12
```

```r
(VC.H5.mod6<-getVC(H5.mod6))
```

```
##             grp            var1            var2      est_SD       est_SD2
## 1  voting.group     (Intercept)            <NA>  0.19774799  3.910427e-02
## 2  voting.group    environ.lvl1            <NA>  0.04074296  1.659989e-03
## 3  voting.group engagement.lvl1            <NA>  0.05485461  3.009028e-03
## 4  voting.group     env.eng.int            <NA>  0.03587717  1.287171e-03
## 5  voting.group     (Intercept)    environ.lvl1  0.05927151  4.775409e-04
## 6  voting.group     (Intercept) engagement.lvl1  0.15361039  1.666272e-03
## 7  voting.group     (Intercept)     env.eng.int -0.07280364 -5.165155e-04
## 8  voting.group    environ.lvl1 engagement.lvl1 -0.01717751 -3.839069e-05
## 9  voting.group    environ.lvl1     env.eng.int -0.17336687 -2.534176e-04
## 10 voting.group engagement.lvl1     env.eng.int -0.44021109 -8.663478e-04
## 11        cntry     (Intercept)            <NA>  0.48293273  2.332240e-01
## 12        cntry    environ.lvl1            <NA>  0.03971959  1.577646e-03
## 13        cntry engagement.lvl1            <NA>  0.04737519  2.244409e-03
## 14        cntry     env.eng.int            <NA>  0.02368430  5.609462e-04
## 15        cntry     (Intercept)    environ.lvl1 -0.06216367 -1.192417e-03
## 16        cntry     (Intercept) engagement.lvl1  0.51179504  1.170937e-02
## 17        cntry     (Intercept)     env.eng.int -0.30965651 -3.541828e-03
## 18        cntry    environ.lvl1 engagement.lvl1  0.11015536  2.072819e-04
## 19        cntry    environ.lvl1     env.eng.int  0.91162729  8.575960e-04
## 20        cntry engagement.lvl1     env.eng.int -0.30542645 -3.427033e-04
## 21     Residual            <NA>            <NA>  1.02567279  1.052005e+00
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.mod6.trends<-emtrends(H5.mod6,specs = c("all.parties.lvl2"),var=c("env.eng.int"))
(H5.mod6.trends.tab<-data.frame(H5.mod6.trends))
```

```
##          all.parties.lvl2 env.eng.int.trend         SE  df    asymp.LCL
## 1  Anti-immigration party      -0.036599505 0.02149166 Inf -0.078722375
## 2            Did not vote      -0.004134442 0.01648642 Inf -0.036447225
## 3              Don't know       0.003360580 0.03807472 Inf -0.071264505
## 4            Invalid vote       0.300733211 0.54566214 Inf -0.768744936
## 5                  NE age      -0.012011093 0.03192660 Inf -0.074586077
## 6              NE citizen       0.095857958 0.03246905 Inf  0.032219791
## 7                NE other      -0.074699048 0.06566682 Inf -0.203403642
## 8               No answer      -0.112330663 0.54089240 Inf -1.172460290
## 9             Other party       0.013466283 0.01126837 Inf -0.008619325
## 10  Pro-environment party       0.010419460 0.02905574 Inf -0.046528739
##      asymp.UCL
## 1  0.005523365
## 2  0.028178342
## 3  0.077985664
## 4  1.370211358
## 5  0.050563891
## 6  0.159496125
## 7  0.054005545
## 8  0.947798964
## 9  0.035551891
## 10 0.067367658
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
## 1  Anti-immigration party             -0.04 0.02 0.0886 0.7972     -0.08
## 2            Did not vote              0.00 0.02 0.8020 1.0000     -0.04
## 3              Don't know              0.00 0.04 0.9297 1.0000     -0.07
## 4            Invalid vote              0.30 0.55 0.5815 1.0000     -0.77
## 5                  NE age             -0.01 0.03 0.7068 1.0000     -0.07
## 6              NE citizen              0.10 0.03 0.0032 0.0315      0.03
## 7                NE other             -0.07 0.07 0.2553 1.0000     -0.20
## 8               No answer             -0.11 0.54 0.8355 1.0000     -1.17
## 9             Other party              0.01 0.01 0.2321 1.0000     -0.01
## 10  Pro-environment party              0.01 0.03 0.7199 1.0000     -0.05
##    asymp.UCL
## 1       0.01
## 2       0.03
## 3       0.08
## 4       1.37
## 5       0.05
## 6       0.16
## 7       0.05
## 8       0.95
## 9       0.04
## 10      0.07
```

```r
write.csv2(H5.mod6.trends.tab,"H5.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.04701896 0.03539014 Inf
##    z.ratio   p.value
## 1 1.328589 0.1839834
```

```r
#contrast for all groups against mean of other groups
contrast(H5.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect -0.06112 0.0885 Inf -0.691  1.0000 
##  Did not vote effect           -0.02505 0.0874 Inf -0.286  1.0000 
##  Don't know effect             -0.01672 0.0938 Inf -0.178  1.0000 
##  Invalid vote effect            0.31370 0.5490 Inf  0.571  1.0000 
##  NE age effect                 -0.03380 0.0915 Inf -0.369  1.0000 
##  NE citizen effect              0.08606 0.0917 Inf  0.938  1.0000 
##  NE other effect               -0.10345 0.1079 Inf -0.959  1.0000 
##  No answer effect              -0.14526 0.5444 Inf -0.267  1.0000 
##  Other party effect            -0.00549 0.0866 Inf -0.063  1.0000 
##  Pro-environment party effect  -0.00887 0.0906 Inf -0.098  1.0000 
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
##                                         contrast     estimate         SE  df
## 1           Other party - Anti-immigration party  0.050065788 0.02311458 Inf
## 2 Pro-environment party - Anti-immigration party  0.047018964 0.03539014 Inf
## 3            Pro-environment party - Other party -0.003046823 0.03025143 Inf
##      z.ratio    p.value
## 1  2.1659834 0.09093737
## 2  1.3285895 0.36796688
## 3 -0.1007167 0.91977537
```

