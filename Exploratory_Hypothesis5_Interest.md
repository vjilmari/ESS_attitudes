---
title: "Exlopratory Hypothesis 5 with Political Interest"
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

# Exploratory Hypothesis 5: The strength of the association between environment and refugee attitudes is stronger among more politically interested individuals.

Omit missing values on political interest item


```r
dat.H5.intr<-dat %>%
  filter(!is.na(polintr))
```

### Model 1: without interactions (only main effects)


```r
H5.exp.intr.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.37
## 2                                                      age  0.00 0.00 34011.21
## 3                                                   gender  0.07 0.01 35529.43
## 4                                                     educ  0.01 0.00 35634.84
## 5                                                    resid -0.06 0.01 35598.12
## 6                            occupClerical support workers -0.04 0.09 35472.39
## 7                    occupCraft and related trades workers -0.06 0.09 35482.31
## 8                              occupElementary occupations  0.02 0.09 35484.65
## 9                                            occupManagers -0.03 0.09 35471.30
## 10                            occupOther: Not in paid work  0.14 0.09 35638.43
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35480.18
## 12                                      occupProfessionals  0.07 0.09 35475.86
## 13                                            occupRetired -0.03 0.10 35473.15
## 14                          occupService and sales workers -0.03 0.09 35478.05
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35487.02
## 16            occupTechnicians and associate professionals -0.03 0.09 35471.93
## 17                                         occupUnemployed -0.01 0.11 35488.15
## 18                                            environ.lvl1  0.12 0.01    19.47
## 19                                            polintr.lvl1  0.06 0.01 35520.18
##        t     p    LL    UL
## 1   0.50 0.616 -0.21  0.36
## 2   3.27 0.001  0.00  0.00
## 3   5.93 0.000  0.05  0.09
## 4   6.12 0.000  0.01  0.02
## 5  -4.98 0.000 -0.08 -0.04
## 6  -0.47 0.638 -0.21  0.13
## 7  -0.69 0.488 -0.23  0.11
## 8   0.26 0.795 -0.15  0.20
## 9  -0.29 0.772 -0.20  0.15
## 10  1.58 0.114 -0.03  0.32
## 11 -0.68 0.497 -0.23  0.11
## 12  0.75 0.454 -0.11  0.24
## 13 -0.30 0.768 -0.22  0.16
## 14 -0.37 0.708 -0.20  0.14
## 15 -0.57 0.569 -0.24  0.13
## 16 -0.38 0.702 -0.20  0.14
## 17 -0.13 0.898 -0.22  0.20
## 18 11.55 0.000  0.10  0.14
## 19  8.76 0.000  0.05  0.08
```

```r
(VC.H5.exp.intr.mod1<-getVC(H5.exp.intr.mod1))
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30858043  0.095221882
## 2 voting.group environ.lvl1         <NA>  0.04247065  0.001803756
## 3 voting.group  (Intercept) environ.lvl1  0.09775183  0.001281098
## 4        cntry  (Intercept)         <NA>  0.49938803  0.249388401
## 5        cntry environ.lvl1         <NA>  0.03885202  0.001509479
## 6        cntry  (Intercept) environ.lvl1 -0.06809065 -0.001321111
## 7     Residual         <NA>         <NA>  1.02785526  1.056486439
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political polintr


```r
H5.exp.intr.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
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
## H5.exp.intr.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.intr.mod1:     age + gender + educ + resid + occup + environ.lvl1 + polintr.lvl1
## H5.exp.intr.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.intr.mod2:     age + gender + educ + resid + occup + environ.lvl1 + polintr.lvl1 + 
## H5.exp.intr.mod2:     environ.lvl1:polintr.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.exp.intr.mod1   26 104091 104311 -52019   104039                       
## H5.exp.intr.mod2   27 104087 104316 -52016   104033 6.1174  1    0.01339 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod2<-getFE(H5.exp.intr.mod2))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.40
## 2                                                      age  0.00 0.00 33996.51
## 3                                                   gender  0.07 0.01 35529.97
## 4                                                     educ  0.01 0.00 35634.15
## 5                                                    resid -0.06 0.01 35598.62
## 6                            occupClerical support workers -0.04 0.09 35472.67
## 7                    occupCraft and related trades workers -0.06 0.09 35482.63
## 8                              occupElementary occupations  0.02 0.09 35484.93
## 9                                            occupManagers -0.03 0.09 35471.55
## 10                            occupOther: Not in paid work  0.14 0.09 35638.81
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35480.48
## 12                                      occupProfessionals  0.07 0.09 35476.12
## 13                                            occupRetired -0.03 0.10 35473.40
## 14                          occupService and sales workers -0.03 0.09 35478.36
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35487.38
## 16            occupTechnicians and associate professionals -0.03 0.09 35472.25
## 17                                         occupUnemployed -0.01 0.11 35488.44
## 18                                            environ.lvl1  0.12 0.01    19.48
## 19                                            polintr.lvl1  0.06 0.01 35520.65
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01 35513.73
##        t     p    LL    UL
## 1   0.49 0.626 -0.22  0.36
## 2   3.24 0.001  0.00  0.00
## 3   5.95 0.000  0.05  0.09
## 4   6.12 0.000  0.01  0.02
## 5  -4.95 0.000 -0.08 -0.04
## 6  -0.46 0.645 -0.21  0.13
## 7  -0.68 0.495 -0.23  0.11
## 8   0.27 0.791 -0.15  0.20
## 9  -0.28 0.776 -0.20  0.15
## 10  1.59 0.111 -0.03  0.32
## 11 -0.67 0.505 -0.23  0.11
## 12  0.75 0.451 -0.11  0.24
## 13 -0.29 0.770 -0.22  0.16
## 14 -0.36 0.717 -0.20  0.14
## 15 -0.55 0.579 -0.23  0.13
## 16 -0.37 0.713 -0.20  0.14
## 17 -0.12 0.908 -0.22  0.20
## 18 11.51 0.000  0.10  0.14
## 19  8.81 0.000  0.05  0.08
## 20  2.47 0.013  0.00  0.02
```

```r
(VC.H5.exp.intr.mod2<-getVC(H5.exp.intr.mod2))
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30807908  0.094912721
## 2 voting.group environ.lvl1         <NA>  0.04244597  0.001801661
## 3 voting.group  (Intercept) environ.lvl1  0.09783966  0.001279421
## 4        cntry  (Intercept)         <NA>  0.49915606  0.249156769
## 5        cntry environ.lvl1         <NA>  0.03879048  0.001504702
## 6        cntry  (Intercept) environ.lvl1 -0.07185366 -0.001391267
## 7     Residual         <NA>         <NA>  1.02777890  1.056329472
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political polintr, allow polintr effect to vary between voting groups and countries


```r
H5.exp.intr.mod3<-lmer(refugees~(environ.lvl1+polintr.lvl1|voting.group)+
                (environ.lvl1+polintr.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polintr.lvl1+
                environ.lvl1:polintr.lvl1
                ,data=dat.H5.intr,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.intr.mod3)
```

```
## [1] FALSE
```

```r
anova(H5.exp.intr.mod2,H5.exp.intr.mod3)
```

```
## Data: dat.H5.intr
## Models:
## H5.exp.intr.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.intr.mod2:     age + gender + educ + resid + occup + environ.lvl1 + polintr.lvl1 + 
## H5.exp.intr.mod2:     environ.lvl1:polintr.lvl1
## H5.exp.intr.mod3: refugees ~ (environ.lvl1 + polintr.lvl1 | voting.group) + (environ.lvl1 + 
## H5.exp.intr.mod3:     polintr.lvl1 | cntry) + age + gender + educ + resid + occup + 
## H5.exp.intr.mod3:     environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.intr.mod2   27 104087 104316 -52016   104033                         
## H5.exp.intr.mod3   33 104035 104315 -51985   103969 63.403  6  9.134e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod3<-getFE(H5.exp.intr.mod3))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.39
## 2                                                      age  0.00 0.00 33995.67
## 3                                                   gender  0.07 0.01 35515.56
## 4                                                     educ  0.01 0.00 35576.48
## 5                                                    resid -0.06 0.01 35562.47
## 6                            occupClerical support workers -0.04 0.09 35436.79
## 7                    occupCraft and related trades workers -0.06 0.09 35444.83
## 8                              occupElementary occupations  0.03 0.09 35444.42
## 9                                            occupManagers -0.02 0.09 35438.98
## 10                            occupOther: Not in paid work  0.15 0.09 35605.45
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35448.35
## 12                                      occupProfessionals  0.07 0.09 35441.87
## 13                                            occupRetired -0.03 0.10 35432.14
## 14                          occupService and sales workers -0.03 0.09 35438.92
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35455.55
## 16            occupTechnicians and associate professionals -0.03 0.09 35436.62
## 17                                         occupUnemployed -0.01 0.11 35458.28
## 18                                            environ.lvl1  0.12 0.01    19.47
## 19                                            polintr.lvl1  0.07 0.01    20.77
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01 35490.32
##        t     p    LL    UL
## 1   0.48 0.632 -0.22  0.35
## 2   3.23 0.001  0.00  0.00
## 3   5.88 0.000  0.05  0.09
## 4   5.93 0.000  0.01  0.01
## 5  -4.97 0.000 -0.08 -0.04
## 6  -0.45 0.650 -0.21  0.13
## 7  -0.66 0.510 -0.23  0.11
## 8   0.28 0.776 -0.15  0.20
## 9  -0.24 0.811 -0.19  0.15
## 10  1.64 0.101 -0.03  0.33
## 11 -0.65 0.512 -0.23  0.12
## 12  0.76 0.448 -0.10  0.24
## 13 -0.26 0.794 -0.22  0.17
## 14 -0.32 0.749 -0.20  0.14
## 15 -0.52 0.600 -0.23  0.13
## 16 -0.35 0.725 -0.20  0.14
## 17 -0.09 0.926 -0.22  0.20
## 18 11.41 0.000  0.10  0.14
## 19  4.75 0.000  0.04  0.09
## 20  2.43 0.015  0.00  0.02
```

```r
(VC.H5.exp.intr.mod3<-getVC(H5.exp.intr.mod3))
```

```
##             grp         var1         var2      est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.30904222  0.0955070913
## 2  voting.group environ.lvl1         <NA>  0.04340784  0.0018842405
## 3  voting.group polintr.lvl1         <NA>  0.07126646  0.0050789082
## 4  voting.group  (Intercept) environ.lvl1  0.03839805  0.0005151042
## 5  voting.group  (Intercept) polintr.lvl1  0.45989618  0.0101289118
## 6  voting.group environ.lvl1 polintr.lvl1  0.05063037  0.0001566262
## 7         cntry  (Intercept)         <NA>  0.49869049  0.2486922002
## 8         cntry environ.lvl1         <NA>  0.03853467  0.0014849210
## 9         cntry polintr.lvl1         <NA>  0.04683406  0.0021934290
## 10        cntry  (Intercept) environ.lvl1 -0.09025503 -0.0017344197
## 11        cntry  (Intercept) polintr.lvl1  0.45309660  0.0105823878
## 12        cntry environ.lvl1 polintr.lvl1  0.02386250  0.0000430655
## 13     Residual         <NA>         <NA>  1.02554419  1.0517408823
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.exp.intr.mod4<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+
                   environ.lvl1:polintr.lvl1|voting.group)+
                (environ.lvl1+polintr.lvl1+
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
## H5.exp.intr.mod3: refugees ~ (environ.lvl1 + polintr.lvl1 | voting.group) + (environ.lvl1 + 
## H5.exp.intr.mod3:     polintr.lvl1 | cntry) + age + gender + educ + resid + occup + 
## H5.exp.intr.mod3:     environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1
## H5.exp.intr.mod4: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod4:     polintr.lvl1 + environ.lvl1:polintr.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.intr.mod3   33 104035 104315 -51985   103969                     
## H5.exp.intr.mod4   41 104043 104391 -51981   103961 7.8403  8     0.4492
```

```r
(FE.H5.exp.intr.mod4<-getFE(H5.exp.intr.mod4))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    49.40
## 2                                                      age  0.00 0.00 33970.60
## 3                                                   gender  0.07 0.01 35497.68
## 4                                                     educ  0.01 0.00 35544.12
## 5                                                    resid -0.06 0.01 35556.74
## 6                            occupClerical support workers -0.04 0.09 35421.40
## 7                    occupCraft and related trades workers -0.06 0.09 35432.01
## 8                              occupElementary occupations  0.02 0.09 35430.94
## 9                                            occupManagers -0.02 0.09 35417.72
## 10                            occupOther: Not in paid work  0.15 0.09 35600.69
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35430.70
## 12                                      occupProfessionals  0.06 0.09 35422.24
## 13                                            occupRetired -0.02 0.10 35422.55
## 14                          occupService and sales workers -0.03 0.09 35424.90
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35439.06
## 16            occupTechnicians and associate professionals -0.03 0.09 35418.90
## 17                                         occupUnemployed -0.01 0.11 35456.40
## 18                                            environ.lvl1  0.12 0.01    19.67
## 19                                            polintr.lvl1  0.07 0.01    20.83
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01    56.55
##        t     p    LL    UL
## 1   0.48 0.630 -0.22  0.35
## 2   3.25 0.001  0.00  0.00
## 3   5.88 0.000  0.05  0.09
## 4   5.94 0.000  0.01  0.01
## 5  -4.97 0.000 -0.08 -0.04
## 6  -0.47 0.642 -0.21  0.13
## 7  -0.67 0.502 -0.23  0.11
## 8   0.28 0.782 -0.15  0.20
## 9  -0.24 0.808 -0.19  0.15
## 10  1.64 0.102 -0.03  0.33
## 11 -0.66 0.506 -0.23  0.11
## 12  0.74 0.457 -0.11  0.24
## 13 -0.25 0.803 -0.22  0.17
## 14 -0.33 0.741 -0.20  0.14
## 15 -0.53 0.595 -0.23  0.13
## 16 -0.36 0.717 -0.20  0.14
## 17 -0.09 0.926 -0.22  0.20
## 18 11.49 0.000  0.10  0.14
## 19  4.75 0.000  0.04  0.09
## 20  1.97 0.054 -0.00  0.03
```

```r
(VC.H5.exp.intr.mod4<-getVC(H5.exp.intr.mod4))
```

```
##             grp                      var1                      var2
## 1  voting.group               (Intercept)                      <NA>
## 2  voting.group              environ.lvl1                      <NA>
## 3  voting.group              polintr.lvl1                      <NA>
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>
## 5  voting.group               (Intercept)              environ.lvl1
## 6  voting.group               (Intercept)              polintr.lvl1
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1
## 8  voting.group              environ.lvl1              polintr.lvl1
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1
## 11        cntry               (Intercept)                      <NA>
## 12        cntry              environ.lvl1                      <NA>
## 13        cntry              polintr.lvl1                      <NA>
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>
## 15        cntry               (Intercept)              environ.lvl1
## 16        cntry               (Intercept)              polintr.lvl1
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1
## 18        cntry              environ.lvl1              polintr.lvl1
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1
## 21     Residual                      <NA>                      <NA>
##          est_SD       est_SD2
## 1   0.308237181  9.501016e-02
## 2   0.043365299  1.880549e-03
## 3   0.071909474  5.170972e-03
## 4   0.036300898  1.317755e-03
## 5   0.040243268  5.379236e-04
## 6   0.455278655  1.009133e-02
## 7   0.156609524  1.752349e-03
## 8   0.047988227  1.496453e-04
## 9  -0.235402468 -3.705704e-04
## 10 -0.157191539 -4.103294e-04
## 11  0.498544030  2.485461e-01
## 12  0.038211064  1.460085e-03
## 13  0.046859579  2.195820e-03
## 14  0.009102557  8.285654e-05
## 15 -0.080477026 -1.533079e-03
## 16  0.456603483  1.066697e-02
## 17  0.037595397  1.706089e-04
## 18  0.039019656  6.986682e-05
## 19  0.955389570  3.323021e-04
## 20  0.331449878  1.413773e-04
## 21  1.024901814  1.050424e+00
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
##                                             0.30075 
##               voting.group.environ.lvl1.(Intercept) 
##                                             0.00170 
##               voting.group.polintr.lvl1.(Intercept) 
##                                             0.03194 
##  voting.group.environ.lvl1:polintr.lvl1.(Intercept) 
##                                             0.00555 
##                           voting.group.environ.lvl1 
##                                             0.04228 
##              voting.group.polintr.lvl1.environ.lvl1 
##                                             0.00208 
## voting.group.environ.lvl1:polintr.lvl1.environ.lvl1 
##                                            -0.00857 
##                           voting.group.polintr.lvl1 
##                                             0.06243 
## voting.group.environ.lvl1:polintr.lvl1.polintr.lvl1 
##                                            -0.00881 
##              voting.group.environ.lvl1:polintr.lvl1 
##                                             0.03275 
##                                   cntry.(Intercept) 
##                                             0.48643 
##                      cntry.environ.lvl1.(Intercept) 
##                                            -0.00300 
##                      cntry.polintr.lvl1.(Intercept) 
##                                             0.02088 
##         cntry.environ.lvl1:polintr.lvl1.(Intercept) 
##                                             0.00033 
##                                  cntry.environ.lvl1 
##                                             0.03716 
##                     cntry.polintr.lvl1.environ.lvl1 
##                                             0.00348 
##        cntry.environ.lvl1:polintr.lvl1.environ.lvl1 
##                                             0.00854 
##                                  cntry.polintr.lvl1 
##                                             0.04053 
##        cntry.environ.lvl1:polintr.lvl1.polintr.lvl1 
##                                             0.00242 
##                     cntry.environ.lvl1:polintr.lvl1 
##                                             0.00000
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and polintr


```r
H5.exp.intr.mod5<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|voting.group)+
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|cntry)+
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
## H5.exp.intr.mod4:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod4:     polintr.lvl1 + environ.lvl1:polintr.lvl1
## H5.exp.intr.mod5: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.intr.mod4   41 104043 104391 -51981   103961                         
## H5.exp.intr.mod5   68 103890 104467 -51877   103754 207.21 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.intr.mod5<-getFE(H5.exp.intr.mod5))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.14    61.98
## 2                                                      age  0.00 0.00 35616.41
## 3                                                   gender  0.07 0.01 35528.36
## 4                                                     educ  0.01 0.00 35592.31
## 5                                                    resid -0.06 0.01 35609.54
## 6                            occupClerical support workers -0.04 0.09 35481.31
## 7                    occupCraft and related trades workers -0.06 0.09 35494.04
## 8                              occupElementary occupations  0.02 0.09 35487.25
## 9                                            occupManagers -0.02 0.09 35477.54
## 10                            occupOther: Not in paid work  0.13 0.09 35545.97
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35494.29
## 12                                      occupProfessionals  0.06 0.09 35482.88
## 13                                            occupRetired -0.03 0.10 35469.24
## 14                          occupService and sales workers -0.03 0.09 35484.58
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35505.47
## 16            occupTechnicians and associate professionals -0.03 0.09 35480.34
## 17                                         occupUnemployed -0.02 0.11 35499.21
## 18                                            environ.lvl1  0.11 0.02   109.59
## 19                                            polintr.lvl1 -0.02 0.03   123.60
## 20                            all.parties.lvl2Did not vote  0.44 0.06   182.55
## 21                              all.parties.lvl2Don't know  0.42 0.07   267.91
## 22                            all.parties.lvl2Invalid vote  0.44 0.36  1125.55
## 23                                  all.parties.lvl2NE age  0.72 0.07   272.47
## 24                              all.parties.lvl2NE citizen  0.86 0.08   268.53
## 25                                all.parties.lvl2NE other  0.71 0.10   653.08
## 26                               all.parties.lvl2No answer  0.57 0.37  1352.13
## 27                             all.parties.lvl2Other party  0.55 0.05   230.27
## 28                   all.parties.lvl2Pro-environment party  0.92 0.06   256.09
## 29                               environ.lvl1:polintr.lvl1  0.01 0.01    54.50
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    87.38
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   340.06
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.05 0.39 24998.32
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.87
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   236.05
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1206.78
## 36                  environ.lvl1:all.parties.lvl2No answer  0.22 0.23  7135.17
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   135.62
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03   199.37
## 39               polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03    95.64
## 40                 polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05   426.05
## 41               polintr.lvl1:all.parties.lvl2Invalid vote  0.20 0.35 10239.75
## 42                     polintr.lvl1:all.parties.lvl2NE age  0.02 0.04   296.63
## 43                 polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04   213.48
## 44                   polintr.lvl1:all.parties.lvl2NE other  0.02 0.08  1261.73
## 45                  polintr.lvl1:all.parties.lvl2No answer -0.24 0.70 31717.96
## 46                polintr.lvl1:all.parties.lvl2Other party  0.10 0.03   141.88
## 47      polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04   250.77
##        t     p    LL    UL
## 1  -3.26 0.002 -0.76 -0.18
## 2   3.63 0.000  0.00  0.00
## 3   5.93 0.000  0.05  0.09
## 4   5.93 0.000  0.01  0.01
## 5  -4.86 0.000 -0.08 -0.03
## 6  -0.48 0.629 -0.22  0.13
## 7  -0.68 0.499 -0.23  0.11
## 8   0.23 0.818 -0.15  0.19
## 9  -0.27 0.789 -0.20  0.15
## 10  1.43 0.153 -0.05  0.31
## 11 -0.68 0.495 -0.23  0.11
## 12  0.71 0.476 -0.11  0.23
## 13 -0.28 0.783 -0.22  0.16
## 14 -0.34 0.731 -0.20  0.14
## 15 -0.56 0.573 -0.23  0.13
## 16 -0.38 0.706 -0.20  0.14
## 17 -0.16 0.872 -0.23  0.19
## 18  5.58 0.000  0.07  0.15
## 19 -0.85 0.399 -0.08  0.03
## 20  6.92 0.000  0.31  0.57
## 21  5.96 0.000  0.28  0.56
## 22  1.24 0.214 -0.26  1.15
## 23 10.15 0.000  0.58  0.86
## 24 11.26 0.000  0.71  1.01
## 25  7.32 0.000  0.52  0.91
## 26  1.52 0.128 -0.16  1.31
## 27 11.33 0.000  0.45  0.64
## 28 14.11 0.000  0.79  1.04
## 29  1.80 0.077 -0.00  0.03
## 30 -0.05 0.956 -0.05  0.04
## 31  0.75 0.456 -0.04  0.09
## 32  0.12 0.901 -0.71  0.81
## 33  1.24 0.217 -0.02  0.10
## 34 -1.20 0.231 -0.11  0.03
## 35  0.58 0.560 -0.08  0.14
## 36  0.96 0.335 -0.23  0.67
## 37  0.79 0.433 -0.02  0.05
## 38  0.30 0.766 -0.05  0.07
## 39  2.47 0.015  0.02  0.14
## 40  2.08 0.038  0.01  0.20
## 41  0.58 0.564 -0.48  0.89
## 42  0.55 0.585 -0.06  0.11
## 43  0.98 0.326 -0.04  0.13
## 44  0.23 0.820 -0.14  0.18
## 45 -0.34 0.735 -1.60  1.13
## 46  3.82 0.000  0.05  0.16
## 47  4.31 0.000  0.09  0.25
```

```r
(VC.H5.exp.intr.mod5<-getVC(H5.exp.intr.mod5))
```

```
##             grp                      var1                      var2
## 1  voting.group               (Intercept)                      <NA>
## 2  voting.group              environ.lvl1                      <NA>
## 3  voting.group              polintr.lvl1                      <NA>
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>
## 5  voting.group               (Intercept)              environ.lvl1
## 6  voting.group               (Intercept)              polintr.lvl1
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1
## 8  voting.group              environ.lvl1              polintr.lvl1
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1
## 11        cntry               (Intercept)                      <NA>
## 12        cntry              environ.lvl1                      <NA>
## 13        cntry              polintr.lvl1                      <NA>
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>
## 15        cntry               (Intercept)              environ.lvl1
## 16        cntry               (Intercept)              polintr.lvl1
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1
## 18        cntry              environ.lvl1              polintr.lvl1
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1
## 21     Residual                      <NA>                      <NA>
##          est_SD       est_SD2
## 1   0.197807826  3.912794e-02
## 2   0.041731515  1.741519e-03
## 3   0.055426253  3.072070e-03
## 4   0.035181558  1.237742e-03
## 5   0.017532777  1.447299e-04
## 6   0.276283508  3.029102e-03
## 7  -0.001200754 -8.356269e-06
## 8  -0.035818963 -8.285003e-05
## 9  -0.167315108 -2.456486e-04
## 10 -0.277168727 -5.404740e-04
## 11  0.481871385  2.322000e-01
## 12  0.039364305  1.549549e-03
## 13  0.049560668  2.456260e-03
## 14  0.009643158  9.299049e-05
## 15 -0.083856128 -1.590628e-03
## 16  0.401380310  9.585711e-03
## 17  0.022251889  1.033992e-04
## 18  0.118916991  2.319977e-04
## 19  0.974694665  3.699904e-04
## 20  0.337502509  1.612996e-04
## 21  1.024870011  1.050359e+00
```

\newpage

#### Look among which voting group there is strongest association between polintr and refugee attitudes


```r
H5.exp.intr.mod5.trends<-emtrends(H5.exp.intr.mod5,specs = c("all.parties.lvl2"),var=c("polintr.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 35702' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 35702)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 35702' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 35702)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.exp.intr.mod5.trends.tab<-data.frame(H5.exp.intr.mod5.trends))
```

```
##          all.parties.lvl2 polintr.lvl1.trend         SE  df    asymp.LCL
## 1  Anti-immigration party      -0.0227014845 0.02685978 Inf -0.075345689
## 2            Did not vote       0.0546378034 0.02260701 Inf  0.010328875
## 3              Don't know       0.0782603751 0.04349991 Inf -0.006997881
## 4            Invalid vote       0.1790907236 0.34924487 Inf -0.505416639
## 5                  NE age       0.0007101869 0.03698999 Inf -0.071788858
## 6              NE citizen       0.0211170667 0.03896882 Inf -0.055260415
## 7                NE other      -0.0039245266 0.07961813 Inf -0.159973185
## 8               No answer      -0.2589513708 0.69734588 Inf -1.625724177
## 9             Other party       0.0811356354 0.01628700 Inf  0.049213706
## 10  Pro-environment party       0.1511936282 0.03394109 Inf  0.084670320
##     asymp.UCL
## 1  0.02994272
## 2  0.09894673
## 3  0.16351863
## 4  0.86359809
## 5  0.07320923
## 6  0.09749455
## 7  0.15212413
## 8  1.10782144
## 9  0.11305756
## 10 0.21771694
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
## 1  Anti-immigration party              -0.02 0.03 0.3980 1.0000     -0.08
## 2            Did not vote               0.05 0.02 0.0157 0.1252      0.01
## 3              Don't know               0.08 0.04 0.0720 0.5040     -0.01
## 4            Invalid vote               0.18 0.35 0.6081 1.0000     -0.51
## 5                  NE age               0.00 0.04 0.9847 1.0000     -0.07
## 6              NE citizen               0.02 0.04 0.5879 1.0000     -0.06
## 7                NE other               0.00 0.08 0.9607 1.0000     -0.16
## 8               No answer              -0.26 0.70 0.7104 1.0000     -1.63
## 9             Other party               0.08 0.02 0.0000 0.0000      0.05
## 10  Pro-environment party               0.15 0.03 0.0000 0.0001      0.08
##    asymp.UCL
## 1       0.03
## 2       0.10
## 3       0.16
## 4       0.86
## 5       0.07
## 6       0.10
## 7       0.15
## 8       1.11
## 9       0.11
## 10      0.22
```

```r
write.csv2(H5.exp.intr.mod5.trends.tab,"H5.exp.intr.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.1738951 0.04032848 Inf
##    z.ratio     p.value
## 1 4.311968 1.61808e-05
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.intr.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect -0.05640 0.0908 Inf -0.621  0.5345 
##  Did not vote effect            0.02953 0.0897 Inf  0.329  0.7419 
##  Don't know effect              0.05578 0.0970 Inf  0.575  0.5652 
##  Invalid vote effect            0.16782 0.3578 Inf  0.469  0.6390 
##  NE age effect                 -0.03039 0.0943 Inf -0.322  0.7473 
##  NE citizen effect             -0.00771 0.0950 Inf -0.081  0.9353 
##  NE other effect               -0.03553 0.1175 Inf -0.302  0.7623 
##  No answer effect              -0.31890 0.6985 Inf -0.457  0.6480 
##  Other party effect             0.05898 0.0883 Inf  0.668  0.5041 
##  Pro-environment party effect   0.13682 0.0932 Inf  1.468  0.1420 
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
## 1           Other party - Anti-immigration party 0.10383712 0.02716097 Inf
## 2 Pro-environment party - Anti-immigration party 0.17389511 0.04032848 Inf
## 3            Pro-environment party - Other party 0.07005799 0.03409876 Inf
##    z.ratio     p.value
## 1 3.823028 1.31823e-04
## 2 4.311968 1.61808e-05
## 3 2.054561 3.99214e-02
```

\newpage

### Model 6: Enter three-way interaction voting group x polintr x environment attitudes


```r
H5.exp.intr.mod6<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|voting.group)+
                (environ.lvl1+polintr.lvl1+environ.lvl1:polintr.lvl1|cntry)+
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
## H5.exp.intr.mod5:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
## H5.exp.intr.mod6: refugees ~ (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod6:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod6:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod6:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2 + 
## H5.exp.intr.mod6:     environ.lvl1:polintr.lvl1:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.intr.mod5   68 103890 104467 -51877   103754                     
## H5.exp.intr.mod6   77 103897 104550 -51872   103743 11.091  9     0.2695
```

```r
(FE.H5.exp.intr.mod6<-getFE(H5.exp.intr.mod6))
```

```
##                                                                Eff   Est   SE
## 1                                                      (Intercept) -0.47 0.14
## 2                                                              age  0.00 0.00
## 3                                                           gender  0.07 0.01
## 4                                                             educ  0.01 0.00
## 5                                                            resid -0.06 0.01
## 6                                    occupClerical support workers -0.04 0.09
## 7                            occupCraft and related trades workers -0.06 0.09
## 8                                      occupElementary occupations  0.02 0.09
## 9                                                    occupManagers -0.03 0.09
## 10                                    occupOther: Not in paid work  0.13 0.09
## 11                occupPlant and machine operators, and assemblers -0.06 0.09
## 12                                              occupProfessionals  0.06 0.09
## 13                                                    occupRetired -0.03 0.10
## 14                                  occupService and sales workers -0.03 0.09
## 15         occupSkilled agricultural, forestry and fishery workers -0.05 0.09
## 16                    occupTechnicians and associate professionals -0.03 0.09
## 17                                                 occupUnemployed -0.02 0.11
## 18                                                    environ.lvl1  0.11 0.02
## 19                                                    polintr.lvl1 -0.02 0.03
## 20                                    all.parties.lvl2Did not vote  0.44 0.06
## 21                                      all.parties.lvl2Don't know  0.42 0.07
## 22                                    all.parties.lvl2Invalid vote  0.44 0.36
## 23                                          all.parties.lvl2NE age  0.72 0.07
## 24                                      all.parties.lvl2NE citizen  0.85 0.08
## 25                                        all.parties.lvl2NE other  0.71 0.10
## 26                                       all.parties.lvl2No answer  0.58 0.37
## 27                                     all.parties.lvl2Other party  0.55 0.05
## 28                           all.parties.lvl2Pro-environment party  0.92 0.07
## 29                                       environ.lvl1:polintr.lvl1 -0.01 0.02
## 30                       environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02
## 31                         environ.lvl1:all.parties.lvl2Don't know  0.02 0.03
## 32                       environ.lvl1:all.parties.lvl2Invalid vote  0.04 0.39
## 33                             environ.lvl1:all.parties.lvl2NE age  0.04 0.03
## 34                         environ.lvl1:all.parties.lvl2NE citizen -0.05 0.03
## 35                           environ.lvl1:all.parties.lvl2NE other  0.04 0.06
## 36                          environ.lvl1:all.parties.lvl2No answer  0.17 0.26
## 37                        environ.lvl1:all.parties.lvl2Other party  0.01 0.02
## 38              environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03
## 39                       polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03
## 40                         polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05
## 41                       polintr.lvl1:all.parties.lvl2Invalid vote  0.21 0.35
## 42                             polintr.lvl1:all.parties.lvl2NE age  0.02 0.04
## 43                         polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04
## 44                           polintr.lvl1:all.parties.lvl2NE other  0.02 0.08
## 45                          polintr.lvl1:all.parties.lvl2No answer -0.26 0.70
## 46                        polintr.lvl1:all.parties.lvl2Other party  0.10 0.03
## 47              polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04
## 48          environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote  0.02 0.02
## 49            environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know  0.03 0.04
## 50          environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote  0.14 0.41
## 51                environ.lvl1:polintr.lvl1:all.parties.lvl2NE age  0.01 0.03
## 52            environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen  0.09 0.03
## 53              environ.lvl1:polintr.lvl1:all.parties.lvl2NE other -0.08 0.06
## 54             environ.lvl1:polintr.lvl1:all.parties.lvl2No answer  0.36 0.76
## 55           environ.lvl1:polintr.lvl1:all.parties.lvl2Other party  0.02 0.02
## 56 environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party -0.00 0.03
##          df     t     p    LL    UL
## 1     62.03 -3.25 0.002 -0.76 -0.18
## 2  35620.83  3.62 0.000  0.00  0.00
## 3  35529.71  5.95 0.000  0.05  0.09
## 4  35596.16  5.91 0.000  0.01  0.01
## 5  35609.27 -4.88 0.000 -0.08 -0.03
## 6  35481.74 -0.50 0.614 -0.22  0.13
## 7  35494.16 -0.69 0.488 -0.23  0.11
## 8  35487.61  0.21 0.832 -0.15  0.19
## 9  35477.83 -0.29 0.774 -0.20  0.15
## 10 35545.94  1.41 0.157 -0.05  0.31
## 11 35494.22 -0.70 0.482 -0.24  0.11
## 12 35482.92  0.69 0.490 -0.11  0.23
## 13 35469.46 -0.30 0.761 -0.22  0.16
## 14 35484.89 -0.37 0.715 -0.20  0.14
## 15 35505.84 -0.58 0.563 -0.24  0.13
## 16 35480.69 -0.39 0.693 -0.21  0.14
## 17 35500.20 -0.16 0.872 -0.23  0.19
## 18   109.85  5.63 0.000  0.07  0.15
## 19   122.70 -0.83 0.410 -0.08  0.03
## 20   182.54  6.91 0.000  0.31  0.57
## 21   267.96  5.93 0.000  0.28  0.56
## 22  1126.05  1.23 0.217 -0.26  1.14
## 23   272.67 10.14 0.000  0.58  0.86
## 24   268.30 11.14 0.000  0.70  1.00
## 25   652.83  7.31 0.000  0.52  0.91
## 26  1357.83  1.54 0.123 -0.16  1.31
## 27   230.38 11.31 0.000  0.45  0.64
## 28   257.35 14.11 0.000  0.79  1.05
## 29   126.69 -0.38 0.708 -0.04  0.03
## 30    87.78 -0.09 0.926 -0.05  0.04
## 31   342.47  0.70 0.484 -0.04  0.09
## 32 25264.47  0.10 0.924 -0.73  0.80
## 33   238.97  1.22 0.223 -0.02  0.10
## 34   236.46 -1.39 0.166 -0.11  0.02
## 35  1205.32  0.63 0.530 -0.07  0.14
## 36 10536.40  0.65 0.514 -0.34  0.67
## 37   136.75  0.72 0.471 -0.02  0.05
## 38   199.54  0.26 0.797 -0.05  0.06
## 39    94.24  2.46 0.016  0.01  0.14
## 40   421.96  2.04 0.042  0.00  0.20
## 41 10518.93  0.60 0.550 -0.48  0.90
## 42   293.07  0.56 0.578 -0.06  0.11
## 43   211.95  0.80 0.423 -0.05  0.12
## 44  1259.27  0.20 0.844 -0.15  0.18
## 45 31722.95 -0.38 0.707 -1.64  1.11
## 46   140.44  3.80 0.000  0.05  0.16
## 47   249.54  4.30 0.000  0.09  0.25
## 48    85.68  0.71 0.477 -0.03  0.06
## 49   448.49  0.86 0.389 -0.04  0.11
## 50 27792.84  0.34 0.736 -0.67  0.95
## 51   364.25  0.22 0.824 -0.06  0.08
## 52   215.95  2.56 0.011  0.02  0.15
## 53  1155.29 -1.37 0.172 -0.21  0.04
## 54 34661.69  0.47 0.637 -1.14  1.86
## 55   124.68  1.10 0.272 -0.02  0.06
## 56   258.63 -0.06 0.956 -0.06  0.06
```

```r
(VC.H5.exp.intr.mod6<-getVC(H5.exp.intr.mod6))
```

```
##             grp                      var1                      var2
## 1  voting.group               (Intercept)                      <NA>
## 2  voting.group              environ.lvl1                      <NA>
## 3  voting.group              polintr.lvl1                      <NA>
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>
## 5  voting.group               (Intercept)              environ.lvl1
## 6  voting.group               (Intercept)              polintr.lvl1
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1
## 8  voting.group              environ.lvl1              polintr.lvl1
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1
## 11        cntry               (Intercept)                      <NA>
## 12        cntry              environ.lvl1                      <NA>
## 13        cntry              polintr.lvl1                      <NA>
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>
## 15        cntry               (Intercept)              environ.lvl1
## 16        cntry               (Intercept)              polintr.lvl1
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1
## 18        cntry              environ.lvl1              polintr.lvl1
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1
## 21     Residual                      <NA>                      <NA>
##          est_SD       est_SD2
## 1   0.198111360  3.924811e-02
## 2   0.041542864  1.725810e-03
## 3   0.055319665  3.060265e-03
## 4   0.030462894  9.279879e-04
## 5   0.009772025  8.042488e-05
## 6   0.278234232  3.049295e-03
## 7  -0.005343558 -3.224861e-05
## 8  -0.038065623 -8.748003e-05
## 9  -0.118625456 -1.501224e-04
## 10 -0.301971679 -5.088818e-04
## 11  0.481662638  2.319989e-01
## 12  0.039339237  1.547576e-03
## 13  0.049623564  2.462498e-03
## 14  0.009958964  9.918096e-05
## 15 -0.084611697 -1.603243e-03
## 16  0.403867631  9.653170e-03
## 17  0.040551719  1.945210e-04
## 18  0.128007153  2.498896e-04
## 19  0.966486879  3.786483e-04
## 20  0.377774008  1.866956e-04
## 21  1.024822113  1.050260e+00
```

#### Refit with manually coded level-1 interaction


```r
dat.H5.intr$env.intr.int<-dat.H5.intr$environ.lvl1*dat.H5.intr$polintr.lvl1

H5.exp.intr.mod6<-lmer(refugees~
                (environ.lvl1+polintr.lvl1+env.intr.int|voting.group)+
                (environ.lvl1+polintr.lvl1+env.intr.int|cntry)+
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
## H5.exp.intr.mod5:     voting.group) + (environ.lvl1 + polintr.lvl1 + environ.lvl1:polintr.lvl1 | 
## H5.exp.intr.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.intr.mod5:     polintr.lvl1 + environ.lvl1:polintr.lvl1 + all.parties.lvl2 + 
## H5.exp.intr.mod5:     environ.lvl1:all.parties.lvl2 + polintr.lvl1:all.parties.lvl2
## H5.exp.intr.mod6: refugees ~ (environ.lvl1 + polintr.lvl1 + env.intr.int | voting.group) + 
## H5.exp.intr.mod6:     (environ.lvl1 + polintr.lvl1 + env.intr.int | cntry) + age + 
## H5.exp.intr.mod6:     gender + educ + resid + occup + environ.lvl1 + polintr.lvl1 + 
## H5.exp.intr.mod6:     env.intr.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.exp.intr.mod6:     polintr.lvl1:all.parties.lvl2 + env.intr.int:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.intr.mod5   68 103890 104467 -51877   103754                     
## H5.exp.intr.mod6   77 103897 104550 -51872   103743 11.091  9     0.2695
```

```r
(FE.H5.exp.intr.mod6<-getFE(H5.exp.intr.mod6))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.14    62.03
## 2                                                      age  0.00 0.00 35620.81
## 3                                                   gender  0.07 0.01 35529.71
## 4                                                     educ  0.01 0.00 35596.15
## 5                                                    resid -0.06 0.01 35609.28
## 6                            occupClerical support workers -0.04 0.09 35481.74
## 7                    occupCraft and related trades workers -0.06 0.09 35494.16
## 8                              occupElementary occupations  0.02 0.09 35487.61
## 9                                            occupManagers -0.03 0.09 35477.83
## 10                            occupOther: Not in paid work  0.13 0.09 35545.93
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35494.22
## 12                                      occupProfessionals  0.06 0.09 35482.93
## 13                                            occupRetired -0.03 0.10 35469.46
## 14                          occupService and sales workers -0.03 0.09 35484.89
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35505.84
## 16            occupTechnicians and associate professionals -0.03 0.09 35480.69
## 17                                         occupUnemployed -0.02 0.11 35500.20
## 18                                            environ.lvl1  0.11 0.02   109.85
## 19                                            polintr.lvl1 -0.02 0.03   122.70
## 20                                            env.intr.int -0.01 0.02   126.69
## 21                            all.parties.lvl2Did not vote  0.44 0.06   182.53
## 22                              all.parties.lvl2Don't know  0.42 0.07   267.95
## 23                            all.parties.lvl2Invalid vote  0.44 0.36  1126.03
## 24                                  all.parties.lvl2NE age  0.72 0.07   272.66
## 25                              all.parties.lvl2NE citizen  0.85 0.08   268.29
## 26                                all.parties.lvl2NE other  0.71 0.10   652.82
## 27                               all.parties.lvl2No answer  0.58 0.37  1357.80
## 28                             all.parties.lvl2Other party  0.55 0.05   230.37
## 29                   all.parties.lvl2Pro-environment party  0.92 0.07   257.35
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    87.78
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.02 0.03   342.46
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.04 0.39 25264.46
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.97
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.05 0.03   236.46
## 35                   environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1205.32
## 36                  environ.lvl1:all.parties.lvl2No answer  0.17 0.26 10536.39
## 37                environ.lvl1:all.parties.lvl2Other party  0.01 0.02   136.75
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03   199.54
## 39               polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03    94.24
## 40                 polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05   421.95
## 41               polintr.lvl1:all.parties.lvl2Invalid vote  0.21 0.35 10518.76
## 42                     polintr.lvl1:all.parties.lvl2NE age  0.02 0.04   293.07
## 43                 polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04   211.95
## 44                   polintr.lvl1:all.parties.lvl2NE other  0.02 0.08  1259.25
## 45                  polintr.lvl1:all.parties.lvl2No answer -0.26 0.70 31722.87
## 46                polintr.lvl1:all.parties.lvl2Other party  0.10 0.03   140.44
## 47      polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04   249.54
## 48               env.intr.int:all.parties.lvl2Did not vote  0.02 0.02    85.68
## 49                 env.intr.int:all.parties.lvl2Don't know  0.03 0.04   448.48
## 50               env.intr.int:all.parties.lvl2Invalid vote  0.14 0.41 27792.75
## 51                     env.intr.int:all.parties.lvl2NE age  0.01 0.03   364.24
## 52                 env.intr.int:all.parties.lvl2NE citizen  0.09 0.03   215.95
## 53                   env.intr.int:all.parties.lvl2NE other -0.08 0.06  1155.27
## 54                  env.intr.int:all.parties.lvl2No answer  0.36 0.76 34661.68
## 55                env.intr.int:all.parties.lvl2Other party  0.02 0.02   124.68
## 56      env.intr.int:all.parties.lvl2Pro-environment party -0.00 0.03   258.62
##        t     p    LL    UL
## 1  -3.25 0.002 -0.76 -0.18
## 2   3.62 0.000  0.00  0.00
## 3   5.95 0.000  0.05  0.09
## 4   5.91 0.000  0.01  0.01
## 5  -4.88 0.000 -0.08 -0.03
## 6  -0.50 0.614 -0.22  0.13
## 7  -0.69 0.488 -0.23  0.11
## 8   0.21 0.832 -0.15  0.19
## 9  -0.29 0.774 -0.20  0.15
## 10  1.41 0.157 -0.05  0.31
## 11 -0.70 0.482 -0.24  0.11
## 12  0.69 0.490 -0.11  0.23
## 13 -0.30 0.761 -0.22  0.16
## 14 -0.37 0.715 -0.20  0.14
## 15 -0.58 0.563 -0.24  0.13
## 16 -0.39 0.693 -0.21  0.14
## 17 -0.16 0.872 -0.23  0.19
## 18  5.63 0.000  0.07  0.15
## 19 -0.83 0.410 -0.08  0.03
## 20 -0.38 0.708 -0.04  0.03
## 21  6.91 0.000  0.31  0.57
## 22  5.93 0.000  0.28  0.56
## 23  1.23 0.217 -0.26  1.14
## 24 10.14 0.000  0.58  0.86
## 25 11.14 0.000  0.70  1.00
## 26  7.31 0.000  0.52  0.91
## 27  1.54 0.123 -0.16  1.31
## 28 11.31 0.000  0.45  0.64
## 29 14.11 0.000  0.79  1.05
## 30 -0.09 0.926 -0.05  0.04
## 31  0.70 0.484 -0.04  0.09
## 32  0.10 0.924 -0.73  0.80
## 33  1.22 0.223 -0.02  0.10
## 34 -1.39 0.166 -0.11  0.02
## 35  0.63 0.530 -0.07  0.14
## 36  0.65 0.514 -0.34  0.67
## 37  0.72 0.471 -0.02  0.05
## 38  0.26 0.797 -0.05  0.06
## 39  2.46 0.016  0.01  0.14
## 40  2.04 0.042  0.00  0.20
## 41  0.60 0.550 -0.48  0.90
## 42  0.56 0.578 -0.06  0.11
## 43  0.80 0.423 -0.05  0.12
## 44  0.20 0.844 -0.15  0.18
## 45 -0.38 0.707 -1.64  1.11
## 46  3.80 0.000  0.05  0.16
## 47  4.30 0.000  0.09  0.25
## 48  0.71 0.477 -0.03  0.06
## 49  0.86 0.389 -0.04  0.11
## 50  0.34 0.736 -0.67  0.95
## 51  0.22 0.824 -0.06  0.08
## 52  2.56 0.011  0.02  0.15
## 53 -1.37 0.172 -0.21  0.04
## 54  0.47 0.637 -1.14  1.86
## 55  1.10 0.272 -0.02  0.06
## 56 -0.06 0.956 -0.06  0.06
```

```r
(VC.H5.exp.intr.mod6<-getVC(H5.exp.intr.mod6))
```

```
##             grp         var1         var2       est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.198111517  3.924817e-02
## 2  voting.group environ.lvl1         <NA>  0.041542908  1.725813e-03
## 3  voting.group polintr.lvl1         <NA>  0.055320255  3.060331e-03
## 4  voting.group env.intr.int         <NA>  0.030463074  9.279989e-04
## 5  voting.group  (Intercept) environ.lvl1  0.009775029  8.044974e-05
## 6  voting.group  (Intercept) polintr.lvl1  0.278234076  3.049328e-03
## 7  voting.group  (Intercept) env.intr.int -0.005342924 -3.224500e-05
## 8  voting.group environ.lvl1 polintr.lvl1 -0.038071087 -8.749361e-05
## 9  voting.group environ.lvl1 env.intr.int -0.118634496 -1.501349e-04
## 10 voting.group polintr.lvl1 env.intr.int -0.301982633 -5.089087e-04
## 11        cntry  (Intercept)         <NA>  0.481658528  2.319949e-01
## 12        cntry environ.lvl1         <NA>  0.039339109  1.547565e-03
## 13        cntry polintr.lvl1         <NA>  0.049623450  2.462487e-03
## 14        cntry env.intr.int         <NA>  0.009958889  9.917946e-05
## 15        cntry  (Intercept) environ.lvl1 -0.084609166 -1.603176e-03
## 16        cntry  (Intercept) polintr.lvl1  0.403865605  9.653017e-03
## 17        cntry  (Intercept) env.intr.int  0.040557267  1.945444e-04
## 18        cntry environ.lvl1 polintr.lvl1  0.128000767  2.498757e-04
## 19        cntry environ.lvl1 env.intr.int  0.966480981  3.786419e-04
## 20        cntry polintr.lvl1 env.intr.int  0.377789797  1.867016e-04
## 21     Residual         <NA>         <NA>  1.024822098  1.050260e+00
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.exp.intr.mod6.trends<-emtrends(H5.exp.intr.mod6,specs = c("all.parties.lvl2"),var=c("env.intr.int"))
(H5.exp.intr.mod6.trends.tab<-data.frame(H5.exp.intr.mod6.trends))
```

```
##          all.parties.lvl2 env.intr.int.trend          SE  df    asymp.LCL
## 1  Anti-immigration party      -0.0069621407 0.018544266 Inf -0.043308234
## 2            Did not vote       0.0095517468 0.014187703 Inf -0.018255640
## 3              Don't know       0.0267672477 0.034578183 Inf -0.041004746
## 4            Invalid vote       0.1320711625 0.412682858 Inf -0.676772377
## 5                  NE age       0.0007638314 0.029563838 Inf -0.057180227
## 6              NE citizen       0.0795677618 0.028392962 Inf  0.023918580
## 7                NE other      -0.0916793524 0.059222250 Inf -0.207752829
## 8               No answer       0.3529258125 0.762887333 Inf -1.142305884
## 9             Other party       0.0154857242 0.008953365 Inf -0.002062548
## 10  Pro-environment party      -0.0086996200 0.025538314 Inf -0.058753796
##     asymp.UCL
## 1  0.02938395
## 2  0.03735913
## 3  0.09453924
## 4  0.94091470
## 5  0.05870789
## 6  0.13521694
## 7  0.02439412
## 8  1.84815751
## 9  0.03303400
## 10 0.04135456
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
## 1  Anti-immigration party              -0.01 0.02 0.7073 1.0000     -0.04
## 2            Did not vote               0.01 0.01 0.5008 1.0000     -0.02
## 3              Don't know               0.03 0.03 0.4389 1.0000     -0.04
## 4            Invalid vote               0.13 0.41 0.7489 1.0000     -0.68
## 5                  NE age               0.00 0.03 0.9794 1.0000     -0.06
## 6              NE citizen               0.08 0.03 0.0051 0.0507      0.02
## 7                NE other              -0.09 0.06 0.1216 0.9729     -0.21
## 8               No answer               0.35 0.76 0.6436 1.0000     -1.14
## 9             Other party               0.02 0.01 0.0837 0.7533      0.00
## 10  Pro-environment party              -0.01 0.03 0.7334 1.0000     -0.06
##    asymp.UCL
## 1       0.03
## 2       0.04
## 3       0.09
## 4       0.94
## 5       0.06
## 6       0.14
## 7       0.02
## 8       1.85
## 9       0.03
## 10      0.04
```

```r
write.csv2(H5.exp.intr.mod6.trends.tab,"H5.exp.intr.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast     estimate        SE  df
## 1 Pro-environment party - Anti-immigration party -0.001737479 0.0314156 Inf
##       z.ratio   p.value
## 1 -0.05530626 0.9558945
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.intr.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.0644 0.0986 Inf -0.653  1.0000 
##  Did not vote effect            -0.0460 0.0979 Inf -0.470  1.0000 
##  Don't know effect              -0.0269 0.1027 Inf -0.262  1.0000 
##  Invalid vote effect             0.0901 0.4214 Inf  0.214  1.0000 
##  NE age effect                  -0.0558 0.1012 Inf -0.551  1.0000 
##  NE citizen effect               0.0318 0.1009 Inf  0.315  1.0000 
##  NE other effect                -0.1585 0.1133 Inf -1.399  1.0000 
##  No answer effect                0.3355 0.7643 Inf  0.439  1.0000 
##  Other party effect             -0.0394 0.0972 Inf -0.406  1.0000 
##  Pro-environment party effect   -0.0663 0.1001 Inf -0.662  1.0000 
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
## 1           Other party - Anti-immigration party  0.022447865 0.02035632 Inf
## 2 Pro-environment party - Anti-immigration party -0.001737479 0.03141560 Inf
## 3            Pro-environment party - Other party -0.024185344 0.02687443 Inf
##       z.ratio   p.value
## 1  1.10274645 0.8104119
## 2 -0.05530626 0.9558945
## 3 -0.89993896 0.8104119
```

