---
title: "Hypothesis 2: Those who voted for pro-environment parties will report higher pro-refugee attitudes than those who voted for anti-immigration parties"
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

# Hypothesis 2: Those who voted for pro-environment parties will report higher pro-refugee attitudes than those who voted for anti-immigration parties

### Model 1: random intercepts + covariates (same as in H1)


```r
H2.mod1<-lmer(refugees~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup
                ,data=dat,REML=F)


(FE.H2.mod1<-getFE(H2.mod1))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.06 0.14    50.42
## 2                                                      age  0.00 0.00 34095.35
## 3                                                   gender  0.06 0.01 35598.09
## 4                                                     educ  0.02 0.00 35706.37
## 5                                                    resid -0.08 0.01 35657.25
## 6                            occupClerical support workers -0.04 0.09 35532.93
## 7                    occupCraft and related trades workers -0.08 0.09 35539.44
## 8                              occupElementary occupations  0.01 0.09 35540.85
## 9                                            occupManagers  0.00 0.09 35533.47
## 10                            occupOther: Not in paid work  0.16 0.09 35696.52
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35539.06
## 12                                      occupProfessionals  0.10 0.09 35535.19
## 13                                            occupRetired -0.03 0.10 35538.37
## 14                          occupService and sales workers -0.04 0.09 35534.90
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35544.87
## 16            occupTechnicians and associate professionals -0.03 0.09 35532.10
## 17                                         occupUnemployed -0.03 0.11 35546.28
##        t     p    LL    UL
## 1   0.44 0.663 -0.22  0.35
## 2   3.30 0.001  0.00  0.00
## 3   4.77 0.000  0.03  0.08
## 4   9.28 0.000  0.01  0.02
## 5  -6.58 0.000 -0.10 -0.05
## 6  -0.40 0.687 -0.21  0.14
## 7  -0.85 0.398 -0.25  0.10
## 8   0.14 0.892 -0.16  0.19
## 9   0.00 0.998 -0.18  0.18
## 10  1.71 0.088 -0.02  0.34
## 11 -0.77 0.442 -0.24  0.11
## 12  1.10 0.273 -0.08  0.27
## 13 -0.28 0.782 -0.22  0.17
## 14 -0.40 0.691 -0.21  0.14
## 15 -0.52 0.601 -0.23  0.14
## 16 -0.29 0.774 -0.20  0.15
## 17 -0.23 0.816 -0.24  0.19
```

```r
(VC.H2.mod1<-getVC(H2.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.3034866 0.09210413
## 2        cntry (Intercept) <NA> 0.4986022 0.24860415
## 3     Residual        <NA> <NA> 1.0417310 1.08520339
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.15    62.98
## 2                                                      age  0.00 0.00 35716.28
## 3                                                   gender  0.06 0.01 35620.34
## 4                                                     educ  0.02 0.00 35722.81
## 5                                                    resid -0.08 0.01 35707.05
## 6                            occupClerical support workers -0.04 0.09 35593.68
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66
## 8                              occupElementary occupations  0.01 0.09 35596.76
## 9                                            occupManagers -0.00 0.09 35595.21
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23
## 12                                      occupProfessionals  0.09 0.09 35595.08
## 13                                            occupRetired -0.03 0.10 35590.71
## 14                          occupService and sales workers -0.04 0.09 35594.22
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11
## 17                                         occupUnemployed -0.03 0.11 35590.31
## 18                            all.parties.lvl2Did not vote  0.45 0.06   181.24
## 19                              all.parties.lvl2Don't know  0.42 0.07   270.53
## 20                            all.parties.lvl2Invalid vote  0.50 0.36  1156.29
## 21                                  all.parties.lvl2NE age  0.73 0.07   274.23
## 22                              all.parties.lvl2NE citizen  0.85 0.08   268.85
## 23                                all.parties.lvl2NE other  0.71 0.10   672.84
## 24                               all.parties.lvl2No answer  0.56 0.38  1449.94
## 25                             all.parties.lvl2Other party  0.54 0.05   230.15
## 26                   all.parties.lvl2Pro-environment party  0.90 0.06   256.45
##        t     p    LL    UL
## 1  -3.25 0.002 -0.76 -0.18
## 2   3.73 0.000  0.00  0.00
## 3   4.74 0.000  0.03  0.08
## 4   9.31 0.000  0.01  0.02
## 5  -6.49 0.000 -0.10 -0.05
## 6  -0.44 0.661 -0.21  0.14
## 7  -0.86 0.390 -0.25  0.10
## 8   0.09 0.928 -0.17  0.18
## 9  -0.03 0.978 -0.18  0.17
## 10  1.49 0.135 -0.04  0.32
## 11 -0.80 0.426 -0.25  0.10
## 12  1.06 0.290 -0.08  0.27
## 13 -0.33 0.741 -0.23  0.16
## 14 -0.42 0.677 -0.21  0.14
## 15 -0.56 0.573 -0.24  0.13
## 16 -0.31 0.753 -0.20  0.15
## 17 -0.30 0.763 -0.25  0.18
## 18  7.15 0.000  0.33  0.57
## 19  6.02 0.000  0.28  0.56
## 20  1.40 0.161 -0.20  1.20
## 21 10.42 0.000  0.59  0.87
## 22 11.25 0.000  0.70  1.00
## 23  7.21 0.000  0.51  0.90
## 24  1.49 0.136 -0.18  1.30
## 25 11.29 0.000  0.45  0.63
## 26 13.92 0.000  0.77  1.02
```

```r
(VC.H2.mod2<-getVC(H2.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1942751 0.03774283
## 2        cntry (Intercept) <NA> 0.4821835 0.23250093
## 3     Residual        <NA> <NA> 1.0417292 1.08519972
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H2.mod1   20 105044 105214 -52502   105004                         
## H2.mod2   29 104880 105126 -52411   104822 182.15  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H2.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
## age               15.095  15.095     1 35716 13.9102 0.0001921 ***
## gender            24.404  24.404     1 35620 22.4878 2.123e-06 ***
## educ              94.051  94.051     1 35723 86.6668 < 2.2e-16 ***
## resid             45.677  45.677     1 35707 42.0906 8.828e-11 ***
## occup            115.467   9.622    12 35609  8.8668 < 2.2e-16 ***
## all.parties.lvl2 281.610  31.290     9   341 28.8334 < 2.2e-16 ***
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
## [1] 0.5902157
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
  cbind.data.frame(group=H2.mod2.mmeans.tab[,"group"],
        emmean=round_tidy(H2.mod2.mmeans.tab[,"emmean"],2),
        CI=paste0("[",
                  round_tidy(H2.mod2.mmeans.tab[,"asymp.LCL"],2),
                  ", ",
                  round_tidy(H2.mod2.mmeans.tab[,"asymp.UCL"],2),
                  "]"),
        p=round_tidy(H2.mod2.mmeans.tab[,"p"],3),
        p.adj=round_tidy(H2.mod2.mmeans.tab[,"adj.p"],3),
        p_less_001=ifelse(H2.mod2.mmeans.tab[,"p"]<.001,"yes","no"),
        p.adj_less_001=ifelse(H2.mod2.mmeans.tab[,"adj.p"]<.001,"yes","no"))
H2.mod2.mmeans.tab
```

```
##                     group emmean             CI     p p.adj p_less_001
## 1  Anti-immigration party  -0.48 [-0.71, -0.25] 0.000 0.000        yes
## 2            Did not vote  -0.03  [-0.26, 0.20] 0.781 1.000         no
## 3              Don't know  -0.06  [-0.30, 0.18] 0.630 1.000         no
## 4            Invalid vote   0.02  [-0.71, 0.74] 0.963 1.000         no
## 5                  NE age   0.25   [0.01, 0.49] 0.040 0.278         no
## 6              NE citizen   0.37   [0.13, 0.61] 0.003 0.024         no
## 7                NE other   0.22  [-0.05, 0.50] 0.109 0.651         no
## 8               No answer   0.08  [-0.68, 0.84] 0.839 1.000         no
## 9             Other party   0.06  [-0.16, 0.27] 0.596 1.000         no
## 10  Pro-environment party   0.41   [0.18, 0.65] 0.000 0.004        yes
##    p.adj_less_001
## 1             yes
## 2              no
## 3              no
## 4              no
## 5              no
## 6              no
## 7              no
## 8              no
## 9              no
## 10             no
```

```r
write.csv2(H2.mod2.mmeans.tab,"H2.mod2.mmeans.tab.csv")


#contrast between anti-immigration and pro-environment
(H2.contrast<-data.frame(pairs(H2.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.8955111 0.06434551 Inf
##    z.ratio      p.value
## 1 13.91723 4.978424e-44
```

```r
#contrast for all groups against mean of other groups
contrast(H2.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect -0.62831 0.0735 Inf -8.551  <.0001 
##  Did not vote effect           -0.12947 0.0749 Inf -1.728  0.4199 
##  Don't know effect             -0.15839 0.0812 Inf -1.950  0.3071 
##  Invalid vote effect           -0.07419 0.3559 Inf -0.208  1.0000 
##  NE age effect                  0.18394 0.0810 Inf  2.271  0.1620 
##  NE citizen effect              0.31827 0.0857 Inf  3.713  0.0016 
##  NE other effect                0.15526 0.1058 Inf  1.468  0.5686 
##  No answer effect              -0.00533 0.3764 Inf -0.014  1.0000 
##  Other party effect            -0.02848 0.0628 Inf -0.453  1.0000 
##  Pro-environment party effect   0.36670 0.0763 Inf  4.805  <.0001 
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
## 1           Other party - Anti-immigration party 0.5398443 0.04779971 Inf
## 2 Pro-environment party - Anti-immigration party 0.8955111 0.06434551 Inf
## 3            Pro-environment party - Other party 0.3556668 0.05168896 Inf
##     z.ratio      p.value
## 1 11.293884 2.813585e-29
## 2 13.917228 1.493527e-43
## 3  6.880904 5.947390e-12
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
## [1] 1.076106
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
## [1] 1.001363
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
## [1] 1.039407
```

```r
(H2.effect.size<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd)
```

```
## Warning in Ops.factor(H2.mod2.mmeans.tab[10, 2], H2.mod2.mmeans.tab[1, 2]): '-'
## not meaningful for factors
```

```
## [1] NA
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
## [1] 1.049983
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
## [1] 1.028043
```

```r
(H2.effect.size.env.other<-(H2.mod2.mmeans.tab[10,2]-
  H2.mod2.mmeans.tab[9,2])/H2.pooled.sd.other.env)
```

```
## Warning in Ops.factor(H2.mod2.mmeans.tab[10, 2], H2.mod2.mmeans.tab[9, 2]): '-'
## not meaningful for factors
```

```
## [1] NA
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
## [1] 1.062005
```

```r
(H2.effect.size.imm.other<-(H2.mod2.mmeans.tab[9,2]-
  H2.mod2.mmeans.tab[1,2])/H2.pooled.sd.other.imm)
```

```
## Warning in Ops.factor(H2.mod2.mmeans.tab[9, 2], H2.mod2.mmeans.tab[1, 2]): '-'
## not meaningful for factors
```

```
## [1] NA
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.02 0.15    64.39
## 2                                                      age  0.00 0.00 35716.28
## 3                                                   gender  0.06 0.01 35620.34
## 4                                                     educ  0.02 0.00 35722.81
## 5                                                    resid -0.08 0.01 35707.05
## 6                            occupClerical support workers -0.04 0.09 35593.68
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66
## 8                              occupElementary occupations  0.01 0.09 35596.76
## 9                                            occupManagers -0.00 0.09 35595.21
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23
## 12                                      occupProfessionals  0.09 0.09 35595.08
## 13                                            occupRetired -0.03 0.10 35590.71
## 14                          occupService and sales workers -0.04 0.09 35594.23
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11
## 17                                         occupUnemployed -0.03 0.11 35590.31
## 18                                       other.party.dummy  0.09 0.05   163.48
## 19                                         dont.know.dummy -0.03 0.07   222.65
## 20                                      invalid.vote.dummy  0.05 0.36  1136.40
## 21                                         no.answer.dummy  0.11 0.38  1416.47
## 22                                  not.eligible.age.dummy  0.28 0.07   219.09
## 23                          not.eligible.citizenship.dummy  0.40 0.08   226.00
## 24                                not.eligible.other.dummy  0.26 0.10   558.14
## 25                                    anti.imm.party.dummy -0.45 0.06   181.24
## 26                                     pro.env.party.dummy  0.45 0.07   206.58
##        t     p    LL    UL
## 1  -0.16 0.875 -0.31  0.27
## 2   3.73 0.000  0.00  0.00
## 3   4.74 0.000  0.03  0.08
## 4   9.31 0.000  0.01  0.02
## 5  -6.49 0.000 -0.10 -0.05
## 6  -0.44 0.661 -0.21  0.14
## 7  -0.86 0.390 -0.25  0.10
## 8   0.09 0.928 -0.17  0.18
## 9  -0.03 0.978 -0.18  0.17
## 10  1.49 0.135 -0.04  0.32
## 11 -0.80 0.426 -0.25  0.10
## 12  1.06 0.290 -0.08  0.27
## 13 -0.33 0.741 -0.23  0.16
## 14 -0.42 0.677 -0.21  0.14
## 15 -0.56 0.573 -0.24  0.13
## 16 -0.31 0.753 -0.20  0.15
## 17 -0.30 0.763 -0.25  0.18
## 18  1.82 0.071 -0.01  0.19
## 19 -0.36 0.717 -0.17  0.12
## 20  0.14 0.889 -0.65  0.75
## 21  0.30 0.767 -0.63  0.85
## 22  3.96 0.000  0.14  0.42
## 23  5.23 0.000  0.25  0.55
## 24  2.59 0.010  0.06  0.45
## 25 -7.15 0.000 -0.57 -0.33
## 26  6.76 0.000  0.32  0.58
```

```r
(VC.H2.mod3<-getVC(H2.mod3))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1942750 0.03774278
## 2        cntry (Intercept) <NA> 0.4821838 0.23250125
## 3     Residual        <NA> <NA> 1.0417292 1.08519972
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
##         npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)    
## H2.mod2   29 104880 105126 -52411   104822                        
## H2.mod3   29 104880 105126 -52411   104822     0  0  < 2.2e-16 ***
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.02 0.15    63.03
## 2                                                      age  0.00 0.00 35701.12
## 3                                                   gender  0.06 0.01 35620.05
## 4                                                     educ  0.02 0.00 35718.65
## 5                                                    resid -0.08 0.01 35708.16
## 6                            occupClerical support workers -0.04 0.09 35595.60
## 7                    occupCraft and related trades workers -0.08 0.09 35602.53
## 8                              occupElementary occupations  0.01 0.09 35598.49
## 9                                            occupManagers -0.00 0.09 35596.94
## 10                            occupOther: Not in paid work  0.14 0.09 35647.42
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35602.28
## 12                                      occupProfessionals  0.09 0.09 35596.84
## 13                                            occupRetired -0.03 0.10 35592.09
## 14                          occupService and sales workers -0.04 0.09 35595.92
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35612.90
## 16            occupTechnicians and associate professionals -0.03 0.09 35594.87
## 17                                         occupUnemployed -0.03 0.11 35592.09
## 18                                       other.party.dummy  0.09 0.05   132.00
## 19                                         dont.know.dummy -0.03 0.07   184.46
## 20                                      invalid.vote.dummy  0.05 0.35  1046.75
## 21                                         no.answer.dummy  0.10 0.37  1321.51
## 22                                  not.eligible.age.dummy  0.28 0.07   181.50
## 23                          not.eligible.citizenship.dummy  0.40 0.07   186.15
## 24                                not.eligible.other.dummy  0.25 0.10   483.87
## 25                                    anti.imm.party.dummy -0.45 0.07    38.57
## 26                                     pro.env.party.dummy  0.44 0.08    30.86
##        t     p    LL    UL
## 1  -0.16 0.873 -0.31  0.27
## 2   3.74 0.000  0.00  0.00
## 3   4.74 0.000  0.03  0.08
## 4   9.31 0.000  0.01  0.02
## 5  -6.51 0.000 -0.10 -0.05
## 6  -0.44 0.661 -0.21  0.14
## 7  -0.86 0.391 -0.25  0.10
## 8   0.09 0.926 -0.17  0.18
## 9  -0.03 0.976 -0.18  0.17
## 10  1.49 0.136 -0.04  0.32
## 11 -0.79 0.427 -0.25  0.10
## 12  1.05 0.292 -0.08  0.27
## 13 -0.33 0.741 -0.23  0.16
## 14 -0.42 0.678 -0.21  0.14
## 15 -0.56 0.573 -0.24  0.13
## 16 -0.32 0.753 -0.20  0.15
## 17 -0.29 0.768 -0.24  0.18
## 18  1.89 0.061 -0.00  0.18
## 19 -0.38 0.701 -0.16  0.11
## 20  0.15 0.883 -0.63  0.74
## 21  0.26 0.796 -0.63  0.82
## 22  4.12 0.000  0.15  0.42
## 23  5.44 0.000  0.26  0.55
## 24  2.64 0.009  0.06  0.44
## 25 -6.83 0.000 -0.58 -0.32
## 26  5.76 0.000  0.28  0.59
```

```r
(VC.H2.mod4<-getVC(H2.mod4))
```

```
##            grp                 var1 var2    est_SD    est_SD2
## 1 voting.group          (Intercept) <NA> 0.1836158 0.03371476
## 2        cntry          (Intercept) <NA> 0.4828848 0.23317769
## 3      cntry.1 anti.imm.party.dummy <NA> 0.1103878 0.01218547
## 4      cntry.2  pro.env.party.dummy <NA> 0.1661464 0.02760461
## 5     Residual                 <NA> <NA> 1.0417661 1.08527658
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H2.mod3   29 104880 105126 -52411   104822                     
## H2.mod4   31 104880 105143 -52409   104818 3.3674  2     0.1857
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.03 0.14    51.13
## 2                                                      age  0.00 0.00 35707.89
## 3                                                   gender  0.06 0.01 35578.58
## 4                                                     educ  0.02 0.00 35682.00
## 5                                                    resid -0.08 0.01 35640.73
## 6                            occupClerical support workers -0.04 0.09 35542.31
## 7                    occupCraft and related trades workers -0.08 0.09 35548.02
## 8                              occupElementary occupations  0.01 0.09 35545.03
## 9                                            occupManagers -0.00 0.09 35543.00
## 10                            occupOther: Not in paid work  0.14 0.09 35580.80
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35546.89
## 12                                      occupProfessionals  0.10 0.09 35544.96
## 13                                            occupRetired -0.03 0.10 35542.87
## 14                          occupService and sales workers -0.04 0.09 35543.35
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35554.44
## 16            occupTechnicians and associate professionals -0.03 0.09 35541.81
## 17                                         occupUnemployed -0.03 0.11 35546.38
## 18                                      did.not.vote.dummy -0.06 0.07   197.35
## 19                                         dont.know.dummy -0.08 0.08   283.64
## 20                                      invalid.vote.dummy  0.01 0.41   684.76
## 21                                         no.answer.dummy  0.04 0.43   809.06
## 22                                  not.eligible.age.dummy  0.23 0.08   286.48
## 23                          not.eligible.citizenship.dummy  0.35 0.08   303.31
## 24                                not.eligible.other.dummy  0.22 0.11   711.72
##        t     p    LL    UL
## 1   0.23 0.823 -0.26  0.32
## 2   3.75 0.000  0.00  0.00
## 3   4.83 0.000  0.03  0.08
## 4   9.46 0.000  0.01  0.02
## 5  -6.49 0.000 -0.10 -0.05
## 6  -0.40 0.686 -0.21  0.14
## 7  -0.86 0.391 -0.25  0.10
## 8   0.10 0.920 -0.17  0.18
## 9  -0.01 0.992 -0.18  0.17
## 10  1.52 0.129 -0.04  0.32
## 11 -0.78 0.433 -0.25  0.11
## 12  1.09 0.276 -0.08  0.27
## 13 -0.31 0.755 -0.22  0.16
## 14 -0.42 0.677 -0.21  0.14
## 15 -0.53 0.596 -0.23  0.13
## 16 -0.29 0.772 -0.20  0.15
## 17 -0.29 0.774 -0.24  0.18
## 18 -0.80 0.425 -0.19  0.08
## 19 -1.10 0.273 -0.23  0.07
## 20  0.01 0.989 -0.81  0.82
## 21  0.08 0.933 -0.81  0.88
## 22  3.01 0.003  0.08  0.38
## 23  4.09 0.000  0.18  0.51
## 24  2.06 0.040  0.01  0.43
```

```r
(VC.H2.mod5<-getVC(H2.mod5))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2848488 0.08113883
## 2        cntry (Intercept) <NA> 0.4994842 0.24948451
## 3     Residual        <NA> <NA> 1.0416978 1.08513430
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H2.total.eff<-(VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"]-
     VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])/
  VC.H2.mod1[VC.H2.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1190533
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.12 0.14    53.92
## 2                                                      age  0.00 0.00 35722.27
## 3                                                   gender  0.06 0.01 35608.88
## 4                                                     educ  0.02 0.00 35718.92
## 5                                                    resid -0.08 0.01 35691.03
## 6                            occupClerical support workers -0.04 0.09 35576.83
## 7                    occupCraft and related trades workers -0.07 0.09 35584.02
## 8                              occupElementary occupations  0.01 0.09 35580.21
## 9                                            occupManagers -0.00 0.09 35577.58
## 10                            occupOther: Not in paid work  0.14 0.09 35625.80
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35583.03
## 12                                      occupProfessionals  0.10 0.09 35579.90
## 13                                            occupRetired -0.03 0.10 35574.87
## 14                          occupService and sales workers -0.04 0.09 35577.46
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35591.80
## 16            occupTechnicians and associate professionals -0.03 0.09 35576.48
## 17                                         occupUnemployed -0.03 0.11 35575.34
## 18                                      did.not.vote.dummy -0.15 0.05   170.98
## 19                                         dont.know.dummy -0.17 0.06   287.04
## 20                                      invalid.vote.dummy -0.06 0.37   987.47
## 21                                         no.answer.dummy -0.04 0.39  1213.22
## 22                                  not.eligible.age.dummy  0.13 0.06   295.41
## 23                          not.eligible.citizenship.dummy  0.25 0.07   288.85
## 24                                not.eligible.other.dummy  0.11 0.09   833.42
## 25                                    anti.imm.party.dummy -0.60 0.05   235.80
##         t     p    LL    UL
## 1    0.86 0.392 -0.16  0.40
## 2    3.59 0.000  0.00  0.00
## 3    4.79 0.000  0.03  0.08
## 4    9.43 0.000  0.01  0.02
## 5   -6.52 0.000 -0.10 -0.05
## 6   -0.41 0.680 -0.21  0.14
## 7   -0.83 0.406 -0.25  0.10
## 8    0.11 0.910 -0.17  0.19
## 9   -0.02 0.988 -0.18  0.17
## 10   1.53 0.127 -0.04  0.32
## 11  -0.77 0.440 -0.25  0.11
## 12   1.09 0.274 -0.08  0.27
## 13  -0.30 0.761 -0.22  0.16
## 14  -0.40 0.690 -0.21  0.14
## 15  -0.54 0.590 -0.23  0.13
## 16  -0.29 0.774 -0.20  0.15
## 17  -0.28 0.781 -0.24  0.18
## 18  -2.67 0.008 -0.25 -0.04
## 19  -2.75 0.006 -0.30 -0.05
## 20  -0.17 0.864 -0.78  0.66
## 21  -0.11 0.910 -0.80  0.72
## 22   2.13 0.034  0.01  0.26
## 23   3.64 0.000  0.12  0.39
## 24   1.21 0.228 -0.07  0.30
## 25 -11.61 0.000 -0.70 -0.49
```

```r
(VC.H2.mod6<-getVC(H2.mod6))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2175942 0.04734725
## 2        cntry (Intercept) <NA> 0.4852287 0.23544693
## 3     Residual        <NA> <NA> 1.0417471 1.08523704
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod6[VC.H2.mod6$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.4164662
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.03 0.14    52.30
## 2                                                      age  0.00 0.00 35717.46
## 3                                                   gender  0.06 0.01 35589.38
## 4                                                     educ  0.02 0.00 35687.43
## 5                                                    resid -0.08 0.01 35660.33
## 6                            occupClerical support workers -0.04 0.09 35558.32
## 7                    occupCraft and related trades workers -0.08 0.09 35564.26
## 8                              occupElementary occupations  0.01 0.09 35560.98
## 9                                            occupManagers -0.00 0.09 35559.65
## 10                            occupOther: Not in paid work  0.14 0.09 35599.92
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35563.39
## 12                                      occupProfessionals  0.09 0.09 35559.75
## 13                                            occupRetired -0.03 0.10 35557.95
## 14                          occupService and sales workers -0.04 0.09 35559.41
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35571.61
## 16            occupTechnicians and associate professionals -0.03 0.09 35557.67
## 17                                         occupUnemployed -0.03 0.11 35559.99
## 18                                      did.not.vote.dummy  0.01 0.06   189.71
## 19                                         dont.know.dummy -0.02 0.07   289.32
## 20                                      invalid.vote.dummy  0.03 0.39   804.27
## 21                                         no.answer.dummy  0.11 0.41   965.65
## 22                                  not.eligible.age.dummy  0.29 0.07   293.85
## 23                          not.eligible.citizenship.dummy  0.41 0.08   302.50
## 24                                not.eligible.other.dummy  0.28 0.10   768.38
## 25                                     pro.env.party.dummy  0.45 0.06   285.94
##        t     p    LL    UL
## 1  -0.19 0.851 -0.31  0.26
## 2   3.88 0.000  0.00  0.00
## 3   4.79 0.000  0.03  0.08
## 4   9.36 0.000  0.01  0.02
## 5  -6.46 0.000 -0.10 -0.05
## 6  -0.43 0.669 -0.21  0.14
## 7  -0.88 0.377 -0.25  0.10
## 8   0.08 0.936 -0.17  0.18
## 9  -0.02 0.984 -0.18  0.17
## 10  1.49 0.136 -0.04  0.32
## 11 -0.81 0.420 -0.25  0.10
## 12  1.06 0.290 -0.08  0.27
## 13 -0.34 0.736 -0.23  0.16
## 14 -0.43 0.665 -0.21  0.14
## 15 -0.55 0.582 -0.24  0.13
## 16 -0.31 0.754 -0.20  0.15
## 17 -0.31 0.757 -0.25  0.18
## 18  0.09 0.930 -0.12  0.13
## 19 -0.31 0.756 -0.16  0.12
## 20  0.07 0.948 -0.74  0.79
## 21  0.27 0.790 -0.70  0.92
## 22  4.13 0.000  0.15  0.43
## 23  5.23 0.000  0.26  0.56
## 24  2.74 0.006  0.08  0.47
## 25  7.32 0.000  0.33  0.57
```

```r
(VC.H2.mod7<-getVC(H2.mod7))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2547649 0.06490517
## 2        cntry (Intercept) <NA> 0.4942343 0.24426750
## 3     Residual        <NA> <NA> 1.0416779 1.08509294
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod7[VC.H2.mod7$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.2000726
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.07 0.14    54.74
## 2                                                      age  0.00 0.00 35716.28
## 3                                                   gender  0.06 0.01 35620.34
## 4                                                     educ  0.02 0.00 35722.81
## 5                                                    resid -0.08 0.01 35707.05
## 6                            occupClerical support workers -0.04 0.09 35593.68
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66
## 8                              occupElementary occupations  0.01 0.09 35596.75
## 9                                            occupManagers -0.00 0.09 35595.21
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23
## 12                                      occupProfessionals  0.09 0.09 35595.08
## 13                                            occupRetired -0.03 0.10 35590.71
## 14                          occupService and sales workers -0.04 0.09 35594.22
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11
## 17                                         occupUnemployed -0.03 0.11 35590.31
## 18                                      did.not.vote.dummy -0.09 0.05   163.48
## 19                                         dont.know.dummy -0.12 0.06   294.72
## 20                                      invalid.vote.dummy -0.04 0.35  1200.25
## 21                                         no.answer.dummy  0.02 0.37  1493.88
## 22                                  not.eligible.age.dummy  0.19 0.06   306.55
## 23                          not.eligible.citizenship.dummy  0.31 0.07   287.33
## 24                                not.eligible.other.dummy  0.17 0.09   897.11
## 25                                    anti.imm.party.dummy -0.54 0.05   230.15
## 26                                     pro.env.party.dummy  0.36 0.05   274.55
##         t     p    LL    UL
## 1    0.49 0.629 -0.21  0.35
## 2    3.73 0.000  0.00  0.00
## 3    4.74 0.000  0.03  0.08
## 4    9.31 0.000  0.01  0.02
## 5   -6.49 0.000 -0.10 -0.05
## 6   -0.44 0.661 -0.21  0.14
## 7   -0.86 0.390 -0.25  0.10
## 8    0.09 0.928 -0.17  0.18
## 9   -0.03 0.978 -0.18  0.17
## 10   1.49 0.135 -0.04  0.32
## 11  -0.80 0.426 -0.25  0.10
## 12   1.06 0.290 -0.08  0.27
## 13  -0.33 0.741 -0.23  0.16
## 14  -0.42 0.677 -0.21  0.14
## 15  -0.56 0.573 -0.24  0.13
## 16  -0.31 0.753 -0.20  0.15
## 17  -0.30 0.763 -0.25  0.18
## 18  -1.82 0.071 -0.19  0.01
## 19  -1.98 0.048 -0.23 -0.00
## 20  -0.12 0.907 -0.73  0.65
## 21   0.06 0.956 -0.71  0.76
## 22   3.23 0.001  0.07  0.31
## 23   4.78 0.000  0.18  0.44
## 24   1.83 0.068 -0.01  0.34
## 25 -11.29 0.000 -0.63 -0.45
## 26   6.88 0.000  0.25  0.46
```

```r
(VC.H2.mod8<-getVC(H2.mod8))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1942750 0.03774279
## 2        cntry (Intercept) <NA> 0.4821848 0.23250222
## 3     Residual        <NA> <NA> 1.0417292 1.08519972
```

```r
(H2.total.eff<-(VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"]-
     VC.H2.mod8[VC.H2.mod8$grp=="voting.group","est_SD2"])/
  VC.H2.mod5[VC.H2.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5348369
```

```r
#how much variance was left at level-2
1-H2.total.eff
```

```
## [1] 0.4651631
```

