---
title: "Hypothesis 4"
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.14    62.06
## 2                                                      age  0.00 0.00 35705.43
## 3                                                   gender  0.06 0.01 35593.59
## 4                                                     educ  0.01 0.00 35719.00
## 5                                                    resid -0.06 0.01 35683.44
## 6                            occupClerical support workers -0.04 0.09 35568.33
## 7                    occupCraft and related trades workers -0.07 0.09 35579.01
## 8                              occupElementary occupations  0.01 0.09 35575.90
## 9                                            occupManagers -0.02 0.09 35568.40
## 10                            occupOther: Not in paid work  0.12 0.09 35624.19
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35576.73
## 12                                      occupProfessionals  0.07 0.09 35571.35
## 13                                            occupRetired -0.04 0.10 35560.28
## 14                          occupService and sales workers -0.04 0.09 35572.75
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35587.51
## 16            occupTechnicians and associate professionals -0.03 0.09 35568.43
## 17                                         occupUnemployed -0.03 0.11 35569.10
## 18                                            environ.lvl1  0.12 0.01    19.43
## 19                            all.parties.lvl2Did not vote  0.45 0.06   182.54
## 20                              all.parties.lvl2Don't know  0.42 0.07   269.12
## 21                            all.parties.lvl2Invalid vote  0.49 0.35  1099.45
## 22                                  all.parties.lvl2NE age  0.74 0.07   272.45
## 23                              all.parties.lvl2NE citizen  0.86 0.08   268.81
## 24                                all.parties.lvl2NE other  0.72 0.10   659.04
## 25                               all.parties.lvl2No answer  0.56 0.37  1373.88
## 26                             all.parties.lvl2Other party  0.54 0.05   230.96
## 27                   all.parties.lvl2Pro-environment party  0.91 0.06   256.35
##        t     p    LL    UL
## 1  -3.22 0.002 -0.75 -0.18
## 2   4.89 0.000  0.00  0.00
## 3   4.67 0.000  0.03  0.08
## 4   7.54 0.000  0.01  0.02
## 5  -5.14 0.000 -0.08 -0.04
## 6  -0.51 0.611 -0.22  0.13
## 7  -0.78 0.434 -0.24  0.10
## 8   0.09 0.931 -0.17  0.18
## 9  -0.24 0.810 -0.19  0.15
## 10  1.32 0.186 -0.06  0.30
## 11 -0.76 0.447 -0.24  0.11
## 12  0.76 0.450 -0.11  0.24
## 13 -0.41 0.682 -0.23  0.15
## 14 -0.46 0.647 -0.21  0.13
## 15 -0.66 0.506 -0.24  0.12
## 16 -0.40 0.690 -0.21  0.14
## 17 -0.30 0.763 -0.24  0.18
## 18 11.50 0.000  0.10  0.15
## 19  7.09 0.000  0.32  0.57
## 20  5.99 0.000  0.28  0.56
## 21  1.38 0.169 -0.21  1.18
## 22 10.51 0.000  0.60  0.88
## 23 11.34 0.000  0.71  1.01
## 24  7.33 0.000  0.52  0.91
## 25  1.50 0.133 -0.17  1.30
## 26 11.28 0.000  0.45  0.64
## 27 14.05 0.000  0.78  1.04
```

```r
(VC.H4.mod1<-getVC(H4.mod1))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.19740553  0.0389689413
## 2 voting.group environ.lvl1         <NA>  0.04295429  0.0018450708
## 3 voting.group  (Intercept) environ.lvl1  0.06604155  0.0005599936
## 4        cntry  (Intercept)         <NA>  0.48262919  0.2329309370
## 5        cntry environ.lvl1         <NA>  0.03997237  0.0015977908
## 6        cntry  (Intercept) environ.lvl1 -0.05083406 -0.0009806823
## 7     Residual         <NA>         <NA>  1.02913757  1.0591241343
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H4.mod1   34 104123 104411 -52027   104055                     
## H4.mod2   43 104133 104498 -52024   104047 7.4587  9     0.5895
```

```r
(FE.H4.mod2<-getFE(H4.mod2))
```

```
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.47 0.14    62.05
## 2                                                      age  0.00 0.00 35703.98
## 3                                                   gender  0.06 0.01 35594.36
## 4                                                     educ  0.01 0.00 35717.58
## 5                                                    resid -0.06 0.01 35680.08
## 6                            occupClerical support workers -0.05 0.09 35565.71
## 7                    occupCraft and related trades workers -0.07 0.09 35577.12
## 8                              occupElementary occupations  0.01 0.09 35573.30
## 9                                            occupManagers -0.02 0.09 35565.33
## 10                            occupOther: Not in paid work  0.12 0.09 35622.51
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35574.13
## 12                                      occupProfessionals  0.07 0.09 35567.54
## 13                                            occupRetired -0.04 0.10 35553.77
## 14                          occupService and sales workers -0.04 0.09 35570.97
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35584.95
## 16            occupTechnicians and associate professionals -0.04 0.09 35565.98
## 17                                         occupUnemployed -0.03 0.11 35567.13
## 18                                            environ.lvl1  0.11 0.02   106.04
## 19                            all.parties.lvl2Did not vote  0.45 0.06   182.70
## 20                              all.parties.lvl2Don't know  0.42 0.07   268.96
## 21                            all.parties.lvl2Invalid vote  0.49 0.35  1099.20
## 22                                  all.parties.lvl2NE age  0.74 0.07   272.33
## 23                              all.parties.lvl2NE citizen  0.86 0.08   268.70
## 24                                all.parties.lvl2NE other  0.72 0.10   658.23
## 25                               all.parties.lvl2No answer  0.57 0.37  1373.19
## 26                             all.parties.lvl2Other party  0.55 0.05   230.94
## 27                   all.parties.lvl2Pro-environment party  0.91 0.06   256.28
## 28               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    85.05
## 29                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   333.74
## 30               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25263.41
## 31                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   232.54
## 32                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   231.02
## 33                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1198.34
## 34                  environ.lvl1:all.parties.lvl2No answer  0.23 0.23  7338.62
## 35                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   132.60
## 36      environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   191.66
##        t     p    LL    UL
## 1  -3.22 0.002 -0.75 -0.18
## 2   4.89 0.000  0.00  0.00
## 3   4.68 0.000  0.03  0.08
## 4   7.53 0.000  0.01  0.02
## 5  -5.15 0.000 -0.08 -0.04
## 6  -0.52 0.603 -0.22  0.13
## 7  -0.79 0.428 -0.24  0.10
## 8   0.07 0.943 -0.17  0.18
## 9  -0.25 0.800 -0.20  0.15
## 10  1.32 0.188 -0.06  0.30
## 11 -0.77 0.442 -0.24  0.11
## 12  0.75 0.456 -0.11  0.24
## 13 -0.41 0.680 -0.23  0.15
## 14 -0.47 0.640 -0.21  0.13
## 15 -0.67 0.503 -0.24  0.12
## 16 -0.41 0.684 -0.21  0.14
## 17 -0.30 0.765 -0.24  0.18
## 18  5.56 0.000  0.07  0.15
## 19  7.07 0.000  0.32  0.57
## 20  6.00 0.000  0.29  0.56
## 21  1.38 0.168 -0.21  1.19
## 22 10.54 0.000  0.61  0.88
## 23 11.31 0.000  0.71  1.01
## 24  7.34 0.000  0.53  0.91
## 25  1.51 0.131 -0.17  1.30
## 26 11.31 0.000  0.45  0.64
## 27 14.07 0.000  0.78  1.04
## 28  0.09 0.926 -0.04  0.05
## 29  0.82 0.413 -0.04  0.10
## 30  0.15 0.878 -0.70  0.82
## 31  1.28 0.203 -0.02  0.10
## 32 -1.12 0.265 -0.10  0.03
## 33  0.59 0.557 -0.08  0.14
## 34  0.99 0.320 -0.22  0.68
## 35  0.98 0.329 -0.02  0.06
## 36  0.92 0.358 -0.03  0.08
```

```r
(VC.H4.mod2<-getVC(H4.mod2))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.19740057  0.0389669851
## 2 voting.group environ.lvl1         <NA>  0.04068572  0.0016553277
## 3 voting.group  (Intercept) environ.lvl1  0.07646619  0.0006141293
## 4        cntry  (Intercept)         <NA>  0.48264308  0.2329443406
## 5        cntry environ.lvl1         <NA>  0.04086829  0.0016702173
## 6        cntry  (Intercept) environ.lvl1 -0.05205339 -0.0010267426
## 7     Residual         <NA>         <NA>  1.02909763  1.0590419302
```

\newpage

#### Marginal effect for pro-environment and anti-immigration voters


```r
H4.mod2.trends<-emtrends(H4.mod2,
                         specs = c("all.parties.lvl2"),
                         var=c("environ.lvl1"))
H4.mod2.trends.tab<-data.frame(H4.mod2.trends)

H4.mod2.trends.tab$p<-
  2*(1-pnorm(abs(H4.mod2.trends.tab$environ.lvl1.trend/
                   H4.mod2.trends.tab$SE)))
H4.mod2.trends.tab$adj.p<-
  p.adjust(H4.mod2.trends.tab$p,method="holm")

H4.mod2.trends.tab<-
  cbind.data.frame(group=H4.mod2.trends.tab[,"all.parties.lvl2"],
                   beta=round_tidy(H4.mod2.trends.tab[,"environ.lvl1.trend"],2),
                   CI=paste0("[",
                             round_tidy(H4.mod2.trends.tab[,"asymp.LCL"],2),
                             ", ",
                             round_tidy(H4.mod2.trends.tab[,"asymp.UCL"],2),
                             "]"),
                   p=round_tidy(H4.mod2.trends.tab[,"p"],3),
                   p.adj=round_tidy(H4.mod2.trends.tab[,"adj.p"],3),
                   p_less_001=ifelse(H4.mod2.trends.tab[,"p"]<.001,"yes","no"),
                   p.adj_less_001=ifelse(H4.mod2.trends.tab[,"adj.p"]<.001,"yes","no"))

H4.mod2.trends.tab
```

```
##                     group beta            CI     p p.adj p_less_001
## 1  Anti-immigration party 0.11  [0.07, 0.15] 0.000 0.000        yes
## 2            Did not vote 0.11  [0.08, 0.14] 0.000 0.000        yes
## 3              Don't know 0.14  [0.08, 0.20] 0.000 0.000        yes
## 4            Invalid vote 0.17 [-0.59, 0.93] 0.664 0.664         no
## 5                  NE age 0.15  [0.10, 0.20] 0.000 0.000        yes
## 6              NE citizen 0.07  [0.01, 0.13] 0.016 0.048         no
## 7                NE other 0.14  [0.04, 0.25] 0.008 0.034         no
## 8               No answer 0.34 [-0.11, 0.79] 0.141 0.282         no
## 9             Other party 0.13  [0.10, 0.15] 0.000 0.000        yes
## 10  Pro-environment party 0.14  [0.09, 0.18] 0.000 0.000        yes
##    p.adj_less_001
## 1             yes
## 2             yes
## 3             yes
## 4              no
## 5             yes
## 6              no
## 7              no
## 8              no
## 9             yes
## 10            yes
```

```r
write.csv2(H4.mod2.trends.tab,"H4.mod2.trends.tab.csv")


#contrast between anti-immigration and pro-environment
(H4.contrast<-data.frame(pairs(H4.mod2.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.02618252 0.02840579 Inf
##     z.ratio   p.value
## 1 0.9217319 0.3566685
```

```r
#contrast for all groups against mean of other groups
contrast(H4.mod2.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect -0.04420 0.0538 Inf -0.821  0.4117 
##  Did not vote effect           -0.04188 0.0529 Inf -0.792  0.4283 
##  Don't know effect             -0.01316 0.0587 Inf -0.224  0.8226 
##  Invalid vote effect            0.02239 0.3905 Inf  0.057  0.9543 
##  NE age effect                 -0.00183 0.0564 Inf -0.032  0.9742 
##  NE citizen effect             -0.08565 0.0583 Inf -1.470  0.1416 
##  NE other effect               -0.00781 0.0734 Inf -0.106  0.9152 
##  No answer effect               0.21035 0.2339 Inf  0.899  0.3685 
##  Other party effect            -0.02310 0.0516 Inf -0.448  0.6545 
##  Pro-environment party effect  -0.01511 0.0556 Inf -0.272  0.7859 
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
## 1           Other party - Anti-immigration party 0.018994304 0.01937022 Inf
## 2 Pro-environment party - Anti-immigration party 0.026182522 0.02840579 Inf
## 3            Pro-environment party - Other party 0.007188218 0.02369191 Inf
##     z.ratio   p.value
## 1 0.9805933 0.3267934
## 2 0.9217319 0.3566685
## 3 0.3034039 0.7615821
```



\newpage

#### Compile Table 2


```r
# Compile Table 2

ref.means<-read.csv2("H2.mod2.mmeans.tab.csv",
                     stringsAsFactors = F)

env.means<-read.csv2("H3.mod2.mmeans.tab.csv",
                     stringsAsFactors = F)

group.trends<-read.csv2("H4.mod2.trends.tab.csv",
                        stringsAsFactors = F)

tab2<-cbind.data.frame(ref.means[,c("group","emmean","CI","p")],
                 env.means[,c("emmean","CI","p")],
                 group.trends[,c("beta","CI","p")])
tab2
```

```
##                     group emmean             CI     p emmean             CI
## 1  Anti-immigration party  -0.48 [-0.71, -0.25] 0.000  -0.22 [-0.36, -0.08]
## 2            Did not vote  -0.03  [-0.26, 0.20] 0.781  -0.12  [-0.27, 0.02]
## 3              Don't know  -0.06  [-0.30, 0.18] 0.630  -0.15  [-0.30, 0.01]
## 4            Invalid vote   0.02  [-0.71, 0.74] 0.963  -0.84 [-1.56, -0.12]
## 5                  NE age   0.25   [0.01, 0.49] 0.040   0.17   [0.01, 0.33]
## 6              NE citizen   0.37   [0.13, 0.61] 0.003   0.02  [-0.15, 0.18]
## 7                NE other   0.22  [-0.05, 0.50] 0.109   0.07  [-0.14, 0.29]
## 8               No answer   0.08  [-0.68, 0.84] 0.839  -0.02  [-0.79, 0.75]
## 9             Other party   0.06  [-0.16, 0.27] 0.596  -0.01  [-0.13, 0.12]
## 10  Pro-environment party   0.41   [0.18, 0.65] 0.000   0.48   [0.33, 0.63]
##        p beta            CI     p
## 1  0.002 0.11  [0.07, 0.15] 0.000
## 2  0.086 0.11  [0.08, 0.14] 0.000
## 3  0.066 0.14  [0.08, 0.20] 0.000
## 4  0.023 0.17 [-0.59, 0.93] 0.664
## 5  0.034 0.15  [0.10, 0.20] 0.000
## 6  0.850 0.07  [0.01, 0.13] 0.016
## 7  0.489 0.14  [0.04, 0.25] 0.008
## 8  0.961 0.34 [-0.11, 0.79] 0.141
## 9  0.911 0.13  [0.10, 0.15] 0.000
## 10 0.000 0.14  [0.09, 0.18] 0.000
```

```r
write.csv2(tab2,"tab2.csv")
```
