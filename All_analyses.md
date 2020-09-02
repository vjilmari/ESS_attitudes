---
title: Association between Attitudes towards Refugees and Attitudes towards the Environment in the European Social Survey
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
    toc_depth: 5
  pdf_document: 
    toc: yes
    toc_depth: 4
fontsize: 12pt
geometry: margin=0.5in
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
## [1] LC_COLLATE=Finnish_Finland.1252  LC_CTYPE=Finnish_Finland.1252    LC_MONETARY=Finnish_Finland.1252
## [4] LC_NUMERIC=C                     LC_TIME=Finnish_Finland.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] finalfit_1.0.2  merTools_0.5.0  arm_1.10-1      MASS_7.3-51.5   metafor_2.4-0   ggplot2_3.3.0   emmeans_1.4.5  
##  [8] psych_1.9.12.31 dplyr_0.8.5     lmerTest_3.1-2  lme4_1.1-23     Matrix_1.2-18  
## 
## loaded via a namespace (and not attached):
##  [1] tidyr_1.0.2         splines_3.6.3       foreach_1.5.0       shiny_1.4.0.2       assertthat_0.2.1    statmod_1.4.34     
##  [7] yaml_2.2.0          numDeriv_2016.8-1.1 pillar_1.4.3        backports_1.1.6     lattice_0.20-38     glue_1.4.0         
## [13] digest_0.6.25       promises_1.1.0      minqa_1.2.4         colorspace_1.4-1    sandwich_2.5-1      htmltools_0.4.0    
## [19] httpuv_1.5.2        pkgconfig_2.0.3     broom_0.5.5         purrr_0.3.3         xtable_1.8-4        mvtnorm_1.1-0      
## [25] scales_1.1.0        later_1.0.0         tibble_3.0.0        generics_0.0.2      ellipsis_0.3.0      TH.data_1.0-10     
## [31] withr_2.1.2         cli_2.0.2           mnormt_1.5-6        survival_3.1-8      magrittr_1.5        crayon_1.3.4       
## [37] mime_0.9            estimability_1.3    evaluate_0.14       mice_3.8.0          fansi_0.4.1         nlme_3.1-144       
## [43] forcats_0.5.0       tools_3.6.3         lifecycle_0.2.0     multcomp_1.4-13     stringr_1.4.0       munsell_0.5.0      
## [49] compiler_3.6.3      rlang_0.4.5         blme_1.0-4          grid_3.6.3          nloptr_1.2.2.1      iterators_1.0.12   
## [55] rmarkdown_2.1       boot_1.3-24         gtable_0.3.0        codetools_0.2-16    abind_1.4-5         R6_2.4.1           
## [61] zoo_1.8-7           knitr_1.28          fastmap_1.0.1       stringi_1.4.6       parallel_3.6.3      Rcpp_1.0.4.6       
## [67] vctrs_0.2.4         tidyselect_1.0.0    xfun_0.13           coda_0.19-3
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
##   AT   BE   CH   CZ   DE   EE   ES   FI   FR   GB   HU   IE   IT   LT   NL   NO   PL   PT   SE   SI 
## 1973 1753 1503 2156 2819 1974 1817 1862 2015 1876 1391 2676 2317 1927 1661 1538 1589 1228 1525 1276
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




\newpage


# Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes

### Model 0: Intercepts only


```r
H1.mod0<-lmer(refugees~(1|voting.group)+(1|cntry),
              data=dat,REML=F)

(FE.H1.mod0<-getFE(H1.mod0))
```

```
##           Eff  Est   SE    df    t     p    LL   UL
## 1 (Intercept) 0.08 0.11 19.92 0.69 0.499 -0.16 0.32
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.06 0.14    50.42  0.44 0.663 -0.22  0.35
## 2                                                      age  0.00 0.00 34095.35  3.30 0.001  0.00  0.00
## 3                                                   gender  0.06 0.01 35598.09  4.77 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35706.37  9.28 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35657.25 -6.58 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35532.93 -0.40 0.687 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35539.44 -0.85 0.398 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35540.85  0.14 0.892 -0.16  0.19
## 9                                            occupManagers  0.00 0.09 35533.47  0.00 0.998 -0.18  0.18
## 10                            occupOther: Not in paid work  0.16 0.09 35696.52  1.71 0.088 -0.02  0.34
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35539.06 -0.77 0.442 -0.24  0.11
## 12                                      occupProfessionals  0.10 0.09 35535.19  1.10 0.273 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35538.37 -0.28 0.782 -0.22  0.17
## 14                          occupService and sales workers -0.04 0.09 35534.90 -0.40 0.691 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35544.87 -0.52 0.601 -0.23  0.14
## 16            occupTechnicians and associate professionals -0.03 0.09 35532.10 -0.29 0.774 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35546.28 -0.23 0.816 -0.24  0.19
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.59  0.51 0.611 -0.21  0.36
## 2                                                      age  0.00 0.00 34270.64  4.39 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35593.07  4.67 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35714.29  7.51 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35651.76 -5.23 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35528.94 -0.46 0.646 -0.21  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35535.31 -0.76 0.446 -0.24  0.11
## 8                              occupElementary occupations  0.01 0.09 35536.62  0.16 0.871 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35529.47 -0.18 0.859 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35691.43  1.54 0.123 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35534.89 -0.72 0.474 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35531.19  0.83 0.407 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35534.45 -0.34 0.734 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35530.91 -0.43 0.668 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35540.59 -0.58 0.564 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35528.10 -0.36 0.718 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35542.75 -0.19 0.852 -0.23  0.19
## 18                                            environ.lvl1  0.13 0.00 35459.04 26.76 0.000  0.12  0.13
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.47  0.53 0.598 -0.21  0.36
## 2                                                      age  0.00 0.00 34211.38  4.49 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35569.60  4.70 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35698.92  7.49 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35634.38 -5.24 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35508.21 -0.48 0.632 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35518.38 -0.77 0.439 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35520.28  0.13 0.900 -0.16  0.18
## 9                                            occupManagers -0.02 0.09 35507.43 -0.22 0.827 -0.19  0.15
## 10                            occupOther: Not in paid work  0.14 0.09 35674.75  1.52 0.127 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35516.07 -0.74 0.460 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35512.07  0.79 0.432 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35508.44 -0.36 0.716 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35513.92 -0.44 0.657 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35523.04 -0.63 0.529 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35508.10 -0.38 0.705 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35525.32 -0.24 0.812 -0.24  0.18
## 18                                            environ.lvl1  0.12 0.01    19.40 11.53 0.000  0.10  0.15
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
## Warning in cor(environ.gmc, refugees, use = "pairwise.complete.obs"): the standard deviation is zero

## Warning in cor(environ.gmc, refugees, use = "pairwise.complete.obs"): the standard deviation is zero
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
## Warning: Column `cntry` joining character vector and factor, coercing into character vector
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
## Warning: Column `voting.group` joining factor and character vector, coercing into character vector
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

# png-file for Figure 1

grDevices::png(file = "forestplot.png",family="sans",units = "mm",
               height=210,width=210,res=300) 

forest.m<-meta::forest(m,overall=T,
                       plotwidth="11cm",
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
##                                                                       COR            95%-CI %W(random)
## AT: BZÖ                                                            0.0000 [-0.6296; 0.6296]        0.3
## AT: FPÖ                                                            0.2587 [ 0.1373; 0.3724]        5.8
## BE: N-VA                                                          -0.0320 [-0.1533; 0.0902]        5.9
## BE: Parti Populaire                                               -0.1021 [-0.9030; 0.8574]        0.1
## BE: Vlaams Belang                                                  0.1874 [-0.1854; 0.5130]        1.2
## CH: Swiss People's Party                                           0.0096 [-0.1531; 0.1718]        4.4
## CZ: ODS                                                            0.1914 [-0.0042; 0.3729]        3.4
## CZ: Usvit                                                          0.3619 [ 0.0379; 0.6171]        1.5
## DE: AfD                                                           -0.0016 [-0.2598; 0.2568]        2.2
## DE: NPD                                                            0.1633 [-0.5197; 0.7190]        0.3
## EE: Eesti Konservatiivne Rahvaerakond                              0.0235 [-0.2127; 0.2570]        2.6
## ES: Partido Popular - PP                                           0.0959 [-0.0145; 0.2041]        6.5
## FI: True Finns                                                     0.1420 [-0.0007; 0.2791]        5.1
## FR: FN (Front National)                                            0.2846 [ 0.1110; 0.4413]        3.9
## FR: MPF (Mouvement pour la France)                                -0.2187 [-0.6029; 0.2478]        0.8
## FR: UMP (Union pour un Mouvement Populaire)                        0.1587 [ 0.0463; 0.2672]        6.3
## GB: Conservative                                                   0.2133 [ 0.1269; 0.2964]        7.6
## GB: Democratic Unionist Party (nir)                                0.1931 [-0.6552; 0.8260]        0.2
## GB: UK Independence Party                                          0.0790 [-0.1134; 0.2658]        3.6
## HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)  0.2420 [ 0.1583; 0.3223]        7.6
## HU: Jobbik (Jobbik Magyarországért Mozgalom)                       0.2087 [-0.0033; 0.4028]        3.1
## IT: Fratelli d'Italia                                              0.0080 [-0.3661; 0.3800]        1.1
## IT: Lega Nord                                                     -0.1013 [-0.3209; 0.1286]        2.8
## IT: Popolo delle Libertà (PdL)                                     0.2035 [ 0.0114; 0.3812]        3.5
## NL: Party for Freedom                                             -0.0299 [-0.2270; 0.1695]        3.4
## NL: Reformed Political Party                                       0.4091 [-0.0151; 0.7085]        0.9
## NO: Progress Party (FRP)                                           0.1942 [ 0.0163; 0.3602]        3.9
## PL: KORWIN                                                        -0.1311 [-0.4493; 0.2167]        1.4
## PL: Kukiz'15                                                       0.0307 [-0.1677; 0.2268]        3.4
## SE: Sverigedemokraterna                                            0.0849 [-0.1221; 0.2847]        3.2
## SI: SDS - Slovenska demokratska stranka                            0.1381 [-0.0343; 0.3024]        4.1
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
##                                                        COR            95%-CI %W(random)
## AT: Grüne                                           0.2685 [ 0.1276; 0.3987]        6.6
## BE: Ecolo                                           0.1091 [-0.1489; 0.3533]        3.7
## BE: Groen!                                         -0.0325 [-0.2532; 0.1914]        4.4
## CH: Green Party                                     0.0892 [-0.1601; 0.3278]        3.9
## CH: Social Democratic Party                         0.3476 [ 0.1866; 0.4904]        5.8
## DE: Bündnis 90/ Die Grünen                          0.1805 [ 0.0545; 0.3007]        7.2
## EE: Erakond Eestimaa Rohelised                      0.4930 [-0.1128; 0.8316]        0.8
## FI: Green League                                    0.1873 [ 0.0351; 0.3311]        6.3
## FR: Autres mouvements écologistes                  -0.0066 [-0.4478; 0.4372]        1.5
## FR: EELV (Europe Ecologie Les Verts)               -0.0099 [-0.2632; 0.2446]        3.7
## GB: Green Party                                     0.3931 [ 0.0840; 0.6332]        2.7
## HU: LMP (Lehet Más A Politika)                     -0.0781 [-0.5258; 0.4035]        1.3
## IE: Green Party                                     0.1208 [-0.2797; 0.4854]        1.9
## IT: Movimento 5 Stelle                              0.3121 [ 0.1884; 0.4260]        7.1
## IT: Sinistra Ecologia e Libertà (SEL)               0.3445 [-0.0586; 0.6510]        1.9
## LT: Lithuanian Greens Party (LZP)                   0.3536 [-0.4676; 0.8472]        0.5
## LT: Lithuanian Peasant and Greens Union (LVZS)      0.1456 [ 0.0199; 0.2668]        7.2
## NL: Green Left                                      0.0306 [-0.2229; 0.2803]        3.8
## NL: Party for the Animals                          -0.2584 [-0.5660; 0.1123]        2.2
## NO: Green Party (MDG)                              -0.0569 [-0.3876; 0.2868]        2.4
## NO: Liberal Party (V)                               0.0484 [-0.2303; 0.3197]        3.3
## NO: Socialist Left Party (SV)                       0.1245 [-0.1336; 0.3669]        3.7
## PT: B.E. - Bloco de Esquerda                        0.1221 [-0.1338; 0.3628]        3.8
## PT: PAN - Pessoas-Animais-Natureza                 -0.0356 [-0.6222; 0.5766]        0.8
## SE: FI (Feministiskt initiativ)                     0.3334 [-0.0453; 0.6283]        2.1
## SE: Miljöpartiet de gröna                           0.4540 [ 0.2779; 0.6006]        4.9
## SE: Vänsterpartiet                                  0.4020 [ 0.1985; 0.5722]        4.4
## SI: ZL - Združena levica (DSD, IDS in Stranka TRS)  0.0579 [-0.3294; 0.4285]        2.0
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

##### Confidence intervals for supplementary tables and Table 1


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
##                                                                 group observed_r    LL   UL     p    r            CI  pno0
## 1                                                             AT: BZÖ       0.00 -0.63 0.63 1.000  .00   [-.63, .63] 1.000
## 4                                                             AT: FPÖ       0.26  0.14 0.37 0.000  .26    [.14, .37]  .000
## 22                                                           BE: N-VA      -0.03 -0.15 0.09 0.608 -.03   [-.15, .09]  .608
## 28                                                BE: Parti Populaire      -0.10 -0.90 0.86 0.885 -.10   [-.90, .86]  .885
## 34                                                  BE: Vlaams Belang       0.19 -0.19 0.51 0.325  .19   [-.19, .51]  .325
## 42                                       CH: Federal Democratic Union       0.87 -1.00 1.00 1.000  .87 [-1.00, 1.00] 1.000
## 51                                           CH: Swiss People's Party       0.01 -0.15 0.17 0.909  .01   [-.15, .17]  .909
## 52                                                  CH: Ticino League         NA    NA   NA    NA <NA>      [NA, NA]  <NA>
## 62                                                            CZ: ODS       0.19  0.00 0.37 0.055  .19    [.00, .37]  .055
## 66                                                          CZ: Usvit       0.36  0.04 0.62 0.029  .36    [.04, .62]  .029
## 67                                                            DE: AfD       0.00 -0.26 0.26 0.990  .00   [-.26, .26]  .990
## 77                                                            DE: NPD       0.16 -0.52 0.72 0.663  .16   [-.52, .72]  .663
## 85                              EE: Eesti Konservatiivne Rahvaerakond       0.02 -0.21 0.26 0.848  .02   [-.21, .26]  .848
## 111                                          ES: Partido Popular - PP       0.10 -0.01 0.20 0.089  .10   [-.01, .20]  .089
## 129                                                    FI: True Finns       0.14  0.00 0.28 0.051  .14    [.00, .28]  .051
## 136                                           FR: FN (Front National)       0.28  0.11 0.44 0.002  .28    [.11, .44]  .002
## 138                                FR: MPF (Mouvement pour la France)      -0.22 -0.60 0.25 0.359 -.22   [-.60, .25]  .359
## 148                       FR: UMP (Union pour un Mouvement Populaire)       0.16  0.05 0.27 0.006  .16    [.05, .27]  .006
## 149                                                  GB: Conservative       0.21  0.13 0.30 0.000  .21    [.13, .30]  .000
## 150                               GB: Democratic Unionist Party (nir)       0.19 -0.66 0.83 0.696  .19   [-.66, .83]  .696
## 163                                         GB: UK Independence Party       0.08 -0.11 0.27 0.422  .08   [-.11, .27]  .422
## 166 HU: Fidesz - KDNP (Fidesz – Magyar Polgári Szövetség Keresztényd)       0.24  0.16 0.32 0.000  .24    [.16, .32]  .000
## 167                      HU: Jobbik (Jobbik Magyarországért Mozgalom)       0.21  0.00 0.40 0.054  .21    [.00, .40]  .054
## 190                                             IT: Fratelli d'Italia       0.01 -0.37 0.38 0.968  .01   [-.37, .38]  .968
## 192                                                     IT: Lega Nord      -0.10 -0.32 0.13 0.388 -.10   [-.32, .13]  .388
## 199                                    IT: Popolo delle Libertà (PdL)       0.20  0.01 0.38 0.038  .20    [.01, .38]  .038
## 234                                             NL: Party for Freedom      -0.03 -0.23 0.17 0.771 -.03   [-.23, .17]  .771
## 237                                      NL: Reformed Political Party       0.41 -0.02 0.71 0.058  .41   [-.02, .71]  .058
## 252                                          NO: Progress Party (FRP)       0.19  0.02 0.36 0.033  .19    [.02, .36]  .033
## 258                                                        PL: KORWIN      -0.13 -0.45 0.22 0.463 -.13   [-.45, .22]  .463
## 259                                                      PL: Kukiz'15       0.03 -0.17 0.23 0.763  .03   [-.17, .23]  .763
## 298                                           SE: Sverigedemokraterna       0.08 -0.12 0.28 0.422  .08   [-.12, .28]  .422
## 311                           SI: SDS - Slovenska demokratska stranka       0.14 -0.03 0.30 0.116  .14   [-.03, .30]  .116
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
##                                                  group observed_r    LL   UL     p    r            CI  pno0
## 5                                            AT: Grüne       0.27  0.13 0.40 0.000  .27    [.13, .40]  .000
## 19                                           BE: Ecolo       0.11 -0.15 0.35 0.408  .11   [-.15, .35]  .408
## 20                                          BE: Groen!      -0.03 -0.25 0.19 0.778 -.03   [-.25, .19]  .778
## 44                                     CH: Green Party       0.09 -0.16 0.33 0.485  .09   [-.16, .33]  .485
## 50                         CH: Social Democratic Party       0.35  0.19 0.49 0.000  .35    [.19, .49]  .000
## 68                          DE: Bündnis 90/ Die Grünen       0.18  0.05 0.30 0.005  .18    [.05, .30]  .005
## 88                      EE: Erakond Eestimaa Rohelised       0.49 -0.11 0.83 0.105  .49   [-.11, .83]  .105
## 118                                   FI: Green League       0.19  0.04 0.33 0.016  .19    [.04, .33]  .016
## 130                  FR: Autres mouvements écologistes      -0.01 -0.45 0.44 0.978 -.01   [-.45, .44]  .978
## 134               FR: EELV (Europe Ecologie Les Verts)      -0.01 -0.26 0.24 0.940 -.01   [-.26, .24]  .940
## 153                                    GB: Green Party       0.39  0.08 0.63 0.014  .39    [.08, .63]  .014
## 168                     HU: LMP (Lehet Más A Politika)      -0.08 -0.53 0.40 0.762 -.08   [-.53, .40]  .762
## 180                                    IE: Green Party       0.12 -0.28 0.49 0.561  .12   [-.28, .49]  .561
## 193                             IT: Movimento 5 Stelle       0.31  0.19 0.43 0.000  .31    [.19, .43]  .000
## 201                   IT: Rivoluzione Civile (Ingroia)       0.76 -1.00 1.00 1.000  .76 [-1.00, 1.00] 1.000
## 203              IT: Sinistra Ecologia e Libertà (SEL)       0.34 -0.06 0.65 0.092  .34   [-.06, .65]  .092
## 211                  LT: Lithuanian Greens Party (LZP)       0.35 -0.47 0.85 0.409  .35   [-.47, .85]  .409
## 212     LT: Lithuanian Peasant and Greens Union (LVZS)       0.15  0.02 0.27 0.023  .15    [.02, .27]  .023
## 227                                     NL: Green Left       0.03 -0.22 0.28 0.815  .03   [-.22, .28]  .815
## 235                          NL: Party for the Animals      -0.26 -0.57 0.11 0.169 -.26   [-.57, .11]  .169
## 245                              NO: Green Party (MDG)      -0.06 -0.39 0.29 0.751 -.06   [-.39, .29]  .751
## 247                              NO: Liberal Party (V)       0.05 -0.23 0.32 0.737  .05   [-.23, .32]  .737
## 254                      NO: Socialist Left Party (SV)       0.12 -0.13 0.37 0.345  .12   [-.13, .37]  .345
## 269                       PT: B.E. - Bloco de Esquerda       0.12 -0.13 0.36 0.350  .12   [-.13, .36]  .350
## 273                         PT: MPT - Partido da Terra         NA    NA   NA    NA <NA>      [NA, NA]  <NA>
## 278                 PT: PAN - Pessoas-Animais-Natureza      -0.04 -0.62 0.58 0.920 -.04   [-.62, .58]  .920
## 286                    SE: FI (Feministiskt initiativ)       0.33 -0.05 0.63 0.083  .33   [-.05, .63]  .083
## 289                          SE: Miljöpartiet de gröna       0.45  0.28 0.60 0.000  .45    [.28, .60]  .000
## 299                                 SE: Vänsterpartiet       0.40  0.20 0.57 0.000  .40    [.20, .57]  .000
## 314 SI: ZL - Združena levica (DSD, IDS in Stranka TRS)       0.06 -0.33 0.43 0.777  .06   [-.33, .43]  .777
```

```r
write.csv2(pro.env.r,"pro.env.r.csv")
      




# For table 1


tab1.mod2<-read.csv2("FE.H1.mod1.csv",
                   stringsAsFactors = F)       

tab1.mod2<-data.frame(Eff=tab1.mod2$Eff,
                      Est=tab1.mod2$Est,
                      CI=paste0("[",tab1.mod2[,"LL"],", ",tab1.mod2[,"UL"],"]"),
                      p=tab1.mod2$p)

tab1.mod3<-read.csv2("FE.H1.mod2.csv",
                     stringsAsFactors = F)

tab1.mod3<-data.frame(Eff=tab1.mod3$Eff,
                      Est=tab1.mod3$Est,
                      CI=paste0("[",tab1.mod3[,"LL"],", ",tab1.mod3[,"UL"],"]"),
                      p=tab1.mod3$p)

tab1.mod4<-read.csv2("FE.H1.mod3.csv",
                     stringsAsFactors = F)

tab1.mod4<-data.frame(Eff=tab1.mod4$Eff,
                      Est=tab1.mod4$Est,
                      CI=paste0("[",tab1.mod4[,"LL"],", ",tab1.mod4[,"UL"],"]"),
                      p=tab1.mod4$p)


tab1<-full_join(x=tab1.mod2,
          y=tab1.mod3,
          by="Eff")
```

```
## Warning: Column `Eff` joining factors with different levels, coercing to character vector
```

```r
tab1<-full_join(x=tab1,
                y=tab1.mod4,
                by="Eff")
```

```
## Warning: Column `Eff` joining character vector and factor, coercing into character vector
```

```r
tab1
```

```
##                                                        Eff Est.x           CI.x   p.x Est.y           CI.y   p.y   Est
## 1                                              (Intercept)  0.06  [-0.22, 0.35] 0.663  0.07  [-0.21, 0.36] 0.611  0.08
## 2                                                      age  0.00   [0.00, 0.00] 0.001  0.00   [0.00, 0.00] 0.000  0.00
## 3                                                   gender  0.06   [0.03, 0.08] 0.000  0.06   [0.03, 0.08] 0.000  0.06
## 4                                                     educ  0.02   [0.01, 0.02] 0.000  0.01   [0.01, 0.02] 0.000  0.01
## 5                                                    resid -0.08 [-0.10, -0.05] 0.000 -0.06 [-0.08, -0.04] 0.000 -0.06
## 6                            occupClerical support workers -0.04  [-0.21, 0.14] 0.687 -0.04  [-0.21, 0.13] 0.646 -0.04
## 7                    occupCraft and related trades workers -0.08  [-0.25, 0.10] 0.398 -0.07  [-0.24, 0.11] 0.446 -0.07
## 8                              occupElementary occupations  0.01  [-0.16, 0.19] 0.892  0.01  [-0.16, 0.19] 0.871  0.01
## 9                                            occupManagers  0.00  [-0.18, 0.18] 0.998 -0.02  [-0.19, 0.16] 0.859 -0.02
## 10                            occupOther: Not in paid work  0.16  [-0.02, 0.34] 0.088  0.14  [-0.04, 0.32] 0.123  0.14
## 11        occupPlant and machine operators, and assemblers -0.07  [-0.24, 0.11] 0.442 -0.06  [-0.24, 0.11] 0.474 -0.07
## 12                                      occupProfessionals  0.10  [-0.08, 0.27] 0.273  0.07  [-0.10, 0.24] 0.407  0.07
## 13                                            occupRetired -0.03  [-0.22, 0.17] 0.782 -0.03  [-0.23, 0.16] 0.734 -0.04
## 14                          occupService and sales workers -0.04  [-0.21, 0.14] 0.691 -0.04  [-0.21, 0.13] 0.668 -0.04
## 15 occupSkilled agricultural, forestry and fishery workers -0.05  [-0.23, 0.14] 0.601 -0.05  [-0.24, 0.13] 0.564 -0.06
## 16            occupTechnicians and associate professionals -0.03  [-0.20, 0.15] 0.774 -0.03  [-0.20, 0.14] 0.718 -0.03
## 17                                         occupUnemployed -0.03  [-0.24, 0.19] 0.816 -0.02  [-0.23, 0.19] 0.852 -0.03
## 18                                            environ.lvl1  <NA>           <NA>  <NA>  0.13   [0.12, 0.13] 0.000  0.12
##                CI     p
## 1   [-0.21, 0.36] 0.598
## 2    [0.00, 0.00] 0.000
## 3    [0.03, 0.08] 0.000
## 4    [0.01, 0.02] 0.000
## 5  [-0.08, -0.04] 0.000
## 6   [-0.22, 0.13] 0.632
## 7   [-0.24, 0.10] 0.439
## 8   [-0.16, 0.18] 0.900
## 9   [-0.19, 0.15] 0.827
## 10  [-0.04, 0.32] 0.127
## 11  [-0.24, 0.11] 0.460
## 12  [-0.10, 0.24] 0.432
## 13  [-0.23, 0.16] 0.716
## 14  [-0.21, 0.13] 0.657
## 15  [-0.24, 0.12] 0.529
## 16  [-0.20, 0.14] 0.705
## 17  [-0.24, 0.18] 0.812
## 18   [0.10, 0.15] 0.000
```

```r
write.csv2(tab1,"tab1.csv")

confint.merMod(H1.mod1,parm="theta_",method="profile",level=.95,oldNames=F)
```

```
## Computing profile confidence intervals ...
```

```
##                                 2.5 %    97.5 %
## sd_(Intercept)|voting.group 0.2739203 0.3370940
## sd_(Intercept)|cntry        0.3716660 0.7102682
## sigma                       1.0341091 1.0494469
```

```r
confint.merMod(H1.mod2,parm="theta_",method="profile",level=.95,oldNames=F)
```

```
## Computing profile confidence intervals ...
```

```
##                                 2.5 %    97.5 %
## sd_(Intercept)|voting.group 0.2790412 0.3428484
## sd_(Intercept)|cntry        0.3719003 0.7110853
## sigma                       1.0237206 1.0389048
```

```r
#don't run the random slopes -model again, CIs below
#confint.merMod(H1.mod3,parm="theta_",method="profile",level=.95,oldNames=F)

#d_(Intercept)|voting.group                0.2793160 0.34313941
#cor_environ.lvl1.(Intercept)|voting.group -0.1898851 0.44957436
#sd_environ.lvl1|voting.group               0.0242823 0.06094131
#sd_(Intercept)|cntry                       0.3720908 0.71144406
#cor_environ.lvl1.(Intercept)|cntry        -0.5288256 0.45582915
#sd_environ.lvl1|cntry                      0.0242841 0.06260749
#sigma                                      1.0215694 1.03677305
```

\newpage

## Alternative (exploratory) approach for Hypothesis 1 with Environment attitudes as dependent variable, and immigrant attitudes as independent

### Model 0: Intercepts only


```r
H1.env.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##           Eff  Est   SE    df    t     p    LL   UL
## 1 (Intercept) 0.05 0.07 20.08 0.78 0.443 -0.09 0.20
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.04 0.12   188.50  -0.32 0.751 -0.27  0.20
## 2                                                      age -0.00 0.00 29829.87  -8.11 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35659.58   1.14 0.253 -0.01  0.04
## 4                                                     educ  0.03 0.00 35342.84  13.39 0.000  0.02  0.03
## 5                                                    resid -0.14 0.01 35731.67 -10.19 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35578.76   0.42 0.672 -0.15  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35588.19  -0.60 0.548 -0.26  0.14
## 8                              occupElementary occupations -0.01 0.10 35591.13  -0.11 0.915 -0.21  0.19
## 9                                            occupManagers  0.13 0.10 35579.35   1.33 0.184 -0.06  0.33
## 10                            occupOther: Not in paid work  0.16 0.10 35719.64   1.55 0.121 -0.04  0.36
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35588.53  -0.41 0.684 -0.24  0.16
## 12                                      occupProfessionals  0.21 0.10 35581.91   2.09 0.037  0.01  0.40
## 13                                            occupRetired  0.06 0.11 35582.94   0.52 0.604 -0.16  0.28
## 14                          occupService and sales workers  0.02 0.10 35580.41   0.25 0.806 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35595.50   0.43 0.664 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35577.81   0.59 0.557 -0.14  0.25
## 17                                         occupUnemployed -0.03 0.12 35586.90  -0.25 0.803 -0.27  0.21
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.04 0.12   182.49 -0.30 0.763 -0.27  0.20
## 2                                                      age -0.00 0.00 30300.28 -8.71 0.000 -0.00 -0.00
## 3                                                   gender  0.01 0.01 35654.32  0.49 0.622 -0.02  0.03
## 4                                                     educ  0.03 0.00 35392.50 12.22 0.000  0.02  0.03
## 5                                                    resid -0.12 0.01 35728.36 -9.40 0.000 -0.15 -0.10
## 6                            occupClerical support workers  0.05 0.10 35574.16  0.48 0.629 -0.15  0.24
## 7                    occupCraft and related trades workers -0.05 0.10 35583.40 -0.49 0.625 -0.24  0.15
## 8                              occupElementary occupations -0.01 0.10 35586.20 -0.12 0.903 -0.21  0.18
## 9                                            occupManagers  0.13 0.10 35574.73  1.34 0.181 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35722.40  1.35 0.178 -0.06  0.34
## 11        occupPlant and machine operators, and assemblers -0.03 0.10 35583.69 -0.30 0.762 -0.23  0.17
## 12                                      occupProfessionals  0.19 0.10 35577.21  1.96 0.050 -0.00  0.38
## 13                                            occupRetired  0.06 0.11 35578.52  0.56 0.573 -0.15  0.28
## 14                          occupService and sales workers  0.03 0.10 35575.86  0.31 0.760 -0.16  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.10 35590.53  0.51 0.610 -0.15  0.26
## 16            occupTechnicians and associate professionals  0.06 0.10 35573.19  0.63 0.527 -0.13  0.25
## 17                                         occupUnemployed -0.02 0.12 35582.73 -0.20 0.838 -0.26  0.21
## 18                                           refugees.lvl1  0.16 0.01 35453.62 26.78 0.000  0.15  0.17
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.04 0.12   181.41 -0.33 0.740 -0.27  0.19
## 2                                                      age -0.00 0.00 30229.24 -8.64 0.000 -0.00 -0.00
## 3                                                   gender  0.01 0.01 35647.27  0.51 0.613 -0.02  0.03
## 4                                                     educ  0.03 0.00 35243.03 12.16 0.000  0.02  0.03
## 5                                                    resid -0.12 0.01 35688.34 -9.42 0.000 -0.15 -0.10
## 6                            occupClerical support workers  0.05 0.10 35558.96  0.52 0.604 -0.14  0.25
## 7                    occupCraft and related trades workers -0.04 0.10 35566.00 -0.45 0.652 -0.24  0.15
## 8                              occupElementary occupations -0.01 0.10 35571.70 -0.10 0.921 -0.20  0.19
## 9                                            occupManagers  0.14 0.10 35560.44  1.37 0.170 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35708.21  1.40 0.163 -0.06  0.34
## 11        occupPlant and machine operators, and assemblers -0.03 0.10 35566.12 -0.27 0.788 -0.22  0.17
## 12                                      occupProfessionals  0.19 0.10 35558.85  1.99 0.047  0.00  0.39
## 13                                            occupRetired  0.06 0.11 35554.14  0.58 0.560 -0.15  0.28
## 14                          occupService and sales workers  0.03 0.10 35560.22  0.34 0.735 -0.16  0.23
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.10 35572.65  0.50 0.615 -0.15  0.26
## 16            occupTechnicians and associate professionals  0.07 0.10 35559.38  0.69 0.489 -0.12  0.26
## 17                                         occupUnemployed -0.02 0.12 35571.95 -0.19 0.848 -0.26  0.21
## 18                                           refugees.lvl1  0.16 0.01    18.18 13.89 0.000  0.13  0.18
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



\newpage



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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.06 0.14    50.42  0.44 0.663 -0.22  0.35
## 2                                                      age  0.00 0.00 34095.35  3.30 0.001  0.00  0.00
## 3                                                   gender  0.06 0.01 35598.09  4.77 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35706.37  9.28 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35657.25 -6.58 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35532.93 -0.40 0.687 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35539.44 -0.85 0.398 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35540.85  0.14 0.892 -0.16  0.19
## 9                                            occupManagers  0.00 0.09 35533.47  0.00 0.998 -0.18  0.18
## 10                            occupOther: Not in paid work  0.16 0.09 35696.52  1.71 0.088 -0.02  0.34
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35539.06 -0.77 0.442 -0.24  0.11
## 12                                      occupProfessionals  0.10 0.09 35535.19  1.10 0.273 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35538.37 -0.28 0.782 -0.22  0.17
## 14                          occupService and sales workers -0.04 0.09 35534.90 -0.40 0.691 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35544.87 -0.52 0.601 -0.23  0.14
## 16            occupTechnicians and associate professionals -0.03 0.09 35532.10 -0.29 0.774 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35546.28 -0.23 0.816 -0.24  0.19
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.15    62.98 -3.25 0.002 -0.76 -0.18
## 2                                                      age  0.00 0.00 35716.28  3.73 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35620.34  4.74 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35722.81  9.31 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35707.05 -6.49 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35593.68 -0.44 0.661 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66 -0.86 0.390 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35596.76  0.09 0.928 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35595.21 -0.03 0.978 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56  1.49 0.135 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23 -0.80 0.426 -0.25  0.10
## 12                                      occupProfessionals  0.09 0.09 35595.08  1.06 0.290 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35590.71 -0.33 0.741 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35594.22 -0.42 0.677 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85 -0.56 0.573 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11 -0.31 0.753 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35590.31 -0.30 0.763 -0.25  0.18
## 18                            all.parties.lvl2Did not vote  0.45 0.06   181.24  7.15 0.000  0.33  0.57
## 19                              all.parties.lvl2Don't know  0.42 0.07   270.53  6.02 0.000  0.28  0.56
## 20                            all.parties.lvl2Invalid vote  0.50 0.36  1156.29  1.40 0.161 -0.20  1.20
## 21                                  all.parties.lvl2NE age  0.73 0.07   274.23 10.42 0.000  0.59  0.87
## 22                              all.parties.lvl2NE citizen  0.85 0.08   268.85 11.25 0.000  0.70  1.00
## 23                                all.parties.lvl2NE other  0.71 0.10   672.84  7.21 0.000  0.51  0.90
## 24                               all.parties.lvl2No answer  0.56 0.38  1449.94  1.49 0.136 -0.18  1.30
## 25                             all.parties.lvl2Other party  0.54 0.05   230.15 11.29 0.000  0.45  0.63
## 26                   all.parties.lvl2Pro-environment party  0.90 0.06   256.45 13.92 0.000  0.77  1.02
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
##                     group emmean             CI     p p.adj p_less_001 p.adj_less_001
## 1  Anti-immigration party  -0.48 [-0.71, -0.25] 0.000 0.000        yes            yes
## 2            Did not vote  -0.03  [-0.26, 0.20] 0.781 1.000         no             no
## 3              Don't know  -0.06  [-0.30, 0.18] 0.630 1.000         no             no
## 4            Invalid vote   0.02  [-0.71, 0.74] 0.963 1.000         no             no
## 5                  NE age   0.25   [0.01, 0.49] 0.040 0.278         no             no
## 6              NE citizen   0.37   [0.13, 0.61] 0.003 0.024         no             no
## 7                NE other   0.22  [-0.05, 0.50] 0.109 0.651         no             no
## 8               No answer   0.08  [-0.68, 0.84] 0.839 1.000         no             no
## 9             Other party   0.06  [-0.16, 0.27] 0.596 1.000         no             no
## 10  Pro-environment party   0.41   [0.18, 0.65] 0.000 0.004        yes             no
```

```r
write.csv2(H2.mod2.mmeans.tab,"H2.mod2.mmeans.tab.csv")


#contrast between anti-immigration and pro-environment
(H2.contrast<-data.frame(pairs(H2.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df  z.ratio      p.value
## 1 Pro-environment party - Anti-immigration party 0.8955111 0.06434551 Inf 13.91723 4.978424e-44
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
##                                         contrast  estimate         SE  df   z.ratio      p.value
## 1           Other party - Anti-immigration party 0.5398443 0.04779971 Inf 11.293884 2.813585e-29
## 2 Pro-environment party - Anti-immigration party 0.8955111 0.06434551 Inf 13.917228 1.493527e-43
## 3            Pro-environment party - Other party 0.3556668 0.05168896 Inf  6.880904 5.947390e-12
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
## Warning in Ops.factor(H2.mod2.mmeans.tab[10, 2], H2.mod2.mmeans.tab[1, 2]): '-' not meaningful for factors
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
## Warning in Ops.factor(H2.mod2.mmeans.tab[10, 2], H2.mod2.mmeans.tab[9, 2]): '-' not meaningful for factors
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
## Warning in Ops.factor(H2.mod2.mmeans.tab[9, 2], H2.mod2.mmeans.tab[1, 2]): '-' not meaningful for factors
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.02 0.15    64.39 -0.16 0.875 -0.31  0.27
## 2                                                      age  0.00 0.00 35716.28  3.73 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35620.34  4.74 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35722.81  9.31 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35707.05 -6.49 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35593.68 -0.44 0.661 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66 -0.86 0.390 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35596.76  0.09 0.928 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35595.21 -0.03 0.978 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56  1.49 0.135 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23 -0.80 0.426 -0.25  0.10
## 12                                      occupProfessionals  0.09 0.09 35595.08  1.06 0.290 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35590.71 -0.33 0.741 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35594.23 -0.42 0.677 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85 -0.56 0.573 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11 -0.31 0.753 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35590.31 -0.30 0.763 -0.25  0.18
## 18                                       other.party.dummy  0.09 0.05   163.48  1.82 0.071 -0.01  0.19
## 19                                         dont.know.dummy -0.03 0.07   222.65 -0.36 0.717 -0.17  0.12
## 20                                      invalid.vote.dummy  0.05 0.36  1136.40  0.14 0.889 -0.65  0.75
## 21                                         no.answer.dummy  0.11 0.38  1416.47  0.30 0.767 -0.63  0.85
## 22                                  not.eligible.age.dummy  0.28 0.07   219.09  3.96 0.000  0.14  0.42
## 23                          not.eligible.citizenship.dummy  0.40 0.08   226.00  5.23 0.000  0.25  0.55
## 24                                not.eligible.other.dummy  0.26 0.10   558.14  2.59 0.010  0.06  0.45
## 25                                    anti.imm.party.dummy -0.45 0.06   181.24 -7.15 0.000 -0.57 -0.33
## 26                                     pro.env.party.dummy  0.45 0.07   206.58  6.76 0.000  0.32  0.58
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.02 0.15    63.03 -0.16 0.873 -0.31  0.27
## 2                                                      age  0.00 0.00 35701.12  3.74 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35620.05  4.74 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35718.65  9.31 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35708.16 -6.51 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35595.60 -0.44 0.661 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35602.53 -0.86 0.391 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35598.49  0.09 0.926 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35596.94 -0.03 0.976 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35647.42  1.49 0.136 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35602.28 -0.79 0.427 -0.25  0.10
## 12                                      occupProfessionals  0.09 0.09 35596.84  1.05 0.292 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35592.09 -0.33 0.741 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35595.92 -0.42 0.678 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35612.90 -0.56 0.573 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35594.87 -0.32 0.753 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35592.09 -0.29 0.768 -0.24  0.18
## 18                                       other.party.dummy  0.09 0.05   132.00  1.89 0.061 -0.00  0.18
## 19                                         dont.know.dummy -0.03 0.07   184.46 -0.38 0.701 -0.16  0.11
## 20                                      invalid.vote.dummy  0.05 0.35  1046.75  0.15 0.883 -0.63  0.74
## 21                                         no.answer.dummy  0.10 0.37  1321.51  0.26 0.796 -0.63  0.82
## 22                                  not.eligible.age.dummy  0.28 0.07   181.50  4.12 0.000  0.15  0.42
## 23                          not.eligible.citizenship.dummy  0.40 0.07   186.15  5.44 0.000  0.26  0.55
## 24                                not.eligible.other.dummy  0.25 0.10   483.87  2.64 0.009  0.06  0.44
## 25                                    anti.imm.party.dummy -0.45 0.07    38.57 -6.83 0.000 -0.58 -0.32
## 26                                     pro.env.party.dummy  0.44 0.08    30.86  5.76 0.000  0.28  0.59
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.03 0.14    51.13  0.23 0.823 -0.26  0.32
## 2                                                      age  0.00 0.00 35707.89  3.75 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35578.58  4.83 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35682.00  9.46 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35640.73 -6.49 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35542.31 -0.40 0.686 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35548.02 -0.86 0.391 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35545.03  0.10 0.920 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35543.00 -0.01 0.992 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35580.80  1.52 0.129 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35546.89 -0.78 0.433 -0.25  0.11
## 12                                      occupProfessionals  0.10 0.09 35544.96  1.09 0.276 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35542.87 -0.31 0.755 -0.22  0.16
## 14                          occupService and sales workers -0.04 0.09 35543.35 -0.42 0.677 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35554.44 -0.53 0.596 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35541.81 -0.29 0.772 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35546.38 -0.29 0.774 -0.24  0.18
## 18                                      did.not.vote.dummy -0.06 0.07   197.35 -0.80 0.425 -0.19  0.08
## 19                                         dont.know.dummy -0.08 0.08   283.64 -1.10 0.273 -0.23  0.07
## 20                                      invalid.vote.dummy  0.01 0.41   684.76  0.01 0.989 -0.81  0.82
## 21                                         no.answer.dummy  0.04 0.43   809.06  0.08 0.933 -0.81  0.88
## 22                                  not.eligible.age.dummy  0.23 0.08   286.48  3.01 0.003  0.08  0.38
## 23                          not.eligible.citizenship.dummy  0.35 0.08   303.31  4.09 0.000  0.18  0.51
## 24                                not.eligible.other.dummy  0.22 0.11   711.72  2.06 0.040  0.01  0.43
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept)  0.12 0.14    53.92   0.86 0.392 -0.16  0.40
## 2                                                      age  0.00 0.00 35722.27   3.59 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35608.88   4.79 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35718.92   9.43 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35691.03  -6.52 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35576.83  -0.41 0.680 -0.21  0.14
## 7                    occupCraft and related trades workers -0.07 0.09 35584.02  -0.83 0.406 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35580.21   0.11 0.910 -0.17  0.19
## 9                                            occupManagers -0.00 0.09 35577.58  -0.02 0.988 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35625.80   1.53 0.127 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35583.03  -0.77 0.440 -0.25  0.11
## 12                                      occupProfessionals  0.10 0.09 35579.90   1.09 0.274 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35574.87  -0.30 0.761 -0.22  0.16
## 14                          occupService and sales workers -0.04 0.09 35577.46  -0.40 0.690 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35591.80  -0.54 0.590 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35576.48  -0.29 0.774 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35575.34  -0.28 0.781 -0.24  0.18
## 18                                      did.not.vote.dummy -0.15 0.05   170.98  -2.67 0.008 -0.25 -0.04
## 19                                         dont.know.dummy -0.17 0.06   287.04  -2.75 0.006 -0.30 -0.05
## 20                                      invalid.vote.dummy -0.06 0.37   987.47  -0.17 0.864 -0.78  0.66
## 21                                         no.answer.dummy -0.04 0.39  1213.22  -0.11 0.910 -0.80  0.72
## 22                                  not.eligible.age.dummy  0.13 0.06   295.41   2.13 0.034  0.01  0.26
## 23                          not.eligible.citizenship.dummy  0.25 0.07   288.85   3.64 0.000  0.12  0.39
## 24                                not.eligible.other.dummy  0.11 0.09   833.42   1.21 0.228 -0.07  0.30
## 25                                    anti.imm.party.dummy -0.60 0.05   235.80 -11.61 0.000 -0.70 -0.49
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.03 0.14    52.30 -0.19 0.851 -0.31  0.26
## 2                                                      age  0.00 0.00 35717.46  3.88 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35589.38  4.79 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35687.43  9.36 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35660.33 -6.46 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35558.32 -0.43 0.669 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35564.26 -0.88 0.377 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35560.98  0.08 0.936 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35559.65 -0.02 0.984 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35599.92  1.49 0.136 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35563.39 -0.81 0.420 -0.25  0.10
## 12                                      occupProfessionals  0.09 0.09 35559.75  1.06 0.290 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35557.95 -0.34 0.736 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35559.41 -0.43 0.665 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35571.61 -0.55 0.582 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35557.67 -0.31 0.754 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35559.99 -0.31 0.757 -0.25  0.18
## 18                                      did.not.vote.dummy  0.01 0.06   189.71  0.09 0.930 -0.12  0.13
## 19                                         dont.know.dummy -0.02 0.07   289.32 -0.31 0.756 -0.16  0.12
## 20                                      invalid.vote.dummy  0.03 0.39   804.27  0.07 0.948 -0.74  0.79
## 21                                         no.answer.dummy  0.11 0.41   965.65  0.27 0.790 -0.70  0.92
## 22                                  not.eligible.age.dummy  0.29 0.07   293.85  4.13 0.000  0.15  0.43
## 23                          not.eligible.citizenship.dummy  0.41 0.08   302.50  5.23 0.000  0.26  0.56
## 24                                not.eligible.other.dummy  0.28 0.10   768.38  2.74 0.006  0.08  0.47
## 25                                     pro.env.party.dummy  0.45 0.06   285.94  7.32 0.000  0.33  0.57
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    54.74   0.49 0.629 -0.21  0.35
## 2                                                      age  0.00 0.00 35716.28   3.73 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35620.34   4.74 0.000  0.03  0.08
## 4                                                     educ  0.02 0.00 35722.81   9.31 0.000  0.01  0.02
## 5                                                    resid -0.08 0.01 35707.05  -6.49 0.000 -0.10 -0.05
## 6                            occupClerical support workers -0.04 0.09 35593.68  -0.44 0.661 -0.21  0.14
## 7                    occupCraft and related trades workers -0.08 0.09 35600.66  -0.86 0.390 -0.25  0.10
## 8                              occupElementary occupations  0.01 0.09 35596.75   0.09 0.928 -0.17  0.18
## 9                                            occupManagers -0.00 0.09 35595.21  -0.03 0.978 -0.18  0.17
## 10                            occupOther: Not in paid work  0.14 0.09 35644.56   1.49 0.135 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35600.23  -0.80 0.426 -0.25  0.10
## 12                                      occupProfessionals  0.09 0.09 35595.08   1.06 0.290 -0.08  0.27
## 13                                            occupRetired -0.03 0.10 35590.71  -0.33 0.741 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35594.22  -0.42 0.677 -0.21  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35609.85  -0.56 0.573 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35593.11  -0.31 0.753 -0.20  0.15
## 17                                         occupUnemployed -0.03 0.11 35590.31  -0.30 0.763 -0.25  0.18
## 18                                      did.not.vote.dummy -0.09 0.05   163.48  -1.82 0.071 -0.19  0.01
## 19                                         dont.know.dummy -0.12 0.06   294.72  -1.98 0.048 -0.23 -0.00
## 20                                      invalid.vote.dummy -0.04 0.35  1200.25  -0.12 0.907 -0.73  0.65
## 21                                         no.answer.dummy  0.02 0.37  1493.88   0.06 0.956 -0.71  0.76
## 22                                  not.eligible.age.dummy  0.19 0.06   306.55   3.23 0.001  0.07  0.31
## 23                          not.eligible.citizenship.dummy  0.31 0.07   287.33   4.78 0.000  0.18  0.44
## 24                                not.eligible.other.dummy  0.17 0.09   897.11   1.83 0.068 -0.01  0.34
## 25                                    anti.imm.party.dummy -0.54 0.05   230.15 -11.29 0.000 -0.63 -0.45
## 26                                     pro.env.party.dummy  0.36 0.05   274.55   6.88 0.000  0.25  0.46
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



\newpage



\newpage

# Hypothesis 3: Those who voted for anti-immigration parties will report lower pro-environment attitudes than those who voted for pro-environment parties.

### Model 0: random intercepts for environment attitudes



```r
H3.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H3.mod0<-getFE(H3.mod0))
```

```
##           Eff  Est   SE    df    t     p    LL   UL
## 1 (Intercept) 0.05 0.07 20.08 0.78 0.443 -0.09 0.20
```

```r
(VC.H3.mod0<-getVC(H3.mod0))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2851622 0.08131749
## 2        cntry (Intercept) <NA> 0.3004455 0.09026750
## 3     Residual        <NA> <NA> 1.1844251 1.40286286
```

```r
#ICC

##voting group

VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## [1] 0.05164826
```

```r
##country

VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]/
  sum(VC.H3.mod0[,"est_SD2"])
```

```
## [1] 0.0573328
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H3.mod0    4 114082 114116 -57037   114074                         
## H3.mod1   20 113129 113299 -56545   113089 984.48 16  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H3.mod1<-getFE(H3.mod1))
```

```
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.04 0.12   188.50  -0.32 0.751 -0.27  0.20
## 2                                                      age -0.00 0.00 29829.87  -8.11 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35659.58   1.14 0.253 -0.01  0.04
## 4                                                     educ  0.03 0.00 35342.84  13.39 0.000  0.02  0.03
## 5                                                    resid -0.14 0.01 35731.67 -10.19 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35578.76   0.42 0.672 -0.15  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35588.19  -0.60 0.548 -0.26  0.14
## 8                              occupElementary occupations -0.01 0.10 35591.13  -0.11 0.915 -0.21  0.19
## 9                                            occupManagers  0.13 0.10 35579.35   1.33 0.184 -0.06  0.33
## 10                            occupOther: Not in paid work  0.16 0.10 35719.64   1.55 0.121 -0.04  0.36
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35588.53  -0.41 0.684 -0.24  0.16
## 12                                      occupProfessionals  0.21 0.10 35581.91   2.09 0.037  0.01  0.40
## 13                                            occupRetired  0.06 0.11 35582.94   0.52 0.604 -0.16  0.28
## 14                          occupService and sales workers  0.02 0.10 35580.41   0.25 0.806 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35595.50   0.43 0.664 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35577.81   0.59 0.557 -0.14  0.25
## 17                                         occupUnemployed -0.03 0.12 35586.90  -0.25 0.803 -0.27  0.21
```

```r
(VC.H3.mod1<-getVC(H3.mod1))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2357854 0.05559476
## 2        cntry (Intercept) <NA> 0.2942470 0.08658130
## 3     Residual        <NA> <NA> 1.1693790 1.36744715
```

```r
#variance explained

##lvl 1: individuals

(VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="Residual","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="Residual","est_SD2"]
```

```
## [1] 0.02524532
```

```r
##lvl 2: voting group

(VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="voting.group","est_SD2"]
```

```
## [1] 0.3163247
```

```r
##lvl 3: country

(VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]-
     VC.H3.mod1[VC.H3.mod1$grp=="cntry","est_SD2"])/
  VC.H3.mod0[VC.H3.mod0$grp=="cntry","est_SD2"]
```

```
## [1] 0.04083635
```

```r
##total

(sum(VC.H3.mod0$est_SD2)-sum(VC.H3.mod1$est_SD2))/
  sum(VC.H3.mod0$est_SD2)
```

```
## [1] 0.04117294
```

```r
#individual contributions of covariates
anova(H3.mod1)
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


### Model 2: Categorical predictor at level-2


```r
H3.mod2<-lmer(environ.gmc~(1|voting.group)+(1|cntry)+
                age+gender+educ+resid+occup+
                all.parties.lvl2
                ,data=dat,REML=F)


(FE.H3.mod2<-getFE(H3.mod2))
```

```
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.27 0.12   280.06  -2.20 0.029 -0.50 -0.03
## 2                                                      age -0.00 0.00 35583.97  -7.51 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35672.30   1.19 0.235 -0.01  0.04
## 4                                                     educ  0.03 0.00 35720.39  13.27 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35734.98 -10.05 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35641.21   0.40 0.689 -0.16  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39  -0.60 0.548 -0.26  0.14
## 8                              occupElementary occupations -0.01 0.10 35645.09  -0.13 0.900 -0.21  0.18
## 9                                            occupManagers  0.13 0.10 35642.84   1.34 0.181 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06   1.32 0.188 -0.07  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76  -0.41 0.682 -0.24  0.16
## 12                                      occupProfessionals  0.20 0.10 35641.83   2.06 0.039  0.01  0.40
## 13                                            occupRetired  0.05 0.11 35635.37   0.49 0.624 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35640.53   0.24 0.812 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17   0.45 0.656 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21   0.58 0.564 -0.14  0.25
## 17                                         occupUnemployed -0.04 0.12 35636.70  -0.30 0.765 -0.27  0.20
## 18                            all.parties.lvl2Did not vote  0.10 0.05   151.88   1.82 0.070 -0.01  0.20
## 19                              all.parties.lvl2Don't know  0.08 0.06   275.08   1.18 0.240 -0.05  0.20
## 20                            all.parties.lvl2Invalid vote -0.62 0.36  2074.65  -1.69 0.092 -1.33  0.10
## 21                                  all.parties.lvl2NE age  0.39 0.06   286.34   6.13 0.000  0.27  0.52
## 22                              all.parties.lvl2NE citizen  0.24 0.07   254.61   3.48 0.001  0.10  0.37
## 23                                all.parties.lvl2NE other  0.30 0.10   835.68   3.08 0.002  0.11  0.49
## 24                               all.parties.lvl2No answer  0.20 0.39  2698.12   0.52 0.603 -0.56  0.97
## 25                             all.parties.lvl2Other party  0.22 0.04   204.52   5.09 0.000  0.13  0.30
## 26                   all.parties.lvl2Pro-environment party  0.70 0.06   241.69  12.15 0.000  0.59  0.82
```

```r
(VC.H3.mod2<-getVC(H3.mod2))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1546593 0.02391949
## 2        cntry (Intercept) <NA> 0.2729290 0.07449022
## 3     Residual        <NA> <NA> 1.1692947 1.36725014
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H3.mod1   20 113129 113299 -56545   113089                         
## H3.mod2   29 112998 113244 -56470   112940 149.01  9  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(H3.mod2)
```

```
## Type III Analysis of Variance Table with Satterthwaite's method
##                   Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
## age               77.159  77.159     1 35584  56.4340 5.948e-14 ***
## gender             1.929   1.929     1 35672   1.4112    0.2349    
## educ             240.908 240.908     1 35720 176.1991 < 2.2e-16 ***
## resid            137.986 137.986     1 35735 100.9224 < 2.2e-16 ***
## occup            194.516  16.210    12 35592  11.8557 < 2.2e-16 ***
## all.parties.lvl2 263.281  29.253     9   336  21.3958 < 2.2e-16 ***
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
## [1] 0.5697528
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
  cbind.data.frame(group=H3.mod2.mmeans.tab[,"group"],
        emmean=round_tidy(H3.mod2.mmeans.tab[,"emmean"],2),
        CI=paste0("[",
                  round_tidy(H3.mod2.mmeans.tab[,"asymp.LCL"],2),
                  ", ",
                  round_tidy(H3.mod2.mmeans.tab[,"asymp.UCL"],2),
                  "]"),
        p=round_tidy(H3.mod2.mmeans.tab[,"p"],3),
        p.adj=round_tidy(H3.mod2.mmeans.tab[,"adj.p"],3),
        p_less_001=ifelse(H3.mod2.mmeans.tab[,"p"]<.001,"yes","no"),
        p.adj_less_001=ifelse(H3.mod2.mmeans.tab[,"adj.p"]<.001,"yes","no"))
H3.mod2.mmeans.tab
```

```
##                     group emmean             CI     p p.adj p_less_001 p.adj_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no             no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no             no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no             no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no             no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no             no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no             no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no             no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no             no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no             no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes            yes
```

```r
write.csv2(H3.mod2.mmeans.tab,"H3.mod2.mmeans.tab.csv")

#contrast between anti-immigration and pro-environment
(H3.contrast<-data.frame(pairs(H3.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df  z.ratio      p.value
## 1 Pro-environment party - Anti-immigration party 0.7023866 0.05780773 Inf 12.15039 5.709007e-34
```

```r
#contrast for all groups against mean of other groups
contrast(H3.mod2.mmeans, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.1782 0.0720 Inf -2.476  0.1064 
##  Did not vote effect            -0.0695 0.0719 Inf -0.967  1.0000 
##  Don't know effect              -0.0949 0.0796 Inf -1.192  1.0000 
##  Invalid vote effect            -0.8625 0.3660 Inf -2.357  0.1291 
##  NE age effect                   0.2565 0.0794 Inf  3.229  0.0112 
##  NE citizen effect               0.0865 0.0831 Inf  1.041  1.0000 
##  NE other effect                 0.1519 0.1074 Inf  1.414  0.9438 
##  No answer effect                0.0472 0.3910 Inf  0.121  1.0000 
##  Other party effect              0.0608 0.0637 Inf  0.954  1.0000 
##  Pro-environment party effect    0.6022 0.0750 Inf  8.026  <.0001 
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
##                                         contrast  estimate         SE  df   z.ratio      p.value
## 1           Other party - Anti-immigration party 0.2151125 0.04227427 Inf  5.088498 3.609110e-07
## 2 Pro-environment party - Anti-immigration party 0.7023866 0.05780773 Inf 12.150394 1.712702e-33
## 3            Pro-environment party - Other party 0.4872741 0.04683923 Inf 10.403119 4.798968e-25
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
## [1] 1.189101
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
## [1] 1.178441
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
## [1] 1.183783
```

```r
(H3.effect.size<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd)
```

```
## Warning in Ops.factor(H3.mod2.mmeans.tab[10, 2], H3.mod2.mmeans.tab[1, 2]): '-' not meaningful for factors
```

```
## [1] NA
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
## [1] 1.193834
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
## [1] 1.186822
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean             CI     p p.adj p_less_001 p.adj_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no             no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no             no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no             no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no             no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no             no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no             no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no             no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no             no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no             no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes            yes
```

```r
(H3.effect.size.env.other<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[9,2])/H3.pooled.sd.other.env)
```

```
## Warning in Ops.factor(H3.mod2.mmeans.tab[10, 2], H3.mod2.mmeans.tab[9, 2]): '-' not meaningful for factors
```

```
## [1] NA
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
## [1] 1.191673
```

```r
H3.mod2.mmeans.tab
```

```
##                     group emmean             CI     p p.adj p_less_001 p.adj_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no             no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no             no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no             no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no             no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no             no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no             no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no             no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no             no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no             no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes            yes
```

```r
(H3.effect.size.imm.other<-(H3.mod2.mmeans.tab[9,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd.other.imm)
```

```
## Warning in Ops.factor(H3.mod2.mmeans.tab[9, 2], H3.mod2.mmeans.tab[1, 2]): '-' not meaningful for factors
```

```
## [1] NA
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.17 0.12   278.89  -1.39 0.167 -0.41  0.07
## 2                                                      age -0.00 0.00 35583.97  -7.51 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35672.29   1.19 0.235 -0.01  0.04
## 4                                                     educ  0.03 0.00 35720.39  13.27 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35734.98 -10.05 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35641.21   0.40 0.689 -0.16  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39  -0.60 0.548 -0.26  0.14
## 8                              occupElementary occupations -0.01 0.10 35645.08  -0.13 0.900 -0.21  0.18
## 9                                            occupManagers  0.13 0.10 35642.84   1.34 0.181 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06   1.32 0.188 -0.07  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76  -0.41 0.682 -0.24  0.16
## 12                                      occupProfessionals  0.20 0.10 35641.83   2.06 0.039  0.01  0.40
## 13                                            occupRetired  0.05 0.11 35635.36   0.49 0.624 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35640.53   0.24 0.812 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17   0.45 0.656 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21   0.58 0.564 -0.14  0.25
## 17                                         occupUnemployed -0.04 0.12 35636.70  -0.30 0.765 -0.27  0.20
## 18                                       other.party.dummy  0.12 0.04   132.01   2.79 0.006  0.03  0.20
## 19                                         dont.know.dummy -0.02 0.06   216.96  -0.36 0.720 -0.15  0.10
## 20                                      invalid.vote.dummy -0.71 0.36  2042.66  -1.96 0.051 -1.43  0.00
## 21                                         no.answer.dummy  0.10 0.39  2638.69   0.27 0.788 -0.66  0.87
## 22                                  not.eligible.age.dummy  0.29 0.06   214.47   4.65 0.000  0.17  0.42
## 23                          not.eligible.citizenship.dummy  0.14 0.07   206.74   2.06 0.041  0.01  0.27
## 24                                not.eligible.other.dummy  0.20 0.10   697.82   2.06 0.039  0.01  0.39
## 25                                    anti.imm.party.dummy -0.10 0.05   151.88  -1.82 0.070 -0.20  0.01
## 26                                     pro.env.party.dummy  0.60 0.06   185.57  10.47 0.000  0.49  0.72
```

```r
(VC.H3.mod3<-getVC(H3.mod3))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1546575 0.02391893
## 2        cntry (Intercept) <NA> 0.2729257 0.07448842
## 3     Residual        <NA> <NA> 1.1692948 1.36725029
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
##         npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)    
## H3.mod2   29 112998 113244 -56470   112940                        
## H3.mod3   29 112998 113244 -56470   112940     0  0  < 2.2e-16 ***
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
## [1] TRUE
```

```r
(FE.H3.mod4<-getFE(H3.mod4))
```

```
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.17 0.12   286.87  -1.41 0.161 -0.41  0.07
## 2                                                      age -0.00 0.00 35536.29  -7.50 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35674.84   1.18 0.239 -0.01  0.04
## 4                                                     educ  0.03 0.00 35701.59  13.27 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35722.10 -10.05 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35644.30   0.41 0.682 -0.16  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35652.54  -0.59 0.553 -0.25  0.14
## 8                              occupElementary occupations -0.01 0.10 35648.15  -0.12 0.907 -0.21  0.19
## 9                                            occupManagers  0.14 0.10 35646.01   1.35 0.178 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35698.11   1.31 0.190 -0.07  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35653.05  -0.40 0.686 -0.24  0.16
## 12                                      occupProfessionals  0.20 0.10 35644.89   2.06 0.039  0.01  0.40
## 13                                            occupRetired  0.06 0.11 35638.46   0.50 0.618 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35643.56   0.24 0.809 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35663.89   0.45 0.649 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35643.20   0.58 0.559 -0.14  0.25
## 17                                         occupUnemployed -0.04 0.12 35640.61  -0.29 0.770 -0.27  0.20
## 18                                       other.party.dummy  0.12 0.04   116.00   2.94 0.004  0.04  0.20
## 19                                         dont.know.dummy -0.02 0.06   198.20  -0.38 0.705 -0.14  0.10
## 20                                      invalid.vote.dummy -0.71 0.36  2073.56  -1.97 0.048 -1.42 -0.00
## 21                                         no.answer.dummy  0.10 0.39  2692.24   0.25 0.799 -0.66  0.86
## 22                                  not.eligible.age.dummy  0.29 0.06   196.18   4.83 0.000  0.17  0.41
## 23                          not.eligible.citizenship.dummy  0.14 0.07   186.70   2.16 0.032  0.01  0.27
## 24                                not.eligible.other.dummy  0.20 0.09   666.71   2.09 0.037  0.01  0.38
## 25                                    anti.imm.party.dummy -0.10 0.05   134.56  -1.86 0.065 -0.20  0.01
## 26                                     pro.env.party.dummy  0.60 0.07    29.76   8.53 0.000  0.46  0.74
```

```r
(VC.H3.mod4<-getVC(H3.mod4))
```

```
##            grp                 var1 var2    est_SD    est_SD2
## 1 voting.group          (Intercept) <NA> 0.1457504 0.02124319
## 2        cntry          (Intercept) <NA> 0.2681739 0.07191726
## 3      cntry.1 anti.imm.party.dummy <NA> 0.0000000 0.00000000
## 4      cntry.2  pro.env.party.dummy <NA> 0.1677282 0.02813274
## 5     Residual                 <NA> <NA> 1.1693366 1.36734801
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
##         npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H3.mod3   29 112998 113244 -56470   112940                       
## H3.mod4   31 112998 113261 -56468   112936 4.6209  2    0.09922 .
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.02 0.12   194.57  -0.21 0.834 -0.26  0.21
## 2                                                      age -0.00 0.00 35715.03  -7.62 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35632.32   1.30 0.193 -0.01  0.04
## 4                                                     educ  0.03 0.00 35736.79  13.50 0.000  0.02  0.03
## 5                                                    resid -0.14 0.01 35717.92 -10.15 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35587.55   0.43 0.665 -0.15  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35595.29  -0.58 0.562 -0.25  0.14
## 8                              occupElementary occupations -0.01 0.10 35591.57  -0.10 0.924 -0.21  0.19
## 9                                            occupManagers  0.13 0.10 35588.02   1.33 0.183 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35642.79   1.35 0.176 -0.06  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35594.52  -0.38 0.701 -0.24  0.16
## 12                                      occupProfessionals  0.21 0.10 35590.81   2.10 0.036  0.01  0.40
## 13                                            occupRetired  0.06 0.11 35584.13   0.51 0.611 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35587.43   0.25 0.802 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35604.00   0.46 0.647 -0.16  0.26
## 16            occupTechnicians and associate professionals  0.06 0.10 35586.88   0.60 0.549 -0.14  0.25
## 17                                         occupUnemployed -0.03 0.12 35584.96  -0.26 0.792 -0.27  0.21
## 18                                      did.not.vote.dummy -0.15 0.06   153.99  -2.63 0.009 -0.26 -0.04
## 19                                         dont.know.dummy -0.17 0.07   280.18  -2.60 0.010 -0.30 -0.04
## 20                                      invalid.vote.dummy -0.83 0.40  1108.73  -2.07 0.038 -1.61 -0.04
## 21                                         no.answer.dummy -0.04 0.42  1381.75  -0.09 0.930 -0.87  0.79
## 22                                  not.eligible.age.dummy  0.15 0.07   290.03   2.29 0.023  0.02  0.28
## 23                          not.eligible.citizenship.dummy  0.00 0.07   274.56   0.02 0.986 -0.14  0.15
## 24                                not.eligible.other.dummy  0.08 0.10   863.54   0.74 0.458 -0.12  0.27
```

```r
(VC.H3.mod5<-getVC(H3.mod5))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2224381 0.04947869
## 2        cntry (Intercept) <NA> 0.2928747 0.08577560
## 3     Residual        <NA> <NA> 1.1693366 1.36734811
```

```r
#see how much variance was explained at level-2

##lvl 2: voting group

(H3.total.eff<-(VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"]-
     VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])/
  VC.H3.mod1[VC.H3.mod1$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1100117
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept)  0.02 0.12   201.45   0.17 0.867 -0.21  0.25
## 2                                                      age -0.00 0.00 35693.16  -7.73 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35634.76   1.25 0.210 -0.01  0.04
## 4                                                     educ  0.03 0.00 35737.10  13.43 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35728.67 -10.12 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35598.90   0.43 0.667 -0.15  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35607.53  -0.56 0.575 -0.25  0.14
## 8                              occupElementary occupations -0.01 0.10 35603.29  -0.09 0.927 -0.21  0.19
## 9                                            occupManagers  0.13 0.10 35599.54   1.34 0.182 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35657.41   1.36 0.175 -0.06  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35606.85  -0.37 0.708 -0.24  0.16
## 12                                      occupProfessionals  0.21 0.10 35601.95   2.10 0.036  0.01  0.40
## 13                                            occupRetired  0.06 0.11 35594.59   0.52 0.606 -0.16  0.28
## 14                          occupService and sales workers  0.03 0.10 35598.77   0.26 0.795 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35616.48   0.46 0.646 -0.16  0.26
## 16            occupTechnicians and associate professionals  0.06 0.10 35598.29   0.60 0.547 -0.13  0.25
## 17                                         occupUnemployed -0.03 0.12 35594.88  -0.27 0.791 -0.27  0.21
## 18                                      did.not.vote.dummy -0.19 0.05   144.47  -3.70 0.000 -0.29 -0.09
## 19                                         dont.know.dummy -0.22 0.06   279.00  -3.46 0.001 -0.34 -0.09
## 20                                      invalid.vote.dummy -0.86 0.39  1274.76  -2.23 0.026 -1.62 -0.10
## 21                                         no.answer.dummy -0.08 0.41  1605.89  -0.20 0.839 -0.89  0.72
## 22                                  not.eligible.age.dummy  0.10 0.06   292.43   1.63 0.103 -0.02  0.23
## 23                          not.eligible.citizenship.dummy -0.05 0.07   265.67  -0.67 0.504 -0.18  0.09
## 24                                not.eligible.other.dummy  0.02 0.10   895.80   0.22 0.830 -0.17  0.21
## 25                                    anti.imm.party.dummy -0.29 0.05   209.84  -5.78 0.000 -0.39 -0.19
```

```r
(VC.H3.mod6<-getVC(H3.mod6))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.2021118 0.04084919
## 2        cntry (Intercept) <NA> 0.2907726 0.08454873
## 3     Residual        <NA> <NA> 1.1693914 1.36747627
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod6[VC.H3.mod6$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.1744084
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.09 0.12   238.38  -0.78 0.437 -0.32  0.14
## 2                                                      age -0.00 0.00 35625.90  -7.41 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35669.67   1.24 0.214 -0.01  0.04
## 4                                                     educ  0.03 0.00 35725.28  13.35 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35738.79 -10.08 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35631.32   0.40 0.686 -0.16  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35639.11  -0.62 0.537 -0.26  0.13
## 8                              occupElementary occupations -0.01 0.10 35635.03  -0.13 0.899 -0.21  0.18
## 9                                            occupManagers  0.13 0.10 35632.76   1.34 0.182 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35686.91   1.32 0.188 -0.07  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35639.22  -0.42 0.677 -0.24  0.16
## 12                                      occupProfessionals  0.20 0.10 35632.21   2.06 0.039  0.01  0.40
## 13                                            occupRetired  0.05 0.11 35625.94   0.48 0.629 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35630.70   0.23 0.818 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35649.61   0.45 0.654 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35630.33   0.58 0.565 -0.14  0.25
## 17                                         occupUnemployed -0.04 0.12 35626.79  -0.30 0.768 -0.27  0.20
## 18                                      did.not.vote.dummy -0.08 0.04   139.08  -1.74 0.084 -0.17  0.01
## 19                                         dont.know.dummy -0.10 0.06   314.98  -1.79 0.074 -0.21  0.01
## 20                                      invalid.vote.dummy -0.80 0.37  1875.75  -2.16 0.031 -1.53 -0.07
## 21                                         no.answer.dummy  0.03 0.40  2402.60   0.07 0.942 -0.75  0.80
## 22                                  not.eligible.age.dummy  0.22 0.06   336.32   3.87 0.000  0.11  0.33
## 23                          not.eligible.citizenship.dummy  0.07 0.06   282.81   1.05 0.296 -0.06  0.19
## 24                                not.eligible.other.dummy  0.13 0.09  1101.05   1.40 0.163 -0.05  0.31
## 25                                     pro.env.party.dummy  0.53 0.05   271.99  10.85 0.000  0.43  0.62
```

```r
(VC.H3.mod7<-getVC(H3.mod7))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1694007 0.02869660
## 2        cntry (Intercept) <NA> 0.2733412 0.07471543
## 3     Residual        <NA> <NA> 1.1692469 1.36713828
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod7[VC.H3.mod7$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.4200209
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
##                                                        Eff   Est   SE       df      t     p    LL    UL
## 1                                              (Intercept) -0.05 0.12   242.61  -0.44 0.663 -0.28  0.18
## 2                                                      age -0.00 0.00 35583.97  -7.51 0.000 -0.00 -0.00
## 3                                                   gender  0.02 0.01 35672.29   1.19 0.235 -0.01  0.04
## 4                                                     educ  0.03 0.00 35720.39  13.27 0.000  0.02  0.03
## 5                                                    resid -0.13 0.01 35734.98 -10.05 0.000 -0.16 -0.11
## 6                            occupClerical support workers  0.04 0.10 35641.21   0.40 0.689 -0.16  0.24
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39  -0.60 0.548 -0.26  0.14
## 8                              occupElementary occupations -0.01 0.10 35645.08  -0.13 0.900 -0.21  0.18
## 9                                            occupManagers  0.13 0.10 35642.84   1.34 0.181 -0.06  0.33
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06   1.32 0.188 -0.07  0.34
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76  -0.41 0.682 -0.24  0.16
## 12                                      occupProfessionals  0.20 0.10 35641.83   2.06 0.039  0.01  0.40
## 13                                            occupRetired  0.05 0.11 35635.36   0.49 0.624 -0.16  0.27
## 14                          occupService and sales workers  0.02 0.10 35640.53   0.24 0.812 -0.17  0.22
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17   0.45 0.656 -0.16  0.25
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21   0.58 0.564 -0.14  0.25
## 17                                         occupUnemployed -0.04 0.12 35636.70  -0.30 0.765 -0.27  0.20
## 18                                      did.not.vote.dummy -0.12 0.04   132.01  -2.79 0.006 -0.20 -0.03
## 19                                         dont.know.dummy -0.14 0.05   316.43  -2.59 0.010 -0.25 -0.03
## 20                                      invalid.vote.dummy -0.83 0.36  2173.94  -2.29 0.022 -1.54 -0.12
## 21                                         no.answer.dummy -0.01 0.39  2799.75  -0.03 0.975 -0.77  0.75
## 22                                  not.eligible.age.dummy  0.18 0.05   342.85   3.22 0.001  0.07  0.28
## 23                          not.eligible.citizenship.dummy  0.02 0.06   276.20   0.39 0.697 -0.09  0.14
## 24                                not.eligible.other.dummy  0.08 0.09  1145.84   0.90 0.367 -0.10  0.26
## 25                                    anti.imm.party.dummy -0.22 0.04   204.52  -5.09 0.000 -0.30 -0.13
## 26                                     pro.env.party.dummy  0.49 0.05   263.65  10.40 0.000  0.40  0.58
```

```r
(VC.H3.mod8<-getVC(H3.mod8))
```

```
##            grp        var1 var2    est_SD    est_SD2
## 1 voting.group (Intercept) <NA> 0.1546575 0.02391895
## 2        cntry (Intercept) <NA> 0.2729257 0.07448844
## 3     Residual        <NA> <NA> 1.1692948 1.36725028
```

```r
(H3.total.eff<-(VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"]-
     VC.H3.mod8[VC.H3.mod8$grp=="voting.group","est_SD2"])/
  VC.H3.mod5[VC.H3.mod5$grp=="voting.group","est_SD2"])
```

```
## [1] 0.5165807
```

```r
#variance not accounted for at level-2
1-H3.total.eff
```

```
## [1] 0.4834193
```


\newpage



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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    62.06 -3.22 0.002 -0.75 -0.18
## 2                                                      age  0.00 0.00 35705.43  4.89 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35593.59  4.67 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35719.00  7.54 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35683.44 -5.14 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35568.33 -0.51 0.611 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35579.01 -0.78 0.434 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35575.90  0.09 0.931 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35568.40 -0.24 0.810 -0.19  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35624.19  1.32 0.186 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35576.73 -0.76 0.447 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35571.35  0.76 0.450 -0.11  0.24
## 13                                            occupRetired -0.04 0.10 35560.28 -0.41 0.682 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35572.75 -0.46 0.647 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35587.51 -0.66 0.506 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35568.43 -0.40 0.690 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35569.10 -0.30 0.763 -0.24  0.18
## 18                                            environ.lvl1  0.12 0.01    19.43 11.50 0.000  0.10  0.15
## 19                            all.parties.lvl2Did not vote  0.45 0.06   182.54  7.09 0.000  0.32  0.57
## 20                              all.parties.lvl2Don't know  0.42 0.07   269.12  5.99 0.000  0.28  0.56
## 21                            all.parties.lvl2Invalid vote  0.49 0.35  1099.45  1.38 0.169 -0.21  1.18
## 22                                  all.parties.lvl2NE age  0.74 0.07   272.45 10.51 0.000  0.60  0.88
## 23                              all.parties.lvl2NE citizen  0.86 0.08   268.81 11.34 0.000  0.71  1.01
## 24                                all.parties.lvl2NE other  0.72 0.10   659.04  7.33 0.000  0.52  0.91
## 25                               all.parties.lvl2No answer  0.56 0.37  1373.88  1.50 0.133 -0.17  1.30
## 26                             all.parties.lvl2Other party  0.54 0.05   230.96 11.28 0.000  0.45  0.64
## 27                   all.parties.lvl2Pro-environment party  0.91 0.06   256.35 14.05 0.000  0.78  1.04
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
## Warning: Some predictor variables are on very different scales: consider rescaling

## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    62.05 -3.22 0.002 -0.75 -0.18
## 2                                                      age  0.00 0.00 35703.98  4.89 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35594.36  4.68 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35717.58  7.53 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35680.08 -5.15 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.05 0.09 35565.71 -0.52 0.603 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35577.12 -0.79 0.428 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35573.30  0.07 0.943 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35565.33 -0.25 0.800 -0.20  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35622.51  1.32 0.188 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35574.13 -0.77 0.442 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35567.54  0.75 0.456 -0.11  0.24
## 13                                            occupRetired -0.04 0.10 35553.77 -0.41 0.680 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35570.97 -0.47 0.640 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35584.95 -0.67 0.503 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.04 0.09 35565.98 -0.41 0.684 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35567.13 -0.30 0.765 -0.24  0.18
## 18                                            environ.lvl1  0.11 0.02   106.04  5.56 0.000  0.07  0.15
## 19                            all.parties.lvl2Did not vote  0.45 0.06   182.70  7.07 0.000  0.32  0.57
## 20                              all.parties.lvl2Don't know  0.42 0.07   268.96  6.00 0.000  0.29  0.56
## 21                            all.parties.lvl2Invalid vote  0.49 0.35  1099.20  1.38 0.168 -0.21  1.19
## 22                                  all.parties.lvl2NE age  0.74 0.07   272.33 10.54 0.000  0.61  0.88
## 23                              all.parties.lvl2NE citizen  0.86 0.08   268.70 11.31 0.000  0.71  1.01
## 24                                all.parties.lvl2NE other  0.72 0.10   658.23  7.34 0.000  0.53  0.91
## 25                               all.parties.lvl2No answer  0.57 0.37  1373.19  1.51 0.131 -0.17  1.30
## 26                             all.parties.lvl2Other party  0.55 0.05   230.94 11.31 0.000  0.45  0.64
## 27                   all.parties.lvl2Pro-environment party  0.91 0.06   256.28 14.07 0.000  0.78  1.04
## 28               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    85.05  0.09 0.926 -0.04  0.05
## 29                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   333.74  0.82 0.413 -0.04  0.10
## 30               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25263.41  0.15 0.878 -0.70  0.82
## 31                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   232.54  1.28 0.203 -0.02  0.10
## 32                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   231.02 -1.12 0.265 -0.10  0.03
## 33                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1198.34  0.59 0.557 -0.08  0.14
## 34                  environ.lvl1:all.parties.lvl2No answer  0.23 0.23  7338.62  0.99 0.320 -0.22  0.68
## 35                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   132.60  0.98 0.329 -0.02  0.06
## 36      environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   191.66  0.92 0.358 -0.03  0.08
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
##                     group beta            CI     p p.adj p_less_001 p.adj_less_001
## 1  Anti-immigration party 0.11  [0.07, 0.15] 0.000 0.000        yes            yes
## 2            Did not vote 0.11  [0.08, 0.14] 0.000 0.000        yes            yes
## 3              Don't know 0.14  [0.08, 0.20] 0.000 0.000        yes            yes
## 4            Invalid vote 0.17 [-0.59, 0.93] 0.664 0.664         no             no
## 5                  NE age 0.15  [0.10, 0.20] 0.000 0.000        yes            yes
## 6              NE citizen 0.07  [0.01, 0.13] 0.016 0.048         no             no
## 7                NE other 0.14  [0.04, 0.25] 0.008 0.034         no             no
## 8               No answer 0.34 [-0.11, 0.79] 0.141 0.282         no             no
## 9             Other party 0.13  [0.10, 0.15] 0.000 0.000        yes            yes
## 10  Pro-environment party 0.14  [0.09, 0.18] 0.000 0.000        yes            yes
```

```r
write.csv2(H4.mod2.trends.tab,"H4.mod2.trends.tab.csv")


#contrast between anti-immigration and pro-environment
(H4.contrast<-data.frame(pairs(H4.mod2.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df   z.ratio   p.value
## 1 Pro-environment party - Anti-immigration party 0.02618252 0.02840579 Inf 0.9217319 0.3566685
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
##                                         contrast    estimate         SE  df   z.ratio   p.value
## 1           Other party - Anti-immigration party 0.018994304 0.01937022 Inf 0.9805933 0.3267934
## 2 Pro-environment party - Anti-immigration party 0.026182522 0.02840579 Inf 0.9217319 0.3566685
## 3            Pro-environment party - Other party 0.007188218 0.02369191 Inf 0.3034039 0.7615821
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
##                     group emmean             CI     p emmean             CI     p beta            CI     p
## 1  Anti-immigration party  -0.48 [-0.71, -0.25] 0.000  -0.22 [-0.36, -0.08] 0.002 0.11  [0.07, 0.15] 0.000
## 2            Did not vote  -0.03  [-0.26, 0.20] 0.781  -0.12  [-0.27, 0.02] 0.086 0.11  [0.08, 0.14] 0.000
## 3              Don't know  -0.06  [-0.30, 0.18] 0.630  -0.15  [-0.30, 0.01] 0.066 0.14  [0.08, 0.20] 0.000
## 4            Invalid vote   0.02  [-0.71, 0.74] 0.963  -0.84 [-1.56, -0.12] 0.023 0.17 [-0.59, 0.93] 0.664
## 5                  NE age   0.25   [0.01, 0.49] 0.040   0.17   [0.01, 0.33] 0.034 0.15  [0.10, 0.20] 0.000
## 6              NE citizen   0.37   [0.13, 0.61] 0.003   0.02  [-0.15, 0.18] 0.850 0.07  [0.01, 0.13] 0.016
## 7                NE other   0.22  [-0.05, 0.50] 0.109   0.07  [-0.14, 0.29] 0.489 0.14  [0.04, 0.25] 0.008
## 8               No answer   0.08  [-0.68, 0.84] 0.839  -0.02  [-0.79, 0.75] 0.961 0.34 [-0.11, 0.79] 0.141
## 9             Other party   0.06  [-0.16, 0.27] 0.596  -0.01  [-0.13, 0.12] 0.911 0.13  [0.10, 0.15] 0.000
## 10  Pro-environment party   0.41   [0.18, 0.65] 0.000   0.48   [0.33, 0.63] 0.000 0.14  [0.09, 0.18] 0.000
```

```r
write.csv2(tab2,"tab2.csv")
```


\newpage



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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.43  0.48 0.634 -0.22  0.35
## 2                                                      age  0.00 0.00 33758.85  2.39 0.017  0.00  0.00
## 3                                                   gender  0.07 0.01 35565.52  5.66 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35675.86  6.63 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35635.00 -4.98 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35509.19 -0.42 0.673 -0.21  0.14
## 7                    occupCraft and related trades workers -0.06 0.09 35519.51 -0.67 0.502 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35521.98  0.27 0.787 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35508.37 -0.20 0.842 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35676.02  1.56 0.118 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.21 -0.64 0.522 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35513.08  0.83 0.409 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35510.01 -0.28 0.776 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35515.16 -0.35 0.725 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35524.18 -0.53 0.598 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35509.06 -0.33 0.738 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35526.75 -0.16 0.876 -0.23  0.19
## 18                                            environ.lvl1  0.12 0.01    19.45 11.43 0.000  0.10  0.14
## 19                                         engagement.lvl1  0.06 0.01 35609.62  7.86 0.000  0.05  0.08
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.44  0.48 0.637 -0.22  0.35
## 2                                                      age  0.00 0.00 33750.64  2.38 0.017  0.00  0.00
## 3                                                   gender  0.07 0.01 35565.68  5.67 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35675.75  6.64 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35635.19 -4.97 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35509.24 -0.42 0.674 -0.21  0.14
## 7                    occupCraft and related trades workers -0.06 0.09 35519.58 -0.67 0.504 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35522.03  0.27 0.786 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35508.41 -0.20 0.843 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35676.21  1.57 0.117 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.29 -0.64 0.525 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35513.13  0.83 0.408 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35510.04 -0.28 0.777 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35515.24 -0.35 0.728 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35524.27 -0.52 0.600 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35509.12 -0.33 0.741 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35526.77 -0.16 0.875 -0.23  0.19
## 18                                            environ.lvl1  0.12 0.01    19.46 11.42 0.000  0.10  0.14
## 19                                         engagement.lvl1  0.06 0.01 35610.30  7.87 0.000  0.05  0.08
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01 35543.53  1.03 0.301 -0.01  0.02
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.06 0.14    49.36  0.46 0.651 -0.22  0.35
## 2                                                      age  0.00 0.00 33343.34  2.59 0.010  0.00  0.00
## 3                                                   gender  0.07 0.01 35558.07  5.68 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35664.79  6.57 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35611.45 -4.94 0.000 -0.08 -0.03
## 6                            occupClerical support workers -0.03 0.09 35484.87 -0.39 0.699 -0.21  0.14
## 7                    occupCraft and related trades workers -0.06 0.09 35493.52 -0.63 0.527 -0.23  0.12
## 8                              occupElementary occupations  0.03 0.09 35493.55  0.31 0.755 -0.15  0.20
## 9                                            occupManagers -0.01 0.09 35483.48 -0.15 0.877 -0.19  0.16
## 10                            occupOther: Not in paid work  0.15 0.09 35654.04  1.62 0.106 -0.03  0.32
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35496.23 -0.59 0.555 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35488.49  0.84 0.402 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35480.11 -0.27 0.785 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35487.34 -0.31 0.758 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.04 0.09 35503.28 -0.48 0.634 -0.23  0.14
## 16            occupTechnicians and associate professionals -0.03 0.09 35485.04 -0.30 0.767 -0.20  0.15
## 17                                         occupUnemployed -0.02 0.11 35505.62 -0.16 0.870 -0.23  0.19
## 18                                            environ.lvl1  0.12 0.01    19.45 11.34 0.000  0.10  0.14
## 19                                         engagement.lvl1  0.06 0.01    21.71  4.38 0.000  0.03  0.09
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01 35519.96  1.05 0.295 -0.01  0.02
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.31  0.46 0.648 -0.22  0.35
## 2                                                      age  0.00 0.00 33320.82  2.60 0.009  0.00  0.00
## 3                                                   gender  0.07 0.01 35527.23  5.66 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35628.81  6.60 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35604.24 -4.95 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35454.99 -0.40 0.688 -0.21  0.14
## 7                    occupCraft and related trades workers -0.06 0.09 35464.39 -0.65 0.518 -0.23  0.12
## 8                              occupElementary occupations  0.03 0.09 35459.54  0.30 0.767 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35443.30 -0.17 0.863 -0.19  0.16
## 10                            occupOther: Not in paid work  0.15 0.09 35629.00  1.62 0.105 -0.03  0.32
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35458.67 -0.60 0.545 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35453.25  0.82 0.413 -0.10  0.24
## 13                                            occupRetired -0.02 0.10 35456.55 -0.25 0.801 -0.22  0.17
## 14                          occupService and sales workers -0.03 0.09 35452.88 -0.32 0.747 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.04 0.09 35474.75 -0.48 0.628 -0.23  0.14
## 16            occupTechnicians and associate professionals -0.03 0.09 35446.84 -0.31 0.755 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35488.30 -0.18 0.855 -0.23  0.19
## 18                                            environ.lvl1  0.12 0.01    19.65 11.47 0.000  0.10  0.14
## 19                                         engagement.lvl1  0.06 0.01    21.72  4.37 0.000  0.03  0.09
## 20                            environ.lvl1:engagement.lvl1  0.01 0.01    27.73  0.86 0.396 -0.01  0.03
```

```r
(VC.H5.mod4<-getVC(H5.mod4))
```

```
##             grp                         var1                         var2      est_SD       est_SD2
## 1  voting.group                  (Intercept)                         <NA>  0.30666977  0.0940463483
## 2  voting.group                 environ.lvl1                         <NA>  0.04316499  0.0018632163
## 3  voting.group              engagement.lvl1                         <NA>  0.06352955  0.0040360037
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>  0.04537679  0.0020590531
## 5  voting.group                  (Intercept)                 environ.lvl1  0.08618317  0.0011408409
## 6  voting.group                  (Intercept)              engagement.lvl1  0.29709444  0.0057881699
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1  0.30736238  0.0042771596
## 8  voting.group                 environ.lvl1              engagement.lvl1  0.04133690  0.0001133562
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1 -0.22779132 -0.0004461723
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.17800013 -0.0005131329
## 11        cntry                  (Intercept)                         <NA>  0.49968388  0.2496839810
## 12        cntry                 environ.lvl1                         <NA>  0.03883750  0.0015083511
## 13        cntry              engagement.lvl1                         <NA>  0.04584691  0.0021019391
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>  0.02252269  0.0005072716
## 15        cntry                  (Intercept)                 environ.lvl1 -0.05677543 -0.0011018106
## 16        cntry                  (Intercept)              engagement.lvl1  0.56550304  0.0129550876
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1 -0.29825728 -0.0033566547
## 18        cntry                 environ.lvl1              engagement.lvl1  0.05592332  0.0000995759
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1  0.90661876  0.0007930420
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.37051502 -0.0003825922
## 21     Residual                         <NA>                         <NA>  1.02577908  1.0522227306
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
##                                  voting.group.(Intercept)                     voting.group.environ.lvl1.(Intercept) 
##                                                   0.29896                                                   0.00363 
##                  voting.group.engagement.lvl1.(Intercept)     voting.group.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                   0.01840                                                   0.01360 
##                                 voting.group.environ.lvl1                 voting.group.engagement.lvl1.environ.lvl1 
##                                                   0.04192                                                   0.00098 
##    voting.group.environ.lvl1:engagement.lvl1.environ.lvl1                              voting.group.engagement.lvl1 
##                                                  -0.01129                                                   0.05913 
## voting.group.environ.lvl1:engagement.lvl1.engagement.lvl1                 voting.group.environ.lvl1:engagement.lvl1 
##                                                  -0.01229                                                   0.03864 
##                                         cntry.(Intercept)                            cntry.environ.lvl1.(Intercept) 
##                                                   0.48713                                                  -0.00215 
##                         cntry.engagement.lvl1.(Intercept)            cntry.environ.lvl1:engagement.lvl1.(Intercept) 
##                                                   0.02528                                                  -0.00655 
##                                        cntry.environ.lvl1                        cntry.engagement.lvl1.environ.lvl1 
##                                                   0.03780                                                   0.00394 
##           cntry.environ.lvl1:engagement.lvl1.environ.lvl1                                     cntry.engagement.lvl1 
##                                                   0.01957                                                   0.03665 
##        cntry.environ.lvl1:engagement.lvl1.engagement.lvl1                        cntry.environ.lvl1:engagement.lvl1 
##                                                  -0.00751                                                   0.00000
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.48 0.14    61.87 -3.30 0.002 -0.77 -0.19
## 2                                                      age  0.00 0.00 35562.27  2.89 0.004  0.00  0.00
## 3                                                   gender  0.07 0.01 35550.62  5.68 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35668.50  6.59 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35653.24 -4.88 0.000 -0.08 -0.03
## 6                            occupClerical support workers -0.04 0.09 35509.60 -0.42 0.673 -0.21  0.14
## 7                    occupCraft and related trades workers -0.06 0.09 35519.06 -0.64 0.521 -0.23  0.12
## 8                              occupElementary occupations  0.02 0.09 35510.18  0.26 0.797 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35496.94 -0.19 0.850 -0.19  0.16
## 10                            occupOther: Not in paid work  0.13 0.09 35568.62  1.43 0.153 -0.05  0.31
## 11        occupPlant and machine operators, and assemblers -0.05 0.09 35515.24 -0.61 0.543 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35509.52  0.80 0.425 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35499.27 -0.27 0.786 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35507.14 -0.33 0.743 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35535.86 -0.51 0.613 -0.23  0.14
## 16            occupTechnicians and associate professionals -0.03 0.09 35502.94 -0.32 0.753 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35531.23 -0.24 0.812 -0.24  0.18
## 18                                            environ.lvl1  0.11 0.02   107.98  5.49 0.000  0.07  0.15
## 19                                         engagement.lvl1  0.01 0.03   142.46  0.41 0.686 -0.04  0.07
## 20                            all.parties.lvl2Did not vote  0.44 0.06   182.83  6.98 0.000  0.32  0.57
## 21                              all.parties.lvl2Don't know  0.42 0.07   268.91  6.00 0.000  0.29  0.56
## 22                            all.parties.lvl2Invalid vote  0.48 0.35  1090.05  1.35 0.176 -0.22  1.17
## 23                                  all.parties.lvl2NE age  0.72 0.07   273.86 10.15 0.000  0.58  0.86
## 24                              all.parties.lvl2NE citizen  0.86 0.08   268.92 11.30 0.000  0.71  1.01
## 25                                all.parties.lvl2NE other  0.71 0.10   657.74  7.29 0.000  0.52  0.90
## 26                               all.parties.lvl2No answer  0.54 0.38  1390.63  1.44 0.150 -0.20  1.28
## 27                             all.parties.lvl2Other party  0.55 0.05   231.13 11.42 0.000  0.46  0.65
## 28                   all.parties.lvl2Pro-environment party  0.92 0.06   256.15 14.13 0.000  0.79  1.04
## 29                            environ.lvl1:engagement.lvl1  0.01 0.01    26.79  0.63 0.535 -0.01  0.02
## 30               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    86.20  0.04 0.967 -0.04  0.05
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   339.50  0.91 0.361 -0.04  0.10
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25331.83  0.15 0.882 -0.70  0.82
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.26  1.34 0.183 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   236.37 -1.06 0.290 -0.10  0.03
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1243.41  0.61 0.545 -0.08  0.14
## 36                  environ.lvl1:all.parties.lvl2No answer  0.24 0.23  7550.15  1.06 0.290 -0.21  0.69
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   134.21  0.92 0.361 -0.02  0.06
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.02 0.03   195.88  0.80 0.423 -0.03  0.08
## 39            engagement.lvl1:all.parties.lvl2Did not vote  0.03 0.03   105.78  0.76 0.450 -0.04  0.09
## 40              engagement.lvl1:all.parties.lvl2Don't know -0.00 0.05   513.61 -0.03 0.979 -0.10  0.10
## 41            engagement.lvl1:all.parties.lvl2Invalid vote  0.11 0.45 18420.88  0.25 0.799 -0.77  1.00
## 42                  engagement.lvl1:all.parties.lvl2NE age -0.01 0.05   354.94 -0.13 0.900 -0.10  0.08
## 43              engagement.lvl1:all.parties.lvl2NE citizen  0.01 0.05   258.49  0.21 0.837 -0.09  0.11
## 44                engagement.lvl1:all.parties.lvl2NE other -0.02 0.09  2570.22 -0.25 0.806 -0.20  0.16
## 45               engagement.lvl1:all.parties.lvl2No answer -0.68 0.72 31757.51 -0.94 0.346 -2.09  0.73
## 46             engagement.lvl1:all.parties.lvl2Other party  0.07 0.03   171.19  2.46 0.015  0.01  0.13
## 47   engagement.lvl1:all.parties.lvl2Pro-environment party  0.10 0.04   296.05  2.33 0.020  0.02  0.19
```

```r
(VC.H5.mod5<-getVC(H5.mod5))
```

```
##             grp                         var1                         var2      est_SD       est_SD2
## 1  voting.group                  (Intercept)                         <NA>  0.19753217  3.901896e-02
## 2  voting.group                 environ.lvl1                         <NA>  0.04094707  1.676662e-03
## 3  voting.group              engagement.lvl1                         <NA>  0.05449208  2.969387e-03
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>  0.04316977  1.863629e-03
## 5  voting.group                  (Intercept)                 environ.lvl1  0.06573845  5.317164e-04
## 6  voting.group                  (Intercept)              engagement.lvl1  0.15485894  1.666892e-03
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1 -0.07651028 -6.524352e-04
## 8  voting.group                 environ.lvl1              engagement.lvl1 -0.01103897 -2.463116e-05
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1 -0.23564127 -4.165373e-04
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.36840660 -8.666436e-04
## 11        cntry                  (Intercept)                         <NA>  0.48255877  2.328630e-01
## 12        cntry                 environ.lvl1                         <NA>  0.03970632  1.576592e-03
## 13        cntry              engagement.lvl1                         <NA>  0.04760744  2.266468e-03
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>  0.02373098  5.631596e-04
## 15        cntry                  (Intercept)                 environ.lvl1 -0.06057173 -1.160593e-03
## 16        cntry                  (Intercept)              engagement.lvl1  0.51021906  1.172146e-02
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1 -0.30364533 -3.477223e-03
## 18        cntry                 environ.lvl1              engagement.lvl1  0.10900597  2.060558e-04
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1  0.89680996  8.450372e-04
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.34142628 -3.857336e-04
## 21     Residual                         <NA>                         <NA>  1.02573145  1.052125e+00
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
##          all.parties.lvl2 engagement.lvl1.trend         SE  df   asymp.LCL  asymp.UCL
## 1  Anti-immigration party           0.011511487 0.02833794 Inf -0.04402985 0.06705282
## 2            Did not vote           0.036598204 0.02300658 Inf -0.00849387 0.08169028
## 3              Don't know           0.010135998 0.04581401 Inf -0.07965782 0.09992981
## 4            Invalid vote           0.126363928 0.45130400 Inf -0.75817565 1.01090351
## 5                  NE age           0.005759236 0.03915849 Inf -0.07098999 0.08250846
## 6              NE citizen           0.021501952 0.04214949 Inf -0.06110953 0.10411343
## 7                NE other          -0.010700005 0.08728812 Inf -0.18178157 0.16038156
## 8               No answer          -0.668262981 0.72111011 Inf -2.08161282 0.74508686
## 9             Other party           0.083001023 0.01669505 Inf  0.05027933 0.11572272
## 10  Pro-environment party           0.111997005 0.03582636 Inf  0.04177864 0.18221537
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
##                     group engagement.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party                  0.01 0.03 0.6846 1.0000     -0.04      0.07
## 2            Did not vote                  0.04 0.02 0.1117 0.8933     -0.01      0.08
## 3              Don't know                  0.01 0.05 0.8249 1.0000     -0.08      0.10
## 4            Invalid vote                  0.13 0.45 0.7795 1.0000     -0.76      1.01
## 5                  NE age                  0.01 0.04 0.8831 1.0000     -0.07      0.08
## 6              NE citizen                  0.02 0.04 0.6100 1.0000     -0.06      0.10
## 7                NE other                 -0.01 0.09 0.9024 1.0000     -0.18      0.16
## 8               No answer                 -0.67 0.72 0.3541 1.0000     -2.08      0.75
## 9             Other party                  0.08 0.02 0.0000 0.0000      0.05      0.12
## 10  Pro-environment party                  0.11 0.04 0.0018 0.0159      0.04      0.18
```

```r
write.csv2(H5.mod5.trends.tab,"H5.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df  z.ratio    p.value
## 1 Pro-environment party - Anti-immigration party 0.1004855 0.04306793 Inf 2.333187 0.01963836
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
##                                         contrast   estimate         SE  df   z.ratio    p.value
## 1           Other party - Anti-immigration party 0.07148954 0.02911709 Inf 2.4552429 0.01407894
## 2 Pro-environment party - Anti-immigration party 0.10048552 0.04306793 Inf 2.3331865 0.01963836
## 3            Pro-environment party - Other party 0.02899598 0.03638382 Inf 0.7969471 0.42548176
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                                   Eff   Est   SE       df     t     p    LL    UL
## 1                                                         (Intercept) -0.47 0.14    61.81 -3.26 0.002 -0.76 -0.18
## 2                                                                 age  0.00 0.00 35561.36  2.89 0.004  0.00  0.00
## 3                                                              gender  0.07 0.01 35553.22  5.72 0.000  0.04  0.09
## 4                                                                educ  0.01 0.00 35669.76  6.60 0.000  0.01  0.02
## 5                                                               resid -0.06 0.01 35654.28 -4.87 0.000 -0.08 -0.03
## 6                                       occupClerical support workers -0.04 0.09 35512.26 -0.45 0.655 -0.21  0.13
## 7                               occupCraft and related trades workers -0.06 0.09 35521.65 -0.67 0.503 -0.23  0.11
## 8                                         occupElementary occupations  0.02 0.09 35512.45  0.23 0.816 -0.15  0.19
## 9                                                       occupManagers -0.02 0.09 35499.42 -0.22 0.830 -0.19  0.15
## 10                                       occupOther: Not in paid work  0.13 0.09 35570.57  1.40 0.160 -0.05  0.31
## 11                   occupPlant and machine operators, and assemblers -0.06 0.09 35517.32 -0.64 0.521 -0.23  0.12
## 12                                                 occupProfessionals  0.07 0.09 35511.77  0.76 0.444 -0.10  0.24
## 13                                                       occupRetired -0.03 0.10 35501.08 -0.31 0.757 -0.22  0.16
## 14                                     occupService and sales workers -0.03 0.09 35509.95 -0.36 0.721 -0.20  0.14
## 15            occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35539.17 -0.53 0.598 -0.23  0.13
## 16                       occupTechnicians and associate professionals -0.03 0.09 35505.94 -0.34 0.731 -0.20  0.14
## 17                                                    occupUnemployed -0.03 0.11 35532.88 -0.24 0.808 -0.24  0.18
## 18                                                       environ.lvl1  0.11 0.02   107.96  5.60 0.000  0.07  0.15
## 19                                                    engagement.lvl1  0.02 0.03   142.26  0.53 0.597 -0.04  0.07
## 20                                       all.parties.lvl2Did not vote  0.44 0.06   182.69  6.95 0.000  0.32  0.57
## 21                                         all.parties.lvl2Don't know  0.42 0.07   268.25  5.95 0.000  0.28  0.56
## 22                                       all.parties.lvl2Invalid vote  0.47 0.35  1087.13  1.34 0.180 -0.22  1.17
## 23                                             all.parties.lvl2NE age  0.72 0.07   273.89 10.12 0.000  0.58  0.86
## 24                                         all.parties.lvl2NE citizen  0.85 0.08   268.13 11.12 0.000  0.70  1.00
## 25                                           all.parties.lvl2NE other  0.71 0.10   656.69  7.25 0.000  0.52  0.90
## 26                                          all.parties.lvl2No answer  0.54 0.38  1390.54  1.44 0.149 -0.20  1.28
## 27                                        all.parties.lvl2Other party  0.55 0.05   230.67 11.33 0.000  0.45  0.64
## 28                              all.parties.lvl2Pro-environment party  0.91 0.06   256.51 14.03 0.000  0.78  1.04
## 29                                       environ.lvl1:engagement.lvl1 -0.04 0.02   144.61 -1.70 0.091 -0.08  0.01
## 30                          environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    86.02 -0.01 0.991 -0.04  0.04
## 31                            environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   338.63  0.86 0.390 -0.04  0.10
## 32                          environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25427.50  0.16 0.870 -0.70  0.83
## 33                                environ.lvl1:all.parties.lvl2NE age  0.04 0.03   237.66  1.31 0.193 -0.02  0.10
## 34                            environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   233.93 -1.31 0.191 -0.11  0.02
## 35                              environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1246.82  0.67 0.505 -0.07  0.15
## 36                             environ.lvl1:all.parties.lvl2No answer  0.24 0.23  7561.32  1.05 0.292 -0.21  0.69
## 37                           environ.lvl1:all.parties.lvl2Other party  0.02 0.02   134.16  0.79 0.432 -0.02  0.05
## 38                 environ.lvl1:all.parties.lvl2Pro-environment party  0.02 0.03   195.40  0.74 0.463 -0.04  0.08
## 39                       engagement.lvl1:all.parties.lvl2Did not vote  0.02 0.03   104.77  0.70 0.487 -0.04  0.09
## 40                         engagement.lvl1:all.parties.lvl2Don't know -0.00 0.05   502.12 -0.09 0.926 -0.11  0.10
## 41                       engagement.lvl1:all.parties.lvl2Invalid vote  0.12 0.45 18649.67  0.28 0.783 -0.76  1.01
## 42                             engagement.lvl1:all.parties.lvl2NE age -0.01 0.05   348.95 -0.15 0.881 -0.10  0.08
## 43                         engagement.lvl1:all.parties.lvl2NE citizen -0.00 0.05   254.20 -0.10 0.919 -0.10  0.09
## 44                           engagement.lvl1:all.parties.lvl2NE other -0.01 0.09  2622.74 -0.12 0.901 -0.19  0.17
## 45                          engagement.lvl1:all.parties.lvl2No answer -0.68 0.72 31703.95 -0.94 0.345 -2.10  0.73
## 46                        engagement.lvl1:all.parties.lvl2Other party  0.07 0.03   168.71  2.29 0.023  0.01  0.12
## 47              engagement.lvl1:all.parties.lvl2Pro-environment party  0.10 0.04   291.53  2.25 0.025  0.01  0.18
## 48          environ.lvl1:engagement.lvl1:all.parties.lvl2Did not vote  0.03 0.03    89.15  1.25 0.216 -0.02  0.08
## 49            environ.lvl1:engagement.lvl1:all.parties.lvl2Don't know  0.04 0.04   450.31  0.93 0.354 -0.04  0.12
## 50          environ.lvl1:engagement.lvl1:all.parties.lvl2Invalid vote  0.34 0.55 31115.62  0.62 0.537 -0.73  1.41
## 51                environ.lvl1:engagement.lvl1:all.parties.lvl2NE age  0.02 0.04   338.20  0.65 0.515 -0.05  0.10
## 52            environ.lvl1:engagement.lvl1:all.parties.lvl2NE citizen  0.13 0.04   228.36  3.47 0.001  0.06  0.21
## 53              environ.lvl1:engagement.lvl1:all.parties.lvl2NE other -0.04 0.07  1447.94 -0.55 0.579 -0.17  0.10
## 54             environ.lvl1:engagement.lvl1:all.parties.lvl2No answer -0.08 0.54 30988.47 -0.14 0.889 -1.14  0.99
## 55           environ.lvl1:engagement.lvl1:all.parties.lvl2Other party  0.05 0.02   146.84  2.17 0.032  0.00  0.10
## 56 environ.lvl1:engagement.lvl1:all.parties.lvl2Pro-environment party  0.05 0.04   281.08  1.33 0.185 -0.02  0.12
```

```r
(VC.H5.mod6<-getVC(H5.mod6))
```

```
##             grp                         var1                         var2      est_SD       est_SD2
## 1  voting.group                  (Intercept)                         <NA>  0.19774809  3.910431e-02
## 2  voting.group                 environ.lvl1                         <NA>  0.04074286  1.659981e-03
## 3  voting.group              engagement.lvl1                         <NA>  0.05485458  3.009024e-03
## 4  voting.group environ.lvl1:engagement.lvl1                         <NA>  0.03587722  1.287175e-03
## 5  voting.group                  (Intercept)                 environ.lvl1  0.05926903  4.775201e-04
## 6  voting.group                  (Intercept)              engagement.lvl1  0.15360858  1.666252e-03
## 7  voting.group                  (Intercept) environ.lvl1:engagement.lvl1 -0.07280733 -5.165427e-04
## 8  voting.group                 environ.lvl1              engagement.lvl1 -0.01718070 -3.839769e-05
## 9  voting.group                 environ.lvl1 environ.lvl1:engagement.lvl1 -0.17335628 -2.534019e-04
## 10 voting.group              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.44021252 -8.663514e-04
## 11        cntry                  (Intercept)                         <NA>  0.48293468  2.332259e-01
## 12        cntry                 environ.lvl1                         <NA>  0.03971972  1.577656e-03
## 13        cntry              engagement.lvl1                         <NA>  0.04737523  2.244412e-03
## 14        cntry environ.lvl1:engagement.lvl1                         <NA>  0.02368416  5.609392e-04
## 15        cntry                  (Intercept)                 environ.lvl1 -0.06217027 -1.192552e-03
## 16        cntry                  (Intercept)              engagement.lvl1  0.51179009  1.170932e-02
## 17        cntry                  (Intercept) environ.lvl1:engagement.lvl1 -0.30965167 -3.541765e-03
## 18        cntry                 environ.lvl1              engagement.lvl1  0.11015106  2.072746e-04
## 19        cntry                 environ.lvl1 environ.lvl1:engagement.lvl1  0.91162995  8.575958e-04
## 20        cntry              engagement.lvl1 environ.lvl1:engagement.lvl1 -0.30542525 -3.427000e-04
## 21     Residual                         <NA>                         <NA>  1.02567279  1.052005e+00
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    61.81 -3.26 0.002 -0.76 -0.18
## 2                                                      age  0.00 0.00 35561.36  2.89 0.004  0.00  0.00
## 3                                                   gender  0.07 0.01 35553.22  5.72 0.000  0.04  0.09
## 4                                                     educ  0.01 0.00 35669.76  6.60 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35654.28 -4.87 0.000 -0.08 -0.03
## 6                            occupClerical support workers -0.04 0.09 35512.26 -0.45 0.655 -0.21  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35521.65 -0.67 0.503 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35512.44  0.23 0.816 -0.15  0.19
## 9                                            occupManagers -0.02 0.09 35499.41 -0.22 0.830 -0.19  0.15
## 10                            occupOther: Not in paid work  0.13 0.09 35570.56  1.40 0.160 -0.05  0.31
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35517.31 -0.64 0.521 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35511.76  0.76 0.444 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35501.08 -0.31 0.757 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35509.95 -0.36 0.721 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35539.17 -0.53 0.598 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35505.93 -0.34 0.731 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35532.88 -0.24 0.808 -0.24  0.18
## 18                                            environ.lvl1  0.11 0.02   107.96  5.60 0.000  0.07  0.15
## 19                                         engagement.lvl1  0.02 0.03   142.26  0.53 0.597 -0.04  0.07
## 20                                             env.eng.int -0.04 0.02   144.61 -1.70 0.091 -0.08  0.01
## 21                            all.parties.lvl2Did not vote  0.44 0.06   182.69  6.95 0.000  0.32  0.57
## 22                              all.parties.lvl2Don't know  0.42 0.07   268.25  5.95 0.000  0.28  0.56
## 23                            all.parties.lvl2Invalid vote  0.47 0.35  1087.13  1.34 0.180 -0.22  1.17
## 24                                  all.parties.lvl2NE age  0.72 0.07   273.89 10.12 0.000  0.58  0.86
## 25                              all.parties.lvl2NE citizen  0.85 0.08   268.13 11.12 0.000  0.70  1.00
## 26                                all.parties.lvl2NE other  0.71 0.10   656.69  7.25 0.000  0.52  0.90
## 27                               all.parties.lvl2No answer  0.54 0.38  1390.55  1.44 0.149 -0.20  1.28
## 28                             all.parties.lvl2Other party  0.55 0.05   230.67 11.33 0.000  0.45  0.64
## 29                   all.parties.lvl2Pro-environment party  0.91 0.06   256.51 14.03 0.000  0.78  1.04
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    86.02 -0.01 0.991 -0.04  0.04
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   338.63  0.86 0.390 -0.04  0.10
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.06 0.39 25427.47  0.16 0.870 -0.70  0.83
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   237.66  1.31 0.193 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   233.93 -1.31 0.191 -0.11  0.02
## 35                   environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1246.82  0.67 0.505 -0.07  0.15
## 36                  environ.lvl1:all.parties.lvl2No answer  0.24 0.23  7561.30  1.05 0.292 -0.21  0.69
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   134.16  0.79 0.432 -0.02  0.05
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.02 0.03   195.40  0.74 0.463 -0.04  0.08
## 39            engagement.lvl1:all.parties.lvl2Did not vote  0.02 0.03   104.77  0.70 0.487 -0.04  0.09
## 40              engagement.lvl1:all.parties.lvl2Don't know -0.00 0.05   502.12 -0.09 0.926 -0.11  0.10
## 41            engagement.lvl1:all.parties.lvl2Invalid vote  0.12 0.45 18649.66  0.28 0.783 -0.76  1.01
## 42                  engagement.lvl1:all.parties.lvl2NE age -0.01 0.05   348.95 -0.15 0.881 -0.10  0.08
## 43              engagement.lvl1:all.parties.lvl2NE citizen -0.00 0.05   254.20 -0.10 0.919 -0.10  0.09
## 44                engagement.lvl1:all.parties.lvl2NE other -0.01 0.09  2622.73 -0.12 0.901 -0.19  0.17
## 45               engagement.lvl1:all.parties.lvl2No answer -0.68 0.72 31703.95 -0.94 0.345 -2.10  0.73
## 46             engagement.lvl1:all.parties.lvl2Other party  0.07 0.03   168.71  2.29 0.023  0.01  0.12
## 47   engagement.lvl1:all.parties.lvl2Pro-environment party  0.10 0.04   291.53  2.25 0.025  0.01  0.18
## 48                env.eng.int:all.parties.lvl2Did not vote  0.03 0.03    89.15  1.25 0.216 -0.02  0.08
## 49                  env.eng.int:all.parties.lvl2Don't know  0.04 0.04   450.31  0.93 0.354 -0.04  0.12
## 50                env.eng.int:all.parties.lvl2Invalid vote  0.34 0.55 31115.63  0.62 0.537 -0.73  1.41
## 51                      env.eng.int:all.parties.lvl2NE age  0.02 0.04   338.21  0.65 0.515 -0.05  0.10
## 52                  env.eng.int:all.parties.lvl2NE citizen  0.13 0.04   228.36  3.47 0.001  0.06  0.21
## 53                    env.eng.int:all.parties.lvl2NE other -0.04 0.07  1447.94 -0.55 0.579 -0.17  0.10
## 54                   env.eng.int:all.parties.lvl2No answer -0.08 0.54 30988.48 -0.14 0.889 -1.14  0.99
## 55                 env.eng.int:all.parties.lvl2Other party  0.05 0.02   146.84  2.17 0.032  0.00  0.10
## 56       env.eng.int:all.parties.lvl2Pro-environment party  0.05 0.04   281.08  1.33 0.185 -0.02  0.12
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
##          all.parties.lvl2 env.eng.int.trend         SE  df    asymp.LCL   asymp.UCL
## 1  Anti-immigration party      -0.036599505 0.02149166 Inf -0.078722375 0.005523365
## 2            Did not vote      -0.004134442 0.01648642 Inf -0.036447225 0.028178342
## 3              Don't know       0.003360580 0.03807472 Inf -0.071264505 0.077985664
## 4            Invalid vote       0.300733211 0.54566214 Inf -0.768744936 1.370211358
## 5                  NE age      -0.012011093 0.03192660 Inf -0.074586077 0.050563891
## 6              NE citizen       0.095857958 0.03246905 Inf  0.032219791 0.159496125
## 7                NE other      -0.074699048 0.06566682 Inf -0.203403642 0.054005545
## 8               No answer      -0.112330663 0.54089240 Inf -1.172460290 0.947798964
## 9             Other party       0.013466283 0.01126837 Inf -0.008619325 0.035551891
## 10  Pro-environment party       0.010419460 0.02905574 Inf -0.046528739 0.067367658
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
##                     group env.eng.int.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party             -0.04 0.02 0.0886 0.7972     -0.08      0.01
## 2            Did not vote              0.00 0.02 0.8020 1.0000     -0.04      0.03
## 3              Don't know              0.00 0.04 0.9297 1.0000     -0.07      0.08
## 4            Invalid vote              0.30 0.55 0.5815 1.0000     -0.77      1.37
## 5                  NE age             -0.01 0.03 0.7068 1.0000     -0.07      0.05
## 6              NE citizen              0.10 0.03 0.0032 0.0315      0.03      0.16
## 7                NE other             -0.07 0.07 0.2553 1.0000     -0.20      0.05
## 8               No answer             -0.11 0.54 0.8355 1.0000     -1.17      0.95
## 9             Other party              0.01 0.01 0.2321 1.0000     -0.01      0.04
## 10  Pro-environment party              0.01 0.03 0.7199 1.0000     -0.05      0.07
```

```r
write.csv2(H5.mod6.trends.tab,"H5.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.contrast<-data.frame(pairs(H5.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate         SE  df  z.ratio   p.value
## 1 Pro-environment party - Anti-immigration party 0.04701896 0.03539014 Inf 1.328589 0.1839834
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
##                                         contrast     estimate         SE  df    z.ratio    p.value
## 1           Other party - Anti-immigration party  0.050065788 0.02311458 Inf  2.1659834 0.09093737
## 2 Pro-environment party - Anti-immigration party  0.047018964 0.03539014 Inf  1.3285895 0.36796688
## 3            Pro-environment party - Other party -0.003046823 0.03025143 Inf -0.1007167 0.91977537
```





\newpage




# Exploratory analysis: Is the association linear?

## Data re-preparations


```r
##Remove two of the smallest voting groups

Rdat<-dat %>%
  filter(all.parties.lvl2!="Invalid vote" & all.parties.lvl2!="No answer")
```

## Model 1 (H1 selected model without random effect correlations)


```r
EX9RR.NCmod1<-lmer(refugees~
                     (environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1,data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod1)
```

```
## [1] FALSE
```

```r
(VC.EX9RR.NCmod1<-getVC(EX9RR.NCmod1))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.30996069 0.096075628
## 2 voting.group.1 environ.lvl1 <NA> 0.04306744 0.001854804
## 3          cntry  (Intercept) <NA> 0.49934698 0.249347405
## 4        cntry.1 environ.lvl1 <NA> 0.04019147 0.001615354
## 5       Residual         <NA> <NA> 1.02930332 1.059465322
```

```r
getFE(EX9RR.NCmod1)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.46  0.53 0.602 -0.21  0.36
## 2                                                      age  0.00 0.00 34257.10  4.47 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35544.08  4.68 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35691.90  7.51 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35610.06 -5.25 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35484.21 -0.48 0.633 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35494.53 -0.77 0.443 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35496.59  0.14 0.890 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35483.77 -0.21 0.834 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35651.66  1.53 0.127 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35492.07 -0.73 0.465 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35488.61  0.80 0.424 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35484.99 -0.35 0.725 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35490.18 -0.44 0.663 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35499.20 -0.61 0.539 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35484.38 -0.37 0.711 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35502.24 -0.23 0.816 -0.24  0.19
## 18                                            environ.lvl1  0.12 0.01    19.45 11.40 0.000  0.10  0.15
```

\newpage

## Model 2 (add squared fixed effect)


```r
EX9RR.NCmod2<-lmer(refugees~(environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2),data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))


anova(EX9RR.NCmod1,EX9RR.NCmod2)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod1: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1
## EX9RR.NCmod2: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod2:     I(environ.lvl1^2)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod1   23 104228 104424 -52091   104182                       
## EX9RR.NCmod2   24 104226 104430 -52089   104178 4.4154  1    0.03562 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX9RR.NCmod2)
```

```
## [1] FALSE
```

```r
(VC.EX9RR.NCmod2<-getVC(EX9RR.NCmod2))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.30979150 0.095970771
## 2 voting.group.1 environ.lvl1 <NA> 0.04341201 0.001884603
## 3          cntry  (Intercept) <NA> 0.49939985 0.249400211
## 4        cntry.1 environ.lvl1 <NA> 0.03921179 0.001537565
## 5       Residual         <NA> <NA> 1.02924080 1.059336634
```

```r
getFE(EX9RR.NCmod2)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.09 0.14    49.58  0.60 0.550 -0.20  0.37
## 2                                                      age  0.00 0.00 34256.05  4.49 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35544.88  4.59 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35691.75  7.54 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35611.25 -5.28 0.000 -0.09 -0.04
## 6                            occupClerical support workers -0.04 0.09 35484.53 -0.47 0.636 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35494.84 -0.76 0.445 -0.24  0.11
## 8                              occupElementary occupations  0.01 0.09 35497.02  0.14 0.885 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35484.32 -0.20 0.842 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35652.15  1.54 0.124 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35492.33 -0.73 0.466 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35489.04  0.80 0.421 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35485.50 -0.34 0.732 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35490.47 -0.43 0.665 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35499.57 -0.61 0.542 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35484.54 -0.37 0.710 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35502.45 -0.23 0.818 -0.24  0.19
## 18                                            environ.lvl1  0.12 0.01    19.50 11.65 0.000  0.10  0.15
## 19                                       I(environ.lvl1^2) -0.01 0.00 24591.48 -2.10 0.035 -0.02 -0.00
```

\newpage

## Model 3 (add fixed cubic effect)


```r
EX9RR.NCmod3<-lmer(refugees~(environ.lvl1||voting.group)+
                     (environ.lvl1||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
            I(environ.lvl1^2)+
              I(environ.lvl1^3),
            data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod3)
```

```
## [1] FALSE
```

```r
#test for cubic term only
anova(EX9RR.NCmod2,EX9RR.NCmod3)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod2: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod2:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod2:     I(environ.lvl1^2)
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod2   24 104226 104430 -52089   104178                       
## EX9RR.NCmod3   25 104224 104436 -52087   104174 3.9124  1    0.04793 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#test for both non-linear terms
anova(EX9RR.NCmod1,EX9RR.NCmod3)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod1: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod1:     cntry) + age + gender + educ + resid + occup + environ.lvl1
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod1   23 104228 104424 -52091   104182                       
## EX9RR.NCmod3   25 104224 104436 -52087   104174 8.3278  2    0.01555 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(VC.EX9RR.NCmod3<-getVC(EX9RR.NCmod3))
```

```
##              grp         var1 var2     est_SD     est_SD2
## 1   voting.group  (Intercept) <NA> 0.31111641 0.096793423
## 2 voting.group.1 environ.lvl1 <NA> 0.04337556 0.001881439
## 3          cntry  (Intercept) <NA> 0.49952815 0.249528372
## 4        cntry.1 environ.lvl1 <NA> 0.03898023 0.001519459
## 5       Residual         <NA> <NA> 1.02915958 1.059169442
```

```r
getFE(EX9RR.NCmod3)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.09 0.14    49.57  0.63 0.530 -0.20  0.38
## 2                                                      age  0.00 0.00 34279.68  4.50 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35544.61  4.57 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35693.02  7.56 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35610.42 -5.27 0.000 -0.09 -0.04
## 6                            occupClerical support workers -0.04 0.09 35483.10 -0.49 0.622 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35493.41 -0.79 0.431 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35495.61  0.12 0.903 -0.16  0.18
## 9                                            occupManagers -0.02 0.09 35483.10 -0.22 0.828 -0.19  0.15
## 10                            occupOther: Not in paid work  0.14 0.09 35651.27  1.52 0.129 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35490.73 -0.75 0.452 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35487.91  0.79 0.431 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35484.20 -0.36 0.717 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35489.32 -0.45 0.652 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35498.32 -0.62 0.533 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35483.37 -0.39 0.693 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35501.65 -0.24 0.808 -0.24  0.18
## 18                                            environ.lvl1  0.11 0.01    49.95  8.03 0.000  0.08  0.13
## 19                                       I(environ.lvl1^2) -0.01 0.00 22734.68 -2.56 0.010 -0.02 -0.00
## 20                                       I(environ.lvl1^3)  0.01 0.00 32672.25  1.98 0.048  0.00  0.01
```

## Model 4 (add random effects for non-linear terms by country)


```r
EX9RR.NCmod4<-lmer(refugees~(environ.lvl1||voting.group)+
                (environ.lvl1+
                   I(environ.lvl1^2)+
                   I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                  I(environ.lvl1^2)+
                  I(environ.lvl1^3),
                data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod4)
```

```
## [1] FALSE
```

```r
anova(EX9RR.NCmod3,EX9RR.NCmod4)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod4: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 + I(environ.lvl1^2) + 
## EX9RR.NCmod4:     I(environ.lvl1^3) || cntry) + age + gender + educ + resid + 
## EX9RR.NCmod4:     occup + environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX9RR.NCmod3   25 104224 104436 -52087   104174                         
## EX9RR.NCmod4   27 104206 104435 -52076   104152 22.138  2  1.559e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(VC.EX9RR.NCmod4<-getVC(EX9RR.NCmod4))
```

```
##              grp              var1 var2     est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.31245182 9.762614e-02
## 2 voting.group.1      environ.lvl1 <NA> 0.04243715 1.800912e-03
## 3          cntry       (Intercept) <NA> 0.48773044 2.378810e-01
## 4        cntry.1      environ.lvl1 <NA> 0.03329338 1.108449e-03
## 5        cntry.2 I(environ.lvl1^2) <NA> 0.02518698 6.343839e-04
## 6        cntry.3 I(environ.lvl1^3) <NA> 0.00926153 8.577593e-05
## 7       Residual              <NA> <NA> 1.02840797 1.057623e+00
```

```r
getFE(EX9RR.NCmod4)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.09 0.14    51.04  0.63 0.534 -0.19  0.37
## 2                                                      age  0.00 0.00 34298.25  4.51 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35539.59  4.63 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35688.76  7.54 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35600.99 -5.25 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35475.45 -0.49 0.625 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35485.47 -0.78 0.435 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35486.43  0.13 0.900 -0.16  0.18
## 9                                            occupManagers -0.02 0.09 35473.91 -0.23 0.822 -0.19  0.15
## 10                            occupOther: Not in paid work  0.14 0.09 35644.31  1.55 0.122 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35481.49 -0.75 0.456 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35477.85  0.79 0.427 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35473.11 -0.37 0.713 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35479.95 -0.44 0.657 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35490.33 -0.62 0.535 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35474.97 -0.38 0.700 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35489.69 -0.23 0.816 -0.24  0.19
## 18                                            environ.lvl1  0.11 0.01    35.69  8.31 0.000  0.08  0.13
## 19                                       I(environ.lvl1^2) -0.01 0.01    20.82 -1.19 0.248 -0.02  0.01
## 20                                       I(environ.lvl1^3)  0.01 0.00    38.41  1.77 0.085 -0.00  0.01
```
## Model 5 (add random effects for non-linear terms by voting group)


```r
EX9RR.NCmod6<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3),data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))

isSingular(EX9RR.NCmod6)
```

```
## [1] FALSE
```

```r
#test against fixed effects only
anova(EX9RR.NCmod3,EX9RR.NCmod6)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod3: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 || 
## EX9RR.NCmod3:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod3:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX9RR.NCmod3   25 104224 104436 -52087   104174                         
## EX9RR.NCmod6   29 104200 104446 -52071   104142 32.135  4  1.795e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#test against random effects by country
anova(EX9RR.NCmod4,EX9RR.NCmod6)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod4: refugees ~ (environ.lvl1 || voting.group) + (environ.lvl1 + I(environ.lvl1^2) + 
## EX9RR.NCmod4:     I(environ.lvl1^3) || cntry) + age + gender + educ + resid + 
## EX9RR.NCmod4:     occup + environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX9RR.NCmod4   27 104206 104435 -52076   104152                        
## EX9RR.NCmod6   29 104200 104446 -52071   104142 9.9974  2   0.006747 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(VC.EX9RR.NCmod6<-getVC(EX9RR.NCmod6))
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.307426195 9.451087e-02
## 2 voting.group.1      environ.lvl1 <NA> 0.014478170 2.096174e-04
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.031444852 9.887787e-04
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.012955276 1.678392e-04
## 5          cntry       (Intercept) <NA> 0.489157753 2.392753e-01
## 6        cntry.1      environ.lvl1 <NA> 0.035967289 1.293646e-03
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.022173861 4.916801e-04
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.008474988 7.182542e-05
## 9       Residual              <NA> <NA> 1.027514662 1.055786e+00
```

```r
getFE(EX9RR.NCmod6)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    50.85  0.60 0.550 -0.20  0.37
## 2                                                      age  0.00 0.00 34169.44  4.52 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35517.87  4.64 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35678.33  7.53 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35581.98 -5.26 0.000 -0.09 -0.04
## 6                            occupClerical support workers -0.04 0.09 35455.22 -0.46 0.648 -0.21  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35465.85 -0.76 0.445 -0.24  0.11
## 8                              occupElementary occupations  0.01 0.09 35465.83  0.14 0.887 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35454.42 -0.20 0.838 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35622.13  1.56 0.119 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35460.83 -0.72 0.474 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35454.69  0.81 0.416 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35453.15 -0.35 0.729 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35459.41 -0.42 0.671 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35456.73 -0.59 0.553 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35453.44 -0.36 0.716 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35462.25 -0.23 0.818 -0.23  0.19
## 18                                            environ.lvl1  0.10 0.01    36.08  8.00 0.000  0.08  0.13
## 19                                       I(environ.lvl1^2) -0.01 0.01    20.33 -1.09 0.290 -0.02  0.01
## 20                                       I(environ.lvl1^3)  0.01 0.00    38.52  2.17 0.036  0.00  0.02
```

\newpage

### Marginal effects for non-linearity at level-1 


```r
EX9RR.NCmod6.trends<-
  emtrends(EX9RR.NCmod6,
           specs = c("environ.lvl1"),
           var=c("environ.lvl1"),
           at=list(environ.lvl1=
                     c(mean(Rdat$environ.lvl1)-1*sd(Rdat$environ.lvl1),
                       mean(Rdat$environ.lvl1)-0*sd(Rdat$environ.lvl1),
                       mean(Rdat$environ.lvl1)+1*sd(Rdat$environ.lvl1))))

(EX9RR.NCmod6.trends.tab<-data.frame(EX9RR.NCmod6.trends))
```

```
##   environ.lvl1 environ.lvl1.trend         SE  df  asymp.LCL asymp.UCL
## 1 -1.175205665          0.1527948 0.02196623 Inf 0.10974177 0.1958478
## 2  0.004159972          0.1014234 0.01270349 Inf 0.07652502 0.1263218
## 3  1.183525609          0.1176032 0.02108788 Inf 0.07627176 0.1589347
```

```r
EX9RR.NCmod6.trends.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod6.trends.tab$environ.lvl1.trend/
                   EX9RR.NCmod6.trends.tab$SE)))
EX9RR.NCmod6.trends.tab$adj.p<-
  p.adjust(EX9RR.NCmod6.trends.tab$p,method="holm")

EX9RR.NCmod6.trends.tab<-
  cbind(group=EX9RR.NCmod6.trends.tab[,1],
        round(EX9RR.NCmod6.trends.tab[,c(2,3)],2),
        round(EX9RR.NCmod6.trends.tab[,c(7,8)],4),
        round(EX9RR.NCmod6.trends.tab[,c(5,6)],2))
EX9RR.NCmod6.trends.tab
```

```
##          group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -1.175205665               0.15 0.02 0     0      0.11      0.20
## 2  0.004159972               0.10 0.01 0     0      0.08      0.13
## 3  1.183525609               0.12 0.02 0     0      0.08      0.16
```

```r
write.csv2(EX9RR.NCmod6.trends.tab,"EX9RR.NCmod6.trends.tab.csv")
```

\newpage

## Model 6 (include voting group main effects at level-2)


```r
EX9RR.NCmod7<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3)+
                     all.parties.lvl2,data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
isSingular(EX9RR.NCmod7)
```

```
## [1] FALSE
```

```r
anova(EX9RR.NCmod6,EX9RR.NCmod7)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod6: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod6:     I(environ.lvl1^2) + I(environ.lvl1^3)
## EX9RR.NCmod7: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod7:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2
##              npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)    
## EX9RR.NCmod6   29 104200 104446 -52071   104142                        
## EX9RR.NCmod7   36 104033 104338 -51980   103961 181.4  7  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod7)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.45 0.14    64.50 -3.16 0.002 -0.74 -0.17
## 2                                                      age  0.00 0.00 35629.81  4.94 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35545.28  4.61 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35686.35  7.58 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35639.40 -5.17 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35519.88 -0.50 0.616 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35530.79 -0.79 0.432 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35525.28  0.09 0.928 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35519.15 -0.24 0.812 -0.19  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35576.97  1.35 0.178 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35525.84 -0.75 0.452 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35518.30  0.77 0.441 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35509.14 -0.41 0.684 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35522.67 -0.45 0.653 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35526.17 -0.64 0.523 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35518.08 -0.40 0.693 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35508.68 -0.30 0.763 -0.24  0.18
## 18                                            environ.lvl1  0.10 0.01    36.24  7.97 0.000  0.08  0.13
## 19                                       I(environ.lvl1^2) -0.01 0.01    20.54 -1.19 0.248 -0.02  0.01
## 20                                       I(environ.lvl1^3)  0.01 0.00    38.86  2.22 0.033  0.00  0.02
## 21                            all.parties.lvl2Did not vote  0.45 0.06   180.45  7.04 0.000  0.32  0.58
## 22                              all.parties.lvl2Don't know  0.42 0.07   264.77  5.84 0.000  0.28  0.56
## 23                                  all.parties.lvl2NE age  0.74 0.07   268.08 10.41 0.000  0.60  0.88
## 24                              all.parties.lvl2NE citizen  0.86 0.08   266.89 11.15 0.000  0.71  1.01
## 25                                all.parties.lvl2NE other  0.71 0.10   654.98  7.24 0.000  0.52  0.91
## 26                             all.parties.lvl2Other party  0.54 0.05   228.58 11.11 0.000  0.44  0.64
## 27                   all.parties.lvl2Pro-environment party  0.91 0.07   254.07 13.95 0.000  0.78  1.04
```

```r
getVC(EX9RR.NCmod7)
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.196817449 0.0387371082
## 2 voting.group.1      environ.lvl1 <NA> 0.017962159 0.0003226391
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.025058248 0.0006279158
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.012426644 0.0001544215
## 5          cntry       (Intercept) <NA> 0.472145737 0.2229215971
## 6        cntry.1      environ.lvl1 <NA> 0.036060765 0.0013003787
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.023060706 0.0005317962
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.008334273 0.0000694601
## 9       Residual              <NA> <NA> 1.027756095 1.0562825916
```

## Model 7 (voting group interactions with linear term)


```r
EX9RR.NCmod8<-lmer(refugees~
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||voting.group)+
                     (environ.lvl1+
                        I(environ.lvl1^2)+
                        I(environ.lvl1^3)||cntry)+
                     age+gender+educ+resid+occup+
                     environ.lvl1+
                     I(environ.lvl1^2)+
                     I(environ.lvl1^3)+
                     all.parties.lvl2+
                     all.parties.lvl2:environ.lvl1,
                   data=Rdat,REML=F,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e8)))
isSingular(EX9RR.NCmod8)
```

```
## [1] FALSE
```

```r
anova(EX9RR.NCmod7,EX9RR.NCmod8)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod7: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod7:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod7:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
##              npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX9RR.NCmod7   36 104033 104338 -51980   103961                     
## EX9RR.NCmod8   43 104040 104405 -51977   103954 6.2321  7     0.5129
```

```r
getFE(EX9RR.NCmod8)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.45 0.14    64.47 -3.16 0.002 -0.74 -0.17
## 2                                                      age  0.00 0.00 35638.53  4.95 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35541.82  4.63 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35685.52  7.57 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35638.50 -5.16 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35515.74 -0.51 0.612 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35527.30 -0.79 0.428 -0.24  0.10
## 8                              occupElementary occupations  0.01 0.09 35521.23  0.08 0.937 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35515.04 -0.25 0.805 -0.20  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35574.34  1.35 0.178 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35522.26 -0.76 0.448 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35512.67  0.76 0.445 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35501.04 -0.41 0.684 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35519.50 -0.46 0.649 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35522.00 -0.64 0.520 -0.24  0.12
## 16            occupTechnicians and associate professionals -0.03 0.09 35514.61 -0.40 0.690 -0.21  0.14
## 17                                         occupUnemployed -0.03 0.11 35508.13 -0.30 0.767 -0.24  0.18
## 18                                            environ.lvl1  0.09 0.02   114.32  4.36 0.000  0.05  0.13
## 19                                       I(environ.lvl1^2) -0.01 0.01    20.71 -1.10 0.283 -0.02  0.01
## 20                                       I(environ.lvl1^3)  0.01 0.00    38.85  2.19 0.034  0.00  0.02
## 21                            all.parties.lvl2Did not vote  0.45 0.06   180.33  7.05 0.000  0.32  0.58
## 22                              all.parties.lvl2Don't know  0.42 0.07   264.74  5.84 0.000  0.28  0.56
## 23                                  all.parties.lvl2NE age  0.74 0.07   268.03 10.41 0.000  0.60  0.88
## 24                              all.parties.lvl2NE citizen  0.86 0.08   266.87 11.15 0.000  0.71  1.01
## 25                                all.parties.lvl2NE other  0.71 0.10   655.71  7.25 0.000  0.52  0.91
## 26                             all.parties.lvl2Other party  0.54 0.05   228.50 11.12 0.000  0.44  0.64
## 27                   all.parties.lvl2Pro-environment party  0.91 0.07   254.10 13.95 0.000  0.78  1.04
## 28               environ.lvl1:all.parties.lvl2Did not vote  0.01 0.02    71.83  0.24 0.811 -0.04  0.05
## 29                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   340.56  0.89 0.376 -0.04  0.10
## 30                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   229.65  1.21 0.227 -0.02  0.09
## 31                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   235.22 -1.19 0.234 -0.10  0.03
## 32                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1280.12  0.46 0.645 -0.08  0.14
## 33                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   115.57  0.92 0.358 -0.02  0.05
## 34      environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   192.73  0.93 0.353 -0.03  0.08
```

```r
getVC(EX9RR.NCmod8)
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.196553680 0.0386333490
## 2 voting.group.1      environ.lvl1 <NA> 0.014829854 0.0002199246
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.025780848 0.0006646521
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.012105853 0.0001465517
## 5          cntry       (Intercept) <NA> 0.472166732 0.2229414229
## 6        cntry.1      environ.lvl1 <NA> 0.037215243 0.0013849743
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.023113385 0.0005342286
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.008449817 0.0000713994
## 9       Residual              <NA> <NA> 1.027708667 1.0561851034
```

\newpage

## Model 8 (interactions with non-linear terms)


```r
EX9RR.NCmod10<-lmer(refugees~
                      (environ.lvl1+
                         I(environ.lvl1^2)+
                         I(environ.lvl1^3)||voting.group)+
                      (environ.lvl1+
                         I(environ.lvl1^2)+
                         I(environ.lvl1^3)||cntry)+
                      age+gender+educ+resid+occup+
                      environ.lvl1+
                      I(environ.lvl1^2)+
                      I(environ.lvl1^3)+
                      all.parties.lvl2+
                      all.parties.lvl2:environ.lvl1+
                      all.parties.lvl2:I(environ.lvl1^2)+
                      all.parties.lvl2:I(environ.lvl1^3),
                    data=Rdat,REML=F,
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e8)))

#test against linear interaction by voting groups
anova(EX9RR.NCmod8,EX9RR.NCmod10)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
## EX9RR.NCmod10: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod10:     all.parties.lvl2:environ.lvl1 + all.parties.lvl2:I(environ.lvl1^2) + 
## EX9RR.NCmod10:     all.parties.lvl2:I(environ.lvl1^3)
##               npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod8    43 104040 104405 -51977   103954                       
## EX9RR.NCmod10   57 104044 104527 -51965   103930 24.497 14    0.03986 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod10)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.40 0.14    66.64 -2.76 0.008 -0.69 -0.11
## 2                                                      age  0.00 0.00 35642.44  4.90 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35547.31  4.62 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35684.40  7.53 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35641.28 -5.18 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.05 0.09 35519.83 -0.53 0.595 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35530.07 -0.82 0.412 -0.24  0.10
## 8                              occupElementary occupations  0.00 0.09 35522.57  0.05 0.959 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35516.46 -0.28 0.779 -0.20  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35575.67  1.31 0.192 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35525.19 -0.79 0.428 -0.24  0.10
## 12                                      occupProfessionals  0.06 0.09 35516.86  0.74 0.458 -0.11  0.24
## 13                                            occupRetired -0.04 0.10 35503.40 -0.43 0.664 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35522.32 -0.48 0.629 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35525.29 -0.68 0.498 -0.25  0.12
## 16            occupTechnicians and associate professionals -0.04 0.09 35518.41 -0.43 0.669 -0.21  0.13
## 17                                         occupUnemployed -0.04 0.11 35512.22 -0.35 0.725 -0.25  0.17
## 18                                            environ.lvl1  0.05 0.03   343.91  1.74 0.083 -0.01  0.11
## 19                                       I(environ.lvl1^2) -0.05 0.02   234.60 -3.33 0.001 -0.09 -0.02
## 20                                       I(environ.lvl1^3)  0.02 0.01   625.02  2.42 0.016  0.00  0.05
## 21                            all.parties.lvl2Did not vote  0.41 0.07   220.90  6.08 0.000  0.28  0.55
## 22                              all.parties.lvl2Don't know  0.34 0.08   388.17  4.32 0.000  0.19  0.50
## 23                                  all.parties.lvl2NE age  0.65 0.08   354.42  8.49 0.000  0.50  0.80
## 24                              all.parties.lvl2NE citizen  0.77 0.08   377.20  9.15 0.000  0.61  0.94
## 25                                all.parties.lvl2NE other  0.64 0.12  1461.99  5.14 0.000  0.39  0.88
## 26                             all.parties.lvl2Other party  0.49 0.05   305.30  9.22 0.000  0.38  0.59
## 27                   all.parties.lvl2Pro-environment party  0.83 0.07   342.20 11.69 0.000  0.69  0.97
## 28               environ.lvl1:all.parties.lvl2Did not vote  0.04 0.04   251.68  1.20 0.230 -0.03  0.12
## 29                 environ.lvl1:all.parties.lvl2Don't know  0.06 0.06  1539.07  1.06 0.288 -0.05  0.18
## 30                     environ.lvl1:all.parties.lvl2NE age -0.00 0.05  1020.29 -0.07 0.943 -0.11  0.10
## 31                 environ.lvl1:all.parties.lvl2NE citizen  0.00 0.06   961.51  0.03 0.973 -0.12  0.12
## 32                   environ.lvl1:all.parties.lvl2NE other  0.02 0.11  7763.51  0.20 0.841 -0.20  0.24
## 33                environ.lvl1:all.parties.lvl2Other party  0.05 0.03   403.27  1.52 0.130 -0.01  0.11
## 34      environ.lvl1:all.parties.lvl2Pro-environment party  0.10 0.05   619.67  2.15 0.032  0.01  0.20
## 35          I(environ.lvl1^2):all.parties.lvl2Did not vote  0.03 0.02   149.19  1.71 0.089 -0.01  0.07
## 36            I(environ.lvl1^2):all.parties.lvl2Don't know  0.07 0.03   573.95  2.29 0.022  0.01  0.13
## 37                I(environ.lvl1^2):all.parties.lvl2NE age  0.08 0.02   336.73  3.39 0.001  0.03  0.13
## 38            I(environ.lvl1^2):all.parties.lvl2NE citizen  0.07 0.03   342.72  2.63 0.009  0.02  0.12
## 39              I(environ.lvl1^2):all.parties.lvl2NE other  0.06 0.05  2545.09  1.28 0.202 -0.03  0.16
## 40           I(environ.lvl1^2):all.parties.lvl2Other party  0.05 0.02   231.68  2.93 0.004  0.02  0.08
## 41 I(environ.lvl1^2):all.parties.lvl2Pro-environment party  0.07 0.02   230.80  2.72 0.007  0.02  0.11
## 42          I(environ.lvl1^3):all.parties.lvl2Did not vote -0.02 0.01   472.82 -1.47 0.144 -0.04  0.01
## 43            I(environ.lvl1^3):all.parties.lvl2Don't know -0.02 0.02  1254.72 -0.89 0.375 -0.06  0.02
## 44                I(environ.lvl1^3):all.parties.lvl2NE age  0.01 0.02  1212.45  0.72 0.474 -0.02  0.05
## 45            I(environ.lvl1^3):all.parties.lvl2NE citizen -0.02 0.02  1196.65 -1.04 0.300 -0.06  0.02
## 46              I(environ.lvl1^3):all.parties.lvl2NE other -0.01 0.03  3207.43 -0.17 0.864 -0.07  0.06
## 47           I(environ.lvl1^3):all.parties.lvl2Other party -0.02 0.01   753.38 -1.45 0.149 -0.04  0.01
## 48 I(environ.lvl1^3):all.parties.lvl2Pro-environment party -0.03 0.01   536.48 -1.93 0.054 -0.06  0.00
```

```r
getVC(EX9RR.NCmod10)
```

```
##              grp              var1 var2      est_SD      est_SD2
## 1   voting.group       (Intercept) <NA> 0.198599864 3.944191e-02
## 2 voting.group.1      environ.lvl1 <NA> 0.017750503 3.150804e-04
## 3 voting.group.2 I(environ.lvl1^2) <NA> 0.019833792 3.933793e-04
## 4 voting.group.3 I(environ.lvl1^3) <NA> 0.011536641 1.330941e-04
## 5          cntry       (Intercept) <NA> 0.471895023 2.226849e-01
## 6        cntry.1      environ.lvl1 <NA> 0.036861258 1.358752e-03
## 7        cntry.2 I(environ.lvl1^2) <NA> 0.024002372 5.761139e-04
## 8        cntry.3 I(environ.lvl1^3) <NA> 0.008735925 7.631638e-05
## 9       Residual              <NA> <NA> 1.027516886 1.055791e+00
```

\newpage

### Marginal trends for each voting group

#### Linear coefficients


```r
EX9RR.NCmod10.linear<-
  emtrends(EX9RR.NCmod10,specs = c("all.parties.lvl2"),var=c("environ.lvl1"))
EX9RR.NCmod10.linear.tab<-data.frame(EX9RR.NCmod10.linear)

EX9RR.NCmod10.linear.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.linear.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.linear.tab$SE)))
EX9RR.NCmod10.linear.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.linear.tab$p,method="holm")

EX9RR.NCmod10.linear.tab<-
  cbind.data.frame(group=EX9RR.NCmod10.linear.tab[,"all.parties.lvl2"],
                   beta=round_tidy(EX9RR.NCmod10.linear.tab[,"environ.lvl1.trend"],2),
                   CI=paste0("[",
                             round_tidy(EX9RR.NCmod10.linear.tab[,"asymp.LCL"],2),
                             ", ",
                             round_tidy(EX9RR.NCmod10.linear.tab[,"asymp.UCL"],2),
                             "]"),
                   p=round_tidy(EX9RR.NCmod10.linear.tab[,"p"],3),
                   p.adj=round_tidy(EX9RR.NCmod10.linear.tab[,"adj.p"],3),
                   p_less_001=ifelse(EX9RR.NCmod10.linear.tab[,"p"]<.001,"yes","no"),
                   p.adj_less_001=ifelse(EX9RR.NCmod10.linear.tab[,"adj.p"]<.001,"yes","no"))

EX9RR.NCmod10.linear.tab
```

```
##                    group beta            CI     p p.adj p_less_001 p.adj_less_001
## 1 Anti-immigration party 0.05 [-0.01, 0.11] 0.088 0.351         no             no
## 2           Did not vote 0.10  [0.05, 0.14] 0.000 0.000        yes            yes
## 3             Don't know 0.12  [0.01, 0.22] 0.028 0.138         no             no
## 4                 NE age 0.05 [-0.04, 0.14] 0.273 0.818         no             no
## 5             NE citizen 0.06 [-0.05, 0.16] 0.305 0.818         no             no
## 6               NE other 0.08 [-0.14, 0.29] 0.485 0.818         no             no
## 7            Other party 0.10  [0.07, 0.14] 0.000 0.000        yes            yes
## 8  Pro-environment party 0.16  [0.08, 0.23] 0.000 0.000        yes            yes
```

```r
write.csv2(EX9RR.NCmod10.linear.tab,"EX9RR.NCmod10.linear.tab.csv")


#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.linear, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast   estimate         SE  df  z.ratio    p.value
## 1           Other party - Anti-immigration party 0.05056946 0.03294233 Inf 1.535090 0.12476163
## 2 Pro-environment party - Anti-immigration party 0.10397477 0.04785163 Inf 2.172857 0.02979105
## 3            Pro-environment party - Other party 0.05340531 0.03968538 Inf 1.345717 0.17839368
```

### Quadratic and cubic coefficients

* Non-linear variables must be manually recoded to new variables for extracting marginal effects


```r
#recoding
Rdat$environ.lvl1.sq<-Rdat$environ.lvl1^2
Rdat$environ.lvl1.cu<-Rdat$environ.lvl1^3

#refitting the model
EX9RR.NCmod10.ma<-lmer(refugees~
                         (environ.lvl1+
                            environ.lvl1.sq+
                            environ.lvl1.cu||voting.group)+
                         (environ.lvl1+
                            environ.lvl1.sq+
                            environ.lvl1.cu||cntry)+
                         age+gender+educ+resid+occup+
                         environ.lvl1+
                         environ.lvl1.sq+
                         environ.lvl1.cu+
                         all.parties.lvl2+
                         all.parties.lvl2:environ.lvl1+
                         all.parties.lvl2:environ.lvl1.sq+
                       all.parties.lvl2:environ.lvl1.cu,
                       data=Rdat,REML=F,
                       control=lmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e8)))
isSingular(EX9RR.NCmod10.ma)
```

```
## [1] FALSE
```

```r
#check if identical
anova(EX9RR.NCmod10,EX9RR.NCmod10.ma)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod10: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod10:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod10:     all.parties.lvl2:environ.lvl1 + all.parties.lvl2:I(environ.lvl1^2) + 
## EX9RR.NCmod10:     all.parties.lvl2:I(environ.lvl1^3)
## EX9RR.NCmod10.ma: refugees ~ (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     voting.group) + (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10.ma:     environ.lvl1.sq + environ.lvl1.cu + all.parties.lvl2 + all.parties.lvl2:environ.lvl1 + 
## EX9RR.NCmod10.ma:     all.parties.lvl2:environ.lvl1.sq + all.parties.lvl2:environ.lvl1.cu
##                  npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)
## EX9RR.NCmod10      57 104044 104527 -51965   103930                    
## EX9RR.NCmod10.ma   57 104044 104527 -51965   103930     0  0          1
```

```r
#
anova(EX9RR.NCmod8,EX9RR.NCmod10.ma)
```

```
## Data: Rdat
## Models:
## EX9RR.NCmod8: refugees ~ (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     voting.group) + (environ.lvl1 + I(environ.lvl1^2) + I(environ.lvl1^3) || 
## EX9RR.NCmod8:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod8:     I(environ.lvl1^2) + I(environ.lvl1^3) + all.parties.lvl2 + 
## EX9RR.NCmod8:     all.parties.lvl2:environ.lvl1
## EX9RR.NCmod10.ma: refugees ~ (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     voting.group) + (environ.lvl1 + environ.lvl1.sq + environ.lvl1.cu || 
## EX9RR.NCmod10.ma:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## EX9RR.NCmod10.ma:     environ.lvl1.sq + environ.lvl1.cu + all.parties.lvl2 + all.parties.lvl2:environ.lvl1 + 
## EX9RR.NCmod10.ma:     all.parties.lvl2:environ.lvl1.sq + all.parties.lvl2:environ.lvl1.cu
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## EX9RR.NCmod8       43 104040 104405 -51977   103954                       
## EX9RR.NCmod10.ma   57 104044 104527 -51965   103930 24.497 14    0.03986 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
getFE(EX9RR.NCmod10.ma)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.40 0.14    66.64 -2.76 0.008 -0.69 -0.11
## 2                                                      age  0.00 0.00 35642.44  4.90 0.000  0.00  0.00
## 3                                                   gender  0.05 0.01 35547.31  4.62 0.000  0.03  0.08
## 4                                                     educ  0.01 0.00 35684.40  7.53 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35641.28 -5.18 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.05 0.09 35519.83 -0.53 0.595 -0.22  0.13
## 7                    occupCraft and related trades workers -0.07 0.09 35530.07 -0.82 0.412 -0.24  0.10
## 8                              occupElementary occupations  0.00 0.09 35522.57  0.05 0.959 -0.17  0.18
## 9                                            occupManagers -0.02 0.09 35516.46 -0.28 0.779 -0.20  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35575.67  1.31 0.192 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.07 0.09 35525.19 -0.79 0.428 -0.24  0.10
## 12                                      occupProfessionals  0.06 0.09 35516.86  0.74 0.458 -0.11  0.24
## 13                                            occupRetired -0.04 0.10 35503.40 -0.43 0.664 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35522.32 -0.48 0.629 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35525.29 -0.68 0.498 -0.25  0.12
## 16            occupTechnicians and associate professionals -0.04 0.09 35518.41 -0.43 0.669 -0.21  0.13
## 17                                         occupUnemployed -0.04 0.11 35512.22 -0.35 0.725 -0.25  0.17
## 18                                            environ.lvl1  0.05 0.03   343.91  1.74 0.083 -0.01  0.11
## 19                                         environ.lvl1.sq -0.05 0.02   234.60 -3.33 0.001 -0.09 -0.02
## 20                                         environ.lvl1.cu  0.02 0.01   625.02  2.42 0.016  0.00  0.05
## 21                            all.parties.lvl2Did not vote  0.41 0.07   220.90  6.08 0.000  0.28  0.55
## 22                              all.parties.lvl2Don't know  0.34 0.08   388.17  4.32 0.000  0.19  0.50
## 23                                  all.parties.lvl2NE age  0.65 0.08   354.42  8.49 0.000  0.50  0.80
## 24                              all.parties.lvl2NE citizen  0.77 0.08   377.20  9.15 0.000  0.61  0.94
## 25                                all.parties.lvl2NE other  0.64 0.12  1461.99  5.14 0.000  0.39  0.88
## 26                             all.parties.lvl2Other party  0.49 0.05   305.30  9.22 0.000  0.38  0.59
## 27                   all.parties.lvl2Pro-environment party  0.83 0.07   342.20 11.69 0.000  0.69  0.97
## 28               environ.lvl1:all.parties.lvl2Did not vote  0.04 0.04   251.68  1.20 0.230 -0.03  0.12
## 29                 environ.lvl1:all.parties.lvl2Don't know  0.06 0.06  1539.07  1.06 0.288 -0.05  0.18
## 30                     environ.lvl1:all.parties.lvl2NE age -0.00 0.05  1020.29 -0.07 0.943 -0.11  0.10
## 31                 environ.lvl1:all.parties.lvl2NE citizen  0.00 0.06   961.51  0.03 0.973 -0.12  0.12
## 32                   environ.lvl1:all.parties.lvl2NE other  0.02 0.11  7763.51  0.20 0.841 -0.20  0.24
## 33                environ.lvl1:all.parties.lvl2Other party  0.05 0.03   403.27  1.52 0.130 -0.01  0.11
## 34      environ.lvl1:all.parties.lvl2Pro-environment party  0.10 0.05   619.67  2.15 0.032  0.01  0.20
## 35            environ.lvl1.sq:all.parties.lvl2Did not vote  0.03 0.02   149.19  1.71 0.089 -0.01  0.07
## 36              environ.lvl1.sq:all.parties.lvl2Don't know  0.07 0.03   573.95  2.29 0.022  0.01  0.13
## 37                  environ.lvl1.sq:all.parties.lvl2NE age  0.08 0.02   336.73  3.39 0.001  0.03  0.13
## 38              environ.lvl1.sq:all.parties.lvl2NE citizen  0.07 0.03   342.72  2.63 0.009  0.02  0.12
## 39                environ.lvl1.sq:all.parties.lvl2NE other  0.06 0.05  2545.09  1.28 0.202 -0.03  0.16
## 40             environ.lvl1.sq:all.parties.lvl2Other party  0.05 0.02   231.68  2.93 0.004  0.02  0.08
## 41   environ.lvl1.sq:all.parties.lvl2Pro-environment party  0.07 0.02   230.80  2.72 0.007  0.02  0.11
## 42            environ.lvl1.cu:all.parties.lvl2Did not vote -0.02 0.01   472.82 -1.47 0.144 -0.04  0.01
## 43              environ.lvl1.cu:all.parties.lvl2Don't know -0.02 0.02  1254.72 -0.89 0.375 -0.06  0.02
## 44                  environ.lvl1.cu:all.parties.lvl2NE age  0.01 0.02  1212.45  0.72 0.474 -0.02  0.05
## 45              environ.lvl1.cu:all.parties.lvl2NE citizen -0.02 0.02  1196.65 -1.04 0.300 -0.06  0.02
## 46                environ.lvl1.cu:all.parties.lvl2NE other -0.01 0.03  3207.43 -0.17 0.864 -0.07  0.06
## 47             environ.lvl1.cu:all.parties.lvl2Other party -0.02 0.01   753.38 -1.45 0.149 -0.04  0.01
## 48   environ.lvl1.cu:all.parties.lvl2Pro-environment party -0.03 0.01   536.48 -1.93 0.054 -0.06  0.00
```

```r
getVC(EX9RR.NCmod10.ma)
```

```
##              grp            var1 var2      est_SD      est_SD2
## 1   voting.group     (Intercept) <NA> 0.198599864 3.944191e-02
## 2 voting.group.1    environ.lvl1 <NA> 0.017750503 3.150804e-04
## 3 voting.group.2 environ.lvl1.sq <NA> 0.019833792 3.933793e-04
## 4 voting.group.3 environ.lvl1.cu <NA> 0.011536641 1.330941e-04
## 5          cntry     (Intercept) <NA> 0.471895023 2.226849e-01
## 6        cntry.1    environ.lvl1 <NA> 0.036861258 1.358752e-03
## 7        cntry.2 environ.lvl1.sq <NA> 0.024002372 5.761139e-04
## 8        cntry.3 environ.lvl1.cu <NA> 0.008735925 7.631638e-05
## 9       Residual            <NA> <NA> 1.027516886 1.055791e+00
```

#### Marginal quadratic effects


```r
EX9RR.NCmod10.quadratic<-
  emtrends(EX9RR.NCmod10.ma,
           specs = c("all.parties.lvl2"),var=c("environ.lvl1.sq"))
EX9RR.NCmod10.quadratic.tab<-data.frame(EX9RR.NCmod10.quadratic)

EX9RR.NCmod10.quadratic.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.quadratic.tab$environ.lvl1.sq.trend/
                   EX9RR.NCmod10.quadratic.tab$SE)))
EX9RR.NCmod10.quadratic.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.quadratic.tab$p,method="holm")

EX9RR.NCmod10.quadratic.tab<-
  cbind.data.frame(group=EX9RR.NCmod10.quadratic.tab[,"all.parties.lvl2"],
                   beta=round_tidy(EX9RR.NCmod10.quadratic.tab[,"environ.lvl1.sq.trend"],2),
                   CI=paste0("[",
                             round_tidy(EX9RR.NCmod10.quadratic.tab[,"asymp.LCL"],2),
                             ", ",
                             round_tidy(EX9RR.NCmod10.quadratic.tab[,"asymp.UCL"],2),
                             "]"),
                   p=round_tidy(EX9RR.NCmod10.quadratic.tab[,"p"],3),
                   p.adj=round_tidy(EX9RR.NCmod10.quadratic.tab[,"adj.p"],3),
                   p_less_001=ifelse(EX9RR.NCmod10.quadratic.tab[,"p"]<.001,"yes","no"),
                   p.adj_less_001=ifelse(EX9RR.NCmod10.quadratic.tab[,"adj.p"]<.001,"yes","no"))

EX9RR.NCmod10.quadratic.tab
```

```
##                    group  beta             CI     p p.adj p_less_001 p.adj_less_001
## 1 Anti-immigration party -0.05 [-0.09, -0.02] 0.001 0.007        yes             no
## 2           Did not vote -0.02  [-0.05, 0.00] 0.075 0.523         no             no
## 3             Don't know  0.01  [-0.04, 0.06] 0.606 1.000         no             no
## 4                 NE age  0.03  [-0.01, 0.07] 0.152 0.913         no             no
## 5             NE citizen  0.02  [-0.03, 0.06] 0.479 1.000         no             no
## 6               NE other  0.01  [-0.08, 0.10] 0.850 1.000         no             no
## 7            Other party -0.01  [-0.02, 0.01] 0.518 1.000         no             no
## 8  Pro-environment party  0.01  [-0.03, 0.05] 0.557 1.000         no             no
```

```r
write.csv2(EX9RR.NCmod10.quadratic.tab,"EX9RR.NCmod10.quadratic.tab.csv")

#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.quadratic, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast   estimate         SE  df   z.ratio     p.value
## 1           Other party - Anti-immigration party 0.04945800 0.01688706 Inf 2.9287517 0.003403261
## 2 Pro-environment party - Anti-immigration party 0.06620059 0.02431888 Inf 2.7221896 0.006485092
## 3            Pro-environment party - Other party 0.01674259 0.01949677 Inf 0.8587367 0.390485802
```

```r
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.quadratic, "del.eff",
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##        all.parties.lvl2_del.eff    estimate         SE  df   z.ratio     p.value
## 1 Anti-immigration party effect -0.05782929 0.01852725 Inf -3.121309 0.001800488
## 2            Other party effect  0.01635770 0.01359399 Inf  1.203304 0.228858667
## 3  Pro-environment party effect  0.04147159 0.02035860 Inf  2.037055 0.041644550
```

\newpage

#### Marginal cubic effects


```r
EX9RR.NCmod10.cubic<-
  emtrends(EX9RR.NCmod10.ma,specs = c("all.parties.lvl2"),var=c("environ.lvl1.cu"))
EX9RR.NCmod10.cubic.tab<-data.frame(EX9RR.NCmod10.cubic)


EX9RR.NCmod10.cubic.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.cubic.tab$environ.lvl1.cu.trend/
                   EX9RR.NCmod10.cubic.tab$SE)))
EX9RR.NCmod10.cubic.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.cubic.tab$p,method="holm")


EX9RR.NCmod10.cubic.tab<-
  cbind.data.frame(group=EX9RR.NCmod10.cubic.tab[,"all.parties.lvl2"],
                   beta=round_tidy(EX9RR.NCmod10.cubic.tab[,"environ.lvl1.cu.trend"],2),
                   CI=paste0("[",
                             round_tidy(EX9RR.NCmod10.cubic.tab[,"asymp.LCL"],2),
                             ", ",
                             round_tidy(EX9RR.NCmod10.cubic.tab[,"asymp.UCL"],2),
                             "]"),
                   p=round_tidy(EX9RR.NCmod10.cubic.tab[,"p"],3),
                   p.adj=round_tidy(EX9RR.NCmod10.cubic.tab[,"adj.p"],3),
                   p_less_001=ifelse(EX9RR.NCmod10.cubic.tab[,"p"]<.001,"yes","no"),
                   p.adj_less_001=ifelse(EX9RR.NCmod10.cubic.tab[,"adj.p"]<.001,"yes","no"))

EX9RR.NCmod10.cubic.tab
```

```
##                    group  beta            CI     p p.adj p_less_001 p.adj_less_001
## 1 Anti-immigration party  0.02  [0.00, 0.05] 0.016 0.110         no             no
## 2           Did not vote  0.01 [-0.01, 0.02] 0.400 1.000         no             no
## 3             Don't know  0.01 [-0.03, 0.04] 0.681 1.000         no             no
## 4                 NE age  0.04  [0.01, 0.07] 0.009 0.071         no             no
## 5             NE citizen  0.01 [-0.03, 0.04] 0.731 1.000         no             no
## 6               NE other  0.02 [-0.04, 0.08] 0.535 1.000         no             no
## 7            Other party  0.01 [-0.00, 0.02] 0.067 0.400         no             no
## 8  Pro-environment party -0.00 [-0.03, 0.02] 0.725 1.000         no             no
```

```r
write.csv2(EX9RR.NCmod10.cubic.tab,"EX9RR.NCmod10.cubic.tab.csv")

#contrast for three voting groups
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.cubic, 
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast    estimate         SE  df   z.ratio    p.value
## 1           Other party - Anti-immigration party -0.01601655 0.01108249 Inf -1.445212 0.14839830
## 2 Pro-environment party - Anti-immigration party -0.02886720 0.01492070 Inf -1.934708 0.05302616
## 3            Pro-environment party - Other party -0.01285065 0.01166814 Inf -1.101346 0.27074613
```

```r
(EX9RR.more.contrasts<-data.frame(pairs(EX9RR.NCmod10.cubic, "del.eff",
                                        exclude=c(2:6), by = NULL,adjust=c("none"),reverse=T)))
```

```
##        all.parties.lvl2_del.eff     estimate          SE  df    z.ratio    p.value
## 1 Anti-immigration party effect  0.022441877 0.011776593 Inf  1.9056340 0.05669771
## 2            Other party effect -0.001582947 0.008592243 Inf -0.1842298 0.85383318
## 3  Pro-environment party effect -0.020858929 0.012193483 Inf -1.7106621 0.08714350
```



\newpage

#### Compile Table 3


```r
# Compile Table 3

lin.coefs<-read.csv2("EX9RR.NCmod10.linear.tab.csv",
                     stringsAsFactors = F)

quad.coefs<-read.csv2("EX9RR.NCmod10.quadratic.tab.csv",
                     stringsAsFactors = F)

cubic.coefs<-read.csv2("EX9RR.NCmod10.cubic.tab.csv",
                      stringsAsFactors = F)

tab3<-cbind.data.frame(Voting_group=lin.coefs[,c("group")],
                 Linear=lin.coefs[,c("beta")],
                 CI=lin.coefs[,c("CI")],
                 p=lin.coefs[,c("p")],
                 Quadratic=quad.coefs[,c("beta")],
                 CI=quad.coefs[,c("CI")],
                 p=quad.coefs[,c("p")],
                 Cubic=cubic.coefs[,c("beta")],
                 CI=cubic.coefs[,c("CI")],
                 p=cubic.coefs[,c("p")])
tab3                 
```

```
##             Voting_group Linear            CI     p Quadratic             CI     p Cubic            CI     p
## 1 Anti-immigration party   0.05 [-0.01, 0.11] 0.088     -0.05 [-0.09, -0.02] 0.001  0.02  [0.00, 0.05] 0.016
## 2           Did not vote   0.10  [0.05, 0.14] 0.000     -0.02  [-0.05, 0.00] 0.075  0.01 [-0.01, 0.02] 0.400
## 3             Don't know   0.12  [0.01, 0.22] 0.028      0.01  [-0.04, 0.06] 0.606  0.01 [-0.03, 0.04] 0.681
## 4                 NE age   0.05 [-0.04, 0.14] 0.273      0.03  [-0.01, 0.07] 0.152  0.04  [0.01, 0.07] 0.009
## 5             NE citizen   0.06 [-0.05, 0.16] 0.305      0.02  [-0.03, 0.06] 0.479  0.01 [-0.03, 0.04] 0.731
## 6               NE other   0.08 [-0.14, 0.29] 0.485      0.01  [-0.08, 0.10] 0.850  0.02 [-0.04, 0.08] 0.535
## 7            Other party   0.10  [0.07, 0.14] 0.000     -0.01  [-0.02, 0.01] 0.518  0.01 [-0.00, 0.02] 0.067
## 8  Pro-environment party   0.16  [0.08, 0.23] 0.000      0.01  [-0.03, 0.05] 0.557 -0.00 [-0.03, 0.02] 0.725
```

```r
write.csv2(tab3,"tab3.csv")
```


\newpage

### Simple slopes for each voting group at -1SD, mean, and +1SD (grand mean used to enhance comparisons between voting groups)

* Because the location of points on grand-mean vary for group-means, the locations of grand mean points for each voting group are calculated separately

#### This script would produce the slopes for group means (not interpreted in the main text)


```r
EX9RR.NCmod10.slopes<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                               at=list(environ.lvl1=c(mean(dat$environ.lvl1)-sd(dat$environ.lvl1),
                                                      mean(dat$environ.lvl1)-0*sd(dat$environ.lvl1),
                                                      mean(dat$environ.lvl1)+sd(dat$environ.lvl1))))
(EX9RR.NCmod10.slopes.tab<-data.frame(EX9RR.NCmod10.slopes))
```

```
##    environ.lvl1       all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL asymp.UCL
## 1  -1.175144455 Anti-immigration party         0.28504089 0.05450375 Inf  0.178215492 0.3918663
## 2   0.004157179 Anti-immigration party         0.05275139 0.03088596 Inf -0.007783985 0.1132868
## 3   1.183458812 Anti-immigration party         0.02882925 0.03872705 Inf -0.047074369 0.1047329
## 4  -1.175144455           Did not vote         0.17687860 0.04063306 Inf  0.097239269 0.2565179
## 5   0.004157179           Did not vote         0.09787830 0.02384802 Inf  0.051137041 0.1446196
## 6   1.183458812           Did not vote         0.07333770 0.03185323 Inf  0.010906513 0.1357689
## 7  -1.175144455             Don't know         0.11605971 0.08371854 Inf -0.048025609 0.2801450
## 8   0.004157179             Don't know         0.11820422 0.05363854 Inf  0.013074621 0.2233338
## 9   1.183458812             Don't know         0.18015552 0.06320766 Inf  0.056270787 0.3040403
## 10 -1.175144455                 NE age         0.13767568 0.05327805 Inf  0.033252627 0.2420987
## 11  0.004157179                 NE age         0.05012206 0.04569817 Inf -0.039444701 0.1396888
## 12  1.183458812                 NE age         0.27498893 0.06211810 Inf  0.153239686 0.3967382
## 13 -1.175144455             NE citizen         0.04074452 0.06776664 Inf -0.092075660 0.1735647
## 14  0.004157179             NE citizen         0.05587648 0.05443604 Inf -0.050816198 0.1625692
## 15  1.183458812             NE citizen         0.11677156 0.06115238 Inf -0.003084893 0.2366280
## 16 -1.175144455               NE other         0.13458712 0.14560598 Inf -0.150795357 0.4199696
## 17  0.004157179               NE other         0.07618990 0.10908074 Inf -0.137604421 0.2899842
## 18  1.183458812               NE other         0.17940556 0.11223628 Inf -0.040573502 0.3993846
## 19 -1.175144455            Other party         0.15296094 0.02599519 Inf  0.102011302 0.2039106
## 20  0.004157179            Other party         0.10332058 0.01618173 Inf  0.071604968 0.1350362
## 21  1.183458812            Other party         0.12839724 0.02452572 Inf  0.080327712 0.1764668
## 22 -1.175144455  Pro-environment party         0.11392498 0.04354684 Inf  0.028574736 0.1992752
## 23  0.004157179  Pro-environment party         0.15672579 0.03824301 Inf  0.081770873 0.2316807
## 24  1.183458812  Pro-environment party         0.16701117 0.06485639 Inf  0.039894987 0.2941273
```

```r
EX9RR.NCmod10.slopes.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.tab$SE)))
EX9RR.NCmod10.slopes.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.tab$p,method="holm")

EX9RR.NCmod10.slopes.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.tab[,1],group=EX9RR.NCmod10.slopes.tab[,2],
        round(EX9RR.NCmod10.slopes.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.tab
```

```
##       env_point                  group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  -1.175144455 Anti-immigration party               0.29 0.05 0.0000 0.0000      0.18      0.39
## 2   0.004157179 Anti-immigration party               0.05 0.03 0.0876 0.7888     -0.01      0.11
## 3   1.183458812 Anti-immigration party               0.03 0.04 0.4566 1.0000     -0.05      0.10
## 4  -1.175144455           Did not vote               0.18 0.04 0.0000 0.0003      0.10      0.26
## 5   0.004157179           Did not vote               0.10 0.02 0.0000 0.0007      0.05      0.14
## 6   1.183458812           Did not vote               0.07 0.03 0.0213 0.2558      0.01      0.14
## 7  -1.175144455             Don't know               0.12 0.08 0.1657 1.0000     -0.05      0.28
## 8   0.004157179             Don't know               0.12 0.05 0.0275 0.3030      0.01      0.22
## 9   1.183458812             Don't know               0.18 0.06 0.0044 0.0699      0.06      0.30
## 10 -1.175144455                 NE age               0.14 0.05 0.0098 0.1367      0.03      0.24
## 11  0.004157179                 NE age               0.05 0.05 0.2727 1.0000     -0.04      0.14
## 12  1.183458812                 NE age               0.27 0.06 0.0000 0.0002      0.15      0.40
## 13 -1.175144455             NE citizen               0.04 0.07 0.5477 1.0000     -0.09      0.17
## 14  0.004157179             NE citizen               0.06 0.05 0.3047 1.0000     -0.05      0.16
## 15  1.183458812             NE citizen               0.12 0.06 0.0562 0.5620      0.00      0.24
## 16 -1.175144455               NE other               0.13 0.15 0.3553 1.0000     -0.15      0.42
## 17  0.004157179               NE other               0.08 0.11 0.4849 1.0000     -0.14      0.29
## 18  1.183458812               NE other               0.18 0.11 0.1099 0.8795     -0.04      0.40
## 19 -1.175144455            Other party               0.15 0.03 0.0000 0.0000      0.10      0.20
## 20  0.004157179            Other party               0.10 0.02 0.0000 0.0000      0.07      0.14
## 21  1.183458812            Other party               0.13 0.02 0.0000 0.0000      0.08      0.18
## 22 -1.175144455  Pro-environment party               0.11 0.04 0.0089 0.1334      0.03      0.20
## 23  0.004157179  Pro-environment party               0.16 0.04 0.0000 0.0007      0.08      0.23
## 24  1.183458812  Pro-environment party               0.17 0.06 0.0100 0.1367      0.04      0.29
```

```r
write.csv2(EX9RR.NCmod10.slopes.tab,"EX9RR.NCmod10.slopes.tab.csv")
```

\newpage

#### Locating the grand-mean points for group-means


```r
all.parties.lvl2.means<-Rdat %>%
  group_by(all.parties.lvl2) %>%
  summarize(env.lvl2.mean=mean(environ.cntrymc),
            env.lvl1.sd=sd(environ.lvl1),
            ref.lvl2.mean=mean(refugees.cntrymc),
            ref.lvl1.sd=sd(refugees.lvl1),
            n=n())

all.parties.lvl2.means<-data.frame(all.parties.lvl2.means)

all.parties.lvl2.means$env.num<-
  (all.parties.lvl2.means$n-1)*
  (all.parties.lvl2.means$env.lvl1.sd^2)

all.parties.lvl2.means$env.pooled.sd<-
  sqrt(sum(all.parties.lvl2.means$env.num)/
         (sum(all.parties.lvl2.means$n)-
            nrow(all.parties.lvl2.means)))

all.parties.lvl2.means$ref.num<-
  (all.parties.lvl2.means$n-1)*
  (all.parties.lvl2.means$ref.lvl1.sd^2)

all.parties.lvl2.means$ref.pooled.sd<-
  sqrt(sum(all.parties.lvl2.means$ref.num)/
         (sum(all.parties.lvl2.means$n)-
            nrow(all.parties.lvl2.means)))
all.parties.lvl2.means
```

```
##         all.parties.lvl2 env.lvl2.mean env.lvl1.sd ref.lvl2.mean ref.lvl1.sd     n    env.num env.pooled.sd    ref.num
## 1 Anti-immigration party   -0.19378366    1.178983   -0.43503661    1.062387  3916  5441.8552      1.179479  4418.7318
## 2           Did not vote   -0.10458783    1.194766   -0.06325198    1.077707  7474 10667.4436      1.179479  8679.5307
## 3             Don't know   -0.14719577    1.124986   -0.07646043    1.045010  1163  1470.6209      1.179479  1268.9574
## 4                 NE age    0.23978598    1.130640    0.23289439    1.012210  1730  2210.2631      1.179479  1771.4817
## 5             NE citizen    0.03985547    1.206684    0.33587859    1.100438  1150  1673.0423      1.179479  1391.3980
## 6               NE other    0.13291121    1.284980    0.19966753    1.061528   249   409.4912      1.179479   279.4566
## 7            Other party    0.01713584    1.180109    0.03516285    1.030356 17905 24934.1409      1.179479 19007.4739
## 8  Pro-environment party    0.50764949    1.161147    0.40214213    0.982227  2129  2869.0999      1.179479  2053.0304
##   ref.pooled.sd
## 1      1.043338
## 2      1.043338
## 3      1.043338
## 4      1.043338
## 5      1.043338
## 6      1.043338
## 7      1.043338
## 8      1.043338
```

##### group.1: Anti-immigration party


```r
x.points.group.1<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[1]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[1]))

x.points.group.1
```

```
## [1] -0.9773814  0.2020980  1.3815774
```

```r
EX9RR.NCmod10.slopes.group.1<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.1)))
EX9RR.NCmod10.slopes.group.1.tab<-data.frame(EX9RR.NCmod10.slopes.group.1)
(EX9RR.NCmod10.slopes.group.1.tab<-
    EX9RR.NCmod10.slopes.group.1.tab[EX9RR.NCmod10.slopes.group.1.tab$all.parties.lvl2=="Anti-immigration party",])
```

```
##   environ.lvl1       all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 1   -0.9773814 Anti-immigration party         0.23154581 0.04109604 Inf  0.15099905 0.3120926
## 2    0.2020980 Anti-immigration party         0.03418445 0.03241705 Inf -0.02935180 0.0977207
## 3    1.3815774 Anti-immigration party         0.04525327 0.04640982 Inf -0.04570832 0.1362149
```

```r
EX9RR.NCmod10.slopes.group.1.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.1.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.1.tab$SE)))
EX9RR.NCmod10.slopes.group.1.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.1.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.1.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.1.tab[,1],group=EX9RR.NCmod10.slopes.group.1.tab[,2],
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.1.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.1.tab
```

```
##    env_point                  group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1 -0.9773814 Anti-immigration party               0.23 0.04 0.0000 0.0000      0.15      0.31
## 2  0.2020980 Anti-immigration party               0.03 0.03 0.2916 0.5833     -0.03      0.10
## 3  1.3815774 Anti-immigration party               0.05 0.05 0.3295 0.5833     -0.05      0.14
```

\newpage

##### group.2: Did not vote -group


```r
x.points.group.2<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[2]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[2]))

x.points.group.2
```

```
## [1] -1.0665772  0.1129022  1.2923816
```

```r
EX9RR.NCmod10.slopes.group.2<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.2)))
EX9RR.NCmod10.slopes.group.2.tab<-data.frame(EX9RR.NCmod10.slopes.group.2)
(EX9RR.NCmod10.slopes.group.2.tab<-
    EX9RR.NCmod10.slopes.group.2.tab[EX9RR.NCmod10.slopes.group.2.tab$all.parties.lvl2=="Did not vote",])
```

```
##   environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 4   -1.0665772     Did not vote         0.16732977 0.03497184 Inf 0.098786211 0.2358733
## 5    0.1129022     Did not vote         0.09333601 0.02453223 Inf 0.045253722 0.1414183
## 6    1.2923816     Did not vote         0.07381838 0.03499545 Inf 0.005228565 0.1424082
```

```r
EX9RR.NCmod10.slopes.group.2.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.2.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.2.tab$SE)))
EX9RR.NCmod10.slopes.group.2.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.2.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.2.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.2.tab[,1],group=EX9RR.NCmod10.slopes.group.2.tab[,2],
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.2.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.2.tab
```

```
##    env_point        group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 4 -1.0665772 Did not vote               0.17 0.03 0.0000 0.0000      0.10      0.24
## 5  0.1129022 Did not vote               0.09 0.02 0.0001 0.0003      0.05      0.14
## 6  1.2923816 Did not vote               0.07 0.03 0.0349 0.0349      0.01      0.14
```

\newpage

##### group.3: Don't know -group


```r
x.points.group.3<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[3]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[3]))

x.points.group.3
```

```
## [1] -1.0239693  0.1555101  1.3349895
```

```r
EX9RR.NCmod10.slopes.group.3<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.3)))
EX9RR.NCmod10.slopes.group.3.tab<-data.frame(EX9RR.NCmod10.slopes.group.3)
(EX9RR.NCmod10.slopes.group.3.tab<-
    EX9RR.NCmod10.slopes.group.3.tab[EX9RR.NCmod10.slopes.group.3.tab$all.parties.lvl2=="Don't know",])
```

```
##   environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 7   -1.0239693       Don't know          0.1129927 0.06771479 Inf -0.01972586 0.2457112
## 8    0.1555101       Don't know          0.1228098 0.05497948 Inf  0.01505204 0.2305676
## 9    1.3349895       Don't know          0.1924518 0.07371813 Inf  0.04796695 0.3369367
```

```r
EX9RR.NCmod10.slopes.group.3.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.3.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.3.tab$SE)))
EX9RR.NCmod10.slopes.group.3.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.3.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.3.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.3.tab[,1],group=EX9RR.NCmod10.slopes.group.3.tab[,2],
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.3.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.3.tab
```

```
##    env_point      group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 7 -1.0239693 Don't know               0.11 0.07 0.0952 0.0952     -0.02      0.25
## 8  0.1555101 Don't know               0.12 0.05 0.0255 0.0510      0.02      0.23
## 9  1.3349895 Don't know               0.19 0.07 0.0090 0.0271      0.05      0.34
```

\newpage

##### group.4: Not eligible Age -group


```r
x.points.group.4<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[4]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[4]))

x.points.group.4
```

```
## [1] -1.4109510 -0.2314716  0.9480078
```

```r
EX9RR.NCmod10.slopes.group.4<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.4)))
EX9RR.NCmod10.slopes.group.4.tab<-data.frame(EX9RR.NCmod10.slopes.group.4)
(EX9RR.NCmod10.slopes.group.4.tab<-
    EX9RR.NCmod10.slopes.group.4.tab[EX9RR.NCmod10.slopes.group.4.tab$all.parties.lvl2=="NE age",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 10   -1.4109510           NE age         0.19266290 0.07209069 Inf  0.05136775 0.3339581
## 11   -0.2314716           NE age         0.04264037 0.04503743 Inf -0.04563137 0.1309121
## 12    0.9480078           NE age         0.20513251 0.04785432 Inf  0.11133976 0.2989253
```

```r
EX9RR.NCmod10.slopes.group.4.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.4.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.4.tab$SE)))
EX9RR.NCmod10.slopes.group.4.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.4.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.4.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.4.tab[,1],group=EX9RR.NCmod10.slopes.group.4.tab[,2],
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.4.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.4.tab
```

```
##     env_point  group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 10 -1.4109510 NE age               0.19 0.07 0.0075 0.0151      0.05      0.33
## 11 -0.2314716 NE age               0.04 0.05 0.3438 0.3438     -0.05      0.13
## 12  0.9480078 NE age               0.21 0.05 0.0000 0.0001      0.11      0.30
```

\newpage

##### group.5: Not eligible citizenship -group


```r
x.points.group.5<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[5]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[5]))

x.points.group.5
```

```
## [1] -1.21102051 -0.03154112  1.14793828
```

```r
EX9RR.NCmod10.slopes.group.5<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.5)))
EX9RR.NCmod10.slopes.group.5.tab<-data.frame(EX9RR.NCmod10.slopes.group.5)
(EX9RR.NCmod10.slopes.group.5.tab<-
    EX9RR.NCmod10.slopes.group.5.tab[EX9RR.NCmod10.slopes.group.5.tab$all.parties.lvl2=="NE citizen",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df    asymp.LCL asymp.UCL
## 13  -1.21102051       NE citizen         0.04100145 0.07103284 Inf -0.098220356 0.1802233
## 14  -0.03154112       NE citizen         0.05474675 0.05419441 Inf -0.051472351 0.1609658
## 15   1.14793828       NE citizen         0.11426897 0.05918868 Inf -0.001738718 0.2302767
```

```r
EX9RR.NCmod10.slopes.group.5.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.5.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.5.tab$SE)))
EX9RR.NCmod10.slopes.group.5.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.5.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.5.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.5.tab[,1],group=EX9RR.NCmod10.slopes.group.5.tab[,2],
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.5.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.5.tab
```

```
##      env_point      group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 13 -1.21102051 NE citizen               0.04 0.07 0.5638 0.6248     -0.10      0.18
## 14 -0.03154112 NE citizen               0.05 0.05 0.3124 0.6248     -0.05      0.16
## 15  1.14793828 NE citizen               0.11 0.06 0.0535 0.1606      0.00      0.23
```

\newpage

##### group.6: Not eligible Other -group


```r
x.points.group.6<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[6]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[6]))

x.points.group.6
```

```
## [1] -1.3040762 -0.1245969  1.0548825
```

```r
EX9RR.NCmod10.slopes.group.6<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.6)))
EX9RR.NCmod10.slopes.group.6.tab<-data.frame(EX9RR.NCmod10.slopes.group.6)
(EX9RR.NCmod10.slopes.group.6.tab<-
    EX9RR.NCmod10.slopes.group.6.tab[EX9RR.NCmod10.slopes.group.6.tab$all.parties.lvl2=="NE other",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend        SE  df  asymp.LCL asymp.UCL
## 16   -1.3040762         NE other          0.1507720 0.1751231 Inf -0.1924630 0.4940069
## 17   -0.1245969         NE other          0.0747065 0.1046609 Inf -0.1304251 0.2798381
## 18    1.0548825         NE other          0.1603026 0.1068783 Inf -0.0491750 0.3697803
```

```r
EX9RR.NCmod10.slopes.group.6.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.6.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.6.tab$SE)))
EX9RR.NCmod10.slopes.group.6.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.6.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.6.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.6.tab[,1],group=EX9RR.NCmod10.slopes.group.6.tab[,2],
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.6.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.6.tab
```

```
##     env_point    group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 16 -1.3040762 NE other               0.15 0.18 0.3893 0.7785     -0.19      0.49
## 17 -0.1245969 NE other               0.07 0.10 0.4754 0.7785     -0.13      0.28
## 18  1.0548825 NE other               0.16 0.11 0.1337 0.4010     -0.05      0.37
```

##### group.7: Other party -group


```r
x.points.group.7<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[7]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[7]))

x.points.group.7
```

```
## [1] -1.188300879 -0.008821484  1.170657911
```

```r
EX9RR.NCmod10.slopes.group.7<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.7)))
EX9RR.NCmod10.slopes.group.7.tab<-data.frame(EX9RR.NCmod10.slopes.group.7)
(EX9RR.NCmod10.slopes.group.7.tab<-
    EX9RR.NCmod10.slopes.group.7.tab[EX9RR.NCmod10.slopes.group.7.tab$all.parties.lvl2=="Other party",])
```

```
##    environ.lvl1 all.parties.lvl2 environ.lvl1.trend         SE  df  asymp.LCL asymp.UCL
## 19 -1.188300879      Other party          0.1539362 0.02639077 Inf 0.10221120 0.2056611
## 20 -0.008821484      Other party          0.1034603 0.01615921 Inf 0.07178879 0.1351317
## 21  1.170657911      Other party          0.1277239 0.02423238 Inf 0.08022934 0.1752185
```

```r
EX9RR.NCmod10.slopes.group.7.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.7.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.7.tab$SE)))
EX9RR.NCmod10.slopes.group.7.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.7.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.7.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.7.tab[,1],group=EX9RR.NCmod10.slopes.group.7.tab[,2],
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.7.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.7.tab
```

```
##       env_point       group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 19 -1.188300879 Other party               0.15 0.03 0     0      0.10      0.21
## 20 -0.008821484 Other party               0.10 0.02 0     0      0.07      0.14
## 21  1.170657911 Other party               0.13 0.02 0     0      0.08      0.18
```

\newpage

##### group.8: Pro-environment party -group


```r
x.points.group.8<-
  c(mean(dat$environ.lvl1)-1*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]),
    mean(dat$environ.lvl1)-0*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]),
    mean(dat$environ.lvl1)+1*all.parties.lvl2.means$env.pooled.sd[8]+
      (mean(dat$environ.lvl1)-all.parties.lvl2.means$env.lvl2.mean[8]))

x.points.group.8
```

```
## [1] -1.6788145 -0.4993351  0.6801443
```

```r
EX9RR.NCmod10.slopes.group.8<-emtrends(EX9RR.NCmod10,specs = c("environ.lvl1","all.parties.lvl2"),var="environ.lvl1",
                                       at=list(environ.lvl1=c(x.points.group.8)))
EX9RR.NCmod10.slopes.group.8.tab<-data.frame(EX9RR.NCmod10.slopes.group.8)
(EX9RR.NCmod10.slopes.group.8.tab<-
    EX9RR.NCmod10.slopes.group.8.tab[EX9RR.NCmod10.slopes.group.8.tab$all.parties.lvl2=="Pro-environment party",])
```

```
##    environ.lvl1      all.parties.lvl2 environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 22   -1.6788145 Pro-environment party         0.08573603 0.07155012 Inf -0.05449963 0.2259717
## 23   -0.4993351 Pro-environment party         0.14243002 0.03896286 Inf  0.06606422 0.2187958
## 24    0.6801443 Pro-environment party         0.16659876 0.03920425 Inf  0.08975984 0.2434377
```

```r
EX9RR.NCmod10.slopes.group.8.tab$p<-
  2*(1-pnorm(abs(EX9RR.NCmod10.slopes.group.8.tab$environ.lvl1.trend/
                   EX9RR.NCmod10.slopes.group.8.tab$SE)))
EX9RR.NCmod10.slopes.group.8.tab$adj.p<-
  p.adjust(EX9RR.NCmod10.slopes.group.8.tab$p,method="holm")

EX9RR.NCmod10.slopes.group.8.tab<-
  cbind(env_point=EX9RR.NCmod10.slopes.group.8.tab[,1],group=EX9RR.NCmod10.slopes.group.8.tab[,2],
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(3,4)],2),
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(8,9)],4),
        round(EX9RR.NCmod10.slopes.group.8.tab[,c(6,7)],2))
EX9RR.NCmod10.slopes.group.8.tab
```

```
##     env_point                 group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 22 -1.6788145 Pro-environment party               0.09 0.07 0.2308 0.2308     -0.05      0.23
## 23 -0.4993351 Pro-environment party               0.14 0.04 0.0003 0.0005      0.07      0.22
## 24  0.6801443 Pro-environment party               0.17 0.04 0.0000 0.0001      0.09      0.24
```

### Plotting the effects

#### Make a separate long format data file for the focal groups + young non-voters


```r
EX9.plot.dat<-Rdat %>%
  filter(anti.imm.party.dummy==1 |
           pro.env.party.dummy==1 |
           other.party.dummy==1 |
           not.eligible.age.dummy==1)

#rename the groups
EX9.plot.dat$party<-ifelse(EX9.plot.dat$anti.imm.party.dummy==1,"Anti-immigration voters",
                           ifelse(EX9.plot.dat$pro.env.party.dummy==1,"Pro-environment voters",
                                  ifelse(EX9.plot.dat$other.party.dummy==1,"Other party voters",
                                         ifelse(EX9.plot.dat$not.eligible.age.dummy==1,"Not eligible to vote: Age",NA))))
```

\newpage

#### Plot the curves


```r
non.lin.plot<-
  ggplot(data=EX9.plot.dat,aes(x=environ.gmc,y=refugees))+
  geom_smooth(method="lm",
              se=T,
              formula=y ~ poly(x, 3),
              size=2,color="black")+
  geom_smooth(method="lm",
              se=F,
              formula=y ~ x,
              size=1,
              color="red",linetype="dashed")+
  xlab("Attitudes towards the Environment")+
  ylab("Attitudes towards Refugees")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.25,
                                        linetype = 'dotted',
                                        colour = "black"))+
  coord_cartesian(xlim=c(-2,2))+
  facet_wrap(~party,ncol=2)+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.background = element_blank())

non.lin.plot
```

![](All_analyses_files/figure-html/unnamed-chunk-92-1.png)<!-- -->

```r
ggsave(plot = non.lin.plot,
       filename="non.lin.plot.png",device = "png",
       units = "cm",width=12,height=18,dpi = 600)
```


\newpage



# Exploratory analyses for moderators

## Does the association vary by age?

### Center the age variable


```r
describe(dat$age)
```

```
##    vars     n  mean   sd median trimmed   mad    min   max range skew kurtosis  se
## X1    1 35740 -0.09 18.4   0.64   -0.17 22.24 -34.36 50.64    85 0.02    -0.91 0.1
```

```r
#already grand-mean centered

#obtain dataframe with country means and add to data

age.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(age.cntry=mean(age,na.rm=T))

dat<-left_join(x=dat,
               y=age.cntry,
               by=c("cntry"))

#center individuals around country means

dat$age.cntrymc<-dat$age-dat$age.cntry

#obtain dataframe with voting group means and add to data

age.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(age.voting.group=mean(age.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=age.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$age.vgmc<-dat$age.cntrymc-dat$age.voting.group

#describe the variable

describe(dat$age.vgmc)
```

```
##    vars     n mean    sd median trimmed   mad    min   max range skew kurtosis   se
## X1    1 35740    0 16.12  -0.02   -0.18 17.63 -41.78 55.52  97.3 0.09    -0.56 0.09
```

```r
#rename as lvl1, lvl2, and lvl3

dat$age.lvl1<-dat$age.vgmc
dat$age.lvl2<-dat$age.voting.group
dat$age.lvl3<-dat$age.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered age)

* Divide age variable by 10 to give interpretation by a decade


```r
dat$age.lvl1.10<-dat$age.lvl1/10

EX3.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+educ+resid+
                 age.lvl1.10+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX3.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX3.mod1<-getVC(EX3.mod1))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30440667  0.0926634196
## 2 voting.group environ.lvl1         <NA>  0.04306793  0.0018548466
## 3 voting.group  (Intercept) environ.lvl1  0.11749126  0.0015403298
## 4        cntry  (Intercept)         <NA>  0.49917568  0.2491763637
## 5        cntry environ.lvl1         <NA>  0.04008438  0.0016067573
## 6        cntry  (Intercept) environ.lvl1 -0.04834386 -0.0009673194
## 7     Residual         <NA>         <NA>  1.02911869  1.0590852861
```

```r
getFE(EX3.mod1)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.52  0.51 0.615 -0.21  0.36
## 2                                                   gender  0.06 0.01 35582.65  4.68 0.000  0.03  0.08
## 3                            occupClerical support workers -0.04 0.09 35512.24 -0.47 0.640 -0.21  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35522.70 -0.76 0.446 -0.24  0.11
## 5                              occupElementary occupations  0.01 0.09 35527.09  0.14 0.892 -0.16  0.19
## 6                                            occupManagers -0.02 0.09 35509.88 -0.22 0.824 -0.19  0.15
## 7                             occupOther: Not in paid work  0.14 0.09 35700.60  1.53 0.126 -0.04  0.32
## 8         occupPlant and machine operators, and assemblers -0.06 0.09 35519.71 -0.73 0.464 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35514.53  0.78 0.434 -0.10  0.24
## 10                                            occupRetired -0.04 0.10 35505.90 -0.40 0.687 -0.23  0.15
## 11                          occupService and sales workers -0.04 0.09 35521.63 -0.43 0.671 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35525.96 -0.63 0.530 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35511.65 -0.37 0.710 -0.20  0.14
## 14                                         occupUnemployed -0.02 0.11 35547.75 -0.20 0.840 -0.23  0.19
## 15                                                    educ  0.01 0.00 35715.13  7.76 0.000  0.01  0.02
## 16                                                   resid -0.06 0.01 35639.20 -5.25 0.000 -0.08 -0.04
## 17                                             age.lvl1.10  0.02 0.00 35528.89  5.25 0.000  0.01  0.03
## 18                                            environ.lvl1  0.12 0.01    19.40 11.52 0.000  0.10  0.15
```


\newpage

### Model 2 (interaction between age and environmental attitudes)


```r
EX3.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+educ+resid+
                 age.lvl1.10+
                environ.lvl1+
                 age.lvl1.10:environ.lvl1,
                 data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX3.mod1,EX3.mod2)
```

```
## Data: dat
## Models:
## EX3.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX3.mod1:     gender + occup + educ + resid + age.lvl1.10 + environ.lvl1
## EX3.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX3.mod2:     gender + occup + educ + resid + age.lvl1.10 + environ.lvl1 + 
## EX3.mod2:     age.lvl1.10:environ.lvl1
##          npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX3.mod1   25 104282 104494 -52116   104232                        
## EX3.mod2   26 104275 104495 -52111   104223 9.3761  1   0.002198 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX3.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX3.mod2<-getVC(EX3.mod2))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30426074  0.0925745987
## 2 voting.group environ.lvl1         <NA>  0.04295373  0.0018450227
## 3 voting.group  (Intercept) environ.lvl1  0.12759011  0.0016674922
## 4        cntry  (Intercept)         <NA>  0.49911415  0.2491149305
## 5        cntry environ.lvl1         <NA>  0.04005143  0.0016041172
## 6        cntry  (Intercept) environ.lvl1 -0.04815853 -0.0009627004
## 7     Residual         <NA>         <NA>  1.02899151  1.0588235363
```

```r
getFE(EX3.mod2)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.52  0.51 0.614 -0.21  0.36
## 2                                                   gender  0.05 0.01 35582.46  4.66 0.000  0.03  0.08
## 3                            occupClerical support workers -0.04 0.09 35512.30 -0.48 0.632 -0.22  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35522.74 -0.77 0.443 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35527.15  0.12 0.902 -0.16  0.18
## 6                                            occupManagers -0.02 0.09 35509.94 -0.23 0.815 -0.19  0.15
## 7                             occupOther: Not in paid work  0.14 0.09 35700.95  1.50 0.134 -0.04  0.31
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35519.77 -0.74 0.458 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35514.57  0.76 0.447 -0.10  0.24
## 10                                            occupRetired -0.04 0.10 35505.92 -0.43 0.665 -0.23  0.15
## 11                          occupService and sales workers -0.04 0.09 35521.68 -0.44 0.661 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35526.03 -0.64 0.525 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35511.69 -0.38 0.702 -0.20  0.14
## 14                                         occupUnemployed -0.02 0.11 35547.78 -0.21 0.837 -0.23  0.19
## 15                                                    educ  0.01 0.00 35716.30  7.84 0.000  0.01  0.02
## 16                                                   resid -0.06 0.01 35639.39 -5.24 0.000 -0.08 -0.04
## 17                                             age.lvl1.10  0.02 0.00 35528.59  5.35 0.000  0.01  0.03
## 18                                            environ.lvl1  0.12 0.01    19.40 11.54 0.000  0.10  0.15
## 19                                age.lvl1.10:environ.lvl1 -0.01 0.00 35544.83 -3.06 0.002 -0.01 -0.00
```


\newpage

#### Marginal effects for ages at -1SD and +1SD


```r
EX3.mod2.trends<-
  emtrends(EX3.mod2,specs = c("age.lvl1.10"),var=c("environ.lvl1"),
           at=list(age.lvl1.10=c(
             
             mean(dat$age.lvl1.10)-sd(dat$age.lvl1.10),
             mean(dat$age.lvl1.10),
             mean(dat$age.lvl1.10)+sd(dat$age.lvl1.10)
             )))

(EX3.mod2.trends.tab<-data.frame(EX3.mod2.trends))
```

```
##     age.lvl1.10 environ.lvl1.trend         SE  df  asymp.LCL asymp.UCL
## 1 -1.611745e+00          0.1381195 0.01171672 Inf 0.11515520 0.1610839
## 2  1.780079e-18          0.1237155 0.01072314 Inf 0.10269849 0.1447324
## 3  1.611745e+00          0.1093114 0.01170135 Inf 0.08637713 0.1322456
```

```r
EX3.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX3.mod2.trends.tab$environ.lvl1.trend/
                   EX3.mod2.trends.tab$SE)))
EX3.mod2.trends.tab$adj.p<-
  p.adjust(EX3.mod2.trends.tab$p,method="holm")

EX3.mod2.trends.tab<-
  cbind(group=round(EX3.mod2.trends.tab[,1],2),
      round(EX3.mod2.trends.tab[,c(2,3)],2),
      round(EX3.mod2.trends.tab[,c(7,8)],4),
      round(EX3.mod2.trends.tab[,c(5,6)],2))
EX3.mod2.trends.tab
```

```
##   group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -1.61               0.14 0.01 0     0      0.12      0.16
## 2  0.00               0.12 0.01 0     0      0.10      0.14
## 3  1.61               0.11 0.01 0     0      0.09      0.13
```

```r
pairs(EX3.mod2.trends,adjust="none")
```

```
##  contrast                                 estimate      SE  df z.ratio p.value
##  -1.61174532093322 - 1.78007855902334e-18   0.0144 0.00470 Inf 3.063   0.0022 
##  -1.61174532093322 - 1.61174532093322       0.0288 0.00941 Inf 3.063   0.0022 
##  1.78007855902334e-18 - 1.61174532093322    0.0144 0.00470 Inf 3.063   0.0022 
## 
## Results are averaged over the levels of: gender, occup, resid 
## Degrees-of-freedom method: asymptotic
```

\newpage

## Does the association vary by sex?

### Center the sex variable


```r
describe(dat$gender)
```

```
##    vars     n mean  sd median trimmed mad  min max range  skew kurtosis se
## X1    1 35740 0.02 0.5    0.5    0.02   0 -0.5 0.5     1 -0.07       -2  0
```

```r
#grand mean center
dat$gender.gmc<-dat$gender-mean(dat$gender,na.rm=T)

#obtain dataframe with country means and add to data

gender.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(gender.cntry=mean(gender.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=gender.cntry,
               by=c("cntry"))

#center individuals around country means

dat$gender.cntrymc<-dat$gender.gmc-dat$gender.cntry

#obtain dataframe with voting group means and add to data

gender.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(gender.voting.group=mean(gender.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=gender.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$gender.vgmc<-dat$gender.cntrymc-dat$gender.voting.group

#describe the variable

describe(dat$gender.vgmc)
```

```
##    vars     n mean   sd median trimmed mad  min  max range  skew kurtosis se
## X1    1 35740    0 0.49   0.34       0 0.4 -0.8 0.86  1.66 -0.07    -1.91  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$gender.lvl1<-dat$gender.vgmc
dat$gender.lvl2<-dat$gender.voting.group
dat$gender.lvl3<-dat$gender.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered sex)


```r
EX1.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+occup+educ+resid+
                 gender.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX1.mod1)
```

```
## [1] FALSE
```

```r
getVC(EX1.mod1)
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31081495  0.096605935
## 2 voting.group environ.lvl1         <NA>  0.04319467  0.001865779
## 3 voting.group  (Intercept) environ.lvl1  0.13043033  0.001751099
## 4        cntry  (Intercept)         <NA>  0.49963845  0.249638579
## 5        cntry environ.lvl1         <NA>  0.03998946  0.001599157
## 6        cntry  (Intercept) environ.lvl1 -0.05059339 -0.001010870
## 7     Residual         <NA>         <NA>  1.02912353  1.059095237
```

```r
getFE(EX1.mod1)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.42  0.53 0.597 -0.21  0.36
## 2                                                      age  0.00 0.00 34236.53  4.50 0.000  0.00  0.00
## 3                            occupClerical support workers -0.04 0.09 35508.20 -0.46 0.646 -0.21  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35517.67 -0.77 0.439 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35520.15  0.14 0.888 -0.16  0.19
## 6                                            occupManagers -0.02 0.09 35506.79 -0.21 0.833 -0.19  0.15
## 7                             occupOther: Not in paid work  0.14 0.09 35673.27  1.54 0.124 -0.04  0.32
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35515.40 -0.74 0.462 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35512.10  0.80 0.423 -0.10  0.24
## 10                                            occupRetired -0.03 0.10 35508.45 -0.35 0.727 -0.23  0.16
## 11                          occupService and sales workers -0.04 0.09 35513.98 -0.42 0.671 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35522.40 -0.62 0.532 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35507.83 -0.37 0.714 -0.20  0.14
## 14                                         occupUnemployed -0.02 0.11 35524.96 -0.23 0.819 -0.23  0.19
## 15                                                    educ  0.01 0.00 35701.13  7.49 0.000  0.01  0.02
## 16                                                   resid -0.06 0.01 35633.28 -5.23 0.000 -0.08 -0.04
## 17                                             gender.lvl1  0.05 0.01 35448.50  4.44 0.000  0.03  0.08
## 18                                            environ.lvl1  0.12 0.01    19.40 11.53 0.000  0.10  0.15
```

\newpage

### Model 2 (interaction between sex and environmental attitudes)


```r
EX1.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+occup+educ+resid+
                 gender.lvl1+
                environ.lvl1+
                 gender.lvl1:environ.lvl1,
                 data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX1.mod1,EX1.mod2)
```

```
## Data: dat
## Models:
## EX1.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX1.mod1:     age + occup + educ + resid + gender.lvl1 + environ.lvl1
## EX1.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX1.mod2:     age + occup + educ + resid + gender.lvl1 + environ.lvl1 + 
## EX1.mod2:     gender.lvl1:environ.lvl1
##          npar    AIC    BIC logLik deviance Chisq Df Pr(>Chisq)
## EX1.mod1   25 104292 104504 -52121   104242                    
## EX1.mod2   26 104294 104514 -52121   104242 0.012  1     0.9128
```

```r
isSingular(EX1.mod2)
```

```
## [1] FALSE
```

```r
getVC(EX1.mod2)
```

```
##            grp         var1         var2      est_SD      est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31082638  0.096613040
## 2 voting.group environ.lvl1         <NA>  0.04319942  0.001866190
## 3 voting.group  (Intercept) environ.lvl1  0.13037807  0.001750654
## 4        cntry  (Intercept)         <NA>  0.49964041  0.249640543
## 5        cntry environ.lvl1         <NA>  0.03999356  0.001599484
## 6        cntry  (Intercept) environ.lvl1 -0.05050727 -0.001009256
## 7     Residual         <NA>         <NA>  1.02912290  1.059093940
```

```r
getFE(EX1.mod2)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.42  0.53 0.598 -0.21  0.36
## 2                                                      age  0.00 0.00 34235.66  4.50 0.000  0.00  0.00
## 3                            occupClerical support workers -0.04 0.09 35508.12 -0.46 0.646 -0.21  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35517.72 -0.77 0.439 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35520.04  0.14 0.888 -0.16  0.19
## 6                                            occupManagers -0.02 0.09 35506.74 -0.21 0.833 -0.19  0.16
## 7                             occupOther: Not in paid work  0.14 0.09 35673.24  1.54 0.124 -0.04  0.32
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35515.40 -0.74 0.462 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35512.07  0.80 0.423 -0.10  0.24
## 10                                            occupRetired -0.03 0.10 35508.36 -0.35 0.728 -0.23  0.16
## 11                          occupService and sales workers -0.04 0.09 35513.90 -0.42 0.671 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35522.38 -0.62 0.532 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35507.78 -0.37 0.714 -0.20  0.14
## 14                                         occupUnemployed -0.02 0.11 35524.83 -0.23 0.819 -0.23  0.19
## 15                                                    educ  0.01 0.00 35700.98  7.48 0.000  0.01  0.02
## 16                                                   resid -0.06 0.01 35633.12 -5.23 0.000 -0.08 -0.04
## 17                                             gender.lvl1  0.05 0.01 35448.16  4.44 0.000  0.03  0.08
## 18                                            environ.lvl1  0.12 0.01    19.42 11.53 0.000  0.10  0.15
## 19                                gender.lvl1:environ.lvl1  0.00 0.01 35516.35  0.11 0.913 -0.02  0.02
```

\newpage


## Does the association vary by education (years)?

### Center the education variable


```r
describe(dat$educ)
```

```
##    vars     n mean   sd median trimmed  mad    min   max range skew kurtosis   se
## X1    1 35740 0.04 3.86  -0.02   -0.01 2.97 -13.02 40.98    54 0.33     1.88 0.02
```

```r
#already grand-mean centered

#obtain dataframe with country means and add to data

educ.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(educ.cntry=mean(educ,na.rm=T))

dat<-left_join(x=dat,
               y=educ.cntry,
               by=c("cntry"))

#center individuals around country means

dat$educ.cntrymc<-dat$educ-dat$educ.cntry

#obtain dataframe with voting group means and add to data

educ.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(educ.voting.group=mean(educ.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=educ.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$educ.vgmc<-dat$educ.cntrymc-dat$educ.voting.group

#describe the variable

describe(dat$educ.vgmc)
```

```
##    vars     n mean   sd median trimmed  mad    min   max range skew kurtosis   se
## X1    1 35740    0 3.56  -0.16   -0.08 3.23 -15.82 39.72 55.54 0.42     2.49 0.02
```

```r
#rename as lvl1, lvl2, and lvl3

dat$educ.lvl1<-dat$educ.vgmc
dat$educ.lvl2<-dat$educ.voting.group
dat$educ.lvl3<-dat$educ.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered education)


```r
EX4.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+age+resid+
                 educ.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX4.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX4.mod1<-getVC(EX4.mod1))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31538141  0.0994654350
## 2 voting.group environ.lvl1         <NA>  0.04319669  0.0018659537
## 3 voting.group  (Intercept) environ.lvl1  0.13786860  0.0018782434
## 4        cntry  (Intercept)         <NA>  0.49968068  0.2496807852
## 5        cntry environ.lvl1         <NA>  0.03998612  0.0015988901
## 6        cntry  (Intercept) environ.lvl1 -0.03730454 -0.0007453556
## 7     Residual         <NA>         <NA>  1.02911789  1.0590836279
```

```r
getFE(EX4.mod1)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.37  0.56 0.579 -0.21  0.37
## 2                                                   gender  0.06 0.01 35565.32  4.69 0.000  0.03  0.08
## 3                            occupClerical support workers -0.04 0.09 35505.86 -0.49 0.622 -0.22  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35518.14 -0.81 0.419 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35523.49  0.08 0.938 -0.17  0.18
## 6                                            occupManagers -0.02 0.09 35505.18 -0.21 0.832 -0.19  0.15
## 7                             occupOther: Not in paid work  0.13 0.09 35680.92  1.46 0.145 -0.05  0.31
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35516.64 -0.78 0.437 -0.24  0.10
## 9                                       occupProfessionals  0.07 0.09 35512.25  0.81 0.417 -0.10  0.24
## 10                                            occupRetired -0.04 0.10 35508.76 -0.39 0.694 -0.23  0.15
## 11                          occupService and sales workers -0.04 0.09 35513.41 -0.47 0.636 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35523.19 -0.67 0.505 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35505.55 -0.38 0.701 -0.21  0.14
## 14                                         occupUnemployed -0.03 0.11 35526.98 -0.27 0.785 -0.24  0.18
## 15                                                     age  0.00 0.00 34776.42  4.34 0.000  0.00  0.00
## 16                                                   resid -0.06 0.01 35638.63 -5.28 0.000 -0.09 -0.04
## 17                                               educ.lvl1  0.01 0.00 35615.19  6.90 0.000  0.01  0.02
## 18                                            environ.lvl1  0.12 0.01    19.39 11.56 0.000  0.10  0.15
```

\newpage

### Model 2 (interaction between education and environment attitudes)


```r
EX4.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+age+resid+
                 educ.lvl1+
                environ.lvl1+
                 environ.lvl1:educ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

anova(EX4.mod1,EX4.mod2)
```

```
## Data: dat
## Models:
## EX4.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX4.mod1:     gender + occup + age + resid + educ.lvl1 + environ.lvl1
## EX4.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX4.mod2:     gender + occup + age + resid + educ.lvl1 + environ.lvl1 + 
## EX4.mod2:     environ.lvl1:educ.lvl1
##          npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## EX4.mod1   25 104298 104510 -52124   104248                         
## EX4.mod2   26 104280 104501 -52114   104228 19.406  1  1.057e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX4.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX4.mod2<-getVC(EX4.mod2))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31474748  0.0990659785
## 2 voting.group environ.lvl1         <NA>  0.04281063  0.0018327499
## 3 voting.group  (Intercept) environ.lvl1  0.14510199  0.0019551822
## 4        cntry  (Intercept)         <NA>  0.49874182  0.2487434013
## 5        cntry environ.lvl1         <NA>  0.04021941  0.0016176010
## 6        cntry  (Intercept) environ.lvl1 -0.03869034 -0.0007760934
## 7     Residual         <NA>         <NA>  1.02886241  1.0585578524
```

```r
getFE(EX4.mod2)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.48  0.56 0.577 -0.21  0.37
## 2                                                   gender  0.05 0.01 35565.51  4.65 0.000  0.03  0.08
## 3                            occupClerical support workers -0.05 0.09 35506.10 -0.52 0.604 -0.22  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35518.38 -0.85 0.397 -0.25  0.10
## 5                              occupElementary occupations  0.00 0.09 35523.75  0.03 0.976 -0.17  0.18
## 6                                            occupManagers -0.02 0.09 35505.43 -0.25 0.804 -0.20  0.15
## 7                             occupOther: Not in paid work  0.13 0.09 35681.55  1.41 0.158 -0.05  0.31
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35516.92 -0.81 0.418 -0.25  0.10
## 9                                       occupProfessionals  0.07 0.09 35512.31  0.76 0.449 -0.11  0.24
## 10                                            occupRetired -0.04 0.10 35509.14 -0.45 0.651 -0.24  0.15
## 11                          occupService and sales workers -0.04 0.09 35513.61 -0.50 0.619 -0.22  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.07 0.09 35523.45 -0.71 0.479 -0.25  0.12
## 13            occupTechnicians and associate professionals -0.04 0.09 35505.75 -0.40 0.686 -0.21  0.14
## 14                                         occupUnemployed -0.03 0.11 35527.26 -0.31 0.760 -0.24  0.18
## 15                                                     age  0.00 0.00 34768.80  4.22 0.000  0.00  0.00
## 16                                                   resid -0.06 0.01 35639.30 -5.23 0.000 -0.08 -0.04
## 17                                               educ.lvl1  0.01 0.00 35614.80  6.75 0.000  0.01  0.02
## 18                                            environ.lvl1  0.12 0.01    19.40 11.49 0.000  0.10  0.15
## 19                                  educ.lvl1:environ.lvl1  0.01 0.00 35550.14  4.41 0.000  0.00  0.01
```


\newpage

#### Marginal effects for different levels of education


```r
EX4.mod2.trends<-
  emtrends(EX4.mod2,specs = c("educ.lvl1"),var=c("environ.lvl1"),
           at=
             list(educ.lvl1=
                    c(mean(dat$educ.lvl1)-sd(dat$educ.lvl1),
                      mean(dat$educ.lvl1),                              mean(dat$educ.lvl1)+sd(dat$educ.lvl1))))
(EX4.mod2.trends.tab<-data.frame(EX4.mod2.trends))
```

```
##       educ.lvl1 environ.lvl1.trend         SE  df  asymp.LCL asymp.UCL
## 1 -3.561414e+00          0.1033404 0.01171741 Inf 0.08037474 0.1263061
## 2  9.056577e-19          0.1234926 0.01075061 Inf 0.10242182 0.1445634
## 3  3.561414e+00          0.1436448 0.01164837 Inf 0.12081444 0.1664752
```

```r
EX4.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX4.mod2.trends.tab$environ.lvl1.trend/
                   EX4.mod2.trends.tab$SE)))
EX4.mod2.trends.tab$adj.p<-
  p.adjust(EX4.mod2.trends.tab$p,method="holm")

EX4.mod2.trends.tab<-
  cbind(group=round(EX4.mod2.trends.tab[,1],2),
        round(EX4.mod2.trends.tab[,c(2,3)],2),
        round(EX4.mod2.trends.tab[,c(7,8)],4),
        round(EX4.mod2.trends.tab[,c(5,6)],2))
EX4.mod2.trends.tab
```

```
##   group environ.lvl1.trend   SE p adj.p asymp.LCL asymp.UCL
## 1 -3.56               0.10 0.01 0     0      0.08      0.13
## 2  0.00               0.12 0.01 0     0      0.10      0.14
## 3  3.56               0.14 0.01 0     0      0.12      0.17
```

\newpage

## Does the association vary by place of residence (urban/rural)?

### Center the residence variable


```r
describe(dat$resid)
```

```
##    vars     n  mean   sd median trimmed mad  min max range skew kurtosis se
## X1    1 35740 -0.11 0.49   -0.5   -0.14   0 -0.5 0.5     1 0.45    -1.79  0
```

```r
dat$resid.gmc<-dat$resid-mean(dat$resid,na.rm=T)


#obtain dataframe with country means and add to data

resid.cntry<-dat %>%
  group_by(cntry) %>%
  summarize(resid.cntry=mean(resid.gmc,na.rm=T))

dat<-left_join(x=dat,
               y=resid.cntry,
               by=c("cntry"))

#center individuals around country means

dat$resid.cntrymc<-dat$resid.gmc-dat$resid.cntry

#obtain dataframe with voting group means and add to data

resid.voting.group<-dat %>%
  group_by(voting.group) %>%
  summarize(resid.voting.group=mean(resid.cntrymc,na.rm=T))

dat<-left_join(x=dat,
               y=resid.voting.group,
               by=c("voting.group"))

#center individuals around voting group means

dat$resid.vgmc<-dat$resid.cntrymc-dat$resid.voting.group

#describe the variable

describe(dat$resid.vgmc)
```

```
##    vars     n mean   sd median trimmed  mad   min  max range skew kurtosis se
## X1    1 35740    0 0.47  -0.25   -0.02 0.34 -0.86 0.97  1.83 0.41    -1.53  0
```

```r
#rename as lvl1, lvl2, and lvl3

dat$resid.lvl1<-dat$resid.vgmc
dat$resid.lvl2<-dat$resid.voting.group
dat$resid.lvl3<-dat$resid.cntry
```

\newpage

### Model 1 (Same as H1 selected model but with centered residence)


```r
EX5.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+age+educ+
                 resid.lvl1+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX5.mod1)
```

```
## [1] FALSE
```

```r
(VC.EX5.mod1<-getVC(EX5.mod1))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31213323  0.0974271515
## 2 voting.group environ.lvl1         <NA>  0.04314186  0.0018612198
## 3 voting.group  (Intercept) environ.lvl1  0.12768151  0.0017193602
## 4        cntry  (Intercept)         <NA>  0.49932897  0.2493294180
## 5        cntry environ.lvl1         <NA>  0.04001042  0.0016008337
## 6        cntry  (Intercept) environ.lvl1 -0.04882110 -0.0009753655
## 7     Residual         <NA>         <NA>  1.02912226  1.0590926168
```

```r
getFE(EX5.mod1)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.44  0.58 0.563 -0.20  0.37
## 2                                                   gender  0.06 0.01 35568.31  4.69 0.000  0.03  0.08
## 3                            occupClerical support workers -0.04 0.09 35507.14 -0.48 0.633 -0.22  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35517.33 -0.78 0.438 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35519.18  0.12 0.902 -0.16  0.18
## 6                                            occupManagers -0.02 0.09 35506.37 -0.22 0.827 -0.19  0.15
## 7                             occupOther: Not in paid work  0.14 0.09 35672.96  1.52 0.128 -0.04  0.32
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35515.02 -0.74 0.458 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35510.98  0.79 0.432 -0.10  0.24
## 10                                            occupRetired -0.04 0.10 35507.41 -0.37 0.715 -0.23  0.16
## 11                          occupService and sales workers -0.04 0.09 35512.89 -0.45 0.656 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35524.34 -0.65 0.516 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35507.07 -0.38 0.704 -0.20  0.14
## 14                                         occupUnemployed -0.03 0.11 35524.45 -0.24 0.811 -0.24  0.18
## 15                                                     age  0.00 0.00 34269.44  4.48 0.000  0.00  0.00
## 16                                                    educ  0.01 0.00 35698.20  7.52 0.000  0.01  0.02
## 17                                              resid.lvl1 -0.06 0.01 35436.62 -4.83 0.000 -0.08 -0.03
## 18                                            environ.lvl1  0.12 0.01    19.40 11.54 0.000  0.10  0.15
```

\newpage

### Model 2 (interaction between residence and environment attitudes)


```r
EX5.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                gender+occup+age+educ+
                 resid.lvl1+
                environ.lvl1+
                 environ.lvl1:resid.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))

anova(EX5.mod1,EX5.mod2)
```

```
## Data: dat
## Models:
## EX5.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX5.mod1:     gender + occup + age + educ + resid.lvl1 + environ.lvl1
## EX5.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX5.mod2:     gender + occup + age + educ + resid.lvl1 + environ.lvl1 + 
## EX5.mod2:     environ.lvl1:resid.lvl1
##          npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## EX5.mod1   25 104293 104506 -52122   104243                     
## EX5.mod2   26 104295 104516 -52122   104243 0.1275  1      0.721
```

```r
isSingular(EX5.mod2)
```

```
## [1] FALSE
```

```r
(VC.EX5.mod2<-getVC(EX5.mod2))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.31212906  0.0974245472
## 2 voting.group environ.lvl1         <NA>  0.04311736  0.0018591063
## 3 voting.group  (Intercept) environ.lvl1  0.12866809  0.0017316382
## 4        cntry  (Intercept)         <NA>  0.49932953  0.2493299812
## 5        cntry environ.lvl1         <NA>  0.04002016  0.0016016135
## 6        cntry  (Intercept) environ.lvl1 -0.04828142 -0.0009648196
## 7     Residual         <NA>         <NA>  1.02912138  1.0590908084
```

```r
getFE(EX5.mod2)
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.08 0.14    49.44  0.58 0.563 -0.20  0.37
## 2                                                   gender  0.06 0.01 35568.19  4.70 0.000  0.03  0.08
## 3                            occupClerical support workers -0.04 0.09 35507.17 -0.48 0.633 -0.22  0.13
## 4                    occupCraft and related trades workers -0.07 0.09 35517.36 -0.78 0.437 -0.24  0.10
## 5                              occupElementary occupations  0.01 0.09 35519.21  0.12 0.902 -0.16  0.18
## 6                                            occupManagers -0.02 0.09 35506.44 -0.22 0.827 -0.19  0.15
## 7                             occupOther: Not in paid work  0.14 0.09 35673.02  1.52 0.128 -0.04  0.32
## 8         occupPlant and machine operators, and assemblers -0.07 0.09 35515.08 -0.74 0.458 -0.24  0.11
## 9                                       occupProfessionals  0.07 0.09 35511.04  0.79 0.432 -0.10  0.24
## 10                                            occupRetired -0.04 0.10 35507.69 -0.37 0.713 -0.23  0.16
## 11                          occupService and sales workers -0.04 0.09 35512.89 -0.44 0.656 -0.21  0.13
## 12 occupSkilled agricultural, forestry and fishery workers -0.06 0.09 35524.50 -0.65 0.515 -0.24  0.12
## 13            occupTechnicians and associate professionals -0.03 0.09 35507.07 -0.38 0.704 -0.20  0.14
## 14                                         occupUnemployed -0.03 0.11 35524.48 -0.24 0.810 -0.24  0.18
## 15                                                     age  0.00 0.00 34269.03  4.48 0.000  0.00  0.00
## 16                                                    educ  0.01 0.00 35698.28  7.51 0.000  0.01  0.02
## 17                                              resid.lvl1 -0.06 0.01 35436.44 -4.84 0.000 -0.08 -0.03
## 18                                            environ.lvl1  0.12 0.01    19.40 11.53 0.000  0.10  0.15
## 19                                 resid.lvl1:environ.lvl1 -0.00 0.01 35520.29 -0.36 0.721 -0.02  0.02
```


\newpage

## Does the association vary by occupational groups

### Model 1 (Same as H1 selected model)


```r
EX2.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))



isSingular(EX2.mod1)
```

```
## [1] FALSE
```

```r
getVC(EX2.mod1)
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

\newpage


### Model 2 (Interaction between environment and occupational groups)


```r
EX2.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                 occup:environ.lvl1,data=dat,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


anova(EX2.mod1,EX2.mod2)
```

```
## Data: dat
## Models:
## EX2.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX2.mod1:     age + gender + educ + resid + occup + environ.lvl1
## EX2.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## EX2.mod2:     age + gender + educ + resid + occup + environ.lvl1 + occup:environ.lvl1
##          npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
## EX2.mod1   25 104289 104502 -52120   104239                        
## EX2.mod2   37 104281 104594 -52103   104207 32.836 12   0.001027 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
isSingular(EX2.mod2)
```

```
## [1] FALSE
```

```r
getVC(EX2.mod2)
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30949953  0.0957899604
## 2 voting.group environ.lvl1         <NA>  0.04210546  0.0017728693
## 3 voting.group  (Intercept) environ.lvl1  0.09193784  0.0011980989
## 4        cntry  (Intercept)         <NA>  0.49892466  0.2489258209
## 5        cntry environ.lvl1         <NA>  0.03929672  0.0015442324
## 6        cntry  (Intercept) environ.lvl1 -0.03786228 -0.0007423317
## 7     Residual         <NA>         <NA>  1.02867966  1.0581818352
```

\newpage

#### Marginal effects for each occupation group


```r
EX2.mod2.trends<-emtrends(EX2.mod2,specs = c("occup"),var=c("environ.lvl1"))
(EX2.mod2.trends.tab<-data.frame(EX2.mod2.trends))
```

```
##                                                 occup environ.lvl1.trend         SE  df   asymp.LCL asymp.UCL
## 1                                        Armed forces        0.002458622 0.07237495 Inf -0.13939367 0.1443109
## 2                            Clerical support workers        0.171992561 0.01883320 Inf  0.13508017 0.2089050
## 3                    Craft and related trades workers        0.119880759 0.01694750 Inf  0.08666427 0.1530972
## 4                              Elementary occupations        0.101679221 0.01877408 Inf  0.06488271 0.1384757
## 5                                            Managers        0.081133954 0.01964255 Inf  0.04263527 0.1196326
## 6                             Other: Not in paid work        0.159899211 0.02207025 Inf  0.11664232 0.2031561
## 7         Plant and machine operators, and assemblers        0.102839704 0.01985327 Inf  0.06392801 0.1417514
## 8                                       Professionals        0.140211289 0.01480285 Inf  0.11119823 0.1692243
## 9                                             Retired        0.131422064 0.03913153 Inf  0.05472567 0.2081185
## 10                          Service and sales workers        0.110732091 0.01505900 Inf  0.08121700 0.1402472
## 11 Skilled agricultural, forestry and fishery workers        0.121043544 0.03187536 Inf  0.05856899 0.1835181
## 12            Technicians and associate professionals        0.130837107 0.01590712 Inf  0.09965973 0.1620145
## 13                                         Unemployed        0.007406265 0.05474532 Inf -0.09989260 0.1147051
```

```r
EX2.mod2.trends.tab$p<-
  2*(1-pnorm(abs(EX2.mod2.trends.tab$environ.lvl1.trend/
                   EX2.mod2.trends.tab$SE)))
EX2.mod2.trends.tab$adj.p<-
  p.adjust(EX2.mod2.trends.tab$p,method="holm")

EX2.mod2.trends.tab<-
  cbind(group=EX2.mod2.trends.tab[,1],
        round(EX2.mod2.trends.tab[,c(2,3)],2),
        round(EX2.mod2.trends.tab[,c(7,8)],4),
        round(EX2.mod2.trends.tab[,c(5,6)],2))
EX2.mod2.trends.tab
```

```
##                                                 group environ.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1                                        Armed forces               0.00 0.07 0.9729 1.0000     -0.14      0.14
## 2                            Clerical support workers               0.17 0.02 0.0000 0.0000      0.14      0.21
## 3                    Craft and related trades workers               0.12 0.02 0.0000 0.0000      0.09      0.15
## 4                              Elementary occupations               0.10 0.02 0.0000 0.0000      0.06      0.14
## 5                                            Managers               0.08 0.02 0.0000 0.0002      0.04      0.12
## 6                             Other: Not in paid work               0.16 0.02 0.0000 0.0000      0.12      0.20
## 7         Plant and machine operators, and assemblers               0.10 0.02 0.0000 0.0000      0.06      0.14
## 8                                       Professionals               0.14 0.01 0.0000 0.0000      0.11      0.17
## 9                                             Retired               0.13 0.04 0.0008 0.0024      0.05      0.21
## 10                          Service and sales workers               0.11 0.02 0.0000 0.0000      0.08      0.14
## 11 Skilled agricultural, forestry and fishery workers               0.12 0.03 0.0001 0.0006      0.06      0.18
## 12            Technicians and associate professionals               0.13 0.02 0.0000 0.0000      0.10      0.16
## 13                                         Unemployed               0.01 0.05 0.8924 1.0000     -0.10      0.11
```

```r
#contrast for all groups against mean of other groups
contrast(EX2.mod2.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                                                  estimate     SE  df z.ratio p.value
##  Armed forces effect                                       -0.11246 0.0721 Inf -1.559  0.9523 
##  Clerical support workers effect                            0.07120 0.0187 Inf  3.806  0.0018 
##  Craft and related trades workers effect                    0.01474 0.0168 Inf  0.876  1.0000 
##  Elementary occupations effect                             -0.00498 0.0186 Inf -0.267  1.0000 
##  Managers effect                                           -0.02723 0.0196 Inf -1.392  1.0000 
##  Other: Not in paid work effect                             0.05810 0.0220 Inf  2.639  0.0999 
##  Plant and machine operators, and assemblers effect        -0.00372 0.0197 Inf -0.188  1.0000 
##  Professionals effect                                       0.03677 0.0148 Inf  2.491  0.1402 
##  Retired effect                                             0.02725 0.0389 Inf  0.700  1.0000 
##  Service and sales workers effect                           0.00483 0.0149 Inf  0.324  1.0000 
##  Skilled agricultural, forestry and fishery workers effect  0.01600 0.0318 Inf  0.504  1.0000 
##  Technicians and associate professionals effect             0.02661 0.0158 Inf  1.683  0.8310 
##  Unemployed effect                                         -0.10710 0.0545 Inf -1.965  0.4940 
## 
## Results are averaged over the levels of: gender, resid 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 13 tests
```







\newpage



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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.37  0.50 0.616 -0.21  0.36
## 2                                                      age  0.00 0.00 34011.21  3.27 0.001  0.00  0.00
## 3                                                   gender  0.07 0.01 35529.43  5.93 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35634.84  6.12 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35598.12 -4.98 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35472.39 -0.47 0.638 -0.21  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35482.31 -0.69 0.488 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35484.65  0.26 0.795 -0.15  0.20
## 9                                            occupManagers -0.03 0.09 35471.30 -0.29 0.772 -0.20  0.15
## 10                            occupOther: Not in paid work  0.14 0.09 35638.43  1.58 0.114 -0.03  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35480.18 -0.68 0.497 -0.23  0.11
## 12                                      occupProfessionals  0.07 0.09 35475.86  0.75 0.454 -0.11  0.24
## 13                                            occupRetired -0.03 0.10 35473.15 -0.30 0.768 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35478.05 -0.37 0.708 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35487.02 -0.57 0.569 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35471.93 -0.38 0.702 -0.20  0.14
## 17                                         occupUnemployed -0.01 0.11 35488.15 -0.13 0.898 -0.22  0.20
## 18                                            environ.lvl1  0.12 0.01    19.47 11.55 0.000  0.10  0.14
## 19                                            polintr.lvl1  0.06 0.01 35520.18  8.76 0.000  0.05  0.08
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.40  0.49 0.626 -0.22  0.36
## 2                                                      age  0.00 0.00 33996.51  3.24 0.001  0.00  0.00
## 3                                                   gender  0.07 0.01 35529.97  5.95 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35634.15  6.12 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35598.62 -4.95 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35472.67 -0.46 0.645 -0.21  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35482.63 -0.68 0.495 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35484.93  0.27 0.791 -0.15  0.20
## 9                                            occupManagers -0.03 0.09 35471.55 -0.28 0.776 -0.20  0.15
## 10                            occupOther: Not in paid work  0.14 0.09 35638.81  1.59 0.111 -0.03  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35480.48 -0.67 0.505 -0.23  0.11
## 12                                      occupProfessionals  0.07 0.09 35476.12  0.75 0.451 -0.11  0.24
## 13                                            occupRetired -0.03 0.10 35473.40 -0.29 0.770 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35478.36 -0.36 0.717 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35487.38 -0.55 0.579 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35472.25 -0.37 0.713 -0.20  0.14
## 17                                         occupUnemployed -0.01 0.11 35488.44 -0.12 0.908 -0.22  0.20
## 18                                            environ.lvl1  0.12 0.01    19.48 11.51 0.000  0.10  0.14
## 19                                            polintr.lvl1  0.06 0.01 35520.65  8.81 0.000  0.05  0.08
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01 35513.73  2.47 0.013  0.00  0.02
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.39  0.48 0.632 -0.22  0.35
## 2                                                      age  0.00 0.00 33995.67  3.23 0.001  0.00  0.00
## 3                                                   gender  0.07 0.01 35515.56  5.88 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35576.48  5.93 0.000  0.01  0.01
## 5                                                    resid -0.06 0.01 35562.47 -4.97 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35436.79 -0.45 0.650 -0.21  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35444.83 -0.66 0.510 -0.23  0.11
## 8                              occupElementary occupations  0.03 0.09 35444.42  0.28 0.776 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35438.98 -0.24 0.811 -0.19  0.15
## 10                            occupOther: Not in paid work  0.15 0.09 35605.45  1.64 0.101 -0.03  0.33
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35448.35 -0.65 0.512 -0.23  0.12
## 12                                      occupProfessionals  0.07 0.09 35441.87  0.76 0.448 -0.10  0.24
## 13                                            occupRetired -0.03 0.10 35432.14 -0.26 0.794 -0.22  0.17
## 14                          occupService and sales workers -0.03 0.09 35438.92 -0.32 0.749 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35455.55 -0.52 0.600 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35436.62 -0.35 0.725 -0.20  0.14
## 17                                         occupUnemployed -0.01 0.11 35458.28 -0.09 0.926 -0.22  0.20
## 18                                            environ.lvl1  0.12 0.01    19.47 11.41 0.000  0.10  0.14
## 19                                            polintr.lvl1  0.07 0.01    20.77  4.75 0.000  0.04  0.09
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01 35490.32  2.43 0.015  0.00  0.02
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.40  0.48 0.630 -0.22  0.35
## 2                                                      age  0.00 0.00 33970.60  3.25 0.001  0.00  0.00
## 3                                                   gender  0.07 0.01 35497.68  5.88 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35544.12  5.94 0.000  0.01  0.01
## 5                                                    resid -0.06 0.01 35556.74 -4.97 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35421.40 -0.47 0.642 -0.21  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35432.01 -0.67 0.502 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35430.94  0.28 0.782 -0.15  0.20
## 9                                            occupManagers -0.02 0.09 35417.72 -0.24 0.808 -0.19  0.15
## 10                            occupOther: Not in paid work  0.15 0.09 35600.69  1.64 0.102 -0.03  0.33
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35430.70 -0.66 0.506 -0.23  0.11
## 12                                      occupProfessionals  0.06 0.09 35422.24  0.74 0.457 -0.11  0.24
## 13                                            occupRetired -0.02 0.10 35422.55 -0.25 0.803 -0.22  0.17
## 14                          occupService and sales workers -0.03 0.09 35424.90 -0.33 0.741 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35439.06 -0.53 0.595 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35418.90 -0.36 0.717 -0.20  0.14
## 17                                         occupUnemployed -0.01 0.11 35456.40 -0.09 0.926 -0.22  0.20
## 18                                            environ.lvl1  0.12 0.01    19.67 11.49 0.000  0.10  0.14
## 19                                            polintr.lvl1  0.07 0.01    20.83  4.75 0.000  0.04  0.09
## 20                               environ.lvl1:polintr.lvl1  0.01 0.01    56.55  1.97 0.054 -0.00  0.03
```

```r
(VC.H5.exp.intr.mod4<-getVC(H5.exp.intr.mod4))
```

```
##             grp                      var1                      var2       est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.308237181  9.501016e-02
## 2  voting.group              environ.lvl1                      <NA>  0.043365299  1.880549e-03
## 3  voting.group              polintr.lvl1                      <NA>  0.071909474  5.170972e-03
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.036300898  1.317755e-03
## 5  voting.group               (Intercept)              environ.lvl1  0.040243268  5.379236e-04
## 6  voting.group               (Intercept)              polintr.lvl1  0.455278655  1.009133e-02
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1  0.156609524  1.752349e-03
## 8  voting.group              environ.lvl1              polintr.lvl1  0.047988227  1.496453e-04
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.235402468 -3.705704e-04
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1 -0.157191539 -4.103294e-04
## 11        cntry               (Intercept)                      <NA>  0.498544030  2.485461e-01
## 12        cntry              environ.lvl1                      <NA>  0.038211064  1.460085e-03
## 13        cntry              polintr.lvl1                      <NA>  0.046859579  2.195820e-03
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>  0.009102557  8.285654e-05
## 15        cntry               (Intercept)              environ.lvl1 -0.080477026 -1.533079e-03
## 16        cntry               (Intercept)              polintr.lvl1  0.456603483  1.066697e-02
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1  0.037595397  1.706089e-04
## 18        cntry              environ.lvl1              polintr.lvl1  0.039019656  6.986682e-05
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.955389570  3.323021e-04
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.331449878  1.413773e-04
## 21     Residual                      <NA>                      <NA>  1.024901814  1.050424e+00
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
##                            voting.group.(Intercept)               voting.group.environ.lvl1.(Intercept) 
##                                             0.30075                                             0.00170 
##               voting.group.polintr.lvl1.(Intercept)  voting.group.environ.lvl1:polintr.lvl1.(Intercept) 
##                                             0.03194                                             0.00555 
##                           voting.group.environ.lvl1              voting.group.polintr.lvl1.environ.lvl1 
##                                             0.04228                                             0.00208 
## voting.group.environ.lvl1:polintr.lvl1.environ.lvl1                           voting.group.polintr.lvl1 
##                                            -0.00857                                             0.06243 
## voting.group.environ.lvl1:polintr.lvl1.polintr.lvl1              voting.group.environ.lvl1:polintr.lvl1 
##                                            -0.00881                                             0.03275 
##                                   cntry.(Intercept)                      cntry.environ.lvl1.(Intercept) 
##                                             0.48643                                            -0.00300 
##                      cntry.polintr.lvl1.(Intercept)         cntry.environ.lvl1:polintr.lvl1.(Intercept) 
##                                             0.02088                                             0.00033 
##                                  cntry.environ.lvl1                     cntry.polintr.lvl1.environ.lvl1 
##                                             0.03716                                             0.00348 
##        cntry.environ.lvl1:polintr.lvl1.environ.lvl1                                  cntry.polintr.lvl1 
##                                             0.00854                                             0.04053 
##        cntry.environ.lvl1:polintr.lvl1.polintr.lvl1                     cntry.environ.lvl1:polintr.lvl1 
##                                             0.00242                                             0.00000
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    61.98 -3.26 0.002 -0.76 -0.18
## 2                                                      age  0.00 0.00 35616.41  3.63 0.000  0.00  0.00
## 3                                                   gender  0.07 0.01 35528.36  5.93 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35592.31  5.93 0.000  0.01  0.01
## 5                                                    resid -0.06 0.01 35609.54 -4.86 0.000 -0.08 -0.03
## 6                            occupClerical support workers -0.04 0.09 35481.31 -0.48 0.629 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35494.04 -0.68 0.499 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35487.25  0.23 0.818 -0.15  0.19
## 9                                            occupManagers -0.02 0.09 35477.54 -0.27 0.789 -0.20  0.15
## 10                            occupOther: Not in paid work  0.13 0.09 35545.97  1.43 0.153 -0.05  0.31
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35494.29 -0.68 0.495 -0.23  0.11
## 12                                      occupProfessionals  0.06 0.09 35482.88  0.71 0.476 -0.11  0.23
## 13                                            occupRetired -0.03 0.10 35469.24 -0.28 0.783 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35484.58 -0.34 0.731 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35505.47 -0.56 0.573 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35480.34 -0.38 0.706 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35499.21 -0.16 0.872 -0.23  0.19
## 18                                            environ.lvl1  0.11 0.02   109.59  5.58 0.000  0.07  0.15
## 19                                            polintr.lvl1 -0.02 0.03   123.60 -0.85 0.399 -0.08  0.03
## 20                            all.parties.lvl2Did not vote  0.44 0.06   182.55  6.92 0.000  0.31  0.57
## 21                              all.parties.lvl2Don't know  0.42 0.07   267.91  5.96 0.000  0.28  0.56
## 22                            all.parties.lvl2Invalid vote  0.44 0.36  1125.55  1.24 0.214 -0.26  1.15
## 23                                  all.parties.lvl2NE age  0.72 0.07   272.47 10.15 0.000  0.58  0.86
## 24                              all.parties.lvl2NE citizen  0.86 0.08   268.53 11.26 0.000  0.71  1.01
## 25                                all.parties.lvl2NE other  0.71 0.10   653.08  7.32 0.000  0.52  0.91
## 26                               all.parties.lvl2No answer  0.57 0.37  1352.13  1.52 0.128 -0.16  1.31
## 27                             all.parties.lvl2Other party  0.55 0.05   230.27 11.33 0.000  0.45  0.64
## 28                   all.parties.lvl2Pro-environment party  0.92 0.06   256.09 14.11 0.000  0.79  1.04
## 29                               environ.lvl1:polintr.lvl1  0.01 0.01    54.50  1.80 0.077 -0.00  0.03
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    87.38 -0.05 0.956 -0.05  0.04
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   340.06  0.75 0.456 -0.04  0.09
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.05 0.39 24998.32  0.12 0.901 -0.71  0.81
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.87  1.24 0.217 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   236.05 -1.20 0.231 -0.11  0.03
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1206.78  0.58 0.560 -0.08  0.14
## 36                  environ.lvl1:all.parties.lvl2No answer  0.22 0.23  7135.17  0.96 0.335 -0.23  0.67
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   135.62  0.79 0.433 -0.02  0.05
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03   199.37  0.30 0.766 -0.05  0.07
## 39               polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03    95.64  2.47 0.015  0.02  0.14
## 40                 polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05   426.05  2.08 0.038  0.01  0.20
## 41               polintr.lvl1:all.parties.lvl2Invalid vote  0.20 0.35 10239.75  0.58 0.564 -0.48  0.89
## 42                     polintr.lvl1:all.parties.lvl2NE age  0.02 0.04   296.63  0.55 0.585 -0.06  0.11
## 43                 polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04   213.48  0.98 0.326 -0.04  0.13
## 44                   polintr.lvl1:all.parties.lvl2NE other  0.02 0.08  1261.73  0.23 0.820 -0.14  0.18
## 45                  polintr.lvl1:all.parties.lvl2No answer -0.24 0.70 31717.96 -0.34 0.735 -1.60  1.13
## 46                polintr.lvl1:all.parties.lvl2Other party  0.10 0.03   141.88  3.82 0.000  0.05  0.16
## 47      polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04   250.77  4.31 0.000  0.09  0.25
```

```r
(VC.H5.exp.intr.mod5<-getVC(H5.exp.intr.mod5))
```

```
##             grp                      var1                      var2       est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.197807826  3.912794e-02
## 2  voting.group              environ.lvl1                      <NA>  0.041731515  1.741519e-03
## 3  voting.group              polintr.lvl1                      <NA>  0.055426253  3.072070e-03
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.035181558  1.237742e-03
## 5  voting.group               (Intercept)              environ.lvl1  0.017532777  1.447299e-04
## 6  voting.group               (Intercept)              polintr.lvl1  0.276283508  3.029102e-03
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1 -0.001200754 -8.356269e-06
## 8  voting.group              environ.lvl1              polintr.lvl1 -0.035818963 -8.285003e-05
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.167315108 -2.456486e-04
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1 -0.277168727 -5.404740e-04
## 11        cntry               (Intercept)                      <NA>  0.481871385  2.322000e-01
## 12        cntry              environ.lvl1                      <NA>  0.039364305  1.549549e-03
## 13        cntry              polintr.lvl1                      <NA>  0.049560668  2.456260e-03
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>  0.009643158  9.299049e-05
## 15        cntry               (Intercept)              environ.lvl1 -0.083856128 -1.590628e-03
## 16        cntry               (Intercept)              polintr.lvl1  0.401380310  9.585711e-03
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1  0.022251889  1.033992e-04
## 18        cntry              environ.lvl1              polintr.lvl1  0.118916991  2.319977e-04
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.974694665  3.699904e-04
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.337502509  1.612996e-04
## 21     Residual                      <NA>                      <NA>  1.024870011  1.050359e+00
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
##          all.parties.lvl2 polintr.lvl1.trend         SE  df    asymp.LCL  asymp.UCL
## 1  Anti-immigration party      -0.0227014845 0.02685978 Inf -0.075345689 0.02994272
## 2            Did not vote       0.0546378034 0.02260701 Inf  0.010328875 0.09894673
## 3              Don't know       0.0782603751 0.04349991 Inf -0.006997881 0.16351863
## 4            Invalid vote       0.1790907236 0.34924487 Inf -0.505416639 0.86359809
## 5                  NE age       0.0007101869 0.03698999 Inf -0.071788858 0.07320923
## 6              NE citizen       0.0211170667 0.03896882 Inf -0.055260415 0.09749455
## 7                NE other      -0.0039245266 0.07961813 Inf -0.159973185 0.15212413
## 8               No answer      -0.2589513708 0.69734588 Inf -1.625724177 1.10782144
## 9             Other party       0.0811356354 0.01628700 Inf  0.049213706 0.11305756
## 10  Pro-environment party       0.1511936282 0.03394109 Inf  0.084670320 0.21771694
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
##                     group polintr.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party              -0.02 0.03 0.3980 1.0000     -0.08      0.03
## 2            Did not vote               0.05 0.02 0.0157 0.1252      0.01      0.10
## 3              Don't know               0.08 0.04 0.0720 0.5040     -0.01      0.16
## 4            Invalid vote               0.18 0.35 0.6081 1.0000     -0.51      0.86
## 5                  NE age               0.00 0.04 0.9847 1.0000     -0.07      0.07
## 6              NE citizen               0.02 0.04 0.5879 1.0000     -0.06      0.10
## 7                NE other               0.00 0.08 0.9607 1.0000     -0.16      0.15
## 8               No answer              -0.26 0.70 0.7104 1.0000     -1.63      1.11
## 9             Other party               0.08 0.02 0.0000 0.0000      0.05      0.11
## 10  Pro-environment party               0.15 0.03 0.0000 0.0001      0.08      0.22
```

```r
write.csv2(H5.exp.intr.mod5.trends.tab,"H5.exp.intr.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df  z.ratio     p.value
## 1 Pro-environment party - Anti-immigration party 0.1738951 0.04032848 Inf 4.311968 1.61808e-05
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
##                                         contrast   estimate         SE  df  z.ratio     p.value
## 1           Other party - Anti-immigration party 0.10383712 0.02716097 Inf 3.823028 1.31823e-04
## 2 Pro-environment party - Anti-immigration party 0.17389511 0.04032848 Inf 4.311968 1.61808e-05
## 3            Pro-environment party - Other party 0.07005799 0.03409876 Inf 2.054561 3.99214e-02
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                                Eff   Est   SE       df     t     p    LL    UL
## 1                                                      (Intercept) -0.47 0.14    62.03 -3.25 0.002 -0.76 -0.18
## 2                                                              age  0.00 0.00 35620.83  3.62 0.000  0.00  0.00
## 3                                                           gender  0.07 0.01 35529.71  5.95 0.000  0.05  0.09
## 4                                                             educ  0.01 0.00 35596.16  5.91 0.000  0.01  0.01
## 5                                                            resid -0.06 0.01 35609.27 -4.88 0.000 -0.08 -0.03
## 6                                    occupClerical support workers -0.04 0.09 35481.74 -0.50 0.614 -0.22  0.13
## 7                            occupCraft and related trades workers -0.06 0.09 35494.16 -0.69 0.488 -0.23  0.11
## 8                                      occupElementary occupations  0.02 0.09 35487.61  0.21 0.832 -0.15  0.19
## 9                                                    occupManagers -0.03 0.09 35477.83 -0.29 0.774 -0.20  0.15
## 10                                    occupOther: Not in paid work  0.13 0.09 35545.94  1.41 0.157 -0.05  0.31
## 11                occupPlant and machine operators, and assemblers -0.06 0.09 35494.22 -0.70 0.482 -0.24  0.11
## 12                                              occupProfessionals  0.06 0.09 35482.92  0.69 0.490 -0.11  0.23
## 13                                                    occupRetired -0.03 0.10 35469.46 -0.30 0.761 -0.22  0.16
## 14                                  occupService and sales workers -0.03 0.09 35484.89 -0.37 0.715 -0.20  0.14
## 15         occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35505.84 -0.58 0.563 -0.24  0.13
## 16                    occupTechnicians and associate professionals -0.03 0.09 35480.69 -0.39 0.693 -0.21  0.14
## 17                                                 occupUnemployed -0.02 0.11 35500.20 -0.16 0.872 -0.23  0.19
## 18                                                    environ.lvl1  0.11 0.02   109.85  5.63 0.000  0.07  0.15
## 19                                                    polintr.lvl1 -0.02 0.03   122.70 -0.83 0.410 -0.08  0.03
## 20                                    all.parties.lvl2Did not vote  0.44 0.06   182.54  6.91 0.000  0.31  0.57
## 21                                      all.parties.lvl2Don't know  0.42 0.07   267.96  5.93 0.000  0.28  0.56
## 22                                    all.parties.lvl2Invalid vote  0.44 0.36  1126.05  1.23 0.217 -0.26  1.14
## 23                                          all.parties.lvl2NE age  0.72 0.07   272.67 10.14 0.000  0.58  0.86
## 24                                      all.parties.lvl2NE citizen  0.85 0.08   268.30 11.14 0.000  0.70  1.00
## 25                                        all.parties.lvl2NE other  0.71 0.10   652.83  7.31 0.000  0.52  0.91
## 26                                       all.parties.lvl2No answer  0.58 0.37  1357.83  1.54 0.123 -0.16  1.31
## 27                                     all.parties.lvl2Other party  0.55 0.05   230.38 11.31 0.000  0.45  0.64
## 28                           all.parties.lvl2Pro-environment party  0.92 0.07   257.35 14.11 0.000  0.79  1.05
## 29                                       environ.lvl1:polintr.lvl1 -0.01 0.02   126.69 -0.38 0.708 -0.04  0.03
## 30                       environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    87.78 -0.09 0.926 -0.05  0.04
## 31                         environ.lvl1:all.parties.lvl2Don't know  0.02 0.03   342.47  0.70 0.484 -0.04  0.09
## 32                       environ.lvl1:all.parties.lvl2Invalid vote  0.04 0.39 25264.47  0.10 0.924 -0.73  0.80
## 33                             environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.97  1.22 0.223 -0.02  0.10
## 34                         environ.lvl1:all.parties.lvl2NE citizen -0.05 0.03   236.46 -1.39 0.166 -0.11  0.02
## 35                           environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1205.32  0.63 0.530 -0.07  0.14
## 36                          environ.lvl1:all.parties.lvl2No answer  0.17 0.26 10536.40  0.65 0.514 -0.34  0.67
## 37                        environ.lvl1:all.parties.lvl2Other party  0.01 0.02   136.75  0.72 0.471 -0.02  0.05
## 38              environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03   199.54  0.26 0.797 -0.05  0.06
## 39                       polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03    94.24  2.46 0.016  0.01  0.14
## 40                         polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05   421.96  2.04 0.042  0.00  0.20
## 41                       polintr.lvl1:all.parties.lvl2Invalid vote  0.21 0.35 10518.93  0.60 0.550 -0.48  0.90
## 42                             polintr.lvl1:all.parties.lvl2NE age  0.02 0.04   293.07  0.56 0.578 -0.06  0.11
## 43                         polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04   211.95  0.80 0.423 -0.05  0.12
## 44                           polintr.lvl1:all.parties.lvl2NE other  0.02 0.08  1259.27  0.20 0.844 -0.15  0.18
## 45                          polintr.lvl1:all.parties.lvl2No answer -0.26 0.70 31722.95 -0.38 0.707 -1.64  1.11
## 46                        polintr.lvl1:all.parties.lvl2Other party  0.10 0.03   140.44  3.80 0.000  0.05  0.16
## 47              polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04   249.54  4.30 0.000  0.09  0.25
## 48          environ.lvl1:polintr.lvl1:all.parties.lvl2Did not vote  0.02 0.02    85.68  0.71 0.477 -0.03  0.06
## 49            environ.lvl1:polintr.lvl1:all.parties.lvl2Don't know  0.03 0.04   448.49  0.86 0.389 -0.04  0.11
## 50          environ.lvl1:polintr.lvl1:all.parties.lvl2Invalid vote  0.14 0.41 27792.84  0.34 0.736 -0.67  0.95
## 51                environ.lvl1:polintr.lvl1:all.parties.lvl2NE age  0.01 0.03   364.25  0.22 0.824 -0.06  0.08
## 52            environ.lvl1:polintr.lvl1:all.parties.lvl2NE citizen  0.09 0.03   215.95  2.56 0.011  0.02  0.15
## 53              environ.lvl1:polintr.lvl1:all.parties.lvl2NE other -0.08 0.06  1155.29 -1.37 0.172 -0.21  0.04
## 54             environ.lvl1:polintr.lvl1:all.parties.lvl2No answer  0.36 0.76 34661.69  0.47 0.637 -1.14  1.86
## 55           environ.lvl1:polintr.lvl1:all.parties.lvl2Other party  0.02 0.02   124.68  1.10 0.272 -0.02  0.06
## 56 environ.lvl1:polintr.lvl1:all.parties.lvl2Pro-environment party -0.00 0.03   258.63 -0.06 0.956 -0.06  0.06
```

```r
(VC.H5.exp.intr.mod6<-getVC(H5.exp.intr.mod6))
```

```
##             grp                      var1                      var2       est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.198111360  3.924811e-02
## 2  voting.group              environ.lvl1                      <NA>  0.041542864  1.725810e-03
## 3  voting.group              polintr.lvl1                      <NA>  0.055319665  3.060265e-03
## 4  voting.group environ.lvl1:polintr.lvl1                      <NA>  0.030462894  9.279879e-04
## 5  voting.group               (Intercept)              environ.lvl1  0.009772025  8.042488e-05
## 6  voting.group               (Intercept)              polintr.lvl1  0.278234232  3.049295e-03
## 7  voting.group               (Intercept) environ.lvl1:polintr.lvl1 -0.005343558 -3.224861e-05
## 8  voting.group              environ.lvl1              polintr.lvl1 -0.038065623 -8.748003e-05
## 9  voting.group              environ.lvl1 environ.lvl1:polintr.lvl1 -0.118625456 -1.501224e-04
## 10 voting.group              polintr.lvl1 environ.lvl1:polintr.lvl1 -0.301971679 -5.088818e-04
## 11        cntry               (Intercept)                      <NA>  0.481662638  2.319989e-01
## 12        cntry              environ.lvl1                      <NA>  0.039339237  1.547576e-03
## 13        cntry              polintr.lvl1                      <NA>  0.049623564  2.462498e-03
## 14        cntry environ.lvl1:polintr.lvl1                      <NA>  0.009958964  9.918096e-05
## 15        cntry               (Intercept)              environ.lvl1 -0.084611697 -1.603243e-03
## 16        cntry               (Intercept)              polintr.lvl1  0.403867631  9.653170e-03
## 17        cntry               (Intercept) environ.lvl1:polintr.lvl1  0.040551719  1.945210e-04
## 18        cntry              environ.lvl1              polintr.lvl1  0.128007153  2.498896e-04
## 19        cntry              environ.lvl1 environ.lvl1:polintr.lvl1  0.966486879  3.786483e-04
## 20        cntry              polintr.lvl1 environ.lvl1:polintr.lvl1  0.377774008  1.866956e-04
## 21     Residual                      <NA>                      <NA>  1.024822113  1.050260e+00
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
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
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
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    62.03 -3.25 0.002 -0.76 -0.18
## 2                                                      age  0.00 0.00 35620.81  3.62 0.000  0.00  0.00
## 3                                                   gender  0.07 0.01 35529.71  5.95 0.000  0.05  0.09
## 4                                                     educ  0.01 0.00 35596.15  5.91 0.000  0.01  0.01
## 5                                                    resid -0.06 0.01 35609.28 -4.88 0.000 -0.08 -0.03
## 6                            occupClerical support workers -0.04 0.09 35481.74 -0.50 0.614 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35494.16 -0.69 0.488 -0.23  0.11
## 8                              occupElementary occupations  0.02 0.09 35487.61  0.21 0.832 -0.15  0.19
## 9                                            occupManagers -0.03 0.09 35477.83 -0.29 0.774 -0.20  0.15
## 10                            occupOther: Not in paid work  0.13 0.09 35545.93  1.41 0.157 -0.05  0.31
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35494.22 -0.70 0.482 -0.24  0.11
## 12                                      occupProfessionals  0.06 0.09 35482.93  0.69 0.490 -0.11  0.23
## 13                                            occupRetired -0.03 0.10 35469.46 -0.30 0.761 -0.22  0.16
## 14                          occupService and sales workers -0.03 0.09 35484.89 -0.37 0.715 -0.20  0.14
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35505.84 -0.58 0.563 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35480.69 -0.39 0.693 -0.21  0.14
## 17                                         occupUnemployed -0.02 0.11 35500.20 -0.16 0.872 -0.23  0.19
## 18                                            environ.lvl1  0.11 0.02   109.85  5.63 0.000  0.07  0.15
## 19                                            polintr.lvl1 -0.02 0.03   122.70 -0.83 0.410 -0.08  0.03
## 20                                            env.intr.int -0.01 0.02   126.69 -0.38 0.708 -0.04  0.03
## 21                            all.parties.lvl2Did not vote  0.44 0.06   182.53  6.91 0.000  0.31  0.57
## 22                              all.parties.lvl2Don't know  0.42 0.07   267.95  5.93 0.000  0.28  0.56
## 23                            all.parties.lvl2Invalid vote  0.44 0.36  1126.03  1.23 0.217 -0.26  1.14
## 24                                  all.parties.lvl2NE age  0.72 0.07   272.66 10.14 0.000  0.58  0.86
## 25                              all.parties.lvl2NE citizen  0.85 0.08   268.29 11.14 0.000  0.70  1.00
## 26                                all.parties.lvl2NE other  0.71 0.10   652.82  7.31 0.000  0.52  0.91
## 27                               all.parties.lvl2No answer  0.58 0.37  1357.80  1.54 0.123 -0.16  1.31
## 28                             all.parties.lvl2Other party  0.55 0.05   230.37 11.31 0.000  0.45  0.64
## 29                   all.parties.lvl2Pro-environment party  0.92 0.07   257.35 14.11 0.000  0.79  1.05
## 30               environ.lvl1:all.parties.lvl2Did not vote -0.00 0.02    87.78 -0.09 0.926 -0.05  0.04
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.02 0.03   342.46  0.70 0.484 -0.04  0.09
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.04 0.39 25264.46  0.10 0.924 -0.73  0.80
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   238.97  1.22 0.223 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.05 0.03   236.46 -1.39 0.166 -0.11  0.02
## 35                   environ.lvl1:all.parties.lvl2NE other  0.04 0.06  1205.32  0.63 0.530 -0.07  0.14
## 36                  environ.lvl1:all.parties.lvl2No answer  0.17 0.26 10536.39  0.65 0.514 -0.34  0.67
## 37                environ.lvl1:all.parties.lvl2Other party  0.01 0.02   136.75  0.72 0.471 -0.02  0.05
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.01 0.03   199.54  0.26 0.797 -0.05  0.06
## 39               polintr.lvl1:all.parties.lvl2Did not vote  0.08 0.03    94.24  2.46 0.016  0.01  0.14
## 40                 polintr.lvl1:all.parties.lvl2Don't know  0.10 0.05   421.95  2.04 0.042  0.00  0.20
## 41               polintr.lvl1:all.parties.lvl2Invalid vote  0.21 0.35 10518.76  0.60 0.550 -0.48  0.90
## 42                     polintr.lvl1:all.parties.lvl2NE age  0.02 0.04   293.07  0.56 0.578 -0.06  0.11
## 43                 polintr.lvl1:all.parties.lvl2NE citizen  0.04 0.04   211.95  0.80 0.423 -0.05  0.12
## 44                   polintr.lvl1:all.parties.lvl2NE other  0.02 0.08  1259.25  0.20 0.844 -0.15  0.18
## 45                  polintr.lvl1:all.parties.lvl2No answer -0.26 0.70 31722.87 -0.38 0.707 -1.64  1.11
## 46                polintr.lvl1:all.parties.lvl2Other party  0.10 0.03   140.44  3.80 0.000  0.05  0.16
## 47      polintr.lvl1:all.parties.lvl2Pro-environment party  0.17 0.04   249.54  4.30 0.000  0.09  0.25
## 48               env.intr.int:all.parties.lvl2Did not vote  0.02 0.02    85.68  0.71 0.477 -0.03  0.06
## 49                 env.intr.int:all.parties.lvl2Don't know  0.03 0.04   448.48  0.86 0.389 -0.04  0.11
## 50               env.intr.int:all.parties.lvl2Invalid vote  0.14 0.41 27792.75  0.34 0.736 -0.67  0.95
## 51                     env.intr.int:all.parties.lvl2NE age  0.01 0.03   364.24  0.22 0.824 -0.06  0.08
## 52                 env.intr.int:all.parties.lvl2NE citizen  0.09 0.03   215.95  2.56 0.011  0.02  0.15
## 53                   env.intr.int:all.parties.lvl2NE other -0.08 0.06  1155.27 -1.37 0.172 -0.21  0.04
## 54                  env.intr.int:all.parties.lvl2No answer  0.36 0.76 34661.68  0.47 0.637 -1.14  1.86
## 55                env.intr.int:all.parties.lvl2Other party  0.02 0.02   124.68  1.10 0.272 -0.02  0.06
## 56      env.intr.int:all.parties.lvl2Pro-environment party -0.00 0.03   258.62 -0.06 0.956 -0.06  0.06
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
##          all.parties.lvl2 env.intr.int.trend          SE  df    asymp.LCL  asymp.UCL
## 1  Anti-immigration party      -0.0069621407 0.018544266 Inf -0.043308234 0.02938395
## 2            Did not vote       0.0095517468 0.014187703 Inf -0.018255640 0.03735913
## 3              Don't know       0.0267672477 0.034578183 Inf -0.041004746 0.09453924
## 4            Invalid vote       0.1320711625 0.412682858 Inf -0.676772377 0.94091470
## 5                  NE age       0.0007638314 0.029563838 Inf -0.057180227 0.05870789
## 6              NE citizen       0.0795677618 0.028392962 Inf  0.023918580 0.13521694
## 7                NE other      -0.0916793524 0.059222250 Inf -0.207752829 0.02439412
## 8               No answer       0.3529258125 0.762887333 Inf -1.142305884 1.84815751
## 9             Other party       0.0154857242 0.008953365 Inf -0.002062548 0.03303400
## 10  Pro-environment party      -0.0086996200 0.025538314 Inf -0.058753796 0.04135456
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
##                     group env.intr.int.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party              -0.01 0.02 0.7073 1.0000     -0.04      0.03
## 2            Did not vote               0.01 0.01 0.5008 1.0000     -0.02      0.04
## 3              Don't know               0.03 0.03 0.4389 1.0000     -0.04      0.09
## 4            Invalid vote               0.13 0.41 0.7489 1.0000     -0.68      0.94
## 5                  NE age               0.00 0.03 0.9794 1.0000     -0.06      0.06
## 6              NE citizen               0.08 0.03 0.0051 0.0507      0.02      0.14
## 7                NE other              -0.09 0.06 0.1216 0.9729     -0.21      0.02
## 8               No answer               0.35 0.76 0.6436 1.0000     -1.14      1.85
## 9             Other party               0.02 0.01 0.0837 0.7533      0.00      0.03
## 10  Pro-environment party              -0.01 0.03 0.7334 1.0000     -0.06      0.04
```

```r
write.csv2(H5.exp.intr.mod6.trends.tab,"H5.exp.intr.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.exp.intr.contrast<-data.frame(pairs(H5.exp.intr.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast     estimate        SE  df     z.ratio   p.value
## 1 Pro-environment party - Anti-immigration party -0.001737479 0.0314156 Inf -0.05530626 0.9558945
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
##                                         contrast     estimate         SE  df     z.ratio   p.value
## 1           Other party - Anti-immigration party  0.022447865 0.02035632 Inf  1.10274645 0.8104119
## 2 Pro-environment party - Anti-immigration party -0.001737479 0.03141560 Inf -0.05530626 0.9558945
## 3            Pro-environment party - Other party -0.024185344 0.02687443 Inf -0.89993896 0.8104119
```


\newpage



# Exploratory Hypothesis 5: The strength of the association between environment and refugee attitudes is stronger among those who consume political news more

Omit missing values on political news consumption item


```r
dat.H5.news<-dat %>%
  filter(!is.na(polnews))
```

### Model 1: without interactions (only main effects)


```r
H5.exp.news.mod1<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))


isSingular(H5.exp.news.mod1)
```

```
## [1] FALSE
```

```r
(FE.H5.exp.news.mod1<-getFE(H5.exp.news.mod1))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.49  0.51 0.612 -0.21  0.36
## 2                                                      age  0.00 0.00 33590.27  3.30 0.001  0.00  0.00
## 3                                                   gender  0.06 0.01 35308.57  5.07 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35436.41  7.31 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35377.31 -5.19 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35251.11 -0.49 0.627 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35261.80 -0.72 0.472 -0.24  0.11
## 8                              occupElementary occupations  0.02 0.09 35263.86  0.18 0.854 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35250.42 -0.19 0.852 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35420.20  1.55 0.122 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35259.33 -0.71 0.479 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35255.30  0.83 0.404 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35250.00 -0.36 0.717 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35257.25 -0.43 0.664 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35266.06 -0.56 0.572 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35251.22 -0.37 0.715 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35270.08 -0.23 0.820 -0.24  0.19
## 18                                            environ.lvl1  0.12 0.01    19.41 11.66 0.000  0.10  0.14
## 19                                            polnews.lvl1  0.03 0.01 35312.08  4.45 0.000  0.01  0.04
```

```r
(VC.H5.exp.news.mod1<-getVC(H5.exp.news.mod1))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30827355  0.0950325843
## 2 voting.group environ.lvl1         <NA>  0.04328811  0.0018738605
## 3 voting.group  (Intercept) environ.lvl1  0.11077227  0.0014782094
## 4        cntry  (Intercept)         <NA>  0.49919789  0.2491985382
## 5        cntry environ.lvl1         <NA>  0.03896499  0.0015182701
## 6        cntry  (Intercept) environ.lvl1 -0.02792210 -0.0005431194
## 7     Residual         <NA>         <NA>  1.02912215  1.0590924056
```

\newpage

### Model 2: Level-1 interaction between environmental attitudes and political polnews


```r
H5.exp.news.mod2<-lmer(refugees~(environ.lvl1|voting.group)+
                (environ.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.news.mod2)
```

```
## [1] FALSE
```

```r
anova(H5.exp.news.mod1,H5.exp.news.mod2)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod1: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.news.mod1:     age + gender + educ + resid + occup + environ.lvl1 + polnews.lvl1
## H5.exp.news.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.news.mod2:     age + gender + educ + resid + occup + environ.lvl1 + polnews.lvl1 + 
## H5.exp.news.mod2:     environ.lvl1:polnews.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.news.mod1   26 103544 103764 -51746   103492                     
## H5.exp.news.mod2   27 103546 103774 -51746   103492 0.1524  1     0.6963
```

```r
(FE.H5.exp.news.mod2<-getFE(H5.exp.news.mod2))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.49  0.51 0.612 -0.21  0.36
## 2                                                      age  0.00 0.00 33587.67  3.30 0.001  0.00  0.00
## 3                                                   gender  0.06 0.01 35308.54  5.06 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35436.51  7.31 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35377.32 -5.19 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35251.12 -0.49 0.628 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35261.80 -0.72 0.472 -0.24  0.11
## 8                              occupElementary occupations  0.02 0.09 35263.86  0.18 0.854 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35250.43 -0.19 0.852 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35420.32  1.54 0.123 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35259.34 -0.71 0.479 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35255.31  0.83 0.404 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35250.02 -0.36 0.717 -0.23  0.16
## 14                          occupService and sales workers -0.04 0.09 35257.26 -0.43 0.664 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35266.06 -0.56 0.573 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35251.23 -0.37 0.715 -0.20  0.14
## 17                                         occupUnemployed -0.02 0.11 35270.03 -0.23 0.821 -0.24  0.19
## 18                                            environ.lvl1  0.12 0.01    19.41 11.66 0.000  0.10  0.14
## 19                                            polnews.lvl1  0.03 0.01 35312.86  4.45 0.000  0.01  0.04
## 20                               environ.lvl1:polnews.lvl1 -0.00 0.00 35279.06 -0.39 0.696 -0.01  0.01
```

```r
(VC.H5.exp.news.mod2<-getVC(H5.exp.news.mod2))
```

```
##            grp         var1         var2      est_SD       est_SD2
## 1 voting.group  (Intercept)         <NA>  0.30828779  0.0950413597
## 2 voting.group environ.lvl1         <NA>  0.04325288  0.0018708114
## 3 voting.group  (Intercept) environ.lvl1  0.11113357  0.0014818921
## 4        cntry  (Intercept)         <NA>  0.49921640  0.2492170167
## 5        cntry environ.lvl1         <NA>  0.03898152  0.0015195590
## 6        cntry  (Intercept) environ.lvl1 -0.02812737 -0.0005473646
## 7     Residual         <NA>         <NA>  1.02912069  1.0590893901
```

\newpage

### Model 3: Level-1 interaction between environmental attitudes and political polnews, allow polnews effect to vary between voting groups and countries


```r
H5.exp.news.mod3<-lmer(refugees~(environ.lvl1+polnews.lvl1|voting.group)+
                (environ.lvl1+polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
isSingular(H5.exp.news.mod3)
```

```
## [1] FALSE
```

```r
anova(H5.exp.news.mod2,H5.exp.news.mod3)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod2: refugees ~ (environ.lvl1 | voting.group) + (environ.lvl1 | cntry) + 
## H5.exp.news.mod2:     age + gender + educ + resid + occup + environ.lvl1 + polnews.lvl1 + 
## H5.exp.news.mod2:     environ.lvl1:polnews.lvl1
## H5.exp.news.mod3: refugees ~ (environ.lvl1 + polnews.lvl1 | voting.group) + (environ.lvl1 + 
## H5.exp.news.mod3:     polnews.lvl1 | cntry) + age + gender + educ + resid + occup + 
## H5.exp.news.mod3:     environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
## H5.exp.news.mod2   27 103546 103774 -51746   103492                       
## H5.exp.news.mod3   33 103543 103823 -51739   103477 14.116  6    0.02836 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod3<-getFE(H5.exp.news.mod3))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.48  0.51 0.611 -0.21  0.36
## 2                                                      age  0.00 0.00 32648.32  3.49 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35306.61  5.09 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35392.24  7.34 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35352.96 -5.10 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35234.96 -0.48 0.633 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35239.84 -0.72 0.472 -0.24  0.11
## 8                              occupElementary occupations  0.02 0.09 35243.26  0.20 0.842 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35232.28 -0.19 0.848 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35406.43  1.54 0.124 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35245.57 -0.69 0.490 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35236.68  0.83 0.409 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35230.83 -0.39 0.699 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35236.06 -0.44 0.662 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35251.34 -0.55 0.583 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35231.31 -0.36 0.719 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35249.83 -0.26 0.796 -0.24  0.18
## 18                                            environ.lvl1  0.12 0.01    19.42 11.65 0.000  0.10  0.14
## 19                                            polnews.lvl1  0.02 0.01    22.50  2.66 0.014  0.01  0.04
## 20                               environ.lvl1:polnews.lvl1 -0.00 0.00 35270.99 -0.35 0.725 -0.01  0.01
```

```r
(VC.H5.exp.news.mod3<-getVC(H5.exp.news.mod3))
```

```
##             grp         var1         var2       est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.308661099  9.527167e-02
## 2  voting.group environ.lvl1         <NA>  0.043203926  1.866579e-03
## 3  voting.group polnews.lvl1         <NA>  0.029900588  8.940451e-04
## 4  voting.group  (Intercept) environ.lvl1  0.113331427  1.511317e-03
## 5  voting.group  (Intercept) polnews.lvl1  0.039458622  3.641695e-04
## 6  voting.group environ.lvl1 polnews.lvl1  0.002452728  3.168490e-06
## 7         cntry  (Intercept)         <NA>  0.499102349  2.491032e-01
## 8         cntry environ.lvl1         <NA>  0.038943792  1.516619e-03
## 9         cntry polnews.lvl1         <NA>  0.026970238  7.273937e-04
## 10        cntry  (Intercept) environ.lvl1 -0.030940455 -6.013877e-04
## 11        cntry  (Intercept) polnews.lvl1  0.436423393  5.874656e-03
## 12        cntry environ.lvl1 polnews.lvl1  0.076178478  8.001203e-05
## 13     Residual         <NA>         <NA>  1.028345791  1.057495e+00
```


\newpage

### Model 4: Allow the level-1 interaction to vary between voting groups and countries


```r
H5.exp.news.mod4<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+
                   environ.lvl1:polnews.lvl1|voting.group)+
                (environ.lvl1+polnews.lvl1+
                   environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Model failed to converge with 1 negative eigenvalue: -1.2e+03
```

```r
isSingular(H5.exp.news.mod4)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod3,H5.exp.news.mod4)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod3: refugees ~ (environ.lvl1 + polnews.lvl1 | voting.group) + (environ.lvl1 + 
## H5.exp.news.mod3:     polnews.lvl1 | cntry) + age + gender + educ + resid + occup + 
## H5.exp.news.mod3:     environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1
## H5.exp.news.mod4: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod4:     polnews.lvl1 + environ.lvl1:polnews.lvl1
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)   
## H5.exp.news.mod3   33 103543 103823 -51739   103477                        
## H5.exp.news.mod4   41 103538 103885 -51728   103456 21.631  8   0.005647 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod4<-getFE(H5.exp.news.mod4))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept)  0.07 0.14    49.45  0.51 0.610 -0.21  0.36
## 2                                                      age  0.00 0.00 32551.27  3.50 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35279.89  5.09 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35331.70  7.38 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35340.92 -5.09 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35216.24 -0.47 0.635 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35219.22 -0.72 0.474 -0.24  0.11
## 8                              occupElementary occupations  0.02 0.09 35218.20  0.19 0.850 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35204.19 -0.20 0.842 -0.19  0.16
## 10                            occupOther: Not in paid work  0.14 0.09 35385.77  1.56 0.119 -0.04  0.32
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35216.81 -0.69 0.487 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35214.37  0.82 0.412 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35205.67 -0.38 0.701 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35211.38 -0.44 0.658 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35233.55 -0.55 0.583 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35207.56 -0.36 0.716 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35215.03 -0.28 0.776 -0.24  0.18
## 18                                            environ.lvl1  0.12 0.01    19.46 11.72 0.000  0.10  0.14
## 19                                            polnews.lvl1  0.02 0.01    22.39  2.65 0.015  0.00  0.04
## 20                               environ.lvl1:polnews.lvl1 -0.00 0.01    19.30 -0.14 0.887 -0.02  0.01
```

```r
(VC.H5.exp.news.mod4<-getVC(H5.exp.news.mod4))
```

```
##             grp                      var1                      var2       est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.308107232  9.493007e-02
## 2  voting.group              environ.lvl1                      <NA>  0.042898183  1.840254e-03
## 3  voting.group              polnews.lvl1                      <NA>  0.030868161  9.528434e-04
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.018365200  3.372806e-04
## 5  voting.group               (Intercept)              environ.lvl1  0.110342172  1.458419e-03
## 6  voting.group               (Intercept)              polnews.lvl1  0.022461296  2.136227e-04
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1  0.615115341  3.480600e-03
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.048684972 -6.446806e-05
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.104949559 -8.268280e-05
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.745042238 -4.223644e-04
## 11        cntry               (Intercept)                      <NA>  0.499211808  2.492124e-01
## 12        cntry              environ.lvl1                      <NA>  0.038514255  1.483348e-03
## 13        cntry              polnews.lvl1                      <NA>  0.026109935  6.817287e-04
## 14        cntry environ.lvl1:polnews.lvl1                      <NA>  0.023597931  5.568624e-04
## 15        cntry               (Intercept)              environ.lvl1 -0.009048442 -1.739723e-04
## 16        cntry               (Intercept)              polnews.lvl1  0.491309760  6.403922e-03
## 17        cntry               (Intercept) environ.lvl1:polnews.lvl1 -0.395609561 -4.660425e-03
## 18        cntry              environ.lvl1              polnews.lvl1  0.089737428  9.024038e-05
## 19        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.732075194  6.653515e-04
## 20        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.338042424 -2.082816e-04
## 21     Residual                      <NA>                      <NA>  1.027788635  1.056349e+00
```

```r
theta <- getME(H5.exp.news.mod4,"theta")

## diagonal elements are identifiable because they are fitted
##  with a lower bound of zero ...
diag.element <- getME(H5.exp.news.mod4,"lower")==0
any(theta[diag.element]<1e-5)
```

```
## [1] TRUE
```

```r
round(theta,5)
```

```
##                            voting.group.(Intercept)               voting.group.environ.lvl1.(Intercept) 
##                                             0.29978                                             0.00461 
##               voting.group.polnews.lvl1.(Intercept)  voting.group.environ.lvl1:polnews.lvl1.(Intercept) 
##                                             0.00067                                             0.01099 
##                           voting.group.environ.lvl1              voting.group.polnews.lvl1.environ.lvl1 
##                                             0.04148                                            -0.00155 
## voting.group.environ.lvl1:polnews.lvl1.environ.lvl1                           voting.group.polnews.lvl1 
##                                            -0.00311                                             0.02999 
## voting.group.environ.lvl1:polnews.lvl1.polnews.lvl1              voting.group.environ.lvl1:polnews.lvl1 
##                                            -0.01374                                             0.00000 
##                                   cntry.(Intercept)                      cntry.environ.lvl1.(Intercept) 
##                                             0.48571                                            -0.00034 
##                      cntry.polnews.lvl1.(Intercept)         cntry.environ.lvl1:polnews.lvl1.(Intercept) 
##                                             0.01248                                            -0.00908 
##                                  cntry.environ.lvl1                     cntry.polnews.lvl1.environ.lvl1 
##                                             0.03747                                             0.00239 
##        cntry.environ.lvl1:polnews.lvl1.environ.lvl1                                  cntry.polnews.lvl1 
##                                             0.01673                                             0.02200 
##        cntry.environ.lvl1:polnews.lvl1.polnews.lvl1                     cntry.environ.lvl1:polnews.lvl1 
##                                            -0.00563                                             0.01154
```

\newpage

### Model 5: Enter voting group variable and both two-way interactions with environment and polnews


```r
H5.exp.news.mod5<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|voting.group)+
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```r
isSingular(H5.exp.news.mod5)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod4,H5.exp.news.mod5)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod4: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod4:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod4:     polnews.lvl1 + environ.lvl1:polnews.lvl1
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
## H5.exp.news.mod4   41 103538 103885 -51728   103456                         
## H5.exp.news.mod5   68 103396 103973 -51630   103260 195.65 27  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
(FE.H5.exp.news.mod5<-getFE(H5.exp.news.mod5))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.48 0.14    62.06 -3.30 0.002 -0.76 -0.19
## 2                                                      age  0.00 0.00 35050.90  3.82 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35303.96  5.07 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35338.83  7.39 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35389.26 -5.04 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35273.69 -0.49 0.623 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35278.10 -0.70 0.483 -0.23  0.11
## 8                              occupElementary occupations  0.01 0.09 35273.68  0.16 0.872 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35262.22 -0.21 0.837 -0.19  0.16
## 10                            occupOther: Not in paid work  0.12 0.09 35332.02  1.37 0.171 -0.05  0.30
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35276.08 -0.69 0.493 -0.23  0.11
## 12                                      occupProfessionals  0.07 0.09 35273.26  0.81 0.417 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35250.31 -0.40 0.686 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35270.09 -0.44 0.660 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35295.96 -0.56 0.574 -0.23  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35267.41 -0.36 0.720 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35260.13 -0.32 0.752 -0.24  0.18
## 18                                            environ.lvl1  0.11 0.02   108.30  5.59 0.000  0.07  0.15
## 19                                            polnews.lvl1  0.02 0.02   141.31  1.25 0.214 -0.01  0.06
## 20                            all.parties.lvl2Did not vote  0.45 0.06   179.10  7.10 0.000  0.32  0.57
## 21                              all.parties.lvl2Don't know  0.43 0.07   269.00  6.05 0.000  0.29  0.57
## 22                            all.parties.lvl2Invalid vote  0.49 0.36  1119.81  1.37 0.172 -0.21  1.19
## 23                                  all.parties.lvl2NE age  0.74 0.07   273.06 10.44 0.000  0.60  0.88
## 24                              all.parties.lvl2NE citizen  0.87 0.08   269.18 11.43 0.000  0.72  1.02
## 25                                all.parties.lvl2NE other  0.72 0.10   658.30  7.36 0.000  0.53  0.91
## 26                               all.parties.lvl2No answer  0.54 0.38  1387.11  1.45 0.148 -0.19  1.28
## 27                             all.parties.lvl2Other party  0.55 0.05   229.57 11.49 0.000  0.46  0.65
## 28                   all.parties.lvl2Pro-environment party  0.92 0.06   255.26 14.16 0.000  0.79  1.05
## 29                               environ.lvl1:polnews.lvl1 -0.00 0.01    19.42 -0.41 0.683 -0.02  0.01
## 30               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    84.50  0.03 0.977 -0.04  0.05
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   333.34  0.92 0.356 -0.04  0.10
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.05 0.39 25029.45  0.13 0.896 -0.71  0.81
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   232.31  1.25 0.212 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.03 0.03   229.81 -1.00 0.319 -0.10  0.03
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1219.70  0.59 0.554 -0.08  0.14
## 36                  environ.lvl1:all.parties.lvl2No answer  0.25 0.23  7350.89  1.06 0.288 -0.21  0.70
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   131.72  0.94 0.351 -0.02  0.06
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   190.28  0.96 0.339 -0.03  0.08
## 39               polnews.lvl1:all.parties.lvl2Did not vote -0.02 0.02    94.21 -0.75 0.457 -0.06  0.03
## 40                 polnews.lvl1:all.parties.lvl2Don't know -0.06 0.04   500.80 -1.67 0.095 -0.13  0.01
## 41               polnews.lvl1:all.parties.lvl2Invalid vote -0.06 0.34 22217.09 -0.17 0.863 -0.73  0.61
## 42                     polnews.lvl1:all.parties.lvl2NE age -0.02 0.03   340.96 -0.56 0.576 -0.08  0.05
## 43                 polnews.lvl1:all.parties.lvl2NE citizen -0.02 0.04   314.07 -0.43 0.669 -0.08  0.05
## 44                   polnews.lvl1:all.parties.lvl2NE other -0.03 0.06  2835.56 -0.44 0.660 -0.15  0.10
## 45                  polnews.lvl1:all.parties.lvl2No answer -0.38 0.42 29170.71 -0.91 0.363 -1.21  0.44
## 46                polnews.lvl1:all.parties.lvl2Other party  0.01 0.02   160.97  0.56 0.577 -0.03  0.05
## 47      polnews.lvl1:all.parties.lvl2Pro-environment party  0.00 0.03   283.10  0.13 0.893 -0.06  0.06
```

```r
(VC.H5.exp.news.mod5<-getVC(H5.exp.news.mod5))
```

```
##             grp                      var1                      var2      est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.19787604  3.915493e-02
## 2  voting.group              environ.lvl1                      <NA>  0.04070099  1.656571e-03
## 3  voting.group              polnews.lvl1                      <NA>  0.02760802  7.622026e-04
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.01552929  2.411588e-04
## 5  voting.group               (Intercept)              environ.lvl1  0.08152100  6.565498e-04
## 6  voting.group               (Intercept)              polnews.lvl1 -0.01322898 -7.226946e-05
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1 -0.33110264 -1.017437e-03
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.03344998 -3.758686e-05
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.23299967 -1.472692e-04
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.90901803 -3.897259e-04
## 11        cntry               (Intercept)                      <NA>  0.48213103  2.324503e-01
## 12        cntry              environ.lvl1                      <NA>  0.03910401  1.529124e-03
## 13        cntry              polnews.lvl1                      <NA>  0.02697269  7.275259e-04
## 14        cntry environ.lvl1:polnews.lvl1                      <NA>  0.02397377  5.747417e-04
## 15        cntry               (Intercept)              environ.lvl1 -0.01261116 -2.377615e-04
## 16        cntry               (Intercept)              polnews.lvl1  0.45249552  5.884419e-03
## 17        cntry               (Intercept) environ.lvl1:polnews.lvl1 -0.40038518 -4.627852e-03
## 18        cntry              environ.lvl1              polnews.lvl1  0.08462077  8.925295e-05
## 19        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.75421565  7.070551e-04
## 20        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.36305362 -2.347639e-04
## 21     Residual                      <NA>                      <NA>  1.02774420  1.056258e+00
```

\newpage

#### Look among which voting group there is strongest association between polnews and refugee attitudes


```r
H5.exp.news.mod5.trends<-emtrends(H5.exp.news.mod5,specs = c("all.parties.lvl2"),var=c("polnews.lvl1"))
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'pbkrtest.limit = 35483' (or larger)
## [or, globally, 'set emm_options(pbkrtest.limit = 35483)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```
## Note: D.f. calculations have been disabled because the number of observations exceeds 3000.
## To enable adjustments, add the argument 'lmerTest.limit = 35483' (or larger)
## [or, globally, 'set emm_options(lmerTest.limit = 35483)' or larger];
## but be warned that this may result in large computation time and memory use.
```

```r
(H5.exp.news.mod5.trends.tab<-data.frame(H5.exp.news.mod5.trends))
```

```
##          all.parties.lvl2 polnews.lvl1.trend         SE  df   asymp.LCL  asymp.UCL
## 1  Anti-immigration party        0.023895565 0.01914989 Inf -0.01363752 0.06142865
## 2            Did not vote        0.007253473 0.01444551 Inf -0.02105920 0.03556615
## 3              Don't know       -0.035982728 0.03146913 Inf -0.09766109 0.02569563
## 4            Invalid vote       -0.034937695 0.34015378 Inf -0.70162686 0.63175147
## 5                  NE age        0.005899655 0.02718764 Inf -0.04738713 0.05918644
## 6              NE citizen        0.008775561 0.03094524 Inf -0.05187600 0.06942712
## 7                NE other       -0.003699292 0.06042415 Inf -0.12212845 0.11472987
## 8               No answer       -0.361024641 0.42242427 Inf -1.18896100 0.46691171
## 9             Other party        0.035171778 0.01076074 Inf  0.01408112 0.05626243
## 10  Pro-environment party        0.027950814 0.02482136 Inf -0.02069816 0.07659979
```

```r
H5.exp.news.mod5.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.news.mod5.trends.tab$polnews.lvl1.trend/
                   H5.exp.news.mod5.trends.tab$SE)))
H5.exp.news.mod5.trends.tab$adj.p<-
  p.adjust(H5.exp.news.mod5.trends.tab$p,method="holm")

H5.exp.news.mod5.trends.tab<-
  cbind(group=H5.exp.news.mod5.trends.tab[,1],
      round(H5.exp.news.mod5.trends.tab[,c(2,3)],2),
      round(H5.exp.news.mod5.trends.tab[,c(7,8)],4),
      round(H5.exp.news.mod5.trends.tab[,c(5,6)],2))
H5.exp.news.mod5.trends.tab
```

```
##                     group polnews.lvl1.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party               0.02 0.02 0.2121 1.0000     -0.01      0.06
## 2            Did not vote               0.01 0.01 0.6156 1.0000     -0.02      0.04
## 3              Don't know              -0.04 0.03 0.2529 1.0000     -0.10      0.03
## 4            Invalid vote              -0.03 0.34 0.9182 1.0000     -0.70      0.63
## 5                  NE age               0.01 0.03 0.8282 1.0000     -0.05      0.06
## 6              NE citizen               0.01 0.03 0.7767 1.0000     -0.05      0.07
## 7                NE other               0.00 0.06 0.9512 1.0000     -0.12      0.11
## 8               No answer              -0.36 0.42 0.3927 1.0000     -1.19      0.47
## 9             Other party               0.04 0.01 0.0011 0.0108      0.01      0.06
## 10  Pro-environment party               0.03 0.02 0.2601 1.0000     -0.02      0.08
```

```r
write.csv2(H5.exp.news.mod5.trends.tab,"H5.exp.news.mod5.trends.tab.csv")

#contrast between anti-immigration and pro-environment
(H5.exp.news.contrast<-data.frame(pairs(H5.exp.news.mod5.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast    estimate         SE  df   z.ratio   p.value
## 1 Pro-environment party - Anti-immigration party 0.004055249 0.03015088 Inf 0.1344985 0.8930084
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.news.mod5.trends, "del.eff", by = NULL,adjust=c("none"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  0.06285 0.0636 Inf  0.988  0.3230 
##  Did not vote effect            0.04436 0.0624 Inf  0.711  0.4769 
##  Don't know effect             -0.00368 0.0683 Inf -0.054  0.9570 
##  Invalid vote effect           -0.00252 0.3435 Inf -0.007  0.9941 
##  NE age effect                  0.04285 0.0664 Inf  0.645  0.5190 
##  NE citizen effect              0.04605 0.0680 Inf  0.677  0.4985 
##  NE other effect                0.03219 0.0854 Inf  0.377  0.7061 
##  No answer effect              -0.36484 0.4242 Inf -0.860  0.3897 
##  Other party effect             0.07538 0.0616 Inf  1.224  0.2211 
##  Pro-environment party effect   0.06736 0.0655 Inf  1.028  0.3040 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic
```

```r
#contrast for three voting groups
(H4.more.contrasts<-data.frame(pairs(H5.exp.news.mod5.trends, 
         exclude=c(2:8), by = NULL,adjust=c("none"),reverse=T)))
```

```
##                                         contrast     estimate         SE  df    z.ratio   p.value
## 1           Other party - Anti-immigration party  0.011276213 0.02017026 Inf  0.5590515 0.5761266
## 2 Pro-environment party - Anti-immigration party  0.004055249 0.03015088 Inf  0.1344985 0.8930084
## 3            Pro-environment party - Other party -0.007220963 0.02556620 Inf -0.2824418 0.7776047
```

\newpage

### Model 6: Enter three-way interaction voting group x polnews x environment attitudes


```r
H5.exp.news.mod6<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|voting.group)+
                (environ.lvl1+polnews.lvl1+environ.lvl1:polnews.lvl1|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                environ.lvl1:polnews.lvl1+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2+
                environ.lvl1:polnews.lvl1:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```r
isSingular(H5.exp.news.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod5,H5.exp.news.mod6)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
## H5.exp.news.mod6: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod6:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod6:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod6:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod6:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2 + 
## H5.exp.news.mod6:     environ.lvl1:polnews.lvl1:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.news.mod5   68 103396 103973 -51630   103260                     
## H5.exp.news.mod6   77 103400 104052 -51623   103246 14.446  9     0.1073
```

```r
(FE.H5.exp.news.mod6<-getFE(H5.exp.news.mod6))
```

```
##                                                                Eff   Est   SE       df     t     p    LL    UL
## 1                                                      (Intercept) -0.47 0.14    61.97 -3.23 0.002 -0.76 -0.18
## 2                                                              age  0.00 0.00 35078.87  3.87 0.000  0.00  0.00
## 3                                                           gender  0.06 0.01 35308.31  5.09 0.000  0.04  0.08
## 4                                                             educ  0.01 0.00 35341.53  7.43 0.000  0.01  0.02
## 5                                                            resid -0.06 0.01 35392.15 -5.02 0.000 -0.08 -0.04
## 6                                    occupClerical support workers -0.04 0.09 35278.42 -0.50 0.614 -0.22  0.13
## 7                            occupCraft and related trades workers -0.06 0.09 35284.07 -0.72 0.472 -0.24  0.11
## 8                                      occupElementary occupations  0.01 0.09 35279.64  0.14 0.885 -0.16  0.19
## 9                                                    occupManagers -0.02 0.09 35268.16 -0.22 0.828 -0.19  0.15
## 10                                    occupOther: Not in paid work  0.12 0.09 35335.83  1.35 0.178 -0.06  0.30
## 11                occupPlant and machine operators, and assemblers -0.06 0.09 35281.32 -0.71 0.479 -0.24  0.11
## 12                                              occupProfessionals  0.07 0.09 35277.47  0.79 0.429 -0.10  0.24
## 13                                                    occupRetired -0.04 0.10 35256.78 -0.43 0.670 -0.23  0.15
## 14                                  occupService and sales workers -0.04 0.09 35276.31 -0.46 0.648 -0.21  0.13
## 15         occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35301.00 -0.57 0.567 -0.24  0.13
## 16                    occupTechnicians and associate professionals -0.03 0.09 35272.76 -0.37 0.708 -0.20  0.14
## 17                                                 occupUnemployed -0.03 0.11 35270.48 -0.32 0.746 -0.25  0.18
## 18                                                    environ.lvl1  0.11 0.02   108.17  5.65 0.000  0.07  0.15
## 19                                                    polnews.lvl1  0.03 0.02   143.52  1.39 0.168 -0.01  0.06
## 20                                    all.parties.lvl2Did not vote  0.45 0.06   182.48  7.01 0.000  0.32  0.57
## 21                                      all.parties.lvl2Don't know  0.42 0.07   268.51  5.95 0.000  0.28  0.56
## 22                                    all.parties.lvl2Invalid vote  0.49 0.36  1122.95  1.37 0.172 -0.21  1.19
## 23                                          all.parties.lvl2NE age  0.73 0.07   273.77 10.37 0.000  0.60  0.87
## 24                                      all.parties.lvl2NE citizen  0.86 0.08   268.93 11.22 0.000  0.71  1.01
## 25                                        all.parties.lvl2NE other  0.71 0.10   656.95  7.29 0.000  0.52  0.91
## 26                                       all.parties.lvl2No answer  0.55 0.38  1407.00  1.45 0.149 -0.20  1.29
## 27                                     all.parties.lvl2Other party  0.55 0.05   230.47 11.30 0.000  0.45  0.64
## 28                           all.parties.lvl2Pro-environment party  0.91 0.06   255.65 13.98 0.000  0.78  1.04
## 29                                       environ.lvl1:polnews.lvl1 -0.04 0.02   190.57 -2.51 0.013 -0.07 -0.01
## 30                       environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    84.04  0.00 0.998 -0.04  0.04
## 31                         environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   333.34  0.91 0.366 -0.04  0.10
## 32                       environ.lvl1:all.parties.lvl2Invalid vote  0.11 0.40 25949.48  0.27 0.788 -0.67  0.89
## 33                             environ.lvl1:all.parties.lvl2NE age  0.04 0.03   233.55  1.25 0.214 -0.02  0.10
## 34                         environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   228.84 -1.11 0.268 -0.10  0.03
## 35                           environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1249.46  0.60 0.548 -0.08  0.14
## 36                          environ.lvl1:all.parties.lvl2No answer  0.23 0.23  7977.79  1.00 0.318 -0.23  0.69
## 37                        environ.lvl1:all.parties.lvl2Other party  0.02 0.02   131.21  0.87 0.387 -0.02  0.05
## 38              environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   190.33  0.90 0.367 -0.03  0.08
## 39                       polnews.lvl1:all.parties.lvl2Did not vote -0.02 0.02    93.48 -0.82 0.415 -0.06  0.03
## 40                         polnews.lvl1:all.parties.lvl2Don't know -0.06 0.04   498.33 -1.74 0.083 -0.13  0.01
## 41                       polnews.lvl1:all.parties.lvl2Invalid vote -0.08 0.34 22208.90 -0.23 0.820 -0.75  0.59
## 42                             polnews.lvl1:all.parties.lvl2NE age -0.02 0.03   336.27 -0.61 0.544 -0.08  0.04
## 43                         polnews.lvl1:all.parties.lvl2NE citizen -0.02 0.04   308.09 -0.67 0.504 -0.09  0.05
## 44                           polnews.lvl1:all.parties.lvl2NE other -0.03 0.06  2838.91 -0.40 0.686 -0.15  0.10
## 45                          polnews.lvl1:all.parties.lvl2No answer -0.37 0.43 29076.52 -0.85 0.394 -1.21  0.48
## 46                        polnews.lvl1:all.parties.lvl2Other party  0.01 0.02   161.97  0.37 0.709 -0.03  0.05
## 47              polnews.lvl1:all.parties.lvl2Pro-environment party  0.00 0.03   284.91  0.03 0.977 -0.06  0.06
## 48          environ.lvl1:polnews.lvl1:all.parties.lvl2Did not vote  0.03 0.02   229.81  1.71 0.089 -0.00  0.07
## 49            environ.lvl1:polnews.lvl1:all.parties.lvl2Don't know  0.02 0.03  1291.66  0.78 0.436 -0.04  0.08
## 50          environ.lvl1:polnews.lvl1:all.parties.lvl2Invalid vote  0.45 0.58 34650.10  0.77 0.444 -0.70  1.59
## 51                environ.lvl1:polnews.lvl1:all.parties.lvl2NE age  0.02 0.03   961.62  0.87 0.386 -0.03  0.08
## 52            environ.lvl1:polnews.lvl1:all.parties.lvl2NE citizen  0.09 0.03   824.81  3.22 0.001  0.04  0.15
## 53              environ.lvl1:polnews.lvl1:all.parties.lvl2NE other  0.02 0.05  4486.28  0.37 0.715 -0.08  0.11
## 54             environ.lvl1:polnews.lvl1:all.parties.lvl2No answer -0.04 0.32 31141.51 -0.13 0.900 -0.66  0.58
## 55           environ.lvl1:polnews.lvl1:all.parties.lvl2Other party  0.04 0.02   390.04  2.63 0.009  0.01  0.08
## 56 environ.lvl1:polnews.lvl1:all.parties.lvl2Pro-environment party  0.06 0.02   683.32  2.27 0.023  0.01  0.11
```

```r
(VC.H5.exp.news.mod6<-getVC(H5.exp.news.mod6))
```

```
##             grp                      var1                      var2      est_SD       est_SD2
## 1  voting.group               (Intercept)                      <NA>  0.19777034  3.911311e-02
## 2  voting.group              environ.lvl1                      <NA>  0.04031794  1.625536e-03
## 3  voting.group              polnews.lvl1                      <NA>  0.02912348  8.481772e-04
## 4  voting.group environ.lvl1:polnews.lvl1                      <NA>  0.01253524  1.571323e-04
## 5  voting.group               (Intercept)              environ.lvl1  0.07569809  6.035933e-04
## 6  voting.group               (Intercept)              polnews.lvl1 -0.01046148 -6.025564e-05
## 7  voting.group               (Intercept) environ.lvl1:polnews.lvl1 -0.26477237 -6.563970e-04
## 8  voting.group              environ.lvl1              polnews.lvl1 -0.03841581 -4.510780e-05
## 9  voting.group              environ.lvl1 environ.lvl1:polnews.lvl1 -0.02683966 -1.356463e-05
## 10 voting.group              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.96052025 -3.506570e-04
## 11        cntry               (Intercept)                      <NA>  0.48279116  2.330873e-01
## 12        cntry              environ.lvl1                      <NA>  0.03914727  1.532509e-03
## 13        cntry              polnews.lvl1                      <NA>  0.02711379  7.351575e-04
## 14        cntry environ.lvl1:polnews.lvl1                      <NA>  0.02400230  5.761102e-04
## 15        cntry               (Intercept)              environ.lvl1 -0.02005174 -3.789770e-04
## 16        cntry               (Intercept)              polnews.lvl1  0.45045941  5.896648e-03
## 17        cntry               (Intercept) environ.lvl1:polnews.lvl1 -0.40919557 -4.741798e-03
## 18        cntry              environ.lvl1              polnews.lvl1  0.07246194  7.691334e-05
## 19        cntry              environ.lvl1 environ.lvl1:polnews.lvl1  0.74793761  7.027804e-04
## 20        cntry              polnews.lvl1 environ.lvl1:polnews.lvl1 -0.32870414 -2.139184e-04
## 21     Residual                      <NA>                      <NA>  1.02754700  1.055853e+00
```

#### Refit with manually coded level-1 interaction


```r
dat.H5.news$env.news.int<-dat.H5.news$environ.lvl1*dat.H5.news$polnews.lvl1

H5.exp.news.mod6<-lmer(refugees~
                (environ.lvl1+polnews.lvl1+env.news.int|voting.group)+
                (environ.lvl1+polnews.lvl1+env.news.int|cntry)+
                age+gender+educ+resid+occup+
                environ.lvl1+
                polnews.lvl1+
                env.news.int+
                all.parties.lvl2+
                environ.lvl1:all.parties.lvl2+
                polnews.lvl1:all.parties.lvl2+
                env.news.int:all.parties.lvl2
                ,data=dat.H5.news,REML=F,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e8)))
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```
## boundary (singular) fit: see ?isSingular
```

```
## Warning: Some predictor variables are on very different scales: consider rescaling
```

```r
isSingular(H5.exp.news.mod6)
```

```
## [1] TRUE
```

```r
anova(H5.exp.news.mod5,H5.exp.news.mod6)
```

```
## Data: dat.H5.news
## Models:
## H5.exp.news.mod5: refugees ~ (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     voting.group) + (environ.lvl1 + polnews.lvl1 + environ.lvl1:polnews.lvl1 | 
## H5.exp.news.mod5:     cntry) + age + gender + educ + resid + occup + environ.lvl1 + 
## H5.exp.news.mod5:     polnews.lvl1 + environ.lvl1:polnews.lvl1 + all.parties.lvl2 + 
## H5.exp.news.mod5:     environ.lvl1:all.parties.lvl2 + polnews.lvl1:all.parties.lvl2
## H5.exp.news.mod6: refugees ~ (environ.lvl1 + polnews.lvl1 + env.news.int | voting.group) + 
## H5.exp.news.mod6:     (environ.lvl1 + polnews.lvl1 + env.news.int | cntry) + age + 
## H5.exp.news.mod6:     gender + educ + resid + occup + environ.lvl1 + polnews.lvl1 + 
## H5.exp.news.mod6:     env.news.int + all.parties.lvl2 + environ.lvl1:all.parties.lvl2 + 
## H5.exp.news.mod6:     polnews.lvl1:all.parties.lvl2 + env.news.int:all.parties.lvl2
##                  npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
## H5.exp.news.mod5   68 103396 103973 -51630   103260                     
## H5.exp.news.mod6   77 103400 104052 -51623   103246 14.446  9     0.1073
```

```r
(FE.H5.exp.news.mod6<-getFE(H5.exp.news.mod6))
```

```
##                                                        Eff   Est   SE       df     t     p    LL    UL
## 1                                              (Intercept) -0.47 0.14    61.96 -3.23 0.002 -0.76 -0.18
## 2                                                      age  0.00 0.00 35078.90  3.87 0.000  0.00  0.00
## 3                                                   gender  0.06 0.01 35308.32  5.09 0.000  0.04  0.08
## 4                                                     educ  0.01 0.00 35341.53  7.43 0.000  0.01  0.02
## 5                                                    resid -0.06 0.01 35392.16 -5.02 0.000 -0.08 -0.04
## 6                            occupClerical support workers -0.04 0.09 35278.52 -0.50 0.614 -0.22  0.13
## 7                    occupCraft and related trades workers -0.06 0.09 35284.18 -0.72 0.472 -0.24  0.11
## 8                              occupElementary occupations  0.01 0.09 35279.78  0.14 0.885 -0.16  0.19
## 9                                            occupManagers -0.02 0.09 35268.30 -0.22 0.828 -0.19  0.15
## 10                            occupOther: Not in paid work  0.12 0.09 35335.92  1.35 0.178 -0.06  0.30
## 11        occupPlant and machine operators, and assemblers -0.06 0.09 35281.45 -0.71 0.479 -0.24  0.11
## 12                                      occupProfessionals  0.07 0.09 35277.61  0.79 0.429 -0.10  0.24
## 13                                            occupRetired -0.04 0.10 35256.92 -0.43 0.670 -0.23  0.15
## 14                          occupService and sales workers -0.04 0.09 35276.44 -0.46 0.648 -0.21  0.13
## 15 occupSkilled agricultural, forestry and fishery workers -0.05 0.09 35301.10 -0.57 0.567 -0.24  0.13
## 16            occupTechnicians and associate professionals -0.03 0.09 35272.90 -0.37 0.708 -0.20  0.14
## 17                                         occupUnemployed -0.03 0.11 35270.69 -0.32 0.746 -0.25  0.18
## 18                                            environ.lvl1  0.11 0.02   108.17  5.65 0.000  0.07  0.15
## 19                                            polnews.lvl1  0.03 0.02   143.52  1.39 0.168 -0.01  0.06
## 20                                            env.news.int -0.04 0.02   190.57 -2.51 0.013 -0.07 -0.01
## 21                            all.parties.lvl2Did not vote  0.45 0.06   182.48  7.01 0.000  0.32  0.57
## 22                              all.parties.lvl2Don't know  0.42 0.07   268.51  5.95 0.000  0.28  0.56
## 23                            all.parties.lvl2Invalid vote  0.49 0.36  1122.96  1.37 0.172 -0.21  1.19
## 24                                  all.parties.lvl2NE age  0.73 0.07   273.77 10.37 0.000  0.60  0.87
## 25                              all.parties.lvl2NE citizen  0.86 0.08   268.93 11.22 0.000  0.71  1.01
## 26                                all.parties.lvl2NE other  0.71 0.10   656.95  7.29 0.000  0.52  0.91
## 27                               all.parties.lvl2No answer  0.55 0.38  1407.01  1.45 0.149 -0.20  1.29
## 28                             all.parties.lvl2Other party  0.55 0.05   230.47 11.30 0.000  0.45  0.64
## 29                   all.parties.lvl2Pro-environment party  0.91 0.06   255.65 13.98 0.000  0.78  1.04
## 30               environ.lvl1:all.parties.lvl2Did not vote  0.00 0.02    84.04  0.00 0.998 -0.04  0.04
## 31                 environ.lvl1:all.parties.lvl2Don't know  0.03 0.03   333.35  0.91 0.366 -0.04  0.10
## 32               environ.lvl1:all.parties.lvl2Invalid vote  0.11 0.40 25949.51  0.27 0.788 -0.67  0.89
## 33                     environ.lvl1:all.parties.lvl2NE age  0.04 0.03   233.55  1.25 0.214 -0.02  0.10
## 34                 environ.lvl1:all.parties.lvl2NE citizen -0.04 0.03   228.85 -1.11 0.268 -0.10  0.03
## 35                   environ.lvl1:all.parties.lvl2NE other  0.03 0.06  1249.48  0.60 0.548 -0.08  0.14
## 36                  environ.lvl1:all.parties.lvl2No answer  0.23 0.23  7977.81  1.00 0.318 -0.23  0.69
## 37                environ.lvl1:all.parties.lvl2Other party  0.02 0.02   131.21  0.87 0.387 -0.02  0.05
## 38      environ.lvl1:all.parties.lvl2Pro-environment party  0.03 0.03   190.33  0.90 0.367 -0.03  0.08
## 39               polnews.lvl1:all.parties.lvl2Did not vote -0.02 0.02    93.48 -0.82 0.415 -0.06  0.03
## 40                 polnews.lvl1:all.parties.lvl2Don't know -0.06 0.04   498.33 -1.74 0.083 -0.13  0.01
## 41               polnews.lvl1:all.parties.lvl2Invalid vote -0.08 0.34 22208.83 -0.23 0.820 -0.75  0.59
## 42                     polnews.lvl1:all.parties.lvl2NE age -0.02 0.03   336.27 -0.61 0.544 -0.08  0.04
## 43                 polnews.lvl1:all.parties.lvl2NE citizen -0.02 0.04   308.09 -0.67 0.504 -0.09  0.05
## 44                   polnews.lvl1:all.parties.lvl2NE other -0.03 0.06  2838.89 -0.40 0.686 -0.15  0.10
## 45                  polnews.lvl1:all.parties.lvl2No answer -0.37 0.43 29076.49 -0.85 0.394 -1.21  0.48
## 46                polnews.lvl1:all.parties.lvl2Other party  0.01 0.02   161.97  0.37 0.709 -0.03  0.05
## 47      polnews.lvl1:all.parties.lvl2Pro-environment party  0.00 0.03   284.91  0.03 0.977 -0.06  0.06
## 48               env.news.int:all.parties.lvl2Did not vote  0.03 0.02   229.83  1.71 0.089 -0.00  0.07
## 49                 env.news.int:all.parties.lvl2Don't know  0.02 0.03  1291.76  0.78 0.436 -0.04  0.08
## 50               env.news.int:all.parties.lvl2Invalid vote  0.45 0.58 34650.11  0.77 0.444 -0.70  1.59
## 51                     env.news.int:all.parties.lvl2NE age  0.02 0.03   961.67  0.87 0.386 -0.03  0.08
## 52                 env.news.int:all.parties.lvl2NE citizen  0.09 0.03   824.91  3.22 0.001  0.04  0.15
## 53                   env.news.int:all.parties.lvl2NE other  0.02 0.05  4486.56  0.37 0.715 -0.08  0.11
## 54                  env.news.int:all.parties.lvl2No answer -0.04 0.32 31141.58 -0.13 0.900 -0.66  0.58
## 55                env.news.int:all.parties.lvl2Other party  0.04 0.02   390.06  2.63 0.009  0.01  0.08
## 56      env.news.int:all.parties.lvl2Pro-environment party  0.06 0.02   683.34  2.27 0.023  0.01  0.11
```

```r
(VC.H5.exp.news.mod6<-getVC(H5.exp.news.mod6))
```

```
##             grp         var1         var2      est_SD       est_SD2
## 1  voting.group  (Intercept)         <NA>  0.19777011  3.911302e-02
## 2  voting.group environ.lvl1         <NA>  0.04031781  1.625526e-03
## 3  voting.group polnews.lvl1         <NA>  0.02912369  8.481895e-04
## 4  voting.group env.news.int         <NA>  0.01253505  1.571275e-04
## 5  voting.group  (Intercept) environ.lvl1  0.07569780  6.035883e-04
## 6  voting.group  (Intercept) polnews.lvl1 -0.01046183 -6.025800e-05
## 7  voting.group  (Intercept) env.news.int -0.26477408 -6.563903e-04
## 8  voting.group environ.lvl1 polnews.lvl1 -0.03840948 -4.510055e-05
## 9  voting.group environ.lvl1 env.news.int -0.02684509 -1.356712e-05
## 10 voting.group polnews.lvl1 env.news.int -0.96051970 -3.506540e-04
## 11        cntry  (Intercept)         <NA>  0.48279638  2.330923e-01
## 12        cntry environ.lvl1         <NA>  0.03914714  1.532499e-03
## 13        cntry polnews.lvl1         <NA>  0.02711396  7.351667e-04
## 14        cntry env.news.int         <NA>  0.02400230  5.761106e-04
## 15        cntry  (Intercept) environ.lvl1 -0.02005803 -3.790988e-04
## 16        cntry  (Intercept) polnews.lvl1  0.45046839  5.896866e-03
## 17        cntry  (Intercept) env.news.int -0.40921099 -4.742029e-03
## 18        cntry environ.lvl1 polnews.lvl1  0.07245464  7.690582e-05
## 19        cntry environ.lvl1 env.news.int  0.74793891  7.027795e-04
## 20        cntry polnews.lvl1 env.news.int -0.32871318 -2.139257e-04
## 21     Residual         <NA>         <NA>  1.02754701  1.055853e+00
```


#### Marginal effect for pro-environment and anti-immigration voters


```r
H5.exp.news.mod6.trends<-emtrends(H5.exp.news.mod6,specs = c("all.parties.lvl2"),var=c("env.news.int"))
(H5.exp.news.mod6.trends.tab<-data.frame(H5.exp.news.mod6.trends))
```

```
##          all.parties.lvl2 env.news.int.trend          SE  df    asymp.LCL    asymp.UCL
## 1  Anti-immigration party       -0.039341965 0.015666246 Inf -0.070047243 -0.008636688
## 2            Did not vote       -0.008904728 0.011350070 Inf -0.031150457  0.013341001
## 3              Don't know       -0.015395924 0.027527224 Inf -0.069348292  0.038556443
## 4            Invalid vote        0.408057711 0.584542696 Inf -0.737624920  1.553740341
## 5                  NE age       -0.016160217 0.022944451 Inf -0.061130515  0.028810082
## 6              NE citizen        0.052434017 0.024900760 Inf  0.003629425  0.101238610
## 7                NE other       -0.021775930 0.046157538 Inf -0.112243041  0.068691181
## 8               No answer       -0.078891894 0.314912820 Inf -0.696109679  0.538325891
## 9             Other party        0.003768864 0.008846732 Inf -0.013570411  0.021108139
## 10  Pro-environment party        0.017392344 0.020815457 Inf -0.023405203  0.058189890
```

```r
H5.exp.news.mod6.trends.tab$p<-
  2*(1-pnorm(abs(H5.exp.news.mod6.trends.tab$env.news.int.trend/
                   H5.exp.news.mod6.trends.tab$SE)))
H5.exp.news.mod6.trends.tab$adj.p<-
  p.adjust(H5.exp.news.mod6.trends.tab$p,method="holm")

H5.exp.news.mod6.trends.tab<-
  cbind(group=H5.exp.news.mod6.trends.tab[,1],
      round(H5.exp.news.mod6.trends.tab[,c(2,3)],2),
      round(H5.exp.news.mod6.trends.tab[,c(7,8)],4),
      round(H5.exp.news.mod6.trends.tab[,c(5,6)],2))
H5.exp.news.mod6.trends.tab
```

```
##                     group env.news.int.trend   SE      p  adj.p asymp.LCL asymp.UCL
## 1  Anti-immigration party              -0.04 0.02 0.0120 0.1203     -0.07     -0.01
## 2            Did not vote              -0.01 0.01 0.4327 1.0000     -0.03      0.01
## 3              Don't know              -0.02 0.03 0.5760 1.0000     -0.07      0.04
## 4            Invalid vote               0.41 0.58 0.4851 1.0000     -0.74      1.55
## 5                  NE age              -0.02 0.02 0.4812 1.0000     -0.06      0.03
## 6              NE citizen               0.05 0.02 0.0352 0.3171      0.00      0.10
## 7                NE other              -0.02 0.05 0.6371 1.0000     -0.11      0.07
## 8               No answer              -0.08 0.31 0.8022 1.0000     -0.70      0.54
## 9             Other party               0.00 0.01 0.6701 1.0000     -0.01      0.02
## 10  Pro-environment party               0.02 0.02 0.4034 1.0000     -0.02      0.06
```

```r
write.csv2(H5.exp.news.mod6.trends.tab,"H5.exp.news.mod6.trends.tab.csv")



#contrast between anti-immigration and pro-environment
(H5.exp.news.contrast<-data.frame(pairs(H5.exp.news.mod6.trends, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast   estimate        SE  df  z.ratio    p.value
## 1 Pro-environment party - Anti-immigration party 0.05673431 0.0249874 Inf 2.270516 0.02317627
```

```r
#contrast for all groups against mean of other groups
contrast(H5.exp.news.mod6.trends, "del.eff", by = NULL,adjust=c("holm"))
```

```
##  contrast                      estimate     SE  df z.ratio p.value
##  Anti-immigration party effect  -0.0772 0.0756 Inf -1.021  1.0000 
##  Did not vote effect            -0.0434 0.0748 Inf -0.579  1.0000 
##  Don't know effect              -0.0506 0.0789 Inf -0.641  1.0000 
##  Invalid vote effect             0.4199 0.5856 Inf  0.717  1.0000 
##  NE age effect                  -0.0514 0.0774 Inf -0.664  1.0000 
##  NE citizen effect               0.0248 0.0780 Inf  0.318  1.0000 
##  NE other effect                -0.0577 0.0870 Inf -0.662  1.0000 
##  No answer effect               -0.1211 0.3216 Inf -0.377  1.0000 
##  Other party effect             -0.0293 0.0745 Inf -0.393  1.0000 
##  Pro-environment party effect   -0.0141 0.0768 Inf -0.184  1.0000 
## 
## Results are averaged over the levels of: gender, resid, occup 
## Degrees-of-freedom method: asymptotic 
## P value adjustment: holm method for 10 tests
```

```r
#contrast for three voting groups
(H5.exp.news.more.contrasts<-data.frame(pairs(H5.exp.news.mod6.trends, 
         exclude=c(2:8), by = NULL,adjust=c("holm"),reverse=T)))
```

```
##                                         contrast   estimate         SE  df   z.ratio    p.value
## 1           Other party - Anti-immigration party 0.04311083 0.01637828 Inf 2.6321952 0.02545053
## 2 Pro-environment party - Anti-immigration party 0.05673431 0.02498740 Inf 2.2705164 0.04635253
## 3            Pro-environment party - Other party 0.01362348 0.02129393 Inf 0.6397825 0.52231404
```

