---
title: "Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_depth: 5
---







\newpage

# Preparations

* This "preparations" part of script is same for all different files





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
## [1] LC_COLLATE=Finnish_Finland.1252  LC_CTYPE=Finnish_Finland.1252    LC_MONETARY=Finnish_Finland.1252
## [4] LC_NUMERIC=C                     LC_TIME=Finnish_Finland.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] merTools_0.5.0  arm_1.10-1      MASS_7.3-51.5   metafor_2.4-0   ggplot2_3.3.0   emmeans_1.4.5   psych_1.9.12.31
##  [8] dplyr_0.8.5     lmerTest_3.1-2  lme4_1.1-23     Matrix_1.2-18  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6        mvtnorm_1.1-0       lattice_0.20-38     tidyr_1.0.2         zoo_1.8-7           foreach_1.5.0      
##  [7] assertthat_0.2.1    digest_0.6.25       mime_0.9            R6_2.4.1            backports_1.1.6     evaluate_0.14      
## [13] coda_0.19-3         pillar_1.4.3        rlang_0.4.5         multcomp_1.4-13     minqa_1.2.4         nloptr_1.2.2.1     
## [19] rmarkdown_2.1       splines_3.6.3       statmod_1.4.34      stringr_1.4.0       munsell_0.5.0       shiny_1.4.0.2      
## [25] broom_0.5.5         httpuv_1.5.2        compiler_3.6.3      numDeriv_2016.8-1.1 xfun_0.13           pkgconfig_2.0.3    
## [31] mnormt_1.5-6        htmltools_0.4.0     tidyselect_1.0.0    tibble_3.0.0        codetools_0.2-16    fansi_0.4.1        
## [37] later_1.0.0         crayon_1.3.4        withr_2.1.2         grid_3.6.3          nlme_3.1-144        xtable_1.8-4       
## [43] gtable_0.3.0        lifecycle_0.2.0     magrittr_1.5        scales_1.1.0        estimability_1.3    cli_2.0.2          
## [49] stringi_1.4.6       promises_1.1.0      generics_0.0.2      ellipsis_0.3.0      vctrs_0.2.4         boot_1.3-24        
## [55] sandwich_2.5-1      blme_1.0-4          TH.data_1.0-10      iterators_1.0.12    tools_3.6.3         glue_1.4.0         
## [61] purrr_0.3.3         fastmap_1.0.1       abind_1.4-5         parallel_3.6.3      survival_3.1-8      yaml_2.2.0         
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
##   AT   BE   CH   CZ   DE   EE   ES   FI   FR   GB   HU   IE   IT   LT   NL   NO   PL   PT   SE   SI 
## 1973 1753 1503 2156 2819 1974 1817 1862 2015 1876 1391 2676 2317 1927 1661 1538 1589 1228 1525 1276
```

#### Voting group


```r
#make voting group variable names unique to each country
dat$voting.group<-paste0(dat$cntry,": ",dat$vote.group.combined)
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

# Hypothesis 1: There will be a positive association between pro-environment and pro-refugee attitudes

### Model 0: Intercepts only


```r
H1.mod0<-lmer(refugees~(1|voting.group)+(1|cntry),
              data=dat,REML=F)

(FE.H1.mod0<-getFE(H1.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL   UL
## (Intercept)     0.08       0.11 19.92    0.69 0.499 -0.16 0.32
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                 0.06       0.14    50.42    0.44 0.663 -0.22  0.35
## age                                                         0.00       0.00 34095.35    3.30 0.001  0.00  0.00
## gender                                                      0.06       0.01 35598.09    4.77 0.000  0.03  0.08
## educ                                                        0.02       0.00 35706.37    9.28 0.000  0.01  0.02
## resid                                                      -0.08       0.01 35657.25   -6.58 0.000 -0.10 -0.05
## occupClerical support workers                              -0.04       0.09 35532.93   -0.40 0.687 -0.21  0.14
## occupCraft and related trades workers                      -0.08       0.09 35539.44   -0.85 0.398 -0.25  0.10
## occupElementary occupations                                 0.01       0.09 35540.85    0.14 0.892 -0.16  0.19
## occupManagers                                               0.00       0.09 35533.47    0.00 0.998 -0.18  0.18
## occupOther: Not in paid work                                0.16       0.09 35696.52    1.71 0.088 -0.02  0.34
## occupPlant and machine operators, and assemblers           -0.07       0.09 35539.06   -0.77 0.442 -0.24  0.11
## occupProfessionals                                          0.10       0.09 35535.19    1.10 0.273 -0.08  0.27
## occupRetired                                               -0.03       0.10 35538.37   -0.28 0.782 -0.22  0.17
## occupService and sales workers                             -0.04       0.09 35534.90   -0.40 0.691 -0.21  0.14
## occupSkilled agricultural, forestry and fishery workers    -0.05       0.09 35544.87   -0.52 0.601 -0.23  0.14
## occupTechnicians and associate professionals               -0.03       0.09 35532.10   -0.29 0.774 -0.20  0.15
## occupUnemployed                                            -0.03       0.11 35546.28   -0.23 0.816 -0.24  0.19
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                 0.07       0.14    49.59    0.51 0.611 -0.21  0.36
## age                                                         0.00       0.00 34270.64    4.39 0.000  0.00  0.00
## gender                                                      0.06       0.01 35593.07    4.67 0.000  0.03  0.08
## educ                                                        0.01       0.00 35714.29    7.51 0.000  0.01  0.02
## resid                                                      -0.06       0.01 35651.76   -5.23 0.000 -0.08 -0.04
## occupClerical support workers                              -0.04       0.09 35528.94   -0.46 0.646 -0.21  0.13
## occupCraft and related trades workers                      -0.07       0.09 35535.31   -0.76 0.446 -0.24  0.11
## occupElementary occupations                                 0.01       0.09 35536.62    0.16 0.871 -0.16  0.19
## occupManagers                                              -0.02       0.09 35529.47   -0.18 0.859 -0.19  0.16
## occupOther: Not in paid work                                0.14       0.09 35691.43    1.54 0.123 -0.04  0.32
## occupPlant and machine operators, and assemblers           -0.06       0.09 35534.89   -0.72 0.474 -0.24  0.11
## occupProfessionals                                          0.07       0.09 35531.19    0.83 0.407 -0.10  0.24
## occupRetired                                               -0.03       0.10 35534.45   -0.34 0.734 -0.23  0.16
## occupService and sales workers                             -0.04       0.09 35530.91   -0.43 0.668 -0.21  0.13
## occupSkilled agricultural, forestry and fishery workers    -0.05       0.09 35540.59   -0.58 0.564 -0.24  0.13
## occupTechnicians and associate professionals               -0.03       0.09 35528.10   -0.36 0.718 -0.20  0.14
## occupUnemployed                                            -0.02       0.11 35542.75   -0.19 0.852 -0.23  0.19
## environ.lvl1                                                0.13       0.00 35459.04   26.76 0.000  0.12  0.13
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                 0.08       0.14    49.47    0.53 0.598 -0.21  0.36
## age                                                         0.00       0.00 34211.38    4.49 0.000  0.00  0.00
## gender                                                      0.06       0.01 35569.60    4.70 0.000  0.03  0.08
## educ                                                        0.01       0.00 35698.92    7.49 0.000  0.01  0.02
## resid                                                      -0.06       0.01 35634.38   -5.24 0.000 -0.08 -0.04
## occupClerical support workers                              -0.04       0.09 35508.21   -0.48 0.632 -0.22  0.13
## occupCraft and related trades workers                      -0.07       0.09 35518.38   -0.77 0.439 -0.24  0.10
## occupElementary occupations                                 0.01       0.09 35520.28    0.13 0.900 -0.16  0.18
## occupManagers                                              -0.02       0.09 35507.43   -0.22 0.827 -0.19  0.15
## occupOther: Not in paid work                                0.14       0.09 35674.75    1.52 0.127 -0.04  0.32
## occupPlant and machine operators, and assemblers           -0.07       0.09 35516.07   -0.74 0.460 -0.24  0.11
## occupProfessionals                                          0.07       0.09 35512.07    0.79 0.432 -0.10  0.24
## occupRetired                                               -0.04       0.10 35508.44   -0.36 0.716 -0.23  0.16
## occupService and sales workers                             -0.04       0.09 35513.92   -0.44 0.657 -0.21  0.13
## occupSkilled agricultural, forestry and fishery workers    -0.06       0.09 35523.04   -0.63 0.529 -0.24  0.12
## occupTechnicians and associate professionals               -0.03       0.09 35508.10   -0.38 0.705 -0.20  0.14
## occupUnemployed                                            -0.03       0.11 35525.32   -0.24 0.812 -0.24  0.18
## environ.lvl1                                                0.12       0.01    19.40   11.53 0.000  0.10  0.15
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
country.effs.dat
```

```
##    cntry slope_with_covs slope_without_covs observed_r  partial_r    n
## 1     AT      0.16388921         0.17073407 0.25596053 0.22287125 1934
## 2     BE      0.13561055         0.14776506 0.17729252 0.14631334 1746
## 3     CH      0.12415313         0.13494604 0.21389146 0.17601223 1469
## 4     CZ      0.13937493         0.14235903 0.17395316 0.16591222 2071
## 5     DE      0.18286129         0.19258682 0.24433787 0.23353039 2791
## 6     EE      0.07728714         0.08881514 0.05001522 0.04719782 1927
## 7     ES      0.10525890         0.11563744 0.15432815 0.13903233 1666
## 8     FI      0.15671522         0.16971104 0.23493567 0.22299268 1840
## 9     FR      0.17072139         0.17996195 0.20552407 0.19741954 1995
## 10    GB      0.14359012         0.15307094 0.20728022 0.18157755 1856
## 11    HU      0.12194688         0.12843216 0.17372436 0.17371641 1341
## 12    IE      0.14056993         0.15218857 0.19521566 0.17388402 2591
## 13    IT      0.11174664         0.11988893 0.14321699 0.13573951 2177
## 14    LT      0.08476625         0.09240974 0.09612563 0.09510338 1759
## 15    NL      0.09383763         0.10475537 0.15260881 0.14788809 1620
## 16    NO      0.13225738         0.14956688 0.26483270 0.22781088 1530
## 17    PL      0.05053810         0.06035134 0.02603599 0.02285898 1495
## 18    PT      0.07681183         0.09485497 0.11228764 0.08639693 1214
## 19    SE      0.13932739         0.15254652 0.25416405 0.22308459 1487
## 20    SI      0.12051011         0.13280999 0.17475936 0.15300594 1231
```

```r
write.csv2(country.effs.dat,"associations_within_countries.csv")

#calculate the meta-analytical estimate for partial correlation
country.effs.dat$partial_r.se<-1/sqrt(country.effs.dat$n-3)
rma.uni(yi=transf.rtoz(country.effs.dat$partial_r),
        sei=country.effs.dat$partial_r.se)
```

```
## 
## Random-Effects Model (k = 20; tau^2 estimator: REML)
## 
## tau^2 (estimated amount of total heterogeneity): 0.0031 (SE = 0.0012)
## tau (square root of estimated tau^2 value):      0.0555
## I^2 (total heterogeneity / total variability):   84.59%
## H^2 (total variability / sampling variability):  6.49
## 
## Test for Heterogeneity:
## Q(df = 19) = 122.1403, p-val < .0001
## 
## Model Results:
## 
## estimate      se     zval    pval   ci.lb   ci.ub 
##   0.1610  0.0135  11.8871  <.0001  0.1344  0.1875  *** 
## 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#back-transform to r
transf.ztor(0.1610)
```

```
## [1] 0.1596232
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

## Alternative (exploratory) approach for Hypothesis 1 with Environment attitudes as dependent variable, and immigrant attitudes as independent

### Model 0: Intercepts only


```r
H1.env.mod0<-lmer(environ.gmc~(1|voting.group)+(1|cntry),data=dat,REML=F)

(FE.H1.env.mod0<-getFE(H1.env.mod0))
```

```
##             Estimate Std..Error    df t.value     p    LL  UL
## (Intercept)     0.05       0.07 20.08    0.78 0.443 -0.09 0.2
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                -0.04       0.12   188.50   -0.32 0.751 -0.27  0.20
## age                                                         0.00       0.00 29829.87   -8.11 0.000  0.00  0.00
## gender                                                      0.02       0.01 35659.58    1.14 0.253 -0.01  0.04
## educ                                                        0.03       0.00 35342.84   13.39 0.000  0.02  0.03
## resid                                                      -0.14       0.01 35731.67  -10.19 0.000 -0.16 -0.11
## occupClerical support workers                               0.04       0.10 35578.76    0.42 0.672 -0.15  0.24
## occupCraft and related trades workers                      -0.06       0.10 35588.19   -0.60 0.548 -0.26  0.14
## occupElementary occupations                                -0.01       0.10 35591.13   -0.11 0.915 -0.21  0.19
## occupManagers                                               0.13       0.10 35579.35    1.33 0.184 -0.06  0.33
## occupOther: Not in paid work                                0.16       0.10 35719.64    1.55 0.121 -0.04  0.36
## occupPlant and machine operators, and assemblers           -0.04       0.10 35588.53   -0.41 0.684 -0.24  0.16
## occupProfessionals                                          0.21       0.10 35581.91    2.09 0.037  0.01  0.40
## occupRetired                                                0.06       0.11 35582.94    0.52 0.604 -0.16  0.28
## occupService and sales workers                              0.02       0.10 35580.41    0.25 0.806 -0.17  0.22
## occupSkilled agricultural, forestry and fishery workers     0.05       0.11 35595.50    0.43 0.664 -0.16  0.25
## occupTechnicians and associate professionals                0.06       0.10 35577.81    0.59 0.557 -0.14  0.25
## occupUnemployed                                            -0.03       0.12 35586.90   -0.25 0.803 -0.27  0.21
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                -0.04       0.12   182.49   -0.30 0.763 -0.27  0.20
## age                                                         0.00       0.00 30300.28   -8.71 0.000  0.00  0.00
## gender                                                      0.01       0.01 35654.32    0.49 0.622 -0.02  0.03
## educ                                                        0.03       0.00 35392.50   12.22 0.000  0.02  0.03
## resid                                                      -0.12       0.01 35728.36   -9.40 0.000 -0.15 -0.10
## occupClerical support workers                               0.05       0.10 35574.16    0.48 0.629 -0.15  0.24
## occupCraft and related trades workers                      -0.05       0.10 35583.40   -0.49 0.625 -0.24  0.15
## occupElementary occupations                                -0.01       0.10 35586.20   -0.12 0.903 -0.21  0.18
## occupManagers                                               0.13       0.10 35574.73    1.34 0.181 -0.06  0.33
## occupOther: Not in paid work                                0.14       0.10 35722.40    1.35 0.178 -0.06  0.34
## occupPlant and machine operators, and assemblers           -0.03       0.10 35583.69   -0.30 0.762 -0.23  0.17
## occupProfessionals                                          0.19       0.10 35577.21    1.96 0.050  0.00  0.38
## occupRetired                                                0.06       0.11 35578.52    0.56 0.573 -0.15  0.28
## occupService and sales workers                              0.03       0.10 35575.86    0.31 0.760 -0.16  0.22
## occupSkilled agricultural, forestry and fishery workers     0.05       0.10 35590.53    0.51 0.610 -0.15  0.26
## occupTechnicians and associate professionals                0.06       0.10 35573.19    0.63 0.527 -0.13  0.25
## occupUnemployed                                            -0.02       0.12 35582.73   -0.20 0.838 -0.26  0.21
## refugees.lvl1                                               0.16       0.01 35453.62   26.78 0.000  0.15  0.17
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
##                                                         Estimate Std..Error       df t.value     p    LL    UL
## (Intercept)                                                -0.04       0.12   181.41   -0.33 0.740 -0.27  0.19
## age                                                         0.00       0.00 30229.24   -8.64 0.000  0.00  0.00
## gender                                                      0.01       0.01 35647.27    0.51 0.613 -0.02  0.03
## educ                                                        0.03       0.00 35243.03   12.16 0.000  0.02  0.03
## resid                                                      -0.12       0.01 35688.34   -9.42 0.000 -0.15 -0.10
## occupClerical support workers                               0.05       0.10 35558.96    0.52 0.604 -0.14  0.25
## occupCraft and related trades workers                      -0.04       0.10 35566.00   -0.45 0.652 -0.24  0.15
## occupElementary occupations                                -0.01       0.10 35571.70   -0.10 0.921 -0.20  0.19
## occupManagers                                               0.14       0.10 35560.44    1.37 0.170 -0.06  0.33
## occupOther: Not in paid work                                0.14       0.10 35708.21    1.40 0.163 -0.06  0.34
## occupPlant and machine operators, and assemblers           -0.03       0.10 35566.12   -0.27 0.788 -0.22  0.17
## occupProfessionals                                          0.19       0.10 35558.85    1.99 0.047  0.00  0.39
## occupRetired                                                0.06       0.11 35554.14    0.58 0.560 -0.15  0.28
## occupService and sales workers                              0.03       0.10 35560.22    0.34 0.735 -0.16  0.23
## occupSkilled agricultural, forestry and fishery workers     0.05       0.10 35572.65    0.50 0.615 -0.15  0.26
## occupTechnicians and associate professionals                0.07       0.10 35559.38    0.69 0.489 -0.12  0.26
## occupUnemployed                                            -0.02       0.12 35571.95   -0.19 0.848 -0.26  0.21
## refugees.lvl1                                               0.16       0.01    18.18   13.89 0.000  0.13  0.18
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

