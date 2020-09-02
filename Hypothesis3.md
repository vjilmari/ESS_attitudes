---
title: "Hypothesis 3"
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.04 0.12   188.50
## 2                                                      age -0.00 0.00 29829.87
## 3                                                   gender  0.02 0.01 35659.58
## 4                                                     educ  0.03 0.00 35342.84
## 5                                                    resid -0.14 0.01 35731.67
## 6                            occupClerical support workers  0.04 0.10 35578.76
## 7                    occupCraft and related trades workers -0.06 0.10 35588.19
## 8                              occupElementary occupations -0.01 0.10 35591.13
## 9                                            occupManagers  0.13 0.10 35579.35
## 10                            occupOther: Not in paid work  0.16 0.10 35719.64
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35588.53
## 12                                      occupProfessionals  0.21 0.10 35581.91
## 13                                            occupRetired  0.06 0.11 35582.94
## 14                          occupService and sales workers  0.02 0.10 35580.41
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35595.50
## 16            occupTechnicians and associate professionals  0.06 0.10 35577.81
## 17                                         occupUnemployed -0.03 0.12 35586.90
##         t     p    LL    UL
## 1   -0.32 0.751 -0.27  0.20
## 2   -8.11 0.000 -0.00 -0.00
## 3    1.14 0.253 -0.01  0.04
## 4   13.39 0.000  0.02  0.03
## 5  -10.19 0.000 -0.16 -0.11
## 6    0.42 0.672 -0.15  0.24
## 7   -0.60 0.548 -0.26  0.14
## 8   -0.11 0.915 -0.21  0.19
## 9    1.33 0.184 -0.06  0.33
## 10   1.55 0.121 -0.04  0.36
## 11  -0.41 0.684 -0.24  0.16
## 12   2.09 0.037  0.01  0.40
## 13   0.52 0.604 -0.16  0.28
## 14   0.25 0.806 -0.17  0.22
## 15   0.43 0.664 -0.16  0.25
## 16   0.59 0.557 -0.14  0.25
## 17  -0.25 0.803 -0.27  0.21
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.27 0.12   280.06
## 2                                                      age -0.00 0.00 35583.97
## 3                                                   gender  0.02 0.01 35672.30
## 4                                                     educ  0.03 0.00 35720.39
## 5                                                    resid -0.13 0.01 35734.98
## 6                            occupClerical support workers  0.04 0.10 35641.21
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39
## 8                              occupElementary occupations -0.01 0.10 35645.09
## 9                                            occupManagers  0.13 0.10 35642.84
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76
## 12                                      occupProfessionals  0.20 0.10 35641.83
## 13                                            occupRetired  0.05 0.11 35635.37
## 14                          occupService and sales workers  0.02 0.10 35640.53
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21
## 17                                         occupUnemployed -0.04 0.12 35636.70
## 18                            all.parties.lvl2Did not vote  0.10 0.05   151.88
## 19                              all.parties.lvl2Don't know  0.08 0.06   275.08
## 20                            all.parties.lvl2Invalid vote -0.62 0.36  2074.65
## 21                                  all.parties.lvl2NE age  0.39 0.06   286.34
## 22                              all.parties.lvl2NE citizen  0.24 0.07   254.61
## 23                                all.parties.lvl2NE other  0.30 0.10   835.68
## 24                               all.parties.lvl2No answer  0.20 0.39  2698.12
## 25                             all.parties.lvl2Other party  0.22 0.04   204.52
## 26                   all.parties.lvl2Pro-environment party  0.70 0.06   241.69
##         t     p    LL    UL
## 1   -2.20 0.029 -0.50 -0.03
## 2   -7.51 0.000 -0.00 -0.00
## 3    1.19 0.235 -0.01  0.04
## 4   13.27 0.000  0.02  0.03
## 5  -10.05 0.000 -0.16 -0.11
## 6    0.40 0.689 -0.16  0.24
## 7   -0.60 0.548 -0.26  0.14
## 8   -0.13 0.900 -0.21  0.18
## 9    1.34 0.181 -0.06  0.33
## 10   1.32 0.188 -0.07  0.34
## 11  -0.41 0.682 -0.24  0.16
## 12   2.06 0.039  0.01  0.40
## 13   0.49 0.624 -0.16  0.27
## 14   0.24 0.812 -0.17  0.22
## 15   0.45 0.656 -0.16  0.25
## 16   0.58 0.564 -0.14  0.25
## 17  -0.30 0.765 -0.27  0.20
## 18   1.82 0.070 -0.01  0.20
## 19   1.18 0.240 -0.05  0.20
## 20  -1.69 0.092 -1.33  0.10
## 21   6.13 0.000  0.27  0.52
## 22   3.48 0.001  0.10  0.37
## 23   3.08 0.002  0.11  0.49
## 24   0.52 0.603 -0.56  0.97
## 25   5.09 0.000  0.13  0.30
## 26  12.15 0.000  0.59  0.82
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
##                     group emmean             CI     p p.adj p_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes
##    p.adj_less_001
## 1              no
## 2              no
## 3              no
## 4              no
## 5              no
## 6              no
## 7              no
## 8              no
## 9              no
## 10            yes
```

```r
write.csv2(H3.mod2.mmeans.tab,"H3.mod2.mmeans.tab.csv")

#contrast between anti-immigration and pro-environment
(H3.contrast<-data.frame(pairs(H3.mod2.mmeans, exclude = c(2:9),reverse=T)))
```

```
##                                         contrast  estimate         SE  df
## 1 Pro-environment party - Anti-immigration party 0.7023866 0.05780773 Inf
##    z.ratio      p.value
## 1 12.15039 5.709007e-34
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
##                                         contrast  estimate         SE  df
## 1           Other party - Anti-immigration party 0.2151125 0.04227427 Inf
## 2 Pro-environment party - Anti-immigration party 0.7023866 0.05780773 Inf
## 3            Pro-environment party - Other party 0.4872741 0.04683923 Inf
##     z.ratio      p.value
## 1  5.088498 3.609110e-07
## 2 12.150394 1.712702e-33
## 3 10.403119 4.798968e-25
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
## Warning in Ops.factor(H3.mod2.mmeans.tab[10, 2], H3.mod2.mmeans.tab[1, 2]): '-'
## not meaningful for factors
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
##                     group emmean             CI     p p.adj p_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes
##    p.adj_less_001
## 1              no
## 2              no
## 3              no
## 4              no
## 5              no
## 6              no
## 7              no
## 8              no
## 9              no
## 10            yes
```

```r
(H3.effect.size.env.other<-(H3.mod2.mmeans.tab[10,2]-
  H3.mod2.mmeans.tab[9,2])/H3.pooled.sd.other.env)
```

```
## Warning in Ops.factor(H3.mod2.mmeans.tab[10, 2], H3.mod2.mmeans.tab[9, 2]): '-'
## not meaningful for factors
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
##                     group emmean             CI     p p.adj p_less_001
## 1  Anti-immigration party  -0.22 [-0.36, -0.08] 0.002 0.020         no
## 2            Did not vote  -0.12  [-0.27, 0.02] 0.086 0.430         no
## 3              Don't know  -0.15  [-0.30, 0.01] 0.066 0.398         no
## 4            Invalid vote  -0.84 [-1.56, -0.12] 0.023 0.183         no
## 5                  NE age   0.17   [0.01, 0.33] 0.034 0.241         no
## 6              NE citizen   0.02  [-0.15, 0.18] 0.850 1.000         no
## 7                NE other   0.07  [-0.14, 0.29] 0.489 1.000         no
## 8               No answer  -0.02  [-0.79, 0.75] 0.961 1.000         no
## 9             Other party  -0.01  [-0.13, 0.12] 0.911 1.000         no
## 10  Pro-environment party   0.48   [0.33, 0.63] 0.000 0.000        yes
##    p.adj_less_001
## 1              no
## 2              no
## 3              no
## 4              no
## 5              no
## 6              no
## 7              no
## 8              no
## 9              no
## 10            yes
```

```r
(H3.effect.size.imm.other<-(H3.mod2.mmeans.tab[9,2]-
  H3.mod2.mmeans.tab[1,2])/H3.pooled.sd.other.imm)
```

```
## Warning in Ops.factor(H3.mod2.mmeans.tab[9, 2], H3.mod2.mmeans.tab[1, 2]): '-'
## not meaningful for factors
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.17 0.12   278.89
## 2                                                      age -0.00 0.00 35583.97
## 3                                                   gender  0.02 0.01 35672.29
## 4                                                     educ  0.03 0.00 35720.39
## 5                                                    resid -0.13 0.01 35734.98
## 6                            occupClerical support workers  0.04 0.10 35641.21
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39
## 8                              occupElementary occupations -0.01 0.10 35645.08
## 9                                            occupManagers  0.13 0.10 35642.84
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76
## 12                                      occupProfessionals  0.20 0.10 35641.83
## 13                                            occupRetired  0.05 0.11 35635.36
## 14                          occupService and sales workers  0.02 0.10 35640.53
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21
## 17                                         occupUnemployed -0.04 0.12 35636.70
## 18                                       other.party.dummy  0.12 0.04   132.01
## 19                                         dont.know.dummy -0.02 0.06   216.96
## 20                                      invalid.vote.dummy -0.71 0.36  2042.66
## 21                                         no.answer.dummy  0.10 0.39  2638.69
## 22                                  not.eligible.age.dummy  0.29 0.06   214.47
## 23                          not.eligible.citizenship.dummy  0.14 0.07   206.74
## 24                                not.eligible.other.dummy  0.20 0.10   697.82
## 25                                    anti.imm.party.dummy -0.10 0.05   151.88
## 26                                     pro.env.party.dummy  0.60 0.06   185.57
##         t     p    LL    UL
## 1   -1.39 0.167 -0.41  0.07
## 2   -7.51 0.000 -0.00 -0.00
## 3    1.19 0.235 -0.01  0.04
## 4   13.27 0.000  0.02  0.03
## 5  -10.05 0.000 -0.16 -0.11
## 6    0.40 0.689 -0.16  0.24
## 7   -0.60 0.548 -0.26  0.14
## 8   -0.13 0.900 -0.21  0.18
## 9    1.34 0.181 -0.06  0.33
## 10   1.32 0.188 -0.07  0.34
## 11  -0.41 0.682 -0.24  0.16
## 12   2.06 0.039  0.01  0.40
## 13   0.49 0.624 -0.16  0.27
## 14   0.24 0.812 -0.17  0.22
## 15   0.45 0.656 -0.16  0.25
## 16   0.58 0.564 -0.14  0.25
## 17  -0.30 0.765 -0.27  0.20
## 18   2.79 0.006  0.03  0.20
## 19  -0.36 0.720 -0.15  0.10
## 20  -1.96 0.051 -1.43  0.00
## 21   0.27 0.788 -0.66  0.87
## 22   4.65 0.000  0.17  0.42
## 23   2.06 0.041  0.01  0.27
## 24   2.06 0.039  0.01  0.39
## 25  -1.82 0.070 -0.20  0.01
## 26  10.47 0.000  0.49  0.72
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.17 0.12   286.87
## 2                                                      age -0.00 0.00 35536.29
## 3                                                   gender  0.02 0.01 35674.84
## 4                                                     educ  0.03 0.00 35701.59
## 5                                                    resid -0.13 0.01 35722.10
## 6                            occupClerical support workers  0.04 0.10 35644.30
## 7                    occupCraft and related trades workers -0.06 0.10 35652.54
## 8                              occupElementary occupations -0.01 0.10 35648.15
## 9                                            occupManagers  0.14 0.10 35646.01
## 10                            occupOther: Not in paid work  0.14 0.10 35698.11
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35653.05
## 12                                      occupProfessionals  0.20 0.10 35644.89
## 13                                            occupRetired  0.06 0.11 35638.46
## 14                          occupService and sales workers  0.02 0.10 35643.56
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35663.89
## 16            occupTechnicians and associate professionals  0.06 0.10 35643.20
## 17                                         occupUnemployed -0.04 0.12 35640.61
## 18                                       other.party.dummy  0.12 0.04   116.00
## 19                                         dont.know.dummy -0.02 0.06   198.20
## 20                                      invalid.vote.dummy -0.71 0.36  2073.56
## 21                                         no.answer.dummy  0.10 0.39  2692.24
## 22                                  not.eligible.age.dummy  0.29 0.06   196.18
## 23                          not.eligible.citizenship.dummy  0.14 0.07   186.70
## 24                                not.eligible.other.dummy  0.20 0.09   666.71
## 25                                    anti.imm.party.dummy -0.10 0.05   134.56
## 26                                     pro.env.party.dummy  0.60 0.07    29.76
##         t     p    LL    UL
## 1   -1.41 0.161 -0.41  0.07
## 2   -7.50 0.000 -0.00 -0.00
## 3    1.18 0.239 -0.01  0.04
## 4   13.27 0.000  0.02  0.03
## 5  -10.05 0.000 -0.16 -0.11
## 6    0.41 0.682 -0.16  0.24
## 7   -0.59 0.553 -0.25  0.14
## 8   -0.12 0.907 -0.21  0.19
## 9    1.35 0.178 -0.06  0.33
## 10   1.31 0.190 -0.07  0.34
## 11  -0.40 0.686 -0.24  0.16
## 12   2.06 0.039  0.01  0.40
## 13   0.50 0.618 -0.16  0.27
## 14   0.24 0.809 -0.17  0.22
## 15   0.45 0.649 -0.16  0.25
## 16   0.58 0.559 -0.14  0.25
## 17  -0.29 0.770 -0.27  0.20
## 18   2.94 0.004  0.04  0.20
## 19  -0.38 0.705 -0.14  0.10
## 20  -1.97 0.048 -1.42 -0.00
## 21   0.25 0.799 -0.66  0.86
## 22   4.83 0.000  0.17  0.41
## 23   2.16 0.032  0.01  0.27
## 24   2.09 0.037  0.01  0.38
## 25  -1.86 0.065 -0.20  0.01
## 26   8.53 0.000  0.46  0.74
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.02 0.12   194.57
## 2                                                      age -0.00 0.00 35715.03
## 3                                                   gender  0.02 0.01 35632.32
## 4                                                     educ  0.03 0.00 35736.79
## 5                                                    resid -0.14 0.01 35717.92
## 6                            occupClerical support workers  0.04 0.10 35587.55
## 7                    occupCraft and related trades workers -0.06 0.10 35595.29
## 8                              occupElementary occupations -0.01 0.10 35591.57
## 9                                            occupManagers  0.13 0.10 35588.02
## 10                            occupOther: Not in paid work  0.14 0.10 35642.79
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35594.52
## 12                                      occupProfessionals  0.21 0.10 35590.81
## 13                                            occupRetired  0.06 0.11 35584.13
## 14                          occupService and sales workers  0.02 0.10 35587.43
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35604.00
## 16            occupTechnicians and associate professionals  0.06 0.10 35586.88
## 17                                         occupUnemployed -0.03 0.12 35584.96
## 18                                      did.not.vote.dummy -0.15 0.06   153.99
## 19                                         dont.know.dummy -0.17 0.07   280.18
## 20                                      invalid.vote.dummy -0.83 0.40  1108.73
## 21                                         no.answer.dummy -0.04 0.42  1381.75
## 22                                  not.eligible.age.dummy  0.15 0.07   290.03
## 23                          not.eligible.citizenship.dummy  0.00 0.07   274.56
## 24                                not.eligible.other.dummy  0.08 0.10   863.54
##         t     p    LL    UL
## 1   -0.21 0.834 -0.26  0.21
## 2   -7.62 0.000 -0.00 -0.00
## 3    1.30 0.193 -0.01  0.04
## 4   13.50 0.000  0.02  0.03
## 5  -10.15 0.000 -0.16 -0.11
## 6    0.43 0.665 -0.15  0.24
## 7   -0.58 0.562 -0.25  0.14
## 8   -0.10 0.924 -0.21  0.19
## 9    1.33 0.183 -0.06  0.33
## 10   1.35 0.176 -0.06  0.34
## 11  -0.38 0.701 -0.24  0.16
## 12   2.10 0.036  0.01  0.40
## 13   0.51 0.611 -0.16  0.27
## 14   0.25 0.802 -0.17  0.22
## 15   0.46 0.647 -0.16  0.26
## 16   0.60 0.549 -0.14  0.25
## 17  -0.26 0.792 -0.27  0.21
## 18  -2.63 0.009 -0.26 -0.04
## 19  -2.60 0.010 -0.30 -0.04
## 20  -2.07 0.038 -1.61 -0.04
## 21  -0.09 0.930 -0.87  0.79
## 22   2.29 0.023  0.02  0.28
## 23   0.02 0.986 -0.14  0.15
## 24   0.74 0.458 -0.12  0.27
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept)  0.02 0.12   201.45
## 2                                                      age -0.00 0.00 35693.16
## 3                                                   gender  0.02 0.01 35634.76
## 4                                                     educ  0.03 0.00 35737.10
## 5                                                    resid -0.13 0.01 35728.67
## 6                            occupClerical support workers  0.04 0.10 35598.90
## 7                    occupCraft and related trades workers -0.06 0.10 35607.53
## 8                              occupElementary occupations -0.01 0.10 35603.29
## 9                                            occupManagers  0.13 0.10 35599.54
## 10                            occupOther: Not in paid work  0.14 0.10 35657.41
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35606.85
## 12                                      occupProfessionals  0.21 0.10 35601.95
## 13                                            occupRetired  0.06 0.11 35594.59
## 14                          occupService and sales workers  0.03 0.10 35598.77
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35616.48
## 16            occupTechnicians and associate professionals  0.06 0.10 35598.29
## 17                                         occupUnemployed -0.03 0.12 35594.88
## 18                                      did.not.vote.dummy -0.19 0.05   144.47
## 19                                         dont.know.dummy -0.22 0.06   279.00
## 20                                      invalid.vote.dummy -0.86 0.39  1274.76
## 21                                         no.answer.dummy -0.08 0.41  1605.89
## 22                                  not.eligible.age.dummy  0.10 0.06   292.43
## 23                          not.eligible.citizenship.dummy -0.05 0.07   265.67
## 24                                not.eligible.other.dummy  0.02 0.10   895.80
## 25                                    anti.imm.party.dummy -0.29 0.05   209.84
##         t     p    LL    UL
## 1    0.17 0.867 -0.21  0.25
## 2   -7.73 0.000 -0.00 -0.00
## 3    1.25 0.210 -0.01  0.04
## 4   13.43 0.000  0.02  0.03
## 5  -10.12 0.000 -0.16 -0.11
## 6    0.43 0.667 -0.15  0.24
## 7   -0.56 0.575 -0.25  0.14
## 8   -0.09 0.927 -0.21  0.19
## 9    1.34 0.182 -0.06  0.33
## 10   1.36 0.175 -0.06  0.34
## 11  -0.37 0.708 -0.24  0.16
## 12   2.10 0.036  0.01  0.40
## 13   0.52 0.606 -0.16  0.28
## 14   0.26 0.795 -0.17  0.22
## 15   0.46 0.646 -0.16  0.26
## 16   0.60 0.547 -0.13  0.25
## 17  -0.27 0.791 -0.27  0.21
## 18  -3.70 0.000 -0.29 -0.09
## 19  -3.46 0.001 -0.34 -0.09
## 20  -2.23 0.026 -1.62 -0.10
## 21  -0.20 0.839 -0.89  0.72
## 22   1.63 0.103 -0.02  0.23
## 23  -0.67 0.504 -0.18  0.09
## 24   0.22 0.830 -0.17  0.21
## 25  -5.78 0.000 -0.39 -0.19
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.09 0.12   238.38
## 2                                                      age -0.00 0.00 35625.90
## 3                                                   gender  0.02 0.01 35669.67
## 4                                                     educ  0.03 0.00 35725.28
## 5                                                    resid -0.13 0.01 35738.79
## 6                            occupClerical support workers  0.04 0.10 35631.32
## 7                    occupCraft and related trades workers -0.06 0.10 35639.11
## 8                              occupElementary occupations -0.01 0.10 35635.03
## 9                                            occupManagers  0.13 0.10 35632.76
## 10                            occupOther: Not in paid work  0.14 0.10 35686.91
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35639.22
## 12                                      occupProfessionals  0.20 0.10 35632.21
## 13                                            occupRetired  0.05 0.11 35625.94
## 14                          occupService and sales workers  0.02 0.10 35630.70
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35649.61
## 16            occupTechnicians and associate professionals  0.06 0.10 35630.33
## 17                                         occupUnemployed -0.04 0.12 35626.79
## 18                                      did.not.vote.dummy -0.08 0.04   139.08
## 19                                         dont.know.dummy -0.10 0.06   314.98
## 20                                      invalid.vote.dummy -0.80 0.37  1875.75
## 21                                         no.answer.dummy  0.03 0.40  2402.60
## 22                                  not.eligible.age.dummy  0.22 0.06   336.32
## 23                          not.eligible.citizenship.dummy  0.07 0.06   282.81
## 24                                not.eligible.other.dummy  0.13 0.09  1101.05
## 25                                     pro.env.party.dummy  0.53 0.05   271.99
##         t     p    LL    UL
## 1   -0.78 0.437 -0.32  0.14
## 2   -7.41 0.000 -0.00 -0.00
## 3    1.24 0.214 -0.01  0.04
## 4   13.35 0.000  0.02  0.03
## 5  -10.08 0.000 -0.16 -0.11
## 6    0.40 0.686 -0.16  0.24
## 7   -0.62 0.537 -0.26  0.13
## 8   -0.13 0.899 -0.21  0.18
## 9    1.34 0.182 -0.06  0.33
## 10   1.32 0.188 -0.07  0.34
## 11  -0.42 0.677 -0.24  0.16
## 12   2.06 0.039  0.01  0.40
## 13   0.48 0.629 -0.16  0.27
## 14   0.23 0.818 -0.17  0.22
## 15   0.45 0.654 -0.16  0.25
## 16   0.58 0.565 -0.14  0.25
## 17  -0.30 0.768 -0.27  0.20
## 18  -1.74 0.084 -0.17  0.01
## 19  -1.79 0.074 -0.21  0.01
## 20  -2.16 0.031 -1.53 -0.07
## 21   0.07 0.942 -0.75  0.80
## 22   3.87 0.000  0.11  0.33
## 23   1.05 0.296 -0.06  0.19
## 24   1.40 0.163 -0.05  0.31
## 25  10.85 0.000  0.43  0.62
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
##                                                        Eff   Est   SE       df
## 1                                              (Intercept) -0.05 0.12   242.61
## 2                                                      age -0.00 0.00 35583.97
## 3                                                   gender  0.02 0.01 35672.29
## 4                                                     educ  0.03 0.00 35720.39
## 5                                                    resid -0.13 0.01 35734.98
## 6                            occupClerical support workers  0.04 0.10 35641.21
## 7                    occupCraft and related trades workers -0.06 0.10 35649.39
## 8                              occupElementary occupations -0.01 0.10 35645.08
## 9                                            occupManagers  0.13 0.10 35642.84
## 10                            occupOther: Not in paid work  0.14 0.10 35697.06
## 11        occupPlant and machine operators, and assemblers -0.04 0.10 35649.76
## 12                                      occupProfessionals  0.20 0.10 35641.83
## 13                                            occupRetired  0.05 0.11 35635.36
## 14                          occupService and sales workers  0.02 0.10 35640.53
## 15 occupSkilled agricultural, forestry and fishery workers  0.05 0.11 35660.17
## 16            occupTechnicians and associate professionals  0.06 0.10 35640.21
## 17                                         occupUnemployed -0.04 0.12 35636.70
## 18                                      did.not.vote.dummy -0.12 0.04   132.01
## 19                                         dont.know.dummy -0.14 0.05   316.43
## 20                                      invalid.vote.dummy -0.83 0.36  2173.94
## 21                                         no.answer.dummy -0.01 0.39  2799.75
## 22                                  not.eligible.age.dummy  0.18 0.05   342.85
## 23                          not.eligible.citizenship.dummy  0.02 0.06   276.20
## 24                                not.eligible.other.dummy  0.08 0.09  1145.84
## 25                                    anti.imm.party.dummy -0.22 0.04   204.52
## 26                                     pro.env.party.dummy  0.49 0.05   263.65
##         t     p    LL    UL
## 1   -0.44 0.663 -0.28  0.18
## 2   -7.51 0.000 -0.00 -0.00
## 3    1.19 0.235 -0.01  0.04
## 4   13.27 0.000  0.02  0.03
## 5  -10.05 0.000 -0.16 -0.11
## 6    0.40 0.689 -0.16  0.24
## 7   -0.60 0.548 -0.26  0.14
## 8   -0.13 0.900 -0.21  0.18
## 9    1.34 0.181 -0.06  0.33
## 10   1.32 0.188 -0.07  0.34
## 11  -0.41 0.682 -0.24  0.16
## 12   2.06 0.039  0.01  0.40
## 13   0.49 0.624 -0.16  0.27
## 14   0.24 0.812 -0.17  0.22
## 15   0.45 0.656 -0.16  0.25
## 16   0.58 0.564 -0.14  0.25
## 17  -0.30 0.765 -0.27  0.20
## 18  -2.79 0.006 -0.20 -0.03
## 19  -2.59 0.010 -0.25 -0.03
## 20  -2.29 0.022 -1.54 -0.12
## 21  -0.03 0.975 -0.77  0.75
## 22   3.22 0.001  0.07  0.28
## 23   0.39 0.697 -0.09  0.14
## 24   0.90 0.367 -0.10  0.26
## 25  -5.09 0.000 -0.30 -0.13
## 26  10.40 0.000  0.40  0.58
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
