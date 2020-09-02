---
title: "Alignment attempt"
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
    toc_depth: 5
---



# Preparations

### Load packages


```r
library(sirt)
```

```
## - sirt 3.9-4 (2020-02-17 12:57:09)
```

```r
library(psych)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(lavaan)
```

```
## This is lavaan 0.6-5
```

```
## lavaan is BETA software! Please report any bugs.
```

```
## 
## Attaching package: 'lavaan'
```

```
## The following object is masked from 'package:psych':
## 
##     cor2cov
```

### Read data


```r
dat<-read.csv2("dat.no.miss.csv",stringsAsFactors = F)

# data for refugee attitudes
imm.vars<-c("gvrfgap.R","rfgfrpc","rfgbfml.R")

# standardization of the indicators (as suggested by Asparouhov & Muthen, 2014)
na.standardize<-function(var){
  (var-mean(var,na.rm=T))/sd(var,na.rm=T)
}

dat$imm.1.z<-na.standardize(dat$gvrfgap.R)
dat$imm.2.z<-na.standardize(dat$rfgfrpc)
dat$imm.3.z<-na.standardize(dat$rfgbfml.R)

imm.z.vars<-c("imm.1.z","imm.2.z","imm.3.z")

describe(dat[,imm.z.vars])
```

```
##         vars     n mean sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 36425    0  1   0.06    0.00 1.25 -1.63 1.74  3.37 -0.05    -1.00
## imm.2.z    2 33198    0  1   0.09    0.01 1.40 -1.81 1.98  3.79  0.07    -0.73
## imm.3.z    3 34927    0  1   0.66    0.05 1.33 -2.04 1.56  3.60 -0.50    -0.67
##           se
## imm.1.z 0.01
## imm.2.z 0.01
## imm.3.z 0.01
```

```r
describeBy(dat[,imm.z.vars],group=dat[,"cntry"])
```

```
## Warning in FUN(newX[, i], ...): no non-missing arguments to min; returning Inf

## Warning in FUN(newX[, i], ...): no non-missing arguments to min; returning Inf
```

```
## Warning in FUN(newX[, i], ...): no non-missing arguments to max; returning -Inf

## Warning in FUN(newX[, i], ...): no non-missing arguments to max; returning -Inf
```

```
## 
##  Descriptive statistics by group 
## group: AT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1950 -0.29 1.04  -0.79   -0.34 1.25 -1.63 1.74  3.37  0.26
## imm.2.z    2 1868 -0.23 1.05   0.09   -0.27 1.40 -1.81 1.98  3.79  0.35
## imm.3.z    3 1918 -0.31 1.14  -0.24   -0.33 1.33 -2.04 1.56  3.60 -0.08
##         kurtosis   se
## imm.1.z    -1.02 0.02
## imm.2.z    -0.59 0.02
## imm.3.z    -1.11 0.03
## ------------------------------------------------------------ 
## group: BE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1750 -0.17 1.04  -0.79   -0.22 1.25 -1.63 1.74  3.37  0.30
## imm.2.z    2 1696  0.25 1.04   0.09    0.26 1.40 -1.81 1.98  3.79 -0.20
## imm.3.z    3 1750 -0.23 1.09  -0.24   -0.23 1.33 -2.04 1.56  3.60 -0.09
##         kurtosis   se
## imm.1.z    -0.98 0.02
## imm.2.z    -0.75 0.03
## imm.3.z    -1.09 0.03
## ------------------------------------------------------------ 
## group: CH
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 1493 0.06 0.91   0.06    0.07 1.25 -1.63 1.74  3.37 -0.06    -0.82
## imm.2.z    2 1385 0.08 0.95   0.09    0.08 1.40 -1.81 1.98  3.79  0.03    -0.71
## imm.3.z    3 1478 0.01 0.99   0.66    0.04 1.33 -2.04 1.56  3.60 -0.44    -0.75
##           se
## imm.1.z 0.02
## imm.2.z 0.03
## imm.3.z 0.03
## ------------------------------------------------------------ 
## group: CZ
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 2140 -0.74 0.91  -0.79   -0.85 1.25 -1.63 1.74  3.37  0.76
## imm.2.z    2 2080 -0.55 1.04  -0.86   -0.64 1.40 -1.81 1.98  3.79  0.56
## imm.3.z    3 2106 -0.30 1.08  -0.24   -0.28 1.33 -2.04 1.56  3.60 -0.21
##         kurtosis   se
## imm.1.z    -0.35 0.02
## imm.2.z    -0.43 0.02
## imm.3.z    -1.01 0.02
## ------------------------------------------------------------ 
## group: DE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 2812 -0.22 0.94  -0.79   -0.25 1.25 -1.63 1.74  3.37  0.32
## imm.2.z    2 2710  0.04 0.96   0.09    0.04 1.40 -1.81 1.98  3.79  0.03
## imm.3.z    3 2802  0.06 0.97   0.66    0.11 1.33 -2.04 1.56  3.60 -0.59
##         kurtosis   se
## imm.1.z    -0.77 0.02
## imm.2.z    -0.75 0.02
## imm.3.z    -0.64 0.02
## ------------------------------------------------------------ 
## group: EE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1958 -0.59 0.84  -0.79   -0.67 1.25 -1.63 1.74  3.37  0.68
## imm.2.z    2 1930 -0.43 0.90  -0.86   -0.47 1.40 -1.81 1.98  3.79  0.40
## imm.3.z    3 1967 -0.35 1.00  -0.24   -0.29 1.33 -2.04 1.56  3.60 -0.20
##         kurtosis   se
## imm.1.z    -0.08 0.02
## imm.2.z    -0.23 0.02
## imm.3.z    -1.07 0.02
## ------------------------------------------------------------ 
## group: ES
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 1747 0.56 0.88   0.90    0.63 1.25 -1.63 1.74  3.37 -0.60    -0.17
## imm.2.z    2 1397 0.28 1.01   0.09    0.26 1.40 -1.81 1.98  3.79 -0.08    -0.81
## imm.3.z    3 1762 0.55 0.85   0.66    0.66 0.00 -2.04 1.56  3.60 -0.95     0.75
##           se
## imm.1.z 0.02
## imm.2.z 0.03
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: FI
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 1850 0.15 0.88   0.06    0.15 1.25 -1.63 1.74  3.37 -0.08    -0.69
## imm.2.z    2 1838 0.13 0.93   0.09    0.15 1.40 -1.81 1.98  3.79 -0.14    -0.67
## imm.3.z    3 1858 0.41 0.82   0.66    0.49 0.00 -2.04 1.56  3.60 -1.00     0.89
##           se
## imm.1.z 0.02
## imm.2.z 0.02
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: FR
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 2000  0.40 1.12   0.90    0.49 1.25 -1.63 1.74  3.37 -0.47
## imm.2.z    2 1913  0.30 1.09   0.09    0.32 1.40 -1.81 1.98  3.79 -0.09
## imm.3.z    3 1995 -0.08 1.18   0.66   -0.04 1.33 -2.04 1.56  3.60 -0.36
##         kurtosis   se
## imm.1.z    -0.94 0.02
## imm.2.z    -0.83 0.02
## imm.3.z    -1.08 0.03
## ------------------------------------------------------------ 
## group: GB
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1867  0.32 0.90   0.06    0.33 1.25 -1.63 1.74  3.37 -0.30
## imm.2.z    2 1821  0.07 0.90   0.09    0.06 1.40 -1.81 1.98  3.79  0.04
## imm.3.z    3 1862 -0.18 0.94  -0.24   -0.12 1.33 -2.04 1.56  3.60 -0.31
##         kurtosis   se
## imm.1.z    -0.69 0.02
## imm.2.z    -0.54 0.02
## imm.3.z    -0.94 0.02
## ------------------------------------------------------------ 
## group: HU
##         vars    n  mean   sd median trimmed  mad   min  max range skew kurtosis
## imm.1.z    1 1391 -0.69 0.95  -0.79   -0.79 1.25 -1.63 1.74  3.37 0.63    -0.71
## imm.2.z    2    0   NaN   NA     NA     NaN   NA   Inf -Inf  -Inf   NA       NA
## imm.3.z    3    0   NaN   NA     NA     NaN   NA   Inf -Inf  -Inf   NA       NA
##           se
## imm.1.z 0.03
## imm.2.z   NA
## imm.3.z   NA
## ------------------------------------------------------------ 
## group: IE
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 2616 0.43 0.90   0.90    0.47 1.25 -1.63 1.74  3.37 -0.65    -0.37
## imm.2.z    2 2400 0.13 0.96   0.09    0.13 1.40 -1.81 1.98  3.79 -0.06    -0.90
## imm.3.z    3 2625 0.05 0.97   0.66    0.11 0.00 -2.04 1.56  3.60 -0.61    -0.63
##           se
## imm.1.z 0.02
## imm.2.z 0.02
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: IT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 2275 -0.16 0.96   0.06   -0.17 1.25 -1.63 1.74  3.37  0.05
## imm.2.z    2 2157 -0.24 0.98   0.09   -0.26 1.40 -1.81 1.98  3.79  0.32
## imm.3.z    3 2264  0.02 0.97   0.66    0.07 1.33 -2.04 1.56  3.60 -0.53
##         kurtosis   se
## imm.1.z    -0.92 0.02
## imm.2.z    -0.48 0.02
## imm.3.z    -0.43 0.02
## ------------------------------------------------------------ 
## group: LT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1854  0.06 0.82   0.06    0.10 1.25 -1.63 1.74  3.37 -0.25
## imm.2.z    2 1779 -0.37 0.91  -0.86   -0.41 1.40 -1.81 1.98  3.79  0.43
## imm.3.z    3 1853  0.02 0.86  -0.24    0.10 1.33 -2.04 1.56  3.60 -0.67
##         kurtosis   se
## imm.1.z    -0.55 0.02
## imm.2.z    -0.09 0.02
## imm.3.z    -0.01 0.02
## ------------------------------------------------------------ 
## group: NL
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1650 -0.48 0.83  -0.79   -0.53 0.00 -1.63 1.74  3.37  0.66
## imm.2.z    2 1470  0.29 0.86   0.09    0.31 1.40 -1.81 1.98  3.79 -0.23
## imm.3.z    3 1644 -0.24 0.94  -0.24   -0.19 1.33 -2.04 1.56  3.60 -0.22
##         kurtosis   se
## imm.1.z    -0.27 0.02
## imm.2.z    -0.59 0.02
## imm.3.z    -1.09 0.02
## ------------------------------------------------------------ 
## group: NO
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 1533 0.42 0.86   0.90    0.45 1.25 -1.63 1.74  3.37 -0.50    -0.39
## imm.2.z    2 1516 0.51 0.90   1.03    0.54 1.40 -1.81 1.98  3.79 -0.44    -0.16
## imm.3.z    3 1533 0.20 0.84   0.66    0.26 0.00 -2.04 1.56  3.60 -0.74     0.06
##           se
## imm.1.z 0.02
## imm.2.z 0.02
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: PL
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1550  0.31 0.82   0.06    0.34 1.25 -1.63 1.74  3.37 -0.44
## imm.2.z    2 1469 -0.28 0.88  -0.86   -0.31 1.40 -1.81 1.98  3.79  0.40
## imm.3.z    3 1533  0.14 0.92   0.66    0.18 1.33 -2.04 1.56  3.60 -0.64
##         kurtosis   se
## imm.1.z    -0.29 0.02
## imm.2.z    -0.38 0.02
## imm.3.z    -0.19 0.02
## ------------------------------------------------------------ 
## group: PT
##         vars    n mean   sd median trimmed mad   min  max range  skew kurtosis
## imm.1.z    1 1219 0.70 0.72   0.90    0.77 0.0 -1.63 1.74  3.37 -0.90     0.89
## imm.2.z    2 1161 0.04 0.87   0.09    0.02 1.4 -1.81 1.98  3.79  0.14    -0.90
## imm.3.z    3 1224 0.34 0.78   0.66    0.42 0.0 -2.04 1.56  3.60 -1.06     0.78
##           se
## imm.1.z 0.02
## imm.2.z 0.03
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: SE
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## imm.1.z    1 1507 0.38 0.80   0.06    0.39 1.25 -1.63 1.74  3.37 -0.32    -0.25
## imm.2.z    2 1401 0.38 0.83   0.09    0.38 1.40 -1.81 1.98  3.79 -0.12    -0.20
## imm.3.z    3 1506 0.14 0.88   0.66    0.16 1.33 -2.04 1.56  3.60 -0.52    -0.25
##           se
## imm.1.z 0.02
## imm.2.z 0.02
## imm.3.z 0.02
## ------------------------------------------------------------ 
## group: SI
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## imm.1.z    1 1263 -0.13 0.87   0.06   -0.12 1.25 -1.63 1.74  3.37  0.13
## imm.2.z    2 1207 -0.04 0.94   0.09   -0.06 1.40 -1.81 1.98  3.79  0.22
## imm.3.z    3 1247 -0.10 0.97  -0.24   -0.05 1.33 -2.04 1.56  3.60 -0.37
##         kurtosis   se
## imm.1.z    -0.77 0.02
## imm.2.z    -0.80 0.03
## imm.3.z    -0.92 0.03
```

```r
#exclude Hungary
dat.imm<-dat %>%
  filter(cntry!="HU") %>%
  dplyr::select(cntry,all_of(imm.z.vars),all_of(imm.vars))
dat.imm<-na.omit(dat.imm)


env.vars<-c("inctxff.R","sbsrnen.R","banhhap.R")

dat$env.1.z<-na.standardize(dat$inctxff.R)
dat$env.2.z<-na.standardize(dat$sbsrnen.R)
dat$env.3.z<-na.standardize(dat$banhhap.R)

env.z.vars<-c("env.1.z","env.2.z","env.3.z")

describe(dat[,env.z.vars])
```

```
##         vars     n mean sd median trimmed  mad   min  max range  skew kurtosis
## env.1.z    1 36131    0  1   0.18   -0.02 1.19 -1.44 1.79  3.22  0.08    -1.10
## env.2.z    2 36509    0  1   0.02    0.16 1.41 -2.84 0.97  3.81 -1.15     0.89
## env.3.z    3 36318    0  1   0.38    0.08 1.27 -2.19 1.23  3.42 -0.58    -0.53
##           se
## env.1.z 0.01
## env.2.z 0.01
## env.3.z 0.01
```

```r
describeBy(dat[,env.z.vars],group=dat[,"cntry"])
```

```
## 
##  Descriptive statistics by group 
## group: AT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1956 -0.01 0.97   0.18   -0.03 1.19 -1.44 1.79  3.22  0.12
## env.2.z    2 1965  0.20 0.85   0.02    0.33 1.41 -2.84 0.97  3.81 -1.32
## env.3.z    3 1945  0.18 0.95   0.38    0.28 1.27 -2.19 1.23  3.42 -0.70
##         kurtosis   se
## env.1.z    -1.01 0.02
## env.2.z     1.99 0.02
## env.3.z    -0.31 0.02
## ------------------------------------------------------------ 
## group: BE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1749 -0.02 0.99   0.18   -0.05 1.19 -1.44 1.79  3.22  0.14
## env.2.z    2 1752 -0.03 0.99   0.02    0.12 1.41 -2.84 0.97  3.81 -1.14
## env.3.z    3 1751  0.15 0.93   0.38    0.24 1.27 -2.19 1.23  3.42 -0.76
##         kurtosis   se
## env.1.z    -1.08 0.02
## env.2.z     0.88 0.02
## env.3.z    -0.12 0.02
## ------------------------------------------------------------ 
## group: CH
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## env.1.z    1 1477 0.35 0.95   0.18    0.37 1.19 -1.44 1.79  3.22 -0.20    -0.98
## env.2.z    2 1492 0.14 0.85   0.02    0.27 1.41 -2.84 0.97  3.81 -1.20     1.50
## env.3.z    3 1491 0.20 0.98   0.38    0.30 1.27 -2.19 1.23  3.42 -0.76    -0.40
##           se
## env.1.z 0.02
## env.2.z 0.02
## env.3.z 0.03
## ------------------------------------------------------------ 
## group: CZ
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 2083 -0.14 1.04   0.18   -0.19 1.19 -1.44 1.79  3.22  0.22
## env.2.z    2 2121 -0.48 1.26   0.02   -0.37 1.41 -2.84 0.97  3.81 -0.64
## env.3.z    3 2128 -0.12 1.17   0.38   -0.03 1.27 -2.19 1.23  3.42 -0.52
##         kurtosis   se
## env.1.z    -1.15 0.02
## env.2.z    -0.77 0.03
## env.3.z    -1.00 0.03
## ------------------------------------------------------------ 
## group: DE
##         vars    n mean   sd median trimmed  mad   min  max range  skew kurtosis
## env.1.z    1 2798 0.15 0.93   0.18    0.16 1.19 -1.44 1.79  3.22 -0.01    -0.97
## env.2.z    2 2810 0.11 0.91   0.02    0.27 1.41 -2.84 0.97  3.81 -1.33     1.68
## env.3.z    3 2811 0.22 1.00   0.38    0.34 1.27 -2.19 1.23  3.42 -0.78    -0.38
##           se
## env.1.z 0.02
## env.2.z 0.02
## env.3.z 0.02
## ------------------------------------------------------------ 
## group: EE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1935 -0.16 0.83   0.18   -0.18 1.19 -1.44 1.79  3.22  0.26
## env.2.z    2 1961 -0.06 0.88   0.02    0.06 0.00 -2.84 0.97  3.81 -1.12
## env.3.z    3 1964 -0.23 0.91  -0.48   -0.21 1.27 -2.19 1.23  3.42 -0.36
##         kurtosis   se
## env.1.z    -0.53 0.02
## env.2.z     1.49 0.02
## env.3.z    -0.54 0.02
## ------------------------------------------------------------ 
## group: ES
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1720 -0.27 0.99  -0.63   -0.33 1.19 -1.44 1.79  3.22  0.43
## env.2.z    2 1763  0.06 1.06   0.02    0.25 1.41 -2.84 0.97  3.81 -1.23
## env.3.z    3 1723  0.12 0.94   0.38    0.22 1.27 -2.19 1.23  3.42 -0.73
##         kurtosis   se
## env.1.z    -0.97 0.02
## env.2.z     0.78 0.03
## env.3.z    -0.11 0.02
## ------------------------------------------------------------ 
## group: FI
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1852  0.47 0.87   0.98    0.50 1.19 -1.44 1.79  3.22 -0.44
## env.2.z    2 1857 -0.05 0.92   0.02    0.09 0.00 -2.84 0.97  3.81 -1.09
## env.3.z    3 1849 -0.01 0.87   0.38    0.04 1.27 -2.19 1.23  3.42 -0.52
##         kurtosis   se
## env.1.z    -0.51 0.02
## env.2.z     1.02 0.02
## env.3.z    -0.22 0.02
## ------------------------------------------------------------ 
## group: FR
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 2010 -0.22 0.95  -0.63   -0.27 1.19 -1.44 1.79  3.22  0.35
## env.2.z    2 2006 -0.09 0.98   0.02    0.04 1.41 -2.84 0.97  3.81 -1.08
## env.3.z    3 2009  0.11 0.96   0.38    0.20 1.27 -2.19 1.23  3.42 -0.71
##         kurtosis   se
## env.1.z    -0.87 0.02
## env.2.z     0.80 0.02
## env.3.z    -0.28 0.02
## ------------------------------------------------------------ 
## group: GB
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1865  0.07 0.94   0.18    0.08 1.19 -1.44 1.79  3.22 -0.04
## env.2.z    2 1868 -0.26 1.02   0.02   -0.16 1.41 -2.84 0.97  3.81 -0.78
## env.3.z    3 1860 -0.10 0.97   0.38   -0.05 1.27 -2.19 1.23  3.42 -0.45
##         kurtosis   se
## env.1.z    -0.99 0.02
## env.2.z     0.00 0.02
## env.3.z    -0.64 0.02
## ------------------------------------------------------------ 
## group: HU
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1341 -0.12 1.02   0.18   -0.16 1.19 -1.44 1.79  3.22  0.14
## env.2.z    2 1377  0.43 0.87   0.97    0.63 0.00 -2.84 0.97  3.81 -1.94
## env.3.z    3 1350 -0.11 1.02   0.38   -0.04 1.27 -2.19 1.23  3.42 -0.44
##         kurtosis   se
## env.1.z    -1.17 0.03
## env.2.z     3.60 0.02
## env.3.z    -0.66 0.03
## ------------------------------------------------------------ 
## group: IE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 2649 -0.16 1.03  -0.63   -0.20 1.19 -1.44 1.79  3.22  0.20
## env.2.z    2 2645 -0.28 1.15   0.02   -0.15 1.41 -2.84 0.97  3.81 -0.83
## env.3.z    3 2619 -0.21 1.07   0.38   -0.15 1.27 -2.19 1.23  3.42 -0.41
##         kurtosis   se
## env.1.z    -1.27 0.02
## env.2.z    -0.24 0.02
## env.3.z    -0.90 0.02
## ------------------------------------------------------------ 
## group: IT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 2210 -0.15 1.00  -0.63   -0.20 1.19 -1.44 1.79  3.22  0.30
## env.2.z    2 2282 -0.07 1.02   0.02    0.08 1.41 -2.84 0.97  3.81 -1.03
## env.3.z    3 2260  0.23 0.86   0.38    0.34 1.27 -2.19 1.23  3.42 -0.84
##         kurtosis   se
## env.1.z    -1.00 0.02
## env.2.z     0.56 0.02
## env.3.z     0.40 0.02
## ------------------------------------------------------------ 
## group: LT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1819 -0.09 1.03   0.18   -0.15 1.19 -1.44 1.79  3.22  0.21
## env.2.z    2 1890 -0.14 0.96   0.02   -0.04 1.41 -2.84 0.97  3.81 -0.74
## env.3.z    3 1872 -0.33 0.98  -0.48   -0.31 1.27 -2.19 1.23  3.42 -0.20
##         kurtosis   se
## env.1.z    -1.06 0.02
## env.2.z     0.21 0.02
## env.3.z    -0.72 0.02
## ------------------------------------------------------------ 
## group: NL
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1631  0.11 1.00   0.18    0.12 1.19 -1.44 1.79  3.22 -0.11
## env.2.z    2 1636  0.23 0.86   0.02    0.39 1.41 -2.84 0.97  3.81 -1.55
## env.3.z    3 1640 -0.12 1.06   0.38   -0.05 1.27 -2.19 1.23  3.42 -0.47
##         kurtosis   se
## env.1.z    -1.16 0.02
## env.2.z     2.75 0.02
## env.3.z    -0.85 0.03
## ------------------------------------------------------------ 
## group: NO
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1535  0.34 1.00   0.18    0.38 1.19 -1.44 1.79  3.22 -0.31
## env.2.z    2 1536  0.24 0.80   0.02    0.36 1.41 -2.84 0.97  3.81 -1.27
## env.3.z    3 1528 -0.18 0.96   0.38   -0.14 1.27 -2.19 1.23  3.42 -0.40
##         kurtosis   se
## env.1.z    -1.02 0.03
## env.2.z     1.92 0.02
## env.3.z    -0.60 0.02
## ------------------------------------------------------------ 
## group: PL
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1531 -0.35 0.84  -0.63   -0.40 1.19 -1.44 1.79  3.22  0.39
## env.2.z    2 1563  0.01 0.91   0.02    0.14 1.41 -2.84 0.97  3.81 -1.09
## env.3.z    3 1554  0.04 0.93   0.38    0.11 1.27 -2.19 1.23  3.42 -0.59
##         kurtosis   se
## env.1.z    -0.62 0.02
## env.2.z     1.12 0.02
## env.3.z    -0.29 0.02
## ------------------------------------------------------------ 
## group: PT
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1223 -0.15 1.06  -0.63   -0.21 1.19 -1.44 1.79  3.22  0.24
## env.2.z    2 1206 -0.20 1.23   0.02   -0.03 1.41 -2.84 0.97  3.81 -0.88
## env.3.z    3 1204  0.25 1.01   0.38    0.38 1.27 -2.19 1.23  3.42 -0.90
##         kurtosis   se
## env.1.z    -1.21 0.03
## env.2.z    -0.35 0.04
## env.3.z    -0.06 0.03
## ------------------------------------------------------------ 
## group: SE
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1505  0.56 0.97   0.98    0.65 1.19 -1.44 1.79  3.22 -0.64
## env.2.z    2 1513  0.26 0.83   0.02    0.40 1.41 -2.84 0.97  3.81 -1.37
## env.3.z    3 1500 -0.19 1.03   0.38   -0.12 1.27 -2.19 1.23  3.42 -0.38
##         kurtosis   se
## env.1.z    -0.53 0.03
## env.2.z     2.09 0.02
## env.3.z    -0.78 0.03
## ------------------------------------------------------------ 
## group: SI
##         vars    n  mean   sd median trimmed  mad   min  max range  skew
## env.1.z    1 1242 -0.11 0.99   0.18   -0.13 1.19 -1.44 1.79  3.22  0.07
## env.2.z    2 1266  0.47 0.78   0.97    0.64 0.00 -2.84 0.97  3.81 -1.97
## env.3.z    3 1260  0.14 1.01   0.38    0.24 1.27 -2.19 1.23  3.42 -0.74
##         kurtosis   se
## env.1.z    -1.19 0.03
## env.2.z     4.43 0.02
## env.3.z    -0.37 0.03
```

```r
dat.env<-dat %>%
  dplyr::select(cntry,all_of(env.z.vars),all_of(env.vars))
dat.env<-na.omit(dat.env)
```

\newpage

# Alignment for refugee attitudes

Approach by Fischer and Karl (2019)

### Configural models


```r
par.imm <- invariance_alignment_cfa_config(dat = dat.imm[,imm.z.vars],
                                       group = dat.imm[,"cntry"])
```

```
## Compute CFA for group 1
## Compute CFA for group 2
## Compute CFA for group 3
## Compute CFA for group 4
```

```
## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
## variances are negative
```

```
## Compute CFA for group 5
## Compute CFA for group 6
## Compute CFA for group 7
## Compute CFA for group 8
## Compute CFA for group 9
## Compute CFA for group 10
## Compute CFA for group 11
## Compute CFA for group 12
## Compute CFA for group 13
```

```
## Warning in lav_model_estimate(lavmodel = lavmodel, lavpartable = lavpartable, :
## lavaan WARNING: the optimizer warns that a solution has NOT been found!
```

```
## Compute CFA for group 14
## Compute CFA for group 15
## Compute CFA for group 16
```

```
## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
## variances are negative
```

```
## Compute CFA for group 17
## Compute CFA for group 18
## Compute CFA for group 19
```

```r
par.imm
```

```
## $nu
##        imm.1.z     imm.2.z      imm.3.z
## AT -0.28552068 -0.22846558 -0.323812213
## BE -0.17008024  0.25219438 -0.233089670
## CH  0.06533838  0.08077002 -0.006555338
## CZ -0.72918204 -0.55619500 -0.305610677
## DE -0.22390999  0.03711166  0.050426379
## EE -0.59343926 -0.43106077 -0.350092186
## ES  0.60013974  0.29165310  0.546240735
## FI  0.15002613  0.12445998  0.411128079
## FR  0.41293860  0.30358733 -0.069160249
## GB  0.31330298  0.06623068 -0.182609363
## IE  0.43334152  0.13420510  0.048148958
## IT -0.16618042 -0.23509948  0.003369280
## LT  0.05884075 -0.38252932 -0.013229983
## NL -0.48063230  0.29072854 -0.244042514
## NO  0.42582401  0.51208729  0.204494779
## PL  0.30274174 -0.27994058  0.121090586
## PT  0.69979362  0.03804292  0.335121115
## SE  0.38419353  0.38071737  0.129420395
## SI -0.11343749 -0.03988538 -0.086666036
## 
## $lambda
##        imm.1.z      imm.2.z     imm.3.z
## AT   0.8344492  0.685604992  0.86647918
## BE   0.8007933  0.470755052  0.80807877
## CH   0.6621644  0.431051329  0.63599550
## CZ   1.0671154  0.241053866  0.32373334
## DE   0.6052193  0.501894679  0.63574573
## EE   0.4382512  0.369404518  0.58887129
## ES   0.7923688  0.430682878  0.46670003
## FI   0.6666981  0.470492420  0.43868680
## FR   0.5721942  0.303261616  1.03853505
## GB   0.6264455  0.384756854  0.62488911
## IE   0.7074223  0.377362850  0.56270773
## IT   0.7924042  0.271051951  0.54335046
## LT -26.7558970 -0.003749287 -0.01299068
## NL   0.4019954  0.398296819  0.53099441
## NO   0.5183304  0.461367544  0.47527862
## PL   0.8737312  0.167346749  0.34060047
## PT   0.4580052  0.312620216  0.38218988
## SE   0.5873763  0.371284510  0.55449604
## SI   0.7228873  0.212564295  0.51214567
## 
## $err_var
##         imm.1.z   imm.2.z   imm.3.z
## AT    0.3882069 0.6504400 0.5472446
## BE    0.4486350 0.8627032 0.5521657
## CH    0.3915550 0.7124709 0.5624090
## CZ   -0.3045766 1.0357144 1.0596240
## DE    0.5237562 0.6643814 0.5420545
## EE    0.5180209 0.6783301 0.6447203
## ES    0.1644722 0.8316737 0.5333963
## FI    0.3282013 0.6492781 0.4885206
## FR    0.9121115 1.1058620 0.3111595
## GB    0.4160337 0.6590392 0.5018051
## IE    0.3249166 0.7827289 0.6221563
## IT    0.2859325 0.8821720 0.6495206
## LT -715.1960485 0.8149666 0.7442800
## NL    0.5231289 0.5819328 0.5995203
## NO    0.4770324 0.5988764 0.4719777
## PL   -0.0876767 0.7340509 0.7308303
## PT    0.2927106 0.6544284 0.4576938
## SE    0.2929494 0.5518801 0.4702758
## SI    0.2347983 0.8445844 0.6602236
## 
## $N
##   AT   BE   CH   CZ   DE   EE   ES   FI   FR   GB   IE   IT   LT   NL   NO   PL 
## 1811 1693 1369 2032 2693 1916 1355 1827 1893 1807 2337 2107 1682 1458 1509 1411 
##   PT   SE   SI 
## 1154 1379 1184 
## 
## $G
## [1] 19
## 
## $I
## [1] 3
## 
## $items
## [1] "imm.1.z" "imm.2.z" "imm.3.z"
## 
## $groups
##  [1] "AT" "BE" "CH" "CZ" "DE" "EE" "ES" "FI" "FR" "GB" "IE" "IT" "LT" "NL" "NO"
## [16] "PL" "PT" "SE" "SI"
## 
## $CALL
## invariance_alignment_cfa_config(dat = dat.imm[, imm.z.vars], 
##     group = dat.imm[, "cntry"])
```

For Lithuania, solution was not found. For Poland and Czech, there are negative variances. Also for France, the model seems to produce Heywood loadings and error variances (>1)

\newpage

### Construct separate configural models for the countries with problems


```r
conf.mod<-"
F.imm=~imm.1.z+imm.2.z+imm.3.z
"
```

\newpage

#### Configural model for Lithuania


```r
dat.imm.LT<-dat.imm %>%
  filter(cntry=="LT") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.LT<-cfa(data=dat.imm.LT,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
```

```
## Warning in lav_model_estimate(lavmodel = lavmodel, lavpartable = lavpartable, :
## lavaan WARNING: the optimizer warns that a solution has NOT been found!
```

```r
summary(imm.conf.LT)
```

```
## lavaan 0.6-5 did NOT end normally after 10000 iterations
## ** WARNING ** Estimates below are most likely unreliable
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          1682
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                    NA
##   Degrees of freedom                                NA
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate   Std.Err  z-value  P(>|z|)
##   F.imm =~                                             
##     imm.1.z          -32.424       NA                  
##     imm.2.z           -0.004       NA                  
##     imm.3.z           -0.015       NA                  
## 
## Variances:
##                    Estimate   Std.Err  z-value  P(>|z|)
##    .imm.1.z        -1050.345       NA                  
##    .imm.2.z            0.999       NA                  
##    .imm.3.z            0.999       NA                  
##     F.imm              1.000
```

\newpage

#### Configural model for Poland


```r
dat.imm.PL<-dat.imm %>%
  filter(cntry=="PL") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.PL<-cfa(data=dat.imm.PL,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
```

```
## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
## variances are negative
```

```r
summary(imm.conf.PL)
```

```
## lavaan 0.6-5 ended normally after 21 iterations
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          1411
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                 0.000
##   Degrees of freedom                                 0
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate  Std.Err  z-value  P(>|z|)
##   F.imm =~                                            
##     imm.1.z           1.063    0.181    5.886    0.000
##     imm.2.z           0.192    0.042    4.587    0.000
##     imm.3.z           0.370    0.068    5.471    0.000
## 
## Variances:
##                    Estimate  Std.Err  z-value  P(>|z|)
##    .imm.1.z          -0.130    0.382   -0.340    0.734
##    .imm.2.z           0.963    0.038   25.127    0.000
##    .imm.3.z           0.862    0.057   15.251    0.000
##     F.imm             1.000
```

\newpage

#### Configural model for Czech Republic


```r
dat.imm.CZ<-dat.imm %>%
  filter(cntry=="CZ") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.CZ<-cfa(data=dat.imm.CZ,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
```

```
## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
## variances are negative
```

```r
summary(imm.conf.CZ)
```

```
## lavaan 0.6-5 ended normally after 23 iterations
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          2032
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                 0.000
##   Degrees of freedom                                 0
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate  Std.Err  z-value  P(>|z|)
##   F.imm =~                                            
##     imm.1.z           1.168    0.171    6.841    0.000
##     imm.2.z           0.230    0.040    5.748    0.000
##     imm.3.z           0.300    0.049    6.146    0.000
## 
## Variances:
##                    Estimate  Std.Err  z-value  P(>|z|)
##    .imm.1.z          -0.365    0.398   -0.917    0.359
##    .imm.2.z           0.946    0.033   28.264    0.000
##    .imm.3.z           0.910    0.039   23.467    0.000
##     F.imm             1.000
```

\newpage

#### Configural model for France


```r
dat.imm.FR<-dat.imm %>%
  filter(cntry=="FR") %>%
  mutate(imm.1.z=na.standardize(gvrfgap.R),
         imm.2.z=na.standardize(rfgfrpc),
         imm.3.z=na.standardize(rfgbfml.R))


imm.conf.FR<-cfa(data=dat.imm.FR,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(imm.conf.FR)
```

```
## lavaan 0.6-5 ended normally after 18 iterations
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          1893
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                 0.000
##   Degrees of freedom                                 0
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate  Std.Err  z-value  P(>|z|)
##   F.imm =~                                            
##     imm.1.z           0.514    0.042   12.195    0.000
##     imm.2.z           0.277    0.030    9.283    0.000
##     imm.3.z           0.881    0.065   13.602    0.000
## 
## Variances:
##                    Estimate  Std.Err  z-value  P(>|z|)
##    .imm.1.z           0.735    0.044   16.603    0.000
##    .imm.2.z           0.923    0.032   28.934    0.000
##    .imm.3.z           0.224    0.110    2.038    0.042
##     F.imm             1.000
```


\newpage


# Alignment for environment attitudes


### Configural models


```r
par.env <- invariance_alignment_cfa_config(dat = dat.env[,env.z.vars],
                                       group = dat.env[,"cntry"])
```

```
## Compute CFA for group 1
## Compute CFA for group 2
## Compute CFA for group 3
## Compute CFA for group 4
## Compute CFA for group 5
## Compute CFA for group 6
## Compute CFA for group 7
## Compute CFA for group 8
## Compute CFA for group 9
## Compute CFA for group 10
## Compute CFA for group 11
```

```
## Warning in lav_model_estimate(lavmodel = lavmodel, lavpartable = lavpartable, :
## lavaan WARNING: the optimizer warns that a solution has NOT been found!
```

```
## Compute CFA for group 12
## Compute CFA for group 13
## Compute CFA for group 14
## Compute CFA for group 15
## Compute CFA for group 16
## Compute CFA for group 17
## Compute CFA for group 18
## Compute CFA for group 19
## Compute CFA for group 20
```

```r
par.env
```

```
## $nu
##        env.1.z     env.2.z      env.3.z
## AT -0.00380100  0.19888606  0.184503447
## BE -0.01780894 -0.03610362  0.148837626
## CH  0.35562809  0.13703068  0.203663682
## CZ -0.13202601 -0.48407184 -0.123536999
## DE  0.15116790  0.11522228  0.224113681
## EE -0.16173882 -0.05409509 -0.219428821
## ES -0.23980207  0.06145417  0.129214174
## FI  0.47267000 -0.04818352 -0.009227077
## FR -0.21551301 -0.09266415  0.108119401
## GB  0.07160703 -0.26288418 -0.100159088
## HU -0.10696530  0.43113728 -0.106775061
## IE -0.14440337 -0.26600599 -0.208756557
## IT -0.14034449 -0.06791908  0.234567775
## LT -0.07871242 -0.16096036 -0.344830664
## NL  0.11769555  0.23075236 -0.131863980
## NO  0.33920229  0.24133236 -0.184151177
## PL -0.33787064  0.01580286  0.032141488
## PT -0.12390489 -0.19353558  0.252706961
## SE  0.56653655  0.26609973 -0.182727192
## SI -0.10412332  0.48297950  0.152794589
## 
## $lambda
##        env.1.z     env.2.z    env.3.z
## AT 0.366378896 0.416570122  0.6019907
## BE 0.511725688 0.547365566  0.3557384
## CH 0.546989687 0.496442011  0.3903070
## CZ 0.554985751 0.733745227  0.6504141
## DE 0.422030164 0.582069644  0.3817274
## EE 0.363889257 0.326850507  0.6572258
## ES 0.345571345 0.548804192  0.5044191
## FI 0.480696324 0.400330553  0.3810049
## FR 0.402922355 0.611851264  0.4036217
## GB 0.437712121 0.685119576  0.4810954
## HU 0.007752346 0.004129361 39.2196109
## IE 0.570436448 0.548172763  0.6867714
## IT 0.400715305 0.693337193  0.4588704
## LT 0.321989539 0.513270709  0.7510091
## NL 0.518470665 0.477207809  0.5366481
## NO 0.567906748 0.424047187  0.4718817
## PL 0.215174623 0.394722347  0.7032612
## PT 0.466262940 0.766183542  0.5120951
## SE 0.556435146 0.364425261  0.4823727
## SI 0.205105231 0.239455075  0.6471300
## 
## $err_var
##      env.1.z   env.2.z       env.3.z
## AT 0.8059733 0.5528352     0.5406610
## BE 0.7251939 0.6829919     0.7451735
## CH 0.5958319 0.4769823     0.8016870
## CZ 0.7737544 1.0701823     0.9604787
## DE 0.6853979 0.4877287     0.8481299
## EE 0.5621019 0.6679707     0.3889200
## ES 0.8722693 0.8227716     0.6286494
## FI 0.5178866 0.6925704     0.6077413
## FR 0.7473996 0.5804562     0.7600567
## GB 0.6924087 0.5783483     0.7168666
## HU 1.0458005 0.7661512 -1537.1226578
## IE 0.7340126 0.9874504     0.6688047
## IT 0.8461694 0.5700815     0.5322671
## LT 0.9670527 0.6667074     0.4055614
## NL 0.7249933 0.5017917     0.8247355
## NO 0.6824385 0.4571760     0.7029381
## PL 0.6534361 0.6781773     0.3748785
## PT 0.9021005 0.9051023     0.7363809
## SE 0.6312888 0.5509960     0.8228221
## SI 0.9464217 0.5335514     0.5685027
## 
## $N
##   AT   BE   CH   CZ   DE   EE   ES   FI   FR   GB   HU   IE   IT   LT   NL   NO 
## 1924 1747 1462 2044 2790 1926 1626 1835 1997 1847 1301 2572 2156 1769 1591 1525 
##   PL   PT   SE   SI 
## 1488 1187 1477 1222 
## 
## $G
## [1] 20
## 
## $I
## [1] 3
## 
## $items
## [1] "env.1.z" "env.2.z" "env.3.z"
## 
## $groups
##  [1] "AT" "BE" "CH" "CZ" "DE" "EE" "ES" "FI" "FR" "GB" "HU" "IE" "IT" "LT" "NL"
## [16] "NO" "PL" "PT" "SE" "SI"
## 
## $CALL
## invariance_alignment_cfa_config(dat = dat.env[, env.z.vars], 
##     group = dat.env[, "cntry"])
```

Configural model does not converge for Hungary. For Czech Republic, one error variance seems to be > 1.


\newpage

### Construct separate configural models for the countries with problems


```r
conf.mod<-"
F.env=~env.1.z+env.2.z+env.3.z
"
```

\newpage

#### Configural model for Hungary


```r
dat.env.HU<-dat.env %>%
  filter(cntry=="HU") %>%
  mutate(env.1.z=na.standardize(inctxff.R),
         env.2.z=na.standardize(sbsrnen.R),
         env.3.z=na.standardize(banhhap.R))


env.conf.HU<-cfa(data=dat.env.HU,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
```

```
## Warning in lav_model_estimate(lavmodel = lavmodel, lavpartable = lavpartable, :
## lavaan WARNING: the optimizer warns that a solution has NOT been found!
```

```r
summary(env.conf.HU)
```

```
## lavaan 0.6-5 did NOT end normally after 10000 iterations
## ** WARNING ** Estimates below are most likely unreliable
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          1301
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                    NA
##   Degrees of freedom                                NA
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate   Std.Err  z-value  P(>|z|)
##   F.env =~                                             
##     env.1.z            0.008       NA                  
##     env.2.z            0.005       NA                  
##     env.3.z           38.450       NA                  
## 
## Variances:
##                    Estimate   Std.Err  z-value  P(>|z|)
##    .env.1.z            0.999       NA                  
##    .env.2.z            0.999       NA                  
##    .env.3.z        -1477.368       NA                  
##     F.env              1.000
```

\newpage

#### Configural model for Czech Republic


```r
dat.env.CZ<-dat.env %>%
  filter(cntry=="CZ") %>%
  mutate(env.1.z=na.standardize(inctxff.R),
         env.2.z=na.standardize(sbsrnen.R),
         env.3.z=na.standardize(banhhap.R))


env.conf.CZ<-cfa(data=dat.env.CZ,
                 model=conf.mod,std.lv=T,auto.fix.first=F)
summary(env.conf.CZ)
```

```
## lavaan 0.6-5 ended normally after 14 iterations
## 
##   Estimator                                         ML
##   Optimization method                           NLMINB
##   Number of free parameters                          6
##                                                       
##   Number of observations                          2044
##                                                       
## Model Test User Model:
##                                                       
##   Test statistic                                 0.000
##   Degrees of freedom                                 0
## 
## Parameter Estimates:
## 
##   Information                                 Expected
##   Information saturated (h1) model          Structured
##   Standard errors                             Standard
## 
## Latent Variables:
##                    Estimate  Std.Err  z-value  P(>|z|)
##   F.env =~                                            
##     env.1.z           0.533    0.031   17.319    0.000
##     env.2.z           0.578    0.032   18.026    0.000
##     env.3.z           0.553    0.031   17.634    0.000
## 
## Variances:
##                    Estimate  Std.Err  z-value  P(>|z|)
##    .env.1.z           0.715    0.033   21.528    0.000
##    .env.2.z           0.665    0.036   18.693    0.000
##    .env.3.z           0.694    0.034   20.319    0.000
##     F.env             1.000
```

