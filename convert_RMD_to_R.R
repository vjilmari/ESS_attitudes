library(knitr)


purl(input="All_analyses.Rmd",
     output="All_analyses.R",documentation = 2)

purl(input="Alignment_attempt.Rmd",
     output="Alignment_attempt.R",documentation = 2)
