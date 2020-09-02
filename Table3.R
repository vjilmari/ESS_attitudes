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

write.csv2(tab3,"tab3.csv")
