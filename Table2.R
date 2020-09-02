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

write.csv2(tab2,"tab2.csv")
