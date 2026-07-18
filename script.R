df<-read.csv("aroma.csv",header = T,sep=";",row.names=1)
library(TeachFreeSortR)
DissTotal<-total_dissim(df)

part1<-dissim_partition(as.matrix(df$S1))
part2<-dissim_partition(as.matrix(df$S2))
p=dim(df)[1]
(sum(part1$V1==part2$V1)-p)/(p*(p-1))

FreeSortR::RandIndex(df$S1,df$S2)
rand_index(df$S1,df$S2)
adjusted_rand_index(df$S1,df$S2)
