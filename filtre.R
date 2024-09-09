args <- commandArgs(trailingOnly = TRUE)

file= args[1]
dir= args[2]
library(dplyr)
filter <- read.table(paste(dir,"/",file,"_TEDseq_ratio_sort.bed",sep=""))
#filter <- read.table(paste(dir,"/",file,"_TEDseq_ratio_sort_average.bed",sep=""))
filter$V2 <- as.numeric(filter$V2)

filter_filtered <- filter(filter, V4=="-" & V2 +20  >= lead(V2, n=2) )
filter_filtered_bis <- filter(filter, V4=="+" & V2 -20  <= lag(V2, n=2) )
df <- rbind(filter_filtered,filter_filtered_bis)
write.table(df,paste(dir,"/",file,"_TEDseq_ratio_filtered.bed",sep=""),quote=FALSE,col.names = FALSE, row.names = FALSE)
