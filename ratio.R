args <- commandArgs(trailingOnly = TRUE)

file= args[1]
dir=args[2]
depth <- read.table(paste(dir,"/",file,"_depth.bed", sep=""))
start <- read.table(paste(dir,"/",file,"_start_site.bed", sep=""))
stop <- read.table(paste(dir,"/",file,"_stop_site.bed", sep=""))
colnames(depth) <- c("col1","col2","col3")
colnames(start) <- c("col1","col2","col3","col4")
colnames(stop) <- c("col1","col2","col3","col4")
mergefile <- merge(depth,start,  by.x = c("col1","col2") , by.y = c("col1","col2"), all = TRUE)
colnames(mergefile) <- c("Chr","start","depth","strand","nbr_start")
mergef <- merge(mergefile,stop, by.x = c("Chr","start") , by.y = c("col1","col2"), all = TRUE)
mergef[is.na(mergef)] <- 0  
mergef_zero <- mergef[!(mergef$strand == 0 & mergef$col3 == 0),]
ratio <- cbind(mergef_zero,(mergef_zero$nbr_start+mergef_zero$col4)/mergef_zero$depth)
write.table(ratio,paste(dir,"/",file,"_ratio.bed",sep=""),quote=FALSE,row.names = FALSE, col.names=FALSE)
