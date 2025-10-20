library(gmodels)
args<-commandArgs(T)

infile = "all.samples.exp.xls"
outpfx = "all.samples"
scale = "no"
log = "no"

all <- read.table(infile,header = T,row.names=1,sep='\t',check.names=F,quote = "",comment.char = "")
all_num = nrow(all)*ncol(all)
if(log == "yes" && length(all[all<1])/all_num>0.2){
	log = "no";
}

if(log == "yes"){
	all[all<1] = 1
	all = log2(all)
}
sds = apply(all,1,sd)
sds[is.na(sds)] = 0 ## rm NA
all = all[sds!=0,]
tmp = all

scale = ifelse(scale=="yes",T,F)

data <- t(tmp)

data.pca <- fast.prcomp(data,retx=T,scale=scale,center=T)
a <- summary(data.pca)
tmp <- a[4]$importance

pc = as.data.frame(data.pca$x)

## sample pc 
pc_data=data.frame(rownames(pc),pc)
pcs = c(as.numeric(sprintf("%.3f",tmp[2,]))*100)
colnames(pc_data) = c("id",paste(colnames(pc),"(",pcs,"%)",sep=""))
write.table(pc_data,file=paste(outpfx,".PC_data.xls",sep=""),sep="\t", quote=FALSE,row.names=FALSE)


## gene pc 
z<-cbind(rownames(all),data.pca$rotation[,1],data.pca$rotation[,2])
colnames(z)<-c("ID","PC1","PC2")
z=data.frame("id"=rownames(all),"PC1"=data.pca$rotation[,1],"PC2"=data.pca$rotation[,2])
write.table(z,file=paste(outpfx,".PC_2Dcomp.xls",sep=""),sep="\t", quote=FALSE,row.names=FALSE)

