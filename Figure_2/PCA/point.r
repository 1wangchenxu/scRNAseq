args <- commandArgs(trailingOnly = FALSE) 
Bin = dirname(sub("--file=", "",args[grep("--file=", args)]))
source(paste(Bin,"ggplot2.r",sep="/"))

args<-commandArgs(T)
infile = 'all.samples.PC_data.xls' 
outpfx = 'all.samples.PCA'
col_list = as.numeric(unlist(strsplit('1,2',','))) ## 1,2 
xlab = NA
ylab = NA
title = 'PCA'
group_info = 'group.list'
vjust=-1.65
psize=0
tsize=0
colors = unlist(strsplit('#D32421,#2771A7',","))
isellipse = 'no'

library(ggplot2)
dat0 = ReadTable(infile)
dat = dat0[,col_list]

if(xlab == "NA") xlab = colnames(dat)[1]
if(ylab == "NA") ylab = colnames(dat)[2]
if(group_info == "none") group_info = "NA"

samples = rownames(dat)
colnames(dat)[1] = "x"
colnames(dat)[2] = "y"

ratio = ifelse(isellipse=="yes",3,1.2)
xmin_max = MaxMin(dat[,1],ratio = ratio)
ymin_max = MaxMin(dat[,2],ratio = ratio)

if(group_info != "NA"){
	classify = read.table(group_info,header = F,sep='\t',check.names=F,quote = "",row.names=1,na.strings="")
	if(ncol(classify) == 0){
		group_info = "NA"
	}
#	head(classify)
}
dat$names = samples
#colnames(dat)
if("group" %in% colnames(dat0)){
	groups = factor(dat0$group,levels = unique(dat0$group))
	legends = 'right'
}else if(group_info == "NA"){
	groups = factor(samples,levels=samples)
	legends = 'none'
}else{
	group_name = classify[samples,1]
	if(length(classify[is.na(group_name),1])>0){
		cat("No samples in group:",samples[is.na(group_name)],"\n")
		q(status=1)
	}
	groups = factor(group_name,levels=unique(group_name))
	legends = 'right'
}

dat$group = groups
length = length(unique(dat$group))

if(!exists('vjust') || is.na(vjust)){
    vjust = -0.5
}
size_info = PointSize(sample.num=length(samples),psize0=psize,tsize0=tsize)
#size_info
psize = size_info$point_size
tsize = size_info$text_size

pca=ggplot(dat,aes(x,y))
if(group_info != "NA"){
	if(length>6){
		pca=pca+geom_point(size=psize,aes(color=group))
	}else{
		pca=pca+geom_point(size=psize,aes(shape=group,color=group))
	}
}else{
	pca=pca+geom_point(size=psize,aes(color=group))
}
mytheme()
pca=pca+labs(x=xlab,y=ylab,title=title)+
	geom_hline(yintercept=0,linetype=4,color="grey") +
	geom_vline(xintercept=0,linetype=4,color="grey")
if(colors[1] != 'none'){
	pca = pca + scale_color_manual(values = colors)
}

pca = pca + guides(color = guide_legend(override.aes = list(size = 4)))
if(isellipse == "yes"){
	pca = pca + stat_ellipse(aes(fill=group),type="norm",geom="polygon",alpha=0.2,color=NA)
}else{
	pca = pca + xlim(xmin_max)+ylim(ymin_max)
}
n=30
if(length > 30){
	pca2=pca+guides(col = guide_legend(nrow=n))
	width=ceiling(length/n)*1.5+10
}else{
	pca2 = pca
	width = 10
}

SaveFig(fig=pca2,outpfx=paste(outpfx,".nosampleid",sep=""),width = width,height=8)

pca=pca+theme(legend.position=legends)
pca=pca+geom_text(aes(label=names),size=tsize,vjust= vjust)
width = 10

SaveFig(fig=pca,outpfx=outpfx,width=width,height=8)
