## .libPaths rm /home

path0 = .libPaths()
path = path0[grep("^/home/",path0,invert=T)]
.libPaths(path)
## theme_get()

#rgb(t(col2rgb("red")),max=255)
#"#FF0000"


## read table 修改某些默认值
ReadTable = function(file,header=TRUE,row.names=1,sep="\t",stringsAsFactors = F,check.names=FALSE,quote="",trans_table = "no"){
	data = read.table(file, header=header, row.names = row.names, sep=sep,stringsAsFactors = stringsAsFactors, check.names=check.names, comment.char = "",quote=quote)
	if(trans_table == "yes") data = data.frame(t(data),check.names=F)
	data
}

## write table 
MyWriteTable = function(data,file,col.names=TRUE){
	write.table(data,file=file,sep="\t",quote=FALSE,row.names=FALSE,col.names=col.names)
}

## default theme 
mytheme =function(axis.text.size=NA,axis.text.x=NA,panel.grid = "no",panel.border="yes",axis.text="yes",legend.title = "yes",axis.ticks="yes",legend.position="yes",Cairo="no", extrafont="yes", font_use="Arial"){

	my_theme = theme_bw()
	panel_border_size = 1
	if(extrafont == "yes"){
		library(extrafont)
		library(extrafontdb)
		library(Rttf2pt1)
		# extrafont::font_import() to import fonts 
		# extrafont::fonts() to check what fonts are imported 
		if( font_use %in% fonts() ) {
			cat("mytheme: use font ",font_use,"\n")
			my_theme <- my_theme + theme(text = element_text(family = font_use))
		}
	}else if(Cairo == "yes"){
		ggsaveR = paste0(Bin,"/ggsave.R")
		library(Cairo)
		Error("Cairo" %in% installed.packages(),"Please install cairo")
		CairoFonts(
			regular = "Arial:style=Regular",
			bold = "Arial:style=Bold",
			italic = "Arial:style=Italic",
			bolditalic = "Arial:style=Bold Italic,BoldItalic",
		)
		panel_border_size = 0.6
		source(ggsaveR)
	}else{
		library(showtext)
		showtext_auto(enable = TRUE)
#	font_add("Arial", regular = "/Bio/User/liuyubin/pipeline/GeneralPlot/fonts/tff/msttcore/arial.ttf")
		regular = paste(Bin,"/../conf/arial.ttf",sep="")
		font_add("Arial", regular = regular)
		my_theme = my_theme + theme(text = element_text(family = "Arial"))
	}
	my_theme = my_theme +  
	theme(

		panel.border = element_rect(color = "#000000", size = panel_border_size),

		axis.text = element_text(color = "#000000", size = 14),
#		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#		axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 16, face = "bold"),
#		axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
#		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
#		axis.ticks = element_line(color = "#000000", size = 0.8),
#		axis.ticks.length = unit(0.15, 'cm'),

#		legend.title = element_blank(),
		legend.text = element_text(size = 12),
#		legend.key = element_blank(),
#		legend.key.height = unit(5, "mm"),
#		legend.key.width = unit(5, "mm"),
#		legend.direction = "vertical",
#		legend.justification = "center",

		plot.title = element_text(color = "#000000", size = 20, face = "bold", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)
	if(!is.na(axis.text.x)){
		my_theme = my_theme + theme(axis.text.x=element_text(size=axis.text.x))
	}
	if(!is.na(axis.text.size)){
		my_theme = my_theme + theme(axis.text=element_text(size=axis.text.size))
	}
	if(panel.grid == "no"){
		my_theme = my_theme + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
	}
	if(panel.border == "no"){ ## 边框
#		my_theme = my_theme + theme(panel.background = element_blank(),panel.grid=element_blank())
		my_theme = my_theme + theme(panel.border = element_blank())
	}
	if(axis.text == "no"){ ## 坐标轴文字
		my_theme = my_theme + theme(axis.text=element_blank(), axis.ticks=element_blank())
	}
	if(legend.title == "no"){ ## legend标题
		my_theme = my_theme + theme(legend.title=element_blank())
	}
	if(legend.position == "none"){ ## legend
		my_theme = my_theme + theme(legend.position="none")
	}
	if(axis.ticks == "no"){ ## 坐标轴刻度
		my_theme = my_theme + theme(axis.ticks = element_blank())
	}

	theme_set(my_theme)
#	my_theme
}


## 字体大小
AxisTextSize = function(sample_name,width=8){
	name = as.character(unique(sample_name))
	name_len = sum(nchar(name))
	name_mean = name_len / length(name)
	
	if(name_mean < 10){
		size = 14
	}else{
		size = 10
	}
	sample_num = length(name)
	if(sample_num >= 40){
		## 
		size = width*0.7*72/sample_num
#		size = width*8*0.9/sample_num
#		size = rel(size)
#		cat(sample_num,size,"\n")
	}
	
	size
}

##柱形图 添加文字大小 
GeomTextSize = function(value,sample_num=10,width = 8,fontsize=4,legend=T){
	name = as.character(value)
	name_len = sum(nchar(name))
	max_name_len = max(nchar(name))
	name_mean = name_len / length(name)
	
	name_font_len = fontsize * 1/72 * name_len * 1 ## size * pt2in * len * 大小写比例,数字为1
	max_name_font_len = fontsize * 1/72 * max_name_len * 1 ## size * pt2in * len * 大小写比例,数字为1
	
	## 作图区域
	draw_width = ifelse(legend,width-3,width-2)

	if(max_name_font_len > draw_width/sample_num){
		##
		fontsize = draw_width*72/sample_num/max_name_len*0.8
#		cat(length(value),fontsize,"\n")
	}
#	cat("GeomTextSize",max_name_font_len,name_font_len,width,draw_width,length(value),sample_num,fontsize,"\n")

	fontsize
}

AxisTextSize2 = function(sample_name,width=8){
	name = as.character(unique(sample_name))
	name_len = sum(nchar(name))
	name_mean = name_len / length(name)
	
	if(name_mean < 10){
		size = 14
	}else{
		size = 10
	}
	sample_num = length(name)
	if(sample_num >= 100){
		if(sample_num >= 1000){
			size = 0.08 
		}else{
			size = switch(as.integer(sample_num/100),
				0.40, #1 - 0.0025 * sample_num,
				0.28,
				0.22,
				0.20,
				0.16,
				0.14,
				0.10,
				0.10,
				0.08)
		}
#		size = width/8*size
		size = rel(size)
	}

	size
}

## 根据名字长度设置角度axis.text
# xlab text names 
xlab_angle = function(sample_name,width = 7,fontsize = 14){
	name = as.character(unique(sample_name))
	name_len = sum(nchar(name))
	name_mean = name_len / length(name)

	name_font_len = fontsize * 1/72 * name_len * 0.6 ## size * pt2in * len * 大小写比例,长度不一

	cut = c(0.6,2)*width  ## 文字 估算占 图片的0.6  放松一些, 0.7 名字会连在一起
#	cat(cut,name_len,name_font_len,"\n")
#	c(width*4*14/fontsize,width*13*14/font)
#	if(name_len < cut[1]){
	if(name_font_len < cut[1]){
		angle = 0   ## 水平
		hjust = 0.5 ## 左右 居中
		vjust = 0   ## 貌似没用
	}else if(name_font_len < cut[2]){
		if(name_mean > 15){ ## 名字过长，角度太大不好看
			angle = 20
			hjust = 1
			vjust = 1
		}else{
			angle = 45
			hjust = 1  ## 左右, 最右
			vjust = 1  ## 上下, 最上
		}
	}else{
		angle = 90  ## 竖直
		hjust = 1   ## 上下, 最上
		vjust = 0.5 ## 左右, 中间
	}
	theme(axis.text.x= element_text(angle=angle,hjust = hjust,vjust = vjust))
}
#theme(axis.text.x= element_text(angle=angle,hjust=hjust))

FigWidth = function(sample_name,width=8){
	sample_num = length(sample_name)
	if(sample_num >= 800){
		width = 30
	}else if(sample_num >= 500){
		width = 25
	}else if(sample_num >= 100){
		width = switch(as.integer(sample_num/100),
			11,14,17,20)
	}else{
		if(sample_num >= 50){
			width = 10
		}else if(sample_num >= 20){
			width = 9
		}
	}

	width
}


## LegendWidth
## legend多的时候，，legend会多列
LegendWidth = function(sample_name,fontsize = 14){
	name = as.character(unique(sample_name))
	max_name_len = max(nchar(name))
#	name_mean = name_len / length(name)
	name_font_len = fontsize * 1/72 * max_name_len * 0.6 ## size * pt2in * len * 大小写比例,长度不一

	col_num = as.integer((length(sample_name)-1)/20)+1

	width = col_num*(name_font_len+0)
#	cat(length(sample_name),name_font_len,col_num,width,"\n")
	
	width
}


## save pdf & png 
# ggplot2 fig, outpfx 
SaveFig = function(fig = fig,outpfx = outpfx,width=7,height=7,eps="no"){
	ggsave(file = paste(outpfx,".png",sep=""),fig,width=width,height=height,limitsize = FALSE)
	ggsave(file = paste(outpfx,".pdf",sep=""),fig,width=width,height=height,limitsize = FALSE)
	if(eps == "yes"){
		ggsave(file = paste(outpfx,".eps",sep=""),fig,width=width,height=height,limitsize = FALSE)
	}

}

## png() pdf() save
SaveFig2 = function(outpfx = outpfx,width=7,height=7,type="pdf"){
	if(type == "pdf"){
		pdf(paste(outpfx,".pdf",sep=""),width=width,height=height)
	}else if(type == "png"){
		png(paste(outpfx,".png",sep=""),width=width,height=height,units = 'in', res = 300)
	}else if(type == "dev.off"){
		dev.off()
	}
}

## like melt, no need library 
Melt = function(dat,remain=NA,rm_col=NA,rr = "no",reverse="no",value_type="numeric",parenthese="no"){

	Error(is.data.frame(dat),"No data frame input for Melt")

	id=rownames(dat)
	if(parenthese == "no") colnames(dat) = gsub("\\(.*","",colnames(dat),perl=TRUE)
	colnames(dat) = gsub("_abundance$","",colnames(dat),perl=TRUE) ## smallRNA 
	colnames(dat) = gsub("_unique$","",colnames(dat),perl=TRUE)  ## smallRNA 
	coln = colnames(dat)
#	print(coln)
#	if(length(coln) < 3){
#		data = data.frame(id = id,value = dat[,1],variable=id)
#		return(data)
#	}

	allcols = 1:ncol(dat)
	if(!is.na(remain)){
		allcols = allcols[-remain]
	}else if(!is.na(rm_col) && rm_col[1] > 0){
		allcols = allcols[-rm_col]
	}
#	print(rm_col)
#	print(allcols)

##  change data format
	melt_col = ifelse(rr == "yes",4,3)

	tmp = rep(0,nrow(dat)*length(allcols)*melt_col)
	dim(tmp) = c(nrow(dat)*length(allcols),melt_col)
	nn = 1

	for (i in allcols){
		tmp[nn:(nn+nrow(dat)-1),1] = id
		tmp[nn:(nn+nrow(dat)-1),2] = dat[,i]
		tmp[nn:(nn+nrow(dat)-1),3] = rep(coln[i],nrow(dat))
		if(rr == "yes"){
			tmp[nn:(nn+nrow(dat)-1),4] = dat[,remain]
		}
		nn = nn + nrow(dat)
	}
	colname = coln[allcols]
	if(reverse != "no") colname = rev(colname)

	data = data.frame(tmp)
	data[,2] = gsub("\\s.*","",data[,2],perl=TRUE)
	if(value_type == "numeric"){
		data[,2] = as.numeric(data[,2])
	}
#	data[,2] = as.numeric(gsub("\\s.*","",data[,2],perl=TRUE))
#	data[,2] = as.numeric(as.character(data[,2])) ## col 2 is number
	data[,3] = factor(data[,3],levels = colname)
	data[,1] = factor(data[,1],levels = id)
	colnames(data) = c("id","value","variable")
	
	if(!is.na(remain)){
		for (i in 1:length(remain)){
			colnames(data)[i+3] = coln[remain[i]]
#			data[coln[i]] = dat[,i]
		}
	}

	data
}

## 根据长度换行
sepline=function(ss,len){
	splits = unlist(strsplit(ss,split=" "));
#	print(splits)
	if(nchar(ss,type="bytes") > len & length(splits)>1){ 
		all = 0
		out = c()
		i = 1
		while(all < nchar(ss,type="bytes")/2){
			out[i] = splits[i]
			all = all+nchar(splits[i],type="bytes") + 1
#			cat(all,nchar(ss,type="bytes"),"\n")
			i = i+1
		}
#		i = i-1
		if(i > length(splits)){
			i = i-1
			out = out[-length(out)]
		}
		ss = paste(out,collapse=" ")
		other = paste(splits[i:length(splits)],collapse=" ")
#		cat("=>",i,length(splits),ss,"\n",other,"\n")
		ss = paste(ss,"\n",other,sep="")
	}
	ss
}

GetSepMaxLen = function(x){
	all = paste("\n",x)
	all_sep = unlist(strsplit(all,"\n"))
	len = nchar(all_sep,type="bytes")

	max(len)
}

## abcde,2 => ab...
CommentWord = function(ss,len){
	if(nchar(ss) > len){
		ss = substr(ss,1,len)
		ss = paste0(ss,"...")
	}
	ss
}

Number = function(x){
	if(x == "NULL"){
		NULL 
	}else{
		as.numeric(x)
	}
}

GetPFCData = function(pfc,class_list){
	pfc_info = c("P"="Biological Process","F"="Molecular Function","C"="Cellular Component","all"="all")
	class_list == pfc_info[pfc]
}

MaxMin = function(x,ratio = 1,ref_min_max = "none"){
	values = unlist(strsplit(ref_min_max,","))
	if(length(values) == 3) ratio = values[3]
	max_v = ifelse(max(x) > 0 , max(x) * ratio, max(x) / ratio)
	min_v = ifelse(min(x) > 0 , min(x) / ratio, min(x) * ratio)
	if(ref_min_max != "none"){
		if(values[1]!="none") min_v = as.numeric(values[1])
		if(values[2]!="none") max_v = as.numeric(values[2])
	}

	c(min_v,max_v)
}

Sample2Colors = function(x,color_list = c("red","green","blue","cyan","yellow","mediumpurple","orange","purple","pink","gray","wheat","brown","darkgreen","greenyellow","black","chocolate")){
	if(color_list[1] == "none"){color_list = c("red","green","blue","cyan","yellow","mediumpurple","orange","purple","pink","gray","wheat","brown","darkgreen","greenyellow","black","chocolate")}
	color_list[x]
}

# scale_x_discrete

## PCA point size 
PointSize = function(sample.num = 1,psize0=0,tsize0=0){
	if(sample.num <= 20) {
		psize = 6
		tsize = 4
	} else if(sample.num <= 40) {
		psize = 5
		tsize = 3.5
	} else if(sample.num <= 80) {
		psize = 3
		tsize = 3
	} else {
		psize = 2
		tsize = 2
	}
	if(psize0 != 0){
		psize = psize0
	}
	if(tsize0 != 0){
		tsize = tsize0
	}

	list(point_size=psize,text_size=tsize)
}

Error = function(TF,word){
	if(!TF){
		cat("Error:",word,"\n")
		q(status=1)
	}
}

MyScale = function(data,scale="none"){ # row  column
	if(scale == "row"){
		if(ncol(data) == 1){
			cat("=>Attention: Can't scale by row with only one row!! So will not scale the data!!\n")
		}else{
			data = t(scale(t(data)))
		}
	}else if(scale == "column" || scale == "col"){
		if(nrow(data) == 1){
			cat("=>Attention: Can't scale by column with only one column!! So will not scale the data!!\n")
		}else{
			data = scale(data)
		}
	}else if(scale == "row0"){
		data = t(scale(t(data)))
		data[is.nan(data)] = 0
	}

	data
}

Sample2Group = function(sample_list,group_info){
	out = data.frame(sample=sample_list,group=rep("others",length(sample_list)))
	for(i in 1:length(sample_list)){
		index = which(group_info[,1] %in% sample_list[i])
		if(length(index) > 0){
			out[i,2] = group_info[index,2]
		}
	}
	out
}

