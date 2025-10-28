args <- commandArgs(trailingOnly = FALSE) 
Bin = dirname(sub("--file=", "",args[grep("--file=", args)]))

args <- commandArgs(TRUE)
infile = 'stat.xls'
outpfx = 'result'
title = ''
xlab = 'cluster'
ylab = 'Percent (%)'
rmcol = 1
colors = as.character(unlist(strsplit('#FF7777,#759EDD,#EDE447,#42B540,#00468B',",")))
type = 'fill' ## fill bar or count bar 
source(paste(Bin,"ggplot2.r",sep="/"))

dat = ReadTable(infile)
id = rownames(dat)

data = Melt(dat,reverse="yes",rm_col = rmcol)

library(ggplot2)
fig = ggplot(data, aes(x = id, y=value,fill = variable))

##
if(type == "fill"){  
	fig = fig + geom_bar(position = "fill", stat = "identity", width=0.8)
}else{ 
	fig = fig + geom_bar(stat = "identity", width=0.8)
}

fig_width = FigWidth(id)
height = 6
asis_x_size = AxisTextSize(id,width=fig_width)
asis_x_size
mytheme(axis.text.x=asis_x_size,legend.title="no")
fig = fig + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0.001, 0.001))
fig	 = fig + labs(x=xlab, y=ylab,title = title) 
if(length(colors) > 0 && !is.na(colors) && colors[1] != "no"){
	fig = fig + scale_fill_manual(values=colors)
}

fig = fig + xlab_angle(data$id)  ## xlab text angle 

SaveFig(fig,outpfx,width = fig_width,height = height)
