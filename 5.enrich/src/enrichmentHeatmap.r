args<-commandArgs(T)
library(pheatmap)

infile = args[1]
type = args[2]
outpfx = args[3]

#infile = "/Bio/Project/PROJECT/kegg_example/new_enrich/src/t/pv2.go.tmp"
#type = "tt"
#outpfx = "tt"


#infile = paste(inpfx,".tmp",sep="")
outpdf = paste(outpfx,".pdf",sep="")


rawdata = read.table(infile, header = T, sep = "\t", quote ="", check.names = F)
#rawdata = rawdata[1:3000,]
mat = data.frame(rawdata[,3:length(colnames(rawdata))])

rownames(mat) = rawdata[,2]
colnames(mat) = colnames(rawdata)[3:length(colnames(rawdata))]
## matx = max(mat)
## mati = min(mat)

matmp = as.matrix(mat)
matmpx = arrayInd(order(matmp,decreasing=TRUE)[1:1],dim(matmp))
matx = mat[matmpx[1,1],matmpx[1,2]]
matmpi = arrayInd(order(matmp,decreasing=FALSE)[1:1],dim(matmp))
mati = mat[matmpi[1,1],matmpi[1,2]]
annot = data.frame(class=rawdata[,1])
rownames(annot)= rownames(mat)

colors_num = 256
sig_colors_num = colors_num/2
nosig_colors_num = colors_num/2
colors = colorRampPalette(c("green","white","red"))(colors_num)
nosig_colors = colorRampPalette(c("white", "green"))(nosig_colors_num)
sig_colors = colorRampPalette(c("red", "#FFB6C1"))(sig_colors_num)

#type

#mycolor = colors 
#pheatmap(mat,cluster_cols=F, cluster_rows=F, color=mycolor, display_numbers=T, annotation_row=annot, number_format="%.3f",filename=outpdf,cellwidth=25,cellheight=12)
#q()


if (matx > mati)
{
	if (type=="rf.go" || type == "rf.kegg"){
		mycolor = colors
		pheatmap(mat,cluster_cols=F, cluster_rows=F, color=mycolor, display_numbers=T, annotation_row=annot, number_format="%.3f",filename=outpdf,cellwidth=25,cellheight=12)
	}else if(matx > 0.05){
		nosig_breaks = seq(0.0500001,matx,length.out=length(nosig_colors))
		sig_breaks = seq(mati,0.05,length.out=length(sig_colors))
		mycolor = c(sig_colors,nosig_colors)
		pheatmap(mat,breaks = c(sig_breaks,nosig_breaks),legend_breaks=c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1),cluster_cols=F, cluster_rows=F, color=mycolor, display_numbers=T, annotation_row=annot, number_format="%.3f",filename=outpdf,cellwidth=25,cellheight=12)
	}else{
		mycolor = colorRampPalette(c("red", "white"))(256)
		pheatmap(mat,cluster_cols=F, cluster_rows=F, color=mycolor, display_numbers=T, annotation_row=annot, number_format="%.3f",filename=outpdf,cellwidth=25,cellheight=12)
	}
}
