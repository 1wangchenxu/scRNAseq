#! /usr/bin/perl
use utf8;
use strict;
use warnings;

die "perl $0 <kegg.path.xls> <the first number> <Q/P>\n" if(@ARGV < 3);

my ($out_perfix) = $ARGV[0] =~ /(\S+)\.xls$/;
my $name_key = "Gene";
if(defined $ARGV[3]){
        $name_key = $ARGV[3];
}

my($n, %ra, %qv);
open FA, $ARGV[0] or die $!;
open TMP, "> $out_perfix.tmp" or die $!;
<FA>;
while(<FA>)
{
	chomp;
	$n ++;
	last if($n > $ARGV[1]);
	my @tmp = split /\t/;
	my $ratio = $tmp[3] / $tmp[4];
	my $path = $tmp[2];
	$path =~ s/'/\'/g;
	$path =~ s/"/\"/g;
	$ra{$path} = $ratio;
	if($ARGV[2] eq "Q"){
		$qv{$path} = "$tmp[6]\t$tmp[3]";
#		$order{$path} = $tmp[6];
	}
	elsif($ARGV[2] eq "P"){
		$qv{$path} = "$tmp[5]\t$tmp[3]";
#		$order{$path} = $tmp[6];
	}
	else{
		die "perl $0 <kegg.path> <the first number> <Q/P>\n";
	}
	print TMP "$path\t$ra{$path}\t$qv{$path}\n";
}
close FA;
close TMP;

open RCMD, "> $out_perfix.r" or die $!;
print RCMD "
library(ggplot2)
mat <- read.table(\"$out_perfix.tmp\", sep = \"\t\", check.names = 0, header = F, quote=\"\")
sepline=function(ss,len){
    if(nchar(ss) > len){ 
        splits = unlist(strsplit(ss,split=\" \"));
        all = 0
        out = ''
        i = 1
        while(all < nchar(ss)/2){
            if(i == 1){
                out = splits[i]
            }else{
                out = paste(out,splits[i],collapse=\" \",sep=\" \")
            }
            all = all+nchar(splits[i]) + 1
            i = i+1
        }
        other = paste(splits[i:length(splits)],collapse=\" \")
        ss = paste(out,\"\\n\",other,sep=\"\")
    }
    ss
}
mat\$V1 = as.character(mat\$V1)
for(i in 1:nrow(mat)){
            mat\$V1[i] = sepline(mat\$V1[i],60)
}

matmp = as.matrix(mat\$V3)
matmpx = arrayInd(order(matmp,decreasing=TRUE)[1:1],dim(matmp))
matx = matmp[matmpx[1,1],matmpx[1,2]]
matmpi = arrayInd(order(matmp,decreasing=FALSE)[1:1],dim(matmp))
mati = matmp[matmpi[1,1],matmpi[1,2]]
porder = factor(mat\$V1,levels=rev(mat\$V1))
p <- ggplot(mat, aes(mat\$V2, porder))
if (matx > mati)
{
#	porder = factor(as.integer(rownames(mat)),labels=mat\$V1)
#    porder = factor(mat\$V1,levels=rev(mat\$V1))
#	p <- ggplot(mat, aes(mat\$V2, porder))
	fig = p + geom_point(aes(size = mat\$V4,colour = mat\$V3)) + scale_colour_continuous(\"$ARGV[2]Value\", low=\"red\", high = \"forestgreen\")+guides(colour = guide_colorbar(order=2),size = guide_legend(order=1))
}else{
	$ARGV[2]Value = as.character(matx)
#	p <- ggplot(mat, aes(mat\$V2, mat\$V1))
	fig = p + geom_point(aes(size = mat\$V4,colour = $ARGV[2]Value))+guides(colour = guide_legend(order=2),size = guide_legend(order=1))
}
fig = fig + scale_size(\"${name_key}Number\") + labs(title = \"Top $ARGV[1] of Pathway Enrichment\", x = \"RichFactor\", y = \"Pathway\") + theme_bw() 
#    theme(axis.text=element_text(size=8),axis.title=element_text(size=10,face=\"bold\"))
ggsave(fig,file = \"$out_perfix.png\", dpi = 300)
ggsave(fig,file = \"$out_perfix.pdf\")
";

`Rscript $out_perfix.r`;
`rm $out_perfix.tmp $out_perfix.r Rplots.pdf -rf`;
