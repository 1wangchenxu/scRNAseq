#!/usr/bin/perl
use warnings;
use strict;
use File::Basename qw(basename);
use FindBin qw($Bin $RealBin);
use Getopt::Long;
use Cwd 'abs_path';
use lib "$RealBin/../";
use Support_Program;
use lib "$RealBin";
use GOVersion;
use Carp qw(carp confess croak);


my %opts = ('k'=>'Genes',"run"=>1,"symbol"=>"none","run_BP"=>1,"run_MF"=>1,"run_CC"=>1);
my @ARGV_bak = @ARGV;
my $src = "$RealBin/src";
my $system_path = "$RealBin/Rscript";
GetOptions(\%opts,"g=s","bg=s","annot=s","op=s", "p=f", "q=f", "ud=s","k=s","a=s","symbol=s","run_BP=s","run_MF=s","run_CC=s","run=s");

#$ENV{'LD_LIBRARY_PATH'}="/Bio/User/luoyue/gcc-4.9.2/lib64:".$ENV{'LD_LIBRARY_PATH'};

my $usage = <<"USAGE";
        Usage:  perl $0 -g <gene.list> -annot <annot> -op <test> [options]
        Notes:  Need R
        version: v2.0

        Options:
                -g      files     gene list
                -annot  files     annotation file
                -symbol files     symbol info
                -op     prefix    output prefix
                -ud     strings   updown table(2 col) or only gene(1 col), default 'nodiff', or 'diff' , 'all'
                -p      float     pvalue, default 0.05
                -q      float     qvalue, default 1

                -k      strings   keyword [Genes]/Proteins

                -BP,-MF,-CC       [1]/0
        eg:
                perl $0 -g gene.list -bg backgroundgene.list -annot annotfile -op test

USAGE

## soft info 
my %apps;
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");

if(defined $opts{a} && defined $opts{bg}){
    print "[Warnings]: found old parameter [bg & pfc], use the old enrich\n";
#    Excu("perl $RealBin/old_version/enrichGO.pl @ARGV_bak");
   Excu("perl $apps{enrichGO_old} @ARGV_bak");
    exit(0);
}

die $usage if(!$opts{g} or !$opts{annot} or !$opts{op});
$opts{p}=$opts{p}?$opts{p}:0.05;
$opts{q}=$opts{q}?$opts{q}:1;
$opts{ud}=$opts{ud}?$opts{ud}:"nodiff";

my @afs;# = qw/P F C/;
push @afs,"P" if($opts{run_BP} == 1);
push @afs,"F" if($opts{run_MF} == 1);
push @afs,"C" if($opts{run_CC} == 1);
die "No any PFC run\n" if(@afs == 0);

$opts{op} = abs_path($opts{op});

my %namespace = (
	"P" => "Biological Process",
	"F" => "Molecular Function",
	"C" => "Cellular Component"
);

my $sep = ",";

##secondLevel 21
my $drawclass = "$RealBin/GOdraw.pl";
my $R_plot = "$RealBin/R_plot.pl";
my $OSGO_conf = "$RealBin/etc/OSGO.$opts{k}.conf";

## get version info 
my ($Rscript,%db_file_info) = GOVersion($opts{annot},\%apps);

## read second level
my %slgo = readSL($db_file_info{sl});

## read diff info  id \t diff
my %gene_ud = readDiff($opts{g},$opts{ud});

## read symbol   id \t symbol
my %symbol = readDiff($opts{symbol},"diff") if(exists $opts{symbol} && $opts{symbol} ne 'none');

open SECL, "> $opts{op}.Level2.xls" or die $!;
open SECL2, "> $opts{op}.Level2.noid.xls" or die $!;
my $title = basename($opts{op}, ".Level2.xls");
my $hh = "GO ID (level1)\tGO Term (level1)\tGO ID (level2)\tGO Term (level2)";
if($opts{ud} eq "nodiff" || $opts{ud} eq 'all'){
    print SECL "$hh\tnumber_of_$title (All)\t$opts{k}_of_$title\n";
    print SECL2 "$hh\tnumber_of_$title (All)\n";
}else{
    print SECL "$hh\tnumber_of_${title} (up)\tnumber_of_${title} (down)\t$opts{k}_of_${title}_up\t$opts{k}_of_${title}_down\n";
    print SECL2 "$hh\tnumber_of_${title} (up)\tnumber_of_${title} (down)\n";
}
my @xls;

Excu("perl $src/annot2osgo.pl $opts{annot} $opts{op}.osgo");
foreach my $t (@afs){
    my %cc = qw/C CC  P BP  F MF/;
    my $type = $cc{$t};
    
    Excu("perl $src/../enrich_xls.pl ${t}enrich -godb $db_file_info{goannot} -godesc $db_file_info{godesc} -bg $opts{annot} -i $opts{g} -symbol $opts{symbol} -outpfx $opts{op}");
    Excu("cut -f1-6 $opts{op}.$t.xls >$opts{op}.$t.xls.tmp && perl $src/../R_plot.pl topGO_graph -Rscript $Rscript -infile $opts{op}.osgo -glist $opts{g} -xls $opts{op}.$t.xls.tmp -type $type -outpfx $opts{op}.$t && rm -f $opts{op}.$t.xls.tmp");

	open PFC, "$opts{op}.$t.xls" or die "$!:$opts{op}.$t.xls";
	my $new_pfc = <PFC>;

    while(<PFC>){
        chomp;
		my @tmp = split /\t/;
        my @genes = split /;/, $tmp[6];
        my @up_gene;
        my @down_gene;
        my @all_genes;

        next if(!exists $slgo{$tmp[0]});
        for my $i(@genes){
            my $gene = $i;
            $gene =~ s/\(.*//;
#            next if($opts{symbol} ne 'none' && !exists $symbol{$gene});
#            $symbol{$i} //= "-";
#            $gene = "$i($symbol{$i})" if(defined $opts{symbol});
            push @all_genes,$i;
            if(!exists $gene_ud{$gene} || $gene_ud{$gene} >0){
                push @up_gene, $i;
            }else{
                push @down_gene,$i;
            }
        }
        if(exists $slgo{$tmp[0]}){
            my $up = @up_gene;
            my $down = @down_gene;
            my $up_gene = @up_gene > 0 ? join "$sep", @up_gene : '--';
            my $down_gene = @down_gene > 0 ? join "$sep", @down_gene : '--';
            if($opts{ud} eq "nodiff" || $opts{ud} eq 'all'){
                print SECL $slgo{$tmp[0]}."\t$up\t$up_gene\n";
                print SECL2 $slgo{$tmp[0]}."\t$up\n";
            }elsif($opts{ud} eq "diff"){
                print SECL $slgo{$tmp[0]}."\t$up\t$down\t$up_gene\t$down_gene\n";
                print SECL2 $slgo{$tmp[0]}."\t$up\t$down\n";
            }
        }
#        $tmp[6] = join(";",@all_genes);
#        $new_pfc .= join("\t",@tmp)."\n";
    }
    close PFC;
#	open OUTPFC, ">$opts{op}.$t.xls" or die $!;
#    print OUTPFC "$new_pfc";
#    close OUTPFC;
    
    push @xls, "$opts{op}.$t.xls";
#    Excu("perl $src/enrich_bar_Gradient.pl $opts{op}.$t.xls $opts{op}.$t 20 GO Q $opts{k}") if($opts{ud} ne 'all');
}

# enrich bar_Gradient
#Excu("perl $src/enrich_bar_Gradient.pl -glist $opts{g} $opts{op}  $opts{op} 20 GO Q $opts{k}") if($opts{ud} ne 'all');

close SECL;
close SECL2;
Excu("rm -f $opts{op}.osgo");
## level2 fig 
my ($dirname,$filename) = $opts{op} =~ /(.*\/)(.*)/;
Excu("perl $R_plot GO_level2 -infile $opts{op}.Level2.xls -conf $OSGO_conf -outdir $dirname -outpfx $filename");
#Excu("cd $dirname;perl $drawclass --prefix $filename $opts{op}.Level2.xls -conf $OSGO_conf;convert -density 300 $opts{op}.GO.level2.bar.svg $opts{op}.GO.level2.bar.png");

## html 
foreach(@xls)
{
	my $i = $_;
	$i =~ s/\.xls$//;
	my ($type) = $i =~ /\.(\w+?)$/;
	my $pi = basename($i);
	open HTML, "> $i.html" or die $!;
	print HTML <<CT;
	<html>
	<head>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>GO Enrichment Analysis</title>
	<style type="text/css">
		a {
			text-decoration: none;
			outline: none;
			hide-focus: expression(this.hideFocus=true);
		}
		a:hover {
			color:#FF0000;
			text-decoration:none;
			outline:none;
			hide-focus: expression(this.hideFocus=true);
		}
		body {
			font-size: 12px;
			font-family: "Microsoft YaHei","微软雅黑","雅黑宋体","新宋体","宋体","Microsoft JhengHei","华文细黑",STHeiti,MingLiu;
			background-color: #FFFFFF;
			padding-left: 8%;
			padding-right: 8%;
		}
		table {
			width: 100%;
			border: 0px;
			border-top: 4px #009933 solid;
			border-bottom: 4px #009933 solid;
			text-align: center;
			border-collapse: collapse;
			caption-side: top;
		}
		th {
			border-bottom: 2px #009933 solid;
			padding-left: 5px;
			padding-right: 5px;
		}
		td {
			padding-left: 5px;
			padding-right: 5px;
		}
		table caption{
			font-weight: bold;
			font-size: 16px;
			color: #009933;
			margin-bottom: 8px;
		}
		
		#bt {
			font-size: 16px;
			position: fixed;
			right: 2%;
			bottom: 5%;
		}
	</style>
	<script type="text/javascript">
	<!--
	function reSize2() {
		try {
			parent.document.getElementsByTagName("iframe")[0].style.height = document.body.scrollHeight + 10;
			parent.parent.document.getElementsByTagName("iframe")[0].style.height = parent.document.body.scrollHeight;
		} catch(e) {}
	}

	preRow = null;
	preColor = null;
	function colorRow(trObj) {
		if (preRow != null) {
			preRow.style.backgroundColor = preColor;
		}
		preRow = trObj;
		preColor = trObj.style.backgroundColor;
		trObj.style.backgroundColor = "FF9900";
	}

	function diffColor(tables) {
		color = ["#FFFFFF", "#CCFF99"];
		for (i = 0; i < tables.length; i++) {
			trObj = tables[i].getElementsByTagName("tr");
			for (j = 1; j < trObj.length; j++) {
				trObj[j].style.backgroundColor = color[j % color.length];
			}
		}
	}

	function markColor(table) {
		trs = table.getElementsByTagName("tr");
			for (i = 1; i < trs.length; i++) {
				if(table.rows[i].cells[6].innerHTML < 0.05){
					//trs[i].style.fontWeight = "500";
					table.rows[i].cells[6].style.color = "#FF0000";
					table.rows[i].cells[6].style.fontWeight = "900";
				}
				if(table.rows[i].cells[5].innerHTML < 0.05){
					table.rows[i].cells[5].style.color = "#FF0000";
				}
			}
	}

	function showPer(tableObj) {
		trObj = tableObj.getElementsByTagName("tr");
		if (trObj.length < 2) {
			return;
		}
		sum1 = trObj[0].cells[3].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
		sum2 = trObj[0].cells[4].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
		/*
		if (trObj[0].cells.length > 4) {
		}
		if (trObj[0].cells.length > 4) {
			trObj[0].cells[2].innerHTML = "DEGs genes with pathway annotation (" + sum1 + ")";
			trObj[0].cells[3].innerHTML = "All genes with pathway annotation (" + sum2 + ")";
		}else{
			trObj[0].cells[2].innerHTML = "All genes with pathway annotation (" + sum1 + ")";
		}
		*/
		for (i = 1; i < trObj.length; i++) {
			trObj[i].cells[3].innerHTML += " (" + (Math.round(trObj[i].cells[3].innerHTML * 10000/ sum1) / 100) + "%)";
			trObj[i].cells[4].innerHTML += " (" + (Math.round(trObj[i].cells[4].innerHTML * 10000/ sum2) / 100) + "%)";
		}
	}


	window.onload = function() {
		setTimeout("reSize2()", 1);
	}
	//-->
	</script>
	
	</head>
		<body>
CT
	
	open FA, $_ or die $!;
	my $head = <FA>; chomp $head;
	my @head = split /\t/, $head;
	my $t1 = "<table><caption>$title GO Enrichment ($namespace{$type})</caption><tr><th>#</th><th>$head[0]</th><th>$head[1]</th><th>$head[2]</th><th>$head[3]</th><th>$head[4]</th><th>$head[5]</th></tr>";
	my $t2 = "<table><caption>$title GO Enrichment ($namespace{$type}) Gene Details</caption><tr><th>#</th><th>$head[0]</th><th>$head[6]</th></tr>";
	my $index = 0;
	while(<FA>)
	{
		chomp;
		my @tmp = split /\t/;
		$index ++;
        last if($index > 5000);
		my $genes = pop @tmp;
#        $genes =~ s/$sep/ /g;
        my @genes = split/[,;]/,$genes;
        @genes = (@genes[1..100],"...") if(@genes > 100);
        $genes = join(" ",@genes);
		$tmp[-1] = sprintf("%.6f", $tmp[-1]);
		$tmp[-2] = sprintf("%.6f", $tmp[-2]);
		$t1 .= "<tr><td>$index</td><td><a href='#gene$index' title='click to view genes' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[1].rows[$index]);'>$tmp[0]</a></td><td>" . (join "</td><td>", @tmp[1..$#tmp]) . "</td></tr>\n";
		$t2 .= "<tr><td>$index</td><td><a href='http://amigo.geneontology.org/amigo/medial_search?q=$tmp[0]' title='click to view GO' target='_blank'>$tmp[0]</a></td><td style=\"text-align: left;\"><a name='gene$index'></a>$genes</td></tr>\n";
	}
	close FA;
	$t1 .= "</table>";
	$t2 .= "</table><script type='text/javascript'>showPer(document.getElementsByTagName('table')[0]);\ndiffColor([document.getElementsByTagName('table')[0], document.getElementsByTagName('table')[1]]);markColor(document.getElementsByTagName('table')[0]);</script>";
	print HTML $t1."\n<br><hr>\n".$t2."\n";
	print HTML <<CT;
		<br><hr>
		<p style="text-align:center;font-weight: bold;font-size: 16px;color: #009933;margin-bottom: 8px;">GO Directed Acycline Graph</p>
		<img src="./$pi.png" style="width:100%">
		<a id="bt" href="#">Back Top</a>
		</body>
	</html>
CT
}

sub version{
    my $file = shift;
    my $dir = shift;

    open FILE,$file or die "$!:$file\n";
    my $version_line = <FILE>;
    chomp($version_line);
    close FILE;
                        ## GO.db     Rscript 
    my %version_info = ("3.3.0" => "$apps{'Rscript.3.3.1.2'}",
                        "3.4.1" => "$apps{'Rscript.3.4.0'}", 
                        "3.5.0" => "$apps{'Rscript.3.4.1'}",
                        "3.6.0" => "$apps{'Rscript.3.5.0'}",
                        "3.8.2" => "$apps{'Rscript.3.6.0'}",
                    );

    my %sl = map{$_=>"$apps{GOinfo}/secondLevel.$_"}keys %version_info;

    my $version = "3.3.0";
    if($version_line =~ /##version (\d.*?)(\(.*)/){
        print "Found version [GO.db $1$2]\n";
        $version = $1;
    }else{
        print "Can't find version info, use default version [GO.db 3.3.0(2016-03-05)]\n";
    }
    if(!exists $version_info{$version} || !exists $sl{$version}){
        my @aa = (%version_info,%sl);
        print STDERR "@aa\n";
        print STDERR "Unknown version [$version]\n";
        exit(1);
    }
    if(!-e $version_info{$version}){
        print STDERR "version [$version] GO.db & Rscript is not exists\n";
        exit(1);
    }
    
    my %hash;
    $hash{sl} = "$apps{GOinfo}/secondLevel.$version";
    $hash{godesc} = "$apps{GOinfo}/goterm.$version";
    $hash{goannot} = "$apps{GOinfo}/goannot.$version";

    return ($version_info{$version},%hash);
}

sub readSL{
    my $sl = shift;
    
    my %slgo;
    open SL, $sl or die "$!:$sl\n";
    while(<SL>){
        chomp;
        my @tmp = split /\t/;
        $slgo{$tmp[1]} = "$tmp[3]\t$tmp[0]\t$tmp[1]\t$tmp[2]";
    }
    close SL;
    
    return %slgo;
}

sub readDiff{
    my $glist = shift;
    my $type = shift;

    my %gene_ud;
    if($type eq "diff"){
        open GENE_LIST, "$glist" or die "$!:$glist";
        while(<GENE_LIST>){
            chomp;
            my @tmp = split /\t/;
            $gene_ud{$tmp[0]} = $tmp[1];
        }
        close GENE_LIST;
    }
    
    return %gene_ud;
}


sub Excu{
    my $cmd = shift;
    print "$cmd\n";

    if($opts{run} == 1){
        my $ret = system($cmd);
        if($ret != 0){
            Error("err$ret\t$cmd");
        }
    }

    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}

