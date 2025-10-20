#!/usr/bin/perl -w

use Getopt::Long;
use FindBin qw($Bin $RealBin);
use File::Basename qw/basename dirname/;

my ($indir, $help, $key_word);
GetOptions("indir:s" => \$indir, "help|?" => \$help, "k:s"=>\$key_word);

if (!defined $indir || defined $help) {
	print STDERR << "USAGE";
description: generate DO html files
usage: perl $0 [options]
options:
	-indir *: input directory, containing *.do.xls files
	-help|?: print help information
USAGE
	exit 1;
}

if (!-d "$indir") {
	print STDERR "directory $indir not exists\n";
	exit 1;
}

#@files = glob("$indir/*.do.xls");
my @files = `find $indir -name "*.do.xls"`;

foreach(@files){
    chomp($_);
	my $i = $_;
	$i =~ s/\.xls$//;
	my $pi = basename($i,"\.do");
	open HTML, "> $i.html" or die $!;
	print HTML <<CT;
	<html>
	<head>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>DO Enrichment Analysis</title>
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
	my $t1 = "<table><caption>$pi DO Enrichment</caption><tr><th>#</th><th>$head[0]</th><th>$head[1]</th><th>$head[3]</th><th>$head[4]</th><th>$head[5]</th><th>$head[6]</th></tr>";
	my $t2 = "<table><caption>$pi DO Gene Details</caption><tr><th>#</th><th>$head[0]</th><th>$head[-1]</th></tr>";
	my $index = 0;
	while(<FA>)
	{
		chomp;
		my @tmp = split /\t/;
		$index ++;
		my $genes = pop @tmp;
		my @genes = split /;/, $genes;
		if (@genes > 100)
		{
			$genes = join(" ",@genes[0 .. 99]);
			$genes .= " (the number of this term's genes more than 100, full information in file)";
		}
		else
		{
			$genes = join(";",@genes);
		}
		$tmp[-1] = sprintf("%.6f", $tmp[-1]);
		$tmp[-2] = sprintf("%.6f", $tmp[-2]);
		$t1 .= "<tr><td>$index</td><td><a href='#gene$index' title='click to view genes' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[1].rows[$index]);'>$tmp[0]</a></td><td>" . (join "</td><td>", @tmp[1,3..$#tmp]) . "</td></tr>";
		$t2 .= "<tr><td>$index</td><td><a href='$tmp[2]' title='click to view DO' target='_blank'>$tmp[0]</a></td><td style=\"text-align: left;\"><a name='gene$index'></a>$genes</td></tr>";
	}
	close FA;
	$t1 .= "</table>";
	$t2 .= "</table><script type='text/javascript'>showPer(document.getElementsByTagName('table')[0]);\ndiffColor([document.getElementsByTagName('table')[0], document.getElementsByTagName('table')[1]]);markColor(document.getElementsByTagName('table')[0]);</script>";
	print HTML $t1."\n<br><hr>\n".$t2."\n";
	print HTML <<CT;
		<p><br /></p><p><br /></p>
		<a id="bt" href="#">Back Top</a>
		</body>
	</html>
CT
}
