
=pod
description: generate pathway html files
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090806
modified date: 20100205, 20100127, 20091204, 20091201, 20091010, 20090814, 20090807
=cut
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin $RealBin);
use File::Basename qw/basename dirname/;

my $version = "1.1";
my ($indir, $help);
my $nomap_die = 1;
my $key_name = "gene"; # metabolites  protein mRNA target_genes
my %type_map = ('gene'=>['genes','K_id'],'metabolite'=>['metabolites','C_id'],'protein'=>['proteins','K_id'],
                'mRNA'=>['mRNAs','K_id'],'target' => ['target_genes','K_id'],
                'source_gene'=>['genes','K_id']);
$type_map{Proteins}  = $type_map{protein};
$type_map{Metabolites} = $type_map{metabolite};
$type_map{Genes} = $type_map{gene};
$type_map{target_gene} = $type_map{target};
#my $key_id = "K_id";
my $ko;
my $run_symbol = "none";
my $diff_type = "nodiff";

GetOptions("indir:s" => \$indir, "help|?" => \$help,"k:s"=> \$key_name,"ko:s"=>\$ko,"nomap_die=s"=>\$nomap_die,"symbol=s"=>\$run_symbol,"diff_type=s");

if (!defined $indir || defined $help) {
	print STDERR << "USAGE";

    description: generate pathway html files
    usage: perl $0 [options]
    options:
        -indir *: input directory, containing *.path files
        -help|?: print help information
        -k     : key name for display, [gene],metabolite,protein,target

USAGE
	exit 1;
}

if (!-d "$indir") {
	print STDERR "directory $indir not exists\n";
	exit 1;
}

#$key_name .= "s";
die "Unknown key name [$key_name]\n"if(!exists $type_map{$key_name});
my $key_id = $type_map{$key_name}[1];
$key_name = $type_map{$key_name}[0];
#$key_id = "C_id" if($key_name=~/metab/i);



my %ko_desc;
my $ko_desc = "$RealBin/ko";
GetKoDesc($ko_desc,\%ko_desc);


#my @files = glob("$indir/*.path.xls");

my @files = `find $indir -name "*.path.xls"`;
#print join("=",@files)."\n";die;

for my $i (0 .. $#files) {
    chomp($files[$i]);
    my $dir = dirname $files[$i];
    my $name = (split /[\\\/]/, $files[$i])[-1];
    $name =~ s/\.path\.xls$//;
    my $htmlFile = $files[$i];
    $htmlFile =~ s/path\.xls$/htm/;

    my $code = <<HTML;
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<title>$name</title>
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
	#backtop{
		font-size: 16px;
		position: fixed;
		bottom: 5%;
		right: 2%;
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
tr=null;
trc=null;
function colorRow(trObj,num,count) {
    if (preRow != null) {
        for (var i = 0; i < trc; i++) {
            preRow.rows[tr].style.backgroundColor = preColor[i];
            tr += 1;
        };

    }
    preRow = trObj;
//    count +=1
    preColor=new Array();
    tr=num;
    trc=count;
    for (var i = 0; i < count; i++) {
        preColor[i] = trObj.rows[num].style.backgroundColor;
        trObj.rows[num].style.backgroundColor = "FF9900";
        num++;
    }
}

function colorRow_old(trObj) {
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
			if(table.rows[i].cells.length > 4){
				if(table.rows[i].cells[5].innerHTML < 0.05){
					//trs[i].style.fontWeight = "500";
					table.rows[i].cells[5].style.color = "#FF0000";
					table.rows[i].cells[5].style.fontWeight = "900";
				}
				if(table.rows[i].cells[4].innerHTML < 0.05){
					table.rows[i].cells[4].style.color = "#FF0000";
				}
			}
		}
}

function showPer(tableObj) {
	trObj = tableObj.getElementsByTagName("tr");
	if (trObj.length < 2) {
		return;
	}
	sum1 = trObj[0].cells[2].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
	if (trObj[0].cells.length > 4) {
		sum2 = trObj[0].cells[3].innerHTML.replace(/^.*\\(([\\d]+)\\).*\$/, "\$1");
	}
	if (trObj[0].cells.length > 4) {
		trObj[0].cells[2].innerHTML = "Candidate $key_name with pathway annotation (" + sum1 + ")";
		trObj[0].cells[3].innerHTML = "All $key_name with pathway annotation (" + sum2 + ")";
	}else{
		trObj[0].cells[2].innerHTML = "All $key_name with pathway annotation (" + sum1 + ")";
	}
	for (i = 1; i < trObj.length; i++) {
		trObj[i].cells[2].innerHTML += " (" + (Math.round(trObj[i].cells[2].innerHTML * 10000/ sum1) / 100) + "%)";
		if (trObj[0].cells.length > 4) {
			trObj[i].cells[3].innerHTML += " (" + (Math.round(trObj[i].cells[3].innerHTML * 10000/ sum2) / 100) + "%)";
		}
	}
}

window.onload = function() {
	setTimeout("reSize2()", 1);
}
//-->
</script>
HTML
	$code .= "</head>\n<body>\n";
	open IN, "< $files[$i]" || die $!;
	chomp(my $content = <IN>);
	my @temp = split /\t/, $content;
	shift @temp; shift @temp;
	my $pre = shift @temp;
	pop @temp;
	my $gene = pop @temp;
    my $gene_list_title = "$key_name";
#	if(scalar(@temp) > 2){## diff 
    if($diff_type eq 'diff'){
		$gene_list_title = "Differentially expressed $key_name";
	}
	$code .= "<table><caption>$name Pathway Enrichment</caption><tr><th>#</th><th>" . substr($pre, 0) . "</th><th>" . (join "</th><th>", @temp) . "</th></tr>";
	my $table2 = "<p><br /></p><table  border frame=box><caption>Pathway Detail</caption><tr><th>#</th><th>" . substr($pre, 0);# . "</th><th>Pathway ID</th><th>$key_id</th><th>$gene_list_title</th><th>$key_id desc</th></tr>";
    my @cols = ("","Pathway ID","$key_id","$gene_list_title","$key_id desc");
#    if($run_symbol ne 'none'){
#        @cols = ("","Pathway ID","$key_id","$gene_list_title","Symbol","$key_id desc");
#    }
    $table2 .= join("</th><th>",@cols)."</tr>";
    my $index = 0;
	my $classnum = 1;
    my $nrow_index_for_color = 1;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my $ldesc;my $lgene;
		$index++;
		@temp = split /\t/, $_;
		shift @temp; shift @temp;
		$pre = shift @temp;
		my $lk_id = pop @temp;
		my @k_ids = split/\Q+\E/,$lk_id;
		my $num = @k_ids;
		$num+=1 if($num > 1);
#		$lk_id =~ s/\Q+\E/<br\/>/g;
		my $gene = pop @temp;
		if(scalar(@temp) > 2){
			$temp[-2] = sprintf("%.6f", $temp[-2]);
			$temp[-3] = sprintf("%.6f", $temp[-3]);
		}
		my $nnum = $index%2+1;
#        my $nrow_index_for_color = $index*2-1;
        
		$code .= "<tr class=\"$nnum\"><td>$index</td><td style=\"text-align: left;\"><a href='#gene$index' title='click to view $key_name' onclick='javascript: colorRow(document.getElementsByTagName(\"table\")[1],$nrow_index_for_color,$num);'>$pre</a></td><td>" . (join "</td><td>", @temp) . "</td></tr>";
        $nrow_index_for_color = $nrow_index_for_color + $num;
        my $map = $temp[-1];
        $map =~ s/ko/map/;
        $table2 .= "<tr class=\"$nnum\">\n";
        $table2 .= "<td rowspan=\"$num\">$index</td><td rowspan=\"$num\" style=\"text-align: left;\">";
#        print "$dir/$name\_map/$map.html";die;
        if (-e "$dir/$name\_map/$map.html") {
            my $ko = $map;
            $ko =~ s/map/ko/;
            $table2 .= "<a href='$name\_map/$map.html?group=$name&id=$ko' title='click to view map' target='_blank'>$pre</a>";
        }else {
#            print STDERR "$indir/$name\_map/$map.html\n";
            my $err = ("#"x80)."\n\nError: no map [$dir/$name\_map/$map.html] in kegg database for [$files[$i]]\n\n".("#"x80)."\n";
            if($nomap_die == 1){
                die "$err";
            }else{
                print STDERR "Error: no map [$dir/$name\_map/$map.html] in kegg database\n";
            }
#			$table2 .= "$pre (no map in kegg database)";
		}
		$table2 .= "</td>\n<td rowspan=\"$num\" style=\"text-align: left;\"><a name='gene$index'></a>$temp[-1]</td>\n"; ## col 2 koid
		$table2 .= GetInfo($lk_id,$gene,\%ko_desc,\$classnum)."</tr>\n";
	}
	$table2 .= "</table>";
	$code .= "</table>\n";
    $code .= "<p>ps: ko011开头的通路为global map，ko012开头的通路为overview map，仅展示kegg例图</p>\n" if($version >= 2);
	$code .= "$table2\n<div id=\"backtop\"><a href=\"#\">Back Top</a></div><script type='text/javascript'>showPer(document.getElementsByTagName('table')[0]);\ndiffColor([document.getElementsByTagName('table')[0], document.getElementsByTagName('table')[1]]);markColor(document.getElementsByTagName('table')[0]);\n";
	$code .=<<HH;
	var tr=document.getElementsByClassName('2');
	for(var i=0;i<tr.length;i++){
		tr[i].style.background='#ccff99';
	}
	var tr=document.getElementsByClassName('1');
	for(var i=0;i<tr.length;i++){
		tr[i].style.background='#ffffff';
	}
</script>
</body>
</html>
HH
	close IN;

	open HTML, "> $htmlFile" || die $!;
	print HTML "$code";
	close HTML;
}

exit 0;



sub GetInfo{
    my $ko = shift;
    my $gene = shift;
    my $ko_desc = shift;
    my $class = shift;

    my @kos = split/\Q+\E/,$ko;
    my @genes = split/;/,$gene;
    if($gene =~ /\);/){
        @genes = split/\);/,$gene;
        @genes = map{"$_)"}@genes;
    }

    my %hash;
    for my $i(0..$#kos){    
        $hash{$i} = $kos[$i];
    }
    
    my $html = '';
    my @keys = (0);
    if(@kos > 1){
        @keys = sort {$hash{$a} cmp $hash{$b} or $a <=> $b} keys %hash;
    }
    for my $y(@keys){
        $$class = GetClass($$class);
        $$ko_desc{$kos[$y]} //= '';
#        if($genes[$y] =~s/\((.*)//){
        my @infos = ($kos[$y],$genes[$y],$$ko_desc{$kos[$y]});
#        if($run_symbol ne 'none'){
#            @infos = ($kos[$y],$genes[$y],$symbol,$$ko_desc{$kos[$y]});
#        }
        $html .= "\t<tr class=\"$$class\">\n" if(@kos > 1);
        for my $info(@infos){
            $html .= "\t\t<td style=\"text-align: left;\">$info</td>\n";
        }
        $html .= "\t</tr>\n" if(@kos > 1);
    }

    return $html;
}


sub GetKoDesc{
	my $file = shift;
	my $hash = shift;

	open FILE,$file or die $!;
	while(<FILE>){
		chomp;
		my @aa = split/\t/,$_;
		$aa[0] =~ s/^ko://;
		$aa[0] =~ s/^cpd://;
		$aa[1] =~ s/.*;// if(/^K/);
#		$aa[1] =~ s/\[.*//;
		$$hash{$aa[0]} = $aa[1];
	}
	return 0;
}

sub GetClass{
	my $num = shift;

	if($num == 1){
		return 2;
	}else{
		return 1;
	}
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    system($cmd);

    return 0;
}
