#! /usr/bin/perl

use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Basename qw/basename dirname/;
use FindBin qw/$Bin $Script/;
use Cwd qw/abs_path/;

#use lib "/Bio/Bin/software/report/lib/"; ## old version
#use lib "/Bio/User/aipeng/bin/report/lib"; ## new version, beta
use lib "/Bio/Bin/software/GDHR/current/lib/"; ## new version
use GDHR;

use Data::Dumper;
use YAML;

#use lib "/Bio/Bin/pipeline/GDreport/v2.0/lib/";
#use enrich;

use lib '/Bio/Bin/pipeline/Toolkits/pipeline_lib/v2.2';
use DIR;
use Refer;
use HelpHtml;


#### Options ####
my ( $old_style, $verbose ) = ( 0, 0 );
GetOptions( "old|o" => \$old_style, "verbose|v" => \$verbose );
die "perl $0 (-old|o) <conf> <upload_dir>\n" unless @ARGV == 2;
my ( $conf_file, $outdir ) = @ARGV;

$outdir //= '.';
$outdir   = abs_path( $outdir );


#### Initals ####
my $conf = &read_conf( $conf_file );
die Dumper $conf if $verbose;

## 文件夹系统
$DIR::is_mkdir = 0;
my $dirs = DIR->new(path => $outdir, base => $outdir); ## old version is "$outdir/src";

## 参考
my $refer = Refer->new();

## 图片注释
my $img_help = HelpHtml->new(title => "图片注释", file => "$outdir/src/doc/img_help.html");


#### Main ####
my $report = GDHR->new( -outdir => $outdir, -pipe => "scRNA-seq组间差异", -name => "index", -nonlazy => $old_style );

&Project_info;
&Introduce;
&Diff_stat;
&Diff_pict;
&Diff_enrich;
&Diff_GSEA;
&Diff_GSVA;
&Catalog;
&Reference;
&Appendix;

$report->write;

system "cp $Bin/../help/*png $outdir/src/image/";
system "cp $Bin/../help/help.html $outdir/src/doc";
#system "cp $Bin/../help/enrichment.html $outdir/src/doc";
system "cp $Bin/../help/scRNAseq_SampleDiffer.method.pdf $outdir/src/doc";
system "cp $Bin/../help/单细胞样本差异分析.pdf $outdir/src/doc";

$img_help->dump();		
#### END ####




################################################
################################################
sub read_conf {
		my ( $file ) = @_;

		my $config = YAML::LoadFile($file);

		return $config;
}

sub Project_info {
		my $infomation = $report->section( id => 'Project_info' );
		$infomation->menu( "项目概述" );

		my $separator = '<span class="project_info_separator">¦</span>';
		my $data = [
				[ '项目编号',   $conf->{project_id} ],
				[ '项目内容',   $conf->{content} ],
				[ '参考基因组', $conf->{ref} ],
				[ '样品名称',   join ( $separator, @{$conf->{sample}} ) ],
		];
		if ( exists $conf->{group} ){
				push @$data, [ "分组方案", join ( $separator, @{$conf->{group}} ) ];
		}
		if ( exists $conf->{group_gene} ){
				push @$data, [ "分组方案（按基因）", join ( $separator, @{$conf->{group_gene}} ) ];
		}
		if ( exists $conf->{Sde} ) {
				push @$data, [ "差异方案", join ( $separator, @{$conf->{Sde}} ) ];
		}

		$infomation->matrix2html( -matrix => $data, -header => 0, -class => "hl_table", -no_order => 1 );
}

sub Introduce {
		my $Introduce = $report->section( id => 'Introduce_' );
		$Introduce->menu( "项目介绍", -help =>'help.html#Introduce' );

		## 2.1 ##
		$Introduce->submenu( "背景介绍" , -help =>'help.html#Introduce_exp');
		$Introduce->add_html("<p style='text-indent:2em'>样本分组信息是在实验设计时就已经确立的细胞间关系，但是在针对细胞亚群的基础分析中并没有得到充分的利用。所以，哪些基因在表型变化中起着关键作用无法得到确认。</p>");
		$Introduce->add_html("<p style='text-indent:2em'>单细胞转录组样本间差异分析针对同时分布于两个及两个以上样本的细胞亚群进行，深入探索同一细胞亚群/细胞类型在不同样本之间的表达差异，为完成转录本变化和表型关联提供新的解析角度。</p>");

		## 2.2 ##
		$Introduce->submenu( "分析流程", -help => 'help.html#Introduce_ana');
		$Introduce->img2html( -file => 'image/h1.pipeline.png', -name => "分析流程示意图", -width=>"80%");
		my $ref = $refer->add( "seurat", "Stuart T, Butler A, Hoffman P and et al. Comprehensive integration of single cell data. [J] Cell, 2019. 177(7):1888-1902.e21." );
		$Introduce->add_html("<p style='text-indent:2em'>利用Seurat$ref对细胞重新标注和差异基因分析，完成对差异基因的统计和分布绘图。然后，使用GO数据库和KEGG数据库对差异基因进行富集分析。</p>");


		$Introduce->break;
}

sub Diff_stat {
		my $Diff_stat = $report->section( id => '_differ_stat' );
		$Diff_stat->menu( "差异基因整体统计" ); #, -help => "help.html#differ_stat" );

		$dirs->child( 'result', is_order => 1 );		
		my $fulldir = $dirs->result->{_path};
		my $reldir = $dirs->result->rel();
#		my $dirname = basename $fulldir;

		my $ref = $refer->add("MAST", "Finak G, McDavid A, Yajima M and et al. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. [J] Genome Biol, 2015. 16:278." );
		$Diff_stat->add_html( "<p style='text-indent:2em'>基于细胞所属的亚群信息和样本信息，使用Seurat软件对细胞进行进行组间差异分析，设置|log2FC|≥$conf->{LOGFC}、任意一组中表达目标基因的细胞比例≥$conf->{PCT}为阈值，使用MAST$ref障碍模型检验差异显著性。通过Seurat的BH法对差异基因的显著性p值进行多重检验矫正，筛选@{[$conf->{'use.q'} =~ /true/i ? '矫正后' : '']}p值≤$conf->{PVALUE}的基因做为显著差异基因，得到显著差异基因结果表。显著差异基因将用于后续分布情况和富集分析。</p>" );
		$Diff_stat->add_html( "<p style='text-indent:2em'>根据每个细胞亚群得到的上调基因数和下调基因数进行统计，获得差异基因统计表并依此绘制差异基因统计图。</p>" );


		my @pngs = map { "$reldir/$_/Stat.$_.png" } @{$conf->{differ_name}};
		$Diff_stat->imgs2html( -files => \@pngs, -names => $conf->{differ_name}, -name => "差异基因统计图@{[help_logo('差异基因统计图')]}");
		$img_help->add("差异基因统计图", "横坐标表示不同的细胞亚群，每个细胞亚群对应一对条形柱，红色代表上调基因，蓝色代表下调基因；纵坐标为基因数目。");

		my @files = map { "$fulldir/$_/Stat.$_.xls" } @{$conf->{differ_name}};
		$Diff_stat->tsvs2html( -files => \@files, -names => $conf->{differ_name}, -name => "差异基因统计表", -help => "help.html#差异基因统计表", -top => 200, -class => "hl_table", -orient => "vertical");

		$Diff_stat->desc("<strong>所有差异基因结果表</strong>");
		my @lists = map { "$reldir/$_/DifferMarker.$_.all.annot.xls" } @{$conf->{differ_name}};
		$Diff_stat->files2list( -files => \@lists, -desc => $conf->{differ_name} );

		$Diff_stat->add_html( "<p style='text-indent:2em'>因为Seurat会先以阈值过滤基因，然后再做显著性检验；所以，对于不满足阈值的基因，我们赋值其矫正后p值=1。最后，根据基因在比较组间的差异倍数和矫正后p值绘制火山图，直观显示基因的差异显著性。</p>" );

#		my $vol_list = YAML::LoadFile("$fulldir/Volcano.list");
		## YAML will turn every charcter into utf-8 format, such as \x{3B3}.
		## It will cause unreadable page in GDHR. So read the file by hand here.
		my $vol_list = {};
		{ 
				if ( -s "$fulldir/Volcano.list" ) {
				my $last = "";
				open my $in, "<", "$fulldir/Volcano.list";
				while (<$in>) {
						chomp;
						if (/:/) {
								if ( /:\s*\S+/ ) {
										my @data = map { s/^['"]|['"]$//gr } split /\s*[:,]\s*/;
										$last = shift @data;
										$vol_list->{$last} = \@data;
								} else {
										$last = $_ =~ s/://r;
								}
						} else {
								s/^\s*-\s*//;
								s/^['"]|['"]$//g;
								push @{$vol_list->{$last}}, $_;
						}
				} 
				close $in;
				} else {
						$vol_list->{order} = [];
				}
		}
#		print Dumper $vol_list;
		$vol_list->{order} = [ $vol_list->{order} ] if ref $vol_list->{order} ne 'ARRAY';
		my @divs = map {
				my $i = $_;
				my @img_div = map { $Diff_stat->_img2html( -file => $_ ) } map { "$reldir/$i/" . ("Volcano.${i}_${_}.png" =~ s#/#_#gr) } @{$vol_list->{$i}};
				my @names = @{$vol_list->{$i}};
				$Diff_stat->any2child_vtab( $Diff_stat->div_pack(@img_div), -names => \@names, -name => '亚群中各组的富集分数热图', -type => 'image' );
		} @{$vol_list->{order}};
		$Diff_stat->any2hontab( $Diff_stat->div_pack(@divs), -names => $vol_list->{order}, -name => "差异基因火山图@{[help_logo('差异基因火山图')]}", -type => 'image' );


		$img_help->add('差异基因火山图', "横坐标表示两个分组间的差异倍数对数值，纵坐标表示两个分组差异的 -log10(FDR或p值)，红色（ group_2 相对于 group_1 表达量上调）和蓝色（表达量下调）的点表示基因的表达量有差异，灰色的点表示没有差异的基因");

		$Diff_stat->tsv2html( -file => "$fulldir/$_/DifferMarker.$_.diff.annot.xls",
				-name => "$_显著差异基因结果表（前20行）",
				-help => "help.html#差异基因结果表",
				-top  => 20) for @{[ (grep { -s $_ } @{$conf->{differ_name}})[0] ]}; ## careful! [for] will change $_.

		$Diff_stat->desc("<strong>显著差异基因结果表</strong>");
		@lists = map { "$reldir/$_/DifferMarker.$_.diff.annot.xls" } @{$conf->{differ_name}};
		$Diff_stat->files2list( -files => \@lists, -desc => $conf->{differ_name} );

		$Diff_stat->break;
}

sub Diff_pict {
		my $Diff_pict = $report->section( id => '_differ_pict' );
		$Diff_pict->menu( "差异基因在细胞亚群分布情况" );

#		$dirs->child( 'result', is_order => 1 );
		my $fulldir = $dirs->result->{_path};
		my $reldir = $dirs->result->rel();
#		my $dirname = basename $fulldir;


		$Diff_pict->add_html( "<p style='text-indent:2em'>为了更直观地展示基因在不同样本及细胞亚群的分布特征，我们将不同基因依据表达量变化情况，通过R包绘制成热图、UMAP映射图和气泡图(这里只展示每个cluster差异倍数前5的基因)。</p>" );

		$Diff_pict->submenu( "差异基因热图分布" );
		$Diff_pict->add_html( "<p style='text-indent:2em'>我们使用z-score法对基因表达量进行归一化处理，然后将归一化后的基因表达量用于热图绘制。如此，避免热图中极端值的出现，使得基因在不同细胞群的表达特征更加明显。</p>" );
		my @pngs = map { "$reldir/$_/Heatmap.$_.png" } @{$conf->{differ_name}};
		$Diff_pict->imgs2html( -files => \@pngs, -names => $conf->{differ_name}, -name => "差异基因表达热图@{['差异基因表达热图']}");
		$img_help->add('差异基因表达热图', "热图每一行代表一个基因，每一列代表一个细胞。图中同一个分组的细胞被归类放在一起，在热图上方以contrast来标注每个细胞所属分组；每个分组内，同一个细胞亚群的细胞被归类放在一起，在热图上方以cluster来标注每个细胞所属细胞亚群。</p><p>绘图主要使用的R包：Seurat、pheatmap");

		$Diff_pict->submenu( "差异基因表达分布" );
		$Diff_pict->add_html( "<p style='text-indent:2em'>我们使用log均一化法对基因表达量进行处理，减小基因在不同细胞间的表达差异性；然后将均一化后的基因表达量用于UMAP映射图的绘制。如此，在保留一定原始表达量特征的情况下减小了极值对数据分布的影响，便于观察基因在不同细胞间的表达特征。</p>" );
		@pngs = map { "$reldir/$_/FeaturePlot.$_.umap.demo.png" } @{$conf->{differ_name}};
		$Diff_pict->imgs2html( -files => \@pngs, -names => $conf->{differ_name}, -name => "差异基因表达分布图@{[help_logo('差异基因表达分布图')]}");
		$img_help->add('差异基因表达分布图', "该图由若干张UMAP映射图组成，每一张UMAP映射图代表一个基因，分属于两个分组的细胞在图中分割为两个子图进行展示。</p><p>绘图主要使用的R包：Seurat、ggplot2" );

		$Diff_pict->submenu( "差异基因气泡图分布" );
		$Diff_pict->add_html( "<p style='text-indent:2em'>气泡图同时展示了基因在细胞亚群中的平均表达量和表达目标基因细胞的比例。为了绘制气泡图，我们对每个基因在每个细胞亚群的未取log值的均一化表达量取均值，然后基于表达量均值在不同亚群间对基因平均表达量做z-score归一化，使用归一化的基因表达量绘制气泡图中气泡的颜色深浅。另外，我们统计表达目标基因的细胞在亚群中占比，用于绘制气泡图中气泡的大小。</p>" );
		@pngs = map { "$reldir/$_/DotPlot.$_.png" } @{$conf->{differ_name}};
		$Diff_pict->imgs2html( -files => \@pngs, -names => $conf->{differ_name}, -name => "差异基因表达气泡图@{[help_logo('差异基因表达气泡图')]}" );
		$img_help->add('差异基因表达气泡图', "气泡图的横坐标为基因名称，纵坐标为细胞的亚群信息和分组信息。气泡大小代表细胞亚群表达横坐标对应基因的细胞比例；颜色代表不同分组信息，对于比较组group_1 vs group_2，红色为group_1，蓝色为group_2；颜色的深浅代表基因在细胞亚群的平均表达量</p><p>绘图主要使用的R包：Seurat、ggplot2");

		$refer->add("PanY", "Pan Y, Lu FC, Fei QL and et al. Single-cell RNA sequencing reveals compartmental remodeling of tumorinfiltrating immune cells induced by anti-CD47 targeting in pancreatic cancer. [J] Journal of Hematology & Oncology, 2019. 12(1):124.");
		$refer->add("Mathys", "Mathys H, Davila-Velderrain J, Peng Z and et al. Single-cell Transcriptomic Analysis of Alzheimer's Disease. [J] Nature, 2019. 570(7761):332-337.");
		$refer->add("Subramanian", "Subramanian A, Sidhom EH, Emani M and et al. Single Cell Census of Human Kidney Organoids Shows Reproducibility and Diminished Off-Target Cells After Transplantation. [J] Nature Communications, 2019. 10(1):5462.");

		$Diff_pict->break;
}

sub Diff_enrich {
		my $Diff_enrich = $report->section( id => '_differ_enirch' );
		$Diff_enrich->menu( "差异基因富集分析" );

		$dirs->child( 'Enrichment', is_order => 1 );
		my $fulldir = $dirs->Enrichment->{_path};
		my $reldir = $dirs->Enrichment->rel();
#		my $dirname = basename $fulldir;
		
		my %list;
		my $array = [ sort { $list{$a} <=> $list{$b} } keys %list ];
		my $run_do = -d "$fulldir/Enrichment/DO" ? 1 : 0;
		my $run_reactome = -d "$fulldir/Enrichment/Reactome" ? 1 : 0;

		use lib $Bin;
		use enrich;

#		enrich::enrich( $Diff_enrich, "$fulldir/Healthy-vs-Asymptomatic", "$reldir/Healthy-vs-Asymptomatic", "$fulldir/Healthy-vs-Asymptomatic/enrich.list", do => 1, reactome => 1, help_dir => "$outdir/src/doc", refer => $refer);
		enrich::multi_enrich( $Diff_enrich, $fulldir, $reldir, "$fulldir/enrich.list", do => 1, reactome => 1, help_dir => "$outdir/src/doc", refer => $refer);

		$Diff_enrich->break;
}

sub Diff_GSEA {
		my $Diff_gsea = $report->section( id => '_differ_gsea' );
		$Diff_gsea->menu( "GSEA分析", -help => "gsea_help.html#GSEA" );

		$dirs->child( 'GSEA', is_order => 1 );
		my $fulldir = $dirs->GSEA->{_path};
		my $reldir = $dirs->GSEA->rel();
		
		use lib $Bin;
		use gsea;

		gsea::multi_gsea($Diff_gsea, $fulldir, $reldir, "$fulldir/enrich.list", -refer => $refer, -help => 'gsea_help.html');

		my $help = HelpHtml->new(title => "GSEA注释", file => "$outdir/src/doc/gsea_help.html");
		$help->add_html(gsea::gsea_help());
		$help->dump;

		$Diff_gsea->break;
}

sub Diff_GSVA { 
		my $GSVA = $report->section( id => 'GSVA_' );
		$GSVA->menu( "GSVA分析", -help => "gsva_help.html#GSVA" );

		$dirs->child( 'GSVA', is_order => 1 );
		my $fulldir = $dirs->GSVA->{_path};
		my $reldir = $dirs->GSVA->rel();

		use lib $Bin;
		use gsva;

		gsva::multi_gsva($GSVA, $fulldir, $reldir, "$fulldir/enrich.list", -refer => $refer);

		my $help = HelpHtml->new(title => "GSVA注释", file => "$outdir/src/doc/gsva_help.html");
		$help->add_html(gsva::gsva_help($refer));
		$help->dump;

		$GSVA->break;
}


sub Catalog {
		my $catalog = $report->section( id => 'Catalog' );
		$catalog->menu( "目录结构" );
		open CATALOG, "<", "$Bin/../help/catalog.html" or die $!;
		my $catalog_code = join "", <CATALOG>;
		close CATALOG;
		$catalog->add_html( $catalog_code );
		$catalog->break;
}

sub Reference {
		my $Reference = $report->section( id => "Reference" );
		$Reference->menu( '参考文献' );


		$Reference->add_html( "<p><ul>" );
		$Reference->add_html( "<li>$_</li>" ) for $refer->getAllRefer;
		$Reference->add_html( "</ul></p>" );
		$Reference->break;
		
		$Reference->add_html( '<script type="text/javascript">
				function close_help(){
						var f = document.getElementById("show_help");
						var b = document.getElementById("bgbox");
						f.style.display = "none";
						document.body.removeChild(b);
				}</script>' ); ## in order to close help page then jump to referrence, when click cite index in help page
}

sub Appendix {
		my $appendix = $report->section( id => 'Appendix' );
		$appendix->menu( "附录" );
#		$appendix->submenu( "分析方法英文文档" );
#		$appendix->add_html( "<div style='text-indent:2em'><p>scRNA-seq分析方法文档（英文）：<a href='../scRNA-seq_method.pdf' target='_blank'>scRNA-seq_method.pdf</a></p></div>" );
		$appendix->submenu( "结果文件查看" );
		$appendix->add_html( "<div style='text-indent:2em'>
						<p>*.xls,*.txt ：结果数据表格文件，文件以制表符（Tab）分隔。unix/Linux/Mac用户使用 less 或 more 命令查看；windows用户使用高级文本编辑器Notepad++ 等查看，也可以用Microsoft Excel打开。</p>
						<p>*.png：结果图像文件，位图，无损压缩。</p>
						<p>*.pdf：结果图像文件，矢量图，可以放大和缩小而不失真，方便用户查看和编辑处理，可使用Adobe Illustrator进行图片编辑，用于文章发表等。</p></div>" );
		$appendix->submenu( "文章引用与致谢" );
		$appendix->add_html( "<div style='text-indent:2em'><p>如果您的研究课题使用了基迪奥的测序和分析服务，我们期望您在论文发表时，在Method部分或Acknowledgements部分引用或提及基迪奥公司。以下语句可供参考：</p></div><ul>
						<li><strong>Method部分</strong>：The cDNA/DNA/Small RNA libraries were sequenced on the Illumina sequencing platform by Genedenovo Biotechnology Co., Ltd (Guangzhou, China).</li>
						<li><strong>Acknowledgements部分</strong>：We are grateful to/thank Guangzhou Genedenovo Biotechnology Co., Ltd for assisting in sequencing and/or bioinformatics analysis.</li></ul>" );

		$appendix->break;
}



sub help_logo {
		my ( $tag, $html ) = @_;
		$html //= "src/doc/img_help.html";

		return qq(<a href="$html#$tag" target="help_page" onclick="show_help();"><img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABoAAAAaCAQAAAADQ4RFAAAABGdBTUEAALGPC/xhBQAAACBjSFJN
AAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAAmJLR0QA/4ePzL8AAAAJcEhZ
cwAACxMAAAsTAQCanBgAAAAHdElNRQfkCQ4SAwcGJVcIAAABEklEQVQ4y52U0XXCMAxFrzn8hw3M
BmGDsAHdIN2AERjBIzACbEA3yAhJJ2AE9SOiYCHRHuwf2/K1ZD/JSXhuKbPW4STfzgapOpnCWC2N
FLLZ9TBsKOaMey80DkTLECKCMNAaiJbrS0QQrjfsFtjwJzJ7a+6QvcuRHZnMjpO9m0JkY+irl+qN
Nc9Q7eeimzs6HV2sLwSjy14QOCLIjHGodROWD+rPbZtWwAdRW6cMXfhWew1vNOtdDPWh5CE0aCZ6
kveLIPIvAD5ZObZpyeRC23QIn2PyLvq6j8ICODunbSRJYuNYzrhpJIhTbnUaOQkbQ+V31SmNVlWK
S+OtIgy1r73YctcgS4j4H4tO//WFpXc+yx9W6OdFyDnIhAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAy
MC0wOS0xNFQxODowMzowNyswODowMK1g45IAAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjAtMDktMTRU
MTg6MDM6MDcrMDg6MDDcPVsuAAAATnRFWHRzb2Z0d2FyZQBJbWFnZU1hZ2ljayA2LjguOC0xMCBR
MTYgeDg2XzY0IDIwMTUtMDctMTkgaHR0cDovL3d3dy5pbWFnZW1hZ2ljay5vcmcFDJw1AAAAGHRF
WHRUaHVtYjo6RG9jdW1lbnQ6OlBhZ2VzADGn/7svAAAAF3RFWHRUaHVtYjo6SW1hZ2U6OkhlaWdo
dAAyNjaCDawAAAAWdEVYdFRodW1iOjpJbWFnZTo6V2lkdGgAMjbOLc0hAAAAGXRFWHRUaHVtYjo6
TWltZXR5cGUAaW1hZ2UvcG5nP7JWTgAAABd0RVh0VGh1bWI6Ok1UaW1lADEzNDAwMDAwMTDVyQWX
AAAAEXRFWHRUaHVtYjo6U2l6ZQA0MDhCQkU+CXwAAABadEVYdFRodW1iOjpVUkkAZmlsZTovLy9o
b21lL3d3d3Jvb3Qvd3d3LmVhc3lpY29uLm5ldC9jZG4taW1nLmVhc3lpY29uLmNuL3NyYy8xMDcz
Ni8xMDczNjY5LnBuZxVTqewAAAAASUVORK5CYII=
" class="help_logo"></a>);
}



