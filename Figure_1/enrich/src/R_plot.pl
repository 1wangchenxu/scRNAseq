#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: plot_pie.pl
#
#        USAGE: ./plot_pie.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: smhong (), smhong@genedenovo.com
# ORGANIZATION: R&D
#      VERSION: 1.0
#      CREATED: 01/09/2019 04:33:40 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use FindBin qw/$RealBin $Script/;
#use File::Basename qw/basename dirname/;
use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;
use Cwd 'abs_path';
use Carp qw(carp confess croak);
use lib "$RealBin/../";
use Support_Program;
#use lib "/Bio/Bin/pipeline/GeneralPlot/v1.0/lib/";
use lib "$RealBin/../../../System_Programs//GeneralPlot/v1.0/lib";
#use lib "/Bio/User/liuyubin/pipeline/GeneralPlot/lib/";
use GeneralPlot;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> [draw_type] -infile -outdir -outpfx -title -col -header
            percent_pie
            stauration
            pq_bar
            vip_bar
            diff_zscore
            dotplot
            pca_point
            Gradient
            enrichmentHeatmap
            pheatmap
            fpkm_distribution
            violin
            bar
            fill_bar
            break_bar
            volcano
            two_pie
            FaLengthDistribution
            Trend_rect_profile
            Trend_profile_fig
            Trend_profile_line
            Trend_profile_table
            GO_level2
            venn
            kegg_PathwayClass
            ReplotGSEA
            cog_kog
            box
            diff_scatter_plot
            histogram
            histogram_facet
            density
            ternary_plot
            busco_figure
            enrichPlot
            enrich_bubble
            random
            ggwordcloud
            ggseqlogo
            ggraph
            ggpubr
            letter_sig

            PCA
            cor
            ggcor
            corrplot
            enrich
            topGO_graph
            plsda_loading
            plsda_permutation
            plsda
            mageck_RRA
            plotPearson
            cluster
            roc
            MSEAp
            MSEAqea
            MSEAora
            ipath
            xGSEAdotplot
            UMAP

            SBV

            fdr
            get_part_col

            example

    e.g: perl $0 <> <>

USAGE

die "$usage" if(@ARGV<1); 

my $src = "$RealBin/src";
my $perl = "$RealBin/perl";
my $python = "$RealBin/python";
my $apple = "$RealBin/apple_lib";
my $sbv = "$RealBin/SBV/bin/sbv.pl";
my $soft_conf = "$RealBin/../Small_pipe.soft.conf";
#my $inkscape = "$RealBin/inkscape";

my %opts = ("Rscript_version"=>"Rscript.3.6.0");
#my $draw_type = shift;
$opts{run} = 1;
$opts{density} = 300;
$opts{new} = 0;
$opts{default} = 1;
GetOptions(\%opts,"i=s","infile=s","out=s","outpfx=s","outdir=s","title=s","col=s","row=s","header=s",
                "exclude_row=s","outfile=s","line=s",
                "title=s","xlab=s","ylab=s","dot=s","ymin=s","ymax=s","color=s","format=s","indir=s",
                "xmin=s","xmax=s",
                "add_ellipse=s",
                "min=s", "max=s", "by=s","name=s","name_list=s",
                "xcol=s","ycol=s","xcut=s","ycut=s",
                "melt=s","rowname=s","slope=s",
                "window=s","all_reads=s","sample=s","type=s","name_key=s","top=s","group=s","rmcol=s",
                "trend_genetable=s","profile_name=s","trend_profiletable=s", # Trend
                "conf=s","permutation=s",
                "cluster_col=s","filter_na=s","display_numbers=s",# pheatmap
                "label_info=s","show_cluster_tree=s","column_dend_height=s","row_dend_width=s","cluster_row=s", "group_color=s",## complex heatmap
                "cell_width=s","cell_height=s", "border_color=s", "border_size=s",# heatmap
                "list=s","glist=s",# venn "header=s","o=s"
                "target=s","cor_cut=s","pvalue_cut=s","name1=s","name2=s","col_or_row=s",  # cor 
                "fontsize=s","method=s","display=s","border=s","showrow=s","showcol=s", #cor
                "logformat=s", ## log2 yes or not 
                "cog_kog=s", ## annot 
                "xls=s","log=s",
                "ref=s",
                "fg_all=s","bg_all=s","diff=s","pfc=s",#enrich 
                "pq=s","pq_cut=s","fc_cut=s","pq_col=s","fc_col=s",## diff scatter plot 
                "updown=s","scale=s",## venn
                "break_pos=s","ratio=s", # break_barplot
                "break_info=s",
                "name_col=s","vip_col=s","vip_cut=s", # loading plot
                "neighbors=s","metric=s","min_dist=s","legend=s","do_label=s","label_size=s","byrow=s",#UMAP
                "ellipticity=s","seed=s",#ggwordcloud2
                "first_pos=s","range=s","units=s","color_info=s","ignore_case=s",##ggseqlogo
                "dot_color=s","dot_alpha=s","error_bar=s",#ggpubr
                "title_size=s","lab_size=s","axis_size=s","text_size=s",## letter_sig
                "back_color=s","order=s","linetype=s",
                "size=s","alpha=s","plot_single=s",
                "filter0=s","show_text=s",
                "lipid=s","class=s",
                "stat=s","soft=s",
                "ncol=s","pvalue=s","cor=s","reverse=s","default=s",
                "new=s","yaml=s","run=s","density=s","Rscript=s","Rscript_version=s");
my $draw_type = shift;

## software 
my %apps;
#Support_Program("$RealBin/../",\apps,1);
Support_Program_conf("$RealBin/../",\%apps,"$soft_conf");
Error("No Rscript version [$opts{Rscript_version}]") if(!exists $apps{$opts{Rscript_version}});
my $Rscript = "$apps{$opts{Rscript_version}}";
print "$Rscript\n";

mkdir $opts{outdir} if(exists $opts{outdir} && !-e $opts{outdir});
my @g_clean_files;
my %draw_pipe = (
    "percent_pie" => \&PercentPie,
    "stauration" => \&Stauration,
    "pq_bar" => \&PQBar,
    vip_bar => \&VIPBar,
    diff_zscore=>\&DiffZscore,
    "Gradient" => \&Gradient,
    "enrichmentHeatmap" => \&EnrichHeatmap,
    "fill_bar" => \&FillBar,
    "pheatmap" => \&Pheatmap,
    "fpkm_distribution" => \&FpkmDistribution,
    "violin" => \&Violin,
    "box" => \&Box,
    "bar" => \&Bar,
    dotplot => \&DotPlot,
    pca_point => \&PcaPoint,
    "break_bar" => \&BreakBar,
    "volcano" => \&Volcano,
    "two_pie" => \&TwoPie,
    "FaLengthDistribution" => \&FaLengthDistribution,
    "Trend_rect_profile" => \&TrendRectProfile,
    "Trend_profile_fig" => \&TrendProfileFig,
    "Trend_profile_line" => \&TrendProfileLine,
    "Trend_profile_table" => \&TrendProfileTable,
    "GO_level2" => \&GOLevel2,
    "venn" => \&Venn,
    "cor" => \&Cor,
    ggcor => \&GGCor,
    corrplot => \&Corrplot,
    kegg_PathwayClass => \&KeggPathwayClass,
    ReplotGSEA => \&ReplotGSEA,
    "PCA" => \&PCA,
    "cog_kog" => \&COG,
    density => \&Density,
    enrich => \&Enrich,
    topGO_graph => \&topGOGraph,
    busco_figure => \&BuscoFigure,
    diff_scatter_plot => \&DiffScatterPlot,
    histogram => \&Histogram,
    histogram_facet => \&HistogramFacet,
    ternary_plot => \&TernaryPlot,
    SBV => \&SBV,
    enrichPlot => \&EnrichPlot,
    enrich_bubble => \&EnrichBubble,
    random => \&Random,
    fdr => \&Fdr,
    get_part_col => \&GetPartCol,
    plsda_loading => \&PlsdaLoading,
    plsda_permutation => \&PlsdaPermutation,
    plsda => \&Plsda,
    mageck_RRA => \&MageckRRA,
    plotPearson => \&PlotPearson,
    cluster => \&Cluster,
    roc => \&ROC,
    MSEAp => \&MSEAp,
    MSEAqea => \&MSEAqea,
    MSEAora => \&MSEAora,
    ipath => \&iPath,
    xGSEAdotplot => \&xGSEAdotplot,
    UMAP => \&UMAP,
    ggwordcloud => \&GGwordcloud,
    ggseqlogo   => \&GGseqlogo,
    ggraph      => \&GGraph,
    ggpubr      => \&GGpubr,
    letter_sig  => \&LetterSig,

    example => \&Example,
);

Error("Unknown type [$draw_type]") if(!exists $draw_pipe{$draw_type});
GetColors($draw_type,\%opts);
my $outpdf = &{$draw_pipe{$draw_type}}(%opts);
Converts($outpdf);

Clean(@g_clean_files) if(@g_clean_files > 0);

sub GetColors{
    my $type = shift;
    my $opts_ref = shift;
   
    return if(!exists $opts_ref->{infile} && !exists $opts_ref->{list}); ## venn -list 
    return if(exists $opts_ref->{color});
    my $color_arr;
    my $color_num = 0;
    my %info = qw/tag default row_col row rm_col_num 0/;
    $info{file} = $opts_ref->{infile} if(exists $opts_ref->{infile});
    if(exists $opts_ref->{rmcol} && $opts_ref->{rmcol} ne '0'){
        $info{rm_col_num} = scalar(split/,/,$opts_ref->{rmcol});
    }
    if($type eq 'fill_bar2'){
#        %infos = qw/type bar tag stack2 row_col col max 10/;
#        $infos{row_col} = "row" if(exists $opts_ref->{row} && $opts_ref->{row} eq 'yes');
        my ($row,$col) = ReadFile($opts_ref->{infile});
        if(exists $opts_ref->{row} && $opts_ref->{row} eq 'yes'){ ## by rownum 
            $color_num = $row - 1;
        }else{
            $color_num = $col - 1;  ## by col_num
        }
        if(exists $opts_ref->{rmcol} && $opts_ref->{rmcol} ne '0'){
            my @rm_col = split/,/,$opts_ref->{rmcol};
            $color_num  = $color_num - scalar(@rm_col);
        }
        my $num = $color_num > 10 ? 0 : $color_num;
#        $color_arr = GeneralPlot::fetch_color(-type=>'all',-num=>$color_num,-tag=>'default2');
        $color_arr = GeneralPlot::fetch_color(-type=>'bar',-num=>$num,-tag=>'stack2');
        # perl /Bio/User/hongshimiao/bin/GeneralPlot/v1.0/bin/general_plot.pl bar -file example/match_annotation.stat.xls  -outprefix t -xlab "Sample" -ylab "Percentage" -title "Tag Match annotation Stat" -barfill T -clean F
    }elsif($type eq 'PCA' || $type eq 'pca_point'){
        my $tag = "default";
        if(exists $opts_ref->{group} && $opts_ref->{group} ne 'NA'){
            my ($row,$col) = ReadFile($opts_ref->{group},col=>2,unique=>1);
            $color_num = $row;
        }else{
            my ($row,$col) = ReadFile($opts_ref->{infile});
            $color_num = $col - 1; ## rowname
            $color_num = $row - 1 if($type eq 'pca_point');
        }
        $color_arr = GeneralPlot::fetch_color(-type=>'dot',-num=>$color_num,-tag=>$tag);
    }elsif($type eq 'pheatmap'){
        $color_arr = GeneralPlot::fetch_color(-type=>'heatmap',-num=>0,-tag=>'default');
        $color_num = 3;
    }elsif($type eq 'bar'){
        my ($row,$col) = ReadFile($opts_ref->{infile});
        if(exists $opts_ref->{row} && $opts_ref->{row} eq 'yes'){ ## by rownum 
            $color_num = $row - 1;
        }else{
            $color_num = $col - 1; ## by colnum
        }
        if(exists $opts_ref->{col} && $opts_ref->{col} ne 'all'){
            $color_num = scalar(split/,/,$opts_ref->{col});
        }
        $color_arr = GeneralPlot::fetch_color(-type=>'bar',-num=>$color_num,-tag=>'default');
    }elsif($type eq 'vip_bar'){
        $info{type} = "bar";
        $info{color_num} = 2;
    }elsif($type eq 'diff_zscore'){
        $info{type} = "bar";
        $info{color_num} = 2;
    }elsif($type eq 'box'){
        $info{type} = "box";
        $info{row_col} = "col";
        $info{rowname} = 1;
    }elsif($type eq 'venn'){
        $info{type} = "venn";
        $info{color_num} = scalar(split/,/,$opts_ref->{list});
    }elsif($type eq 'percent_pie'){
        $info{type} = "pie";
        $info{rowname} = 1;
        $info{reverse} = 1;
    }elsif($type eq 'ggcor'){
        $info{type} = "heatmap";
        $info{color_num} = 0;
    }elsif($type eq 'corrplot'){
        $info{type} = "heatmap";
        $info{color_num} = 0;
    }elsif($type eq 'cor'){
        $info{tag} = "corr";
        $info{type} = "heatmap";
        $info{color_num} = 0;
    }elsif($type eq 'fill_bar'){
        $info{type} = "bar";
        $info{row_col} = exists $opts_ref->{row} && $opts_ref->{row} eq 'yes' ? 'row' : 'col';
        $info{tag} = "stack";
        $info{reverse} = 1;
    }elsif($type eq 'violin'){ ## GetColor in Violin
    }elsif($type eq 'volcano'){
        $info{type} = "dot";
        $info{color_num} = "volcano";        
    }elsif($type eq 'fpkm_distribution'){
    }elsif($type eq 'two_pie'){
        $info{type} = "twopie";
        $info{row_col} = "row";
    }elsif($type eq 'histogram_facet'){
        $info{type} = "hist";
        $info{unique} = 1;
        $info{col} = 1;
        $info{row_col} = "row";
    }

    if(exists $info{type}){
        GetColor($opts_ref,%info);
    }else{
        if($color_num > 0){
            my $rel_color_num = @$color_arr;
#            $color_arr = [@$color_arr[$rel_color_num-$color_num..$rel_color_num-1]];
            $color_arr = [@$color_arr[0..$color_num-1]];
            print "Get [$color_num] color from [GeneralPlot::fetch_color] [@$color_arr]\n";
            $opts_ref->{color} = join(",",@$color_arr);
        }
    }

    return;
}

sub GetColor{
    my $opts_ref = shift;
    my %info = @_;
    
    $info{tag} //= "default";
    $info{row_col} //= "col";
    $info{rowname} //= 1;
    $info{reverse} //= 0;
    $info{rm_col_num} //= 0;
    Error("No type info to GetColor, please add type=>xxx") if(!exists $info{type});
    
    $opts_ref->{GetColor}++;
    my $color_num = 0;
    if(exists $info{color_num}){
        $color_num = $info{color_num};
    }else{
        my ($row,$col) = ReadFile($info{file});
        if(exists $info{unique}){
            Error("No col_num info to ReadFile for unique, please add col=>xxx") if(!exists $info{col});
            ($row,$col) = ReadFile($info{file},unique=>$info{unique},col=>$info{col});
            $info{row_col} = "row";
        }
        $color_num = $info{row_col} eq 'col' ? $col : $row;
        $color_num -= $info{rowname};  ## remove rowname num 
        $color_num -= $info{rm_col_num};
    }
    my $color_arr = GeneralPlot::fetch_color(-type=>$info{type},-num=>$color_num,-tag=>$info{tag});
#    print "$color_num\n";
    
    if($color_num eq 'volcano' || $color_num >= 0){
        my $rel_color_num = @$color_arr;
#        $color_arr = [@$color_arr[$rel_color_num-$color_num..$rel_color_num-1]];
#        $color_arr = [@$color_arr[0..$color_num-1]] if($color_num > 0);
        my $add_info = exists $info{file} ? " by file [$info{file}]" : "";
        print "Get [$color_num] color from [GeneralPlot::fetch_color] [@$color_arr]$add_info\n";
        @$color_arr = reverse(@$color_arr) if($info{reverse} == 1);
        $opts_ref->{color} = join(",",@$color_arr);
    }

    return;
}

sub ReadFile{
    my $file = shift;
    my %opts = @_;
    
    $opts{unique} //= 0;
    $opts{col} //= 1;
    $opts{col} -= 1;

    open my $in_fh,"$file" or die "$!:$file";
    my $col;
    my $row = 0;
    my %e;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        $col = @aa;
        if($opts{unique} == 1){
            Error("No col info [$opts{col}] in [@aa] [$file]") if(!defined $aa[$opts{col}]);
            my $info = $aa[$opts{col}];
            if(!exists $e{$info}){
                $row++;
                $e{$info} = 1;
            }
        }else{
            $row++;
        }
    }
    close $in_fh;

    return($row,$col);
}

sub PercentPie{ ## 百分比饼图
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx,col,title","outdir","header","exclude_row","row","color");
    $opts{exclude_row} //= "";
    $opts{row} //= "";
    $opts{color} //= "NA";
    $opts{header} //= "TRUE";
    $opts{outpfx} = "$opts{outdir}/$opts{outpfx}" if(exists $opts{outdir});

    Excu("$Rscript $src/percent_pie.r $opts{infile} $opts{outpfx} $opts{col} $opts{header} \"$opts{title}\" \"$opts{exclude_row}\" \"$opts{row}\" \"$opts{color}\"");


    return "$opts{outpfx}.pdf";
}

sub Stauration{ ## 随机曲线
    my %opts = @_;

    Excu("$Rscript $src/stauration_analysis.R $opts{infile} $opts{window} $opts{sample} $opts{all_reads} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub PQBar{
    my %opts = @_;

    my $Rscript = $apps{'Rscript.3.6.0'};
    $opts{name} //= "KEGG"; #GO
    $opts{top} //= 20;
    $opts{name_key} //= "Gene"; 
    $opts{ylab} //= $opts{name} eq 'KEGG' ? "Pathway" : "GOterm";
    $opts{title} //= "Top $opts{top} of $opts{name} Enrichment";
    $opts{pfc} //= "all"; # for go pfc
    $opts{color} //= "none";

    #Excu("Rscript $RealBin/pq_bar.r $tmp_xls $outpfx.barplot $enrich_type \"$name_key\" $first_num");
    Excu(qq($Rscript $src/pq_bar.r $opts{infile} $opts{outpfx} $opts{top} "$opts{name_key}" "$opts{ylab}" "$opts{title}" $opts{pfc} "$opts{color}"));

    return "$opts{outpfx}.pdf";
}

sub VIPBar{
    my %opts = @_;

    $opts{top} //= 15;
    my $Rscript = $apps{'Rscript.3.6.0'};
    ParaError2(\%opts,"infile,outpfx,xcol,vip_col,name_col,name_list","top","color");
    $opts{top} //= 15;
    
    CheckFile(\%opts,'infile',line=>1);

    Excu(qq($Rscript $src/vip_bar.r $opts{infile} $opts{outpfx} $opts{top} $opts{xcol} $opts{vip_col}  $opts{name_col} "$opts{color}" $opts{name_list}));

    return "";
    return "$opts{outpfx}.pdf";
}

sub DiffZscore{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx,col,name_list","color","top","name");
    $opts{top} //= 20;

    CheckFile(\%opts,'infile',line=>1);
    Excu(qq($Rscript $src/diff_zscore.r $opts{infile} $opts{outpfx} $opts{top} $opts{col} "$opts{color}" $opts{name} $opts{name_list}));

    return "$opts{outpfx}.pdf";
}

sub Gradient{
    my %opts = @_;

    $opts{name} //= "KEGG"; #GO
    $opts{name_key} //= "Gene";
    $opts{ylab} //= $opts{name} eq 'KEGG' ? "Pathway" : "GOterm";
    $opts{title} //= "Top 20 of $opts{name} Enrichment";
    $opts{reverse} //= 0;
    $opts{pfc} //= "all"; ## for GO PFC

    my $Rscript = $apps{'Rscript.3.6.0'};
    #Excu("Rscript $RealBin/Gradient.r $tmp_xls $outpfx.gradient $enrich_type \"$name_key\" $first_num");
    Excu(qq($Rscript $src/Gradient.r $opts{infile} $opts{outpfx} $opts{top} "$opts{name_key}" "$opts{ylab}" "$opts{title}" $opts{reverse} $opts{pfc}));

    return "$opts{outpfx}.pdf";
}

sub EnrichHeatmap{
    my %opts = @_;

    ParaError2(\%opts,"infile,type,outpfx");

    CheckFile(\%opts,'infile',line=>1);
#    return if(CheckExcu($opts{infile}=>"line1"));

    Excu("$Rscript $src/enrichmentHeatmap.r $opts{infile} $opts{type} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub FillBar{ ## stack graph 
    my %opts = @_;

    $opts{title} //= "title";
    $opts{xlab} //= "xlab";
    $opts{ylab} //= "ylab";
    $opts{rmcol} //= 0;
    $opts{color} //= "no";
    $opts{type} //= "fill";  ## fill bar or count bar 
    
    ParaError2(\%opts,"infile,outpfx","title","xlab","ylab","rmcol","color","type");

    Excu(qq($Rscript $src/fill_bar.r $opts{infile} $opts{outpfx} "$opts{title}" "$opts{xlab}" "$opts{ylab}" "$opts{rmcol}" "$opts{color}" "$opts{type}"));

    return "$opts{outpfx}.pdf";
}

sub Pheatmap{
    my %opts = @_; 

    $opts{col} //= "NA";
    $opts{col} = ColInfo($opts{col}); ## 1,3-7
    $opts{color} //= "steelblue,white,red";# "green,black,red"
    $opts{type} //= "normal";
    $opts{name} //= "none";
    $opts{cluster_col} //= "yes";
    $opts{filter_na} //= "no";
    $opts{label_info} //= "none";
    $opts{group} //= "none";
    $opts{show_cluster_tree} //= "all",
    $opts{column_dend_height} //= 10;
    $opts{row_dend_width} //= 10;
    $opts{cluster_row} //= "TRUE";
    $opts{scale} //= "row";
    $opts{fontsize} //= 10;
    $opts{display_numbers} //= "no";
    $opts{title} //= "";
    $opts{group_color} //= "default";
    $opts{cell_width} //= "NA";
    $opts{cell_height} //= "NA";
    $opts{border_color} //= "NA";
    $opts{border_size} //= "NA";
    $opts{legend} //= "TRUE";

    ParaError2(\%opts,"infile,outdir,outpfx","col","type");
    CheckFile(\%opts,'infile',line=>1);
#    return if(CheckExcu($opts{infile}=>"line1"));

#    if($opts{default} == 1){
#        Excu("$apps{perl} $apps{general_plot} cluster -file $opts{infile} -outprefix $opts{outdir}/$opts{outpfx}");
#    }else{
#    print "=>$opts{type}\n";
    if($opts{type} eq 'complex_heatmap'){
        $Rscript = $apps{'Rscript.4.1.2'};
        Excu("$Rscript $src/complex_heatmap.r $opts{infile} $opts{outdir} $opts{outpfx} $opts{col} \"$opts{color}\" $opts{cluster_col} $opts{name} $opts{filter_na} $opts{label_info} $opts{group} $opts{show_cluster_tree} $opts{column_dend_height} $opts{row_dend_width} $opts{cluster_row} $opts{scale} $opts{fontsize} $opts{display_numbers} \"$opts{title}\" \"$opts{group_color}\" $opts{cell_width} $opts{cell_height} \"$opts{border_color}\" $opts{border_size}");
    }elsif($opts{type} eq 'dia'){
        Excu("$Rscript $src/pheatmap.dia.r $opts{infile} $opts{outdir} $opts{outpfx}");
    }else{
        Excu("$Rscript $src/pheatmap.r $opts{infile} $opts{outdir} $opts{outpfx} $opts{col} \"$opts{color}\" $opts{cluster_col} $opts{name} $opts{filter_na} $opts{cluster_row} $opts{scale} $opts{cell_width} $opts{cell_height} $opts{fontsize} $opts{border_color} $opts{legend}");
    }

    return "";
#    return "$opts{outdir}/$opts{outpfx}.pdf";
#    }
}

sub FpkmDistribution{
    my %opts = @_;

    $opts{title} //= "FPKM distribution of all samples";
    $opts{xlab} //= "log10(FPKM) of Gene";
#    $opts{color} //= "none";

    ParaError2(\%opts,"infile,outpfx","title","xlab","color");
    my $tmp_file = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",name_key=>"fpkm,rpkm,tpm,rpm",clean=>1);
    GetColor(\%opts,file=>$tmp_file,type=>"hist",row_col=>'col') if(!exists $opts{color});
#    $tmp_file = $opts{infile};
#    FormatData(\%opts);
    Excu("$Rscript $src/fpkm_distribution.r $tmp_file $opts{outpfx} \"$opts{title}\" \"$opts{xlab}\" \"$opts{color}\"");

    return "$opts{outpfx}.pdf";
}

sub Violin{
    my %opts = @_;

    ## default parameters 
#    $opts{color} //= "rainbow";
    $opts{ymin} //= -999999;
    $opts{ymax} //= 999999;
    $opts{group} //= "none";
    $opts{xlab} //= "Sample";
    $opts{title} //= "Sample Expression Violin Plot";
    $opts{format} //= "log10+1";
    $opts{rowname} //= 1;
    $opts{dot} //= "no";
    $opts{filter0} //= 0;
    $opts{name_key} //= "_fpkm"; ## s/_fpkm//
    ParaError2(\%opts,"infile,outdir,outpfx","title","xlab","ylab","color","ymin","ymax","format","rowname","dot","filter0","name_key");
    my $outpfx = "$opts{outdir}/$opts{outpfx}";

    my $name = $opts{name_key} eq '_rpkm' ? "rpkm" : "fpkm";
    $name = "rpm" if($opts{name_key} eq '_rpm');
    $name = "tpm" if($opts{name_key} eq '_tpm');
    if($opts{format} eq 'log10'){
        $opts{ylab} //= "log10($name)";
    }elsif($opts{format} eq 'no'){
        $opts{ylab} //= "$name";
    }else{
        $opts{ylab} //= "log10($name+1)";
    }


    Excu("perl $perl/violin.pl -i $opts{infile} -outfile $outpfx.4draw.xls -filter0 $opts{filter0} -n $opts{group} -rowname $opts{rowname} -name_key $opts{name_key}");
    GetColor(\%opts,file=>"$outpfx.4draw.xls",unique=>1,col=>3,type=>'box') if(!exists $opts{color});
    $opts{color} //= "rainbow";
    ## Rscript infile outpfx  title xlab ylab color \"$ymin\ \"$ymax\" $dot $format"
    Excu(qq($Rscript $src/violin.r $outpfx.4draw.xls $outpfx "$opts{title}" "$opts{xlab}" "$opts{ylab}" "$opts{color}" "$opts{ymin}" "$opts{ymax}" $opts{dot} "$opts{format}"));
#    Excu("Rscript $src/violin.r $outpfx.4draw.xls $outpfx \"$opts{title}\" \"$opts{xlab}\" \"$opts{ylab}\"  ");

    return "$outpfx.pdf";
}

sub Box{
    my %opts = @_;
    
    $opts{title} //= "";
    $opts{xlab} //= "";
    $opts{ylab} //= "";
    $opts{melt} //= "no";
    $opts{type} //= "all";
    $opts{format} //= "none";
    $opts{color} //= "none";

    ParaError2(\%opts,"infile,outpfx",qw/title xlab ylab melt type format color/);
    if($opts{type} eq 'DIA'){
        Excu("perl $perl/boxplot.pl -i $opts{infile} -outpfx $opts{outpfx}");
        $opts{infile} = "$opts{outpfx}.4draw";
        GetColor(\%opts,file=>$opts{infile},unique=>1,col=>2,type=>'box');
    }

    Excu(qq($Rscript $src/box.r $opts{infile} $opts{outpfx} "$opts{title}" "$opts{xlab}" "$opts{ylab}" $opts{melt} $opts{format} "$opts{color}"));

    return "$opts{outpfx}.pdf";
}

sub Bar{
    my %opts = @_;

    $opts{title} //= "DiffExp Gene Statistics";
    $opts{ylab} //= "Number of Genes";
    $opts{xlab} //= "";
    $opts{top} //= "all";
    $opts{col} //= "all";
    $opts{color} //= "#568FC5";#blue  #F8766D "red";
    $opts{row} //= "no"; ## by row, colname => xlab
    $opts{rowname} //= 1;
    $opts{stat} //= "identity";
    $opts{show_text} //= "yes";
    $opts{break_info} //= "none";

    ParaError2(\%opts,'infile,outpfx','title','ylab','xlab','top','col','color','row',"stat");

    Excu("$Rscript $src/bar.r $opts{infile} $opts{outpfx} \"$opts{title}\" \"$opts{ylab}\" \"$opts{xlab}\" $opts{top} $opts{col} \"$opts{color}\" $opts{row} $opts{rowname} $opts{stat} $opts{show_text} $opts{break_info}");

    return "$opts{outpfx}.pdf";
}

sub DotPlot{
    my %opts = @_;

    $opts{title} //= "Title";
    $opts{ylab} //= "Ylab";
    $opts{xlab} //= "Xlab";
    $opts{col} //= "1,2";
    $opts{color} //= "#568FC5";#blue  #F8766D "red";
    $opts{rowname} //= "NULL"; ## doplot no rowname
    $opts{ymax} //= "none";
    $opts{format} //= "none";

    ParaError2(\%opts,'infile,outpfx',qw/title ylab xlab col color rowname/);
    
    Excu("$Rscript $src/dotplot.r $opts{infile} $opts{outpfx} \"$opts{title}\" \"$opts{ylab}\" \"$opts{xlab}\" $opts{col} \"$opts{color}\" $opts{rowname} $opts{ymax} $opts{format}");

    return "$opts{outpfx}.pdf";
}

sub BreakBar{
    my %opts = @_;

    $opts{title} //= "Break Bar";
    $opts{ylab} //= "ylab";
    $opts{xlab} //= "xlab";
    $opts{color} //= "red";
    $opts{group} //= "none";
    $opts{ratio} //= "0.65"; ## 下半部比例

    ParaError2(\%opts,'infile,outpfx,break_pos','title','ylab','xlab','color','group','ratio');
    CheckFile(\%opts,'infile',line=>1);
#    return if(CheckExcu($opts{infile}=>"line1"));
    my $mean = $opts{infile};
    my $sd = "none";
    if(exists $opts{group} && $opts{group} ne 'none'){
        Excu("perl $perl/pre_break_braplot.pl $opts{infile} $opts{group} $opts{outpfx}");
        ($mean,$sd) = ("$opts{outpfx}.mean.xls","$opts{outpfx}.sd.xls");
    }
    Excu("$Rscript $src/break_barplot.r $mean $sd $opts{outpfx} $opts{break_pos} $opts{ratio} \"$opts{color}\" \"$opts{title}\" \"$opts{xlab}\" \"$opts{ylab}\"");
    
    return "$opts{outpfx}.pdf";
}

sub Volcano{
    my %opts = @_;

    $opts{logformat} //= "no";
    
    ParaError2(\%opts,"infile,outpfx,xcol,ycol,xcut,ycut,title",'logformat',"vip_col","vip_cut","xlab","ylab","size","xmax","ymax");
    $opts{color} //= "#F8766D,#00000032,#619CFF"; # up nosig down
    $opts{vip_col} //= 0;
    $opts{vip_cut} //= 0;
    $opts{xlab} //= "log2(fc)";
    $opts{ylab} //= "auto"; # auto use p/q
    $opts{size} //= 1;
    $opts{xmax} //= "auto";
    $opts{ymax} //= "auto";

    Excu(qq($Rscript $src/volcano.r $opts{infile} $opts{outpfx} $opts{xcol} $opts{ycol} $opts{vip_col} $opts{xcut} $opts{ycut} $opts{vip_cut} "$opts{title}" "$opts{xlab}" "$opts{ylab}" $opts{logformat} "$opts{color}" $opts{size} $opts{xmax} $opts{ymax}));

    return "$opts{outpfx}.pdf";
}

sub TwoPie{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx","color");
    $opts{color} //= "none";

    Excu("$Rscript $src/two_pie.r $opts{infile} $opts{outpfx} \"$opts{color}\"");

    return "$opts{outpfx}.pdf";
}

sub FaLengthDistribution{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx,min,max,by",'title','col',"xlab","ylab");

    $opts{title} //= "Length Distribution";
    $opts{xlab} //= "Length(bp)";
    $opts{ylab} //= "Frequence(#)";
    $opts{col} //= 2;
    Excu("$Rscript $src/length_distribution.r $opts{infile} $opts{outpfx} $opts{min} $opts{max} $opts{by} \"$opts{title}\" \"$opts{xlab}\" \"$opts{ylab}\"  $opts{col}");

    return "$opts{outpfx}.pdf";
}

sub TrendRectProfile{
    my %opts = @_;

    $opts{profile_name} //= "profile";
    # perl /state/partition1/workspace/User/hongshimiao/pipeline/web_src/Small_pipe/trend/src/draw_rect.profileGens.pl -i /state/partition1/workspace/User/hongshimiao/bin/web_src/temp/test/diff_trend/Trend_analysis/profile_stat.xls -o /state/partition1/workspace/User/hongshimiao/bin/web_src/temp/test/diff_trend/Trend_analysis/profile_stat -t Genes_in_Profiele -h yes 
    Excu("perl $perl/draw_rect.profileGens.pl -i $opts{infile} -o $opts{outpfx} -p $opts{profile_name} -t $opts{title}  -h $opts{header}");

    return "$opts{outpfx}.svg";
}

sub TrendProfileFig{
    my %opts = @_;

    $opts{profile_name} //= "profile";
    my $ncol = exists $opts{ncol} ? "-N $opts{ncol}" : "";
    $opts{name_key} //= "gene";
    $opts{pvalue} //=  0.05;

    Excu("perl $apple/profiles_fig.pl $ncol -g $opts{name_key} -v $opts{pvalue} -z $apps{inkscape} -npt -o $opts{outdir} $opts{infile} $opts{profile_name}");

    return;
}

sub TrendProfileLine{
    my %opts = @_;

    $opts{name_key} //= "gene";
    $opts{profile_name} //= "profile";
    $opts{cor} //= "no";
    
    ParaError2(\%opts,"infile,outpfx,trend_genetable","name_key","profile_name","cor");

    Excu("perl $perl/trend_scale_data.pl $opts{infile} $opts{trend_genetable} $opts{outpfx}.draw $opts{profile_name} $opts{name_key}");
    Excu("$Rscript $src/trend_profile_line.r $opts{outpfx}.draw $opts{outpfx} $opts{cor}");

    return "$opts{outpfx}.pdf";
}

sub TrendProfileTable{
    my %opts = @_;

    $opts{name_key} //= "Genes";
    $opts{profile_name} //= "profile";
    $opts{method} //= "log";
    $opts{pvalue} //= 0.05;
    ParaError2(\%opts,"trend_genetable,trend_profiletable,outdir","name_key","profile_name","cor");
    mkdir $opts{outdir} if(!-e $opts{outdir});
    my $outdir = MakeDir("$opts{outdir}/all_$opts{profile_name}s");
    $opts{outpfx} = "$outdir/$opts{profile_name}";
    Excu("perl $perl/draw_trend_analysis.pl $opts{trend_genetable} $opts{outpfx}");
    my %ylab = ("log"=>"log2(V(i)/V(0))","n"=>"V(i) - V(0)","i"=>"V(i)");
    Excu("$Rscript $src/trend_profiletable.r $opts{trend_profiletable} $opts{outpfx} $opts{outpfx} \"$ylab{$opts{method}}\" $opts{pvalue} $opts{profile_name} $opts{name_key}");

    open IN,"$opts{trend_profiletable}" or die $!;
    <IN>;
    my @pdfs;
    while(<IN>){
        chomp;
        my @aa = split/\t/,$_;
        push @pdfs,"$opts{outpfx}$aa[0].pdf";
    }
    close IN;
    
    return join(",",@pdfs);
}

sub GOLevel2{
    my %opts = @_;
    
    die "\n Usage: perl $0 $draw_type -infile -conf -outdir -outpfx\n\n" if(ParaError(\%opts,'outdir','infile','outpfx','conf'));
    CheckFile(\%opts,'infile',line=>1);

    Excu("cd $opts{outdir} && perl $apple/GOdraw.pl --prefix $opts{outpfx} $opts{infile} -conf $opts{conf}");

    return "$opts{outdir}/$opts{outpfx}.GO.level2.bar.svg";
}

sub Venn{
    my %opts = @_;

    $opts{header} //= "F"; 
    $opts{glist} //= "yes";
    $opts{updown} //= "no";
    $opts{scale} //= "no";
    $opts{color} //= "none";
    $opts{soft} //= "venn_diff_perl";
    my $name_list = exists $opts{name} ? "-name \"$opts{name}\"" : "";

    #
    $opts{type} //= "long";
    $opts{outpfx} //= "Venn";
    $opts{show_text} //= "counts"; # counts, percent, all
    $opts{linetype}  //= 1; # 0 => none, 1=>实线, 2=>虚线
    $opts{title}  //= "";
    $opts{alpha}  //= 0.6;

    ParaError2(\%opts,"list,outdir",'header','glist','name',"updown");

    if($opts{soft} eq 'venn_diff_perl'){
        Excu("perl $perl/venn_diff.pl -l $opts{list} $name_list -h $opts{header} -outlist $opts{glist} -updown $opts{updown} -scale $opts{scale} -colors \"$opts{color}\" -outdir $opts{outdir}");
        return "$opts{outdir}/Venn.svg";
    }else{
        my $Rscript = $apps{'Rscript.4.1.2'};
        Excu("perl $perl/venn_stat.pl $opts{list} $opts{type} $opts{outdir}/venn_stat.xls");
        Excu("$Rscript $src/venn_eulerr.r $opts{outdir}/venn_stat.xls $opts{outdir}/$opts{outpfx} $opts{show_text} $opts{linetype} \"$opts{title}\" \"$opts{color}\" $opts{alpha}");
    }

}

sub KeggPathwayClass{
    my %opts = @_;

    $opts{name_key} //= "Genes";
    
    die "\n Usage: perl $0 kegg_PathwayClass -name_key -infile -outpfx \n" if(ParaError(\%opts,'infile','outpfx'));
    CheckFile(\%opts,'infile',line=>1);
#    return if(CheckExcu($opts{infile}=>'line1'));

    Excu("perl $apple/KEGGdraw.pl -xlab \"Number of $opts{name_key}\" --prefix $opts{outpfx} $opts{infile} ");
#    Excu("perl $perl/draw_PathwayClass.pl -k $opts{name_key} -f $opts{infile} -p $opts{outpfx}");

    return "$opts{outpfx}.svg";
}

sub Cor{
    my %opts = @_;

    $opts{cor_cut} //= "0";
    $opts{pvalue_cut} //= "1";
    $opts{method} //= "pearson";
    $opts{name1} //= "Var1";
    $opts{name2} //= "Var2";
    $opts{type} //= "corr";
    $opts{col_or_row} //= "row";
    $opts{color} //= "white,forestgreen";
    $opts{fontsize} //= 14;
    $opts{display} //= "T";
    $opts{border} //= "grey60";
    $opts{showrow} //= "T";
    $opts{showcol} //= "T";
    $opts{cluster_row} //= "F";
    $opts{cluster_col} //= "F";
    $opts{title} //= "NA";
    $opts{soft} //= "ggplot2";

    ParaError2(\%opts,"infile,outpfx",'cor_cut','pvalue_cut','method','name1','name2','type','col_or_row','target');
    my $Rscript_cor = "$apps{Rscript}"; ## 21 3.6.0 no psych
    
    CheckFile(\%opts,'infile',line=>1); 
    $opts{infile} = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",clean=>1);

    my $pheatmap_file = "xxxx";
#    $pheatmap_file = "$opts{outpfx}.pearson.xls";
    if(exists $opts{target}){
        CheckFile(\%opts,'target',line=>2); 
        Excu("$Rscript_cor $src/cor.r $opts{infile} $opts{target} $opts{outpfx} $opts{col_or_row} $opts{method} $opts{name1} $opts{name2} $opts{cor_cut} $opts{pvalue_cut} $opts{type}");
        $pheatmap_file = "$opts{outpfx}.cor.matrix.xls";
    }else{
        Excu("$Rscript_cor $src/cor.r $opts{infile} NA $opts{outpfx}.pearson.xls $opts{col_or_row} $opts{method} $opts{name1} $opts{name2} $opts{cor_cut} $opts{pvalue_cut} $opts{type}");
        $pheatmap_file = "$opts{outpfx}.pearson.xls";
    }

    if($opts{cor_cut} > 0 && exists $opts{target}){
        return "";
    }else{
        $pheatmap_file = abs_path($pheatmap_file);
        my ($outdir,$outpfx) = ("./","$opts{outpfx}");
        if($opts{outpfx} =~ /(.*)\/(.*)/){
            $outdir = $1;
            $outpfx = $2;
        }
        $outpfx = exists $opts{target} ? "$outpfx.cor_heatmap" : "$outpfx.pearson";
        if($opts{col_or_row} eq 'row' && $opts{GetColor} >= 1){
            GetColor(\%opts,type=>'heatmap',color_num=>0);
        }
        ##  outdir outpfx scale cluster_row cluster_col 
        if($opts{soft} eq 'ggplot2'){
            Excu("$Rscript $src/ggplot2_pheatmap.r $pheatmap_file $outdir/$outpfx \"$opts{color}\" \"$opts{title}\"");
        }else{
            Excu("$Rscript $src/pheatmap_display.r $pheatmap_file $outdir $outpfx none \"$opts{color}\" $opts{fontsize} $opts{display} $opts{border} $opts{showrow} $opts{showcol} FALSE FALSE \"$opts{title}\"");
        }
#        Excu("$Rscript $src/ggcor_heatmap.r $pheatmap_file $outdir $outpfx \"$opts{color}\"");
        return "$outdir/$outpfx.pdf";
    }

#    return "$outpfx.pdf";   
}

sub GGCor{
    my %opts = @_;
   
    $opts{type} //= "full";
    $opts{color} //= "white,forestgreen";
    $opts{col_or_row} //= "row";
    $opts{name_col} //= "no";
    ParaError2(\%opts,"infile,outpfx","color",'cor_cut','pvalue_cut','method',"col","type","col_or_row","name_col");
    
    $opts{infile} = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",clean=>1);
    Excu("$Rscript $src/ggcor.r $opts{infile} $opts{outpfx} \"$opts{color}\" $opts{type} $opts{col_or_row} $opts{name_col}");
#    Excu("$Rscript $src/ggcor.r $opts{infile} $opts{outpfx} \"$opts{color}\"");

    return "$opts{outpfx}.pdf";
}

sub Corrplot{
    my %opts = @_;

    my $Rscript = $apps{"Rscript.3.6.3"};
    $opts{color} //= "white,forestgreen";
    $opts{type} //= "full";
    $opts{col_or_row} //= "row";
    $opts{name_col} //= "no";
    $opts{line} //= 0;

    ParaError2(\%opts,"infile,outpfx","color",'cor_cut','pvalue_cut','method',"col","type","col_or_row","name_col","line");
    $opts{infile} = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",clean=>1);
    CheckFile(\%opts,'infile',line=>1); 
    Excu("$Rscript $src/corrplot.r $opts{infile} $opts{outpfx} \"$opts{color}\" $opts{type} $opts{col_or_row} $opts{name_col} $opts{line}");

    return "$opts{outpfx}.pdf";
}

sub PCA{
    my %opts = @_;

    $opts{scale} //= "no";
    $opts{group} //= "NA";
    $opts{type} //= "PCA";
    $opts{color} //= "none";
    $opts{hjust} //= "-1.65";
    $opts{log} //= "no";
    $opts{scale} = "yes" if($opts{scale} eq "T");
    $opts{soft} //= "gmodels";
    
    ParaError2(\%opts,"infile,outpfx",'scale','group',"log","soft");
    CheckFile(\%opts,'infile',line=>2);
#    return if(CheckExcu($opts{infile}=>"line1"));

    if($opts{soft} eq 'gmodels'){
        Excu("$Rscript $src/PCA.r $opts{infile} $opts{outpfx} $opts{scale} $opts{log}");
    }else{
        Excu("$apps{'Rscript.3.3.1'} $src/pls-da.pca.r $opts{infile} $opts{group} $opts{outpfx} $opts{scale}");
    }
    my $file_list = PcaPoint(%opts,infile=>"$opts{outpfx}.PC_data.xls",col=>"1,2",title=>"PCA");

    return $file_list;
}

sub PcaPoint{
    my %opts = @_;
    
    ParaError2(\%opts,"infile,outpfx","color","hjust","col","title","group","add_ellipse");

    $opts{color} //= "none";
    $opts{hjust} //= "-1.65";
    $opts{col} //= "1,2";
    $opts{title} //= "PCA";
    $opts{type} //= "PCA";
    $opts{group} //= "none";
    $opts{add_ellipse} //= "no";

    ##   file outpfx  col_list xlab ylab  title group_info  vjust  point_size text_size 
    Excu("$Rscript $src/point.r $opts{infile} $opts{outpfx}.PCA $opts{col} NA NA $opts{title} $opts{group} $opts{hjust} 0 0 \"$opts{color}\" $opts{add_ellipse}") if($opts{type} =~ /PCA/);
    Excu("$Rscript $src/3d_point.r $opts{infile} $opts{outpfx}.PCA3D 1,2,3 NA NA NA $opts{title} $opts{group} -0.5 3 6 \"$opts{color}\"") if($opts{type} =~ /3D/);

    my @files;
    push @files,("$opts{outpfx}.PCA.pdf","$opts{outpfx}.PCA.nosampleid.pdf") if($opts{type} =~ /PCA/);
    push @files,"$opts{outpfx}.PCA3D.pdf" if($opts{type} =~ /3D/);
#    print "@files\n";

    return join(",",@files);
}

sub COG{
    my %opts = @_;

    $opts{name_key} //= "Gene";
    $opts{cog_kog} //= "COG";

    my $tmp_file = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",name_key=>"none",col=>"2,3,4",clean=>1);
    Excu("$Rscript $src/cog_kog.r $tmp_file $opts{outpfx} \"$opts{name_key}\" \"$opts{cog_kog}\"");

    return "$opts{outpfx}.pdf";
}

sub Density{
    my %opts = @_;

    ParaError2(\%opts,'infile,outpfx');
    Excu("$Rscript $src/density.r $opts{infile} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub Enrich{
    my %opts = @_;

    ParaError2(\%opts,'infile,outpfx,fg_all,bg_all,col');   

    Excu("$Rscript $src/enrich.r $opts{infile} $opts{outpfx}.xls $opts{fg_all} $opts{bg_all} $opts{col}");

    return;
}

sub topGOGraph{
    my %opts = @_;

    ParaError2(\%opts,'infile,outpfx,Rscript,xls,glist,type');
    CheckFile(\%opts,'xls',line=>1);
#    return if(CheckExcu($opts{xls}=>"line1"));

    Excu("$opts{Rscript} $src/topGO_graph.r $opts{infile} $opts{xls} $opts{glist} $opts{type} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub DiffScatterPlot{
    my %opts = @_;

    ParaError2(\%opts,"infile,xlab,ylab,pq_cut,fc_cut,pq_col,fc_col,outpfx");

    Excu("perl $perl/diff_scatter_plot.pl $opts{infile} $opts{xlab} $opts{ylab} $opts{pq_cut} $opts{fc_cut} $opts{pq_col} $opts{fc_col} $opts{outpfx}");

    return "$opts{outpfx}.svg";
}

sub Histogram{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx,col,xlab,ylab");

    Excu("$Rscript $src/histogram.r $opts{infile} $opts{outpfx} $opts{col} \"$opts{xlab}\" \"$opts{ylab}\"");

    return "$opts{outpfx}.pdf";
}

sub HistogramFacet{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx","col","xlab","ylab","title","color");
    $opts{color} //= "none";

    $opts{xlab} = "ratio of editing";
    $opts{ylab} = "Frequency";
    $opts{title} = "Frequency distribution of all samples";

    Excu("$Rscript $src/histogram_facet.r $opts{infile} $opts{outpfx} $opts{col} \"$opts{title}\" \"$opts{xlab}\" \"$opts{ylab}\" \"$opts{color}\"");

    return "$opts{outpfx}.pdf";
}

sub ReplotGSEA{
    my %opts = @_;

    $opts{slope} //= 0.05;
#    $opts{type} //= "xls"; #old is edb
    ParaError2(\%opts,"infile,outpfx,name,title","slope");
    
    if(-e "$opts{infile}.all.xls.gz"){ ## use all.xls.gz
        print "$opts{infile}.all.xls.gz\n";
        Excu("perl $src/pre_RepoltGSEA.pl $opts{infile}.all.xls.gz $opts{name} $opts{outpfx}.tmp");
        my $rank_file = "$opts{infile}/ranked_gene_list.xls";
        if(!-e $rank_file){
            my @files = glob("$opts{infile}/ranked_gene_list*xls");
            $rank_file = $files[0];
        }
        Excu("$Rscript $src/ReplotGSEA_xls.r $opts{outpfx}.tmp $opts{infile}.xls $rank_file $opts{name} \"$opts{title}\" $opts{outpfx}.pdf $opts{slope}");
        Excu("rm -rf $opts{outpfx}.tmp");
    }else{
        Excu("$Rscript $src/ReplotGSEA.r $opts{infile} $opts{name} \"$opts{title}\" $opts{outpfx}.pdf $opts{slope}");
    }

    return "$opts{outpfx}.pdf";
}

sub TernaryPlot{
    my %opts = @_;

    $opts{outpfx} //= "Ternary";
    ParaError2(\%opts,"infile,group,outdir","yaml","outpfx");
    my $yaml = exists $opts{yaml} ? "-y $opts{yaml}" : ""; 
#    my $Rscript = $apps{"Rscript.3.6"};
    my $Rscript = $apps{"Rscript.3.3.1"};
    Excu("$Rscript $src/ternary-plot.R -i $opts{infile} -m $opts{group} -o $opts{outdir} -p $opts{outpfx}  $yaml");

    return "$opts{outdir}/$opts{outpfx}_plot.pdf";
}

sub BuscoFigure{
    my %opts = @_;
    ParaError2(\%opts,"i,outpfx");
    my $Rscript = $apps{'Rscript.3.6.0'};
    Excu("$Rscript $src/busco_figure.R $opts{i} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub SBV{
    my %opts = @_;

    ParaError2(\%opts,"infile,type,outpfx,outdir","conf");
    my $conf = exists $opts{conf} ? "-conf $opts{conf}" : "";
    my $infile = $opts{infile} ne 'none' ? "$opts{infile}" : "";
    Excu("cd $opts{outdir}; perl $sbv $opts{type} $conf -out $opts{outpfx}.svg $infile");

    return "$opts{outpfx}.svg";
}

sub EnrichPlot{
    my %opts = @_;

    $opts{type} //= "KO";
    $opts{pq} //= "q";
    $opts{top} //= 20;
    $opts{name_key} //= "Genes";
    $opts{name} //= "DEGs";
    
    ParaError2(\%opts,"infile,outpfx","diff","pq[q]","type","top[20]");
    if($opts{new} == 1){
        CheckFile(\%opts,'infile',line=>1);
#        return if(CheckExcu($opts{infile}=>'line1'));
        Excu("perl $perl/enrichPlot.pl -t $opts{top} -o $opts{outpfx}.svg -V $apps{sbv} -q $opts{pq} -N $opts{name_key} -T $opts{name} $opts{infile}");
        return "$opts{outpfx}.svg";
    }
    if($opts{type} eq 'GO'){
        Excu("$apps{python2} $python/merge.py $opts{infile} $opts{outpfx}.xls");
        $opts{infile} = "$opts{outpfx}.xls";
    }
    CheckFile(\%opts,'infile',line=>1);
#    return if(CheckExcu($opts{infile}=>'line1'));
    my $diff_info = exists $opts{diff} ? "-d $opts{diff}" : "";
    Excu("perl $apps{enrichPlot} $diff_info -o $opts{outpfx}.svg $opts{infile}");

    if($opts{type} eq 'GO'){
        Clean("$opts{outpfx}.xls");
    }

    return "$opts{outpfx}.svg";
}

sub EnrichBubble{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx","pq_cut","pq","type","top");
    $opts{pq_cut} //= 0.05;
    $opts{pq} //= "Qvalue";
    $opts{type} //= "KO";
    $opts{top} //= 20;

#    $opts{outpfx} = "$opts{outdir}/$opts{name}_bubble";
#    Excu("$apps{python2} $python/draw_Bubble.py $opts{infile} $opts{diff} $opts{outpfx} $opts{type} $opts{pq_cut} $opts{pq}"); 
#    Excu("$Rscript $src/enrich_bubble.r $opts{outpfx}.xls $opts{outpfx} $opts{pq_cut} $opts{pq}");
    Excu("$Rscript $src/enrich_bubble.r $opts{infile} $opts{outpfx} $opts{pq_cut} $opts{pq} $opts{top}");

    my $pdf = "$opts{outpfx}.pdf";
    if($opts{type} eq 'KO' || $opts{type} eq 'GO'){
        $pdf = "$opts{outpfx}.pdf,$opts{outpfx}_sp.pdf"
    }
    return $pdf;
}

sub Random{
    my %opts = @_;

#    ParaError2(\%opts,"infile,outpfx");
#    Excu("$Rscript $src/random.r $opts{infile} $opts{outpfx}");
    $opts{type} //= "sample";
    $opts{format} //= "log10";
    ParaError2(\%opts,"infile,outpfx","type","format");
    if($opts{type} eq 'sample'){
        Excu("$Rscript $src/random.r $opts{infile} $opts{outpfx}");
    }else{
        Excu("$Rscript $src/random_ggplot2.r $opts{infile} $opts{outpfx} random pos value yes $opts{format}");
    }

    return "$opts{outpfx}.pdf";
}

sub PlsdaLoading{
    my %opts = @_;
        
    $opts{name_col} //= 1;
    ParaError2(\%opts,"infile,outpfx");
    Excu(qq($Rscript $src/plsda_loading_plot.r $opts{infile} $opts{outpfx} 2,3 "#F8766D" "NULL" none $opts{name_col}));

    return "$opts{outpfx}.pdf";
}

sub PlsdaPermutation{
    my %opts = @_;

    $opts{color} //= "green,blue";
    ParaError2(\%opts,"infile,outpfx","color");
    
    Excu(qq($Rscript $src/plsda_permutation_plot.r $opts{infile} $opts{outpfx} "$opts{color}"));

    return "$opts{outpfx}.pdf";
}

sub Plsda{
    my %opts = @_; 

    $opts{scale} //= "pareto";
    $opts{type} //= "OPLS-DA";
    $opts{permutation} //= 200;
    
    if($opts{type} eq 'OPLS-DA'){
        $opts{type} = "1,1,1"; # PCA PLS-DA OPLS-DA
    }elsif($opts{type} eq 'PLS-DA'){
        $opts{type} = "1,1,0";
    }
    $opts{scale} = "$opts{scale},$opts{scale},$opts{scale}" if($opts{scale} !~ /,/);

    ParaError2(\%opts,"infile,outpfx,group","scale[pareto]","type[OPLS-DA]","permutation[200]");
    $Rscript = $apps{'Rscript.3.3.1'};

    Excu(qq($Rscript $src/pls-da.r $opts{infile} $opts{group} $opts{outpfx} $opts{scale} $opts{type} $opts{permutation}));

    my @pdfs = map{"$opts{outpfx}.$_"}qw/loading_plot.pdf permutation.pdf plsda.nosampleid.pdf plsda.pdf/;
    if($opts{type} eq 'OPLS-DA' || $opts{type} =~ /,1$/){
        push @pdfs,map{"$opts{outpfx}.$_"}qw/oplsda.nosampleid.pdf oplsda.pdf/;
    }
    return join(",",@pdfs);
}

sub MageckRRA{
    my %opts = @_;

    $opts{xlab} //= "Genes";
    ParaError2(\%opts,"infile,out,title,col","xlab");
    # col 3=>neg|score 4=>neg|p-value 9=>pos|score 10=>pos|p-value

    Excu(qq($Rscript $src/mageck_RRA_score.r $opts{infile} $opts{out} $opts{col} "$opts{title}" "$opts{xlab}"));

    return $opts{out};
}
sub PlotPearson{
    my %opts = @_;

    $opts{name_key} //= "FPKM";

    ParaError2(\%opts,"infile,outdir,sample","name_key");
    Excu(qq($Rscript $src/plotPearson.r $opts{infile} $opts{outdir} $opts{sample} $opts{name_key}));

    return;
}

sub Cluster{
    my %opts = @_;
    
    ParaError2(\%opts,"infile,outpfx","name_key");
    
#    return if(CheckExcu($opts{infile}=>"col3"));
    CheckFile(\%opts,'infile',col=>3);
    
    Excu(qq($Rscript $src/cluster.r $opts{infile} $opts{outpfx}));

    return "$opts{outpfx}.pdf";
}

sub ROC{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx","col","group","diff","top");
    $opts{top} //= 10;
    $opts{clean} //= 0;
    $opts{plot_single} //= "yes";

    my $roc_file = "$opts{outpfx}.tmp";
    CheckFile(\%opts,'infile',line=>1);
    if(exists $opts{col}){
        Excu(qq(perl $perl/pre_roc.pl $opts{infile} $opts{col} $opts{top} >$roc_file));
    }elsif(exists $opts{group} && exists $opts{diff}){
        Excu(qq(perl $perl/pre_roc.pl $opts{infile} none $opts{top} -group $opts{group} -diff $opts{diff} >$roc_file));
    }
    GetColor(\%opts,file=>$roc_file,type=>"hist",row_col=>'col') if(!exists $opts{color});
    Excu(qq($Rscript $src/roc.r -i $roc_file -a -t -s -c "$opts{color}" -p $opts{plot_single} -o $opts{outpfx}));
    push @g_clean_files,"$roc_file" if($opts{clean} == 1);

    return "$opts{outpfx}.roc.pdf";
}

sub MSEAp{
    my %opts = @_;

    $Rscript = $apps{'Rscript.3.6.3_MSEA'}; 
    $opts{color} //= "cm.colors";
    ParaError2(\%opts,"infile,outpfx,glist","color");
    if(exists $opts{ref}){
        Excu("perl $perl/format_MSEAp.pl $opts{infile} $opts{glist} $opts{ref} >$opts{outpfx}.txt");
        $opts{glist} = "$opts{outpfx}.txt";
    }else{
        Excu("perl $perl/format_MSEAp2.pl  $opts{infile} $opts{glist} >$opts{outpfx}.txt");
        $opts{glist} = "$opts{outpfx}.txt";
    }
    Excu(qq($Rscript $src/MSEAp.r $opts{infile} $opts{glist} "$opts{color}" $opts{outpfx}));

    return "$opts{outpfx}.pdf";
}

sub MSEAqea{
    my %opts = @_;

    $opts{lipid} //= "F";
    $opts{class} //= "smpdb_pathway";

    CheckFile(\%opts,'infile',line=>1);
    $Rscript = $apps{'Rscript.3.6.3_MSEA'};
    ParaError2(\%opts,"infile,outdir","lipid","class");
    Error("infile is not exists [$opts{infiel}]") if(!-e $opts{infile});
    Excu(qq($Rscript $src/MSEA_qea.r $opts{infile} name $opts{outdir} $opts{lipid} $opts{class}));

    return;
}

sub MSEAora{
    my %opts = @_;

    $Rscript = $apps{'Rscript.3.6.3_MSEA'};
    ParaError2(\%opts,"infile,outdir");
    Excu(qq($Rscript $src/MSEA_ora.r $opts{infile} name $opts{outdir}));

    Excu("rename oradpi72 ora $opts{outdir}/msea_oradpi72*");

    return;
}

sub iPath{
    my %opts = @_;

    $Rscript = $apps{'Rscript.3.6.3'};
    ParaError2(\%opts,"infile,outpfx","ref");

    if(exists $opts{ref}){
        Excu("perl $perl/format_ipath.pl $opts{infile} $opts{ref} >$opts{outpfx}.txt");
        $opts{infile} = "$opts{outpfx}.txt";
    }

    Excu(qq($Rscript $src/ipath.r $opts{infile} $opts{outpfx}.svg));
    Excu("sed -i s/26/29/g -i $opts{outpfx}.svg");
    Excu("sed -i s/16/18/g -i $opts{outpfx}.svg");
    Excu("$apps{inkscape} -d 300 -e $opts{outpfx}.png $opts{outpfx}.svg");

    return;
}

sub xGSEAdotplot{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx,list","name");

    if(!exists $opts{name}){
        $opts{name} = $opts{list};
        $opts{name} =~ s/.*\///;
        $opts{name} =~ s/\..*//;
    }
    $opts{name} = GetGOName($opts{infile},$opts{name}) if($opts{name} =~ /^GO/);# s/GO_/GO:/; ## GO_0005581.xl
    Excu("$Rscript $src/xGSEAdotplot.r $opts{infile} $opts{name} $opts{list} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub GetGOName{
    my $infile = shift;
    my $name = shift;

    my %name_cc;
    $name_cc{$name} = 1;
    $name =~ s/GO_/GO:/;
    $name_cc{$name} = 1;

    open my $in_fh,"$infile" or die "$!:$infile\n";
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        if(exists $name_cc{$aa[0]}){
            return $aa[0];
        }
    }
    close $in_fh;

    return $name;
}

sub UMAP{
    my %opts = @_;

    $opts{neighbors} //= 15;
    $opts{metric} //= "euclidean";
    $opts{min_dist} //= 0.1;
    $opts{size} //= 1;
    $opts{alpha} //= 1;
    $opts{title} //= "UMAP Plot";
    $opts{legend} //= "T";
    $opts{do_label} //= "T";
    $opts{label_size} //= 5;
    $opts{byrow} //= "T";
    $Rscript = $apps{'Rscript.3.6.3'};
    my $color = $opts{color} ne 'NA' ? "--point_color \"$opts{color}\"" : "";

    ParaError2(\%opts,"infile,outdir,group","color","neighbors","metric","min_dist","size","alpha","title","legend","do_label","label_size","byrow");
    Excu("$Rscript $src/OmicUmapTools.r --expr $opts{infile} --group $opts{group} --neighbors $opts{neighbors} --metric $opts{metric} --min_dist $opts{min_dist} --point_size $opts{size} --point_alpha $opts{alpha} $color --title \"$opts{title}\" --legend $opts{legend} --do_label $opts{do_label} --label_size $opts{label_size} --byrow $opts{byrow} --outdir $opts{outdir}");

    return "$opts{outdir}/UMAP.pdf";
}

sub GGwordcloud{
    my %opts = @_;

    $opts{color} //= "#0077c1,#00a99e,#6bc72b,#ff5a20,#ff1620,#752995";#"random-dark";# random-light
    $opts{color} = "#0077c1,#00a99e,#6bc72b,#ff5a20,#ff1620,#752995" if($opts{color} eq 'default');
    $opts{top} //= "all";
    $opts{type} //= "auto"; # horizon
    $opts{ellipticity} //= 0.65;
    $opts{seed} //= 10;

    $Rscript = $apps{'Rscript.4.1.2'};
    ParaError2(\%opts,"infile,outpfx","color","top","top");
    
    Excu("$Rscript $src/ggwordcloud.r $opts{infile} $opts{outpfx} $opts{top} $opts{type} \"$opts{color}\" $opts{ellipticity} $opts{seed}");

    return "$opts{outpfx}.pdf";
}

sub GGseqlogo{
    my %opts = @_;

    $opts{first_pos} //= 1;
    $opts{range} //= 'all';
    $opts{ignore_case} //= "yes";
    $opts{units} //= "bits"; # prob
    $opts{color} //= "auto";
    $opts{color_info} //= "none";
    $opts{sequence_type} //= "auto";

#    if($opts{color} eq 'auto' && $opts{sequence_type} eq 'protein'){
#        $opts{color} = 
#        $

    $Rscript = $apps{'Rscript.3.6.3'};
    ParaError2(\%opts,"infile,outpfx","first_pos","range","ignore_case","units","color");
    ($opts{color},$opts{color_info}) = FormatColor($opts{color},$opts{color_info}) if($opts{color_info} ne 'none');
    Excu("perl $perl/seq2logoMatrix.pl -first_pos $opts{first_pos} -range $opts{range} -ignore_case $opts{ignore_case} $opts{infile} $opts{outpfx}.pfsm.matrix.xls ");
    Excu("$Rscript $src/ggseqlogo.r $opts{outpfx}.pfsm.matrix.xls $opts{outpfx} $opts{units} \"$opts{color}\" $opts{color_info} ");

    return "$opts{outpfx}.pdf";
}

sub GGraph{
    my %opts = @_;

    $opts{color} //= "none";
    $opts{name} //= "auto";
    $opts{back_color} //= "black";
    $opts{type} //= "ratius"; # area
    $opts{order} //= "none"; # none ascending descending

    $Rscript = $apps{'Rscript.3.6.0'};
    ParaError2(\%opts,"infile,outpfx","color","name");

    Excu("perl $perl/prepare_circlepacker.pl -color \"$opts{color}\" -name $opts{name} $opts{infile} $opts{outpfx}.edges.xls $opts{outpfx}.nodes.xls");
    Excu("$Rscript $src/circlepacker.r $opts{outpfx}.edges.xls $opts{outpfx}.nodes.xls $opts{outpfx} \"$opts{back_color}\" $opts{type} $opts{order}");

    return;
}

sub GGpubr{
    my %opts = @_;

    $opts{type} //= "wilcox.test";
    $opts{title} //= "ggpubr plot";
    $opts{xlab} //= "types";
    $opts{ylab} //= "values";
    $opts{col_or_row} //= "0"; # 0=>col 1=>row
    $opts{dot} //= "none";
    $opts{dot_color} //= "none";
    $opts{dot_alpha} //= 1;
    $opts{error_bar} //= "sd"; #se sd ci
    
    $Rscript = $apps{'Rscript.3.6.0'};
    ParaError2(\%opts,"infile,outdir","type","title","xlab","ylab");
    Excu("$Rscript $src/ggpubr_plot.R $opts{infile} $opts{type} $opts{outdir} \"$opts{title}\" \"$opts{xlab}\" \"$opts{ylab}\" $opts{col_or_row} $opts{dot} \"$opts{dot_color}\" $opts{dot_alpha} $opts{error_bar}");

    return ;
}

sub LetterSig{
    my %opts = @_;

    $opts{title} //= "";
    $opts{col_or_row} //= "0"; # 0=>col 1=>row
    $opts{title_size} //= 15;
    $opts{lab_size} //= 15;
    $opts{axis_size} //= 10;
    $opts{text_size} //= 10;
    $opts{dot} //= "none";
    $opts{dot_color} //= "none";
    $opts{dot_alpha} //= 1;
    $opts{xlab} //= "";
    $opts{ylab} //= "";
    $opts{error_bar} //= "sd"; #se sd ci

    $Rscript = $apps{'Rscript.3.6.0'};
    ParaError2(\%opts,"infile,outdir","title","col_or_row","title","title_size","lab_size","axis_size","text_size");

    Excu("$Rscript $src/cld_plot.R -i $opts{infile} -t \"$opts{title}\" --xlab \"$opts{xlab}\" --ylab \"$opts{ylab}\" --row $opts{col_or_row} --ts $opts{title_size}  --ls $opts{lab_size}  --as $opts{axis_size}  --size $opts{text_size} --dot $opts{dot} --dot_color \"$opts{dot_color}\" --dot_alpha $opts{dot_alpha} --error_bar $opts{error_bar} -o  $opts{outdir}");

    return;
}

sub FormatColor{
    my $color = shift;
    my $color_info = shift;

    my @colors = split/,/,$color;
    my @color_infos = split/;/,$color_info;

    my @o_colors;
    my @o_infos;
    for my $i(0..$#color_infos){
        my @infos = split/,/,$color_infos[$i];
        push @o_infos,@infos;
        my $num = scalar(@infos);
        push @o_colors,($colors[$i])x$num;
    }

    my $o_color = join(",",@o_colors);
    my $o_info = join(",",@o_infos);

    return ($o_color,$o_info);
}

sub Fdr{
    my %opts = @_;

    ParaError2(\%opts,"infile,outfile","type");
    $opts{type} //= "BH";

    Excu("$Rscript $src/fdr.r $opts{infile} $opts{type} $opts{outfile}");

    return;
}

sub GetPartCol{
    my %opts = @_;

    ParaError2(\%opts,"infile,outfile","name_key","clean","col");
    $opts{clean} //= 0;
    $opts{name_key} //= "none";
    $opts{col} //= "none";
    $opts{line} //= "none";

    return $opts{infile} if($opts{col} eq 'none' && $opts{line} eq 'none' && $opts{name_key} eq 'none');

    Excu("perl $perl/get_exp.pl $opts{infile} $opts{outfile} -name_key $opts{name_key} -col $opts{col} -line $opts{line}");
    push @g_clean_files,"$opts{outfile}" if($opts{clean} == 1);

    return $opts{outfile};
}

sub Example{
    my %opts = @_;

    ParaError2(\%opts,"type,outdir","version=s");
    my @files;    
    my $file_dir = "$RealBin/example";
    $file_dir = "/Bio/Bin/pipeline/Small_pipe/R_plot_test/";
    my $type = $opts{type};
    my $outdir = $opts{outdir};
    my $version= exists $opts{version} ? $opts{version} : $outdir;
    my $g_plot = "perl /Bio/User/liuyubin/pipeline/GeneralPlot/bin/general_plot.pl";
    $g_plot = "perl /Bio/Bin/pipeline/System_Programs/GeneralPlot/v1.0/bin/general_plot.pl";
#    $g_plot = "perl /home/hongshimiao/GeneralPlot/bin/general_plot.pl";
    chdir($file_dir);
    if($type eq 'all' || $type eq 'pq_bar'){
        Excu(qq(perl $0 pq_bar -infile pq_bar/qvalue1.xls -name_key "Proteins" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/pq_bar1.barplot.$version));
        Excu(qq(perl $0 pq_bar -infile pq_bar/qvalue2.xls -name_key "Proteins" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/pq_bar2.barplot.$version));
        Excu(qq(perl $0 pq_bar -infile pq_bar/qvalue3.xls -name_key "Proteins" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/pq_bar3.barplot.$version));
        Excu(qq(perl $0 pq_bar -infile pq_bar/qvalue4.xls -name_key "Proteins" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/pq_bar4.barplot.$version));
    }
    if($type eq 'all' || $type eq 'Gradient'){
        Excu(qq(perl $0 Gradient -infile pq_bar/qvalue1.xls -reverse 0 -name_key "Genes" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/Gradient1.$version));
        Excu(qq(perl $0 Gradient -infile pq_bar/qvalue2.xls -reverse 0 -name_key "Genes" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/Gradient2.$version));
        Excu(qq(perl $0 Gradient -infile pq_bar/qvalue3.xls -reverse 0 -name_key "Genes" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/Gradient3.$version));
        Excu(qq(perl $0 Gradient -infile pq_bar/qvalue4.xls -reverse 0 -name_key "Genes" -top 20 -ylab GOterm -title "Top 20 of GO Enrichment" -outpfx $outdir/Gradient4.$version));
    }
    if($type eq 'all' || $type eq 'fill_bar'){
        my $color = "";
        $color = "-color \"#EE1C27,#A3218E,#F58225,#FEF104,#2ECDB3\"" if($version eq 'v1.7');
        Excu(qq(perl $0 fill_bar -infile fill_bar/all.filter_stat.xls -title "Reads Filter Stat" -xlab "Sample" -ylab "Percentage" -rmcol "1" $color -type fill -outpfx $outdir/fill_bar.$version));
        $color = "";
        $color = "-color \"#619CFF,#00BA38,#F8766D\"" if($version eq 'v1.7');
        Excu(qq(perl $0 fill_bar -infile fill_bar/all.region_stat.xls -title "Reads Align Region" -xlab "Sample" -ylab "Percentage" $color -outpfx $outdir/fill_bar2.$version));
        Excu("$g_plot bar -file fill_bar/all.region_stat.xls -barfill T -outprefix $outdir/fill_bar2.g.$version");

        Excu(qq(perl $0 fill_bar -infile fill_bar/match_annotation.stat.xls -type fill -xlab "Sample" -ylab "Percentage" -title "Tag Match annotation Stat" -outpfx $outdir/fill_bar3.$version));
        Excu("$g_plot bar -file fill_bar/match_annotation.stat.xls -barfill T -outprefix $outdir/fill_bar3.g.$version");
#        Excu(qq($g_plot bar -file fill_bar/all.filter_stat.xls -barfill T -outprefix $outdir/fill_bar.2.$version));
#        Excu(qq($g_plot bar -file fill_bar/all.region_stat.xls -barfill T -outprefix $outdir/fill_bar2.2.$version));
    }
    if($type eq 'all' || $type eq 'PCA'){
        Excu(qq(perl $0 PCA -scale yes -type PCA,3D -group pca/all_mirna.group -infile pca/all_mirna.exp.filter.xls -outpfx $outdir/pca1.$version));
        Excu(qq($g_plot pca -file pca/all_mirna.exp.filter.xls -group pca/all_mirna.group -label T -outprefix $outdir/pca1.g.$version));
        Excu(qq(perl $0 PCA -scale yes -type PCA,3D -infile pca/all_mirna.exp.filter.xls -outpfx $outdir/pca2.$version));
        Excu(qq(perl $0 PCA -scale yes -type PCA,3D -group pca/all.exp.xls.group -infile pca/all.exp.xls -outpfx $outdir/pca3.$version));
        Excu(qq($g_plot pca -file pca/all.exp.xls -group pca/all.exp.xls.group -label T -outprefix $outdir/pca3.g.$version));
        
        Excu(qq(perl $0 PCA -scale yes -type PCA,3D -infile pca/all_mirna.exp.2.xls -outpfx $outdir/pca4.$version)); ## PCA3 too small
    }
    if($type eq 'all' || $type eq 'pheatmap'){
        Excu(qq(perl $0 pheatmap -infile heatmap/all_mirna.exp.filter.xls -outdir $outdir -outpfx heatmap1.$version));
        Excu("$g_plot cluster -file heatmap/all_mirna.exp.filter.xls -outprefix $outdir/heatmap1.g.$version");
    }

    if($type eq 'all' || $type eq 'bar'){
        Excu(qq(perl $0 bar -infile bar/diff.stat.xls  -ylab "Number of miRNAs" -title "DiffExp miRNA Statistics"  -outpfx $outdir/bar1.$version));
        Excu("$g_plot bar -file bar/diff.stat.xls -bardodge T -outprefix $outdir/bar1.g.$version");

        Excu(qq(perl $0 bar -infile bar/circ_candidates.info -stat count -col 2 -outpfx $outdir/bar2.$version));
    
        Excu(qq(perl $0 bar -infile bar/Nr.species.stat.xls -ylab "Number of Unigenes" -title "" -top 10 -outpfx $outdir/bar3.$version));
    }

    if($type eq 'all' || $type eq 'box'){
        Excu(qq(perl $0 box -infile box/all_Peptides.xls -type DIA -title Boxplot -xlab Samples -ylab "log10(quantity)" -format log10 -outpfx $outdir/box1.$version));
        Excu(qq(perl $0 box -infile box/all.genes.expression.xls -title Boxplot -xlab Samples -ylab "log10(fpkm)" -format log10 -melt yes -outpfx $outdir/box2.$version));
        Excu("$g_plot box -file box/all.genes.expression.xls -outprefix $outdir/box2.g.$version -scale log10"); # -namefix _fpkm -filter0 T
    }
    if($type eq 'all' || $type eq 'percent_pie'){
        Excu(qq(perl $0 percent_pie -infile percent_pie/Gre2.coverage.xls -header T -col 2 -title "Distribution of Genes' Coverage" -outdir $outdir -outpfx percent_pie1.$version));
        Excu(qq(perl $0 percent_pie -infile percent_pie/protein_cov.xls -header T -col 2 -title "Distribution of Genes' Coverage" -outdir $outdir -outpfx percent_pie2.$version));
        Excu("$g_plot pie -file percent_pie/protein_cov.xls -outprefix $outdir/percent_pie2.g.$version");
    }

    if($type eq 'all' || $type eq 'venn'){   
        Excu(qq(perl $0 venn -list venn/g1.list,venn/g2.list -name GO,KEGG -header F -outdir $outdir/venn2.$version));
        Excu("$g_plot venn -file venn/venn2.dat.txt -outprefix $outdir/venn2.g.$version");
        
        Excu(qq(perl $0 venn -list venn/g1.list,venn/g2.list,venn/g3.list -name GO,KEGG,COG -header F -outdir $outdir/venn3.$version));
        Excu("$g_plot venn -file venn/venn3.dat.txt -outprefix $outdir/venn3.g.$version");
        
        Excu(qq(perl $0 venn -list venn/g1.list,venn/g2.list,venn/g3.list,venn/g4.list -name GO,KEGG,COG,NR -header F -outdir $outdir/venn4.$version));
        Excu("$g_plot venn -file venn/venn4.dat.txt -outprefix $outdir/venn4.g.$version");
    }

    if($type eq 'all' || $type eq 'cor'){
        my $i = 0;
        for my $exp("cor/all.samples.exp.xls","cor/all.samples.exp.18.xls","cor/all.samples.exp.36.xls"){
            $i++;
            Excu(qq(perl $0 cor -infile $exp -col_or_row col -outpfx $outdir/cor$i.1.$version));
            Excu(qq(perl $0 ggcor -infile $exp -outpfx $outdir/cor$i.2.$version));
            Excu(qq(perl $0 cor -infile $exp -col_or_row col -soft pheatmap -outpfx $outdir/cor$i.3.$version));
            Excu("$g_plot correlation -file $exp -outprefix $outdir/cor$i.1.g.$version");
            Excu("$g_plot correlation -file $exp -label T -cortest T -postype upper -outprefix $outdir/cor$i.2.g.$version");
            Excu("$g_plot correlation_old -file $exp -label T -cortest T -postype upper -outprefix $outdir/cor$i.3.g.$version");
        }
    }

    if($type eq 'all' || $type eq 'cluster'){
        Excu(qq(perl $0 cluster -infile cluster/all.samples.exp.xls -outpfx $outdir/cluster1.$version));
    }

    if($type eq 'all' || $type eq 'violin'){
        Excu(qq(perl $0 violin -infile violin/all.isoforms.expression.1000.xls -format log10 -name_key _fpkm -outdir $outdir -filter0 1 -outpfx violin.$version));
        Excu("$g_plot violin -file violin/all.isoforms.expression.1000.xls -filter0 T -namefix _fpkm -scale log10 -xlab sample -ylab \"log10(fpkm)\" -title \"Violin Plot\" -outprefix $outdir/violin.g.version");
    }

    if($type eq 'volcano' || $type eq 'all'){
        Excu(qq(perl $0 volcano -infile volcano/CSH-vs-DPH.xls -xcol 2 -ycol 3 -xcut 2 -ycut 0.05 -title "Volcano Plot" -outpfx $outdir/volcano1.$version ));
        Excu("$g_plot volcano -file volcano/CSH-vs-DPH.xls -outdir $outdir -outname volcano1.g.$version");
    }

    if($type eq 'all' || $type eq 'fpkm_distribution'){
        Excu(qq(perl $0 fpkm_distribution -infile fpkm_distribution/all.isoforms.expression.1000.xls -outpfx $outdir/fpkm_distribution.$version));
        Excu("$g_plot density -file fpkm_distribution/all.isoforms.expression.1000.fpkm.xls -format matrix -colnpg T -clean F -inline T -scale log10 -filter0 T -outdir $outdir  -outname fpkm_distribution.g.$version");
    }

    if($type eq 'all' || $type eq 'two_pie'){
        Excu(qq(perl $0 two_pie -infile two_pie/all.trans_stat.xls -outpfx $outdir/two_pie1.$version));
        Excu("$g_plot twopie -file two_pie/twopie.1.xls  -outprefix $outdir/two_pie1.g.$version -label T");
        
        Excu(qq(perl $0 two_pie -infile two_pie/all.trans_stat.2.xls -outpfx $outdir/two_pie2.$version));
    }

    if($type eq 'all' || $type eq 'histogram_facet'){
        Excu(qq(perl $0 histogram_facet -infile histogram_facet/editing.xls -col frequency -outpfx $outdir/histogram_facet1.$version));
        Excu("$g_plot histgram -file histogram_facet/editing.xls -outprefix $outdir/histogram_facet1.g.$version -colname frequency -facet T");
        Excu(qq(perl $0 histogram_facet -infile histogram_facet/editing.36.xls -col frequency -outpfx $outdir/histogram_facet2.$version));
        Excu("$g_plot histgram -file histogram_facet/editing.36.xls -outprefix $outdir/histogram_facet2.g.$version -colname frequency -facet T");
    }

    if($type eq 'all' || $type eq 'ReplotGSEA'){
        Excu(qq(perl $0 ReplotGSEA -infile GSEA -name KO01210 -title "Enrichment plot:KO01210 2-Oxocarboxylic acid metabolism" -outpfx $outdir/ReplotGSEA1.$version));
    }

    if($type eq 'all' || $type eq 'plsda'){
        Excu(qq(perl $0 plsda -infile plsda/NC-vs-mp1p.tmp -group plsda/NC-vs-mp1p.group.txt -outpfx $outdir/plsda1.$version -scale pareto -type OPLS-DA));
    }

    if($type eq 'all' || $type eq 'plsda_loading'){
        Excu(qq(perl $0 plsda_loading -infile plsda_loading_plot/plsda1.v1.9.loading_plot.xls -outpfx $outdir/plsda_loading_plot1.$version));
    }

    if($type eq 'enrichPlot' || $type eq 'all'){
        Excu(qq(perl $0 enrichPlot -infile bar_Gradient/reactome_bar_Gradient.xls  -pq q -outpfx $outdir/enrichPlot1.$version -name_key Genes -name DGEs -new 1));
    }

#    print @{GeneralPlot::fetch_color(-type=>"venn")};

    return "";
}


sub Converts{
    my $file_list = shift;

    return if(!defined $file_list);
    my @files = split/,/,$file_list;
    for my $file(@files){
        Convert($file);
    }
}

sub Convert{
    my $file = shift;

    return if(!defined $file);
    if(!-e $file){
        print "[$file] is not exist\n";
    }
    return if(!-e $file);

    my $outfile = $file;
    if($outfile =~ /.pdf$/){
        $outfile =~ s/.pdf$/.png/;
        Excu("convert -density $opts{density} $file $outfile");
    }elsif($outfile =~ /.svg$/){
        $outfile =~ s/.svg$/.png/;
        Excu("$apps{inkscape} --export-background=white -d 300 -f $file --export-png=$outfile");
    }

    return;
}

sub ParaError{
    my $opts_ref = shift;
    my @keys = @_;

    my $err = "";
    for my $key(@keys){
        $err .= "Error: No info of [$key]\n" if(!exists $opts_ref->{$key});
    }
    
    if($err eq ''){
        return 0;
    }else{
        print STDERR "$err";
        return 1;
    }   
}

sub ParaError22222{
    my $opts_ref = shift;
    my $need_list = shift;
    my @others = @_;

    my @needs = split/,/,$need_list;
    my $err = "";
    for my $key(@needs){
        $err .= "Error: No info of [$key]\n" if(!exists $opts_ref->{$key});
    }
    if($err eq ''){
        return 0;
    }else{
        print STDERR "$err";
        Error("\n\tUsage: perl $0 $draw_type -[@needs] optional : -[@others]");
    }
}

sub FormatData{
    my $opts_ref = shift;

    if($opts_ref->{melt} eq 'yes'){
        Melt($opts_ref->{infile},"$opts_ref->{outpfx}.tmp");
    }


    return;
}

sub Melt{
    my $infile = shift;
    my $outfile = shift;
    my %opts = @_;

    open IN,$infile or die $!;
    my @head;
    my @key_col = (0);
    my @value_col;
    open OUT,">$outfile" or die $!;
    print OUT "id\ttype\tvalue\n";
    while(<IN>){
        chomp;
        my @aa = split/\t/,$_;
        if($. == 1){
            @head = @aa;
            @value_col = (1..$#aa);
        }else{
            my $key = join("\t",@aa[@key_col]);
            for my $i(@value_col){
                print OUT join("\t",$key,$head[$i],$aa[$i])."\n";
            }
        }
    }
    close IN;
    close OUT;


    return;
}


sub ColInfo{
    my $col = shift;

    return $col if($col eq 'NA');

    my @cols;
    for my $info(split/,/,$col){
        if($info =~ /(\d+)-(\d+)/){
            push @cols,($1..$2);
        }else{
            push @cols,$info;
        }
    }

    return join(",",@cols);
}

sub CheckExcu{
    my %info = @_;

    my $err;
    for my $file(sort keys %info){
        if($info{$file} =~ /line(\d+)/){
            if(CheckFile($file,line=>$1,die=>0)==0){
                $err .= "file [$file] is empty\n";
            }
        }
        if($info{$file} =~ /col(\d+)/){
            if(CheckFile($file,col=>$1,die=>0) == 0){
                $err .= "file [$file] col_num is too less, must more than [$1]\n";
            }
        }
    }
    
    if(defined $err){
        print STDERR "Warings:$err";
        return 1;
    }else{
        return 0;
    }
}

sub Clean{
    my @files = @_;

    for my $file(@files){
        Excu("rm -f $file");
    }

    return;
}

sub CheckFile{
    my $opts_ref = shift;
    my $key = shift;
    my %opts = @_;

    my $file = $opts_ref->{$key};
    $opts{die} //= 1;

    my $ok = 0;
    my $error;
    if(exists $opts{line} || exists $opts{col}){
        open FILE,$file or die "$!:$file";
        my $i = 0;
        my $col = 0;
        while(<FILE>){
            $i++;
            if($. == 1 && exists $opts{col}){
                my @aa = split/\t/,$_;
                $col = @aa;
                print "@aa\t$opts{col}\n";
                if(@aa > $opts{col}){
                    $ok = 1;
                }else{
                    $ok = 0;
                }
            }
            if(exists $opts{line}){
                if($i > $opts{line}){
                    $ok = 1;
                }else{
                    $ok = 0;
                }
            }
        }
        close FILE;
        $error = "[$key] file [$file] line [$i] <= [$opts{line}]" if(exists $opts{line} && $i <= $opts{line});
        $error = "[$key] file [$file] col [$col] <= [$opts{col}]" if(exists $opts{col} && $col <= $opts{col});
    }
    if($ok == 0){
        print STDERR "Warnings: $error\n";
        exit(0) if($opts{die} == 1);
    }

    return $ok;
}

sub MakeDir2{
    my $dir = shift;

    mkdir $dir if(!-e $dir);

    return $dir;
}

sub Excu{
    my $cmd = shift;
    
    print "$cmd\n";

    if($opts{run} == 1){
        my $ret = system($cmd);
        Error("[$ret] [$cmd]") if($ret > 0);
    }

    return;
}

sub Error{
     
    confess "FATAL ERROR: @_\n";

}
