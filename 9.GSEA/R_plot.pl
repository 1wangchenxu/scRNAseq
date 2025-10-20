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

            PCA
            cor
            ggcor
            enrich
            topGO_graph
            plsda_loading
            plsda
            mageck_RRA
            plotPearson
            cluster

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
                "exclude_row=s","outfile=s",
                "title=s","xlab=s","ylab=s","dot=s","ymin=s","ymax=s","color=s","format=s","indir=s",
                "min=s", "max=s", "by=s","name=s",
                "xcol=s","ycol=s","xcut=s","ycut=s",
                "melt=s","rowname=s",
                "window=s","all_reads=s","sample=s","type=s","name_key=s","top=s","group=s","rmcol=s",
                "trend_genetable=s","profile_name=s","trend_profiletable=s", # Trend
                "conf=s",
                "list=s","glist=s",# venn "header=s","o=s"
                "target=s","cor_cut=s","pvalue_cut=s","name1=s","name2=s","col_or_row=s",  # cor 
                "fontsize=s","method=s","display=s","border=s","showrow=s","showcol=s", #cor
                "logformat=s", ## log2 yes or not 
                "cog_kog=s", ## annot 
                "xls=s","log=s",
                "slope=s",
                "fg_all=s","bg_all=s","diff=s","pfc=s",#enrich 
                "pq=s","pq_cut=s","fc_cut=s","pq_col=s","fc_col=s",## diff scatter plot 
                "updown=s","scale=s",## venn
                "break_pos=s","ratio=s", # break_barplot
                "name_col=s", # loading plot
                "filter0=s","show_text=s",
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
    plsda => \&Plsda,
    mageck_RRA => \&MageckRRA,
    plotPearson => \&PlotPearson,
    cluster => \&Cluster,
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
        if(exists $opts_ref->{group}){
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

    open my $in_fh,"$file" or die "$!:$file\n";
    my $col;
    my $row = 0;
    my %e;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        $col = @aa;
        if($opts{unique} == 1){
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


    #Excu("Rscript $RealBin/pq_bar.r $tmp_xls $outpfx.barplot $enrich_type \"$name_key\" $first_num");
    Excu(qq($Rscript $src/pq_bar.r $opts{infile} $opts{outpfx} $opts{top} "$opts{name_key}" "$opts{ylab}" "$opts{title}" $opts{pfc}));

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

    return if(CheckExcu($opts{infile}=>"line1"));

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
    
    ParaError2(\%opts,"infile,outdir,outpfx","col","type");
    return if(CheckExcu($opts{infile}=>"line1"));

#    if($opts{default} == 1){
#        Excu("$apps{perl} $apps{general_plot} cluster -file $opts{infile} -outprefix $opts{outdir}/$opts{outpfx}");
#    }else{
    if($opts{type} eq 'dia'){
        Excu("$Rscript $src/pheatmap.dia.r $opts{infile} $opts{outdir} $opts{outpfx}");
    }else{
        Excu("$Rscript $src/pheatmap.r $opts{infile} $opts{outdir} $opts{outpfx} $opts{col} \"$opts{color}\"");
    }

    return "";
#    return "$opts{outdir}/$opts{outpfx}.pdf";
#    }
}

sub FpkmDistribution{
    my %opts = @_;

    $opts{title} //= "FPKM distribution of all samples";
    $opts{xlab} //= "log10(FPKM) of Gene";
    $opts{color} //= "none";

    ParaError2(\%opts,"infile,outpfx","title","xlab","color");
    my $tmp_file = GetPartCol(%opts,outfile=>"$opts{outpfx}.tmp",name_key=>"fpkm,rpkm,tpm,rpm",clean=>1);
    GetColor(\%opts,file=>$tmp_file,type=>"hist",row_col=>'col');
#    $tmp_file = $opts{infile};
#    FormatData(\%opts);
    Excu("$Rscript $src/fpkm_distribution.r $tmp_file $opts{outpfx} \"$opts{title}\" \"$opts{xlab}\" \"$opts{color}\"");

    return "$opts{outpfx}.pdf";
}

sub Violin{
    my %opts = @_;

    ## default parameters 
    $opts{color} //= "rainbow";
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
    if($opts{format} eq 'log10'){
        $opts{ylab} //= "log10($name)";
    }elsif($opts{format} eq 'no'){
        $opts{ylab} //= "$name";
    }else{
        $opts{ylab} //= "log10($name+1)";
    }


    Excu("perl $perl/violin.pl -i $opts{infile} -outfile $outpfx.4draw.xls -filter0 $opts{filter0} -n $opts{group} -color \"$opts{color}\" -rowname $opts{rowname} -name_key $opts{name_key}");
    GetColor(\%opts,file=>"$outpfx.4draw.xls",unique=>1,col=>3,type=>'box');
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

    ParaError2(\%opts,'infile,outpfx','title','ylab','xlab','top','col','color','row',"stat");

    Excu("$Rscript $src/bar.r $opts{infile} $opts{outpfx} \"$opts{title}\" \"$opts{ylab}\" \"$opts{xlab}\" $opts{top} $opts{col} \"$opts{color}\" $opts{row} $opts{rowname} $opts{stat} $opts{show_text}");

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

    ParaError2(\%opts,'infile,outpfx',qw/title ylab xlab col color rowname/);
    
    Excu("$Rscript $src/dotplot.r $opts{infile} $opts{outpfx} \"$opts{title}\" \"$opts{ylab}\" \"$opts{xlab}\" $opts{col} \"$opts{color}\" $opts{rowname} $opts{ymax}");

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
    return if(CheckExcu($opts{infile}=>"line1"));
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
    
    ParaError2(\%opts,"infile,outpfx,xcol,ycol,xcut,ycut,title",'logformat');
    $opts{color} //= "#F8766D,#00000032,#619CFF";


    Excu(qq($Rscript $src/volcano.r $opts{infile} $opts{outpfx} $opts{xcol} $opts{ycol} $opts{xcut} $opts{ycut} "$opts{title}" $opts{logformat} "$opts{color}"));

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
    if(CheckFile($opts{infile},"line"=>1) ==  0){
        print STDERR "Warnings: input file [$opts{infile}] is empty or only head line\n";
        exit(0);
    }

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
    my $name_list = exists $opts{name} ? "-name \"$opts{name}\"" : "";

    ParaError2(\%opts,"list,outdir",'header','glist','name',"updown");

    Excu("perl $perl/venn_diff.pl -l $opts{list} $name_list -h $opts{header} -outlist $opts{glist} -updown $opts{updown} -scale $opts{scale} -colors \"$opts{color}\" -outdir $opts{outdir}");

    return "$opts{outdir}/Venn.svg";
}

sub KeggPathwayClass{
    my %opts = @_;

    $opts{name_key} //= "Genes";
    
    die "\n Usage: perl $0 kegg_PathwayClass -name_key -infile -outpfx \n" if(ParaError(\%opts,'infile','outpfx'));
    return if(CheckExcu($opts{infile}=>'line1'));

    Excu("perl $apple/KEGGdraw.pl -xlab \"Number of $opts{name_key}\" --prefix $opts{outpfx} $opts{infile} ");
#    Excu("perl $perl/draw_PathwayClass.pl -k $opts{name_key} -f $opts{infile} -p $opts{outpfx}");

    return "$opts{outpfx}.svg";
}

sub Cor{
    my %opts = @_;

    $opts{cor_cut} //= "0";
    $opts{pvalue_cut} //= "0";
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

    my $pheatmap_file = "xxxx";
#    $pheatmap_file = "$opts{outpfx}.pearson.xls";
    if(exists $opts{target}){
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
        ##  outdir outpfx scale cluster_row cluster_col 
        if($opts{soft} eq 'ggplot2'){
            Excu("$Rscript $src/ggplot2_pheatmap.r $pheatmap_file $opts{outpfx}.pearson \"$opts{color}\" \"$opts{title}\"");
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
    
    $opts{color} //= "white,forestgreen";
    ParaError2(\%opts,"infile,outpfx","color",'cor_cut','pvalue_cut','method');
    
    Excu("$Rscript $src/ggcor.r $opts{infile} $opts{outpfx} \"$opts{color}\"");

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
    
    ParaError2(\%opts,"infile,outpfx",'scale','group',"log");
    return if(CheckExcu($opts{infile}=>"line1"));

    Excu("$Rscript $src/PCA.r $opts{infile} $opts{outpfx} $opts{scale} $opts{log}");
    my $file_list = PcaPoint(%opts,infile=>"$opts{outpfx}.PC_data.xls",col=>"1,2",title=>"PCA");

    return $file_list;
}

sub PcaPoint{
    my %opts = @_;
    
    ParaError2(\%opts,"infile,outpfx","color","hjust","col","title","group");

    $opts{color} //= "none";
    $opts{hjust} //= "-1.65";
    $opts{col} //= "1,2";
    $opts{title} //= "PCA";
    $opts{type} //= "PCA";
    $opts{group} //= "none";

    ##   file outpfx  col_list xlab ylab  title group_info  vjust  point_size text_size 
    Excu("$Rscript $src/point.r $opts{infile} $opts{outpfx}.PCA $opts{col} NA NA $opts{title} $opts{group} $opts{hjust} 0 0 \"$opts{color}\"") if($opts{type} =~ /PCA/);
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

    return if(CheckExcu($opts{xls}=>"line1"));

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
    ParaError2(\%opts,"infile,outpfx,name,title","slope");
    
    Excu("$Rscript $src/ReplotGSEA.r $opts{infile} $opts{name} \"$opts{title}\" $opts{outpfx}.pdf $opts{slope}");

    return "$opts{outpfx}.pdf";
}

sub TernaryPlot{
    my %opts = @_;

    ParaError2(\%opts,"infile,group,outdir","yaml");
    my $yaml = exists $opts{yaml} ? "-y $opts{yaml}" : ""; 
    my $Rscript = $apps{"Rscript.3.6"};
    Excu("$Rscript $src/ternary-plot.R -i $opts{infile} -m $opts{group} -o $opts{outdir} $yaml");

    return "$opts{outdir}/Ternary_plot.pdf";
}

sub BuscoFigure{
    my %opts = @_;
    ParaError2(\%opts,"i,outpfx");
    Excu("Rscript $src/busco_figure.R $opts{i} $opts{outpfx}");

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
    
    ParaError2(\%opts,"infile,outpfx","diff","pq[q]","type","top[20]");
    if($opts{new} == 1){
        return if(CheckExcu($opts{infile}=>'line1'));
        Excu("perl $perl/enrichPlot.pl -t $opts{top} -o $opts{outpfx}.svg -V $apps{sbv} -q $opts{pq} $opts{infile}");
        return "$opts{outpfx}.svg";
    }
    if($opts{type} eq 'GO'){
        Excu("$apps{python2} $python/merge.py $opts{infile} $opts{outpfx}.xls");
        $opts{infile} = "$opts{outpfx}.xls";
    }
    return if(CheckExcu($opts{infile}=>'line1'));
    my $diff_info = exists $opts{diff} ? "-d $opts{diff}" : "";
    Excu("perl $apps{enrichPlot} $diff_info -o $opts{outpfx}.svg $opts{infile}");

    if($opts{type} eq 'GO'){
        Clean("$opts{outpfx}.xls");
    }

    return "$opts{outpfx}.svg";
}

sub EnrichBubble{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx","pq_cut","pq","type");
    $opts{pq_cut} //= 0.05;
    $opts{pq} //= "Qvalue";
    $opts{type} //= "KO";

#    $opts{outpfx} = "$opts{outdir}/$opts{name}_bubble";
#    Excu("$apps{python2} $python/draw_Bubble.py $opts{infile} $opts{diff} $opts{outpfx} $opts{type} $opts{pq_cut} $opts{pq}"); 
#    Excu("$Rscript $src/enrich_bubble.r $opts{outpfx}.xls $opts{outpfx} $opts{pq_cut} $opts{pq}");
    Excu("$Rscript $src/enrich_bubble.r $opts{infile} $opts{outpfx} $opts{pq_cut} $opts{pq}");

    my $pdf = "$opts{outpfx}.pdf";
    if($opts{type} eq 'KO' || $opts{type} eq 'GO'){
        $pdf = "$opts{outpfx}.pdf,$opts{outpfx}_sp.pdf"
    }
    return $pdf;
}

sub Random{
    my %opts = @_;

    ParaError2(\%opts,"infile,outpfx");
    Excu("$Rscript $src/random.r $opts{infile} $opts{outpfx}");

    return "$opts{outpfx}.pdf";
}

sub PlsdaLoading{
    my %opts = @_;
        
    $opts{name_col} //= 1;
    ParaError2(\%opts,"infile,outpfx");
    Excu(qq($Rscript $src/plsda_loading_plot.r $opts{infile} $opts{outpfx} 2,3 "#F8766D" "NULL" none $opts{name_col}));

    return "$opts{outpfx}.pdf";
}

sub Plsda{
    my %opts = @_; 

    $opts{scale} //= "pareto";
    $opts{type} //= "OPLS-DA";
    ParaError2(\%opts,"infile,outpfx,group","scale[pareto]","type[OPLS-DA]");
    $Rscript = $apps{'Rscript.3.3.1'};

    Excu(qq($Rscript $src/pls-da.r $opts{infile} $opts{group} $opts{outpfx} $opts{scale} $opts{type}));

    my @pdfs = map{"$opts{outpfx}.$_"}qw/loading_plot.pdf permutation.pdf plsda.nosampleid.pdf plsda.pdf/;
    if($opts{type} eq 'OPLS-DA'){
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
    
    return if(CheckExcu($opts{infile}=>"col3"));
    
    Excu(qq($Rscript $src/cluster.r $opts{infile} $opts{outpfx}));

    return "$opts{outpfx}.pdf";
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
        
    Excu("perl $perl/get_exp.pl $opts{infile} $opts{outfile} -name_key $opts{name_key} -col $opts{col}");
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
        Excu(qq(perl $0 $type -infile plsda_loading_plot/plsda1.v1.9.loading_plot.xls -outpfx $outdir/plsda_loading_plot1.$version));
    }

    if($type eq 'enrichPlot' || $type eq 'all'){
        Excu(qq(perl $0 $type -infile bar_Gradient/reactome_bar_Gradient.xls  -pq q -outpfx $outdir/enrichPlot1.$version -new 1));
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
            if(CheckFile($file,line=>$1)==0){
                $err .= "file [$file] is empty\n";
            }
        }
        if($info{$file} =~ /col(\d+)/){
            if(CheckFile($file,col=>$1) == 0){
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
    my $file = shift;
    my %opts = @_;

    my $ok = 0;
    if(exists $opts{line} || exists $opts{col}){
        open FILE,$file or die $!;
        my $i = 0;
        while(<FILE>){
            $i++;
            if($. == 1 && exists $opts{col}){
                my @aa = split/\t/,$_;
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
