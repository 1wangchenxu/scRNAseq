#!/usr/bin/perl -w

=pod
description: draw KEGG Map, referenced from kegg_parser.pl
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090728
modified date: 20100824, 20100406, 20100205, 20091224, 20091204, 20091127, 20091124, 20091123, 20091119, 20091010, 20090918, 20090825, 20090806, 20090731
=cut
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use FindBin qw($RealBin $Bin);
use GD;
use Data::Dumper;
use lib "$RealBin/lib";
use version;
use Carp qw(carp confess croak);

use lib "$RealBin/../";
use Support_Program;

my $no_version_check = 1;
my $DEBUG = 0;
my ($ko,  $diff, $outdir, $help, $header, $nomap,$gene_type,$gene_type_map);
my $out_name;
my $diff_type = "nodiff";
my $color_list = "red,green";
my $version = '3.0';
my $coord = 0;
GetOptions("ko:s" => \$ko, "diff:s" => \$diff,
           "outdir:s" => \$outdir, "header!"=>\$header, 
           "outname:s" => \$out_name,"diff_type:s" => \$diff_type,
           "help|?" => \$help, "nomap!" => \$nomap,
           "no_version_check:s"=>\$no_version_check,
           "genetype=s"=>\$gene_type,"genetypemap=s"=>\$gene_type_map,
           "coord=s"=>\$coord, 
           "color=s"=>\$color_list,"DEBUG" => \$DEBUG,
        );

our %color;# = ('red'=>[255,0,0],'green'=>[0,255,0]);
my @colors = split/,/,$color_list;
my %updown2color = ("up"=>$colors[0],"down"=>$colors[1]);
        
$outdir ||= ".";


if (!defined $ko || defined $help) {
    &usage();
}
if(!-s $ko){
    print STDERR "kopath file is not exists!\n";
    &usage();
}

## soft 
my %apps;
#Support_Program("$RealBin/../",\%apps,1);
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");


## conf files 
my $kegg_version = GetVersion($ko);
#my %conf_file = VersionInfo($kegg_version,"$RealBin/../../Database/kegg");
my %conf_file = VersionInfo($kegg_version,"$apps{kegg_db}");
my $kegg_report_map_dir = "$apps{kegg_db}/report_kegg_map";

my $map_js_dir = "$outdir/map_js" if($coord == 1);
$outdir = "$outdir/$out_name\_map" if(defined $out_name);

if(-l $outdir){
    `rm $outdir`;
}elsif(-d $outdir){
    `rm -f $outdir/map*png`;
    `rm -f $outdir/map*html`;
#    `rm -rf $outdir/map_js`
    `rmdir $outdir`;
}else{
    system("mkdir -p $outdir") if($coord == 0);
}

# main
versions($ko,$version) if($no_version_check == 0);
my @overview_list = qw/01110 01100 01120 01130 01200 01210 01212 01230 01220/;
@overview_list = ();
my %overview_list = map{$_=>1} map{$_,"map$_","ko$_"} @overview_list;

# get diff ratio
my %ratios;
%ratios =  ReadDiff($diff) if(defined $diff && $diff_type ne 'nodiff');

## add gene_type infos
my %g_conf_info;
my (%gene_type,%gene_type_map);
my $gene_type_flag = 0;
if(defined $gene_type){
    %gene_type = ReadGeneType($gene_type);
    $gene_type_flag = 1;
    %gene_type_map = split/,/,$gene_type_map
}

# read kopath 
my %kopath = ReadKopath($ko);
#print Dumper(%kopath);

# read conf 
my %conf;
ReadConf($conf_file{map_conf},\%kopath,\%conf);
#print Dumper(%conf);

if($coord == 1){
    print "$map_js_dir\t$coord\n";
    mkdir $map_js_dir if(!-e $map_js_dir);
    print "$kegg_report_map_dir/map_v$kegg_version\n";
    WriteCoordsTempData(\%conf, $ko, $map_js_dir);
    Excu("ln -sf $kegg_report_map_dir/map_v$kegg_version $outdir");    
    exit(0);
}
mkdir $outdir if(!-e $outdir);

# draw map 
DrawMap(\%conf,$conf_file{map_dir},$outdir);

# write html
WriteHtml(\%conf,$conf_file{map_dir},$outdir);

CheckKopath(\%kopath,\%g_conf_info);

exit(0);

sub WriteCoordsTempData {

    my $conf_ref = shift;
    my $kopath_file = shift;
    my $outdir = shift;

    my $id = $kopath_file;
    $id =~ s/.*\///;
    $id =~ s/\.kopath$//;
    print "$outdir/coord\n";
    mkdir("$outdir/coord") if(!-d "$outdir/coord");

    open my $ofh, "> $outdir/coord/$id.coord" or die "$!:$outdir/coord/$id.coord";
    for my $path (sort keys %{$conf_ref}) {
#        print Dumper(%{$conf_ref->{$path}});
        for my $pos (sort keys %{$conf_ref->{$path}}) {
            print $ofh join("\t", $id, "ko$path", $conf_ref->{$path}->{$pos}->{idlist},
            $conf_ref->{$path}->{$pos}->{type},$conf_ref->{$path}->{$pos}->{word}), "\n";
        }
    }
    close $ofh;
    
#    mkdir("$outdir/js") if(!-d "$outdir/js");
#    my $coord_to_js = "$RealBin/generate_enrich_js.pl";
#    Excu("perl $coord_to_js $outdir/coord ko_coord $outdir/js");

    return 0;
}

sub WriteHtml{
    my $conf_ref = shift;
    my $mapdir = shift;
    my $outdir = shift;

    # generate html
    for my $path(keys %{$conf_ref}){
#        next if(exists $overview_list{$path});
        my $map = "map$path";
        open HTML,"$mapdir/$map.html" or die $!;
        open OUT,">","$outdir/$map.html" or die $!;
    
        my $js =<<JS;
    <script type="text/javascript">
        var area =document.getElementsByTagName("area");
        for(var i =0;i < area.length;i++){
            if(!area[i].getAttribute("onmouseover")){
                area[i].addEventListener('mouseover',function(event){
                    var text= this.getAttribute("title");
                    var e = event || window.event; 
                    var scrollX = document.documentElement.scrollLeft || document.body.scrollLeft;
                    var scrollY = document.documentElement.scrollTop || document.body.scrollTop;
                    var x = e.pageX || e.clientX + scrollX;
                    var y = e.pageY || e.clientY + scrollY;
                    \$('body').append('<div id="showDiv'+i+'" style="display:block;width:500px;border:soild 1px;background:#ccc;padding:20px 5px;position:absolute;z-index:1000;top:'+(y+1)+';left:'+(x+1)+';">'+text+'</div>');
                });

                area[i].addEventListener('mouseout',function(event){
                    \$("#showDiv"+i).hover(function(){},function(){\$(this).remove();});
                    var div = document.getElementById("showDiv"+i);
                    var x=event.clientX;
                    var y=event.clientY;
                    var divx1 = div.offsetLeft;
                    var divy1 = div.offsetTop;
                    var divx2 = div.offsetLeft + div.offsetWidth;
                    var divy2 = div.offsetTop + div.offsetHeight;
                    if( x < divx1 || x > divx2 || y < divy1 || y > divy2){
                        \$("#showDiv"+i).remove();
                    }
                });
            }
        }
    </script>
JS
        local($_);
        while(<HTML>){
            if(/^\s*<area/){
                my @aa = split/\t/,$_;
                my ($pos) = $aa[1] =~ /coords=(\S+)$/;
                ($pos) = $_ =~ /coords=\"(\S+)\"/ if(!defined $pos);
                if($_ !~ /title/){
                    print OUT "$_";
                    next;
                }
                my ($former,$knums,$later) = $_  =~ /(^.*?title=\")(.*?)(\"[\d\D]*)/;
                if(exists $conf_ref->{$path}->{$pos}){
#                    $former =~ s/title=/data=/;
                    print OUT"$former$knums\n$conf_ref->{$path}->{$pos}->{word}\n$later";
                }
                else{
                    print OUT "$_";
                } 
            }elsif(/<\/html>/){
                print OUT "$js\n$_";
            }elsif(/dhtml.js/){
                print OUT "$_<script language=\"javascript\" src=\"http://code.jquery.com/jquery-1.8.3.min.js\"></script>\n";
            }else{
                print OUT"$_";
            }
        }
        close HTML;
    }

    return;
}

sub DrawMap{
    my $conf_ref = shift;
    my $mapdir = shift;
    my $outdir = shift;

#    my %color = ('red'=>[255,0,0],'green'=>[0,255,0]);

    for my $path(sort keys %{$conf_ref}){
#        next if($path ne '00030');
        my $map = "map$path.png";
#        print "$map\n";
        if(exists $overview_list{$path}){
#            print "cp $mapDir/$map $outdir/$map\n";
            `cp $mapdir/$map $outdir/$map`; 
            next;
        }
        open PNG, "< $mapdir/$map" or die $!;
#        print STDERR "$mapdir/$map\n";
        my $im = GD::Image->new(*PNG);
        open RES, "> $outdir/$map" or die $!;
        binmode RES;
        print RES $im->png;
        close PNG;
        close RES;

        open PNG, "< $outdir/$map" or die $!;
        $im = GD::Image->new(*PNG);
        my %color_pan = ("red"=>[(255,0,0)],"green"=>[(0,255,0)],"blue"=>[(0,0,255)],"pink"=>[255,192,203]);
        my @colors = split/,/,$color_list; ## red,green
        for my $ud(qw/up down/){
            $color{$ud} = $im->colorAllocate(@{$color_pan{$updown2color{$ud}}});
        }
#        print "$colors[0]\n";
#        print "@{$color_pan{$colors[0]}}\n";
#        $color{up} = $im->colorAllocate(@{$color_pan{$colors[0]}});
#        $color{down} = $im->colorAllocate(@{$color_pan{$colors[1]}});
#        DrawRect("100,200,200,300","up",$im);        
        for my $pos(sort keys %{$conf_ref->{$path}}){
            my %pos = %{$conf_ref->{$path}->{$pos}};
#            print Dumper(%pos);
            next if(!exists $pos{shape});
            Check($path,$pos{shape},$pos{type},$pos);
            if($pos{shape} eq 'rect'){
                DrawRect($pos,$pos{type},$im);
            }elsif($pos{shape} eq 'line'){
                DrawLine($pos,$pos{type},$im);
            }elsif($pos{shape} eq 'circle'){
                DrawCircle($pos,$pos{type},$im);
            }elsif($pos{shape} eq 'poly'){
                my @pos = split/,/,$pos;
                @pos = RePos(@pos) if($path =~ /^011/);
                my $i = 0;
#                print "=>$path\t$conf_ref->{$path}->{$pos}->{word}\t$pos\n";
                while($i < @pos-4){
                    $pos = join(",",@pos[$i..$i+7]);
                    $pos = Poly2Line($pos);
                    DrawLine($pos,$pos{type},$im);
                    $i+=4;
                }
            }
#            print Dumper(%pos);die;
        }
        close PNG;

        open RES, "> $outdir/$map" or die $!;
        binmode RES;
        print RES $im->png;
        close RES;
#        die;
    }
}

sub RePos{
    my @pos = @_;

    my @out;
    my $ratio = 0.35;
    $ratio = 0.3 if($kegg_version > 4.0);

    for my $pos(@pos){
        push @out,$pos/$ratio;
    }

    return @out;
}


sub DrawRect{
    my $pos = shift;
    my $type = shift;
    my $im = shift;

#    my %color = ('red'=>[255,0,0],'green'=>[0,255,0]);   
#    my $red = $im->colorAllocate(255, 0, 0);    
#    $im->rectangle(100,100,500,500,$red);
    my @pos = split/,/,$pos;
    my ($x1,$y1,$x2,$y2) = @pos;
    my @pos2 = ($pos[0]+1,$pos[1]+1,$pos[2]-1,$pos[3]-1);

    #Check("Rect",$type,@pos);
    if($type =~ /up/){ # red
        $im->rectangle(@pos, $color{up});
        $im->rectangle(@pos2, $color{up});
    }
    elsif($type =~ /down/){ # green
        $im->rectangle(@pos, $color{down});
        $im->rectangle(@pos2, $color{down});
    }
    if($type =~ /up/ && $type =~ /down/){## half red, half green
        $im->line($x1,($y1+$y2)/2,$x1,$y2,$color{down}); ## left
        $im->line($x1,$y2,$x2,$y2,$color{down});  ## bottom
        $im->line($x2,($y1+$y2)/2,$x2,$y2,$color{down}); ## right
        
        ## 2 line
        $im->line($x1+1,($y1+$y2)/2,$x1+1,$y2,$color{down}); ## left
        $im->line($x1,$y2-1,$x2,$y2-1,$color{down});  ## bottom
        $im->line($x2-1,($y1+$y2)/2,$x2-1,$y2,$color{down}); ## right
    }
    
    return;
}

sub Poly2Line{
    my $pos = shift;

    my @pos = split/,/,$pos;

    my ($x1,$y1,$x2,$y2,$x3,$y3,$x4,$y4) = split/,/,$pos;  

    my $sx = ($x1+$x4)/2;
    my $sy = ($y1+$y4)/2;
    my $ex = ($x2+$x3)/2;
    my $ey = ($y2+$y3)/2;

    return "$sx,$sy,$ex,$ey";
}


sub DrawLine{
    my $pos = shift;
    my $type = shift;
    my $im = shift;

#    my %color = ('red'=>[255,0,0],'green'=>[0,255,0]);   

    my @pos = split/,/,$pos;
    #Check("Line",$type,@pos);
    my ($x1,$y1,$x2,$y2) = @pos;
    my ($x_move,$y_move) = (0,0);
    if($y1 == $y2){
        $y_move = 1;
    }else{
        $x_move = 1;
    }
    my @pos_1 = ($x1-$x_move,$y1-$y_move,$x2-$x_move,$y2-$y_move);
    my @pos1 = ($x1+$x_move,$y1+$y_move,$x2+$x_move,$y2+$y_move);
    my @pos2 = ($x1+$x_move*2,$y1+$y_move*2,$x2+$x_move*2,$y2+$y_move*2);
    if($type =~ /up/){ # red
        $im->line(@pos,$color{up});
        $im->line(@pos_1,$color{up});
        $im->line(@pos1,$color{up});
        $im->line(@pos2,$color{up});
    }elsif($type =~ /down/){ # green
        $im->line(@pos,$color{down});
        $im->line(@pos_1,$color{down});
        $im->line(@pos1,$color{down});
        $im->line(@pos2,$color{down});
    }
    if($type =~ /up/ && $type =~ /down/){ ## half red half green
        $im->line(@pos1,$color{down});
        $im->line(@pos2,$color{down});
    }

    return;
}

sub DrawCircle{
    my $pos = shift;
    my $type = shift;
    my $im = shift;

#    my %color = ('red'=>[255,0,0],'green'=>[0,255,0]);   

    my @pos = split/,/,$pos;
    #Check("Circle",$type,@pos);
    $pos[-1] *= 2.5;
    push @pos,$pos[-1];
    if($type =~ /up/){ # red
         $im->filledArc(@pos,0,360,$color{up});
    }
    elsif($type =~ /down/){ # green
         $im->filledArc(@pos,0,360,$color{down});
    }
    if($type =~ /up/ && $type =~ /down/){ # half red  half green
        $im->filledArc(@pos,0,180,$color{down});
    }

    return;
}

sub ReadConf{
    my $conf_file  = shift;
    my $kopath_ref = shift;
    my $conf_ref = shift;

    local($_);
    open CONF,$conf_file or die $!;
    while(<CONF>){
        chomp;
        my ($path,$shape,$pos,$idlist) = split/\t/,$_;
        my @ids = split/\Q+\E/,$idlist;
        for my $id(@ids){
            $g_conf_info{$id}{$path} = 1;
        }

        next if(!exists $kopath_ref->{$path});
        my ($type,$word) = GetType($idlist,$kopath_ref->{$path});
        next if($type eq '');
        $conf_ref->{$path}->{$pos}->{type} = $type;
        $conf_ref->{$path}->{$pos}->{word} = $word;
        $conf_ref->{$path}->{$pos}->{shape} = $shape;
        $conf_ref->{$path}->{$pos}->{idlist} = $idlist;
    }
    close CONF;
#    print Dumper($conf_ref);

    return;
}

sub GetType{
    my $idlist = shift;
    my $path_ref = shift;

    my @ids = split/\Q+\E/,$idlist;

    my $type = '';
    my $word = '';

    my $ratio_flag = 0;
    if(keys %ratios > 0){
        $ratio_flag = 1;
    }
    my %hash;
    my %hash2;
    
    for my $id(@ids){
        next if(!exists $path_ref->{$id});
        for my $gene(@{$path_ref->{$id}}){
            if($ratio_flag){
                next if(!exists $ratios{$gene});
                $hash{$gene} = $ratios{$gene};
                if($ratios{$gene}>0){
                    $hash2{up}++;
                }
                else{
                    $hash2{down}++;
                }
            }elsif($gene_type_flag){
                next if(!exists $gene_type{$gene});
                for my $g_info(@{$gene_type{$gene}}){
                    my ($type,$log2fc) = split/,/,$g_info; ## type mRNA or protein ..
                    Error("No gene_type_map info of [$type]") if(!exists $gene_type_map{$type});
                    my $ud_type = $gene_type_map{$type}; ## change to up or down
                    $hash2{$ud_type}++;
                    my $gene_info = "$gene ($log2fc)";
                    $hash{$type}{$gene_info}++;
                }
            }else{
                $hash{$gene}++;
            }
        }
    }
    if($ratio_flag){
        my @keys = keys %hash2;
        
        my @genes = keys %hash;
        if(@keys > 1){ ## up down
            $type = "updown";
        }
        elsif(@keys == 1){
            $type = $keys[0];
        }
        if(@genes > 0){
            for my $gene(sort {$hash{$a} <=> $hash{$b} or $a cmp $b} @genes){
                $word .= " $gene ($ratios{$gene}),";
            }
            chop($word);
        }
    }elsif($gene_type_flag){
        my @keys = keys %hash2;
        if(@keys > 1){ ## up down
            $type = "updown";
        }elsif(@keys == 1){
            $type = $keys[0];
        }
        $word = "<br/>";
        for my $type(sort keys %hash){ 
            my $ud_type = $gene_type_map{$type};
            Error("No color of [$ud_type]") if(!exists $updown2color{$ud_type});
            my $color = $updown2color{$ud_type};
            $word .= "<font color='$color'>$type:".join(",",sort keys %{$hash{$type}})."</font><br/>";
        }
    }else{
        my @keys = keys %hash;
        if(@keys > 0){
            $word = join(",",sort keys %hash);
            $type = "up";
        }
    }

    return ($type,$word);
}


sub ReadDiff{
    my $diff = shift;
    open DIFF, "< $diff" or die "$!:$diff";
    <DIFF> if (defined $header);
    my %ratios;
    while (<DIFF>) {
        chomp;
        my @temp = split /\s+/, $_;
        next if(@temp == 1);
        $ratios{$temp[0]} = sprintf("%.1f", $temp[1]);
    }
    close DIFF;
    
    return %ratios;
}

sub ReadGeneType{
    my $file = shift;

    local($_);
    open FILE,"$file" or die $!;
    my %hash;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        $aa[1] = sprintf("%.1f", $aa[1]);
        push @{$hash{$aa[0]}},"$aa[2],$aa[1]";
#        $hash{"$aa[0] ($aa[1])"} = $aa[2];
    }
    close FILE;

    return %hash;
}

sub ReadKopath{
    my $ko = shift;
    
    open KO, "< $ko" or die $!;
    my %kopath;
    local($_);
    while (<KO>) {
        chomp;
        next if ($_ eq "" or index($_, "#") == 0);
        my @temp = split /\t/, $_;
        next if (@temp < 3 || $temp[1] =~ /^\s*$/);
        my @kos = split/,/,$temp[2];

        for my $ko(@kos){
            push @{$kopath{$ko}{$temp[1]}},$temp[0];
            push @{$kopath{$ko}{$temp[3]}},$temp[0] if(@temp == 4);
        }
    }
    close KO;

    return %kopath;
}

sub CheckKopath{
    my $kopath_ref = shift;
    my $conf_ref = shift;

    my $err;
    for my $path(sort keys %{$kopath_ref}){
        next if($path eq '-');
        for my $k(sort keys %{$kopath_ref->{$path}}){
            if($k =~ /^[KC]\d+$/){
                if(!exists $conf_ref->{$k}{$path}){
                    $err .= "No $k in $path in conf\n";
#                    print STDERR "$conf_file{map_conf}\n";
                }
            }
        }
    }
    print STDERR "$err$conf_file{map_conf}\n" if(defined $err);
    return;
}

sub Check{
    my @infos = @_;

    return if($DEBUG==0);
    print join("\t",@infos)."\n";

    return;
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    system($cmd);
    return ;
}

sub Error{

    confess("Error @_\n");

}

sub usage {

    print STDERR << "USAGE";

    description: draw KEGG Map, referenced from kegg_parser.pl
    usage: perl $0 [options]
    options:
        -ko *: .kopath file
        -diff *: DiffGeneExpFilter.xls, result of solexa-mRNAtag_pipeline.pl
        -header : DiffGeneExpFilter.xls contains header
        -nomap : output nomap ko and KO and geneid
        -outdir: output directory, default is current directory "."
        -help|?: help information

        -genetype
        -genetypemap
        -color
    e.g.:
        perl $0 -ko ko -diff DiffGeneExpFiltr.xls -outdir .

USAGE

    exit 1;
}
