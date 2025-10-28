#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use threads;
use FindBin qw($Bin $RealBin);
use Cwd 'abs_path';
use Carp qw(carp confess croak);
use lib "$RealBin/..";
use Support_Program;

my $soft_conf = "$RealBin/../Small_pipe.soft.conf";
my $version = "v1.1";

my ($indir,$glist_file,$kopath,$bgl,$go,$help,$format,$debug,$outdir,$symbol,$ref);
my $type = "nodiff";
my $key = "gene";
my $nomap_die = 1;
my %key_word = ('gene'=>{'s'=>'Genes','a'=>'gene','u'=>'Gene'},
                'protein'=>{'s'=>'Proteins','a'=>'protein','u'=>'Protein'},
                'metabolite' => {'s' => 'Metabolites', 'a' => 'metabolite', 'u' =>'Metabolite'},
                'target_genes'=>{'s'=>'Target_Genes','a'=>'target_gene','u'=>'Target_Gene'},
                'source_gene' =>{'s'=>'Source_Genes','a'=>'source_gene','u'=>'Source_Gene'},);
$key_word{target_gene} = $key_word{target_genes};
$key_word{mirnas} = $key_word{target_genes};
$key_word{circrnas} = $key_word{source_gene};

my $run = 1;
my $thread = 3;
my $no_thread = 1;
my $run_enrichment = 1;
my @updown_types = qw/All Up Down/;
my ($pq,$top) = ("q",20); ## for draw
my ($do,$reactome);
my ($run_ko,$run_go,$run_do,$run_reactome) = qw/1 1 1 1/;
my ($run_MF,$run_BP,$run_CC) = qw/1 1 1/;
my ($run_bar_Gradient) = (1);
my $summary = 0;
my $err = 0;
my $draw_type = 'all';
my $omic_report = 0;
GetOptions(
    "i=s"=>\$indir,
    "glist=s" => \$glist_file,
    "ref=s" => \$ref,
    "annot=s"=>\$go,
    "ko=s"=>\$kopath,
    "do=s" => \$do,
    "reactome=s" => \$reactome,
    "bgl=s"=>\$bgl,
    "symbol=s" => \$symbol,
    "type=s"=>\$type,
    "o=s" =>\$outdir,
    "run=s" => \$run,
    "h"=>\$help,
    "format"=>\$format,
    "thread=s" => \$thread,
    "enrichment=s" => \$run_enrichment,
    "keyword=s" => \$key,
    "run_ko=s" => \$run_ko,
    "run_do=s" => \$run_do,
    "run_go=s" => \$run_go,
    "run_reactome=s" => \$run_reactome,
    "run_BP=s" => \$run_BP,
    "run_MF=s" => \$run_MF,
    "run_CC=s" => \$run_CC,
    "pq=s"     => \$pq,
    "top=s"    => \$top,
    "summary=s" => \$summary,
    "run_bar_Gradient=s" => \$run_bar_Gradient,
    "draw_type=s" => \$draw_type,
    "nomap_die=s" => \$nomap_die,
    "omic_report=s" => \$omic_report,
    "d"=>\$debug,
    );

if (defined $help or (!defined $indir && !defined $glist_file)){# or ((!defined $go or !defined $bgl) && !defined $kopath && !defined $do && !defined $reactome)){
    &usage();
}

InitOpts(run=>$run,'multi_type'=>'threads',cpu=>$thread);

my %apps;
#Support_Program_conf("$RealBin",\%apps,$soft_conf);
#$Bin = "/Bio/Project/PROJECT/kegg_example/new_enrich";
## pipe start
my $shell = "$RealBin/src";
$indir = abs_path($indir) if(defined $indir);
$kopath = abs_path($kopath) if(defined $kopath);
$bgl = abs_path($bgl) if(defined $bgl);
#$go = abs_path($go) if(defined $go);

$key = lc($key);
if(exists $key_word{$key}){
    %key_word = %{$key_word{$key}};
#    print %key_word;
}
else{
    Info("key [$key] is unrecognizable!");die;
}

$outdir //= $indir;
$outdir = abs_path($outdir);

mkdir $outdir if(!-e $outdir);

## software
my $src = "$RealBin/src";
my $R_plot = "perl $src/R_plot.pl";
my $keggpath = "perl $src/keggpath.pl";
#my $keggGradient = "perl $src/keggGradient.pl";
my $keggGradient = "perl $src/enrich_bar_Gradient.pl -key_word $key_word{u} -top $top -pq $pq -fig $run_bar_Gradient -draw_type $draw_type";
my $draw_PathwayClass = "perl  $src/draw_PathwayClass.pl";
my $keggMap = "perl $src/keggMap.pl";
my $keggMap_nodiff = "perl $src/keggMap_nodiff.pl";
my $genPathHTML = "perl $src/genPathHTML.pl";
my $enrichmentHeatmap = "perl $src/enrichmentHeatmap.pl";
my $enrichGO = "perl $src/enrichGO.pl";
my $doenrich = "perl $src/do_enrich.pl";
my $reactome_enrich = "perl $src/reactome_enrich.pl";
my $summary_perl = "perl $src/enrich_summary.pl";
my $generate_enrich_js = "perl $src/generate_enrich_js.pl";

## ref 
AllRef($ref) if(defined $ref);
my $g_symbol_para = defined $symbol ? "-symbol $symbol" : "";

## prepare glist 
my %glist_list;
if(defined $glist_file){
    GetGlist('none',$type,$outdir,\%glist_list,'-glist_file'=>$glist_file);
    $run_enrichment = 0;
}else{
    GetGlist($indir,$type,$outdir,\%glist_list);
}

## kegg enrich 
if(defined $kopath && -s $kopath){
     RunKeggEnrich(\%glist_list,$type,"$outdir");
}

## go enrich 
if(defined $go){
    RunGoEnrich(\%glist_list,$type,"$outdir");
}

#print "$do\n$reactome\n";
if(defined $do && -s $do){
    RunDoEnrich(\%glist_list,$type,"$outdir");
}

if(defined $reactome && -s $reactome){
    RunReactomeEnrich(\%glist_list,$type,"$outdir");
}


if($err > 0){
    Info("ERROR: All [$err] errors, please check");
    exit(1);
}else{
    Info("All Done");
    exit(0);
}

sub AllRef{
    my $ref = shift;

    $kopath = "$ref.kopath";
    if(-s "$ref.annot"){
        $go = "$ref.annot";
    }elsif(-s "$ref.bgl"){
        $go = "$ref";
        $bgl = "$ref.bgl";
    }else{
        $run_go = 0;
    }
    $do = "$ref.do";
    $reactome = "$ref.reactome";
    $symbol = "$ref.gen2sym" if(-e "$ref.gen2sym");

    if(!-e $kopath && $run_go == 0 && !-e $reactome && !-e $do){
        Error("No any ref info of [$ref] found\n");
    }

    return;
}

sub RunKeggEnrich{
    my $file_ref = shift;
    my $type = shift;
    my $outdir = shift;
 
    if($run_ko == 0){
        return;
    }
    MakeEnrichDir($file_ref,$outdir,"KO",$type);
    
    my @cmds;
    ## kegg path 
    Info("run kegg path ");
    @cmds = map{"$keggpath PATH -f $_ -b $kopath -k $key_word{s} $g_symbol_para -o $file_ref->{$_}->{kegg}"} keys %{$file_ref};
    Excu2(\@cmds);
    
    ## keggGradient
    Info("run keggGradient");
#    @cmds = map{"$keggGradient $file_ref->{$_}->{kegg}.path.xls 20 Q $key_word{u}"} keys %{$file_ref};
    @cmds = map{"$keggGradient -glist $_ $file_ref->{$_}->{kegg}.path.xls $file_ref->{$_}->{kegg} KEGG"} keys %{$file_ref};
    Excu2(\@cmds);
    
    ## PathwayClass
    Info("run PathwayClass");
#    @cmds = map{"$draw_PathwayClass  -f $file_ref->{$_}->{kegg}.path.xls -k $key_word{s} -p $file_ref->{$_}->{kegg}.path.xls"}  keys %{$file_ref};
    @cmds = map{"$R_plot kegg_PathwayClass -infile $file_ref->{$_}->{kegg}.path.xls -name_key $key_word{s} -outpfx $file_ref->{$_}->{kegg}.path.xls"}  keys %{$file_ref};
    Excu2(\@cmds);

    ## kegg bubble
#    if($type eq 'diff'){
#        Info("run kegg bubble");
#        @cmds = map{"$R_plot enrich_bubble -infile $file_ref->{$_}->{kegg}.path.xls -diff $_ -type KO -outpfx $file_ref->{$_}->{kegg}.bubble"}keys %{$file_ref};
#    }
#    Excu($thread,@cmds);
    
    ## keggMap
    Info("run keggMap");
    my $keggMap_cmd = "$keggMap -diff_type $type -coord $omic_report ";

    @cmds = map{"$keggMap_cmd -outdir $file_ref->{$_}->{kegg_dir} -ko $file_ref->{$_}->{kegg}.kopath -diff $_ -outname $file_ref->{$_}->{filename}"}  keys %{$file_ref};

#    if($type eq 'nodiff'){
#        @cmds = map{"$keggMap -ko $file_ref->{$_}->{kegg}.kopath -outdir $file_ref->{$_}->{kegg}_map"} keys %{$file_ref};
#    }else{
#        @cmds = map{"$keggMap -ko $file_ref->{$_}->{kegg}.kopath -diff $_ -outdir $file_ref->{$_}->{kegg}_map"}  keys %{$file_ref};
#    }
    Excu2(\@cmds);

    ## js 
    if($omic_report == 1){
        Info("generate_enrich_js");
        @cmds = ("$generate_enrich_js $outdir/KO/map_js/coord ko_coord $outdir/KO/map_js");
        Excu2(\@cmds);
    }

    ## genPathHTML
    Info("run genPathHTML");
    @cmds = ("$genPathHTML -k $key_word{a} -indir $outdir/KO -diff_type $type -nomap_die $nomap_die");
    if($type eq 'updown'){
        @cmds = map{"$genPathHTML -indir $outdir/$_/KO"}@updown_types;
    }
    Excu2(\@cmds);

    ## enrichmentHeatmap
    Enrichment($outdir,"KO",$type);

    ## summary
    Info("run kegg summary");
    if($summary == 1){
        @cmds = ("$summary_perl -indir $outdir/KO -name_pfx .path.xls -output $outdir/KO/all_pathway.xls");
        if($type eq 'updown'){
            @cmds = map{"$summary_perl -indir $outdir/$_/KO -name_pfx .path.xls -output $outdir/$_/KO/all_pathway.xls"} @updown_types;
        }
        Excu2(\@cmds);
    }
#    Excu(1,"$summary_perl -indir $outdir/$_/KO -name_pfx .path.xls -output $outdir/KO/all_pathway.xls") if($summary == 1);

    return;
}

sub RunGoEnrich{
    my $file_ref = shift;
    my $type = shift;
    my $outdir = shift;

    if($run_go == 0){
        return;
    }

    MakeEnrichDir($file_ref,$outdir,"GO",$type);
    
    my @cmds;
    ## go enrich
    Info("run enrichGO");
#    if(defined $go && -e $go){  ## annot file 
    if(defined $go && -e $go && $go =~ /annot$/){  ## annot file
        @cmds = map{"$enrichGO -g $_ -annot $go -op $file_ref->{$_}->{go} -k $key_word{s} $g_symbol_para -ud $file_ref->{$_}->{ud} -run_BP $run_BP -run_MF $run_MF -run_CC $run_CC"} keys %{$file_ref};
    }elsif(defined $bgl && defined $go){
        @cmds = map{"$enrichGO -g $_ -bg $bgl -a $go -op $file_ref->{$_}->{go} -k $key_word{s} $g_symbol_para -ud $file_ref->{$_}->{ud}"} keys %{$file_ref};
    }else{
        Info("go is not defined or not exist!");
        return;
    }
    Excu($thread,@cmds);

    ## go bubble
    Info("run go gradient bubble pqbar");
    @cmds = map{"$keggGradient -glist $_ $file_ref->{$_}->{go} $file_ref->{$_}->{go} GO"} keys %{$file_ref};
    Excu($thread,@cmds);

    ## enrichment 
    Enrichment($outdir,"GO",$type);

    ## summary
    Excu(1,"$summary_perl -indir $outdir/GO -name_pfx .Level2.xls -output $outdir/GO/all_level2.xls") if($summary == 1);

    return;
}

sub RunDoEnrich{
    my $file_ref = shift;
    my $type = shift;
    my $outdir = shift;

    if($run_do == 0){
        return;
    }

    MakeEnrichDir($file_ref,$outdir,"DO",$type);
    
    my @cmds;
    ## Do enrich 
    Info("run DO enrich");
    @cmds = map{"$doenrich ENRICH -i $_ -bg $do -k $key_word{s}  $g_symbol_para -outpfx $file_ref->{$_}->{do}"} keys %{$file_ref};
    Excu($thread,@cmds);

    ## bar gradiant 
    @cmds = map{"$keggGradient -glist $_ $file_ref->{$_}->{do}.do.xls $file_ref->{$_}->{do} DO"} keys %{$file_ref};
    Excu($thread,@cmds);
    
    ## do bubble
#    if($type eq 'diff'){
#        Info("run do bubble");
#        @cmds = map{"$R_plot enrich_bubble -infile $file_ref->{$_}->{do}.do.xls -diff $_ -type DO -outpfx $file_ref->{$_}->{do}.bubble"}keys %{$file_ref};
#    }
#    Excu($thread,@cmds);
   
    ## html 
    Info("run do HTML");
    @cmds = ("perl $src/DO_HTML.pl -k $key_word{a} -indir $outdir/DO ");
    if($type eq 'updown'){
        @cmds = map{"perl $src/DO_HTML.pl -k $key_word{a} -indir $outdir/$_/DO"}@updown_types;
    }
    Excu($thread,@cmds);


    ## enrichment 
#    Enrichment($outdir,"DO",$type);
    
    ## summary
    Excu(1,"$summary_perl -indir $outdir/ -name_pfx .do.xls -output $outdir/DO/all_do.xls") if($summary == 1);

    return;
}

sub RunReactomeEnrich{
    my $file_ref = shift;
    my $type = shift;
    my $outdir = shift;

    if($run_reactome == 0){
        return;
    }

    MakeEnrichDir($file_ref,$outdir,"Reactome",$type);
    
    my @cmds;
    ## enrich 
    @cmds = map{"$reactome_enrich ENRICH -i $_ -bg $reactome -k $key_word{s} $g_symbol_para -outpfx $file_ref->{$_}->{Reactome}"} keys %{$file_ref};
    Excu($thread,@cmds);

    ## bar gradiant 
    @cmds = map{"$keggGradient -glist $_ $file_ref->{$_}->{Reactome}.reactome.xls $file_ref->{$_}->{Reactome} reactome"} keys %{$file_ref};
    Excu($thread,@cmds);
    
    ## reactome bubble
#    if($type eq 'diff'){
#        Info("run Reactome bubble");
#        @cmds = map{"$R_plot enrich_bubble -infile $file_ref->{$_}->{Reactome}.reactome.xls -diff $_ -type Reactome -outpfx $file_ref->{$_}->{Reactome}.bubble"}keys %{$file_ref};
#    }
#    Excu($thread,@cmds);
    
    ## html 
    Info("run Reactome HTML");
    @cmds = ("perl $src/Reactome_HTML.pl -k $key_word{a} -indir $outdir/Reactome");
    if($type eq 'updown'){
        @cmds = map{"perl $src/Reactome_HTML.pl -k $key_word{a} -indir $outdir/$_/Reactome"}@updown_types;
    }
    Excu($thread,@cmds);
   
    
    
    ## enrichment 
#    Enrichment($outdir,"Reactome",$type);

    ## summary
    Excu(1,"$summary_perl -indir $outdir/ -name_pfx .reactome.xls -output $outdir/Reactome/all_reactome.xls") if($summary == 1);

    return;
}

sub MakeEnrichDir{
    my $file_ref = shift;
    my $outdir = shift;
    my $enrich_type = shift; # KO GO DO Reactome 
    my $ud_type = shift;

    my @dirs = ("$outdir/$enrich_type");
    if($ud_type eq 'updown'){
        @dirs = map{"$outdir/$_/$enrich_type"} @updown_types;
    }
    for my $dir(@dirs){
        MakeDir($dir);
        map{MakeDir("$dir/$file_ref->{$_}->{filename}")} keys %{$file_ref} if($format);
    }

    return;
}

sub Enrichment{
    my $outdir = shift;
    my $enrich_type = shift;
    my $ud_type = shift;


    return if(!$run_enrichment);
    Info("run $enrich_type enrichmentHeatmap");
    my @cmds = ("$enrichmentHeatmap $outdir/$enrich_type $enrich_type");
    if($ud_type eq 'updown'){
        @cmds = map{"$enrichmentHeatmap $outdir/$_/$enrich_type $enrich_type"} @updown_types;
    }
    Excu($thread,@cmds);

    return;
}

sub GetGlist{
    my $indir = shift;
    my $type = shift;
    my $outdir = shift;
    my $file_ref = shift;
    my %opts = @_;

    my @files;
    my $file_num = 0;
    if(exists $opts{'-glist_file'}){
        @files = split/,/,$opts{'-glist_file'};
        $file_num = @files;
    }else{
        @files = glob("$indir/*.glist");
        $file_num = @files;
        Info("Get glist from dir: $indir, found [$file_num] glist");
    }

    my @cmds;
    if($file_num > 0){
        if($type eq 'updown'){
            for my $i(@updown_types){
                mkdir "$outdir/$i" if(!-e "$outdir/$i");
            }
        }
        for my $file(@files){
            my $file_name = $file;
            $file_name =~ s/.*\///;
            $file_name =~ s/.glist$//;

            my @types = ('');
            if($type eq 'updown'){
                @types = @updown_types;
            }
            for my $i(@types){
                my $t_file = "$outdir/$i/$file_name.glist";
                if($i eq 'All' || $i eq ''){
                    $t_file = $file;
                }
                $file_ref->{$t_file}->{filename} = $file_name;
                $file_ref->{$t_file}->{kegg} = "$outdir/$i/KO/$file_name";
                $file_ref->{$t_file}->{kegg} = "$outdir/$i/KO/$file_name/$file_name" if($format);
                $file_ref->{$t_file}->{kegg_dir} = "$outdir/$i/KO";
                $file_ref->{$t_file}->{kegg_dir} = "$outdir/$i/KO/$file_name" if($format);
                $file_ref->{$t_file}->{go} = "$outdir/$i/GO/$file_name";
                $file_ref->{$t_file}->{go} = "$outdir/$i/GO/$file_name/$file_name" if($format);
                $file_ref->{$t_file}->{do} = "$outdir/$i/DO/$file_name";
                $file_ref->{$t_file}->{do} = "$outdir/$i/DO/$file_name/$file_name" if($format);
                $file_ref->{$t_file}->{Reactome} = "$outdir/$i/Reactome/$file_name";
                $file_ref->{$t_file}->{Reactome} = "$outdir/$i/Reactome/$file_name/$file_name" if($format);
                if($type eq 'diff' || $i eq 'All'){
                    $file_ref->{$t_file}->{ud} = "diff";
                }
                else{
                    $file_ref->{$t_file}->{ud} = "nodiff";
                }
                if($i eq $updown_types[1]){ ## up
                    push @cmds,"awk -F\$'\\t' -vOFS='\\t' '\$2>0{print \$1,\$2}' $file >$t_file";
                }elsif($i eq $updown_types[2]){ ## down
                    push @cmds,"awk -F\$'\\t' -vOFS='\\t' '\$2<0{print \$1,\$2}' $file >$t_file"
                }
            }
        }

    }else{
        Info("No glist input, exit now!");
        die;
    }
    Excu($thread,@cmds) if(@cmds > 0);

    for my $glist(sort keys %{$file_ref}){
        if(!-s $glist){
            Info("$file_ref->{$glist}->{filename} glist file [$glist] is empty!delete it!");
            delete $file_ref->{$glist};
        }
    }
    if(keys %$file_ref == 0){
        Info("Warings: All glists are empty");
        exit(0);
    }
#    print keys %$file_ref;
#    die;
    
    return;
}


sub Excu{
    my $thread  = shift;
    my @cmds = @_;

    Excu2(\@cmds);

    return;
}

sub Error{

    confess("Error:@_\n");
    $err++;
    return;
}

sub Info{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
    $year += 1900;
    $mon += 1;

    my $ftime = sprintf("%d\-%02d\-%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
    print STDERR "[$ftime]: @_\n";

    return ;
}

sub usage{

my $usage = <<"USAGE";

    version: $version

    Usage:  perl $0 -i <indir> -ko <all.kopath> -bgl <all.bgl> -annot <dir/all>  -type <str>
    Options:
        -i      dir       dir contains genelist(format: *.glist)
        -glist  file_list glist file list, sep by ","
        
        -ref    refpfx    xxx.kopath, xxx.annot, xxx.gen2sym, xxx.do xxx.rectome

        or 

        -ko     file      ko file(kopath)
        -bgl    file      go bgl file
        -type   strings   diff enrich or nodiff enrich [nodiff], [updown]
        -symbol file      gene symbol file 
        new go ref : 
            -annot  file      xxx.annot file 
        
        old go ref : 
            -annot  file      fullpath's dir and pfx of go file(dir/*.[CPF])
            -bgl    file      go bgl file 
        -do     file      do file 
        -reactome file    rectome file 

        -o      strings   output dir, default: indir
        -thread int       thread for run [3]
        -enrichment int   run enrichment or not [1]/0 

        -key    strings   [gene] protein,metabolite,mRNA target_genes

        -format           format dir(one sample one dir)
                            default:  indir/KO/name
                            format :  indir/KO/name/name 
        -run_BP,run_MF,run_CC     GO PFC
        -run              [1]/0
                           run_do, run_go, run_ko, run_reactome
        -nomap_die        [1]/0    die when no map found

        -h                help
    eg:
    perl $0 -i indir -annot all -bgl all.bgl -ko all.kopath -type diff 

USAGE

    print "$usage";
    exit(0);

    return;
}

