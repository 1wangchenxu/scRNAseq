#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use Getopt::Std;
use FindBin qw($Bin $RealBin);
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename qw/basename dirname/;
use lib "$RealBin/lib";
use version;

use lib "$RealBin/../";
use Support_Program;

my $src = "$RealBin";

my $type = shift;

my %opts = ('k'=>'Genes','t'=>'geneid','o'=>'out',v=>"5.0");  #  for output title 
my %type_map = ('gene' => ['Genes','K_IDs'],'Metabolites'=>['Metabolites','C_IDs'],
                'Proteins'=>['Proteins','K_IDs'],'mRNA'=>['mRNAs','K_IDs'],
                'target' => ['target_genes','K_IDs'],'Source_Genes'=>['Genes','K_IDs']);
$type_map{protein}  = $type_map{Proteins};
$type_map{metabolite} = $type_map{Metabolites};
$type_map{Genes} = $type_map{gene};
$type_map{target_gene} = $type_map{target};
$type_map{Target_Genes} = $type_map{target};
#getopts('a:b:e:f:i:k:o:p:r:t:v:hn',\%opts);
GetOptions(\%opts,"a=s","b=s","e=s","f=s","i=s","k=s","o=s","p=s","r=s","t=s","v=s","h","n","symbol=s");
Error2("No type [$opts{k}] info") if(!exists $type_map{$opts{k}});

#$Bin = "/Bio/Database/Database/kegg/latest_kegg/shell";
## set soft path fishInWinter and Rscript
my %apps;
#Support_Program("$RealBin/../",\%apps,1);
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");

## select version 
my $version = $opts{v};
#my %conf_file = VersionInfo($opts{v},"$src/../../Database/kegg"); ## get db path 
my %conf_file = VersionInfo($version,"$apps{kegg_db}"); ## get db path 
my $title_file = $conf_file{map_title};

## check files 
Check($conf_file{map_class});

## help inofs 
&Help();


## outfile, name 
my $out_kopath="$opts{o}.kopath";
my $out_path="$opts{o}.path.xls";
my $name = (split /[\\\/]/, $opts{o})[-1];
$name =~ s/\.path$//;

## edge species infos
my %edge_path = ReadEdge($conf_file{edge},$opts{e}) if(defined $opts{e});

## title info 
my %title;
my %no_human_disease = read_title($title_file,\%title);
my %symbol = ReadIdmap($opts{symbol}) if(exists $opts{symbol});

if(exists $opts{r} && $opts{r} !~ /animal/i){
    $opts{n} = 1;
}
if($type eq 'OLD' && ($opts{e} !~ /animal/i && $opts{e} !~ /^vertebrates/i && $opts{e} !~ /^all$/i)){
    $opts{n} = 1;
}

my %hash;
my  $kopath;
#id 2 kopath file
if ($type eq 'ID'){
    #read ref
    &read_ref($opts{r},\%hash,$opts{t});
    open IN,$opts{i} or die $!;
    open OUT,">","$out_kopath";
    print OUT"###version $version\n";
    my $a = 0;
    my ($id,$ko);
    my %exists;
    while(<IN>){
        chomp;
        next if (/^#/);
        next if (/^\s*$/);
        my @aa = split/\t/,$_;
        next if(exists $exists{$aa[0]});
        $exists{$aa[0]} = 1;
        if (@aa == 2){ # only line with 2 columns will be used
            ## id 2 ko and path  (id,ref hash, org pfx())
            my $out = &id2ko_path($_,\%hash,$opts{p});
            print OUT"$out\n";
        }else{
            print OUT "$_\n";
        }
    }
    close IN;
    close OUT;
}

if ($type eq 'OLD') {
	my %hash2;
	&read_ko_tab($conf_file{kotab},\%hash2);
	open IN,$opts{i} or die $!;
	open OUT,">","$out_kopath" or die "$!:$out_kopath";
	print OUT"###version $version\n";
	my %exists;
	while(<IN>){
		chomp;
		next if (/^#/);
		next if (/^\s*$/);
		my @aa = split/\t/,$_;
		next if(exists $exists{$aa[0]});
		$exists{$aa[0]} = 1;
		if (@aa >= 2) {
            $aa[1] =~ s/\|.*//;
		    my $out;
		    if(exists $hash2{$aa[1]}){ ## KO => path 
                my $path = filter_edge_path(\%edge_path,$hash2{$aa[1]});
#                print Dumper(%no_human_disease);die;
                $path = filter_edge_path(\%no_human_disease,$path) if($opts{n}); ## filter human disease
                if($path eq ''){
                    $out = "$aa[0]";
                }else{
                    $out = "$aa[0]\t$aa[1]\t$path";
                }
            }else{
                $out = "$aa[0]";
            }
		    print OUT"$out\n";
		}else{
			print OUT"$_\n";
		}
	}
	close IN;
	close OUT;
}

## filter pathway near edge species
#&filter_edge_path($type,\%opts,$edge) if(defined $opts{e});

#kopath file 2 path
$out_kopath = "$opts{f}" if ($type eq 'PATH');
my (%kos,%genes,%bg_kos,%bg_genes,%num_bg,%num_fg,%pq,$sum_fg,$sum_bg);
open PATH,">","$out_path" or die "$!: $out_path";
if (defined $opts{b}){
    $out_kopath = "$opts{o}.kopath";
    if (abs_path($opts{f}) ne abs_path($out_kopath)){
        if(`head -n 1 $opts{b}` =~ /###version/){
            `head -n 1 $opts{b} >$out_kopath`;
            `perl $apps{fish} -bf table -ff table $opts{f} $opts{b} >>$out_kopath`;
        }else{
            `perl $apps{fish} -bf table -ff table $opts{f} $opts{b} >$out_kopath`;
        }
        #`sed 1i"###version $version" -i $out_kopath`;
    }else{
        print "output kopath file: $out_kopath is the same as input fgfile. And will be used. please check\n";
    }
    $sum_fg = &gene_ko_path($out_kopath,\%kos,\%genes,\%num_fg);
    $sum_bg = &gene_ko_path($opts{b},\%bg_kos,\%bg_genes,\%num_bg);
    print PATH"KEGG_A_class\tKEGG_B_class\tPathway\t$name ($sum_fg)\tAll ($sum_bg)\tPvalue\tQvalue\tPathway ID\t$type_map{$opts{k}}[0]\t$type_map{$opts{k}}[1]\n";
    if($sum_fg > 0){
        &get_pq(\%kos,\%num_fg,\%num_bg,\%pq,$sum_fg,$sum_bg);
        for (sort { $pq{$a}{'p'} <=> $pq{$b}{'p'}  } sort keys %num_fg){
            $kos{$_} =~ s/\Q+\E$//;$genes{$_} =~ s/;$//;
            die "No title for [$_] in config file\n" if(!exists $title{$_});
#            print STDERR "No title for [$_] in config file\n" if(!exists $title{$_});
            print PATH"$title{$_}\t$num_fg{$_}\t$num_bg{$_}\t$pq{$_}{'p'}\t$pq{$_}{'q'}\t$_\t$genes{$_}\t$kos{$_}\n";
        }
    }else{
        print STDERR "Warnings: fg gene num is 0\n";
    }
}else{
    $sum_fg = &gene_ko_path($out_kopath,\%kos,\%genes,\%num_fg);
    print PATH"KEGG_A_class\tKEGG_B_class\tPathway\tCount ($sum_fg)\tPathway ID\t$type_map{$opts{k}}[0]\t$type_map{$opts{k}}[1]\n";
    for (sort { $num_fg{$b} <=> $num_fg{$a} } sort keys %num_fg){
        $kos{$_} =~ s/\Q+\E$//;$genes{$_} =~ s/;$//;
        die "No title for [$_] in config file\n" if(!exists $title{$_});
#        print STDERR "No title for [$_] in config file\n" if(!exists $title{$_});
        print PATH"$title{$_}\t$num_fg{$_}\t$_\t$genes{$_}\t$kos{$_}\n";
    }
}

close PATH;

&add_next(\%opts) if (defined $opts{a});

sub ReadEdge{
    my $edge_file = shift;
    my $sp_name = shift;
    $sp_name = lc($sp_name);

    open my $edge_fh,$edge_file or die "$!:$edge_file\n";
    my $edge_info = "";
    while(<$edge_fh>){
        chomp;
        next if(/^\s*$/);
        next if(/^#/);
        my @aa = split/\t/,$_;
        if(lc($aa[0]) eq $sp_name){ ## get sp info only 
            $edge_info = $aa[1];
            last;
        }
    }
    close $edge_fh;

    if($edge_info eq ''){ ## no info get 
        print STDERR "The edge species [$sp_name] you input is not exists, please check the ref file: $conf_file{edge_xls}\n";
        exit(1);
    }

    my %edge_path = map{$_=>1} (split/,/,$edge_info);
    
#    print Dumper(%edge_path);

    return %edge_path; 
}

sub filter_edge_path{
    my $edge_path_ref = shift;
    my $path_list = shift;


    my @paths = split/,/,$path_list;

    my @out;
    for my $path(@paths){
        next if(!exists $edge_path_ref->{$path});
        push @out,$path;
    }

    my $out = @out > 0 ? join(",",@out) : ""; 

    return $out;
}


sub get_pq{#calcute pq value,
#    my ($h_kos,$num,$num_all,$pq,$sum,$sum_all) = @_;
    my ($h_kos,$termdiff,$all_diff,$pq,$term,$sum_all) = @_;
    open PT, "> $opts{o}.ptemp.R" or die $!;
    foreach (sort keys %{$h_kos}){
#	print PT"phyper($$num{$_}-1,$$num_all{$_},$sum_all-$$num_all{$_},$sum,lower.tail=F)\n";
#	         phyper(m-1,M,N-M,n,lower.tail=F)
        print PT"phyper($$termdiff{$_}-1,$$all_diff{$_},$sum_bg-$$all_diff{$_},$sum_fg,lower.tail=F)\n";
    }
    close PT;
    #my $Rscript="/Bio/bin/Rscript" if (-x "/Bio/bin/Rscript");
    my @pValues = split /\n/, `$apps{Rscript} $opts{o}.ptemp.R | grep -v "WARNING" | awk '{print \$2}' 2> /dev/null`;
    my @qValues;
    if($#pValues == 0){
        @qValues = @pValues;
    }else{
        my $p = join(", ", @pValues);
        open QT, "> $opts{o}.qtemp.R" or die $!;
        print QT <<QT;
library(qvalue)
p <- c($p)
q <- qvalue(p, lambda = 0)
q[3]
QT
        close QT;
        @qValues = split /\n/, `$apps{Rscript} $opts{o}.qtemp.R 2> /dev/null | grep -v "WARNING"  | awk 'NR > 1 {for (x = 2; x <= NF; x++) print \$x}'`;
        `$apps{Rscript} $opts{o}.qtemp.R`;
    }
    my $i = 0;
    foreach (sort keys %{$h_kos}){
        $$pq{$_}{'p'}=$pValues[$i];
        $$pq{$_}{'q'}=$qValues[$i];
        $i++;
    }
    `rm $opts{o}.ptemp.R` if($#pValues > 0);
    `rm $opts{o}.qtemp.R`;
}

## 1. input koapth file
## 2. return hash : koid => K_id;K_id;...    eg.ko01100
## 3. return hash : koid => gene+gene+...
## 4. return hash : koid => count
sub gene_ko_path{
    my ($file,$h_kos,$h_genes,$h_count) = @_;

    open IN,$file or die "$!: $file";
    my %count;
    while (<IN>){
        chomp;
        next if(/^#/);
        my @aa = split/\t/,$_;
        next unless(@aa > 2);
        next if($aa[2] eq '-');
        my %exist1;
        foreach my $ko(split/,/,$aa[2]){
            next if(exists $exist1{$ko}); ## rm dup
            next if($ko !~ /^\d+$/);
            $exist1{$ko} = 1;
            $ko = "ko$ko"; ## 
            $$h_kos{$ko} .= "$aa[1]+";
            if(exists $opts{symbol}){
                $symbol{$aa[0]} //= "-";
                $h_genes->{$ko} .= "$aa[0]($symbol{$aa[0]});";
            }else{
                $$h_genes{$ko} .= "$aa[0];";
            }
            $$h_count{$ko}++;
        }
        $count{$aa[0]}++;
    }
    close IN;

    my $num = scalar(keys %count);
    return $num;
}

sub read_ko_tab{
	my $file = shift;
	my $hash = shift;
    print "read_ko_tab [$file]\n";
	open FILE,$file or die $!;
	while(<FILE>){
		chomp;
		my @aa = split/\t/,$_;
		$$hash{$aa[0]} = $aa[1];
	}
	close FILE;
}

sub read_title{# title file, return hash
    my $file = shift;
    my $hash = shift;

    local($_);
    my %no_human_disease;
    open IN,$file or die $!;
    while (<IN>){
        chomp;
        s/\s*$//;
        my @aa = split/\t/,$_,2;
        if($aa[1] !~ /^Human Diseases/){
            $no_human_disease{$aa[0]} = 1;
        }
        $aa[0] = "ko$aa[0]";
        $$hash{$aa[0]} = $aa[1];
    }
    close IN;

    return %no_human_disease;
}

sub id2ko_path{ #line info contains id, hash contains ref inf, pfx just use pfx:id as kegid
    my ($line,$hash,$sp) = @_;
    my @aa = split/\t/,$line;

    my $out;

    my ($id,$ko,$path,$keg) = ($aa[0],"","","");
    if ($aa[1] =~ /.*\|(.*\:.*?)\|/){#blastout
        $aa[1] = $1;
    }elsif ($sp) {
        $aa[1] = "$sp:$aa[1]";
    }
    $aa[1] =~ s/\s.*//;
    $aa[1] =~ s/\.\d+//;#pepid with version like XP_0012313.1 
    if (exists $$hash{$aa[1]}){
        $ko = "$$hash{$aa[1]}{'kos'}";
        $path = "$$hash{$aa[1]}{'genes'}";
        $keg = $$hash{$aa[1]}{'keg'};
    }

    $path = filter_edge_path(\%edge_path,$path) if(defined $opts{e} && $path ne '');
    $path = filter_edge_path(\%no_human_disease,$path) if($opts{n}); ## filter human disease
    if($path ne ''){
        $out = join("\t",$id,$ko,$path,$keg);
    }else{
        $out = "$aa[0]";
    }

    return $out;
}


sub read_ref{# ref file, hash address return, defined the col read
    my ($file,$hash,$type) = @_;

    if($file =~ /.gz$/){
        open REF,"gzip -cd $file|" or die "$!:$file\n";
    }else{
        open REF,$file or die "$!:$file\n";
    }
    my %hash = ('keggid' => 0, 'geneid' =>3, 'gi' =>4, 'pepid' =>5);
    my @bb;
    if ($type eq 'none'){
	@bb = (0);
    }else{
	@bb = (0,$hash{"$type"});
    }
    while (<REF>){
	chomp;
	my @aa = split/\t/,$_;
	foreach (@bb){
	    $$hash{$aa[$_]}{'kos'} = $aa[1];
	    $$hash{$aa[$_]}{'genes'} = $aa[2];
	    $$hash{$aa[$_]}{'keg'} = $aa[0];
	    #$$hash{$aa[$_]} = $aa[1];
	    #$$hash2{$aa[$_]} = $aa[2];
	    #$$kegid{$aa[$_]} = $aa[0];
	}
    }
    delete $$hash{'-'} if (exists $$hash{'-'});
    close REF;
}


sub add_next{
	my $opts = shift;
	
#    Excu("perl $src/keggGradient.pl $$opts{o}.path.xls 20 Q");
	my $outdir = dirname($opts->{o});

#    my $diff_para = exists $opts->{f} ? "-diff $opts->{f} " : "";    
    my $diff_para = exists $opts->{f} ? "" : "";    
    my $ko = "$$opts{o}.kopath";
    if(exists $opts->{f} && !exists $opts{b}){
        $ko = $opts{f}
    }
    Excu("perl $src/keggMap.pl -ko $ko $diff_para -outdir $$opts{o}_map");
    Excu("perl $src/genPathHTML.pl -indir $outdir -k $opts->{k}");
#    Excu("perl $src/draw_PathwayClass.pl -f $opts->{o}.path.xls -p $opts->{o}.path.xls -k $opts->{k}");
    Excu("perl $src/R_plot.pl kegg_PathwayClass -infile $opts->{o}.path.xls -name_key $opts->{k} -outpfx $opts->{o}.path.xls");
    
#	    if (defined $$opts{b}){
#		`perl $RealBin/keggMap_nodiff.pl -ko $$opts{o}.kopath -outdir $$opts{o}_map`;
#	    }elsif(defined $$opts{f}){
#		`perl $RealBin/keggMap_nodiff.pl -ko $$opts{f} -outdir $$opts{o}_map`;
#	    }else{
#		`perl $RealBin/keggMap_nodiff.pl -ko $$opts{o}.kopath -outdir $$opts{o}_map`;
#	    }
    return;
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";

    system($cmd);

    return;
}

sub StatKNum{
    my $file = shift;

    local($_);
    open FILE,$file or die "$!: $file";
    my $all = 0;
    my %exist;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        if(@aa > 2){
            if($aa[1] =~ /^[CK][0-9]{5}$/ && $aa[2] =~ /\d+/){
                next if(exists $exist{$aa[0]}); ## rm dup id 
                $exist{$aa[0]} = 1;
                $all++;
            }
        }
    }
    close FILE;

    return $all;
}

sub ReadIdmap{
    my $file = shift;
    
    local($_);
    my %hash;
    open my $file_fh,$file or die $!;
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/,$_;
        $hash{$aa[0]} = $aa[1];
    }
    close $file_fh;

    return %hash;
}

sub Check{
    my $map_class_dir = shift;


    ## ref file: map_class, for find ko and path
    if (defined $opts{r}){
        my $tmp = $opts{r};
        $opts{r} = "$map_class_dir/$opts{r}.tab";
        $opts{r} = "$map_class_dir/$tmp.tab.gz" if(-e "$map_class_dir/$tmp.tab.gz");
        die "No [$tmp] file [$opts{r}]\n" if(!-e $opts{r});
    }
    if ( !-e $title_file ){
        print "=>No such file:$title_file\n         it uses for add pathway title in .path.xls file\n";
        exit 1;
    }
    
    return;
}


sub help_all{

    my $usage_all=<<USAGE;

    description:get ko and path from id or get path from ko

    usage:perl $0 TYPE  -i input -o output(.kopath .path)
              ID    get gene K number and path by (gene\\tid)(geneid,keggid,blastout)
              PATH  use (gene\\tko\\tpath) file to get path file
              OLD   old method: use KO find path

    example: perl $0 ID/PATH

USAGE
    
    print $usage_all;
    
    exit 0;
}

sub help_path{

    my $usage_path=<<USAGE;
    
    usage:perl $0 PATH -f fg -b bg -o output
           -f .kopath file / genelist.
        optional:
           -b<str> background file 
           -e<str> filter the (kopath)pathway with the Near edge species(See: $RealBin/../near_edge_species.xls)
           -a<str> [diff/nodiff] enrich process
           -k<str> keys for display [Genes],Proteins

    example:perl $0 PATH -f test1.kopath -o test2
            perl $0 PATH -f test.100 -b test1.kopath -o test3

USAGE
    print $usage_path;
    exit 0;
}
sub help_id{

    my $usage_id=<<USAGE;

    usage:perl $0 ID -i input -r ref -o output(.kopath .path)

             -i<str>  blast2ko result or gene\\id file
             -r<str>  ref db:animal,plant,fungi,micro,bacteria,archaea,all
                      ref:keggid\\tKO\\tko\\tncbigeneid\\tgi\\tncbiproteinid\\t

         optional:
             -p<str> short name of sp. (like:hsa).=>only when using "pfx:id" to find K number and path
             -t<str> choose id type to find KO/path,(geneid,gi,pepid,none,default:geneid)
             -n      no Human disease output
             -e<str> filter the pathway with the Near edge species(See: $conf_file{edge_xls})
             -a<str> [diff/nodiff] enrich process
             -k<str> keys for display [Genes],Proteins

             -v<str> version 

    example:perl $0 ID -i input -r animal -o test

USAGE
    print $usage_id;
    
    exit 0;
}

sub help_old{

    my $usage_old=<<USAGE;

    usage:perl $0 OLD -i input -e species type  -o output(.kopath .path)

             -i<str>  gene KO file
             -e<str>  the Near edge species(See: $conf_file{edge_xls})
             -o<str>  pfx of output
             -n       no Human disease output
             -k<str> keys for display [Genes],Proteins

    example:perl $0 OLD -i input -e Animal -o test

USAGE
    
    print $usage_old;
    
    exit 0;
}

sub Help{

    &help_all() if ($opts{h});
    &help_all() unless (defined $type);
    &help_all() unless ($type eq 'ID' || $type eq 'PATH' || $type eq 'OLD');
    &help_id() if ($type eq 'ID' and (!defined $opts{i} or !defined $opts{r} or defined $opts{b} or defined $opts{f}));
    &help_path() if($type eq 'PATH' and (!defined $opts{f}));
    &help_old() if($type eq 'OLD' and (!defined $opts{i} or !defined $opts{e}));
}
