#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin $Script/;
use Cwd qw/abs_path/;
use File::Basename qw/basename dirname/;
use lib "$RealBin/../";
use Carp qw(carp confess croak);
use Support_Program;

my $src = "$RealBin";
my $R_plot = "$src/R_plot.pl";
sub usage {
		my $use = <<USAGE;
usage:
        perl $0 ( [options] ) -c <con1&treat1,con2&treat2...> -exp <file> -cls <file>
        
        -diff            <FILE>  diff info, contains col [ "Pvalue","FDR" ]
        -dataset|s       <STR>   select one or more of "h,c1,...,c7", default is 'all', separate with comma;
        -minSet          <INT>   Gene sets smaller than this number are EXLCUDED from the analysis, Default: 15;
        -maxSet          <INT>   Gene sets larger than this number are EXLCUDED from the analysis, Default: 500;
        -permutation|p   <INT>   permutation test times, default 1000;
        -outdir|o        <DIR>   output dir, default ./;
        -contrast|c    * <STR>   set contrasts, format: -c "control&treatment", multi contrasts separate with comma; if compare to all others, use 'RSET', like "group1&REST";
        -report_num      <INT>   Plot Top N results, default 20;
        -rnd_seed        <INT>   Seed to use for randomization (a long number), default: 'timestamp';
        -exp           * <FILE>  expression file, see '-h' for more details;
        -cls           * <FILE>  grouping file, see '-h' for more details;
        -all_xls                 output all_xls 
        -help|h                  full help information.

Notice:
        All paths contain only [a-zA-Z0-9_.], other characters may cause GSEA breakdown.

USAGE
        my $full = <<HELP;
format:
    Expression file(exp):
        [Gene Symbol]  [Description]  [group1-1]  [group1-2] ... [group2-1]  [group2-2] ...

        1. make sure your [Gene Symbol] match the names in $RealBin/check.name.txt ( can use $RealBin/check.name.pl )
        2. [Description] is useless, just a place holder;
        3. other columns are expression value;
        4. separate with tab.

    Grouping file(cls):
        (line 1)[sample number] [group number] 1
        (line 2)# [group1 name] [group2 name] ...
        (line 3)[group1 name] [group1 name] ... [group2 name] [group2 name] ...

        1. separate with space;
        2. line 1 have three number, first is [sample number] which match the number of expression columns in Expression file; second is [group number]; and the last is always '1';
        3. line 2 start with '#', then follow your [groups name];
        4. line 3 is the project of expression columns of Expression file, describing how they grouping.

        Read more format details in http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

HELP
    print $use;
    print $full if $_[0];
    exit;
}

&usage(0) if(@ARGV<1);

my ( $contrasts, $sets, $perm, $cls, $minSet, $maxSet, $report_num, $rnd_seed );
my ($diff, $outdir, $outpfx, $idmap, $enrich_gmt);
my ($compare,$group,$exp,$gmt);
my $symbol = "none";
my $run_sge = "yes";
my $run = 1;
my ($run_ko,$run_go,$run_reactome,$run_do, $run_symbol) = (1,1,0,0,0);
my $all_xls = 0;
my $queue = "all.q,avx.q";
GetOptions( 
    "outpfx=s" => \$outpfx,
    "idmap=s" => \$idmap,
    "diff=s" => \$diff,
    "contrast|con|c=s"=>\$contrasts, 
    "dataset|set|s=s"=>\$sets, 
    "permutation|perm|p=s"=>\$perm,
    "minSet|min=i"=>\$minSet,
    "maxSet|max=i"=>\$maxSet,
    "outdir|out|o=s"=>\$outdir, 
    "report_num=i"=>\$report_num,
    "rnd_seed=i"=>\$rnd_seed,
    "exp=s"=>\$exp,
    "compare=s" => \$compare,
    "group=s" => \$group,
    "gmt=s" => \$gmt,
    "enrich_gmt=s" => \$enrich_gmt,
    "symbol=s" => \$symbol,
    "cls=s"=>\$cls,
    "all_xls=s" => \$all_xls,
    "run=s" => \$run,
    "run_ko=s" => \$run_ko,
    "run_go=s" => \$run_go,
    "run_reactome=s" => \$run_reactome,
    "run_do=s" => \$run_do,
    'run_symbol=s' => \$run_symbol,
    "run_sge:s"=> \$run_sge,
    "queue=s" => \$queue,
    "help|h"=>sub{&usage(1)}
);

#&usage(0) unless $contrasts and $exp and $cls;
my %apps;
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");
#Support_Program("$RealBin/../",\%apps,1);

$minSet //= 15;
$maxSet //= 500;
$perm   //= 1000;
$outdir //= '.';
$report_num //= 20;
$rnd_seed //= 'timestamp';
mkdir $outdir if(!-e $outdir);
$outdir   = abs_path($outdir);
#$exp      = abs_path($exp);
#$cls      = abs_path($cls);

#my @sets = ( !$sets or $sets =~ /all/i ) ? ( 'h','c1','c2','c3','c4','c5','c6','c7' ) : split /,/,$sets;
my @sets = ( !$sets or $sets =~ /all/i ) ? "msigdb" : split /,/,$sets;
print "$sets\t@sets\n";

my $gmt_dir = $apps{gsea_gmt};#"/Bio/User/hongshimiao/pipeline/db/GSEA/20190107/msigdb_v6.2_GMTs/";
my %gmt = map{$_ => "$gmt_dir/$_.all.v6.2.symbols.gmt"}qw/h c1 c2 c3 c4 c5 c6 c7/;
$gmt{msigdb} = "$gmt_dir/msigdb.v6.2.symbols.gmt";
# print OUT "java -cp /home/xushuyang/software/GSEA/gsea2-2.2.4.jar -Xmx1024m xtools.gsea.Gsea -gmx /home/xushuyang/software/GSEA/msigdb_v6.0_files/msigdb_v6.0_GMTs/$set.all.v6.0.symbols.gmt -nperm $perm -collapse false -set_min $minSet -set_max $maxSet -permute gene_set -plot_top_x $report_num -rpt_label ${treat}_versus_${con}.$set.$perm -rnd_seed $rnd_seed -res $exp -cls $cls#${treat}_versus_${con} -out $outdir\n";
my $gsea = "$apps{java} -cp $apps{gsea}";#/Bio/User/hongshimiao/software/GSEA/gsea2-2.2.4.jar";
my $gsea_para = "-collapse false -nperm $perm -set_min $minSet -set_max $maxSet -plot_top_x $report_num -rnd_seed $rnd_seed";

my ($symbol_para,%gmt_list) = GmtInfos($outdir);

if(defined $diff){
    DiffGSEA($diff,$outdir,$outpfx);
}

if(defined $compare && defined $group && $exp){
    AllExpGSEA($exp,$compare,$group,$outdir,$symbol_para,%gmt_list);
}

sub GmtInfos{
    my $outdir = shift;
 
    ($gmt,$symbol) = EnrichGmt($enrich_gmt) if(defined $enrich_gmt);
    my @gmts = split/,/,$gmt if(defined $gmt);
    
    my $version = "v6.2";
    if(defined $sets){
        my %cc = ("C1(positional)"=>"c1.all","C2(curated) chemical and genetic perturbations"=>"c2.cgp",
        "C2(curated) BioCarta"=>"c2.cp.biocarta","C2(curated) KEGG"=>"c2.cp.kegg",
        "C2(curated) Reactome"=>"c2.cp.reactome","C2(curated) canonical pathways"=>"c2.cp",   
        "C3(motif) miRNA Target"=>"c3.mir","C3(motif) transcription factor targets"=>"c3.tft",
        "C4(computational) cancer gene neighborhoods"=>"c4.cgn","C4(computational) cancer modules"=>"c4.cm",
        "C5(GO) cellular components"=>"c5.cc","C5(GO) molecular functions"=>"c5.mf",
        "C5(GO) biological processes"=>"c5.bp",
        "C6(oncogenic signatures)"=>"c6.all","C7(immunologic signatures)"=>"c7.all",
        "H(hallmark)"=>"h.all",'GSEA MSigDB'=>"msigdb");
        my $gmt_dir = $apps{gsea_gmt};#"/Bio/User/hongshimiao/pipeline/db/GSEA/20190107/msigdb_v6.2_GMTs/";
#        $sets = "h.all,c1.all,c2.cgp,c3.mir,c3.tft,c4.cgn,c4.cm,c5.cc,c5.mf,c5.bp,c6.all,c7" if($sets eq 'all');
        $sets = join(",",keys %cc) if($sets eq 'all');
        my @sets =  map{$cc{$_} || $_} split /,/,$sets;
        push @gmts,map{"$gmt_dir/$_.$version.symbols.gmt"}@sets;
#        print "@gmts\n";die;
    }

    my %gmt_list;
    for my $g(@gmts){
        print "$g\n";
        $g = abs_path($g);
        my $gn = basename($g);
        $gn =~ s/$version.symbols.//;
        Error("gmt file [$g] is not exists") if(!-e $g);
        Excu("ln -sf $g $outdir/$gn",1) if($g ne "$outdir/$gn");
#        next if(CheckGmtSize($g,$exp) == 0);
        $gmt_list{$gn}{origin_path} = $g;
        $gmt_list{$gn}{name} = $gn;
        $gmt_list{$gn}{name} =~ s/.gmt$//;
        push @{$gmt_list{gmt_all}},$gn;
    }
    my $symbol_para = "";
    if(defined $symbol && $symbol ne 'none'){
        Error("symbol file [$symbol] is not exists") if(!-e $symbol);
        Excu("ln -sf $symbol $outdir/all.symbol.chip",1);
        $symbol_para = "-chip all.symbol.chip";
    }

    return ($symbol_para,%gmt_list);
}

sub AllExpGSEA{
    my $exp = shift;
    my $compare = shift;
    my $group = shift;
    my $outdir = shift;
    my $symbol_para = shift;
    my %gmt_list = @_;

    FormatExp($exp,$group,$outdir);

    open COM,$compare or die $!;
    open OUT,">$outdir/run.sh" or die $!;

    my @result_dirs;
    my @gsea_xls_shell;
    my @replot_shell;
    my @all_xls_shell;
    while(<COM>){
        chomp;
        my @names = split/\t/,$_;
        s/-/_/g;
        my @aa = split/\t/,$_;

        for my $gmt(@{$gmt_list{gmt_all}}){
            next if(CheckGmtSize("$outdir/$gmt","$outdir/all.exp.gct") == 0);
            my $out_name =  "$names[0]-vs-$names[1].$gmt_list{$gmt}{name}.Gsea";
            push @result_dirs,"$outdir/$out_name";
            push @gsea_xls_shell,"perl $src/gsea_xls.pl $gmt_list{$gmt}{origin_path} $symbol $outdir/$out_name.xls $outdir/$out_name/gsea_report_for_*xls";
            push @replot_shell,"perl $src/ReplotGSEA.pl $outdir/$out_name $outdir/$out_name.xls $outdir/$out_name $R_plot";
            push @all_xls_shell,"perl $src/GSEA_data.pl $outdir/$out_name/edb all $outdir/$out_name.all.xls.gz";
            if(-d "$outdir/$out_name"){
                my $old = "$outdir/$out_name".time();
                Excu("mv $outdir/$out_name $old");
            }

            my $max_mem = $gmt_list{$gmt}{name} eq 'go' ? 1024*15: 1024*5;
            $max_mem = "-Xmx${max_mem}m";
            print OUT "cd $outdir; $gsea $max_mem xtools.gsea.Gsea -gmx $gmt $gsea_para -permute gene_set  -cls all.cls#$aa[1]_versus_$aa[0]  -rpt_label $aa[0]_vs_$aa[1].$gmt_list{$gmt}{name} $symbol_para -res all.exp.gct -out ./; mv $aa[0]_vs_$aa[1].$gmt_list{$gmt}{name}.Gsea.* $out_name\n";
        }
    }
    close COM;
    close OUT;

#    Excu("$apps{qsub_sge} --convert no --maxjob 10 --jobprefix GSEA $outdir/run.sh");
    my $qsub_sge = $run_sge eq 'pbs' ? "perl $apps{qsub_pbs} --convert no --maxjob 4" : "perl $apps{qsub_sge} --convert no --queue $queue --maxjob 10 --reqsub";
    if($run_sge eq 'sh'){
        Excu("sh $outdir/run.sh");
    }elsif($run_sge ne 'no'){
        Excu("$qsub_sge --jobprefix GSEA $outdir/run.sh");
    }

    for my $shell(@gsea_xls_shell){
        Excu($shell);
#    for my $result_dir(@result_dirs){
#        Excu("perl $src/gsea_xls.pl $result_dir.xls $result_dir/gsea_report_for_*xls");
    }
    for my $shell(@replot_shell){
        Excu($shell);
    }

    for my $shell(@all_xls_shell){
        Excu($shell);
    }

    return;
}

sub EnrichGmt{
    my $enrich_gmt = shift;

    my @gmts;
    push @gmts,"$enrich_gmt/kegg.gmt" if($run_ko == 1);
    push @gmts,"$enrich_gmt/go.gmt" if($run_go == 1);
    push @gmts,"$enrich_gmt/do.gmt" if($run_do == 1);
    push @gmts,"$enrich_gmt/reactome.gmt" if($run_reactome == 1);

    my $symbol = "$enrich_gmt/all.chip";
    my $gmt = join(",",@gmts);

    return ($gmt,$symbol);
}

sub CheckGmtSize{
    my $gmt = shift;
    my $exp = shift;
    print "=====>$exp\n";
    my $out = 0;
    
    open EXP,$exp or die $!;
    my %id_list;
    while(<EXP>){
        chomp;
        my @aa = split/\t/,$_;
        $id_list{$aa[0]} = 1;
    }
    close EXP;

    my $info;
    open my $gmt_fh,$gmt or die "$!:$gmt";
    while(<$gmt_fh>){
        chomp;
        my @aa = split/\t/,$_;
        my $id = $aa[0];
        @aa = FilterID(\%id_list,@aa[2..$#aa]);
#        if(@aa >= $minSet+2 && @aa <= $maxSet+2){
        if(@aa >= $minSet && @aa <= $maxSet){
            $out = 1;
        }
        $info .= "$id\t".(scalar(@aa))."\n";
    }
    close $gmt_fh;
    print STDERR "No Set in [$minSet, $maxSet] of [$gmt]\n$info\n" if($out == 0);

    return $out;
}

sub FilterID{
    my $id_list_ref = shift;
    my @ids = @_;

    my @out;

    for my $id(@ids){
        next if(!exists $id_list_ref->{$id});
        push @out,$id;
    }

    return @out;
}

sub FormatExp{
    my $exp = shift;
    my $group = shift;
    my $outdir = shift;

    my %group;
    open GROUP,$group or die $!;
    my %count;
    while(<GROUP>){
        chomp;
        my @aa = split;
        $aa[1] =~ s/-/_/g;
        $group{$aa[0]} = $aa[1];
        $count{sample}{$aa[0]}++;
        $count{group}{$aa[1]}++;
    }
    close GROUP;
    
    open CLS,">$outdir/all.cls" or die $!;
    my $sample_num = scalar(keys %{$count{sample}});
    my $group_num = scalar(keys %{$count{group}});
    print CLS "$sample_num $group_num 1\n";
    my @groups;
    my @samples;

    my %symbol;
    if($run_symbol == 1){
        open SYMBOL,$symbol or die $!;
        while(<SYMBOL>){
            chomp;
            my @aa = split/\t/,$_;
            $symbol{$aa[0]} = uc($aa[1]);
        }
        close SYMBOL;
    }

    open EXP,$exp or die $!;
    my @cols;
    open OUT,">$outdir/all.exp.gct" or die $!;
    my $gene_num = -1;
    my $tmp;
    my %e_group;
    my %e_name;
    while(<EXP>){
        chomp;
        my @aa = split;
        my $filter_na = 0;
        if($. == 1){
            my @cls_line3;
            for my $i(0..$#aa){
                if($aa[$i] =~ s/_fpkm// || $aa[$i] =~ s/_rpkm// || exists $group{$aa[$i]}){
                    my $group = $group{$aa[$i]};
                    push @samples,$group;
                    push @groups,$group if(!exists $e_group{$group});
                    $e_group{$group} = 1;
                    push @cols,$i;
                }
            }
        }else{
            for my $i(@cols){
                if($aa[$i] eq 'NA'){
                    $filter_na = 1;
                }
            }
            next if($filter_na == 1);
        }
        if($run_symbol == 1){
            next if(!exists $symbol{$aa[0]} || $symbol{$aa[0]} eq '-');
#            next if(!exists $symbol{$aa[0]});
            $aa[0] = $symbol{$aa[0]};
        }
        next if(exists $e_name{$aa[0]});
        $e_name{$aa[0]} = 1;
        $gene_num++;
        $tmp .= join("\t",$aa[0],$aa[0],@aa[@cols])."\n";
    }
    close EXP;
    print CLS "# @groups\n";
    print CLS "@samples\n";
    print OUT "#1.2\n";
    print OUT "$gene_num\t".scalar(@cols)."\n";
    print OUT "$tmp";

    close OUT;
    close CLS;

    return;
}


sub DiffGSEA{
    my $diff = shift;
    my $outdir = shift;
    my $outpfx = shift;
    
    mkdir "$outdir" if(!-e $outdir);

    $outpfx =~ s/-/_/g;
    if($outpfx =~ /\W/){
        die "[$outpfx] can only contants [a-zA-Z0-9_.]\n";
    }

    ## pre diff, diff file is ordered by pvalue 
    my ($pvalue_col,$qvalue_col) = GetCol($diff,"Pvalue|PValue","FDR");
    Excu("sort -k${pvalue_col}g $diff | awk 'NR>1{print \$1\"\\t\"1-\$$qvalue_col}' >$outdir/$outpfx.rnk");
#    Excu("awk 'NR>1{print \$1\"\\t\"-NR}' $diff >$outdir/$outpfx.rnk");

    ## run GSEA
    for my $set(@sets){
        my $gmt = $gmt{$set};
        if(defined $idmap){
            Excu("perl $src/reform_gmt.pl $gmt $idmap $outdir/$outpfx.$set.gmt");
            $gmt = "$outpfx.$set.gmt";
        }
        Excu("cd $outdir && $gsea xtools.gsea.GseaPreranked $gsea_para -rnk $outpfx.rnk -gmx $gmt -rpt_label $outpfx"); 
        my $result_dir = "$outdir/GSEA_$outpfx";
        Excu("rm -rf $result_dir") if(-e $result_dir);
        Excu("mv $outdir/*/$outpfx.GseaPreranked* $result_dir");
        Excu("rm -rf $outdir/$outpfx.$set.gmt") if(defined $idmap);
    }

    return;
}

sub GetCol{
    my $file = shift;
    my @colnames = @_;

    open FILE,$file or die $!;
    my $head = <FILE>;
    close FILE;
    chomp($head);
    my @head = split/\t/,$head;
    my %head = map{$head[$_] => $_+1} 0..$#head;

    my @out;
    for my $name_list(@colnames){
        my @names = split/\|/,$name_list;
        my $ok = 0;
        for my $name(@names){
            if(exists $head{$name}){
                push @out,$head{$name};
                $ok = 1;
            }
        }

        die "No col name [$name_list]\n" if($ok == 0);
    }

    return @out;
}

sub Excu{
    my $cmd = shift;
    my $all_run = shift;

    $all_run //= 0;

    print "$cmd\n";
    if($run == 1 || $all_run == 1){
        my $ret = system($cmd);
        Error("Error:$ret [$cmd]") if($ret != 0);
    }

    return;
}

sub Error{

    confess("Error:@_\n");
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
