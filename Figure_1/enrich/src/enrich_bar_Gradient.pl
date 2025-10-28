#! /usr/bin/perl
use utf8;
use strict;
use warnings;
use Getopt::Long;
use FindBin'$RealBin';
use Carp qw(carp confess croak);
use Data::Dumper;

die "    Usage: draw pq barplot &&  gradient plot
    perl $0 <kegg.path.xls/ [PFC].xls> <outpfx>  <the first number> <enrich_type: [KEGG/GO]  <Q/P> <name_key>\n" if(@ARGV < 3);

## file_type pqbar_type1 : id<tab>pvalue<tab>per
my %opts = (file_type=>'pipe',draw_type=>'all',"rev"=>0,"run"=>1,pfc=>1,fig=>1,
            "pq_cut"=>0.05,"key_word"=>"Genes","top"=>20,"pq"=>"q","name"=>"DGEs");
GetOptions(\%opts,"file_type=s","draw_type=s","rev=s","glist=s","run=s","pq_cut=s","top=s","pq=s","key_word=s","pfc=s","fig=s","name=s");

my $infiles = shift;
my $outpfx = shift;
my $enrich_type = shift;

my $first_num = $opts{top};
my $pq_type = lc($opts{pq});
my $name_key = $opts{key_word};

my %pfc = ("P"=>"Biological Process","F"=>"Molecular Function","C"=>"Cellular Component");
my %glist = ReadGlist($opts{glist});

my $tmp_xls;
if($opts{file_type} eq 'circle_type2' || $opts{file_type} eq 'pqbar_type1'){
    $tmp_xls = $infiles;
}elsif($opts{file_type} eq 'gradient_type1'){
    open my $in_fh,"$infiles" or die "$!:$infiles\n";
    open OUT,">$outpfx.bar_Gradient.xls" or die $!;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        if($. == 1){
            @aa = ("id","num","${pq_type}value","per","ratio");
        }else{
            @aa = ($aa[0],$aa[1],$aa[3],printf("%.3f",$aa[1]/$aa[2]*100),sprintf("%.3f",$aa[1]/$aa[2]));
        }
        print OUT join("\t",@aa)."\n";
    }
    close $in_fh;
    close OUT;
    $tmp_xls = "$outpfx.bar_Gradient.xls";
#    Excu("cp $tmp_xls $outpfx.bar_Gradient.xls");
}else{
    my @files = split/,/,$infiles;  ## read multi files, like PFC draw one, no use now
    @files = map{"$infiles.$_.xls"}qw/P F C/ if($enrich_type eq 'GO' && $opts{pfc} == 1);
#    my %info = GetFileInfo($files[0],$enrich_type,$opts{file_type},$first_num);
    my @infos;
    for my $file(@files){
        push @infos,GetFileInfo($file,$enrich_type,$opts{file_type},$first_num);
    }
    my %info = @infos;
#    print Dumper(%info);

    $tmp_xls = "$outpfx.bar_Gradient.xls"; ## "$outpfx.tmp"
    open TMP,">$tmp_xls" or die $!;
    #                              class ID  Descrption  Pvalue  Qvalue  Up  Down
    #                              Cellular Processes    ko04530 Tight junction  0.01082009  0.1825742   1   1
    print TMP "id\tnum\t${pq_type}value\tper\tratio\tclass\tID\tDescrption\tPvalue\tQvalue\tUp\tDown\tfg_num\tbg_num\n";
    my $num = 0;
    my @ids = sort {$info{$a}{q} <=> $info{$b}{q} or $info{$a}{p} <=> $info{$b}{p} or $a cmp $b} keys %info;
    my ($min_big_0_pvalue,$min_big_0_qvalue) = MinPQ(%info);
    if(@files == 1){
        @ids = sort {$info{$a}{order} <=> $info{$b}{order}} keys %info;
    }
    for my $id(@ids){
        $num++;
        my ($up,$down) = GetUD(\%glist,$info{$id}{genes}) if(exists $info{$id}{genes});
        if(exists $info{$id}{up} && exists $info{$id}{down}){
            $up = $info{$id}{up};
            $down = $info{$id}{down};
        }
        
        $info{$id}{p} = $min_big_0_pvalue / 10 if($info{$id}{p} == 0);
        $info{$id}{q} = $min_big_0_qvalue / 10 if($info{$id}{q} == 0);
#        CheckKey($info{$id},qw/num per ratio class ID desc p q num all/,$pq_type);
#    last if($num > $first_num);
        print TMP join("\t",$id,$info{$id}{num},$info{$id}{$pq_type},$info{$id}{per},$info{$id}{ratio},
                           $info{$id}{class},$info{$id}{ID},$info{$id}{desc},$info{$id}{p},$info{$id}{q},
                           $up,$down,$info{$id}{num},$info{$id}{all}
                    )."\n";
    }
    close TMP;
}
exit(0) if($opts{fig} == 0);

my %ylab = qw/KEGG Pathway  GO GOterm  DO DOterm reactome Reactome/;
Error("No ylab type of [$enrich_type]") if(!exists $ylab{$enrich_type});
my $title = "Top $first_num of $enrich_type Enrichment";
my @types = ($opts{pfc} == 1 && $enrich_type eq 'GO') ? qw/P F C/ : ("all");
for my $type(@types){
#    print "$type";
    my $outpfx_tmp = $type eq 'all' ? $outpfx : "$outpfx.$type";
#    print "$outpfx_tmp\t$enrich_type\n";
    Excu("perl $RealBin/R_plot.pl pq_bar -infile $tmp_xls -pfc $type -outpfx $outpfx_tmp.barplot -name_key \"$name_key\" -top $first_num -ylab $ylab{$enrich_type} -title \"$title\"") if($opts{draw_type} eq 'all' || $opts{draw_type} =~ 'pqbar');
#$title = "Tor $first_num of $enrich_type Enrichment";
    Excu("perl $RealBin/R_plot.pl Gradient -infile $tmp_xls -pfc $type -outpfx $outpfx_tmp.gradient -reverse $opts{rev} -name_key \"$name_key\" -top $first_num -ylab $ylab{$enrich_type} -title \"$title\"") if($opts{draw_type} eq 'all' || $opts{draw_type} =~ 'gradient');
}
my %pq = qw/p Pvalue q Qvalue/;
Excu("perl $RealBin/R_plot.pl enrich_bubble -infile $tmp_xls -type $enrich_type -outpfx $outpfx.bubble -pq $pq{$pq_type} -pq_cut $opts{pq_cut}") if(exists $glist{diff} && ($opts{draw_type} eq 'all' || $opts{draw_type} =~ 'bubble'));

my %max_cc = qw/KEGG 30  GO 25  DO 25 reactome 20/;
my $circ_top_num = $first_num > $max_cc{$enrich_type} ? $max_cc{$enrich_type} : $first_num;
Excu("perl $RealBin/R_plot.pl enrichPlot -infile $tmp_xls -pq $pq_type -top $circ_top_num -outpfx $outpfx.circular -name_key \"$name_key\" -new 1") if($opts{draw_type} eq 'all' ||  $opts{draw_type} =~ 'circular');

sub MinPQ{
    my %info = @_;

    my @keys = sort {$info{$a}{p} <=> $info{$b}{p}} keys %info;
    my @big0 = grep{$info{$_}{p}>0}@keys;
#    print "$big0[0]\t$info{$big0[0]}{p}\t$keys[0]\t$info{$keys[0]}{p}\n";
    @keys = sort {$info{$a}{q} <=> $info{$b}{q}} keys %info;
    my @big0_q = grep{$info{$_}{q}>0}@keys;
#    print "$big0_q[0]\t$info{$big0_q[0]}{q}\t$keys[0]\t$info{$keys[0]}{q}\n";die;

    return ($info{$big0[0]}{p},$info{$big0_q[0]}{q});
}

sub GetFileInfo{
    my $file = shift;
    my $type = shift;
    my $file_type = shift;
    my $first_num = shift;

    my %hash;
    local($_);
    open my $file_fh,$file or die "$!:$file";
    my %col;
    my @keys = qw/num all p q class genes ID desc/;
    my %class_cc = (biological_process  => "Biological Process",
                    molecular_function  => "Molecular Function",
                    cellular_component  => "Cellular Component");
    if($file_type eq 'gradient_type1'){
        die "No config [$file_type]\n";
        # Pathway   OS-vs-QS (2187) All (8483)  Pvalue
        %col = qw/num 1  p 3 q 3 all 2/;
        @{$col{id}} = (0);
    }elsif($file_type eq 'pqbar_type2' || $file_type eq 'gradient_type2'){
#        die "No config [$file_type]\n";
        # Pathway   OS-vs-QS (2187) All (8483)  Pvalue  Qvalue  Pathway ID  Genes   K_IDs
        %col = qw/num 1  p 3 q 4 all 2 genes -1 ID -1 desc -1 class -1/;
        @{$col{id}} = (0);
    }elsif($file_type eq 'circle_type1'){
        %col = qw/p 3 q 3 all 2  genes -1 ID 0 desc -1 class 1 up 4 down 5/;
        @{$col{id}} = (0);
    }elsif($file_type eq 'diff_bubble_type1'){ ## class ID  Descrption  Up  Down    Pvalue
        %col = qw/class 0 ID 1 all -1 desc 2 up 3 down 4 p 5 q 5 pvalue 5 qvalue 5/;
        @{$col{id}} = (1);
        push @keys,qw/pvalue qvalue/;
        $glist{diff} = 1;
    }elsif($type eq 'KEGG'){ # col 0 base
        %col = qw/num 3 p 5 q 6 all 4 class 0 genes 8 ID 7 desc 2/;
        @{$col{id}} = (2);
    }elsif($type eq 'GO'){
        %col = qw/num 2 p 4 q 5 all 3 genes 6 ID 0 desc 1/;
        @{$col{id}} = (0,1);
    }elsif($type eq 'DO'){
        %col = qw/num 3 p 5 q 6 all 4 genes 7 ID 0 desc 1/;
        @{$col{id}} = (0,1);
    }elsif($type eq 'reactome'){
        %col = qw/num 3 p 5 q 6 all 4 genes 7 ID 0 desc 1/;
        @{$col{id}} = (0,1);
    }else{
        die "Unknown enrich type [$type]\n";
    }
    my $all = 1;
    my $i = -1;
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/,$_;
#        last if($i == $first_num);
        if($i == -1){
            $i++;
            ($all) = $aa[$col{num}] =~ /\((\d+)\)/ if(exists $col{num});
            next;
        }
        $i++;
        my $id = join(" ",@aa[@{$col{id}}]);
        for my $t(@keys){
            next if(!exists $col{$t});
            if($col{$t} eq '-1'){
                $hash{$id}{$t} = "--";
            }else{
                $hash{$id}{$t} = $aa[$col{$t}];
            }
        }
#        $hash{$id}{per} = exists $col{per} ? $aa[$col{per}] : 
#                          ( (defined $all && exists $col{num}) ? sprintf("%.3f",$aa[$col{num}]/$all*100) : "--");
        $hash{$id}{ratio} = exists $col{num} ? sprintf("%.3f",$aa[$col{num}]/$aa[$col{all}]) : "--";
        $hash{$id}{per} = exists $col{per} ? $aa[$col{per}] : 
                          ( (defined $all && exists $col{num}) ? sprintf("%.3f",$aa[$col{num}]/$all*100) : $hash{$id}{ratio});
        $hash{$id}{ratio} = exists $col{num} ? sprintf("%.3f",$aa[$col{num}]/$aa[$col{all}]) : "--";
        $hash{$id}{num} //= "--";
        $hash{$id}{order} = $i;

        ## up down 
        if(exists $col{up} &&  exists $col{down}){
            $hash{$id}{up} = $aa[$col{up}];
            if($#aa >= $col{down}){
                $hash{$id}{down} = $aa[$col{down}];
                $hash{$id}{num} = $hash{$id}{up} + $hash{$id}{down} if($hash{$id}{num} eq '--');
            }else{
                $hash{$id}{num} = $hash{$id}{up} if($hash{$id}{num} eq '--');
            }
        }

        # class 
        $hash{$id}{class} = "DO" if($enrich_type eq 'DO');
        $hash{$id}{class} = "Reactome" if($enrich_type eq 'reactome');
        if($file_type ne 'diff_bubble_type1'){
            if($enrich_type eq 'GO' && $file =~ /([PFC]).xls$/){
                $hash{$id}{class} = $pfc{$1};
            }elsif($enrich_type eq 'GO'){
                $hash{$id}{class} = "GO";
            }
        }
        my $class0 = lc($hash{$id}{class});
        $class0 =~ s/\s/_/g;
        if(exists $class_cc{$class0}){
            $hash{$id}{class} = $class_cc{$class0};
        }
    }
    close $file_fh;
#    print Dumper(%hash);

    return %hash;
}

sub GetUD{
    my $glist_ref = shift;
    my $genes = shift;
    
    my @genes = split/;/,$genes;
    my %hash;
    for my $gene(@genes){
        $gene =~ s/\(.*//;
        next if(!exists $glist_ref->{$gene});
        $hash{$glist_ref->{$gene}}++;
    }
    $hash{up} //= 0;
    $hash{down} //= 0;

    return ($hash{up},$hash{down});
}

sub ReadGlist{
    my $glist = shift;

    return if(!defined $glist);
    
    my %glist;
    open GLIST,$glist or die $!;
    while(<GLIST>){
        chomp;
        my @aa = split/\t/,$_;
        if(@aa > 1){
            next if($aa[1] =~ /log/i);
            $glist{$aa[0]} = $aa[1] > 0 ? "up" : "down";
            $glist{diff} = 1;
        }
    }
    close GLIST;

    return %glist;
}

sub Excu{
    my $cmd  = shift;
    print "$cmd\n";

    if($opts{run} == 1){
        my $ret = system($cmd);
        Error("Error $ret : [$cmd]") if($ret != 0);
    }

    return;
}

sub CheckKey{
    my $hash_ref = shift;
    my @keys = @_;
 
    my @err;
    for my $key(@keys){
        next if(exists $hash_ref->{$key});
        push @err,$key;
    }

    Error("No keys [@err] info") if(@err > 0);

    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}

