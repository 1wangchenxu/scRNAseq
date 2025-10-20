#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: reactome_enrich.pl
#
#        USAGE: ./reactome_enrich.pl  
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
#      CREATED: 12/27/2018 11:14:25 AM
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
use lib "$RealBin/../";
use Support_Program;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 

            -type  Ensembl,NCBI
        
    e.g: perl $0 <run_type> <>
                    ID
                   PATH 
                  ENRICH 

USAGE

die "$usage" if(@ARGV<1); 

#my $db_dir = "$RealBin/db/20181227";
my $src = "$RealBin";

my %opts = ("type" => "Ensembl","idmap"=>"none","k"=>"Genes","DEBUG"=>0,version=>"v2");
my $run_type = shift;
GetOptions(\%opts,"bg=s","i=s","type=s","idmap=s","outpfx=s","k=s","symbol=s","DEBUG=s","version=s");

## help 
die "$usage" if(!exists $opts{i} || !exists $opts{outpfx});

## soft 
my %apps;
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");
$opts{version} = ReadVersion(\%opts,$run_type);
print "reactome version $opts{version}\n";
my $db_dir = $apps{"reactome_db_dir_$opts{version}"};

## get back path info
my %bg_info;
if($run_type eq 'ID'){
    FormatDB($opts{type},$opts{idmap},$opts{i},"$opts{outpfx}.reactome");
    %bg_info = ReadIdmap("$opts{outpfx}.reactome");
}elsif($run_type eq 'PATH'){
    %bg_info = ReadIdmap("$opts{i}");
}else{
    %bg_info = ReadIdmap("$opts{bg}");
}

## read title 
my %g_title = ReadIdmap("$db_dir/reactome.title.txt");
my %symbol = ReadIdmap("$opts{symbol}") if(exists $opts{symbol});

##
my ($bg_all,%bg_stat_info) = StatPath(\%bg_info);
if($run_type eq 'ID' || $run_type eq 'PATH'){
    OutputPath("$opts{outpfx}.reactome.xls",\%bg_stat_info,'bg_all'=>$bg_all,"k"=>$opts{k});
}else{
    my @fg = GetGene($opts{i});
    my ($fg_all,%fg_stat_info) = StatPath(\%bg_info,\@fg);
    OutputPath("$opts{outpfx}.reactome.tmp",\%bg_stat_info,'fg'=>\%fg_stat_info,'bg_all'=>$bg_all,'fg_all'=>$fg_all,'name'=>$opts{outpfx},"k"=>$opts{k});
    Excu("Rscript $src/enrich.r $opts{outpfx}.reactome.tmp $opts{outpfx}.reactome.xls $fg_all $bg_all");
    Excu("rm -f $opts{outpfx}.reactome.tmp") if($opts{DEBUG}==0);
}



sub GetGene{
    my $file = shift;

    open my $file_fh,$file or die $!;
    my @out;
    my %e;
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/,$_;
        next if(exists $e{$aa[0]});
        $e{$aa[0]} = 1;
        push @out,$aa[0];
    }
    close $file_fh;

    return @out;
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    system($cmd);

    return;
}


sub FormatDB{
    my $type = shift;
    my $idmap = shift;
    my $idlist = shift;
    my $outfile = shift;

    my %idmap = ReadIdmap($idmap);  ## id <tab> dbid
    my %db = ReadIdmap("$db_dir/$type.txt"); ## id <tab> path info

    open ID,$idlist or die "$!:$idlist";
    open OUT,">$outfile" or die $!;
    print OUT "## version $opts{version}\n";
    while(<ID>){
        chomp;
        my @aa = split/\t/,$_;
        my $id = $aa[0];
        if($idmap ne 'none'){
            next if(!exists $idmap{$aa[0]});
            $id = $idmap{$aa[0]};
        }
        $id =~ s/^ncbi_//;
        next if(!exists $db{$id});
        print OUT "$aa[0]\t$db{$id}\n";
    }
    close ID;
    close OUT;


    return;
}

sub OutputPath{
    my $outfile = shift;
    my $bg_ref = shift;

    my %opts = @_;
#    my $fg_ref = shift;
    my $fg_flag = 1 if(exists $opts{fg});
    my $fg_ref = $opts{fg} if(exists $opts{fg});
    
    my @path = exists $opts{fg} ? keys %{$opts{fg}} : keys %{$bg_ref};

    if(!$fg_flag){
        @path = sort {$bg_ref->{$b}{num} <=> $bg_ref->{$a}{num} or $a cmp $b} @path;
    }
    my $url = "https://reactome.org/PathwayBrowser/#";

    open OUT,">$outfile" or die $!;
    my $title0 = "Reactome\tReactome_Name\tURL";
    if($fg_flag){
        $opts{name} =~ s/.*\///;
        print OUT join("\t",$title0,"$opts{name} ($opts{fg_all})","All ($opts{bg_all})","$opts{k}")."\n";
    }else{
        print OUT join("\t",$title0,"All ($opts{bg_all})","$opts{k}")."\n";
    }

    for my $path(@path){
        next if(!exists $g_title{$path});
#        $g_title{$path} //= "NA";
        if(defined $fg_ref){
            print OUT join("\t",$path,$g_title{$path},"$url/$path",$fg_ref->{$path}{num},$bg_ref->{$path}{num},
                        &GeneSymbol(@{$fg_ref->{$path}{id}}))."\n";
#                        join(";",sort @{$fg_ref->{$path}{id}}))."\n";
        }else{
            print OUT join("\t",$path,$g_title{$path},"$url/$path",$bg_ref->{$path}{num},
                        &GeneSymbol(@{$bg_ref->{$path}{id}}))."\n";
#                join(";",sort @{$bg_ref->{$path}{id}}))."\n";
        }
    }
    close OUT;
    

    return;
}

sub GeneSymbol{
    my @genes = @_;

    @genes = sort @genes;
    if(exists $opts{symbol}){
        my @tmp;
        for my $gene(@genes){
            $symbol{$gene} //= "-";
            push @tmp,"$gene($symbol{$gene})";
        }
        @genes = @tmp;
    }

    return join(";",@genes);
}

sub StatPath{
    my $bg_ref = shift;
    my $fg_arr = shift;

    my @keys = defined $fg_arr ? @{$fg_arr} : keys %$bg_ref;

    my %stat;
    my $all = 0;
    for my $key(@keys){
        next if(!exists $bg_ref->{$key});
        $all++;
        my @path = split/,/,$bg_ref->{$key};
        for my $path(@path){
            $stat{$path}{num}++;
            push @{$stat{$path}{id}},$key;
        }
    }
    
    return ($all,%stat);
}

sub ReadIdmap{
    my $file = shift;
    
    my %hash;
    return %hash if($file eq 'none');

    open my $file_fh,$file or die "$!:$file";
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/;
        next if(@aa < 2);
        $hash{$aa[0]} = $aa[1];
#        push @{$hash{$aa[1]}},$aa[0];
    }
    close $file_fh;

    return %hash;
}

sub ReadDb{
    my $file = shift;

    my %hash;
#    for my $type(qw/Reactome.txt ReactomeReactions.txt Reactome_All_Levels.txt Reactome_PE_All_Levels.txt Reactome_PE_Pathway.txt Reactome_PE_Reactions.txt/){
#        my $file = "$dir/${name}2$type";
    open my $in_fh,"$file" or die $!;
    while(<$in_fh>){
        my @aa = split/\t/,$_;
        $hash{$aa[0]} = $aa[1]
    }
    close $in_fh;

    return %hash;
}

sub ReadVersion{
    my $opts_ref = shift;
    my $run_type = shift;

    my $version;
    my $file = "none";
    if($run_type eq 'ID'){
        $version = $opts_ref->{version};
    }elsif($run_type eq 'PATH'){
        $file = "$opts{i}";
    }else{
        $file = "$opts{bg}";
    }
    if($file ne 'none'){
        open my $in_fh,"$file" or die "$!:$file\n";
        my $head = <$in_fh>;
        close $in_fh;

        if($head =~ /## version (\S+)/){
            $version = $1;
        }else{
            $version = "v1";
        }
    }

    return $version;
}
