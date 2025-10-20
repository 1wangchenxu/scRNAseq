#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: annot2gmt.pl
#
#        USAGE: ./annot2gmt.pl  
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
#      CREATED: 01/18/2019 09:51:02 AM
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
use  Support_Program;
use lib "$RealBin";
use GOVersion;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 
        -kegg 
        -go 
	-do
	-reactome
        -outfile

    e.g: perl $0 <> <>

USAGE

die "$usage" if(@ARGV<1); 

my %opts;
GetOptions(\%opts,"kegg=s","go=s","do=s","reactome=s","outfile=s","symbol=s");

my %apps;
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");

#Support_Program("$RealBin/../",\%apps,1);

my $out;
my $annot;

my %g_symbol_info = ReadSymbol($opts{symbol}) if(exists $opts{symbol});

open OUT,">$opts{outfile}" or die $!;
if(exists $opts{kegg}){
    ($out,$annot) = GetKegg($opts{kegg});
}

if(exists $opts{go}){
    ($out,$annot) = GetGO($opts{go});
}

if(exists $opts{do}){
    ($out,$annot) = GetDOReac($opts{do});
}

if(exists $opts{reactome}){
    ($out,$annot) = GetDOReac($opts{reactome})
}

print OUT "$out";
close OUT;
if(defined $annot && $annot ne ''){
    open ANNOT,">$opts{outfile}.annot.list" or die $!;
    print ANNOT "$annot";
    close ANNOT;
}

sub ReadSymbol{
    my $symbol = shift;

    open SB,$symbol or die $!;
    my %hash;
    while(<SB>){
        chomp;
        my @aa = split/\t/,$_;
        next if(@aa < 2);
        $hash{$aa[0]} = uc($aa[1]);
    }
    close SB;

    return %hash;
}

sub GetKegg{
    my $file = shift;

    open my $file_fh,$file or die $!;
#    <$file_fh>;
    # KEGG_A_class  KEGG_B_class    Pathway Count (6722)    Pathway ID  Genes   K_IDs
    my $out = "";
    my $annot = "";
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/,$_;
        my $genes = $aa[5];
        if($. == 1){
            $annot .= "head\t$aa[0]\t$aa[1]\n";
            next;
        }
        $genes = Gene2Symbol($genes);
        next if($genes eq '');
        $out .= "$aa[4]\t$aa[2]\t$genes\n";
        $annot .= "$aa[4]\t$aa[0]\t$aa[1]\n";
    }
    close $file_fh;

    return ($out,$annot);
}

sub GetDOReac{
    my $file = shift;
    open my $file_fh,$file or die $!;

    my $out = "";
    my $annot = "";
    while(<$file_fh>){
        chomp;
        my @aa = split/\t/,$_;
        my $genes = $aa[4];
        if($. == 1){
            $annot .= "head\t$aa[2]\n";
            next;
        }
        $genes = Gene2Symbol($genes);
        next if($genes eq '');
        $out .= "$aa[0]\t$aa[1]\t$genes\n";
        $annot .= "$aa[0]\t$aa[2]\n";
    }
    close $file_fh;
    return ($out,$annot);
}
sub GetGO{
    my $file = shift;

#    my $db = $apps{GOannot};#"/Bio/User/hongshimiao/bin/Small_pipe/GO/src/goannot.3.3.1";
#    my $goterm = $apps{GOterm};#"/Bio/User/hongshimiao/bin/Small_pipe/GO/src/goterm.3.3.1";
    my ($Rscript,%go_info) = GOVersion($file,\%apps);
    my $db = $go_info{goannot};
    my $goterm = $go_info{godesc};
    my $second_level = $go_info{sl};

    my %slgo = readSL($second_level,1,"0,2");
    my %slgo_all;
#    print %go_info;die;
    open DB,$db or die $!;
    my %db;
    while(<DB>){
        chomp;
        my @aa = split;
        my @bb = split/,/,$aa[1];
        @{$db{$aa[0]}} = @bb;
        for my $go(@bb){
            if(exists $slgo{$go}){
                $slgo_all{$aa[0]} = $slgo{$go};
            }
        }
    }
    close DB;

    open TERM,$goterm or die $!;
    my %goterm;
    while(<TERM>){
        chomp;
        my @aa = split/\t/,$_;
        $goterm{$aa[0]} = $aa[2];
    }
    close TERM;

    open my $file_fh,$file or die $!;
    my %go;
    while(<$file_fh>){
        chomp;
        my @aa = split;
        next if(@aa < 2);
        next if(!exists $db{$aa[1]});
        for my $go(@{$db{$aa[1]}}){
            $go{$go}{$aa[0]} = 1;
        }
    }
    close $file_fh;

#    print "$second_level\n";

    my $out = "";
    my $annot = "head\tGO Term(level1)\tGO Term(level2)\n";
    my %cc = ("biological_process" => "Biological Process",
              "cellular_component" => "Cellular Component", 
              "molecular_function" => "Molecular Function"
            );
    for my $go(sort keys %go){
        next if(!exists $goterm{$go});
        my $genes = join(";",sort keys %{$go{$go}});
        $genes = Gene2Symbol($genes);
        next if($genes eq '');
        $out .= join("\t",$go,$goterm{$go},$genes)."\n";

        $slgo_all{$go} = "$cc{$goterm{$go}}\t--" if(exists $cc{$goterm{$go}});
        $annot .= join("\t",$go,$slgo_all{$go})."\n";
    }

    return ($out,$annot);
}

sub readSL{
    my $sl = shift;
    my $key_col = shift; ## 0 base
    my $need_col = shift;

    my @need_cols = split/,/,$need_col;


    my %slgo;
    open SL, $sl or die "$!:$sl\n";
    while(<SL>){
        chomp;
        my @tmp = split /\t/;
        $slgo{$tmp[$key_col]} = join("\t",@tmp[@need_cols]);
#        $slgo{$tmp[1]} = "$tmp[3]\t$tmp[0]\t$tmp[1]\t$tmp[2]";
    }
    close SL;

    return %slgo;
}

sub Gene2Symbol{
    my $gene_list = shift;
    my $symbol_ref = \%g_symbol_info;

    my @genes = split/;/,$gene_list;
    

    my @out;
    my %e;
    if(exists $opts{symbol}){
        for my $gene(@genes){
            next if(!exists $symbol_ref->{$gene});
            my $symbol = $symbol_ref->{$gene};
            next if($symbol eq '--');
            next if(exists $e{$symbol});
            $e{$symbol} = 1;
            push @out,$symbol;
        }
    }else{
        @out = @genes;
    }

    return join("\t",@out);
}
