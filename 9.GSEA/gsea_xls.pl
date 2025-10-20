#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gsea_xls.pl
#
#        USAGE: ./gsea_xls.pl  
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
#      CREATED: 01/22/2019 03:52:05 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use FindBin qw/$Bin $Script/;
#use File::Basename qw/basename dirname/;
#use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 



    e.g: perl $0 <> <>

USAGE

die "$usage" if(@ARGV<1); 

my $gmt = shift;
my $symbol = shift;
my $outfile = shift;
my @files = @ARGV;

my %symbol = ReadSymbol($symbol) if($symbol ne 'none');
open GMT,$gmt or die $!;
my %gmt;
my $annot_flag = (-e "$gmt.annot.list") ? 1 : 0;
my %annot = ReadAnnot("$gmt.annot.list") if($annot_flag);
my %genes;
while(<GMT>){
    chomp;
    my @aa = split/\t/,$_;
    my $id = $aa[0];
    $aa[0] = uc($aa[0]);
    $gmt{$aa[0]} = $aa[1];
    $gmt{$aa[0]} .= "\t$annot{$id}" if($annot_flag);
    my @genes;
    for my $gene(@aa[2..$#aa]){
        if($symbol ne 'none'){
            $symbol{$gene} //= "--";
            push @genes,"$gene($symbol{$gene})";
        }else{
            push @genes,$gene;
        }
    }
    $genes{$aa[0]} = join(";",@genes);
}
close GMT;
$genes{id} = "Genes";
$gmt{id} = "name";
$gmt{id} .= "\t$annot{head}" if($annot_flag);
#$gmt{id} .= "\tGenes";

open OUT,">$outfile" or die $!;
for my $file(@files){
    my @cols;
    open FILE,$file or die $!;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        if($. == 1){
            @cols = GetCol(@aa);
            $aa[0] = "id";
            next if($file ne $files[0]);
        }
        for my $i(@cols){
            if($aa[$i] eq ''){
                $aa[$i] = "0.01";
#                print "$aa[0]\t$aa[$i]\t$cols[1]\t$aa[$cols[1]]\n";
                $aa[$i] = -$aa[$i] if($aa[$cols[1]] < 0);
            }

        }
        print OUT join("\t",$aa[0],$gmt{$aa[0]},@aa[@cols],$genes{$aa[0]})."\n";
#        print OUT join("\t",@aa[@cols])."\n";
    }
    close FILE;
}
close OUT;

sub GetCol{
    my @all = @_;

    my %need = qw/SIZE 1 ES 1 NES 1 /;
    $need{"NOM p-val"} = 1;
    $need{"FDR q-val"} = 1;
    my @out;
    for my $i(0..$#all){
        if(exists $need{$all[$i]}){
            push @out,$i;
        }
    }

    return @out;
}

sub ReadAnnot{
    my $infile = shift;

    open FILE,$infile or die $!;
    my %hash;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_,2;
        $hash{$aa[0]} = $aa[1];
    }
    close FILE;

    return %hash;
}

sub ReadSymbol{
    my $infile = shift;
    
    open FILE,$infile or die $!;
    my %hash;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        $hash{$aa[0]} = $aa[1];
    }
    close FILE;

    return %hash;
}

