#!/usr/bin/perl
use warnings;
use strict;
use File::Basename qw(basename);
use FindBin qw($Bin $RealBin);
use Getopt::Long;
use Cwd 'abs_path';
use lib "$RealBin/../";
use Support_Program;
use Carp qw(carp confess croak);

my %opts =(ud=>'all');
GetOptions(\%opts,"ud=s","op=s","g=s","k=s");

my $second_level_file = shift;
my $pfc = shift;
my $inpfx = shift;


my @afs = split/,/,$pfc;# = qw/P F C/;

## read second level
my %slgo = readSL($second_level_file);

## read diff info  id \t diff
my %gene_ud = readDiff($opts{g},$opts{ud});

open SECL, "> $opts{op}.Level2.xls" or die $!;
my $title = basename($opts{op}, ".Level2.xls");
my $hh = "GO ID (level1)\tGO Term (level1)\tGO ID (level2)\tGO Term (level2)";
if($opts{ud} eq "nodiff" || $opts{ud} eq 'all'){
    print SECL "$hh\tnumber_of_$title (All)\t$opts{k}_of_$title\n";
}else{
    print SECL "$hh\tnumber_of_${title} (up)\tnumber_of_${title} (down)\t$opts{k}_of_${title}_up\t$opts{k}_of_${title}_down\n";
}
my @xls;

foreach my $t (@afs){
    my %cc = qw/C CC  P BP  F MF/;
    my $type = $cc{$t};
    
    open PFC, "$opts{op}.$t.xls" or die $!;
    my $new_pfc = <PFC>;

    while(<PFC>){
        chomp;
        my @tmp = split /\t/;
        my @genes = split /;/, $tmp[6];
        my @up_gene;
        my @down_gene;
        my @all_genes;

        next if(!exists $slgo{$tmp[0]});
        for my $i(@genes){
            my $gene = $i;
            $gene =~ s/\(.*//;
#            $symbol{$i} //= "-";
#            $gene = "$i($symbol{$i})" if(defined $opts{symbol});
            push @all_genes,$i;
            if(!exists $gene_ud{$gene} || $gene_ud{$gene} >0){
                push @up_gene, $i;
            }else{
                push @down_gene,$i;
            }
        }
        if(exists $slgo{$tmp[0]}){
            my $up = @up_gene;
            my $down = @down_gene;
            my $up_gene = @up_gene > 0 ? join "$sep", @up_gene : '--';
            my $down_gene = @down_gene > 0 ? join "$sep", @down_gene : '--';
            if($opts{ud} eq "nodiff" || $opts{ud} eq 'all'){
                print SECL $slgo{$tmp[0]}."\t$up\t$up_gene\n";
            }elsif($opts{ud} eq "diff"){
                print SECL $slgo{$tmp[0]}."\t$up\t$down\t$up_gene\t$down_gene\n";
            }
        }
#        $tmp[6] = join(";",@all_genes);
#        $new_pfc .= join("\t",@tmp)."\n";
    }
    close PFC;
}
close SECL;

sub readSL{
    my $sl = shift;
    
    my %slgo;
    open SL, $sl or die "$!:$sl\n";
    while(<SL>){
        chomp;
        my @tmp = split /\t/;
        $slgo{$tmp[1]} = "$tmp[3]\t$tmp[0]\t$tmp[1]\t$tmp[2]";
    }
    close SL;
    
    return %slgo;
}

sub readDiff{
    my $glist = shift;
    my $type = shift;

    my %gene_ud;
    if($type eq "diff"){
        open GENE_LIST, "$glist" or die "$!:$glist";
        while(<GENE_LIST>){
            chomp;
            my @tmp = split /\t/;
            $gene_ud{$tmp[0]} = $tmp[1];
        }
        close GENE_LIST;
    }
    
    return %gene_ud;
}


sub Excu{
    my $cmd = shift;
    print "$cmd\n";

    if($opts{run} == 1){
        my $ret = system($cmd);
        if($ret != 0){
            Error("err$ret\t$cmd");
        }
    }

    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}

