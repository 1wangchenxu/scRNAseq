#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: ReplotGSEA.pl
#
#        USAGE: ./ReplotGSEA.pl  
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
#      CREATED: 06/10/2019 06:13:43 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use FindBin qw/$Bin $Script/;
#use File::Basename qw/basename dirname/;
use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;
use Carp qw(carp confess croak);

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 



    e.g: perl $0 <indir> <list> <outdir> <R_plot> <top 20/all> <run or not>
                -id

USAGE

die "$usage" if(@ARGV<1); 

my %opts = (slope=>0.05);
GetOptions(\%opts,"id=s","slope=s");

my $indir = shift;
my $list = shift;
my $outdir = shift;
my $R_plot = shift;
my $top = shift;
my $run = shift;

mkdir $outdir if(!-e $outdir);

$top //= 20;
$run //= 1;

open LIST,$list or die $!;
my $pos_top = 0;
my $neg_top = 0;
my ($id_col,$desc_col,$nes_col);
my %png_name;
my @rename_files = ("$indir/pos_snapshot.html","$indir/neg_snapshot.html");
while(<LIST>){
    chomp;
    my @aa = split/\t/,$_;
    if($. == 1){
        $id_col = GetCol(\@aa,"NAME",'id');
        $desc_col = GetCol(\@aa,"Description",'name');
        $nes_col = GetCol(\@aa,"NES");
        next;
    }
    my $name = $aa[$id_col];
    if(exists $opts{id}){
        next if($name ne $opts{id});
    }else{
        next if($aa[$nes_col] eq '-');
        if($aa[$nes_col] > 0){
            $pos_top++;
            next if($top ne 'all' && $pos_top > $top);
        }else{
            $neg_top++;
            next if($top ne 'all' && $neg_top > $top);
        }
    }
    if(length($aa[$desc_col]) > 100){
        $aa[$desc_col] = substr($aa[$desc_col],0,100);
    }
    push @rename_files,"$indir/$name.html";
    $name =~ s/:/_/;
    $aa[$desc_col] =~ s/\///g;
    my $outpfx = "$outdir/$aa[$id_col]_$aa[$desc_col].enplot";
    $outpfx =~ s/\s/_/g;
    $outpfx =~ s/[\"\'\(\)]//g;
    $outpfx =~ s/://g;
    $outpfx =~ s/,//g;
    $png_name{$name} = "$outpfx.png";
    $png_name{$name} =~ s/.*\///;
    $aa[$desc_col] =~ s/\"/\\"/g;
    Excu("perl $R_plot ReplotGSEA -infile $indir -name $aa[$id_col] -title \"Enrichment plot:$aa[$id_col] $aa[$desc_col]\" -slope $opts{slope} -outpfx $outpfx");

}
close LIST;
exit(0) if(exists $opts{id});

## re html 
for my $file(@rename_files){
    my $name = $file;
    $name =~ s/.*\///;
    open FILE,$file or die "$!:$file";
    my $out;
    while(<FILE>){
        chomp;
        while(/(enplot_(\S+)_\d+.png)/){
            my ($old_name,$term) = ($1,$2);
            die "No png name of [$term] [$old_name] [$file]\n" if(!exists $png_name{$term});
            s/$old_name/$png_name{$term}/g;
        }
        $out .= "$_\n";
    }
    close FILE;
    open OUT,">$outdir/$name" or die "$!:$outdir/$name";
    print OUT "$out";
    close OUT;
}

## rename GO: => GO_
my @files = glob("$indir/gsea_report_for_*html");
push @files,("$indir/pos_snapshot.html","$indir/neg_snapshot.html");
for my $file(@files){
    print "$file\n";
    my $name = $file;
    $name =~ s/.*\///;
    open FILE,$file or die "$!:$file\n";
    my $out;
    while(<FILE>){
        chomp;
        while(/(GO:(\d+.html))/){
            my $old = $1;
            my $new = "GO_$2";
#            print "$old\t$new\n";
            s/$old/$new/;
        }
        while(/(DOID:(\d+.html))/){
            my $old = $1;
            my $new = "DOID_$2";
            s/$old/$new/;
        }
        $out .= "$_\n";
    }
    close FILE;

    open OUT,">$outdir/$name" or die "$!:$outdir/$name";
    print OUT "$out";
    close OUT;
}

my @go_htmls = glob("$indir/GO:*.html");
my @do_htmls = glob("$indir/DOID:*.html");
if (@go_htmls > 0 ){
    Rename(\@go_htmls,$indir,$outdir,"GO");
}
if (@do_htmls > 0 ){
    Rename(\@do_htmls,$indir,$outdir,"DO");
}

sub Rename{
    my $htmls = shift;
    my $indir = shift;
    my $outdir = shift;
    my $ty = shift;

    for my $htmls(@$htmls){
        if($indir eq $outdir){
            Excu("sed \"s/png\'>/png\' width=\\\"400\\\">/\" -i $htmls");
        }else{
            my $name = $htmls;
            $name =~ s/.*\///;
            Excu("sed \"s/png\'>/png\' width=\\\"400\\\">/\" $htmls >$outdir/$name");
        }
    }

    if($ty eq "GO"){
        Excu("rename GO: GO_ $outdir/GO:*");
    }else{
        Excu("rename DOID: DOID_ $outdir/DOID:*");
    }
    return;
}


sub Excu{
    my $cmd = shift;

    print "$cmd\n";

    system($cmd) if($run == 1);
    return;
}

sub GetCol{
    my $col_arr = shift;
    my @names = @_;

    my %col_info = map{$col_arr->[$_]=>$_} 0..$#$col_arr;
    my $col;
    for my $name(@names){
        next if(!exists $col_info{$name});
        $col = $col_info{$name};
        last;
    }

    Error("No col info of [@names] in [@$col_arr]") if(!defined $col);
    return $col;
}

sub Error{

    confess "FATAL ERROR: @_\n";
    die;
}

