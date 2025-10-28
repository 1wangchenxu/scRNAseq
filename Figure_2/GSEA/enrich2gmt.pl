#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: enrich2gmt.pl
#
#        USAGE: ./enrich2gmt.pl  
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
#      CREATED: 12/11/2019 02:43:10 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use FindBin qw/$Bin $Script/;
#use File::Basename qw/basename dirname/;
#use Getopt::Long;
#use lib $Bin;
#use Data::Dumper;
use utf8;

my $usage =<<USAGE;

    Desc: use for xxx

    perl $0 <> 



    e.g: perl $0 <enrich_ref> <outdir> <type> <run>

USAGE

die "$usage" if(@ARGV<1); 
my $annot2gmt = "perl $Bin/annot2gmt.pl";
my $pre_symbol = "perl $Bin/pre_symbol.pl";

my $enrich_ref = shift;
my $outdir = shift;
my $type = shift;
my $run = shift;

mkdir $outdir if(!-e $outdir);

$type //= "id";
$run //= 1;

my @key1 = qw/kegg      go      do       reactome/;
my @key2 = qw/.path.xls .annot  .do.xls  .reactome.xls/;

my $symbol_info = "";
if(-e "$enrich_ref.gen2sym"){
    Excu("$pre_symbol $enrich_ref.gen2sym $outdir/all.chip");
    $annot2gmt = "$annot2gmt -symbol $enrich_ref.gen2sym" if($type eq 'symbol');
}

for my $i(0..$#key1){
    next if(!-e "$enrich_ref$key2[$i]");
    Excu("$annot2gmt -$key1[$i] $enrich_ref$key2[$i] -outfile $outdir/$key1[$i].gmt");
}

if(-e "$enrich_ref.reactome"){
    Excu("perl $Bin/reactome_enrich.pl PATH -i $enrich_ref.reactome -outpfx $outdir/all");
    Excu("$annot2gmt -reactome $outdir/all.reactome.xls -outfile $outdir/reactome.gmt");
}

if(-e "$enrich_ref.do"){
    Excu("perl $Bin/do_enrich.pl PATH -i $enrich_ref.do -outpfx $outdir/all");
    Excu("$annot2gmt -reactome $outdir/all.do.xls -outfile $outdir/do.gmt");
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
#    my $run = 1;

    system($cmd) if($run == 1);
    return;
}

