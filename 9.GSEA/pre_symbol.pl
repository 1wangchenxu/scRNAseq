#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: pre_symbol.pl
#
#        USAGE: ./pre_symbol.pl  
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
#      CREATED: 09/01/2019 10:26:14 AM
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



    e.g: perl $0 <symbol> <out>

USAGE

die "$usage" if(@ARGV<1); 

my $symbol = shift;
my $outfile = shift;

open IN,$symbol or die $!;
open OUT,">$outfile" or die "$!:$outfile";
print OUT "Probe Set ID\tGene Symbol\tGene Title\n";
while(<IN>){
    chomp;
    my @aa = split/\t/,$_;
    next if($aa[0] eq 'geneID');
    print OUT "$_\tNA\n";
}
close IN;
close OUT;
