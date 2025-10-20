#! /usr/bin/perl

use strict;
use warnings;
use autodie;
#use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;

die "perl $0 <indir> <names> ...\n" unless @ARGV > 1;
my ( $indir, @names ) = @ARGV;

my @types = ( "GO", "KO", "DO", "Reactome" );
my %types_name = ( "GO" => "GO", "KO" => "KEGG", "DO" => "DO", "Reactome" => "Reactome");

for my $name ( @names ) {
		open IN, "$indir/$name/enrich.list";
		while (<IN>){
				chomp;
				my ( $cls ) = split /\t/;
				for my $type ( @types ) {
						my @files = glob "$indir/$name/$type/$cls*";
						if ( @files ) {
								print join ( "\t", $name, $cls, $types_name{$type} ), "\n";
						}
				}
		}
		close IN;
}

