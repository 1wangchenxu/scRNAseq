#! /usr/bin/perl

use strict;
use warnings;
use autodie;
use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;
use Cwd qw/abs_path/;

my ( $verbose );
GetOptions( "verbose|v" => \$verbose );
my ( $indir, $contrasts, $outdir ) = @ARGV;
die "perl $0 <indir> <contrast> <outdir>\n" unless @ARGV == 3;

$indir = abs_path($indir);
$outdir //= '.';
$outdir = abs_path($outdir);


my @types = qw/kegg go do reactome/;


for my $contrast ( split /,/, $contrasts ) {
		open IN, "$indir/$contrast/enrich.list";
		while (my $name = <IN>){
				chomp($name);
				for my $type ( @types ) {
						my $result = "$indir/$contrast/$name/$contrast.$type.Gsea";
						if ( -d $result ) {
								my $cmd = "mkdir -p '$outdir/$contrast/$name'; ln -sf '$result' '$result.xls' '$outdir/$contrast/$name'";
								print join( "\t", $contrast, $name, $type, "$contrast/$name/$contrast.$type.Gsea" ), "\n";
								system $cmd;
						} else {
#								warn "No Gsea [$type] results in [$indir/$contrast/$name]\n" if $verbose;
						}
				}
		}
		close IN;
}

