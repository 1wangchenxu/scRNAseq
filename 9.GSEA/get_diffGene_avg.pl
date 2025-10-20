#! /usr/bin/perl

use strict;
use warnings;
use autodie;
#use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;

die "perl $0 <list_file> <indir> <outdir>\n" unless @ARGV == 3;
my ( $file, $indir, $outdir ) = @ARGV;

my %list;
open OUT, '>', "$outdir/enrich.list";
open IN, $file;
<IN>;
while (<IN>){
		chomp;
		my ( $contrast, $cluster, $gene ) = split /\t/;
		print OUT "$cluster\n" unless exists $list{$cluster};
		$list{$cluster}{$gene} = 1;
}
close IN;
close OUT;

for my $cluster ( sort keys %list ) {
		open my $in,  '<', "$indir/AllGene.avg_exp.$cluster.xls";
		open my $out, '>', "$outdir/DifferGene.avg_exp.$cluster.xls";
		print { $out } scalar <$in>;
		while (<$in>){
				my ( $id ) = split /\t/;
				if ( exists $list{$cluster}{$id} ) {
						print { $out } $_;
				}
		}
		close $in;
		close $out;
}

