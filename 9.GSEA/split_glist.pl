#! /usr/bin/perl

use strict;
use warnings;
use autodie;
use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;


my ( $c_col, $g_col, $l_col ) = (1, 2, 10);
GetOptions( 
				"c_col|c=i" => \$c_col,
				"g_col|g=i" => \$g_col,
				"l_col|l=i" => \$l_col,
		  );
die "perl $0 ( -c_col <int> -g_col <int> -l_col <int> ) <file> <outdir>\n" unless @ARGV == 2;
my ( $file, $outdir ) = @ARGV;

my @order;
my %glist;
open IN, "<", $file;
<IN>;
while (<IN>){
		chomp;
		my ( $cls, $id, $fc ) = (split /\t/)[$c_col - 1, $g_col - 1, $l_col - 1];
		$cls =~ s#/#_#g;
		$cls =~ s/[ \(\)]/_/g;
		if ( ! exists $glist{$cls} ) {
				push @order, $cls;
		}
		$glist{$cls}{$id} = $fc;
}
close IN;

`mkdir $outdir` unless -d $outdir;
open my $fh_list, ">", "$outdir/enrich.list";
for my $cls ( @order ) {
		open my $fh, '>', "$outdir/Cluster_$cls.glist";
		for my $id ( sort keys %{$glist{$cls}} ) {
				print $fh join( "\t", $id, $glist{$cls}{$id} ), "\n";
		}
		close $fh;

		print $fh_list "Cluster_$cls\n";
}
close $fh_list;

