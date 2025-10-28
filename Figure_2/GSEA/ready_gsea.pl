#! /usr/bin/perl

use strict;
use warnings;
use autodie;
#use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;

die "perl $0 <file> <refdir> <outdir>\n" unless @ARGV >= 2;
my ( $file, $refdir, $outdir ) = @ARGV;
$outdir //= '.';

my ( $ctrl, $treat ) = $file =~ /DifferMarker\.(.+)-vs-(.+)\.(diff|all)(\.annot)?\.xls/;

my $perl = "/Bio/bin/perl";
my $run_GSEA = "/Bio/Bin/pipeline//System_Programs/Small_pipe/v2.0/GSEA/run_GSEA.pl";

open OUT, '>', "$outdir/group_diff.list";
print OUT join( "\t", $ctrl, $treat ), "\n";
close OUT;

open OUT, '>', "$outdir/group.list";
print OUT join( "\t", $ctrl,  $ctrl  ), "\n";
print OUT join( "\t", $treat, $treat ), "\n";
close OUT;


open OUT, '>', "$outdir/enrich.list";
my %fh;
open IN, $file;
my $header = join ( "\t", (split /\t/, scalar <IN>)[1,3,4] );
while (<IN>){
		chomp;
		my ( $cls, $id, undef, $ctr, $tr ) = split /\t/;
		$cls =~ s#/#_#g;
		$cls =~ s/[ \(\)]/_/g;
		if ( ! exists $fh{$cls} ) {				
				print OUT $cls, "\n";
				open $fh{$cls}, ">", "$outdir/$cls.xls";
				print { $fh{$cls} } $header, "\n";
		}
		print { $fh{$cls} } join( "\t", $id, $ctr, $tr ), "\n";
}
close IN;

my $run;
for my $type ( qw/go ko reactome do/ ){
		$run .= " -run_$type ";
		my $file = $type eq 'ko' ? "$refdir/kegg.gmt" : "$refdir/$type.gmt";
		$run .= -s $file ? 1 : 0;
}

`$perl -ne 'chomp; print "$perl $run_GSEA $run -run_sge yes -compare $outdir/group_diff.list -group $outdir/group.list -enrich_gmt $refdir -outdir $outdir/\$_ -exp $outdir/\$_.xls\\n";' $outdir/enrich.list > $outdir/ALLSTEP.sh`;


