#! /usr/bin/perl

use strict;
use warnings;
use autodie;
use Getopt::Long;
#use Bio::SeqIO;
#use List::Util qw/sum min max/;
#use File::Basename qw/basename dirname/;
#use FindBin qw/$Bin $Script/;


my ( $bc, $ac, $type ) = ( 1, 1, "append" );
GetOptions(
				"bc|b=i" => \$bc,
				"ac|a=i" => \$ac,
				"type|t=s" => \$type,
		  );
die "perl $0 (-bc <int> -ac <int> -type [append|insert]) <base_file> <add_file>\n\t-bc, -ac : also can be negative int\n" unless @ARGV == 2;
my ( $b_file, $a_file ) = @ARGV;

open IN, $a_file;
my %a_head = get_header( $ac, scalar <IN> );
my %annot = map { chomp; my @data = split /\t/; my $id = splice @data, $ac < 0 ? $ac : $ac - 1, 1; ( $id, join( "\t", @data ) ) } <IN>;
close IN;

open IN, $b_file;
chomp(my @b_head = split /\t/, <IN>);
$bc = scalar(@b_head) + $bc + 1 if $bc < 0;
my $add_pos = $type eq "append" ? scalar(@b_head)
			: $type eq "insert" ? $bc
			: die;
splice @b_head, $add_pos, 0, $a_head{header};
print join( "\t", @b_head ), "\n";
while (<IN>){
		chomp;
		my @data = split /\t/;
		splice @data, $add_pos, 0, exists $annot{$data[$bc - 1]} ? $annot{$data[$bc - 1]} : $a_head{null};
		print join( "\t", @data ), "\n";
}
close IN;

sub get_header {
		my ( $col, $string ) = @_;
		chomp(my @data = split /\t/, $string);
		splice @data, $col < 0 ? $col : $col - 1, 1;
		my %head = ( "header" => join( "\t", @data ),
					"null" => join( "\t", ( '-' ) x scalar @data )
				   );
		return %head;
}
