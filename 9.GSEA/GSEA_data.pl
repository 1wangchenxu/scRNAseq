#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: GSEA_data.pl
#
#        USAGE: ./GSEA_data.pl  
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
#      CREATED: 04/01/2020 03:02:09 PM
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



    e.g: perl $0 <edb> <id>

USAGE

die "$usage" if(@ARGV<1); 

my $gsea_edb = shift;
my $need_id = shift;
my $outfile = shift;

my $edb = "$gsea_edb/results.edb";
die "no [results.edb] in [$gsea_edb]" if(!-e "$gsea_edb/results.edb");
my @rnks = glob("$gsea_edb/../ranked_gene_list*xls");

print "$edb\n$rnks[0]\n";

my %rnk = ReadRank($rnks[0]);

open EDB,$edb or die $!;
my $head = join("\t","ID","NAME","PROBE","DESCRIPTION<br>(from dataset)","GENE SYMBOL","GENE_TITLE","RANK IN GENE LIST","RANK METRIC SCOR","RUNNING ES","CORE ENRICHMENT")."\n";

open OUT,">$outfile" or die $!;
print OUT "$head";

while(<EDB>){
    chomp;
    if(/<DTG/){
        my %info = GetInfo($_);
        my ($id) = $info{GENESET} =~ /gene_sets.gmt#(\S+)/;
        if($need_id eq 'all' || $id eq $need_id){
            my @gene_rnks = split/\s+/,$info{HIT_INDICES};
            my @rnd_es = split/\s+/,$info{RND_ES};
            my @es_profile = split/\s+/,$info{ES_PROFILE};
            my ($last_x,$last_y) = (0,0);
            for my $i(0..$#gene_rnks){
                my $rnk = $gene_rnks[$i];
                my $gene_id = $rnk{$gene_rnks[$i]};
                my $score = $rnk{score}{$rnk};
                my $enrichment = $score > 0 ? "Yes" : "No";
                my $es_profile_exp = CalExpectES($last_x,$last_y,$rnk);
                print OUT join("\t",$id,"row_$i\_0",$gene_id,$gene_id,$rnk{symbol}{$rnk},$rnk{title}{$rnk},$gene_rnks[$i],$score,$es_profile_exp,$enrichment)."\n";
                print OUT join("\t",$id,"row_$i",$gene_id,$gene_id,$rnk{symbol}{$rnk},$rnk{title}{$rnk},$gene_rnks[$i],$score,$es_profile[$i],$enrichment)."\n";
                ($last_x,$last_y) = ($rnk,$es_profile[$i]);
            }
            my $num = @gene_rnks;
#            print OUT join("\t","row_$num","add0","add0","add0","add0","$rnk{max}","0","0","No")."\n";
        }
    }
}
close EDB;
close OUT;# if($need_id ne 'all');

sub CalExpectES{
    my $last_x = shift;
    my $last_y = shift;
    my $x = shift;

    my $y = $last_y + ($x-$last_x)*-0.00005;

    return $y;
}

sub GetInfo{
    my $line = shift;

    my %info;
    while($line =~ s/(\S+)=\"(.*?)\"//){
        $info{$1} = $2;
    }

    return %info;
}

sub ReadRank{
    my $file = shift;

    open IN,$file or die $!;
    my %rnk;
    my $i = -1;
    while(<IN>){
        chomp;
        my @aa = split/\t/,$_;
        next if($aa[0] eq 'NAME');
        $i++;
        $rnk{$i} = $aa[0];
        $rnk{symbol}{$i} = $aa[2];
        $rnk{score}{$i} = $aa[4];
        $rnk{title}{$i} = $aa[3];
    }
    close IN;
    $rnk{max} = $i;

    return %rnk;
}

