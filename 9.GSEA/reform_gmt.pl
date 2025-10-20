#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin $Script/;
use Cwd qw/abs_path/;

die "Usage: perl $0 <gmt> <idmap> <outfile>\n" if(@ARGV<1);

my $gmt = shift;
my $idmap = shift;
my $outfile = shift;

ReformGMT($gmt,$idmap,$outfile);

sub ReformGMT{
    my $gmt = shift;
    my $idmap = shift;
    my $new_gmt = shift;

    local($_);
    open IDMAP,$idmap or die "$!: $idmap\n";
    my %idmap;
    while(<IDMAP>){
        chomp;
        my @aa = split/\t/,$_;
        push @{$idmap{$aa[1]}},$aa[0];
    }
    close IDMAP;

    open GMT,$gmt or die "$!: $gmt\n";
    open OUT,">$new_gmt" or die "$!: $new_gmt\n";;
    while(<GMT>){
        chomp;
        my @aa = split/\t/,$_;
        my %gene;
        for my $i(3..$#aa){
            next if(!exists $idmap{$aa[$i]});
            for my $gene(@{$idmap{$aa[$i]}}){
                $gene{$gene}++;
            }
        }
        next if(!%gene);
        my @genes = sort keys %gene;
        print OUT join("\t",$aa[0],$aa[1],@genes)."\n";
    }
    close GMT;
    close OUT;

    return;
}
