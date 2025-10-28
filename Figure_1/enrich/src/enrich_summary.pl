#!/usr/bin/perl -w
=head1

	perl path_sta.pl -i indir -o output.xls


=cut

use strict;
use Getopt::Long;
use Carp qw(carp confess croak);

die "Usage: perl $0 -indir -name_pfx -output\n" if(@ARGV < 1);

my %opts;
GetOptions(\%opts,"map_title=s","indir=s","output=s","name_pfx=s");

if(!-e $opts{indir}){
    print STDERR "WARNINGS: No any enrich result\n";
    exit(0);
}
#my @infiles = `ls -v $opts{indir}/*$opts{name_pfx}`;
my @infiles = `find $opts{indir} -name "*$opts{name_pfx}"`;# `find $indir -name "$file_pat"`;
if(@infiles == 0 && -e $opts{indir}){
    print STDERR "WARNINGS: Can't find any file of [*$opts{name_pfx}] in dir [$opts{indir}]\n";
    exit(0);
}
Error("Can't find any file of [*$opts{name_pfx}] in dir [$opts{indir}]") if(@infiles == 0);

my %path_info;
for my $file(sort @infiles){
    chomp($file);
    my $name = $file;
    $name =~ s/.*\///;
    $name =~ s/$opts{name_pfx}$//;
    push @{$path_info{name_list}},$name;
    ReadInfo($file,\%path_info,$name);
}

Output($opts{output},\%path_info);


sub ReadInfo{
    my $file = shift;
    my $info_ref = shift;
    my $name = shift;

    #KEGG_A_class   KEGG_B_class    Pathway profile0 (11)   All (5632)  Pvalue  Qvalue  Pathway ID  Genes   K_IDs
    open FILE,$file or die $!;
    my @title_col;
    my @info_col;
    my %hash;
    while(<FILE>){
        chomp;
        $_ =~ s/\t$/\t--/ if(/\t$/);
        my @aa = split/\t/,$_;
        if($. == 1){
            $hash{type} = CheckFile(@aa);
            if($hash{type} eq "kegg"){
                # KEGG_A_class   KEGG_B_class    Pathway profile0 (11)   All (5632)  Pvalue  Qvalue  Pathway ID  Genes   K_IDs
                @title_col = (2,7,0,1);
                @info_col = (3,5,6);
                $info_ref->{$name}{all} = GetNum($aa[3]);
                $info_ref->{head} = "Pathway\tPathway_ID\tKEGG_A_class\tKEGG_B_class" if(!exists $info_ref->{head});
                $info_ref->{head} .= "\t$name($info_ref->{$name}{all})\t$name\_Pvalue\t$name\_Qvalue";
                $info_ref->{all_num} = 1;
                $info_ref->{default} = "0\t1\t1"; ## num pvalue qvalue
            }elsif($hash{type} eq 'kegg_old'){
                # #Pathway AN-VS-AR (765) All-Unigene (6784) Pvalue Qvalue Pathway ID Genes KOs 
                @title_col = (0,5);
                @info_col = (1);
                $info_ref->{$name}{all} = GetNum($aa[1]);
                $info_ref->{head} = "Pathway\tPathway_ID" if(!exists $info_ref->{head});
                $info_ref->{head} .= "\t$name($info_ref->{$name}{all})";
                $info_ref->{all_num} = 1;
                $info_ref->{default} = "0\t1\t1"; ## num pvalue qvalue
            }elsif($hash{type} eq 'go_old' || $hash{type} eq 'godiff_old'){
                # Ontology Class number_of_AN-VS-AR_up number_of_AN-VS-AR_down genes_of_AN-VS-AR_up genes_of_AN-VS-AR_down
                @title_col = (0,1);
                $info_ref->{head} = "Ontology\tClass" if(!exists $info_ref->{head});
                if($hash{type} eq 'godiff_old'){
                    @info_col = (2,3,4,5);                    
                    $info_ref->{default} = join("\t",0,"--",0,"--"); ## up up_glist down down_glist
                    $info_ref->{head} .= "\tGene Number($name\_up)\tGene Number($name\_down)\tGene List($name\_up)\tGene List($name\_down)";
                }else{
                    @info_col = (2,3);
                    $info_ref->{default} = join("\t",0,"--"); ## gene_number  gene_list
                    $info_ref->{head} .= "\tGene Number($name)\tGene List($name)";
                }
            }elsif($hash{type} =~ 'go'){
                @title_col = (0,1,2,3);
                $info_ref->{head} = "GO ID (level1)\tGO Term (level1)\tGO ID (level2)\tGO Term (level2)" if(!exists $info_ref->{head});
                if($hash{type} eq 'godiff'){
                    @info_col = (4,5,6,7);
                    $info_ref->{default} = join("\t",0,0,"--","--");
                    $info_ref->{head} .= "\tGene Number($name\_up)\tGene Number($name\_down)\tGene List($name\_up)\tGene List($name\_down)";
                }else{
                    @info_col = (4,5);
                    $info_ref->{default} = join("\t",0,"--");
                    $info_ref->{head} .= "\tGene Number($name)\tGene List($name)";
                }
            }elsif($hash{type} =~ /reactome/i){
                ## Reactome Reactome_Name   URL A-vs-B (436)    All (10734) pvalue  qvalue  Genes_of_all 
                @title_col = (0,1,2); ##
                @info_col = (3,5,6);
                $info_ref->{$name}{all} = GetNum($aa[3]);
                $info_ref->{head} = "Reactome\tReactome_Name\tURL" if(!exists $info_ref->{head});
                $info_ref->{head} .= "\t$name($info_ref->{$name}{all})\t$name\_Pvalue\t$name\_Qvalue";
                $info_ref->{all_num} = 1;
                $info_ref->{default} = "0\t1\t1"; ## number  pvalue qvalue 
            }elsif($hash{type} =~ /do/i){
                # doid  DO_Name URL A-vs-B (436)    All (6229)  pvalue  qvalue  Genes 
                @title_col = (0,1,2);
                @info_col = (3,5,6);
                $info_ref->{$name}{all} = GetNum($aa[3]);
                $info_ref->{head} = "doid\tDO_Name\tURL" if(!exists $info_ref->{head});
                $info_ref->{head} .= "\t$name($info_ref->{$name}{all})\t$name\_Pvalue\t$name\_Qvalue";
                $info_ref->{all_num} = 1;
                $info_ref->{default} = "0\t1\t1"; ## number  pvalue qvalue 
            }
        }else{
            my $key = join("\t",@aa[@title_col]);
            my $num = scalar(@aa);
            Error("col [$num] <= need max col [$info_col[-1]], [@aa] [@info_col] in [$file]") if(@aa <= $info_col[-1]);
            $info_ref->{$name}{$key} = join("\t",@aa[@info_col]);
            $info_ref->{path_list}{$key} = 1;
        }
    }
    close FILE;

    return;
}

sub InitGOInfo{
    my $info_ref = shift;
    my $head = shift;
    my $type = shift;
    my $name = shift;

    $info_ref->{head} = $head if(!exists $info_ref->{head});
    if($type eq 'godiff'){
        $info_ref->{head} .= "\tGene Number($name\_up)\tGene Number($name\_down)\tGene List($name\_up)    \tGene List($name\_down)";
    }

}

sub CheckFile{
    my @head = @_;

    #KEGG_A_class   KEGG_B_class    Pathway profile0 (11)   All (5632)  Pvalue  Qvalue  Pathway ID  Genes   K_IDs
    if($head[0] eq 'KEGG_A_class' && $head[5] eq 'Pvalue' && ($head[9] eq 'K_IDs' || $head[9] eq 'C_IDs')){
        return 'kegg';
    }
    # #Pathway AN-VS-AR (765) All-Unigene (6784) Pvalue Qvalue Pathway ID Genes KOs
    if($head[0] eq '#Pathway'){
        return 'kegg_old';
    }
    # GO ID (level1)    GO Term (level1)    GO ID (level2)  GO Term (level2)    number_of_all (All) Genes_of_all
    if($head[0] eq 'GO ID (level1)' && $head[3] eq 'GO Term (level2)'){
        if(@head == 6){
            return "go";
        }elsif(@head == 8){
            return "godiff";
        }
    }
    # Ontology Class number_of_AN-VS-AR_up number_of_AN-VS-AR_down genes_of_AN-VS-AR_up genes_of_AN-VS-AR_down
    if($head[0] eq 'Ontology'){
        if(@head == 6){
            return 'godiff_old';
        }else{
            return 'go_old';
        }
    }
    ## Reactome Reactome_Name   URL A-vs-B (436)    All (10734) pvalue  qvalue  Genes 
    if($head[0] eq 'Reactome' && $head[5] eq 'pvalue'){
        return 'reactome';
    }

    ## doid DO_Name URL A-vs-B (436)    All (6229)  pvalue  qvalue  Genes 
    if($head[0] eq 'doid' && $head[5] eq 'pvalue'){ 
        return 'do';
    }

    Error("Unknown type [@head]");

    return "none";
}

sub GetNum{
    my $col = shift;

    $col =~ s/.*\(//;
    $col =~ s/\).*//;
#    print "$col\n";

    return $col;
}

sub Output{
    my $outfile = shift;
    my $info_ref = shift;

    my @name_list = @{$info_ref->{name_list}};

    open OUT,">$outfile" or die $!;
    print OUT "$info_ref->{head}\n";
#    print "@name_list\n";
#    if(exists $info_ref->{all_num}){
#        print OUT join("\t",$info_ref->{head},map{"$_($info_ref->{$_}{all})"}@name_list)."\n";
#    }else{
#        print OUT join("\t",$info_ref->{head},map{my $a="Gene Number($_)","Gene List($_)"}@name_list)."\n";
#    }
    for my $path(sort keys %{$info_ref->{path_list}}){
        print OUT join("\t",$path,map{$info_ref->{$_}{$path} //= $info_ref->{default};$info_ref->{$_}{$path}} @name_list);
        print OUT "\n";
    }
    close OUT;

    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}

