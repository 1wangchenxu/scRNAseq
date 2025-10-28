#!perl
use warnings;
use strict;
use FindBin qw($RealBin);
use File::Basename qw(basename);

my $src = "$RealBin";
die "perl $0 <input dir> <go or kegg, do, reactome>\n" if @ARGV != 2;

my $indir = shift;
my $enrich_type = shift;


my $file_pat;
my %col;
my $head;
if($enrich_type eq "go" || $enrich_type eq 'GO'){
    $file_pat = "*.[PFC].xls";
    %col = qw/num 2 all 3 p 4 q 5/;
    $head = "Type\tTerm";
}elsif($enrich_type eq 'kegg' || $enrich_type eq 'KO'){
    $file_pat = "*.path.xls";
    %col = qw/num 3 all 4 p 5 q 6/;
    $head = "B_Pathway\tC_Pathway";
}elsif($enrich_type eq 'DO'){
    $file_pat = "*.do.xls";
    %col = qw/num 3 all 4 p 5 q 6/;
}elsif($enrich_type eq 'Reactome'){
    $file_pat = "*.reactome.xls";
    %col = qw/num 3 all 4 p 5 q 6/;
}
my @files = `find $indir -name "$file_pat"`;
@files = sort @files;

my %go_pfc = (".P.xls" => "Biological Process", ".F.xls"=> "Molecular Function", ".C.xls" => "Cellular Component");

my (%info, %all, %name);
my @sample_names;
my %e;
for my $file(@files){
    chomp($file);
    my $name = basename($file);
#    my ($name,$pat) = $name_tmp =~ /(\S+)({$file_pat})/;
    my ($pat) = $name =~ /(.[PFC].xls)$/;
    $name =~ s/\.\w+\.xls$//;
#    my ($name,$pat) = $name_tmp =~ /(\S+)(.path.xls)/;
    
    open FILE,$file or die "$!:$file\n";
    my $id;
    <FILE>;
    push @sample_names,$name if(!exists $e{$name});
    $e{$name} = 1;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        if($enrich_type eq "go" || $enrich_type eq 'GO'){
            $id = "$go_pfc{$pat}\t$aa[0] $aa[1]"
        }elsif($enrich_type eq 'kegg' || $enrich_type eq 'KO'){
            $id = "$aa[1]\t$aa[7] $aa[2]";
        }else{
            $id = "aa\t$aa[0] $aa[1]";
        }
        $info{$id}{$name}{num} = $aa[$col{num}];
        $info{$id}{$name}{pv} = $aa[$col{p}];
        $info{$id}{$name}{qv} = $aa[$col{q}];
        $info{$id}{$name}{rf} = $aa[$col{num}] / $aa[$col{all}];
    }
    close FILE;

}

my %p_sig;
my %out_name = qw/go go  GO go  KO kegg  kegg kegg  DO DO  Reactome Reactome/;
die "Unknown [$enrich_type]\n" if(!exists $out_name{$enrich_type});
for my $t(qw/pv qv rf/){
    my $outpfx = "$indir/$t.$out_name{$enrich_type}";
    open OUT,">$outpfx.xls" or die "$!:$outpfx.xls\n";
    print OUT join("\t",$head,@sample_names)."\n";
    for my $id(sort keys %info){
        my @num = map{$info{$id}{$_}{$t}//="NA";$info{$id}{$_}{$t}}@sample_names;
        if($t ne 'rf'){
            my @aa = grep{$_ ne 'NA' && $_ <=0.05} @num;
            if(@aa > 0){
                $p_sig{$id} = 1 if($t eq 'pv');
            }else{
                next;
            }
        }else{
            next if(!exists $p_sig{$id});
        }
        print OUT join("\t",$id,@num)."\n";
    }
    close OUT;
    
    Excu("perl $src/R_plot.pl enrichmentHeatmap -infile $outpfx.xls -type $t -outpfx $outpfx");
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    system($cmd);

    return;
}
