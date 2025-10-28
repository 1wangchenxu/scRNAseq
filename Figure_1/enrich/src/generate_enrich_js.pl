#!/Bio/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Data::Dumper;
# Program:       generate_enrich_js.pl
# Author:        Liu yubin
# Date:          Mon 05 Jul 2021 02:10:32 PM CST
#-----------------------------------------------------------------------------
die "perl $0 <data_dir> <type: ko_data ko_coord go_data> <outdir>\n" unless(@ARGV);
my $in_dir  = shift;
my $type    = shift;    # ko_data ko_coord go_data 
my $out_dir = shift;

#-----------------------------------------------------------------------------
if(!defined $type or $type !~ /^ko_data$|^go_data$|^ko_coord$/) {
    die "[Error]: not right format for tpye: [$type] .\n";
}
if($type eq 'ko_data') {
    ## assume the outdir is ./enrich_src
    my $data_hs = fetch_ko_data($in_dir);
    generate_ko_table_js($data_hs, $out_dir);
}
elsif($type eq 'go_data') {
    ## assume the outdir is ./enrich_src
    my $data_hs = fetch_go_data($in_dir);
    generate_go_table_js($data_hs, $out_dir);
}
elsif($type eq 'ko_coord') {
    ## assume the input the is ./map/coord and the outdir is ./map/js
    my $data_hs = fetch_ko_coord($in_dir);
    generate_coord_js($data_hs, $out_dir);
}

#-----------------------------------------------------------------------------
## generate goTableData.js
sub fetch_go_data {

    my $dir = shift;

    if(!-d $dir) {
        die "[Error]: the go_data input_dir isnot exist. \n";
    }

    my %hash;
    my @files = `find $dir -name *.P.xls -or -name *.F.xls -or -name *.C.xls`;
    for my $file (@files) {
        chomp($file);
        my $dir_name = dirname($file);
        my $bas_name = basename($file);
        $bas_name =~ s/\.xls$//;

        my ($type) = $bas_name =~ /\.(\w+?)$/;
        $hash{$bas_name}{type} = $type;

        open IN, "< $file" or die "$!";
        while (<IN>) {
            chomp;
            my @arr = split /\t/, $_;
            if($. == 1) {
                if($_ !~ /^GO ID/) {
                    die "[Error]: not right format for go.[PFC].xls \n";
                }
                $hash{$bas_name}{colname} = $_;
            }
            else {
                push @{$hash{$bas_name}{data}}, $_;
            }
        }
        close IN;
    }

    return \%hash;
}

sub generate_go_table_js {

    my $hash = shift;
    my $outdir = shift;

    my %namespace = (
        "P" => "Biological Process",
        "F" => "Molecular Function",
        "C" => "Cellular Component",
    );

    my $id_tmp = "";
    for my $id (sort keys %$hash) {
        my $caption1 = "\"$id GO Enrichment ($namespace{$hash->{$id}->{type}})\"";
        my $caption2 = "\"$id GO Enrichment ($namespace{$hash->{$id}->{type}}) Details\"";
        my $picTitle = "\"GO Directed Acycline Graph\"";
        my $picSrc   = "\"./go_img/$id.png\"";

        my @columns = split /\t/, $hash->{$id}->{colname};
        @columns = map { "\"$_\"" } @columns;
        my $columns = "[" . join(",", @columns). "]";

        my @data = ();
        for my $line (@{$hash->{$id}->{data}}) {
            my @tmp = split /\t/, $line;
            @tmp = map { "\"$_\"" } @tmp;
            push @data, "[" . join(",", @tmp) . "]";
        }
        my $data = "[" . join(",", @data) . "]";

        $id_tmp .= <<ID;
    '$id': {
        'caption1': $caption1,
        'caption2': $caption2,
        'picTitle': $picTitle,
        'picSrc': $picSrc,
        'columns': $columns,
        'data': $data
    },
ID
    }
    $id_tmp =~ s/},\n$/}/;

    my $js_tmp = <<JS;
getTableData({
$id_tmp
})
JS

    if(defined $outdir) {
        if(!-d $outdir) {
            die "[Error]: no out_dir exist. \n";
        }

        `mkdir -p $outdir/js`  unless(-d "$outdir/js");
        `mkdir -p $outdir/css` unless(-d "$outdir/css");
        open my $ofh, "> $outdir/js/goTableData.js" or die "$!";
        print $ofh $js_tmp;
        close $ofh;

        if(!-e "$outdir/go.html") {
            `cp -rp $RealBin/enrich_src/go.html $outdir`;
        }
        if(!-e "$outdir/js/tableLoading.min.js") {
            `cp -rp $RealBin/enrich_src/js/tableLoading.min.js $outdir/js`;
        }
        if(!-e "$outdir/css/index.min.css") {
            `cp -rp $RealBin/enrich_src/css/index.min.css $outdir/css`;
        }
    }
    else {
        print "$js_tmp";
    }

    return 0;
}

#-----------------------------------------------------------------------------
## generate koTableData.js
sub fetch_ko_data {

    my $dir = shift;

    if(!-d $dir) {
        die "[Error]: the ko_data input_dir isnot exist. \n";
    }

    my %hash;
    my @files = `find $dir -name "*.path.xls"`;
    for my $file (@files) {
        chomp($file);
        my $d_name = dirname($file);
        my $f_name = (split /[\\\/]/, $file)[-1];
        $f_name =~ s/\.path\.xls$//;

        open IN, "< $file" or die "$!";
        while (<IN>) {
            chomp;
            my @arr = split /\t/, $_;
            if($. == 1) {
                if($_ !~ /^KEGG_A_class/) {
                    die "[Error]: not right format for path.xls \n";
                }
                $hash{$f_name}{colname} = $_;
            }
            else {
                push @{$hash{$f_name}{data}}, $_;
            }
        }
        close IN;
    }

    return \%hash;
}

sub generate_ko_table_js {

    my $hash = shift;
    my $outdir = shift;

    my $id_tmp = "";
    for my $id (sort keys %$hash) {
        my $caption1 = "\"$id Pathway Enrichment\"";
        my $caption2 = "\"$id Pathway Details\"";

        my @columns = split /\t/, $hash->{$id}->{colname};
        @columns = map { "\"$_\"" } @columns;
        my $columns = "[" . join(",", @columns). "]";

        my @data = ();
        for my $line (@{$hash->{$id}->{data}}) {
            my @tmp = split /\t/, $line;
            @tmp = map { "\"$_\"" } @tmp;
            push @data, "[" . join(",", @tmp) . "]";
        }
        my $data = "[" . join(",", @data) . "]";

        $id_tmp .= <<ID;
    '$id': {
        'caption1': $caption1,
        'caption2': $caption2,
        'columns': $columns,
        'data': $data
    },
ID
    }
    $id_tmp =~ s/},\n$/}/;

    my $js_tmp = <<JS;
getTableData({
$id_tmp
})
JS

    if(defined $outdir) {
        if(!-d $outdir) {
            die "[Error]: no out_dir exist. \n";
        }

        `mkdir -p $outdir/js`  unless(-d "$outdir/js");
        `mkdir -p $outdir/css` unless(-d "$outdir/css");
        open my $ofh, "> $outdir/js/koTableData.js" or die "$!";
        print $ofh $js_tmp;
        close $ofh;

        if(!-e "$outdir/ko.html") {
            `cp -rp $RealBin/enrich_src/ko.html $outdir`;
        }
        if(!-e "$outdir/js/tableLoading.min.js") {
            `cp -rp $RealBin/enrich_src/js/tableLoading.min.js $outdir/js`;
        }
        if(!-e "$outdir/js/koDesc.js") {
            `cp -rp $RealBin/enrich_src/js/koDesc.js $outdir/js`;
        }
        if(!-e "$outdir/css/index.min.css") {
            `cp -rp $RealBin/enrich_src/css/index.min.css $outdir/css`;
        }
    }
    else {
        print "$js_tmp";
    }

    return 0;
}

#-----------------------------------------------------------------------------
## generate ko_map coordinates.js
sub fetch_ko_coord {

    my $dir = shift;

    if(!-d $dir) {
        die "[Error]: no coord_dir exist. \n";
    }

    my %coord_hs = ();
#   my @files = glob("$dir/*.coord");
    my @files = `find $dir -name "*.coord"`;
    for my $file (@files) {
        chomp($file);
        open my $fh, "< $file" or die "$!";
        while (<$fh>) {
            chomp;
            my @arr = split /\t/, $_;
            my @k_ids = split /,/, $arr[2];
            my @types = split /,/, $arr[3];
            for my $i (0 .. $#k_ids) {
                $coord_hs{$arr[0]}{$arr[1]}{$k_ids[$i]} = $types[$i];
            }
            $coord_hs{word}{$arr[0]}{$arr[1]}{$arr[2]} = $arr[4];
        }
        close $fh;
    }

    return \%coord_hs;
}

sub generate_coord_js {

    my $hash = shift;
    my $outdir = shift;

    my %color = (
        up     => '#FF0000',
        down   => '#00FF00',
        updown => '#0000FF',
    );

    my $id_tmp = "";
    for my $id (sort keys %$hash) {
        my $path_tmp = "";
        next if($id eq 'word');
        for my $path (sort keys %{$hash->{$id}}) {
            my $kid_tmp = "";
            my @kids = sort keys %{$hash->{$id}->{$path}};

            for my $kid (@kids) {
                die "No [$kid] [$path] [$id]" if(!exists $hash->{$id}->{$path}->{$kid});
                die "No color [$kid] [$path] [$id] [$hash->{$id}->{$path}->{$kid}]" if(!exists $color{$hash->{$id}->{$path}->{$kid}});
                $kid_tmp .= <<KID;
            {
                id: '$kid',
                color: '$color{$hash->{$id}->{$path}->{$kid}}'
            },
KID
            }
            $kid_tmp =~ s/\}\,\n$/\}/;
            $path_tmp .= <<PATH;
        '$path': [
$kid_tmp
        ],
PATH
        }
        $path_tmp =~ s/\]\,\n$/\]/;
        $id_tmp .= <<ID;
    '$id': {
$path_tmp
    },
ID
    }
    $id_tmp =~ s/\}\,\n$/\}/;

    my $word_info = GetWordInfo($hash);
    my $js_tmp = <<JS;
getCoor({
$id_tmp
}$word_info)
JS

    if(defined $outdir) {
        if(!-d $outdir) {
            die "[Error]: no out_dir exist. \n";
        }

        open my $ofh, "> $outdir/coordinates.js" or die "$!";
        print $ofh "$js_tmp";
        close $ofh;

#        if(!-e "$outdir/canvasToImg.min.js") {
#            `cp -rp $RealBin/enrich_src/js/canvasToImg.min.js $outdir`;
#        }
    }
    else {
        print "$js_tmp";
    }

    return 0;
}

sub GetWordInfo{
    my $hash = shift;

    my @glists = keys %{$hash};
   
    my $all_info = "";
    for my $id(@glists){
        next if($id eq 'word');
        my $word_info = "";
        for my $path (sort keys %{$hash->{$id}}) {
            my @kids = sort keys %{$hash->{word}{$id}{$path}};
#            print Dumper($hash);
            my $k_infos = join(",\n",map{"\t\t\t'$_': '$hash->{word}{$id}{$path}{$_}'"}@kids);
#            print "$path\t$k_infos\n";
            $word_info .= "\t\t'$path':{\n$k_infos\n\t\t},\n";
        }
        $word_info =~ s/,\n$/\n/;
        $all_info .= "\t'$id':{\n$word_info\t},\n";
    }

    $all_info =~ s/,\n$/\n/;
    $all_info = ",{\n$all_info}" if($all_info ne '');

    return $all_info;
}

