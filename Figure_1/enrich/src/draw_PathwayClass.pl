#!/usr/bin/perl -w
use strict;

=head1 name

	 draw_PathwayClass.pl

=head1 descripyion
	-f   file of pathway_class
	-t   type of font position : right / left.The default is left.
	-p   prefix of output file

=head1 example
	perl draw_PathwayClass.pl -f xx -t xx -p xx

=cut

use Getopt::Long;
my ($file,$type,$prefix);
my $key_word = "Genes";
GetOptions(
	"f=s"=>\$file,
	"t=s"=>\$type,
	"p=s"=>\$prefix,
    "k=s"=>\$key_word,
);
die ` pod2text $0` unless ($file&&$prefix);
$type = "left" unless ($type);


my $col_code = "#3bb387;#db8253;#F08080;#759fe8;#9F79EE;#A2CD5A";
my $class_all = "Metabolism;Genetic Information Processing;Environmental Information Processing;Cellular Processes;Organismal Systems;Human Diseases";
my @col = split /;/,$col_code;
my @class = split /;/,$class_all;

open IN,"<$file" or die "$! : $file";
my %A;
my %B;
my %C;
my %D;
my %hash;
my %Hash;
my $max = 0;
my $total;
my $mid = <IN>;
my $num =1;
chomp $mid;
my @Mid = split /\t/,$mid;
$Mid[3] =~ /(.*)\((.*)\)/;
$total = $2;
while ( my $line = <IN> ){
	chomp $line;
	my @array = split /\t/,$line;
	$A{$array[0]}++;
	$B{A}{$array[1]} = $array[0];
	$B{gene}{$array[1]} .= "$array[-2];";

#	$B{$array[1]} = $array[0];
#	$C{$array[2]} = $array[3];
#	$hash{$array[0]}{$array[1]}{$array[2]} = $array[3];
#	if ( $array[3] > $max ){
#		$max = $array[3];
#	}
#	if ( $D{$array[1]} ){
#		$D{$array[1]}+=$array[3];
#	}else{
#		$D{$array[1]} = $array[3];
#	}
}
close IN;

for my $i(keys %{$B{gene}}) {
	chop($B{gene}{$i});
	my @aa = split/;/,$B{gene}{$i};
	my %aaa; #= map {$_=> 1} @aa ;
	my $num = 0;
	my $aaaaa = 0;
	for(@aa){
		$aaaaa++;
#		print "$i\t$_\n";
		next if(exists $aaa{$_});
		$num++;
		$aaa{$_} =1;
	}
	$D{$B{A}{$i}}{$i} = $num;
#	print %aaa;
	$max = $max>$D{$B{A}{$i}}{$i}?$max:$D{$B{A}{$i}}{$i};
#	print "$i\t$aaaaa\t$num\n";
}

#for my $ii ( keys %D ){
#	if ( $D{$ii} > $max ){
#		$max = $D{$ii};
#	}
#}


my @ARRAY = keys %{$B{A}} ;
my @array = keys %A;

my $num_2 = @ARRAY; #class_c
my $num_3 = @array; #class_b
my $Height = $num_3*26+$num_2*26+300;
my $height = $Height-170;
$Height = $Height-25;
my $Width;
my $width;
my $scale;
#if ( $max >= 400 ){
#	$scale = int($max/400)+1;
#}else{
#	$scale = int (400/$max)-1;
#	if ( $scale == 0 ){
#		$scale=1;
#	}
#}
my $max_1=(int($max/5)+1)*5;
$scale = 400/$max_1;

open OUT,">$prefix.svg" or die "$!: $prefix.svg";
if ( $type eq "left" ){
	$Width=680+450+20;
	$width=$Width-100;
	print OUT "<?xml version=\"1.0\" encoding=\"utf-8\" ?><!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	print OUT "<svg width=\"$Width\" height=\"$Height\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	print OUT "<line x1=\"580\" y1=\"100\" x2=\"580\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"580\" y1=\"100\" x2=\"$width\" y2=\"100\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"580\" y1=\"$height\" x2=\"$width\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"$width\" y1=\"100\" x2=\"$width\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
}else{
#	$Width = (int($max/50)+9)*10+850;
#	$width = $Width-750;
	$Width=680+450+30+20;
	$width=570;
	print OUT "<?xml version=\"1.0\" encoding=\"utf-8\" ?><!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	print OUT "<svg width=\"$Width\" height=\"$Height\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	print OUT "<line x1=\"100\" y1=\"100\" x2=\"100\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"100\" y1=\"100\" x2=\"$width\" y2=\"100\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"100\" y1=\"$height\" x2=\"$width\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
	print OUT "<line x1=\"$width\" y1=\"100\" x2=\"$width\" y2=\"$height\" style=\"stroke: #000; stroke-width: 1;\" />\n";
}
my $y1 = 90;
for my $k1 ( @class ){
	if ( exists $A{$k1} ){
		my $x1 = $width+30;
		my $x2 = $width+10;
		my $y2;
		$y1+=26;
		$y2=$y1+17;
		my $Col = pop@col;
		unshift @col,$Col;
		if ( $type eq "left" ){
			print OUT "<text x=\"565\" y=\"$y2\" style=\"fill: black; text-anchor: end;\" font-size=\"27\" font-family=\"Arial\"><tspan font-weight=\"bold\">$k1</tspan></text>\n";
		}else{
			print OUT "<text x=\"$x2\" y=\"$y2\" style=\"fill: black; text-anchor: start;\" font-size=\"27\" font-family=\"Arial\"><tspan font-weight=\"bold\">$k1</tspan></text>\n";
		}
		my %Num;
		if ( $type eq "left" ){
#			for my $k2 ( keys %{$hash{$k1}} ){
#				$Num{$k2} = 0;
#				for my $k3 ( keys %{$hash{$k1}{$k2}} ){
#					$Num{$k2} += $hash{$k1}{$k2}{$k3};
#				}
#			}
#			for my $k2 ( sort {$Num{$b} <=> $Num{$a}} keys %Num ){
			for my $k2( sort {$D{$k1}{$b} <=> $D{$k1}{$a}} keys %{$D{$k1}}) {
				$y1+=26;
				$y2=$y1+17;
#				for my $k3 ( keys %{$hash{$k1}{$k2}} ){
#					$y2 = $y1+17;
#				}
				my $x3;
#				if ( $max >= 400 ){
#					$x3 = $Num{$k2}/$scale;
#				}else{
#					$x3 = $Num{$k2}*$scale;
#				}
				$x3=$D{$k1}{$k2}*$scale;
				my $x4 = 590+$x3;
				print OUT "<rect x=\"585\" y=\"$y1\" width=\"$x3\" height=\"20\" style=\"fill:$Col;\" />\n";
				print OUT "<text x=\"560\" y=\"$y2\" style=\"fill: $Col; text-anchor: end;\" font-size=\"25\" font-family=\"Arial\"><tspan font-weight=\"bold\">$k2</tspan></text>\n";
				print OUT "<text x=\"$x4\" y=\"$y2\" style=\"fill: black; text-anchor: start;\" font-size=\"23\" font-family=\"Arial\"><tspan font-weight=\"bold\">$D{$k1}{$k2}</tspan></text>\n";
			}
		}else{
#			for my $k2 ( keys %{$hash{$k1}} ){
			for my $k2( sort {$D{$k1}{$b} <=> $D{$k1}{$a}} keys %{$D{$k1}}) {
				$y1+=26;
				$y2=$y1+17;
				my $Num=$D{$k1}{$k2};
#				for my $k3 ( keys %{$hash{$k1}{$k2}} ){
#					$y2 = $y1+17;
#					$Num = $Num+$hash{$k1}{$k2}{$k3};
#				}
				my $x3;
				my $x4;
#				if ( $max >= 400 ){
#					$x3 = $width-$Num/$scale-5;
#					$x4 = $Num/$scale;
#				}else{
#					$x3 = $width-$Num*$scale-5;
#					$x4 = $Num*$scale;
#				}
				$x3 = $width-$Num*$scale;
				$x4 = $Num*$scale;
				my $x5 = $x3-5;
				print OUT "<rect x=\"$x3\" y=\"$y1\" width=\"$x4\" height=\"20\" style=\"fill:$Col;\" />\n";
				print OUT "<text x=\"$x1\" y=\"$y2\" style=\"fill: $Col; text-anchor: start;\" font-size=\"25\" font-family=\"Arial\"><tspan font-weight=\"bold\">$k2</tspan></text>\n";
				print OUT "<text x=\"$x5\" y=\"$y2\" style=\"fill: black; text-anchor: end;\" font-size=\"23\" font-family=\"Arial\"><tspan font-weight=\"bold\">$Num</tspan></text>\n";
			}
		}
	}
	for my $I ( 0 .. 5 ){
#		my $ll;
#		if ( $max >= 400 ){
#			$ll = 80*$scale;
#			$ll = sprintf "%.0f",$ll;
#		}else{
#			$ll = 80/$scale;
#			$ll = sprintf "%.0f",$ll;
#		}
		if ( $type eq "left" ){
			my $x = 580+$I*80+5;
			my $y = $height+8;
			my $yy = $y+25;
#			my $OUT=$I*$ll;
			my $OUT=$I*$max_1/5;
#			$OUT = sprintf "%.0f",$OUT;
			print OUT "<line x1=\"$x\" y1=\"$height\" x2=\"$x\" y2=\"$y\" style=\"stroke: #000; stroke-width: 2;\" />\n";
			print OUT "<text x=\"$x\" y=\"$yy\" style=\"fill: black; text-anchor: middle;\" font-size=\"17\" font-family=\"Arial\"><tspan font-weight=\"bold\">$OUT</tspan></text>\n";
		}else{
			my $x = $width-$I*80-5;
			my $y = $height+8;
			my $yy = $y+25;
#			my $OUT=$I*$ll;
			my $OUT=$I*$max_1/5;
#			$OUT = sprintf "%.0f",$OUT;
			print OUT "<line x1=\"$x\" y1=\"$height\" x2=\"$x\" y2=\"$y\" style=\"stroke: #000; stroke-width: 2;\" />\n";
			print OUT "<text x=\"$x\" y=\"$yy\" style=\"fill: black; text-anchor: middle;\" font-size=\"20\" font-family=\"Arial\"><tspan font-weight=\"bold\">$OUT</tspan></text>\n";
		}
	}
	my $y3 = $height+75;
	my $x6 = $Width/2;
	if ( $type eq "left" ){
		my $x5 = 290+$width/2;
		print OUT "<text x=\"$x6\" y=\"75\" style=\"fill: black; text-anchor: middle;\" font-size=\"35\" font-family=\"Arial\"><tspan font-weight=\"bold\">KEGG pathway annotation</tspan></text>\n";
		print OUT "<text x=\"$x5\" y=\"$y3\" style=\"fill: black; text-anchor: middle;\" font-size=\"25\" font-family=\"Arial\"><tspan font-weight=\"bold\">Number of $key_word</tspan></text>\n";
	}else{
		my $x5 = 115+$width/2;
		print OUT "<text x=\"$x6\" y=\"75\" style=\"fill: black; text-anchor: middle;\" font-size=\"35\" font-family=\"Arial\"><tspan font-weight=\"bold\">KEGG pathway annotation</tspan></text>\n";
		print OUT "<text x=\"$x5\" y=\"$y3\" style=\"fill: black; text-anchor: end;\" font-size=\"20\" font-family=\"Arial\"><tspan font-weight=\"bold\">Number of $key_word</tspan></text>\n";
	}
}
print OUT "</svg>";
close OUT;

#system "/usr/bin/java -jar /workspace/Web/htdocs/gdcloudtest/Public/src/kogo/bin/batik-1.7/batik-rasterizer.jar -m image/png $prefix\.svg -dpi 300";
#system "/home/miaoxin/Toolkit/svg2xxx/svg2xxx -t png -dpi 300 $prefix\.svg";
system "convert -density 300 $prefix.svg $prefix.png";
