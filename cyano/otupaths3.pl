# this program reads an distance matrix in phylip format and distributes the sequences into OTUs.
# it uses taxonomic thresholds
# results are given as a taxonomic path
# first argument = distance matrix in phylip format

use strict;

if ($ARGV[0]=~m/\./){die"please remove dots from file name\n";}

# read distance matrix file in phylip format and run mothur
open (MOT, ">$ARGV[0].mothur");
print MOT "cluster(phylip=$ARGV[0], cutoff=0.500, method=furthest, precision=1000)";
close MOT;
system ("mothur", "$ARGV[0].mothur");
#exit;

# extract OTUs deliveries at distinct cutoffs.
open (IN1, "<$ARGV[0]fn.list");my@IN1=<IN1>;
print "\nreading file $ARGV[0]fn.list\n";
my%CUTOFFS=("genus"=>"0.055","family"=>"0.135","order"=>"0.180","class"=>"0.215","phylum"=>"0.250");
my%GENUS;my%FAMILY;my%ORDER;my%CLASS;my%PHYLUM;
while ( my ($key, $value) = each(%CUTOFFS) ) {
	my$a;my$e;my%OTU=();
	for ($a=1;$a<=$#IN1;$a++){
		$e=$a-1;
		chomp $IN1[$a];
		my@COLSA=split(/\t/,$IN1[$a]);
		if ($COLSA[0] == $value){
			print "\n=>file is $COLSA[0], cutoff is $value\n";
			%OTU=&printotulist($ARGV[0],$key,$IN1[$a],$COLSA[0]);
			last;
		}	
		if ($COLSA[0] > $value){
			print "\n=>file is $COLSA[0], cutoff is $value\n";
			print "cutoff = $value not found. It selects the one before because mothur ommits cutoffs providing identical results\n";
			my@COLSE=split(/\t/,$IN1[$e]);			
			chomp $IN1[$e];
			%OTU=&printotulist($ARGV[0],$key,$IN1[$e],$COLSE[0]);
			last;
		}
		if ($a==$#IN1){
			if ($COLSA[0] < $value){
				print "\n=>file is $COLSA[0], cutoff is $value. No distinct OTUs found at cutoff $value = All entries belong to the same $key\n";
				my@TMP=split(/\t/,$IN1[$a]); #separate first two columns and join the rest into a third (and last) column.
				my$tmp1=shift@TMP;
				my$tmp2=shift@TMP;
				my$tmp3=join(",",@TMP);
				$IN1[$a]="$tmp1\t$tmp2\t$tmp3";
				%OTU=&printotulist($ARGV[0],$key,$IN1[$a],$COLSA[0]);
			}
		}	
	}
	if ($key eq "genus"){
		%GENUS=%OTU;
		my$size= keys %GENUS; print "$size entries\n";
	} 
	if ($key eq "family"){
		%FAMILY=%OTU;
		my$size= keys %FAMILY; print "$size entries\n";
	} 
	if ($key eq "order"){
		%ORDER=%OTU;
		my$size= keys %ORDER; print "$size entries\n";
	} 
	if ($key eq "class"){
		%CLASS=%OTU;
		my$size= keys %CLASS; print "$size entries\n";
	} 
	if ($key eq "phylum"){
		%PHYLUM=%OTU;
		my$size= keys %PHYLUM; print "$size entries\n";
	}
}
close (IN1);

# joining OTU deliveries into a single hash
my%FINAL;
foreach my$keygen (keys %GENUS){
	my$valuegen = $GENUS{$keygen};
	foreach my$keyfam (keys %FAMILY){
		my$valuefam = $FAMILY{$keyfam};
		if ($keygen eq $keyfam){
			foreach my$keyord (keys %ORDER){
				my$valueord = $ORDER{$keyord};
				if ($keygen eq $keyord){
					foreach my$keyclas (keys %CLASS){
						my$valueclas = $CLASS{$keyclas};
						if ($keygen eq $keyclas){
							foreach my$keyphy (keys %PHYLUM){
								my$valuephy = $PHYLUM{$keyphy};
								if ($keygen eq $keyphy){
									my$finalchain="gen $valuegen\;fam $valuefam\;ord $valueord\;clas $valueclas\;phyl $valuephy";
									$FINAL{$keygen}=$finalchain;
									last;
								}
							}
						}
					}
				}
			}
		}
	}
}

# print result to file
open (OUT2,">$ARGV[0].otupaths");
foreach my$keyfinal (keys %FINAL){
	my$valuefinal = $FINAL{$keyfinal};
	print OUT2 "$keyfinal\t$valuefinal\n";	
}
close OUT2;
print "\ndone\n";


### sub print OTU lists
sub printotulist {
	my%HASH;
	open (OUT1, ">$_[0].otu.$_[1]");
	print "building otu list with cutoff $_[3]\n";
	my@LIST = split (/\t/,$_[2]);
	my$c=0;
	for (my$b=2;$b<=$#LIST;$b++){ #dismiss first and second columns
		$c++;
		my@TEMP=split(/,/,$LIST[$b]);
		for (my$d=0;$d<=$#TEMP;$d++){
			print OUT1 "$TEMP[$d]\tOTUdata=OTU$c\n";
			$HASH{$TEMP[$d]}=$c;
		}
	}
	close OUT1;
	return(%HASH);
}
