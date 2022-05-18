# PROGRAM TO TRIM FASTQ FILES RETAIN ONLY EVERY 4TH LINE WITH FIRST n NUCLEOTIDES ONLY

use strict;
use warnings;

my $fil=$ARGV[0]; # Input File
my $kmer=$ARGV[1]; # K-mer Random Sequence

open(f2,$fil);
while (<f2>)
{
	chomp;
	my @ele=split(/\t/,$_);
	my $file_in=$ele[0].".fastq";
	my $file_out=$ele[1].".txt";
	
	print $file_in."\t".$file_out."\n";
	
	open(f3,$file_in);
	open(f3o,'>'.$file_out);
	my $i=0;
	while (<f3>)
	{
		$i=$i+1;
		if (($i%4)==2)
		{
			print f3o substr($_,0,$kmer)."\n";
		}
	}
	close(f3);
	close(f3o);
}
close(f2);




    
