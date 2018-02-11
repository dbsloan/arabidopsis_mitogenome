#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 variant_file fasta_file flank_len list_of_read_files(space delimited) > output\n\n";

my $var_file = shift or die ($usage);
my $fasta_file = shift or die ($usage);
my $flank = shift or die ($usage);

my @read_file_list;
while (my $read_file = shift){
	push (@read_file_list, $read_file);
}

my %fasta = fasta2hash($fasta_file);
my @var_lines = file_to_array($var_file);

my $first_line = shift (@var_lines);

chomp $first_line;
print $first_line;

my @sfl = split (/\t/, $first_line);

my $genome1;
if ($sfl[0] =~ /(.*)\spos$/){
	$genome1 = $1
}else{
	die ("\nERROR: could not extract genome name from $sfl[0]\n\n");
}
 
my $genome2; 
if ($sfl[1] =~ /(.*)\spos$/){
	$genome2 = $1
}else{
	die ("\nERROR: could not extract genome name from $sfl[1]\n\n");
}

exists ($fasta{$genome1}) or die ("\nERROR: could not find $genome1 in $fasta_file\n\n"); 
exists ($fasta{$genome2}) or die ("\nERROR: could not find $genome2 in $fasta_file\n\n"); 

print "\t$genome1 kmer\t$genome2 kmer";

foreach (@read_file_list){
	print "\t$genome1 $_\t$genome2 $_";
}

print "\n";

foreach (@var_lines){
	chomp $_;
	
	print $_;
	
	my @sl = split (/\t/, $_);
	my $start1 = $sl[0] - $flank;
	my $start2 = $sl[1] - $flank;

	my $seq1 = substr ($fasta{$genome1}, $start1 - 1, 2*$flank + 1);	
	my $rc1 = revcom ($seq1);	
	my $seq2 = substr ($fasta{$genome2}, $start2 - 1, 2*$flank + 1);	
	my $rc2 = revcom ($seq2);
	
	print "\t$seq1\t$seq2";
	
	foreach my $reads (@read_file_list){
		my $seq1_count = `grep $seq1 $reads | wc -l`;
		my $rc1_count = `grep $rc1 $reads | wc -l`;
		my $seq2_count = `grep $seq2 $reads | wc -l`;
		my $rc2_count = `grep $rc2 $reads | wc -l`;
		
		my $sum1 = $seq1_count + $rc1_count; 
		my $sum2 = $seq2_count + $rc2_count;
		
		print "\t$sum1\t$sum2"; 
	}	
	
	print "\n";
	
}
