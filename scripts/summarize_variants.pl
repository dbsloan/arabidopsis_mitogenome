#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_fasta > output\n\n";

my $fasta_file = shift or die ($usage);

my %fasta = fasta2hash($fasta_file);


my @keys = sort keys %fasta;
my $seq1pos = 0;
my $seq1name = $keys[0];
my $seq2pos = 0;
my $seq2name = $keys[1];
my $seq1 = $fasta{$seq1name};
my $seq2 = $fasta{$seq2name};
my $aln_len = length ($seq1);
$aln_len == length ($seq2) or die ("\nERROR: sequence lengths do not match\n\n");

my $seq1gap = 0;
my $seq2gap = 0;
my $seq1gaplen = 0;
my $seq2gaplen = 0;
my $seq1gap_partner_start = 0;
my $seq2gap_partner_start = 0;

print "$seq1name pos\t$seq2name pos\tType\tLength\t$seq1name allele\t$seq2name allele\n";

for (my $i = 0; $i < $aln_len; ++$i){
	substr($seq1, $i, 1) eq '-' or ++$seq1pos;
	substr($seq2, $i, 1) eq '-' or ++$seq2pos;
	
	if ($seq1gap){
		substr($seq1, $i, 1) eq '-' and ++$seq1gaplen and next;
		print $seq1pos - 1, "\t$seq1gap_partner_start\tDel\t$seq1gaplen\t\t", uc(substr($seq2, $i - $seq1gaplen, $seq1gaplen)), "\n";
		$seq1gaplen = 0;
		$seq1gap = 0;
		$seq1gap_partner_start = 0
	}

	if ($seq2gap){
		substr($seq2, $i, 1) eq '-' and ++$seq2gaplen and next;
		print "$seq2gap_partner_start\t", $seq2pos - 1, "\tIns\t$seq2gaplen\t", uc(substr($seq1, $i - $seq2gaplen, $seq2gaplen)), "\t\n";
		$seq2gaplen = 0;
		$seq2gap = 0;
		$seq2gap_partner_start = 0;
	}

	substr($seq1, $i, 1) eq substr($seq2, $i, 1) and next;
	
	if (substr($seq1, $i, 1) eq '-'){
		$seq1gap = 1;
		$seq1gaplen = 1;
		$seq1gap_partner_start = $seq2pos;
		next;
	}

	if (substr($seq2, $i, 1) eq '-'){
		$seq2gap = 1;
		$seq2gaplen = 1;
		$seq2gap_partner_start = $seq1pos;
		next;
	}

	print "$seq1pos\t$seq2pos\tSNP\t1\t", uc(substr($seq1, $i, 1)), "\t" ,uc(substr($seq2, $i, 1)), "\n";
	
}

($seq1gap or $seq2gap) and print STDERR "\nWARNING: Open gap at end of alignment not recorded in output\n\n";
