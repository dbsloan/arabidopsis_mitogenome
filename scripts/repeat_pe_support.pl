#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 repeat_file sam_lines_file flank_dist seq_name\n\n";

my $rep_file = shift or die ($usage);
my $sam_file = shift or die ($usage);
my $flank = shift or die ($usage);
my $seq_name = shift or die ($usage);


my @rep_lines = file_to_array ($rep_file);
my $first_line = shift @rep_lines;

chomp $first_line;
print "$first_line\tRef1\tRef2\tRec1\tRec2\n";

my %reps;
my %reps_PE;

for (my $i = 1; $i <= scalar(@rep_lines); ++$i){
	chomp $rep_lines[$i - 1];
	my @sl = split(/\t/, $rep_lines[$i - 1]);
	for (my $j = 0; $j < scalar(@sl); ++$j){
		$reps{$i}[$j] = $sl[$j];
	}
	$reps_PE{$i}[0]=0;
	$reps_PE{$i}[1]=0;
	$reps_PE{$i}[2]=0;
	$reps_PE{$i}[3]=0;	
}

my $rep_count = scalar (keys %reps);

my @sam_lines = file_to_array ($sam_file);

my $i = 0;
while ($i < scalar (@sam_lines)){
	my @sl1 = split (/\t/, $sam_lines[$i]); 
	my @sl2 = split (/\t/, $sam_lines[$i+1]);
	unless ($sl1[0] eq $sl2[0]){
		++$i;
		next;
	}
	
	unless ($sl1[5] =~ /^\d+M$/ and $sl2[5] =~ /^\d+M$/){
		$i += 2;
		next;
	}
	
	unless ($sl1[2] eq $seq_name and $sl2[2] eq $seq_name){
		$i += 2;
		next;
	}
	
	my $len1 = substr($sl1[5], 0, -1);
	my $len2 = substr($sl2[5], 0, -1);
	
	if ($sl1[1] == 97 || $sl1[1] == 99){ #first read is forward, second is reverse
		
		my $supp_repeats = 0;
		for (my $j = 1; $j <= $rep_count; ++$j){
				
			my $copy1L;
			my $copy1R;
			my $copy2L;
			my $copy2R;

			$sl1[3] >= $reps{$j}[0] - $flank and $sl1[3] < $reps{$j}[0] and $copy1L = 1;
			$sl1[3] >= $reps{$j}[2] - $flank and $sl1[3] < $reps{$j}[2] and $copy2L = 1;
			$sl2[3] + $len2 - 1 > $reps{$j}[1] and $sl2[3] + $len2 - 1 <= $reps{$j}[1] + $flank and $copy1R = 1;
			$sl2[3] + $len2 - 1 > $reps{$j}[3] and $sl2[3] + $len2 - 1 <= $reps{$j}[3] + $flank and $copy2R = 1;			
			
			$copy1L and $copy2L and print STDERR ("WARNING: Both left flanks for repeat $j for read pair $sl1[0]\n") and next;			
			$copy1R and $copy2R and print STDERR ("WARNING: Both right flanks for repeat $j for read pair $sl1[0]\n") and next;
			
			$copy1L and $copy1R and ++$reps_PE{$j}[0] and ++$supp_repeats;
			$copy2L and $copy2R and ++$reps_PE{$j}[1] and ++$supp_repeats;
			
			if ($reps{$j}[4] == 1){
				$copy1L and $copy2R and ++$reps_PE{$j}[2] and ++$supp_repeats;			
				$copy2L and $copy1R and ++$reps_PE{$j}[3] and ++$supp_repeats;
			}
		}
		$supp_repeats > 1 and print STDERR ("WARNING: read pair $sl1[0] counted for more than one repeat pair\n");
			
	}elsif ($sl1[1] == 81 || $sl1[1] == 83){ #first read is reverse, second is forward

		my $supp_repeats = 0;
		for (my $j = 1; $j <= $rep_count; ++$j){


			my $copy1L;
			my $copy1R;
			my $copy2L;
			my $copy2R;
			$sl2[3] >= $reps{$j}[0] - $flank and $sl2[3] < $reps{$j}[0] and $copy1L = 1;
			$sl2[3] >= $reps{$j}[2] - $flank and $sl2[3] < $reps{$j}[2] and $copy2L = 1;
			$sl1[3] + $len2 - 1 > $reps{$j}[1] and $sl1[3] + $len2 - 1 <= $reps{$j}[1] + $flank and $copy1R = 1;
			$sl1[3] + $len2 - 1 > $reps{$j}[3] and $sl1[3] + $len2 - 1 <= $reps{$j}[3] + $flank and $copy2R = 1;			
			
			$copy1L and $copy2L and print STDERR ("WARNING: Both left flanks for repeat $j for read pair $sl1[0]\n") and next;			
			$copy1R and $copy2R and print STDERR ("WARNING: Both right flanks for repeat $j for read pair $sl1[0]\n") and next;
			
			$copy1L and $copy1R and ++$reps_PE{$j}[0] and ++$supp_repeats;
			$copy2L and $copy2R and ++$reps_PE{$j}[1] and ++$supp_repeats;

			if ($reps{$j}[4] == 1){
				$copy1L and $copy2R and ++$reps_PE{$j}[2] and ++$supp_repeats;			
				$copy2L and $copy1R and ++$reps_PE{$j}[3] and ++$supp_repeats;			
			}


		}		
		$supp_repeats > 1 and print STDERR ("WARNING: read pair $sl1[0] counted for more than one repeat pair\n");
	}elsif ($sl1[1] == 65 || $sl1[1] == 67){ #both reads are forward

		my $supp_repeats = 0;
		for (my $j = 1; $j <= $rep_count; ++$j){

			my $copy1L = 0;
			my $copy2L = 0;
		
			$sl1[3] >= $reps{$j}[0] - $flank and $sl1[3] < $reps{$j}[0] and ++$copy1L;
			$sl2[3] >= $reps{$j}[0] - $flank and $sl2[3] < $reps{$j}[0] and ++$copy1L;

			$sl1[3] >= $reps{$j}[2] - $flank and $sl1[3] < $reps{$j}[2] and ++$copy2L;
			$sl2[3] >= $reps{$j}[2] - $flank and $sl2[3] < $reps{$j}[2] and ++$copy2L;
			
			$copy1L > 1 and print STDERR ("WARNING: Both reads map to first left flank for repeat $j for read pair $sl1[0]\n") and next;	 
			$copy2L > 1 and print STDERR ("WARNING: Both reads map to second left flank for repeat $j for read pair $sl1[0]\n") and next;	 

			if ($reps{$j}[4] == -1){	
				$copy1L and $copy2L and ++$reps_PE{$j}[2] and ++$supp_repeats;
			}
		}
		$supp_repeats > 1 and print STDERR ("WARNING: read pair $sl1[0] counted for more than one repeat pair\n");
	}elsif ($sl1[1] == 113 || $sl1[1] == 115){ #both reads are reverse

		my $supp_repeats = 0;
		for (my $j = 1; $j <= $rep_count; ++$j){

			my $copy1R = 0;
			my $copy2R = 0;
		
			$sl1[3] + $len1 - 1 > $reps{$j}[1] and $sl1[3] + $len1 - 1 <= $reps{$j}[1] + $flank and ++$copy1R;
			$sl2[3] + $len2 - 1 > $reps{$j}[1] and $sl2[3] + $len2 - 1 <= $reps{$j}[1] + $flank and ++$copy1R;

			$sl1[3] + $len1 - 1 > $reps{$j}[3] and $sl1[3] + $len1 - 1 <= $reps{$j}[3] + $flank and ++$copy2R;
			$sl2[3] + $len2 - 1 > $reps{$j}[3] and $sl2[3] + $len2 - 1 <= $reps{$j}[3] + $flank and ++$copy2R;

			
			$copy1R > 1 and print STDERR ("WARNING: Both reads map to first right flank for repeat $j for read pair $sl1[0]\n") and next;	 
			$copy2R > 1 and print STDERR ("WARNING: Both reads map to second right flank for repeat $j for read pair $sl1[0]\n") and next;	 

			if ($reps{$j}[4] == -1){	
				$copy1R and $copy2R and ++$reps_PE{$j}[3] and ++$supp_repeats;
			}
		}
		$supp_repeats > 1 and print STDERR ("WARNING: read pair $sl1[0] counted for more than one repeat pair\n");
	}
	$i += 2;
}


for (my $j = 1; $j <= $rep_count; ++$j){
	print "$rep_lines[$j-1]\t$reps_PE{$j}[0]\t$reps_PE{$j}[1]\t$reps_PE{$j}[2]\t$reps_PE{$j}[3]\n";
}
