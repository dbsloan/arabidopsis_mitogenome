# arabidopsis_mitogenome


## Overview:

The Perl scripts included in this repository were used to validate a new assembly of the Arabidopsis thaliana Col-0 mitochondrial genome. A related manuscript is currently [available on bioRxiv](https://www.biorxiv.org/content/early/2018/01/18/249086) and has been submitted for formal publication. The scripts require the [sloan.pm Perl module](https://github.com/dbsloan/perl_modules). None of the scripts require BioPerl, but some of the (unused) functions in the sloan.pm module do. If you don't have BioPerl installed, simply deleted those functions from the module before attempting to run. Example inputs for some of the scripts are provided in the input_files subdirectory. Scripts have been run on Mac OSX 10.13 and/or Linux CentOS 6 operating systems.

## summarize_variants.pl

This script takes an aligned fasta file containing two sequences as input and reports a tabular summary of sequence differences between the two (tab-delimited text).

Example usage:

perl  summarize_variants.pl  postalign.fas  >  output.txt


## kmer_grep.pl

This scripts takes the variant file output from summarize_variants.pl and the corresponding unaligned fasta file (not the aligned file that was used with summarize_variants.pl) as inputs. The user most also specify at least one fastq read file and a flanking sequence length. The script will extract k-mers surrounding each variant in the variant file and perform a read count for each k-mer pair on each specified fastq file. The k-mer length will be twice the specified flanking sequence length plus one (e.g., entering a flanking sequence length of 12 will give a k-mer length of 25). This is an inefficient script that relies on the unix grep command. It should only be used on small datasets.


Example usage:

perl  kmer_grep.pl  variant_table.txt  prealign.fas  12  fastq1  fastq2  fastq3...  >  output.txt


## repeat_pe_support.pl

This script takes a file summarizing repeat pairs in a genome (see input_files/At_repeats.txt for example format) and a SAM file of mapped reads to that genome as inputs. The user must also specify the amount of flanking sequence to search for adjacent mapping next to each repeat and the name of the genome (same name used in reference mapping). It returns counts of paired-end reads that support (i.e., span in a consistent orientation) each repeat and counts of paired-end reads that support the expected pair of repeats that would result from recombination between the two sequences,

Example usage:

perl  repeat_pe_support.pl  At_repeats.txt  input_sam_file  500  genome_seq_name  >  output.txt

Note that this script has only been used on reduced SAM files that have been pared down with grep by only taking lines that contain the name of target genome of interest. Performance on larger input SAM files is unknown.
