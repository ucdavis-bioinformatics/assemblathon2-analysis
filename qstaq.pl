#!/usr/bin/perl
#
# qstaq.pl
#
# A script to multiplex BLAST searches
#
# Author: Ian Korf, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.
use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_d $opt_g $opt_h $opt_s $opt_w $opt_x $opt_m $opt_n);
getopts('dgh:s:w:x:m:n:');

my $MINSCORE  = 20;
my $MAXCONTIG = 100000;
my $WORDSIZE  = 13;
my $HSPMAX    = 0;
my $MISMATCH  = -1;
my $MATCH     = 1;

die "
QSTAQ - Multiplex BLASTN (fixed read length)
usage: qstaq.pl [options] <blastn db> <fasta file>
options:
  -d         dust sequence       [default off]
  -g         allow gaps (-3)     [default no gaps]
  -h  <int>  hsp max             [$HSPMAX]
  -s  <int>  min score           [default $MINSCORE]
  -m  <int>  match score         [default $MATCH]
  -n  <int>  mismatch score      [default $MISMATCH]
  -w  <int>  wordsize            [default $WORDSIZE]
  -x  <int>  multiplex length    [default $MAXCONTIG]
" unless @ARGV == 2;

my ($DB, $QUERY) = @ARGV;
my $DUST   = $opt_d ? 'filter=dust' : '';
my $GAPS   = $opt_g ? 'Q=3 R=3' : 'nogap';
$MINSCORE  = $opt_s if $opt_s;
$WORDSIZE  = $opt_w if $opt_w;
$MAXCONTIG = $opt_x if $opt_x;
$HSPMAX    = $opt_h if $opt_h;
$MATCH     = $opt_m if $opt_m;
$MISMATCH  = $opt_n if $opt_n;

###############################################################################
# Part 1: Stack the query sequences
###############################################################################
my $tmp_dna = "tmp.qstaq.$$.dna";
open(OUT, ">$tmp_dna") or die;
open(IN, $QUERY) or die;
my $contig_length = 0;
my $def_count     = 0;
my $read_count    = 0;
my $read_length;
my @id;

print OUT ">sequence-$def_count\n";
while (<IN>) {
	if (/^>(\S+)/) {
		my $id = $1;
		my $seq = <IN>;
		if ($contig_length > $MAXCONTIG) {
			$def_count++;
			print OUT ">sequence-$def_count\n";
			$contig_length = 0;
			$read_count = 0;
		}
		$read_length = length($seq) if not defined $read_length;
		die "reads must all be equal lengths\n" if length($seq) != $read_length;
		print OUT $seq, "-\n";
		push @{$id[$def_count]}, $id;
		$contig_length += $read_length;
		$read_count++;
		
	} else {die "unexpected file format"}
}

close IN;
close OUT;

###############################################################################
# Part 2: BLASTN search
###############################################################################
my $p1 = "M=$MATCH N=$MISMATCH kap mformat=3 hspmax=0 B=2147483647 V=0";
my $p2 = "$GAPS W=$WORDSIZE S=$MINSCORE";
$p2 .= " hspmax=$HSPMAX -warnings" if $HSPMAX;
my $tmp_rep = "tmp.qstaq.$$.blast";
system("blastn $DB $tmp_dna $p1 $p2 > $tmp_rep") == 0 or die;

###############################################################################
# Part 3: Convert output
###############################################################################
open(IN, $tmp_rep) or die;
while (<IN>) {
	next unless /^\w/;
	my @f = split;
	my ($contig_num) = $f[0] =~ /^sequence-(\d+)/;
	my $seq_num = int($f[17] / $read_length);
	my $offset = $seq_num * $read_length;
	$f[17] -= $offset;
	$f[18] -= $offset;
	$f[0] = $id[$contig_num][$seq_num];
	print join("\t", @f), "\n";
}
close IN;

###############################################################################
# Part 4: Clean up
###############################################################################

END {
	unlink $tmp_dna if defined $tmp_dna and -e $tmp_dna;
	unlink $tmp_rep if defined $tmp_rep and -e $tmp_rep;
}

