#!/usr/bin/perl -w

# AUTHOR: Joseph Fass, Nikhil Joshi, Michael Lewis
# LAST REVISED: March 2011
#
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2011 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;
use Getopt::Std;

my $usage = "\n\nusage: $0 [options] <reference (fasta)> <contigs (fasta)>\n\n".
            "Aligns assembled contigs to reference using lastz, and generates various measures of ".
            "assembly quality, for comparison purposes.  Measures generated are: NX, longest contig, ".
            "coverage, fidelity, continuity, and parsimony.\n\n".
            "-s #          use samtools binary #\n".
            "-l #          minimum contig size to be considered in measures (default 100 bp)\n".
            "-n #          1/2 (?) the expected assembly size in bp, for NX measure (default 75000000)\n".
            "-c #          consensus quality cutoff for qualifying variant positions (default 20)\n".
            "-x #          maximum length in bp, among reference sequences, contigs, etc.; for graphing (default 50000000)\n".
            "-y #          maximum cumulative length in bp of any sequence set; for graphing (default 200000000)\n\n";
our($opt_b,$opt_s,$opt_l,$opt_n,$opt_c,$opt_x,$opt_y);
getopts('s:l:n:c:x:y:') or die $usage;
if (!defined($opt_s)) {$opt_s = "/share/apps/samtools-0.1.18/samtools"}
if (!defined($opt_l) or !($opt_l =~ m/^[0-9]+$/)) {$opt_l = 100}
if (!defined($opt_n) or !($opt_n =~ m/^[0-9]+$/)) {$opt_n = 75000000}
if (!defined($opt_c) or !($opt_c =~ m/^[0-9]+$/)) {$opt_c = 20}
if (!defined($opt_x) or !($opt_x =~ m/^[0-9]+$/)) {$opt_x = 50000000}
if (!defined($opt_y) or !($opt_y =~ m/^[0-9]+$/)) {$opt_y = 200000000}

# reference and contig filenames
my $ref = shift or die $usage;
my $ctg = shift or die $usage;
# my $newctg = "compass.contigs.lt$opt_l";  # new contig filename, short sequences removed
system "cp $ref compassRun.TEMP.ref";  # only for Assemblathon 2 fosmid evaluation - ref is small!
system "cp $ctg compassRun.TEMP.contigs.original";
print "\n";

## apply contig length cutoff; initial mapping and pileup
# create new fasta of contigs longer than minimum size
my $header;  my %seqHash;  my $line;  my $key;  my $seqLength;
my @lengthList;  my $nonACGTflag=0;
print "<compass> ... processing contigs\n";
open FID, "<compassRun.TEMP.contigs.original";
while ($line = <FID>) {
	if ($line =~ m/^>/) { $header = $line }
	else {
        chomp $line;
        if (!($line =~ m/^[ATCGNatcgn]+$/)) {  # if not nuc's from beg to end
            $nonACGTflag++;
            $line =~ s/[^ATCGNatcgn]/N/g;  # change non-nuc chars to N
        }
        $seqHash{$header} .= $line
    }
}
close FID;
print "NON-NUCLEOTIDE CHARACTERS IN CONTIGS! ... on $nonACGTflag lines\n" if $nonACGTflag;
# print new contig fasta file of contigs longer than cutoff
my $newID = 1;  # having problems with some contig names and lastz ... could happen to ref too!
open FID, ">compassRun.TEMP.contigs";
foreach $key (keys %seqHash) {
	$seqLength = length($seqHash{$key});
	if ( $seqLength >= $opt_l) {
		print FID ">$newID\n".$seqHash{$key}."\n";
		$newID++;
		push @lengthList, $seqLength;
	}
}
close FID;
undef %seqHash;
## Contig calculations: find longest contig, # of contigs, NX, total length
## output sorted contig length list
@lengthList = reverse sort {$a <=> $b} @lengthList;
open FID, ">compassRun.OUT.contigLengths";
my $longestContig = $lengthList[0];
my $contigCount = $#lengthList + 1;
my $contigLength = 0;
my $NX;
# start summing contig lengths to find N(X)
while ($contigLength < $opt_n and $seqLength = shift @lengthList) {
	print FID "$seqLength\n";
	$NX = $seqLength;
	$contigLength += $seqLength;
}
if ($contigLength < $opt_n) { $NX = "N/A" }
# finish summing contig lengths to find total contig length
while ($seqLength = shift @lengthList) {
	print FID "$seqLength\n";
	$contigLength += $seqLength;
}
close FID;
print "<compass> ... done processing contigs! (1/5 steps)\n";
# align new set of contigs against reference, with default parameters
#print "<compass> ... running alignment with lastz\n";
system "/share/apps/lastz/lastz compassRun.TEMP.ref[multiple] compassRun.TEMP.contigs[multiple] --ambiguous=n --ambiguous=iupac --notransition --step=20 --match=1,5 --chain --identity=98 --format=sam 2>> compassRun.TEMP.log > compassRun.TEMP.contigsVSref.sam";
system "$opt_s view -uS compassRun.TEMP.contigsVSref.sam 2>> compassRun.TEMP.log | $opt_s sort - compassRun.TEMP.contigsVSref 2>> compassRun.TEMP.log";
print "<compass> ... done running alignment! (2/5 steps)\n";
# pileup
print "<compass> ... running pileup with 'samtools mpileup'\n";
system "$opt_s faidx compassRun.TEMP.ref";
system "$opt_s mpileup -f compassRun.TEMP.ref compassRun.TEMP.contigsVSref.bam 2> /dev/null > compassRun.TEMP.contigsVSref.fullPileup";
print "<compass> ... done running pileup! (3/5 steps)\n";


## Reference calculations: total reference length, output cumulative length data
print "<compass> ... processing reference\n";
open FID, "<compassRun.TEMP.ref";
my $refLength = 0;  my $refCount = 0;  $nonACGTflag=0;
while ($line = <FID>) {
	if ($line =~ m/^>/) {
		$header = $line;
		$refCount++;
	}
	else {
		chomp $line;
        if (!($line =~ m/^[ATCGNatcgn]+$/)) {  # if not nuc's from beg to end
	    $nonACGTflag++;
            $line =~ s/[^ATCGNatcgn]+//g;  # squeeze out non-nuc's
        }
		$seqHash{$header} .= $line;
		$refLength += length($line);
	}
}
close FID;
print "NON-NUCLEOTIDE CHARACTERS IN REFERENCE! ... on $nonACGTflag lines\n" if $nonACGTflag;
# output sorted reference sequence length list
@lengthList = ();
foreach $key (keys %seqHash) {
	push @lengthList, length($seqHash{$key});
}
@lengthList = reverse sort {$a <=> $b} @lengthList;
open FID, ">compassRun.OUT.refLengths";
while ($seqLength = shift @lengthList) { print FID "$seqLength\n" }
close FID;
print "<compass> ... done with reference! (4/5 steps)\n";


## Calculate continuity, fidelity, parsimony
## and alignment / perfect alignment lengths
## and coverage island lengths
print "<compass> ... processing pileup\n";
open FID, "<compassRun.TEMP.contigsVSref.fullPileup";
my @coverageIsland;  my $islandSize = 0;
my $refName = "N/A";  my $currentCoord = -1;  my $startCoord = -1;
my $coveredPositions = 0;  my @field;  my $pile;  my @nt;  my $nt;
my $i;  my @seq;  my @mappedContig;  my @unbrokenLength;
while ($line = <FID>) {
	@field = split(/\t/, $line);
	if ($field[0] ne $refName or $field[1]-$currentCoord > 1) {  # jumped to new region; end old island and start new one
        push @coverageIsland, ($currentCoord - $startCoord + 1);  # add previous island's length to list (will be 0 on first loop)
        $refName = $field[0];  # set to current ref sequence name
        $currentCoord = $field[1];  # set to new position
        $startCoord = $currentCoord;  # starting new 0-width island
	}
	else {  # current position is immediately after previous position, so increment size of current growing island
        $currentCoord = $field[1];
	}
	$coveredPositions++;  # tally # of covered positions
	# process pileup to separate sequences from each other
	$pile = $field[4];
	$pile =~ s/\^./\^/g;  # remove quality characters after sequence start flag
	$pile =~ s/-(\d+)(??{".{$1}"})//g;  # remove deletion flags and nt's
	$pile =~ s/\+(\d+)(??{".{$1}"})/&/g;  # remove insertion flags and nt's, replace with insertion indicator
	$pile =~ s/\^[\.,]/\(/g;  # convert 2 characters indicating read start to a single '(', if nt matches ref
	$pile =~ s/[\.,]\$/\)/g;  # convert 2 characters indicating read end to a single ')', if nt matches ref
	$pile =~ s/\^[ATCGNatcgn]/\{/g;  # convert 2 read start characters to a single '{', if nt doesn't match ref
	$pile =~ s/[ATCGNatcgn]\$/\}/g;  # convert 2 read end characters to a single '}', if nt doesn't match ref
	@nt = split(//,$pile);
	$i = 0;
	while ($nt = shift @nt) {
		if ($nt eq "(") { $seq[$i] = "." }  # start with generic match
		elsif ($nt eq "{") { $seq[$i] = "X" }  # start with generic mismatch
		elsif ($nt eq ")") {
			$seq[$i] .= ".";
			push @mappedContig, splice(@seq,$i--,1);  # terminate with generic match, remove from list, and print
			push @unbrokenLength, length($mappedContig[$#mappedContig]);
		}
		elsif ($nt eq "}") {
			$seq[$i] .= "X";
			push @mappedContig, splice(@seq,$i--,1);  # terminate with generic mismatch, remove from list, print
			push @unbrokenLength, length($mappedContig[$#mappedContig]);
		}
		else { $seq[$i] .= $nt }
		$i++;
	}
}
close FID;
push @coverageIsland, ($currentCoord - $startCoord + 1);  # add last island length to list
## calculate coverage statistic (name??)
@coverageIsland = reverse sort {$a <=> $b} @coverageIsland;
open FID, ">compassRun.OUT.coverageIslandLengths";
# start summing to find coverage statistic (name??)
my $islandStat;  my $totalIslandLength = 0;
while ($totalIslandLength < $opt_n and $seqLength = shift @coverageIsland) {
    print FID "$seqLength\n";
    $totalIslandLength += $seqLength;
    $islandStat = $seqLength;
}
if ($totalIslandLength < $opt_n) { $islandStat = "N/A" }
# print "$islandStat\n";
# finish summing to find sum of coverage island lengths and finish file of lengths
while ($seqLength = shift @coverageIsland) {
    print FID "$seqLength\n";
    $totalIslandLength += $seqLength;
}
close FID;
my $islandCoverage;
if ($refLength == 0) { $islandCoverage = "N/A" }
else { $islandCoverage = sprintf("%.5e",$totalIslandLength/$refLength) }
##
## calculate coverage fraction
my $coverage;
if ($refLength == 0) { $coverage = "N/A" }
else { $coverage = sprintf("%.5e",$coveredPositions/$refLength) }
##
## write out continuity data
@unbrokenLength = reverse sort {$a <=> $b} @unbrokenLength;
# print join("\n",@unbrokenLength)."\n";
# my $unbrokenCount = $#unbrokenLength;
my $continuity;  my $totalUnbroken = 0;
open FID, ">compassRun.OUT.unbrokenLengths";
# start summing to find continuity
while ($totalUnbroken < $opt_n and $seqLength = shift @unbrokenLength) {
	print FID "$seqLength\n";
	$totalUnbroken += $seqLength;
	$continuity = $seqLength;
}
if ($totalUnbroken < $opt_n) { $continuity = "N/A" }
# print "$continuity\n";
# finish summing to find sum of alignment lengths and finish file of lengths
while ($seqLength = shift @unbrokenLength) {
	print FID "$seqLength\n";
	$totalUnbroken += $seqLength;
}
close FID;
##
## write out fidelity data
my $sequence = join("X",@mappedContig);
$sequence =~ s/&/\.\*/g;  # create splittable insertion character, retaining preceding nt
my @perfectMapper = split(/[A-Za-z\{\}\*]+/,$sequence);
# print join("\n",@perfectMapper)."\n";
my @perfectLength;
for ($i=0; $i<=$#perfectMapper; $i++) {
	push @perfectLength, length($perfectMapper[$i]);
}
@perfectLength = reverse sort {$a <=> $b} @perfectLength;
my $perfectCount = $#perfectLength;
my $fidelity;  my $totalPerfect = 0;
open FID, ">compassRun.OUT.perfectLengths";
while ($totalPerfect < $opt_n and $seqLength = shift @perfectLength) {
	print FID "$seqLength\n";
	$totalPerfect += $seqLength;
	$fidelity = $seqLength;
}
if ($totalPerfect < $opt_n) { $fidelity = "N/A" }
# print "$fidelity\n";
while ($seqLength = shift @perfectLength) {
	print FID "$seqLength\n";
	$totalPerfect += $seqLength;
}
close FID;
print "<compass> ... done with pileup! (5/5 steps)\n";
## calculate validity
my $validity;
if ($totalUnbroken < 1) { $validity = "N/A" }
else { $validity = sprintf("%.5e", $totalUnbroken/$contigLength) }
## calculate multiplicity
my $multiplicity;
if ($coveredPositions < 1) { $multiplicity = "N/A" }
else { $multiplicity = sprintf("%.5e", $totalUnbroken/$coveredPositions) }
## calculate parsimony
my $parsimony;
if ($coveredPositions < 1) { $parsimony = "N/A" }
else { $parsimony = sprintf("%.5e", $contigLength/$coveredPositions) }
print "<compass> ... DONE ... length data in four files named 'compassRun.OUT.~~~'\n\n";


## OUTPUT SUMMARY STATS
# columnar output:
print "By-column output:\n\n";
print "Reference\t$ref\n".
      "Assembly\t$ctg\n".
      "Minimum Contig Length (bp)\t$opt_l\n".
      "Reference Size (bp)\t$refLength\n".
      "# Reference Sequences\t$refCount\n".
      "Assembly Size (bp)\t$contigLength\n".
      "# Assembled Sequences\t$contigCount\n".
      "Longest Contig (bp)\t$longestContig\n".
      "N($opt_n) (bp)\t$NX\n".
      "Coverage (by counting pileup positions)\t$coverage\n".
      "Coverage (by summing \"island\" lengths)\t$islandCoverage\n".
      "Unnamed Coverage Statistic (N($opt_n) of \"island\" lengths)\t$islandStat\n".
      "Continuity($opt_n)\t$continuity\n".
      "Fidelity($opt_n)\t$fidelity\n".
      "Validity\t$validity\n".
      "Multiplicity\t$multiplicity\n".
      "Parsimony\t$parsimony\n\n";
# row output:
print "By-row output:\n\n";
print "Reference\tAssembly\tMinimum Contig Length (bp)\tReference Size (bp)\t".
      "# Reference Sequences\tAssembly Size (bp)\t# Assembled Sequences\t".
      "Longest Contig (bp)\tN($opt_n) (bp)\tCoverage (by counting pileup positions)\t".
      "Coverage (by summing \"island\" lengths)\tUnnamed Coverage Statistic (N($opt_n) of \"island\" lengths)\t".
      "Continuity\tFidelity($opt_n)\tValidity\tMultiplicity\tParsimony\n";
print "$ref\t$ctg\t$opt_l\t$refLength\t$refCount\t$contigLength\t$contigCount\t".
      "$longestContig\t$NX\t$coverage\t$islandCoverage\t$islandStat\t$continuity\t$fidelity\t$validity\t$multiplicity\t$parsimony\n\n";

# R script for graphing:
print "<compass> ... graphing results\n";
# system "Rscript compass_graph.R &> $fnGoo.log";
system "Rscript ./compass_graph.R $opt_x $opt_y &> compassRun.TEMP.log";
print "Done; look for (basename_of_this_directory).graph.tiff, with a cumulative length / density graph from this analysis.\n\n";

print "# IGNORE: sanity check ... pileup positions = $coveredPositions; sum of island sizes = $totalIslandLength\n\n";






