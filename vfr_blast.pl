#!/usr/bin/perl
#
# vfr_blast.pl
#
# A script to match pairs of sequences at increasing distances from a known genome to an assembly
#
# Authors: Ian Korf, Ken Yu, and Keith Bradnam: Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict; use warnings;
use FAlite; use DataBrowser;
use Getopt::Std;
use vars qw($opt_c $opt_d $opt_l $opt_s $opt_o);
getopts('cd:l:s:o');

my $SEED  = 1;
my $DISTANCE = 800;
my $LENGTH = 100;
my $CSV = 0;

die "
usage: $0.pl [options] <assembly blast database> <validated sequence file> <species snake|bird>
options:
  -d <int> distance between reads [$DISTANCE]
  -l <int> read length [$LENGTH]
  -o saves blast output file [default OFF]
  -c write data to CSV file [default OFF]
" unless @ARGV == 3;

my ($BLAST_DB, $INPUT, $SPECIES) = @ARGV;

$DISTANCE = $opt_d if $opt_d;
$LENGTH   = $opt_l if $opt_l;
$SEED     = $opt_s if $opt_s;
$CSV      = 1      if $opt_c;

die "bad seed" unless $SEED == int $SEED and $SEED > 0 and $SEED < 10;
srand($SEED);

# format BLAST databases if not already done
unless (-s "$INPUT.xni") {system("xdformat -n -I $INPUT") == 0 or die}
#unless (-s "$ASSEMBLY.xni")  {system("xdformat -n -I $ASSEMBLY")  == 0 or die}

# find sequence lengths
open(my $fh, "<",  $INPUT) or die "Can't read from $INPUT";
my $fasta = new FAlite($fh);
my $total_length = 0;
my %length;
while (my $entry = $fasta->nextEntry) {
	my ($name) = $entry->def =~ /^>(\S+)/;
	my $len = length($entry->seq);
	$length{$name} = $len;
	$total_length += $len;
}
close($fh);

print STDERR scalar keys %length, " contigs in validated sequences of $total_length bp\n";


# generate $LENGTH bp paired fragment files (if necessary)
my %generated;
my $frags = "${SPECIES}_fragments.L${LENGTH}.D$DISTANCE";

# if fragments file already exists, we don't need to re-generate it, but we do need to count how many 
# fragments there were
my $total_tag_count = 0;
my $total_pair_count = 0;

if (-s $frags){
	$total_tag_count = `grep -c \">\" $frags`;
	chomp($total_tag_count);
} else{

	print STDERR "generating read pairs of length $LENGTH bp, at $DISTANCE bp apart, with seed $SEED\n";
	open(my $out, ">$frags") or die "Can't write to $frags";
	
	my $counter = 0;
	foreach my $name (keys %length) {
	#	print "$name $length{$name}\n";
		my $span = $LENGTH + $DISTANCE + $LENGTH;
		for (my $i = 0; $i < $length{$name} - $span; $i += $span) {
			my $pos1 = $i + 1;
			my $end1 = $pos1 + $LENGTH - 1;
			my $pos2 = $pos1 + $LENGTH + $DISTANCE;
			my $end2 = $pos2 + $LENGTH -1;
			
	#			print "\t$pos1-$end1\t$pos2-$end2\n";
			my ($def1, @seq1) = `xdget -n -a $pos1 -b $end1 $INPUT $name`;
			my ($def2, @seq2) = `xdget -n -a $pos2 -b $end2 $INPUT $name`;
			$def1 =~ s/\s//g;
			chomp @seq1;
			chomp @seq2;
	# 		$generated{$r}++;
	
			$counter++;		
			print $out ">L-$counter $name\n", @seq1, "\n";
			print $out ">R-$counter $name\n", @seq2, "\n";
		}
	}
	close $out;

	$total_tag_count = $counter * 2;
}
$total_pair_count = $total_tag_count / 2 ;

print STDERR "$total_tag_count validated fosmid region tags (VFRTs) generated ($total_pair_count pairs)\n\n";



#############
# Run BLAST
#############

my $minscore = 90; # 95% identity for 100 nt fragment

my $blast_file = "$INPUT.L$LENGTH.D$DISTANCE.blast.out";
my $command = "qstaq.pl -h 0 -s $minscore $BLAST_DB $frags > $blast_file";

unless (-e "$blast_file"){ 
	system("qstaq.pl -h 0 -s $minscore $BLAST_DB $frags > $blast_file") == 0 or die "Can't run $command $!";
} 



##########################
# Process BLAST output
##########################

my %hit;

# track stats at the level of assembly or scaffold
# could probably combine both in one hash
my %assembly_stats;
my %scaffold_stats;
my %location_stats;
my %strand_stats;

my $blast;
open($blast, "<$blast_file") or die "can't open $blast_file";

while (<$blast>) {
	my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, 
		$pct, $ppos, $qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se) = split;

	# want at least 95 nt matching
	next unless $len >= 95;

	my ($scaffold, $assembly) = $sid =~ m/^(\d+)_${SPECIES}_(\d+[CE])$/;
	
	# keep track of how many VFRTs match each assembly
	$assembly_stats{$assembly}{$qid}++;

	# this will be potentially overwritten if a tag matches
	# multiple scaffolds, but we will only be interested in the unique ones
	$scaffold_stats{$assembly}{$qid} = $scaffold;
		
	# track strand of match
	$strand_stats{$assembly}{$qid} = $qf;
	
	# track the start coord of each tag
	$location_stats{$assembly}{$qid} = $ss;
	
	next;
	my ($side, $num) = split("-", $qid);
	push @{$hit{$num}{$side}}, {
		parent => $sid,
		start  => $ss,
		end    => $se,
		strand => $qf,
	}
}

close($blast);

my $out;
my $csv_file = "$SPECIES.csv";
if($CSV){
	open($out, ">", $csv_file) or die "Can't write to $csv_file";
}

if($CSV){

	print $out "Assembly,Number of $total_tag_count VFRTs that match,%Matching VFRTs,";
	print $out "Number of unique VFRTs,%Unique VFRTs,";
	print $out "Number of unique VFRTs pairs matching same scaffold,%Unique pairs matching same scaffold,";
	print $out "Number of unique VFRTs pairs matching same scaffold on same strand,%Unique pairs matching same scaffold on same strand,";
	print $out "Number of unique VFRTs pairs matching same scaffold on same strand at expected distance,%Unique pairs matching same scaffold on same strand at expected distance,";
	
	print $out "\n";
}


foreach my $assembly (sort keys %assembly_stats){

	print "$assembly\n===\n";
	print $out "$assembly," if ($CSV);

	# open a file to record distances between tags
	my $distance_file_name = "${assembly}_${SPECIES}_distances.L${LENGTH}.D$DISTANCE.dat";

	open(my $distance_file, ">", "$distance_file_name") or die "Can't write to $distance_file_name";
		
	my $matched_tags = scalar keys %{$assembly_stats{$assembly}};
	my $percent = sprintf("%.1f", $matched_tags / $total_tag_count * 100);
	print "$matched_tags ($percent%) - VFRTs that match assembly\n";
	print $out "$matched_tags,$percent%," if ($CSV);

	
	# will keep track of maximum number of tags per assembly (duplications/repeats in assembly?)
	my $max_tag_id;
	my $max_tag_count = 0;
	
	# now count how many tags are unique
	my $unique = 0;
	foreach my $tag (sort keys %{$assembly_stats{$assembly}}){
		$unique++ if ($assembly_stats{$assembly}{$tag} == 1);
	}
	$percent = sprintf("%.1f", $unique / $total_tag_count * 100);
	print "$unique ($percent%) - VFRTs that match uniquely to assembly\n";
	print $out "$unique,$percent%," if ($CSV);
	
	
	# now see how many pairs of VFRTs — that have unique matches — match same scaffold
	# and also how many are on same strand
	my $match_same_scaffold = 0;
	my $same_strand_count = 0;
	my $correct_distance = 0;
	my $distance_errors;
	
	for (my $i = 1; $i <= ($total_pair_count); $i++){
		my $ltag = "L-$i";
		my $rtag = "R-$i";
		
		# first check that both tags had a match
		next unless (exists $scaffold_stats{$assembly}{$ltag});	
		next unless (exists $scaffold_stats{$assembly}{$rtag});	
	
		# keep track of maximum number of hits per tag
		my $ltag_count = $assembly_stats{$assembly}{$ltag};
		my $rtag_count = $assembly_stats{$assembly}{$rtag};

		if ($ltag_count > $max_tag_count){
			$max_tag_count = $ltag_count;
			$max_tag_id = $ltag;
		}
		if ($rtag_count > $max_tag_count){
			$max_tag_count = $rtag_count;
			$max_tag_id = $rtag;
		}

		# only continue now if tag counts are both unique
		next unless ($ltag_count == 1 && $rtag_count == 1);
		
		# do they both match the same scaffold		
		if ($scaffold_stats{$assembly}{$ltag} eq $scaffold_stats{$assembly}{$rtag}){
			$match_same_scaffold++;
		}

		# only continue now if tags both match same scaffold
		next unless ($scaffold_stats{$assembly}{$ltag} eq $scaffold_stats{$assembly}{$rtag});
		
	
		# check for same orientation
		if ($strand_stats{$assembly}{$ltag} eq $strand_stats{$assembly}{$rtag}){
			$same_strand_count++;
		}

		# only continue now if tags are on same strand
		next unless ($strand_stats{$assembly}{$ltag} eq $strand_stats{$assembly}{$rtag});

		# check for distance
		my ($lstart, $rstart) = ($location_stats{$assembly}{$ltag}, $location_stats{$assembly}{$rtag});
		my $distance = $rstart - $lstart;
		$distance = $lstart - $rstart if ($lstart > $rstart);

		print $distance_file "$distance\n";

#		my $tolerance_low  = ($DISTANCE + $LENGTH) * 0.95;
#		my $tolerance_high = ($DISTANCE + $LENGTH) * 1.05;
#		# The distance between the pair of fragments in the assembly must be 95-105% of the distance between fragments in the known genome
#		if ($distance >= $tolerance_low && $distance <= $tolerance_high){
#			$correct_distance++;
#		} else{
#			$distance_errors .= "\tINCORRECT DISTANCE: $ltag & $rtag: $distance bp between tags\n";
#		}

		if ($distance == $DISTANCE + $LENGTH){
			$correct_distance++;
		} else{
	 		$distance_errors .= "\tINCORRECT DISTANCE: $ltag & $rtag: $distance bp between tags\n";
		}


	} 
	$percent = sprintf("%.1f", ($match_same_scaffold / $total_pair_count) * 100);
	print "$match_same_scaffold ($percent%) - unique pairs that match same scaffold\n";
	print $out "$match_same_scaffold,$percent%," if ($CSV);

	$percent = sprintf("%.1f", ($same_strand_count / $total_pair_count) * 100);
	print "$same_strand_count ($percent%) - unique pairs, match same scaffold, in same orientation\n";
	print $out "$same_strand_count,$percent%," if ($CSV);

	$percent = sprintf("%.1f", ($correct_distance / $total_pair_count) * 100);
	my $expected_distance = $DISTANCE + $LENGTH;
	print "$correct_distance ($percent%) - unique pairs, match same scaffold, in same orientation, at expected distance (~$expected_distance bp)\n";
	print $out "$correct_distance,$percent%," if ($CSV);

	print "$distance_errors" if ($distance_errors);

	print "Most hits to any one tag: $max_tag_id ($max_tag_count hits)\n";
	print "\n";	
	print $out "\n" if ($CSV);
	close($distance_file);
}

close($out) if ($CSV);
exit;

__END__


my $tolerance_low  = $r * 0.95;
my $tolerance_high = $r * 1.05;
my $count = 0;

OUTER: foreach my $num (keys %hit) {
	my $left  = $hit{$num}{L};
	my $right = $hit{$num}{R};
	next unless defined $left and defined $right;
	foreach my $hsp1 (@$left) {
		foreach my $hsp2 (@$right) {
			next if $hsp1->{parent} ne $hsp2->{parent}; # both fragments must match same contig/scaffold
			next if $hsp1->{strand} ne $hsp2->{strand}; # both fragments must be in same orientation

			# calculate distance between fragments, but also need to check if we have a reverse strand match
			my $distance = abs($hsp1->{end}+1 - $hsp2->{start});
			$distance    = abs($hsp2->{end}+1 - $hsp1->{start}) if ($hsp1->{start} > $hsp2->{start});
			
			# The distance between the pair of fragments in the assembly must be 95-105% of the distance between fragments in the known genome
			if ($distance >= $tolerance_low && $distance <= $tolerance_high){
				$count++;
				next OUTER;			
			}
		}
	}
}
printf "%d\t%.4f\n", $r, $count / $generated{$r};

