assemblathon2-analysis
======================

This repo contains a motley assortment of unpublished scripts and commands used by Ian Korf, Keith Bradnam, and Joe Fass in the analysis of Assemblathon 2 competition entries (assemblies).
Any published software that was used in the analysis is cited in the Assemblathon 2 paper and available elsewhere (e.g. [CEGMA](http://korflab.ucdavis.edu/Datasets/cegma/)).

1. compass.pl and compass_graph.R - used to align competition assemblies to the fosmid assemblies and produce the length sets described in the manuscript, as well as coverage, validity, multiplicity, and parsimony metrics (also described in the manuscript). These scripts require lastz, samtools, and R (with package ggplot2), and can be blamed o- ..., were written by Joe Fass and Vince Buffalo.
2. assemblathon_stats.pl - used to calculate many of the basic contig- and scaffold-level statistics (requires FAlite.pm)
3. vfr_blast.pl - used to generate tag sequences from validated fosmid regions (VFRs) and BLAST against assemblies. Requires qstaq.pl and FAlite.pm. Uses WU-BLAST or AB-BLAST (not NCBI BLAST)
4 blasting for fosmid read coverage?
5.
 

