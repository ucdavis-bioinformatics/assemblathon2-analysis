## Compass Graphics
# AUTHOR: Joseph Fass, Vince Buffalo
# LAST REVISED: March 2011
#
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2011 The Regents of University of California, Davis Campus.
# All rights reserved.


require(ggplot2)


wd <- getwd()
args <- commandArgs()
print(args)


# Define files here (should not change)
data.files <- c(reference='compassRun.OUT.refLengths',
				contigs='compassRun.OUT.contigLengths',
                alignments='compassRun.OUT.unbrokenLengths',
                perfectAlignments='compassRun.OUT.perfectLengths')
# data.dir <- 'compassRun.OUT'
data.dir <- getwd();
files <- sapply(data.files, function(file) paste(data.dir, file, sep='/'))
data <- as.list(sapply(files, function(file) read.table(file, header=FALSE)))
# function for calculating cumulative lengths
makeCumLines <- function(x) {
  l <- sort(x, decreasing=TRUE)
  y <- cumsum(l)
  return(data.frame(length=l, cum.length=y))
}
#
all <- lapply(data, makeCumLines)
all <- do.call(rbind, all)
tmp <- strsplit(rownames(all), '.', fixed=TRUE)
sequenceSet <- unlist(lapply(tmp, function(x) x[[1]]))
all <- cbind(sequenceSet, all)
# add length differences to plot horizontal segments
tmp <- split(all$length, list(all$sequenceSet))
tmp <- lapply(tmp, function(x) c(diff(x), 0))
tmp <- unsplit(tmp, all$sequenceSet)
# tmp <- unlist(tmp)
all$length.diff <- tmp

# reorder sequenceSet levels for desired legend order
all$sequenceSet <- factor(all$sequenceSet, levels=unique(all$sequenceSet))

## plotting

# color palette (modified from colorblind-friendly colors here: http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)
CBpalette <- c("#000000","#DD2000","#E69F00","#009EAA")

# # cumulative length
# tiff(file=paste(basename(wd),'.CLP.tiff',sep=''), width=800, height=640, title=basename(wd))
# # create plot
# p <- ggplot(data=all, aes(x=length, y=cum.length, color=sequenceSet))
# # define custom palette
# p <- p + scale_color_manual(name="sequenceSet", values=CBpalette)
# # restrict axes (to enable uniform ranges across multiple graphs)
# p <- p + xlim(1,as.integer(args[6])) + ylim(0,as.integer(args[7]))
# # use log-log scale for x-axis, since many contigs are << genome size
# p <- p + scale_x_log10()
# # p <- p + scale_y_log10()
# # add reference sequence endpoints
# p <- p + geom_point(data=subset(all, sequenceSet=="reference"), aes(x=length, y=cum.length), color="#C0C0C0", size=4)
# p <- p + geom_point(data=subset(all, sequenceSet=="reference"), aes(x=length, y=cum.length), color="#D0D0D0", size=3)
# p <- p + geom_point(data=subset(all, sequenceSet=="reference"), aes(x=length, y=cum.length), color="white", size=2)
# # add dark vertical lines
# p <- p + geom_segment(aes(y=(cum.length-length), yend=cum.length, x=length, xend=length), size=1, alpha=1)
# # add pale connecting lines
# p <- p + geom_segment(aes(x=length+length.diff, xend=length, y=cum.length, yend=cum.length), size=0.5, alpha=0.2)
# # add rugs, reference last so it appears on top of other rugs
# p <- p + geom_rug(data=subset(all, sequenceSet!="reference"), aes(x=length, y=NULL), size=0.1, alpha=0.1)
# p <- p + geom_rug(data=subset(all, sequenceSet=="reference"), aes(x=length, y=NULL), size=0.25, alpha=1)
# # adjust theme coloring and text size
# p <- p + theme_gray(base_size=24)
# # vertical x-axis labels, and add a title
# graphTitle <- basename(wd)
# p <- p + opts(axis.text.x=theme_text(angle=-90, size=14, hjust=0), axis.text.y=theme_text(size=14), title=graphTitle, plot.title=theme_text(size=18))
# print(p)
# dev.off()
# 
# # histogram
# tiff(file=paste(basename(wd),'.histogram.tiff',sep=''), width=800, height=640, title=basename(wd))
# # create plot
# p <- ggplot(data=all, aes(x=length, ..density.., fill=sequenceSet, color=sequenceSet))
# # define custom palette
# p <- p + scale_color_manual(name="sequenceSet", values=CBpalette)
# # define custom palette
# p <- p + scale_fill_manual(name="sequenceSet", values=CBpalette)
# # plot using histogram
# # p <- p + geom_bar(alpha=0.2, binwidth=0.1, position="dodge")
# p <- p + geom_bar(data=subset(all, sequenceSet!="reference"), alpha=0.2, binwidth=0.1, position="dodge")
# p <- p + geom_bar(data=subset(all, sequenceSet=="reference"), alpha=0.2, binwidth=0.1)
# # restrict axes (to enable uniform ranges across multiple graphs)
# p <- p + xlim(1,as.integer(args[6]))
# # use log-log scale for x-axis, since many contigs are << genome size
# p <- p + scale_x_log10()
# # and y-axis, for large counts
# # p <- p + scale_y_log10()
# # polish text
# p <- p + opts(axis.text.x=theme_text(angle=-90, size=14, hjust=0), axis.text.y=theme_text(size=14), title=graphTitle, plot.title=theme_text(size=18))
# # add rugs, reference last so it appears on top of other rugs
# p <- p + geom_rug(data=subset(all, sequenceSet!="reference"), aes(x=length, y=NULL), size=0.1, alpha=0.1)
# p <- p + geom_rug(data=subset(all, sequenceSet=="reference"), aes(x=length, y=NULL), size=0.25, alpha=1)
# print(p)
# dev.off()


## produce combined cumulative length plot *and* length density plot
tiff(file=paste(basename(wd),'.tiff',sep=''), width=800, height=640, title=basename(wd))
# create fake paneled data
all.cum <- cbind(all, panel="Cumulative Length (bp)")
all.hist <- cbind(all, panel="Length Density (per ~ bp)")
allnew <- rbind(all.cum,all.hist)
# create plot
p <- ggplot(data=all, mapping=aes(x=length, y=cum.length, color=sequenceSet))
p <- p + scale_color_manual(name="sequenceSet", values=CBpalette)
p <- p + scale_fill_manual(name="sequenceSet", values=CBpalette)
p <- p + facet_grid(panel ~ ., scale="free")
# p <- p + xlim(1,as.integer(args[6])) + ylim(0,as.integer(args[7]))  # this overrides free scaling
# p <- p + xlim(1,80000000) + ylim(0,225000000)
p <- p + scale_x_log10()
# cumulative length plot panel
p <- p + layer(data=subset(all.cum, sequenceSet=="reference"), geom="point", color="#C0C0C0", size=4)
p <- p + layer(data=subset(all.cum, sequenceSet=="reference"), geom="point", color="#E0E0E0", size=3)
p <- p + layer(data=subset(all.cum, sequenceSet=="reference"), geom="point", color="white", size=2)
p <- p + layer(data=all.cum, geom="segment", mapping=aes(x=length+length.diff, xend=length, y=cum.length, yend=cum.length), size=0.5, alpha=0.2)
p <- p + layer(data=all.cum, geom="segment", mapping=aes(y=(cum.length-length), yend=cum.length, x=length, xend=length), size=1.5, alpha=1)
# histogram panel
# p <- p + layer(data=subset(all.hist, sequenceSet =="reference"), geom="bar", geom_params=list(alpha=0.2, size=0.25), mapping=aes(x=length, ..density.., fill=sequenceSet, color=sequenceSet), stat="bin", stat_params=list(binwidth=0.1))
p <- p + layer(data=subset(all.hist, sequenceSet !="reference"), geom="density", geom_params=list(alpha=0.2, size=0.5), mapping=aes(x=length, ..density.., fill=sequenceSet, color=sequenceSet), stat="bin", stat_params=list(binwidth=0.1))
# p <- p + ylab("")
bs <- 18
p <- p + theme_grey(base_size=bs) + labs(x="Length (bp)", y="")
p <- p + opts(axis.text.x=theme_text(size=0.75*bs), title=basename(wd))
print(p)
dev.off()
