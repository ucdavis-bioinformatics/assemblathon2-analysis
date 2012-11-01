
library(Rsamtools)
library(GenomicRanges)
library(ggbio)
library(lattice)

# general naming
# normalize around "bird, fish, snake" naming and order
# assemblies named like "fish 12E"
assemblies <- dir(pattern="^[bs].+_.*[EC]$")
bird <- grepl("^b", assemblies)
snake <- grepl("^s", assemblies)
alnFiles <- paste(assemblies, "/compassRun.TEMP.contigsVSref.bam", sep="")
ctgFiles <- paste(assemblies, "/compassRun.OUT.contigLengths", sep="")
unbFiles <- paste(assemblies, "/compassRun.OUT.unbrokenLengths", sep="")
perFiles <- paste(assemblies, "/compassRun.OUT.perfectLengths", sep="")
#assembliesRenamed <- assemblies
#assembliesRenamed <- sub("snake", "boa", assembliesRenamed)
#assembliesRenamed <- sub("bird", "parrot", assembliesRenamed)


## Assemblies --> fosmids (COMPASS)

## single-# stats ...
ss <- read.table("single_stats.txt")
colnames(ss) <- c("assembly", "value", "metric")
ss$value[ss$metric=="Parsimony"] <- ss$value[ss$metric=="Parsimony"]^(-1)  # fix major error in compass.pl, inverting parsimony formula
ss$metric <- ordered(ss$metric, levels=c("Multiplicity","Parsimony","Coverage","Validity"))
#ss$assembly <- sub("snake", "boa", ss$assembly)
#ss$assembly <- sub("bird", "parrot", ss$assembly)
trellis.par.set(canonical.theme(color = FALSE))
ss$organism <- sub("_.*", "", ss$assembly)
ss$assembly <- sub(".*_", "", ss$assembly)
tempss <- subset(ss, assembly != "8C")  # sloppy workaround for ylim issue on bird plot
tempss$assembly <- ordered(tempss$assembly, levels=c("1C","2C","3E","4C","5C","6C","7C","9C","10C","11C","12C","13E","14E","15C"))  # bird assemblies, sorted; removed "8C" for tiny bird submission
png("birdNumbers.png")
barchart(data=tempss, subset = organism=="bird", log10(value) ~ assembly | metric, scales=list(x=list(rot=90), y=list(relation='free')), main="Bird Assemblies", col="black")
dev.off()
rm(tempss)
tempss <- subset(ss, assembly != "12C")  # sloppy workaround for ylim issue on snake plot
tempss$assembly <- ordered(tempss$assembly, levels=c("1C","2C","3C","4C","5C","6C","7C","8C","9C","10C","11C"))  # snake assemblies, sorted; removed "12C" for poorly performing snake submission
png("snakeNumbers.png")
barchart(data=tempss, subset = organism=="snake", log10(value) ~ assembly | metric, scales=list(x=list(rot=90), y=list(relation='free')), main="Snake Assemblies", col="black")
dev.off()
rm(tempss)


## lengths ...
# getting length data from COMPASS files
scafLengths <- lapply(ctgFiles, function(f) lengths <- scan(f))
names(scafLengths) <- assemblies
alignLengths <- lapply(unbFiles, function(f) lengths <- scan(f))
names(alignLengths) <- assemblies
#plengths <- lapply(perFiles, function(f) lengths <- scan(f))
#names(plengths) <- assemblies
scafData <- do.call(rbind, lapply(seq_along(scafLengths), function(i, lnames) { data.frame(length=rev(scafLengths[[i]]), cum.length=rev(cumsum(scafLengths[[i]])), assembly=lnames[i]) }, lnames=names(scafLengths)) )
scafData$lengthOf <- c("scaffolds")
alignData <- do.call(rbind, lapply(seq_along(alignLengths), function(i, lnames) { data.frame(length=rev(alignLengths[[i]]), cum.length=rev(cumsum(alignLengths[[i]])), assembly=lnames[i]) }, lnames=names(alignLengths)) )
alignData$lengthOf <- c("alignments")
plotData <- rbind(scafData, alignData)
plotData$organism <- sub("_.*", "", plotData$assembly)
plotData$assembly <- sub(".*_", "", plotData$assembly)
# plot cumulative c,ulengths for bird and snake, using lattice:
# bird ...
trellis.par.set(canonical.theme(color = TRUE))
birdData <- subset(plotData, organism=="bird" & assembly %in% c("1C","2C","3E","4C","5C","6C","7C","9C","10C","11C","12C","13E","14E","15C"))
birdData$lengthOf <- ordered(birdData$lengthOf, levels=c("scaffolds","alignments"))
birdData$assembly <- ordered(birdData$assembly, levels=c("1C","2C","3E","4C","5C","6C","7C","9C","10C","11C","12C","13E","14E","15C"))  # 8C removed
png("birdCLP.png", width=1360, height=700, units="px", pointsize=12)
trellis.par.set(fontsize=list(text=14,points=12), superpose.line=list(col=rep(c("purple4","mediumpurple1","lightgreen","darkgreen"),times=4,length.out=15), lwd = 2, lty=rep(c(1,2,3,4),each=4,length.out=15)))
xyplot(data=birdData, log10(cum.length) ~ log10(length) | lengthOf, group=assembly, type='a', auto.key=list(space='right', points=FALSE, lines=TRUE), main="Bird Assemblies", scales=list(relation='free'))
dev.off()
# snake ... 
trellis.par.set(canonical.theme(color = TRUE))
snakeData <- subset(plotData, organism=="snake" & assembly %in% c("1C","2C","3C","4C","5C","6C","7C","8C","9C","10C","11C"))
snakeData$lengthOf <- ordered(snakeData$lengthOf, levels=c("scaffolds","alignments"))
snakeData$assembly <- ordered(snakeData$assembly, levels=c("1C","2C","3C","4C","5C","6C","7C","8C","9C","10C","11C"))  # 12C removed
png("snakeCLP.png", width=1360, height=700, units="px", pointsize=12)
trellis.par.set(fontsize=list(text=14,points=12), superpose.line = list(col=rep(c("purple4","mediumpurple1","lightgreen","darkgreen"),times=4,length.out=16), lwd = 2, lty=rep(c(1,2,3,4),each=4,length.out=16)))
xyplot(data=snakeData, log10(cum.length) ~ log10(length) | lengthOf, group=assembly, type='a', auto.key=list(space='right', points=FALSE, lines=TRUE), main="Snake Assemblies", scales=list(relation='free'))
dev.off()




# contigs:
plotData <- do.call(rbind, lapply(seq_along(clengths), function(i, lnames) { data.frame(length=rev(clengths[[i]]), cum.length=rev(cumsum(clengths[[i]])), assembly=lnames[i]) }, lnames=names(clengths)) )
svg("contigs.svg", width=14, height=7)
trellis.par.set(superpose.line = list(col=c("purple2","forestgreen","darkorange2"), lwd = 1.5, lty=c(1,2,3,4)))
xyplot(data=plotData, log10(cum.length) ~ log10(length) | sub("_.*","",assembly), group=assembly, type='a', auto.key=list(space='right', points=FALSE, lines=TRUE))
dev.off
# alignment blocks:
plotData <- do.call(rbind, lapply(seq_along(ulengths), function(i, lnames) { data.frame(length=rev(ulengths[[i]]), cum.length=rev(cumsum(ulengths[[i]])), assembly=lnames[i]) }, lnames=names(ulengths)) )
svg("unbroken.svg", width=14, height=7)
trellis.par.set(superpose.line = list(col=c("purple2","forestgreen","darkorange2"), lwd = 1.5, lty=c(1,2,3,4))); xyplot(data=plotData, log10(cum.length) ~ log10(length) | sub("_.*","",assembly), group=assembly, type='a', auto.key=list(space='right', points=FALSE, lines=TRUE))
dev.off()
# perfect alignment blocks:
#plotData <- do.call(rbind, lapply(seq_along(plengths), function(i, lnames) { data.frame(length=rev(plengths[[i]]), cum.length=rev(cumsum(plengths[[i]])), assembly=lnames[i]) }, lnames=names(plengths)) )
#svg("perfect.svg", width=14, height=7)
#trellis.par.set(superpose.line = list(col=c("purple2","forestgreen","darkorange2"), lwd = 1.5, lty=c(1,2,3,4))); xyplot(data=plotData, log10(cum.length) ~ log10(length) | sub("_.*","",assembly), group=assembly, type='a', auto.key=list(space='right', points=FALSE, lines=TRUE))
#dev.off()




## alignments ...
alignments <- lapply(alnFiles, function(f) readBamGappedAlignments(f)) 
alignmentsGRL <- GRangesList(lapply(alignments, granges))
names(alignmentsGRL) <- assembliesRenamed
autoplot(alignmentsGRL[parrot], xlab=NULL, ylab=NULL,
         main="Parrot assemblies aligned to Parrot fosmid assemblies",
         geom='line', scales='free', indName=NULL, type="none", group.selfish=TRUE)
#autoplot(object, ..., xlab, ylab, main, indName = "grl_name",
#                               geom = NULL, stat = NULL, type = c("none", "sashimi"),
#                               coverage.col = "gray50", coverage.fill = coverage.col,
#                               group.selfish = FALSE, arch.offset = 1.3)



## Competition reads --> fosmids
alnParrot <- readBamGappedAlignments("../reads/Parrot/illumina_uk_qseq/mapToAll/s_6.to.All.bam", use.name=TRUE)
alnParrotGRL <- grglist(alnParrot)  # ?
alnBoa <- readBamGappedAlignments("../reads/Boa/short_inserts/mapToAll/lane1ToAll.bam", use.names=TRUE)
# for quicker testing:
#alnTest <- alnParrot[seqnames(alnParrot) %in% c("NODE_1_length_34498_cov_482.640869", "NODE_7_length_34648_cov_668.417358")]
#alnTestGRL <- grglist(alnTest)
p <- autoplot(alnParrot) + scale_y_log10(); print(p)
## For object BamFile
     ## S4 method for signature 'BamFile'
     autoplot(object, ..., which,
                         xlab, ylab, main, bsgenome, geom = "line", stat = "coverage",
                        method = c("estimate", "raw"), resize.extra = 10, show.coverage = TRUE)



