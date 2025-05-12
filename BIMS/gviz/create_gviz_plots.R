#load libraryx
library(GenomicAlignments)
library(GenomicRanges)
library(Gviz) 

# minimap no secondary alignments 
# https://github.com/lh3/minimap2/issues/98
# https://github.com/ivanek/Gviz/issues/21

.import.bam.alignments.ignore.secondary <- function(file, selection) {
  indNames <- c(sub("\\.bam$", ".bai", file), paste(file, "bai", sep = "."))
  index <- NULL
  for (i in indNames) {
    if (file.exists(i)) {
      index <- i
      break
    }
  }
  if (is.null(index)) {
    stop(
      "Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
      "library(Rsamtools)\n\tindexBam(\"", file, "\")"
    )
  }
  pairedEnd <- parent.env(environment())[["._isPaired"]]
  if (is.null(pairedEnd)) {
    pairedEnd <- TRUE
  }
  flag <- parent.env(environment())[["._flag"]]
  if (is.null(flag)) {
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
  }
  flag <- bamFlagAND(flag, scanBamFlag(isSupplementaryAlignment = FALSE))
  flag <- bamFlagAND(flag, scanBamFlag(isSecondaryAlignment = FALSE))
  bf <- BamFile(file, index = index, asMates = pairedEnd)
  param <- ScanBamParam(which = selection, what = scanBamWhat(), tag = "MD", flag = flag)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  md <- if (is.null(reads$tag$MD)) rep(as.character(NA), length(reads$pos)) else reads$tag$MD
  if (length(reads$pos)) {
    layed_seq <- sequenceLayer(reads$seq, reads$cigar)
    region <- unlist(bamWhich(param), use.names = FALSE)
    ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
    names(ans) <- seq_along(reads$qname)
  } else {
    ans <- DNAStringSet()
  }
  return(GRanges(
    seqnames = if (is.null(reads$rname)) character() else reads$rname,
    strand = if (is.null(reads$strand)) character() else reads$strand,
    ranges = IRanges(start = reads$pos, width = reads$qwidth),
    id = if (is.null(reads$qname)) character() else reads$qname,
    cigar = if (is.null(reads$cigar)) character() else reads$cigar,
    mapq = if (is.null(reads$mapq)) integer() else reads$mapq,
    flag = if (is.null(reads$flag)) integer() else reads$flag,
    md = md, seq = ans,
    isize = if (is.null(reads$isize)) integer() else reads$isize,
    groupid = if (pairedEnd) if (is.null(reads$groupid)) integer() else reads$groupid else seq_along(reads$pos),
    status = if (pairedEnd) {
      if (is.null(reads$mate_status)) factor(levels = c("mated", "ambiguous", "unmated")) else reads$mate_status
    } else {
      rep(
        factor("unmated", levels = c("mated", "ambiguous", "unmated")),
        length(reads$pos)
      )
    }
  ))
}


##### BIM 10 

fasta_file <- "../ATCC_700802_trycycler_plus_phages_reference.fasta"

bam <- "../BAMS/EF_original_aln_sorted.bam"
bam2 <- "../BAMS/BIM_10_aln_sorted.bam"




options(ucscChromosomeNames=FALSE)
seqset <- readDNAStringSet(fasta_file)

sTrack <- SequenceTrack(seqset,
                        genome = "chromosome0001",
                        chromosome = "chromosome0001")



# https://github.com/ivanek/Gviz/issues/16
alTrack <- AlignmentsTrack(
  bam, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= FALSE,
  name="EF Control", 
  type=c("coverage", "pileup"))


# https://github.com/ivanek/Gviz/issues/16
alTrack2 <- AlignmentsTrack(
  bam2, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= FALSE,
  name="BIM10", 
  type=c("coverage", "pileup"))



# change the gff minorly to get rid of some stuff and change contig_1 to 420_T0
geneTrack <- GeneRegionTrack(
  range = "ATCC_700802.gff3",
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  start = 2085000,
  end = 2095000,
  shape = "arrow", 
  col= "black",
  fill= c("#d0d0d0","#d0d0d0", "#d0d0d0" ,"#0bd3d3", "#d0d0d0"),
  cex = 15.5,         # Increase the size of the feature labels
  #cex.title = 1.5    # Increase the size of the track title
)

geneTrack@range[2270:2275] 

geneTrack@range$transcript[2270:2275] <- c("IS256",
                                           "hydrolase",
                                "Lysozyme",
                                "LicD",
                                "glycosyltransferase",
                                "wcaJ")
library(grDevices)
grDevices::dev.size("in")

pdf("BIM10.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack2),
            from = 2085000, 
            to = 2090000, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()


### BIM 2010


bam3 <- "../BAMS/BIM_2010_aln_sorted.bam"


# https://github.com/ivanek/Gviz/issues/16
alTrack3 <- AlignmentsTrack(
  bam3, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="BIM2010", 
  type=c("coverage", "pileup"))


# change the gff minorly to get rid of some stuff and change contig_1 to 420_T0
geneTrack <- GeneRegionTrack(
  range = "ATCC_700802.gff3",
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  start = 2085000,
  end = 2095000,
  shape = "arrow", 
  col= "black",
  fill= c("#d0d0d0","#d0d0d0", "#d0d0d0" ,"#0bd3d3", "#d0d0d0"),
  cex = 15.5,         # Increase the size of the feature labels
  #cex.title = 1.5    # Increase the size of the track title
)

geneTrack@range[2270:2275] 

geneTrack@range$transcript[2270:2275] <- c("IS256",
                                           "hydrolase",
                                           "Lysozyme",
                                           "LicD",
                                           "glycosyltransferase",
                                           "wcaJ")

pdf("BIM2010.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack,alTrack, alTrack3),
            from = 2085000, 
            to = 2090000, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()




pdf("BIM2010_zoomed.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack3),
            from = 2087773, 
            to = 2087833, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



### BIM 1610


bam4 <- "../BAMS/BIM_1610_aln_sorted.bam"


# https://github.com/ivanek/Gviz/issues/16
alTrack4 <- AlignmentsTrack(
  bam4, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="BIM1610", 
  type=c("coverage", "pileup"))



pdf("BIM1610.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack,alTrack, alTrack4),
            from = 2085000, 
            to = 2090000, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()




pdf("BIM1610_zoomed.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack4),
            from = 2087773, 
            to = 2087833, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



### BIM 16


bam5 <- "../BAMS/BIM_16_aln_sorted.bam"


# https://github.com/ivanek/Gviz/issues/16
alTrack5 <- AlignmentsTrack(
  bam5, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="BIM16", 
  type=c("coverage", "pileup"))



pdf("BIM16.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack,alTrack, alTrack5),
            from = 2085000, 
            to = 2090000, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()




pdf("BIM16_zoomed.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack5),
            from = 2087773, 
            to = 2087833, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()




### BIM 20 resequencing

bam6 <- "../BAMS/BIM_20_200125_aln_sorted.bam"

# https://github.com/ivanek/Gviz/issues/16
alTrack6 <- AlignmentsTrack(
  bam6, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="BIM20", 
  type=c("coverage", "pileup"))



pdf("BIM20_licd_zoomed.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack6),
            from = 2087773, 
            to = 2087833, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()

pdf("BIM20_23S_rrna_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack6),
            from = 251380, 
            to = 251415, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



pdf("BIM20_FabZ_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack6),
            from = 270150, 
            to = 270165, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



pdf("BIM20_wcaG_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack6),
            from = 2073980, 
            to = 2073995, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



### BIM 2016 resequencing

bam7 <- "../BAMS/BIM_2016_aln_sorted.bam"

# https://github.com/ivanek/Gviz/issues/16
alTrack7 <- AlignmentsTrack(
  bam7, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "chromosome0001",
  chromosome = "chromosome0001",
  isPaired = FALSE,
  fill.coverage= "#133e7c",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion = "#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="BIM2016", 
  type=c("coverage", "pileup"))



pdf("BIM2016_licd_zoomed.pdf", width = 10, height = 8)
plotTracks( c(sTrack, geneTrack, alTrack, alTrack7),
            from = 2087773, 
            to = 2087833, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()

pdf("BIM2016_23S_rrna_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack7),
            from = 251380, 
            to = 251415, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



pdf("BIM2016_FabZ_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack7),
            from = 270150, 
            to = 270165, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()



pdf("BIM2016_wcaG_zoomed.pdf", width = 10, height = 8)

plotTracks( c(sTrack, geneTrack, alTrack, alTrack7),
            from = 2073980, 
            to = 2073995, 
            chromosome = "chromosome0001", 
            genome = "chromosome0001", 
            transcriptAnnotation = "transcript",
            background.title="black",
            col.axis="white",
            col.title= "white")

dev.off()




