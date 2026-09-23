## make-vignette-figures.R
##
## Draws the figures that the peaks vignette includes as images instead of
## building them while it knits.
##
## peakUpset.png: the UpSet plot of the consensus peaks. On the R-devel
## builders ComplexHeatmap stops with "node stack overflow" while drawing it
## inside knitr, and it draws fine on the release builders and in a plain R
## session. The vignette still shows the plotPeakUpset() call, the image is
## what that call returns.
##
## The consensus is built with exactly the settings of the vignette, so the
## figure has to be drawn again whenever those settings or the example data
## change.
##
## Run once from the package root, with RegionSetDE installed.
##
## Author: Sebastian Gregoricchio


# ---- Parameters -------------------------------------------------------------
outputDir <- file.path("vignettes", "figures")
figureWidth <- 7
figureHeight <- 4
figureResolution <- 150

dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

library(RegionSetDE)


# ---- Consensus, as in the vignette ------------------------------------------
sampleSheet <- loadExampleData("peakSheet")

blacklist <- loadBlacklist("hg38", seqlevelsStyle = "Ensembl")
greylist <- makeGreylist(sampleSheet, binSize = 10000)
artefactRegions <- c(GenomicRanges::granges(blacklist), GenomicRanges::granges(greylist))

consensus <- loadConsensusPeaks(sampleSheet,
                                groupBy = "condition",
                                excludeRegions = artefactRegions,
                                seqlevelsStyle = "Ensembl")


# ---- UpSet plot ---------------------------------------------------------------
png(filename = file.path(outputDir, "peakUpset.png"),
    width = figureWidth,
    height = figureHeight,
    units = "in",
    res = figureResolution)
ComplexHeatmap::draw(plotPeakUpset(consensus))
dev.off()


# ---- Session ------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path("inst", "scripts", "make-vignette-figures-sessionInfo.txt"))
