## Session cache for the derived objects. The fit is rebuilt from the packaged
## counts rather than shipped, so it is worth computing only once per session.
.exampleCache <- new.env(parent = emptyenv())


#' @title loadExampleData
#'
#' @description Loads one of the example objects installed with RegionSetDE. The objects cover the whole path from a collection of genomic regions to a fitted model, and are the ones used throughout the vignette and the examples.
#'
#' @param dataset String indicating which object to load, one among \code{"counts"}, \code{"fit"}, \code{"regions"}, \code{"exclusionRegions"}, \code{"sampleSheet"} and \code{"buildMetadata"}. The string \code{"blacklist"} is accepted as a synonym of \code{"exclusionRegions"}. Default: \code{"counts"}.
#' @param verbose Logical value to indicate whether the loading message must be printed. Default: \code{TRUE}.
#'
#' @return The requested object. \code{"counts"} returns a \code{RegionSetDE.counts} object, unnormalised and unfiltered, with the background bins already stored in its metadata. \code{"fit"} returns a \code{RegionSetDE.fit} object ready for \code{\link{testRegions}} and \code{\link{testRegionSets}}. \code{"regions"} returns a data.frame with one row per region and its set membership, \code{"exclusionRegions"} a \code{GRanges}, \code{"sampleSheet"} a data.frame describing every library of the source dataset, and \code{"buildMetadata"} a list with the parameters used to generate the others.
#'
#' @details The example data comes from the liver ChIP-seq libraries of the EURATRANS project, distributed by the \code{chromstaRData} package and aligned to rn4. The contrast is H3K4me3 in the spontaneously hypertensive (SHR) rat against the Brown Norway (BN) strain, two biological replicates each, restricted to chromosome 12.
#'
#' The regions are one kilobase windows split into four sets ordered by the amount of H3K4me3 they are expected to carry: promoters overlapping a CpG island, promoters without one, positions inside gene bodies away from any transcription start site, and intergenic positions that serve as the low-signal control. The sets are disjoint, and a window claimed by an earlier set is never reused by a later one.
#'
#' Only the counts are stored on disk. Asking for \code{"fit"} runs \code{\link{normalizeCounts}}, \code{\link{filterRegions}} and \code{\link{fitRegions}} on them, with the parameters recorded in \code{"buildMetadata"}, and keeps the result for the rest of the session. A fitted model carries the internals of its engine, and those change between releases of \code{edgeR} and \code{limma}, so a serialised fit would break as soon as the engine moved underneath it.
#'
#' rn4 has no curated blacklist, so \code{"exclusionRegions"} is assembled from the UCSC assembly gap track and from bins carrying implausible coverage in the input libraries. It is not an ENCODE-grade exclusion list and should not be reused outside these examples.
#'
#' The script that generated the stored objects is installed with the package, at \code{system.file("scripts", "make-data.R", package = "RegionSetDE")}.
#'
#' @examples
#' counts <- loadExampleData("counts")
#' counts
#'
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' topRegions(results, n = 5, FDR = 1)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadRegions}}, \code{\link{splitLoadRegions}}, \code{\link{applyBlacklist}}, \code{\link{countReads}}, \code{\link{countBackground}}, \code{\link{fitRegions}}
#'
#' @importFrom stats as.formula
#'
#' @export loadExampleData

loadExampleData <-
  function(dataset = "counts",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    dataset <- unname(c("counts" = "counts",
                        "fit" = "fit",
                        "regions" = "regions",
                        "exclusionregions" = "exclusionRegions",
                        "blacklist" = "exclusionRegions",
                        "samplesheet" = "sampleSheet",
                        "buildmetadata" = "buildMetadata")[tolower(dataset[1])])

    if (is.na(dataset)) {
      stop("The 'dataset' parameter must be one among 'counts', 'fit', 'regions', 'exclusionRegions', 'sampleSheet' or 'buildMetadata'.", call. = FALSE)
    }

    #------------------------------------#
    # Objects read straight from the disk #
    #------------------------------------#
    readStoredObject <- function(objectName) {
      filePath <- system.file("extdata", paste0("euratrans_", objectName, ".rds"), package = "RegionSetDE")

      if (filePath == "" | !file.exists(filePath)) {
        stop("The example object '", objectName, "' is not installed with RegionSetDE, please reinstall the package.", call. = FALSE)
      }

      return(readRDS(filePath))
    }

    if (dataset != "fit") {
      exampleObject <- readStoredObject(dataset)

      if (isTRUE(verbose)) {
        message("Loaded the '", dataset, "' example object (", class(exampleObject)[1], ").")
      }

      return(exampleObject)
    }

    #---------------------------------------------------#
    # The fit is derived, not stored: engine internals   #
    # change between releases and a serialised one ages  #
    #---------------------------------------------------#
    if (!is.null(.exampleCache$fit)) {
      if (isTRUE(verbose)) {
        message("Loaded the 'fit' example object from the session cache.")
      }

      return(.exampleCache$fit)
    }

    buildMetadata <- readStoredObject("buildMetadata")
    exampleCounts <- readStoredObject("counts")

    if (isTRUE(verbose)) {
      message("Building the 'fit' example object from the packaged counts.")
    }

    exampleCounts <- normalizeCounts(counts = exampleCounts,
                                     method = buildMetadata$normalizationMethod,
                                     verbose = FALSE)

    exampleCounts <- filterRegions(counts = exampleCounts,
                                   method = buildMetadata$filterMethod,
                                   foldChange = buildMetadata$filterFoldChange,
                                   verbose = FALSE)

    # Built against baseenv() so that the formula does not drag this frame along
    designFormula <- stats::as.formula(buildMetadata$fitDesign, env = baseenv())

    exampleFit <- fitRegions(counts = exampleCounts,
                             design = designFormula,
                             engine = buildMetadata$fitEngine,
                             verbose = FALSE)

    .exampleCache$fit <- exampleFit

    return(exampleFit)
  } # END function
