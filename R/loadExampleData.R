#' @title loadExampleData
#'
#' @description Loads one of the example objects installed with RegionSetDE. The objects cover the whole path from a collection of genomic regions to a fitted model, and are the ones used throughout the vignette and the examples.
#'
#' @param dataset String indicating which object to load, one among \code{"counts"}, \code{"fit"}, \code{"regions"}, \code{"exclusionRegions"}, \code{"sampleSheet"} and \code{"buildMetadata"}. The string \code{"blacklist"} is accepted as a synonym of \code{"exclusionRegions"}. Default: \code{"counts"}.
#' @param verbose Logical value to indicate whether the loading message must be printed. Default: \code{TRUE}.
#'
#' @return The requested object. \code{"counts"} returns a \code{RegionSetDE.counts} object, unnormalised and unfiltered, with the background bins already stored in its metadata. \code{"fit"} returns a \code{RegionSetDE.fit} object obtained from the same counts after background normalisation and filtering, ready for \code{\link{testRegions}} and \code{\link{testRegionSets}}. \code{"regions"} returns a data.frame with one row per region and its set membership, \code{"exclusionRegions"} a \code{GRanges}, \code{"sampleSheet"} a data.frame describing every library of the source dataset, and \code{"buildMetadata"} a list with the parameters used to generate the others.
#'
#' @details The example data comes from the liver ChIP-seq libraries of the EURATRANS project, distributed by the \code{chromstaRData} package and aligned to rn4. The contrast is H3K4me3 in the spontaneously hypertensive (SHR) rat against the Brown Norway (BN) strain, two biological replicates each, restricted to chromosome 12.
#'
#' The regions are one kilobase windows split into four sets ordered by the amount of H3K4me3 they are expected to carry: promoters overlapping a CpG island, promoters without one, positions inside gene bodies away from any transcription start site, and intergenic positions that serve as the low-signal control. The sets are disjoint, and a window claimed by an earlier set is never reused by a later one.
#'
#' The counts are shipped before normalisation so that the vignette can run \code{\link{normalizeCounts}} and \code{\link{filterRegions}} itself. The fit is shipped after both, since estimating the dispersion and matching the universe is the slow step and every example that needs a model would otherwise repeat it.
#'
#' rn4 has no curated blacklist, so \code{"exclusionRegions"} is assembled from the UCSC assembly gap track and from bins carrying implausible coverage in the input libraries. It is not an ENCODE-grade exclusion list and should not be reused outside these examples.
#'
#' The script that generated every object is installed with the package, at \code{system.file("scripts", "make-data.R", package = "RegionSetDE")}.
#'
#' @examples
#' counts <- loadExampleData("counts")
#' counts
#'
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' topRegions(results, n = 5)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadRegions}}, \code{\link{splitLoadRegions}}, \code{\link{applyBlacklist}}, \code{\link{countReads}}, \code{\link{countBackground}}, \code{\link{fitRegions}}
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

    #-----------------------#
    # Locate the object     #
    #-----------------------#
    filePath <- system.file("extdata", paste0("euratrans_", dataset, ".rds"), package = "RegionSetDE")

    if (filePath == "" | !file.exists(filePath)) {
      stop(paste0("The example object '", dataset, "' is not installed with RegionSetDE, please reinstall the package."), call. = FALSE)
    }

    #-----------------------#
    # Read it back          #
    #-----------------------#
    exampleObject <- readRDS(filePath)

    if (isTRUE(verbose)) {
      message(paste0("Loaded the '", dataset, "' example object (", class(exampleObject)[1], ")."))
    }

    return(exampleObject)
  } # END function
