#' Load the packaged NF-YA example objects
#'
#' Loads one of the precomputed objects shipped in \code{inst/extdata}. These
#' were built from the NF-YA ChIP-seq libraries distributed by
#' \pkg{chipseqDBData}, restricted to chr17 and chr19 of mm10. See
#' \code{inst/scripts/make-data.R} for how each object was generated.
#'
#' @param dataset Character, the object to load. One of \code{"readCounts"},
#'   \code{"backgroundCounts"}, \code{"regions"}, \code{"regionSets"},
#'   \code{"blacklist"}, \code{"sampleSheet"} or \code{"buildMetadata"}.
#'
#' @return The requested object.
#'
#' @examples
#' regionSets <- loadExampleData("regionSets")
#' head(regionSets)
#'
#' @author Sebastian Gregoricchio
#'
#' @export
loadExampleData <- function(dataset = c("readCounts", "backgroundCounts",
                                        "regions", "regionSets", "blacklist",
                                        "sampleSheet", "buildMetadata")) {
  # Resolve the requested object to a file shipped with the package.
  dataset <- match.arg(dataset)
  filePath <- system.file("extdata", paste0("nfya_", dataset, ".rds"),
                          package = "RegionSetDE")

  if (filePath == "") {
    stop("Could not find the '", dataset, "' example object. Reinstall RegionSetDE.")
  }

  readRDS(filePath)
}
