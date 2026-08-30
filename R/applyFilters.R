#' @title .filterRegionSets
#'
#' @description Internal worker shared by \code{applyBlacklist} and \code{applyWhitelist}. Filters each region set against a reference set of regions, either discarding or retaining the regions that overlap it.
#'
#' @param regionSets Named list of \code{GRanges} to filter.
#' @param filterRegions \code{GRanges} used as reference, already reduced.
#' @param keepOverlapping Logical value: \code{TRUE} retains the overlapping regions (whitelist), \code{FALSE} discards them (blacklist).
#' @param overlapType String passed to \code{GenomicRanges::findOverlaps}.
#' @param minOverlapBp Numeric value with the minimum number of overlapping bases.
#' @param minOverlapFraction Numeric value with the minimum fraction of the region that must overlap.
#' @param trimRegions Logical value to indicate whether the regions must be clipped instead of discarded.
#' @param ignoreStrand Logical value to indicate whether the strand must be ignored.
#'
#' @return A named list of filtered \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps setdiff
#' @importFrom IRanges pintersect
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' @importFrom BiocGenerics width
#'
#' @keywords internal
#' @noRd

.filterRegionSets <-
  function(regionSets,
           filterRegions,
           keepOverlapping,
           overlapType,
           minOverlapBp,
           minOverlapFraction,
           trimRegions,
           ignoreStrand) {

    filteredSets <-
      lapply(X = regionSets,
             FUN = function(gr) {
               hits <- GenomicRanges::findOverlaps(query = gr, subject = filterRegions,
                                                   type = overlapType, ignore.strand = ignoreStrand)

               # Covered bases are summed per region, the reference set is reduced so nothing is counted twice
               coveredBp <- rep(0, length(gr))
               if (length(hits) > 0) {
                 bpPerHit <- BiocGenerics::width(IRanges::pintersect(gr[S4Vectors::queryHits(hits)],
                                                                     filterRegions[S4Vectors::subjectHits(hits)],
                                                                     ignore.strand = TRUE))
                 aggregatedBp <- tapply(bpPerHit, S4Vectors::queryHits(hits), sum)
                 coveredBp[as.numeric(names(aggregatedBp))] <- as.numeric(aggregatedBp)
               }

               isOverlapping <- coveredBp > 0 &
                 coveredBp >= minOverlapBp &
                 (coveredBp / BiocGenerics::width(gr)) >= minOverlapFraction

               ### Discard mode, the regions are kept entire
               if (trimRegions == FALSE) {
                 if (keepOverlapping == TRUE) {return(gr[isOverlapping])} else {return(gr[!isOverlapping])}
               }

               ### Trim mode, the regions are clipped at the boundaries of the reference set
               if (keepOverlapping == TRUE) {
                 # One piece per hit, each carrying the metadata of the region it comes from
                 if (length(hits) == 0) {return(gr[0])}
                 trimmed <- IRanges::pintersect(gr[S4Vectors::queryHits(hits)],
                                                filterRegions[S4Vectors::subjectHits(hits)],
                                                ignore.strand = TRUE)
                 S4Vectors::mcols(trimmed) <- S4Vectors::mcols(gr)[S4Vectors::queryHits(hits), , drop = FALSE]
               } else {
                 trimmed <- GenomicRanges::setdiff(gr, filterRegions, ignore.strand = ignoreStrand)
                 if (length(trimmed) == 0) {return(trimmed)}
                 # Each piece falls inside one original region, from which the metadata are recovered
                 parentIdx <- GenomicRanges::findOverlaps(trimmed, gr, select = "first", ignore.strand = ignoreStrand)
                 S4Vectors::mcols(trimmed) <- S4Vectors::mcols(gr)[parentIdx, , drop = FALSE]
               }

               return(trimmed)
             })

    names(filteredSets) <- names(regionSets)
    return(filteredSets)
  } # END function




#' @title .applyRegionFilter
#'
#' @description Internal function handling the input/output classes, the loading of the reference regions and the logging shared by \code{applyBlacklist} and \code{applyWhitelist}.
#'
#' @param regionSet Object to filter.
#' @param filterSet Reference regions, as a path, a \code{GRanges} or a data.frame.
#' @param keepOverlapping Logical value: \code{TRUE} for a whitelist, \code{FALSE} for a blacklist.
#' @param overlapType String passed to \code{GenomicRanges::findOverlaps}.
#' @param minOverlapBp Numeric value with the minimum number of overlapping bases.
#' @param minOverlapFraction Numeric value with the minimum overlapping fraction.
#' @param trimRegions Logical value to indicate whether the regions must be clipped.
#' @param ignoreStrand Logical value to indicate whether the strand must be ignored.
#' @param emptySets String indicating how to handle the emptied sets.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return An object of the same class as \code{regionSet}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges reduce
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom methods is as validObject
#'
#' @keywords internal
#' @noRd

.applyRegionFilter <-
  function(regionSet,
           filterSet,
           keepOverlapping,
           overlapType,
           minOverlapBp,
           minOverlapFraction,
           trimRegions,
           ignoreStrand,
           emptySets,
           verbose) {

    filterLabel <- ifelse(keepOverlapping == TRUE, "whitelist", "blacklist")

    #------------------------#
    # Check of the arguments #
    #------------------------#
    overlapType <- tolower(overlapType[1])
    if (!(overlapType %in% c("any", "start", "end", "within", "equal"))) {
      stop("The 'overlapType' parameter must be one among 'any', 'start', 'end', 'within' or 'equal'.", call. = FALSE)
    }

    emptySets <- tolower(emptySets[1])
    if (!(emptySets %in% c("stop", "remove", "keep"))) {
      stop("The 'emptySets' parameter must be one among 'stop', 'remove' or 'keep'.", call. = FALSE)
    }

    if (minOverlapFraction < 0 | minOverlapFraction > 1) {
      stop("The 'minOverlapFraction' parameter must be a value between 0 and 1.", call. = FALSE)
    }

    #-------------------------#
    # Uniform the input class #
    #-------------------------#
    # The input class is recorded so that the same class can be returned at the end
    inputClass <- class(regionSet)[1]

    if (methods::is(regionSet, "RegionSetDE")) {
      regionSets <- as.list(regionSet@regions)
      targetStyle <- regionSet@seqlevels.style
      targetAssembly <- regionSet@genome.assembly
    } else if (methods::is(regionSet, "RegionSetDE.provenance")) {
      # Removing regions after the counting would leave the library sizes inconsistent with the rows
      stop("The blacklist and the whitelist must be applied before the counting, on the RegionSetDE object returned by 'loadRegions'.", call. = FALSE)
    } else if (methods::is(regionSet, "GRangesList")) {
      regionSets <- as.list(regionSet)
      targetStyle <- NULL
      targetAssembly <- NULL
    } else if (methods::is(regionSet, "GRanges")) {
      regionSets <- list(regions = regionSet)
      targetStyle <- NULL
      targetAssembly <- NULL
    } else if (is.list(regionSet)) {
      regionSets <- regionSet
      targetStyle <- NULL
      targetAssembly <- NULL
    } else {
      stop("The 'regionSet' parameter must be a RegionSetDE object, a GRangesList, a list of GRanges or a GRanges.", call. = FALSE)
    }

    if (is.null(names(regionSets))) {names(regionSets) <- paste0("regions_", seq_along(regionSets))}

    #----------------------------#
    # Load the reference regions #
    #----------------------------#
    # Reusing loadRegions keeps the chromosome style aligned with the sets, otherwise the overlaps silently return zero
    filterRegions <- loadRegions(regions = list(filterSet),
                                 keepMetadata = FALSE,
                                 seqlevelsStyle = ifelse(is.null(targetStyle) | isTRUE(is.na(targetStyle)), "UCSC", targetStyle),
                                 genomeAssembly = targetAssembly,
                                 outputFormat = "list",
                                 verbose = FALSE)[[1]]

    # The reference set is reduced so that the covered bases are not counted twice
    filterRegions <- IRanges::reduce(filterRegions, ignore.strand = TRUE)

    sharedSeqlevels <- vapply(regionSets,
                              function(gr) {any(GenomeInfoDb::seqlevels(gr) %in% GenomeInfoDb::seqlevels(filterRegions))},
                              logical(1))
    if (all(sharedSeqlevels == FALSE)) {
      stop(paste0("The ", filterLabel, " shares no chromosome name with the region sets, please check the chromosome naming styles."), call. = FALSE)
    }

    #--------------------#
    # Filter the regions #
    #--------------------#
    if (verbose == TRUE) {message(paste0("Applying the ", filterLabel, " (", length(filterRegions), " regions):"))}

    nBefore <- vapply(regionSets, length, numeric(1))

    filteredSets <- .filterRegionSets(regionSets = regionSets,
                                      filterRegions = filterRegions,
                                      keepOverlapping = keepOverlapping,
                                      overlapType = overlapType,
                                      minOverlapBp = minOverlapBp,
                                      minOverlapFraction = minOverlapFraction,
                                      trimRegions = trimRegions,
                                      ignoreStrand = ignoreStrand)

    nAfter <- vapply(filteredSets, length, numeric(1))

    if (verbose == TRUE) {
      for (i in seq_along(filteredSets)) {
        message(paste0("  ", names(filteredSets)[i], ": ", nAfter[i], "/", nBefore[i], " regions retained (",
                       round(100 * (nBefore[i] - nAfter[i]) / max(nBefore[i], 1), 1), "% removed)"))
      }
    }

    #-----------------------#
    # Handle the empty sets #
    #-----------------------#
    # A set emptied by the filter breaks any downstream comparison, so it is fatal unless stated otherwise
    if (any(nAfter == 0) & emptySets != "keep") {
      emptyNames <- names(filteredSets)[nAfter == 0]

      if (emptySets == "stop") {
        stop(paste0("The ", filterLabel, " left the following region sets without any region: ",
                    paste(emptyNames, collapse = ", "),
                    ". Set 'emptySets' to 'remove' or 'keep' to proceed."), call. = FALSE)
      } else {
        warning(paste0("Empty region sets removed: ", paste(emptyNames, collapse = ", "), "."), call. = FALSE)
        filteredSets <- filteredSets[nAfter > 0]
      }
    }

    #--------------------#
    # Export the regions #
    #--------------------#
    if (inputClass == "GRanges") {return(filteredSets[[1]])}
    if (inputClass == "list") {return(filteredSets)}
    if (methods::is(regionSet, "GRangesList")) {return(methods::as(filteredSets, "GRangesList"))}

    ### RegionSetDE object, the filter and its counts are stored alongside the regions
    filteringStep <- data.frame(step = filterLabel,
                                region.set = names(regionSets),
                                n.before = nBefore,
                                n.after = nAfter,
                                n.removed = nBefore - nAfter,
                                row.names = NULL,
                                stringsAsFactors = FALSE)

    regionSet@regions <- methods::as(filteredSets, "GRangesList")
    regionSet@filtering.log <- rbind(regionSet@filtering.log, filteringStep)
    regionSet@parameters[[filterLabel]] <- list(overlapType = overlapType,
                                                minOverlapBp = minOverlapBp,
                                                minOverlapFraction = minOverlapFraction,
                                                trimRegions = trimRegions,
                                                ignoreStrand = ignoreStrand)

    if (keepOverlapping == TRUE) {regionSet@whitelist <- filterRegions} else {regionSet@blacklist <- filterRegions}

    methods::validObject(regionSet)
    return(regionSet)
  } # END function




#' @title applyBlacklist
#'
#' @description Removes from every region set the regions overlapping a blacklist, such as the ENCODE blacklisted regions. The blacklist can be provided as a BED-like file, a \code{GRanges} or a data.frame.
#'
#' @param regionSet A \code{RegionSetDE} object, a \code{GRangesList}, a named list of \code{GRanges} or a single \code{GRanges}.
#' @param blacklist String indicating the path to a BED-like file, a \code{GRanges} or a data.frame with the regions to exclude.
#' @param overlapType String indicating the type of overlap required to blacklist a region, one among \code{"any"}, \code{"within"}, \code{"start"}, \code{"end"} or \code{"equal"}. Default: \code{"any"}.
#' @param minOverlapBp Numeric value indicating the minimum number of bases that must overlap the blacklist for a region to be removed. Default: \code{1}.
#' @param minOverlapFraction Numeric value between 0 and 1 indicating the minimum fraction of a region that must overlap the blacklist for it to be removed. Default: \code{0}, any overlap is sufficient.
#' @param trimRegions Logical value to indicate whether the blacklisted portion must be subtracted from the regions instead of removing them entirely. Notice that trimming collapses the regions overlapping each other within the same set. Default: \code{FALSE}.
#' @param ignoreStrand Logical value to indicate whether the strand must be ignored when computing the overlaps. Default: \code{TRUE}.
#' @param emptySets String indicating how to handle the sets left without any region, one among \code{"stop"}, \code{"remove"} or \code{"keep"}. Default: \code{"stop"}.
#' @param verbose Logical value to indicate whether the filtering messages must be printed. Default: \code{TRUE}.
#'
#' @return An object of the same class as the input, with the blacklisted regions removed. For a \code{RegionSetDE} object the blacklist and the filtering counts are stored in the corresponding slots.
#'
#' @examples
#' \dontrun{
#' regions <- applyBlacklist(regions, blacklist = "hg38-blacklist.v2.bed")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{applyWhitelist}}, \code{\link{loadRegions}}
#'
#' @export applyBlacklist

applyBlacklist <-
  function(regionSet,
           blacklist,
           overlapType = "any",
           minOverlapBp = 1,
           minOverlapFraction = 0,
           trimRegions = FALSE,
           ignoreStrand = TRUE,
           emptySets = "stop",
           verbose = TRUE) {

    return(.applyRegionFilter(regionSet = regionSet,
                              filterSet = blacklist,
                              keepOverlapping = FALSE,
                              overlapType = overlapType,
                              minOverlapBp = minOverlapBp,
                              minOverlapFraction = minOverlapFraction,
                              trimRegions = trimRegions,
                              ignoreStrand = ignoreStrand,
                              emptySets = emptySets,
                              verbose = verbose))
  } # END function




#' @title applyWhitelist
#'
#' @description Restricts every region set to the regions overlapping a whitelist, for instance a set of accessible or mappable regions. The whitelist can be provided as a BED-like file, a \code{GRanges} or a data.frame.
#'
#' @param regionSet A \code{RegionSetDE} object, a \code{GRangesList}, a named list of \code{GRanges} or a single \code{GRanges}.
#' @param whitelist String indicating the path to a BED-like file, a \code{GRanges} or a data.frame with the regions to retain.
#' @param overlapType String indicating the type of overlap required to retain a region, one among \code{"any"}, \code{"within"}, \code{"start"}, \code{"end"} or \code{"equal"}. Default: \code{"any"}.
#' @param minOverlapBp Numeric value indicating the minimum number of bases that must overlap the whitelist for a region to be retained. Default: \code{1}.
#' @param minOverlapFraction Numeric value between 0 and 1 indicating the minimum fraction of a region that must overlap the whitelist for it to be retained. Default: \code{0}, any overlap is sufficient.
#' @param trimRegions Logical value to indicate whether the regions must be clipped at the whitelist boundaries instead of being retained entirely. A region spanning two whitelisted blocks is split accordingly. Default: \code{FALSE}.
#' @param ignoreStrand Logical value to indicate whether the strand must be ignored when computing the overlaps. Default: \code{TRUE}.
#' @param emptySets String indicating how to handle the sets left without any region, one among \code{"stop"}, \code{"remove"} or \code{"keep"}. Default: \code{"stop"}.
#' @param verbose Logical value to indicate whether the filtering messages must be printed. Default: \code{TRUE}.
#'
#' @return An object of the same class as the input, restricted to the whitelisted regions. For a \code{RegionSetDE} object the whitelist and the filtering counts are stored in the corresponding slots.
#'
#' @examples
#' \dontrun{
#' regions <- applyWhitelist(regions, whitelist = "hg38-accessible-regions.bed")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{applyBlacklist}}, \code{\link{loadRegions}}
#'
#' @export applyWhitelist

applyWhitelist <-
  function(regionSet,
           whitelist,
           overlapType = "any",
           minOverlapBp = 1,
           minOverlapFraction = 0,
           trimRegions = FALSE,
           ignoreStrand = TRUE,
           emptySets = "stop",
           verbose = TRUE) {

    return(.applyRegionFilter(regionSet = regionSet,
                              filterSet = whitelist,
                              keepOverlapping = TRUE,
                              overlapType = overlapType,
                              minOverlapBp = minOverlapBp,
                              minOverlapFraction = minOverlapFraction,
                              trimRegions = trimRegions,
                              ignoreStrand = ignoreStrand,
                              emptySets = emptySets,
                              verbose = verbose))
  } # END function
