# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

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




#' @title .poolStoredList
#'
#' @description Puts the list a filter has just applied together with the one the object already stores, so that the slot describes every call rather than the last one only. Blacklists add up, while whitelists restrict each other.
#'
#' @param storedList \code{GRanges} already stored in the object, or \code{NULL}.
#' @param newList \code{GRanges} with the list just applied, already reduced.
#' @param sharedOnly Logical value: \code{FALSE} returns the union of the two lists, as for a blacklist, \code{TRUE} the stretches they have in common, as for a whitelist. Default: \code{FALSE}.
#'
#' @return A \code{GRanges} without overlaps, \code{newList} itself when nothing was stored.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges granges intersect
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges ranges reduce
#'
#' @keywords internal
#' @noRd

.poolStoredList <-
  function(storedList,
           newList,
           sharedOnly = FALSE) {

    if (is.null(storedList)) {
      return(newList)
    }

    combineLists <- function(firstList, secondList) {
      if (isTRUE(sharedOnly)) {
        return(GenomicRanges::intersect(firstList, secondList, ignore.strand = TRUE))
      }
      return(IRanges::reduce(c(firstList, secondList), ignore.strand = TRUE))
    }

    plainRanges <- function(gr) {
      return(GenomicRanges::GRanges(seqnames = as.character(GenomeInfoDb::seqnames(gr)), ranges = IRanges::ranges(gr)))
    }

    # The sequence information is kept as long as the two lists agree on it. Lists loaded at different times
    # may differ in the chromosome lengths or in the assembly tag, and are then combined on the names alone
    pooledList <- suppressWarnings(tryCatch(expr = combineLists(GenomicRanges::granges(storedList), GenomicRanges::granges(newList)),
                                            error = function(e) {return(combineLists(plainRanges(storedList), plainRanges(newList)))}))

    return(pooledList)
  } # END function




#' @title .applyRegionFilter
#'
#' @description Internal function handling the input/output classes, the loading of the reference regions and the logging shared by \code{applyBlacklist}, \code{applyGreylist} and \code{applyWhitelist}.
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
#' @param filterLabel String with the name of the filter in the messages and in the \code{filtering.log}, one among \code{"blacklist"}, \code{"greylist"} and \code{"whitelist"}. Default: \code{NULL}, \code{"whitelist"} or \code{"blacklist"} following \code{keepOverlapping}.
#'
#' @return An object of the same class as \code{regionSet}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges reduce
#' @importFrom GenomeInfoDb seqlevels genome
#' @importFrom BiocGenerics width
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
           verbose,
           filterLabel = NULL) {

    if (is.null(filterLabel)) {filterLabel <- ifelse(keepOverlapping == TRUE, "whitelist", "blacklist")}

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
      stop("The ", filterLabel, " must be applied before the counting, on the RegionSetDE object returned by 'loadRegions'.", call. = FALSE)
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
    # A list built for another assembly matches chromosome names and nothing else, which is silent and wrong
    filterAssembly <- if (methods::is(filterSet, "GRanges")) {unique(as.character(GenomeInfoDb::genome(filterSet)))} else {character(0)}
    filterAssembly <- filterAssembly[!is.na(filterAssembly) & filterAssembly != ""]
    regionAssembly <- if (is.null(targetAssembly)) {character(0)} else {targetAssembly[!is.na(targetAssembly) & targetAssembly != ""]}

    # GRCh38 and hg38 name the same assembly, the aliases are resolved before the comparison
    if (length(filterAssembly) == 1 & length(regionAssembly) == 1 && .resolveGenomeName(filterAssembly) != .resolveGenomeName(regionAssembly)) {
      stop("The ", filterLabel, " was built for ", filterAssembly, " and the regions for ", regionAssembly,
           ". Overlapping them would match the chromosome names and nothing else. Clear the assembly with genome() on one of the two to force it.", call. = FALSE)
    }

    # Reusing loadRegions keeps the chromosome style aligned with the sets, otherwise the overlaps silently return zero
    filterRegions <- loadRegions(regions = list(filterSet),
                                 keepMetadata = FALSE,
                                 seqlevelsStyle = if (is.null(targetStyle) || is.na(targetStyle)) {NULL} else {targetStyle},
                                 genomeAssembly = targetAssembly,
                                 outputFormat = "list",
                                 verbose = FALSE)[[1]]

    # Sets loaded without a style, or given as plain ranges, declare none: the list then follows the names they
    # actually carry, chromosome by chromosome. A list that cannot be reconciled is left to the check below
    regionSeqlevels <- unique(unlist(lapply(regionSets, GenomeInfoDb::seqlevels), use.names = FALSE))
    filterRegions <- tryCatch(expr = .matchSeqlevels(x = filterRegions, targetSeqlevels = regionSeqlevels, verbose = FALSE),
                              error = function(e) {return(filterRegions)})

    # The reference set is reduced so that the covered bases are not counted twice
    filterRegions <- IRanges::reduce(filterRegions, ignore.strand = TRUE)

    sharedSeqlevels <- vapply(regionSets,
                              function(gr) {any(GenomeInfoDb::seqlevels(gr) %in% GenomeInfoDb::seqlevels(filterRegions))},
                              logical(1))
    if (all(sharedSeqlevels == FALSE)) {
      stop("The ", filterLabel, " shares no chromosome name with the region sets, please check the chromosome naming styles.", call. = FALSE)
    }

    #--------------------#
    # Filter the regions #
    #--------------------#
    if (verbose == TRUE) {message("Applying the ", filterLabel, " (", length(filterRegions), " regions):")}

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
        message("  ", names(filteredSets)[i], ": ", nAfter[i], "/", nBefore[i], " regions retained (",
                round(100 * (nBefore[i] - nAfter[i]) / max(nBefore[i], 1), 1), "% removed)")
      }
    }

    #-----------------------#
    # Handle the empty sets #
    #-----------------------#
    # A set emptied by the filter breaks any downstream comparison, so it is fatal unless stated otherwise
    if (any(nAfter == 0) & emptySets != "keep") {
      emptyNames <- names(filteredSets)[nAfter == 0]

      if (emptySets == "stop") {
        stop("The ", filterLabel, " left the following region sets without any region: ",
             paste(emptyNames, collapse = ", "),
             ". Set 'emptySets' to 'remove' or 'keep' to proceed.", call. = FALSE)
      } else {
        warning("Empty region sets removed: ", paste(emptyNames, collapse = ", "), ".", call. = FALSE)
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

    # A second list adds to the one already stored, replacing it would lose track of the first.
    # The greylist has a slot of its own and leaves the stored blacklist alone
    if (keepOverlapping == TRUE) {
      regionSet@whitelist <- .poolStoredList(storedList = regionSet@whitelist, newList = filterRegions, sharedOnly = TRUE)
    } else if (filterLabel == "blacklist") {
      regionSet@blacklist <- .poolStoredList(storedList = regionSet@blacklist, newList = filterRegions)
    } else {
      regionSet@greylist <- .poolStoredList(storedList = regionSet@greylist, newList = filterRegions)
      regionSet@parameters[[filterLabel]]$n.regions <- length(regionSet@greylist)
      regionSet@parameters[[filterLabel]]$covered.bp <- sum(as.numeric(BiocGenerics::width(regionSet@greylist)))
    }

    methods::validObject(regionSet)
    return(regionSet)
  } # END function




#' @title .applyCountingLists
#'
#' @description Applies a blacklist and a greylist to the region sets right before they are counted, for \code{countReads} and \code{countBigwig}. The regions overlapping a list are removed as \code{applyBlacklist} and \code{applyGreylist} remove them, the steps are recorded, and the sets left without any region are dropped.
#'
#' @param regionSet \code{RegionSetDE} object, \code{GRangesList} or named list of \code{GRanges}.
#' @param blacklist A \code{GRanges}, a path to a BED-like file, a data.frame, or a list of them, or \code{NULL}.
#' @param greylist Same forms as \code{blacklist}, or \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with \code{regionSet}, the filtered regions in the class they came in; \code{provenance}, the list returned by \code{.provenanceSlots} with the lists and the steps recorded, the \code{blacklist} and the \code{greylist} merged with the ones stored in the object; and \code{lists}, the \code{blacklist} and the \code{greylist} of the call as \code{GRanges}, \code{NULL} for a list that was not given.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocGenerics width
#' @importFrom methods is
#'
#' @keywords internal
#' @noRd

.applyCountingLists <-
  function(regionSet,
           blacklist = NULL,
           greylist = NULL,
           verbose = TRUE) {

    isRegionSet <- methods::is(regionSet, "RegionSetDE")
    provenance <- .provenanceSlots(regionSet)
    listRanges <- list(blacklist = NULL, greylist = NULL)

    setSizes <- function(x) {
      regionList <- if (methods::is(x, "RegionSetDE")) {as.list(x@regions)} else {as.list(x)}
      return(vapply(regionList, length, numeric(1)))
    }

    for (listLabel in c("blacklist", "greylist")) {
      listInput <- if (listLabel == "blacklist") {blacklist} else {greylist}

      if (is.null(listInput)) {next}

      # An empty list removes nothing, and it could not be loaded as a set of regions
      if (methods::is(listInput, "GRanges") && length(listInput) == 0) {
        if (isTRUE(verbose)) {message("The ", listLabel, " holds no region, nothing is removed.")}
        next
      }

      # The filter would name the unnamed sets by itself, and the counting has to refuse them as it always did
      if (is.null(names(setSizes(regionSet)))) {
        stop("All the region sets must be named.", call. = FALSE)
      }

      .checkListAssembly(listInput = listInput, genomeAssembly = provenance$genome.assembly, listLabel = listLabel, regionLabel = "regions")
      listRanges[[listLabel]] <- .loadExclusionRegions(excludeRegions = listInput, seqlevelsStyle = NULL, listLabel = listLabel)

      sizeBefore <- setSizes(regionSet)

      # The emptied sets are kept for now, they are handled once both lists have gone through
      regionSet <- .applyRegionFilter(regionSet = regionSet,
                                      filterSet = listRanges[[listLabel]],
                                      keepOverlapping = FALSE,
                                      overlapType = "any",
                                      minOverlapBp = 1,
                                      minOverlapFraction = 0,
                                      trimRegions = FALSE,
                                      ignoreStrand = TRUE,
                                      emptySets = "keep",
                                      verbose = verbose,
                                      filterLabel = listLabel)

      if (isRegionSet) {next}

      # Plain ranges carry no history: the step is written here, the way the RegionSetDE object writes its own
      sizeAfter <- setSizes(regionSet)

      provenance$filtering.log <- rbind(provenance$filtering.log,
                                        data.frame(step = listLabel,
                                                   region.set = names(sizeBefore),
                                                   n.before = sizeBefore,
                                                   n.after = sizeAfter,
                                                   n.removed = sizeBefore - sizeAfter,
                                                   row.names = NULL,
                                                   stringsAsFactors = FALSE))

      if (listLabel == "blacklist") {
        provenance$blacklist <- listRanges$blacklist
      } else {
        provenance$greylist <- listRanges$greylist
        provenance$parameters$greylist <- list(n.regions = length(listRanges$greylist),
                                               covered.bp = sum(as.numeric(BiocGenerics::width(listRanges$greylist))))
      }
    }

    if (is.null(listRanges$blacklist) & is.null(listRanges$greylist)) {
      return(list(regionSet = regionSet, provenance = provenance, lists = listRanges))
    }

    #-----------------------#
    # Handle the empty sets #
    #-----------------------#
    # A set without regions has no row to count, it leaves with a warning rather than stopping the whole counting
    sizeAfter <- setSizes(regionSet)

    if (all(sizeAfter == 0)) {
      stop("The blacklist and the greylist left no region to count.", call. = FALSE)
    }

    if (any(sizeAfter == 0)) {
      warning("Region sets left without any region by the blacklist or the greylist, and not counted: ",
              paste(names(sizeAfter)[sizeAfter == 0], collapse = ", "), ".", call. = FALSE)

      if (isRegionSet) {
        regionSet@regions <- regionSet@regions[sizeAfter > 0]
      } else {
        regionSet <- regionSet[sizeAfter > 0]
      }
    }

    if (isRegionSet) {
      provenance <- .provenanceSlots(regionSet)
    }

    return(list(regionSet = regionSet, provenance = provenance, lists = listRanges))
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
#' @details A blacklist applied to an object that already stores one adds to it. The \code{blacklist} slot then holds the two lists merged, and the \code{filtering.log} one step per call, so nothing is lost when an assay-specific list follows the ENCODE one. \code{\link{countReads}} reads the stored blacklist to leave its reads out of the counts and of the library sizes.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#' exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(
#'   GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'   splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
#'
#' # The example regions are shipped before the filtering, so this removes rows
#' regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)
#' regions
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{applyWhitelist}}, \code{\link{applyGreylist}}, \code{\link{loadBlacklist}}, \code{\link{countReads}}, \code{\link{loadRegions}}
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




#' @title applyGreylist
#'
#' @description Removes from every region set the regions overlapping a greylist, typically the one built from the input libraries by \code{\link{makeGreylist}}, or one exported by another tool as a BED-like file. It works as \code{\link{applyBlacklist}} does, but the list goes to the \code{greylist} slot and the step is recorded as a greylist, while the blacklist already stored in the object is left as it is.
#'
#' @param regionSet A \code{RegionSetDE} object, a \code{GRangesList}, a named list of \code{GRanges} or a single \code{GRanges}.
#' @param greylist \code{GRanges} returned by \code{\link{makeGreylist}}, string indicating the path to a BED-like file, or data.frame with the regions to exclude.
#' @param overlapType String indicating the type of overlap required to greylist a region, one among \code{"any"}, \code{"within"}, \code{"start"}, \code{"end"} or \code{"equal"}. Default: \code{"any"}.
#' @param minOverlapBp Numeric value indicating the minimum number of bases that must overlap the greylist for a region to be removed. Default: \code{1}.
#' @param minOverlapFraction Numeric value between 0 and 1 indicating the minimum fraction of a region that must overlap the greylist for it to be removed. Default: \code{0}, any overlap is sufficient.
#' @param trimRegions Logical value to indicate whether the greylisted portion must be subtracted from the regions instead of removing them entirely. Notice that trimming collapses the regions overlapping each other within the same set. Default: \code{FALSE}.
#' @param ignoreStrand Logical value to indicate whether the strand must be ignored when computing the overlaps. Default: \code{TRUE}.
#' @param emptySets String indicating how to handle the sets left without any region, one among \code{"stop"}, \code{"remove"} or \code{"keep"}. Default: \code{"stop"}.
#' @param verbose Logical value to indicate whether the filtering messages must be printed. Default: \code{TRUE}.
#'
#' @return An object of the same class as \code{regionSet}. For a \code{RegionSetDE} object the greylist is stored in the \code{greylist} slot, merged with the one already there, the step is added to the \code{filtering.log} as \code{"greylist"}, and the number of regions of the stored greylist and the bases they cover go to \code{parameters$greylist}.
#'
#' @details A greylist removes what the inputs flag as artefacts, and for broad marks some of what it flags is genuine signal: heterochromatin marks such as H3K9me3 sit on satellites and repeats, where inputs pile up too. \code{trimRegions = TRUE} cuts the greylisted stretch out of a broad domain and keeps the rest of it, and \code{minOverlapFraction} removes only the regions mostly covered by the greylist. The messages report how many regions each set keeps, which is the number to look at before going further.
#'
#' \code{\link{countReads}} reads the stored greylist to leave its reads out of the counts and of the library sizes, as it does with the blacklist. The \code{greylist} argument of \code{\link{countReads}} does the work of this function at the counting step, so the two need not be combined.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'                             splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
#'
#' # The exclusion list shipped with the package stands in for a greylist
#' greylist <- loadExampleData("exclusionRegions", verbose = FALSE)
#'
#' greylisted <- applyGreylist(regions, greylist = greylist)
#' filteringLog(greylisted)
#'
#' # Broad domains lose the greylisted stretch only
#' trimmed <- applyGreylist(regions, greylist = greylist, trimRegions = TRUE, verbose = FALSE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{makeGreylist}}, \code{\link{applyBlacklist}}, \code{\link{applyWhitelist}}, \code{\link{countReads}}
#'
#' @export applyGreylist

applyGreylist <-
  function(regionSet,
           greylist,
           overlapType = "any",
           minOverlapBp = 1,
           minOverlapFraction = 0,
           trimRegions = FALSE,
           ignoreStrand = TRUE,
           emptySets = "stop",
           verbose = TRUE) {

    return(.applyRegionFilter(regionSet = regionSet,
                              filterSet = greylist,
                              keepOverlapping = FALSE,
                              overlapType = overlapType,
                              minOverlapBp = minOverlapBp,
                              minOverlapFraction = minOverlapFraction,
                              trimRegions = trimRegions,
                              ignoreStrand = ignoreStrand,
                              emptySets = emptySets,
                              verbose = verbose,
                              filterLabel = "greylist"))
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
#' @details A second whitelist restricts the first one: the regions left overlap both lists. The \code{whitelist} slot then holds the stretches the two lists have in common, and the \code{filtering.log} one step per call.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(
#'   GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'   splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
#'
#' # Restrict the analysis to the first half of the chromosome
#' whitelist <- GenomicRanges::GRanges(
#'   seqnames = "chr12",
#'   ranges = IRanges::IRanges(start = 1, end = 25e6))
#'
#' regions <- applyWhitelist(regions, whitelist = whitelist, verbose = FALSE)
#' regions
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
