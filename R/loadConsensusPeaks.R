# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title loadConsensusPeaks
#'
#' @description Builds the regions of an analysis from the peaks called on each sample. The peaks are combined into a consensus within each group of samples through \code{consensusRegions}, and the group consensus are pooled into a total one, which becomes the region set counted downstream. Regions supplied by the user can split the total consensus into sets, or take its place.
#'
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, with at least the \code{sample} and \code{peaks} columns.
#' @param groupBy String with the column of the sample sheet defining the groups, for instance \code{"condition"}. Default: \code{NULL}, all the samples form a single group.
#' @param excludeRegions Regions whose peaks must be dropped before any consensus is built, typically a blacklist and the greylist returned by \code{\link{makeGreylist}}. Either a \code{GRanges}, a path to a BED-like file, a data.frame, or a list of them. Default: \code{NULL}.
#' @param regionSets Regions of interest of the user, in any form accepted by \code{\link{loadRegions}}: a named list of paths, \code{GRanges} or data.frames. Default: \code{NULL}, the total consensus is the only set.
#' @param regionMode String indicating what \code{regionSets} do, either \code{"split"}, where every region of the total consensus goes to the first set it overlaps, or \code{"replace"}, where the sets of the user are the regions and the peaks only annotate them. Ignored without \code{regionSets}. Default: \code{"split"}.
#' @param unassignedSet String with the name of the set collecting the consensus regions that overlap none of \code{regionSets} in \code{"split"} mode. \code{NULL} drops them. Default: \code{"other"}.
#' @param seqlevelsStyle String indicating the chromosome naming style, one among \code{"UCSC"}, \code{"Ensembl"} or \code{"NCBI"}, or \code{NULL} to keep the names as they are. Default: \code{"UCSC"}.
#' @param genomeAssembly String indicating the genome assembly to store with the regions, e.g. \code{"hg38"}. Default: \code{NULL}.
#' @param nThreads Number of threads used by the consensus. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#' @param ... Further arguments passed to \code{consensusRegions::runConsensus}, for instance \code{minReplicates}, \code{combinedThreshold}, \code{calibrate} or \code{weightMethod}.
#'
#' @return A \code{RegionSetDE} object. Every region carries \code{peak.<group>}, telling whether the consensus of that group overlaps it, \code{peak.groups}, the number of groups with a consensus peak on it, and \code{peak.samples}, the number of samples with a peak on it. The \code{consensus} slot, read with \code{\link{consensusData}}, keeps the consensus of every group, the \code{consensusRegions} object behind it, the peaks of every sample after the exclusion, the total consensus, the table of the samples and the sample sheet itself, which \code{\link{countReads}} and \code{\link{countBigwig}} read when no file is given to them.
#'
#' @details The consensus is built group by group and then pooled, rather than once over all the samples. A single requirement over all the libraries, such as a peak in at least two of them, favours the larger group, while the same requirement applied within each group treats them alike, and the union keeps the regions found in one group only. A group with a single sample has no replicate to agree with, and its peaks stand for the group as they are.
#'
#' The peaks are read by \code{consensusRegions}, which keeps the standard chromosomes only, as \code{GenomeInfoDb::keepStandardChromosomes} defines them: scaffolds, patches and unplaced contigs are dropped, and so are chromosomes whose names it does not recognise.
#'
#' The peaks overlapping \code{excludeRegions} are removed before the consensus, not after it. An artefact lying next to a genuine peak would otherwise merge with it, and removing the merged region afterwards would take the genuine peak away as well.
#'
#' With \code{regionMode = "split"} the consensus regions are assigned to the sets of \code{regionSets} in the order they are given, a region overlapping two sets going to the first one. The sets then describe classes of peaks, such as promoter and distal ones, and \code{\link{testRegionSets}} can compare them. With \code{regionMode = "replace"} the counted regions are those of \code{regionSets}, and the peaks only tell which of them are occupied in each group.
#'
#' The occupancy columns describe where the peaks were called, and they are a poor basis for a set of regions to test. A set of the regions found in one group only was selected on the signal of that group, so a test of its change between the groups answers a question already settled by the selection.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'   peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
#'                           pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE)
#'
#'   # Three replicates of one condition and two of another, peaks only
#'   sampleSheet <- data.frame(sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
#'                             bam = c("A_1.bam", "A_2.bam", "A_3.bam", "B_1.bam", "B_2.bam"),
#'                             peaks = peakFiles[c(1, 2, 3, 1, 2)],
#'                             condition = c("A", "A", "A", "B", "B"))
#'
#'   sampleSheet <- loadSampleSheet(sampleSheet, checkFiles = FALSE, verbose = FALSE)
#'
#'   regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition")
#'   regions
#'
#'   head(regionRanges(regions)$consensus, 3)
#'
#'   # The same peaks split by regions of interest of the user
#'   firstHalf <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 5e5))
#'   splitRegions <- loadConsensusPeaks(sampleSheet, groupBy = "condition",
#'                                      regionSets = list(firstHalf = firstHalf),
#'                                      verbose = FALSE)
#'   lengths(regionRanges(splitRegions))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotPeakUpset}}, \code{\link{consensusData}}, \code{\link{loadSampleSheet}}, \code{\link{makeGreylist}}
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges reduce overlapsAny
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is validObject
#'
#' @export loadConsensusPeaks

loadConsensusPeaks <-
  function(sampleSheet,
           groupBy = NULL,
           excludeRegions = NULL,
           regionSets = NULL,
           regionMode = "split",
           unassignedSet = "other",
           seqlevelsStyle = "UCSC",
           genomeAssembly = NULL,
           nThreads = 1,
           verbose = TRUE,
           ...) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!requireNamespace("consensusRegions", quietly = TRUE)) {
      stop("The 'consensusRegions' package is needed to build the consensus peaks.", call. = FALSE)
    }

    if (is.character(sampleSheet) & length(sampleSheet) == 1) {
      sampleSheet <- loadSampleSheet(sampleSheet, verbose = verbose)
    }

    if (!is.data.frame(sampleSheet) || !all(c("sample", "peaks") %in% colnames(sampleSheet))) {
      stop("The 'sampleSheet' parameter must be a table with the 'sample' and 'peaks' columns, as returned by loadSampleSheet().", call. = FALSE)
    }

    if (!is.null(groupBy)) {
      if (!is.character(groupBy) | length(groupBy) != 1 || !(groupBy %in% colnames(sampleSheet))) {
        stop("The 'groupBy' parameter must name a column of the sample sheet.", call. = FALSE)
      }
    }

    regionMode <- tolower(as.character(regionMode[1]))
    if (!(regionMode %in% c("split", "replace"))) {
      stop("The 'regionMode' parameter must be either 'split' or 'replace'.", call. = FALSE)
    }

    if (!is.null(unassignedSet) && (!is.character(unassignedSet) | length(unassignedSet) != 1)) {
      stop("The 'unassignedSet' parameter must be a single string, or NULL.", call. = FALSE)
    }

    consensusArguments <- list(...)
    reservedArguments <- intersect(names(consensusArguments), c("peaks", "sampleNames"))
    if (length(reservedArguments) > 0) {
      stop("The following arguments of runConsensus are set from the sample sheet and cannot be given: ",
           paste(reservedArguments, collapse = ", "), ".", call. = FALSE)
    }

    #-------------------------------#
    # Samples and groups            #
    #-------------------------------#
    groupValues <- if (is.null(groupBy)) {rep("all", nrow(sampleSheet))} else {as.character(sampleSheet[[groupBy]])}

    sampleTable <- data.frame(sample = as.character(sampleSheet$sample),
                              group = groupValues,
                              peak.file = as.character(sampleSheet$peaks),
                              bam = if ("bam" %in% colnames(sampleSheet)) {as.character(sampleSheet$bam)} else {NA_character_},
                              stringsAsFactors = FALSE)

    if (anyNA(sampleTable$group)) {
      stop("The following samples have no value in the '", groupBy, "' column: ",
           paste(sampleTable$sample[is.na(sampleTable$group)], collapse = ", "), ".", call. = FALSE)
    }

    withoutPeaks <- dplyr::filter(sampleTable, is.na(.data$peak.file))
    if (isTRUE(verbose) & nrow(withoutPeaks) > 0) {
      message("The following samples have no peak file and stay out of the consensus: ", paste(withoutPeaks$sample, collapse = ", "), ".")
    }

    sampleTable <- dplyr::filter(sampleTable, !is.na(.data$peak.file))

    if (nrow(sampleTable) < 2) {
      stop("At least two samples with a peak file are needed to build a consensus.", call. = FALSE)
    }

    # The levels of a factor are the order the user chose, otherwise the groups follow the sheet
    groupLevels <- if (!is.null(groupBy) && is.factor(sampleSheet[[groupBy]])) {
      intersect(levels(sampleSheet[[groupBy]]), sampleTable$group)
    } else {
      unique(sampleTable$group)
    }

    # Recentring works within a group, and the union merges the overlapping windows of different groups back
    if (isTRUE(consensusArguments$recentre) & length(groupLevels) > 1) {
      warning("Recentred regions of different groups overlap and are merged back by the union, so the total consensus will not have a fixed width.", call. = FALSE)
    }

    #-------------------------------#
    # Read and clean the peaks      #
    #-------------------------------#
    peakList <- consensusRegions::readPeakSets(peaks = sampleTable$peak.file,
                                               sampleNames = sampleTable$sample,
                                               seqlevelsStyle = seqlevelsStyle,
                                               verbose = FALSE)

    sampleTable$n.peaks <- as.integer(lengths(peakList))
    sampleTable$n.excluded <- 0L
    excludedRanges <- NULL

    # An artefact next to a genuine peak would merge with it, so the exclusion comes before the consensus
    if (!is.null(excludeRegions)) {
      excludedRanges <- .loadExclusionRegions(excludeRegions = excludeRegions, seqlevelsStyle = seqlevelsStyle)

      peakList <- GenomicRanges::GRangesList(lapply(as.list(peakList),
                                                    function(peakRanges) {
                                                      peakRanges[!IRanges::overlapsAny(peakRanges, excludedRanges, ignore.strand = TRUE)]
                                                    }))

      sampleTable$n.excluded <- sampleTable$n.peaks - as.integer(lengths(peakList))

      if (isTRUE(verbose)) {
        message("Peaks overlapping the excluded regions removed: ",
                paste(sampleTable$sample, paste(sampleTable$n.excluded, sampleTable$n.peaks, sep = "/"), collapse = ", "), ".")
      }
    }

    #-------------------------------#
    # Consensus of every group      #
    #-------------------------------#
    groupConsensus <- list()
    consensusObjects <- list()

    for (groupName in groupLevels) {
      groupSamples <- sampleTable$sample[sampleTable$group == groupName]

      # A single sample has no replicate to agree with, its peaks stand for the group
      if (length(groupSamples) == 1) {
        groupConsensus[[groupName]] <- IRanges::reduce(peakList[[groupSamples]], ignore.strand = TRUE)
        consensusObjects[groupName] <- list(NULL)

        if (isTRUE(verbose)) {
          message("Group ", groupName, ": a single sample, its ", length(groupConsensus[[groupName]]), " peaks are taken as they are.")
        }
        next
      }

      groupArguments <- consensusArguments
      if (is.null(groupArguments$BPPARAM)) {groupArguments$BPPARAM <- nThreads}
      if (is.null(groupArguments$verbose)) {groupArguments$verbose <- FALSE}

      # runConsensus renames the chromosomes to UCSC on its own, which would leave the consensus in
      # one style and the peaks it came from in another, and every overlap between the two empty
      if (!("seqlevelsStyle" %in% names(groupArguments))) {
        groupArguments["seqlevelsStyle"] <- list(seqlevelsStyle)
      }

      # Weights from the libraries need the BAM files, which the sheet already lists
      weightMethod <- groupArguments$weightMethod
      if (!is.null(weightMethod) && weightMethod[1] %in% c("frip", "librarySize") &
          is.null(groupArguments$bamFiles) & is.null(groupArguments$frip) & is.null(groupArguments$librarySize)) {
        groupArguments$bamFiles <- sampleTable$bam[match(groupSamples, sampleTable$sample)]
      }

      consensusObject <- do.call(what = consensusRegions::runConsensus,
                                 args = c(list(peaks = peakList[groupSamples], sampleNames = groupSamples), groupArguments))

      consensusRanges <- consensusRegions::consensusRanges(consensusObject)

      # The statistics of one group mean nothing for the others, only the coordinates are pooled
      S4Vectors::mcols(consensusRanges) <- NULL
      groupConsensus[[groupName]] <- consensusRanges
      consensusObjects[[groupName]] <- consensusObject

      if (isTRUE(verbose)) {
        message("Group ", groupName, ": ", length(groupSamples), " samples, ", length(consensusRanges), " consensus regions.")
      }
    }

    #-------------------------------#
    # Total consensus               #
    #-------------------------------#
    totalConsensus <- IRanges::reduce(unlist(GenomicRanges::GRangesList(unname(groupConsensus)), use.names = FALSE), ignore.strand = TRUE)

    if (length(totalConsensus) == 0) {
      stop("No consensus region was found in any group.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message("Total consensus: ", length(totalConsensus), " regions from ", length(groupLevels), " groups.")
    }

    #-------------------------------#
    # Regions of the user           #
    #-------------------------------#
    analysisMode <- "consensus"
    regionList <- list(consensus = totalConsensus)

    if (!is.null(regionSets)) {
      userList <- loadRegions(regions = regionSets,
                              keepMetadata = TRUE,
                              seqlevelsStyle = seqlevelsStyle,
                              genomeAssembly = genomeAssembly,
                              outputFormat = "list",
                              verbose = FALSE)

      analysisMode <- regionMode

      regionList <- if (regionMode == "split") {
        .splitConsensus(totalConsensus = totalConsensus, userList = userList, unassignedSet = unassignedSet, verbose = verbose)
      } else {
        userList
      }
    }

    #-------------------------------#
    # Occupancy of every region     #
    #-------------------------------#
    regionList <- lapply(regionList,
                         function(regionRanges) {
                           occupancyTable <- .peakOccupancy(regionRanges = regionRanges,
                                                            groupConsensus = groupConsensus,
                                                            peakList = peakList)
                           S4Vectors::mcols(regionRanges) <- cbind(S4Vectors::mcols(regionRanges), occupancyTable)
                           return(regionRanges)
                         })

    regionSet <- loadRegions(regions = regionList,
                             keepMetadata = TRUE,
                             seqlevelsStyle = seqlevelsStyle,
                             genomeAssembly = genomeAssembly,
                             verbose = FALSE)

    #-------------------------------#
    # Keep the consensus data       #
    #-------------------------------#
    regionSet@consensus <- list(groups = groupConsensus,
                                objects = consensusObjects,
                                total = totalConsensus,
                                peaks = peakList,
                                samples = sampleTable,
                                sheet = sampleSheet,
                                groupBy = groupBy,
                                mode = analysisMode,
                                excluded = excludedRanges)

    # Only plain values go to the parameters, a BiocParallel object has no place in an exported record
    regionSet@parameters$loadConsensusPeaks <- list(groupBy = groupBy,
                                                    regionMode = analysisMode,
                                                    unassignedSet = unassignedSet,
                                                    seqlevelsStyle = seqlevelsStyle,
                                                    n.samples = nrow(sampleTable),
                                                    n.groups = length(groupLevels),
                                                    n.excluded.regions = if (is.null(excludedRanges)) {0L} else {length(excludedRanges)},
                                                    consensusArguments = Filter(is.atomic, consensusArguments))

    methods::validObject(regionSet)
    return(regionSet)
  } # END function




#' @title .loadExclusionRegions
#'
#' @description Loads one or several exclusion lists and pools them into a single set of ranges, in the chromosome style of the peaks.
#'
#' @param excludeRegions A \code{GRanges}, a path, a data.frame, or a list of them.
#' @param seqlevelsStyle String with the chromosome naming style, or \code{NULL}.
#'
#' @return A \code{GRanges} without overlaps.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges reduce
#' @importFrom methods is
#'
#' @keywords internal

.loadExclusionRegions <-
  function(excludeRegions,
           seqlevelsStyle = "UCSC") {

    # A single list of regions is wrapped, a list of lists is taken as it is
    if (methods::is(excludeRegions, "GRanges") | is.data.frame(excludeRegions) | is.character(excludeRegions)) {
      excludeRegions <- list(excludeRegions)
    }

    if (!is.list(excludeRegions)) {
      stop("The 'excludeRegions' parameter must be a GRanges, a path, a data.frame, or a list of them.", call. = FALSE)
    }

    names(excludeRegions) <- paste0("exclusion_", seq_along(excludeRegions))

    exclusionList <- loadRegions(regions = excludeRegions,
                                 keepMetadata = FALSE,
                                 duplicatedSets = "remove",
                                 seqlevelsStyle = seqlevelsStyle,
                                 outputFormat = "list",
                                 verbose = FALSE)

    return(IRanges::reduce(unlist(GenomicRanges::GRangesList(unname(exclusionList)), use.names = FALSE), ignore.strand = TRUE))
  } # END function




#' @title .splitConsensus
#'
#' @description Assigns every region of the total consensus to the first set of the user it overlaps, the others to the unassigned set.
#'
#' @param totalConsensus \code{GRanges} with the total consensus.
#' @param userList Named list of \code{GRanges} with the regions of the user.
#' @param unassignedSet String with the name of the set collecting the regions overlapping no set, or \code{NULL} to drop them.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A named list of \code{GRanges}, the empty sets left out.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges overlapsAny
#'
#' @keywords internal

.splitConsensus <-
  function(totalConsensus,
           userList,
           unassignedSet = "other",
           verbose = TRUE) {

    if (!is.null(unassignedSet) && unassignedSet %in% names(userList)) {
      stop("The 'unassignedSet' name '", unassignedSet, "' is already used by one of the region sets.", call. = FALSE)
    }

    # The sets are visited in the order given, so a region overlapping two of them goes to the first
    setAssignment <- rep(NA_character_, length(totalConsensus))

    for (setName in names(userList)) {
      overlapping <- is.na(setAssignment) & IRanges::overlapsAny(totalConsensus, userList[[setName]], ignore.strand = TRUE)
      setAssignment[overlapping] <- setName
    }

    if (!is.null(unassignedSet)) {
      setAssignment[is.na(setAssignment)] <- unassignedSet
    }

    setLevels <- c(names(userList), unassignedSet)
    splitList <- lapply(setLevels, function(setName) {totalConsensus[which(setAssignment == setName)]})
    names(splitList) <- setLevels

    emptySets <- setLevels[lengths(splitList) == 0]
    if (isTRUE(verbose) & length(emptySets) > 0) {
      message("No consensus region falls in the following sets, which are left out: ", paste(emptySets, collapse = ", "), ".")
    }

    splitList <- splitList[lengths(splitList) > 0]

    if (length(splitList) == 0) {
      stop("No consensus region overlaps the region sets, and 'unassignedSet' is NULL.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message("Consensus regions per set: ", paste(names(splitList), lengths(splitList), sep = " ", collapse = ", "), ".")
    }

    return(splitList)
  } # END function




#' @title .peakOccupancy
#'
#' @description Tells, for every region, which groups have a consensus peak on it and how many samples have a peak on it.
#'
#' @param regionRanges \code{GRanges} with the regions.
#' @param groupConsensus Named list of \code{GRanges}, the consensus of every group.
#' @param peakList Named \code{GRangesList} with the peaks of every sample.
#'
#' @return A \code{DataFrame} with one logical \code{peak.<group>} column per group, \code{peak.groups} and \code{peak.samples}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom utils head
#'
#' @keywords internal

.peakOccupancy <-
  function(regionRanges,
           groupConsensus,
           peakList) {

    regionNumber <- length(regionRanges)

    # Two naming styles meeting here would give an empty overlap for every sample rather than an error
    sharedSeqlevels <- intersect(GenomeInfoDb::seqlevels(regionRanges),
                                 unique(unlist(lapply(as.list(peakList), GenomeInfoDb::seqlevels), use.names = FALSE)))

    if (length(sharedSeqlevels) == 0) {
      stop("The regions and the peaks of the samples share no chromosome, their naming styles differ: ",
           paste(utils::head(GenomeInfoDb::seqlevels(regionRanges), 3), collapse = ", "), " against ",
           paste(utils::head(GenomeInfoDb::seqlevels(peakList[[1]]), 3), collapse = ", "), ".", call. = FALSE)
    }

    # matrix() keeps the shape when there is a single region, which vapply would flatten
    groupMatrix <- matrix(vapply(groupConsensus,
                                 function(groupRanges) {IRanges::overlapsAny(regionRanges, groupRanges, ignore.strand = TRUE)},
                                 logical(regionNumber)),
                          nrow = regionNumber)

    sampleMatrix <- matrix(vapply(as.list(peakList),
                                  function(peakRanges) {IRanges::overlapsAny(regionRanges, peakRanges, ignore.strand = TRUE)},
                                  logical(regionNumber)),
                           nrow = regionNumber)

    colnames(groupMatrix) <- paste0("peak.", make.names(names(groupConsensus)))

    occupancyTable <- S4Vectors::DataFrame(groupMatrix, check.names = FALSE)
    occupancyTable$peak.groups <- as.integer(rowSums(groupMatrix))
    occupancyTable$peak.samples <- as.integer(rowSums(sampleMatrix))

    return(occupancyTable)
  } # END function




#' @title consensusData
#'
#' @description Returns the consensus data kept by a \code{RegionSetDE} object built from peaks with \code{\link{loadConsensusPeaks}}.
#'
#' @param object \code{RegionSetDE} object.
#'
#' @return A list with \code{groups}, the consensus regions of every group; \code{objects}, the \code{consensusRegions} object behind each of them (\code{NULL} for the groups with a single sample); \code{total}, the pooled consensus; \code{peaks}, the peaks of every sample after the exclusion; \code{samples}, a data.frame with the group, the peak file and the number of peaks kept and excluded for every sample; \code{sheet}, the sample sheet the consensus was built from; \code{groupBy}; \code{mode}, one among \code{"consensus"}, \code{"split"} and \code{"replace"}; and \code{excluded}, the regions the peaks were cleaned against.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'   peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
#'                           pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE)
#'
#'   sampleSheet <- loadSampleSheet(data.frame(sample = c("A_1", "A_2", "B_1", "B_2"),
#'                                             bam = c("A_1.bam", "A_2.bam", "B_1.bam", "B_2.bam"),
#'                                             peaks = peakFiles[c(1, 2, 2, 3)],
#'                                             condition = c("A", "A", "B", "B")),
#'                                  checkFiles = FALSE, verbose = FALSE)
#'
#'   regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)
#'
#'   consensusData(regions)$samples
#'   lengths(consensusData(regions)$groups)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadConsensusPeaks}}, \code{\link{plotPeakUpset}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export

setGeneric(name = "consensusData", def = function(object) {standardGeneric("consensusData")})


#' @rdname consensusData
#' @export

setMethod(f = "consensusData",
          signature = "RegionSetDE",
          definition = function(object) {
            if (length(object@consensus) == 0) {
              stop("The regions were not built from peaks, there is no consensus to return.", call. = FALSE)
            }
            return(object@consensus)
          })
