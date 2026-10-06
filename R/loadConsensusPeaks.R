# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title loadConsensusPeaks
#'
#' @description Builds the regions of an analysis from the peaks called on each sample. The peaks are combined into a consensus within each group of samples through \code{consensusRegions}, and the group consensus are pooled into a total one, which becomes the region set counted downstream. Regions supplied by the user can split the total consensus into sets, or take its place.
#'
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, with at least the \code{sample} and \code{peaks} columns.
#' @param groupBy String with the column of the sample sheet defining the groups, for instance \code{"condition"}. Default: \code{NULL}, all the samples form a single group.
#' @param blacklist Regions of the assembly whose peaks must be dropped before any consensus is built, typically the list returned by \code{\link{loadBlacklist}}. Either a \code{GRanges}, a path to a BED-like file, a data.frame, or a list of them, which are pooled. The list is stored in the \code{blacklist} slot of the object, as \code{\link{applyBlacklist}} does. Default: \code{NULL}.
#' @param greylist Regions of this experiment whose peaks must be dropped before any consensus is built, typically the list built from the inputs by \code{\link{makeGreylist}}. Accepts the same forms as \code{blacklist}. Default: \code{NULL}.
#' @param regionSets Regions of interest of the user, in any form accepted by \code{\link{loadRegions}}: a named list of paths, \code{GRanges} or data.frames. Default: \code{NULL}, the total consensus is the only set.
#' @param regionMode String indicating what \code{regionSets} do, either \code{"split"}, where every region of the total consensus goes to the first set it overlaps, or \code{"replace"}, where the sets of the user are the regions and the peaks only annotate them. Ignored without \code{regionSets}. Default: \code{"split"}.
#' @param unassignedSet String with the name of the set collecting the consensus regions that overlap none of \code{regionSets} in \code{"split"} mode. \code{NULL} drops them. Default: \code{"other"}.
#' @param seqlevelsStyle String indicating the chromosome naming style, one among \code{"UCSC"}, \code{"Ensembl"} or \code{"NCBI"}, or \code{NULL} to keep the names as they are. Default: \code{"UCSC"}.
#' @param genomeAssembly String indicating the genome assembly to store with the regions, e.g. \code{"hg38"}. Default: \code{NULL}.
#' @param nThreads Number of threads. The groups with more than one sample are built side by side, one worker per group, and the threads left over go to the calibration of the threshold inside each group when \code{calibrate = TRUE} is passed on to \code{consensusRegions::runConsensus}. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#' @param ... Further arguments passed to \code{consensusRegions::runConsensus}, for instance \code{minReplicates}, \code{combinedThreshold}, \code{calibrate} or \code{weightMethod}.
#'
#' @return A \code{RegionSetDE} object. Every region carries \code{peak.<group>}, telling whether the consensus of that group overlaps it, \code{peak.groups}, the number of groups with a consensus peak on it, and \code{peak.samples}, the number of samples with a peak on it. The \code{consensus} slot, read with \code{\link{consensusData}}, keeps the consensus of every group, the \code{consensusRegions} object behind it, the peaks of every sample after the blacklist and the greylist, the peaks they removed, the two lists, the total consensus, the table of the samples and the sample sheet itself, which \code{\link{countReads}} and \code{\link{countBigwig}} read when no file is given to them. The blacklist goes to the \code{blacklist} slot, the size of the greylist to \code{parameters$greylist}, and the peaks each of them removed from every sample to the \code{filtering.log}, read with \code{\link{filteringLog}}.
#'
#' @details The consensus is built group by group and then pooled, rather than once over all the samples. A single requirement over all the libraries, such as a peak in at least two of them, favours the larger group, while the same requirement applied within each group treats them alike, and the union keeps the regions found in one group only. A group with a single sample has no replicate to agree with, and its peaks stand for the group as they are.
#'
#' The peaks are read by \code{consensusRegions}, which keeps the standard chromosomes only, as \code{GenomeInfoDb::keepStandardChromosomes} defines them: scaffolds, patches and unplaced contigs are dropped, and so are chromosomes whose names it does not recognise.
#'
#' The peaks overlapping \code{blacklist} and \code{greylist} are removed before the consensus, not after it. An artefact lying next to a genuine peak would otherwise merge with it, and removing the merged region afterwards would take the genuine peak away as well. The blacklist is applied first and the greylist to the peaks left, so a peak lying on both is counted as blacklisted. The two lists are kept apart because they say different things: a blacklist describes the assembly and is the same for every experiment, a greylist describes the inputs of this one. Both are recorded the way \code{\link{applyBlacklist}} and \code{\link{applyGreylist}} record them, so the counts, the fit and the results built from the object carry the blacklist with them and \code{\link{exportResults}} writes it down. \code{consensusData(x)$removed} holds the peaks taken out, with the sample they came from and the list that removed them, and \code{consensusData(x)$samples} their number per sample. A blacklist built for another assembly than \code{genomeAssembly} is refused, as \code{\link{applyBlacklist}} refuses it.
#'
#' With \code{regionMode = "split"} the consensus regions are assigned to the sets of \code{regionSets} in the order they are given, a region overlapping two sets going to the first one. The sets then describe classes of peaks, such as promoter and distal ones, and \code{\link{testRegionSets}} can compare them. With \code{regionMode = "replace"} the counted regions are those of \code{regionSets}, and the peaks only tell which of them are occupied in each group.
#'
#' The occupancy columns describe where the peaks were called, and they are a poor basis for a set of regions to test. A set of the regions found in one group only was selected on the signal of that group, so a test of its change between the groups answers a question already settled by the selection.
#'
#' The consensus of a group does not depend on the other groups, so with \code{nThreads} above one the groups run in parallel. \code{consensusRegions} builds a consensus on a single thread and only parallelises the permutations of the calibration, so without \code{calibrate} the groups are the one place where more threads save time: three groups on three threads take about as long as the largest of them. With \code{calibrate = TRUE} the threads are shared, \code{nThreads} divided by the number of groups for the permutations of each group. A \code{BPPARAM} passed in \code{...} is handed to every group as it is and the groups then run one after the other, so that two levels of workers are never stacked on each other by accident. The random numbers of the calibration come from \code{BiocParallel}, which gives each group a stream of its own, so a \code{set.seed()} before the call returns the same consensus whatever the number of threads.
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
#' @seealso \code{\link{plotPeakUpset}}, \code{\link{consensusData}}, \code{\link{loadSampleSheet}}, \code{\link{loadBlacklist}}, \code{\link{makeGreylist}}
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges reduce overlapsAny IRanges
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocParallel bplapply
#' @importFrom methods is validObject
#'
#' @export loadConsensusPeaks

loadConsensusPeaks <-
  function(sampleSheet,
           groupBy = NULL,
           blacklist = NULL,
           greylist = NULL,
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

    if ("excludeRegions" %in% names(consensusArguments)) {
      stop("The 'excludeRegions' parameter has been replaced by 'blacklist' and 'greylist', which take the two lists separately.", call. = FALSE)
    }

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

    #-------------------------------#
    # Blacklist, then greylist      #
    #-------------------------------#
    # An artefact next to a genuine peak would merge with it, so both lists act before the consensus.
    # The greylist sees only the peaks the blacklist left, so a peak on both is counted once, as blacklisted.
    listRanges <- list(blacklist = NULL, greylist = NULL)
    removedList <- list()
    filteringSteps <- list()

    for (listLabel in c("blacklist", "greylist")) {
      listInput <- if (listLabel == "blacklist") {blacklist} else {greylist}
      sampleTable[[paste0("n.", listLabel)]] <- 0L

      if (is.null(listInput)) {next}

      .checkListAssembly(listInput = listInput, genomeAssembly = genomeAssembly, listLabel = listLabel)
      listRanges[[listLabel]] <- .loadExclusionRegions(excludeRegions = listInput, seqlevelsStyle = seqlevelsStyle, listLabel = listLabel)

      # With seqlevelsStyle = NULL the peaks keep the names of their files, and the list has to follow them
      listRanges[[listLabel]] <- tryCatch(expr = .matchSeqlevels(x = listRanges[[listLabel]],
                                                                 targetSeqlevels = GenomeInfoDb::seqlevels(peakList),
                                                                 verbose = FALSE),
                                          error = function(e) {return(listRanges[[listLabel]])})

      peaksBefore <- as.integer(lengths(peakList))
      overlapList <- lapply(as.list(peakList), function(peakRanges) {IRanges::overlapsAny(peakRanges, listRanges[[listLabel]], ignore.strand = TRUE)})

      # The peaks taken out are kept, with the sample they came from and the list that removed them
      removedList[[listLabel]] <- lapply(names(peakList), function(sampleName) {
        removedPeaks <- peakList[[sampleName]][overlapList[[sampleName]]]
        S4Vectors::mcols(removedPeaks)$sample <- rep(sampleName, length(removedPeaks))
        S4Vectors::mcols(removedPeaks)$removed.by <- rep(listLabel, length(removedPeaks))
        return(removedPeaks)
      })

      peakList <- GenomicRanges::GRangesList(lapply(names(peakList), function(sampleName) {peakList[[sampleName]][!overlapList[[sampleName]]]}))
      names(peakList) <- sampleTable$sample
      peaksAfter <- as.integer(lengths(peakList))

      sampleTable[[paste0("n.", listLabel)]] <- peaksBefore - peaksAfter

      filteringSteps[[listLabel]] <- data.frame(step = listLabel,
                                                region.set = paste(sampleTable$sample, "peaks"),
                                                n.before = peaksBefore,
                                                n.after = peaksAfter,
                                                n.removed = peaksBefore - peaksAfter,
                                                stringsAsFactors = FALSE)

      if (isTRUE(verbose)) {
        message("Peaks overlapping the ", listLabel, " (", length(listRanges[[listLabel]]), " regions) removed: ",
                paste(sampleTable$sample, paste(peaksBefore - peaksAfter, peaksBefore, sep = "/"), collapse = ", "), ".")
      }
    }

    sampleTable$n.excluded <- sampleTable$n.blacklist + sampleTable$n.greylist

    removedPeaks <- if (length(removedList) == 0) {
      GenomicRanges::GRanges()
    } else {
      unlist(GenomicRanges::GRangesList(unlist(removedList, recursive = FALSE, use.names = FALSE)), use.names = FALSE)
    }

    #-------------------------------#
    # Consensus of every group      #
    #-------------------------------#
    groupConsensus <- list()
    consensusObjects <- list()

    groupSizes <- vapply(groupLevels, function(groupName) {sum(sampleTable$group == groupName)}, integer(1))
    multiGroups <- groupLevels[groupSizes > 1]

    # A single sample has no replicate to agree with, its peaks stand for the group
    for (groupName in groupLevels[groupSizes == 1]) {
      groupSamples <- sampleTable$sample[sampleTable$group == groupName]
      groupConsensus[[groupName]] <- IRanges::reduce(peakList[[groupSamples]], ignore.strand = TRUE)
      consensusObjects[groupName] <- list(NULL)
    }

    # How the threads are shared between the groups and the calibration inside each of them
    workerPlan <- .consensusWorkers(nThreads = nThreads,
                                    nGroups = length(multiGroups),
                                    calibrate = isTRUE(consensusArguments$calibrate),
                                    userBPPARAM = consensusArguments$BPPARAM)

    groupArgumentList <-
      lapply(multiGroups,
             function(groupName) {
               groupSamples <- sampleTable$sample[sampleTable$group == groupName]

               groupArguments <- consensusArguments
               if (is.null(groupArguments$BPPARAM)) {groupArguments$BPPARAM <- workerPlan$inner}
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

               return(c(list(peaks = peakList[groupSamples], sampleNames = groupSamples), groupArguments))
             })
    names(groupArgumentList) <- multiGroups

    if (isTRUE(verbose) & length(multiGroups) > 1 & workerPlan$outer > 1) {
      message("Building the consensus of ", length(multiGroups), " groups on ", workerPlan$outer, " workers",
              if (workerPlan$inner > 1) {paste(",", workerPlan$inner, "threads each for the calibration")} else {""}, ".")
    }

    # The groups do not depend on each other, so they run side by side
    groupObjects <- BiocParallel::bplapply(X = groupArgumentList,
                                           FUN = .runGroupConsensus,
                                           BPPARAM = .makeParallelParam(nThreads = workerPlan$outer,
                                                                        tasks = length(multiGroups)))

    for (groupName in multiGroups) {
      consensusRanges <- consensusRegions::consensusRanges(groupObjects[[groupName]])

      # The statistics of one group mean nothing for the others, only the coordinates are pooled
      S4Vectors::mcols(consensusRanges) <- NULL
      groupConsensus[[groupName]] <- consensusRanges
      consensusObjects[[groupName]] <- groupObjects[[groupName]]
    }

    # Back in the order of the groups, whatever order they were computed in
    groupConsensus <- groupConsensus[groupLevels]
    consensusObjects <- consensusObjects[groupLevels]

    if (isTRUE(verbose)) {
      for (groupName in groupLevels) {
        if (groupSizes[[groupName]] == 1) {
          message("Group ", groupName, ": a single sample, its ", length(groupConsensus[[groupName]]), " peaks are taken as they are.")
        } else {
          message("Group ", groupName, ": ", groupSizes[[groupName]], " samples, ", length(groupConsensus[[groupName]]), " consensus regions.")
        }
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

      # With seqlevelsStyle = NULL the regions of the user may be written in another style than the peaks.
      # Regions that cannot be reconciled are left as they are, and the occupancy step says so
      userList <- lapply(userList,
                         function(userRanges) {
                           tryCatch(expr = .matchSeqlevels(x = userRanges,
                                                           targetSeqlevels = GenomeInfoDb::seqlevels(totalConsensus),
                                                           verbose = FALSE),
                                    error = function(e) {return(userRanges)})
                         })

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
                                removed = removedPeaks,
                                samples = sampleTable,
                                sheet = sampleSheet,
                                groupBy = groupBy,
                                mode = analysisMode,
                                blacklist = listRanges$blacklist,
                                greylist = listRanges$greylist)

    # Recorded as applyBlacklist and applyGreylist record them, so the lists follow the object downstream
    if (!is.null(listRanges$blacklist)) {
      regionSet@blacklist <- listRanges$blacklist
    }

    if (!is.null(listRanges$greylist)) {
      regionSet@parameters$greylist <- list(n.regions = length(listRanges$greylist),
                                            covered.bp = sum(as.numeric(BiocGenerics::width(listRanges$greylist))),
                                            applied.to = "peaks")
    }

    if (length(filteringSteps) > 0) {
      regionSet@filtering.log <- rbind(regionSet@filtering.log, do.call(what = rbind, args = unname(filteringSteps)))
    }

    # Only plain values go to the parameters, a BiocParallel object has no place in an exported record
    regionSet@parameters$loadConsensusPeaks <- list(groupBy = groupBy,
                                                    regionMode = analysisMode,
                                                    unassignedSet = unassignedSet,
                                                    seqlevelsStyle = seqlevelsStyle,
                                                    n.samples = nrow(sampleTable),
                                                    n.groups = length(groupLevels),
                                                    n.blacklist.regions = if (is.null(listRanges$blacklist)) {0L} else {length(listRanges$blacklist)},
                                                    n.greylist.regions = if (is.null(listRanges$greylist)) {0L} else {length(listRanges$greylist)},
                                                    consensusArguments = Filter(is.atomic, consensusArguments))

    methods::validObject(regionSet)
    return(regionSet)
  } # END function




#' @title .consensusWorkers
#'
#' @description Shares the threads between the groups built in parallel and the calibration run inside each group.
#'
#' @param nThreads Number of threads asked for.
#' @param nGroups Number of groups needing a consensus, the groups of a single sample left out.
#' @param calibrate Logical value to indicate whether the threshold is calibrated inside each group.
#' @param userBPPARAM The \code{BPPARAM} passed by the user to \code{consensusRegions::runConsensus}, or \code{NULL}.
#'
#' @return A list with \code{outer}, the number of groups run at once, and \code{inner}, the number of threads each group gets for its permutations.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.consensusWorkers <-
  function(nThreads,
           nGroups,
           calibrate = FALSE,
           userBPPARAM = NULL) {

    nThreads <- as.integer(nThreads[1])
    if (is.na(nThreads) | nThreads < 1) {
      stop("The 'nThreads' parameter must be a positive integer.", call. = FALSE)
    }

    # A back end chosen by the user goes to every group untouched, and the groups wait for each other
    if (!is.null(userBPPARAM)) {
      return(list(outer = 1L, inner = userBPPARAM))
    }

    outerWorkers <- as.integer(max(1L, min(nThreads, nGroups)))

    # Without calibration a consensus runs on one thread, so the spare threads would only sit idle
    innerWorkers <- if (isTRUE(calibrate)) {as.integer(max(1L, nThreads %/% outerWorkers))} else {1L}

    return(list(outer = outerWorkers, inner = innerWorkers))
  } # END function




#' @title .runGroupConsensus
#'
#' @description Builds the consensus of one group through \code{consensusRegions::runConsensus}, called by the workers of \code{\link{loadConsensusPeaks}}.
#'
#' @param groupArguments List with the arguments of \code{consensusRegions::runConsensus} for the group, the peaks included.
#'
#' @return A \code{consensusRegions} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.runGroupConsensus <-
  function(groupArguments) {

    return(do.call(what = consensusRegions::runConsensus, args = groupArguments))
  } # END function




#' @title .loadExclusionRegions
#'
#' @description Loads one or several exclusion lists and pools them into a single set of ranges, in the chromosome style of the peaks.
#'
#' @param excludeRegions A \code{GRanges}, a path, a data.frame, or a list of them.
#' @param seqlevelsStyle String with the chromosome naming style, or \code{NULL}.
#' @param listLabel String with the name of the parameter the list came from, used in the messages. Default: \code{"blacklist"}.
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
           seqlevelsStyle = "UCSC",
           listLabel = "blacklist") {

    # A single list of regions is wrapped, a list of lists is taken as it is
    if (methods::is(excludeRegions, "GRanges") | is.data.frame(excludeRegions) | is.character(excludeRegions)) {
      excludeRegions <- list(excludeRegions)
    }

    if (!is.list(excludeRegions)) {
      stop("The '", listLabel, "' parameter must be a GRanges, a path, a data.frame, or a list of them.", call. = FALSE)
    }

    names(excludeRegions) <- paste0(listLabel, "_", seq_along(excludeRegions))

    exclusionList <- loadRegions(regions = excludeRegions,
                                 keepMetadata = FALSE,
                                 duplicatedSets = "remove",
                                 seqlevelsStyle = seqlevelsStyle,
                                 outputFormat = "list",
                                 verbose = FALSE)

    return(IRanges::reduce(unlist(GenomicRanges::GRangesList(unname(exclusionList)), use.names = FALSE), ignore.strand = TRUE))
  } # END function




#' @title .checkListAssembly
#'
#' @description Refuses a list of regions built for another assembly than the one declared for the peaks, which would overlap them on the chromosome names alone.
#'
#' @param listInput The list as given: a \code{GRanges}, a path, a data.frame, or a list of them.
#' @param genomeAssembly String with the assembly of the peaks, or \code{NULL}.
#' @param listLabel String with the name of the parameter the list came from.
#'
#' @return Nothing, it stops when the two assemblies differ.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb genome
#' @importFrom methods is
#'
#' @keywords internal

.checkListAssembly <-
  function(listInput,
           genomeAssembly,
           listLabel) {

    if (is.null(genomeAssembly) || is.na(genomeAssembly[1]) || genomeAssembly[1] == "") {
      return(invisible(TRUE))
    }

    # Only a GRanges carries its assembly, a file or a table cannot be checked
    listElements <- if (methods::is(listInput, "GRanges")) {list(listInput)} else if (is.list(listInput) & !is.data.frame(listInput)) {listInput} else {list()}

    for (listElement in listElements) {
      if (!methods::is(listElement, "GRanges")) {next}

      listAssembly <- unique(as.character(GenomeInfoDb::genome(listElement)))
      listAssembly <- listAssembly[!is.na(listAssembly) & listAssembly != ""]

      # GRCh38 and hg38 name the same assembly, the aliases are resolved before the comparison
      if (length(listAssembly) == 1 && .resolveGenomeName(listAssembly) != .resolveGenomeName(as.character(genomeAssembly[1]))) {
        stop("The ", listLabel, " was built for ", listAssembly, " and the peaks are declared as ", genomeAssembly[1],
             ". Overlapping them would match the chromosome names and nothing else. Clear the assembly with genome() on the list to force it.", call. = FALSE)
      }
    }

    return(invisible(TRUE))
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
#' @return A list with \code{groups}, the consensus regions of every group; \code{objects}, the \code{consensusRegions} object behind each of them (\code{NULL} for the groups with a single sample); \code{total}, the pooled consensus; \code{peaks}, the peaks of every sample after the blacklist and the greylist; \code{removed}, the peaks they took out, with the \code{sample} they came from and the list that removed them in \code{removed.by}; \code{samples}, a data.frame with the group, the peak file, the number of peaks read (\code{n.peaks}), removed by the blacklist (\code{n.blacklist}), by the greylist (\code{n.greylist}) and by both (\code{n.excluded}) for every sample; \code{sheet}, the sample sheet the consensus was built from; \code{groupBy}; \code{mode}, one among \code{"consensus"}, \code{"split"} and \code{"replace"}; and \code{blacklist} and \code{greylist}, the regions the peaks were cleaned against, \code{NULL} when a list was not given.
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
#' @seealso \code{\link{loadConsensusPeaks}}, \code{\link{consensusGroupList}}, \code{\link{plotPeakUpset}}
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




#' @title consensusGroupList
#'
#' @description Returns the consensus of every group of samples built by \code{\link{loadConsensusPeaks}}, as a named list with one \code{GRanges} per group, ready to be annotated or compared group by group with tools taking a list of peak sets, such as \code{ChIPseeker}.
#'
#' @param object \code{RegionSetDE} object returned by \code{\link{loadConsensusPeaks}}.
#' @param groups Character vector with the groups returned, in the order wanted. Default: \code{NULL}, every group, in the order of the consensus.
#' @param seqlevelsStyle String with the chromosome naming style of the output, one among \code{"UCSC"}, \code{"Ensembl"} and \code{"NCBI"}, for instance \code{"UCSC"} to match a \code{TxDb} of the UCSC annotation. Default: \code{NULL}, the style of the object.
#' @param asGRangesList Logical value to indicate whether a \code{GRangesList} must be returned instead of a list. Default: \code{FALSE}.
#'
#' @return A named list of \code{GRanges}, or a \code{GRangesList}, with one element per group holding its consensus regions, named after the group.
#'
#' @details The consensus of a group is the one built within that group, before the groups are pooled, so a region found in two groups is in both elements and a region found in one group only is in that one alone. The regions carry no metadata column, since the statistics of the consensus of one group mean nothing for the others; the consensus object of every group, with its statistics, is in \code{consensusData(object)$objects}. A group with a single sample holds the peaks of that sample, merged where they overlap, and the peaks removed by the \code{blacklist} and the \code{greylist} are absent from every group.
#'
#' The list goes as it is to the functions of \code{ChIPseeker} that take several peak sets, for instance \code{lapply(consensusGroupList(x, seqlevelsStyle = "UCSC"), ChIPseeker::annotatePeak, TxDb = txdb)} followed by \code{ChIPseeker::plotAnnoBar()}. The naming style of the chromosomes has to match the one of the annotation, which is what \code{seqlevelsStyle} is for.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'   sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#'   consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
#'
#'   groupList <- consensusGroupList(consensus, seqlevelsStyle = "UCSC")
#'   lengths(groupList)
#'   groupList$R1881_24h
#'
#'   # Two groups only, in the order given
#'   names(consensusGroupList(consensus, groups = c("R1881_24h", "DMSO")))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{consensusData}}, \code{\link{loadConsensusPeaks}}, \code{\link{plotPeakUpset}}
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom methods is
#'
#' @export consensusGroupList

consensusGroupList <-
  function(object,
           groups = NULL,
           seqlevelsStyle = NULL,
           asGRangesList = FALSE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(object, "RegionSetDE") || length(object@consensus) == 0) {
      stop("The 'object' parameter must be the RegionSetDE object returned by loadConsensusPeaks().", call. = FALSE)
    }

    if (!is.logical(asGRangesList) | length(asGRangesList) != 1 | anyNA(asGRangesList)) {
      stop("The 'asGRangesList' parameter must be TRUE or FALSE.", call. = FALSE)
    }

    groupList <- object@consensus$groups

    if (!is.null(groups)) {
      absentGroups <- setdiff(groups, names(groupList))
      if (length(absentGroups) > 0) {
        stop("The following groups are absent from the consensus: ", paste(absentGroups, collapse = ", "),
             ". The groups are: ", paste(names(groupList), collapse = ", "), ".", call. = FALSE)
      }
      groupList <- groupList[groups]
    }

    #------------------------#
    # Naming style           #
    #------------------------#
    # An annotation of another style would find no gene next to any region, so the names are converted on request
    if (!is.null(seqlevelsStyle)) {
      groupList <- lapply(groupList, .styleSeqlevels, seqlevelsStyle = seqlevelsStyle)
    }

    if (isTRUE(asGRangesList)) {
      return(GenomicRanges::GRangesList(groupList))
    }

    return(groupList)
  } # END function
