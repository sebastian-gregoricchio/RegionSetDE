# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title loadConsensusPeaks
#'
#' @description Builds the regions of an analysis from the peaks called on each sample. The peaks are combined into a consensus within each group of samples through \code{consensusRegions}, and the group consensus are pooled into a total one, which becomes the region set counted downstream. Regions supplied by the user can split the total consensus into sets, or take its place.
#'
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, with at least the \code{sample} and \code{peaks} columns.
#' @param groupBy String with the column of the sample sheet defining the groups, for instance \code{"condition"}. Default: \code{NULL}, all the samples form a single group.
#' @param technicalBy String with the column of the sample sheet naming the biological sample each library comes from, for instance \code{"bioRep"}. The libraries of a group sharing a value are technical replicates: their peaks are joined first into one peak set per biological sample, which then enters the consensus of the group as a single replicate. A label found in more than one group, such as \code{R1} in every condition, names a different sample in each of them. Default: \code{NULL}, every library is a biological replicate.
#' @param technicalArguments Named list of arguments of \code{consensusRegions::runConsensus}, or of \code{consensusRegions::buildConsensus}, for the step joining the technical replicates, for instance \code{list(replicateType = "biological", minReplicates = 2)} to keep the peaks found in two libraries out of three. They take precedence over the ones this step inherits from the consensus of the groups. The calibration, the recentring and their arguments are reserved to the consensus of the groups. Ignored without \code{technicalBy}. Default: \code{NULL}, the arguments of the consensus of the groups, with \code{replicateType = "technical"} and equal weights.
#' @param blacklist Regions of the assembly whose peaks must be dropped before any consensus is built, typically the list returned by \code{\link{loadBlacklist}}. Either a \code{GRanges}, a path to a BED-like file, a data.frame, or a list of them, which are pooled. The list is stored in the \code{blacklist} slot (read with \code{\link{blacklistRanges}}) of the object, as \code{\link{applyBlacklist}} does. Default: \code{NULL}.
#' @param greylist Regions of this experiment whose peaks must be dropped before any consensus is built, typically the list built from the inputs by \code{\link{makeGreylist}}. Accepts the same forms as \code{blacklist}. Default: \code{NULL}.
#' @param regionSets Regions of interest of the user, in any form accepted by \code{\link{loadRegions}}: a named list of paths, \code{GRanges} or data.frames. Default: \code{NULL}, the total consensus is the only set.
#' @param regionMode String indicating what \code{regionSets} do, either \code{"split"}, where every region of the total consensus goes to the first set it overlaps, or \code{"replace"}, where the sets of the user are the regions and the peaks only annotate them. Ignored without \code{regionSets}. Default: \code{"split"}.
#' @param unassignedSet String with the name of the set collecting the consensus regions that overlap none of \code{regionSets} in \code{"split"} mode. \code{NULL} drops them. Default: \code{"other"}.
#' @param seqlevelsStyle String indicating the chromosome naming style, one among \code{"UCSC"}, \code{"Ensembl"} or \code{"NCBI"}, or \code{NULL} to keep the names as they are. Default: \code{"UCSC"}.
#' @param genomeAssembly String indicating the genome assembly to store with the regions, e.g. \code{"hg38"}. Default: \code{NULL}.
#' @param BPPARAM A \code{BiocParallelParam} object, used as it is, or a number of cores as a shorthand for one. A number of cores is shared out: the groups with more than one sample are built side by side, one worker per group, and with \code{calibrate = TRUE} the cores left over go to the calibration inside each group. A \code{BiocParallelParam} object goes to the calibration of every group when \code{calibrate = TRUE}, the groups being built one after the other, and runs the groups side by side otherwise, so that two levels of workers are never stacked. Default: \code{NULL}, a \code{BiocParallel::SerialParam()}, so that nothing runs in parallel unasked.
#' @param seed Number seeding the random numbers of the calibration, when \code{calibrate = TRUE}. It goes to \code{runConsensus} for every group, which hands it to the \code{BiocParallel} back end of the permutations as its \code{RNGseed}, and it is stored in \code{parameters$loadConsensusPeaks$seed}, so that the consensus can be built again identically, on any number of cores. Needs consensusRegions 0.99.2 or later. Default: \code{NULL}, no seed, and a calibrated consensus may change from one call to the next.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#' @param weightMethod String with the weights of the replicates within a group, passed to \code{consensusRegions::computeReplicateWeights}: \code{"equal"}, \code{"frip"}, \code{"librarySize"} or \code{"intrinsic"}. \code{"equal"} counts every replicate the same and is the reproducible choice. With \code{"frip"} or \code{"librarySize"} the BAM files of the sample sheet are used. Default: \code{"equal"}.
#' @param calibrate Logical value to indicate whether the threshold on the combined p-value must be estimated within each group by \code{consensusRegions::calibrateThreshold}, against shuffled peaks, instead of being fixed by \code{combinedThreshold}. Default: \code{FALSE}.
#' @param nPermutations Number of shuffles of the calibration. Ignored when \code{calibrate = FALSE}. Default: \code{50}.
#' @param targetFDR Numeric value with the empirical false discovery rate the calibration aims for. Ignored when \code{calibrate = FALSE}. Default: \code{0.05}.
#' @param recentre Logical value to indicate whether the consensus of each group must be reduced to windows of fixed width around the summits, as \code{consensusRegions::recentrePeaks} does. It suits transcription factors and ATAC-seq, not broad domains, and with several groups the union merges overlapping windows back into wider regions. Default: \code{FALSE}.
#' @param width Width in base pairs of the recentred windows. Ignored when \code{recentre = FALSE}. Default: \code{400}.
#' @param combinationMethod String with the way the p-values of overlapping peaks are combined, one among \code{"stouffer"}, \code{"fisher"}, \code{"lancaster"} and \code{"rankProduct"}. Default: \code{"stouffer"}.
#' @param replicateType String with the kind of replicates in a group: \code{"biological"} asks for \code{minSupport} supporting replicates, \code{"technical"} for all of them. Default: \code{"biological"}.
#' @param stringencyThreshold Numeric value with the p-value below which a peak is stringent. Default: \code{1e-8}.
#' @param weakThreshold Numeric value with the p-value above which a peak is taken as background and plays no further part. Default: \code{1e-4}.
#' @param combinedThreshold Numeric value with the threshold on the combined p-value. Its scale depends on \code{combinationMethod}, so a value chosen for one scheme does not mean the same under another. Default: \code{NULL}, \code{stringencyThreshold} is used.
#' @param alpha Numeric value with the level of the Benjamini-Hochberg step within each replicate. Default: \code{0.05}.
#' @param minSupport Number of other replicates of the group that must hold an overlapping peak. Give this or \code{minReplicates}, not both. Default: \code{NULL}, one supporting replicate is enough.
#' @param minReplicates Number of replicates of the group that must hold the peak, counting the one it came from, as a count (\code{2}), a proportion (\code{0.66}) or a percentage (\code{"66\%"}), the last two rounded up. \code{minReplicates = 2} is the same requirement as \code{minSupport = 1}. Default: \code{NULL}.
#' @param minOverlap Minimum overlap in base pairs between two peaks for one to support the other. Default: \code{1}.
#' @param minOverlapFraction Minimum overlap as a fraction of the shorter of the two peaks, required on top of \code{minOverlap}. Default: \code{NULL}, no such requirement.
#' @param mergeMethod String with the way the confirmed peaks are merged: \code{"reduce"} merges anything that touches, \code{"iterative"} keeps the most significant peak, removes what overlaps it and repeats. Default: \code{"reduce"}.
#' @param maxConsensusWidth Maximum width in base pairs of a merged region. Wider regions are rebuilt with the iterative rule. Default: \code{NULL}, no limit.
#' @param ... Any other argument of \code{consensusRegions::runConsensus} or \code{consensusRegions::buildConsensus}, for instance \code{adjustmentFamily} or \code{minSupportWeight}. \code{peaks} and \code{sampleNames} come from the sample sheet and cannot be given.
#'
#' @return A \code{RegionSetDE} object. Every region carries \code{peak.<group>}, telling whether the consensus of that group overlaps it, \code{peak.groups}, the number of groups with a consensus peak on it, and \code{peak.samples}, the number of samples with a peak on it. With \code{technicalBy}, \code{peak.samples} counts the libraries and \code{peak.biological} the biological samples. The \code{consensus} slot (read with \code{\link{consensusData}}), read with \code{\link{consensusData}}, keeps the consensus of every group, the \code{consensusRegions} object behind it, the peaks of every sample after the blacklist and the greylist, the peaks they removed, the two lists, the total consensus, the table of the samples and the sample sheet itself, which \code{\link{countReads}} and \code{\link{countBigwig}} read when no file is given to them. The blacklist goes to the \code{blacklist} slot (read with \code{\link{blacklistRanges}}), the greylist to the \code{greylist} slot (read with \code{\link{greylistRanges}}) with its size in \code{parameters$greylist}, and the peaks each of them removed from every sample to the \code{filtering.log}, read with \code{\link{filteringLog}}. Where the random numbers of the calibration came from is recorded in \code{parameters$loadConsensusPeaks}, see Details.
#'
#' @details The consensus is built group by group and then pooled, rather than once over all the samples. A single requirement over all the libraries, such as a peak in at least two of them, favours the larger group, while the same requirement applied within each group treats them alike, and the union keeps the regions found in one group only. A group with a single sample has no replicate to agree with, and its peaks stand for the group as they are.
#'
#' Technical replicates, given through \code{technicalBy}, are joined before the consensus of the group. The combination of \code{consensusRegions} takes every replicate as an independent piece of evidence, and two libraries of the same biological sample are not: put side by side with the other samples of the group, they would count the same sample twice, and a peak found in the two libraries of one sample alone would meet a requirement of two supporting replicates. Here the libraries of each biological sample are joined first, by default keeping the peaks found in all of them, and the peak set of every biological sample enters the consensus of its group as one replicate. That peak set carries the combined p-value of its libraries, which with the Stouffer method is close to what a peak caller would report on the pooled reads, and the summit of its most significant peak, so that \code{recentre} still has a summit to work from. The calibration, when asked for, runs on the consensus of the groups only. A biological sample with a single library keeps its peaks as they are. The consensus of every biological sample is kept in \code{consensusData(x)$technical}. Weights taken from the BAM files of the sheet, with \code{weightMethod = "frip"} or \code{"librarySize"}, describe libraries rather than biological samples, and they are refused with \code{technicalBy} unless the values are given through \code{...}.
#'
#' The peaks are read by \code{consensusRegions}, which keeps the standard chromosomes only, as \code{GenomeInfoDb::keepStandardChromosomes} defines them: scaffolds, patches and unplaced contigs are dropped, and so are chromosomes whose names it does not recognise.
#'
#' The peaks overlapping \code{blacklist} and \code{greylist} are removed before the consensus, not after it. An artefact lying next to a genuine peak would otherwise merge with it, and removing the merged region afterwards would take the genuine peak away as well. The blacklist is applied first and the greylist to the peaks left, so a peak lying on both is counted as blacklisted. The two lists are kept apart because they say different things: a blacklist describes the assembly and is the same for every experiment, a greylist describes the inputs of this one. Both are recorded the way \code{\link{applyBlacklist}} and \code{\link{applyGreylist}} record them, so the counts, the fit and the results built from the object carry the two lists with them, and \code{\link{countReads}} leaves their reads out of the library sizes. \code{consensusData(x)$removed} holds the peaks taken out, with the sample they came from and the list that removed them, and \code{consensusData(x)$samples} their number per sample. A blacklist built for another assembly than \code{genomeAssembly} is refused, as \code{\link{applyBlacklist}} refuses it.
#'
#' With \code{regionMode = "split"} the consensus regions are assigned to the sets of \code{regionSets} in the order they are given, a region overlapping two sets going to the first one. The sets then describe classes of peaks, such as promoter and distal ones, and \code{\link{testRegionSets}} can compare them. With \code{regionMode = "replace"} the counted regions are those of \code{regionSets}, and the peaks only tell which of them are occupied in each group.
#'
#' The occupancy columns describe where the peaks were called, and they are a poor basis for a set of regions to test. A set of the regions found in one group only was selected on the signal of that group, so a test of its change between the groups answers a question already settled by the selection.
#'
#' The consensus of a group does not depend on the other groups, so with more than one core the groups run in parallel. \code{consensusRegions} builds a consensus on a single core and only parallelises the permutations of the calibration, so without \code{calibrate} the groups are the one place where more cores save time: three groups on three cores take about as long as the largest of them. With \code{calibrate = TRUE} a number of cores is shared, divided by the number of groups for the permutations of each group. A \code{BiocParallelParam} object is never split: with \code{calibrate = TRUE} it is handed to every group as it is and the groups run one after the other, without calibration it runs the groups.
#'
#' The calibration is the only step drawing random numbers: it places the peaks at random to see how often they overlap by chance. The draws come from \code{BiocParallel}, and without a seed of its own a \code{BiocParallel} back end takes them from a stream it keeps for the whole session, which a \code{set.seed()} before the call does not reset: the calibrated threshold, and with it the consensus, can then change from one call to the next. \code{seed} is passed on to \code{consensusRegions::runConsensus} in every group, which gives it to the back end of the permutations as its \code{RNGseed}: the same seed returns the same consensus whatever the number of cores, and the random numbers of the session are left as they were. Neither package calls \code{set.seed()}. \code{parameters$loadConsensusPeaks} keeps the \code{seed} and its \code{seed.source}: \code{"seed"}, \code{"BPPARAM"} when a \code{BiocParallelParam} object given as \code{BPPARAM} carries its own \code{RNGseed}, which is then the one recorded, \code{"unseeded"} when the calibration ran without a seed, or \code{"none"} without calibration. A \code{BiocParallelParam} object is used as it is, so \code{seed} is ignored with it, with a warning.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'     peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
#'         pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE
#'     )
#'
#'     # Three replicates of one condition and two of another, peaks only
#'     sampleSheet <- data.frame(
#'         sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
#'         bam = c("A_1.bam", "A_2.bam", "A_3.bam", "B_1.bam", "B_2.bam"),
#'         peaks = peakFiles[c(1, 2, 3, 1, 2)],
#'         condition = c("A", "A", "A", "B", "B")
#'     )
#'
#'     sampleSheet <- loadSampleSheet(sampleSheet, checkFiles = FALSE, verbose = FALSE)
#'
#'     regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition")
#'     regions
#'
#'     head(regionRanges(regions)$consensus, 3)
#'
#'     # The same peaks split by regions of interest of the user
#'     firstHalf <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 5e5))
#'     splitRegions <- loadConsensusPeaks(sampleSheet,
#'         groupBy = "condition",
#'         regionSets = list(firstHalf = firstHalf),
#'         verbose = FALSE
#'     )
#'     lengths(regionRanges(splitRegions))
#'
#'     # Two libraries per biological sample, joined before the consensus of each condition
#'     libraryNames <- c("A_1a", "A_1b", "A_2a", "A_2b", "B_1a", "B_1b", "B_2a", "B_2b")
#'     technicalSheet <- loadSampleSheet(
#'         data.frame(
#'             sample = libraryNames,
#'             bam = paste0(libraryNames, ".bam"),
#'             peaks = peakFiles[c(1, 2, 2, 3, 1, 3, 2, 3)],
#'             condition = rep(c("A", "B"), each = 4),
#'             bioRep = rep(c("R1", "R1", "R2", "R2"), times = 2)
#'         ),
#'         checkFiles = FALSE, verbose = FALSE
#'     )
#'
#'     technicalRegions <- loadConsensusPeaks(technicalSheet,
#'         groupBy = "condition",
#'         technicalBy = "bioRep",
#'         verbose = FALSE
#'     )
#'     names(consensusData(technicalRegions)$technical)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotPeakUpset}}, \code{\link{consensusData}}, \code{\link{loadSampleSheet}}, \code{\link{loadBlacklist}}, \code{\link{makeGreylist}}
#'
#' @importFrom dplyr filter distinct count mutate group_by summarise first arrange %>%
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges reduce overlapsAny IRanges
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocParallel bplapply bpRNGseed
#' @importFrom methods is validObject
#' @importFrom utils packageVersion
#'
#' @export loadConsensusPeaks

loadConsensusPeaks <-
    function(sampleSheet,
             groupBy = NULL,
             technicalBy = NULL,
             technicalArguments = NULL,
             blacklist = NULL,
             greylist = NULL,
             regionSets = NULL,
             regionMode = "split",
             unassignedSet = "other",
             seqlevelsStyle = "UCSC",
             genomeAssembly = NULL,
             BPPARAM = NULL,
             seed = NULL,
             weightMethod = "equal",
             calibrate = FALSE,
             nPermutations = 50L,
             targetFDR = 0.05,
             recentre = FALSE,
             width = 400L,
             combinationMethod = "stouffer",
             replicateType = "biological",
             stringencyThreshold = 1e-8,
             weakThreshold = 1e-4,
             combinedThreshold = NULL,
             alpha = 0.05,
             minSupport = NULL,
             minReplicates = NULL,
             minOverlap = 1L,
             minOverlapFraction = NULL,
             mergeMethod = "reduce",
             maxConsensusWidth = NULL,
             verbose = TRUE,
             ...) {
        #------------------------#
        # Check of the arguments #
        #------------------------#
        checkedArguments <- .checkLoadConsensusPeaksArguments(
            sampleSheet = sampleSheet,
            groupBy = groupBy,
            technicalBy = technicalBy,
            technicalArguments = technicalArguments,
            regionMode = regionMode,
            unassignedSet = unassignedSet,
            seed = seed,
            weightMethod = weightMethod,
            calibrate = calibrate,
            nPermutations = nPermutations,
            targetFDR = targetFDR,
            recentre = recentre,
            width = width,
            combinationMethod = combinationMethod,
            replicateType = replicateType,
            stringencyThreshold = stringencyThreshold,
            weakThreshold = weakThreshold,
            combinedThreshold = combinedThreshold,
            alpha = alpha,
            minSupport = minSupport,
            minReplicates = minReplicates,
            minOverlap = minOverlap,
            minOverlapFraction = minOverlapFraction,
            mergeMethod = mergeMethod,
            maxConsensusWidth = maxConsensusWidth,
            verbose = verbose,
            ...
        )
        consensusArguments <- checkedArguments$consensusArguments
        technicalArguments <- checkedArguments$technicalArguments
        regionMode <- checkedArguments$regionMode
        sampleSheet <- checkedArguments$sampleSheet

        #-------------------------------#
        # Samples and groups            #
        #-------------------------------#
        sampleTable <- .consensusSampleTable(
            sampleSheet = sampleSheet,
            groupBy = groupBy,
            technicalBy = technicalBy,
            verbose = verbose
        )

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
        peakList <- consensusRegions::readPeakSets(
            peaks = sampleTable$peak.file,
            sampleNames = sampleTable$sample,
            seqlevelsStyle = seqlevelsStyle,
            verbose = FALSE
        )

        sampleTable$n.peaks <- as.integer(lengths(peakList))

        #-------------------------------#
        # Blacklist, then greylist      #
        #-------------------------------#
        # An artefact next to a genuine peak would merge with it, so both lists act before the consensus.
        # The greylist sees only the peaks the blacklist left, so a peak on both is counted once, as blacklisted.
        exclusionStep <- .excludeListedPeaks(
            blacklist = blacklist,
            greylist = greylist,
            seqlevelsStyle = seqlevelsStyle,
            genomeAssembly = genomeAssembly,
            verbose = verbose,
            peakList = peakList,
            sampleTable = sampleTable
        )
        filteringSteps <- exclusionStep$filteringSteps
        listRanges <- exclusionStep$listRanges
        peakList <- exclusionStep$peakList
        removedPeaks <- exclusionStep$removedPeaks
        sampleTable <- exclusionStep$sampleTable

        #-------------------------------#
        # Join the technical replicates #
        #-------------------------------#
        # Without technicalBy every library is its own biological sample and its peaks go on unchanged
        technicalStep <- .biologicalPeakSets(
            technicalBy = technicalBy,
            technicalArguments = technicalArguments,
            seqlevelsStyle = seqlevelsStyle,
            BPPARAM = BPPARAM,
            verbose = verbose,
            peakList = peakList,
            sampleTable = sampleTable
        )
        biologicalPeaks <- technicalStep$biologicalPeaks
        biologicalTable <- technicalStep$biologicalTable
        technicalObjects <- technicalStep$technicalObjects

        #-------------------------------#
        # Consensus of every group      #
        #-------------------------------#
        consensusStep <- .groupConsensusList(
            technicalBy = technicalBy,
            seqlevelsStyle = seqlevelsStyle,
            BPPARAM = BPPARAM,
            seed = seed,
            verbose = verbose,
            biologicalPeaks = biologicalPeaks,
            biologicalTable = biologicalTable,
            consensusArguments = consensusArguments,
            groupLevels = groupLevels,
            sampleTable = sampleTable
        )
        consensusObjects <- consensusStep$consensusObjects
        groupConsensus <- consensusStep$groupConsensus
        recordedSeed <- consensusStep$recordedSeed
        seedSource <- consensusStep$seedSource

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
        regionStep <- .consensusRegionList(
            technicalBy = technicalBy,
            regionSets = regionSets,
            regionMode = regionMode,
            unassignedSet = unassignedSet,
            seqlevelsStyle = seqlevelsStyle,
            genomeAssembly = genomeAssembly,
            verbose = verbose,
            biologicalPeaks = biologicalPeaks,
            groupConsensus = groupConsensus,
            peakList = peakList,
            totalConsensus = totalConsensus
        )
        analysisMode <- regionStep$analysisMode
        regionList <- regionStep$regionList

        regionSet <- loadRegions(
            regions = regionList,
            keepMetadata = TRUE,
            seqlevelsStyle = seqlevelsStyle,
            genomeAssembly = genomeAssembly,
            verbose = FALSE
        )

        #-------------------------------#
        # Keep the consensus data       #
        #-------------------------------#
        consensusData(regionSet) <- list(
            groups = groupConsensus,
            objects = consensusObjects,
            technical = technicalObjects,
            total = totalConsensus,
            peaks = peakList,
            removed = removedPeaks,
            samples = sampleTable,
            sheet = sampleSheet,
            groupBy = groupBy,
            technicalBy = technicalBy,
            mode = analysisMode,
            blacklist = listRanges$blacklist,
            greylist = listRanges$greylist
        )

        # Recorded as applyBlacklist and applyGreylist record them, so the lists follow the object downstream
        if (!is.null(listRanges$blacklist)) {
            blacklistRanges(regionSet) <- listRanges$blacklist
        }

        if (!is.null(listRanges$greylist)) {
            greylistRanges(regionSet) <- listRanges$greylist
            parameterLog(regionSet)$greylist <- list(
                n.regions = length(listRanges$greylist),
                covered.bp = sum(as.numeric(BiocGenerics::width(listRanges$greylist))),
                applied.to = "peaks"
            )
        }

        if (length(filteringSteps) > 0) {
            filteringLog(regionSet) <- rbind(filteringLog(regionSet), do.call(what = rbind, args = unname(filteringSteps)))
        }

        # Only plain values go to the parameters, a BiocParallel object has no place in an exported record
        parameterLog(regionSet)$loadConsensusPeaks <- list(
            groupBy = groupBy,
            technicalBy = technicalBy,
            regionMode = analysisMode,
            unassignedSet = unassignedSet,
            seqlevelsStyle = seqlevelsStyle,
            n.samples = nrow(sampleTable),
            n.biological = nrow(biologicalTable),
            n.groups = length(groupLevels),
            n.blacklist.regions = if (is.null(listRanges$blacklist)) {
                0L
            } else {
                length(listRanges$blacklist)
            },
            n.greylist.regions = if (is.null(listRanges$greylist)) {
                0L
            } else {
                length(listRanges$greylist)
            },
            seed = recordedSeed,
            seed.source = seedSource,
            consensusArguments = Filter(function(x) {
                !is.null(x) && is.atomic(x)
            }, consensusArguments),
            technicalArguments = if (is.null(technicalBy)) {
                NULL
            } else {
                Filter(function(x) {
                    !is.null(x) && is.atomic(x)
                }, technicalArguments)
            }
        )

        methods::validObject(regionSet)
        return(regionSet)
    } # END function


#' @title .consensusSampleTable
#'
#' @description Builds the table of the libraries of \code{\link{loadConsensusPeaks}}, one row per library with its peak file, group and biological sample, and checks it.
#'
#' @param sampleSheet Sample sheet as checked by \code{.checkLoadConsensusPeaksArguments}.
#' @param groupBy See \code{\link{loadConsensusPeaks}}.
#' @param technicalBy See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#'
#' @return A data.frame with one row per library holding a peak file.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr count distinct filter mutate
#'
#' @noRd

.consensusSampleTable <-
    function(sampleSheet,
             groupBy,
             technicalBy,
             verbose) {
        groupValues <- if (is.null(groupBy)) {
            rep("all", nrow(sampleSheet))
        } else {
            as.character(sampleSheet[[groupBy]])
        }

        sampleTable <- data.frame(
            sample = as.character(sampleSheet$sample),
            group = groupValues,
            biological = if (is.null(technicalBy)) {
                as.character(sampleSheet$sample)
            } else {
                as.character(sampleSheet[[technicalBy]])
            },
            peak.file = as.character(sampleSheet$peaks),
            bam = if ("bam" %in% colnames(sampleSheet)) {
                as.character(sampleSheet$bam)
            } else {
                NA_character_
            },
            stringsAsFactors = FALSE
        )

        if (anyNA(sampleTable$group)) {
            stop("The following samples have no value in the '", groupBy, "' column: ",
                paste(sampleTable$sample[is.na(sampleTable$group)], collapse = ", "), ".",
                call. = FALSE
            )
        }

        if (anyNA(sampleTable$biological)) {
            stop("The following samples have no value in the '", technicalBy, "' column: ",
                paste(sampleTable$sample[is.na(sampleTable$biological)], collapse = ", "), ".",
                call. = FALSE
            )
        }

        withoutPeaks <- dplyr::filter(sampleTable, is.na(.data$peak.file))
        if (isTRUE(verbose) & nrow(withoutPeaks) > 0) {
            message("The following samples have no peak file and stay out of the consensus: ", paste(withoutPeaks$sample, collapse = ", "), ".")
        }

        sampleTable <- dplyr::filter(sampleTable, !is.na(.data$peak.file))

        if (nrow(sampleTable) < 2) {
            stop("At least two samples with a peak file are needed to build a consensus.", call. = FALSE)
        }

        # A label repeated across the groups, such as R1 in every condition, names a different sample in each of them
        if (!is.null(technicalBy)) {
            sharedLabels <- sampleTable %>%
                dplyr::distinct(.data$group, .data$biological) %>%
                dplyr::count(.data$biological) %>%
                dplyr::filter(.data$n > 1)

            if (nrow(sharedLabels) > 0) {
                sampleTable <- dplyr::mutate(sampleTable, biological = paste(.data$group, .data$biological, sep = "_"))
                if (isTRUE(verbose)) {
                    message(
                        "The labels of '", technicalBy, "' repeat across the groups, the biological samples are named after both: ",
                        paste(unique(sampleTable$biological), collapse = ", "), "."
                    )
                }
            }
        }

        return(sampleTable)
    }


#' @title .excludeListedPeaks
#'
#' @description Removes from the peaks of every library those overlapping the blacklist, then those overlapping the greylist, and records what each list removed.
#'
#' @param blacklist See \code{\link{loadConsensusPeaks}}.
#' @param greylist See \code{\link{loadConsensusPeaks}}.
#' @param seqlevelsStyle See \code{\link{loadConsensusPeaks}}.
#' @param genomeAssembly See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#' @param peakList List with the peaks of every library.
#' @param sampleTable Data.frame with one row per library.
#'
#' @return A list with \code{peakList} and \code{sampleTable} updated, \code{listRanges}, the two lists as used, \code{removedPeaks}, the peaks removed with their sample and list, and \code{filteringSteps}, the entries of the filtering log.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors mcols
#'
#' @noRd

.excludeListedPeaks <-
    function(blacklist,
             greylist,
             seqlevelsStyle,
             genomeAssembly,
             verbose,
             peakList,
             sampleTable) {
        listRanges <- list(blacklist = NULL, greylist = NULL)
        removedList <- list()
        filteringSteps <- list()

        for (listLabel in c("blacklist", "greylist")) {
            listInput <- if (listLabel == "blacklist") {
                blacklist
            } else {
                greylist
            }
            sampleTable[[paste0("n.", listLabel)]] <- 0L

            if (is.null(listInput)) {
                next
            }

            .checkListAssembly(listInput = listInput, genomeAssembly = genomeAssembly, listLabel = listLabel)
            listRanges[[listLabel]] <- .loadExclusionRegions(excludeRegions = listInput, seqlevelsStyle = seqlevelsStyle, listLabel = listLabel)

            # With seqlevelsStyle = NULL the peaks keep the names of their files, and the list has to follow them
            listRanges[[listLabel]] <- tryCatch(
                expr = .matchSeqlevels(
                    x = listRanges[[listLabel]],
                    targetSeqlevels = GenomeInfoDb::seqlevels(peakList),
                    verbose = FALSE
                ),
                error = function(e) {
                    return(listRanges[[listLabel]])
                }
            )

            peaksBefore <- as.integer(lengths(peakList))
            overlapList <- lapply(as.list(peakList), function(peakRanges) {
                IRanges::overlapsAny(peakRanges, listRanges[[listLabel]], ignore.strand = TRUE)
            })

            # The peaks taken out are kept, with the sample they came from and the list that removed them
            removedList[[listLabel]] <- lapply(names(peakList), function(sampleName) {
                removedPeaks <- peakList[[sampleName]][overlapList[[sampleName]]]
                S4Vectors::mcols(removedPeaks)$sample <- rep(sampleName, length(removedPeaks))
                S4Vectors::mcols(removedPeaks)$removed.by <- rep(listLabel, length(removedPeaks))
                return(removedPeaks)
            })

            peakList <- GenomicRanges::GRangesList(lapply(names(peakList), function(sampleName) {
                peakList[[sampleName]][!overlapList[[sampleName]]]
            }))
            names(peakList) <- sampleTable$sample
            peaksAfter <- as.integer(lengths(peakList))

            sampleTable[[paste0("n.", listLabel)]] <- peaksBefore - peaksAfter

            filteringSteps[[listLabel]] <- data.frame(
                step = listLabel,
                region.set = paste(sampleTable$sample, "peaks"),
                n.before = peaksBefore,
                n.after = peaksAfter,
                n.removed = peaksBefore - peaksAfter,
                stringsAsFactors = FALSE
            )

            if (isTRUE(verbose)) {
                message(
                    "Peaks overlapping the ", listLabel, " (", length(listRanges[[listLabel]]), " regions) removed: ",
                    paste(sampleTable$sample, paste(peaksBefore - peaksAfter, peaksBefore, sep = "/"), collapse = ", "), "."
                )
            }
        }

        sampleTable$n.excluded <- sampleTable$n.blacklist + sampleTable$n.greylist

        removedPeaks <- if (length(removedList) == 0) {
            GenomicRanges::GRanges()
        } else {
            unlist(GenomicRanges::GRangesList(unlist(removedList, recursive = FALSE, use.names = FALSE)), use.names = FALSE)
        }

        return(list(
            filteringSteps = filteringSteps,
            listRanges = listRanges,
            peakList = peakList,
            removedPeaks = removedPeaks,
            sampleTable = sampleTable
        ))
    }


#' @title .biologicalPeakSets
#'
#' @description Joins the libraries of every biological sample of \code{\link{loadConsensusPeaks}} into one peak set, when \code{technicalBy} is given. Without it every library stands for its own sample.
#'
#' @param technicalBy See \code{\link{loadConsensusPeaks}}.
#' @param technicalArguments See \code{\link{loadConsensusPeaks}}.
#' @param seqlevelsStyle See \code{\link{loadConsensusPeaks}}.
#' @param BPPARAM See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#' @param peakList List with the peaks of every library, after the exclusion.
#' @param sampleTable Data.frame with one row per library.
#'
#' @return A list with \code{biologicalTable}, one row per biological sample, \code{biologicalPeaks}, its peak sets, and \code{technicalObjects}, the consensus objects of the joining step.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr arrange filter first group_by summarise
#'
#' @noRd

.biologicalPeakSets <-
    function(technicalBy,
             technicalArguments,
             seqlevelsStyle,
             BPPARAM,
             verbose,
             peakList,
             sampleTable) {
        biologicalTable <- sampleTable %>%
            dplyr::group_by(.data$biological) %>%
            dplyr::summarise(group = dplyr::first(.data$group), .groups = "drop") %>%
            dplyr::arrange(match(.data$biological, unique(sampleTable$biological))) %>%
            as.data.frame()

        technicalObjects <- list()
        biologicalPeaks <- peakList

        if (!is.null(technicalBy)) {
            technicalList <- .joinTechnicalReplicates(
                sampleTable = sampleTable,
                peakList = peakList,
                technicalArguments = technicalArguments,
                seqlevelsStyle = seqlevelsStyle,
                BPPARAM = BPPARAM,
                verbose = verbose
            )
            technicalObjects <- technicalList$objects
            biologicalPeaks <- technicalList$peaks

            # An empty peak set would leave the group without a statistic to combine, so the sample stays out
            emptySamples <- names(biologicalPeaks)[lengths(biologicalPeaks) == 0]
            if (length(emptySamples) > 0) {
                warning("The libraries of the following biological samples share no peak, and the samples stay out of the consensus of their group: ",
                    paste(emptySamples, collapse = ", "), ".",
                    call. = FALSE
                )
                biologicalTable <- dplyr::filter(biologicalTable, !(.data$biological %in% emptySamples))
                biologicalPeaks <- biologicalPeaks[biologicalTable$biological]
            }
        }

        return(list(
            biologicalPeaks = biologicalPeaks,
            biologicalTable = biologicalTable,
            technicalObjects = technicalObjects
        ))
    }


#' @title .groupConsensusList
#'
#' @description Builds the consensus of every group of \code{\link{loadConsensusPeaks}} through \code{consensusRegions}, in parallel over the groups, and resolves where the seed of the calibration comes from.
#'
#' @param technicalBy See \code{\link{loadConsensusPeaks}}.
#' @param seqlevelsStyle See \code{\link{loadConsensusPeaks}}.
#' @param BPPARAM See \code{\link{loadConsensusPeaks}}.
#' @param seed See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#' @param biologicalPeaks List with the peaks of every biological sample, after the exclusion.
#' @param biologicalTable Data.frame with one row per biological sample and its group.
#' @param consensusArguments List with the arguments passed on to \code{consensusRegions::runConsensus()}.
#' @param groupLevels Character vector with the groups, in the order they are reported.
#' @param sampleTable Data.frame with one row per library.
#'
#' @return A list with \code{groupConsensus} and \code{consensusObjects}, named by group, \code{seedSource} and \code{recordedSeed}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel bpRNGseed bplapply
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges reduce
#' @importFrom S4Vectors mcols
#' @importFrom consensusRegions consensusRanges
#' @importFrom dplyr filter
#' @importFrom methods is
#' @importFrom utils packageVersion
#'
#' @noRd

.groupConsensusList <-
    function(technicalBy,
             seqlevelsStyle,
             BPPARAM,
             seed,
             verbose,
             biologicalPeaks,
             biologicalTable,
             consensusArguments,
             groupLevels,
             sampleTable) {
        groupConsensus <- list()
        consensusObjects <- list()

        groupSizes <- vapply(groupLevels, function(groupName) {
            sum(biologicalTable$group == groupName)
        }, integer(1))
        multiGroups <- groupLevels[groupSizes > 1]

        # Only the calibration draws random numbers, and the record says where they came from.
        # A back end of the user is handed over untouched, so its own seed, or the lack of one, decides
        userBPPARAM <- if (methods::is(BPPARAM, "BiocParallelParam")) {
            BPPARAM
        } else {
            NULL
        }
        calibrated <- isTRUE(consensusArguments$calibrate) & length(multiGroups) > 0
        userRNGseed <- if (calibrated & !is.null(userBPPARAM)) {
            BiocParallel::bpRNGseed(userBPPARAM)
        } else {
            NULL
        }

        if (calibrated & !is.null(seed) & !is.null(userBPPARAM)) {
            warning("The 'seed' parameter is ignored when 'BPPARAM' is a BiocParallelParam object: give the seed to the back end as its RNGseed.", call. = FALSE)
        }

        seedSource <- if (!calibrated) {
            "none"
        } else if (!is.null(userRNGseed)) {
            "BPPARAM"
        } else if (!is.null(seed) & is.null(userBPPARAM)) {
            "seed"
        } else {
            "unseeded"
        }

        recordedSeed <- switch(seedSource,
            "seed" = seed,
            "BPPARAM" = userRNGseed,
            NULL
        )

        # Earlier versions of consensusRegions would take the seed for an argument of buildConsensus, and stop there
        if (seedSource == "seed" && utils::packageVersion("consensusRegions") < "0.99.2") {
            stop("The 'seed' parameter needs consensusRegions 0.99.2 or later, which passes it on to the calibration. ",
                "Update consensusRegions, or leave 'seed' out.",
                call. = FALSE
            )
        }

        if (isTRUE(verbose) & seedSource == "unseeded") {
            message(
                "The calibration draws random positions and no seed was given: the consensus may change from one call to the next. ",
                "Set 'seed' to build it again identically."
            )
        }

        # A group left without biological samples has no consensus
        for (groupName in groupLevels[groupSizes == 0]) {
            groupConsensus[[groupName]] <- GenomicRanges::GRanges()
            consensusObjects[groupName] <- list(NULL)
        }

        # A single sample has no replicate to agree with, its peaks stand for the group
        for (groupName in groupLevels[groupSizes == 1]) {
            groupSamples <- dplyr::filter(biologicalTable, .data$group == groupName)$biological
            groupConsensus[[groupName]] <- IRanges::reduce(biologicalPeaks[[groupSamples]], ignore.strand = TRUE)
            consensusObjects[groupName] <- list(NULL)
        }

        # How the cores are shared between the groups and the calibration inside each of them
        workerPlan <- .consensusWorkers(
            BPPARAM = BPPARAM,
            nGroups = length(multiGroups),
            calibrate = isTRUE(consensusArguments$calibrate)
        )

        groupArgumentList <-
            lapply(
                multiGroups,
                function(groupName) {
                    groupSamples <- dplyr::filter(biologicalTable, .data$group == groupName)$biological

                    groupArguments <- consensusArguments
                    groupArguments$BPPARAM <- workerPlan$inner
                    if (is.null(groupArguments$verbose)) {
                        groupArguments$verbose <- FALSE
                    }

                    # Every group calibrates with the same seed, so that its consensus depends neither on the number of
                    # cores nor on the other groups. consensusRegions hands it to its back end as the RNGseed
                    if (seedSource == "seed") {
                        groupArguments$seed <- seed
                    }

                    # runConsensus renames the chromosomes to UCSC on its own, which would leave the consensus in
                    # one style and the peaks it came from in another, and every overlap between the two empty
                    if (!("seqlevelsStyle" %in% names(groupArguments))) {
                        groupArguments["seqlevelsStyle"] <- list(seqlevelsStyle)
                    }

                    # Weights from the libraries need the BAM files, which the sheet already lists. With technicalBy this is refused upstream
                    weightMethod <- groupArguments$weightMethod
                    if (!is.null(weightMethod) && weightMethod[1] %in% c("frip", "librarySize") &
                        is.null(groupArguments$bamFiles) & is.null(groupArguments$frip) & is.null(groupArguments$librarySize)) {
                        groupArguments$bamFiles <- sampleTable$bam[match(groupSamples, sampleTable$sample)]
                    }

                    return(c(list(peaks = biologicalPeaks[groupSamples], sampleNames = groupSamples), groupArguments))
                }
            )
        names(groupArgumentList) <- multiGroups

        if (isTRUE(verbose) & length(multiGroups) > 1 & workerPlan$nOuter > 1) {
            message(
                "Building the consensus of ", length(multiGroups), " groups on ", workerPlan$nOuter, " workers",
                if (workerPlan$nInner > 1) {
                    paste(",", workerPlan$nInner, "workers each for the calibration")
                } else {
                    ""
                }, "."
            )
        }

        # The groups do not depend on each other, so they run side by side
        groupObjects <- BiocParallel::bplapply(
            X = groupArgumentList,
            FUN = .runGroupConsensus,
            BPPARAM = .resolveBPPARAM(
                BPPARAM = workerPlan$outer,
                tasks = length(multiGroups)
            )
        )

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
            sampleLabel <- if (is.null(technicalBy)) {
                "sample"
            } else {
                "biological sample"
            }
            for (groupName in groupLevels) {
                libraryNote <- if (is.null(technicalBy)) {
                    ""
                } else {
                    paste0(" (", sum(sampleTable$group == groupName), " libraries)")
                }
                if (groupSizes[[groupName]] == 1) {
                    message("Group ", groupName, ": a single ", sampleLabel, libraryNote, ", its ", length(groupConsensus[[groupName]]), " peaks are taken as they are.")
                } else {
                    message("Group ", groupName, ": ", groupSizes[[groupName]], " ", sampleLabel, "s", libraryNote, ", ", length(groupConsensus[[groupName]]), " consensus regions.")
                }
            }
        }

        return(list(
            consensusObjects = consensusObjects,
            groupConsensus = groupConsensus,
            recordedSeed = recordedSeed,
            seedSource = seedSource
        ))
    }


#' @title .consensusRegionList
#'
#' @description Turns the total consensus into the list of region sets of \code{\link{loadConsensusPeaks}}, alone or split or replaced by the regions of the user, and carries the occupancy columns onto every set.
#'
#' @param technicalBy See \code{\link{loadConsensusPeaks}}.
#' @param regionSets See \code{\link{loadConsensusPeaks}}.
#' @param regionMode See \code{\link{loadConsensusPeaks}}.
#' @param unassignedSet See \code{\link{loadConsensusPeaks}}.
#' @param seqlevelsStyle See \code{\link{loadConsensusPeaks}}.
#' @param genomeAssembly See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#' @param biologicalPeaks List with the peaks of every biological sample, after the exclusion.
#' @param groupConsensus Named list with the consensus of every group.
#' @param peakList List with the peaks of every library, after the exclusion.
#' @param totalConsensus \code{GRanges} with the pooled consensus of the groups.
#'
#' @return A list with \code{regionList}, the named list of region sets, and \code{analysisMode}, one of \code{"consensus"}, \code{"split"} or \code{"replace"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors mcols
#'
#' @noRd

.consensusRegionList <-
    function(technicalBy,
             regionSets,
             regionMode,
             unassignedSet,
             seqlevelsStyle,
             genomeAssembly,
             verbose,
             biologicalPeaks,
             groupConsensus,
             peakList,
             totalConsensus) {
        analysisMode <- "consensus"
        regionList <- list(consensus = totalConsensus)

        if (!is.null(regionSets)) {
            userList <- loadRegions(
                regions = regionSets,
                keepMetadata = TRUE,
                seqlevelsStyle = seqlevelsStyle,
                genomeAssembly = genomeAssembly,
                outputFormat = "list",
                verbose = FALSE
            )

            # With seqlevelsStyle = NULL the regions of the user may be written in another style than the peaks.
            # Regions that cannot be reconciled are left as they are, and the occupancy step says so
            userList <- lapply(
                userList,
                function(userRanges) {
                    tryCatch(
                        expr = .matchSeqlevels(
                            x = userRanges,
                            targetSeqlevels = GenomeInfoDb::seqlevels(totalConsensus),
                            verbose = FALSE
                        ),
                        error = function(e) {
                            return(userRanges)
                        }
                    )
                }
            )

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
        regionList <- lapply(
            regionList,
            function(regionRanges) {
                occupancyTable <- .peakOccupancy(
                    regionRanges = regionRanges,
                    groupConsensus = groupConsensus,
                    peakList = peakList,
                    biologicalPeaks = if (is.null(technicalBy)) {
                        NULL
                    } else {
                        biologicalPeaks
                    }
                )
                S4Vectors::mcols(regionRanges) <- cbind(S4Vectors::mcols(regionRanges), occupancyTable)
                return(regionRanges)
            }
        )

        return(list(
            analysisMode = analysisMode,
            regionList = regionList
        ))
    }


#' @title .checkLoadConsensusPeaksArguments
#'
#' @description Checks the arguments of \code{\link{loadConsensusPeaks}} and brings them to the form the rest of the function works with.
#'
#' @param sampleSheet See \code{\link{loadConsensusPeaks}}.
#' @param groupBy See \code{\link{loadConsensusPeaks}}.
#' @param technicalBy See \code{\link{loadConsensusPeaks}}.
#' @param technicalArguments See \code{\link{loadConsensusPeaks}}.
#' @param regionMode See \code{\link{loadConsensusPeaks}}.
#' @param unassignedSet See \code{\link{loadConsensusPeaks}}.
#' @param seed See \code{\link{loadConsensusPeaks}}.
#' @param weightMethod See \code{\link{loadConsensusPeaks}}.
#' @param calibrate See \code{\link{loadConsensusPeaks}}.
#' @param nPermutations See \code{\link{loadConsensusPeaks}}.
#' @param targetFDR See \code{\link{loadConsensusPeaks}}.
#' @param recentre See \code{\link{loadConsensusPeaks}}.
#' @param width See \code{\link{loadConsensusPeaks}}.
#' @param combinationMethod See \code{\link{loadConsensusPeaks}}.
#' @param replicateType See \code{\link{loadConsensusPeaks}}.
#' @param stringencyThreshold See \code{\link{loadConsensusPeaks}}.
#' @param weakThreshold See \code{\link{loadConsensusPeaks}}.
#' @param combinedThreshold See \code{\link{loadConsensusPeaks}}.
#' @param alpha See \code{\link{loadConsensusPeaks}}.
#' @param minSupport See \code{\link{loadConsensusPeaks}}.
#' @param minReplicates See \code{\link{loadConsensusPeaks}}.
#' @param minOverlap See \code{\link{loadConsensusPeaks}}.
#' @param minOverlapFraction See \code{\link{loadConsensusPeaks}}.
#' @param mergeMethod See \code{\link{loadConsensusPeaks}}.
#' @param maxConsensusWidth See \code{\link{loadConsensusPeaks}}.
#' @param verbose See \code{\link{loadConsensusPeaks}}.
#' @param ... See \code{\link{loadConsensusPeaks}}.
#'
#' @return A list with the arguments as checked and completed: \code{consensusArguments}, \code{technicalArguments}, the arguments of the step joining the technical replicates (\code{NULL} without \code{technicalBy}), \code{regionMode}, \code{sampleSheet}.
#'
#' @author Sebastian Gregoricchio
#'
#' @noRd

.checkLoadConsensusPeaksArguments <-
    function(sampleSheet,
             groupBy,
             technicalBy,
             technicalArguments,
             regionMode,
             unassignedSet,
             seed,
             weightMethod,
             calibrate,
             nPermutations,
             targetFDR,
             recentre,
             width,
             combinationMethod,
             replicateType,
             stringencyThreshold,
             weakThreshold,
             combinedThreshold,
             alpha,
             minSupport,
             minReplicates,
             minOverlap,
             minOverlapFraction,
             mergeMethod,
             maxConsensusWidth,
             verbose,
             ...) {
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

        if (!is.null(technicalBy)) {
            if (!is.character(technicalBy) | length(technicalBy) != 1 || !(technicalBy %in% colnames(sampleSheet))) {
                stop("The 'technicalBy' parameter must name a column of the sample sheet.", call. = FALSE)
            }
            if (technicalBy %in% c("sample", groupBy)) {
                stop("The 'technicalBy' parameter must name the biological sample of every library, a column other than 'sample' and 'groupBy'.", call. = FALSE)
            }
        }

        if (!is.null(technicalArguments) && (!is.list(technicalArguments) || is.null(names(technicalArguments)) || any(names(technicalArguments) == ""))) {
            stop("The 'technicalArguments' parameter must be a named list, or NULL.", call. = FALSE)
        }

        regionMode <- tolower(as.character(regionMode[1]))
        if (!(regionMode %in% c("split", "replace"))) {
            stop("The 'regionMode' parameter must be either 'split' or 'replace'.", call. = FALSE)
        }

        if (!is.null(unassignedSet) && (!is.character(unassignedSet) | length(unassignedSet) != 1)) {
            stop("The 'unassignedSet' parameter must be a single string, or NULL.", call. = FALSE)
        }

        # The arguments of consensusRegions named here go to every group, next to any other one given through '...'
        consensusArguments <- c(
            list(
                weightMethod = weightMethod,
                calibrate = calibrate,
                nPermutations = nPermutations,
                targetFDR = targetFDR,
                recentre = recentre,
                width = width,
                combinationMethod = combinationMethod,
                replicateType = replicateType,
                stringencyThreshold = stringencyThreshold,
                weakThreshold = weakThreshold,
                combinedThreshold = combinedThreshold,
                alpha = alpha,
                minSupport = minSupport,
                minReplicates = minReplicates,
                minOverlap = minOverlap,
                minOverlapFraction = minOverlapFraction,
                mergeMethod = mergeMethod,
                maxConsensusWidth = maxConsensusWidth
            ),
            list(...)
        )

        if ("excludeRegions" %in% names(consensusArguments)) {
            stop("The 'excludeRegions' parameter has been replaced by 'blacklist' and 'greylist', which take the two lists separately.", call. = FALSE)
        }

        if (!is.null(seed) && (!is.numeric(seed) | length(seed) != 1 || !is.finite(seed))) {
            stop("The 'seed' parameter must be a single number, or NULL.", call. = FALSE)
        }

        reservedArguments <- intersect(names(consensusArguments), c("peaks", "sampleNames"))
        if (length(reservedArguments) > 0) {
            stop("The following arguments of runConsensus are set from the sample sheet and cannot be given: ",
                paste(reservedArguments, collapse = ", "), ".",
                call. = FALSE
            )
        }

        #-------------------------------#
        # Technical replicates          #
        #-------------------------------#
        if (!is.null(technicalBy)) {
            # The weights of the sheet describe libraries, the consensus of a group is built on biological samples
            weightSource <- c("bamFiles", "frip", "librarySize")
            if (!is.null(weightMethod) && weightMethod[1] %in% c("frip", "librarySize") & !any(weightSource %in% names(consensusArguments))) {
                stop("With 'technicalBy' the groups hold biological samples, and the BAM files of the sheet describe libraries. ",
                    "Give 'frip' or 'librarySize' per biological sample through '...', or use weightMethod = 'equal' or 'intrinsic'.",
                    call. = FALSE
                )
            }

            reservedTechnical <- intersect(
                names(technicalArguments),
                c(
                    "peaks", "sampleNames", "calibrate", "calibration", "nPermutations", "targetFDR", "recentre", "width",
                    "seed", "BPPARAM", "outputFile", "excludeRegions", "seqlevelsStyle"
                )
            )
            if (length(reservedTechnical) > 0) {
                stop("The following arguments cannot be given to the step joining the technical replicates: ",
                    paste(reservedTechnical, collapse = ", "), ".",
                    call. = FALSE
                )
            }

            # Inherited from the groups, apart from what only makes sense between biological samples
            technicalDefaults <- consensusArguments[setdiff(
                names(consensusArguments),
                c(
                    "calibrate", "nPermutations", "targetFDR", "recentre", "width", "minSupport", "minReplicates",
                    "minSupportWeight", "bamFiles", "frip", "librarySize", "seqlevelsStyle"
                )
            )]
            technicalDefaults$weightMethod <- "equal"
            technicalDefaults$replicateType <- "technical"
            technicalDefaults[names(technicalArguments)] <- technicalArguments
            technicalArguments <- technicalDefaults
        } else {
            technicalArguments <- NULL
        }

        return(list(
            consensusArguments = consensusArguments,
            technicalArguments = technicalArguments,
            regionMode = regionMode,
            sampleSheet = sampleSheet
        ))
    } # END function


#' @title .consensusWorkers
#'
#' @description Shares the cores between the groups built in parallel and the calibration run inside each group.
#'
#' @param BPPARAM A \code{BiocParallelParam} object, a number of cores, or \code{NULL} for one core.
#' @param nGroups Number of groups needing a consensus, the groups of a single sample left out.
#' @param calibrate Logical value to indicate whether the threshold is calibrated inside each group.
#'
#' @return A list with \code{outer}, the back end or the number of cores running the groups, \code{inner}, the back end or the number of cores each group gets for its permutations, and \code{nOuter} and \code{nInner}, the number of workers of each.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#' @importFrom BiocParallel bpnworkers
#'
#' @noRd

.consensusWorkers <-
    function(BPPARAM,
             nGroups,
             calibrate = FALSE) {
        # A back end of the user is never split: it goes to the calibration when there is one, to the groups otherwise
        if (methods::is(BPPARAM, "BiocParallelParam")) {
            nWorkers <- BiocParallel::bpnworkers(BPPARAM)
            if (isTRUE(calibrate)) {
                return(list(outer = 1L, inner = BPPARAM, nOuter = 1L, nInner = nWorkers))
            }
            return(list(outer = BPPARAM, inner = 1L, nOuter = nWorkers, nInner = 1L))
        }

        nCores <- if (is.null(BPPARAM)) {
            1L
        } else {
            BPPARAM
        }
        if (!is.numeric(nCores) || length(nCores) != 1 || !is.finite(nCores) || nCores < 1) {
            stop("The 'BPPARAM' parameter must be a BiocParallelParam object, a number of cores, or NULL.", call. = FALSE)
        }
        nCores <- as.integer(nCores)

        outerWorkers <- as.integer(max(1L, min(nCores, nGroups)))

        # Without calibration a consensus runs on one core, so the spare cores would only sit idle
        innerWorkers <- if (isTRUE(calibrate)) {
            as.integer(max(1L, nCores %/% outerWorkers))
        } else {
            1L
        }

        return(list(outer = outerWorkers, inner = innerWorkers, nOuter = outerWorkers, nInner = innerWorkers))
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
#' @noRd

.runGroupConsensus <-
    function(groupArguments) {
        return(do.call(what = consensusRegions::runConsensus, args = groupArguments))
    } # END function


#' @title .joinTechnicalReplicates
#'
#' @description Joins the libraries of every biological sample into one peak set through \code{consensusRegions::runConsensus}, before the consensus of the groups of \code{\link{loadConsensusPeaks}}.
#'
#' @param sampleTable Data.frame with one row per library and the \code{sample}, \code{group} and \code{biological} columns.
#' @param peakList Named \code{GRangesList} with the peaks of every library, after the blacklist and the greylist.
#' @param technicalArguments List with the arguments of \code{consensusRegions::runConsensus} for this step.
#' @param seqlevelsStyle String with the chromosome naming style of the peaks, or \code{NULL}.
#' @param BPPARAM A \code{BiocParallelParam} object, a number of cores, or \code{NULL}. Default: \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with \code{peaks}, a named \code{GRangesList} with one peak set per biological sample, and \code{objects}, a named list with the \code{consensusRegions} object behind each of them, \code{NULL} for the biological samples with a single library.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRangesList
#' @importFrom BiocParallel bplapply
#'
#' @noRd

.joinTechnicalReplicates <-
    function(sampleTable,
             peakList,
             technicalArguments,
             seqlevelsStyle,
             BPPARAM = NULL,
             verbose = TRUE) {
        biologicalNames <- unique(sampleTable$biological)
        libraryNumber <- vapply(biologicalNames, function(biologicalName) {
            sum(sampleTable$biological == biologicalName)
        }, integer(1))
        multiSamples <- biologicalNames[libraryNumber > 1]

        if (length(multiSamples) == 0) {
            warning("Every biological sample of 'technicalBy' has a single library, there is no technical replicate to join.", call. = FALSE)
        }

        #-------------------------------#
        # Arguments of every sample     #
        #-------------------------------#
        technicalArguments["seqlevelsStyle"] <- list(seqlevelsStyle)
        if (is.null(technicalArguments$verbose)) {
            technicalArguments$verbose <- FALSE
        }

        argumentList <- lapply(multiSamples, function(biologicalName) {
            libraryNames <- dplyr::filter(sampleTable, .data$biological == biologicalName)$sample
            return(c(list(peaks = peakList[libraryNames], sampleNames = libraryNames), technicalArguments))
        })
        names(argumentList) <- multiSamples

        # Nothing is calibrated here, so every core can take a biological sample
        workerPlan <- .consensusWorkers(BPPARAM = BPPARAM, nGroups = length(multiSamples), calibrate = FALSE)

        technicalObjects <- BiocParallel::bplapply(
            X = argumentList,
            FUN = .runGroupConsensus,
            BPPARAM = .resolveBPPARAM(BPPARAM = workerPlan$outer, tasks = length(multiSamples))
        )

        #-------------------------------#
        # One peak set per sample       #
        #-------------------------------#
        biologicalPeaks <- lapply(biologicalNames, function(biologicalName) {
            libraryNames <- dplyr::filter(sampleTable, .data$biological == biologicalName)$sample
            if (length(libraryNames) == 1) {
                return(.asBiologicalPeaks(peakRanges = peakList[[libraryNames]], libraryPeaks = NULL))
            }
            return(.asBiologicalPeaks(
                peakRanges = consensusRegions::consensusRanges(technicalObjects[[biologicalName]]),
                libraryPeaks = peakList[libraryNames]
            ))
        })
        names(biologicalPeaks) <- biologicalNames

        if (isTRUE(verbose)) {
            message(
                "Technical replicates joined into ", length(biologicalNames), " biological samples: ",
                paste(biologicalNames, " (", libraryNumber, " libraries, ", lengths(biologicalPeaks), " peaks)", sep = "", collapse = ", "), "."
            )
        }

        objectList <- lapply(biologicalNames, function(biologicalName) {
            if (biologicalName %in% multiSamples) {
                return(technicalObjects[[biologicalName]])
            }
            return(NULL)
        })
        names(objectList) <- biologicalNames

        return(list(peaks = GenomicRanges::GRangesList(biologicalPeaks), objects = objectList))
    } # END function


#' @title .asBiologicalPeaks
#'
#' @description Brings the peaks of a biological sample to the columns the consensus of a group reads: \code{pValue}, the significance on the -log10 scale, and \code{peak}, the summit as an offset from the start, as in a narrowPeak file.
#'
#' @param peakRanges \code{GRanges} with the peaks of a single library, carrying \code{negLog10P}, or with the consensus of several libraries, carrying \code{combinedNegLog10P}.
#' @param libraryPeaks Named \code{GRangesList} with the peaks of the libraries behind a consensus, used to place its summits, or \code{NULL} for a single library.
#'
#' @return A \code{GRanges} with the \code{pValue} and \code{peak} columns only. \code{peak} is -1 where the peaks carry no summit.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom BiocGenerics start width
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
#' @importFrom dplyr tibble group_by slice_max ungroup %>%
#' @importFrom rlang .data
#'
#' @noRd

.asBiologicalPeaks <-
    function(peakRanges,
             libraryPeaks = NULL) {
        metadataTable <- S4Vectors::mcols(peakRanges)

        #-------------------------------#
        # A single library              #
        #-------------------------------#
        if (is.null(libraryPeaks)) {
            summitOffset <- if ("peak" %in% colnames(metadataTable)) {
                as.integer(metadataTable$peak)
            } else {
                rep(-1L, length(peakRanges))
            }
            S4Vectors::mcols(peakRanges) <- NULL
            S4Vectors::mcols(peakRanges)$pValue <- as.numeric(metadataTable$negLog10P)
            S4Vectors::mcols(peakRanges)$peak <- summitOffset
            return(peakRanges)
        }

        #-------------------------------#
        # A consensus of libraries      #
        #-------------------------------#
        combinedSignificance <- if ("combinedNegLog10P" %in% colnames(metadataTable)) {
            as.numeric(metadataTable$combinedNegLog10P)
        } else {
            rep(NA_real_, length(peakRanges))
        }

        S4Vectors::mcols(peakRanges) <- NULL
        S4Vectors::mcols(peakRanges)$pValue <- combinedSignificance
        S4Vectors::mcols(peakRanges)$peak <- rep(-1L, length(peakRanges))

        pooledPeaks <- unlist(libraryPeaks, use.names = FALSE)
        if (length(peakRanges) == 0 | !("peak" %in% colnames(S4Vectors::mcols(pooledPeaks)))) {
            return(peakRanges)
        }

        # The summit of a region is the one of its most significant library peak
        summitOffset <- as.integer(S4Vectors::mcols(pooledPeaks)$peak)
        hasSummit <- !is.na(summitOffset) & summitOffset >= 0
        summitPosition <- ifelse(hasSummit,
            BiocGenerics::start(pooledPeaks) + summitOffset,
            BiocGenerics::start(pooledPeaks) + BiocGenerics::width(pooledPeaks) %/% 2L
        )

        hits <- GenomicRanges::findOverlaps(peakRanges, pooledPeaks, ignore.strand = TRUE)
        if (length(hits) == 0) {
            return(peakRanges)
        }

        bestHits <- dplyr::tibble(
            region = S4Vectors::queryHits(hits),
            negLog10P = as.numeric(S4Vectors::mcols(pooledPeaks)$negLog10P)[S4Vectors::subjectHits(hits)],
            summit = summitPosition[S4Vectors::subjectHits(hits)]
        ) %>%
            dplyr::group_by(.data$region) %>%
            dplyr::slice_max(order_by = .data$negLog10P, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup()

        # Kept inside the region, a summit from a library peak reaching past its edge would fall outside
        regionStart <- BiocGenerics::start(peakRanges)[bestHits$region]
        regionWidth <- BiocGenerics::width(peakRanges)[bestHits$region]
        S4Vectors::mcols(peakRanges)$peak[bestHits$region] <- as.integer(pmin(pmax(bestHits$summit - regionStart, 0L), regionWidth - 1L))

        return(peakRanges)
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
#' @noRd

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

        exclusionList <- loadRegions(
            regions = excludeRegions,
            keepMetadata = FALSE,
            duplicatedSets = "remove",
            seqlevelsStyle = seqlevelsStyle,
            outputFormat = "list",
            verbose = FALSE
        )

        return(IRanges::reduce(unlist(GenomicRanges::GRangesList(unname(exclusionList)), use.names = FALSE), ignore.strand = TRUE))
    } # END function


#' @title .checkListAssembly
#'
#' @description Refuses a list of regions built for another assembly than the one declared for the peaks, which would overlap them on the chromosome names alone.
#'
#' @param listInput The list as given: a \code{GRanges}, a path, a data.frame, or a list of them.
#' @param genomeAssembly String with the assembly of the peaks, or \code{NULL}.
#' @param listLabel String with the name of the parameter the list came from.
#' @param regionLabel String naming what the list is applied to in the message. Default: \code{"peaks"}.
#'
#' @return Nothing, it stops when the two assemblies differ.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb genome
#' @importFrom methods is
#'
#' @noRd

.checkListAssembly <-
    function(listInput,
             genomeAssembly,
             listLabel,
             regionLabel = "peaks") {
        if (is.null(genomeAssembly) || is.na(genomeAssembly[1]) || genomeAssembly[1] == "") {
            return(invisible(TRUE))
        }

        # Only a GRanges carries its assembly, a file or a table cannot be checked
        listElements <- if (methods::is(listInput, "GRanges")) {
            list(listInput)
        } else if (is.list(listInput) & !is.data.frame(listInput)) {
            listInput
        } else {
            list()
        }

        for (listElement in listElements) {
            if (!methods::is(listElement, "GRanges")) {
                next
            }

            listAssembly <- unique(as.character(GenomeInfoDb::genome(listElement)))
            listAssembly <- listAssembly[!is.na(listAssembly) & listAssembly != ""]

            # GRCh38 and hg38 name the same assembly, the aliases are resolved before the comparison
            if (length(listAssembly) == 1 && .resolveGenomeName(listAssembly) != .resolveGenomeName(as.character(genomeAssembly[1]))) {
                stop("The ", listLabel, " was built for ", listAssembly, " and the ", regionLabel, " are declared as ", genomeAssembly[1],
                    ". Overlapping them would match the chromosome names and nothing else. Clear the assembly with genome() on the list to force it.",
                    call. = FALSE
                )
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
#' @noRd

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
        splitList <- lapply(setLevels, function(setName) {
            totalConsensus[which(setAssignment == setName)]
        })
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
#' @param biologicalPeaks Named \code{GRangesList} with the peak set of every biological sample, or \code{NULL} when every library is one. Default: \code{NULL}.
#'
#' @return A \code{DataFrame} with one logical \code{peak.<group>} column per group, \code{peak.groups} and \code{peak.samples}, and \code{peak.biological} when \code{biologicalPeaks} is given.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom utils head
#'
#' @noRd

.peakOccupancy <-
    function(regionRanges,
             groupConsensus,
             peakList,
             biologicalPeaks = NULL) {
        regionNumber <- length(regionRanges)

        # Two naming styles meeting here would give an empty overlap for every sample rather than an error
        sharedSeqlevels <- intersect(
            GenomeInfoDb::seqlevels(regionRanges),
            unique(unlist(lapply(as.list(peakList), GenomeInfoDb::seqlevels), use.names = FALSE))
        )

        if (length(sharedSeqlevels) == 0) {
            stop("The regions and the peaks of the samples share no chromosome, their naming styles differ: ",
                paste(utils::head(GenomeInfoDb::seqlevels(regionRanges), 3), collapse = ", "), " against ",
                paste(utils::head(GenomeInfoDb::seqlevels(peakList[[1]]), 3), collapse = ", "), ".",
                call. = FALSE
            )
        }

        # matrix() keeps the shape when there is a single region, which vapply would flatten
        groupMatrix <- matrix(
            vapply(
                groupConsensus,
                function(groupRanges) {
                    IRanges::overlapsAny(regionRanges, groupRanges, ignore.strand = TRUE)
                },
                logical(regionNumber)
            ),
            nrow = regionNumber
        )

        sampleMatrix <- matrix(
            vapply(
                as.list(peakList),
                function(peakRanges) {
                    IRanges::overlapsAny(regionRanges, peakRanges, ignore.strand = TRUE)
                },
                logical(regionNumber)
            ),
            nrow = regionNumber
        )

        colnames(groupMatrix) <- paste0("peak.", make.names(names(groupConsensus)))

        occupancyTable <- S4Vectors::DataFrame(groupMatrix, check.names = FALSE)
        occupancyTable$peak.groups <- as.integer(rowSums(groupMatrix))
        occupancyTable$peak.samples <- as.integer(rowSums(sampleMatrix))

        # With technical replicates the libraries and the biological samples are two different counts
        if (!is.null(biologicalPeaks)) {
            biologicalMatrix <- matrix(
                vapply(
                    as.list(biologicalPeaks),
                    function(peakRanges) {
                        IRanges::overlapsAny(regionRanges, peakRanges, ignore.strand = TRUE)
                    },
                    logical(regionNumber)
                ),
                nrow = regionNumber
            )
            occupancyTable$peak.biological <- as.integer(rowSums(biologicalMatrix))
        }

        return(occupancyTable)
    } # END function


#' @title consensusData
#'
#' @description Returns the consensus data kept by a \code{RegionSetDE} object built from peaks with \code{\link{loadConsensusPeaks}}.
#'
#' @param object \code{RegionSetDE} object.
#'
#' @return A list with \code{groups}, the consensus regions of every group; \code{objects}, the \code{consensusRegions} object behind each of them (\code{NULL} for the groups with a single sample); \code{technical}, the \code{consensusRegions} object joining the libraries of every biological sample (\code{NULL} for the samples with a single library), an empty list without \code{technicalBy}; \code{total}, the pooled consensus; \code{peaks}, the peaks of every sample after the blacklist and the greylist; \code{removed}, the peaks they took out, with the \code{sample} they came from and the list that removed them in \code{removed.by}; \code{samples}, a data.frame with the group, the biological sample (the library itself without \code{technicalBy}), the peak file, the number of peaks read (\code{n.peaks}), removed by the blacklist (\code{n.blacklist}), by the greylist (\code{n.greylist}) and by both (\code{n.excluded}) for every sample; \code{sheet}, the sample sheet the consensus was built from; \code{groupBy}; \code{technicalBy}; \code{mode}, one among \code{"consensus"}, \code{"split"} and \code{"replace"}; and \code{blacklist} and \code{greylist}, the regions the peaks were cleaned against, \code{NULL} when a list was not given.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'     peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
#'         pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE
#'     )
#'
#'     sampleSheet <- loadSampleSheet(
#'         data.frame(
#'             sample = c("A_1", "A_2", "B_1", "B_2"),
#'             bam = c("A_1.bam", "A_2.bam", "B_1.bam", "B_2.bam"),
#'             peaks = peakFiles[c(1, 2, 2, 3)],
#'             condition = c("A", "A", "B", "B")
#'         ),
#'         checkFiles = FALSE, verbose = FALSE
#'     )
#'
#'     regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)
#'
#'     consensusData(regions)$samples
#'     lengths(consensusData(regions)$groups)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadConsensusPeaks}}, \code{\link{consensusGroupList}}, \code{\link{plotPeakUpset}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export

setGeneric(name = "consensusData", def = function(object) {
    standardGeneric("consensusData")
})


#' @rdname consensusData
#' @export

setMethod(
    f = "consensusData",
    signature = "RegionSetDE",
    definition = function(object) {
        if (length(object@consensus) == 0) {
            stop("The regions were not built from peaks, there is no consensus to return.", call. = FALSE)
        }
        return(object@consensus)
    }
)


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
#'     sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#'     consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
#'
#'     groupList <- consensusGroupList(consensus, seqlevelsStyle = "UCSC")
#'     lengths(groupList)
#'     groupList$R1881_24h
#'
#'     # Two groups only, in the order given
#'     names(consensusGroupList(consensus, groups = c("R1881_24h", "DMSO")))
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
        if (length(.storedConsensus(object)) == 0) {
            stop("The 'object' parameter must be the RegionSetDE object returned by loadConsensusPeaks().", call. = FALSE)
        }

        if (!is.logical(asGRangesList) | length(asGRangesList) != 1 | anyNA(asGRangesList)) {
            stop("The 'asGRangesList' parameter must be TRUE or FALSE.", call. = FALSE)
        }

        groupList <- consensusData(object)$groups

        if (!is.null(groups)) {
            absentGroups <- setdiff(groups, names(groupList))
            if (length(absentGroups) > 0) {
                stop("The following groups are absent from the consensus: ", paste(absentGroups, collapse = ", "),
                    ". The groups are: ", paste(names(groupList), collapse = ", "), ".",
                    call. = FALSE
                )
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
