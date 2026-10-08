# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title computeProfiles
#'
#' @description Computes the signal of every sample in bins around the centre of a group of regions, by default the regions found to change in a contrast split into those going up and those going down. The signal is read from the BAM files the counts came from, with the same read filters and scaled by the normalisation of the object, or from bigWig files. The result feeds \code{\link{plotProfile}}.
#'
#' @param object \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object carrying its counts, whose differential regions are profiled, or a \code{RegionSetDE.counts} or \code{RegionSetDE.fit} object, whose region sets are.
#' @param regions Regions to profile instead of those picked from \code{object}: a \code{GRanges}, a named list of \code{GRanges} or a \code{GRangesList}, one group of rows per element, or a \code{RegionSetDE} object, one group per set. Default: \code{NULL}.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param samples Samples profiled: a character vector with their names, a numeric vector with their positions, or a logical vector with one value per sample. The scaling factors of the whole analysis are kept. Default: \code{NULL}, all of them.
#' @param direction String with the differential regions profiled, one among \code{"both"}, which draws the regions going up and those going down as two groups of rows, \code{"up"} and \code{"down"}. Only for results. Default: \code{"both"}.
#' @param FDR Numeric value with the adjusted p-value cut-off defining a differential region. Default: \code{NULL}, the one used by the test.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off defining a differential region. Default: \code{NULL}, the one used by the test.
#' @param maxRegions Numeric value with the maximum number of regions per group of rows: the most significant ones for results, the strongest ones for the region sets of a counts object, and regions spread evenly along the list for regions given through \code{regions}. It is applied after \code{blacklist} and \code{whitelist}. Default: \code{1000}.
#' @param blacklist Regions whose signal must stay out of the figure: a row is left out when its drawn window, \code{distance} on each side of the centre, overlaps them. A \code{GRanges}, a path to a BED-like file, a data.frame, a list of them, or \code{TRUE} for the blacklist stored in the object by \code{\link{applyBlacklist}} or \code{\link{loadConsensusPeaks}}. Default: \code{NULL}.
#' @param whitelist Regions the rows are restricted to: a row is kept when its region overlaps them. Accepts the same forms as \code{blacklist}, except \code{TRUE}. Default: \code{NULL}.
#' @param groupBy String with the name of a \code{colData} column whose groups are averaged into one profile each, for instance \code{"condition"}. Default: \code{NULL}, one profile per sample.
#' @param signal String with the files the signal is read from: \code{"bam"}, the BAM files of \code{\link{countReads}}, \code{"bigwig"}, the bigWig files of the sample sheet or of \code{\link{countBigwig}}, or \code{"auto"}, the BAM files when they are known and the bigWig files otherwise. Default: \code{"auto"}.
#' @param signalFiles Character vector with one BAM or bigWig file per sample, in the order of the samples, overriding the ones recorded in the object. Default: \code{NULL}.
#' @param distance Numeric value with the number of base pairs drawn on each side of the centre. Default: \code{1500}.
#' @param binWidth Numeric value with the width of the bins the signal is averaged over, in base pairs. Default: \code{50}.
#' @param centre String with the point the regions are aligned on: \code{"summit"}, the \code{summit} column written by \code{\link{countReads}} with \code{summits}, falling back to the midpoint where there is none, or \code{"midpoint"}. Default: \code{"summit"}.
#' @param useOffsets Logical value to indicate whether the signal read from BAM files must be scaled by the scaling factors of \code{\link{normalizeCounts}}, \code{FALSE} scaling it by the library sizes alone. Ignored for bigWig files, whose signal is taken as it is. Default: \code{TRUE}.
#' @param nThreads Number of threads, one file per thread. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with \code{profiles}, a named list with one matrix per sample, or per group with \code{groupBy}, holding one row per region and one column per bin; \code{regions}, a \code{GRanges} with the window of every row, its \code{row.group} and its \code{region.key}; \code{bins}, the distance of the centre of every bin from the centre of the regions; \code{samples}, the annotation of the samples; and \code{parameters}.
#'
#' @details The windows span \code{distance} base pairs on each side of the centre and are cut into bins of \code{binWidth}, the signal of a bin being the mean coverage over its bases. Regions on the minus strand are read from right to left, so that a set of promoters shows the transcription start site facing the same way. Windows running past the end of a chromosome are left out, and how many were is reported.
#'
#' From BAM files the coverage is built from the fragments as \code{\link{countReads}} counts them: paired-end fragments from their pairs, single-end reads extended to the fragment length of their sample, with the same mapping quality, duplicate and discarded region filters. It is then divided by the scaling factor of the sample and expressed per million fragments of the average library, so that the profiles of different samples compare the way the normalised counts do. A normalisation with one offset per region, as \code{method = "loess"} gives, has no factor to divide by, and the library sizes are used instead. This is the counterpart of \code{DiffBind::dba.plotProfile}, which also reads the reads of the analysis and applies its normalisation.
#'
#' A bigWig file is read as it is, so it has to be normalised already, and to the same scale for every sample, for its profiles to be compared.
#'
#' A region cleaned of blacklisted regions can still sit next to one, and the window drawn around it then shows the artefact. \code{blacklist} therefore acts on the whole window rather than on the region, while \code{whitelist} asks the region itself to overlap the list, which is how a profile is restricted to promoters, enhancers or any other class of sites. How many rows each of them left out is reported.
#'
#' @examples
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#'
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
#'
#' peakProfiles <- computeProfiles(counts, groupBy = "condition", distance = 1000, verbose = FALSE)
#' names(peakProfiles$profiles)
#' dim(peakProfiles$profiles$DMSO)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotProfile}}, \code{\link{countReads}}, \code{\link{testRegions}}
#'
#' @importFrom SummarizedExperiment colData rowData rowRanges assay
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
#'
#' @export computeProfiles

computeProfiles <-
  function(object,
           regions = NULL,
           contrast = NULL,
           set = NULL,
           samples = NULL,
           direction = "both",
           FDR = NULL,
           log2FC = NULL,
           maxRegions = 1000,
           blacklist = NULL,
           whitelist = NULL,
           groupBy = NULL,
           signal = "auto",
           signalFiles = NULL,
           distance = 1500,
           binWidth = 50,
           centre = "summit",
           useOffsets = TRUE,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    resolvedObject <- .resolveCounts(object = object, counts = NULL, contrast = contrast)
    counts <- resolvedObject$counts
    results <- resolvedObject$results

    if (!(direction[1] %in% c("both", "up", "down"))) {
      stop("The 'direction' parameter must be one of 'both', 'up', 'down'.", call. = FALSE)
    }
    direction <- direction[1]

    if (!(centre[1] %in% c("summit", "midpoint"))) {
      stop("The 'centre' parameter must be either 'summit' or 'midpoint'.", call. = FALSE)
    }
    centre <- centre[1]

    if (!is.numeric(distance) | length(distance) != 1 || distance < 1) {
      stop("The 'distance' parameter must be a positive number of base pairs.", call. = FALSE)
    }

    if (!is.numeric(binWidth) | length(binWidth) != 1 || binWidth < 1) {
      stop("The 'binWidth' parameter must be a positive number of base pairs.", call. = FALSE)
    }

    # Whole bins on each side, so that the centre falls on the border between two of them
    binWidth <- as.integer(round(binWidth))
    binsPerSide <- max(1L, as.integer(round(distance / binWidth)))
    distance <- binsPerSide * binWidth

    if (!is.numeric(maxRegions) | length(maxRegions) != 1 || maxRegions < 1) {
      stop("The 'maxRegions' parameter must be a positive number.", call. = FALSE)
    }

    # The files given here follow the samples of the object, and are subset with them
    sampleIndex <- .sampleIndex(sampleNames = colnames(counts), samples = samples)
    if (!is.null(signalFiles) && length(signalFiles) == ncol(counts)) {
      signalFiles <- signalFiles[sampleIndex]
    }

    # The unit of the coverage is set on every sample of the object, so a subset keeps the values it has in the whole
    allCounts <- counts
    counts <- .subsetSamples(counts = counts, samples = samples)

    sampleTable <- .sampleAnnotation(counts = counts)

    if (!is.null(groupBy) && !(groupBy %in% colnames(sampleTable))) {
      stop("The column '", groupBy, "' is absent from the colData.", call. = FALSE)
    }

    #------------------------#
    # Regions to profile     #
    #------------------------#
    # Every candidate first, so that the rows filtered out are replaced by the next ones rather than lost
    profileRegions <- .profileRegions(counts = counts,
                                      results = results,
                                      regions = regions,
                                      set = set,
                                      direction = direction,
                                      FDR = FDR,
                                      log2FC = log2FC,
                                      maxRegions = Inf,
                                      centre = centre)

    #------------------------#
    # Blacklist, whitelist   #
    #------------------------#
    # The blacklist looks at the whole window drawn, the whitelist at the region itself
    if (isTRUE(blacklist)) {
      if (is.null(allCounts@blacklist)) {
        stop("blacklist = TRUE uses the blacklist stored in the object, and there is none: pass the list itself.", call. = FALSE)
      }
      blacklist <- allCounts@blacklist
    }

    profileSeqlevels <- unique(as.character(GenomeInfoDb::seqnames(profileRegions)))

    if (!is.null(blacklist)) {
      blacklistRanges <- .matchSeqlevels(x = .loadExclusionRegions(excludeRegions = blacklist, seqlevelsStyle = NULL, listLabel = "blacklist"),
                                         targetSeqlevels = profileSeqlevels, verbose = FALSE)
      profileWindows <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(profileRegions),
                                               ranges = IRanges::IRanges(start = profileRegions$centre.position - distance, width = 2L * distance))
      blacklistedRows <- IRanges::overlapsAny(profileWindows, blacklistRanges, ignore.strand = TRUE)

      if (isTRUE(verbose)) {
        message(sum(blacklistedRows), " regions whose window overlaps the blacklist were left out.")
      }
      profileRegions <- profileRegions[!blacklistedRows]
    }

    if (!is.null(whitelist)) {
      whitelistRanges <- .matchSeqlevels(x = .loadExclusionRegions(excludeRegions = whitelist, seqlevelsStyle = NULL, listLabel = "whitelist"),
                                         targetSeqlevels = profileSeqlevels, verbose = FALSE)
      whitelistedRows <- IRanges::overlapsAny(profileRegions, whitelistRanges, ignore.strand = TRUE)

      if (isTRUE(verbose)) {
        message(sum(!whitelistedRows), " regions outside the whitelist were left out.")
      }
      profileRegions <- profileRegions[whitelistedRows]
    }

    if (length(profileRegions) == 0) {
      stop("No region is left to profile after the blacklist and the whitelist.", call. = FALSE)
    }

    #------------------------#
    # Rows per group         #
    #------------------------#
    # The rows come ordered by priority, regions given by hand are thinned evenly along the list instead
    keptRows <- unlist(lapply(split(seq_along(profileRegions), factor(profileRegions$row.group, levels = unique(profileRegions$row.group))),
                              function(groupRows) {
                                if (length(groupRows) <= maxRegions) {return(groupRows)}
                                if (is.null(regions)) {return(groupRows[seq_len(maxRegions)])}
                                return(groupRows[.thinIndex(n = length(groupRows), maxPoints = maxRegions)])
                              }),
                       use.names = FALSE)
    profileRegions <- profileRegions[keptRows]

    #------------------------#
    # Files to read          #
    #------------------------#
    signalInfo <- .profileSignalFiles(counts = counts, signal = signal, signalFiles = signalFiles)

    if (isTRUE(verbose)) {
      message("Profiling ", length(profileRegions), " regions (",
              paste(names(table(profileRegions$row.group)), table(profileRegions$row.group), sep = " ", collapse = ", "),
              ") over ", length(signalInfo$files), " ", if (signalInfo$type == "bam") {"BAM"} else {"bigWig"}, " files, ",
              distance, " bp on each side in bins of ", binWidth, " bp...")
    }

    #------------------------#
    # Binned signal          #
    #------------------------#
    windowRanges <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(profileRegions),
                                           ranges = IRanges::IRanges(start = profileRegions$centre.position - distance,
                                                                     width = 2L * distance),
                                           strand = BiocGenerics::strand(profileRegions))
    S4Vectors::mcols(windowRanges) <- S4Vectors::mcols(profileRegions)

    binnedList <- .binnedSignal(files = signalInfo$files,
                                type = signalInfo$type,
                                windows = windowRanges,
                                binWidth = binWidth,
                                counts = counts,
                                nThreads = nThreads)

    keptWindows <- binnedList$kept
    if (isTRUE(verbose) & sum(!keptWindows) > 0) {
      message(sum(!keptWindows), " windows run past the end of their chromosome and were left out.")
    }

    if (sum(keptWindows) == 0) {
      stop("No window lies within the chromosomes of the signal files.", call. = FALSE)
    }

    windowRanges <- windowRanges[keptWindows]
    profileList <- lapply(binnedList$matrices, function(profileMatrix) {profileMatrix[keptWindows, , drop = FALSE]})
    names(profileList) <- sampleTable$sample

    #------------------------#
    # Scale the BAM signal   #
    #------------------------#
    signalLabel <- "Signal (bigWig)"

    if (signalInfo$type == "bam") {
      sampleScale <- .profileScaling(counts = allCounts, useOffsets = useOffsets, verbose = verbose)[sampleIndex]
      profileList <- lapply(seq_along(profileList), function(i) {profileList[[i]] / sampleScale[i]})
      names(profileList) <- sampleTable$sample
      signalLabel <- if (isTRUE(useOffsets)) {"Normalised coverage per million fragments"} else {"Coverage per million fragments"}
    }

    #------------------------#
    # Average the groups     #
    #------------------------#
    if (!is.null(groupBy)) {
      groupValues <- as.character(sampleTable[[groupBy]])
      groupLevels <- if (is.factor(sampleTable[[groupBy]])) {intersect(levels(sampleTable[[groupBy]]), groupValues)} else {unique(groupValues)}

      profileList <- lapply(groupLevels, function(groupName) {Reduce(f = `+`, x = profileList[groupValues == groupName]) / sum(groupValues == groupName)})
      names(profileList) <- groupLevels
    }

    binCentres <- seq(from = -distance + binWidth / 2, by = binWidth, length.out = 2L * binsPerSide)
    profileList <- lapply(profileList, function(profileMatrix) {
      dimnames(profileMatrix) <- list(windowRanges$region.key, as.character(binCentres))
      return(profileMatrix)
    })

    return(list(profiles = profileList,
                regions = windowRanges,
                bins = binCentres,
                samples = sampleTable,
                parameters = list(signal = signalInfo$type,
                                  signal.label = signalLabel,
                                  groupBy = groupBy,
                                  direction = direction,
                                  distance = distance,
                                  binWidth = binWidth,
                                  centre = centre,
                                  useOffsets = useOffsets,
                                  contrast = if (is.null(results)) {NULL} else {contrastName(results)})))
  } # END function




#' @title plotProfile
#'
#' @description Draws the signal around the centre of a group of regions, by default the regions going up and those going down in a contrast, as DiffBind does with \code{dba.plotProfile}. The \code{"heatmap"} style draws one heatmap per sample or group, one row per region and the mean profile on top; the \code{"lines"} style draws the mean profiles alone, with their standard error.
#'
#' @param object Any object accepted by \code{\link{computeProfiles}}, or the list it returns.
#' @param style String with the drawing, either \code{"heatmap"} or \code{"lines"}. Default: \code{"heatmap"}.
#' @param colours Character vector with the colours of the heatmap scale, from low to high. Default: \code{NULL}, \code{viridisLite::mako(100, direction = -1)}.
#' @param groupColours Character vector with the colours of the groups of rows in the heatmap, or of the samples or groups in the \code{"lines"} style, named after them or given in their order. Default: \code{NULL}, the palette of the package, with the regions going up in red and those going down in blue as in \code{\link{plotVolcano}}.
#' @param limits Numeric vector of length two with the range of the colour scale of the heatmap. Default: \code{NULL}, from zero to the 99th percentile of the values.
#' @param title String with the title of the \code{"lines"} style. A \code{HeatmapList} takes its title when it is drawn, as \code{ComplexHeatmap::draw(x, column_title = "title")}. Default: \code{NULL}, the contrast when there is one.
#' @param fontSize Numeric value with the font size of the heatmap. Default: \code{9}.
#' @param baseSize Numeric value with the base font size of the \code{"lines"} style. Default: \code{12}.
#' @param legendPosition String with the position of the legend of the \code{"lines"} style. Default: \code{"right"}.
#' @param ... Arguments passed to \code{\link{computeProfiles}} when \code{object} is not already a set of profiles, such as \code{groupBy}, \code{contrast}, \code{distance} or \code{signalFiles}. They are ignored when \code{object} is the list returned by \code{\link{computeProfiles}}, which already holds them.
#'
#' @return A \code{HeatmapList} built by \code{ComplexHeatmap}, drawn when printed, for the \code{"heatmap"} style, or a \code{ggplot} object for the \code{"lines"} style. \code{ComplexHeatmap::draw} opens the layout of the heatmaps, for instance \code{draw(x, column_title = "R1881 against DMSO", heatmap_legend_side = "bottom")}.
#'
#' @details Given a results, counts or fit object, \code{plotProfile()} calls \code{\link{computeProfiles}} and draws what it returns. The reading of the files is the slow part, so profiles drawn in both styles, or with different colours, are better computed once with \code{\link{computeProfiles}} and handed over as they are.
#'
#' The rows of each group are sorted by their mean signal over all the samples, strongest on top, so that the same row holds the same region in every heatmap, and the eye can run across the panels to see a region gain or lose signal. The colour scale is shared by the panels for the same reason. The profile on top of each heatmap is the mean of the rows of every group, drawn on a common axis.
#'
#' A change found by the test should show here as a difference in height between the conditions, centred on the summit. A difference spread evenly across the whole window points to background rather than binding, and a profile whose peak sits away from the centre says the regions were not aligned on the signal, which recentring with \code{countReads(summits = )} fixes.
#'
#' @examples
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#'
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, countInput = FALSE, verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
#' fit <- fitRegions(counts, design = ~ condition, verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)
#'
#' # The regions going up and down, one heatmap per condition
#' plotProfile(results, groupBy = "condition", distance = 1000, verbose = FALSE)
#'
#' # Computed once, drawn twice: the BAM files are read a single time
#' profileData <- computeProfiles(results, groupBy = "condition", distance = 1000, verbose = FALSE)
#' plotProfile(profileData)
#' plotProfile(profileData, style = "lines")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{computeProfiles}}, \code{\link{plotTopHeatmap}}, \code{\link{countReads}}
#'
#' @importFrom grid gpar unit grid.text
#' @importFrom stats quantile sd
#' @importFrom dplyr bind_rows
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_vline facet_wrap labs scale_colour_manual scale_fill_manual theme
#'
#' @export plotProfile

plotProfile <-
  function(object,
           style = "heatmap",
           colours = NULL,
           groupColours = NULL,
           limits = NULL,
           title = NULL,
           fontSize = 9,
           baseSize = 12,
           legendPosition = "right",
           ...) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!(style[1] %in% c("heatmap", "lines"))) {
      stop("The 'style' parameter must be either 'heatmap' or 'lines'.", call. = FALSE)
    }
    style <- style[1]

    precomputed <- is.list(object) && all(c("profiles", "regions", "bins") %in% names(object))
    profileData <- if (isTRUE(precomputed)) {object} else {computeProfiles(object = object, ...)}

    profileList <- profileData$profiles
    rowGroups <- factor(profileData$regions$row.group, levels = unique(profileData$regions$row.group))
    binCentres <- profileData$bins
    signalLabel <- profileData$parameters$signal.label
    distance <- profileData$parameters$distance

    if (is.null(title) && !is.null(profileData$parameters$contrast)) {
      title <- profileData$parameters$contrast
    }

    #------------------------#
    # Mean of every group    #
    #------------------------#
    # Mean and standard error of the rows of each group, for every sample or group of samples
    meanTable <- dplyr::bind_rows(lapply(names(profileList), function(profileName) {
      dplyr::bind_rows(lapply(levels(rowGroups), function(groupName) {
        groupMatrix <- profileList[[profileName]][rowGroups == groupName, , drop = FALSE]
        data.frame(profile = profileName,
                   row.group = groupName,
                   position = binCentres,
                   mean = colMeans(groupMatrix),
                   se = apply(groupMatrix, MARGIN = 2, FUN = stats::sd) / sqrt(nrow(groupMatrix)),
                   stringsAsFactors = FALSE)
      }))
    }))
    meanTable$se[!is.finite(meanTable$se)] <- 0
    meanTable$profile <- factor(meanTable$profile, levels = names(profileList))
    meanTable$row.group <- factor(meanTable$row.group, levels = levels(rowGroups))

    #------------------------#
    # Lines                  #
    #------------------------#
    if (style == "lines") {
      lineColours <- .upsetGroupColours(groupLevels = names(profileList), groupColours = groupColours)

      profilePlot <-
        ggplot2::ggplot(data = meanTable, mapping = ggplot2::aes(x = .data$position, y = .data$mean, colour = .data$profile, fill = .data$profile)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey80", linetype = 2) +
        ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = .data$mean - .data$se, ymax = .data$mean + .data$se), colour = NA, alpha = 0.2) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::facet_wrap(~ row.group) +
        ggplot2::scale_colour_manual(values = lineColours, name = NULL) +
        ggplot2::scale_fill_manual(values = lineColours, name = NULL) +
        ggplot2::labs(x = "Distance from the centre (bp)", y = signalLabel, title = title) +
        .regionSetTheme(legendPosition = legendPosition, baseSize = baseSize) +
        ggplot2::theme(panel.spacing = grid::unit(1.5, "lines"))

      return(profilePlot)
    }

    #------------------------#
    # Heatmaps               #
    #------------------------#
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) | !requireNamespace("circlize", quietly = TRUE)) {
      stop("The 'ComplexHeatmap' and 'circlize' packages are needed to draw the profile heatmaps.", call. = FALSE)
    }

    if (is.null(colours)) {
      if (!requireNamespace("viridisLite", quietly = TRUE)) {
        stop("The 'viridisLite' package is needed for the default palette, or pass one through 'colours'.", call. = FALSE)
      }
      colours <- viridisLite::mako(n = 100, direction = -1)
    }

    # The same row order in every panel, strongest regions on top within each group
    rowSignal <- Reduce(f = `+`, x = lapply(profileList, rowMeans)) / length(profileList)
    rowOrder <- order(as.integer(rowGroups), -rowSignal)

    allValues <- unlist(lapply(profileList, as.numeric), use.names = FALSE)
    scaleLimits <- if (is.null(limits)) {c(0, as.numeric(stats::quantile(allValues, probs = 0.99, na.rm = TRUE)))} else {limits}
    if (anyNA(scaleLimits)) {scaleLimits[is.na(scaleLimits)] <- range(allValues, na.rm = TRUE)[is.na(scaleLimits)]}
    if (diff(scaleLimits) <= 0) {scaleLimits <- c(scaleLimits[1], scaleLimits[1] + 1)}

    colourFunction <- circlize::colorRamp2(breaks = seq(scaleLimits[1], scaleLimits[2], length.out = length(colours)), colors = colours)
    # Up and down take the colours they have in the volcano and MA plots
    rowGroupColours <- if (is.null(groupColours) & all(levels(rowGroups) %in% c("up", "down"))) {
      .diffStatusColours(colours = NULL)[levels(rowGroups)]
    } else {
      .upsetGroupColours(groupLevels = levels(rowGroups), groupColours = groupColours)
    }

    # The ends and the centre only, the ends justified inwards so that neighbouring panels do not collide
    endLabels <- paste0(c("-", "+"), format(distance, big.mark = ",", trim = TRUE), " bp")
    axisFontSize <- fontSize - 1

    axisAnnotation <- ComplexHeatmap::HeatmapAnnotation(
      distance = ComplexHeatmap::AnnotationFunction(fun = function(index, k, n) {
        grid::grid.text(endLabels[1], x = grid::unit(0, "npc"), y = grid::unit(1, "npc"), just = c("left", "top"), gp = grid::gpar(fontsize = axisFontSize))
        grid::grid.text("centre", x = grid::unit(0.5, "npc"), y = grid::unit(1, "npc"), just = c("centre", "top"), gp = grid::gpar(fontsize = axisFontSize))
        grid::grid.text(endLabels[2], x = grid::unit(1, "npc"), y = grid::unit(1, "npc"), just = c("right", "top"), gp = grid::gpar(fontsize = axisFontSize))
      },
      which = "column",
      height = grid::unit(4, "mm"),
      var_import = list(endLabels = endLabels, axisFontSize = axisFontSize),
      subsettable = FALSE),
      show_annotation_name = FALSE)

    profileMaximum <- max(meanTable$mean + meanTable$se, na.rm = TRUE)
    rowTitles <- paste0(levels(rowGroups), "\n(", as.numeric(table(rowGroups)), ")")

    heatmapList <- NULL

    for (i in seq_along(profileList)) {
      profileName <- names(profileList)[i]
      groupMeans <- vapply(levels(rowGroups), function(groupName) {colMeans(profileList[[profileName]][rowGroups == groupName, , drop = FALSE])}, numeric(length(binCentres)))
      groupMeans <- matrix(groupMeans, nrow = length(binCentres))

      topAnnotation <- ComplexHeatmap::HeatmapAnnotation(profile = ComplexHeatmap::anno_lines(groupMeans,
                                                                                              gp = grid::gpar(col = rowGroupColours, lwd = 1.5),
                                                                                              ylim = c(0, profileMaximum),
                                                                                              axis = i == 1,
                                                                                              height = grid::unit(2, "cm")),
                                                         show_annotation_name = FALSE)

      panelHeatmap <- ComplexHeatmap::Heatmap(matrix = profileList[[profileName]][rowOrder, , drop = FALSE],
                                              name = paste0("profile_", i),
                                              col = colourFunction,
                                              cluster_rows = FALSE,
                                              cluster_columns = FALSE,
                                              row_split = rowGroups[rowOrder],
                                              row_title = if (i == 1) {rowTitles} else {NULL},
                                              row_title_gp = grid::gpar(fontsize = fontSize, col = rowGroupColours),
                                              show_row_names = FALSE,
                                              show_column_names = FALSE,
                                              bottom_annotation = axisAnnotation,
                                              column_title = profileName,
                                              column_title_gp = grid::gpar(fontsize = fontSize + 1, fontface = "bold"),
                                              top_annotation = topAnnotation,
                                              show_heatmap_legend = i == 1,
                                              heatmap_legend_param = list(title = gsub(" per ", "\nper ", signalLabel),
                                                                          title_gp = grid::gpar(fontsize = fontSize),
                                                                          labels_gp = grid::gpar(fontsize = fontSize - 1)),
                                              border = TRUE)

      heatmapList <- if (is.null(heatmapList)) {panelHeatmap} else {heatmapList + panelHeatmap}
    }

    if (length(profileList) == 1) {heatmapList <- heatmapList + NULL}

    return(heatmapList)
  } # END function




#' @title .profileRegions
#'
#' @description Picks the regions of a profile and the point each of them is aligned on: the differential regions of a result, the region sets of a counts object, or regions given by the user.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param results \code{RegionSetDE.results} object, or \code{NULL}.
#' @param regions Regions given by the user, or \code{NULL}.
#' @param set Character vector with the region sets used, or \code{NULL}.
#' @param direction String, one among \code{"both"}, \code{"up"} and \code{"down"}.
#' @param FDR Numeric value with the adjusted p-value cut-off, or \code{NULL}.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off, or \code{NULL}.
#' @param maxRegions Numeric value with the maximum number of regions per group of rows.
#' @param centre String, either \code{"summit"} or \code{"midpoint"}.
#'
#' @return A \code{GRanges} with one element per row, carrying \code{row.group}, \code{region.key} and \code{centre.position}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment rowData rowRanges assay
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
#'
#' @keywords internal

.profileRegions <-
  function(counts,
           results,
           regions,
           set,
           direction,
           FDR,
           log2FC,
           maxRegions,
           centre) {

    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))
    rowTable$region.key <- paste(rowTable$region.set, rowTable$region.id, sep = "|")
    hasSummit <- "summit" %in% colnames(rowTable)

    #-------------------------------#
    # Regions given by the user     #
    #-------------------------------#
    if (!is.null(regions)) {
      if (methods::is(regions, "RegionSetDE")) {regions <- regions@regions}
      if (methods::is(regions, "GRanges")) {regions <- list(regions = regions)}
      if (methods::is(regions, "GRangesList")) {regions <- as.list(regions)}

      if (!is.list(regions) || is.null(names(regions)) || any(names(regions) == "") ||
          !all(vapply(regions, function(x) {methods::is(x, "GRanges")}, logical(1)))) {
        stop("The 'regions' parameter must be a GRanges, a named list of GRanges, a GRangesList or a RegionSetDE object.", call. = FALSE)
      }

      groupList <- lapply(names(regions), function(groupName) {
        groupRanges <- regions[[groupName]]
        if (length(groupRanges) > maxRegions) {groupRanges <- groupRanges[.thinIndex(n = length(groupRanges), maxPoints = maxRegions)]}

        # A summit column carried by the regions is used like the one of the counts
        summitPosition <- if (centre == "summit" && "summit" %in% colnames(S4Vectors::mcols(groupRanges))) {as.integer(groupRanges$summit)} else {rep(NA_integer_, length(groupRanges))}
        midpoint <- as.integer((BiocGenerics::start(groupRanges) + BiocGenerics::end(groupRanges)) %/% 2L)

        GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(groupRanges),
                               ranges = IRanges::IRanges(start = BiocGenerics::start(groupRanges), end = BiocGenerics::end(groupRanges)),
                               strand = BiocGenerics::strand(groupRanges),
                               row.group = groupName,
                               region.key = paste0(groupName, "|", as.character(GenomeInfoDb::seqnames(groupRanges)), ":", BiocGenerics::start(groupRanges), "-", BiocGenerics::end(groupRanges)),
                               centre.position = ifelse(is.na(summitPosition), midpoint, summitPosition))
      })

      return(.stackProfileRegions(groupList))
    }

    #-------------------------------#
    # Differential regions          #
    #-------------------------------#
    if (!is.null(results)) {
      directionList <- if (direction == "both") {c("up", "down")} else {direction}
      setNames <- if (is.null(set)) {unique(as.character(results@results$region.set))} else {set}

      groupList <- lapply(directionList, function(directionName) {
        # A direction with nothing passing drops out, the error below speaks when no direction is left
        topTable <- try(suppressWarnings(topRegions(results = results, n = maxRegions, set = setNames, sortBy = "FDR",
                                                    direction = directionName, FDR = FDR, log2FC = log2FC)),
                        silent = TRUE)

        if (inherits(topTable, "try-error") || nrow(topTable) == 0) {return(NULL)}

        topTable$region.key <- paste(topTable$region.set, topTable$region.id, sep = "|")
        midpoint <- as.integer((topTable$start + topTable$end) %/% 2L)

        # The summit sits in the counts, one row per region when they were not tiled
        summitPosition <- if (centre == "summit" & hasSummit) {rowTable$summit[match(topTable$region.key, rowTable$region.key)]} else {rep(NA_integer_, nrow(topTable))}
        regionStrand <- as.character(BiocGenerics::strand(SummarizedExperiment::rowRanges(counts)))[match(topTable$region.key, rowTable$region.key)]
        regionStrand[is.na(regionStrand)] <- "*"

        GenomicRanges::GRanges(seqnames = topTable$seqnames,
                               ranges = IRanges::IRanges(start = topTable$start, end = topTable$end),
                               strand = regionStrand,
                               row.group = directionName,
                               region.key = topTable$region.key,
                               centre.position = ifelse(is.na(summitPosition), midpoint, as.integer(summitPosition)))
      })

      groupList <- groupList[!vapply(groupList, is.null, logical(1))]

      if (length(groupList) == 0) {
        stop("No region passes the thresholds in the direction asked for, there is nothing to profile.", call. = FALSE)
      }

      return(.stackProfileRegions(groupList))
    }

    #-------------------------------#
    # Region sets of the counts     #
    #-------------------------------#
    if (identical(counts@counting.level, "tile")) {
      stop("A tiled object has no single centre per region: profile the results of a test, or pass the regions through 'regions'.", call. = FALSE)
    }

    setNames <- if (is.null(set)) {unique(as.character(rowTable$region.set))} else {set}
    absentSets <- setdiff(setNames, rowTable$region.set)
    if (length(absentSets) > 0) {
      stop("The following sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
    }

    rowStrength <- rowMeans(as.matrix(SummarizedExperiment::assay(counts, "counts")))
    rowRangesObject <- SummarizedExperiment::rowRanges(counts)

    groupList <- lapply(setNames, function(setName) {
      setIndex <- which(rowTable$region.set == setName)

      # The strongest regions of the set first, the ones drawn when there are too many
      setIndex <- setIndex[order(rowStrength[setIndex], decreasing = TRUE)]
      setIndex <- setIndex[seq_len(min(length(setIndex), maxRegions))]

      setRanges <- rowRangesObject[setIndex]
      midpoint <- as.integer((BiocGenerics::start(setRanges) + BiocGenerics::end(setRanges)) %/% 2L)
      summitPosition <- if (centre == "summit" & hasSummit) {as.integer(rowTable$summit[setIndex])} else {rep(NA_integer_, length(setIndex))}

      GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(setRanges),
                             ranges = IRanges::IRanges(start = BiocGenerics::start(setRanges), end = BiocGenerics::end(setRanges)),
                             strand = BiocGenerics::strand(setRanges),
                             row.group = setName,
                             region.key = rowTable$region.key[setIndex],
                             centre.position = ifelse(is.na(summitPosition), midpoint, summitPosition))
    })

    return(.stackProfileRegions(groupList))
  } # END function




#' @title .stackProfileRegions
#'
#' @description Stacks the groups of rows of a profile into one \code{GRanges}, dropping their seqinfo so that groups from different sources can be pasted together.
#'
#' @param groupList List of \code{GRanges}, one per group of rows.
#'
#' @return A \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand
#' @importFrom IRanges ranges
#' @importFrom S4Vectors mcols
#'
#' @keywords internal

.stackProfileRegions <-
  function(groupList) {

    groupList <- lapply(groupList, function(groupRanges) {
      GenomicRanges::GRanges(seqnames = as.character(GenomeInfoDb::seqnames(groupRanges)),
                             ranges = IRanges::ranges(groupRanges),
                             strand = as.character(BiocGenerics::strand(groupRanges)),
                             S4Vectors::mcols(groupRanges))
    })

    stackedRanges <- do.call(what = c, args = unname(groupList))

    if (length(stackedRanges) == 0) {
      stop("There is no region to profile.", call. = FALSE)
    }

    return(stackedRanges)
  } # END function




#' @title .profileSignalFiles
#'
#' @description Works out which files the signal of a profile is read from, and of which type they are.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param signal String, one among \code{"auto"}, \code{"bam"} and \code{"bigwig"}.
#' @param signalFiles Character vector with one file per sample, or \code{NULL}.
#'
#' @return A list with the \code{files} and their \code{type}, either \code{"bam"} or \code{"bigwig"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal

.profileSignalFiles <-
  function(counts,
           signal,
           signalFiles) {

    signal <- tolower(as.character(signal[1]))
    if (!(signal %in% c("auto", "bam", "bigwig"))) {
      stop("The 'signal' parameter must be one of 'auto', 'bam', 'bigwig'.", call. = FALSE)
    }

    sampleTable <- as.data.frame(SummarizedExperiment::colData(counts))
    bamFiles <- counts@parameters$countReads$bamFiles
    bigwigFiles <- if ("bigwig" %in% colnames(sampleTable)) {as.character(sampleTable$bigwig)} else {counts@parameters$countBigwig$bigwigFiles}

    #-------------------------------#
    # Files given here              #
    #-------------------------------#
    if (!is.null(signalFiles)) {
      if (!is.character(signalFiles) | length(signalFiles) != ncol(counts)) {
        stop("The 'signalFiles' parameter must hold one file per sample, in the order of the samples.", call. = FALSE)
      }

      fileType <- if (signal != "auto") {signal} else if (all(grepl("\\.bam$", signalFiles, ignore.case = TRUE))) {"bam"} else {"bigwig"}
      files <- signalFiles
    } else {
      fileType <- if (signal == "auto") {if (!is.null(bamFiles)) {"bam"} else {"bigwig"}} else {signal}
      files <- if (fileType == "bam") {bamFiles} else {bigwigFiles}

      if (is.null(files) || all(is.na(files))) {
        stop("No ", if (fileType == "bam") {"BAM"} else {"bigWig"}, " file is recorded in the object, pass them through 'signalFiles'.", call. = FALSE)
      }
    }

    #-------------------------------#
    # Files usable at all           #
    #-------------------------------#
    remoteFiles <- grepl("^[A-Za-z][A-Za-z0-9+.-]*://", files)
    missingFiles <- files[!remoteFiles & (is.na(files) | !file.exists(files))]

    if (length(missingFiles) > 0) {
      stop("The following signal files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    if (fileType == "bam") {
      unindexedFiles <- files[!.hasBamIndex(files)]
      if (length(unindexedFiles) > 0) {
        stop("The following BAM files are not indexed: ", paste(basename(unindexedFiles), collapse = ", "), ".", call. = FALSE)
      }
    }

    return(list(files = files, type = fileType))
  } # END function




#' @title .binnedSignal
#'
#' @description Reads the signal of every file over a set of windows and averages it in bins of fixed width, one matrix per file.
#'
#' @param files Character vector with the BAM or bigWig files, one per sample.
#' @param type String, either \code{"bam"} or \code{"bigwig"}.
#' @param windows \code{GRanges} with the windows, all of the same width, a multiple of \code{binWidth}.
#' @param binWidth Integer value with the width of the bins.
#' @param counts \code{RegionSetDE.counts} object, whose counting parameters are used for BAM files.
#' @param nThreads Number of threads, one file per thread.
#'
#' @return A list with \code{matrices}, one matrix per file with one row per window and one column per bin, and \code{kept}, a logical vector telling which windows lie within their chromosome.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqlengths seqlevels
#' @importFrom BiocGenerics start end width strand
#' @importFrom IRanges IRanges ranges Views viewMeans
#' @importFrom BiocParallel bplapply
#' @importFrom rtracklayer BigWigFile import
#'
#' @keywords internal

.binnedSignal <-
  function(files,
           type,
           windows,
           binWidth,
           counts,
           nThreads = 1) {

    #-------------------------------#
    # Chromosomes of the files      #
    #-------------------------------#
    # The files may name the chromosomes in different styles, those of the first file stand for all of them
    if (type == "bam") {
      chromosomeLengths <- .bamChromosomeMap(bamFiles = files)$lengths
    } else {
      chromosomeLengths <- GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(files[1]))
    }

    # The windows are renamed to the style of the files for the reading only
    fileWindows <- .matchSeqlevels(x = windows, targetSeqlevels = names(chromosomeLengths), fileName = files[1], verbose = FALSE)
    windowChromosomes <- as.character(GenomeInfoDb::seqnames(fileWindows))

    keptWindows <- windowChromosomes %in% names(chromosomeLengths) &
      BiocGenerics::start(fileWindows) >= 1 &
      BiocGenerics::end(fileWindows) <= as.numeric(chromosomeLengths[windowChromosomes])
    keptWindows[is.na(keptWindows)] <- FALSE

    fileWindows <- fileWindows[keptWindows]
    windowChromosomes <- windowChromosomes[keptWindows]

    binNumber <- as.integer(BiocGenerics::width(fileWindows)[1] %/% binWidth)
    onMinus <- as.character(BiocGenerics::strand(fileWindows)) == "-"

    # Every bin of every window, the windows one after the other
    binStarts <- as.vector(t(outer(BiocGenerics::start(fileWindows), (seq_len(binNumber) - 1L) * binWidth, FUN = "+")))
    binChromosomes <- rep(windowChromosomes, each = binNumber)

    #-------------------------------#
    # Counting parameters           #
    #-------------------------------#
    countingParameters <- counts@parameters$countReads

    # The reads left out of the counting stay out of the profiles as well
    discardRegions <- if (type == "bam") {.storedDiscardRegions(counts = counts, targetSeqlevels = names(chromosomeLengths))} else {NULL}
    pairedEnd <- if (is.null(countingParameters$pairedEnd)) {rep(FALSE, length(files))} else {rep_len(countingParameters$pairedEnd, length(files))}
    fragmentLength <- if (is.null(countingParameters$fragmentLength)) {rep(150, length(files))} else {rep_len(countingParameters$fragmentLength, length(files))}

    #-------------------------------#
    # One file per thread           #
    #-------------------------------#
    matrixList <-
      BiocParallel::bplapply(seq_along(files),
                             function(fileIndex) {
                               fileCoverage <- if (type == "bam") {
                                 .windowCoverage(bamFile = files[fileIndex],
                                                 windows = fileWindows,
                                                 chromosomeLengths = chromosomeLengths,
                                                 isPairedEnd = pairedEnd[fileIndex],
                                                 fragmentLength = fragmentLength[fileIndex],
                                                 maxFragmentLength = if (is.null(countingParameters$maxFragmentLength)) {1000} else {countingParameters$maxFragmentLength[1]},
                                                 minMapq = if (is.null(countingParameters$minMapq)) {20} else {countingParameters$minMapq},
                                                 removeDuplicates = if (is.null(countingParameters$removeDuplicates)) {TRUE} else {countingParameters$removeDuplicates},
                                                 discardRegions = discardRegions)$coverage
                               } else {
                                 # A bigWig is asked under its own chromosome names, and answers under them
                                 bigwigFile <- rtracklayer::BigWigFile(files[fileIndex])
                                 fileChromosomes <- .translateChromosomeNames(chromosomeNames = names(chromosomeLengths),
                                                                              targetSeqlevels = GenomeInfoDb::seqlevels(bigwigFile))
                                 names(fileChromosomes) <- names(chromosomeLengths)

                                 readableWindows <- !is.na(fileChromosomes[windowChromosomes])

                                 if (any(readableWindows)) {
                                   bigwigWindows <- GenomicRanges::GRanges(seqnames = unname(fileChromosomes[windowChromosomes[readableWindows]]),
                                                                           ranges = IRanges::ranges(fileWindows[readableWindows]))

                                   bigwigCoverage <- rtracklayer::import(bigwigFile, which = bigwigWindows, as = "RleList")

                                   # Back under the names of the windows, the chromosomes without an equivalent left out
                                   commonNames <- names(fileChromosomes)[match(names(bigwigCoverage), fileChromosomes)]
                                   bigwigCoverage <- bigwigCoverage[!is.na(commonNames)]
                                   names(bigwigCoverage) <- commonNames[!is.na(commonNames)]
                                   bigwigCoverage
                                 } else {
                                   list()
                                 }
                               }

                               binMeans <- numeric(length(binStarts))

                               for (chromosome in unique(windowChromosomes)) {
                                 binIndex <- which(binChromosomes == chromosome)
                                 chromosomeCoverage <- if (chromosome %in% names(fileCoverage)) {fileCoverage[[chromosome]]} else {NULL}
                                 if (is.null(chromosomeCoverage)) {next}

                                 binMeans[binIndex] <- IRanges::viewMeans(IRanges::Views(chromosomeCoverage, IRanges::IRanges(start = binStarts[binIndex], width = binWidth)))
                               }

                               binMeans[!is.finite(binMeans)] <- 0
                               profileMatrix <- matrix(binMeans, ncol = binNumber, byrow = TRUE)

                               # The minus strand is read from right to left, so that every region faces the same way
                               if (any(onMinus)) {profileMatrix[onMinus, ] <- profileMatrix[onMinus, rev(seq_len(binNumber)), drop = FALSE]}

                               return(profileMatrix)
                             },
                             BPPARAM = .makeParallelParam(nThreads = nThreads, tasks = length(files)))

    return(list(matrices = matrixList, kept = keptWindows))
  } # END function




#' @title .profileScaling
#'
#' @description Returns the divisor turning the coverage of every sample into coverage per million fragments of the average library, scaled by the normalisation of the object when there is one.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param useOffsets Logical value indicating whether the scaling factors of the normalisation must be used.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A numeric vector with one divisor per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal

.profileScaling <-
  function(counts,
           useOffsets,
           verbose) {

    sampleTable <- as.data.frame(SummarizedExperiment::colData(counts))
    librarySizes <- as.numeric(sampleTable$library.size)

    if (is.null(sampleTable$library.size) || any(is.na(librarySizes) | librarySizes <= 0)) {
      stop("The object carries no library size to scale the coverage with.", call. = FALSE)
    }

    # A library of average depth is the unit, and one million of its fragments the scale
    depthScale <- mean(librarySizes) / 1e6
    scalingFactors <- sampleTable$scaling.factor

    if (isTRUE(useOffsets) & !is.null(scalingFactors) && !any(is.na(scalingFactors))) {
      return(as.numeric(scalingFactors) * depthScale)
    }

    if (isTRUE(useOffsets) & isTRUE(verbose)) {
      message("No scaling factor is stored in the object, the coverage is scaled by the library sizes alone.")
    }

    return(librarySizes / mean(librarySizes) * depthScale)
  } # END function
