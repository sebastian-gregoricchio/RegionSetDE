# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .provenanceSlots
#'
#' @description Collects the provenance slots of a \code{RegionSetDE} object. Regions arriving as a plain \code{GRangesList} carry no history, so the empty defaults are returned instead.
#'
#' @param regionSet Object passed to the counting functions.
#'
#' @return A named list with the slots shared by the \code{RegionSetDE.provenance} classes.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.provenanceSlots <-
  function(regionSet) {
    if (methods::is(regionSet, "RegionSetDE.provenance")) {
      return(list(blacklist = regionSet@blacklist,
                  whitelist = regionSet@whitelist,
                  genome.assembly = regionSet@genome.assembly,
                  seqlevels.style = regionSet@seqlevels.style,
                  filtering.log = regionSet@filtering.log,
                  parameters = regionSet@parameters))
    }

    return(list(blacklist = NULL,
                whitelist = NULL,
                genome.assembly = NULL,
                seqlevels.style = NA_character_,
                filtering.log = data.frame(step = character(0), region.set = character(0),
                                           n.before = numeric(0), n.after = numeric(0), n.removed = numeric(0)),
                parameters = list()))
  } # END function




#' @title .tileRegionSets
#'
#' @description Cuts each region into adjacent tiles of fixed width, propagating the metadata columns of the parent region to all its tiles.
#'
#' @param regions \code{GRanges} with the flattened region sets.
#' @param tileWidth Numeric value with the width of the tiles, in base pairs.
#' @param partialTiles Logical value: \code{TRUE} keeps the trailing tile even when shorter than \code{tileWidth}, \code{FALSE} discards it. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} with one element per tile and an extra \code{tile.id} metadata column.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges slidingWindows
#' @importFrom S4Vectors mcols mcols<- elementNROWS
#' @importFrom BiocGenerics width
#' @importFrom dplyr filter n_distinct
#' @importFrom rlang .data
#'
#' @keywords internal

.tileRegionSets <-
  function(regions,
           tileWidth,
           partialTiles = TRUE,
           verbose = TRUE) {

    tileWidth <- as.integer(tileWidth[1])
    if (is.na(tileWidth) | tileWidth < 1) {
      stop("The 'tileWidth' parameter must be a positive integer.", call. = FALSE)
    }

    # A step equal to the width returns adjacent, non overlapping tiles
    tileList <- IRanges::slidingWindows(x = regions, width = tileWidth, step = tileWidth)
    tilesPerRegion <- S4Vectors::elementNROWS(tileList)

    tiles <- unlist(tileList, use.names = FALSE)
    S4Vectors::mcols(tiles) <- S4Vectors::mcols(regions)[rep(seq_along(regions), times = tilesPerRegion), , drop = FALSE]
    S4Vectors::mcols(tiles)$tile.id <- unlist(lapply(tilesPerRegion, seq_len), use.names = FALSE)

    # A shorter trailing tile collects less signal than its siblings and distorts any per-tile comparison
    if (isFALSE(partialTiles)) {
      tileTable <- data.frame(region.key = paste(S4Vectors::mcols(tiles)$region.set, S4Vectors::mcols(tiles)$region.id, sep = "|"),
                              tile.width = BiocGenerics::width(tiles),
                              stringsAsFactors = FALSE)

      keptTiles <- dplyr::filter(tileTable, .data$tile.width == tileWidth)
      lostRegions <- dplyr::n_distinct(tileTable$region.key) - dplyr::n_distinct(keptTiles$region.key)

      if (lostRegions > 0 & isTRUE(verbose)) {
        warning(lostRegions, " regions are narrower than 'tileWidth' and have been removed together with their partial tiles.", call. = FALSE)
      }

      tiles <- tiles[BiocGenerics::width(tiles) == tileWidth]
    }

    return(tiles)
  } # END function




#' @title .flattenRegionSets
#'
#' @description Turns a collection of region sets into a single \code{GRanges}, one element per row of the future counts matrix. The set name and the region identifier are stored in the metadata columns, whatever else the regions carried is harmonised across the sets and kept alongside them, and the regions are optionally cut into tiles.
#'
#' @param regionSet \code{RegionSetDE} object, \code{GRangesList} or named list of \code{GRanges}.
#' @param tileWidth Numeric value with the width of the tiles. Default: \code{NULL}, one row per region.
#' @param keepMetadata Logical value to indicate whether the metadata columns of the regions must be carried over. Default: \code{TRUE}.
#' @param regionId String with the name of a metadata column holding the region identifiers. Default: \code{NULL}, the names of the ranges, and their coordinates when they are unnamed.
#' @param partialTiles Logical value indicating whether the trailing shorter tile must be kept. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A named \code{GRanges} with the \code{region.set}, \code{region.id} and \code{tile.id} metadata columns, followed by whatever the regions carried.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols mcols<- DataFrame
#' @importFrom BiocGenerics start end
#' @importFrom methods is as
#'
#' @keywords internal

.flattenRegionSets <-
  function(regionSet,
           tileWidth = NULL,
           partialTiles = TRUE,
           keepMetadata = TRUE,
           regionId = NULL,
           verbose = TRUE) {

    #------------------------------#
    # Uniform the input to a list  #
    #------------------------------#
    if (methods::is(regionSet, "RegionSetDE.counts")) {
      stop("The counts have already been computed for this object.", call. = FALSE)
    } else if (methods::is(regionSet, "RegionSetDE")) {
      regionList <- as.list(regionSet@regions)
    } else if (methods::is(regionSet, "GRangesList")) {
      regionList <- as.list(regionSet)
    } else if (is.list(regionSet) & all(vapply(regionSet, function(x) {methods::is(x, "GRanges")}, logical(1)))) {
      regionList <- regionSet
    } else {
      stop("The 'regionSet' parameter must be a RegionSetDE object, a GRangesList, a named list of GRanges or a GRanges.", call. = FALSE)
    }

    if (is.null(names(regionList)) | any(is.na(names(regionList))) | any(names(regionList) == "")) {
      stop("All the region sets must be named.", call. = FALSE)
    }

    #--------------------------------#
    # Stack the sets, one after the  #
    # other, harmonising what they   #
    # carry                          #
    #--------------------------------#
    # Sets loaded from different files rarely share the same columns, and the stacking needs them to
    annotationPlan <- .metadataPlan(regionList = regionList, keepMetadata = keepMetadata, verbose = verbose)

    flatList <-
      lapply(names(regionList),
             function(setName) {
               gr <- regionList[[setName]]

               identifiers <- .regionIdentifiers(gr = gr, regionId = regionId, setName = setName)
               annotationTable <- .alignMetadata(gr = gr, plan = annotationPlan)

               S4Vectors::mcols(gr) <- cbind(S4Vectors::DataFrame(region.set = setName, region.id = identifiers),
                                             annotationTable)
               names(gr) <- NULL
               return(gr)
             })

    allRegions <- do.call(what = c, args = flatList)

    #-----------------#
    # Optional tiling #
    #-----------------#
    if (!is.null(tileWidth)) {
      allRegions <- .tileRegionSets(regions = allRegions, tileWidth = tileWidth, partialTiles = partialTiles, verbose = verbose)
      names(allRegions) <- paste0(S4Vectors::mcols(allRegions)$region.set, "|", S4Vectors::mcols(allRegions)$region.id, "|tile", S4Vectors::mcols(allRegions)$tile.id)
    } else {
      S4Vectors::mcols(allRegions)$tile.id <- NA_integer_
      names(allRegions) <- paste0(S4Vectors::mcols(allRegions)$region.set, "|", S4Vectors::mcols(allRegions)$region.id)
    }

    # The three identifiers first, the annotation of the regions behind them
    identifierColumns <- c("region.set", "region.id", "tile.id")
    S4Vectors::mcols(allRegions) <- S4Vectors::mcols(allRegions)[, c(identifierColumns,
                                                                     setdiff(colnames(S4Vectors::mcols(allRegions)), identifierColumns)),
                                                                 drop = FALSE]

    # The row names index the object from here on, a collision would silently mix up two regions
    if (any(duplicated(names(allRegions)))) {
      stop("Some regions share the same identifier within a set, run 'loadRegions' with 'removeDuplicatedRegions = TRUE'.", call. = FALSE)
    }

    return(allRegions)
  } # END function




#' @title .buildSampleTable
#'
#' @description Assembles the sample table used as \code{colData}, deriving the sample names from the file names when they are not provided and attaching the user metadata.
#'
#' @param files Character vector with the paths of the signal files.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, derived from the file names.
#' @param sampleMetadata Data.frame with the sample annotation. Default: \code{NULL}.
#' @param fileColumn String with the name of the column storing the file paths. Default: \code{"file"}.
#' @param extensionPattern Regular expression removed from the file names to build the sample names. Default: \code{"\\.[^.]*$"}.
#'
#' @return A data.frame with one row per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr left_join bind_cols
#'
#' @keywords internal

.buildSampleTable <-
  function(files,
           sampleNames = NULL,
           sampleMetadata = NULL,
           fileColumn = "file",
           extensionPattern = "\\.[^.]*$") {

    if (is.null(sampleNames)) {
      sampleNames <- sub(extensionPattern, "", basename(files), ignore.case = TRUE)
    }

    if (length(sampleNames) != length(files)) {
      stop("The 'sampleNames' parameter must have the same length as the number of files.", call. = FALSE)
    }

    # The names become the column names of the assay, duplicates would make the samples unaddressable
    if (any(duplicated(sampleNames))) {
      stop("The sample names must be unique.", call. = FALSE)
    }

    sampleTable <- data.frame(sample = sampleNames, file = files, stringsAsFactors = FALSE)
    colnames(sampleTable)[2] <- fileColumn

    if (!is.null(sampleMetadata)) {
      sampleMetadata <- as.data.frame(sampleMetadata, stringsAsFactors = FALSE)

      if ("sample" %in% colnames(sampleMetadata)) {
        # Joining on the sample name tolerates a metadata table in a different order
        absentSamples <- setdiff(sampleTable$sample, sampleMetadata$sample)
        if (length(absentSamples) > 0) {
          stop("The following samples are missing from 'sampleMetadata': ", paste(absentSamples, collapse = ", "), ".", call. = FALSE)
        }
        sampleTable <- dplyr::left_join(sampleTable, sampleMetadata, by = "sample")
      } else {
        if (nrow(sampleMetadata) != length(files)) {
          stop("Without a 'sample' column, 'sampleMetadata' must have one row per file, in the same order.", call. = FALSE)
        }
        sampleTable <- dplyr::bind_cols(sampleTable, sampleMetadata)
      }
    }

    return(sampleTable)
  } # END function




#' @title .makeParallelParam
#'
#' @description Builds the \code{BiocParallel} back end matching the number of requested threads and the operating system.
#'
#' @param nThreads Number of threads. Default: \code{1}.
#' @param tasks Number of tasks the work is split into, see \code{\link[BiocParallel]{MulticoreParam}}. Setting it to the number of jobs hands the jobs out one at a time, as the threads become free. Default: \code{0}, one task per thread.
#' @param progressBar Logical value to indicate whether the back end must draw its progress bar, which advances by one step for every task that comes back. Default: \code{FALSE}.
#'
#' @return A \code{BiocParallelParam} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#'
#' @keywords internal

.makeParallelParam <-
  function(nThreads = 1,
           tasks = 0L,
           progressBar = FALSE) {
    nThreads <- as.integer(nThreads[1])
    progressBar <- isTRUE(progressBar)

    if (is.na(nThreads) | nThreads < 1) {
      stop("The 'nThreads' parameter must be a positive integer.", call. = FALSE)
    }

    if (nThreads == 1) {
      return(BiocParallel::SerialParam(progressbar = progressBar))
    }

    # Windows has no forking, sockets give the same result at a higher start-up cost
    if (.Platform$OS.type == "windows") {
      return(BiocParallel::SnowParam(workers = nThreads, tasks = as.integer(tasks), progressbar = progressBar))
    }

    # Forked workers already see the options of the session, sending them along with every task only costs time
    return(BiocParallel::MulticoreParam(workers = nThreads, tasks = as.integer(tasks), exportglobals = FALSE, progressbar = progressBar))
  } # END function




#' @title .elapsedTime
#'
#' @description Writes the time gone by since a starting point in the unit that reads best: seconds, minutes or hours.
#'
#' @param startTime Value returned by \code{Sys.time()} when the work started.
#'
#' @return A string such as \code{"42 s"}, \code{"3.5 min"} or \code{"1.2 h"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.elapsedTime <-
  function(startTime) {

    elapsedSeconds <- as.numeric(difftime(Sys.time(), startTime, units = "secs"))

    if (elapsedSeconds < 60) {
      return(paste(round(elapsedSeconds), "s"))
    }

    if (elapsedSeconds < 3600) {
      return(paste(round(elapsedSeconds / 60, 1), "min"))
    }

    return(paste(round(elapsedSeconds / 3600, 1), "h"))
  } # END function





#' @title .registeredChromosomeNames
#'
#' @description Converts chromosome names to the naming style of another set of names through the tables of \code{GenomeInfoDb}, which cover the chromosomes of the species it lists and their mitochondrion.
#'
#' @param chromosomeNames Character vector with the chromosome names to convert.
#' @param targetSeqlevels Character vector with chromosome names written in the style to convert to.
#'
#' @return A character vector as long as \code{chromosomeNames}, with \code{NA} for the names the tables do not hold, scaffolds and custom contigs for instance, or when the style of \code{targetSeqlevels} is not recognised.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevelsStyle mapSeqlevels
#'
#' @keywords internal

.registeredChromosomeNames <-
  function(chromosomeNames,
           targetSeqlevels) {

    unknownNames <- rep(NA_character_, length(chromosomeNames))

    registeredNames <-
      tryCatch({
        mappedNames <- GenomeInfoDb::mapSeqlevels(seqlevels = chromosomeNames,
                                                  style = GenomeInfoDb::seqlevelsStyle(targetSeqlevels)[1])
        if (is.matrix(mappedNames)) {mappedNames <- mappedNames[1, ]}
        as.character(mappedNames)
      },
      error = function(e) {unknownNames},
      warning = function(w) {unknownNames})

    if (length(registeredNames) != length(chromosomeNames)) {
      return(unknownNames)
    }

    return(registeredNames)
  } # END function




#' @title .translateChromosomeNames
#'
#' @description Finds, for every chromosome name, the name the same chromosome goes under in another set of names. A name already there is kept. For the others the conversions of \code{GenomeInfoDb} are tried first, then the \code{chr} prefix is added or dropped, and the mitochondrion is looked for under the four names it is given (\code{chrM}, \code{MT}, \code{chrMT}, \code{M}). Each name is settled on its own, so a set mixing two styles is translated as well as a set written in one.
#'
#' @param chromosomeNames Character vector with the chromosome names to translate.
#' @param targetSeqlevels Character vector with the chromosome names to translate into.
#'
#' @return A character vector as long as \code{chromosomeNames}, with the matching name of \code{targetSeqlevels} and \code{NA} where there is none.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.translateChromosomeNames <-
  function(chromosomeNames,
           targetSeqlevels) {

    chromosomeNames <- as.character(chromosomeNames)
    translatedNames <- as.character(ifelse(chromosomeNames %in% targetSeqlevels, chromosomeNames, NA_character_))
    missingIndex <- which(is.na(translatedNames))

    if (length(missingIndex) == 0 | length(targetSeqlevels) == 0) {
      return(translatedNames)
    }

    missingNames <- chromosomeNames[missingIndex]

    #-------------------------------#
    # Candidates for every name     #
    #-------------------------------#
    # GenomeInfoDb knows the registered styles, the mitochondrion of the species it lists included
    registeredNames <- .registeredChromosomeNames(chromosomeNames = missingNames, targetSeqlevels = targetSeqlevels)

    # Scaffolds and custom contigs are in no table, for them the prefix is the only thing to try
    prefixNames <- ifelse(grepl("^chr", missingNames), sub("^chr", "", missingNames), paste0("chr", missingNames))

    # The mitochondrion is the one chromosome the styles disagree on beyond the prefix
    mitochondrialNames <- c("chrM", "MT", "chrMT", "M")
    candidateList <- c(list(registeredNames, prefixNames),
                       lapply(mitochondrialNames, function(mitochondrialName) {
                         ifelse(missingNames %in% mitochondrialNames, mitochondrialName, NA_character_)
                       }))

    #-------------------------------#
    # First candidate that exists   #
    #-------------------------------#
    foundNames <- rep(NA_character_, length(missingNames))

    for (candidateNames in candidateList) {
      usable <- is.na(foundNames) & !is.na(candidateNames) & candidateNames %in% targetSeqlevels
      foundNames[usable] <- candidateNames[usable]
    }

    translatedNames[missingIndex] <- foundNames

    return(translatedNames)
  } # END function




#' @title .matchChromosomeNames
#'
#' @description Renames a plain vector of chromosome names into the naming style of a signal file, the way \code{.matchSeqlevels} does for a set of ranges. What it is for is \code{excludeChromosomes}: a name that matches nothing is not an error, it simply excludes nothing, and since the chromosomes left out decide the library sizes the silence would be paid for by the normalisation.
#'
#' @param chromosomeNames Character vector with the chromosome names given by the user.
#' @param targetSeqlevels Character vector with the chromosome names of the files.
#' @param argumentName String with the name of the argument, used in the warning. Default: \code{"excludeChromosomes"}.
#'
#' @return The names, converted where a conversion was needed.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils head
#'
#' @keywords internal

.matchChromosomeNames <-
  function(chromosomeNames,
           targetSeqlevels,
           argumentName = "excludeChromosomes") {

    if (length(chromosomeNames) == 0) {
      return(chromosomeNames)
    }

    # Only the names that miss are touched, so a list mixing the two styles still works
    translatedNames <- .translateChromosomeNames(chromosomeNames = chromosomeNames, targetSeqlevels = targetSeqlevels)
    chromosomeNames <- ifelse(is.na(translatedNames), as.character(chromosomeNames), translatedNames)

    stillMissing <- setdiff(chromosomeNames, targetSeqlevels)

    if (length(stillMissing) > 0) {
      warning("The following names in '", argumentName, "' match no chromosome of the files and exclude nothing: ",
              paste(stillMissing, collapse = ", "), ". The files use names such as ",
              paste(utils::head(targetSeqlevels, 3), collapse = ", "), ".", call. = FALSE)
    }

    return(chromosomeNames)
  } # END function




#' @title .matchSeqlevels
#'
#' @description Renames the chromosomes of a set of ranges so that they follow the names of a signal file, or of another set of ranges. The chromosomes are settled one by one: those already written as in the target are left alone, and the others take the name the target gives them. Those the target does not have under any name, scaffolds for the most part, follow its style, so that the object does not come back half renamed. Only the copy used for the counting is renamed, so the object returned to the user keeps the style of the regions it was built from.
#'
#' @param x \code{GRanges}, or any object accepting \code{seqlevels}, to be renamed.
#' @param targetSeqlevels Character vector with the chromosome names to align to, usually read from the header of a signal file.
#' @param fileName String with the file path, used in the messages. Default: \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The input object with the renamed chromosomes.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges ranges
#' @importFrom BiocGenerics strand
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
#' @importFrom utils head
#'
#' @keywords internal

.matchSeqlevels <-
  function(x,
           targetSeqlevels,
           fileName = NULL,
           verbose = TRUE) {

    currentSeqlevels <- GenomeInfoDb::seqlevels(x)

    if (length(currentSeqlevels) == 0 | length(targetSeqlevels) == 0) {
      return(x)
    }

    translatedSeqlevels <- .translateChromosomeNames(chromosomeNames = currentSeqlevels, targetSeqlevels = targetSeqlevels)

    # Not a single chromosome in common under any name: the two sides are not on the same assembly
    if (all(is.na(translatedSeqlevels))) {
      stop("The chromosome names of ", if (is.null(fileName)) {"the signal file"} else {c("'", basename(fileName), "'")},
           " cannot be reconciled with the ones of the regions: the file uses ", paste(utils::head(targetSeqlevels, 3), collapse = ", "),
           " while the regions use ", paste(utils::head(currentSeqlevels, 3), collapse = ", "), ".", call. = FALSE)
    }

    # The chromosomes the target does not have do no harm. They follow its style all the same, GenomeInfoDb
    # first and the prefix by hand for the rest, the mitochondrion being renamed before the prefix is taken off it
    unmatchedSeqlevels <- is.na(translatedSeqlevels)
    styledSeqlevels <- currentSeqlevels

    if (any(unmatchedSeqlevels)) {
      registeredSeqlevels <- .registeredChromosomeNames(chromosomeNames = currentSeqlevels, targetSeqlevels = targetSeqlevels)

      manualSeqlevels <-
        if (any(grepl("^chr", targetSeqlevels))) {
          ifelse(grepl("^chr", currentSeqlevels), currentSeqlevels, paste0("chr", sub("^MT$", "M", currentSeqlevels)))
        } else {
          ifelse(grepl("^chr", currentSeqlevels),
                 ifelse(currentSeqlevels == "chrM", "MT", sub("^chr", "", currentSeqlevels)),
                 currentSeqlevels)
        }

      styledSeqlevels <- ifelse(is.na(registeredSeqlevels), manualSeqlevels, registeredSeqlevels)
    }

    renamedSeqlevels <- ifelse(unmatchedSeqlevels, styledSeqlevels, translatedSeqlevels)

    # A styled name must not land on a chromosome that is there already, it is then left as it was
    clashingSeqlevels <- unmatchedSeqlevels &
      (renamedSeqlevels %in% renamedSeqlevels[!unmatchedSeqlevels] | duplicated(ifelse(unmatchedSeqlevels, renamedSeqlevels, NA_character_)) |
         (renamedSeqlevels != currentSeqlevels & renamedSeqlevels %in% currentSeqlevels))
    renamedSeqlevels[clashingSeqlevels] <- currentSeqlevels[clashingSeqlevels]

    changedSeqlevels <- renamedSeqlevels != currentSeqlevels

    if (!any(changedSeqlevels)) {
      return(x)
    }

    if (any(duplicated(renamedSeqlevels))) {
      # Ranges written half as chr1 and half as 1 end up on one chromosome, and the two levels have to merge
      if (!methods::is(x, "GRanges")) {
        stop("The conversion of the chromosome names produced duplicated entries, harmonise the styles before counting.", call. = FALSE)
      }

      rangeSeqnames <- renamedSeqlevels[match(as.character(GenomeInfoDb::seqnames(x)), currentSeqlevels)]

      mergedRanges <- GenomicRanges::GRanges(seqnames = factor(rangeSeqnames, levels = unique(renamedSeqlevels)),
                                             ranges = IRanges::ranges(x),
                                             strand = BiocGenerics::strand(x))
      S4Vectors::mcols(mergedRanges) <- S4Vectors::mcols(x)
      x <- mergedRanges
    } else {
      GenomeInfoDb::seqlevels(x) <- renamedSeqlevels
    }

    if (isTRUE(verbose)) {
      message("The chromosome names have been converted from ", paste(utils::head(currentSeqlevels[changedSeqlevels], 2), collapse = ", "),
              " to ", paste(utils::head(renamedSeqlevels[changedSeqlevels], 2), collapse = ", "), " to match the signal files.")
    }

    return(x)
  } # END function




#' @title .asRegionSet
#'
#' @description Lets the counting functions take a single \code{GRanges} where they expect region sets. The ranges go through \code{\link{loadRegions}} and come back as a \code{RegionSetDE} object with one set, so that everything downstream finds the object it expects. Anything else is returned as it is.
#'
#' @param regionSet Object passed to the counting functions.
#' @param setName String with the name given to the set. Default: \code{"regions"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE} object when \code{regionSet} is a \code{GRanges}, \code{regionSet} itself otherwise.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb genome
#' @importFrom methods is
#' @importFrom stats setNames
#'
#' @keywords internal

.asRegionSet <-
  function(regionSet,
           setName = "regions",
           verbose = TRUE) {

    if (!methods::is(regionSet, "GRanges")) {
      return(regionSet)
    }

    if (length(regionSet) == 0) {
      stop("The GRanges given as 'regionSet' holds no region.", call. = FALSE)
    }

    # The assembly written in the ranges is the only thing known about them, it follows the object
    regionAssembly <- unique(as.character(GenomeInfoDb::genome(regionSet)))
    regionAssembly <- regionAssembly[!is.na(regionAssembly) & regionAssembly != ""]

    # The order and the chromosome names of the ranges are kept, so that the rows can be matched back to them
    loadedSet <- loadRegions(regions = stats::setNames(list(regionSet), setName),
                             sortRegions = FALSE,
                             seqlevelsStyle = NULL,
                             genomeAssembly = if (length(regionAssembly) == 1) {regionAssembly} else {NULL},
                             verbose = FALSE)

    if (isTRUE(verbose)) {
      removedRegions <- length(regionSet) - length(loadedSet@regions[[1]])

      message("The GRanges given as 'regionSet' has been loaded with loadRegions() as a single set named '", setName, "' (",
              length(loadedSet@regions[[1]]), " regions",
              if (removedRegions > 0) {c(", ", removedRegions, " duplicated ones removed")} else {""}, ").")
    }

    return(loadedSet)
  } # END function




#' @title .newCountsObject
#'
#' @description Assembles a \code{RegionSetDE.counts} object from a matrix of values, the regions and the sample table, carrying over the provenance of the region sets.
#'
#' @param countMatrix Matrix with the values, one row per region and one column per sample.
#' @param regions \code{GRanges} used to compute the counts, in the same order as the matrix rows.
#' @param sampleTable Data.frame with the sample annotation, in the same order as the matrix columns.
#' @param provenance List returned by \code{.provenanceSlots}.
#' @param countingLevel String indicating whether the rows are regions (\code{"region"}) or tiles of a region (\code{"tile"}). Default: \code{"region"}.
#' @param newParameters List with the arguments of the calling function, appended to the stored parameters.
#' @param metadataList List stored in the \code{metadata} of the object. Default: \code{list()}.
#'
#' @return A \code{RegionSetDE.counts} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#'
#' @keywords internal

.newCountsObject <-
  function(countMatrix,
           regions,
           sampleTable,
           provenance,
           countingLevel = "region",
           newParameters = list(),
           metadataList = list()) {

    rownames(countMatrix) <- names(regions)
    colnames(countMatrix) <- sampleTable$sample

    countsExperiment <-
      SummarizedExperiment::SummarizedExperiment(assays = list(counts = countMatrix),
                                                 rowRanges = regions,
                                                 colData = S4Vectors::DataFrame(sampleTable, row.names = sampleTable$sample),
                                                 metadata = metadataList)

    return(methods::new("RegionSetDE.counts",
                        countsExperiment,
                        counting.level = countingLevel,
                        blacklist = provenance$blacklist,
                        whitelist = provenance$whitelist,
                        genome.assembly = provenance$genome.assembly,
                        seqlevels.style = provenance$seqlevels.style,
                        filtering.log = provenance$filtering.log,
                        parameters = c(provenance$parameters, newParameters)))
  } # END function




#' @title .metadataPlan
#'
#' @description Works out which metadata columns of a collection of region sets can be stacked into one table, and in which type, so that sets loaded from different files can be put one after the other.
#'
#' @param regionList Named list of \code{GRanges}.
#' @param keepMetadata Logical value to indicate whether the metadata must be carried over at all.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with, for every column kept, the name it takes in the output and the type it is stored as.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#'
#' @keywords internal

.metadataPlan <-
  function(regionList,
           keepMetadata = TRUE,
           verbose = TRUE) {

    if (isFALSE(keepMetadata)) {
      return(list())
    }

    identifierColumns <- c("region.set", "region.id", "tile.id")

    #-------------------------------#
    # What every set carries        #
    #-------------------------------#
    columnClasses <- list()
    droppedColumns <- character(0)

    for (setName in names(regionList)) {
      setMetadata <- S4Vectors::mcols(regionList[[setName]])

      if (is.null(setMetadata) | ncol(setMetadata) == 0) {
        next
      }

      for (columnName in colnames(setMetadata)) {
        # A list column cannot be stacked into a flat table and would break the binding rather than the row
        if (!is.atomic(setMetadata[[columnName]])) {
          droppedColumns <- unique(c(droppedColumns, columnName))
          next
        }

        columnClasses[[columnName]] <- unique(c(columnClasses[[columnName]], class(setMetadata[[columnName]])[1]))
      }
    }

    if (length(droppedColumns) > 0 & isTRUE(verbose)) {
      message("The following region columns are not atomic and have been left out: ",
              paste(droppedColumns, collapse = ", "), ".")
    }

    if (length(columnClasses) == 0) {
      return(list())
    }

    #-------------------------------#
    # One name and one type each    #
    #-------------------------------#
    plan <- list()
    renamedColumns <- character(0)
    coercedColumns <- character(0)

    for (columnName in names(columnClasses)) {
      outputName <- columnName

      # The package writes these three itself, so a column of the same name has to step aside
      if (columnName %in% identifierColumns) {
        outputName <- paste0(columnName, ".original")
        renamedColumns <- c(renamedColumns, columnName)
      }

      # One set calling a column a number and another calling it a word leaves character as the only common ground
      outputClass <- columnClasses[[columnName]]
      if (length(outputClass) > 1) {
        outputClass <- "character"
        coercedColumns <- c(coercedColumns, columnName)
      }

      plan[[columnName]] <- list(name = outputName, class = outputClass)
    }

    if (isTRUE(verbose)) {
      if (length(renamedColumns) > 0) {
        message("The following region columns share a name with an identifier and carry the suffix '.original': ",
                paste(renamedColumns, collapse = ", "), ".")
      }
      if (length(coercedColumns) > 0) {
        message("The following region columns hold different types across the sets and have been read as text: ",
                paste(coercedColumns, collapse = ", "), ".")
      }
    }

    return(plan)
  } # END function




#' @title .alignMetadata
#'
#' @description Builds the metadata table of one region set in the shape every set has to share, filling with \code{NA} the columns that set does not carry.
#'
#' @param gr \code{GRanges} of one region set.
#' @param plan List returned by \code{.metadataPlan}.
#'
#' @return A \code{DataFrame} with one row per region and one column per entry of the plan.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols DataFrame
#'
#' @keywords internal

.alignMetadata <-
  function(gr,
           plan) {

    if (length(plan) == 0) {
      return(S4Vectors::DataFrame(row.names = seq_along(gr)))
    }

    setMetadata <- S4Vectors::mcols(gr)

    columnList <-
      lapply(names(plan),
             function(columnName) {
               targetClass <- plan[[columnName]]$class

               # A set that never had this column contributes a column of the right type full of NA
               if (is.null(setMetadata) || !(columnName %in% colnames(setMetadata))) {
                 return(as.vector(rep(NA, length(gr)), mode = targetClass))
               }

               columnValues <- setMetadata[[columnName]]

               if (!identical(class(columnValues)[1], targetClass)) {
                 columnValues <- as.vector(columnValues, mode = targetClass)
               }

               return(columnValues)
             })

    names(columnList) <- vapply(plan, function(x) {x$name}, character(1))

    return(S4Vectors::DataFrame(columnList, row.names = NULL, check.names = FALSE))
  } # END function




#' @title .regionIdentifiers
#'
#' @description Returns the identifier of every region of a set: the names of the ranges, a metadata column when one is asked for, or the coordinates when there is nothing else.
#'
#' @param gr \code{GRanges} of one region set.
#' @param regionId String with the name of a metadata column, or \code{NULL}.
#' @param setName String with the name of the set, used in the messages.
#'
#' @return A character vector with one identifier per region.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#'
#' @keywords internal

.regionIdentifiers <-
  function(gr,
           regionId = NULL,
           setName = "") {

    coordinateId <- paste0(as.character(GenomeInfoDb::seqnames(gr)), ":",
                           BiocGenerics::start(gr), "-", BiocGenerics::end(gr))

    #-------------------------------#
    # A column asked for by name    #
    #-------------------------------#
    if (!is.null(regionId)) {
      setMetadata <- S4Vectors::mcols(gr)

      if (is.null(setMetadata) || !(regionId %in% colnames(setMetadata))) {
        stop("The column '", regionId, "' is absent from the set '", setName, "'.", call. = FALSE)
      }

      identifiers <- as.character(setMetadata[[regionId]])

      # The identifier indexes the object from here on, and two regions sharing one would be mixed up
      if (any(is.na(identifiers)) | any(duplicated(identifiers))) {
        stop("The column '", regionId, "' does not hold a unique identifier for every region of the set '",
             setName, "'.", call. = FALSE)
      }

      return(identifiers)
    }

    if (is.null(names(gr))) {
      return(coordinateId)
    }

    return(names(gr))
  } # END function
