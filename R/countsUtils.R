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
        warning(paste0(lostRegions, " regions are narrower than 'tileWidth' and have been removed together with their partial tiles."), call. = FALSE)
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
      stop("The 'regionSet' parameter must be a RegionSetDE object, a GRangesList or a named list of GRanges.", call. = FALSE)
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
          stop(paste0("The following samples are missing from 'sampleMetadata': ", paste(absentSamples, collapse = ", "), "."), call. = FALSE)
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
#'
#' @return A \code{BiocParallelParam} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#'
#' @keywords internal

.makeParallelParam <-
  function(nThreads = 1) {
    nThreads <- as.integer(nThreads[1])

    if (is.na(nThreads) | nThreads < 1) {
      stop("The 'nThreads' parameter must be a positive integer.", call. = FALSE)
    }

    if (nThreads == 1) {
      return(BiocParallel::SerialParam())
    }

    # Windows has no forking, sockets give the same result at a higher start-up cost
    if (.Platform$OS.type == "windows") {
      return(BiocParallel::SnowParam(workers = nThreads))
    }

    return(BiocParallel::MulticoreParam(workers = nThreads))
  } # END function




#' @title .matchSeqlevels
#'
#' @description Renames the chromosomes of a set of ranges so that they follow the naming style of a signal file, leaving them untouched when the two already agree. Only the copy used for the counting is renamed, so the object returned to the user keeps the style of the regions it was built from.
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
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlevelsStyle mapSeqlevels
#' @importFrom utils head
#'
#' @keywords internal

.matchSeqlevels <-
  function(x,
           targetSeqlevels,
           fileName = NULL,
           verbose = TRUE) {

    currentSeqlevels <- GenomeInfoDb::seqlevels(x)

    # Sharing at least one chromosome is enough, the extra scaffolds of either side do no harm
    if (length(intersect(currentSeqlevels, targetSeqlevels)) > 0) {
      return(x)
    }

    # GenomeInfoDb resolves the registered styles, the prefix is added or dropped by hand for everything else
    renamedSeqlevels <-
      tryCatch({
        as.character(GenomeInfoDb::mapSeqlevels(seqlevels = currentSeqlevels,
                                                style = GenomeInfoDb::seqlevelsStyle(targetSeqlevels)[1]))
      },
      error = function(e) {rep(NA_character_, length(currentSeqlevels))},
      warning = function(w) {rep(NA_character_, length(currentSeqlevels))})

    manualSeqlevels <-
      if (any(grepl("^chr", targetSeqlevels))) {
        ifelse(grepl("^chr", currentSeqlevels), currentSeqlevels, paste0("chr", sub("^MT$", "M", currentSeqlevels)))
      } else {
        ifelse(grepl("^chr", currentSeqlevels), sub("^chrM$", "MT", sub("^chr", "", currentSeqlevels)), currentSeqlevels)
      }

    renamedSeqlevels[is.na(renamedSeqlevels)] <- manualSeqlevels[is.na(renamedSeqlevels)]

    if (any(duplicated(renamedSeqlevels))) {
      stop("The conversion of the chromosome names produced duplicated entries, harmonise the styles before counting.", call. = FALSE)
    }

    if (length(intersect(renamedSeqlevels, targetSeqlevels)) == 0) {
      stop(paste0("The chromosome names of ", ifelse(is.null(fileName), "the signal file", paste0("'", basename(fileName), "'")),
                  " cannot be reconciled with the ones of the regions: the file uses ", paste(utils::head(targetSeqlevels, 3), collapse = ", "),
                  " while the regions use ", paste(utils::head(currentSeqlevels, 3), collapse = ", "), "."), call. = FALSE)
    }

    GenomeInfoDb::seqlevels(x) <- renamedSeqlevels

    if (isTRUE(verbose)) {
      message(paste0("The chromosome names have been converted from ", paste(utils::head(currentSeqlevels, 2), collapse = ", "),
                     " to ", paste(utils::head(renamedSeqlevels, 2), collapse = ", "), " to match the signal files."))
    }

    return(x)
  } # END function




#' @title .newCountsObject
#'
#' @description Assembles a \code{RegionSetDE.counts} object from a matrix of values, the regions and the sample table, carrying over the provenance of the region sets.
#'
#' @param countMatrix Matrix with the values, one row per region and one column per sample.
#' @param regions \code{GRanges} used to compute the counts, in the same order as the matrix rows.
#' @param sampleTable Data.frame with the sample annotation, in the same order as the matrix columns.
#' @param provenance List returned by \code{.provenanceSlots}.
#' @param countingLevel String indicating whether the rows are regions or sets. Default: \code{"region"}.
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
