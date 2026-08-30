#' @title exportResults
#'
#' @description Writes a result to disk: the table as a compressed TSV, the coordinates as a BED that a genome browser will colour by direction, and every parameter the analysis was run with as a flat file next to them.
#'
#' @param results \code{RegionSetDE.results}, \code{RegionSetDE.setResults}, or either of the two list classes holding several contrasts.
#' @param path String with the directory the files are written to, created when it does not exist. Default: \code{"."}.
#' @param prefix String prepended to every file name. Default: \code{NULL}, the name of the contrast.
#' @param contrast String with the name of the contrast to write, or its position, when \code{results} holds several of them. Default: \code{NULL}, every one of them.
#' @param set Character vector with the names of the region sets to write. Default: \code{NULL}, all of them.
#' @param onlyChanging Logical value to indicate whether the files must hold only the regions labelled \code{"up"} or \code{"down"}. Default: \code{FALSE}.
#' @param splitByDirection Logical value to indicate whether a separate BED must be written for each direction, on top of the combined one. Default: \code{FALSE}.
#' @param bedScore String with the column mapped onto the BED score, one of \code{"FDR"}, \code{"log2FC"} and \code{"none"}. Default: \code{"FDR"}.
#' @param colourByStatus Logical value to indicate whether the BED must carry an item colour per direction, which makes it nine columns rather than six. Default: \code{TRUE}.
#' @param writeTiles Logical value to indicate whether the tile level table must be written as well, when the counts were tiled. Default: \code{FALSE}.
#' @param provenance Logical value to indicate whether the parameters of the analysis must be written alongside. Default: \code{TRUE}.
#' @param compress Logical value to indicate whether the tables must be gzipped. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return Invisibly, a character vector with the paths written.
#'
#' @details The provenance file is the part worth keeping. It flattens everything the object carries in its \code{parameters} slot into one parameter per line: the counting level, the normalisation method and its factors, the filter and its threshold, the engine, the design, the contrast, the multiple testing correction, and on a fit with no replicates the dispersion and where it came from. Six months later that file is the difference between a result that can be reproduced and one that can only be repeated.
#'
#' The BED is written with zero-based starts, as the format requires, so its coordinates are one lower than the ones in the TSV. With \code{colourByStatus} it carries nine columns and an item colour per direction, which is what makes a browser show the increasing and decreasing regions apart without a second file. The score maps \code{-log10(FDR)} onto the range the format allows and saturates at the top of it, so a score of 1000 means "at least that", not "exactly that".
#'
#' A set level result has no coordinates of its own and only the table is written.
#'
#' @examples
#' \dontrun{
#' exportResults(res, path = "results/H3K27ac")
#'
#' # One set of files per contrast
#' exportResults(resList, path = "results", onlyChanging = TRUE, splitByDirection = TRUE)
#'
#' exportResults(setRes, path = "results", prefix = "sets")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{topRegions}}, \code{\link{resultRanges}}
#'
#' @importFrom utils write.table
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export exportResults

exportResults <-
  function(results,
           path = ".",
           prefix = NULL,
           contrast = NULL,
           set = NULL,
           onlyChanging = FALSE,
           splitByDirection = FALSE,
           bedScore = "FDR",
           colourByStatus = TRUE,
           writeTiles = FALSE,
           provenance = TRUE,
           compress = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!(bedScore %in% c("FDR", "log2FC", "none"))) {
      stop("The 'bedScore' parameter must be one of 'FDR', 'log2FC', 'none'.", call. = FALSE)
    }

    #-------------------------------#
    # Several contrasts at once     #
    #-------------------------------#
    isList <- methods::is(results, "RegionSetDE.resultsList") | methods::is(results, "RegionSetDE.setResultsList")

    if (isTRUE(isList) & is.null(contrast)) {
      writtenPaths <-
        unlist(lapply(names(results@results),
                      function(contrastName) {
                        return(exportResults(results = results@results[[contrastName]],
                                             path = path,
                                             prefix = if (is.null(prefix)) {contrastName} else {paste(prefix, contrastName, sep = "_")},
                                             set = set, onlyChanging = onlyChanging, splitByDirection = splitByDirection,
                                             bedScore = bedScore, colourByStatus = colourByStatus, writeTiles = writeTiles,
                                             provenance = provenance, compress = compress, verbose = verbose))
                      }))
      return(invisible(writtenPaths))
    }

    results <- .pickResults(results = results, contrast = contrast)

    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }

    # A contrast name carries operators and spaces, neither of which belongs in a file name
    filePrefix <- if (is.null(prefix)) {gsub("[^A-Za-z0-9._-]+", "_", results@contrast)} else {prefix}
    tableSuffix <- if (isTRUE(compress)) {".tsv.gz"} else {".tsv"}
    writtenPaths <- character(0)

    #-------------------------------#
    # A set level result            #
    #-------------------------------#
    if (methods::is(results, "RegionSetDE.setResults")) {
      setPath <- file.path(path, paste0(filePrefix, "_sets", tableSuffix))
      .writeTable(table = results@results, path = setPath, compress = compress)
      writtenPaths <- c(writtenPaths, setPath)

    } else {
      #-------------------------------#
      # Region level tables           #
      #-------------------------------#
      resultTable <- results@results

      if (!is.null(set)) {
        absentSets <- setdiff(set, unique(resultTable$region.set))
        if (length(absentSets) > 0) {
          stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
        }
        resultTable <- dplyr::filter(resultTable, .data$region.set %in% set)
      }

      if (isTRUE(onlyChanging)) {
        resultTable <- dplyr::filter(resultTable, .data$diff.status != "null")
      }

      if (nrow(resultTable) == 0) {
        stop("The selection leaves no region to write.", call. = FALSE)
      }

      tablePath <- file.path(path, paste0(filePrefix, "_regions", tableSuffix))
      .writeTable(table = resultTable, path = tablePath, compress = compress)
      writtenPaths <- c(writtenPaths, tablePath)

      #-------------------------------#
      # Coordinates                   #
      #-------------------------------#
      bedPath <- file.path(path, paste0(filePrefix, "_regions.bed"))
      .writeBed(resultTable = resultTable, path = bedPath, bedScore = bedScore, colourByStatus = colourByStatus)
      writtenPaths <- c(writtenPaths, bedPath)

      if (isTRUE(splitByDirection)) {
        for (directionName in c("up", "down")) {
          directionTable <- dplyr::filter(resultTable, .data$diff.status == directionName)

          if (nrow(directionTable) > 0) {
            directionPath <- file.path(path, paste0(filePrefix, "_regions_", directionName, ".bed"))
            .writeBed(resultTable = directionTable, path = directionPath,
                      bedScore = bedScore, colourByStatus = colourByStatus)
            writtenPaths <- c(writtenPaths, directionPath)
          }
        }
      }

      if (isTRUE(writeTiles) & nrow(results@tiles) > 0) {
        tilePath <- file.path(path, paste0(filePrefix, "_tiles", tableSuffix))
        .writeTable(table = results@tiles, path = tilePath, compress = compress)
        writtenPaths <- c(writtenPaths, tilePath)
      }
    }

    #-------------------------------#
    # Provenance                    #
    #-------------------------------#
    if (isTRUE(provenance)) {
      parameterTable <- .flattenParameters(x = c(list(contrast = results@contrast,
                                                      engine = results@engine,
                                                      genome.assembly = results@genome.assembly,
                                                      seqlevels.style = results@seqlevels.style,
                                                      thresholds = results@thresholds),
                                                 results@parameters))

      parameterPath <- file.path(path, paste0(filePrefix, "_parameters.tsv"))
      .writeTable(table = parameterTable, path = parameterPath, compress = FALSE)
      writtenPaths <- c(writtenPaths, parameterPath)
    }

    if (isTRUE(verbose)) {
      message(paste0("Written to ", normalizePath(path), ":\n  ", paste(basename(writtenPaths), collapse = "\n  ")))
    }

    return(invisible(writtenPaths))
  } # END function




#' @title .writeTable
#'
#' @description Writes a table as a tab separated file, gzipped when asked for.
#'
#' @param table Data.frame to write.
#' @param path String with the destination.
#' @param compress Logical value indicating whether the file must be gzipped.
#'
#' @return Invisibly the path.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils write.table
#'
#' @keywords internal

.writeTable <-
  function(table,
           path,
           compress = TRUE) {

    fileConnection <- if (isTRUE(compress)) {gzfile(path, open = "wt")} else {file(path, open = "wt")}
    on.exit(close(fileConnection), add = TRUE)

    utils::write.table(x = as.data.frame(table), file = fileConnection,
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    return(invisible(path))
  } # END function




#' @title .writeBed
#'
#' @description Writes the coordinates of a result table as a BED file, coloured by the direction of the change when asked for.
#'
#' @param resultTable Data.frame with the \code{seqnames}, \code{start} and \code{end} columns.
#' @param path String with the destination.
#' @param bedScore String with the column mapped onto the score.
#' @param colourByStatus Logical value indicating whether an item colour must be written.
#'
#' @return Invisibly the path.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils write.table
#'
#' @keywords internal

.writeBed <-
  function(resultTable,
           path,
           bedScore = "FDR",
           colourByStatus = TRUE) {

    #-------------------------------#
    # Score, on the scale BED wants #
    #-------------------------------#
    # The format takes integers up to 1000, so the score saturates rather than running off the top
    scoreValues <- switch(bedScore,
                          "FDR" = round(-log10(pmax(resultTable$FDR, 1e-300)) * 100),
                          "log2FC" = round(abs(resultTable$log2FC) * 100),
                          "none" = rep(0, nrow(resultTable)))

    scoreValues <- pmin(pmax(scoreValues, 0), 1000)
    scoreValues[!is.finite(scoreValues)] <- 0

    bedTable <- data.frame(chrom = resultTable$seqnames,
                           chromStart = resultTable$start - 1,
                           chromEnd = resultTable$end,
                           name = paste(resultTable$region.set, resultTable$region.id, sep = "|"),
                           score = scoreValues,
                           strand = ".",
                           stringsAsFactors = FALSE)

    #-------------------------------#
    # Colour per direction          #
    #-------------------------------#
    if (isTRUE(colourByStatus) & "diff.status" %in% colnames(resultTable)) {
      statusColours <- c("down" = "33,102,172", "null" = "180,180,180", "up" = "178,24,43")

      bedTable$thickStart <- bedTable$chromStart
      bedTable$thickEnd <- bedTable$chromEnd
      bedTable$itemRgb <- statusColours[as.character(resultTable$diff.status)]
      bedTable$itemRgb[is.na(bedTable$itemRgb)] <- statusColours[["null"]]
    }

    utils::write.table(x = bedTable, file = path, sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = FALSE)

    return(invisible(path))
  } # END function




#' @title .flattenParameters
#'
#' @description Turns the nested list of parameters an object carries into one line per value, so that the record of an analysis can be read, grepped and compared against another.
#'
#' @param x List, or any value held inside one.
#' @param prefix String with the path to the current value, built as the recursion descends.
#' @param maxLength Numeric value with the number of characters a value is truncated to. Default: \code{200}.
#'
#' @return A data.frame with the \code{parameter} and \code{value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.flattenParameters <-
  function(x,
           prefix = "",
           maxLength = 200) {

    if (is.null(x)) {
      return(data.frame(parameter = prefix, value = "NULL", stringsAsFactors = FALSE))
    }

    #-------------------------------#
    # Descend into the lists        #
    #-------------------------------#
    if (is.list(x) & !is.data.frame(x)) {
      if (length(x) == 0) {
        return(data.frame(parameter = prefix, value = "empty", stringsAsFactors = FALSE))
      }

      elementNames <- names(x)
      if (is.null(elementNames)) {
        elementNames <- as.character(seq_along(x))
      }

      return(do.call(what = rbind,
                     args = lapply(seq_along(x),
                                   function(i) {
                                     return(.flattenParameters(x = x[[i]],
                                                               prefix = paste(c(prefix, elementNames[i])[nchar(c(prefix, elementNames[i])) > 0],
                                                                              collapse = "."),
                                                               maxLength = maxLength))
                                   })))
    }

    # A matrix or a table is recorded by its shape, writing it out would drown the file it lives in
    if (is.matrix(x) | is.data.frame(x)) {
      return(data.frame(parameter = prefix,
                        value = paste0("<", class(x)[1], " ", nrow(x), "x", ncol(x), ": ",
                                       paste(colnames(x), collapse = ", "), ">"),
                        stringsAsFactors = FALSE))
    }

    flatValue <- paste(as.character(x), collapse = ", ")
    if (nchar(flatValue) > maxLength) {
      flatValue <- paste0(substr(flatValue, 1, maxLength), " [...]")
    }

    return(data.frame(parameter = prefix, value = flatValue, stringsAsFactors = FALSE))
  } # END function
