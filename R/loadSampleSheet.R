# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title loadSampleSheet
#'
#' @description Reads a table describing the libraries of an experiment, one row per sample, and returns it with the file paths resolved and checked. Beside the sample name and its signal file, a row can point to the peaks called on that sample and to the input library it was sequenced against. Every other column is kept as sample annotation, with no restriction on names or number.
#'
#' @param sampleSheet String with the path to a comma or tab separated file, or a data.frame.
#' @param columns Named character vector telling which column of the table holds each standard field, for instance \code{c(sample = "library", bam = "file")}. The standard fields are \code{sample}, \code{bam}, \code{bigwig}, \code{peaks}, \code{input} and \code{input.id}. A field not given here is looked for under its own name and then under the one used by DiffBind (\code{SampleID}, \code{bamReads}, \code{Peaks}, \code{bamControl}, \code{ControlID}), ignoring the case. Default: \code{NULL}.
#' @param basePath String with the directory the relative paths of the table refer to. Default: \code{NULL}, the directory holding the table when it is read from a file, the working directory otherwise.
#' @param checkFiles Logical value to indicate whether the existence of the files, and the index of the BAM files, must be checked. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A data.frame with one row per sample. The standard fields found in the table come first, under their standard names and with absolute paths, followed by every other column as it was. Samples without an input or without peaks carry \code{NA} there. When the table has inputs, \code{input.id} names each distinct input library, taken from the table when it has such a column and from the file names otherwise.
#'
#' @details Only \code{sample} and one of \code{bam} and \code{bigwig} are required. The other fields are optional, and so are their values: a sample can come without peaks or without an input, and one input can serve several samples. \code{\link{makeGreylist}} reads each distinct input file once, whatever the number of samples pointing to it.
#'
#' A DiffBind sample sheet is read as it is. \code{SampleID}, \code{bamReads}, \code{Peaks}, \code{bamControl} and \code{ControlID} take the standard names, while \code{Tissue}, \code{Factor}, \code{Condition}, \code{Treatment}, \code{Replicate} and any other column stay as annotation, available to the design of \code{\link{fitRegions}} under their own names.
#'
#' Relative paths are resolved against \code{basePath} and returned as absolute paths, so the table stays valid when the working directory changes. Addresses starting with a protocol, such as \code{https://}, are left untouched and are not checked, since \code{\link{countBigwig}} can read remote bigWig files.
#'
#' The table goes straight into the counting, as the source of the files, of the sample names and of the annotation: \code{countReads(regions, bamFiles = sampleSheet$bam, sampleNames = sampleSheet$sample, sampleMetadata = sampleSheet)}.
#'
#' @examples
#' # Two samples sharing an input and a third one without any, paths relative to a project folder
#' sheetTable <- data.frame(sample = c("treated_1", "control_1", "control_2"),
#'                          bam = c("reads/treated_1.bam", "reads/control_1.bam", "reads/control_2.bam"),
#'                          peaks = c("peaks/treated_1.narrowPeak", "peaks/control_1.narrowPeak", NA),
#'                          input = c("reads/input_A.bam", "reads/input_A.bam", NA),
#'                          condition = c("treated", "control", "control"))
#'
#' sampleSheet <- loadSampleSheet(sheetTable, basePath = "/data/project", checkFiles = FALSE)
#' sampleSheet
#'
#' # A DiffBind sample sheet is read as it is
#' diffbindTable <- data.frame(SampleID = c("MCF7_1", "MCF7_2"),
#'                             Condition = c("resistant", "responsive"),
#'                             Replicate = c(1, 1),
#'                             bamReads = c("reads/MCF7_1.bam", "reads/MCF7_2.bam"),
#'                             ControlID = c("MCF7_input", "MCF7_input"),
#'                             bamControl = c("reads/MCF7_input.bam", "reads/MCF7_input.bam"))
#'
#' loadSampleSheet(diffbindTable, basePath = "/data/project", checkFiles = FALSE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{makeGreylist}}, \code{\link{countReads}}, \code{\link{countBigwig}}
#'
#' @importFrom dplyr relocate all_of
#'
#' @export loadSampleSheet

loadSampleSheet <-
  function(sampleSheet,
           columns = NULL,
           basePath = NULL,
           checkFiles = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!is.logical(checkFiles) | length(checkFiles) != 1 | anyNA(checkFiles)) {
      stop("The 'checkFiles' parameter must be a single logical value.", call. = FALSE)
    }

    if (!is.null(basePath) & (!is.character(basePath) | length(basePath) != 1)) {
      stop("The 'basePath' parameter must be a single string.", call. = FALSE)
    }

    #-------------------------------#
    # Read the table                #
    #-------------------------------#
    if (is.character(sampleSheet) & length(sampleSheet) == 1) {
      if (!file.exists(sampleSheet)) {
        stop("The sample sheet '", sampleSheet, "' does not exist.", call. = FALSE)
      }
      sheetTable <- .readSheetFile(path = sampleSheet)

      # Paths written in the sheet are usually meant from where the sheet sits
      if (is.null(basePath)) {basePath <- dirname(normalizePath(sampleSheet))}
    } else if (is.data.frame(sampleSheet)) {
      sheetTable <- as.data.frame(sampleSheet, stringsAsFactors = FALSE)
      if (is.null(basePath)) {basePath <- getwd()}
    } else {
      stop("The 'sampleSheet' parameter must be the path to a table or a data.frame.", call. = FALSE)
    }

    if (nrow(sheetTable) == 0) {
      stop("The sample sheet has no rows.", call. = FALSE)
    }

    basePath <- normalizePath(path.expand(basePath), mustWork = FALSE)

    #-------------------------------#
    # Standard fields               #
    #-------------------------------#
    fieldTable <- .matchSheetColumns(columnNames = colnames(sheetTable), columns = columns)

    for (i in seq_len(nrow(fieldTable))) {
      colnames(sheetTable)[colnames(sheetTable) == fieldTable$column[i]] <- fieldTable$field[i]
    }

    renamedTable <- fieldTable[fieldTable$column != fieldTable$field, , drop = FALSE]
    if (isTRUE(verbose) & nrow(renamedTable) > 0) {
      message("Columns read as standard fields: ", paste(renamedTable$column, "as", renamedTable$field, collapse = ", "), ".")
    }

    if (!("sample" %in% fieldTable$field)) {
      stop("No sample column was found, call it 'sample' or point to it through 'columns'.", call. = FALSE)
    }

    # A sheet of peaks alone builds a consensus, it just has nothing to count afterwards
    if (!any(c("bam", "bigwig", "peaks") %in% fieldTable$field)) {
      stop("The sample sheet needs a 'bam', a 'bigwig' or a 'peaks' column, or one pointed to through 'columns'.", call. = FALSE)
    }

    if (isTRUE(verbose) & !any(c("bam", "bigwig") %in% fieldTable$field)) {
      message("The sheet carries peaks and no signal file, which is enough for a consensus and not for the counting.")
    }

    # Empty cells and blanks mean the same as a missing value in any of these fields
    for (field in fieldTable$field) {
      fieldValues <- trimws(as.character(sheetTable[[field]]))
      fieldValues[!is.na(fieldValues) & fieldValues == ""] <- NA_character_
      sheetTable[[field]] <- fieldValues
    }

    #-------------------------------#
    # Samples                       #
    #-------------------------------#
    if (anyNA(sheetTable$sample)) {
      stop("The following rows have no sample name: ", paste(which(is.na(sheetTable$sample)), collapse = ", "), ".", call. = FALSE)
    }

    if (anyDuplicated(sheetTable$sample) > 0) {
      stop("The sample names must be unique, repeated: ",
           paste(unique(sheetTable$sample[duplicated(sheetTable$sample)]), collapse = ", "), ".", call. = FALSE)
    }

    signalFields <- intersect(c("bam", "bigwig"), fieldTable$field)
    noSignal <- Reduce(`&`, lapply(signalFields, function(field) {is.na(sheetTable[[field]])}))

    if (any(noSignal)) {
      stop("The following samples have no signal file: ", paste(sheetTable$sample[noSignal], collapse = ", "), ".", call. = FALSE)
    }

    #-------------------------------#
    # Paths                         #
    #-------------------------------#
    pathFields <- intersect(c("bam", "bigwig", "peaks", "input"), fieldTable$field)

    for (field in pathFields) {
      sheetTable[[field]] <- .resolveSheetPaths(paths = sheetTable[[field]], basePath = basePath)
    }

    # The same library listed twice is almost always a copy-paste slip in the sheet
    for (field in signalFields) {
      repeatedFiles <- unique(sheetTable[[field]][!is.na(sheetTable[[field]]) & duplicated(sheetTable[[field]])])
      if (length(repeatedFiles) > 0) {
        warning("The following ", field, " files are listed for more than one sample: ",
                paste(basename(repeatedFiles), collapse = ", "), ".", call. = FALSE)
      }
    }

    #-------------------------------#
    # Inputs                        #
    #-------------------------------#
    if ("input" %in% colnames(sheetTable)) {
      sheetTable$input.id <- .inputIdentifiers(inputPaths = sheetTable$input,
                                               inputIds = if ("input.id" %in% colnames(sheetTable)) {sheetTable$input.id} else {NULL})
    } else if ("input.id" %in% colnames(sheetTable)) {
      warning("The sample sheet names inputs through 'input.id' without any 'input' file, the column is kept as annotation only.", call. = FALSE)
    }

    #-------------------------------#
    # Files on disk                 #
    #-------------------------------#
    if (isTRUE(checkFiles)) {
      .checkSheetFiles(sheetTable = sheetTable, pathFields = pathFields)
    }

    #-------------------------------#
    # Order of the columns          #
    #-------------------------------#
    standardColumns <- intersect(c("sample", "bam", "bigwig", "peaks", "input", "input.id"), colnames(sheetTable))
    sheetTable <- dplyr::relocate(sheetTable, dplyr::all_of(standardColumns))
    rownames(sheetTable) <- NULL

    if (isTRUE(verbose)) {
      sampleNumber <- nrow(sheetTable)
      peakSummary <- if ("peaks" %in% colnames(sheetTable)) {paste(sum(!is.na(sheetTable$peaks)), "of", sampleNumber)} else {"none"}
      inputSummary <- if ("input" %in% colnames(sheetTable)) {
        paste(length(unique(sheetTable$input[!is.na(sheetTable$input)])), "distinct, none for",
              sum(is.na(sheetTable$input)), "of", sampleNumber, "samples")
      } else {
        "none"
      }

      message("Sample sheet with ", sampleNumber, " samples. Peaks: ", peakSummary, ". Inputs: ", inputSummary, ".")
    }

    return(sheetTable)
  } # END function




#' @title .readSheetFile
#'
#' @description Reads a sample sheet from a comma or tab separated file, the separator being taken from the extension or, failing that, from the first line.
#'
#' @param path String with the path to the file.
#'
#' @return A data.frame, with the column names as written in the file.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils read.table
#'
#' @keywords internal

.readSheetFile <-
  function(path) {

    extension <- tolower(sub("^.*\\.", "", basename(sub("\\.gz$", "", path, ignore.case = TRUE))))

    if (extension %in% c("xls", "xlsx")) {
      stop("Spreadsheet files are not read directly: save the sheet as .csv or .tsv, or read it into a data.frame first.", call. = FALSE)
    }

    # A csv is split on commas and a tsv on tabs, anything else on whichever the header line uses
    firstLine <- readLines(path, n = 1, warn = FALSE)

    separator <- if (extension == "csv") {
      ","
    } else if (extension %in% c("tsv", "tab")) {
      "\t"
    } else if (grepl("\t", firstLine)) {
      "\t"
    } else {
      ","
    }

    sheetTable <- utils::read.table(path,
                                    header = TRUE,
                                    sep = separator,
                                    quote = "\"",
                                    comment.char = "",
                                    na.strings = c("", "NA"),
                                    strip.white = TRUE,
                                    check.names = FALSE,
                                    stringsAsFactors = FALSE)

    return(sheetTable)
  } # END function




#' @title .matchSheetColumns
#'
#' @description Finds the columns of a sample sheet holding the standard fields, from the names given by the user first, then from the standard names and from the DiffBind ones, ignoring the case.
#'
#' @param columnNames Character vector with the column names of the table.
#' @param columns Named character vector given by the user, or \code{NULL}.
#'
#' @return A data.frame with the \code{field} and the \code{column} holding it, one row per field found.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.matchSheetColumns <-
  function(columnNames,
           columns = NULL) {

    aliasList <- list(sample = c("sample", "sampleid"),
                      bam = c("bam", "bamreads"),
                      bigwig = c("bigwig"),
                      peaks = c("peaks"),
                      input = c("input", "bamcontrol"),
                      input.id = c("input.id", "controlid"))

    #-------------------------------#
    # Columns named by the user     #
    #-------------------------------#
    if (!is.null(columns)) {
      if (!is.character(columns) | is.null(names(columns)) | any(names(columns) == "")) {
        stop("The 'columns' parameter must be a named character vector, e.g. c(sample = \"library\").", call. = FALSE)
      }

      unknownFields <- setdiff(names(columns), names(aliasList))
      if (length(unknownFields) > 0) {
        stop("The following names of 'columns' are not standard fields: ", paste(unknownFields, collapse = ", "),
             ". The fields are: ", paste(names(aliasList), collapse = ", "), ".", call. = FALSE)
      }

      absentColumns <- setdiff(columns, columnNames)
      if (length(absentColumns) > 0) {
        stop("The following columns given in 'columns' are absent from the sample sheet: ", paste(absentColumns, collapse = ", "), ".", call. = FALSE)
      }
    }

    #-------------------------------#
    # Every other field by alias    #
    #-------------------------------#
    fieldTable <- data.frame(field = character(0), column = character(0), stringsAsFactors = FALSE)

    for (field in names(aliasList)) {
      if (field %in% names(columns)) {
        matchedColumn <- columns[[field]]
      } else {
        matchedColumn <- columnNames[tolower(columnNames) %in% aliasList[[field]]]
        if (length(matchedColumn) == 0) {next}

        # Two candidate columns, e.g. 'sample' and 'SampleID', would have to be chosen between blindly
        if (length(matchedColumn) > 1) {
          stop("Several columns could hold the '", field, "' field: ", paste(matchedColumn, collapse = ", "),
               ". Point to the right one through 'columns'.", call. = FALSE)
        }
      }

      fieldTable <- rbind(fieldTable, data.frame(field = field, column = matchedColumn, stringsAsFactors = FALSE))
    }

    if (anyDuplicated(fieldTable$column) > 0) {
      stop("The same column cannot hold two fields: ", paste(unique(fieldTable$column[duplicated(fieldTable$column)]), collapse = ", "), ".", call. = FALSE)
    }

    # A renamed column would clash with another column already carrying the standard name
    clashingFields <- fieldTable$field[fieldTable$field != fieldTable$column & fieldTable$field %in% columnNames]
    if (length(clashingFields) > 0) {
      stop("The sample sheet already has columns named ", paste(clashingFields, collapse = ", "),
           ", which would be duplicated by the renaming. Rename them first.", call. = FALSE)
    }

    return(fieldTable)
  } # END function




#' @title .resolveSheetPaths
#'
#' @description Turns the relative paths of a sample sheet into absolute ones, leaving alone the missing values, the absolute paths and the addresses starting with a protocol.
#'
#' @param paths Character vector with the paths.
#' @param basePath String with the directory the relative paths refer to.
#'
#' @return A character vector with the resolved paths.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveSheetPaths <-
  function(paths,
           basePath) {

    isRemote <- !is.na(paths) & grepl("^[A-Za-z][A-Za-z0-9+.-]*://", paths)
    isAbsolute <- !is.na(paths) & grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", paths)
    isRelative <- !is.na(paths) & !isRemote & !isAbsolute

    paths[isAbsolute] <- path.expand(paths[isAbsolute])
    paths[isRelative] <- file.path(basePath, paths[isRelative])

    return(paths)
  } # END function




#' @title .inputIdentifiers
#'
#' @description Names every distinct input library of a sample sheet, checking that names and files correspond one to one when the names come with the table.
#'
#' @param inputPaths Character vector with the input file of every sample, \code{NA} for the samples without one.
#' @param inputIds Character vector with the names given in the table, or \code{NULL}.
#'
#' @return A character vector with the name of the input of every sample, \code{NA} where there is no input.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr distinct filter group_by summarise n_distinct
#' @importFrom rlang .data
#'
#' @keywords internal

.inputIdentifiers <-
  function(inputPaths,
           inputIds = NULL) {

    hasInput <- !is.na(inputPaths)

    #-------------------------------#
    # Names taken from the files    #
    #-------------------------------#
    if (is.null(inputIds)) {
      distinctPaths <- unique(inputPaths[hasInput])

      # Two inputs stored under the same file name in different folders still need two names
      distinctIds <- make.unique(sub("\\.bam$", "", basename(distinctPaths), ignore.case = TRUE), sep = "_")

      return(distinctIds[match(inputPaths, distinctPaths)])
    }

    #-------------------------------#
    # Names given with the table    #
    #-------------------------------#
    inputIds[!hasInput] <- NA_character_

    if (anyNA(inputIds[hasInput])) {
      stop("Some samples have an input file but no 'input.id': fill the column, or drop it to name the inputs after their files.", call. = FALSE)
    }

    pairTable <- dplyr::distinct(data.frame(path = inputPaths[hasInput], id = inputIds[hasInput], stringsAsFactors = FALSE))

    # One name for two files, or two names for one file, would make the inputs impossible to tell apart
    idsPerPath <- dplyr::filter(dplyr::summarise(dplyr::group_by(pairTable, .data$path), n = dplyr::n_distinct(.data$id), .groups = "drop"), .data$n > 1)
    pathsPerId <- dplyr::filter(dplyr::summarise(dplyr::group_by(pairTable, .data$id), n = dplyr::n_distinct(.data$path), .groups = "drop"), .data$n > 1)

    if (nrow(idsPerPath) > 0 | nrow(pathsPerId) > 0) {
      stop("Each input must have one name and each name one input file. Check: ",
           paste(c(basename(idsPerPath$path), pathsPerId$id), collapse = ", "), ".", call. = FALSE)
    }

    return(inputIds)
  } # END function




#' @title .checkSheetFiles
#'
#' @description Checks that the files of a sample sheet exist, and that the BAM files are indexed, reporting every problem at once.
#'
#' @param sheetTable Data.frame with the standard fields already resolved.
#' @param pathFields Character vector with the fields holding paths.
#'
#' @return Nothing, it stops when a file is missing or a BAM file has no index.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.checkSheetFiles <-
  function(sheetTable,
           pathFields) {

    problemList <- character(0)

    for (field in pathFields) {
      fieldPaths <- unique(sheetTable[[field]][!is.na(sheetTable[[field]])])

      # Remote addresses cannot be checked from here and are left to the function reading them
      fieldPaths <- fieldPaths[!grepl("^[A-Za-z][A-Za-z0-9+.-]*://", fieldPaths)]

      missingPaths <- fieldPaths[!file.exists(fieldPaths)]
      if (length(missingPaths) > 0) {
        problemList <- c(problemList, paste("missing", field, "files:", paste(missingPaths, collapse = ", ")))
      }

      if (field %in% c("bam", "input")) {
        presentPaths <- setdiff(fieldPaths, missingPaths)
        unindexedPaths <- presentPaths[!.hasBamIndex(presentPaths)]
        if (length(unindexedPaths) > 0) {
          problemList <- c(problemList, paste("unindexed", field, "files:", paste(unindexedPaths, collapse = ", ")))
        }
      }
    }

    if (length(problemList) > 0) {
      stop("The sample sheet points to files that cannot be used, set checkFiles = FALSE to skip this check.\n  ",
           paste(problemList, collapse = "\n  "), call. = FALSE)
    }

    return(invisible(TRUE))
  } # END function




#' @title .hasBamIndex
#'
#' @description Tells whether BAM files come with an index, under any of the names \code{.bamIndexPath} looks for.
#'
#' @param bamFiles Character vector with the paths of the BAM files.
#'
#' @return A logical vector, one value per file.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.hasBamIndex <-
  function(bamFiles) {

    return(!is.na(vapply(bamFiles, .bamIndexPath, character(1), USE.NAMES = FALSE)))
  } # END function




#' @title .sheetCountingInput
#'
#' @description Takes the signal files, the sample names and the annotation of the counting from a sample sheet, or from the sample sheet a consensus was built from when no file is given at all.
#'
#' @param regionSet Object handed to the counting function.
#' @param sampleSheet Data.frame returned by \code{loadSampleSheet}, the path to a sample sheet, or \code{NULL}.
#' @param files Character vector with the files given directly, or \code{NULL}.
#' @param sampleNames Character vector with the sample names given directly, or \code{NULL}.
#' @param sampleMetadata Data.frame with the annotation given directly, or \code{NULL}.
#' @param fileField String with the column of the sheet holding the files, either \code{"bam"} or \code{"bigwig"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with the \code{files}, the \code{sampleNames} and the \code{sampleMetadata} to count with.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.sheetCountingInput <-
  function(regionSet,
           sampleSheet,
           files,
           sampleNames,
           sampleMetadata,
           fileField,
           verbose = TRUE) {

    # Files given directly keep working as they always did
    if (is.null(sampleSheet) & !is.null(files)) {
      return(list(files = files, sampleNames = sampleNames, sampleMetadata = sampleMetadata))
    }

    if (!is.null(sampleSheet) & !is.null(files)) {
      stop("Give the ", fileField, " files either directly or through 'sampleSheet', not both.", call. = FALSE)
    }

    #-------------------------------#
    # Which sheet                   #
    #-------------------------------#
    # Regions built from peaks remember the sheet they came from
    if (is.null(sampleSheet)) {
      if (methods::is(regionSet, "RegionSetDE") && length(regionSet@consensus) > 0 && !is.null(regionSet@consensus$sheet)) {
        sampleSheet <- regionSet@consensus$sheet

        if (isTRUE(verbose)) {
          message("The files and the annotation are taken from the sample sheet the consensus was built from.")
        }
      } else {
        stop("No ", fileField, " file was given: pass them directly, or through 'sampleSheet'.", call. = FALSE)
      }
    }

    if (is.character(sampleSheet) & length(sampleSheet) == 1) {
      sampleSheet <- loadSampleSheet(sampleSheet, verbose = verbose)
    }

    if (!is.data.frame(sampleSheet) || !all(c("sample", fileField) %in% colnames(sampleSheet))) {
      stop("The sample sheet needs the 'sample' and '", fileField, "' columns.", call. = FALSE)
    }

    if (!is.null(sampleNames) | !is.null(sampleMetadata)) {
      stop("With a sample sheet the names and the annotation come from the sheet, leave 'sampleNames' and 'sampleMetadata' empty.", call. = FALSE)
    }

    withoutFile <- sampleSheet$sample[is.na(sampleSheet[[fileField]])]
    if (length(withoutFile) > 0) {
      stop("The following samples of the sheet have no ", fileField, " file: ", paste(withoutFile, collapse = ", "), ".", call. = FALSE)
    }

    # The counting writes the path under its own column, a second copy would only duplicate it
    sheetMetadata <- as.data.frame(sampleSheet, stringsAsFactors = FALSE)
    sheetMetadata <- sheetMetadata[, setdiff(colnames(sheetMetadata), fileField), drop = FALSE]

    return(list(files = as.character(sampleSheet[[fileField]]),
                sampleNames = as.character(sampleSheet$sample),
                sampleMetadata = sheetMetadata))
  } # END function
