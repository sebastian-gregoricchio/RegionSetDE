#' @title loadRegions
#'
#' @description Imports a collection of genomic regions provided as BED-like files and/or GRanges objects and returns them as a single \code{RegionSetDE} object. Files, ranges and data.frames can be mixed within the same list.
#'
#' @param regions List in which each element is either a string indicating the path to a BED-like file (bed/narrowPeak/broadPeak, gzipped files allowed), a \code{GRanges} object, or a data.frame containing chromosome/start/end columns. A single path or a single \code{GRanges} is accepted as well.
#' @param regionNames Character vector with the names to assign to each element of the list. Default: \code{NULL}, meaning that the names of the input list are used and, when missing, the file base names.
#' @param keepMetadata Logical value to indicate whether the metadata columns must be kept. Default: \code{TRUE}.
#' @param sortRegions Logical value to indicate whether the regions must be sorted by coordinate. Default: \code{TRUE}.
#' @param reduceRegions Logical value to indicate whether the overlapping regions within the same element must be collapsed (\code{IRanges::reduce}). Notice that metadata are lost upon reduction. Default: \code{FALSE}.
#' @param removeDuplicatedRegions Logical value to indicate whether the regions with identical coordinates within the same element must be collapsed to a single entry. Only the metadata of the first occurrence are kept. Ignored when \code{reduceRegions = TRUE}. Default: \code{TRUE}.
#' @param duplicatedSets String indicating how to handle two or more elements containing exactly the same regions: \code{"stop"} interrupts the loading, \code{"remove"} keeps only the first occurrence, \code{"keep"} disables the check. Notice that identical sets are perfectly correlated and inflate the number of tests in the multiple-testing correction. Default: \code{"stop"}.
#' @param seqlevelsStyle String indicating the chromosome naming style to apply to all the elements, one among \code{"UCSC"} (chr1), \code{"Ensembl"} (1) or \code{"NCBI"}. When set to \code{NULL} the names are left as they are and the loading is interrupted if the elements use different styles. Default: \code{"UCSC"}.
#' @param genomeAssembly String indicating the genome assembly to store in the ranges seqinfo, e.g. \code{"hg38"}. Default: \code{NULL}, no assembly is assigned.
#' @param skipInvalid Logical value to indicate whether the elements that cannot be loaded (missing files, malformed tables, unsupported classes) must be skipped with a warning rather than interrupting the loading. Default: \code{FALSE}.
#' @param outputFormat String indicating the class of the returned object, one among \code{"RegionSetDE"}, \code{"GRangesList"} or \code{"list"}. Default: \code{"RegionSetDE"}.
#' @param verbose Logical value to indicate whether the loading messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE} object or, depending on \code{outputFormat}, a named \code{GRangesList} or a named list of \code{GRanges}.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regionList <- split(regionTable[, c("seqnames", "start", "end")],
#'                     regionTable$setName)
#'
#' regions <- loadRegions(regionList, genomeAssembly = "rn4", verbose = FALSE)
#' regions
#'
#' \dontrun{
#' regions <- loadRegions(list(promoters = "peaks/promoters.bed",
#'                             enhancers = grEnhancers,
#'                             "peaks/CTCF.narrowPeak"),
#'                        genomeAssembly = "hg38")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{splitLoadRegions}}, \code{\link{applyBlacklist}}, \code{\link{applyWhitelist}}
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<- seqnames genome genome<-
#' @importFrom IRanges reduce
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics sort unique duplicated start end
#' @importFrom methods is as new
#'
#' @export loadRegions

loadRegions <-
  function(regions,
           regionNames = NULL,
           keepMetadata = TRUE,
           sortRegions = TRUE,
           reduceRegions = FALSE,
           removeDuplicatedRegions = TRUE,
           duplicatedSets = "stop",
           seqlevelsStyle = "UCSC",
           genomeAssembly = NULL,
           skipInvalid = FALSE,
           outputFormat = "RegionSetDE",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    duplicatedSets <- tolower(duplicatedSets[1])
    if (!(duplicatedSets %in% c("stop", "remove", "keep"))) {
      stop("The 'duplicatedSets' parameter must be one among 'stop', 'remove' or 'keep'.", call. = FALSE)
    }

    outputFormat <- unname(c("regionsetde" = "RegionSetDE", "grangeslist" = "GRangesList", "list" = "list")[tolower(outputFormat[1])])
    if (is.na(outputFormat)) {
      stop("The 'outputFormat' parameter must be one among 'RegionSetDE', 'GRangesList' or 'list'.", call. = FALSE)
    }

    if (!is.null(seqlevelsStyle)) {
      seqlevelsStyle <- unname(c("ucsc" = "UCSC", "ensembl" = "Ensembl", "ncbi" = "NCBI")[tolower(seqlevelsStyle[1])])
      if (is.na(seqlevelsStyle)) {
        stop("The 'seqlevelsStyle' parameter must be one among 'UCSC', 'Ensembl' or 'NCBI'.", call. = FALSE)
      }
    }

    #---------------------------#
    # Uniform the input to list #
    #---------------------------#
    # Single objects and GRangesList are wrapped/unwrapped so that everything downstream is a plain list
    if (methods::is(regions, "RegionSetDE")) {
      regions <- as.list(regions@regions)
    } else if (methods::is(regions, "GRanges") | is.data.frame(regions)) {
      regions <- list(regions)
    } else if (methods::is(regions, "GRangesList")) {
      regions <- as.list(regions)
    } else if (is.character(regions)) {
      regions <- as.list(regions)
    } else if (!is.list(regions)) {
      stop("The 'regions' parameter must be a (named) list of BED file paths and/or GRanges objects.", call. = FALSE)
    }

    if (length(regions) == 0) {
      stop("The 'regions' list is empty.", call. = FALSE)
    }

    #-------------------------#
    # Define the region names #
    #-------------------------#
    # Fallback names: file base name for the paths, positional label for the objects
    defaultNames <-
      vapply(X = seq_along(regions),
             FUN.VALUE = character(1),
             FUN = function(i) {
               if (is.character(regions[[i]])) {
                 gsub("\\.(bed|narrowPeak|broadPeak|gtf|gff|gff3|txt|tsv)(\\.gz|\\.bz2)?$", "",
                      basename(regions[[i]]), ignore.case = TRUE)
               } else {
                 paste0("regions_", i)
               }})

    # User-provided names win over the list names, which in turn win over the fallbacks
    if (!is.null(regionNames)) {
      if (length(regionNames) != length(regions)) {
        stop(paste0("The number of 'regionNames' (", length(regionNames),
                    ") differs from the number of elements in 'regions' (", length(regions), ")."), call. = FALSE)
      }
      if (any(duplicated(regionNames))) {
        stop(paste0("The 'regionNames' contain duplicated values: ",
                    paste(unique(regionNames[duplicated(regionNames)]), collapse = ", "), "."), call. = FALSE)
      }
      names(regions) <- regionNames
    } else if (is.null(names(regions))) {
      names(regions) <- defaultNames
    }

    # Partially named lists leave empty slots
    emptyNames <- which(is.na(names(regions)) | names(regions) == "")
    names(regions)[emptyNames] <- defaultNames[emptyNames]

    # Files sharing a base name across directories are renamed, but the user must know about it
    if (any(duplicated(names(regions)))) {
      warning(paste0("Duplicated region names (", paste(unique(names(regions)[duplicated(names(regions))]), collapse = ", "),
                     ") have been made unique by adding a numeric suffix."), call. = FALSE)
      names(regions) <- make.unique(names(regions), sep = "_")
    }

    if (any(duplicated(names(regions)))) {
      stop("The region names could not be made unique, please provide explicit 'regionNames'.", call. = FALSE)
    }

    #-----------------#
    # Failure handler #
    #-----------------#
    # Loading problems are fatal unless the user explicitly asks to skip the faulty elements
    handleFailure <-
      function(errorMessage) {
        if (skipInvalid == TRUE) {
          warning(paste0(errorMessage, " The element will be skipped."), call. = FALSE)
          return(NULL)
        } else {
          stop(errorMessage, call. = FALSE)
        }
      }

    #----------------------------#
    # Chromosome name conversion #
    #----------------------------#
    # Scaffolds and custom contigs make the seqlevelsStyle conversion fail, in that case only the chr prefix is edited
    harmonizeSeqlevels <-
      function(gr, style) {
        newLevels <- GenomeInfoDb::seqlevels(gr)

        if (style == "UCSC") {
          newLevels <- ifelse(grepl("^chr", newLevels), newLevels, paste0("chr", newLevels))
          newLevels <- gsub("^chrMT$", "chrM", newLevels)
        } else {
          newLevels <- gsub("^chr", "", newLevels)
          newLevels <- gsub("^M$", "MT", newLevels)
        }

        if (any(duplicated(newLevels))) {
          stop("The chromosome names could not be converted to the requested style without collisions.", call. = FALSE)
        }

        GenomeInfoDb::seqlevels(gr) <- newLevels
        return(gr)
      }

    #------------------#
    # Load the regions #
    #------------------#
    if (verbose == TRUE) {message("Loading regions:")}

    # Each element is converted to GRanges depending on its class
    regionList <-
      lapply(X = seq_along(regions),
             FUN = function(i) {
               element <- regions[[i]]

               if (methods::is(element, "GRanges")) {
                 gr <- element
               } else if (is.data.frame(element)) {
                 gr <- GenomicRanges::makeGRangesFromDataFrame(element, keep.extra.columns = TRUE)
               } else if (is.character(element) & length(element) == 1) {
                 # The reader interrupts on a faulty file, the error is re-routed to honour 'skipInvalid'
                 gr <- tryCatch(expr = .readRegionFile(filePath = element, header = FALSE, asGRanges = TRUE),
                                error = function(e) {return(handleFailure(conditionMessage(e)))})
               } else {
                 gr <- handleFailure(paste0("The element '", names(regions)[i],
                                            "' is not a file path, a data.frame or a GRanges."))
               }

               # Reached only when 'skipInvalid' is active, otherwise the loading has already been interrupted
               if (is.null(gr)) {return(NULL)}

               ### Chromosome naming
               # Harmonized here so that the duplication check compares the sets on the same coordinate space
               if (!is.null(seqlevelsStyle)) {
                 gr <- tryCatch(expr = {GenomeInfoDb::seqlevelsStyle(gr) <- seqlevelsStyle; gr},
                                error = function(e) {return(harmonizeSeqlevels(gr, seqlevelsStyle))})
               }

               ### Post-processing of the ranges
               if (reduceRegions == TRUE) {
                 gr <- IRanges::reduce(gr)
               } else if (removeDuplicatedRegions == TRUE) {
                 # Repeated coordinates within a file would be counted twice by any overlap-based test
                 nDuplicated <- sum(BiocGenerics::duplicated(gr))
                 if (nDuplicated > 0) {
                   gr <- BiocGenerics::unique(gr)
                   if (verbose == TRUE) {message(paste0("  ", names(regions)[i], ": ", nDuplicated, " duplicated regions removed"))}
                 }
               }

               if (keepMetadata == FALSE | reduceRegions == TRUE) {S4Vectors::mcols(gr) <- NULL} # reduction makes the metadata meaningless
               if (sortRegions == TRUE) {gr <- BiocGenerics::sort(gr, ignore.strand = TRUE)}
               if (!is.null(genomeAssembly)) {GenomeInfoDb::genome(gr) <- genomeAssembly}

               if (verbose == TRUE) {message(paste0("  ", names(regions)[i], ": ", length(gr), " regions"))}

               return(gr)
             })

    # Names are re-assigned before dropping the skipped elements, to keep the correspondence with the input
    names(regionList) <- names(regions)
    regionList <- regionList[!vapply(regionList, is.null, logical(1))]

    if (length(regionList) == 0) {
      stop("None of the elements provided could be loaded.", call. = FALSE)
    }

    #-------------------------#
    # Check the naming styles #
    #-------------------------#
    # A mixture of UCSC and Ensembl names returns zero overlaps downstream without raising any error
    if (is.null(seqlevelsStyle)) {
      hasChrPrefix <- vapply(regionList, function(gr) {any(grepl("^chr", GenomeInfoDb::seqlevels(gr)))}, logical(1))

      if (length(unique(hasChrPrefix)) > 1) {
        stop(paste0("The region sets mix UCSC ('chr1') and Ensembl/NCBI ('1') chromosome names: ",
                    paste0(names(regionList), " [", ifelse(hasChrPrefix, "UCSC", "Ensembl/NCBI"), "]", collapse = ", "),
                    ". Use the 'seqlevelsStyle' parameter to harmonize them."), call. = FALSE)
      }
    }

    #---------------------------#
    # Check the duplicated sets #
    #---------------------------#
    # Sets are compared on sorted coordinates, so that the order of the regions does not matter
    identicalRanges <-
      function(a, b) {
        a <- BiocGenerics::sort(a, ignore.strand = TRUE)
        b <- BiocGenerics::sort(b, ignore.strand = TRUE)

        return(identical(as.character(GenomeInfoDb::seqnames(a)), as.character(GenomeInfoDb::seqnames(b))) &
                 identical(BiocGenerics::start(a), BiocGenerics::start(b)) &
                 identical(BiocGenerics::end(a), BiocGenerics::end(b)))
      }

    if (duplicatedSets != "keep" & length(regionList) > 1) {
      # Only the sets with the same number of regions can be identical, the pairwise comparison is restricted to those
      setLengths <- vapply(regionList, length, numeric(1))
      duplicatedIdx <- integer(0)
      duplicatedMsg <- character(0)

      for (i in seq_along(regionList)[-1]) {
        candidates <- setdiff(which(setLengths[seq_len(i - 1)] == setLengths[i]), duplicatedIdx)

        for (j in candidates) {
          if (identicalRanges(regionList[[i]], regionList[[j]]) == TRUE) {
            duplicatedIdx <- c(duplicatedIdx, i)
            duplicatedMsg <- c(duplicatedMsg, paste0("'", names(regionList)[i], "' (identical to '", names(regionList)[j], "')"))
            break
          }
        }
      }

      if (length(duplicatedIdx) > 0) {
        if (duplicatedSets == "stop") {
          stop(paste0("Some region sets contain exactly the same regions: ", paste(duplicatedMsg, collapse = ", "),
                      ". Redundant sets are perfectly correlated and bias the multiple-testing correction. ",
                      "Set 'duplicatedSets' to 'remove' or 'keep' to proceed."), call. = FALSE)
        } else {
          warning(paste0("Redundant region sets removed: ", paste(duplicatedMsg, collapse = ", "), "."), call. = FALSE)
          regionList <- regionList[-duplicatedIdx]
        }
      }
    }

    #--------------------#
    # Export the regions #
    #--------------------#
    if (outputFormat == "list") {return(regionList)}

    # The coercion fails when the elements carry different metadata columns
    grangesList <-
      tryCatch(expr = methods::as(regionList, "GRangesList"),
               error = function(e) {
                 warning("The elements carry different metadata columns and could not be combined, the metadata have been dropped.", call. = FALSE)
                 return(methods::as(lapply(regionList, function(gr) {S4Vectors::mcols(gr) <- NULL; return(gr)}), "GRangesList"))
               })

    if (outputFormat == "GRangesList") {return(grangesList)}

    ### RegionSetDE object, the loading parameters travel with the regions
    return(methods::new("RegionSetDE",
                        regions = grangesList,
                        genome.assembly = genomeAssembly,
                        seqlevels.style = ifelse(is.null(seqlevelsStyle), NA_character_, seqlevelsStyle),
                        filtering.log = data.frame(step = character(0), region.set = character(0),
                                                   n.before = numeric(0), n.after = numeric(0), n.removed = numeric(0),
                                                   stringsAsFactors = FALSE),
                        parameters = list(loadRegions = list(keepMetadata = keepMetadata,
                                                             sortRegions = sortRegions,
                                                             reduceRegions = reduceRegions,
                                                             removeDuplicatedRegions = removeDuplicatedRegions,
                                                             duplicatedSets = duplicatedSets,
                                                             seqlevelsStyle = seqlevelsStyle,
                                                             genomeAssembly = genomeAssembly))))
  } # END function
