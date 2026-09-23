# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title loadBlacklist
#'
#' @description Returns a blacklist shipped with the package, ready for \code{\link{applyBlacklist}} or for the \code{excludeRegions} argument of \code{\link{loadConsensusPeaks}}. The files travel with the package, so nothing is downloaded and nothing depends on a hub being reachable.
#'
#' @param genome String with the genome assembly, such as \code{"hg38"}, \code{"mm10"} or \code{"hs1"} for T2T-CHM13v2.0. The usual aliases are understood, \code{"GRCh38"} and \code{"T2T"} for instance. Run \code{\link{availableRegionLists}} for the whole list.
#' @param assay String with the assay the list was built for, \code{"cutrun"} or \code{"cuttag"}. Default: \code{NULL}, the ENCODE blacklist, which is not tied to an assay.
#' @param source String with the laboratory or project the list comes from, as \code{\link{availableRegionLists}} prints it, needed only when a genome carries more than one list for the same assay. Default: \code{NULL}.
#' @param seqlevelsStyle String with the chromosome naming style of the output, one among \code{"UCSC"} (chr1), \code{"Ensembl"} (1) and \code{"NCBI"}, or \code{NULL} to keep the names of the file. Default: \code{"UCSC"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} with the regions of the list and their name, carrying the assembly in its \code{genome} and the source, the version and the reference in its metadata.
#'
#' @details Without an assay the list is the one meant for ChIP-seq and ATAC-seq, the regions of anomalous coverage found across many experiments: the ENCODE blacklist version 2 for the assemblies it covers, and for T2T-CHM13v2.0, which ENCODE never covered, the set the excluderanges authors built by running the same software on that assembly. Naming an assay returns instead the high signal regions of the CUT&RUN greenlist paper, built from hundreds of CUT&RUN or CUT&Tag libraries, one list per assay. They are not interchangeable: the first kind is about the genome, the other two are about what a given protocol does to it.
#'
#' The files are the published ones, converted to gzipped BED and nothing else. \code{inst/extdata/regionLists/SOURCES.md} records where each of them comes from, the version, the date it was taken and the licence.
#'
#' @examples
#' blacklist <- loadBlacklist("hg38")
#' head(blacklist, 3)
#'
#' # The regions the CUT&RUN protocol piles reads on, which the ENCODE list does not cover
#' cutrunBlacklist <- loadBlacklist("hg38", assay = "cutrun")
#'
#' # T2T-CHM13v2.0, under any of the names it goes by
#' t2tBlacklist <- loadBlacklist("T2T")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadGreenlist}}, \code{\link{availableRegionLists}}, \code{\link{applyBlacklist}}
#'
#' @export loadBlacklist

loadBlacklist <-
  function(genome,
           assay = NULL,
           source = NULL,
           seqlevelsStyle = "UCSC",
           verbose = TRUE) {

    return(.loadRegionList(type = "blacklist",
                           genome = genome,
                           assay = assay,
                           source = source,
                           seqlevelsStyle = seqlevelsStyle,
                           verbose = verbose))
  } # END function




#' @title loadGreenlist
#'
#' @description Returns a CUT&RUN or CUT&Tag greenlist shipped with the package, the regions whose background is consistent enough between experiments to normalise on. It goes to \code{\link{countGreenlist}}, which counts the libraries over it before \code{\link{normalizeCounts}} turns those counts into scaling factors.
#'
#' @param genome String with the genome assembly, \code{"hg38"} or \code{"mm39"} for the published lists. The usual aliases are understood, \code{"GRCh38"} for instance.
#' @param assay String with the assay the list was built for, \code{"cutrun"} or \code{"cuttag"}.
#' @param source String with the laboratory or project the list comes from, as \code{\link{availableRegionLists}} prints it, needed only when a genome carries more than one list for the same assay. Default: \code{NULL}.
#' @param seqlevelsStyle String with the chromosome naming style of the output, one among \code{"UCSC"} (chr1), \code{"Ensembl"} (1) and \code{"NCBI"}, or \code{NULL} to keep the names of the file. Default: \code{"UCSC"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} with the regions of the greenlist and their name, carrying the assembly in its \code{genome} and the source, the version and the reference in its metadata.
#'
#' @details The greenlist was built by measuring, over hundreds of public libraries, which 1 kb bins carry background of a consistent magnitude, keeping the most consistent of them and discarding anything within 5 kb of a gene so that genuine signal is not called noise. The reads landing there follow the amount of material sequenced rather than the factor being mapped, which is what makes them usable as an internal reference when no spike-in was added.
#'
#' The list is specific to the protocol. CUT&RUN and CUT&Tag have their own, and using one for the other means normalising on regions whose background was never shown to be consistent in that assay. There is no greenlist for mm10, only for mm39, so an mm10 analysis needs the mm39 list lifted over, and the lifted file should be kept beside the chain it came from rather than passed off as the published one.
#'
#' @examples
#' greenlist <- loadGreenlist("hg38", assay = "cutrun")
#' length(greenlist)
#' sum(BiocGenerics::width(greenlist))
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countGreenlist}}, \code{\link{normalizeCounts}}, \code{\link{loadBlacklist}}, \code{\link{availableRegionLists}}
#'
#' @export loadGreenlist

loadGreenlist <-
  function(genome,
           assay,
           source = NULL,
           seqlevelsStyle = "UCSC",
           verbose = TRUE) {

    if (missing(assay) || is.null(assay)) {
      stop("Name the assay the greenlist was built for, 'cutrun' or 'cuttag'. They are not interchangeable.", call. = FALSE)
    }

    return(.loadRegionList(type = "greenlist",
                           genome = genome,
                           assay = assay,
                           source = source,
                           seqlevelsStyle = seqlevelsStyle,
                           verbose = verbose))
  } # END function




#' @title availableRegionLists
#'
#' @description Lists the blacklists and greenlists shipped with the package, with the assembly, the assay, the source and the size of each of them.
#'
#' @param type String with the lists returned, either \code{"blacklist"} or \code{"greenlist"}. Default: \code{NULL}, both.
#' @param genome String with an assembly the lists are restricted to. Default: \code{NULL}, all of them.
#'
#' @return A data.frame with one row per list: the \code{type}, the \code{genome}, the \code{assay} it was built for or \code{"any"}, the \code{source}, the \code{version}, the number of regions, the base pairs they cover and the reference.
#'
#' @examples
#' availableRegionLists()
#'
#' availableRegionLists(type = "greenlist")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadBlacklist}}, \code{\link{loadGreenlist}}
#'
#' @importFrom dplyr filter arrange
#' @importFrom rlang .data
#'
#' @export availableRegionLists

availableRegionLists <-
  function(type = NULL,
           genome = NULL) {

    listIndex <- .regionListIndex()

    if (!is.null(type)) {
      listType <- tolower(type[1])
      if (!(listType %in% c("blacklist", "greenlist"))) {
        stop("The 'type' parameter must be either 'blacklist' or 'greenlist'.", call. = FALSE)
      }
      listIndex <- dplyr::filter(listIndex, .data$type == listType)
    }

    if (!is.null(genome)) {
      genomeName <- .resolveGenomeName(genome = genome)
      listIndex <- dplyr::filter(listIndex, .data$genome == genomeName)
    }

    listIndex <- dplyr::arrange(listIndex, .data$type, .data$genome, .data$assay)

    return(listIndex[, setdiff(colnames(listIndex), "file")])
  } # END function




#' @title .regionListIndex
#'
#' @description Reads the table indexing the region lists shipped with the package.
#'
#' @return A data.frame with one row per file.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils read.delim
#'
#' @keywords internal

.regionListIndex <-
  function() {

    indexPath <- system.file("extdata", "regionLists", "regionLists.tsv", package = "RegionSetDE")

    if (indexPath == "" | !file.exists(indexPath)) {
      stop("The table of the region lists is missing from the installation of the package.", call. = FALSE)
    }

    return(utils::read.delim(indexPath, stringsAsFactors = FALSE))
  } # END function




#' @title .resolveGenomeName
#'
#' @description Turns the name of an assembly into the one the files are indexed under, so that GRCh38 and hg38 reach the same list.
#'
#' @param genome String with the assembly given by the user.
#'
#' @return A string with the name used in the index, the input itself when it is not a known alias.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveGenomeName <-
  function(genome) {

    if (!is.character(genome) | length(genome) != 1) {
      stop("The 'genome' parameter must be a single assembly name, such as 'hg38'.", call. = FALSE)
    }

    # The assembly names of Ensembl and of the consortia, mapped onto the UCSC ones the files use
    genomeAliases <- c("grch37" = "hg19", "hg19" = "hg19",
                       "grch38" = "hg38", "hg38" = "hg38",
                       "t2t" = "hs1", "chm13" = "hs1", "chm13v2" = "hs1", "chm13v20" = "hs1",
                       "t2tchm13" = "hs1", "t2tchm13v20" = "hs1", "hs1" = "hs1",
                       "grcm38" = "mm10", "mm10" = "mm10",
                       "grcm39" = "mm39", "mm39" = "mm39",
                       "bdgp5" = "dm3", "dm3" = "dm3",
                       "bdgp6" = "dm6", "dm6" = "dm6",
                       "ws220" = "ce10", "ce10" = "ce10",
                       "wbcel235" = "ce11", "ce11" = "ce11")

    # 'T2T-CHM13v2.0' and 'chm13v2' name the same assembly, the punctuation is dropped before the lookup
    matchedGenome <- genomeAliases[gsub("[^a-z0-9]", "", tolower(genome))]

    return(unname(ifelse(is.na(matchedGenome), genome, matchedGenome)))
  } # END function




#' @title .loadRegionList
#'
#' @description Picks one of the region lists shipped with the package and reads it into a \code{GRanges}.
#'
#' @param type String with the type of list, \code{"blacklist"} or \code{"greenlist"}.
#' @param genome String with the genome assembly.
#' @param assay String with the assay the list was built for, or \code{NULL} for the lists tied to none.
#' @param source String with the source of the list, or \code{NULL} when the genome and the assay are enough to tell them apart.
#' @param seqlevelsStyle String with the chromosome naming style of the output, or \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A \code{GRanges} with the regions of the list.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb genome genome<-
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#' @importFrom BiocGenerics width
#'
#' @keywords internal

.loadRegionList <-
  function(type,
           genome,
           assay = NULL,
           source = NULL,
           seqlevelsStyle = "UCSC",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    listIndex <- .regionListIndex()
    listType <- type[1]
    genomeName <- .resolveGenomeName(genome = genome)

    # 'CUT&RUN', 'CUTnRUN' and 'cut_run' all mean the same list
    assayName <- if (is.null(assay)) {"any"} else {gsub("[^a-z0-9]", "", tolower(assay[1]))}

    typeIndex <- dplyr::filter(listIndex, .data$type == listType)
    genomeIndex <- dplyr::filter(typeIndex, .data$genome == genomeName)

    if (nrow(genomeIndex) == 0) {
      stop("No ", type, " is shipped for the genome '", genomeName, "'. Available: ",
           paste(sort(unique(typeIndex$genome)), collapse = ", "),
           ". See availableRegionLists().", call. = FALSE)
    }

    pickedList <- dplyr::filter(genomeIndex, .data$assay == assayName)

    if (nrow(pickedList) == 0) {
      stop("No ", type, " is shipped for the genome '", genomeName, "' and the assay '", assayName,
           "'. Available for this genome: ", paste(sort(unique(genomeIndex$assay)), collapse = ", "),
           ". See availableRegionLists().", call. = FALSE)
    }

    # A genome can carry several lists of one kind, in which case the source tells them apart
    if (!is.null(source)) {
      sourceName <- tolower(source[1])
      pickedList <- dplyr::filter(pickedList, tolower(.data$source) == sourceName)

      if (nrow(pickedList) == 0) {
        stop("No ", type, " of '", source[1], "' is shipped for the genome '", genomeName, "'. Available: ",
             paste(sort(unique(genomeIndex$source)), collapse = ", "), ".", call. = FALSE)
      }
    }

    if (nrow(pickedList) > 1) {
      stop("Several ", type, "s match the genome '", genomeName, "' and the assay '", assayName,
           "'. Name one through 'source': ", paste(sort(unique(pickedList$source)), collapse = ", "), ".", call. = FALSE)
    }

    #-------------------------------#
    # Read the file                 #
    #-------------------------------#
    listRanges <- .readRegionFile(filePath = system.file("extdata", "regionLists", pickedList$file, package = "RegionSetDE"),
                                  header = FALSE,
                                  asGRanges = TRUE)

    # The fourth column of these files holds the kind of region, which is worth keeping under a readable name
    columnNames <- colnames(S4Vectors::mcols(listRanges))
    columnNames[columnNames == "V4"] <- "name"
    colnames(S4Vectors::mcols(listRanges)) <- columnNames

    S4Vectors::mcols(listRanges) <- S4Vectors::mcols(listRanges)[, intersect("name", columnNames), drop = FALSE]

    if (!is.null(seqlevelsStyle)) {
      listRanges <- .styleSeqlevels(x = listRanges, seqlevelsStyle = seqlevelsStyle)
    }

    GenomeInfoDb::genome(listRanges) <- genomeName

    S4Vectors::metadata(listRanges) <- list(type = type,
                                            genome = genomeName,
                                            assay = assayName,
                                            source = pickedList$source,
                                            version = pickedList$version,
                                            reference = pickedList$reference)

    if (isTRUE(verbose)) {
      message("The ", pickedList$source, " ", type, " ", pickedList$version, " for ", genomeName,
              ": ", length(listRanges), " regions covering ",
              format(round(sum(as.numeric(BiocGenerics::width(listRanges))) / 1e6, 1), nsmall = 1), " Mb.")
    }

    return(listRanges)
  } # END function




#' @title .styleSeqlevels
#'
#' @description Renames the chromosomes of a \code{GRanges} to a naming style, editing the prefix by hand when the conversion of \code{GenomeInfoDb} does not cover the contigs at hand.
#'
#' @param x \code{GRanges} object.
#' @param seqlevelsStyle String with the style, one among \code{"UCSC"}, \code{"Ensembl"} and \code{"NCBI"}.
#'
#' @return The object with its chromosomes renamed.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#'
#' @keywords internal

.styleSeqlevels <-
  function(x,
           seqlevelsStyle = "UCSC") {

    seqlevelsStyle <- unname(c("ucsc" = "UCSC", "ensembl" = "Ensembl", "ncbi" = "NCBI")[tolower(seqlevelsStyle[1])])

    if (is.na(seqlevelsStyle)) {
      stop("The 'seqlevelsStyle' parameter must be one among 'UCSC', 'Ensembl' or 'NCBI'.", call. = FALSE)
    }

    # Scaffolds and custom contigs make the conversion of GenomeInfoDb fail, there only the prefix is edited
    manualConversion <-
      function(gr) {
        newLevels <- GenomeInfoDb::seqlevels(gr)

        if (seqlevelsStyle == "UCSC") {
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

    return(tryCatch(expr = {GenomeInfoDb::seqlevelsStyle(x) <- seqlevelsStyle; x},
                    error = function(e) {manualConversion(x)},
                    warning = function(w) {manualConversion(x)}))
  } # END function
