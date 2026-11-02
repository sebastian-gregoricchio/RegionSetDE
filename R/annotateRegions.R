# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .regionKeys
#'
#' @description Builds the coordinate key used to match a link table to the rows of an object. Coordinates rather than identifiers, because a link table computed elsewhere knows where a region is and not what the package decided to call it. The same coordinates appearing in two sets return the same key twice, so one link reaches both rows, which is what a region shared between two sets should do.
#'
#' @param x \code{GRanges} with the regions, or a data.frame carrying \code{seqnames}, \code{start} and \code{end}.
#'
#' @return A character vector with one key per region, written as \code{chr:start-end}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom methods is
#'
#' @keywords internal

.regionKeys <-
  function(x) {

    ### Integers pasted by hand, format() would turn the large coordinates into scientific notation
    if (methods::is(x, "GRanges")) {
      return(paste0(as.character(GenomeInfoDb::seqnames(x)), ":",
                    format(BiocGenerics::start(x), scientific = FALSE, trim = TRUE), "-",
                    format(BiocGenerics::end(x), scientific = FALSE, trim = TRUE)))
    }

    if (!all(c("seqnames", "start", "end") %in% colnames(x))) {
      stop("The coordinate key needs the 'seqnames', 'start' and 'end' columns.", call. = FALSE)
    }

    return(paste0(as.character(x$seqnames), ":",
                  format(as.integer(x$start), scientific = FALSE, trim = TRUE), "-",
                  format(as.integer(x$end), scientific = FALSE, trim = TRUE)))
  } # END function




#' @title .linksFromGRanges
#'
#' @description Turns a \code{GRanges} whose metadata columns name the genes into the long link table the package stores. Convenient because this is the shape \code{ChIPseeker::annotatePeak} and most contact callers return.
#'
#' @param linkRanges \code{GRanges} with at least a gene identifier among its metadata columns.
#'
#' @return A data.frame in the long link format.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom dplyr rename any_of %>%
#'
#' @keywords internal

.linksFromGRanges <-
  function(linkRanges) {

    linkTable <- data.frame(seqnames = as.character(GenomeInfoDb::seqnames(linkRanges)),
                            start = BiocGenerics::start(linkRanges),
                            end = BiocGenerics::end(linkRanges),
                            stringsAsFactors = FALSE)

    metadataTable <- as.data.frame(S4Vectors::mcols(linkRanges), stringsAsFactors = FALSE)

    ### Names used by ChIPseeker and by the contact callers, mapped onto the ones of the package
    metadataTable <-
      metadataTable %>%
      dplyr::rename(dplyr::any_of(c("gene.id" = "geneId",
                                    "gene.id" = "gene",
                                    "gene.symbol" = "SYMBOL",
                                    "gene.symbol" = "symbol",
                                    "distance.to.tss" = "distanceToTSS",
                                    "distance.to.tss" = "distanceToTss",
                                    "link.weight" = "score",
                                    "link.type" = "annotation")))

    return(cbind(linkTable, metadataTable))
  } # END function




#' @title .validateLinkTable
#'
#' @description Checks a region-to-gene link table and brings it onto the chromosome names of the regions. The two directions of mismatch are not treated alike: links landing on regions the object does not hold are dropped, since a contact map is routinely computed on a wider set of peaks than the one being tested, while regions the table says nothing about are kept with a missing gene, since dropping them would quietly shrink the universe of every enrichment computed afterwards.
#'
#' @param linkTable Data.frame or \code{GRanges} with the links.
#' @param regions \code{GRanges} with the regions of the object.
#' @param weightSource String naming where the weights come from, recorded alongside them. Default: \code{"supplied"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A data.frame with the validated links, carrying the \code{region} column used for the matching.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter mutate group_by summarise n select all_of %>%
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom methods is
#' @importFrom stats na.omit setNames
#' @importFrom utils head
#'
#' @keywords internal

.validateLinkTable <-
  function(linkTable,
           regions,
           weightSource = "supplied",
           verbose = TRUE) {

    #-------------------------------#
    # Shape of the input            #
    #-------------------------------#
    if (methods::is(linkTable, "GRanges")) {
      linkTable <- .linksFromGRanges(linkTable)
    }

    if (!is.data.frame(linkTable)) {
      stop("The 'linkTable' parameter must be a data.frame or a GRanges.", call. = FALSE)
    }

    if (nrow(linkTable) == 0) {
      stop("The 'linkTable' parameter holds no link.", call. = FALSE)
    }

    if (!("gene.id" %in% colnames(linkTable))) {
      stop("The 'linkTable' parameter requires a 'gene.id' column.", call. = FALSE)
    }

    ### Three ways of naming a region, because the sources differ: contact callers emit coordinates,
    ### a table written by hand usually carries the identifiers the package assigned
    byCoordinates <- all(c("seqnames", "start", "end") %in% colnames(linkTable))
    byRange <- "region" %in% colnames(linkTable)
    byIdentifier <- "region.id" %in% colnames(linkTable)

    if (!any(byCoordinates, byRange, byIdentifier)) {
      stop("The 'linkTable' parameter must name the regions through 'seqnames', 'start' and 'end', through a 'region' column, or through a 'region.id' column.", call. = FALSE)
    }

    linkTable <- dplyr::mutate(linkTable, gene.id = as.character(.data$gene.id))

    #-------------------------------#
    # Weights as they were given    #
    #-------------------------------#
    # A missing weight is ambiguous between "no weight was measured" and "the weight is zero", and
    # guessing between the two is worse than refusing the table
    if ("link.weight" %in% colnames(linkTable)) {
      linkTable <- dplyr::mutate(linkTable, link.weight = as.numeric(.data$link.weight))

      unusableWeight <- is.na(linkTable$link.weight) | linkTable$link.weight < 0

      if (any(unusableWeight)) {
        stop("The 'linkTable' parameter holds ", sum(unusableWeight),
             " links with a negative or missing 'link.weight'. Drop those links, or set their weight to zero explicitly.", call. = FALSE)
      }

      linkTable$weight.source <- as.character(weightSource[1])
    }

    #-------------------------------#
    # Key of every region           #
    #-------------------------------#
    regionKeys <- .regionKeys(regions)

    if (byCoordinates) {
      # The chromosomes of the table follow the names of the regions, one by one, so that a table
      # computed on Ensembl names reaches regions written the UCSC way. A name that matches nothing
      # is left as it was and reported below with the links it carries, rather than warned about twice
      translatedNames <- .translateChromosomeNames(chromosomeNames = as.character(linkTable$seqnames),
                                                  targetSeqlevels = GenomeInfoDb::seqlevels(regions))

      linkTable$seqnames <- ifelse(is.na(translatedNames), as.character(linkTable$seqnames), translatedNames)
      linkTable$region <- .regionKeys(linkTable)

    } else if (byRange) {
      linkTable$region <- as.character(linkTable$region)

    } else {
      # The identifiers of the object are not unique across the sets, one link then reaches every set holding it
      identifierKeys <- stats::setNames(regionKeys, as.character(.regionIdentifiers(regions)))
      linkTable$region <- unname(identifierKeys[as.character(linkTable$region.id)])
    }

    #-------------------------------#
    # Duplicated pairs              #
    #-------------------------------#
    # The same region linked twice to the same gene is an upstream problem rather than something to collapse,
    # since a duplicated link would count that gene twice in every aggregation
    duplicatedPairs <-
      linkTable %>%
      dplyr::filter(!is.na(.data$gene.id)) %>%
      dplyr::group_by(.data$region, .data$gene.id) %>%
      dplyr::summarise(n.links = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(.data$n.links > 1)

    if (nrow(duplicatedPairs) > 0) {
      stop("The 'linkTable' parameter holds ", nrow(duplicatedPairs), " duplicated region-gene pairs, the first being ",
           duplicatedPairs$region[1], " linked to ", duplicatedPairs$gene.id[1], ".", call. = FALSE)
    }

    #-------------------------------#
    # The two kinds of mismatch     #
    #-------------------------------#
    unmatchedLinks <- is.na(linkTable$region) | !(linkTable$region %in% regionKeys)

    if (all(unmatchedLinks)) {
      stop("Not one link of the 'linkTable' parameter falls on a region of the object. The table names regions such as ",
           paste(utils::head(stats::na.omit(linkTable$region), 2), collapse = ", "), " while the object holds ",
           paste(utils::head(regionKeys, 2), collapse = ", "), ".", call. = FALSE)
    }

    if (any(unmatchedLinks)) {
      if (isTRUE(verbose)) {
        message(sum(unmatchedLinks), " links fall on regions absent from the object and have been dropped.")
      }
      linkTable <- linkTable[!unmatchedLinks, , drop = FALSE]
    }

    unlinkedRegions <- sum(!(regionKeys %in% linkTable$region))

    if (unlinkedRegions > 0 & isTRUE(verbose)) {
      message(unlinkedRegions, " regions carry no link and keep a missing 'gene.id'. They stay in the object, and in the universe of any enrichment computed from it.")
    }

    #-------------------------------#
    # Column order                  #
    #-------------------------------#
    linkColumns <- c("region", "gene.id", "gene.symbol", "distance.to.tss", "link.weight", "weight.source", "link.type")
    linkTable <- dplyr::select(linkTable, dplyr::all_of(intersect(linkColumns, colnames(linkTable))))

    return(linkTable)
  } # END function




#' @title .regionIdentifiers
#'
#' @description Returns the identifier of every region, the \code{region.id} metadata column when the regions carry one and their names otherwise.
#'
#' @param regions \code{GRanges} with the regions of the object.
#'
#' @return A character vector with one identifier per region.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#'
#' @keywords internal

.regionIdentifiers <-
  function(regions) {

    if ("region.id" %in% colnames(S4Vectors::mcols(regions))) {
      return(as.character(S4Vectors::mcols(regions)$region.id))
    }

    if (!is.null(names(regions))) {
      return(as.character(names(regions)))
    }

    return(.regionKeys(regions))
  } # END function




#' @title .collapseLinks
#'
#' @description Picks one link per region, for the columns written next to the regions. The closest transcription start site when the distances are known, the heaviest link when the weights are, and the first one otherwise. It is a view for reading and for labelling a plot, never the quantity anything is computed from: that is the full table, which \code{\link{regionLinks}} returns.
#'
#' @param linkTable Data.frame with the validated links.
#' @param regionKeys Character vector with the key of every region of the object.
#'
#' @return A data.frame with one row per region, in the order of \code{regionKeys}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr arrange filter mutate %>%
#' @importFrom rlang .data
#'
#' @keywords internal

.collapseLinks <-
  function(linkTable,
           regionKeys) {

    emptyTable <- data.frame(gene.id = rep(NA_character_, length(regionKeys)),
                             gene.symbol = rep(NA_character_, length(regionKeys)),
                             distance.to.tss = rep(NA_integer_, length(regionKeys)),
                             n.links = rep(0L, length(regionKeys)),
                             stringsAsFactors = FALSE)

    assignedLinks <- dplyr::filter(linkTable, !is.na(.data$gene.id))

    if (nrow(assignedLinks) == 0) {
      return(emptyTable)
    }

    ### The ordering decides which link represents the region, nearest first and heaviest after it
    assignedLinks <-
      assignedLinks %>%
      dplyr::mutate(sort.distance = if ("distance.to.tss" %in% colnames(assignedLinks)) {abs(.data$distance.to.tss)} else {0},
                    sort.weight = if ("link.weight" %in% colnames(assignedLinks)) {-.data$link.weight} else {0}) %>%
      dplyr::arrange(.data$region, .data$sort.distance, .data$sort.weight)

    # Sorted already, so the first row of every region is its representative. A grouped slice would
    # give the same answer and take far longer on a contact table holding millions of links
    bestLinks <- assignedLinks[!duplicated(assignedLinks$region), , drop = FALSE]

    linkCounts <- table(assignedLinks$region)
    matchedRows <- match(regionKeys, bestLinks$region)

    collapsedTable <-
      data.frame(gene.id = bestLinks$gene.id[matchedRows],
                 gene.symbol = if ("gene.symbol" %in% colnames(bestLinks)) {bestLinks$gene.symbol[matchedRows]} else {NA_character_},
                 distance.to.tss = if ("distance.to.tss" %in% colnames(bestLinks)) {as.integer(bestLinks$distance.to.tss[matchedRows])} else {NA_integer_},
                 n.links = as.integer(ifelse(is.na(linkCounts[regionKeys]), 0L, linkCounts[regionKeys])),
                 stringsAsFactors = FALSE)

    return(collapsedTable)
  } # END function




#' @title .nearestGeneLinks
#'
#' @description Builds the links of the regions by proximity to a transcription start site, reading the sites from a \code{TxDb} or an \code{EnsDb}. Regions whose nearest site sits further than \code{maxDistance} are written into the table with a missing gene rather than left out of it, so that an unassigned region stays distinguishable from a region nobody looked at.
#'
#' @param regions \code{GRanges} with the regions of the object.
#' @param txdb \code{TxDb} or \code{EnsDb} object with the gene models. Default: \code{NULL}, the sites come from \code{geneRanges}.
#' @param geneRanges \code{GRanges} with the genes, or their start sites, carrying a gene identifier. Default: \code{NULL}, the genes come from \code{txdb}.
#' @param maxDistance Numeric value with the distance beyond which no gene is assigned, in base pairs. Default: \code{NULL}, every region takes its nearest site.
#' @param flankDistance Numeric value: every gene whose start site falls within it is linked, rather than the nearest one alone. Default: \code{NULL}, the nearest site.
#' @param tssRegion Numeric vector of length two with the window around a start site, signed and read on the strand of the gene, inside which a link counts as a promoter one. Default: \code{c(-3000, 3000)}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A data.frame in the long link format, with a \code{link.type} column holding \code{"promoter"} or \code{"distal"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges resize distanceToNearest findOverlaps
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom BiocGenerics start end strand
#'
#' @keywords internal

.nearestGeneLinks <-
  function(regions,
           txdb = NULL,
           geneRanges = NULL,
           maxDistance = NULL,
           flankDistance = NULL,
           tssRegion = c(-3000, 3000),
           verbose = TRUE) {

    ### The genes come from a TxDb, or from ranges the caller already has
    if (is.null(geneRanges)) {
      if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
        stop("Reading the genes from a 'txdb' needs the 'GenomicFeatures' package. Install it, or pass the genes as ranges through 'genes'.", call. = FALSE)
      }

      geneRanges <- tryCatch(GenomicFeatures::genes(txdb),
                             error = function(e) {
                               stop("The genes could not be read from the 'txdb' parameter: ", conditionMessage(e), call. = FALSE)
                             })
    }

    # The identifier is a metadata column or the names of the ranges, whichever the caller used
    geneIdentifiers <-
      if ("gene.id" %in% colnames(S4Vectors::mcols(geneRanges))) {
        as.character(S4Vectors::mcols(geneRanges)$gene.id)
      } else if (!is.null(names(geneRanges))) {
        as.character(names(geneRanges))
      } else {
        stop("The genes carry no identifier: give them names, or a 'gene.id' metadata column.", call. = FALSE)
      }

    geneLabels <- if ("gene.symbol" %in% colnames(S4Vectors::mcols(geneRanges))) {
      as.character(S4Vectors::mcols(geneRanges)$gene.symbol)
    } else {
      NULL
    }

    ### One start site per gene, taken on its strand
    startSites <- GenomicRanges::resize(geneRanges, width = 1, fix = "start")
    names(startSites) <- geneIdentifiers

    # The annotation and the regions rarely agree on the chromosome names
    startSites <- .matchSeqlevels(x = startSites, targetSeqlevels = GenomeInfoDb::seqlevels(regions), verbose = verbose)

    ### Either every gene whose start site falls in the flank, which is what makes a region
    ### over a bidirectional promoter keep both of them, or the single nearest one
    if (!is.null(flankDistance)) {
      flankDistance <- as.numeric(flankDistance[1])

      if (is.na(flankDistance) | flankDistance < 0) {
        stop("The 'flankDistance' parameter must be a positive number, or NULL.", call. = FALSE)
      }

      flankHits <- GenomicRanges::findOverlaps(query = regions,
                                               subject = startSites,
                                               maxgap = as.integer(flankDistance),
                                               ignore.strand = TRUE)

      queryIndex <- S4Vectors::queryHits(flankHits)
      subjectIndex <- S4Vectors::subjectHits(flankHits)

      if (length(queryIndex) == 0) {
        stop("No gene start site falls within ", format(flankDistance, big.mark = ","),
             " bp of any region. Widen 'flankDistance', or check that the regions and the annotation share an assembly.", call. = FALSE)
      }

    } else {
      nearestSites <- GenomicRanges::distanceToNearest(regions, startSites, ignore.strand = TRUE, select = "arbitrary")

      if (length(nearestSites) == 0) {
        stop("No region could be matched to a transcription start site, check that the regions and the annotation share an assembly.", call. = FALSE)
      }

      queryIndex <- S4Vectors::queryHits(nearestSites)
      subjectIndex <- S4Vectors::subjectHits(nearestSites)
    }

    # Signed distance, negative upstream of the site, the convention ChIPseeker follows
    siteStart <- BiocGenerics::start(startSites)[subjectIndex]
    regionCentre <- (BiocGenerics::start(regions)[queryIndex] + BiocGenerics::end(regions)[queryIndex]) %/% 2L
    onMinusStrand <- as.character(BiocGenerics::strand(startSites))[subjectIndex] == "-"
    signedDistance <- as.integer(ifelse(onMinusStrand, siteStart - regionCentre, regionCentre - siteStart))

    linkTable <- data.frame(region = .regionKeys(regions)[queryIndex],
                            gene.id = geneIdentifiers[subjectIndex],
                            distance.to.tss = signedDistance,
                            stringsAsFactors = FALSE)

    if (!is.null(geneLabels)) {
      linkTable$gene.symbol <- geneLabels[subjectIndex]
    }

    ### A link is a promoter one when the region sits inside the window around the start site,
    ### read on the strand of the gene, and distal otherwise. It is what splits the two analyses
    linkTable$link.type <- ifelse(signedDistance >= tssRegion[1] & signedDistance <= tssRegion[2],
                                  "promoter", "distal")

    ### Beyond maxDistance the nearest gene is a coordinate accident, the region is recorded as looked at
    ### and left unassigned, which keeps it in the universe while keeping it out of every gene set
    if (!is.null(maxDistance)) {
      maxDistance <- as.numeric(maxDistance[1])

      if (is.na(maxDistance) | maxDistance < 0) {
        stop("The 'maxDistance' parameter must be a positive number, or NULL.", call. = FALSE)
      }

      tooFar <- abs(linkTable$distance.to.tss) > maxDistance

      if (any(tooFar) & isTRUE(verbose)) {
        message(sum(tooFar), " regions sit further than ", format(maxDistance, big.mark = ","),
                " bp from any start site and are recorded as unassigned.")
      }

      linkTable$gene.id[tooFar] <- NA_character_
      linkTable$distance.to.tss[tooFar] <- NA_integer_
      linkTable$link.type[tooFar] <- NA_character_
    }

    ### Regions the flank reached nothing on are written back unassigned, so that they stay in the
    ### universe rather than disappearing from the object they were counted in
    if (!is.null(flankDistance)) {
      unreached <- setdiff(seq_along(regions), unique(queryIndex))

      if (length(unreached) > 0) {
        if (isTRUE(verbose)) {
          message(length(unreached), " regions have no gene start site within ",
                  format(flankDistance, big.mark = ","), " bp and are recorded as unassigned.")
        }

        emptyRows <- linkTable[rep(1L, length(unreached)), , drop = FALSE]
        emptyRows$region <- .regionKeys(regions)[unreached]
        emptyRows$gene.id <- NA_character_
        emptyRows$distance.to.tss <- NA_integer_
        emptyRows$link.type <- NA_character_
        if ("gene.symbol" %in% colnames(emptyRows)) {emptyRows$gene.symbol <- NA_character_}

        linkTable <- rbind(linkTable, emptyRows)
      }
    }

    return(linkTable)
  } # END function




#' @title .addGeneSymbols
#'
#' @description Adds the gene symbols to a link table, read from an \code{OrgDb}. Cosmetic: the symbols are for reading and for labelling, nothing is ever matched on them, since they are neither stable between releases nor unique.
#'
#' @param linkTable Data.frame with the links.
#' @param orgdb \code{OrgDb} object.
#' @param keyType String with the key type of \code{gene.id} in the \code{OrgDb}. Default: \code{"ENTREZID"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The link table with a \code{gene.symbol} column.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats na.omit
#'
#' @keywords internal

.addGeneSymbols <-
  function(linkTable,
           orgdb,
           keyType = "ENTREZID",
           verbose = TRUE) {

    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      if (isTRUE(verbose)) {
        message("The 'AnnotationDbi' package is not installed, the gene symbols are left out.")
      }
      return(linkTable)
    }

    assignedGenes <- unique(stats::na.omit(linkTable$gene.id))

    symbolMap <- tryCatch(suppressMessages(AnnotationDbi::mapIds(orgdb, keys = assignedGenes, column = "SYMBOL",
                                                                 keytype = keyType, multiVals = "first")),
                          error = function(e) {
                            if (isTRUE(verbose)) {
                              message("The gene symbols could not be read from the 'orgdb' parameter: ", conditionMessage(e))
                            }
                            return(NULL)
                          })

    if (is.null(symbolMap)) {
      return(linkTable)
    }

    linkTable$gene.symbol <- unname(symbolMap[linkTable$gene.id])

    return(linkTable)
  } # END function




#' @title annotateRegions
#'
#' @description Links the regions of an object to genes and stores the links inside it. The links are the weakest step of any functional analysis built on chromatin data, so they are kept as data rather than resolved on the fly: the full region-to-gene table travels with the object, one row per link, and a region linked to three genes keeps all three.
#'
#' Sources are accepted in three forms. Proximity to a transcription start site, read either from a \code{TxDb} or an \code{EnsDb} through \code{txdb}, or from ranges you already have through \code{genes}, which is what a GTF or a species with no \code{TxDb} leaves you with. A table of links of your own through \code{linkTable}, which is how the output of a contact-based method arrives, Hi-C loops, ABC scores, universalEPI predictions or anything else. Or a \code{GRanges} of regions carrying the gene identifiers in its metadata columns, also through \code{linkTable}, which is the shape \code{ChIPseeker::annotatePeak} returns.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object. Not a \code{RegionSetDE} object: the links live beside the counted rows, so the regions have to be counted first.
#' @param txdb \code{TxDb} or \code{EnsDb} object with the gene models, used when the links are built by proximity. Default: \code{NULL}.
#' @param genes \code{GRanges} with the genes, or their start sites, carrying their identifier in a \code{gene.id} metadata column or in the names of the ranges, and optionally a \code{gene.symbol} column. The links are built by proximity to them, the way \code{txdb} does it, which is the route to take for an annotation read from a GTF or for a species with no \code{TxDb}. Default: \code{NULL}.
#' @param orgdb \code{OrgDb} object used to add the gene symbols. Default: \code{NULL}, no symbols.
#' @param linkTable Data.frame or \code{GRanges} with links of your own. Default: \code{NULL}, the links are built from \code{txdb}.
#' @param maxDistance Numeric value with the distance beyond which no gene is assigned, in base pairs, when the nearest gene is taken. Regions further than this are recorded with a missing \code{gene.id} rather than dropped. Default: \code{NULL}, every region takes its nearest site whatever the distance.
#' @param flankDistance Numeric value: every gene whose start site falls within it is linked to the region, rather than the nearest one alone. This is what makes the proximity route many-to-many, and it is how a region over a bidirectional promoter keeps both of its genes. Regions reaching nothing are recorded unassigned. Cannot be used with \code{maxDistance}, which caps the nearest-gene route instead. Default: \code{NULL}, the nearest site.
#' @param tssRegion Numeric vector of length two with the window around a start site, signed and read on the strand of the gene, inside which a link is called a promoter one. It fills the \code{link.type} column and decides nothing else. Default: \code{c(-3000, 3000)}, the window \code{ChIPseeker} uses.
#' @param keyType String with the key type the gene identifiers follow in \code{orgdb}. Default: \code{"ENTREZID"}.
#' @param weightSource String naming where the weights of \code{linkTable} come from, stored next to them. Default: \code{"supplied"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The input object with the links stored in it. The full table is in \code{metadata(x)$links} and is read with \code{\link{regionLinks}}, while \code{gene.id}, \code{gene.symbol}, \code{distance.to.tss} and \code{n.links} are written next to the regions as one representative link per region, for reading and for labelling plots.
#'
#' @details A region-to-gene link is a prediction, and the quality of every gene-level result rests on it rather than on the statistics applied afterwards. Nearest transcription start site is wrong for a substantial share of distal elements, which skip over genes and reach promoters megabases away. That is the reason this function stores a table instead of a column: the proximity links can be replaced by measured or predicted contacts without recomputing anything else, and the two can be compared on the same object.
#'
#' The links can be many-to-many, and whether they are is a choice. A supplied \code{linkTable} carries whatever links it was built with, so a region overlapping a bidirectional promoter keeps both genes and a region inside an intron keeps the host gene alongside whatever else its contacts reach. The proximity route takes the single nearest start site by default, one link per region, and \code{flankDistance} is what turns it many-to-many: every gene whose start site falls within that distance is linked, each with its own signed distance. Either way nothing decides which link is right, and the number of links is not the number of regions nor the number of genes.
#'
#' \code{tssRegion} fills \code{link.type} with \code{"promoter"} or \code{"distal"}, the window being read on the strand of the gene so that a region upstream of a minus-strand gene sits at a higher coordinate than its start site. The column decides nothing by itself; it is there so that \code{\link{geneSetsToRegionSets}} and \code{\link{aggregateToGenes}} can run the promoter-proximal and the distal analysis separately, which is usually worth doing, since much of the apparent disagreement within a gene is promoters staying put while enhancers move. What this function does not do is classify a region into a genomic feature, promoter against intron against exon, with a priority order; that is what \code{ChIPseeker::annotatePeak} is for, and its output goes into \code{linkTable} directly.
#'
#' Weights are optional and they carry no assumed scale. An ABC score, a Hi-C loop score, a predicted contact probability and a distance decay are all accepted as they are, with \code{weightSource} recording which of them it was, and nothing is normalised. Supplying no weight at all is the normal case and is fully supported: the links are then unweighted, and every function reading them treats the links of a region as equally plausible. What a weight cannot do is reach the set-level tests, since neither \code{camera} nor \code{fry} accepts a per-region weight. Where it does reach is \code{\link{geneSetsToRegionSets}}, which uses it as a confidence filter through its \code{minWeight} argument.
#'
#' Two kinds of mismatch between the table and the object are handled differently on purpose. Links on regions the object does not hold are dropped, because a contact map is routinely computed over a wider peak set than the one being analysed. Regions the table says nothing about are kept, with a missing gene, because dropping them would shrink the universe of every enrichment run afterwards without anything in the output saying so. The counts of both are reported and recorded in the parameters.
#'
#' \code{maxDistance} makes a third case explicit. A region whose nearest start site is 800 kb away has a nearest gene in the arithmetic sense only, and assigning it is worse than leaving it alone. Those regions are written into the table with a missing \code{gene.id}, which keeps them in the universe, keeps them out of every gene set, and leaves the decision visible in the object rather than buried in the choice of which rows were loaded.
#'
#' Gene identifiers are matched exactly, as strings, and must follow the same namespace as the gene sets they will be tested against. Symbols are stored when an \code{OrgDb} is supplied and are never used for matching: they are not stable between annotation releases and not unique within one.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' # Links of your own, the shape a contact-based method returns. Weights are optional.
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:6, 1:3]
#' linkTable <- data.frame(regionTable,
#'                         gene.id = c("11287", "11287", "66968", "73824", "73824", "20683"),
#'                         link.weight = c(0.42, 0.08, 0.31, 0.57, 0.12, 0.66))
#'
#' counts <- annotateRegions(counts, linkTable = linkTable,
#'                           weightSource = "ABC score", verbose = FALSE)
#'
#' regionLinks(counts)
#' head(SummarizedExperiment::rowData(counts)$gene.id)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{regionLinks}}, \code{\link{geneSetsToRegionSets}}, \code{\link{testRegionSets}}
#'
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
#' @importFrom S4Vectors mcols mcols<- metadata metadata<-
#' @importFrom methods is new
#' @importFrom stats na.omit
#'
#' @export annotateRegions

annotateRegions <-
  function(object,
           txdb = NULL,
           genes = NULL,
           orgdb = NULL,
           linkTable = NULL,
           maxDistance = NULL,
           flankDistance = NULL,
           tssRegion = c(-3000, 3000),
           keyType = "ENTREZID",
           weightSource = "supplied",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .countsOfObject(object)

    givenSources <- c("txdb", "genes", "linkTable")[c(!is.null(txdb), !is.null(genes), !is.null(linkTable))]

    if (length(givenSources) == 0) {
      stop("One of the 'txdb', 'genes' or 'linkTable' parameters must be given, there is nothing to build the links from otherwise.", call. = FALSE)
    }

    if (length(givenSources) > 1) {
      stop("The ", paste(givenSources, collapse = " and "), " parameters cannot be used together, the links come from one source at a time.", call. = FALSE)
    }

    if (!is.null(genes) & !methods::is(genes, "GRanges")) {
      stop("The 'genes' parameter must be a GRanges.", call. = FALSE)
    }

    # One caps the nearest gene, the other replaces it with every gene around, so together they
    # would describe two different sets of links
    if (!is.null(maxDistance) & !is.null(flankDistance)) {
      stop("The 'maxDistance' and 'flankDistance' parameters cannot be used together: the first caps the nearest gene, the second links every gene within reach.", call. = FALSE)
    }

    tssRegion <- sort(as.numeric(tssRegion))
    if (length(tssRegion) != 2 | any(is.na(tssRegion))) {
      stop("The 'tssRegion' parameter must be two numbers, the window around a start site.", call. = FALSE)
    }

    regions <- SummarizedExperiment::rowRanges(counts)

    #-------------------------------#
    # Links of the regions          #
    #-------------------------------#
    if (!is.null(linkTable)) {
      linkSource <- "table"
      links <- .validateLinkTable(linkTable = linkTable, regions = regions,
                                  weightSource = weightSource, verbose = verbose)
    } else {
      linkSource <- if (is.null(genes)) {"nearestTSS"} else {"nearestGene"}
      links <- .nearestGeneLinks(regions = regions, txdb = txdb, geneRanges = genes,
                                 maxDistance = maxDistance, flankDistance = flankDistance,
                                 tssRegion = tssRegion, verbose = verbose)
      links <- .validateLinkTable(linkTable = links, regions = regions,
                                  weightSource = weightSource, verbose = FALSE)
    }

    if (!is.null(orgdb)) {
      links <- .addGeneSymbols(linkTable = links, orgdb = orgdb, keyType = keyType, verbose = verbose)
    }

    #-------------------------------#
    # Storage                       #
    #-------------------------------#
    regionKeys <- .regionKeys(regions)

    # The full table is the truth, the columns next to the regions are a view of it
    collapsedLinks <- .collapseLinks(linkTable = links, regionKeys = regionKeys)

    S4Vectors::mcols(regions)$gene.id <- collapsedLinks$gene.id
    S4Vectors::mcols(regions)$gene.symbol <- collapsedLinks$gene.symbol
    S4Vectors::mcols(regions)$distance.to.tss <- collapsedLinks$distance.to.tss
    S4Vectors::mcols(regions)$n.links <- collapsedLinks$n.links

    SummarizedExperiment::rowRanges(counts) <- regions
    S4Vectors::metadata(counts)$links <- links

    counts@parameters$annotateRegions <- list(source = linkSource,
                                              max.distance = if (is.null(maxDistance)) {NA_real_} else {as.numeric(maxDistance[1])},
                                              flank.distance = if (is.null(flankDistance)) {NA_real_} else {as.numeric(flankDistance[1])},
                                              tss.region = if (is.null(linkTable)) {tssRegion} else {NULL},
                                              weighted = "link.weight" %in% colnames(links),
                                              weight.source = if ("link.weight" %in% colnames(links)) {as.character(weightSource[1])} else {NA_character_},
                                              n.links = nrow(links),
                                              n.regions.linked = sum(!is.na(collapsedLinks$gene.id)),
                                              n.regions.unlinked = sum(is.na(collapsedLinks$gene.id)),
                                              n.genes = length(unique(stats::na.omit(links$gene.id))))

    if (isTRUE(verbose)) {
      message(nrow(links), " links stored, reaching ", length(unique(stats::na.omit(links$gene.id))), " genes over ",
              sum(!is.na(collapsedLinks$gene.id)), " of the ", length(regionKeys), " regions.")
    }

    return(.replaceCounts(object = object, counts = counts))
  } # END function




#' @title .countsOfObject
#'
#' @description Returns the counts object carried by any of the classes holding one, so that the annotation reads the rows wherever it is called from.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object. Not a \code{RegionSetDE} object: the links live beside the counted rows, so the regions have to be counted first.
#'
#' @return A \code{RegionSetDE.counts} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.countsOfObject <-
  function(object) {

    if (methods::is(object, "RegionSetDE.counts")) {
      return(object)
    }

    # RegionSetDE.setScores is left out on purpose: it holds no counts slot to read the rows from
    if (methods::is(object, "RegionSetDE.fit") | methods::is(object, "RegionSetDE.results") |
        methods::is(object, "RegionSetDE.setResults")) {

      if (nrow(object@counts) == 0) {
        stop("The object carries no counts, which is where the links are stored. Annotate the counts, or the fit, instead.", call. = FALSE)
      }

      return(object@counts)
    }

    stop("The 'object' parameter must be a RegionSetDE.counts, RegionSetDE.fit, RegionSetDE.results or RegionSetDE.setResults object.", call. = FALSE)
  } # END function




#' @title .replaceCounts
#'
#' @description Writes a counts object back where it came from, the counterpart of \code{\link{.countsOfObject}}.
#'
#' @param object Object the counts were taken from.
#' @param counts \code{RegionSetDE.counts} object to store.
#'
#' @return The input object carrying the new counts.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom methods is
#'
#' @keywords internal

.replaceCounts <-
  function(object,
           counts) {

    if (methods::is(object, "RegionSetDE.counts")) {
      return(counts)
    }

    object@counts <- counts

    # The results keep their own copy of the ranges, which would otherwise disagree with the counts.
    # The slot is reached only once the class is known, since & evaluates both of its sides
    if (methods::is(object, "RegionSetDE.results")) {
      if (length(object@regions) > 0) {
        object@regions <- .carryLinkColumns(regions = object@regions,
                                            annotatedRegions = SummarizedExperiment::rowRanges(counts))
      }
    }

    return(object)
  } # END function




#' @title .carryLinkColumns
#'
#' @description Copies the representative link columns from the annotated rows onto another set of ranges describing the same regions, matched on the coordinates.
#'
#' @param regions \code{GRanges} to write the columns onto.
#' @param annotatedRegions \code{GRanges} carrying the link columns.
#'
#' @return The input ranges with the link columns added.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols mcols<-
#'
#' @keywords internal

.carryLinkColumns <-
  function(regions,
           annotatedRegions) {

    matchedRows <- match(.regionKeys(regions), .regionKeys(annotatedRegions))

    for (linkColumn in c("gene.id", "gene.symbol", "distance.to.tss", "n.links")) {
      if (linkColumn %in% colnames(S4Vectors::mcols(annotatedRegions))) {
        S4Vectors::mcols(regions)[[linkColumn]] <- S4Vectors::mcols(annotatedRegions)[[linkColumn]][matchedRows]
      }
    }

    return(regions)
  } # END function




#' @title regionLinks
#'
#' @description Returns the region-to-gene links stored in an object by \code{\link{annotateRegions}}, one row per link. Also sets them, which is how a contact table computed elsewhere is attached to an object that never went through a \code{TxDb}.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object. Not a \code{RegionSetDE} object: the links live beside the counted rows, so the regions have to be counted first.
#' @param value Data.frame or \code{GRanges} with the links, for the replacement form.
#'
#' @return A data.frame with the columns \code{region} and \code{gene.id}, together with whatever else the links carried: \code{gene.symbol}, \code{distance.to.tss}, \code{link.weight}, \code{weight.source} and \code{link.type}. \code{NULL} when the object has not been annotated.
#'
#' @details A region linked to no gene is kept in the table with a missing \code{gene.id}, which is how an unassigned region stays distinguishable from a region nothing was computed for. Count the genes with \code{na.omit} rather than on the rows.
#'
#' The replacement form runs the same checks \code{\link{annotateRegions}} does and refuses a table that would be ambiguous, so attaching links by hand is as safe as building them.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:4, 1:3]
#' regionLinks(counts) <- data.frame(regionTable, gene.id = c("11287", "66968", "73824", "20683"))
#'
#' regionLinks(counts)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{annotateRegions}}, \code{\link{geneSetsToRegionSets}}
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom S4Vectors metadata
#'
#' @export regionLinks

setGeneric(name = "regionLinks", def = function(object) {standardGeneric("regionLinks")})

#' @rdname regionLinks
#' @export
setMethod(f = "regionLinks",
          signature = "RegionSetDE.provenance",
          definition = function(object) {
            return(S4Vectors::metadata(.countsOfObject(object))$links)
          })

#' @rdname regionLinks
#' @export
setGeneric(name = "regionLinks<-", def = function(object, value) {standardGeneric("regionLinks<-")})

#' @rdname regionLinks
#' @export
setMethod(f = "regionLinks<-",
          signature = "RegionSetDE.provenance",
          definition = function(object, value) {
            return(annotateRegions(object, linkTable = value, verbose = FALSE))
          })
