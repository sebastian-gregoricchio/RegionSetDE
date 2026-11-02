# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .setKeys
#'
#' @description Builds the key identifying a region row inside a set, the set name and the region identifier pasted together. It is what \code{\link{geneSetsToRegionSets}} returns and what \code{\link{testRegionSets}} resolves a supplied set against, and it survives the collapsing of the tiles, which row positions do not.
#'
#' @param x Data.frame or \code{GRanges} carrying the \code{region.set} and \code{region.id} columns.
#'
#' @return A character vector with one key per row.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#'
#' @keywords internal

.setKeys <-
  function(x) {

    keyTable <- if (methods::is(x, "GRanges")) {as.data.frame(S4Vectors::mcols(x), stringsAsFactors = FALSE)} else {x}

    if (!all(c("region.set", "region.id") %in% colnames(keyTable))) {
      stop("The region rows carry no 'region.set' and 'region.id' columns, which is what identifies a region inside a set.", call. = FALSE)
    }

    return(paste(as.character(keyTable$region.set), as.character(keyTable$region.id), sep = "|"))
  } # END function




#' @title .asGeneSetList
#'
#' @description Accepts the two shapes gene sets arrive in, a named list of gene identifiers or a two-column table of the kind \code{clusterProfiler} calls a \code{TERM2GENE}, and returns the list.
#'
#' @param geneSets Named list of character vectors, or a data.frame with one column naming the set and one naming the gene.
#'
#' @return A named list of character vectors.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#'
#' @keywords internal

.asGeneSetList <-
  function(geneSets) {

    if (is.data.frame(geneSets)) {
      if (ncol(geneSets) < 2) {
        stop("The 'geneSets' parameter given as a table needs two columns, one naming the set and one naming the gene.", call. = FALSE)
      }

      # The named columns win, the first two otherwise, which is the TERM2GENE order
      setColumn <- if ("set.name" %in% colnames(geneSets)) {"set.name"} else {colnames(geneSets)[1]}
      geneColumn <- if ("gene.id" %in% colnames(geneSets)) {"gene.id"} else {setdiff(colnames(geneSets), setColumn)[1]}

      setTable <-
        data.frame(set.name = as.character(geneSets[[setColumn]]),
                   gene.id = as.character(geneSets[[geneColumn]]),
                   stringsAsFactors = FALSE) %>%
        dplyr::filter(!is.na(.data$set.name), !is.na(.data$gene.id))

      return(split(setTable$gene.id, setTable$set.name))
    }

    if (!is.list(geneSets)) {
      stop("The 'geneSets' parameter must be a named list of gene identifiers, or a two-column table.", call. = FALSE)
    }

    if (is.null(names(geneSets)) | any(is.na(names(geneSets))) | any(names(geneSets) == "")) {
      stop("Every element of the 'geneSets' parameter must be named, the names are what the sets are reported under.", call. = FALSE)
    }

    return(lapply(geneSets, function(setGenes) {unique(as.character(setGenes[!is.na(setGenes)]))}))
  } # END function




#' @title geneSetsToRegionSets
#'
#' @description Turns gene sets into region sets, by mapping every gene of a set through the links stored by \code{\link{annotateRegions}}, so that a pathway becomes the regions assigned to its genes and can be handed straight to \code{\link{testRegionSets}}.
#'
#' This is what makes pathway-level testing of differential regions go through \code{camera} and \code{fry} rather than through a ranked gene list. The directionality problem that ordinarily has to be solved before a gene list can be written, one gene carrying several regions moving in different directions, is not solved here: it is left to the test, which reports a self-contained and a mixed p-value and so answers it rather than deciding it in advance.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or \code{RegionSetDE.results} object carrying the links.
#' @param geneSets Named list of gene identifiers, one element per set, or a two-column table naming the set and the gene.
#' @param minWeight Numeric value with the smallest \code{link.weight} a link must carry to be used. Default: \code{NULL}, every link is used whatever its weight.
#' @param linkType Character vector with the kinds of link to follow, among \code{"promoter"} and \code{"distal"} as \code{\link{annotateRegions}} labelled them with its \code{tssRegion}. This is how the promoter-proximal and the distal analysis are run separately. Default: \code{NULL}, every link.
#' @param minRegions Numeric value with the smallest number of regions a set must reach to be returned. Note that this counts regions and not genes. Default: \code{2}, the fewest a set level test can work with.
#' @param maxRegions Numeric value with the largest number of regions a set may reach. Default: \code{NULL}, no upper bound.
#' @param regionSets Character vector with the names of the region sets the regions are taken from, which is how the analysis is restricted to one class of element, promoters for instance. Default: \code{NULL}, every set of the object.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A named list of character vectors, one element per gene set, holding the keys of the regions reached by that set. The keys are the set name and the region identifier pasted with a vertical bar, which is the form \code{\link{testRegionSets}} accepts through its \code{regionSets} argument. A \code{geneSetSummary} attribute carries one row per set with the number of genes it holds, the number of them the links reach, and the number of regions.
#'
#' @details Gene sets overlap and region sets of this package do not. A promoter belonging to fifty pathways is one row of the object that fifty sets reach, which is why this returns keys rather than a \code{GRangesList}: a list of ranges would have to repeat the row fifty times, and the counting would then be wrong in fifty ways. \code{\link{testRegionSets}} resolves the keys against its own rows, after any collapsing of the tiles, so the sets stay correct whatever the counting level.
#'
#' The comparison universe follows from that. With sets supplied this way, each one is compared against every other row of the object rather than against the other sets, which is the right comparison for a pathway: the question is whether the regions of this pathway moved more than the regions of the rest of the genome that was tested, not whether they moved more than the regions of the other pathways. It also means the universe is the annotated and tested regions and nothing wider, which is the background any enrichment of this kind should have and the one most often got wrong.
#'
#' What this does not fix is the bias that comes before it. A gene with a long regulatory locus carries more regions than a gene with a short one, so a set rich in such genes reaches more rows, and \code{camera} corrects for the correlation between those rows without correcting for how many of them there are. Where that matters, compare the result against a count-based method built for it, \code{chipenrich} or \code{polyenrich} for instance, run on the regions rather than on the genes.
#'
#' \code{minWeight} is the one place a link weight is used, and it is used as a filter rather than as a weight: neither \code{camera} nor \code{fry} accepts a per-region weight, so a contact probability can decide whether a link is trusted enough to be followed and cannot make one link count more than another. Asking for it on links that carry no weight stops, rather than silently keeping everything, since a filter that quietly does nothing is worse than one that is not there.
#'
#' Splitting promoter-proximal from distal links is usually worth doing, and there are two ways to reach it. \code{linkType} filters the links themselves, on the labels \code{\link{annotateRegions}} wrote from its \code{tssRegion}, which is the route when the regions are one undifferentiated set of peaks. \code{regionSets} selects whole region sets instead, which is the route when they were already split at load time, by \code{loadConsensusPeaks(regionSets = , regionMode = "split")} for instance. Either way the point is the same: promoter-proximal and distal regions of one gene behave differently often enough that mixing them is where much of the apparent disagreement within a gene comes from.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:8, 1:3]
#' linkTable <- data.frame(regionTable,
#'                         gene.id = c("geneA", "geneA", "geneB", "geneB",
#'                                     "geneC", "geneC", "geneD", "geneD"))
#'
#' counts <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)
#'
#' geneSets <- list(pathwayOne = c("geneA", "geneB"),
#'                  pathwayTwo = c("geneC", "geneD", "geneE"))
#'
#' regionSets <- geneSetsToRegionSets(counts, geneSets = geneSets, verbose = FALSE)
#' lengths(regionSets)
#' attr(regionSets, "geneSetSummary")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{annotateRegions}}, \code{\link{regionLinks}}, \code{\link{testRegionSets}}
#'
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors mcols
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export geneSetsToRegionSets

geneSetsToRegionSets <-
  function(object,
           geneSets,
           minWeight = NULL,
           linkType = NULL,
           minRegions = 2,
           maxRegions = NULL,
           regionSets = NULL,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .countsOfObject(object)
    links <- regionLinks(counts)

    if (is.null(links)) {
      stop("The object carries no region to gene links, run annotateRegions() on it first.", call. = FALSE)
    }

    geneSets <- .asGeneSetList(geneSets)

    if (length(geneSets) == 0) {
      stop("The 'geneSets' parameter holds no set.", call. = FALSE)
    }

    minRegions <- as.integer(minRegions[1])
    if (is.na(minRegions) | minRegions < 2) {
      stop("The 'minRegions' parameter must be at least 2, a set level test cannot run on fewer regions.", call. = FALSE)
    }

    #-------------------------------#
    # Links kept                    #
    #-------------------------------#
    # A region recorded with no gene stays in the object, and in the universe, while reaching no set
    links <- dplyr::filter(links, !is.na(.data$gene.id))

    if (!is.null(linkType)) {
      if (!("link.type" %in% colnames(links))) {
        stop("The 'linkType' parameter was given but the links carry no 'link.type' column. Annotate the regions through 'txdb' or 'genes', which labels them, or leave 'linkType' at NULL.", call. = FALSE)
      }

      storedTypes <- unique(stats::na.omit(links$link.type))
      absentTypes <- setdiff(linkType, storedTypes)

      if (length(absentTypes) > 0) {
        stop("The links carry no ", paste(absentTypes, collapse = " or "), " link. They hold: ",
             paste(storedTypes, collapse = ", "), ".", call. = FALSE)
      }

      links <- dplyr::filter(links, .data$link.type %in% linkType)
    }

    if (!is.null(minWeight)) {
      if (!("link.weight" %in% colnames(links))) {
        stop("The 'minWeight' parameter was given but the links carry no 'link.weight' column, so there is nothing to filter on. Annotate the regions with weighted links, or leave 'minWeight' at NULL.", call. = FALSE)
      }

      minWeight <- as.numeric(minWeight[1])
      if (is.na(minWeight)) {
        stop("The 'minWeight' parameter must be a number, or NULL.", call. = FALSE)
      }

      linksBefore <- nrow(links)
      links <- dplyr::filter(links, .data$link.weight >= minWeight)

      if (nrow(links) == 0) {
        stop("No link carries a weight of at least ", minWeight, ", the highest one is ", max(regionLinks(counts)$link.weight, na.rm = TRUE), ".", call. = FALSE)
      }

      if (isTRUE(verbose)) {
        message(linksBefore - nrow(links), " of the ", linksBefore, " links fall below a weight of ", minWeight, " and are not followed.")
      }
    }

    #-------------------------------#
    # Regions reached by a link     #
    #-------------------------------#
    regions <- SummarizedExperiment::rowRanges(counts)
    rowTable <- data.frame(region = .regionKeys(regions),
                           set.key = .setKeys(regions),
                           region.set = as.character(S4Vectors::mcols(regions)$region.set),
                           stringsAsFactors = FALSE)

    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, unique(rowTable$region.set))
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }

      rowTable <- dplyr::filter(rowTable, .data$region.set %in% regionSets)
    }

    # The same coordinates can be held by several sets, so a link reaches every row carrying them
    # rather than the first of them, which is what match() would have given
    keysByRegion <- split(rowTable$set.key, rowTable$region)

    geneToKeys <- split(links$region, links$gene.id)
    geneToKeys <- lapply(geneToKeys,
                         function(linkedRegions) {
                           unique(unlist(keysByRegion[intersect(linkedRegions, names(keysByRegion))], use.names = FALSE))
                         })
    geneToKeys <- geneToKeys[lengths(geneToKeys) > 0]

    #-------------------------------#
    # One region set per gene set   #
    #-------------------------------#
    setKeyList <- lapply(geneSets, function(setGenes) {unique(unlist(geneToKeys[intersect(setGenes, names(geneToKeys))], use.names = FALSE))})

    summaryTable <-
      data.frame(set.name = names(geneSets),
                 n.genes = lengths(geneSets),
                 n.genes.linked = vapply(geneSets, function(setGenes) {length(intersect(setGenes, names(geneToKeys)))}, integer(1)),
                 n.regions = lengths(setKeyList),
                 stringsAsFactors = FALSE,
                 row.names = NULL)

    #-------------------------------#
    # Sets of a usable size         #
    #-------------------------------#
    keptSets <- summaryTable$n.regions >= minRegions

    if (!is.null(maxRegions)) {
      maxRegions <- as.integer(maxRegions[1])
      if (is.na(maxRegions) | maxRegions < minRegions) {
        stop("The 'maxRegions' parameter must be at least as large as 'minRegions'.", call. = FALSE)
      }
      keptSets <- keptSets & summaryTable$n.regions <= maxRegions
    }

    if (!any(keptSets)) {
      stop("Not one gene set reaches ", minRegions, " regions. The links reach ", length(geneToKeys),
           " genes in all, check that the identifiers of 'geneSets' follow the same namespace as the links.", call. = FALSE)
    }

    setKeyList <- setKeyList[keptSets]
    summaryTable <- summaryTable[keptSets, , drop = FALSE]

    if (isTRUE(verbose)) {
      message(length(setKeyList), " of the ", length(geneSets), " gene sets reach between ",
              min(summaryTable$n.regions), " and ", max(summaryTable$n.regions), " regions and are returned.")
    }

    attr(setKeyList, "geneSetSummary") <- summaryTable

    return(setKeyList)
  } # END function




#' @title .setIndexList
#'
#' @description Resolves the sets a set level test runs on, from either the \code{region.set} column of the rows or a named list of region keys supplied by the user. The two differ in one way that matters: the sets of the column are a partition, one row belonging to one set, while supplied sets overlap freely, which is what gene sets do.
#'
#' @param regionStats Data.frame with the per-region statistics, carrying \code{region.set} and \code{region.id}.
#' @param regionSets Character vector with the names of the sets to test, or a named list of region keys. Default: \code{NULL}, every set of the column.
#'
#' @return A named list with the row positions of every set.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @keywords internal

.setIndexList <-
  function(regionStats,
           regionSets = NULL) {

    #-------------------------------#
    # Sets of the object            #
    #-------------------------------#
    if (is.null(regionSets) | is.character(regionSets)) {
      setNames <- unique(as.character(regionStats$region.set))

      if (!is.null(regionSets)) {
        absentSets <- setdiff(regionSets, setNames)
        if (length(absentSets) > 0) {
          stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
        }
        setNames <- as.character(regionSets)
      }

      return(stats::setNames(lapply(setNames, function(setName) {which(regionStats$region.set == setName)}), setNames))
    }

    #-------------------------------#
    # Sets supplied by the user     #
    #-------------------------------#
    if (!is.list(regionSets)) {
      stop("The 'regionSets' parameter must be a character vector with the names of the sets, or a named list of region keys.", call. = FALSE)
    }

    if (is.null(names(regionSets)) | any(is.na(names(regionSets))) | any(names(regionSets) == "")) {
      stop("Every element of the 'regionSets' parameter must be named, the names are what the sets are reported under.", call. = FALSE)
    }

    rowKeys <- .setKeys(regionStats)

    setIndexList <- lapply(regionSets, function(setKeys) {which(rowKeys %in% as.character(setKeys))})
    names(setIndexList) <- names(regionSets)

    # A set that reaches nothing is a namespace mismatch rather than an empty set, and it is worth saying so once
    emptySets <- lengths(setIndexList) == 0

    if (all(emptySets)) {
      stop("Not one region of the supplied sets matches a row of the object. The sets name regions such as ",
           paste(utils::head(unlist(regionSets, use.names = FALSE), 2), collapse = ", "), " while the rows are named ",
           paste(utils::head(rowKeys, 2), collapse = ", "), ". Build the sets with geneSetsToRegionSets() on this object.", call. = FALSE)
    }

    if (any(emptySets)) {
      warning(sum(emptySets), " of the supplied sets match no row of the object and are dropped: ",
              paste(utils::head(names(setIndexList)[emptySets], 3), collapse = ", "), ".", call. = FALSE)
      setIndexList <- setIndexList[!emptySets]
    }

    return(setIndexList)
  } # END function
