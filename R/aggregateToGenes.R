# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .resolveLinkWeights
#'
#' @description Decides what weights the regions of a gene, from what the links actually carry. A weight measured or predicted for the link itself comes first, the decay of the distance to the start site second, and nothing at all last, which is a legitimate answer rather than a failure.
#'
#' @param linkTable Data.frame with the links of the regions kept.
#' @param decay Numeric value with the distance, in base pairs, at which the exponential decay falls to \code{1/e}.
#'
#' @return A list with \code{weight}, a numeric vector parallel to \code{linkTable}, and \code{source}, one among \code{"link.weight"}, \code{"distance"} and \code{"none"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveLinkWeights <-
  function(linkTable,
           decay = 5000) {

    ### A weight that was measured or predicted for the link beats anything derived from the coordinates
    if ("link.weight" %in% colnames(linkTable)) {
      return(list(weight = as.numeric(linkTable$link.weight), source = "link.weight"))
    }

    ### Failing that, the distance to the start site, which is a guess about contact and not a measure of it
    if ("distance.to.tss" %in% colnames(linkTable)) {
      decayedWeight <- exp(-abs(as.numeric(linkTable$distance.to.tss)) / decay)
      decayedWeight[is.na(decayedWeight)] <- 0

      return(list(weight = decayedWeight, source = "distance"))
    }

    return(list(weight = rep(1, nrow(linkTable)), source = "none"))
  } # END function




#' @title aggregateToGenes
#'
#' @description Turns a per-region result into a per-gene table, through the links stored by \code{\link{annotateRegions}}. One row per gene, with the regions that reached it counted in each direction and the disagreement between them reported rather than resolved.
#'
#' This is the table a conventional over-representation analysis needs, and the ranked list a gene-level GSEA needs, for the cases where the set level route of \code{\link{geneSetsToRegionSets}} and \code{\link{testRegionSets}} is not what is wanted. It is the weaker of the two routes and it is honest about why.
#'
#' @param results \code{RegionSetDE.results} object returned by \code{\link{testRegions}}, carrying links.
#' @param method Character vector with what to compute, among \code{"simes"}, which combines the p-values of the regions of a gene, and \code{"potential"}, which sums their fold changes weighted by how plausible their links are. Default: \code{"simes"}.
#' @param decay Numeric value with the distance, in base pairs, at which the weight of a link falls to \code{1/e}, used by \code{"potential"} when the links carry no weight of their own. Default: \code{5000}.
#' @param linkType Character vector with the kinds of link to follow, among \code{"promoter"} and \code{"distal"} as \code{\link{annotateRegions}} labelled them with its \code{tssRegion}. This is how the promoter-proximal and the distal table are built separately. Default: \code{NULL}, every link.
#' @param regionSets Character vector with the names of the region sets the regions are taken from, which is the other way of restricting the table to one class of element, when the sets were already split at load time. Default: \code{NULL}, every set of the object.
#' @param adjustMethod String with the multiple testing correction applied across the genes. Default: \code{"BH"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.geneResults} object. Its table, read with \code{\link{resultsTable}}, holds one row per gene:
#' \itemize{
#'   \item \code{gene.id}, \code{gene.symbol}: the gene, and its symbol when the links carry one.
#'   \item \code{n.regions}: how many regions of the object reached it.
#'   \item \code{n.up}, \code{n.down}: how many of those are labelled \code{"up"} and \code{"down"} in the region table.
#'   \item \code{direction}: \code{"up"} or \code{"down"} when the regions that changed agree, \code{"mixed"} when they do not, \code{"unchanged"} when none of them changed.
#'   \item \code{log2FC}: the fold change of the region carrying the combined p-value, with \code{"simes"}.
#'   \item \code{p.value}, \code{FDR}: the Simes combination over the regions of the gene, corrected across the genes, with \code{"simes"}.
#'   \item \code{score}, \code{score.magnitude}: with \code{"potential"}, the weighted sum of the fold changes of the regions of the gene, \code{sum(w * log2FC)}, and the same sum on their absolute values, \code{sum(w * abs(log2FC))}.
#' }
#'
#' @details The reason this function exists alongside \code{\link{geneSetsToRegionSets}}, and the reason it is the second choice, is worth stating. Aggregating to genes throws away the ranking of the regions and replaces several measurements with one, which is exactly the step the set level route avoids. Use it when an external tool needs a gene table, when a reviewer asks for a conventional over-representation result, or when the gene rather than the region is genuinely the unit of the question. For asking which pathways responded, test the sets.
#'
#' \code{n.up} and \code{n.down} are read from the \code{diff.status} column of the region table rather than recomputed, so they honour the \code{FDR} and \code{log2FC} the regions were labelled with, including after \code{\link{updateThresholds}} has changed them. The gene table therefore cannot disagree with the region table about what counts as up.
#'
#' \code{direction} is the answer to the question that makes this whole exercise awkward. A gene reached by one region that opened and another that closed is \code{"mixed"}, and nothing here decides which of the two wins, because there is no defensible way to decide it: cancelling them loses a real signal, and picking the stronger one is arbitrary. Treat the mixed genes as their own category rather than dropping them. They are where a rearranged regulatory landscape shows up, and a gene with three regions up and three down is a different biological statement from a gene with nothing happening, which a signed score cannot tell apart.
#'
#' \code{"simes"} asks whether any region of the gene changed. It is relatively relaxed, valid under the positive dependence that regions of one locus have, and corrects for the number of regions a gene carries, which a minimum p-value does not. The fold change reported beside it belongs to the region that carried the p-value, not to an average, so a gene with one strong region among many quiet ones is described by that region.
#'
#' \code{"potential"} is a score and not a test: it has no error model and no p-value of its own, which is why the Simes columns are worth keeping alongside it. For a gene reached by regions \code{i} carrying weights \code{w[i]}, \code{score} is \code{sum(w * log2FC)} and \code{score.magnitude} is \code{sum(w * abs(log2FC))}. The first is signed, so a gene gaining at one region and losing at another lands near zero; the second ignores sign and does not, so the same gene scores high. Genes high in the second and near zero in the first are the rearranged ones, the same signal \code{direction} reports as a label.
#'
#' \code{score} is the ranked gene-level statistic a GSEA takes, so \code{fgsea::fgseaMultilevel(pathways = geneSets, stats = setNames(score, gene.id))} is the whole of that route; ranking on \code{score.magnitude} asks which sets were remodelled rather than which shifted, and \code{sign(score) * -log10(p.value)} ranks on significance rather than on effect. No wrapper for it is provided, since the call is one line and wrapping it would add a dependency without adding a decision. What does need stating in the methods is that the genes of the table are not independent, so the permutation behind any such test is anti-conservative.
#'
#' The weights come from the links. A weight measured or predicted for the link itself is used as it is, whatever it measures, which is what \code{weightSource} in \code{\link{annotateRegions}} records. Failing that, the distance to the start site is turned into \code{exp(-d / decay)}, which is a guess about contact rather than a measure of it. Failing both, every region of a gene counts the same. Which of the three applied is in the \code{weighting} slot and is printed by \code{show}.
#'
#' Two things this does not fix, both worth a line in the methods. Genes are not independent here: a region linked to three genes contributes to all three, so the rows of the table share data and any test run on them afterwards is anti-conservative to an unquantified degree. And the number of regions a gene carries depends on the length of its regulatory locus, so a gene set rich in long-locus genes is easier to hit; a count-based method built for that bias, \code{chipenrich} or \code{polyenrich}, is the comparison to run when it matters.
#'
#' A tiled result that was not combined back to regions is refused, since its rows are tiles and a gene would then be weighted by how finely its regions were cut. Test it with \code{combine = TRUE} first.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # Links of your own, four regions per gene
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fit@counts))[, 1:3]
#' geneNumber <- ceiling(nrow(regionTable) / 4)
#'
#' annotatedFit <- annotateRegions(fit,
#'                                 linkTable = data.frame(regionTable,
#'                                                        gene.id = rep(paste0("gene", seq_len(geneNumber)),
#'                                                                      each = 4)[seq_len(nrow(regionTable))]),
#'                                 verbose = FALSE)
#'
#' results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' geneResults <- aggregateToGenes(results, method = c("simes", "potential"), verbose = FALSE)
#' geneResults
#'
#' geneTable <- resultsTable(geneResults)
#' table(geneTable$direction)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{annotateRegions}}, \code{\link{geneSetsToRegionSets}}, \code{\link{testRegionSets}}, \code{\link{updateThresholds}}
#'
#' @importFrom csaw combineTests
#' @importFrom dplyr filter arrange desc %>%
#' @importFrom rlang .data
#' @importFrom methods is new
#' @importFrom stats p.adjust na.omit
#'
#' @export aggregateToGenes

aggregateToGenes <-
  function(results,
           method = "simes",
           decay = 5000,
           linkType = NULL,
           regionSets = NULL,
           adjustMethod = "BH",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (methods::is(results, "RegionSetDE.resultsList")) {
      stop("The 'results' parameter holds several contrasts. Take the one to aggregate out of it first, with results[[\"", names(results@results)[1], "\"]].", call. = FALSE)
    }

    if (!methods::is(results, "RegionSetDE.results")) {
      stop("The 'results' parameter must be a RegionSetDE.results object returned by testRegions().", call. = FALSE)
    }

    method <- unique(tolower(as.character(method)))
    if (length(method) == 0 | !all(method %in% c("simes", "potential"))) {
      stop("The 'method' parameter must contain only 'simes' and 'potential'.", call. = FALSE)
    }

    decay <- as.numeric(decay[1])
    if (is.na(decay) | decay <= 0) {
      stop("The 'decay' parameter must be a positive number.", call. = FALSE)
    }

    # Rows that are tiles would weight a gene by how finely its regions were cut rather than by how many it has
    if (identical(results@counting.level, "tile") & isFALSE(results@combination$applied)) {
      stop("The rows of this result are tiles rather than regions. Run testRegions() with 'combine = TRUE' before aggregating to genes.", call. = FALSE)
    }

    links <- regionLinks(results)
    if (is.null(links)) {
      stop("The object carries no region to gene links, run annotateRegions() on it first.", call. = FALSE)
    }

    #-------------------------------#
    # Regions entering the table    #
    #-------------------------------#
    regionTable <- results@results

    if (!all(c("region.set", "log2FC", "p.value", "diff.status") %in% colnames(regionTable))) {
      stop("The region table lacks the columns the aggregation reads: region.set, log2FC, p.value and diff.status.", call. = FALSE)
    }

    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, unique(as.character(regionTable$region.set)))
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      regionTable <- dplyr::filter(regionTable, .data$region.set %in% regionSets)
    }

    regionTable$region <- .regionKeys(regionTable)

    #-------------------------------#
    # One row per link              #
    #-------------------------------#
    # A region linked to three genes reaches all three, so the rows below are links and not regions
    usableLinks <- dplyr::filter(links, !is.na(.data$gene.id), .data$region %in% regionTable$region)

    if (!is.null(linkType)) {
      if (!("link.type" %in% colnames(usableLinks))) {
        stop("The 'linkType' parameter was given but the links carry no 'link.type' column. Annotate the regions through 'txdb' or 'genes', which labels them, or leave 'linkType' at NULL.", call. = FALSE)
      }

      storedTypes <- unique(stats::na.omit(usableLinks$link.type))
      absentTypes <- setdiff(linkType, storedTypes)

      if (length(absentTypes) > 0) {
        stop("The links carry no ", paste(absentTypes, collapse = " or "), " link. They hold: ",
             paste(storedTypes, collapse = ", "), ".", call. = FALSE)
      }

      usableLinks <- dplyr::filter(usableLinks, .data$link.type %in% linkType)
    }

    if (nrow(usableLinks) == 0) {
      stop("Not one link falls on a region of this result. Check that the links were built on the same object, and that 'regionSets' has not excluded all of them.", call. = FALSE)
    }

    matchedRows <- match(usableLinks$region, regionTable$region)

    linkStatistics <-
      data.frame(gene.id = usableLinks$gene.id,
                 gene.symbol = if ("gene.symbol" %in% colnames(usableLinks)) {usableLinks$gene.symbol} else {NA_character_},
                 log2FC = regionTable$log2FC[matchedRows],
                 p.value = regionTable$p.value[matchedRows],
                 diff.status = as.character(regionTable$diff.status[matchedRows]),
                 stringsAsFactors = FALSE)

    #-------------------------------#
    # Simes over the regions        #
    #-------------------------------#
    # combineTests reads the grouping as a factor, its level order is what the output rows follow
    geneGroups <- factor(linkStatistics$gene.id, levels = unique(linkStatistics$gene.id))

    combinedTable <- csaw::combineTests(ids = geneGroups,
                                        tab = data.frame(logFC = linkStatistics$log2FC,
                                                         PValue = linkStatistics$p.value,
                                                         stringsAsFactors = FALSE),
                                        pval.col = "PValue",
                                        fc.col = "logFC")

    combinedTable <- as.data.frame(combinedTable)
    geneIds <- rownames(combinedTable)

    #-------------------------------#
    # Directions of the regions     #
    #-------------------------------#
    # Read from diff.status rather than recomputed, so the gene table and the region table cannot
    # disagree about what counts as up, updateThresholds() included
    upCounts <- rowsum(x = as.numeric(linkStatistics$diff.status == "up"), group = as.character(geneGroups))[geneIds, 1]
    downCounts <- rowsum(x = as.numeric(linkStatistics$diff.status == "down"), group = as.character(geneGroups))[geneIds, 1]

    geneDirection <- ifelse(upCounts > 0 & downCounts > 0, "mixed",
                            ifelse(upCounts > 0, "up",
                                   ifelse(downCounts > 0, "down", "unchanged")))

    symbolByGene <- linkStatistics$gene.symbol[match(geneIds, linkStatistics$gene.id)]

    geneTable <- data.frame(gene.id = geneIds,
                            gene.symbol = symbolByGene,
                            n.regions = as.integer(combinedTable$num.tests),
                            n.up = as.integer(upCounts),
                            n.down = as.integer(downCounts),
                            direction = geneDirection,
                            stringsAsFactors = FALSE,
                            row.names = NULL)

    if ("simes" %in% method) {
      geneTable$log2FC <- linkStatistics$log2FC[combinedTable$rep.test]
      geneTable$p.value <- combinedTable$PValue
      geneTable$FDR <- stats::p.adjust(combinedTable$PValue, method = adjustMethod)
    }

    #-------------------------------#
    # Regulatory potential          #
    #-------------------------------#
    weighting <- NA_character_

    if ("potential" %in% method) {
      linkWeights <- .resolveLinkWeights(linkTable = usableLinks, decay = decay)
      weighting <- linkWeights$source

      weightedSum <- rowsum(x = linkWeights$weight * linkStatistics$log2FC, group = as.character(geneGroups))[geneIds, 1]
      magnitudeSum <- rowsum(x = linkWeights$weight * abs(linkStatistics$log2FC), group = as.character(geneGroups))[geneIds, 1]

      geneTable$score <- as.numeric(weightedSum)
      geneTable$score.magnitude <- as.numeric(magnitudeSum)
    }

    #-------------------------------#
    # Order and assemble            #
    #-------------------------------#
    geneTable <-
      if ("simes" %in% method) {
        dplyr::arrange(geneTable, .data$p.value)
      } else {
        dplyr::arrange(geneTable, dplyr::desc(abs(.data$score)))
      }

    rownames(geneTable) <- NULL

    if (isTRUE(verbose)) {
      message(nrow(geneTable), " genes aggregated from ", nrow(usableLinks), " links over ",
              length(unique(usableLinks$region)), " regions, of which ",
              sum(geneTable$direction == "mixed"), " are reached by regions moving both ways.")
    }

    geneResults <- methods::new(Class = "RegionSetDE.geneResults",
                                results = geneTable,
                                contrast = results@contrast,
                                method = method,
                                links = usableLinks,
                                weighting = if (is.na(weighting)) {"none"} else {weighting},
                                engine = results@engine,
                                thresholds = c(results@thresholds, list(adjust.method = adjustMethod, decay = decay)),
                                blacklist = results@blacklist,
                                greylist = results@greylist,
                                whitelist = results@whitelist,
                                genome.assembly = results@genome.assembly,
                                seqlevels.style = results@seqlevels.style,
                                filtering.log = results@filtering.log,
                                parameters = c(results@parameters,
                                               list(aggregateToGenes = list(method = method,
                                                                            weighting = weighting,
                                                                            decay = decay,
                                                                            linkType = linkType,
                                                                            regionSets = regionSets,
                                                                            adjustMethod = adjustMethod,
                                                                            n.genes = nrow(geneTable),
                                                                            n.links = nrow(usableLinks)))))

    return(geneResults)
  } # END function
