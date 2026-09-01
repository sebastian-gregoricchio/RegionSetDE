#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClassUnion
setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))





#' @importFrom methods setClassUnion
setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))


#' @title RegionSetDE.provenance class
#'
#' @description Virtual class collecting the filters and the parameters shared by all the RegionSetDE objects, so that the origin of the regions survives every step of the analysis. Not meant to be instantiated directly.
#'
#' @slot blacklist \code{GRanges} with the regions removed by \code{\link{applyBlacklist}}, \code{NULL} when no blacklist has been applied.
#' @slot whitelist \code{GRanges} with the regions used by \code{\link{applyWhitelist}} to restrict the sets, \code{NULL} when no whitelist has been applied.
#' @slot genome.assembly String with the genome assembly of the regions, \code{NULL} when not declared.
#' @slot seqlevels.style String with the chromosome naming style shared by all the sets.
#' @slot filtering.log Data.frame collecting the number of regions before and after each filtering step.
#' @slot parameters List with the arguments used at each step.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.provenance",
         representation = representation("VIRTUAL",
                                         blacklist = "GRangesOrNULL",
                                         whitelist = "GRangesOrNULL",
                                         genome.assembly = "characterOrNULL",
                                         seqlevels.style = "character",
                                         filtering.log = "data.frame",
                                         parameters = "list"),
         prototype = prototype(blacklist = NULL,
                               whitelist = NULL,
                               genome.assembly = NULL,
                               seqlevels.style = NA_character_,
                               filtering.log = data.frame(),
                               parameters = list()))


#' @title RegionSetDE class
#'
#' @description S4 class collecting a group of genomic region sets together with the filters applied to them. The regions are stored as a \code{GRangesList}, so that any Bioconductor operation remains available through the \code{regions} slot.
#'
#' @slot regions \code{GRangesList} containing the region sets, one element per set.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation
#'
#' @export
setClass(Class = "RegionSetDE",
         contains = "RegionSetDE.provenance",
         representation = representation(regions = "GRangesList"))


#' @title RegionSetDE.counts class
#'
#' @description S4 class storing the read counts computed over a collection of region sets. It extends \code{RangedSummarizedExperiment}, therefore \code{assay}, \code{colData}, \code{rowRanges} and the subsetting operators behave as usual, while the filters applied upstream remain accessible in the inherited provenance slots.
#'
#' @slot counting.level String indicating whether the rows correspond to the single regions (\code{"region"}) or to the whole sets (\code{"set"}).
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.counts",
         contains = c("RangedSummarizedExperiment", "RegionSetDE.provenance"),
         representation = representation(counting.level = "character"),
         prototype = prototype(counting.level = NA_character_))


#' @importFrom methods setValidity
#' @importFrom GenomeInfoDb seqlevels
setValidity(Class = "RegionSetDE",
            method = function(object) {
              issues <- character(0)

              # Every set must be addressable by name, the whole package indexes on them
              if (length(object@regions) == 0) {
                issues <- c(issues, "The 'regions' slot is empty.")
              }
              if (is.null(names(object@regions)) | any(is.na(names(object@regions))) | any(names(object@regions) == "")) {
                issues <- c(issues, "All the elements of the 'regions' slot must be named.")
              }
              if (any(duplicated(names(object@regions)))) {
                issues <- c(issues, "The names of the 'regions' slot must be unique.")
              }

              # Mixed chromosome styles return zero overlaps without raising any error
              if (length(object@regions) > 1) {
                hasChrPrefix <- vapply(as.list(object@regions),
                                       function(gr) {any(grepl("^chr", GenomeInfoDb::seqlevels(gr)))},
                                       logical(1))
                if (length(unique(hasChrPrefix)) > 1) {
                  issues <- c(issues, "The region sets mix UCSC and Ensembl/NCBI chromosome names.")
                }
              }

              if (length(issues) > 0) {return(issues)} else {return(TRUE)}
            })


#' @importFrom methods setValidity
#' @importFrom SummarizedExperiment rowData assayNames
setValidity(Class = "RegionSetDE.counts",
            method = function(object) {
              issues <- character(0)

              if (!("counts" %in% SummarizedExperiment::assayNames(object))) {
                issues <- c(issues, "A 'counts' assay is required.")
              }

              # The set membership is what links the counts back to the region sets
              if (!("region.set" %in% colnames(SummarizedExperiment::rowData(object)))) {
                issues <- c(issues, "The rowData must contain a 'region.set' column.")
              }

              if (!(object@counting.level %in% c("region", "set"))) {
                issues <- c(issues, "The 'counting.level' slot must be either 'region' or 'set'.")
              }

              if (length(issues) > 0) {return(issues)} else {return(TRUE)}
            })


#' @title show method for RegionSetDE objects
#'
#' @description Prints a summary of the region sets stored in a \code{RegionSetDE} object.
#'
#' @param object A \code{RegionSetDE} object.
#'
#' @return Prints a summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod show
#' @importFrom BiocGenerics width
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE",
          definition = function(object) {
            cat("### RegionSetDE object ###\n")
            cat(paste0("Genome assembly:   ", ifelse(is.null(object@genome.assembly), "not declared", object@genome.assembly), "\n"))
            cat(paste0("Chromosome style:  ", ifelse(is.na(object@seqlevels.style), "not declared", object@seqlevels.style), "\n"))
            cat(paste0("Region sets:       ", length(object@regions), "\n\n"))

            # Set names are padded so that the counts stay aligned
            setNames <- format(names(object@regions))
            for (i in seq_along(object@regions)) {
              cat(paste0("  ", setNames[i], "  ",
                         format(length(object@regions[[i]]), big.mark = ",", trim = TRUE), " regions  (",
                         format(sum(BiocGenerics::width(object@regions[[i]])), big.mark = ",", trim = TRUE), " bp)\n"))
            }

            cat(paste0("\nBlacklist:  ", ifelse(is.null(object@blacklist), "not applied",
                                                paste0("applied (", format(length(object@blacklist), big.mark = ",", trim = TRUE), " regions)")), "\n"))
            cat(paste0("Whitelist:  ", ifelse(is.null(object@whitelist), "not applied",
                                              paste0("applied (", format(length(object@whitelist), big.mark = ",", trim = TRUE), " regions)")), "\n"))

            if (nrow(object@filtering.log) > 0) {
              cat(paste0("\nFiltering steps: ", paste(unique(object@filtering.log$step), collapse = ", "), "\n"))
              cat("(see the 'filtering.log' slot for the details)\n")
            }
          })




##########################################################################################
###    UNIVERSE
##########################################################################################


#' @title RegionSetDE.universe class
#'
#' @description S4 class holding, for every region set, the rows that make up the universe of its competitive test: the set itself together with the rows it is compared against. The positions refer to the object the universe was built from, so the same universe serves every contrast run on that object.
#'
#' @slot index List with one element per region set, each a vector of row positions forming its universe.
#' @slot type String with the way the universe was built, one of \code{"otherSets"} and \code{"supplied"}.
#' @slot matching Character vector with the covariates the comparison rows were matched on.
#' @slot diagnostics Data.frame with, for every set, the median width and abundance of the set and of the rows it is compared against.
#' @slot n.rows Numeric value with the number of rows of the object the positions refer to.
#'
#' @details The universe answers the question "more than what". A competitive test asks whether the regions of a set respond more strongly than the rows around them, and the answer depends entirely on which rows those are, so the choice is stored rather than assumed. This has nothing to do with \code{\link{countBackground}}, which counts genome bins for the normalisation.
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{makeSetUniverse}}, \code{\link{testRegionSets}}
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.universe",
         representation = representation(index = "list",
                                         type = "character",
                                         matching = "character",
                                         diagnostics = "data.frame",
                                         n.rows = "numeric"),
         prototype = prototype(index = list(),
                               type = NA_character_,
                               matching = character(0),
                               diagnostics = data.frame(),
                               n.rows = 0))




#' @importFrom methods setValidity
setValidity(Class = "RegionSetDE.universe",
            method = function(object) {
              issues <- character(0)

              # The positions are meaningless once they point outside the object they were built on
              if (any(vapply(object@index, function(x) {any(x < 1) | any(x > object@n.rows)}, logical(1)))) {
                issues <- c(issues, "The 'index' slot points at rows outside the object.")
              }

              if (length(object@index) > 0) {
                if (is.null(names(object@index)) | any(names(object@index) == "")) {
                  issues <- c(issues, "All the elements of the 'index' slot must be named.")
                }
              }

              if (length(issues) == 0) {TRUE} else {issues}
            })




#' @title show method for RegionSetDE.universe
#'
#' @description Prints a summary of a comparison universe.
#'
#' @param object \code{RegionSetDE.universe} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.universe",
          definition = function(object) {
            cat("An object of class 'RegionSetDE.universe'\n")
            cat("  type            :", object@type, "\n")
            cat("  matched on      :", if (length(object@matching) == 0) {"nothing"} else {paste(object@matching, collapse = ", ")}, "\n")
            cat("  sets            :", length(object@index), "\n")

            if (nrow(object@diagnostics) > 0) {
              cat("\n")
              print(object@diagnostics, row.names = FALSE)
            }
            invisible(NULL)
          })


# ====================================================================================================
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ====================================================================================================


########################################################################################################
######  ANALYSES per INVIDUAL REGION
########################################################################################################


#' @title RegionSetDE.fit class
#'
#' @description S4 class storing a model fitted on a \code{RegionSetDE.counts} object. The object returned by the engine is kept as it is inside the \code{fit} slot, so that any function of \code{edgeR}, \code{limma}, \code{variancePartition} or \code{DESeq2} remains usable on it, while the design, the offsets and the samples that produced it travel with it. The counts are carried along as well: the per-set tests read the fit and never the counts directly, which is what guarantees that the two levels share the same dispersion and the same normalisation.
#'
#' @slot fit List with a single element, \code{object}, containing the fit returned by the engine, together with the auxiliary objects needed to test it.
#' @slot engine String with the engine used, one of \code{"edgeR"}, \code{"voom"}, \code{"dream"} and \code{"deseq2"}.
#' @slot design Matrix of the fixed effects, with one row per sample used.
#' @slot design.formula List holding the formula when the design was declared as such, empty otherwise.
#' @slot blocking List with the blocking variable, the random effect formula and the consensus correlation, when they apply.
#' @slot dispersion List with the dispersion estimates of the engine, summarised.
#' @slot counts \code{RegionSetDE.counts} object on which the model has been fitted, restricted to the samples used.
#' @slot universe \code{RegionSetDE.universe} object with the rows every set will be compared against at the set level. It depends on the rows and not on the contrast, so it is built once here and reused by every test run on this fit. Empty when none could be built.
#' @slot samples Character vector with the names of the samples used.
#' @slot counting.level String indicating whether the rows are regions (\code{"region"}) or tiles (\code{"tile"}).
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.fit",
         contains = "RegionSetDE.provenance",
         representation = representation(fit = "list",
                                         engine = "character",
                                         design = "matrix",
                                         design.formula = "list",
                                         blocking = "list",
                                         dispersion = "list",
                                         counts = "RegionSetDE.counts",
                                         universe = "RegionSetDE.universe",
                                         samples = "character",
                                         counting.level = "character"),
         prototype = prototype(fit = list(),
                               engine = NA_character_,
                               design = matrix(nrow = 0, ncol = 0),
                               design.formula = list(),
                               blocking = list(),
                               dispersion = list(),
                               samples = character(0),
                               counting.level = NA_character_))




#' @title RegionSetDE.results class
#'
#' @description S4 class holding the outcome of a per-region contrast. The table in the \code{results} slot has one row per region and follows the coordinates in the \code{regions} slot, so the two can be bound together without any matching step. When the counts were tiled, the tile level statistics are kept in the \code{tiles} slot next to the combined ones, since a region that changes only over part of its length is easier to read from the tiles than from the summary.
#'
#' @slot results Data.frame with one row per region.
#' @slot tiles Data.frame with one row per tile, empty when the counts were not tiled.
#' @slot regions \code{GRanges} with the coordinates of the rows of \code{results}.
#' @slot contrast String describing the contrast that was tested.
#' @slot contrast.vector Numeric vector with the coefficients of the contrast over the columns of the design.
#' @slot engine String with the engine that produced the statistics.
#' @slot counting.level String indicating whether the model was fitted on regions or on tiles.
#' @slot combination List with the method used to combine the tiles and whether it was applied.
#' @slot thresholds List with the FDR and log2 fold change cut-offs used to fill the \code{diff.status} column.
#' @slot counts \code{RegionSetDE.counts} object the contrast was computed on, so that the values behind a result can be drawn without carrying a second object around. Empty when \code{testRegions} was called with \code{carryCounts = FALSE}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#' @importClassesFrom GenomicRanges GRanges
#'
#' @export
setClass(Class = "RegionSetDE.results",
         contains = "RegionSetDE.provenance",
         representation = representation(results = "data.frame",
                                         tiles = "data.frame",
                                         regions = "GRanges",
                                         contrast = "character",
                                         contrast.vector = "numeric",
                                         engine = "character",
                                         counting.level = "character",
                                         combination = "list",
                                         thresholds = "list",
                                         counts = "RegionSetDE.counts"),
         prototype = prototype(results = data.frame(),
                               tiles = data.frame(),
                               contrast = NA_character_,
                               contrast.vector = numeric(0),
                               engine = NA_character_,
                               counting.level = NA_character_,
                               combination = list(),
                               thresholds = list()))




#' @importFrom methods setValidity
setValidity(Class = "RegionSetDE.fit",
            method = function(object) {
              issues <- character(0)

              if (!(object@engine %in% c("edgeR", "voom", "dream", "deseq2"))) {
                issues <- c(issues, "The 'engine' slot must be one of edgeR, voom, dream, deseq2.")
              }

              # The contrast resolution indexes the design by name, an unnamed column cannot be addressed
              if (ncol(object@design) > 0 & is.null(colnames(object@design))) {
                issues <- c(issues, "The columns of the 'design' slot must be named.")
              }

              if (nrow(object@design) != length(object@samples)) {
                issues <- c(issues, "The 'design' slot must have one row per sample.")
              }

              if (length(issues) == 0) {TRUE} else {issues}
            })


#' @importFrom methods setValidity
setValidity(Class = "RegionSetDE.results",
            method = function(object) {
              issues <- character(0)

              if (nrow(object@results) != length(object@regions)) {
                issues <- c(issues, "The 'results' and 'regions' slots must have the same length.")
              }

              requiredColumns <- c("region.set", "region.id", "log2FC", "p.value", "FDR")
              if (nrow(object@results) > 0 & !all(requiredColumns %in% colnames(object@results))) {
                issues <- c(issues, paste0("The 'results' slot must contain the columns: ", paste(requiredColumns, collapse = ", "), "."))
              }

              if (length(issues) == 0) {TRUE} else {issues}
            })




#' @title show method for RegionSetDE.fit
#'
#' @description Prints a summary of a fitted model.
#'
#' @param object \code{RegionSetDE.fit} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.fit",
          definition = function(object) {
            cat("An object of class 'RegionSetDE.fit'\n")
            cat("  engine          :", object@engine, "\n")
            cat("  rows            :", nrow(object@counts), paste0("(", object@counting.level, " level)"), "\n")
            cat("  samples         :", length(object@samples), paste0("(", paste(utils::head(object@samples, 4), collapse = ", "),
                                                                      if (length(object@samples) > 4) {", ..."} else {""}, ")"), "\n")
            cat("  coefficients    :", paste(colnames(object@design), collapse = ", "), "\n")

            if (length(object@blocking) > 0) {
              cat("  blocking        :", paste(names(object@blocking), collapse = ", "), "\n")
            }
            if (length(object@universe@index) > 0) {
              cat("  set universe    :", object@universe@type,
                  if (length(object@universe@matching) > 0) {paste0("(matched on ", paste(object@universe@matching, collapse = " and "), ")")} else {""},
                  "\n")
            }
            if (!is.null(object@dispersion$common)) {
              cat("  common disp.    :", signif(object@dispersion$common, 3),
                  paste0("(BCV ", signif(sqrt(object@dispersion$common), 3), ")"),
                  if (isTRUE(object@dispersion$fixed)) {"fixed"} else {""}, "\n")
            }
            if (isTRUE(object@dispersion$no.replicates)) {
              cat("  replicates      : none, dispersion from", object@dispersion$source, "\n")
            }
            invisible(NULL)
          })




#' @title show method for RegionSetDE.results
#'
#' @description Prints a summary of a per-region contrast.
#'
#' @param object \code{RegionSetDE.results} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#' @importFrom dplyr count filter
#' @importFrom rlang .data
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.results",
          definition = function(object) {
            cat("An object of class 'RegionSetDE.results'\n")
            cat("  contrast        :", object@contrast, "\n")
            cat("  engine          :", object@engine, "\n")
            cat("  regions         :", nrow(object@results), "\n")

            if (nrow(object@tiles) > 0) {
              cat("  tiles           :", nrow(object@tiles), paste0("(combined by ", object@combination$method, ")"), "\n")
            }

            if (ncol(object@counts) > 0) {
              cat("  counts carried  :", ncol(object@counts), "samples\n")
            }

            if (nrow(object@results) > 0) {
              statusCount <- dplyr::count(object@results, .data$region.set, .data$diff.status)
              upCount <- dplyr::filter(statusCount, .data$diff.status == "up")
              downCount <- dplyr::filter(statusCount, .data$diff.status == "down")

              cat("  thresholds      : FDR <", object@thresholds$FDR, "| |log2FC| >", object@thresholds$log2FC, "\n")
              cat("  changing regions:\n")

              for (setName in unique(object@results$region.set)) {
                upNumber <- sum(upCount$n[upCount$region.set == setName])
                downNumber <- sum(downCount$n[downCount$region.set == setName])
                cat("    ", setName, ": ", upNumber, " up, ", downNumber, " down\n", sep = "")
              }
            }
            invisible(NULL)
          })




##################################################################################################
###    ANALYSES by SETs
##################################################################################################


#' @title RegionSetDE.setResults class
#'
#' @description S4 class holding the outcome of a set level test. The table in the \code{results} slot has one row per region set and carries the effect size with its confidence interval next to the p-values, since with tens of thousands of regions in a set the p-value stops discriminating long before the effect size does.
#'
#' @slot results Data.frame with one row per region set, or per pair of sets for a set contrast.
#' @slot regionStats Data.frame with the per-region statistics the set level test was computed on.
#' @slot contrast String describing the contrast that was tested.
#' @slot contrast.groups List with the \code{column} of the \code{colData} the contrast separates and the two \code{groups} it compares, the first one being the level the fold change is positive for. Empty when the contrast is not a difference between two levels of one variable.
#' @slot test String with the kind of test, either \code{"set"} or \code{"setContrast"}.
#' @slot methods Character vector with the tests that were run.
#' @slot universe \code{RegionSetDE.universe} object with the rows every set was compared against. Empty when the test was self-contained only.
#' @slot engine String with the engine that produced the per-region statistics.
#' @slot thresholds List with the FDR cut-off and the correction method.
#' @slot counts \code{RegionSetDE.counts} object the test was computed on, so that the signal behind a set level claim can be drawn without carrying a second object around. Empty when \code{carryCounts} was \code{FALSE}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.setResults",
         contains = "RegionSetDE.provenance",
         representation = representation(results = "data.frame",
                                         regionStats = "data.frame",
                                         contrast = "character",
                                         contrast.groups = "list",
                                         test = "character",
                                         methods = "character",
                                         universe = "RegionSetDE.universe",
                                         engine = "character",
                                         thresholds = "list",
                                         counts = "RegionSetDE.counts"),
         prototype = prototype(results = data.frame(),
                               regionStats = data.frame(),
                               contrast = NA_character_,
                               contrast.groups = list(),
                               test = NA_character_,
                               methods = character(0),
                               engine = NA_character_,
                               thresholds = list()))




#' @title show method for RegionSetDE.setResults
#'
#' @description Prints a summary of a set level test.
#'
#' @param object \code{RegionSetDE.setResults} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.setResults",
          definition = function(object) {
            cat("An object of class 'RegionSetDE.setResults'\n")
            cat("  test            :", object@test, "\n")
            cat("  contrast        :", object@contrast, "\n")
            cat("  engine          :", object@engine, "\n")
            cat("  methods         :", paste(object@methods, collapse = ", "), "\n")

            if (!is.na(object@universe@type)) {
              cat("  universe        :", object@universe@type,
                  if (length(object@universe@matching) > 0) {paste0("(matched on ", paste(object@universe@matching, collapse = " and "), ")")} else {""},
                  "\n")
            }

            if (nrow(object@results) > 0) {
              shownColumns <- intersect(c("region.set", "set.1", "set.2", "n.regions", "mean.log2FC",
                                          "delta.log2FC", "CI.lower", "CI.upper", "camera.FDR", "fry.FDR"),
                                        colnames(object@results))
              printTable <- object@results[, shownColumns, drop = FALSE]
              numericColumns <- vapply(printTable, is.numeric, logical(1))
              printTable[numericColumns] <- lapply(printTable[numericColumns], signif, digits = 3)

              cat("\n")
              print(printTable, row.names = FALSE)
            }
            invisible(NULL)
          })




#' @title RegionSetDE.setResultsList class
#'
#' @description S4 class holding the set level outcome of several contrasts run on the same fit. Every function that takes a \code{RegionSetDE.setResults} object also takes this one, together with a \code{contrast} argument naming which of them to use.
#'
#' @slot results Named list of \code{RegionSetDE.setResults} objects.
#' @slot contrasts Character vector with the names of the contrasts, in the order they were run.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.setResultsList",
         representation = representation(results = "list",
                                         contrasts = "character"),
         prototype = prototype(results = list(),
                               contrasts = character(0)))




#' @importFrom methods setValidity is
setValidity(Class = "RegionSetDE.setResultsList",
            method = function(object) {
              issues <- character(0)

              if (length(object@results) == 0) {
                issues <- c(issues, "The \'results\' slot is empty.")
              }

              if (!all(vapply(object@results, function(x) {methods::is(x, "RegionSetDE.setResults")}, logical(1)))) {
                issues <- c(issues, "All the elements of the \'results\' slot must be RegionSetDE.setResults objects.")
              }

              if (any(duplicated(object@contrasts))) {
                issues <- c(issues, "The names of the contrasts must be unique.")
              }

              if (length(issues) == 0) {TRUE} else {issues}
            })




#' @title show method for RegionSetDE.setResultsList
#'
#' @description Prints the contrasts held by the object and the sets tested in each of them.
#'
#' @param object \code{RegionSetDE.setResultsList} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.setResultsList",
          definition = function(object) {
            cat("An object of class \'RegionSetDE.setResultsList\'\n")
            cat("  contrasts       :", length(object@results), "\n\n")

            summaryTable <-
              do.call(what = rbind,
                      args = lapply(names(object@results),
                                    function(contrastName) {
                                      setResults <- object@results[[contrastName]]
                                      significanceColumn <- if ("camera.FDR" %in% colnames(setResults@results)) {"camera.FDR"} else {"fry.FDR"}

                                      return(data.frame(name = contrastName,
                                                        contrast = setResults@contrast,
                                                        n.sets = nrow(setResults@results),
                                                        n.significant = sum(setResults@results[[significanceColumn]] < setResults@thresholds$FDR),
                                                        stringsAsFactors = FALSE))
                                    }))

            print(summaryTable, row.names = FALSE)
            invisible(NULL)
          })




#' @rdname resultsList-accessors
#' @export
setMethod(f = "[[",
          signature = "RegionSetDE.setResultsList",
          definition = function(x, i, ...) {
            return(x@results[[i]])
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "$",
          signature = "RegionSetDE.setResultsList",
          definition = function(x, name) {
            return(x@results[[name]])
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "names",
          signature = "RegionSetDE.setResultsList",
          definition = function(x) {
            return(names(x@results))
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "length",
          signature = "RegionSetDE.setResultsList",
          definition = function(x) {
            return(length(x@results))
          })






# ----------------------------------------------------------------------------------------
# MULTI CONTRAST



#' @title RegionSetDE.resultsList class
#'
#' @description S4 class holding the outcome of several contrasts run on the same fit. Since all of them come from one model, the design, the offsets and the dispersion are shared and the contrasts can be compared with each other directly. Every function that takes a \code{RegionSetDE.results} object also takes this one, together with a \code{contrast} argument naming which of them to use.
#'
#' @slot results Named list of \code{RegionSetDE.results} objects.
#' @slot contrasts Character vector with the names of the contrasts, in the order they were run.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @export
setClass(Class = "RegionSetDE.resultsList",
         representation = representation(results = "list",
                                         contrasts = "character"),
         prototype = prototype(results = list(),
                               contrasts = character(0)))




#' @importFrom methods setValidity is
setValidity(Class = "RegionSetDE.resultsList",
            method = function(object) {
              issues <- character(0)

              if (length(object@results) == 0) {
                issues <- c(issues, "The 'results' slot is empty.")
              }

              if (!all(vapply(object@results, function(x) {methods::is(x, "RegionSetDE.results")}, logical(1)))) {
                issues <- c(issues, "All the elements of the 'results' slot must be RegionSetDE.results objects.")
              }

              # The name is how a contrast is addressed downstream, a duplicate makes one of them unreachable
              if (any(duplicated(object@contrasts))) {
                issues <- c(issues, "The names of the contrasts must be unique.")
              }

              if (length(issues) == 0) {TRUE} else {issues}
            })




#' @title show method for RegionSetDE.resultsList
#'
#' @description Prints the contrasts held by the object and how many regions change in each of them.
#'
#' @param object \code{RegionSetDE.resultsList} object.
#'
#' @return Prints the summary to the console.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(f = "show",
          signature = "RegionSetDE.resultsList",
          definition = function(object) {
            cat("An object of class 'RegionSetDE.resultsList'\n")
            cat("  contrasts       :", length(object@results), "\n\n")

            summaryTable <-
              do.call(what = rbind,
                      args = lapply(names(object@results),
                                    function(contrastName) {
                                      resultTable <- object@results[[contrastName]]@results
                                      return(data.frame(name = contrastName,
                                                        contrast = object@results[[contrastName]]@contrast,
                                                        n.regions = nrow(resultTable),
                                                        up = sum(resultTable$diff.status == "up"),
                                                        down = sum(resultTable$diff.status == "down"),
                                                        stringsAsFactors = FALSE))
                                    }))

            print(summaryTable, row.names = FALSE)
            invisible(NULL)
          })




#' @title Accessors of RegionSetDE.resultsList
#'
#' @description Extract one contrast, or the names of the contrasts, from a \code{RegionSetDE.resultsList} object.
#'
#' @param x \code{RegionSetDE.resultsList} object.
#' @param i String with the name of a contrast, or its position.
#' @param j Not used, present because the \code{[[} generic carries it.
#' @param name String with the name of a contrast.
#' @param ... Not used.
#'
#' @return A \code{RegionSetDE.results} object, or a character vector for \code{names}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @rdname resultsList-accessors
#' @export
setMethod(f = "[[",
          signature = "RegionSetDE.resultsList",
          definition = function(x, i, ...) {
            return(x@results[[i]])
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "$",
          signature = "RegionSetDE.resultsList",
          definition = function(x, name) {
            return(x@results[[name]])
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "names",
          signature = "RegionSetDE.resultsList",
          definition = function(x) {
            return(names(x@results))
          })

#' @rdname resultsList-accessors
#' @export
setMethod(f = "length",
          signature = "RegionSetDE.resultsList",
          definition = function(x) {
            return(length(x@results))
          })




#' @title .pickResults
#'
#' @description Returns the single \code{RegionSetDE.results} object a function has to work on, whether it was handed one directly or one contrast out of a list.
#'
#' @param results \code{RegionSetDE.results}, \code{RegionSetDE.setResults}, or one of the two list classes holding several of them.
#' @param contrast String with the name of the contrast, or its position. Default: \code{NULL}.
#'
#' @return A \code{RegionSetDE.results} or a \code{RegionSetDE.setResults} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.pickResults <-
  function(results,
           contrast = NULL) {

    if (methods::is(results, "RegionSetDE.results") | methods::is(results, "RegionSetDE.setResults")) {
      if (!is.null(contrast)) {
        warning("The object holds a single contrast, the 'contrast' parameter has been ignored.", call. = FALSE)
      }
      return(results)
    }

    if (!methods::is(results, "RegionSetDE.resultsList") & !methods::is(results, "RegionSetDE.setResultsList")) {
      stop("The object must be a results object of the package, or a list of them.", call. = FALSE)
    }

    # One contrast in the list leaves nothing to choose, asking for a name would only be noise
    if (is.null(contrast)) {
      if (length(results@results) == 1) {
        return(results@results[[1]])
      }
      stop("The object holds several contrasts, name one with 'contrast': ",
           paste(names(results@results), collapse = ", "), ".", call. = FALSE)
    }

    if (is.numeric(contrast)) {
      contrastIndex <- as.integer(contrast[1])
      if (is.na(contrastIndex) | contrastIndex < 1 | contrastIndex > length(results@results)) {
        stop("The 'contrast' position is outside the range of the object.", call. = FALSE)
      }
      return(results@results[[contrastIndex]])
    }

    if (!(contrast[1] %in% names(results@results))) {
      stop("The contrast '", contrast[1], "' is absent from the object. Available: ",
           paste(names(results@results), collapse = ", "), ".", call. = FALSE)
    }

    return(results@results[[contrast[1]]])
  } # END function
