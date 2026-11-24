# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.


#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods setClassUnion
setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))


#' @title RegionSetDE.provenance class
#'
#' @description Virtual class collecting the filters and the parameters shared by all the RegionSetDE objects, so that the origin of the regions survives every step of the analysis. Not meant to be instantiated directly.
#'
#' @slot blacklist \code{GRanges} with the regions removed by \code{\link{applyBlacklist}}, \code{NULL} when no blacklist has been applied.
#' @slot greylist \code{GRanges} with the regions removed by \code{\link{applyGreylist}}, by the \code{greylist} argument of \code{\link{countReads}} or by \code{\link{loadConsensusPeaks}}, \code{NULL} when no greylist has been applied. Several greylists applied one after the other are merged.
#' @slot whitelist \code{GRanges} with the regions used by \code{\link{applyWhitelist}} to restrict the sets, \code{NULL} when no whitelist has been applied.
#' @slot genome.assembly String with the genome assembly of the regions, \code{NULL} when not declared.
#' @slot seqlevels.style String with the chromosome naming style shared by all the sets.
#' @slot filtering.log \code{DataFrame} collecting the number of regions before and after each filtering step. Read it with \code{\link{filteringLog}}.
#' @slot parameters List with the arguments used at each step.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' # Every object of the package carries the provenance slots
#' is(counts, "RegionSetDE.provenance")
#' genomeAssembly(counts)
#' filteringLog(counts)
#'
#' @export
setClass(
    Class = "RegionSetDE.provenance",
    representation = representation("VIRTUAL",
        blacklist = "GRangesOrNULL",
        greylist = "GRangesOrNULL",
        whitelist = "GRangesOrNULL",
        genome.assembly = "characterOrNULL",
        seqlevels.style = "character",
        filtering.log = "DataFrame",
        parameters = "list"
    ),
    prototype = prototype(
        blacklist = NULL,
        greylist = NULL,
        whitelist = NULL,
        genome.assembly = NULL,
        seqlevels.style = NA_character_,
        filtering.log = S4Vectors::DataFrame(),
        parameters = list()
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.provenance",
    method = function(object) {
        issues <- character(0)

        # The steps are found again by their name, an unnamed entry could never be read back
        if (length(object@parameters) > 0 && (is.null(names(object@parameters)) || any(names(object@parameters) == ""))) {
            issues <- c(issues, "All the elements of the 'parameters' slot must be named after the step that wrote them.")
        }

        if (length(object@seqlevels.style) != 1) {
            issues <- c(issues, "The 'seqlevels.style' slot must hold a single string.")
        }

        if (!is.null(object@genome.assembly) && length(object@genome.assembly) != 1) {
            issues <- c(issues, "The 'genome.assembly' slot must hold a single string, or be NULL.")
        }

        # Every filter writes its rows under the name of its step
        if (nrow(object@filtering.log) > 0 && !("step" %in% colnames(object@filtering.log))) {
            issues <- c(issues, "The 'filtering.log' slot must have a 'step' column.")
        }

        if (length(issues) > 0) {
            return(issues)
        } else {
            return(TRUE)
        }
    }
)


#' @title RegionSetDE class
#'
#' @description S4 class collecting a group of genomic region sets together with the filters applied to them. The regions are stored as a \code{GRangesList}, so that any Bioconductor operation remains available through the \code{regions} slot.
#'
#' @slot regions \code{GRangesList} containing the region sets, one element per set.
#' @slot consensus List with the consensus data when the regions were built from peaks by \code{\link{loadConsensusPeaks}}: the consensus of every group, the peaks of every sample after the exclusion, the total consensus, the table of the samples and the sample sheet. Empty for regions loaded any other way. Read it with \code{\link{consensusData}}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regions <- loadRegions(split(regionTable[, c("seqnames", "start", "end")], regionTable$setName),
#'     genomeAssembly = "rn4", verbose = FALSE
#' )
#' regions
#' regionSetNames(regions)
#'
#' @export
setClass(
    Class = "RegionSetDE",
    contains = "RegionSetDE.provenance",
    representation = representation(
        regions = "GRangesList",
        consensus = "list"
    ),
    prototype = prototype(consensus = list())
)


#' @title RegionSetDE.counts class
#'
#' @description S4 class storing the read counts computed over a collection of region sets. It extends \code{RangedSummarizedExperiment}, therefore \code{assay}, \code{colData}, \code{rowRanges} and the subsetting operators behave as usual, while the filters applied upstream remain accessible in the inherited provenance slots.
#'
#' @slot counting.level String indicating whether the rows are regions (\code{"region"}) or tiles of a region (\code{"tile"}). It is the same vocabulary the \code{RegionSetDE.fit} and \code{RegionSetDE.results} classes use, since the value travels from here into both of them.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts
#' head(sampleInfo(counts))
#'
#' @export
setClass(
    Class = "RegionSetDE.counts",
    contains = c("RangedSummarizedExperiment", "RegionSetDE.provenance"),
    representation = representation(counting.level = "character"),
    prototype = prototype(counting.level = NA_character_)
)


#' @importFrom methods setValidity
#' @importFrom GenomeInfoDb seqlevels
setValidity(
    Class = "RegionSetDE",
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

        # The consensus is read by name downstream, a list without these parts is not one
        if (length(object@consensus) > 0) {
            missingParts <- setdiff(c("groups", "total", "samples", "mode"), names(object@consensus))
            if (length(missingParts) > 0) {
                issues <- c(issues, paste("The 'consensus' slot lacks:", paste(missingParts, collapse = ", ")))
            }
        }

        # Mixed chromosome styles return zero overlaps without raising any error
        if (length(object@regions) > 1) {
            hasChrPrefix <- vapply(
                as.list(object@regions),
                function(gr) {
                    any(grepl("^chr", GenomeInfoDb::seqlevels(gr)))
                },
                logical(1)
            )
            if (length(unique(hasChrPrefix)) > 1) {
                issues <- c(issues, "The region sets mix UCSC and Ensembl/NCBI chromosome names.")
            }
        }

        if (length(issues) > 0) {
            return(issues)
        } else {
            return(TRUE)
        }
    }
)


#' @importFrom methods setValidity
#' @importFrom SummarizedExperiment rowData assayNames
setValidity(
    Class = "RegionSetDE.counts",
    method = function(object) {
        issues <- character(0)

        if (!("counts" %in% SummarizedExperiment::assayNames(object))) {
            issues <- c(issues, "A 'counts' assay is required.")
        }

        # The set membership is what links the counts back to the region sets
        if (!("region.set" %in% colnames(SummarizedExperiment::rowData(object)))) {
            issues <- c(issues, "The rowData must contain a 'region.set' column.")
        }

        if (!(object@counting.level %in% c("region", "tile"))) {
            issues <- c(issues, "The 'counting.level' slot must be either 'region' or 'tile'.")
        }

        if (length(issues) > 0) {
            return(issues)
        } else {
            return(TRUE)
        }
    }
)


#' @rdname RegionSetDE-class
#'
#' @param object A \code{RegionSetDE} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod show
#' @importFrom BiocGenerics width
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE",
    definition = function(object) {
        cat("### RegionSetDE object ###\n")
        cat(paste0("Genome assembly:   ", ifelse(is.null(object@genome.assembly), "not declared", object@genome.assembly), "\n"))
        cat(paste0("Chromosome style:  ", ifelse(is.na(object@seqlevels.style), "not declared", object@seqlevels.style), "\n"))
        cat(paste0("Region sets:       ", length(object@regions), "\n\n"))

        # Set names are padded so that the counts stay aligned
        setNames <- format(names(object@regions))
        for (i in seq_along(object@regions)) {
            cat(paste0(
                "  ", setNames[i], "  ",
                format(length(object@regions[[i]]), big.mark = ",", trim = TRUE), " regions  (",
                format(sum(BiocGenerics::width(object@regions[[i]])), big.mark = ",", trim = TRUE), " bp)\n"
            ))
        }

        cat(paste0("\nBlacklist:  ", ifelse(is.null(object@blacklist), "not applied",
            paste0("applied (", format(length(object@blacklist), big.mark = ",", trim = TRUE), " regions)")
        ), "\n"))
        cat(paste0("Greylist:   ", ifelse(is.null(object@greylist), "not applied",
            paste0("applied (", format(length(object@greylist), big.mark = ",", trim = TRUE), " regions)")
        ), "\n"))
        cat(paste0("Whitelist:  ", ifelse(is.null(object@whitelist), "not applied",
            paste0("applied (", format(length(object@whitelist), big.mark = ",", trim = TRUE), " regions)")
        ), "\n"))

        if (nrow(object@filtering.log) > 0) {
            cat(paste0("\nFiltering steps: ", paste(unique(object@filtering.log$step), collapse = ", "), "\n"))
            cat("(see filteringLog() for the details)\n")
        }

        # Regions built from peaks say so, together with the groups behind them
        if (length(object@consensus) > 0) {
            cat(paste0(
                "\nConsensus:  ", nrow(object@consensus$samples), " samples in ", length(object@consensus$groups), " groups (",
                object@consensus$mode, " mode), ",
                format(length(object@consensus$total), big.mark = ",", trim = TRUE), " regions in the total consensus\n"
            ))
            cat("(see consensusData() for the details)\n")
        }
    }
)


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
#' @slot comparison.sets Character vector with the names of the region sets the comparison rows were drawn from. A competitive p-value is relative to these and to nothing else, so they are stored beside the result rather than left to be reconstructed.
#' @slot diagnostics \code{DataFrame} with, for every set, the median width and abundance of the set and of the rows it is compared against.
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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' universe <- makeSetUniverse(fit, verbose = FALSE)
#' universe
#'
#' @export
setClass(
    Class = "RegionSetDE.universe",
    representation = representation(
        index = "list",
        type = "character",
        matching = "character",
        comparison.sets = "character",
        diagnostics = "DataFrame",
        n.rows = "numeric"
    ),
    prototype = prototype(
        index = list(),
        type = NA_character_,
        matching = character(0),
        comparison.sets = character(0),
        diagnostics = S4Vectors::DataFrame(),
        n.rows = 0
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.universe",
    method = function(object) {
        issues <- character(0)

        # The positions are meaningless once they point outside the object they were built on
        if (any(vapply(object@index, function(x) {
            any(x < 1) | any(x > object@n.rows)
        }, logical(1)))) {
            issues <- c(issues, "The 'index' slot points at rows outside the object.")
        }

        if (length(object@index) > 0) {
            if (is.null(names(object@index)) | any(names(object@index) == "")) {
                issues <- c(issues, "All the elements of the 'index' slot must be named.")
            }
        }

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @rdname RegionSetDE.universe-class
#'
#' @param object \code{RegionSetDE.universe} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.universe",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.universe'\n")
        cat("  type            :", object@type, "\n")
        cat("  matched on      :", if (length(object@matching) == 0) {
            "nothing"
        } else {
            paste(object@matching, collapse = ", ")
        }, "\n")
        cat("  compared against:", if (length(object@comparison.sets) == 0) {
            "not recorded"
        } else {
            paste(object@comparison.sets, collapse = ", ")
        }, "\n")
        cat("  sets            :", length(object@index), "\n")

        if (nrow(object@diagnostics) > 0) {
            cat("\n")
            print(as.data.frame(object@diagnostics, optional = TRUE), row.names = FALSE)
        }
        invisible(NULL)
    }
)


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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' fit
#' fitEngine(fit)
#'
#' @export
setClass(
    Class = "RegionSetDE.fit",
    contains = "RegionSetDE.provenance",
    representation = representation(
        fit = "list",
        engine = "character",
        design = "matrix",
        design.formula = "list",
        blocking = "list",
        dispersion = "list",
        counts = "RegionSetDE.counts",
        universe = "RegionSetDE.universe",
        samples = "character",
        counting.level = "character"
    ),
    prototype = prototype(
        fit = list(),
        engine = NA_character_,
        design = matrix(nrow = 0, ncol = 0),
        design.formula = list(),
        blocking = list(),
        dispersion = list(),
        samples = character(0),
        counting.level = NA_character_
    )
)


#' @title RegionSetDE.results class
#'
#' @description S4 class holding the outcome of a per-region contrast. The table in the \code{results} slot has one row per region and follows the coordinates in the \code{regions} slot, so the two can be bound together without any matching step. When the counts were tiled, the tile level statistics are kept in the \code{tiles} slot next to the combined ones, since a region that changes only over part of its length is easier to read from the tiles than from the summary.
#'
#' @slot results \code{DataFrame} with one row per region. Read it with \code{\link{resultsTable}}.
#' @slot tiles \code{DataFrame} with one row per tile, empty when the counts were not tiled. Read it with \code{\link{tileTable}}.
#' @slot regions \code{GRanges} with the coordinates of the rows of \code{results}.
#' @slot contrast String describing the contrast that was tested.
#' @slot contrast.vector Numeric vector with the coefficients of the contrast over the columns of the design.
#' @slot contrast.groups List with the \code{column} of the \code{colData} the contrast separates and the two \code{groups} it compares, the first one being the level the fold change is positive for, and \code{n.samples}, the number of samples in each of the two. Empty when the contrast is not a difference between two levels of one variable. \code{\link{plotRegion}} reads it to draw the model brackets, and \code{\link{contrastInfo}} reports it.
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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' results
#' head(resultsTable(results), 3)
#'
#' @export
setClass(
    Class = "RegionSetDE.results",
    contains = "RegionSetDE.provenance",
    representation = representation(
        results = "DataFrame",
        tiles = "DataFrame",
        regions = "GRanges",
        contrast = "character",
        contrast.vector = "numeric",
        contrast.groups = "list",
        engine = "character",
        counting.level = "character",
        combination = "list",
        thresholds = "list",
        counts = "RegionSetDE.counts"
    ),
    prototype = prototype(
        results = S4Vectors::DataFrame(),
        tiles = S4Vectors::DataFrame(),
        contrast = NA_character_,
        contrast.vector = numeric(0),
        contrast.groups = list(),
        engine = NA_character_,
        counting.level = NA_character_,
        combination = list(),
        thresholds = list()
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.fit",
    method = function(object) {
        issues <- character(0)

        if (length(object@engine) != 1 || !(object@engine %in% c("edgeR", "voom", "dream", "deseq2", "limma"))) {
            issues <- c(issues, "The 'engine' slot must be one of edgeR, voom, dream, deseq2, limma.")
        }

        # Every test reads the model of the engine from the 'object' element
        if (length(object@fit) > 0 && !("object" %in% names(object@fit))) {
            issues <- c(issues, "The 'fit' slot must hold the model of the engine in its 'object' element.")
        }

        if (length(object@design.formula) > 0 && !methods::is(object@design.formula$formula, "formula")) {
            issues <- c(issues, "The 'design.formula' slot must be empty, or hold a formula in its 'formula' element.")
        }

        for (slotName in c("blocking", "dispersion")) {
            slotValue <- methods::slot(object, slotName)
            if (length(slotValue) > 0 && (is.null(names(slotValue)) || any(names(slotValue) == ""))) {
                issues <- c(issues, paste0("All the elements of the '", slotName, "' slot must be named."))
            }
        }

        # The contrast resolution indexes the design by name, an unnamed column cannot be addressed
        if (ncol(object@design) > 0 & is.null(colnames(object@design))) {
            issues <- c(issues, "The columns of the 'design' slot must be named.")
        }

        if (nrow(object@design) != length(object@samples)) {
            issues <- c(issues, "The 'design' slot must have one row per sample.")
        }

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.results",
    method = function(object) {
        issues <- character(0)

        if (nrow(object@results) != length(object@regions)) {
            issues <- c(issues, "The 'results' and 'regions' slots must have the same length.")
        }

        requiredColumns <- c("region.set", "region.id", "log2FC", "p.value", "FDR")
        if (nrow(object@results) > 0 & !all(requiredColumns %in% colnames(object@results))) {
            issues <- c(issues, paste0("The 'results' slot must contain the columns: ", paste(requiredColumns, collapse = ", "), "."))
        }

        issues <- c(issues, .checkContrastLists(object))

        # The tile level table only exists together with the way its tiles were combined
        if (nrow(object@tiles) > 0 && !("method" %in% names(object@combination))) {
            issues <- c(issues, "The 'combination' slot must name the 'method' that combined the tiles.")
        }

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @rdname RegionSetDE.fit-class
#'
#' @param object \code{RegionSetDE.fit} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.fit'\n")
        cat("  engine          :", object@engine, "\n")
        cat("  rows            :", nrow(object@counts), paste0("(", object@counting.level, " level)"), "\n")
        cat("  samples         :", length(object@samples), paste0(
            "(", paste(utils::head(object@samples, 4), collapse = ", "),
            if (length(object@samples) > 4) {
                ", ..."
            } else {
                ""
            }, ")"
        ), "\n")
        cat("  coefficients    :", paste(colnames(object@design), collapse = ", "), "\n")

        if (length(object@blocking) > 0) {
            cat("  blocking        :", paste(names(object@blocking), collapse = ", "), "\n")
        }
        if (length(object@universe@index) > 0) {
            cat(
                "  set universe    :", object@universe@type,
                if (length(object@universe@matching) > 0) {
                    paste0("(matched on ", paste(object@universe@matching, collapse = " and "), ")")
                } else {
                    ""
                },
                "\n"
            )
        }
        if (!is.null(object@dispersion$common)) {
            cat(
                "  common disp.    :", signif(object@dispersion$common, 3),
                paste0("(BCV ", signif(sqrt(object@dispersion$common), 3), ")"),
                if (isTRUE(object@dispersion$fixed)) {
                    "fixed"
                } else {
                    ""
                }, "\n"
            )
        }
        if (isTRUE(object@dispersion$no.replicates)) {
            cat("  replicates      : none, dispersion from", object@dispersion$source, "\n")
        }
        invisible(NULL)
    }
)


#' @rdname RegionSetDE.results-class
#'
#' @param object \code{RegionSetDE.results} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr count filter
#' @importFrom rlang .data
#'
#' @export
setMethod(
    f = "show",
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
            statusCount <- dplyr::count(as.data.frame(object@results, optional = TRUE), .data$region.set, .data$diff.status)
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
    }
)


##################################################################################################
###    ANALYSES by SETs
##################################################################################################


#' @title RegionSetDE.setResults class
#'
#' @description S4 class holding the outcome of a set level test. The table in the \code{results} slot has one row per region set and carries the effect size with its confidence interval next to the p-values, since with tens of thousands of regions in a set the p-value stops discriminating long before the effect size does.
#'
#' @slot results \code{DataFrame} with one row per region set, or per pair of sets for a set contrast. Read it with \code{\link{resultsTable}}.
#' @slot regionStats \code{DataFrame} with the per-region statistics the set level test was computed on.
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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' setResults <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' setResults
#'
#' @export
setClass(
    Class = "RegionSetDE.setResults",
    contains = "RegionSetDE.provenance",
    representation = representation(
        results = "DataFrame",
        regionStats = "DataFrame",
        contrast = "character",
        contrast.groups = "list",
        test = "character",
        methods = "character",
        universe = "RegionSetDE.universe",
        engine = "character",
        thresholds = "list",
        counts = "RegionSetDE.counts"
    ),
    prototype = prototype(
        results = S4Vectors::DataFrame(),
        regionStats = S4Vectors::DataFrame(),
        contrast = NA_character_,
        contrast.groups = list(),
        test = NA_character_,
        methods = character(0),
        engine = NA_character_,
        thresholds = list()
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.setResults",
    method = function(object) {
        issues <- character(0)

        if (length(object@test) != 1 || !(is.na(object@test) | object@test %in% c("set", "setContrast"))) {
            issues <- c(issues, "The 'test' slot must be either 'set' or 'setContrast'.")
        }

        issues <- c(issues, .checkContrastLists(object))

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @rdname RegionSetDE.setResults-class
#'
#' @param object \code{RegionSetDE.setResults} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.setResults",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.setResults'\n")
        cat("  test            :", object@test, "\n")
        cat("  contrast        :", object@contrast, "\n")
        cat("  engine          :", object@engine, "\n")
        cat("  methods         :", paste(object@methods, collapse = ", "), "\n")

        if (!is.na(object@universe@type)) {
            cat(
                "  universe        :", object@universe@type,
                if (length(object@universe@matching) > 0) {
                    paste0("(matched on ", paste(object@universe@matching, collapse = " and "), ")")
                } else {
                    ""
                },
                "\n"
            )
        }

        if (nrow(object@results) > 0) {
            shownColumns <- intersect(
                c(
                    "region.set", "set.1", "set.2", "n.regions", "mean.log2FC",
                    "delta.log2FC", "CI.lower", "CI.upper", "CI.type",
                    "camera.FDR", "camera.mixed.FDR", "fry.FDR", "inter.region.cor.universe"
                ),
                colnames(object@results)
            )
            printTable <- as.data.frame(object@results[, shownColumns, drop = FALSE], optional = TRUE)
            numericColumns <- vapply(printTable, is.numeric, logical(1))
            printTable[numericColumns] <- lapply(printTable[numericColumns], signif, digits = 3)

            cat("\n")
            print(printTable, row.names = FALSE)
        }
        invisible(NULL)
    }
)


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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # A named list of contrasts gives one set level result per contrast
#' setResultsList <- testRegionSets(fit,
#'     contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
#'     verbose = FALSE
#' )
#' setResultsList
#' names(setResultsList)
#'
#' @export
setClass(
    Class = "RegionSetDE.setResultsList",
    representation = representation(
        results = "list",
        contrasts = "character"
    ),
    prototype = prototype(
        results = list(),
        contrasts = character(0)
    )
)


#' @importFrom methods setValidity is
setValidity(
    Class = "RegionSetDE.setResultsList",
    method = function(object) {
        issues <- character(0)

        if (length(object@results) == 0) {
            issues <- c(issues, "The \'results\' slot is empty.")
        }

        if (!all(vapply(object@results, function(x) {
            methods::is(x, "RegionSetDE.setResults")
        }, logical(1)))) {
            issues <- c(issues, "All the elements of the \'results\' slot must be RegionSetDE.setResults objects.")
        }

        if (any(duplicated(object@contrasts))) {
            issues <- c(issues, "The names of the contrasts must be unique.")
        }

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @rdname RegionSetDE.setResultsList-class
#'
#' @param object \code{RegionSetDE.setResultsList} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.setResultsList",
    definition = function(object) {
        cat("An object of class \'RegionSetDE.setResultsList\'\n")
        cat("  contrasts       :", length(object@results), "\n\n")

        summaryTable <-
            do.call(
                what = rbind,
                args = lapply(
                    names(object@results),
                    function(contrastName) {
                        setResults <- object@results[[contrastName]]
                        significanceColumn <- if ("camera.FDR" %in% colnames(setResults@results)) {
                            "camera.FDR"
                        } else {
                            "fry.FDR"
                        }

                        return(data.frame(
                            name = contrastName,
                            contrast = setResults@contrast,
                            n.sets = nrow(setResults@results),
                            n.significant = sum(setResults@results[[significanceColumn]] < setResults@thresholds$FDR),
                            stringsAsFactors = FALSE
                        ))
                    }
                )
            )

        print(summaryTable, row.names = FALSE)
        invisible(NULL)
    }
)


#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "[[",
    signature = "RegionSetDE.setResultsList",
    definition = function(x, i, ...) {
        return(x@results[[i]])
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "$",
    signature = "RegionSetDE.setResultsList",
    definition = function(x, name) {
        return(x@results[[name]])
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "names",
    signature = "RegionSetDE.setResultsList",
    definition = function(x) {
        return(names(x@results))
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "length",
    signature = "RegionSetDE.setResultsList",
    definition = function(x) {
        return(length(x@results))
    }
)


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
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # A named list of contrasts gives one result per contrast
#' resultsList <- testRegions(fit,
#'     contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
#'     verbose = FALSE
#' )
#' resultsList
#' resultsList[["shr"]]
#'
#' @export
setClass(
    Class = "RegionSetDE.resultsList",
    representation = representation(
        results = "list",
        contrasts = "character"
    ),
    prototype = prototype(
        results = list(),
        contrasts = character(0)
    )
)


#' @importFrom methods setValidity is
setValidity(
    Class = "RegionSetDE.resultsList",
    method = function(object) {
        issues <- character(0)

        if (length(object@results) == 0) {
            issues <- c(issues, "The 'results' slot is empty.")
        }

        if (!all(vapply(object@results, function(x) {
            methods::is(x, "RegionSetDE.results")
        }, logical(1)))) {
            issues <- c(issues, "All the elements of the 'results' slot must be RegionSetDE.results objects.")
        }

        # The name is how a contrast is addressed downstream, a duplicate makes one of them unreachable
        if (any(duplicated(object@contrasts))) {
            issues <- c(issues, "The names of the contrasts must be unique.")
        }

        if (length(issues) == 0) {
            TRUE
        } else {
            issues
        }
    }
)


#' @rdname RegionSetDE.resultsList-class
#'
#' @param object \code{RegionSetDE.resultsList} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.resultsList",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.resultsList'\n")
        cat("  contrasts       :", length(object@results), "\n\n")

        summaryTable <-
            do.call(
                what = rbind,
                args = lapply(
                    names(object@results),
                    function(contrastName) {
                        resultTable <- object@results[[contrastName]]@results
                        return(data.frame(
                            name = contrastName,
                            contrast = object@results[[contrastName]]@contrast,
                            n.regions = nrow(resultTable),
                            up = sum(resultTable$diff.status == "up"),
                            down = sum(resultTable$diff.status == "down"),
                            stringsAsFactors = FALSE
                        ))
                    }
                )
            )

        print(summaryTable, row.names = FALSE)
        invisible(NULL)
    }
)


#' @title Accessors of RegionSetDE.resultsList
#'
#' @description Extract one contrast, all of them, or the names of the contrasts, from a \code{RegionSetDE.resultsList} or a \code{RegionSetDE.setResultsList} object.
#'
#' @param x \code{RegionSetDE.resultsList} or \code{RegionSetDE.setResultsList} object.
#' @param i String with the name of a contrast, or its position.
#' @param j Not used, present because the \code{[[} generic carries it.
#' @param name String with the name of a contrast.
#' @param ... Not used.
#'
#' @return A \code{RegionSetDE.results} object, or a \code{RegionSetDE.setResults} object for the list of set level results. \code{names} returns a character vector, \code{length} the number of contrasts, and \code{as.list} a named list with every result, one per contrast.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setMethod
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' resultsList <- testRegions(fit,
#'     contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
#'     verbose = FALSE
#' )
#' names(resultsList)
#' length(resultsList)
#' resultsList$bn
#'
#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "[[",
    signature = "RegionSetDE.resultsList",
    definition = function(x, i, ...) {
        return(x@results[[i]])
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "$",
    signature = "RegionSetDE.resultsList",
    definition = function(x, name) {
        return(x@results[[name]])
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "names",
    signature = "RegionSetDE.resultsList",
    definition = function(x) {
        return(names(x@results))
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "length",
    signature = "RegionSetDE.resultsList",
    definition = function(x) {
        return(length(x@results))
    }
)


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
#' @noRd

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
            if (length(results) == 1) {
                return(results[[1]])
            }
            stop("The object holds several contrasts, name one with 'contrast': ",
                paste(names(results), collapse = ", "), ".",
                call. = FALSE
            )
        }

        if (is.numeric(contrast)) {
            contrastIndex <- as.integer(contrast[1])
            if (is.na(contrastIndex) | contrastIndex < 1 | contrastIndex > length(results)) {
                stop("The 'contrast' position is outside the range of the object.", call. = FALSE)
            }
            return(results[[contrastIndex]])
        }

        if (!(contrast[1] %in% names(results))) {
            stop("The contrast '", contrast[1], "' is absent from the object. Available: ",
                paste(names(results), collapse = ", "), ".",
                call. = FALSE
            )
        }

        return(results[[contrast[1]]])
    } # END function


#' @title RegionSetDE.setScores class
#'
#' @description S4 class storing the signal scores computed by \code{\link{scoreRegionSets}}, one per region set per library, together with the paired comparisons between the sets. No contrast lies behind these numbers and no model was fitted on them, which is what separates this class from \code{RegionSetDE.setResults}: the replication is the libraries, and the comparison runs between sets rather than between conditions.
#'
#' @slot scores \code{DataFrame} with one row per library per set, carrying the summarised signal, the reference it was divided by and the \code{score} on a log2 scale.
#' @slot comparisons \code{DataFrame} with one row per pair of sets, carrying the mean paired difference between their scores, its confidence interval, the paired t-statistic and the adjusted p-value.
#' @slot sample.metadata \code{DataFrame} with the \code{colData} of the counts object the scores came from, so that the libraries can be grouped and coloured downstream without handing that object over again.
#' @slot reference String naming what the signal of a set was divided by, one among \code{"background"}, \code{"regions"} and \code{"none"}.
#' @slot assay String with the name of the assay the signal was read from.
#' @slot summary String with how the regions of a set were summarised, either \code{"mean"} or \code{"median"}.
#' @slot per.basepair Logical value indicating whether the signal was divided by the width of the region before being summarised.
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{scoreRegionSets}}
#'
#' @importFrom methods setClass representation prototype
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' setScores <- scoreRegionSets(counts, verbose = FALSE)
#' setScores
#'
#' @export
setClass(
    Class = "RegionSetDE.setScores",
    contains = "RegionSetDE.provenance",
    representation = representation(
        scores = "DataFrame",
        comparisons = "DataFrame",
        sample.metadata = "DataFrame",
        reference = "character",
        assay = "character",
        summary = "character",
        per.basepair = "logical"
    ),
    prototype = prototype(
        scores = S4Vectors::DataFrame(),
        comparisons = S4Vectors::DataFrame(),
        sample.metadata = S4Vectors::DataFrame(),
        reference = NA_character_,
        assay = NA_character_,
        summary = NA_character_,
        per.basepair = NA
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.setScores",
    method = function(object) {
        issues <- character(0)

        requiredScoreColumns <- c("sample", "region.set", "n.regions", "set.signal", "reference.signal", "score")
        missingScoreColumns <- setdiff(requiredScoreColumns, colnames(object@scores))

        if (nrow(object@scores) > 0 & length(missingScoreColumns) > 0) {
            issues <- c(issues, paste0("The 'scores' slot lacks the columns: ", paste(missingScoreColumns, collapse = ", "), "."))
        }

        requiredComparisonColumns <- c("set.1", "set.2", "n.libraries", "mean.delta.score", "p.value", "FDR")
        missingComparisonColumns <- setdiff(requiredComparisonColumns, colnames(object@comparisons))

        if (nrow(object@comparisons) > 0 & length(missingComparisonColumns) > 0) {
            issues <- c(issues, paste0("The 'comparisons' slot lacks the columns: ", paste(missingComparisonColumns, collapse = ", "), "."))
        }

        # A comparison naming a set that was never scored has nothing behind it to read or to draw
        if (length(issues) == 0 & nrow(object@comparisons) > 0) {
            unknownSets <- setdiff(c(object@comparisons$set.1, object@comparisons$set.2), object@scores$region.set)

            if (length(unknownSets) > 0) {
                issues <- c(issues, paste0("The 'comparisons' slot names sets absent from 'scores': ", paste(unknownSets, collapse = ", "), "."))
            }
        }

        if (length(issues) > 0) {
            return(issues)
        } else {
            return(TRUE)
        }
    }
)


#' @rdname RegionSetDE.setScores-class
#'
#' @param object \code{RegionSetDE.setScores} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' scoreRegionSets(counts, verbose = FALSE)
#'
#' @importFrom methods setMethod
#' @importFrom dplyr n_distinct
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.setScores",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.setScores'\n")
        cat("  region sets     :", paste(sort(unique(object@scores$region.set)), collapse = ", "), "\n")
        cat("  libraries       :", dplyr::n_distinct(object@scores$sample), "\n")
        cat("  reference       :", object@reference, "\n")
        cat("  signal          :", paste0(
            object@summary, " of '", object@assay, "'",
            if (isTRUE(object@per.basepair)) {
                " per base pair"
            } else {
                " per region"
            }
        ), "\n")
        cat("  comparisons     :", nrow(object@comparisons), "\n")

        if (nrow(object@comparisons) > 0) {
            shownColumns <- intersect(
                c(
                    "set.1", "set.2", "n.libraries", "mean.delta.score",
                    "CI.lower", "CI.upper", "p.value", "FDR"
                ),
                colnames(object@comparisons)
            )
            printTable <- as.data.frame(object@comparisons[, shownColumns, drop = FALSE], optional = TRUE)
            numericColumns <- vapply(printTable, is.numeric, logical(1))
            printTable[numericColumns] <- lapply(printTable[numericColumns], signif, digits = 3)

            cat("\n")
            print(printTable, row.names = FALSE)
        }

        cat("\nThe libraries are the replication, the composition of the sets is not controlled for.\n")

        invisible(NULL)
    }
)


##########################################################################################
###    UPDATE OF STORED OBJECTS
##########################################################################################


#' @title .updateProvenance
#'
#' @description Brings an object saved by an earlier version of the package to the current definition of its class. An object saved without the \code{greylist} slot gets the greylist its consensus was cleaned against, when there is one, and \code{NULL} otherwise. The tables saved as data.frames are turned into the \code{DataFrame} their slots now hold. The counts and the universe carried by a fit or a result are updated together with it.
#'
#' @param object An object extending \code{RegionSetDE.provenance}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{FALSE}.
#'
#' @return The object with the slots it was missing.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods .hasSlot slot<- is
#' @importFrom BiocGenerics updateObject
#'
#' @noRd
#' @noRd

.updateProvenance <-
    function(object,
             verbose = FALSE) {
        # The slot is written directly: the object does not have it yet, and the usual assignment would refuse it
        if (!methods::.hasSlot(object, "greylist")) {
            storedGreylist <- NULL

            if (methods::.hasSlot(object, "consensus") && length(object@consensus) > 0) {
                storedGreylist <- object@consensus$greylist
            }

            methods::slot(object, "greylist", check = FALSE) <- storedGreylist

            if (isTRUE(verbose)) {
                message(
                    "[updateObject] The '", class(object)[1], "' object gained the 'greylist' slot",
                    if (is.null(storedGreylist)) {
                        ", left empty."
                    } else {
                        ", filled with the greylist of its consensus."
                    }
                )
            }
        }

        object <- .updateTableSlots(object = object, verbose = verbose)

        # A fit or a result carries the counts it was computed on, and its universe, saved together with it
        if (methods::.hasSlot(object, "counts") && methods::is(object@counts, "RegionSetDE.counts")) {
            object@counts <- BiocGenerics::updateObject(object@counts, verbose = verbose)
        }

        if (methods::.hasSlot(object, "universe") && methods::is(object@universe, "RegionSetDE.universe")) {
            object@universe <- BiocGenerics::updateObject(object@universe, verbose = verbose)
        }

        return(object)
    } # END function


#' @title .updateTableSlots
#'
#' @description Turns the tables an earlier version of the package saved as data.frames into the \code{DataFrame} their slots hold now.
#'
#' @param object An object of one of the classes of the package.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{FALSE}.
#'
#' @return The object, with every table slot holding a \code{DataFrame}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods getSlots .hasSlot slot slot<- as
#'
#' @noRd
#' @noRd

.updateTableSlots <-
    function(object,
             verbose = FALSE) {
        slotClasses <- methods::getSlots(class(object))
        tableSlots <- names(slotClasses)[slotClasses == "DataFrame"]
        convertedSlots <- character(0)

        for (slotName in tableSlots) {
            # The slot is read and written without checks, its content does not match the definition yet
            if (methods::.hasSlot(object, slotName)) {
                slotValue <- methods::slot(object, slotName)
                if (is.data.frame(slotValue)) {
                    methods::slot(object, slotName, check = FALSE) <- methods::as(slotValue, "DataFrame")
                    convertedSlots <- c(convertedSlots, slotName)
                }
            }
        }

        if (isTRUE(verbose) & length(convertedSlots) > 0) {
            message(
                "[updateObject] The '", class(object)[1], "' object now holds its tables as DataFrame: ",
                paste(convertedSlots, collapse = ", "), "."
            )
        }

        return(object)
    } # END function


#' @title updateObject methods for the RegionSetDE classes
#'
#' @description Bring an object saved by an earlier version of RegionSetDE to the current definition of its class, as \code{BiocGenerics::updateObject} does for the Bioconductor classes. An object saved before the \code{greylist} slot existed gets it, filled with the greylist of the consensus when the regions were built by \code{\link{loadConsensusPeaks}}, and empty otherwise. Tables saved as data.frames become the \code{S4Vectors::DataFrame} their slots hold now. A fit or a result updates the counts it carries as well, and a list of results every result in it.
#'
#' @param object An object of one of the RegionSetDE classes.
#' @param ... Passed on to the method of \code{SummarizedExperiment} for a counts object, ignored otherwise.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{FALSE}.
#'
#' @return The object, updated to the current definition of its class. An object already up to date comes back as it was.
#'
#' @details A greylist applied with \code{\link{applyGreylist}} before the slot existed left only its size in the object, so the update cannot recover it: apply it again, or pass it to the \code{greylist} argument of \code{\link{countReads}}.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- BiocGenerics::updateObject(counts, verbose = TRUE)
#'
#' @author Sebastian Gregoricchio
#'
#' @rdname updateObject-RegionSetDE
#'
#' @importFrom BiocGenerics updateObject
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod callNextMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.counts",
    definition = function(object, ..., verbose = FALSE) {
        object <- .updateProvenance(object = object, verbose = verbose)

        # The SummarizedExperiment part has an update of its own, run once the slot is in place
        return(methods::callNextMethod(object, ..., verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.fit",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.results",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.setResults",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.setScores",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.universe",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateTableSlots(object = object, verbose = verbose))
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#' @importFrom BiocGenerics updateObject
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.resultsList",
    definition = function(object, ..., verbose = FALSE) {
        object@results <- lapply(object@results, BiocGenerics::updateObject, verbose = verbose)
        return(object)
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#' @importFrom BiocGenerics updateObject
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.setResultsList",
    definition = function(object, ..., verbose = FALSE) {
        object@results <- lapply(object@results, BiocGenerics::updateObject, verbose = verbose)
        return(object)
    }
)


#' @title RegionSetDE.geneResults class
#'
#' @description S4 class holding the per-gene table built by \code{\link{aggregateToGenes}} from a per-region result. One row per gene, with the regions that reached it counted and their directions reported rather than resolved, and the filters applied upstream still readable in the inherited provenance slots.
#'
#' @slot results \code{DataFrame} with one row per gene. Read it with \code{\link{resultsTable}}.
#' @slot contrast String with the contrast the regions were tested on.
#' @slot method String with the aggregation used, \code{"simes"} or \code{"potential"}.
#' @slot links \code{DataFrame} with the region to gene links the aggregation was built on, one row per link.
#' @slot weighting String saying what weighted the regions of a gene, one among \code{"link.weight"}, \code{"distance"} and \code{"none"}.
#' @slot engine String with the engine the regions were fitted with.
#' @slot thresholds List with the cut-offs the directions were read at.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setClass representation prototype
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # Links of your own, four regions per gene
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fitCounts(fit)))[, 1:3]
#' geneNumber <- ceiling(nrow(regionTable) / 4)
#' geneIds <- rep(paste0("gene", seq_len(geneNumber)), each = 4)[seq_len(nrow(regionTable))]
#'
#' annotatedFit <- annotateRegions(fit, linkTable = data.frame(regionTable, gene.id = geneIds), verbose = FALSE)
#' results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' geneResults <- aggregateToGenes(results, method = "simes", verbose = FALSE)
#' geneResults
#'
#' @export
setClass(
    Class = "RegionSetDE.geneResults",
    contains = "RegionSetDE.provenance",
    representation = representation(
        results = "DataFrame",
        contrast = "character",
        method = "character",
        links = "DataFrame",
        weighting = "character",
        engine = "character",
        thresholds = "list"
    ),
    prototype = prototype(
        results = S4Vectors::DataFrame(),
        contrast = NA_character_,
        method = "simes",
        links = S4Vectors::DataFrame(),
        weighting = NA_character_,
        engine = NA_character_,
        thresholds = list()
    )
)


#' @importFrom methods setValidity
setValidity(
    Class = "RegionSetDE.geneResults",
    method = function(object) {
        issues <- character(0)

        if (length(object@method) == 0 | !all(object@method %in% c("simes", "potential"))) {
            issues <- c(issues, "The 'method' slot must hold only simes and potential.")
        }

        if (!(object@weighting %in% c("link.weight", "distance", "none"))) {
            issues <- c(issues, "The 'weighting' slot must be one among link.weight, distance and none.")
        }

        issues <- c(issues, .checkContrastLists(object))

        # A gene appearing twice would be counted twice by anything reading the table
        if (nrow(object@results) > 0 & anyDuplicated(object@results$gene.id) > 0) {
            issues <- c(issues, "The 'results' slot holds the same gene more than once.")
        }

        if (length(issues) == 0) {
            return(TRUE)
        }
        return(issues)
    }
)


#' @rdname RegionSetDE.geneResults-class
#'
#' @param object \code{RegionSetDE.geneResults} object.
#'
#' @return \code{show()} prints a summary of the object to the console.
#'
#' @importFrom methods setMethod
#' @importFrom utils head
#'
#' @export
setMethod(
    f = "show",
    signature = "RegionSetDE.geneResults",
    definition = function(object) {
        cat("An object of class 'RegionSetDE.geneResults'\n")
        cat("  contrast        :", object@contrast, "\n")
        cat("  engine          :", object@engine, "\n")
        cat("  aggregation     :", paste(object@method, collapse = ", "), "\n")
        cat("  weighting       :", object@weighting, "\n")
        cat("  genes           :", nrow(object@results), "\n")
        cat("  links           :", nrow(object@links), "\n")

        if (nrow(object@results) > 0 & "direction" %in% colnames(object@results)) {
            directionCounts <- table(factor(object@results$direction, levels = c("up", "down", "mixed", "unchanged")))
            cat("  direction       :", paste(names(directionCounts), directionCounts, sep = " ", collapse = ", "), "\n")
        }

        if (nrow(object@results) > 0) {
            shownColumns <- intersect(
                c(
                    "gene.id", "gene.symbol", "n.regions", "n.up", "n.down",
                    "direction", "score", "p.value", "FDR"
                ),
                colnames(object@results)
            )
            printTable <- utils::head(as.data.frame(object@results[, shownColumns, drop = FALSE], optional = TRUE), 6)
            numericColumns <- vapply(printTable, is.numeric, logical(1))
            printTable[numericColumns] <- lapply(printTable[numericColumns], signif, digits = 3)

            cat("\n")
            print(printTable, row.names = FALSE)
        }
        invisible(NULL)
    }
)


#' @rdname updateObject-RegionSetDE
#'
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "updateObject",
    signature = "RegionSetDE.geneResults",
    definition = function(object, ..., verbose = FALSE) {
        return(.updateProvenance(object = object, verbose = verbose))
    }
)


##########################################################################################
###    HELPERS OF THE CLASSES
##########################################################################################


#' @title .checkContrastLists
#'
#' @description Checks the list slots shared by the result classes: the cut-offs in \code{thresholds} and the two groups of the contrast in \code{contrast.groups}. Used by the validity methods.
#'
#' @param object \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or \code{RegionSetDE.geneResults} object.
#'
#' @return A character vector with the problems found, empty when there is none.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods .hasSlot
#'
#' @noRd

.checkContrastLists <-
    function(object) {
        issues <- character(0)

        # The labels of the results were set at these cut-offs, and they are read back by name
        if (length(object@thresholds) > 0 && (is.null(names(object@thresholds)) || any(names(object@thresholds) == ""))) {
            issues <- c(issues, "All the elements of the 'thresholds' slot must be named.")
        }

        if (nrow(object@results) > 0 && length(object@thresholds) > 0 && !("FDR" %in% names(object@thresholds))) {
            issues <- c(issues, "The 'thresholds' slot must hold the 'FDR' cut-off the results were labelled at.")
        }

        # A contrast between two groups names its column and both groups, or nothing at all
        if (methods::.hasSlot(object, "contrast.groups") && length(object@contrast.groups) > 0) {
            contrastGroups <- object@contrast.groups
            if (!all(c("column", "groups") %in% names(contrastGroups)) || length(contrastGroups$groups) != 2) {
                issues <- c(issues, "The 'contrast.groups' slot must be empty, or name the 'column' and the two 'groups' of the contrast.")
            }
        }

        return(issues)
    } # END function


#' @title .asSlotValue
#'
#' @description Brings a value to the class of the slot it goes in. A data.frame meant for a \code{DataFrame} slot is converted, so that the functions of the package can keep building their tables with dplyr.
#'
#' @param className String with the class of the object.
#' @param slotName String with the name of the slot.
#' @param value The value to store.
#'
#' @return The value, converted when the slot needs it.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods getSlots as
#'
#' @noRd

.asSlotValue <-
    function(className,
             slotName,
             value) {
        slotClass <- methods::getSlots(className)[slotName]

        if (!is.na(slotClass) && slotClass == "DataFrame" && is.data.frame(value)) {
            value <- methods::as(value, "DataFrame")
        }

        return(value)
    } # END function


#' @title .newObject
#'
#' @description Creates an object of one of the classes of the package, as \code{methods::new} does, converting the tables given as data.frames for the slots holding a \code{DataFrame}.
#'
#' @param Class String with the class of the object.
#' @param ... Values of the slots, named after them. An unnamed value, such as the \code{RangedSummarizedExperiment} a counts object extends, is passed on as it is.
#'
#' @return The new object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods new
#'
#' @noRd

.newObject <-
    function(Class,
             ...) {
        slotValues <- list(...)
        valueNames <- names(slotValues)

        for (i in seq_along(slotValues)) {
            if (!is.null(valueNames) && valueNames[i] != "") {
                slotValues[i] <- list(.asSlotValue(className = Class, slotName = valueNames[i], value = slotValues[[i]]))
            }
        }

        return(do.call(what = methods::new, args = c(list(Class = Class), slotValues)))
    } # END function


#' @title .replaceSlots
#'
#' @description Writes new values into the slots of an object, the way \code{object@slot <- value} would, with the class of every value checked against its slot. The tables given as data.frames are converted for the slots holding a \code{DataFrame}. Validity is left to the caller, since several slots may have to change before the object is consistent again.
#'
#' @param object An object of one of the classes of the package.
#' @param ... New values, named after their slots.
#'
#' @return The object with the new values.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods slot<-
#'
#' @noRd

.replaceSlots <-
    function(object,
             ...) {
        slotValues <- list(...)

        for (slotName in names(slotValues)) {
            methods::slot(object, slotName) <- .asSlotValue(className = class(object), slotName = slotName, value = slotValues[[slotName]])
        }

        return(object)
    } # END function
