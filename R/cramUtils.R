# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .isCramFile
#'
#' @description Tells which alignment files are CRAM files, from their extension. Every other file is read as a BAM file.
#'
#' @param files Character vector with the paths of the alignment files.
#'
#' @return A logical vector with one value per file.
#'
#' @author Sebastian Gregoricchio
#'
#' @noRd

.isCramFile <-
    function(files) {
        return(grepl("\\.cram$", files, ignore.case = TRUE))
    }


#' @title .alignmentTargets
#'
#' @description Returns the chromosomes of an alignment file with their lengths, as \code{Rsamtools::scanBamHeader} gives them for a BAM file. Rsamtools cannot open CRAM files, whose header is read by the C code of the package.
#'
#' @param bamFile String with the path of the BAM or CRAM file.
#'
#' @return A named integer vector with the length of every chromosome.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools scanBamHeader
#'
#' @noRd

.alignmentTargets <-
    function(bamFile) {
        if (.isCramFile(bamFile)) {
            return(.Call(C_cramTargets, path.expand(bamFile)))
        }

        return(Rsamtools::scanBamHeader(bamFile)[[1]]$targets)
    }


#' @title .alignmentHeadFlags
#'
#' @description Returns the flags of the first records of an alignment file, read from its head without the index.
#'
#' @param bamFile String with the path of the BAM or CRAM file.
#' @param nRecords Numeric value with the number of records to read.
#'
#' @return An integer vector with the flag of every record read, empty for a file without records.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools BamFile scanBam ScanBamParam
#'
#' @noRd

.alignmentHeadFlags <-
    function(bamFile,
             nRecords) {
        if (.isCramFile(bamFile)) {
            return(.Call(C_cramHeadFlags, path.expand(bamFile), as.integer(nRecords)))
        }

        # No index is needed to read a file from its first record
        return(Rsamtools::scanBam(Rsamtools::BamFile(file = bamFile, yieldSize = as.integer(nRecords)),
            param = Rsamtools::ScanBamParam(what = "flag")
        )[[1]]$flag)
    }


#' @title .alignmentIndexStats
#'
#' @description Counts the mapped and unmapped records of every chromosome of an alignment file, as \code{samtools idxstats} does. The counts of a BAM file come from its index. A CRAM index holds no counts, so the whole file is read once, decoding only the flag and the chromosome of every record: a matter of minutes for a large file.
#'
#' @param bamFile String with the path of the BAM or CRAM file.
#'
#' @return A data.frame with the columns \code{seqnames}, \code{seqlength}, \code{mapped} and \code{unmapped}, with one row per chromosome and one, named \code{*}, for the reads placed on none.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools idxstatsBam
#'
#' @noRd

.alignmentIndexStats <-
    function(bamFile) {
        if (.isCramFile(bamFile)) {
            return(as.data.frame(.Call(C_cramIndexStats, path.expand(bamFile)), stringsAsFactors = FALSE))
        }

        return(Rsamtools::idxstatsBam(.bamWithIndex(bamFile)))
    }


#' @title .scanAlignments
#'
#' @description Reads the records of a BAM or CRAM file over a set of stretches of genome and returns them as \code{Rsamtools::scanBam} does, one list per stretch. Both formats are read by the C code of the package on htslib, so that they give the same reads and the same fragments. The filters are those of the counting: primary mapped reads only, no duplicates when they are removed, a minimum mapping quality and, for paired-end data, proper pairs only.
#'
#' CRAM files are read without the reference genome. CRAM stores the bases of a read as their differences from the reference, while the read name, the flag, the position, the mapping quality and the CIGAR are stored on their own, and only these are decoded. The template length (TLEN) is not used, in either format: htslib rebuilds it from the reference as soon as a CRAM container holds reads of several chromosomes, and aligners do not all compute it alike. The two mates of a pair are matched by read name instead, and the fragment runs from the 5' end of the forward mate to the 5' end of the reverse mate, as in csaw. Mates on the same strand, or with the forward mate starting past the end of the reverse one, make no fragment. Each stretch is read further on both sides by the longest fragment accepted, so that the mate of a read near the border is found, and the mapping quality of the mate comes from its own record.
#'
#' @param bamFile String with the path of the BAM or CRAM file.
#' @param which \code{GRanges} with the stretches of genome to read, named after the chromosomes of the file.
#' @param isPairedEnd Logical value, \code{TRUE} for a paired-end file.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment, which is also how far beyond the stretches the file is read to find the mates.
#'
#' @return A list with one element per stretch, named \code{chromosome:start-end}. For a paired-end file, the \code{pos} of the first mate of every proper pair making a fragment, the \code{mpos} of its mate, the \code{fragment.start}, the fragment length \code{isize} and the mapping quality of the mate in \code{tag$MQ}. For a single-end file, the \code{pos}, the \code{strand} and the \code{width} on the genome of every read.
#'
#' @details Only the fragments whose two mates are both found are returned: a fragment longer than \code{maxFragmentLength} may be missing, and is dropped by the counting anyway.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#'
#' @noRd

.scanAlignments <-
    function(bamFile,
             which,
             isPairedEnd,
             minMapq,
             removeDuplicates,
             maxFragmentLength) {
        indexPath <- .bamIndexPath(bamFile)
        if (is.na(indexPath)) {
            stop("The file '", basename(bamFile), "' is not indexed.", call. = FALSE)
        }

        checkMapq <- isTRUE(!is.na(minMapq[1]) & minMapq[1] > 0)

        # Beyond a gigabase the mates are looked for over the whole chromosome
        padding <- if (is.finite(maxFragmentLength[1])) {
            as.integer(min(ceiling(maxFragmentLength[1]), 1e9))
        } else {
            1e9L
        }

        whichChromosomes <- as.character(GenomeInfoDb::seqnames(which))
        whichStarts <- as.integer(BiocGenerics::start(which))
        whichEnds <- as.integer(BiocGenerics::end(which))

        readList <- .Call(
            C_readAlignmentRegions,
            path.expand(bamFile),
            path.expand(indexPath),
            whichChromosomes,
            whichStarts,
            whichEnds,
            isTRUE(isPairedEnd),
            if (checkMapq) {
                as.integer(minMapq[1])
            } else {
                0L
            },
            isTRUE(removeDuplicates),
            padding
        )

        # A reader compiled from older sources gives no fragment start, and every fragment would be lost without a word
        if (isTRUE(isPairedEnd) & length(readList) > 0) {
            if (is.null(readList[[1]]$fragment.start)) {
                stop("The compiled reader of RegionSetDE is older than its R code. Reinstall the package from a clean source folder, without the old src/*.o and src/*.so files.", call. = FALSE)
            }
        }

        # The names scanBam gives, which the callers use to find each stretch
        names(readList) <- paste0(whichChromosomes, ":", whichStarts, "-", whichEnds)

        return(readList)
    } # END function
