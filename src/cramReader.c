/* Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio. */

/*
 * Reading of BAM and CRAM files, CRAM without the reference genome.
 *
 * Both formats go through the same code, so that they give the same reads and
 * the same fragments. CRAM stores the bases of a read as their differences from
 * the reference, so
 * htslib needs the reference to rebuild the sequence, the MD and NM tags and,
 * within a container holding reads of several chromosomes, the template length
 * (TLEN). Read name, flag, chromosome, position, mapping quality and CIGAR are
 * stored on their own. Only these are decoded here, so no reference is ever
 * looked for, neither on disk nor on the EBI server. Paired-end fragments are
 * rebuilt from the two mates, matched by read name: from the 5' end of the
 * forward mate to the 5' end of the reverse mate, as csaw defines them, so that
 * mates running past each other give the fragment and not the span of the reads.
 *
 * Every call runs inside R_UnwindProtect, so that the file, the index and the
 * buffers are released when an error or an interrupt leaves the call early.
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <R.h>
#include <Rinternals.h>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/khash.h"

/* Records read between two checks for an interrupt */
#define INTERRUPT_CHECK 1048576

/* Vectors filled while a region is read */
#define BUFFER_NUMBER 5

/* A mate waiting for the other one, under its read name */
typedef struct {
    int position;
    int end;
    int mapq;
    int flag;
} pendingMate;

KHASH_MAP_INIT_STR(mateTable, pendingMate)

/* Integer vector growing while the records are read */
typedef struct {
    int *values;
    R_xlen_t length;
    R_xlen_t capacity;
} intBuffer;

/* Everything a call opens, and its arguments */
typedef struct {
    const char *path;
    const char *indexPath;
    htsFile *file;
    sam_hdr_t *header;
    hts_idx_t *index;
    hts_itr_t *iterator;
    bam1_t *record;
    khash_t(mateTable) *mates;
    intBuffer buffers[BUFFER_NUMBER];
    SEXP chromosomes;
    SEXP starts;
    SEXP ends;
    int pairedEnd;
    int minMapq;
    int removeDuplicates;
    int padding;
    int nRecords;
} cramReader;


/*---------------------------*/
/* Buffers                   */
/*---------------------------*/

static void bufferPush(intBuffer *buffer, int value)
{
    if (buffer->length == buffer->capacity) {
        R_xlen_t capacity = buffer->capacity == 0 ? 4096 : 2 * buffer->capacity;
        int *values = (int *) realloc(buffer->values, (size_t) capacity * sizeof(int));

        /* The old block is still owned by the buffer and freed by the cleanup */
        if (values == NULL) {
            Rf_error("Not enough memory to read the CRAM file.");
        }

        buffer->values = values;
        buffer->capacity = capacity;
    }

    buffer->values[buffer->length++] = value;
}


static SEXP bufferToVector(const intBuffer *buffer)
{
    SEXP vector = PROTECT(Rf_allocVector(INTSXP, buffer->length));

    if (buffer->length > 0) {
        memcpy(INTEGER(vector), buffer->values, (size_t) buffer->length * sizeof(int));
    }

    UNPROTECT(1);
    return vector;
}


/* The strand as scanBam gives it, a factor with the levels +, - and * */
static SEXP bufferToStrand(const intBuffer *buffer)
{
    SEXP strand = PROTECT(bufferToVector(buffer));
    SEXP levels = PROTECT(Rf_allocVector(STRSXP, 3));

    SET_STRING_ELT(levels, 0, Rf_mkChar("+"));
    SET_STRING_ELT(levels, 1, Rf_mkChar("-"));
    SET_STRING_ELT(levels, 2, Rf_mkChar("*"));
    Rf_setAttrib(strand, R_LevelsSymbol, levels);
    Rf_setAttrib(strand, R_ClassSymbol, Rf_mkString("factor"));

    UNPROTECT(2);
    return strand;
}


/* Bounds of the fragment of two mates: from the 5' end of the forward mate to the
   5' end of the reverse one. Mates on the same strand, or with the forward mate
   starting past the end of the reverse one, make no fragment and give 0 */
static int fragmentBounds(const pendingMate *mate1, const pendingMate *mate2, int *start, int *end)
{
    const pendingMate *forwardMate, *reverseMate;
    int reverse1 = (mate1->flag & BAM_FREVERSE) != 0;
    int reverse2 = (mate2->flag & BAM_FREVERSE) != 0;

    if (reverse1 == reverse2) {
        return 0;
    }

    forwardMate = reverse1 ? mate2 : mate1;
    reverseMate = reverse1 ? mate1 : mate2;

    if (forwardMate->position > reverseMate->end) {
        return 0;
    }

    *start = forwardMate->position;
    *end = reverseMate->end;
    return 1;
}


/* Copy of a read name kept as key of the mate table (strdup is not standard C before C23) */
static char *copyName(const char *name)
{
    size_t length = strlen(name) + 1;
    char *copy = (char *) malloc(length);

    if (copy != NULL) {
        memcpy(copy, name, length);
    }

    return copy;
}


static void resetBuffers(cramReader *reader)
{
    int i;

    for (i = 0; i < BUFFER_NUMBER; i++) {
        reader->buffers[i].length = 0;
    }
}


/*---------------------------*/
/* Opening and releasing     */
/*---------------------------*/

static void clearMates(khash_t(mateTable) *mates)
{
    khint_t k;

    if (mates == NULL) {
        return;
    }

    for (k = kh_begin(mates); k != kh_end(mates); ++k) {
        if (kh_exist(mates, k)) {
            free((char *) kh_key(mates, k));
        }
    }

    kh_clear(mateTable, mates);
}


static void releaseReader(void *data, Rboolean jump)
{
    cramReader *reader = (cramReader *) data;
    int i;

    (void) jump;

    if (reader->mates != NULL) {
        clearMates(reader->mates);
        kh_destroy(mateTable, reader->mates);
        reader->mates = NULL;
    }
    if (reader->iterator != NULL) {
        hts_itr_destroy(reader->iterator);
        reader->iterator = NULL;
    }
    if (reader->index != NULL) {
        hts_idx_destroy(reader->index);
        reader->index = NULL;
    }
    if (reader->record != NULL) {
        bam_destroy1(reader->record);
        reader->record = NULL;
    }
    if (reader->header != NULL) {
        sam_hdr_destroy(reader->header);
        reader->header = NULL;
    }
    if (reader->file != NULL) {
        hts_close(reader->file);
        reader->file = NULL;
    }
    for (i = 0; i < BUFFER_NUMBER; i++) {
        free(reader->buffers[i].values);
        reader->buffers[i].values = NULL;
        reader->buffers[i].length = 0;
        reader->buffers[i].capacity = 0;
    }
}


static void openAlignments(cramReader *reader, int requiredFields, int withIndex, int onlyCram)
{
    enum htsExactFormat fileFormat;

    reader->file = hts_open(reader->path, "r");
    if (reader->file == NULL) {
        Rf_error("Cannot open '%s'.", reader->path);
    }

    fileFormat = hts_get_format(reader->file)->format;
    if (fileFormat != bam && fileFormat != cram) {
        Rf_error("'%s' is not a BAM or CRAM file.", reader->path);
    }
    if (onlyCram && fileFormat != cram) {
        Rf_error("'%s' is not a CRAM file.", reader->path);
    }

    /* Of a CRAM file only the fields stored apart from the reference are decoded */
    if (fileFormat == cram &&
        (hts_set_opt(reader->file, CRAM_OPT_REQUIRED_FIELDS, requiredFields) != 0 ||
         hts_set_opt(reader->file, CRAM_OPT_DECODE_MD, 0) != 0)) {
        Rf_error("Cannot set the decoding options of '%s'.", reader->path);
    }

    reader->header = sam_hdr_read(reader->file);
    if (reader->header == NULL) {
        Rf_error("Cannot read the header of '%s'.", reader->path);
    }

    reader->record = bam_init1();
    if (reader->record == NULL) {
        Rf_error("Not enough memory to read '%s'.", reader->path);
    }

    if (withIndex) {
        reader->index = sam_index_load2(reader->file, reader->path, reader->indexPath);
        if (reader->index == NULL) {
            Rf_error("Cannot load the index '%s' of '%s'.", reader->indexPath, reader->path);
        }
    }
}


static const char *stringArgument(SEXP value, const char *name)
{
    if (!Rf_isString(value) || XLENGTH(value) != 1 || STRING_ELT(value, 0) == NA_STRING) {
        Rf_error("'%s' must be a single string.", name);
    }

    return Rf_translateChar(STRING_ELT(value, 0));
}


static int integerArgument(SEXP value, const char *name)
{
    int number = Rf_asInteger(value);

    if (number == NA_INTEGER) {
        Rf_error("'%s' must be a single integer.", name);
    }

    return number;
}


/*---------------------------*/
/* Reads of one region       */
/*---------------------------*/

/* Single-end: every primary mapped read overlapping the region, with its strand
   and the width it covers on the genome */
static SEXP readSingleRegion(cramReader *reader, int tid, int start, int end)
{
    intBuffer *positions = &reader->buffers[0];
    intBuffer *strands = &reader->buffers[1];
    intBuffer *widths = &reader->buffers[2];
    long long recordNumber = 0;
    int status;

    resetBuffers(reader);

    reader->iterator = sam_itr_queryi(reader->index, tid, (hts_pos_t) start - 1, (hts_pos_t) end);
    if (reader->iterator == NULL) {
        Rf_error("Cannot query '%s' at %s:%d-%d.", reader->path, sam_hdr_tid2name(reader->header, tid), start, end);
    }

    while ((status = sam_itr_next(reader->file, reader->iterator, reader->record)) >= 0) {
        const bam1_core_t *core = &reader->record->core;

        if (++recordNumber % INTERRUPT_CHECK == 0) {
            R_CheckUserInterrupt();
        }

        if (core->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            continue;
        }
        if (reader->removeDuplicates && (core->flag & BAM_FDUP)) {
            continue;
        }
        if (core->qual < reader->minMapq) {
            continue;
        }

        bufferPush(positions, (int) core->pos + 1);
        bufferPush(strands, (core->flag & BAM_FREVERSE) ? 2 : 1);
        bufferPush(widths, (int) bam_cigar2rlen(core->n_cigar, bam_get_cigar(reader->record)));
    }

    if (status < -1) {
        Rf_error("Cannot decode the records of '%s' at %s:%d-%d.", reader->path, sam_hdr_tid2name(reader->header, tid), start, end);
    }

    hts_itr_destroy(reader->iterator);
    reader->iterator = NULL;

    SEXP region = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP fieldNames = PROTECT(Rf_allocVector(STRSXP, 3));

    SET_VECTOR_ELT(region, 0, bufferToVector(positions));
    SET_VECTOR_ELT(region, 1, bufferToStrand(strands));
    SET_VECTOR_ELT(region, 2, bufferToVector(widths));
    SET_STRING_ELT(fieldNames, 0, Rf_mkChar("pos"));
    SET_STRING_ELT(fieldNames, 1, Rf_mkChar("strand"));
    SET_STRING_ELT(fieldNames, 2, Rf_mkChar("width"));
    Rf_setAttrib(region, R_NamesSymbol, fieldNames);

    UNPROTECT(2);
    return region;
}


/* Paired-end: the first mate of every proper pair overlapping the region, with
   the position of its mate, the start and the length of the fragment and the
   mapping quality of the mate. The region is widened by the padding so that the
   mate of a read near the border is read as well. */
static SEXP readPairedRegion(cramReader *reader, int tid, int start, int end)
{
    intBuffer *positions = &reader->buffers[0];
    intBuffer *matePositions = &reader->buffers[1];
    intBuffer *fragmentSizes = &reader->buffers[2];
    intBuffer *fragmentStarts = &reader->buffers[3];
    intBuffer *mateMapqs = &reader->buffers[4];
    hts_pos_t windowStart = (hts_pos_t) start - 1 - reader->padding;
    hts_pos_t windowEnd = (hts_pos_t) end + reader->padding;
    long long recordNumber = 0;
    int status;

    resetBuffers(reader);

    if (windowStart < 0) {
        windowStart = 0;
    }

    reader->iterator = sam_itr_queryi(reader->index, tid, windowStart, windowEnd);
    if (reader->iterator == NULL) {
        Rf_error("Cannot query '%s' at %s:%d-%d.", reader->path, sam_hdr_tid2name(reader->header, tid), start, end);
    }

    while ((status = sam_itr_next(reader->file, reader->iterator, reader->record)) >= 0) {
        const bam1_core_t *core = &reader->record->core;
        pendingMate current, firstMate, secondMate;
        khint_t k;
        int isFirst, absent, fragmentStart, fragmentEnd;

        if (++recordNumber % INTERRUPT_CHECK == 0) {
            R_CheckUserInterrupt();
        }

        if (!(core->flag & BAM_FPAIRED) || !(core->flag & BAM_FPROPER_PAIR)) {
            continue;
        }
        if (core->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            continue;
        }

        isFirst = (core->flag & BAM_FREAD1) != 0;
        if (isFirst == ((core->flag & BAM_FREAD2) != 0)) {
            continue;
        }

        current.position = (int) core->pos + 1;
        current.end = (int) (core->pos + bam_cigar2rlen(core->n_cigar, bam_get_cigar(reader->record)));
        current.mapq = core->qual;
        current.flag = core->flag;

        /* The first mate of the pair to arrive waits for the other one */
        k = kh_get(mateTable, reader->mates, bam_get_qname(reader->record));

        if (k == kh_end(reader->mates)) {
            char *key = copyName(bam_get_qname(reader->record));

            if (key == NULL) {
                Rf_error("Not enough memory to read '%s'.", reader->path);
            }

            k = kh_put(mateTable, reader->mates, key, &absent);
            if (absent < 0) {
                free(key);
                Rf_error("Not enough memory to read '%s'.", reader->path);
            }

            kh_val(reader->mates, k) = current;
            continue;
        }

        /* Two records of the same mate under one name: the later one waits instead */
        if (((kh_val(reader->mates, k).flag & BAM_FREAD1) != 0) == isFirst) {
            kh_val(reader->mates, k) = current;
            continue;
        }

        firstMate = isFirst ? current : kh_val(reader->mates, k);
        secondMate = isFirst ? kh_val(reader->mates, k) : current;

        free((char *) kh_key(reader->mates, k));
        kh_del(mateTable, reader->mates, k);

        /* The filters of the counting apply to the first mate, the mapping quality of the mate is returned */
        if (reader->removeDuplicates && (firstMate.flag & BAM_FDUP)) {
            continue;
        }
        if (firstMate.mapq < reader->minMapq) {
            continue;
        }

        /* Only the first mates overlapping the region itself, as scanBam would return them */
        if (firstMate.position > end || (firstMate.end < firstMate.position ? firstMate.position : firstMate.end) < start) {
            continue;
        }

        if (!fragmentBounds(&firstMate, &secondMate, &fragmentStart, &fragmentEnd)) {
            continue;
        }

        bufferPush(positions, firstMate.position);
        bufferPush(matePositions, secondMate.position);
        bufferPush(fragmentSizes, fragmentEnd - fragmentStart + 1);
        bufferPush(fragmentStarts, fragmentStart);
        bufferPush(mateMapqs, secondMate.mapq);
    }

    if (status < -1) {
        Rf_error("Cannot decode the records of '%s' at %s:%d-%d.", reader->path, sam_hdr_tid2name(reader->header, tid), start, end);
    }

    hts_itr_destroy(reader->iterator);
    reader->iterator = NULL;

    /* Mates whose partner lies beyond the widened region belong to fragments longer than the padding */
    clearMates(reader->mates);

    SEXP region = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP fieldNames = PROTECT(Rf_allocVector(STRSXP, 5));
    SEXP tagList = PROTECT(Rf_allocVector(VECSXP, 1));

    SET_VECTOR_ELT(tagList, 0, bufferToVector(mateMapqs));
    Rf_setAttrib(tagList, R_NamesSymbol, Rf_mkString("MQ"));

    SET_VECTOR_ELT(region, 0, bufferToVector(positions));
    SET_VECTOR_ELT(region, 1, bufferToVector(matePositions));
    SET_VECTOR_ELT(region, 2, bufferToVector(fragmentSizes));
    SET_VECTOR_ELT(region, 3, bufferToVector(fragmentStarts));
    SET_VECTOR_ELT(region, 4, tagList);
    SET_STRING_ELT(fieldNames, 0, Rf_mkChar("pos"));
    SET_STRING_ELT(fieldNames, 1, Rf_mkChar("mpos"));
    SET_STRING_ELT(fieldNames, 2, Rf_mkChar("isize"));
    SET_STRING_ELT(fieldNames, 3, Rf_mkChar("fragment.start"));
    SET_STRING_ELT(fieldNames, 4, Rf_mkChar("tag"));
    Rf_setAttrib(region, R_NamesSymbol, fieldNames);

    UNPROTECT(3);
    return region;
}


/*---------------------------*/
/* Bodies of the calls       */
/*---------------------------*/

static SEXP readRegionsBody(void *data)
{
    cramReader *reader = (cramReader *) data;
    int requiredFields = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR;
    R_xlen_t i, regionNumber = XLENGTH(reader->chromosomes);

    /* The read name is what pairs the mates */
    if (reader->pairedEnd) {
        requiredFields |= SAM_QNAME;
    }

    openAlignments(reader, requiredFields, 1, 0);

    reader->mates = kh_init(mateTable);
    if (reader->mates == NULL) {
        Rf_error("Not enough memory to read '%s'.", reader->path);
    }

    SEXP regionList = PROTECT(Rf_allocVector(VECSXP, regionNumber));

    for (i = 0; i < regionNumber; i++) {
        const char *chromosome = Rf_translateChar(STRING_ELT(reader->chromosomes, i));
        int tid = sam_hdr_name2tid(reader->header, chromosome);
        int start = INTEGER(reader->starts)[i];
        int end = INTEGER(reader->ends)[i];

        if (tid < 0) {
            Rf_error("The chromosome '%s' is not in the header of '%s'.", chromosome, reader->path);
        }
        if (start == NA_INTEGER || end == NA_INTEGER || start < 1 || end < start) {
            Rf_error("Invalid region %s:%d-%d.", chromosome, start, end);
        }

        SET_VECTOR_ELT(regionList, i, reader->pairedEnd ? readPairedRegion(reader, tid, start, end) : readSingleRegion(reader, tid, start, end));
    }

    UNPROTECT(1);
    return regionList;
}


static SEXP readTargetsBody(void *data)
{
    cramReader *reader = (cramReader *) data;
    int i, targetNumber;

    openAlignments(reader, SAM_FLAG, 0, 1);
    targetNumber = sam_hdr_nref(reader->header);

    SEXP targetLengths = PROTECT(Rf_allocVector(INTSXP, targetNumber));
    SEXP targetNames = PROTECT(Rf_allocVector(STRSXP, targetNumber));

    for (i = 0; i < targetNumber; i++) {
        hts_pos_t targetLength = sam_hdr_tid2len(reader->header, i);

        INTEGER(targetLengths)[i] = targetLength > INT_MAX ? NA_INTEGER : (int) targetLength;
        SET_STRING_ELT(targetNames, i, Rf_mkChar(sam_hdr_tid2name(reader->header, i)));
    }

    Rf_setAttrib(targetLengths, R_NamesSymbol, targetNames);

    UNPROTECT(2);
    return targetLengths;
}


static SEXP readHeadFlagsBody(void *data)
{
    cramReader *reader = (cramReader *) data;
    int status = 0, recordNumber = 0;

    openAlignments(reader, SAM_FLAG, 0, 1);

    while (recordNumber < reader->nRecords && (status = sam_read1(reader->file, reader->header, reader->record)) >= 0) {
        bufferPush(&reader->buffers[0], reader->record->core.flag);
        recordNumber++;
    }

    if (status < -1) {
        Rf_error("Cannot decode the records of '%s'.", reader->path);
    }

    return bufferToVector(&reader->buffers[0]);
}


/* The counts samtools idxstats gives, from a pass over the whole file, since a CRAM index holds none */
static SEXP readIndexStatsBody(void *data)
{
    cramReader *reader = (cramReader *) data;
    long long recordNumber = 0;
    int i, status, targetNumber;
    double *mappedReads, *unmappedReads;

    openAlignments(reader, SAM_FLAG | SAM_RNAME, 0, 1);
    targetNumber = sam_hdr_nref(reader->header);

    /* The last slot holds the reads placed on no chromosome */
    mappedReads = (double *) R_alloc((size_t) targetNumber + 1, sizeof(double));
    unmappedReads = (double *) R_alloc((size_t) targetNumber + 1, sizeof(double));
    memset(mappedReads, 0, ((size_t) targetNumber + 1) * sizeof(double));
    memset(unmappedReads, 0, ((size_t) targetNumber + 1) * sizeof(double));

    while ((status = sam_read1(reader->file, reader->header, reader->record)) >= 0) {
        int tid = reader->record->core.tid;
        int slot = (tid >= 0 && tid < targetNumber) ? tid : targetNumber;

        if (++recordNumber % INTERRUPT_CHECK == 0) {
            R_CheckUserInterrupt();
        }

        if (reader->record->core.flag & BAM_FUNMAP) {
            unmappedReads[slot] += 1;
        } else {
            mappedReads[slot] += 1;
        }
    }

    if (status < -1) {
        Rf_error("Cannot decode the records of '%s'.", reader->path);
    }

    SEXP statList = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP fieldNames = PROTECT(Rf_allocVector(STRSXP, 4));
    SEXP targetNames = PROTECT(Rf_allocVector(STRSXP, targetNumber + 1));
    SEXP targetLengths = PROTECT(Rf_allocVector(INTSXP, targetNumber + 1));
    SEXP mapped = PROTECT(Rf_allocVector(REALSXP, targetNumber + 1));
    SEXP unmapped = PROTECT(Rf_allocVector(REALSXP, targetNumber + 1));

    for (i = 0; i < targetNumber; i++) {
        hts_pos_t targetLength = sam_hdr_tid2len(reader->header, i);

        SET_STRING_ELT(targetNames, i, Rf_mkChar(sam_hdr_tid2name(reader->header, i)));
        INTEGER(targetLengths)[i] = targetLength > INT_MAX ? NA_INTEGER : (int) targetLength;
    }
    SET_STRING_ELT(targetNames, targetNumber, Rf_mkChar("*"));
    INTEGER(targetLengths)[targetNumber] = 0;
    memcpy(REAL(mapped), mappedReads, ((size_t) targetNumber + 1) * sizeof(double));
    memcpy(REAL(unmapped), unmappedReads, ((size_t) targetNumber + 1) * sizeof(double));

    SET_VECTOR_ELT(statList, 0, targetNames);
    SET_VECTOR_ELT(statList, 1, targetLengths);
    SET_VECTOR_ELT(statList, 2, mapped);
    SET_VECTOR_ELT(statList, 3, unmapped);
    SET_STRING_ELT(fieldNames, 0, Rf_mkChar("seqnames"));
    SET_STRING_ELT(fieldNames, 1, Rf_mkChar("seqlength"));
    SET_STRING_ELT(fieldNames, 2, Rf_mkChar("mapped"));
    SET_STRING_ELT(fieldNames, 3, Rf_mkChar("unmapped"));
    Rf_setAttrib(statList, R_NamesSymbol, fieldNames);

    UNPROTECT(6);
    return statList;
}


/*---------------------------*/
/* Entry points              */
/*---------------------------*/

static SEXP runReader(SEXP (*body)(void *), cramReader *reader)
{
    SEXP continuation = PROTECT(R_MakeUnwindCont());
    SEXP result = R_UnwindProtect(body, reader, releaseReader, reader, continuation);

    UNPROTECT(1);
    return result;
}


SEXP C_readAlignmentRegions(SEXP path, SEXP indexPath, SEXP chromosomes, SEXP starts, SEXP ends,
                       SEXP pairedEnd, SEXP minMapq, SEXP removeDuplicates, SEXP padding)
{
    cramReader reader;

    memset(&reader, 0, sizeof(reader));
    reader.path = stringArgument(path, "path");
    reader.indexPath = stringArgument(indexPath, "indexPath");

    if (!Rf_isString(chromosomes) || !Rf_isInteger(starts) || !Rf_isInteger(ends) ||
        XLENGTH(starts) != XLENGTH(chromosomes) || XLENGTH(ends) != XLENGTH(chromosomes)) {
        Rf_error("The regions must be given as a character vector of chromosomes and two integer vectors of the same length.");
    }

    reader.chromosomes = chromosomes;
    reader.starts = starts;
    reader.ends = ends;
    reader.pairedEnd = Rf_asLogical(pairedEnd) == TRUE;
    reader.minMapq = integerArgument(minMapq, "minMapq");
    reader.removeDuplicates = Rf_asLogical(removeDuplicates) == TRUE;
    reader.padding = integerArgument(padding, "padding");

    if (reader.padding < 0) {
        reader.padding = 0;
    }

    return runReader(readRegionsBody, &reader);
}


SEXP C_cramTargets(SEXP path)
{
    cramReader reader;

    memset(&reader, 0, sizeof(reader));
    reader.path = stringArgument(path, "path");

    return runReader(readTargetsBody, &reader);
}


SEXP C_cramHeadFlags(SEXP path, SEXP nRecords)
{
    cramReader reader;

    memset(&reader, 0, sizeof(reader));
    reader.path = stringArgument(path, "path");
    reader.nRecords = integerArgument(nRecords, "nRecords");

    return runReader(readHeadFlagsBody, &reader);
}


SEXP C_cramIndexStats(SEXP path)
{
    cramReader reader;

    memset(&reader, 0, sizeof(reader));
    reader.path = stringArgument(path, "path");

    return runReader(readIndexStatsBody, &reader);
}
