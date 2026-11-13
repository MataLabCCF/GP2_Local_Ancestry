/*
 * extract_tracts.c
 *
 * C rewrite of extract_tracts.py
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <zlib.h>

/* ── Constants ─────────────────────────────────────────────────────────── */

#define VERSION "1.2.0-c"
#define MAX_LINE       (64 * 1024 * 1024)  /* 64 MiB per VCF line – large cohorts need this */
#define MAX_ANCS       64                   /* max supported ancestries */
#define MAX_PATH_LEN   4096
#define MAX_INDIV      1000000              /* generous cap on individuals */

/* ── Debug helper macros ───────────────────────────────────────────────── */
/* Set DEBUG=1 at compile time (e.g. gcc -DDEBUG=1 ...) for verbose output */
#ifndef DEBUG
#define DEBUG 0
#endif

#define DBG(fmt, ...)                                               \
    do {                                                            \
        if (DEBUG)                                                  \
            fprintf(stderr, "[DEBUG %s:%d] " fmt "\n",              \
                    __FILE__, __LINE__, ##__VA_ARGS__);              \
    } while (0)

#define INFO(fmt, ...)                                              \
    fprintf(stderr, "[INFO] " fmt "\n", ##__VA_ARGS__)

#define WARN(fmt, ...)                                              \
    fprintf(stderr, "[WARN] " fmt "\n", ##__VA_ARGS__)

#define ERR(fmt, ...)                                               \
    fprintf(stderr, "[ERROR] " fmt "\n", ##__VA_ARGS__)

/* ── StrBuf: growable string buffer (avoids repeated strcat scanning) ── */

typedef struct {
    char  *buf;
    size_t len;
    size_t cap;
} StrBuf;

static void strbuf_init(StrBuf *sb, size_t cap) {
    sb->buf = (char *)malloc(cap);
    if (!sb->buf) { ERR("strbuf malloc failed"); exit(1); }
    sb->buf[0] = '\0';
    sb->len = 0;
    sb->cap = cap;
}

static void strbuf_clear(StrBuf *sb) {
    sb->len = 0;
    sb->buf[0] = '\0';
}

static void strbuf_append(StrBuf *sb, const char *s, size_t slen) {
    if (sb->len + slen + 1 > sb->cap) {
        sb->cap = (sb->len + slen + 1) * 2;
        sb->buf = (char *)realloc(sb->buf, sb->cap);
        if (!sb->buf) { ERR("strbuf realloc failed"); exit(1); }
    }
    memcpy(sb->buf + sb->len, s, slen);
    sb->len += slen;
    sb->buf[sb->len] = '\0';
}

static void strbuf_free(StrBuf *sb) {
    free(sb->buf);
    sb->buf = NULL;
    sb->len = sb->cap = 0;
}

/* ── Utility: open a file for writing, optionally gzip-compressed ──── */

typedef struct {
    FILE  *fp;
    gzFile gz;
    int    is_gz;
} OutFile;

static OutFile outfile_open(const char *path, int compressed) {
    OutFile o = {NULL, NULL, compressed};
    if (compressed) {
        /* Level 9, as Python's gzip.open in extract_tracts.py. zlib's default (6)
         * searches too few matches to reach the identical previous hapcount row,
         * so hapcount files come out ~17x larger. */
        o.gz = gzopen(path, "wb9");
        if (!o.gz) { ERR("Cannot open %s for writing (gz)", path); exit(1); }
    } else {
        o.fp = fopen(path, "w");
        if (!o.fp) { ERR("Cannot open %s for writing", path); exit(1); }
    }
    return o;
}

static void outfile_write(OutFile *o, const char *data, int len) {
    if (o->is_gz) gzwrite(o->gz, data, len);
    else          fwrite(data, 1, len, o->fp);
}

static void outfile_close(OutFile *o) {
    if (o->is_gz && o->gz) gzclose(o->gz);
    else if (o->fp)        fclose(o->fp);
    o->fp = NULL;
    o->gz = NULL;
}

/* ── Utility: read a line from plain or gzip file ─────────────────── */

typedef struct {
    FILE  *fp;
    gzFile gz;
    int    is_gz;
} InFile;

static InFile infile_open(const char *path, int compressed) {
    InFile f = {NULL, NULL, compressed};
    if (compressed) {
        f.gz = gzopen(path, "rb");
        if (!f.gz) { ERR("Cannot open %s for reading (gz)", path); exit(1); }
    } else {
        f.fp = fopen(path, "r");
        if (!f.fp) { ERR("Cannot open %s for reading", path); exit(1); }
    }
    return f;
}

static char *infile_gets(InFile *f, char *buf, int size) {
    if (f->is_gz) return gzgets(f->gz, buf, size);
    else          return fgets(buf, size, f->fp);
}

static void infile_close(InFile *f) {
    if (f->is_gz && f->gz) gzclose(f->gz);
    else if (f->fp)        fclose(f->fp);
}

/* ── Parse VCF filepath to extract path, prefix, and compression ───── */

typedef struct {
    char path[MAX_PATH_LEN];
    char prefix[MAX_PATH_LEN];
    int  zipped;
} VcfInfo;

static VcfInfo parse_vcf_filepath(const char *vcf_file_path) {
    VcfInfo info;
    memset(&info, 0, sizeof(info));

    const char *slash = strrchr(vcf_file_path, '/');
    const char *filename;
    if (slash) {
        size_t dirlen = (size_t)(slash - vcf_file_path);
        if (dirlen >= MAX_PATH_LEN) dirlen = MAX_PATH_LEN - 1;
        memcpy(info.path, vcf_file_path, dirlen);
        info.path[dirlen] = '\0';
        filename = slash + 1;
    } else {
        info.path[0] = '.';
        info.path[1] = '\0';
        filename = vcf_file_path;
    }

    size_t flen = strlen(filename);
    if (flen > 7 && strcmp(filename + flen - 7, ".vcf.gz") == 0) {
        info.zipped = 1;
        size_t plen = flen - 7;
        if (plen >= MAX_PATH_LEN) plen = MAX_PATH_LEN - 1;
        memcpy(info.prefix, filename, plen);
        info.prefix[plen] = '\0';
    } else if (flen > 4 && strcmp(filename + flen - 4, ".vcf") == 0) {
        info.zipped = 0;
        size_t plen = flen - 4;
        if (plen >= MAX_PATH_LEN) plen = MAX_PATH_LEN - 1;
        memcpy(info.prefix, filename, plen);
        info.prefix[plen] = '\0';
    } else {
        ERR("Unexpected file extension: %s. VCF file must end in .vcf or .vcf.gz", filename);
        exit(1);
    }

    DBG("parse_vcf_filepath: path='%s' prefix='%s' zipped=%d",
        info.path, info.prefix, info.zipped);
    return info;
}

/* ── Tokenizer ─────────────────────────────────────────────────────── */

/* Split `str` on '\t' in place.  Returns number of fields.
 * If max_fields > 0, last field gets the unsplit remainder. */
static int split_tabs(char *str, char **fields, int max_fields) {
    int n = 0;
    char *p = str;
    while (p && *p) {
        if (max_fields > 0 && n == max_fields - 1) {
            fields[n++] = p;
            break;
        }
        fields[n++] = p;
        char *tab = strchr(p, '\t');
        if (tab) { *tab = '\0'; p = tab + 1; }
        else break;
    }
    return n;
}

/* ── Main extraction logic ─────────────────────────────────────────── */

static void extract_tracts(const char *vcf_path, const char *msp_path,
                           int num_ancs, const char *output_dir,
                           int output_vcf, int compress_output,
                           int emit_freqs)
{
    int i, j;

    /* ── Validate inputs ─────────────────────────────────────────── */
    {
        FILE *f = fopen(vcf_path, "r");
        if (!f) { ERR("VCF path '%s' does not exist", vcf_path); exit(1); }
        fclose(f);
    }
    {
        FILE *f = fopen(msp_path, "r");
        if (!f) { ERR("MSP path '%s' does not exist", msp_path); exit(1); }
        fclose(f);
    }
    if (num_ancs < 1 || num_ancs > MAX_ANCS) {
        ERR("num_ancs must be between 1 and %d (got %d)", MAX_ANCS, num_ancs);
        exit(1);
    }

    VcfInfo vi = parse_vcf_filepath(vcf_path);

    INFO("# Running script version      : %s", VERSION);
    INFO("# VCF File                    : %s", vcf_path);
    INFO("# Prefix of output file names : %s", vi.prefix);
    INFO("# VCF File is compressed?     : %s", vi.zipped ? "True" : "False");
    INFO("# Number of Ancestries in VCF : %d", num_ancs);
    INFO("# Output Directory            : %s", output_dir ? output_dir : vi.path);
    if (compress_output)
        INFO("# Files will be gzip compressed (not bgzip)");

    /* ── Build output base path ──────────────────────────────────── */
    char output_base[MAX_PATH_LEN];
    if (output_dir)
        snprintf(output_base, MAX_PATH_LEN, "%s/%s", output_dir, vi.prefix);
    else
        snprintf(output_base, MAX_PATH_LEN, "%s/%s", vi.path, vi.prefix);

    const char *ext = compress_output ? ".gz" : "";

    /* ── Open output files (heap-allocated arrays) ───────────────── */
    OutFile *dos_files    = (OutFile *)calloc(num_ancs, sizeof(OutFile));
    OutFile *ancdos_files = (OutFile *)calloc(num_ancs, sizeof(OutFile));
    OutFile *vcf_files    = output_vcf ? (OutFile *)calloc(num_ancs, sizeof(OutFile)) : NULL;

    INFO("Creating output files for %d ancestries", num_ancs);
    for (i = 0; i < num_ancs; i++) {
        char fname[MAX_PATH_LEN];

        snprintf(fname, MAX_PATH_LEN, "%s.anc%d.dosage.txt%s", output_base, i, ext);
        DBG("Opening dosage file: %s", fname);
        dos_files[i] = outfile_open(fname, compress_output);

        snprintf(fname, MAX_PATH_LEN, "%s.anc%d.hapcount.txt%s", output_base, i, ext);
        DBG("Opening hapcount file: %s", fname);
        ancdos_files[i] = outfile_open(fname, compress_output);

        if (output_vcf) {
            snprintf(fname, MAX_PATH_LEN, "%s.anc%d.vcf%s", output_base, i, ext);
            DBG("Opening ancestry VCF file: %s", fname);
            vcf_files[i] = outfile_open(fname, compress_output);
        }
    }

    /* Per-ancestry frequency sidecar (AC/AN/MAC/MAF per variant) */
    OutFile freqs_file = {NULL, NULL, compress_output};
    char    freqs_path[MAX_PATH_LEN] = "";
    if (emit_freqs) {
        snprintf(freqs_path, MAX_PATH_LEN, "%s.anc_freqs.tsv%s", output_base, ext);
        DBG("Opening freqs sidecar: %s", freqs_path);
        freqs_file = outfile_open(freqs_path, compress_output);
        INFO("# Per-ancestry frequency sidecar: %s", freqs_path);
    }

    /* ── Open input files ────────────────────────────────────────── */
    InFile vcf_in = infile_open(vcf_path, vi.zipped);
    FILE  *msp_fp = fopen(msp_path, "r");
    if (!msp_fp) { ERR("Cannot open MSP file: %s", msp_path); exit(1); }

    /* ── Pre-read MSP header to extract sample names ────────────── */
    /*
     * MSP header line 2 looks like:
     *   #chm  spos  epos  sgpos  egpos  n snps  SAMPLE1.0  SAMPLE1.1  SAMPLE2.0  SAMPLE2.1
     */
    char **msp_sample_names = NULL;   /* unique sample names from MSP */
    int    n_msp_samples    = 0;      /* number of unique samples in MSP */
    {
        char *hdr_buf = (char *)malloc(MAX_LINE);
        if (!hdr_buf) { ERR("Failed to allocate MSP header buffer"); exit(1); }

        while (fgets(hdr_buf, MAX_LINE, msp_fp)) {
            size_t hlen = strlen(hdr_buf);
            while (hlen > 0 && (hdr_buf[hlen-1] == '\n' || hdr_buf[hlen-1] == '\r'))
                hdr_buf[--hlen] = '\0';

            if (hdr_buf[0] != '#') break;  /* past headers, stop */

            /* Look for the column header line (starts with #chm or has spos) */
            if (strstr(hdr_buf, "spos") == NULL) continue;

            /* Parse: skip first 6 fields, then haplotype columns */
            char *saved = strdup(hdr_buf);
            char *p = saved;
            while (*p == '#' || *p == ' ') p++;

            char *msp_hdr_fields[7];
            int mhf = split_tabs(p, msp_hdr_fields, 7);
            if (mhf < 7) {
                WARN("MSP header has fewer than 7 fields, cannot extract sample names");
                free(saved);
                break;
            }

            /* Field 6 has all haplotype columns tab-separated */
            char **hap_cols = (char **)malloc(sizeof(char *) * MAX_INDIV * 2);
            int n_hap_cols = 0;
            {
                char *hp = msp_hdr_fields[6];
                while (hp && *hp) {
                    hap_cols[n_hap_cols++] = hp;
                    char *t = strchr(hp, '\t');
                    if (t) { *t = '\0'; hp = t + 1; }
                    else break;
                    if (n_hap_cols >= MAX_INDIV * 2) break;
                }
            }

            /* Extract unique sample names (strip .0 / .1 suffix) */
            msp_sample_names = (char **)malloc(sizeof(char *) * (n_hap_cols / 2 + 1));
            for (int h = 0; h < n_hap_cols; h++) {
                char *name = strdup(hap_cols[h]);
                /* Strip trailing .0 or .1 */
                size_t nlen = strlen(name);
                if (nlen >= 2 && name[nlen-2] == '.' &&
                    (name[nlen-1] == '0' || name[nlen-1] == '1')) {
                    name[nlen-2] = '\0';
                }
                /* Only add if not same as previous (haplotypes come in pairs) */
                if (n_msp_samples == 0 || strcmp(name, msp_sample_names[n_msp_samples - 1]) != 0) {
                    msp_sample_names[n_msp_samples++] = name;
                } else {
                    free(name);
                }
            }

            free(hap_cols);
            free(saved);
            break;
        }

        /* Rewind MSP file so the main loop reads from the start */
        rewind(msp_fp);
        free(hdr_buf);

        if (n_msp_samples > 0) {
            INFO("Found %d samples in MSP file", n_msp_samples);
            DBG("First MSP sample: %s, last: %s",
                msp_sample_names[0], msp_sample_names[n_msp_samples - 1]);
        } else {
            WARN("Could not parse sample names from MSP header; all VCF samples will be kept");
        }
    }

    /*
     * Sample filtering arrays (populated when #CHROM line is processed):
     *   keep[i]        = 1 if VCF sample i is also in MSP, 0 otherwise
     *   vcf_to_msp[i]  = MSP sample index for VCF sample i (only valid if keep[i])
     */
    int *keep       = NULL;
    int *vcf_to_msp = NULL;
    int  n_vcf_samples = 0;
    int  n_kept        = 0;

    /* ── Heap-allocate all working buffers ────────────────────────── */
    char *vcf_line       = (char *)malloc(MAX_LINE);
    char *msp_line       = (char *)malloc(MAX_LINE);
    char *msp_line_saved = (char *)malloc(MAX_LINE);
    if (!vcf_line || !msp_line || !msp_line_saved) {
        ERR("Failed to allocate line buffers (need ~%d MB)", MAX_LINE * 3 / (1024*1024));
        exit(1);
    }

    /* VCF header accumulator */
    char  *vcf_header     = NULL;
    size_t vcf_header_len = 0;
    size_t vcf_header_cap = 0;

    /* Current MSP window */
    char win_chrom[256] = "";
    long win_spos = 0, win_epos = 0;

    /* MSP calls — pointers into msp_line_saved */
    char **calls  = (char **)malloc(sizeof(char *) * MAX_INDIV * 2);
    int    n_calls = 0;

    /* Per-ancestry string builders */
    StrBuf *dos_sb    = (StrBuf *)malloc(sizeof(StrBuf) * num_ancs);
    StrBuf *ancdos_sb = (StrBuf *)malloc(sizeof(StrBuf) * num_ancs);
    StrBuf *vcf_sb    = output_vcf ? (StrBuf *)malloc(sizeof(StrBuf) * num_ancs) : NULL;
    for (i = 0; i < num_ancs; i++) {
        strbuf_init(&dos_sb[i],    65536);
        strbuf_init(&ancdos_sb[i], 65536);
        if (output_vcf) strbuf_init(&vcf_sb[i], 65536);
    }

    /* Genos pointer array */
    char **genos = (char **)malloc(sizeof(char *) * MAX_INDIV);

    /* Per-variant per-ancestry frequency accumulators */
    long *ac_per_anc = (long *)calloc(num_ancs, sizeof(long));
    long *an_per_anc = (long *)calloc(num_ancs, sizeof(long));

    if (!calls || !dos_sb || !ancdos_sb || !genos || !ac_per_anc || !an_per_anc) {
        ERR("Failed to allocate working buffers"); exit(1);
    }

    long vcf_variant_count = 0;
    long skipped_count     = 0;

    INFO("Iterating through VCF file");

    /* ── Main VCF loop ───────────────────────────────────────────── */
    while (infile_gets(&vcf_in, vcf_line, MAX_LINE)) {
        /* Strip trailing newline/carriage return */
        size_t vlen = strlen(vcf_line);
        while (vlen > 0 && (vcf_line[vlen-1] == '\n' || vcf_line[vlen-1] == '\r'))
            vcf_line[--vlen] = '\0';

        /* ── Meta-header lines (##) ──────────────────────────────── */
        if (vlen >= 2 && vcf_line[0] == '#' && vcf_line[1] == '#') {
            if (output_vcf) {
                size_t need = vcf_header_len + vlen + 2;
                if (need > vcf_header_cap) {
                    vcf_header_cap = (vcf_header_cap == 0) ? (1 << 20) : vcf_header_cap * 2;
                    if (vcf_header_cap < need) vcf_header_cap = need;
                    vcf_header = (char *)realloc(vcf_header, vcf_header_cap);
                }
                memcpy(vcf_header + vcf_header_len, vcf_line, vlen);
                vcf_header_len += vlen;
                vcf_header[vcf_header_len++] = '\n';
                vcf_header[vcf_header_len]   = '\0';
            }
            continue;
        }

        /* ── #CHROM header line ──────────────────────────────────── */
        if (vcf_line[0] == '#') {
            DBG("Processing #CHROM header line");

            /*
             * Parse VCF sample names from field 9 onward.
             * Build keep[] and vcf_to_msp[] arrays to filter samples
             * that are in VCF but not in MSP.
             */
            char *hdr_copy = strdup(vcf_line);
            char *hp = hdr_copy;
            while (*hp == '#' || *hp == ' ') hp++;

            char *hfields[10];
            int hn = split_tabs(hp, hfields, 10);
            if (hn < 10) {
                ERR("VCF #CHROM line has fewer than 10 tab-separated fields (%d found)", hn);
                exit(1);
            }

            /* Parse individual sample names from field 9 */
            char **vcf_sample_names = (char **)malloc(sizeof(char *) * MAX_INDIV);
            n_vcf_samples = 0;
            {
                /* We need a separate copy since hfields[9] is in hdr_copy which we'll free */
                char *samples_copy = strdup(hfields[9]);
                char *sp = samples_copy;
                while (sp && *sp) {
                    vcf_sample_names[n_vcf_samples] = sp;
                    n_vcf_samples++;
                    char *t = strchr(sp, '\t');
                    if (t) { *t = '\0'; sp = t + 1; }
                    else break;
                    if (n_vcf_samples >= MAX_INDIV) break;
                }
                /* Duplicate each name so they survive after samples_copy is freed */
                for (i = 0; i < n_vcf_samples; i++)
                    vcf_sample_names[i] = strdup(vcf_sample_names[i]);
                free(samples_copy);
            }

            INFO("VCF has %d samples", n_vcf_samples);

            /* Build keep[] and vcf_to_msp[] */
            keep       = (int *)calloc(n_vcf_samples, sizeof(int));
            vcf_to_msp = (int *)calloc(n_vcf_samples, sizeof(int));
            n_kept = 0;

            if (n_msp_samples > 0) {
                /* Match VCF samples against MSP samples */
                for (i = 0; i < n_vcf_samples; i++) {
                    int found = 0;
                    for (int m = 0; m < n_msp_samples; m++) {
                        if (strcmp(vcf_sample_names[i], msp_sample_names[m]) == 0) {
                            keep[i] = 1;
                            vcf_to_msp[i] = m;
                            found = 1;
                            n_kept++;
                            break;
                        }
                    }
                    if (!found) {
                        // DBG("  VCF sample '%s' (index %d) NOT in MSP, dropping", vcf_sample_names[i], i);
                    }
                }
                INFO("Keeping %d of %d VCF samples (dropping %d not in MSP)",
                     n_kept, n_vcf_samples, n_vcf_samples - n_kept);
            } else {
                /* No MSP sample names parsed — keep all (original behavior) */
                // WARN("No MSP sample names available; keeping all %d VCF samples", n_vcf_samples);
                for (i = 0; i < n_vcf_samples; i++) {
                    keep[i] = 1;
                    vcf_to_msp[i] = i;
                    n_kept++;
                }
            }

            /* Build filtered header with only kept samples */
            StrBuf anc_hdr;
            strbuf_init(&anc_hdr, 65536);

            /* First 5 fixed columns */
            char *fixed5[5];
            {
                /* Re-parse from original vcf_line for clean field values */
                char *hdr2 = strdup(vcf_line);
                char *hp2 = hdr2;
                while (*hp2 == '#' || *hp2 == ' ') hp2++;
                char *hf2[10];
                split_tabs(hp2, hf2, 10);

                char prefix_buf[4096];
                int plen = snprintf(prefix_buf, sizeof(prefix_buf), "%s\t%s\t%s\t%s\t%s",
                                    hf2[0], hf2[1], hf2[2], hf2[3], hf2[4]);
                strbuf_append(&anc_hdr, prefix_buf, plen);
                (void)fixed5; /* unused, just keeping variable for clarity */
                free(hdr2);
            }

            /* Append only kept sample names */
            for (i = 0; i < n_vcf_samples; i++) {
                if (keep[i]) {
                    strbuf_append(&anc_hdr, "\t", 1);
                    strbuf_append(&anc_hdr, vcf_sample_names[i], strlen(vcf_sample_names[i]));
                }
            }
            strbuf_append(&anc_hdr, "\n", 1);

            DBG("Filtered anc_header length = %zu (%d samples)", anc_hdr.len, n_kept);

            /* Build filtered VCF header for output_vcf (only kept samples) */
            if (output_vcf) {
                /* Rebuild #CHROM line with only kept samples, append to vcf_header */
                char *hdr3 = strdup(vcf_line);
                char *hp3 = hdr3;
                /* Don't strip # for VCF header */
                char *hf3[10];
                split_tabs(hp3, hf3, 10);

                StrBuf vcf_chrom_line;
                strbuf_init(&vcf_chrom_line, 65536);
                /* Fixed 9 columns */
                {
                    char buf9[8192];
                    int blen = snprintf(buf9, sizeof(buf9), "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                                        hf3[0], hf3[1], hf3[2], hf3[3], hf3[4],
                                        hf3[5], hf3[6], hf3[7], hf3[8]);
                    strbuf_append(&vcf_chrom_line, buf9, blen);
                }
                for (i = 0; i < n_vcf_samples; i++) {
                    if (keep[i]) {
                        strbuf_append(&vcf_chrom_line, "\t", 1);
                        strbuf_append(&vcf_chrom_line, vcf_sample_names[i], strlen(vcf_sample_names[i]));
                    }
                }
                strbuf_append(&vcf_chrom_line, "\n", 1);

                /* Append to vcf_header */
                size_t need = vcf_header_len + vcf_chrom_line.len + 1;
                if (need > vcf_header_cap) {
                    vcf_header_cap = need * 2;
                    vcf_header = (char *)realloc(vcf_header, vcf_header_cap);
                }
                memcpy(vcf_header + vcf_header_len, vcf_chrom_line.buf, vcf_chrom_line.len);
                vcf_header_len += vcf_chrom_line.len;
                vcf_header[vcf_header_len] = '\0';

                strbuf_free(&vcf_chrom_line);
                free(hdr3);
            }

            for (i = 0; i < num_ancs; i++) {
                outfile_write(&dos_files[i],    anc_hdr.buf, (int)anc_hdr.len);
                outfile_write(&ancdos_files[i], anc_hdr.buf, (int)anc_hdr.len);
                if (output_vcf)
                    outfile_write(&vcf_files[i], vcf_header, (int)vcf_header_len);
            }

            /* Write freqs sidecar header: CHROM POS ID REF ALT [AC AN MAC MAF]_anc{j} */
            if (emit_freqs) {
                StrBuf freqs_hdr;
                strbuf_init(&freqs_hdr, 4096);
                strbuf_append(&freqs_hdr, "CHROM\tPOS\tID\tREF\tALT", 20);
                for (i = 0; i < num_ancs; i++) {
                    char col[128];
                    int clen = snprintf(col, sizeof(col),
                        "\tAC_anc%d\tAN_anc%d\tMAC_anc%d\tMAF_anc%d",
                        i, i, i, i);
                    strbuf_append(&freqs_hdr, col, clen);
                }
                strbuf_append(&freqs_hdr, "\n", 1);
                outfile_write(&freqs_file, freqs_hdr.buf, (int)freqs_hdr.len);
                strbuf_free(&freqs_hdr);
            }

            strbuf_free(&anc_hdr);
            /* Free sample name copies */
            for (i = 0; i < n_vcf_samples; i++) free(vcf_sample_names[i]);
            free(vcf_sample_names);
            free(hdr_copy);
            continue;
        }

        /* ── Data lines ──────────────────────────────────────────── */
        vcf_variant_count++;

        char *row[10];
        char *line_copy = strdup(vcf_line);
        int nf = split_tabs(line_copy, row, 10);
        if (nf < 10) {
            WARN("VCF data line has fewer than 10 fields, skipping (variant #%ld)", vcf_variant_count);
            free(line_copy);
            continue;
        }

        char chrom[256];
        strncpy(chrom, row[0], 255); chrom[255] = '\0';
        long pos = atol(row[1]);

        DBG("Variant #%ld: chrom=%s pos=%ld", vcf_variant_count, chrom, pos);

        /* dos_anc_out = fields 0-4 tab-joined */
        char dos_anc_out[4096];
        int dos_anc_out_len = snprintf(dos_anc_out, sizeof(dos_anc_out),
                                       "%s\t%s\t%s\t%s\t%s",
                                       row[0], row[1], row[2], row[3], row[4]);

        /* vcf_out_prefix = fields 0-8 with FORMAT="GT" */
        char vcf_out_prefix[8192];
        int vcf_out_prefix_len = snprintf(vcf_out_prefix, sizeof(vcf_out_prefix),
                 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT",
                 row[0], row[1], row[2], row[3], row[4],
                 row[5], row[6], row[7]);

        /* Parse genotypes from field 9 */
        int n_genos = 0;
        {
            char *p = row[9];
            while (p && *p) {
                genos[n_genos++] = p;
                char *t = strchr(p, '\t');
                if (t) { *t = '\0'; p = t + 1; }
                else break;
                if (n_genos >= MAX_INDIV) break;
            }
        }

        // DBG("  n_genos (individuals) = %d", n_genos);

        /* Initialize output buffers for each ancestry */
        int skip_line = 0;

        for (i = 0; i < num_ancs; i++) {
            strbuf_clear(&dos_sb[i]);
            strbuf_append(&dos_sb[i], dos_anc_out, dos_anc_out_len);
            strbuf_clear(&ancdos_sb[i]);
            strbuf_append(&ancdos_sb[i], dos_anc_out, dos_anc_out_len);
            if (output_vcf) {
                strbuf_clear(&vcf_sb[i]);
                strbuf_append(&vcf_sb[i], vcf_out_prefix, vcf_out_prefix_len);
            }
            ac_per_anc[i] = 0;
            an_per_anc[i] = 0;
        }

        /* ── Advance MSP window until it covers this position ────── */
        while (!(strcmp(chrom, win_chrom) == 0 && win_spos <= pos && pos < win_epos)) {
            /* If same chrom but window is already past pos, this VCF site has no window */
            if (strcmp(chrom, win_chrom) == 0 && win_spos > pos) {
                skip_line = 1;
                DBG("  VCF pos %ld before MSP window start %ld, skipping", pos, win_spos);
                break;
            }

            if (!fgets(msp_line, MAX_LINE, msp_fp)) {
                DBG("  Reached end of MSP file");
                break;
            }
            size_t mlen = strlen(msp_line);
            while (mlen > 0 && (msp_line[mlen-1] == '\n' || msp_line[mlen-1] == '\r'))
                msp_line[--mlen] = '\0';

            if (msp_line[0] == '#') continue;
            if (mlen == 0) continue;

            memcpy(msp_line_saved, msp_line, mlen + 1);

            /* chm, spos, epos, sgpos, egpos, nsnps, calls... */
            char *msp_fields[7];
            int mn = split_tabs(msp_line_saved, msp_fields, 7);
            if (mn < 7) {
                WARN("MSP line has fewer than 7 fields (%d), skipping", mn);
                continue;
            }

            strncpy(win_chrom, msp_fields[0], 255); win_chrom[255] = '\0';
            win_spos = atol(msp_fields[1]);
            win_epos = atol(msp_fields[2]);

            DBG("  New MSP window: chrom=%s spos=%ld epos=%ld",
                win_chrom, win_spos, win_epos);

            /* Parse calls from field 6 */
            n_calls = 0;
            {
                char *p = msp_fields[6];
                while (p && *p) {
                    calls[n_calls++] = p;
                    char *t = strchr(p, '\t');
                    if (t) { *t = '\0'; p = t + 1; }
                    else break;
                    if (n_calls >= MAX_INDIV * 2) break;
                }
            }
            DBG("  MSP calls parsed: %d (expect 2 * n_indiv)", n_calls);

            if (strcmp(chrom, win_chrom) == 0 && win_spos > pos) {
                skip_line = 1;
                DBG("  After advance, pos %ld still before window start %ld", pos, win_spos);
                break;
            }
        }

        if (skip_line) {
            // INFO("VCF position, %ld is not in an msp window, skipping site", pos);
            skipped_count++;
            free(line_copy);
            continue;
        }

        /* ── Per-individual, per-ancestry processing ─────────────── */
        for (i = 0; i < n_genos; i++) {
            /* Skip samples not present in MSP */
            if (keep && !keep[i]) continue;

            /* Parse genotype: take part before ':', split on '|' */
            char geno_buf[64];
            strncpy(geno_buf, genos[i], 63); geno_buf[63] = '\0';
            char *colon = strchr(geno_buf, ':');
            if (colon) *colon = '\0';

            char *pipe_pos = strchr(geno_buf, '|');
            char geno_a, geno_b;
            if (pipe_pos) {
                *pipe_pos = '\0';
                geno_a = geno_buf[0];
                geno_b = *(pipe_pos + 1);
            } else {
                geno_a = geno_buf[0];
                geno_b = (strlen(geno_buf) > 1) ? geno_buf[1] : '0';
            }

            /* Use MSP index for this sample (remapped from VCF index) */
            int msp_idx = (vcf_to_msp) ? vcf_to_msp[i] : i;
            int call_a_idx = 2 * msp_idx;
            int call_b_idx = 2 * msp_idx + 1;

            if (call_a_idx >= n_calls || call_b_idx >= n_calls) {
                if (vcf_variant_count <= 3) {
                    WARN("Individual %d (msp_idx %d): call index %d/%d out of range (n_calls=%d)",
                         i, msp_idx, call_a_idx, call_b_idx, n_calls);
                }
                break;
            }

            int call_a = atoi(calls[call_a_idx]);
            int call_b = atoi(calls[call_b_idx]);

            for (j = 0; j < num_ancs; j++) {
                int count     = 0;
                int anc_count = 0;

                int hap_a_here = (call_a == j);
                int hap_b_here = (call_b == j);

                /* Only count haplotypes with a valid (non-missing) allele.
                 * '.' = missing; excluded from both AC (count) and AN (anc_count). */
                if (hap_a_here && geno_a != '.') {
                    anc_count++;
                    if (geno_a == '1') count++;
                }
                if (hap_b_here && geno_b != '.') {
                    anc_count++;
                    if (geno_b == '1') count++;
                }

                ac_per_anc[j] += count;
                an_per_anc[j] += anc_count;

                char tmp[32];
                int tlen;

                tlen = snprintf(tmp, sizeof(tmp), "\t%d", count);
                strbuf_append(&dos_sb[j], tmp, tlen);

                tlen = snprintf(tmp, sizeof(tmp), "\t%d", anc_count);
                strbuf_append(&ancdos_sb[j], tmp, tlen);

                if (output_vcf) {
                    if (!hap_a_here && !hap_b_here) {
                        /* Neither haplotype on this ancestry = truly missing */
                        tlen = snprintf(tmp, sizeof(tmp), "\t.|.");
                    } else {
                        /* Haplotype present: use actual genotype (0 or 1)
                         * Haplotype absent:  use 0 (not missing, just not on this ancestry) */
                        char ga = hap_a_here ? geno_a : '0';
                        char gb = hap_b_here ? geno_b : '0';
                        tlen = snprintf(tmp, sizeof(tmp), "\t%c/%c", ga, gb);
                    }
                    strbuf_append(&vcf_sb[j], tmp, tlen);
                }
            }
        }

        /* ── Write output lines ──────────────────────────────────── */
        for (j = 0; j < num_ancs; j++) {
            strbuf_append(&dos_sb[j],    "\n", 1);
            strbuf_append(&ancdos_sb[j], "\n", 1);
            outfile_write(&dos_files[j],    dos_sb[j].buf,    (int)dos_sb[j].len);
            outfile_write(&ancdos_files[j], ancdos_sb[j].buf, (int)ancdos_sb[j].len);

            if (output_vcf) {
                strbuf_append(&vcf_sb[j], "\n", 1);
                outfile_write(&vcf_files[j], vcf_sb[j].buf, (int)vcf_sb[j].len);
            }
        }

        /* ── Write per-ancestry frequency sidecar row ──────────────
         * Columns per ancestry j:
         *   AC  = alt allele count on ancestry j (== Tractor's AC)
         *   AN  = total chromosomes (haplotypes) assigned to ancestry j
         *   MAC = min(AC, AN - AC)  (minor allele count, allele-agnostic)
         *   MAF = MAC / AN          ("NA" when AN == 0)
         */
        if (emit_freqs) {
            char freqs_row[16384];
            int  fr_off = 0;
            fr_off += snprintf(freqs_row + fr_off, sizeof(freqs_row) - fr_off,
                               "%s\t%s\t%s\t%s\t%s",
                               row[0], row[1], row[2], row[3], row[4]);
            for (j = 0; j < num_ancs; j++) {
                long AC  = ac_per_anc[j];
                long AN  = an_per_anc[j];
                long MAC = (AC < AN - AC) ? AC : (AN - AC);
                if (AN > 0) {
                    double MAF = (double)MAC / (double)AN;
                    fr_off += snprintf(freqs_row + fr_off, sizeof(freqs_row) - fr_off,
                                       "\t%ld\t%ld\t%ld\t%.6f", AC, AN, MAC, MAF);
                } else {
                    fr_off += snprintf(freqs_row + fr_off, sizeof(freqs_row) - fr_off,
                                       "\t%ld\t%ld\t%ld\tNA", AC, AN, MAC);
                }
            }
            freqs_row[fr_off++] = '\n';
            outfile_write(&freqs_file, freqs_row, fr_off);
        }

        if (vcf_variant_count % 100000 == 0) {
            INFO("Processed %ld variants so far...", vcf_variant_count);
        }

        free(line_copy);
    }

    /* ── Summary & cleanup ───────────────────────────────────────── */
    INFO("Finished extracting tracts per %d ancestries", num_ancs);
    INFO("Total variants processed: %ld", vcf_variant_count);
    INFO("Total variants skipped:   %ld", skipped_count);

    infile_close(&vcf_in);
    fclose(msp_fp);

    for (i = 0; i < num_ancs; i++) {
        outfile_close(&dos_files[i]);
        outfile_close(&ancdos_files[i]);
        if (output_vcf) outfile_close(&vcf_files[i]);
        strbuf_free(&dos_sb[i]);
        strbuf_free(&ancdos_sb[i]);
        if (output_vcf) strbuf_free(&vcf_sb[i]);
    }

    if (emit_freqs) outfile_close(&freqs_file);

    if (vcf_variant_count == 0) {
        WARN("No variants passed filtering in '%s' — output files contain headers only, removing", vcf_path);

        char fname[MAX_PATH_LEN];
        for (i = 0; i < num_ancs; i++) {
            snprintf(fname, MAX_PATH_LEN, "%s.anc%d.dosage.txt%s", output_base, i, ext);
            if (remove(fname) != 0)
                WARN("Could not remove %s: %s", fname, strerror(errno));
            else
                INFO("Removed empty output: %s", fname);

            snprintf(fname, MAX_PATH_LEN, "%s.anc%d.hapcount.txt%s", output_base, i, ext);
            if (remove(fname) != 0)
                WARN("Could not remove %s: %s", fname, strerror(errno));
            else
                INFO("Removed empty output: %s", fname);

            if (output_vcf) {
                snprintf(fname, MAX_PATH_LEN, "%s.anc%d.vcf%s", output_base, i, ext);
                if (remove(fname) != 0)
                    WARN("Could not remove %s: %s", fname, strerror(errno));
                else
                    INFO("Removed empty output: %s", fname);
            }
        }

        if (emit_freqs && freqs_path[0] != '\0') {
            if (remove(freqs_path) != 0)
                WARN("Could not remove %s: %s", freqs_path, strerror(errno));
            else
                INFO("Removed empty output: %s", freqs_path);
        }
    }
    
    free(dos_files);
    free(ancdos_files);
    if (vcf_files)  free(vcf_files);
    free(dos_sb);
    free(ancdos_sb);
    if (vcf_sb)     free(vcf_sb);
    free(genos);
    free(vcf_line);
    free(msp_line);
    free(msp_line_saved);
    free(calls);
    free(ac_per_anc);
    free(an_per_anc);
    if (vcf_header) free(vcf_header);
    if (keep)       free(keep);
    if (vcf_to_msp) free(vcf_to_msp);
    if (msp_sample_names) {
        for (i = 0; i < n_msp_samples; i++) free(msp_sample_names[i]);
        free(msp_sample_names);
    }
}

/* ── CLI ─────────────────────────────────────────────────────────────── */

static void usage(const char *progname) {
    fprintf(stderr,
        "Usage: %s --vcf <file> --msp <file> --num-ancs <int>\n"
        "          [--output-dir <dir>] [--output-vcf] [--compress-output] [--no-freqs]\n"
        "\n"
        "Options:\n"
        "  --vcf              Path to phased VCF file (*.vcf or *.vcf.gz)\n"
        "  --msp              Path to ancestry calls MSP file\n"
        "  --num-ancs         Number of ancestral populations\n"
        "  --output-dir       Output directory (default: same as VCF)\n"
        "  --output-vcf       Generate ancestry-specific VCF files\n"
        "  --compress-output  gzip compress output files\n"
        "  --no-freqs         Suppress the <prefix>.anc_freqs.tsv sidecar\n"
        "                     (default: emit per-variant AC/AN/MAC/MAF per ancestry)\n"
        "\n"
        "Compile with -DDEBUG=1 for verbose debug output.\n",
        progname);
}

int main(int argc, char **argv) {
    const char *vcf_path   = NULL;
    const char *msp_path   = NULL;
    int         num_ancs   = -1;
    const char *output_dir = NULL;
    int         output_vcf_flag = 0;
    int         compress_flag   = 0;
    int         emit_freqs_flag = 1;  /* sidecar on by default */

    static struct option long_opts[] = {
        {"vcf",             required_argument, NULL, 'v'},
        {"msp",             required_argument, NULL, 'm'},
        {"num-ancs",        required_argument, NULL, 'n'},
        {"output-dir",      required_argument, NULL, 'd'},
        {"output-vcf",      no_argument,       NULL, 'o'},
        {"compress-output", no_argument,       NULL, 'z'},
        {"no-freqs",        no_argument,       NULL, 'F'},
        {"help",            no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "v:m:n:d:ozFh", long_opts, NULL)) != -1) {
        switch (c) {
            case 'v': vcf_path   = optarg; break;
            case 'm': msp_path   = optarg; break;
            case 'n': num_ancs   = atoi(optarg); break;
            case 'd': output_dir = optarg; break;
            case 'o': output_vcf_flag = 1; break;
            case 'z': compress_flag   = 1; break;
            case 'F': emit_freqs_flag = 0; break;
            case 'h':
            default:
                usage(argv[0]);
                return (c == 'h') ? 0 : 1;
        }
    }

    if (!vcf_path || !msp_path || num_ancs < 1) {
        ERR("Missing required arguments (--vcf, --msp, --num-ancs)");
        usage(argv[0]);
        return 1;
    }

    DBG("CLI parsed: vcf=%s msp=%s num_ancs=%d output_dir=%s output_vcf=%d compress=%d emit_freqs=%d",
        vcf_path, msp_path, num_ancs,
        output_dir ? output_dir : "(default)",
        output_vcf_flag, compress_flag, emit_freqs_flag);

    extract_tracts(vcf_path, msp_path, num_ancs, output_dir,
                   output_vcf_flag, compress_flag, emit_freqs_flag);

    return 0;
}