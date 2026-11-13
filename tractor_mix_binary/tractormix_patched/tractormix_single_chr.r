#!/usr/bin/env Rscript
#
# Tractor-Mix association for ONE chromosome, against a GMMAT null model from
# 02_fit_null_binary.r. Submitted per piece (or chromosome) and scan by 07_submit_tractormix.py.
#
# Input: extract_tracts output (05, split by 06) for one piece or chromosome, one file per ancestry k,
# given as templates in which {anc} is replaced by 0, 1, ..., n_anc - 1:
#   --dosage    <...>.anc{anc}.dosage.txt     ALT alleles on ancestry-k haplotypes
#   --hapcount  <...>.anc{anc}.hapcount.txt   ancestry-k haplotypes (cond only)
# Plain or .gz (gz needs the R.utils package for data.table::fread).
#
# Tests (--test):
#   uncond  TractorMix.score       Y* ~ X + G_anc0 + ... + G_anc{K-1}    (primary)
#   cond    TractorMix.score_cond  Y* ~ X + LA + G_anc0 + ... + G_anc{K-1}
#           with LA = K-1 hapcounts. Which ancestry is left out (--la-drop,
#           default the last) does not change the results: the K hapcounts sum
#           to 2, so with the intercept in X any K-1 of them span the same space.
#
# Output (Tractor-Mix's own columns):
#   CHR POS ID REF ALT   variant, copied from the first dosage file
#   Chi2, P              joint K-df test that the effect is 0 in EVERY ancestry;
#                        NOT a heterogeneity test. NA if any ancestry was dropped
#   Eff_anc{k}, SE_anc{k}, Pval_anc{k}
#                        effect per ALT allele on an ancestry-k haplotype (Y*
#                        units), estimated jointly across ancestries
#   AC_count{k}          ALT alleles on ancestry-k haplotypes, null-model samples
#   include_anc{k}       1 if AC_count{k} > --ac-threshold; 0 = ancestry dropped
# Tractor-Mix rounds Eff, SE and Chi2 to 5 decimals.
#
# CHECKS BEFORE RUNNING
# ---------------------
# Tractor-Mix assumes the following without checking, so they are checked here:
#   - every input file has the same header (same samples, same order) and the
#     same variant rows (CHROM POS ID REF ALT, same order). Tractor-Mix takes the
#     variant annotation from the first dosage file only, compares the headers of
#     the dosage files only, and reads all files chunk by chunk in step.
#   - every null-model sample is in the files (Tractor-Mix does stop on this, but
#     only after setting up).
#   - chunking: Tractor-Mix fails on a chunk holding exactly one variant (a
#     one-row matrix drops to a vector and abind(along = 3) errors), so
#     --chunk-size is lowered when the variant count leaves a remainder of 1. A
#     file holding a single variant cannot be run.
#
# Tractor-Mix appends results chunk by chunk, so it writes to <out>.tmp, which is
# renamed to <out> only after its row count and variant order are verified. An
# existing <out> is therefore always complete.
#
# --cores is Tractor-Mix's n_core: it starts n_core - 1 worker processes when
# n_core > 1. Each worker receives the null model (P is N x N) and the current
# chunk, so memory grows with --cores and --chunk-size.
#
# Usage:
#   Rscript tractormix_patched/tractormix_single_chr.r \
#     --null-model gmmat_null/gmmat_null_N0.rds \
#     --dosage     tracts/chr22.anc{anc}.dosage.txt \
#     --hapcount   tracts/chr22.anc{anc}.hapcount.txt \
#     --test       uncond \
#     --out        assoc_tm/N0_uncond_chr22.tsv \
#     --cores      4          # default: $SLURM_CPUS_PER_TASK, else 1
#
suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0) stop("Arguments must be --key value pairs.")
  if (length(a)) {
    for (i in seq(1, length(a), by = 2)) {
      key <- gsub("-", "_", sub("^--", "", a[i]))
      if (!key %in% names(defaults)) {
        stop("Unknown argument ", a[i], ". Known: ",
             paste0("--", gsub("_", "-", names(defaults)), collapse = ", "))
      }
      defaults[[key]] <- a[i + 1]
    }
  }
  defaults
}

opt <- parse_args(list(
  null_model     = "gmmat_null/gmmat_null_N0.rds",
  dosage         = "tracts/chr22.anc{anc}.dosage.txt",
  hapcount       = "tracts/chr22.anc{anc}.hapcount.txt",
  test           = "uncond",
  out            = "assoc_tm/N0_uncond_chr22.tsv",
  n_anc          = "3",
  anc_labels     = "NAT,AFR,EUR",
  la_drop        = "",          # ancestry index left out of LA (cond); default: last
  ac_threshold   = "50",
  chunk_size     = "2048",
  cores          = Sys.getenv("SLURM_CPUS_PER_TASK", "1"),
  tractormix_dir = "."
))

msg <- function(...) cat(..., "\n", sep = "")
msg(">>> ", opt$test, "  ", opt$dosage, "  ", format(Sys.time()))

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if (!opt$test %in% c("uncond", "cond")) stop("--test must be 'uncond' or 'cond'.")
K <- as.integer(opt$n_anc)
if (is.na(K) || K < 2) stop("--n-anc must be an integer >= 2.")
labels <- strsplit(opt$anc_labels, ",", fixed = TRUE)[[1]]
if (length(labels) != K) stop("--anc-labels has ", length(labels), " labels for --n-anc ", K, ".")
ac_threshold <- as.integer(opt$ac_threshold)
if (is.na(ac_threshold) || ac_threshold < 0) stop("--ac-threshold must be a non-negative integer.")
chunk_size <- as.integer(opt$chunk_size)
if (is.na(chunk_size) || chunk_size < 2) stop("--chunk-size must be an integer >= 2.")
cores <- as.integer(opt$cores)
if (is.na(cores) || cores < 1) cores <- 1L
la_drop <- if (nzchar(opt$la_drop)) as.integer(opt$la_drop) else K - 1L
if (is.na(la_drop) || !la_drop %in% 0:(K - 1)) stop("--la-drop must be in 0..", K - 1, ".")

anc_files <- function(template, what) {
  if (!grepl("{anc}", template, fixed = TRUE)) {
    stop("--", what, " must contain {anc}, e.g. tracts/chr22.anc{anc}.", what, ".txt")
  }
  vapply(0:(K - 1), function(k) gsub("{anc}", k, template, fixed = TRUE), character(1))
}
dos_files <- anc_files(opt$dosage, "dosage")
la_files  <- if (opt$test == "cond") anc_files(opt$hapcount, "hapcount")[-(la_drop + 1)] else character(0)
inputs    <- c(dos_files, la_files)

tm_script <- file.path(opt$tractormix_dir,
                       if (opt$test == "uncond") "TractorMix.score.R" else "TractorMix.score_cond.R")
absent <- c(opt$null_model, tm_script, inputs)[!file.exists(c(opt$null_model, tm_script, inputs))]
if (length(absent)) stop("Not found: ", paste(absent, collapse = ", "))
if (any(endsWith(inputs, ".gz")) && !requireNamespace("R.utils", quietly = TRUE)) {
  stop("Compressed inputs need the R.utils package (data.table::fread reads .gz through it).")
}

msg("dosage files:   ", paste(basename(dos_files), collapse = ", "))
if (length(la_files)) {
  msg("hapcount files: ", paste(basename(la_files), collapse = ", "),
      "  (", labels[la_drop + 1], " left out)")
}
msg("ancestries: ", paste0("anc", 0:(K - 1), "=", labels, collapse = ", "))

# ---------------------------------------------------------------------------
# Input consistency across files
# ---------------------------------------------------------------------------
headers <- lapply(inputs, function(f) names(fread(f, sep = "\t", nrows = 0)))
bad_hdr <- inputs[!vapply(headers, identical, logical(1), headers[[1]])]
if (length(bad_hdr)) {
  stop("Header (samples or their order) differs from ", basename(inputs[1]), " in: ",
       paste(basename(bad_hdr), collapse = ", "))
}
lead <- headers[[1]][1:5]
if (!identical(lead, c("CHROM", "POS", "ID", "REF", "ALT"))) {
  msg("note: leading columns are ", paste(lead, collapse = " "),
      "; expected CHROM POS ID REF ALT (extract_tracts layout)")
}
samples <- headers[[1]][-(1:5)]
if (anyDuplicated(samples)) {
  stop("Duplicate sample IDs in the header, e.g. ", samples[duplicated(samples)][1])
}

meta <- lapply(inputs, function(f) fread(f, sep = "\t", select = 1:5, colClasses = "character"))
n_var <- nrow(meta[[1]])
same_rows <- vapply(meta, function(m) {
  nrow(m) == n_var && all(vapply(1:5, function(j) identical(m[[j]], meta[[1]][[j]]), logical(1)))
}, logical(1))
if (!all(same_rows)) {
  stop("Variant rows (CHROM POS ID REF ALT, in order) differ from ", basename(inputs[1]),
       " in: ", paste(basename(inputs[!same_rows]), collapse = ", "))
}
msg("variants: ", n_var, " | samples in files: ", length(samples))
if (n_var == 0) stop("No variants in the input files.")
if (n_var == 1) {
  stop("A single variant cannot be run: Tractor-Mix fails on a one-variant chunk.")
}

# ---------------------------------------------------------------------------
# Null model
# ---------------------------------------------------------------------------
obj <- readRDS(opt$null_model)
if (!inherits(obj, "glmmkin")) stop(opt$null_model, " is not a glmmkin object.")
null_ids <- as.character(obj$id_include)
not_in_files <- setdiff(null_ids, samples)
if (length(not_in_files)) {
  stop(length(not_in_files), " null-model sample(s) absent from the input files, e.g. ",
       paste(head(not_in_files, 5), collapse = ", "))
}
msg("null model: ", opt$null_model, "  (", length(null_ids), " samples; ",
    if (is.null(obj$P)) "sparse path: Sigma_i / Sigma_iX" else "dense path: P", ")")
extra <- length(samples) - length(null_ids)
if (extra > 0) msg("note: ", extra, " sample(s) in the files but not in the null model; ignored")

# ---------------------------------------------------------------------------
# Chunk size
# ---------------------------------------------------------------------------
if (n_var %% chunk_size == 1) {
  cs <- chunk_size
  while (cs >= 2 && n_var %% cs == 1) cs <- cs - 1L
  if (cs < 2) stop("No chunk size >= 2 avoids a one-variant chunk for ", n_var, " variants.")
  msg("chunk size ", chunk_size, " would leave a one-variant last chunk; using ", cs)
  chunk_size <- cs
}
msg("chunk size ", chunk_size, " (", ceiling(n_var / chunk_size), " chunks) | n_core ", cores,
    " | AC threshold ", ac_threshold)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
tmp <- paste0(opt$out, ".tmp")
unlink(tmp)

source(tm_script)
t0 <- Sys.time()
if (opt$test == "uncond") {
  TractorMix.score(obj = obj, infiles = dos_files, AC_threshold = ac_threshold,
                   outfiles = tmp, n_core = cores, chunk_size = chunk_size)
} else {
  TractorMix.score_cond(obj = obj, infiles_geno = dos_files, infiles_la = la_files,
                        AC_threshold = ac_threshold, outfiles = tmp,
                        n_core = cores, chunk_size = chunk_size)
}
msg("\nTractor-Mix finished in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))

# ---------------------------------------------------------------------------
# Verify, summarize, finalize
# ---------------------------------------------------------------------------
res <- fread(tmp, sep = "\t", colClasses = list(character = c("CHR", "POS", "ID", "REF", "ALT")))
if (nrow(res) != n_var) {
  stop("Output has ", nrow(res), " rows for ", n_var, " variants; left as ", tmp)
}
if (!identical(res$ID, meta[[1]]$ID) || !identical(res$POS, meta[[1]]$POS)) {
  stop("Output variant order differs from the input; left as ", tmp)
}

k <- 0:(K - 1)
inc <- as.matrix(res[, paste0("include_anc", k), with = FALSE])
n_inc <- rowSums(inc == 1, na.rm = TRUE)
all_in <- n_inc == K
msg("\nvariants with all ancestries included (AC > ", ac_threshold, "): ", sum(all_in),
    " | some dropped: ", sum(n_inc > 0 & n_inc < K), " | none: ", sum(n_inc == 0))
msg("all included but no result (singular variance): ", sum(all_in & is.na(res$Chi2)))

lambda <- function(chi2, df) median(chi2, na.rm = TRUE) / qchisq(0.5, df)
if (any(all_in & !is.na(res$Chi2))) {
  msg("joint test (", K, " df), this chromosome only: lambda = ",
      sprintf("%.4f", lambda(res$Chi2[all_in], K)),
      ", min P = ", format(min(res$P[all_in], na.rm = TRUE), digits = 4))
}
for (i in seq_along(k)) {
  p <- res[[paste0("Pval_anc", k[i])]]
  if (any(!is.na(p))) {
    msg(sprintf("  anc%d %-4s  %8d tested  lambda = %.4f  min p = %s", k[i], labels[i],
                sum(!is.na(p)), lambda(qchisq(p, 1, lower.tail = FALSE), 1),
                format(min(p, na.rm = TRUE), digits = 4)))
  }
}

if (!file.rename(tmp, opt$out)) stop("Could not rename ", tmp, " to ", opt$out)
msg("wrote ", opt$out, "  (", nrow(res), " variants)")
msg("<<< done  ", format(Sys.time()))
