#!/usr/bin/env Rscript
#
# Logistic mixed null model for the binary Tractor-Mix run, fitted with GMMAT (we could use PC-Relate but it is slower and we have way too many samples, if we ended up with weird results could be worth to do it)
#
# GRM: the dense GCTA GRM from 01_make_grm_gcta.sh, passed to GMMAT as a base R matrix. GMMAT then returns P, and Tractor-Mix uses its dense code branch
#
# Covariates: listed with their type in model_covariates.txt (--model), the same
# file 00_define_samples.py reads. Add a covariate there, not here.
#
# Coding:
#   STATUS        2 = case, 1 = control, recoded to 1 / 0 for the model
#   categorical   a factor; the first value in sorted order is the reference,
#                 so the coding (1/2, 0/1, M/F) does not change the fit
#   quantitative  a number, one linear coefficient
#   PCs           GCTA_PC1.. or PC1.., whichever the file has (--pc_prefix if both)
#
# The script stops if:
#   - the GRM samples are not exactly the samples in analysis_samples.txt
#   - an analysis sample is not in the covariate file or has a missing value
#   - STATUS is not 1/2
#   - a categorical covariate has a single value, or a value with no cases or no controls
#   - a quantitative covariate or PC is not a number, or has the same value for everyone
#   - the model does not converge, or GMMAT leaves any sample out
#
# Outputs:
#   <out>                    the null model (--null-model for Tractor-Mix)
#   <out minus .rds>_fixef.tsv  fixed effects (log-odds, SE, OR, Wald p)
#
# Usage:
#   Rscript 02_fit_null_binary.r \
#     --covar   covariates.tsv \
#     --model   model_covariates.txt \
#     --samples samples/analysis_samples.txt \
#     --grm     grm/grm \
#     --out     gmmat_null/gmmat_null_N0.rds
#
suppressPackageStartupMessages(library(GMMAT))

parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0) stop("Arguments must be --key value pairs.")
  for (i in seq_len(length(a) / 2) * 2 - 1) {
    key <- gsub("-", "_", sub("^--", "", a[i]))
    if (!key %in% names(defaults)) {
      stop("Unknown argument ", a[i], ". Known: ", paste0("--", names(defaults), collapse = ", "))
    }
    defaults[[key]] <- a[i + 1]
  }
  defaults
}

opt <- parse_args(list(
  covar     = "covariates.tsv",
  model     = "model_covariates.txt",
  samples   = "samples/analysis_samples.txt",
  grm       = "grm/grm",
  out       = "gmmat_null/gmmat_null_N0.rds",
  pc_prefix = ""
))

show <- function(x, n = 10) {
  paste0(paste(head(x, n), collapse = ", "), if (length(x) > n) paste0(" ... (+", length(x) - n, " more)"))
}

# model_covariates.txt: one "column type" pair per line, plus "PCs <number>".
# Same rules as read_model() in 00_define_samples.py.
read_model <- function(path) {
  lines <- trimws(sub("#.*", "", readLines(path)))
  lines <- lines[nzchar(lines)]
  categorical <- character(0)
  quantitative <- character(0)
  n_pcs <- 0
  for (line in lines) {
    fields <- strsplit(line, "[[:space:]]+")[[1]]
    if (length(fields) != 2) stop(path, ": expected a column name and a type, got '", line, "'")
    name <- fields[1]
    kind <- fields[2]
    if (name %in% c(categorical, quantitative)) stop(path, ": ", name, " is listed twice")
    if (name == "PCs") {
      if (!grepl("^[0-9]+$", kind)) stop(path, ": PCs needs a number, got '", kind, "'")
      n_pcs <- as.integer(kind)
    } else if (name %in% c("IID", "STATUS")) {
      stop(path, ": ", name, " is not a covariate")
    } else if (kind == "categorical") {
      categorical <- c(categorical, name)
    } else if (kind == "quantitative") {
      quantitative <- c(quantitative, name)
    } else {
      stop(path, ": type of ", name, " must be categorical or quantitative, got '", kind, "'")
    }
  }
  if (sum(grepl("^PCs[[:space:]]", lines)) > 1) stop(path, ": PCs is given twice")
  list(categorical = categorical, quantitative = quantitative, n_pcs = n_pcs)
}

# ---------------------------------------------------------------------------
# GRM
# ---------------------------------------------------------------------------
# .grm.bin holds the lower triangle with the diagonal, row by row, as 4-byte
# floats: (1,1), (2,1), (2,2), (3,1), ... The order of the samples is the order
# of the lines in .grm.id. Read row by row, the lower triangle is the upper
# triangle read column by column, which is how R fills upper.tri().
read_gcta_grm <- function(prefix) {
  id_file  <- paste0(prefix, ".grm.id")
  bin_file <- paste0(prefix, ".grm.bin")
  ids <- read.table(id_file, colClasses = "character")[[2]]
  n <- length(ids)
  if (anyDuplicated(ids)) stop("Duplicate IIDs in ", id_file, ": ", show(ids[duplicated(ids)]))

  n_values <- n * (n + 1) / 2
  if (file.size(bin_file) != 4 * n_values) {
    stop(bin_file, " has ", file.size(bin_file), " bytes; ", n, " samples in ", id_file,
         " need ", 4 * n_values, ". The two files do not belong together.")
  }
  values <- readBin(bin_file, what = "numeric", n = n_values, size = 4, endian = "little")

  grm <- matrix(0, n, n)
  grm[upper.tri(grm, diag = TRUE)] <- values
  grm[lower.tri(grm)] <- t(grm)[lower.tri(grm)]
  dimnames(grm) <- list(ids, ids)
  grm
}

samples <- readLines(opt$samples)
samples <- samples[nzchar(samples)]
cat("analysis samples:", length(samples), "\n")

grm <- read_gcta_grm(opt$grm)
cat("GRM:", nrow(grm), "x", ncol(grm), "\n")
only_grm  <- setdiff(rownames(grm), samples)
only_list <- setdiff(samples, rownames(grm))
if (length(only_grm) || length(only_list)) {
  stop("GRM samples differ from ", opt$samples, "\n",
       "  only in the GRM:  ", show(only_grm), "\n",
       "  only in the list: ", show(only_list))
}
off <- grm[lower.tri(grm)]
cat("GRM diagonal:     ", paste(names(summary(diag(grm))), signif(summary(diag(grm)), 4), collapse = "  "), "\n")
cat("GRM off-diagonal: ", paste(names(summary(off)), signif(summary(off), 4), collapse = "  "), "\n")
cat("pairs above 0.125 / 0.25 / 0.4:", sum(off > 0.125), "/", sum(off > 0.25), "/", sum(off > 0.4), "\n")
rm(off)

# ---------------------------------------------------------------------------
# Covariates, in GRM order
# ---------------------------------------------------------------------------
model <- read_model(opt$model)

covar <- read.delim(opt$covar, colClasses = "character", check.names = FALSE,
                    na.strings = c("", "NA", "na", "nan", "NaN", "."))
for (col in c("IID", "STATUS", model$categorical, model$quantitative)) {
  if (!col %in% names(covar)) stop(opt$covar, " has no ", col, " column.")
}
if (anyDuplicated(covar$IID)) stop("Duplicate IIDs in ", opt$covar, ": ", show(covar$IID[duplicated(covar$IID)]))

pc_cols <- character(0)
if (model$n_pcs > 0) {
  prefixes <- if (nzchar(opt$pc_prefix)) opt$pc_prefix else c("GCTA_PC", "PC")
  pc_sets <- lapply(prefixes, function(p) paste0(p, seq_len(model$n_pcs)))
  found <- Filter(function(cols) all(cols %in% names(covar)), pc_sets)
  if (length(found) == 0) stop("No complete set of ", model$n_pcs, " PCs named ", paste0(prefixes, "1..", collapse = " or "))
  if (length(found) > 1) stop("Both GCTA_PC1.. and PC1.. are in the covariate file; choose with --pc_prefix.")
  pc_cols <- found[[1]]
}
numeric_cols <- c(model$quantitative, pc_cols)

cat("\nmodel (", opt$model, "):\n", sep = "")
for (col in model$categorical)  cat("  ", col, "categorical\n")
for (col in model$quantitative) cat("  ", col, "quantitative\n")
cat("   PCs", if (length(pc_cols)) paste0(pc_cols[1], "..", tail(pc_cols, 1)) else "none", "\n\n")

not_in_covar <- setdiff(rownames(grm), covar$IID)
if (length(not_in_covar)) stop(length(not_in_covar), " GRM samples are not in the covariate file: ", show(not_in_covar))

covar <- covar[match(rownames(grm), covar$IID), c("IID", "STATUS", model$categorical, numeric_cols)]
stopifnot(identical(covar$IID, rownames(grm)))

# Missing values: stop and list them rather than let GMMAT drop the rows.
na_cells <- which(is.na(covar), arr.ind = TRUE)
if (nrow(na_cells)) {
  stop(nrow(na_cells), " missing value(s): ",
       show(paste0(covar$IID[na_cells[, 1]], " (", names(covar)[na_cells[, 2]], ")"), 50))
}

bad_status <- !covar$STATUS %in% c("1", "2")
if (any(bad_status)) {
  stop("STATUS must be 1 (control) or 2 (case): ", show(paste0(covar$IID[bad_status], " (", covar$STATUS[bad_status], ")"), 50))
}

dat <- data.frame(IID = covar$IID,
                  CASE = as.numeric(covar$STATUS == "2"),
                  stringsAsFactors = FALSE)
cat(sum(dat$CASE == 1), " cases, ", sum(dat$CASE == 0), " controls (",
    sprintf("%.1f", 100 * mean(dat$CASE)), "% cases)\n", sep = "")

# Categorical: a factor needs at least two values, and every value needs both
# cases and controls, otherwise its coefficient cannot be estimated.
for (col in model$categorical) {
  dat[[col]] <- factor(covar[[col]])
  if (nlevels(dat[[col]]) < 2) stop(col, " has a single value (", levels(dat[[col]]), "); it cannot be a covariate.")
  counts <- table(dat[[col]], factor(dat$CASE, levels = c(1, 0), labels = c("cases", "controls")))
  names(dimnames(counts)) <- c(col, "STATUS")
  cat("\n")
  print(counts)
  cat("reference value:", levels(dat[[col]])[1], "\n")
  empty <- rownames(counts)[counts[, "cases"] == 0 | counts[, "controls"] == 0]
  if (length(empty)) stop(col, " value(s) with no cases or no controls: ", paste(empty, collapse = ", "))
}

# Quantitative and PCs: must be numbers, and must vary.
for (col in numeric_cols) {
  dat[[col]] <- suppressWarnings(as.numeric(covar[[col]]))
  if (anyNA(dat[[col]])) stop(col, " has values that are not numbers, e.g. ", show(covar[[col]][is.na(dat[[col]])]))
  if (length(unique(dat[[col]])) == 1) stop(col, " has the same value (", dat[[col]][1], ") for every sample.")
}
for (col in model$quantitative) {
  cat(sprintf("\n%s: mean %.4g in cases, %.4g in controls; range %.4g to %.4g\n", col,
              mean(dat[[col]][dat$CASE == 1]), mean(dat[[col]][dat$CASE == 0]),
              min(dat[[col]]), max(dat[[col]])))
}

# ---------------------------------------------------------------------------
# Null model
# ---------------------------------------------------------------------------
model_terms <- c(model$categorical, numeric_cols)
if (length(model_terms) == 0) model_terms <- "1"
model_formula <- reformulate(model_terms, response = "CASE")
cat("\nmodel:", deparse(model_formula, width.cutoff = 500), "+ random effect (GRM)\n\n")

null <- glmmkin(fixed = model_formula, data = dat, kins = grm, id = "IID",
                family = binomial(link = "logit"),
                method = "REML", method.optim = "AI", verbose = TRUE)
                #these two arguments method and method.optim are the defaults, just typed out

if (!isTRUE(null$converged)) stop("The null model did not converge; nothing saved.")

used <- as.character(null$id_include)
if (length(used) != nrow(dat) || !setequal(used, dat$IID)) {
  stop("GMMAT used ", length(used), " samples, expected ", nrow(dat), ". Left out: ", show(setdiff(dat$IID, used)))
}
if (!grepl("binomial", paste(deparse(null$call), collapse = " "))) {
  stop("'binomial' is not in the stored call; Tractor-Mix would not recognise the family.")
}
if (is.null(null$P)) stop("GMMAT did not return P; was the GRM passed as a Matrix object?")

tau <- null$theta[["kins1"]]
cat("\nvariance component of the GRM (tau):", signif(tau, 5), "\n")
if (tau < 1e-4) {
  cat("NOTE: tau is ~0. The GRM explains almost no variance on the liability scale,",
      "so the model is close to a plain logistic regression.\n")
}

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
est <- null$coefficients
se  <- sqrt(diag(as.matrix(null$cov)))
fixef <- data.frame(term = names(est), Est = est, SE = se, OR = exp(est),
                    P = 2 * pnorm(-abs(est / se)), row.names = NULL)
cat("\nfixed effects (log-odds):\n")
print(fixef, digits = 4, row.names = FALSE)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(null, opt$out)
fixef_file <- paste0(sub("\\.rds$", "", opt$out), "_fixef.tsv")
write.table(fixef, fixef_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwrote", opt$out, "\nwrote", fixef_file, "\n")
