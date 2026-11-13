#!/usr/bin/env python3
"""
Define the analysis samples for the binary Tractor-Mix run.

A sample is kept when it is in all of these:
  - the covariate file (IID, STATUS and the covariates in model_covariates.txt)
  - the typed genotype file (.psam), used for the GRM
  - the gnomix .msp of every chromosome
  - the VCF of every chromosome (optional, but recommended)

Everything downstream (GRM, null model, VCF subset, msp) is built from the list
written here, so they all hold the same samples.

The covariates come from model_covariates.txt (--model), the same file
02_fit_null_binary.r reads. Add covariates there to be consistent (and not depend on hard coding always).

The script stops, listing the samples involved, when a kept sample has:
  - a missing value in STATUS or in any model covariate
  - a STATUS other than 1 (control) or 2 (case)
  - a quantitative covariate or PC that is not a number
It also stops when the .msp (or VCF) sample sets differ between chromosomes.

PCs are read as GCTA_PC1..GCTA_PC<n> or PC1..PC<n>, whichever the file has.
If it has both, choose one with --pc-prefix.

Outputs in --outdir:
  analysis_samples.txt  one IID per line, in .psam order. This is the format
                        plink2 --keep (01), bcftools -S (03) and 04_prep_msp.py --samples read.
  sample_report.tsv     every sample seen in any source, where it was found,
                        and whether it was kept

Usage:
  python3 00_define_samples.py \\
      --covar        covariates.tsv \\
      --model        model_covariates.txt \\
      --psam         typed/merged_plink2.onlyTyped.allvarids.psam \\
      --msp-template 'gnomix/LPD_chrom{chrom}/query_results.msp' \\
      --vcf-template 'imputed/chr{chrom}.dose.vcf.gz' \\
      --outdir       samples
"""

import argparse
import gzip
import os
import sys
from collections import Counter

MISSING = {"", "NA", "na", "nan", "NaN", "."}


def parse_chrs(spec):
    chroms = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-")
            chroms.extend(range(int(lo), int(hi) + 1))
        elif part:
            chroms.append(int(part))
    return sorted(set(chroms))


def show(ids, n=10):
    ids = sorted(ids)
    text = ", ".join(ids[:n])
    if len(ids) > n:
        text += f" ... (+{len(ids) - n} more)"
    return text


def stop(message):
    sys.exit(f"ERROR: {message}")


def find_duplicates(ids):
    return [i for i, n in Counter(ids).items() if n > 1]


def read_model(path):
    """Read model_covariates.txt: {column: type} and the number of PCs."""
    covariates = {}
    n_pcs = 0
    seen_pcs = False
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.split("#")[0].strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 2:
                stop(f"{path} line {lineno}: expected a column name and a type, got '{line}'")
            name, kind = fields
            if name == "PCs":
                if seen_pcs:
                    stop(f"{path} line {lineno}: PCs is given twice")
                if not kind.isdigit():
                    stop(f"{path} line {lineno}: PCs needs a number, got '{kind}'")
                n_pcs = int(kind)
                seen_pcs = True
            elif name in ("IID", "STATUS"):
                stop(f"{path} line {lineno}: {name} is not a covariate")
            elif kind not in ("categorical", "quantitative"):
                stop(f"{path} line {lineno}: type of {name} must be categorical or quantitative, got '{kind}'")
            elif name in covariates:
                stop(f"{path} line {lineno}: {name} is listed twice")
            else:
                covariates[name] = kind
    return covariates, n_pcs


def pick_pc_columns(header, n_pcs, prefix):
    """Return the PC column names to use, e.g. GCTA_PC1..GCTA_PC20."""
    if prefix:
        candidates = [prefix]
    else:
        candidates = ["GCTA_PC", "PC"]

    complete = []
    for p in candidates:
        cols = [f"{p}{i}" for i in range(1, n_pcs + 1)]
        if all(c in header for c in cols):
            complete.append(cols)

    if len(complete) == 1:
        return complete[0]
    if len(complete) > 1:
        stop("the covariate file has both GCTA_PC1.. and PC1.. columns; "
             "choose one with --pc-prefix GCTA_PC or --pc-prefix PC")
    stop(f"could not find {n_pcs} PC columns named {' or '.join(c + '1..' for c in candidates)}. "
         f"Columns in the file: {', '.join(header)}")


def read_covar(path, covariates, n_pcs, pc_prefix):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for col in ["IID", "STATUS"] + list(covariates):
            if col not in header:
                stop(f"{path} has no {col} column")
        pc_cols = pick_pc_columns(header, n_pcs, pc_prefix) if n_pcs > 0 else []
        keep_cols = ["STATUS"] + list(covariates) + pc_cols

        rows = {}
        order = []
        for lineno, line in enumerate(fh, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                stop(f"{path} line {lineno} has {len(fields)} fields, the header has {len(header)}")
            record = dict(zip(header, fields))
            iid = record["IID"]
            order.append(iid)
            rows[iid] = {c: record[c].strip() for c in keep_cols}

    dup = find_duplicates(order)
    if dup:
        stop(f"duplicate IIDs in {path}: {show(dup)}")
    return rows, pc_cols


def read_psam(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split()
        if not header[0].startswith("#"):
            stop(f"{path} does not start with a #FID or #IID header line")
        header[0] = header[0][1:]
        if "IID" not in header:
            stop(f"{path} has no IID column")
        col = header.index("IID")
        ids = [line.split()[col] for line in fh if line.strip()]
    dup = find_duplicates(ids)
    if dup:
        stop(f"duplicate IIDs in {path}: {show(dup)}")
    return ids


def read_msp_samples(path):
    """Sample IDs from the '#chm ...' header line (columns <id>.0 and <id>.1)."""
    with open(path) as fh:
        code_line = fh.readline().rstrip("\n")
        header = fh.readline().rstrip("\n").split("\t")
    samples = []
    for name in header[6:]:
        sample, dot, hap = name.rpartition(".")
        if not dot or hap not in ("0", "1"):
            stop(f"{path}: column '{name}' is not <sample>.0 or <sample>.1")
        if hap == "0":
            samples.append(sample)
    return code_line, samples


def read_vcf_samples(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("##"):
                break
    stop(f"{path}: no #CHROM header line")


def same_set_on_every_chromosome(template, chroms, reader, what):
    """Read the sample set of every chromosome and require them to be identical."""
    first = None
    for chrom in chroms:
        path = template.replace("{chrom}", str(chrom))
        if not os.path.exists(path):
            stop(f"{what} for chr{chrom} not found: {path}")
        samples = set(reader(path))
        if first is None:
            first, first_chrom = samples, chrom
        elif samples != first:
            stop(f"the {what} samples of chr{chrom} differ from chr{first_chrom}.\n"
                 f"  only in chr{chrom}: {show(samples - first)}\n"
                 f"  only in chr{first_chrom}: {show(first - samples)}")
    return first


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--covar", required=True, help="Tab-delimited covariate file with IID, STATUS and the model covariates.")
    ap.add_argument("--model", default="model_covariates.txt", help="Covariates and their types (see the file).")
    ap.add_argument("--psam", required=True, help=".psam of the typed-variant pfile used for the GRM.")
    ap.add_argument("--msp-template", required=True, help="gnomix .msp per chromosome, with {chrom}.")
    ap.add_argument("--vcf-template", default="", help="VCF per chromosome, with {chrom} (optional).")
    ap.add_argument("--chrs", default="1-22")
    ap.add_argument("--pc-prefix", default="", help="GCTA_PC or PC. Default: whichever the file has.")
    ap.add_argument("--outdir", default="samples")
    args = ap.parse_args()

    chroms = parse_chrs(args.chrs)

    covariates, n_pcs = read_model(args.model)
    covar, pc_cols = read_covar(args.covar, covariates, n_pcs, args.pc_prefix)
    numeric_cols = [c for c, kind in covariates.items() if kind == "quantitative"] + pc_cols
    print(f"model ({args.model}):")
    for col, kind in covariates.items():
        print(f"  {col:<12} {kind}")
    print(f"  PCs          {pc_cols[0]}..{pc_cols[-1]}" if pc_cols else "  PCs          none")
    print(f"covariate file: {len(covar)} samples")

    psam = read_psam(args.psam)
    print(f"typed genotypes (.psam): {len(psam)} samples")

    code_lines = set()

    def msp_reader(path):
        code_line, samples = read_msp_samples(path)
        code_lines.add(code_line)
        return samples

    msp = same_set_on_every_chromosome(args.msp_template, chroms, msp_reader, ".msp")
    if len(code_lines) > 1:
        stop("the ancestry code line differs between chromosomes: " + " | ".join(sorted(code_lines)))
    print(f"gnomix .msp: {len(msp)} samples on every chromosome; codes: {code_lines.pop()}")

    sources = {"covar": set(covar), "psam": set(psam), "msp": msp}
    if args.vcf_template:
        sources["vcf"] = same_set_on_every_chromosome(args.vcf_template, chroms, read_vcf_samples, "VCF")
        print(f"VCF: {len(sources['vcf'])} samples on every chromosome")

    # Kept samples, in .psam order.
    kept = [s for s in psam if all(s in ids for ids in sources.values())]
    if not kept:
        stop("no sample is present in every source")

    # Missing values, STATUS coding and numeric covariates, among the kept samples only.
    missing = []
    bad_status = []
    not_numbers = []
    for s in kept:
        row = covar[s]
        for col, value in row.items():
            if value in MISSING:
                missing.append(f"{s} ({col})")
        if row["STATUS"] not in MISSING and row["STATUS"] not in ("1", "2"):
            bad_status.append(f"{s} ({row['STATUS']})")
        for col in numeric_cols:
            if row[col] not in MISSING:
                try:
                    float(row[col])
                except ValueError:
                    not_numbers.append(f"{s} ({col}={row[col]})")

    problems = []
    if missing:
        problems.append(f"{len(missing)} missing value(s): {show(missing, 50)}")
    if bad_status:
        problems.append(f"{len(bad_status)} STATUS value(s) other than 1/2: {show(bad_status, 50)}")
    if not_numbers:
        problems.append(f"{len(not_numbers)} quantitative value(s) that are not numbers: {show(not_numbers, 50)}")
    if problems:
        stop("covariate problems among samples that would otherwise be analysed:\n  "
             + "\n  ".join(problems))

    # Report.
    os.makedirs(args.outdir, exist_ok=True)
    samples_file = os.path.join(args.outdir, "analysis_samples.txt")
    with open(samples_file, "w") as fh:
        fh.write("\n".join(kept) + "\n")

    kept_set = set(kept)
    everyone = list(psam) + sorted(set().union(*sources.values()) - set(psam))
    report_file = os.path.join(args.outdir, "sample_report.tsv")
    with open(report_file, "w") as fh:
        names = list(sources)
        fh.write("IID\t" + "\t".join(f"in_{n}" for n in names) + "\tkept\n")
        for s in everyone:
            flags = ["1" if s in sources[n] else "0" for n in names]
            fh.write(f"{s}\t" + "\t".join(flags) + f"\t{int(s in kept_set)}\n")

    print("\nnot kept, by missing source (a sample can miss several):")
    for name, ids in sources.items():
        others = set().union(*[v for k, v in sources.items() if k != name])
        print(f"  not in {name:<6} {len(others - ids)}")

    n_case = sum(covar[s]["STATUS"] == "2" for s in kept)
    n_ctrl = len(kept) - n_case
    print(f"\nkept: {len(kept)} samples, {n_case} cases (STATUS 2), {n_ctrl} controls (STATUS 1)")

    # Categorical covariates: a value with no cases or no controls cannot be
    # estimated, and 02_fit_null_binary.r will stop on it. Flag it here already.
    for col, kind in covariates.items():
        if kind != "categorical":
            continue
        counts = Counter((covar[s][col], covar[s]["STATUS"]) for s in kept)
        print(f"{col} by STATUS among kept samples:")
        for value in sorted({k[0] for k in counts}):
            flag = ""
            if counts[(value, "2")] == 0 or counts[(value, "1")] == 0:
                flag = "   <- no cases or no controls; 02_fit_null_binary.r will stop"
            print(f"  {col}={value}: {counts[(value, '2')]} cases, {counts[(value, '1')]} controls{flag}")

    print(f"\nwrote {samples_file}\nwrote {report_file}")


if __name__ == "__main__":
    main()
