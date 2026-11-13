#!/usr/bin/env python3
"""
Reorder gnomix .msp files to the sample order of the phased VCFs, for Tractor's
extract_tracts (c), and check that every analysis sample has genotypes and local
ancestry.

WHY THE ORDER MATTERS
---------------------
extract_tracts (c) skips every '#' line of the .msp, the header included, and
pairs haplotype columns with VCF samples by POSITION. It never reads sample
names. A .msp whose columns are ordered differently from the VCF gives each
person someone else's local ancestry, with no error. After this script, columns
2k-1 and 2k of the .msp belong to sample k of the VCF, verified by name.

CHECKS (any failure stops the script)
------
Per chromosome:
  - VCF samples are unique and are exactly the --samples list (the null-model
    samples: cases, controls, complete covariates, present in the GRM), in both directions.
  - Every VCF sample has both haplotype columns (<id>.0 and <id>.1) in the .msp.
  - The first --check-variants VCF records have GT as the first FORMAT key and
    phased ('|'), non-missing genotypes: extract_tracts (c) reads only GT and
    splits it on '|'.
  - The .msp header is chm, spos, epos, sgpos, egpos, n snps, then <id>.0/<id>.1
    columns; extract_tracts (c) takes haplotype calls from column index 6.
  - Every .msp row has the full column count, the expected chromosome, and only
    ancestry codes declared on the first line.
Across chromosomes:
  - The ancestry code line is identical, so anc0/anc1/anc2 mean the same ancestry
    in every chromosome's extract_tracts (c) output.

Within a sample, haplotype order (.0 then .1) is preserved: swapping it would
break the phase pairing with the VCF.

Also reports mean local ancestry proportions of the analysis samples, weighted
by window length in bp, per chromosome and genome-wide. The Tractor-Mix authors
advise modelling only ancestries with >=10% average ancestry, but in reality we can go as low as 2-5% in our experience, 
but the MAC filter should be a little bit more stringent, specially since TRACTOR-Mix suffers whith high case control imbalance per anc background.

Outputs:
    <outdir>/chr{chrom}.query_results.ordered.msp
    <outdir>/prep_msp_report.tsv

Usage:
    python3 04_prep_msp.py \\
    --vcf-template 'vcf_analysis/chr{chrom}.analysis.vcf.gz' \
    --msp-template 'gnomix/LPD_chrom{chrom}/query_results.msp' \
    --samples      samples/analysis_samples.txt \
    --outdir       msp_ordered
"""

import argparse
import gzip
import os
import re
import sys
from collections import Counter

# ---------------------------------------------------------------------------

# gnomix and RFMix2 share this layout; "n snps" contains a space, so the .msp
# must be split on tabs only.
FIXED_COLS = ["chm", "spos", "epos", "sgpos", "egpos", "n snps"]
HAP_RE = re.compile(r"^(.+)\.([01])$")


def parse_chrs(spec: str) -> list[int]:
    out: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        elif part:
            out.append(int(part))
    return sorted(set(out))


def norm_chrom(c: str) -> str:
    return c.removeprefix("chr")


def examples(items, n: int = 5) -> str:
    items = sorted(items)
    return ", ".join(items[:n]) + (f" ... (+{len(items) - n})" if len(items) > n else "")


def duplicates(ids: list[str]) -> list[str]:
    return sorted(k for k, v in Counter(ids).items() if v > 1)


def read_samples(path: str) -> list[str]:
    with open(path) as fh:
        ids = [line.strip() for line in fh if line.strip()]
    if not ids:
        sys.exit(f"ERROR: {path} has no sample IDs")
    dup = duplicates(ids)
    if dup:
        sys.exit(f"ERROR: duplicate IDs in {path}: {examples(dup)}")
    return ids


# ---------------------------------------------------------------------------
# VCF: sample order, and GT layout on the first records
# ---------------------------------------------------------------------------
def read_vcf(path: str, n_check: int) -> tuple[list[str], str]:
    """Return (samples in file order, CHROM of the first record)."""
    opener = gzip.open if path.endswith(".gz") else open
    samples: list[str] | None = None
    chrom = None
    checked = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if samples is None:
                sys.exit(f"ERROR: {path}: variant record before the #CHROM header")

            f = line.rstrip("\n").split("\t")
            where = f"{path} {f[0]}:{f[1]}"
            if chrom is None:
                chrom = f[0]
            if len(f) - 9 != len(samples):
                sys.exit(f"ERROR: {where} has {len(f) - 9} genotypes for {len(samples)} samples")
            if f[8].split(":")[0] != "GT":
                sys.exit(
                    f"ERROR: {where}: FORMAT is '{f[8]}'. extract_tracts (c) reads the first "
                    "FORMAT key as GT; a dosage-only VCF cannot be used."
                )
            bad = []
            for s, g in zip(samples, f[9:]):
                gt = g.split(":", 1)[0]
                if "|" not in gt or "." in gt:
                    bad.append(f"{s}={gt}")
            if bad:
                sys.exit(
                    f"ERROR: {where}: {len(bad)} unphased or missing GT, e.g. {examples(bad, 3)}"
                )
            checked += 1
            if checked >= n_check:
                break

    if samples is None:
        sys.exit(f"ERROR: {path}: no #CHROM header")
    if checked == 0:
        sys.exit(f"ERROR: {path}: no variant records")
    return samples, chrom


# ---------------------------------------------------------------------------
# .msp: validate, reorder to the VCF, tally ancestry
# ---------------------------------------------------------------------------
def reorder_msp(
    msp: str, out: str, vcf_samples: list[str], chrom: str, expected_codes: str | None
) -> tuple[str, dict[str, str], Counter, int]:
    """
    Write `msp` to `out` with haplotype columns in `vcf_samples` order.
    Return (code line, {code: ancestry}, bp-weighted haplotype counts by code,
    number of windows).
    """
    with open(msp) as fh:
        code_line = fh.readline().rstrip("\n")
        codes = dict(re.findall(r"([A-Za-z0-9_]+)=(\d+)", code_line))
        if not code_line.startswith("#") or not codes:
            sys.exit(f"ERROR: {msp}: first line is not an ancestry code line: {code_line[:120]}")
        if expected_codes is not None and code_line != expected_codes:
            sys.exit(
                f"ERROR: {msp}: ancestry codes differ from the previous chromosome.\n"
                f"       here:     {code_line}\n       previous: {expected_codes}\n"
                "       anc0/anc1/... would mean different ancestries across chromosomes."
            )
        code_to_anc = {v: k for k, v in codes.items()}

        header_line = fh.readline().rstrip("\n")
        if not header_line.startswith("#"):
            sys.exit(f"ERROR: {msp}: second line is not a '#chm ...' header")
        header = header_line[1:].split("\t")
        if header[:6] != FIXED_COLS:
            sys.exit(
                f"ERROR: {msp}: leading columns are {header[:6]}, expected {FIXED_COLS}. "
                "extract_tracts (c) takes haplotype calls from column index 6."
            )

        col_of: dict[str, int] = {}
        for j, name in enumerate(header[6:], start=6):
            if not HAP_RE.match(name):
                sys.exit(f"ERROR: {msp}: column {j + 1} '{name}' is not <sample>.0 or <sample>.1")
            if name in col_of:
                sys.exit(f"ERROR: {msp}: duplicate haplotype column '{name}'")
            col_of[name] = j

        no_la = [s for s in vcf_samples if f"{s}.0" not in col_of or f"{s}.1" not in col_of]
        if no_la:
            sys.exit(
                f"ERROR: {len(no_la)} VCF sample(s) lack local ancestry in {msp}: {examples(no_la)}"
            )

        order = [col_of[f"{s}.{h}"] for s in vcf_samples for h in (0, 1)]
        n_cols = len(header)
        weighted: Counter = Counter()
        n_rows = 0

        tmp = out + ".tmp"
        with open(tmp, "w") as oh:
            oh.write(code_line + "\n")
            oh.write("#" + "\t".join(header[:6] + [header[j] for j in order]) + "\n")
            for lineno, line in enumerate(fh, start=3):
                f = line.rstrip("\n").split("\t")
                if len(f) != n_cols:
                    sys.exit(f"ERROR: {msp}:{lineno} has {len(f)} columns, header has {n_cols}")
                if norm_chrom(f[0]) != norm_chrom(chrom):
                    sys.exit(f"ERROR: {msp}:{lineno} is chromosome '{f[0]}', expected {chrom}")
                calls = [f[j] for j in order]
                tally = Counter(calls)
                unknown = set(tally) - set(code_to_anc)
                if unknown:
                    sys.exit(f"ERROR: {msp}:{lineno} has undeclared ancestry code(s) {sorted(unknown)}")
                # Ancestry proportions are weighted by window length in bp.
                #
                # gnomix windows hold a fixed number of SNPs, not a fixed number of bp,
                # so their length varies: ~100 kb where SNPs are dense, ~2 Mb where they
                # are sparse. Each window has one call per haplotype, and that call
                # covers the whole window. Counting windows would let SNP density. Example, one haplotype, two windows:
                #     window A, 100 kb, AFR;  window B, 1,900 kb, NAT
                #     counting windows:  50% AFR, 50% NAT
                #     weighting by bp:    5% AFR, 95% NAT  (100 / 2,000 and 1,900 / 2,000)
                # The second is what "ancestry proportion" means: the fraction of the DNA
                # that came from each ancestry.
                #
                # So n haplotypes with a given code in this window add n * width to that
                # code: the number of haplotype bp labelled with that ancestry here.
                # Gaps between windows (centromere, stretches without typed SNPs) are
                # not covered by any call and add nothing. sgpos/egpos (cM) could be
                # used instead of bp; both measure "fraction of the genome" and
                # generally agree closely.
                width = float(f[2]) - float(f[1])
                for code, n in tally.items():
                    weighted[code] += n * width
                oh.write("\t".join(f[:6] + calls) + "\n")
                n_rows += 1

    if n_rows == 0:
        os.remove(tmp)
        sys.exit(f"ERROR: {msp}: no windows")
    os.replace(tmp, out)
    return code_line, code_to_anc, weighted, n_rows


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-template", required=True,
                    help="Phased, VCF per chromosome, with {chrom}.")
    ap.add_argument("--msp-template", required=True,
                    help="gnomix .msp per chromosome, with {chrom}.")
    ap.add_argument("--samples", required=True,
                    help="Analysis sample IDs, one per line (samples/analysis_samples.txt from step 0).")
    ap.add_argument("--chrs", default="1-22", help="e.g. '1-22' or '22' or '1,2,19'")
    ap.add_argument("--outdir", default="msp_ordered")
    ap.add_argument("--check-variants", type=int, default=200,
                    help="VCF records checked for GT layout and phasing.")
    args = ap.parse_args()

    samples = read_samples(args.samples)
    sample_set = set(samples)
    print(f"analysis samples: {len(samples)} ({args.samples})")

    chroms = parse_chrs(args.chrs)
    for chrom in chroms:
        for path in (args.vcf_template.format(chrom=chrom), args.msp_template.format(chrom=chrom)):
            if not os.path.exists(path):
                sys.exit(f"ERROR: chr{chrom}: {path} not found")
    os.makedirs(args.outdir, exist_ok=True)

    ref_codes = None
    code_to_anc: dict[str, str] = {}
    genome: Counter = Counter()
    report = []
    vcf_orders: set[tuple[str, ...]] = set()

    for chrom in chroms:
        vcf = args.vcf_template.format(chrom=chrom)
        msp = args.msp_template.format(chrom=chrom)
        out = os.path.join(args.outdir, f"chr{chrom}.query_results.ordered.msp")

        vcf_samples, vcf_chrom = read_vcf(vcf, args.check_variants)
        dup = duplicates(vcf_samples)
        if dup:
            sys.exit(f"ERROR: duplicate sample IDs in {vcf}: {examples(dup)}")
        only_vcf = set(vcf_samples) - sample_set
        only_null = sample_set - set(vcf_samples)
        if only_vcf or only_null:
            sys.exit(
                f"ERROR: chr{chrom}: {vcf} does not hold exactly the null-model samples.\n"
                f"       in VCF, not in null model ({len(only_vcf)}): {examples(only_vcf)}\n"
                f"       in null model, not in VCF ({len(only_null)}): {examples(only_null)}"
            )
        if norm_chrom(vcf_chrom) != str(chrom):
            sys.exit(f"ERROR: {vcf}: first record is on '{vcf_chrom}', expected chromosome {chrom}")
        vcf_orders.add(tuple(vcf_samples))

        ref_codes, code_to_anc, weighted, n_win = reorder_msp(msp, out, vcf_samples,
                                                              vcf_chrom, ref_codes)
        # Per chromosome: (bp labelled ancestry a) / (all labelled bp), over every
        # haplotype of every analysis sample. All samples share the same windows,
        # so pooling haplotypes gives the same number as averaging each person's
        # own proportion: it is the average person's ancestry on this chromosome.
        genome.update(weighted)
        total = sum(weighted.values())
        props = {code_to_anc[c]: weighted[c] / total for c in code_to_anc}
        report.append((chrom, len(vcf_samples), n_win, props))
        print(f"chr{chrom}: {len(vcf_samples)} samples, {n_win} windows -> {out}  "
              + "  ".join(f"{a}={p:.3f}" for a, p in props.items()))

    if len(vcf_orders) > 1:
        print("note: VCF sample order differs between chromosomes. Each .msp was matched to "
              "its own VCF, so this is safe; Tractor-Mix matches dosage files to the null "
              "model by ID.")

    ancs = [code_to_anc[c] for c in sorted(code_to_anc, key=int)]
    report_path = os.path.join(args.outdir, "prep_msp_report.tsv")
    with open(report_path, "w") as fh:
        fh.write("chrom\tn_samples\tn_windows\t" + "\t".join(f"prop_{a}" for a in ancs) + "\n")
        for chrom, n, n_win, props in report:
            fh.write(f"{chrom}\t{n}\t{n_win}\t" + "\t".join(f"{props[a]:.5f}" for a in ancs) + "\n")

    # Genome-wide: the bp totals of all chromosomes are added before dividing, so
    # each chromosome counts in proportion to its length (chr1 more than chr22).
    total = sum(genome.values())
    print(f"\nancestry codes: {ref_codes}")
    print(f"mean local ancestry, analysis samples, chromosomes {args.chrs} (bp-weighted):")
    for code in sorted(code_to_anc, key=int):
        p = genome[code] / total
        flag = "   <10%: Tractor-Mix authors put 10% as the lower end. In our experience, we can go as low as around 4-5%, considering bumping up the MAC filter" if p < 0.10 else ""
        print(f"  anc{code} {code_to_anc[code]:<6} {p:.3f}{flag}")
    print(f"\nwrote {report_path}")
    print("Run extract_tracts (c) on each ordered .msp with its own VCF, "
          f"--num-ancs {len(code_to_anc)}.")


if __name__ == "__main__":
    main()
