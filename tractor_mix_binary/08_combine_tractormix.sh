#!/bin/bash
#
# Combine the per-piece Tractor-Mix results of ONE scan
# (<SCAN>_chr<N>_part<P>.tsv from 07_submit_tractormix.py) into one genome-wide
# file, in chromosome and piece order.
#
# Every check is fatal; nothing is written unless all pass:
#   - each chromosome has a result for every piece listed by
#     06_filter_split_tracts.sh (<prefix>.filter_summary.tsv: pieces) and none
#     beyond them (a stale result from an earlier split);
#   - every result file has the same header, and every row as many fields;
#   - the result rows are EXACTLY the kept variants: CHR POS ID REF ALT are
#     compared as strings, row by row, with the filtered
#     <prefix>.anc_freqs.tsv.gz, which lists the kept variants in order. This
#     catches a missing, repeated, reordered or altered row (e.g. an allele T
#     read by R as TRUE), which a row count alone would not.
# Rows are copied verbatim: $0 is never modified, so awk does not reformat any
# number. Leftover <result>.tmp files (runs that failed their own check) are
# reported but do not count as results.
#
# Usage:
#   bash 08_combine_tractormix.sh RESULTS_DIR SCAN TRACTS_DIR [OUT]

# TRACTS_DIR is the OUTDIR of 06_filter_split_tracts.sh (chr<N>/ subdirectories).
# CHRS: chromosomes to combine (default 1-22); a chromosome without results is
# an error, so set CHRS to combine a subset.
# OUT defaults to RESULTS_DIR/<SCAN>_genome.tsv.gz; per-chromosome row counts
# go to OUT with .counts.tsv in place of .tsv.gz.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: [CHRS=\"1 2 ...\"] $0 RESULTS_DIR SCAN TRACTS_DIR [OUT]" >&2
    exit 1
fi
results=$1
scan=$2
tracts=$3
out=${4:-$results/${scan}_genome.tsv.gz}
counts="${out%.tsv.gz}.counts.tsv"
chrs=${CHRS:-$(seq 1 22)}
shopt -s nullglob

# ---------------------------------------------------------------------------
# Completeness, per chromosome: build the awk operands as
#   chrom=<N> freqs=<kept variants> <part0001> <part0002> ...
# (awk applies each var=value before reading the files after it).
# ---------------------------------------------------------------------------
operands=()
problems=0
for chrom in $chrs; do
    summary=("$tracts/chr$chrom"/*.filter_summary.tsv)
    if [[ ${#summary[@]} -ne 1 ]]; then
        echo "chr$chrom: expected one *.filter_summary.tsv in $tracts/chr$chrom, found ${#summary[@]}" >&2
        problems=1
        continue
    fi
    base=$(basename "${summary[0]}" .filter_summary.tsv)
    freqs="$tracts/chr$chrom/$base.anc_freqs.tsv.gz"
    if [[ ! -f $freqs ]]; then
        echo "chr$chrom: not found: $freqs" >&2
        problems=1
        continue
    fi
    n_pieces=$(awk -F'\t' '$1 == "pieces" { print $2 }' "${summary[0]}")

    files=()
    n_missing=0
    for ((p = 1; p <= n_pieces; p++)); do
        f=$(printf '%s/%s_chr%s_part%04d.tsv' "$results" "$scan" "$chrom" "$p")
        if [[ -f $f ]]; then
            files+=("$f")
        else
            [[ $n_missing -lt 5 ]] && echo "chr$chrom: missing $(basename "$f")" >&2
            n_missing=$((n_missing + 1))
        fi
    done
    if [[ $n_missing -gt 0 ]]; then
        echo "chr$chrom: $n_missing of $n_pieces pieces have no result" >&2
        problems=1
    fi
    for f in "$results/${scan}_chr${chrom}_part"*.tsv; do
        p=${f##*_part}
        p=${p%.tsv}
        if ! [[ $p =~ ^[0-9]+$ ]] || ((10#$p > n_pieces)); then
            echo "chr$chrom: $(basename "$f") is not one of the $n_pieces pieces in $tracts/chr$chrom (stale?)" >&2
            problems=1
        fi
    done
    tmps=("$results/${scan}_chr${chrom}_part"*.tsv.tmp)
    if [[ ${#tmps[@]} -gt 0 ]]; then
        echo "chr$chrom: note: ${#tmps[@]} leftover .tmp file(s), e.g. $(basename "${tmps[0]}")" >&2
    fi

    if [[ $n_missing -eq 0 ]]; then
        operands+=("chrom=$chrom" "freqs=$freqs" "${files[@]}")
    fi
done
if [[ $problems -ne 0 ]]; then
    echo "ERROR: results incomplete; nothing written" >&2
    exit 1
fi
if [[ ${#operands[@]} -eq 0 ]]; then
    echo "ERROR: no chromosomes to combine" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Concatenate, checking every row against the kept variants.
# ---------------------------------------------------------------------------
trap 'rm -f "$out.tmp" "$counts.tmp"' EXIT

mkdir -p "$(dirname "$out")"

awk -F'\t' -v counts="$counts.tmp" '
function fail(msg) {
    print "ERROR: " msg > "/dev/stderr"
    failed = 1
    exit 1
}
function kept_key(line,   v) {
    split(line, v, "\t")
    return v[1] ":" v[2] ":" v[3] ":" v[4] ":" v[5]
}
function end_chrom(   rest) {
    if (fcmd == "") return
    if ((fcmd | getline rest) > 0)
        fail("chr" cur ": results end after " n_chr " variants; next kept variant is " kept_key(rest))
    close(fcmd)
    print "chr" cur "\t" n_chr > counts
    fcmd = ""
}
FNR == 1 {
    if (header == "") {
        header = $0
        nf = NF
        if ($1 != "CHR" || $2 != "POS" || $3 != "ID" || $4 != "REF" || $5 != "ALT")
            fail("header of " FILENAME " does not start with CHR POS ID REF ALT")
        print
    } else if ($0 != header) {
        fail("header of " FILENAME " differs from the first file")
    }
    if (chrom != cur) {
        end_chrom()
        cur = chrom
        n_chr = 0
        fcmd = "gzip -dc \"" freqs "\""
        if ((fcmd | getline line) <= 0) fail("cannot read " freqs)
    }
    next
}
{
    if (NF != nf) fail(FILENAME " line " FNR ": " NF " fields, header has " nf)
    if ((fcmd | getline line) <= 0)
        fail(FILENAME " line " FNR ": more results than kept variants on chr" cur)
    got = $1 ":" $2 ":" $3 ":" $4 ":" $5
    want = kept_key(line)
    if (got != want) fail(FILENAME " line " FNR ": " got ", but the kept variant is " want)
    print
    n_chr++
    n_all++
}
END {
    if (failed) exit 1
    end_chrom()
    print "total\t" n_all > counts
}
' "${operands[@]}" | gzip > "$out.tmp"

mv "$counts.tmp" "$counts"
mv "$out.tmp" "$out"
trap - EXIT
echo "wrote $out"
cat "$counts"
