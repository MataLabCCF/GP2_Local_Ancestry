#!/bin/bash
#
# Split the combined Tractor-Mix results of ONE scan (08_combine_tractormix.sh output)
# into analysis-ready sumstats.
#
# Outputs (<prefix> = OUTPREFIX, default GENOME_TSV_GZ without _genome.tsv.gz):
#   <prefix>_joint.tsv.gz        joint test present (all ancestries tested) and every
#                                ancestry usable; wide, one line per variant
#   <prefix>_partial.tsv.gz      same columns, every other variant with >= 1 usable
#                                ancestry (joint NA, or an ancestry not usable)
#   <prefix>_joint_partial.tsv.gz    joint + partial in one file, in genomic order, with a
#                                final "class" column (joint / partial). Every variant
#                                with >= 1 usable ancestry; ignore the column and it is
#                                simply the two concatenated.
#   <prefix>_anc<k>_<LABEL>.tsv.gz   long sumstats per ancestry, including variants
#                                whose joint Chi2 is NA
#   <prefix>_dropped.tsv.gz      the wide row of every variant with NO usable ancestry,
#                                plus a drop_reason column (not_tested / singular /
#                                below_mac / mixed). Kept for inspection: an all-NA row
#                                loses even the ancestries whose own MAC is fine, because
#                                Tractor-Mix inverts one variance matrix per variant, covering all ancestries
#                                at once. When one ancestry is degenerate — typically the ALT allele sitting on every haplotype of
#                                that ancestry, so its dosage column equals its hapcount and conditioning on local ancestry leaves
#                                 nothing — the patched TractorMix.score_cond returns the all-NA row:.
#                                joint + partial + dropped = every input row.
#   <prefix>_split_summary.tsv   counts per chromosome and total
#
# An ancestry is USABLE when include_anc<k> == 1 (Tractor-Mix tested it: AC > its
# --ac-threshold, ALT alleles only), Eff/SE are not NA (not singular), and
# mac_pass_anc<k> == 1, i.e. MAC_anc<k> > MAC_THRESHOLD from the sidecar. That MAC is
# min(AC, AN - AC) among the ancestry's haplotypes, so the filter is symmetric in
# REF/ALT, unlike Tractor-Mix's own ALT-only rule.
#
# Added columns per ancestry, all from the sidecar except AF: AN_anc<k>, MAC_anc<k>
# (the value the filter tests) and MAF_anc<k>, plus AF_alt_anc<k> = AC/AN (the
# effect-allele frequency) and mac_pass_anc<k>.
#
# Rows are matched to the sidecar POSITIONALLY (both list the same variants in the same
# order); CHR, POS, ID, REF and ALT are then compared field by field as strings to prove
# the two streams are in step. Nothing is hashed or sorted, so memory is constant and an
# ID of "." or one containing colons is harmless.
#
# Usage:
#   bash 09_split_sumstats_by_anc.sh GENOME_TSV_GZ TRACTS_DIR [OUTPREFIX]

# TRACTS_DIR is the OUTDIR of 06_filter_split_tracts.sh (chr<N>/ subdirectories holding
# <prefix>.anc_freqs.tsv.gz).
#
# Env: ANC_LABELS (default NAT,AFR,EUR), MAC_THRESHOLD (20), SE_MAX (1000),
#      ALLOW_AC_MISMATCH=1 (AC check becomes a warning),
#      ALLOW_DEGENERATE=1 (a written SE > SE_MAX becomes a warning).

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 GENOME_TSV_GZ TRACTS_DIR [OUTPREFIX]" >&2
    exit 1
fi
genome=$1
tracts=$2
prefix=${3:-}
if [[ -z $prefix ]]; then
    prefix=${genome%.tsv.gz}
    prefix=${prefix%_genome}
fi
labels=${ANC_LABELS:-NAT,AFR,EUR}
mac_threshold=${MAC_THRESHOLD:-20}
se_max=${SE_MAX:-1000}
allow_ac=${ALLOW_AC_MISMATCH:-0}
allow_degenerate=${ALLOW_DEGENERATE:-0}

[[ -f $genome ]] || { echo "ERROR: not found: $genome" >&2; exit 1; }
[[ -d $tracts ]] || { echo "ERROR: not a directory: $tracts" >&2; exit 1; }

shopt -s nullglob
work=$(mktemp -d)
summary="${prefix}_split_summary.tsv"
cleanup() { rm -rf "$work"; rm -f "${prefix}_"*.tsv.gz.tmp "$summary.tmp"; }
trap cleanup EXIT

# chrom -> sidecar, read by awk at BEGIN. Only chromosomes present in the results are
# required, so a missing one is reported by awk, not here.
map="$work/sidecars.tsv"
: > "$map"
for c in {1..22}; do
    files=("$tracts/chr$c"/*.anc_freqs.tsv.gz)
    if [[ ${#files[@]} -eq 1 ]]; then
        printf 'chr%s\t%s\n' "$c" "${files[0]}" >> "$map"
    elif [[ ${#files[@]} -gt 1 ]]; then
        echo "ERROR: ${#files[@]} *.anc_freqs.tsv.gz files in $tracts/chr$c" >&2
        exit 1
    fi
done
[[ -s $map ]] || { echo "ERROR: no chr<N>/*.anc_freqs.tsv.gz under $tracts" >&2; exit 1; }

gzip -dc "$genome" | awk -F'\t' -v OFS='\t' \
    -v map="$map" -v prefix="$prefix" -v labels="$labels" -v mac_threshold="$mac_threshold" \
    -v se_max="$se_max" -v allow_ac="$allow_ac" -v allow_degenerate="$allow_degenerate" \
    -v summary="$summary.tmp" '
function fail(msg) {
    print "ERROR: " msg > "/dev/stderr"
    failed = 1
    exit 1
}
# Columns of a header line, by name.
function index_header(line, into,   n, i, f) {
    n = split(line, f, "\t")
    for (i = 1; i <= n; i++) into[f[i]] = i
    return n
}
function open_chrom(chrom,   hdr, k) {
    if (!(chrom in side)) fail("no sidecar for " chrom " under the tracts directory")
    if (chrom in seen) fail(chrom " appears again after another chromosome; results must be contiguous")
    seen[chrom] = 1
    cur = chrom
    scmd = "gzip -dc \"" side[chrom] "\""
    if ((scmd | getline hdr) <= 0) fail("cannot read " side[chrom])
    delete sidx
    index_header(hdr, sidx)
    for (k = 0; k < n_anc; k++) {
        if (!(("AC_anc" k) in sidx) || !(("AN_anc" k) in sidx) ||
            !(("MAC_anc" k) in sidx) || !(("MAF_anc" k) in sidx))
            fail("sidecar " side[chrom] " has no AC/AN/MAC/MAF columns for anc" k)
    }
}
function close_chrom(   rest, k) {
    if (scmd == "") return
    if ((scmd | getline rest) > 0)
        fail(cur ": the sidecar has more variants than the results")
    close(scmd)
    scmd = ""
    line = cur "\t" c_in "\t" c_joint "\t" c_partial "\t" c_dropped \
               "\t" c_dnot "\t" c_dsing "\t" c_dmac "\t" c_dmix
    for (k = 0; k < n_anc; k++)
        line = line "\t" c_keep[k] "\t" c_nottested[k] "\t" c_singular[k] "\t" c_belowmac[k] \
                    "\t" c_pre[k] "\t" c_post[k]
    print line > summary
    t_in += c_in; t_joint += c_joint; t_partial += c_partial; t_dropped += c_dropped
    t_dnot += c_dnot; t_dsing += c_dsing; t_dmac += c_dmac; t_dmix += c_dmix
    for (k = 0; k < n_anc; k++) {
        t_keep[k] += c_keep[k]; t_nottested[k] += c_nottested[k]
        t_singular[k] += c_singular[k]; t_belowmac[k] += c_belowmac[k]
        t_pre[k] += c_pre[k]; t_post[k] += c_post[k]
        c_keep[k] = c_nottested[k] = c_singular[k] = c_belowmac[k] = c_pre[k] = c_post[k] = 0
    }
    c_in = c_joint = c_partial = c_dropped = 0
    c_dnot = c_dsing = c_dmac = c_dmix = 0
}
BEGIN {
    n_lab = split(labels, lab, ",")
    while ((getline line < map) > 0) {
        split(line, a, "\t")
        side[a[1]] = a[2]
    }
    close(map)
}
NR == 1 {
    n_col = index_header($0, idx)
    for (k = 0; ("Eff_anc" k) in idx; k++) n_anc = k + 1
    if (n_anc < 2) fail("no Eff_anc* columns in the results header")
    if (n_anc != n_lab)
        fail("ANC_LABELS has " n_lab " labels for " n_anc " ancestries in " FILENAME)
    for (k = 0; k < n_anc; k++) {
        if (!(("SE_anc" k) in idx) || !(("Pval_anc" k) in idx) ||
            !(("AC_count" k) in idx) || !(("include_anc" k) in idx))
            fail("results header lacks SE/Pval/AC_count/include for anc" k)
    }
    for (i = 1; i <= 5; i++) {
        want = (i == 1 ? "CHR" : (i == 2 ? "POS" : (i == 3 ? "ID" : (i == 4 ? "REF" : "ALT"))))
        if (idx[want] != i) fail("results column " i " is not " want)
    }
    if (!("Chi2" in idx) || !("P" in idx)) fail("results header lacks Chi2 / P")

    # Output streams and their headers.
    wide = "CHR\tPOS\tID\tREF\tALT\tChi2\tP"
    for (k = 0; k < n_anc; k++)
        wide = wide "\tEff_anc" k "\tSE_anc" k "\tPval_anc" k "\tAC_count" k \
                    "\tAN_anc" k "\tAF_alt_anc" k "\tMAC_anc" k "\tMAF_anc" k \
                    "\tinclude_anc" k "\tmac_pass_anc" k
    joint_cmd = "gzip -9 > \"" prefix "_joint.tsv.gz.tmp\""
    partial_cmd = "gzip -9 > \"" prefix "_partial.tsv.gz.tmp\""
    dropped_cmd = "gzip -9 > \"" prefix "_dropped.tsv.gz.tmp\""
    both_cmd = "gzip -9 > \"" prefix "_joint_partial.tsv.gz.tmp\""
    print wide | joint_cmd
    print wide | partial_cmd
    print wide "\tdrop_reason" | dropped_cmd
    print wide "\tclass" | both_cmd
    for (k = 0; k < n_anc; k++) {
        anc_cmd[k] = "gzip -9 > \"" prefix "_anc" k "_" lab[k + 1] ".tsv.gz.tmp\""
        print "CHR\tPOS\tID\tREF\tALT\tEff_anc" k "\tSE_anc" k "\tPval_anc" k "\tAC_count" k \
              "\tAN_anc" k "\tAF_alt_anc" k "\tMAC_anc" k "\tMAF_anc" k | anc_cmd[k]
    }
    line = "chrom\trows_in\tjoint\tpartial\tdropped_no_usable" \
           "\tdropped_not_tested\tdropped_singular\tdropped_below_mac\tdropped_mixed"
    for (k = 0; k < n_anc; k++)
        line = line "\tkept_" lab[k + 1] "\tnot_tested_" lab[k + 1] "\tsingular_" lab[k + 1] \
                    "\tbelow_mac_" lab[k + 1] "\tse_pre_" lab[k + 1] "\tse_post_" lab[k + 1]
    print line > summary
    next
}
{
    if (NF != n_col) fail("line " NR ": " NF " fields, header has " n_col)
    if ($1 != cur) { close_chrom(); open_chrom($1) }
    c_in++

    # --- Check 1: the two streams are in step (field by field, as strings) ---
    if ((scmd | getline sline) <= 0)
        fail(cur ": more result rows than variants in " side[cur] " (at line " NR ")")
    split(sline, s, "\t")
    for (i = 1; i <= 5; i++) {
        if (($i "") != (s[i] "")) {
            want = (i == 1 ? "CHR" : (i == 2 ? "POS" : (i == 3 ? "ID" : (i == 4 ? "REF" : "ALT"))))
            fail("line " NR " (" $1 " " $2 "): " want " is " $i " in the results but " s[i] \
                 " in " side[cur])
        }
    }

    n_usable = 0; row_singular = 0; row_belowmac = 0
    for (k = 0; k < n_anc; k++) {
        eff = $(idx["Eff_anc" k]); se = $(idx["SE_anc" k])
        ac_tm = $(idx["AC_count" k]); inc = $(idx["include_anc" k])
        ac = s[sidx["AC_anc" k]]; an = s[sidx["AN_anc" k]]
        mac = s[sidx["MAC_anc" k]]; maf = s[sidx["MAF_anc" k]]

        # --- Check 2: Tractor-Mix and the sidecar count the same samples ---
        if (ac_tm + 0 != ac + 0) {
            n_ac_bad++
            if (n_ac_bad <= 5)
                print "AC mismatch: " $1 " " $2 " " $3 " anc" k ": AC_count" k "=" ac_tm \
                      " but sidecar AC_anc" k "=" ac > "/dev/stderr"
            if (!allow_ac)
                fail("AC_count" k " != sidecar AC_anc" k " at line " NR \
                     "; the results and the sidecar are not the same samples" \
                     " (ALLOW_AC_MISMATCH=1 to continue)")
        }

        # --- Check 3a: degenerate SE before any filtering ---
        if (se != "NA" && se + 0 > se_max + 0) {
            c_pre[k]++
            if (t_pre[k] + c_pre[k] <= 3)
                print "note: SE_anc" k " = " se " at " $1 " " $2 " " $3 \
                      " (unpatched TractorMix.score_cond writes these)" > "/dev/stderr"
        }

        af[k] = (an + 0 > 0) ? sprintf("%.6g", ac / an) : "NA"
        pass[k] = (mac + 0 > mac_threshold + 0) ? 1 : 0
        an_v[k] = an; mac_v[k] = mac; maf_v[k] = maf; inc_v[k] = inc

        usable[k] = (inc + 0 == 1 && eff != "NA" && se != "NA" && pass[k])
        if (usable[k]) {
            n_usable++
            c_keep[k]++
            print $1, $2, $3, $4, $5, eff, se, $(idx["Pval_anc" k]), ac_tm, an, af[k], mac, maf | anc_cmd[k]
            # --- Check 3b: a degenerate SE that reached an output file ---
            if (se + 0 > se_max + 0) c_post[k]++
        } else if (inc + 0 != 1) {
            c_nottested[k]++
        } else if (eff == "NA" || se == "NA") {
            c_singular[k]++
            row_singular++
        } else {
            c_belowmac[k]++
            row_belowmac++
        }
    }

    row = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $(idx["Chi2"]) OFS $(idx["P"])
    for (k = 0; k < n_anc; k++)
        row = row OFS $(idx["Eff_anc" k]) OFS $(idx["SE_anc" k]) OFS $(idx["Pval_anc" k]) \
                  OFS $(idx["AC_count" k]) OFS an_v[k] OFS af[k] OFS mac_v[k] OFS maf_v[k] \
                  OFS inc_v[k] OFS pass[k]

    if (n_usable == 0) {
        # Nothing survived. Keep the row with why: "singular" means Tractor-Mix tested an
        # ancestry but returned NA, which takes every other ancestry of that variant with
        # it (one variance matrix per variant), however healthy their own MAC.
        if (row_singular == 0 && row_belowmac == 0)      { reason = "not_tested"; c_dnot++ }
        else if (row_singular > 0 && row_belowmac > 0)   { reason = "mixed";      c_dmix++ }
        else if (row_singular > 0)                       { reason = "singular";   c_dsing++ }
        else                                             { reason = "below_mac";  c_dmac++ }
        print row OFS reason | dropped_cmd
        c_dropped++
        next
    }

    # The composite file gets the same row, tagged, so it stays in genomic order.
    if ($(idx["Chi2"]) != "NA" && $(idx["P"]) != "NA" && n_usable == n_anc) {
        print row | joint_cmd
        print row OFS "joint" | both_cmd
        c_joint++
    } else {
        print row | partial_cmd
        print row OFS "partial" | both_cmd
        c_partial++
    }
}
END {
    if (failed) exit 1
    close_chrom()
    line = "total\t" t_in "\t" t_joint "\t" t_partial "\t" t_dropped \
           "\t" t_dnot "\t" t_dsing "\t" t_dmac "\t" t_dmix
    for (k = 0; k < n_anc; k++)
        line = line "\t" t_keep[k] "\t" t_nottested[k] "\t" t_singular[k] "\t" t_belowmac[k] \
                    "\t" t_pre[k] "\t" t_post[k]
    print line > summary
    close(summary)
    if (close(joint_cmd) != 0 || close(partial_cmd) != 0 || close(dropped_cmd) != 0 ||
        close(both_cmd) != 0)
        fail("gzip failed on a wide output")
    for (k = 0; k < n_anc; k++)
        if (close(anc_cmd[k]) != 0) fail("gzip failed on the anc" k " output")

    pre_total = post_total = 0
    for (k = 0; k < n_anc; k++) { pre_total += t_pre[k]; post_total += t_post[k] }
    if (pre_total > 0)
        print "WARNING: SE > " se_max " in " pre_total " input row(s) before the MAC filter;" \
              " those parts predate the patched TractorMix.score_cond" > "/dev/stderr"
    if (post_total > 0) {
        print "SE > " se_max " in " post_total " row(s) WRITTEN to the ancestry files" > "/dev/stderr"
        if (!allow_degenerate)
            fail("degenerate results reached the sumstats; rerun those parts with the patched" \
                 " TractorMix.score_cond (ALLOW_DEGENERATE=1 to keep them)")
    }
    if (t_in == 0) fail("no data rows in the results file")
}
'

for f in "${prefix}_"*.tsv.gz.tmp; do
    mv "$f" "${f%.tmp}"
done
mv "$summary.tmp" "$summary"
trap - EXIT
rm -rf "$work"

echo "wrote:"
ls -1 "${prefix}_joint.tsv.gz" "${prefix}_partial.tsv.gz" "${prefix}_joint_partial.tsv.gz" \
      "${prefix}_dropped.tsv.gz" "${prefix}_anc"*.tsv.gz "$summary"
echo
awk -F'\t' 'NR == 1 || $1 == "total"' "$summary" | column -t
