#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --partition=defq
#SBATCH --job-name=filter_tracts
#SBATCH -o logs/filter_tracts_chr%a.out
#SBATCH -e logs/filter_tracts_chr%a.err
#
# Filter and split ONE chromosome of extract_tracts output for Tractor-Mix.
#
# Why: Tractor-Mix reads its inputs with data.table::fread(skip =, nrows =
# chunk_size), and fread decompresses a whole .gz file (or scans a plain file
# from the top) on every call. The cost grows with chunks x file size, which
# does not finish for whole-chromosome files of all imputed variants.
#
# What:
#   1. Keep a variant if MAC > MAC_THRESHOLD in at least one ancestry, taken
#      from <prefix>.anc_freqs.tsv.gz. Those counts cover every MSP sample, so
#      they are >= the counts on the null-model samples: nothing Tractor-Mix
#      would test with --ac-threshold MAC_THRESHOLD is dropped. MAC rather than
#      Tractor-Mix's ALT count also drops variants whose ALT allele is ~fixed.
#   2. Write the kept rows of every dosage and hapcount file in pieces of
#      PIECE_SIZE variants, all files in parallel:
#        <prefix>.part0001.anc0.dosage.txt.gz, <prefix>.part0001.anc0.hapcount.txt.gz, ...
#      A lone last variant joins the previous piece (Tractor-Mix cannot run a
#      one-variant file). Each kept row's CHROM:POS:ID:REF:ALT is checked
#      against anc_freqs, and every file's row count is checked at the end.
#   3. Work in OUTDIR/chr<N>.tmp and rename it to OUTDIR/chr<N> only when every
#      file succeeded, so an existing OUTDIR/chr<N> is complete (and skipped).
#
# Output is gzip -9: zlib's default level 6 gave some problems.
#
# Usage (the chromosome comes from $CHROM, else the array index):
#   mkdir -p logs
#   sbatch --array=1-22 06_filter_split_tracts.sh INDIR OUTDIR [MAC_THRESHOLD=20] [PIECE_SIZE=16384]
#   CHROM=22 bash 06_filter_split_tracts.sh INDIR OUTDIR 20 16384
# INDIR holds chr<N>/<prefix>.anc_freqs.tsv.gz, <prefix>.anc{k}.dosage.txt.gz and
# <prefix>.anc{k}.hapcount.txt.gz (the extract_tracts --output-dir layout).
# Keep PIECE_SIZE a multiple of the Tractor-Mix --chunk-size (2048).

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: [CHROM=N] $0 INDIR OUTDIR [MAC_THRESHOLD=20] [PIECE_SIZE=16384]" >&2
    exit 1
fi
indir=$1
outdir=$2
mac_threshold=${3:-20}
piece_size=${4:-16384}
chrom=${CHROM:-${SLURM_ARRAY_TASK_ID:-}}
if [[ -z $chrom ]]; then
    echo "ERROR: set CHROM or submit as an array job" >&2
    exit 1
fi
if ! [[ $mac_threshold =~ ^[0-9]+$ && $piece_size =~ ^[0-9]+$ && $piece_size -ge 2 ]]; then
    echo "ERROR: MAC_THRESHOLD must be an integer >= 0 and PIECE_SIZE an integer >= 2" >&2
    exit 1
fi

src="$indir/chr$chrom"
dest="$outdir/chr$chrom"
work="$outdir/chr$chrom.tmp"
log() { echo "[$(date '+%F %T')] chr$chrom: $*"; }

if [[ -d $dest ]]; then
    log "$dest exists, skipping"
    exit 0
fi

shopt -s nullglob
freqs_files=("$src"/*.anc_freqs.tsv.gz)
if [[ ${#freqs_files[@]} -ne 1 ]]; then
    echo "ERROR: expected one *.anc_freqs.tsv.gz in $src, found ${#freqs_files[@]}" >&2
    exit 1
fi
freqs=${freqs_files[0]}
base=$(basename "$freqs" .anc_freqs.tsv.gz)

# head closes the pipe early; the braces keep gzip's SIGPIPE from tripping pipefail.
n_anc=$({ gzip -dc "$freqs" || true; } | head -n 1 | tr '\t' '\n' | grep -c '^MAC_anc') || true
if [[ $n_anc -lt 1 ]]; then
    echo "ERROR: no MAC_anc columns in the header of $freqs" >&2
    exit 1
fi
inputs=()
for ((k = 0; k < n_anc; k++)); do
    for what in dosage hapcount; do
        f="$src/$base.anc$k.$what.txt.gz"
        [[ -f $f ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
        inputs+=("$f")
    done
done

log "$n_anc ancestries | keep MAC > $mac_threshold in any ancestry | pieces of $piece_size"
rm -rf "$work"
mkdir -p "$work"
mask="$work/mask.txt"

# ---------------------------------------------------------------------------
# 1. Mask (one line per variant: CHROM:POS:ID:REF:ALT if kept, else 0) and the
#    filtered anc_freqs.
# ---------------------------------------------------------------------------
gzip -dc "$freqs" | awk -F'\t' -v t="$mac_threshold" -v mask="$mask" -v counts="$work/counts" '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i ~ /^MAC_anc[0-9]+$/) mac[++n_mac] = i
    print
    next
}
{
    keep = 0
    for (j = 1; j <= n_mac; j++) if ($(mac[j]) + 0 > t + 0) { keep = 1; break }
    if (keep) { print $1 ":" $2 ":" $3 ":" $4 ":" $5 > mask; print; kept++ }
    else print "0" > mask
}
END { print NR - 1, kept + 0 > counts }
' | gzip -9 > "$work/$base.anc_freqs.tsv.gz"

read -r n_total n_kept < "$work/counts"
log "variants: $n_total | kept: $n_kept"
if [[ $n_kept -eq 0 ]]; then
    echo "ERROR: no variant has MAC > $mac_threshold in any ancestry" >&2
    exit 1
fi
if [[ $n_kept -eq 1 ]]; then
    log "WARNING: a single kept variant; Tractor-Mix cannot run a one-variant file"
fi

# ---------------------------------------------------------------------------
# 2. Filter and split every dosage / hapcount file, in parallel.
# ---------------------------------------------------------------------------
split_file() {  # $1 input .gz, $2 label (e.g. anc0.dosage)
    gzip -dc "$1" | awk -v mask="$mask" -v kept="$n_kept" -v piece="$piece_size" \
                    -v prefix="$work/$base" -v label="$2" '
    NR == 1 { header = $0; next }
    {
        if ((getline m < mask) <= 0) { err = "more rows than anc_freqs"; exit }
        if (m == "0") next
        # Only the leading 5 fields are compared; $0 is never split into fields.
        match($0, /^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*/)
        key = substr($0, 1, RLENGTH)
        gsub(/\t/, ":", key)
        if (key != m) { err = "variant row " NR - 1 " is " key ", anc_freqs has " m; exit }
        # New piece every `piece` variants, except that a lone last variant stays
        # in the previous piece.
        if (n % piece == 0 && (n == 0 || n != kept - 1)) {
            if (cmd != "" && close(cmd) != 0) { err = "gzip failed on " cmd; exit }
            cmd = sprintf("gzip -9 > \"%s.part%04d.%s.txt.gz\"", prefix, ++part, label)
            print header | cmd
        }
        print | cmd
        n++
    }
    END {
        if (err == "" && cmd != "" && close(cmd) != 0) err = "gzip failed on " cmd
        if (err == "" && (getline m < mask) > 0) err = "fewer rows than anc_freqs"
        if (err == "" && n != kept) err = n " variants written, expected " kept
        if (err != "") { print "ERROR: " label ": " err > "/dev/stderr"; exit 1 }
        print label ": " n " variants in " part " pieces"
    }'
}

pids=()
labels=()
for f in "${inputs[@]}"; do
    label=$(basename "$f" .txt.gz)
    label=${label#"$base."}
    split_file "$f" "$label" > "$work/$label.log" 2>&1 &
    pids+=($!)
    labels+=("$label")
done
failed=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        failed=1
    fi
    log "$(cat "$work/${labels[$i]}.log")"
done
if [[ $failed -ne 0 ]]; then
    echo "ERROR: at least one file failed; partial output left in $work" >&2
    exit 1
fi

# Every file must have the same pieces.
n_pieces=
for label in "${labels[@]}"; do
    files=("$work/$base".part*."$label".txt.gz)
    if [[ -z $n_pieces ]]; then
        n_pieces=${#files[@]}
    elif [[ ${#files[@]} -ne $n_pieces ]]; then
        echo "ERROR: $label has ${#files[@]} pieces, ${labels[0]} has $n_pieces" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# 3. Finalize.
# ---------------------------------------------------------------------------
rm -f "$mask" "$work/counts" "$work"/*.log
{
    echo -e "source\t$src"
    echo -e "mac_threshold\t$mac_threshold"
    echo -e "piece_size\t$piece_size"
    echo -e "variants_total\t$n_total"
    echo -e "variants_kept\t$n_kept"
    echo -e "pieces\t$n_pieces"
} > "$work/$base.filter_summary.tsv"
mv "$work" "$dest"
log "done: $n_kept variants in $n_pieces pieces -> $dest"
