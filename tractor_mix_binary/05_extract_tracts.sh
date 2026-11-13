#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --partition=defq
#SBATCH --job-name=extract_tracts
#SBATCH -o logs/extract_tracts_chr%a.out
#SBATCH -e logs/extract_tracts_chr%a.err
#
# Run extract_tracts on one chromosome.
#
# The number of ancestries is read from the ancestry code line of the ordered
# .msp (e.g. "#Subpopulation order/codes: NAT=0 AFR=1 EUR=2"), so the same
# script works for 3 or more ancestries. anc0, anc1, ... in the output follow
# the codes on that line.
#
# Usage:
#   mkdir -p logs
#   sbatch --array=1-22 05_extract_tracts.sh VCF_DIR MSP_DIR OUTDIR
#     VCF_DIR  output of 03_subset_vcf.sh (chr<N>.analysis.vcf.gz)
#     MSP_DIR  output of 04_prep_msp.py (chr<N>.query_results.ordered.msp)
#     OUTDIR   results go to OUTDIR/chr<N>/, the layout 06_filter_split_tracts.sh reads
#   The extract_tracts binary can be set with EXTRACT_TRACTS=... (default ./extract_tracts)
#   One chromosome without SLURM: CHROM=22 bash 05_extract_tracts.sh ...

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 VCF_DIR MSP_DIR OUTDIR" >&2
    exit 1
fi
vcf_dir=$1
msp_dir=$2
outdir=$3
chrom=${CHROM:-${SLURM_ARRAY_TASK_ID:-}}
extract_tracts=${EXTRACT_TRACTS:-./extract_tracts}
[[ -n $chrom ]] || { echo "ERROR: set CHROM or submit as an array job" >&2; exit 1; }

vcf="$vcf_dir/chr$chrom.analysis.vcf.gz"
msp="$msp_dir/chr$chrom.query_results.ordered.msp"
for f in "$vcf" "$msp"; do
    [[ -f $f ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

code_line=$(head -n 1 "$msp")
n_anc=$(( $(echo "$code_line" | grep -o '[A-Za-z0-9_]*=[0-9]*' | wc -l) ))
if [[ $n_anc -lt 2 ]]; then
    echo "ERROR: could not read the ancestries from the first line of $msp: $code_line" >&2
    exit 1
fi
echo "chr$chrom: $n_anc ancestries ($code_line)"

mkdir -p "$outdir/chr$chrom"
"$extract_tracts" --vcf "$vcf" --msp "$msp" --num-ancs "$n_anc" \
    --compress-output --output-dir "$outdir/chr$chrom"
