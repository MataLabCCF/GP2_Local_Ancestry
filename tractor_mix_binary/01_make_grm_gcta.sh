#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=defq
#SBATCH --job-name=grm_gcta
#SBATCH -o grm_gcta.out
#SBATCH -e grm_gcta.err
#
# Dense GCTA GRM for the binary Tractor-Mix run, from the typed variants
# (TOPMed TYPED, already QC'd) of the analysis samples.
#
#   1. Keep the analysis samples and autosomes, then MAF >= 0.01. plink2 applies
#      --keep before --maf, so allele frequencies come from the analysis samples.
#   2. LD pruning: --indep-pairwise 200 50 0.2 (same settings as the Mata lab's PCA script).
#   3. gcta --make-grm on the pruned variants.
#
# Usage:
#   sbatch 01_make_grm_gcta.sh GENO_PFILE SAMPLES OUTDIR
#     GENO_PFILE  typed-variant pfile prefix 
#     SAMPLES     analysis_samples.txt from 00_define_samples.py
#     OUTDIR      output folder
#   Programs can be set with PLINK2=... GCTA=... (defaults: plink2, gcta64)
#
# Output:
#   OUTDIR/grm.grm.bin, grm.grm.N.bin, grm.grm.id   read by 02_fit_null_binary.r
#   OUTDIR/grm_snps.bed/.bim/.fam                   the variants in the GRM

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 GENO_PFILE SAMPLES OUTDIR" >&2
    exit 1
fi
geno=$1
samples=$2
outdir=$3
plink2=${PLINK2:-plink2}
gcta=${GCTA:-gcta64}
threads=${SLURM_CPUS_PER_TASK:-1}

[[ -f $geno.pgen ]] || { echo "ERROR: not found: $geno.pgen" >&2; exit 1; }
[[ -f $samples ]] || { echo "ERROR: not found: $samples" >&2; exit 1; }
mkdir -p "$outdir"

n_samples=$(grep -c . "$samples")
echo "analysis samples: $n_samples"

# 1-2. Filters and LD pruning. The same filters are repeated in the second call
# so the pruned list is applied to exactly the same samples and variants.
filters=(--pfile "$geno" --keep "$samples" --chr 1-22 --maf 0.01 --threads "$threads")

"$plink2" "${filters[@]}" --indep-pairwise 200 50 0.2 --out "$outdir/prune"
"$plink2" "${filters[@]}" --extract "$outdir/prune.prune.in" --make-bed --out "$outdir/grm_snps"

# plink2 skips --keep IDs it cannot find, so check that every sample made it.
n_fam=$(( $(wc -l < "$outdir/grm_snps.fam") ))
if [[ $n_fam -ne $n_samples ]]; then
    echo "ERROR: $n_fam samples in grm_snps.fam, $n_samples in $samples" >&2
    exit 1
fi

# Duplicate variant IDs would be extracted twice and counted twice in the GRM.
dups=$(cut -f2 "$outdir/grm_snps.bim" | sort | uniq -d)
if [[ -n $dups ]]; then
    echo "WARNING: $(echo "$dups" | wc -l) duplicate variant IDs in grm_snps.bim, e.g.:" >&2
    echo "$dups" | head -n 10 >&2
fi
echo "variants in the GRM: $(( $(wc -l < "$outdir/grm_snps.bim") ))"

# 3. Dense GRM.
"$gcta" --bfile "$outdir/grm_snps" --autosome --maf 0.01 --make-grm \
    --thread-num "$threads" --out "$outdir/grm"

n_grm=$(( $(wc -l < "$outdir/grm.grm.id") ))
if [[ $n_grm -ne $n_samples ]]; then
    echo "ERROR: $n_grm samples in grm.grm.id, $n_samples in $samples" >&2
    exit 1
fi
echo "done: $outdir/grm.grm.bin ($n_grm samples)"
