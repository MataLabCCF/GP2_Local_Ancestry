#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=subset_vcf_chr%a
#SBATCH --output=logs/subset_vcf_chr%a.out      
#SBATCH --error=logs/subset_vcf_chr%a.err
#SBATCH --array=1-22                      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=defq

#
# Make the phased, imputed VCF of one chromosome hold exactly the analysis samples.
#
# Why: 04_prep_msp.py requires the VCF to hold exactly the null-model samples,
# and 09_split_sumstats_by_anc.sh checks Tractor-Mix's allele counts against the
# extract_tracts sidecar, which counts every sample in the VCF. With the VCF
# matching the analysis samples, both hold.
#
# Only the set of samples matters here, not their order: 04_prep_msp.py puts the
# msp in the VCF's order, and Tractor-Mix matches samples to the null by ID.
# So if the VCF header already lists exactly the analysis samples, the VCF is
# linked (symlink) instead of rewritten. Otherwise it is subset with bcftools,
# which reads and rewrites the whole file and can take hours on large chromosomes.
#
# Since the subset rewrites the file anyway, it also drops variants with
# MAC < MIN_MAC (default 21, i.e. MAC > 20) in the analysis samples, which removes
# the monomorphic and very rare ones. This loses nothing: 06_filter_split_tracts.sh
# later keeps a variant only if MAC > 20 in at least one ancestry, and that
# implies MAC > 20 overall. The MAC is counted from GT after subsetting (the two
# bcftools calls are piped for that reason), the same GT that extract_tracts reads.
# A linked VCF is not filtered; 06_filter_split_tracts.sh still applies the real filter.
#
# bcftools stops if a sample in the list is not in the VCF (no --force-samples).
# and the list comes form the covariate file which would be used to fit the null model
#
# Usage:
#   mkdir -p logs
#   sbatch --array=1-22 03_subset_vcf.sh VCF_TEMPLATE SAMPLES OUTDIR
#     VCF_TEMPLATE  per-chromosome VCF with {chrom}, e.g. 'imputed/chr{chrom}.dose.vcf.gz'
#     SAMPLES       analysis_samples.txt from 00_define_samples.py
#   Output: OUTDIR/chr<N>.analysis.vcf.gz (+ .tbi), either a new file or a
#           symlink to the original VCF
#   One chromosome without SLURM: CHROM=22 bash 03_subset_vcf.sh ...

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 VCF_TEMPLATE SAMPLES OUTDIR" >&2
    exit 1
fi
template=$1
samples=$2
outdir=$3
chrom=${CHROM:-${SLURM_ARRAY_TASK_ID:-}}
threads=${SLURM_CPUS_PER_TASK:-1}
min_mac=${MIN_MAC:-21}
[[ -n $chrom ]] || { echo "ERROR: set CHROM or submit as an array job" >&2; exit 1; }

vcf=${template//\{chrom\}/$chrom}
out="$outdir/chr$chrom.analysis.vcf.gz"
[[ -f $vcf ]] || { echo "ERROR: not found: $vcf" >&2; exit 1; }
mkdir -p "$outdir"

# Remove the output of an earlier run first. If it was a symlink, its index is
# a symlink to the original's index, and bcftools index would overwrite that.
rm -f "$out" "$out.tbi" "$out.csi"

if cmp -s <(bcftools query -l "$vcf" | sort) <(grep . "$samples" | sort); then
    echo "chr$chrom: VCF already has exactly the analysis samples, linking instead of subsetting"
    vcf_abs="$(cd "$(dirname "$vcf")" && pwd)/$(basename "$vcf")"
    ln -s "$vcf_abs" "$out"
    if [[ -f $vcf_abs.tbi ]]; then
        ln -s "$vcf_abs.tbi" "$out.tbi"
    elif [[ -f $vcf_abs.csi ]]; then
        ln -s "$vcf_abs.csi" "$out.csi"
    fi
else
    echo "chr$chrom: VCF samples differ from the analysis samples, subsetting and keeping MAC >= $min_mac"
    bcftools view -S "$samples" -Ou "$vcf" \
        | bcftools view -c "$min_mac:minor" --threads "$threads" -Oz -o "$out.tmp.vcf.gz"
    mv "$out.tmp.vcf.gz" "$out"
    bcftools index -t "$out"
    echo "chr$chrom: $(bcftools index -n "$out") variants kept"
fi

n_out=$(( $(bcftools query -l "$out" | wc -l) ))
n_list=$(grep -c . "$samples")
if [[ $n_out -ne $n_list ]]; then
    echo "ERROR: $out has $n_out samples, $samples has $n_list" >&2
    exit 1
fi
echo "chr$chrom: $n_out samples -> $out"
