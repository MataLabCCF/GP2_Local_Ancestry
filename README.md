# GP2 Local Ancestry

*This repository contains the code and overall pipeline that was designed in order to implement local ancestry inference in GP2's release 10 data and to run a local-ancestry-aware genome-wide association study (plus downstream analysis).*

---
## Reference panel assembly 

---

## QC of non-European populations altogether
We must first know which samples are labeled to belong to the EUR label assigned by genotools.

```
##get the list of EUR samples 
awk 'NR>1 {print $2}' /config/workspace/gp2_tier2_eu_release10/raw_genotypes/EUR/EUR_release10_vwb.psam > fullEUR_POPlist.txt


##make sure the LA ref samples are in the all EUR list 
grep -Ff /config/workspace/ws_files/LA_noAJ/Best1/reference4ProjectedPCA/runProjectedPCA/plotProjectedPCA/parentals/EUR.txt \
fullEUR_POPlist.txt | wc -l 
```
And then get the list of samples from EUR that do not belong to the reference panel. Because the one that belong to the ref panel must go inside the QC with all the other non-EUR samples. Because ultimately they need to go into topmed for phasing. Since the reference panel for local ancestry inference, and the target data must be in the same phase. 

```
##get the list of EUR samples WITHOUT the members of the LA reference panel 
awk 'NR==FNR {remove[$1]; next} !($1 in remove)' \
/config/workspace/ws_files/LA_noAJ/Best1/reference4ProjectedPCA/runProjectedPCA/plotProjectedPCA/parentals/EUR.txt \
fullEUR_POPlist.txt \
> fullEUR_POPlist_noLAref.txt
```
Now, the exclusion of the EUR samples (that are also not part of the reference sample), from the merged dataset
```
plink2   --exclude ultimate_list_snps2remove_nba_qc.txt \
  --make-pgen \
  --memory 50000 \
  --out GP2_r10_final_with_gp2id_ready_genotools_noEUR_exceptLAref \
  --pfile /config/workspace/gp2_makariousmb/R10/GP2_r10_final_with_gp2id_ready_genotools \
  --remove fullEUR_POPlist_noLAref.txt \
  --threads 8
```
The plink2 files named GP2_r10_final_with_gp2id_ready_genotools* are the ones that normally would be QCed with genotools (to partitioning them to the 11 genotools's populations). However, for this particular analysis, the working dataset must have enough representation of the ancestries that are included in the association model. And if we use the 11 populations derived from genotools, we cannot achieve enough representation of other ancestries within each cluster (not even at continental level). 

That is why we must perform the QC with the samples altogether, and thereafter --> imputation with topmed --> local ancestry inference --> GWAS --> post GWAS.

Since we are excluding the EUR samples (as already classified by genotools for GP2 default QC pipeline), the analysis would be, in principle, more computationally affordable. 

For the QC we are going to use a pipeline designed for admixed populations

---

# Tractor-Mix, binary phenotype

Case/control Tractor-Mix run. The null model is a logistic mixed model fitted
directly with GMMAT, using a dense GCTA GRM
built from typed variants.

Run the steps in numerical order, from this folder.

| Step | Script | Reads | Writes |
|---|---|---|---|
| 0 | `00_define_samples.py` | covariates, .psam, msp, VCF | `samples/analysis_samples.txt` |
| 1 | `01_make_grm_gcta.sh` | typed pfile, step 0 | `gcta_grm/grm.*` |
| 2 | `02_fit_null_binary.r` | covariates, step 0, step 1 | `gmmat_null/gmmat_null_N0.rds` |
| 3 | `03_subset_vcf.sh` | imputed VCFs, step 0 | `vcf_files/chr<N>.analysis.vcf.gz` |
| 4 | `04_prep_msp.py` | step 3, msp, step 0 | `msp_ordered/chr<N>.query_results.ordered.msp` |
| 5 | `05_extract_tracts.sh` (+ `extract_tracts.c`) | steps 3 and 4 | `extracted_tracts/chr<N>/` |
| 6 | `06_filter_split_tracts.sh` | step 5 | `tracts_filtered/chr<N>/...part<P>.anc<k>.*.txt.gz` |
| 7 | `07_submit_tractormix.py` (+ `tractormix_patched/`) | steps 2 and 6 | `assoc_tm/N0_cond_chr<N>_part<P>.tsv` |
| 8 | `08_combine_tractormix.sh` | steps 6 and 7 | `combined_parts/N0_cond_genome.tsv.gz` |
| 9 | `09_split_sumstats_by_anc.sh` | steps 6 and 8 | `combined_parts/N0_cond_*.tsv.gz` + summary |

## More ancestries: what to change

This particular example run used 3 ancestries (NAT, AFR, EUR). With K ancestries, **only steps 7
and 9 need different arguments.** Every other step reads K from the data.

| Step | Change needed? | What to do |
|---|---|---|
| 0–3 | No | Ancestry plays no role yet. |
| 4 `04_prep_msp.py` | No, but **read its output** | It prints the msp code line, e.g. `NAT=0 AFR=1 EUR=2`. That line defines what anc0, anc1, ... mean and the order of every label list below. It also flags ancestries below 10% mean local ancestry. |
| 5 `05_extract_tracts.sh` | No | Counts K from the msp code line and passes `--num-ancs K`. |
| 6 `06_filter_split_tracts.sh` | No (optional speed-up) | Counts K from the `MAC_anc*` columns. It splits all 2K files (hap and dosage files) at once, so raise `#SBATCH --cpus-per-task` (8 now) to about 2K. |
| **7 `07_submit_tractormix.py`** | **Yes** | Set `--n-anc K` and `--anc-labels` with K labels **in the step 4 order**. If you leave the defaults (`3`, `NAT,AFR,EUR`), the script finds anc0–anc2, runs, and silently ignores every ancestry after the third. The conditional scan reads 2K − 1 files per piece (K dosage + K − 1 hapcount), so check `--mem` too. |
| 8 `08_combine_tractormix.sh` | No | Doesn't depend on ancestry. |
| **9 `09_split_sumstats_by_anc.sh`** | **Yes** | Set `ANC_LABELS=` to the same K labels in the same order. The script stops if the label count doesn't match the results. It can't catch the right count in the wrong order, which mislabels the per-ancestry files. |

The same two settings apply with 3 ancestries in another order or with other
names: always copy the order from the step 4 code line.

What changes in the results with more (and smaller) ancestries:
- Tractor-Mix reports the joint K-df test only when every ancestry passes
  `--ac-threshold`, so more variants end up in `_partial` and fewer in `_joint`.
- Small ancestries give many `not_tested` / `below_mac` rows in the step 9
  summary. The Tractor-Mix authors suggest modelling only ancestries with at
  least 10% mean ancestry. Lower (around 4–5%) has worked for us, with a
  stricter MAC filter.

## Sample consistency

Samples can be missing from some files (e.g. genotyped but no covariates, or the
reverse), and files can list samples in different orders. Where each risk is handled:

| Risk | Where it's handled |
|---|---|
| Samples with covariates but no genotypes, msp or VCF | `00_define_samples.py` keeps only samples present in all of them; the GRM and the null use that set |
| Samples with genotypes but not in the null | `03_subset_vcf.sh` makes each VCF hold exactly that set, so no extra people reach extract_tracts or the MAC counts |
| msp and VCF in different sample orders | `04_prep_msp.py` reorders the msp to the VCF's order. This matters because extract_tracts pairs them by position. |
| Dosage files and null in different orders | Tractor-Mix matches them by sample ID (`match(obj$id_include, sample_id)`), so this order doesn't matter |

Because only the set matters after step 0, step 3 does not rewrite a VCF whose
samples already match: it links to it instead with a symlink

## Covariates: `model_covariates.txt`

The model's covariates are listed once, in `model_covariates.txt`. You will have to edit this file manually. Both
`00_define_samples.py` (missing-value and number checks) and
`02_fit_null_binary.r` (the model) read this file, so they always agree.

```
SEX     categorical
PCs     20
# AGE   quantitative   (example of a quantitative covariate)
```

- `categorical`: entered as a factor. Each value gets a coefficient, except the
  first in sorted order, which is the reference. The log prints the reference
  and a cases/controls table. The script stops if there is only one value, or
  if any value has no cases or no controls.
- `quantitative`: entered as a number, with one linear coefficient. The script
  stops if a value is not a number, or if every sample has the same value.
- `PCs N`: uses the first N PCs (`GCTA_PC1..` or `PC1..`). `PCs 0` means no PCs.

To add a covariate, add a line to the file and rerun from step 0. A new covariate
can bring new missing values, which changes the sample set.

## Run order

These are the commands of the test run, with the long cluster paths replaced by
short placeholders. **Replace only these four inputs with your own files:**

| Placeholder | What it is |
|---|---|
| `covariates.tsv` | tab-delimited: IID, STATUS and the columns in `model_covariates.txt` |
| `typed/merged_typed` | plink2 pfile prefix of the typed variants (`.pgen/.pvar/.psam`) |
| `imputed/chr{chrom}.dose.vcf.gz` | phased, imputed VCF per chromosome; keep the literal `{chrom}` |
| `gnomix/LPD_chrom{chrom}/query_results.msp` | gnomix (or RFMix2) .msp per chromosome; keep the literal `{chrom}` |

The other names (`samples`, `gcta_grm`, `gmmat_null`, `vcf_files`,
`msp_ordered`, `extracted_tracts`, `tracts_filtered`, `assoc_tm`,
`combined_parts`) are output folders. You can rename them, but each one is
an input to a later step, so rename it there too.

Keep the single quotes around templates with `{chrom}`, `{part}`, `{anc}` or
`{model}`: the scripts fill them in, not the shell.

The commands below cover all 22 chromosomes. The test run used only chromosomes
21 and 22 (I was testing them in LARGE-PD): `--array=21-22` in steps 3, 5 and 6, `--chrs 21,22` in step 7 and
`CHRS="21 22"` in step 8. Steps 0 and 4 ran on all chromosomes (their default).

**0. Analysis samples**

Before the first step, clone the repository and create the conda environment
from `environment.yml`. It holds the versions used in the test run (Python,
R and the R packages, bcftools, plink2 and GCTA). GMMAT is the one exception,
installed from CRAN right after:
```bash
git clone https://github.com/MataLabCCF/GP2_Local_Ancestry
cd GP2_Local_Ancestry/tractor_mix_binary

conda env create -f environment.yml     # creates the environment "genesis_pcair" (happens to be the one I already had with the R packages needed)
conda activate genesis_pcair

# GMMAT (step 2) has no conda package; install it from CRAN into the environment.
# It compiles C++ code with the compilers that come with the environment.
Rscript -e 'install.packages("GMMAT", repos = "https://cloud.r-project.org")'
```
To give it another name, use `conda env create -n <name> -f environment.yml`,
and pass `--conda-env <name>` in step 7. The file was exported on a Linux
x86-64 cluster, so it creates on Linux, not on macOS or Windows.

Check that everything the scripts load is there:
```bash
Rscript -e 'for (p in c("GMMAT", "data.table", "R.utils", "Matrix", "doParallel", "foreach", "dplyr", "abind"))
  cat(sprintf("%-11s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING"))'
python3 -c 'import numpy; print("numpy", numpy.__version__)'
bcftools --version | head -n 1
plink2 --version
gcta64 2>&1 | grep -m 1 version
```
If an R package shows `MISSING`, install it into the active environment: GMMAT
with the `install.packages` line above, the others with e.g.
`conda install -c conda-forge r-doparallel r-abind r-r.utils`.

**SLURM settings.** The scripts are set up for a partition named `defq`
(`#SBATCH --partition=defq` in `01`, `03`, `05` and `06`). If your cluster uses
another one, override it when submitting, without editing the scripts: e.g.
`sbatch --partition=<name> --array=1-22 03_subset_vcf.sh ...`, and
`--partition <name>` in step 7 (`--partition ''` uses the cluster default).
Job emails (END, FAIL) go to SLURM's default address for the submitting user.
To send them elsewhere, add `--mail-user=<email>` to `sbatch`, or
`--mail-user <email>` in step 7.

Each step 7 job activates the environment by sourcing your conda's `conda.sh`.
The script finds it from the conda you are using when you submit, so submit
step 7 from a shell where conda works. If it can't find it, it stops and asks for
`--conda-sh <conda base>/etc/profile.d/conda.sh`.

Then define the samples:
```bash
python3 00_define_samples.py \
    --covar        covariates.tsv \
    --psam         typed/merged_typed.psam \
    --msp-template 'gnomix/LPD_chrom{chrom}/query_results.msp' \
    --vcf-template 'imputed/chr{chrom}.dose.vcf.gz' \
    --model        model_covariates.txt \
    --outdir       samples 2>&1 | tee 00_define_samples.log
```
Writes `samples/analysis_samples.txt` (one IID per line) and
`samples/sample_report.tsv` (every sample seen, where, and whether it was kept).

**1. Dense GCTA GRM**
```bash
bash 01_make_grm_gcta.sh \
    typed/merged_typed samples/analysis_samples.txt gcta_grm 2>&1 | tee 01_make_grm_gcta.log
```
The environment provides `gcta64` and `plink2`. To use other copies, put
`GCTA=/path/to/gcta` or `PLINK2=/path/to/plink2` in front of the command, as the
test run did with `GCTA=`. The script also has `#SBATCH` headers, so `sbatch` works
in place of `bash`.

**2. Null model**
```bash
Rscript 02_fit_null_binary.r \
    --covar   covariates.tsv \
    --model   model_covariates.txt \
    --samples samples/analysis_samples.txt \
    --grm     gcta_grm/grm \
    --out     gmmat_null/gmmat_null_N0.rds 2>&1 | tee 02_fit_null_binary.log
```
The null is named `..._N0.rds` only because step 7 finds it through the
`{model}` template. Here N0 is simply the binary null.

**3. VCF subset to the analysis samples**
```bash
mkdir -p logs
sbatch --array=1-22 03_subset_vcf.sh \
    'imputed/chr{chrom}.dose.vcf.gz' samples/analysis_samples.txt vcf_files
```
For each chromosome, the script compares the VCF header with the sample list.
If they hold the same samples (in any order), `vcf_files/chr<N>.analysis.vcf.gz`
is a symlink to the original VCF and nothing is rewritten. Otherwise the VCF is
subset with `bcftools view -S`, which can take hours on large chromosomes. The log
says which of the two happened.

When the VCF is subset, variants with MAC < 21 in the analysis samples are
dropped in the same pass (`MIN_MAC=...` changes the threshold). This can't
remove a variant that step 6 would keep, because MAC > 20 in any ancestry
implies MAC > 20 overall. We add this here because we seize the opportunity of opening already the files to subset it to the analysis samples. A linked VCF is not filtered. Step 6 applies the real
filter either way in an ancestry specific manner.

**4. Reorder the msp to the VCF sample order**
```bash
python3 04_prep_msp.py \
    --vcf-template 'vcf_files/chr{chrom}.analysis.vcf.gz' \
    --msp-template 'gnomix/LPD_chrom{chrom}/query_results.msp' \
    --samples      samples/analysis_samples.txt \
    --outdir       msp_ordered 2>&1 | tee 04_prep_msp.log
```
Point `--vcf-template` at the step 3 output, which is the VCF that step 5 reads.
The test run passed the original imputed VCFs here. That worked only because they
already held exactly the analysis samples, so step 3 had linked them.

At the end it prints the ancestry code line and the mean local ancestry per
ancestry, e.g.:

**This order defines anc0, anc1, ... from here on.** Use it for `--anc-labels` in
step 7 and `ANC_LABELS` in step 9.

**5. extract_tracts**

Compile the C version of Tractor's `extract_tracts.py` once (same output, much
faster). It needs gcc and zlib:
```bash
gcc -O2 -o extract_tracts ./extract_tracts.c -lz

sbatch --array=1-22 05_extract_tracts.sh vcf_files/ msp_ordered/ extracted_tracts 2>&1 | tee 05_extract_tracts.log
```
The script looks for the binary as `./extract_tracts` in the folder you submit
from. If it is elsewhere, submit with `EXTRACT_TRACTS=/path/to/extract_tracts sbatch ...`.
The number of ancestries comes from the msp code line, so nothing changes with more than 3.

**6. Filter (MAC > 20 in any ancestry) and split into pieces**
```bash
sbatch --array=1-22 ./06_filter_split_tracts.sh extracted_tracts tracts_filtered
```
This keeps a variant when its MAC is above 20 in at least one ancestry. For more ancestries we can consider to bump it up to 50, for example. It then
splits each chromosome into pieces of 16,384 variants, so step 7 can run many
jobs per chromosome. Optional positional arguments: `MAC_THRESHOLD`
(default 20) and `PIECE_SIZE` (default 16384; keep it a multiple of the
Tractor-Mix `--chunk-size`, 2048). `tracts_filtered/chr<N>/*.filter_summary.tsv`
records the variants kept and the number of pieces. Step 8 uses it to check
that every piece has a result.

This filter is per ancestry. Step 3's filter is on the whole sample and runs
only when the VCF had to be subset, so for linked VCFs this is the first MAC
filter.

**7. Tractor-Mix**
```bash
python3 ./07_submit_tractormix.py \
    --scans N0_cond \
    --null-template     'gmmat_null/gmmat_null_{model}.rds' \
    --dosage-template   'tracts_filtered/chr{chrom}/chr{chrom}.analysis.part{part}.anc{anc}.dosage.txt.gz' \
    --hapcount-template 'tracts_filtered/chr{chrom}/chr{chrom}.analysis.part{part}.anc{anc}.hapcount.txt.gz' \
    --rscript           tractormix_patched/tractormix_single_chr.r \
    --tractormix-dir    tractormix_patched/ \
    --n-anc 3 --anc-labels NAT,AFR,EUR \
    --chrs 1-22 \
    --parts-per-job 1 \
    --ac-threshold 20 \
    --cpus 8 \
    --mem 16G \
    --outdir assoc_tm 2>&1 | tee 07_submit_tractormix.log
```
- **`--n-anc` and `--anc-labels`: set them from the step 4 code line** (see
  [More than 3 ancestries](#more-than-3-ancestries-what-to-change)).
- `--rscript` and `--tractormix-dir` point to `tractormix_patched/`. That is
  also their default, so you can leave both out when running from this folder.
  Each job log
  (`logs/tm_N0_cond_chr<N>_part<P>.out`) should show
  `TractorMix.score_cond Script Version : 0.0.1-patched-singular`, which
  confirms the patched file was used.
- Only the conditional scan (`N0_cond`) is set up. `N0_uncond` would also need
  the original `TractorMix.score.R` in `--tractormix-dir`, and that file is not
  patched (and honestly we care mostly about the conditional model because it etells us the effect of a given variant in a LA background independently of the number of ancestral haplotypes in the region). 
- `--parts-per-job 1` submits one job per piece. A larger value runs that many
  pieces one after another in each job, which means fewer and longer jobs. 1 part per job is as much parallelization as possible. 
- **Memory:** `--mem` applies to chromosomes 6–22 and `--mem-big` (default 64G) to
  chromosomes 1–5. With one job per piece, `--mem-big` doesn't really matter.
  Every piece has the same number of variants, so a chr1 job needs the same
  memory as a chr22 job, and `--mem 16G` was enough (or maybe even lower?). You can set `--mem-big` to
  the same value as `--mem` to avoid asking for more memory than the job needs.
- The script skips pieces whose output already exists. If some jobs fail,
  rerun the same command and only the missing pieces are submitted.
  `--dry-run` writes the job scripts to `logs/` without submitting them.

**8. Combine the pieces into one genome-wide file**
```bash
bash 08_combine_tractormix.sh \
    assoc_tm N0_cond tracts_filtered/ combined_parts/N0_cond_genome.tsv.gz 2>&1 | tee 08_combine_tractormix.log
```
By default it combines chromosomes 1–22. For a subset, put `CHRS` in front of
the command, as in the test run: `CHRS="21 22" bash 08_combine_tractormix.sh ...`. The script stops and writes nothing if any piece is missing,
or if a result row doesn't match the kept variants of step 6.

**9. Split into analysis-ready sumstats**
```bash
ANC_LABELS=NAT,AFR,EUR bash 09_split_sumstats_by_anc.sh \
    combined_parts/N0_cond_genome.tsv.gz tracts_filtered/ 2>&1 | tee 09_split_sumstats.log
```
**`ANC_LABELS` must match `--anc-labels` from step 7.** NAT,AFR,EUR is the
default, so the test run left it out. The script processes every chromosome in
the step 8 file, so `CHRS` has no effect here. Outputs go next to the input, with
the `_genome` suffix dropped from the prefix:

| File | Contents |
|---|---|
| `N0_cond_joint.tsv.gz` | joint test available and every ancestry usable; one wide row per variant |
| `N0_cond_partial.tsv.gz` | at least one ancestry usable, but not all (joint test missing) (mainly for logging purposes) |
| `N0_cond_joint_partial.tsv.gz` | both of the above, in genomic order, with a `class` column |
| `N0_cond_anc<k>_<LABEL>.tsv.gz` | long sumstats for one ancestry, its usable variants only |
| `N0_cond_dropped.tsv.gz` | variants with no usable ancestry, with a `drop_reason` |
| `N0_cond_split_summary.tsv` | counts per chromosome and in total |

An ancestry is usable for a variant when Tractor-Mix tested it, returned an
effect and SE, and its MAC in that ancestry is above 20 (`MAC_THRESHOLD=`).

### Reading the step 9 screen output

At the end, the script lists the files it wrote and prints the header and the
`total` row of `N0_cond_split_summary.tsv`. The per-chromosome rows are only
in the file.

- `rows_in`: variants that came out of Tractor-Mix. It equals
  `joint + partial + dropped_no_usable`.
- `joint` / `partial` / `dropped_no_usable`: how many went to each file.
- `dropped_*` (not_tested, singular, below_mac, mixed): why the dropped
  variants had no usable ancestry. `singular` means the variance matrix couldn't
  be inverted, which removes every ancestry of that variant at once.
- For each ancestry, `kept_X` is the number of rows in its `_anc` file, and
  `not_tested_X`, `singular_X` and `below_mac_X` say why the rest were left out.
  These four add up to `rows_in`.
- `se_pre_X` / `se_post_X`: variants with an extreme SE (> 1000) before and
  after filtering. `se_post` must be 0. If it isn't, the script stops, because
  those pieces were run with the unpatched Tractor-Mix.





## AI use disclosure

Generative AI (Claude, Anthropic, via Claude Code) was used to help write and
document the code in this repository, including wrapper and SLURM submission
scripts, input checks, command-line examples and comments. The study design,
methods, parameter choices, and scientific interpretation are the author's own. The authors have
reviewed, edited and tested the code, and take full responsibility for the
results and for any inconsistencies in the code.



