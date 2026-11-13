#!/usr/bin/env python3
"""
Submit Tractor-Mix jobs: one SLURM job per chromosome per scan, or per group of
pieces when the inputs were split by 06_filter_split_tracts.sh.

writes a job script per task into logs/ and
sbatch's it, skipping outputs that already exist so the whole thing is
resumable. tractormix_patched/tractormix_single_chr.r writes <out>.tmp and renames it only once
verified, so an existing output is always complete.

Templates: {chrom} and {part} are filled here, {anc} by the R script
(0 .. n_anc-1), and {model} (N0 / N1) in the null-model template. {part} is
optional. With it, the pieces of each chromosome are found by globbing the anc0
dosage file; each piece is its own Tractor-Mix run, written to
<scan>_chr<N>_part<P>.tsv, and --parts-per-job pieces run one after another in
each job.

Usage:
    python3 07_submit_tractormix.py --scans N0_cond \\
        --dosage-template   'filt/chr{chrom}/chr{chrom}_x.part{part}.anc{anc}.dosage.txt.gz' \\
        --hapcount-template 'filt/chr{chrom}/chr{chrom}_x.part{part}.anc{anc}.hapcount.txt.gz'

        you can also submit individual chromosomes or a list of them
        potentially can also run the unconditional test, but for the purposes of this project
        we are sticking with the conditional (and the unconditional script is not patched yet)
        """

import argparse
import glob
import os
import shlex
import subprocess
import sys

# ---------------------------------------------------------------------------
SCANS = {
    "N0_uncond": ("N0", "uncond"),
    "N0_cond":   ("N0", "cond"),
    "N1_uncond": ("N1", "uncond"),
    "N1_cond":   ("N1", "cond"),
}
DEFAULT_SCANS = ["N0_uncond", "N0_cond"]

NULL_TEMPLATE = "gmmat_null/gmmat_null_{model}.rds"
DOSAGE_TEMPLATE = "tracts/chr{chrom}.anc{anc}.dosage.txt"
HAPCOUNT_TEMPLATE = "tracts/chr{chrom}.anc{anc}.hapcount.txt"
R_SCRIPT = "tractormix_patched/tractormix_single_chr.r"
TRACTORMIX_DIR = "tractormix_patched"
TM_SCRIPTS = {"uncond": "TractorMix.score.R", "cond": "TractorMix.score_cond.R"}
CONDA_ENV = "genesis_pcair"   # the name in environment.yml
PARTITION = "defq"

# Chromosomes 1-5 carry the most variants; give them more headroom.
BIG_CHRS = {1, 2, 3, 4, 5}


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


def fill_chrom(template: str, chrom: int) -> str:
    # str.replace, not str.format: {anc} must survive for the R script.
    return template.replace("{chrom}", str(chrom))


def find_parts(template: str, chrom: int) -> list[str]:
    """Piece labels on disk for this chromosome, from its anc0 files, sorted."""
    pattern = fill_chrom(template, chrom).replace("{anc}", "0")
    pre, post = pattern.split("{part}")
    hits = glob.glob(glob.escape(pre) + "*" + glob.escape(post))
    return sorted(h[len(pre):len(h) - len(post)] for h in hits)


def find_conda_sh(explicit: str) -> str:
    """The conda.sh each job sources so that `conda activate` works in it.

    --conda-sh if given; otherwise derived from the conda install running here:
    $CONDA_EXE (set in any shell where conda is initialized) or `conda info --base`.
    """
    candidates = []
    if explicit:
        candidates.append(explicit)
    else:
        conda_exe = os.environ.get("CONDA_EXE", "")
        if conda_exe:
            base = os.path.dirname(os.path.dirname(conda_exe))
            candidates.append(os.path.join(base, "etc", "profile.d", "conda.sh"))
        try:
            base = subprocess.run(["conda", "info", "--base"], capture_output=True,
                                  text=True, check=True).stdout.strip()
            candidates.append(os.path.join(base, "etc", "profile.d", "conda.sh"))
        except (OSError, subprocess.CalledProcessError):
            pass
    for path in candidates:
        if os.path.isfile(path):
            return os.path.abspath(path)
    sys.exit("ERROR: cannot find conda.sh"
             + (f" at {explicit}" if explicit else " from $CONDA_EXE or `conda info --base`")
             + ".\n  Pass --conda-sh <conda base>/etc/profile.d/conda.sh")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrs", default="1-22", help="e.g. '1-22' or '19' or '1,2,19'")
    ap.add_argument("--scans", nargs="+", default=DEFAULT_SCANS, choices=list(SCANS))
    ap.add_argument("--null-template", default=NULL_TEMPLATE)
    ap.add_argument("--dosage-template", default=DOSAGE_TEMPLATE)
    ap.add_argument("--hapcount-template", default=HAPCOUNT_TEMPLATE)
    ap.add_argument("--outdir", default="assoc_tm")
    ap.add_argument("--logdir", default="logs")
    ap.add_argument("--rscript", default=R_SCRIPT)
    ap.add_argument("--tractormix-dir", default=TRACTORMIX_DIR)
    ap.add_argument("--n-anc", type=int, default=3)
    ap.add_argument("--anc-labels", default="NAT,AFR,EUR")
    ap.add_argument("--ac-threshold", type=int, default=50)
    ap.add_argument("--chunk-size", type=int, default=2048)
    ap.add_argument("--cpus", type=int, default=4,
                    help="Tractor-Mix n_core: n_core - 1 worker processes.")
    ap.add_argument("--mem", default="32G", help="SLURM memory, e.g. 16G; a bare number is taken as GB.")
    ap.add_argument("--mem-big", default="64G", help="Memory for chromosomes 1-5 (same format).")
    ap.add_argument("--parts-per-job", type=int, default=25,
                    help="With {part} in the templates: pieces run one after another per job.")
    # ap.add_argument("--time", default="48:00:00")
    ap.add_argument("--conda-env", default=CONDA_ENV)
    ap.add_argument("--conda-sh", default="",
                    help="conda.sh sourced by each job. Default: found from the conda in use.")
    ap.add_argument("--partition", default=PARTITION,
                    help="SLURM partition; '' leaves it to the cluster default.")
    ap.add_argument("--mail-user", default="",
                    help="Email for END/FAIL notices. Default: SLURM's default (the submitting user).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Write the job scripts but do not sbatch them.")
    ap.add_argument("--force", action="store_true",
                    help="Resubmit even when the output already exists.")
    args = ap.parse_args()

    for name, tpl in (("--dosage-template", args.dosage_template),
                      ("--hapcount-template", args.hapcount_template)):
        if "{chrom}" not in tpl or "{anc}" not in tpl:
            sys.exit(f"ERROR: {name} must contain both {{chrom}} and {{anc}}: {tpl}")
        if tpl.count("{part}") > 1:
            sys.exit(f"ERROR: {name} contains {{part}} more than once: {tpl}")
    has_parts = "{part}" in args.dosage_template
    if has_parts != ("{part}" in args.hapcount_template):
        sys.exit("ERROR: use {part} in both --dosage-template and --hapcount-template, or in neither")
    if args.parts_per_job < 1:
        sys.exit("ERROR: --parts-per-job must be >= 1")
    if "{model}" not in args.null_template:
        sys.exit(f"ERROR: --null-template must contain {{model}}: {args.null_template}")
    # SLURM reads a bare --mem number as megabytes: --mem 16 would request 16M.
    for attr in ("mem", "mem_big"):
        value = getattr(args, attr)
        if value.isdigit():
            print(f"note: --{attr.replace('_', '-')} {value} has no unit; using {value}G")
            setattr(args, attr, f"{value}G")
    if len(args.anc_labels.split(",")) != args.n_anc:
        sys.exit(f"ERROR: --anc-labels has {len(args.anc_labels.split(','))} labels "
                 f"for --n-anc {args.n_anc}")

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.logdir, exist_ok=True)

    required = [args.rscript]
    for scan in args.scans:
        model, test = SCANS[scan]
        required.append(args.null_template.replace("{model}", model))
        required.append(os.path.join(args.tractormix_dir, TM_SCRIPTS[test]))
    absent = sorted({p for p in required if not os.path.exists(p)})
    if absent:
        sys.exit("ERROR: not found:\n  " + "\n  ".join(absent) +
                 "\n  (the null model comes from 02_fit_null_binary.r)")
    conda_sh = find_conda_sh(args.conda_sh)

    jobs = runs = skipped = missing = 0

    # Descending so the big chromosomes queue first.
    for chrom in sorted(parse_chrs(args.chrs), reverse=True):
        parts: list[str | None] = find_parts(args.dosage_template, chrom) if has_parts else [None]
        if not parts:
            print(f"chr{chrom}: no pieces match {fill_chrom(args.dosage_template, chrom)}; skipping")
            missing += 1
            continue

        for scan in args.scans:
            model, test = SCANS[scan]
            pending = []  # (tag, dosage, hapcount, out) for each run still to do
            n_exist = n_absent = 0
            first_absent = ""
            for part in parts:
                dosage = fill_chrom(args.dosage_template, chrom)
                hapcount = fill_chrom(args.hapcount_template, chrom)
                tag = f"chr{chrom}"
                if part is not None:
                    dosage = dosage.replace("{part}", part)
                    hapcount = hapcount.replace("{part}", part)
                    tag += f"_part{part}"

                # The conditional test reads K-1 hapcount files; the R script
                # leaves out the last ancestry by default.
                needed = [dosage.replace("{anc}", str(k)) for k in range(args.n_anc)]
                if test == "cond":
                    needed += [hapcount.replace("{anc}", str(k)) for k in range(args.n_anc - 1)]
                absent = [p for p in needed if not os.path.exists(p)]
                if absent:
                    n_absent += 1
                    first_absent = first_absent or absent[0]
                    continue

                out = os.path.join(args.outdir, f"{scan}_{tag}.tsv")
                if os.path.exists(out) and not args.force:
                    n_exist += 1
                    continue
                pending.append((tag, dosage, hapcount, out))

            if n_absent:
                print(f"chr{chrom} {scan}: {n_absent} run(s) with input file(s) not found, "
                      f"e.g. {first_absent}; skipping")
                missing += n_absent
            if n_exist:
                print(f"chr{chrom} {scan}: {n_exist} output(s) exist, skipping")
                skipped += n_exist

            mem = args.mem_big if chrom in BIG_CHRS else args.mem
            step = args.parts_per_job if has_parts else 1
            for start in range(0, len(pending), step):
                group = pending[start:start + step]
                job = f"tm_{scan}_{group[0][0]}"
                if len(group) > 1:
                    job += "-" + group[-1][0].rsplit("part", 1)[1]
                pbs = os.path.join(args.logdir, f"{job}.pbs")

                with open(pbs, "w") as fh:
                    fh.write("#!/bin/sh\n")
                    fh.write("#SBATCH --mail-type=END,FAIL\n")
                    if args.mail_user:
                        fh.write(f"#SBATCH --mail-user={args.mail_user}\n")
                    fh.write("#SBATCH --ntasks=1\n")
                    fh.write(f"#SBATCH --cpus-per-task={args.cpus}\n")
                    fh.write(f"#SBATCH --mem={mem}\n")
                    # fh.write(f"#SBATCH --time={args.time}\n")
                    if args.partition:
                        fh.write(f"#SBATCH --partition={args.partition}\n")
                    fh.write(f"#SBATCH --job-name={job}\n")
                    fh.write(f"#SBATCH -o {args.logdir}/{job}.out\n")
                    fh.write(f"#SBATCH -e {args.logdir}/{job}.err\n\n")
                    # "." rather than "source": the job runs under /bin/sh.
                    fh.write(f". {shlex.quote(conda_sh)}\n")
                    fh.write(f"conda activate {args.conda_env}\n\n")
                    # Tractor-Mix parallelizes over worker processes; a multithreaded
                    # BLAS in each worker would oversubscribe the allocation.
                    fh.write("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n\n")
                    # A failed run does not stop the rest of the group; the job
                    # still exits non-zero so SLURM reports FAIL.
                    fh.write("status=0\n\n")
                    for _, dosage, hapcount, out in group:
                        fh.write(
                            f"Rscript {shlex.quote(args.rscript)} \\\n"
                            f"    --null-model {shlex.quote(args.null_template.replace('{model}', model))} \\\n"
                            f"    --dosage {shlex.quote(dosage)} \\\n"
                            f"    --hapcount {shlex.quote(hapcount)} \\\n"
                            f"    --test {test} \\\n"
                            f"    --out {shlex.quote(out)} \\\n"
                            f"    --n-anc {args.n_anc} \\\n"
                            f"    --anc-labels {shlex.quote(args.anc_labels)} \\\n"
                            f"    --ac-threshold {args.ac_threshold} \\\n"
                            f"    --chunk-size {args.chunk_size} \\\n"
                            f"    --cores {args.cpus} \\\n"
                            f"    --tractormix-dir {shlex.quote(args.tractormix_dir)} \\\n"
                            f"    || status=1\n\n"
                        )
                    fh.write("exit $status\n")

                if args.dry_run:
                    print(f"{job}: wrote {pbs} ({len(group)} run(s), dry run)")
                else:
                    subprocess.run(["sbatch", pbs], check=True)
                    print(f"{job}: submitted ({len(group)} run(s))")
                jobs += 1
                runs += len(group)

    print(f"\njobs {jobs} ({runs} runs), existing outputs skipped {skipped}, "
          f"missing inputs {missing}")
    if not args.dry_run and jobs:
        print("\nWhen the jobs finish, combine the pieces with 08_combine_tractormix.sh "
              "and split the result with 09_split_sumstats_by_anc.sh.")

if __name__ == "__main__":
    main()
