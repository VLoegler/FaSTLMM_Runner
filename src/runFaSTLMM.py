#!/usr/bin/env python
# -*- coding: utf8 -*-
# ---------------------------------------------------------------------------
# Created By  : vloegler
# Creation Date: 2022/11/04
# Modified On : 2026/04/01
# Version     : 4.0
# ---------------------------------------------------------------------------
"""
Run Genome Wide Association Study using FaST-LMM (v0.6.4).
Performs association and permutation tests to identify significant SNPs.

Output:
    <outdir>/<phen_name>.first_assoc.txt.gz  -- FaST-LMM association results
    <outdir>/<phen_name>.threshold.txt        -- Permutation-based p-value threshold
    <outdir>/<phen_name>.signif_snps.txt      -- SNPs below threshold
    <outdir>/<phen_name>.log                  -- Human-readable log
"""

import argparse
import glob
import logging
import math
import os
import shutil
import time
from pathlib import Path
from typing import Set, Tuple

import numpy as np
import pandas as pd
from fastlmm.association import single_snp
from pysnptools.snpreader import SnpData
from pysnptools.util.mapreduce1.runner import LocalMultiProc

BATCH_SIZE = 101  # FaST-LMM max phenotypes per call (1 real + 100 permutations)

logger = logging.getLogger("runFaSTLMM")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GWAS using FaST-LMM with optional permutation testing."
    )
    parser.add_argument(
        "-g", "--genotype", required=True, help="Genotype matrix (Plink .bed prefix)"
    )
    parser.add_argument(
        "-k",
        "--kinship",
        type=str,
        default="",
        help="Kinship matrix (Plink .bed prefix)",
    )
    parser.add_argument(
        "-p", "--phenotypes", required=True, nargs="+", help="Phenotype file(s)"
    )
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument(
        "-c", "--covariance", type=str, default="", help="Covariance matrix (optional)"
    )
    parser.add_argument(
        "-n",
        "--nb-permutations",
        type=int,
        default=100,
        help="Number of permutations (default: 100, 0 to disable)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads"
    )
    parser.add_argument(
        "--uncompressed",
        action="store_true",
        help="Write first association results as plain .txt instead of .txt.gz",
    )
    parser.add_argument(
        "--double-id",
        action="store_true",
        help="Set both FID and IID to the sample name. "
        "Without this flag, sample names are split on '_' into FID/IID "
        "(mimics plink 1.9 default behavior).",
    )
    return parser.parse_args()


def setup_logging(log_file: Path) -> logging.FileHandler:
    """Configure logging to both stderr and a file. Returns the file handler."""
    logger.setLevel(logging.DEBUG)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

    fh = logging.FileHandler(str(log_file), mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )

    logger.addHandler(console)
    logger.addHandler(fh)
    return fh


def build_iid_array(strain_index: pd.Index, *, double_id: bool) -> np.ndarray:
    if double_id:
        pairs = [(s, s) for s in strain_index]
    else:
        pairs = [(s.split("_", 1) if "_" in s else (s, s)) for s in strain_index]
    return np.array(pairs, dtype="<U23")


def read_fam_iids(bed_file: str) -> Set[str]:
    """Read IIDs from the .fam file associated with a plink bed prefix."""
    bed_path = bed_file if bed_file.endswith(".bed") else bed_file + ".bed"
    fam_path = bed_path.replace(".bed", ".fam")
    iids = set()
    with open(fam_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                iids.add(parts[1])
    return iids


def log_sample_overlap(phen_iids: Set[str], geno_iids: Set[str]) -> None:
    """Warn about samples present in only one of phenotype / genotype."""
    in_phen_only = sorted(phen_iids - geno_iids)
    in_geno_only = sorted(geno_iids - phen_iids)
    common = phen_iids & geno_iids

    logger.info(
        "%d samples in phenotype file, %d in genotype file, %d in common",
        len(phen_iids),
        len(geno_iids),
        len(common),
    )

    if in_phen_only:
        logger.warning(
            "%d sample(s) in phenotype but NOT in genotypes (will be ignored): %s",
            len(in_phen_only),
            ", ".join(in_phen_only),
        )
    if in_geno_only:
        logger.warning(
            "%d sample(s) in genotypes but NOT in phenotype (will be ignored): %s",
            len(in_geno_only),
            ", ".join(in_geno_only),
        )
    if not common:
        logger.error("No samples in common between phenotype and genotypes!")


def shuffle_columns(matrix: np.ndarray) -> None:
    for col in range(matrix.shape[1]):
        np.random.shuffle(matrix[:, col])


def build_phenotype_data(df: pd.DataFrame, nperm: int) -> Tuple[np.ndarray, np.ndarray]:
    phen_name = df.columns[0]
    if nperm == 0:
        return np.asarray(df), np.array([phen_name])

    values = np.asarray(df.values)
    permut_matrix = np.repeat(values[:, [0]], repeats=nperm, axis=1)
    shuffle_columns(permut_matrix)
    merged = np.hstack((values, permut_matrix))
    headers = np.array([phen_name] + ["perm%d" % i for i in range(1, nperm + 1)])
    return merged, headers


def run_single_snp(
    bed_file: str,
    kinship_file: str,
    covar_file: str,
    pheno: SnpData,
    runner: LocalMultiProc,
) -> pd.DataFrame:
    kwargs = dict(
        test_snps=bed_file,
        pheno=pheno,
        count_A1=True,
        runner=runner,
        show_snp_fract_var_exp=True,
    )  # type: dict
    if kinship_file:
        kwargs["G0"] = kinship_file
    if covar_file:
        kwargs["covar"] = covar_file
    return single_snp(**kwargs)


def run_gwas(
    bed_file: str,
    kinship_file: str,
    pheno_file: str,
    covar_file: str,
    outdir: str,
    nperm: int,
    threads: int,
    *,
    double_id: bool,
    uncompressed: bool = False,
) -> None:
    # Resolve all paths before os.chdir() changes the working directory.    
    bed_file = str(Path(bed_file).resolve())
    if kinship_file:
        kinship_file = str(Path(kinship_file).resolve())
    if covar_file:
        covar_file = str(Path(covar_file).resolve())
    pheno_file = str(Path(pheno_file).resolve())

    df = pd.read_csv(pheno_file, sep="\t", index_col=0)
    phen_name = df.columns[0]

    phen_dir = (Path(outdir) / phen_name).resolve()
    if phen_dir.exists():
        raise FileExistsError(
            "Output directory '%s' already exists. Provide a new one." % phen_dir
        )
    phen_dir.mkdir(parents=True)

    # Set up logging to file + console for this phenotype
    file_handler = setup_logging(phen_dir / (phen_name + ".log"))

    logger.info("=== GWAS for phenotype: %s ===", phen_name)
    logger.info("Phenotype file: %s", pheno_file)
    logger.info("Genotype file:  %s", bed_file)
    if kinship_file:
        logger.info("Kinship file:   %s", kinship_file)
    if covar_file:
        logger.info("Covariance file: %s", covar_file)
    logger.info("Permutations: %d", nperm)
    logger.info("Threads: %d", threads)
    logger.info(
        "ID mode: %s",
        "double-id (FID=IID=name)" if double_id else "split on '_' (plink 1.9 default)",
    )

    iid = build_iid_array(df.index, double_id=double_id)

    # Check sample overlap between phenotype and genotype
    phen_iids = set(iid[:, 1])
    geno_iids = read_fam_iids(bed_file)
    log_sample_overlap(phen_iids, geno_iids)

    merged_values, headers = build_phenotype_data(df, nperm)

    os.chdir(phen_dir)

    runner = LocalMultiProc(threads)
    best_pvalues = []  # type: List[float]
    first_association = None  # type: Optional[pd.DataFrame]
    n_total = nperm + 1
    n_batches = math.ceil(n_total / BATCH_SIZE)

    for batch_idx in range(n_batches):
        start = batch_idx * BATCH_SIZE
        end = min(start + BATCH_SIZE, n_total)

        if nperm > 0:
            batch_values = np.c_[merged_values[:, start:end]]
        else:
            batch_values = np.c_[merged_values[:]]

        pheno_snp = SnpData(iid=iid, sid=headers[start:end], val=batch_values)

        t0 = time.time()
        results_df = run_single_snp(
            bed_file, kinship_file, covar_file, pheno_snp, runner
        )
        elapsed = time.time() - t0
        logger.info(
            "Batch %d/%d (%d phenotypes): association took %.1f seconds",
            batch_idx + 1,
            n_batches,
            end - start,
            elapsed,
        )

        if nperm > 0:
            batch_headers = set(headers[start:end])
            if phen_name in batch_headers:
                first_association = results_df[results_df.Pheno == phen_name]
                permutations = results_df[results_df.Pheno != phen_name]
            else:
                permutations = results_df

            top_per_perm = permutations.drop_duplicates(subset=["Pheno"])
            top_per_perm = top_per_perm[top_per_perm.Pheno != phen_name]
            best_pvalues.extend(top_per_perm.PValue)
        else:
            first_association = results_df

    if first_association is not None:
        suffix = ".first_assoc.txt" if uncompressed else ".first_assoc.txt.gz"
        out_assoc = phen_dir / (phen_name + suffix)
        first_association.to_csv(out_assoc, index=False, sep="\t")
        logger.info(
            "Association results written to %s (%d SNPs tested)",
            out_assoc.name,
            len(first_association),
        )

    if nperm > 0:
        threshold = np.percentile(best_pvalues, 5)
        (phen_dir / (phen_name + ".threshold.txt")).write_text("x\n%s\n" % threshold)
        logger.info(
            "Permutation threshold (5th percentile of %d best p-values): %.2e",
            len(best_pvalues),
            threshold,
        )

        if first_association is not None:
            signif = first_association[first_association.PValue < threshold]
            signif.to_csv(
                phen_dir / (phen_name + ".signif_snps.txt"), index=False, sep="\t"
            )
            logger.info("Significant SNPs: %d", len(signif))

    # Clean up FaST-LMM temp files
    runs_dir = phen_dir / "runs"
    if runs_dir.exists():
        shutil.rmtree(runs_dir)
    for work_dir in glob.glob(str(phen_dir / ".work*")):
        shutil.rmtree(work_dir)

    logger.info("=== Done: %s ===", phen_name)

    # Remove handlers so next phenotype gets a fresh log file
    logger.removeHandler(file_handler)
    file_handler.close()


def main() -> None:
    args = parse_args()

    # Suppress FaST-LMM's own verbose logging
    logging.getLogger("fastlmm").setLevel(logging.ERROR)

    kinship = args.kinship
    if kinship and not kinship.endswith(".bed"):
        kinship += ".bed"

    for pheno in args.phenotypes:
        t0 = time.time()
        run_gwas(
            bed_file=args.genotype,
            kinship_file=kinship,
            pheno_file=pheno,
            covar_file=args.covariance,
            outdir=args.outdir,
            nperm=args.nb_permutations,
            threads=args.threads,
            double_id=args.double_id,
            uncompressed=args.uncompressed,
        )
        logger.info("Total wall time: %.1f seconds", time.time() - t0)


if __name__ == "__main__":
    main()
