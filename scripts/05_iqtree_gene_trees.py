import os
import sys
import argparse
import logging
import subprocess
import re
from concurrent.futures import ProcessPoolExecutor
from typing import Iterable, List, Optional, Tuple

from tqdm import tqdm

logger = logging.getLogger("gene_trees")

FASTA_SUFFIXES = (".faa", ".fa", ".fasta", ".fas", ".fsa", ".aln")

DEFAULT_MSA_BASE = "allignment"
DEFAULT_OUT_BASE = os.path.join("trees", "gene_trees")

DEFAULT_WORKERS = 4
DEFAULT_THREADS_PER_JOB = 1
DEFAULT_MODEL = "MFP"


def ensure_dirs_and_check_content(dir_paths: Iterable[str], force: bool = False) -> None:
    """Create missing directories. By default, abort if any directory is non empty."""
    for path in dir_paths:
        if not path:
            continue
        if os.path.isdir(path) and os.listdir(path) and not force:
            logger.error("Directory '%s' is not empty. Aborting to avoid overwrite.", path)
            sys.exit(1)
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logger.info("Created directory: %s", path)


def list_msa_files(msa_dir: str) -> List[str]:
    if not os.path.isdir(msa_dir):
        return []
    files = [f for f in os.listdir(msa_dir) if f.lower().endswith(FASTA_SUFFIXES)]
    files.sort()
    return files


def shutil_which(cmd: str) -> Optional[str]:
    """Tiny 'which' replacement."""
    paths = os.environ.get("PATH", "").split(os.pathsep)
    for p in paths:
        full = os.path.join(p, cmd)
        if os.path.isfile(full) and os.access(full, os.X_OK):
            return full
    return None


def detect_iqtree_binary(user_bin: Optional[str] = None) -> str:
    if user_bin:
        return user_bin

    for cand in ("iqtree2", "iqtree"):
        if shutil_which(cand):
            return cand

    logger.error("IQ-TREE not found in PATH (expected iqtree2 or iqtree)")
    sys.exit(2)


def average_bootstrap_support(newick: str) -> float:
    """Mean support from node labels in Newick.

    This is a heuristic. It looks for integers right after ')' and before ':'.
    """
    vals = re.findall(r"\)([0-9]+(?:\.[0-9]+)?):", newick)
    if not vals:
        return 0.0
    nums = [float(v) for v in vals]
    return sum(nums) / len(nums)


def run_iqtree_on_alignment(msa_file: str, msa_dir: str, out_dir: str, iqtree_bin: str, threads: int, model: str, ufboot: int) -> Tuple[str, str]:
    """Run IQ-TREE for one alignment."""
    msa_path = os.path.join(msa_dir, msa_file)
    if not os.path.exists(msa_path):
        return (msa_file, "Missing input")

    stem = os.path.splitext(msa_file)[0]
    out_prefix = os.path.join(out_dir, stem)

    expected = [out_prefix + ".treefile"]
    if ufboot > 0:
        expected.append(out_prefix + ".contree")

    if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in expected):
        return (msa_file, "OK (Skipped, exists)")
    
    cmd = [
        iqtree_bin,
        "-s", msa_path,
        "-T", str(max(1, int(threads))),
        "-pre", out_prefix,
        "-m", model,
        "-quiet",
    ]

    if ufboot > 0:
        if ufboot < 1000:
            raise ValueError("UFBoot (-B) should have at least 1000 replicates.")
        cmd += ["-B", str(ufboot)]


    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return (msa_file, "OK")
    except subprocess.CalledProcessError as e:
        return (msa_file, f"IQ-TREE failed (code {e.returncode}), see {out_prefix}.log")
    except FileNotFoundError:
        return (msa_file, "IQ-TREE not found")
    except Exception as e:
        return (msa_file, f"Unexpected error: {e}")


def merge_tree_outputs(out_dir: str, output_name: str, ext: str, filter_by_support: bool, support_threshold: Optional[float] = None) -> Tuple[int, int, str]:
    """MMerge one tree per file into a single file, one Newick per line.
    ext:
      ".treefile" -> srandard ML trees
      ".contree"  -> bootstrap (supports)
    Returns:
      (written, skipped, used_extension)
    """
    files = [f for f in os.listdir(out_dir) if f.endswith(ext)]
    files.sort()
    if not files:
        return (0, 0)
    written = 0
    skipped = 0
    out_path = os.path.join(out_dir, output_name)

    with open(out_path, "w", encoding="utf-8") as out:
        for fn in files:
            p = os.path.join(out_dir, fn)
            with open(p, "r", encoding="utf-8") as f:
                tree = f.readline().strip()
            if not tree:
                skipped += 1
                continue

            if filter_by_support and support_threshold is not None:
                mean_sup = average_bootstrap_support(tree)
                if mean_sup < float(support_threshold):
                    skipped += 1
                    continue

            out.write(tree + "\n")
            written += 1
    return (written, skipped)


def run_batch(msa_dir: str, out_dir: str, iqtree_bin: str, workers: int, threads_per_job: int, model: str, ufboot: int) -> None:
    msa_files = list_msa_files(msa_dir)
    if not msa_files:
        logger.warning("No MSA files found in: %s", msa_dir)
        return

    os.makedirs(out_dir, exist_ok=True)

    logger.info("Input MSA dir: %s", msa_dir)
    logger.info("Output dir: %s", out_dir)
    logger.info("MSA files: %d", len(msa_files))
    logger.info("Workers: %d | Threads per IQ-TREE job: %d", workers, threads_per_job)
    logger.info("Model: %s | UFBoot: %d", model, ufboot)

    ok = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(tqdm(
            ex.map(
                run_iqtree_on_alignment,
                msa_files,
                [msa_dir] * len(msa_files),
                [out_dir] * len(msa_files),
                [iqtree_bin] * len(msa_files),
                [threads_per_job] * len(msa_files),
                [model] * len(msa_files),
                [ufboot] * len(msa_files),
            ),
            total=len(msa_files),
            desc=f"IQ-TREE ({os.path.basename(out_dir)})",
            unit="file",
        ))

    for fn, status in results:
        if status == "OK":
            ok += 1
        else:
            failed += 1
            logger.warning("Tree failed for %s: %s", fn, status)

    logger.info("Tree building done. OK: %d | Failed: %d", ok, failed)


def default_dirs_for_tag(tag: str) -> Tuple[str, str]:
    msa_dir = os.path.join(DEFAULT_MSA_BASE, tag)
    out_dir = os.path.join(DEFAULT_OUT_BASE, tag)
    return msa_dir, out_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute gene family trees from MSAs using IQ-TREE and merge them.")
    parser.add_argument("--tag", choices=["orthologs", "paralogs"], required=True,
                        help="Which dataset to process.")
    parser.add_argument("--msa_dir", default=None,
                        help="Override input MSA dir. Default: allignment/<tag>.")
    parser.add_argument("--out_dir", default=None,
                        help="Override output dir. Default: trees/gene_trees/<tag>.")
    parser.add_argument("--ufboot", type=int, default=0,
                    help="Ultrafast bootstrap replicates for IQ-TREE (-B). 0 disables.")
    parser.add_argument("--support_threshold", type=float, default=0.0,
                        help="Mean bootstrap support threshold for filtering.")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help="Parallel processes.")
    parser.add_argument("--threads_per_job", type=int, default=DEFAULT_THREADS_PER_JOB,
                        help="Threads per IQ-TREE process.")
    parser.add_argument("--force", action="store_true",
                        help="Allow writing into non empty output directories.")
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    msa_def, out_def = default_dirs_for_tag(args.tag)
    msa_dir = args.msa_dir or msa_def
    out_dir = args.out_dir or out_def

    iqtree_bin = detect_iqtree_binary()

    ensure_dirs_and_check_content([out_dir], force=args.force)

    run_batch(
        msa_dir=msa_dir,
        out_dir=out_dir,
        iqtree_bin=iqtree_bin,
        workers=int(args.workers),
        threads_per_job=int(args.threads_per_job),
        model=DEFAULT_MODEL,
        ufboot=int(args.ufboot),
    )

    # Merge trees into one file
    w1, s1 = merge_tree_outputs(
        out_dir=out_dir,
        output_name="all_trees.nwk",
        ext=".treefile",
        filter_by_support=False,
    )
    logger.info("Merged base trees. Written: %d | Skipped: %d", w1, s1)

    if args.ufboot > 0:
        w2, s2 = merge_tree_outputs(
            out_dir=out_dir,
            output_name="all_trees_bootstrap.nwk",
            ext=".contree",        
            filter_by_support=False
        )

        logger.info("Merged bootstrap (contree) trees. Written: %d | Skipped: %d", w2, s2)

        if args.support_threshold > 0:
            wf, sf = merge_tree_outputs(
                out_dir=out_dir,
                output_name="all_trees_filtered.nwk",
                ext=".contree",    
                filter_by_support=True,
                support_threshold=float(args.support_threshold) ,
            )
            logger.info("Merged filtered trees. Written: %d | Skipped: %d", wf, sf)
            

    elif args.support_threshold > 0:
        logger.warning("filtering ignored because ufboot=0 (no supports).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())