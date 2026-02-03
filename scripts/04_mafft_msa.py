import os
import sys
import argparse
import logging
import subprocess
from concurrent.futures import ProcessPoolExecutor
from typing import Iterable, List, Tuple

from tqdm import tqdm

logger = logging.getLogger("mafft_msa")

FASTA_SUFFIXES = (".faa", ".fa", ".fasta", ".fsa", ".fas", ".aln")

DEFAULT_FAMILIES_DIR = "families"
DEFAULT_OUT_BASE = "allignment"

DEFAULT_WORKERS = 4
DEFAULT_MAFFT_THREADS = 1


def ensure_dirs_and_check_content(dir_paths: Iterable[str], force: bool = False) -> None:
    """Create missing directories. By default, abort if any target dir is non empty."""
    for path in dir_paths:
        if not path:
            continue

        if os.path.isdir(path) and os.listdir(path) and not force:
            logger.error("Directory '%s' is not empty. Aborting to avoid overwrite.", path)
            sys.exit(1)

        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logger.info("Created directory: %s", path)


def list_fasta_files(input_dir: str) -> List[str]:
    """Return sorted FASTA like filenames from input_dir."""
    if not os.path.isdir(input_dir):
        return []
    files = [f for f in os.listdir(input_dir) if f.lower().endswith(FASTA_SUFFIXES)]
    files.sort()
    return files


def build_mafft_command(input_path: str, threads: int) -> List[str]:
    cmd = ["mafft", "--auto", "--quiet", "--thread", str(max(1, int(threads))), input_path]
    return cmd


def align_one_file(filename: str, input_dir: str, output_dir: str, mafft_threads: int) -> Tuple[str, str]:
    in_path = os.path.join(input_dir, filename)
    out_name = os.path.splitext(filename)[0] + ".aln.faa"
    out_path = os.path.join(output_dir, out_name)

    if not os.path.exists(in_path):
        return (filename, "Missing input file")

    cmd = build_mafft_command(in_path, threads=mafft_threads)

    try:
        with open(out_path, "w", encoding="utf-8") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.DEVNULL, text=True)
        return (filename, "OK")
    except FileNotFoundError:
        return (filename, "MAFFT not found in PATH")
    except subprocess.CalledProcessError as e:
        return (filename, f"MAFFT failed (code {e.returncode})")
    except Exception as e:
        return (filename, f"Unexpected error: {e}")


def run_mafft_batch(input_dir: str, output_dir: str, workers: int, mafft_threads: int) -> None:
    fasta_files = list_fasta_files(input_dir)
    if not fasta_files:
        logger.warning("No FASTA files found in: %s", input_dir)
        return

    os.makedirs(output_dir, exist_ok=True)

    logger.info("Input: %s", input_dir)
    logger.info("Output: %s", output_dir)
    logger.info("Files to align: %d", len(fasta_files))
    logger.info("Workers: %d | MAFFT threads per job: %d", workers, mafft_threads)

    ok = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(tqdm(
            ex.map(
                align_one_file,
                fasta_files,
                [input_dir] * len(fasta_files),
                [output_dir] * len(fasta_files),
                [mafft_threads] * len(fasta_files),
            ),
            total=len(fasta_files),
            desc=f"MAFFT ({os.path.basename(output_dir)})",
            unit="file",
        ))

    for fn, status in results:
        if status == "OK":
            ok += 1
        else:
            failed += 1
            logger.warning("Alignment issue for %s: %s", fn, status)

    logger.info("MAFFT finished. OK: %d | Failed: %d", ok, failed)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run MAFFT on families (Step IV).")
    parser.add_argument("--tag", choices=["orthologs", "paralogs"], required=True,
                        help="Which family set to align.")

    parser.add_argument("--input_dir", default=None,
                        help="Override input directory. Default: families/<tag>")
    parser.add_argument("--output_dir", default=None,
                        help="Override output directory. Default: alignment/<tag>")

    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--mafft_threads", type=int, default=DEFAULT_MAFFT_THREADS)
    parser.add_argument("--force", action="store_true", help="Allow writing into non empty output dir")
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

    in_dir = args.input_dir or os.path.join(DEFAULT_FAMILIES_DIR, args.tag)
    out_dir = args.output_dir or os.path.join(DEFAULT_OUT_BASE, args.tag)

    if not os.path.isdir(in_dir):
        logger.error("Input directory does not exist: %s", in_dir)
        return 2

    ensure_dirs_and_check_content([out_dir], force=args.force)

    run_mafft_batch(
        input_dir=in_dir,
        output_dir=out_dir,
        workers=max(1, int(args.workers)),
        mafft_threads=max(1, int(args.mafft_threads)),
    )

    logger.info("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
