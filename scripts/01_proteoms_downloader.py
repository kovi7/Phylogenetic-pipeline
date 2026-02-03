import os
import sys
import zipfile
import subprocess
import argparse
import logging
from typing import List, Iterable

from tqdm import tqdm

# ----Config----
logger = logging.getLogger("proteomes_downloader")
OUTPUT_BASE_DIR = "data"
ARCHIVE_DIR = os.path.join(OUTPUT_BASE_DIR, "proteome_archives")
FASTA_DIR = os.path.join(OUTPUT_BASE_DIR, "proteome_fasta")

# ---- Helpers ----
def ensure_dirs_and_check_content(dir_paths: Iterable[str]) -> None:
    """Ensure directories exist and exit if any directory already contains files.

    This safety check avoids mixing outputs between runs.
    """
    for path in dir_paths:
        if not path:
            continue

        if os.path.isdir(path) and os.listdir(path):
            logger.error("Directory '%s' is not empty. Aborting to avoid overwrite.", path)
            sys.exit(1)

        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logger.info("Created directory: %s", path)


def load_accession_lines(file: str) -> List[str]:
    """Load accession lines from a text file.

    Empty lines and comment lines (starting with '#') are ignored.
    """
    with open(file, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


# ----Functions----
def download_proteomes_from_ncbi(file: str, output_dir: str) -> None:
    """Download proteome ZIP archives from NCBI using the Datasets CLI.

    Expected input format per line:
        GCF_...; optional_label
    """
    lines = load_accession_lines(file)
    logger.info("Loaded %d accessions from %s", len(lines), file)

    for line in tqdm(lines, desc="Downloading proteomes"):
        acc_id, sep, label = line.partition(";")
        acc_id = acc_id.strip()
        label = label.strip() if sep else "Missing organism label"

        output_file = os.path.join(output_dir, f"{acc_id}.zip")
        command = [
            "datasets",
            "download",
            "genome",
            "accession",
            acc_id,
            "--include",
            "protein",
            "--filename",
            output_file,
        ]

        try:
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            logger.error("Download failed for %s (%s): %s", acc_id, label, e)

    logger.info("Downloading finished")


def extract_protein_fastas(zip_dir: str, fasta_dir: str) -> None:
    """Extract protein FASTA (protein.faa) from each NCBI Datasets ZIP archive.

    Output files are written as:
        <accession>.faa

    This keeps genome_id consistent across the whole pipeline.
    """
    zip_files = [f for f in os.listdir(zip_dir) if f.endswith(".zip")]
    logger.info("ZIP archives to process: %d", len(zip_files))

    for filename in tqdm(zip_files, desc="Extracting protein.faa"):
        zip_path = os.path.join(zip_dir, filename)
        try:
            with zipfile.ZipFile(zip_path, "r") as z:
                protein_member = None
                for member in z.namelist():
                    if member.endswith("/protein.faa"):
                        protein_member = member
                        break

                if not protein_member:
                    logger.warning("protein.faa missing in %s", filename)
                    continue

                acc_id = os.path.splitext(filename)[0]
                out_faa = os.path.join(fasta_dir, acc_id + ".faa")

                with z.open(protein_member) as source, open(out_faa, "wb") as target:
                    target.write(source.read())

        except zipfile.BadZipFile:
            logger.error("Bad ZIP file: %s", filename)

    logger.info("Extraction finished")


# ----Main----
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download and extract proteomes (NCBI Datasets CLI).")
    parser.add_argument("--input_file", "-i", required=True, help="File with accession IDs and optional organism labels.")
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

    if not os.path.exists(args.input_file):
        logger.error("Input file does not exist: %s", args.input_file)
        sys.exit(1)

    ensure_dirs_and_check_content([ARCHIVE_DIR, FASTA_DIR])

    logger.info("Starting download")
    download_proteomes_from_ncbi(args.input_file, ARCHIVE_DIR)

    logger.info("Extracting protein FASTAs")
    extract_protein_fastas(ARCHIVE_DIR, FASTA_DIR)

    logger.info("Done")
