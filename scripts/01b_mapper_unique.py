import os
import sys
import argparse
import logging
import csv
from typing import Dict, List, Iterable, Tuple

from Bio import SeqIO
from tqdm import tqdm

# ----Config----
logger = logging.getLogger("proteome_mapping")
OUTPUT_BASE_DIR = "data"
FASTA_DIR = os.path.join(OUTPUT_BASE_DIR, "proteome_fasta")
MAPS_DIR = os.path.join(OUTPUT_BASE_DIR, "maps")
MERGED_FASTA_DIR = os.path.join(OUTPUT_BASE_DIR, "merged_fasta")


# ---- Helpers ----
def ensure_dirs_and_check_content(dir_paths: Iterable[str], force: bool = False) -> None:
    """
    Create missing directories. By default, stop if a directory already has files.
    This helps avoid accidentally mixing outputs from different runs.
    """
    for path in dir_paths:
        if not path:
            continue

        if os.path.isdir(path) and os.listdir(path) and not force:
            logger.error("Directory '%s' is not empty. Aborting to avoid overwrite.", path)
            sys.exit(1)

        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logger.info("Created directory: %s", path)


def load_accession_lines(file: str) -> List[str]:
    """
    Read accession lines from a text file.

    Notes:
    - Blank lines are skipped
    - Lines starting with '#' are treated as comments
    """
    with open(file, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]


def unique_protein_id(genome_id: str, protein_id: str) -> str:
    """
    Build a globally unique protein identifier.

    RefSeq protein IDs (e.g., WP_...) can be shared across many genomes, so we
    prefix them with the assembly/genome ID (derived from the *.faa filename).
    """
    return f"{genome_id}__{protein_id}"


def ensure_file_not_exists(path: str, force: bool = False) -> None:
    """
    Abort if an output file already exists, unless --force is used.
    """
    if os.path.exists(path) and not force:
        logger.error("Output file already exists: %s (use --force to overwrite)", path)
        sys.exit(1)


# ---- Functions ----
def build_accession_label_map_tsv(file: str, out_dir: str, out_name: str = "names_map.tsv", force:bool=False) -> Dict[str, str]:
    """
    Build accession -> label mapping and save it as TSV.

    Output columns:
        accession <tab> label
    """
    logger.info("Creating accession-to-label mapping (TSV)")
    mapping: Dict[str, str] = {}

    out_tsv = os.path.join(out_dir, out_name)
    ensure_file_not_exists(out_tsv, force=force)

    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "label"])

        for line in load_accession_lines(file):
            acc, sep, label = line.partition(";")
            acc = acc.strip()
            label = label.strip() if sep else "Missing label"
            mapping[acc] = label
            w.writerow([acc, label])

    logger.info("Saved accession label map: %s", out_tsv)
    return mapping


def write_mapping_and_merged_fasta(protein_files_dir: str, genome_id_to_name: Dict[str, str], out_dir: str, out_map_name: str = "protein2genome.tsv", merged_fasta_path: str = os.path.join(MERGED_FASTA_DIR, "merged_proteins.faa"), force: bool = False,) -> Tuple[int, int]:
    """
    Write two consistent outputs in one pass over the proteomes:
    1) protein2genome.tsv  (unique protein_id -> genome_id, genome_label)
    2) merged_proteins.faa (FASTA with the same unique protein_id)

    The key point:
    - protein_id is written as: genome_id__original_protein_id

    Returns:
        (n_proteins_written, n_fasta_files)
    """
    logger.info("Building protein-to-genome mapping and merged FASTA with unique IDs")

    fasta_files = [f for f in os.listdir(protein_files_dir) if f.endswith(".faa")]
    if not fasta_files:
        logger.error("No .faa files found in: %s", protein_files_dir)
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(merged_fasta_path), exist_ok=True)

    out_tsv = os.path.join(out_dir, out_map_name)
    ensure_file_not_exists(out_tsv, force=force)
    ensure_file_not_exists(merged_fasta_path, force=force)

    n_proteins = 0

    with open(out_tsv, "w", encoding="utf-8", newline="") as map_f, \
            open(merged_fasta_path, "w", encoding="utf-8") as merged_f:

        w = csv.writer(map_f, delimiter="\t")
        # keep raw protein ID for debugging / traceability
        w.writerow(["protein_id", "genome_id", "genome_label", "protein_id_raw"])

        for file in tqdm(sorted(fasta_files), desc="Parsing FASTA"):
            genome_id = os.path.splitext(file)[0]
            genome_label = genome_id_to_name.get(genome_id, "Missing genome label")
            file_path = os.path.join(protein_files_dir, file)

            for rec in SeqIO.parse(file_path, "fasta"):
                raw_id = rec.id
                pid = unique_protein_id(genome_id, raw_id)

                # mapping row
                w.writerow([pid, genome_id, genome_label, raw_id])

                # merged FASTA (unique headers -> no duplicate keys later)
                merged_f.write(f">{pid}\n{str(rec.seq)}\n")

                n_proteins += 1

    logger.info("Saved protein2genome TSV: %s", out_tsv)
    logger.info("Saved merged FASTA: %s", merged_fasta_path)
    logger.info("Proteins written: %d", n_proteins)

    return n_proteins, len(fasta_files)


# ---- Main ----
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build mapping files from extracted proteomes (TSV) and create a merged FASTA with unique protein IDs."
    )
    parser.add_argument("--input_file", "-i", required=True, help="File with accession IDs and organism labels.")
    parser.add_argument("--fasta_dir", default=FASTA_DIR, help="Directory with extracted *.faa proteomes.")
    parser.add_argument("--maps_dir", default=MAPS_DIR, help="Directory where mapping TSV files will be written.")
    parser.add_argument(
        "--merged_fasta",
        default=os.path.join(MERGED_FASTA_DIR, "merged_proteins.faa"),
        help="Output path for merged FASTA with unique protein IDs."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting existing output files/directories."
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level."
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    if not os.path.exists(args.input_file):
        logger.error("Input file does not exist: %s", args.input_file)
        sys.exit(1)

    if not os.path.isdir(args.fasta_dir) or not os.listdir(args.fasta_dir):
        logger.error("FASTA directory is missing or empty: %s", args.fasta_dir)
        sys.exit(1)

    ensure_dirs_and_check_content([args.maps_dir], force=args.force)

    logger.info("Building names map")
    id2name = build_accession_label_map_tsv(args.input_file, args.maps_dir, force=args.force)

    logger.info("Building protein2genome + merged FASTA")
    write_mapping_and_merged_fasta(
        protein_files_dir=args.fasta_dir,
        genome_id_to_name=id2name,
        out_dir=args.maps_dir,
        out_map_name="protein2genome.tsv",
        merged_fasta_path=args.merged_fasta,
        force=args.force,
    )

    logger.info("Mapping finished")
