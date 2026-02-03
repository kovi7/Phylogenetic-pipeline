import os
import sys
import argparse
import logging
import csv
import re
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable, DefaultDict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from tqdm import tqdm

# ---- Config ----
logger = logging.getLogger("families_builder")

OUTPUT_BASE_DIR = "data"
MAPS_DIR = os.path.join(OUTPUT_BASE_DIR, "maps")
MERGED_FASTA_DIR = os.path.join(OUTPUT_BASE_DIR, "merged_fasta")

DEFAULT_CLUSTER_TSV = os.path.join("clustering", "clu_cluster.tsv")
DEFAULT_P2G = os.path.join(MAPS_DIR, "protein2genome.tsv")
DEFAULT_COMBINED = os.path.join(MERGED_FASTA_DIR, "merged_proteins.faa")

DEFAULT_OUT_DIR = "families"

# How to pick 1 protein per genome when building ortholog families
# Use strict mode to require exactly one sequence per genome
ORTHO_PICK_MODE = "strict"  # closest | longest | strict

# Used only by ORTHO_PICK_MODE="closest"
MAFFT_BIN = "mafft"
MAFFT_THREADS = 1

# Globals for workers
worker_seq_index = None
worker_p2g = None

# ---- Helpers ----
def ensure_dirs_and_check_content(dir_paths: Iterable[str], force: bool = False) -> None:
    """Ensure directories exist. Abort if a directory is non-empty (unless --force)."""
    for path in dir_paths:
        if not path: continue
        if os.path.isdir(path) and os.listdir(path) and not force:
            logger.error("Directory '%s' is not empty. Aborting to avoid overwrite.", path)
            sys.exit(1)
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logger.info("Created directory: %s", path)


def sanitize_label(label: str) -> str:
    """Convert a label into a safe FASTA header."""
    label = label.strip()
    label = re.sub(r"\s+", "_", label)
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
    return label

def load_map_tsv(p2g_tsv: str) -> Dict[str, Tuple[str, str]]:
    """Load mapping protein_id -> (genome_id, genome_label) from TSV."""
    mapping: Dict[str, Tuple[str, str]] = {}
    with open(p2g_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row["protein_id"].strip()
            gid = row["genome_id"].strip()
            glab = row["genome_label"].strip()
            if pid: mapping[pid] = (gid, glab)
    return mapping

def parse_mmseqs_clusters(cluster_tsv: str) -> Dict[str, List[str]]:
    """Parse MMseqs2 cluster TSV (rep -> member) into dict cluster_id -> members."""
    clusters: DefaultDict[str, List[str]] = defaultdict(list)
    with open(cluster_tsv, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            rep, member = (line.split("\t") + [""])[:2]
            if not rep or not member:
                continue
            clusters[rep].append(member)

    for rep in list(clusters.keys()):
        if rep not in clusters[rep]:
            clusters[rep].append(rep)

    logger.info("Parsed clusters: %d", len(clusters))
    return dict(clusters)


def pick_longest_protein(seq_index, protein_ids: List[str]) -> str:
    """Pick the longest sequence among candidate protein IDs."""
    best_id = protein_ids[0]
    best_len = -1
    for pid in protein_ids:
        if pid not in seq_index:
            continue
        L = len(seq_index[pid].seq)
        if L > best_len:
            best_len = L
            best_id = pid
    return best_id


def run_mafft_alignment(seq_records: List[SeqIO.SeqRecord], mafft_bin: str, threads: int):
    """Run MAFFT on provided SeqRecords and return an alignment, or None."""
    if not seq_records or len(seq_records) < 3:
        return None

    with tempfile.TemporaryDirectory(prefix="mafft_tmp_") as tmpdir:
        in_faa = os.path.join(tmpdir, "in.faa")
        out_faa = os.path.join(tmpdir, "out.faa")

        with open(in_faa, "w", encoding="utf-8") as f:
            SeqIO.write(seq_records, f, "fasta")

        cmd = [
            mafft_bin,
            "--quiet",
            "--auto",
            "--thread", str(max(1, int(threads))),
            in_faa,
        ]

        try:
            with open(out_faa, "w", encoding="utf-8") as out:
                subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.DEVNULL, text=True)
        except FileNotFoundError:
            logger.error("MAFFT not found in PATH (expected '%s').", mafft_bin)
            return None
        except subprocess.CalledProcessError:
            logger.warning("MAFFT failed for one cluster, falling back to longest-seq pruning.")
            return None

        try:
            return AlignIO.read(out_faa, "fasta")
        except Exception:
            logger.warning("Failed to parse MAFFT output, falling back to longest-seq pruning.")
            return None


def pick_medoid_per_genome(members: List[str], p2g: Dict[str, Tuple[str, str]], seq_index, mafft_bin: str, mafft_threads: int) -> Optional[Dict[str, str]]:
    """Choose 1 protein per genome, by closest-to-others criterion."""
    pid2gid: Dict[str, str] = {}
    seq_records: List[SeqIO.SeqRecord] = []

    for pid in members:
        if pid not in p2g or pid not in seq_index:
            continue
        gid, _ = p2g[pid]
        pid2gid[pid] = gid
        rec = seq_index[pid]
        rec.id = pid
        rec.name = pid
        rec.description = ""
        seq_records.append(rec)

    if len(seq_records) < 3:
        return None

    aln = run_mafft_alignment(seq_records, mafft_bin=mafft_bin, threads=mafft_threads)
    if aln is None:
        return None

    calc = DistanceCalculator("identity")
    dm = calc.get_distance(aln)
    names = list(dm.names)

    by_genome: DefaultDict[str, List[str]] = defaultdict(list)
    for pid, gid in pid2gid.items():
        by_genome[gid].append(pid)

    chosen: Dict[str, str] = {}
    for gid, pids in by_genome.items():
        if len(pids) == 1:
            chosen[gid] = pids[0]
            continue

        other = [n for n in names if pid2gid.get(n) != gid]
        if not other:
            chosen[gid] = pick_longest_protein(seq_index, pids)
            continue

        best_pid = None
        best_mean = None

        for pid in pids:
            if pid not in names:
                continue
            dists = [dm[pid, o] for o in other if o in names]
            if not dists:
                continue
            mean_d = sum(dists) / float(len(dists))

            if best_mean is None or mean_d < best_mean:
                best_mean = mean_d
                best_pid = pid

        chosen[gid] = best_pid if best_pid is not None else pick_longest_protein(seq_index, pids)

    return chosen


def build_ortholog_family(members: List[str], p2g: Dict[str, Tuple[str, str]], seq_index, mode: str) -> Optional[Dict[str, Tuple[str, str]]]:
    """Convert a MMseqs cluster into a 1-per-genome family."""
    by_genome: DefaultDict[str, List[str]] = defaultdict(list)

    for pid in members:
        if pid not in p2g:
            continue
        gid, _ = p2g[pid]
        by_genome[gid].append(pid)

    if not by_genome:
        return None

    chosen: Dict[str, Tuple[str, str]] = {}

    if mode == "strict":
        for gid, pid_list in by_genome.items():
            if len(pid_list) != 1:
                return None
            selected = pid_list[0]
            _, glab = p2g[selected]
            chosen[gid] = (glab, selected)
        return chosen

    if mode == "closest":
        picked = pick_medoid_per_genome(members, p2g, seq_index, mafft_bin=MAFFT_BIN, mafft_threads=MAFFT_THREADS)
        if picked is not None:
            for gid, selected in picked.items():
                _, glab = p2g[selected]
                chosen[gid] = (glab, selected)
            return chosen
        # fallback
        mode = "longest"

    for gid, pid_list in by_genome.items():
        selected = pid_list[0] if len(pid_list) == 1 else pick_longest_protein(seq_index, pid_list)
        _, glab = p2g[selected]
        chosen[gid] = (glab, selected)

    return chosen

# ---- Worker Function ----
def init_worker(p2g_path, fasta_path):
    """Initialize global variables in worker process."""
    global worker_p2g, worker_seq_index
    worker_p2g = load_map_tsv(p2g_path)
    worker_seq_index = SeqIO.index(fasta_path, "fasta")


def process_single_cluster(args):
    """
    args: (cluster_id, members, family_index, min_genomes, req_all, n_all, tag, out_ortho, out_para)
    Returns: stats_row (list)
    """
    cluster_id, members, i, min_genomes, req_all, n_all, tag, out_ortho, out_para = args
    
    global worker_p2g, worker_seq_index
    p2g = worker_p2g
    seq_index = worker_seq_index

    by_genome = defaultdict(int)
    unmapped = 0
    for pid in members:
        if pid not in p2g:
            unmapped += 1
            continue
        gid, _ = p2g[pid]
        by_genome[gid] += 1

    n_genomes = len(by_genome)
    n_multi = sum(1 for c in by_genome.values() if c > 1)
    
    note_parts = []
    if unmapped: note_parts.append(f"unmapped={unmapped}")

    # separate keep logic for paralogs vs orthologs
    keep_for_paralogs = n_genomes >= min_genomes
    keep_for_orthologs = keep_for_paralogs
    if req_all and n_genomes != n_all:
        keep_for_orthologs = False
        note_parts.append("missing_genomes")

    kept_o = "NO"; kept_p = "NO"

    # create family if either orthologs or paralogs should be written
    if keep_for_paralogs or keep_for_orthologs:
        family_name = f"F{i:06d}"

        # 1. ORTHOLOGS (apply require_all_genomes here only)
        if tag in {"orthologs", "both"} and keep_for_orthologs:
            try:
                chosen = build_ortholog_family(members, p2g, seq_index, mode=ORTHO_PICK_MODE)
                if chosen:
                    out_faa = os.path.join(out_ortho, family_name + ".faa")
                    seen = set()
                    with open(out_faa, "w", encoding="utf-8") as out:
                        for gid, (glab, pid) in sorted(chosen.items(), key=lambda x: x[0]):
                            if pid not in seq_index: continue
                            header = sanitize_label(glab) if glab else gid
                            if header in seen: header = f"{header}__{gid}"
                            seen.add(header)
                            out.write(f">{header}\n{str(seq_index[pid].seq)}\n")
                    kept_o = "YES"
                else:
                    note_parts.append("cannot_make_1per_genome")
            except Exception as e:
                note_parts.append(f"err_ortho:{str(e)}")

        # 2. PARALOGS 
        if tag in {"paralogs", "both"} and keep_for_paralogs:
            try:
                out_faa = os.path.join(out_para, family_name + ".faa")
                with open(out_faa, "w", encoding="utf-8") as out:
                    for pid in members:
                        if pid in p2g and pid in seq_index:
                            gid, glab = p2g[pid]
                            header = sanitize_label(glab) if glab else gid
                            out.write(f">{header}__{pid}\n{str(seq_index[pid].seq)}\n")
                kept_p = "YES"
            except Exception as e:
                note_parts.append(f"err_para:{str(e)}")

    return [cluster_id, len(members), n_genomes, n_multi, kept_o, kept_p, ";".join(note_parts)]

# ---- Main ----
def main() -> int:
    parser = argparse.ArgumentParser(description="Step III: prepare gene families from MMseqs2 clusters.")
    parser.add_argument("--tag", choices=["orthologs", "paralogs", "both"], default="orthologs",
                        help="Which family set to write.")
    parser.add_argument("--clusters_tsv", "-c", default=DEFAULT_CLUSTER_TSV,
                        help="MMseqs2 *_cluster.tsv file (rep -> member).")
    parser.add_argument("--out_dir", "-o", default=DEFAULT_OUT_DIR, help="Output base directory (default: families).")
    parser.add_argument("--min_genomes", type=int, default=None, help="Keep families with at least this many genomes. If not provided and writing orthologs, defaults to all genomes (strict 1-1 like genomics_project).")
    parser.add_argument("--require_all_genomes", action="store_true",
                        help="Keep only families containing all genomes (useful for consensus).")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes.")
    parser.add_argument("--force", action="store_true", help="Allow writing into non-empty output directories.")
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

    if not os.path.exists(args.clusters_tsv):
        logger.error("Clusters TSV does not exist: %s", args.clusters_tsv)
        return 2
    if not os.path.exists(DEFAULT_P2G):
        logger.error("protein2genome TSV does not exist: %s", DEFAULT_P2G)
        return 2
    if not os.path.exists(DEFAULT_COMBINED):
        logger.error("Merged FASTA does not exist: %s", DEFAULT_COMBINED)
        return 2

    out_ortho = os.path.join(args.out_dir, "orthologs")
    out_para = os.path.join(args.out_dir, "paralogs")
    out_dirs = []
    if args.tag in {"orthologs", "both"}: out_dirs.append(out_ortho)
    if args.tag in {"paralogs", "both"}: out_dirs.append(out_para)

    ensure_dirs_and_check_content(out_dirs, force=args.force)
    os.makedirs(args.out_dir, exist_ok=True)

    # Pre-scan mapping for genome count
    logger.info("Loading protein-to-genome mapping")
    temp_p2g = load_map_tsv(DEFAULT_P2G)
    genomes_all = {gid for (gid, _) in temp_p2g.values()}
    n_all = len(genomes_all)
    logger.info("Detected %d genomes in mapping", n_all)

    del temp_p2g

    logger.info("Indexing merged FASTA")
    seq_index = SeqIO.index(DEFAULT_COMBINED, "fasta")

    logger.info("Reading MMseqs clusters")
    cluster_map = parse_mmseqs_clusters(args.clusters_tsv)

    if args.min_genomes is None:
        args.min_genomes = n_all if args.tag == "orthologs" else 25
        logger.info("--min_genomes not provided; defaulting to %d (tag=%s)", args.min_genomes, args.tag)

    tasks = []
    for i, (cid, mem) in enumerate(cluster_map.items(), start=1):
        tasks.append((cid, mem, i, args.min_genomes, args.require_all_genomes, n_all, args.tag, out_ortho, out_para))

    stats_path = os.path.join(args.out_dir, "families_stats.tsv")
    kept_ortho_count = 0
    kept_para_count = 0
    logger.info(f"Processing {len(tasks)} clusters with {args.threads} threads...")

    with open(stats_path, "w", encoding="utf-8", newline="") as sf:
        w = csv.writer(sf, delimiter="\t")
        w.writerow(["cluster_id", "n_members", "n_genomes", "n_multi_copy_genomes", "kept_ortholog", "kept_paralog", "note"])
        
        with ProcessPoolExecutor(max_workers=args.threads, initializer=init_worker, initargs=(DEFAULT_P2G, DEFAULT_COMBINED)) as executor:
            futures = [executor.submit(process_single_cluster, t) for t in tasks]
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Building families"):
                res = future.result()
                w.writerow(res)
                if res[4] == "YES": kept_ortho_count += 1
                if res[5] == "YES": kept_para_count += 1

    logger.info("Stats saved: %s", stats_path)
    logger.info("Ortholog families: %d | Paralog families: %d", kept_ortho_count, kept_para_count)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
