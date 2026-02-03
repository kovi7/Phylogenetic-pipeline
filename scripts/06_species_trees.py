import os
import argparse
import logging
import subprocess
import glob
import re
import shutil
from typing import Tuple, List, Set

logger = logging.getLogger("species_trees")

DEFAULT_GENE_TREES_BASE = os.path.join("trees", "gene_trees")
DEFAULT_OUT_DIR = os.path.join("trees", "species_trees")
DEFAULT_CONS_R = os.path.join("scripts", "06b_consensus_species_tree.R")


DEFAULT_STRIP_REGEX = r"__[A-Za-z0-9_.]+"

def ensure_dir(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)

def pick_gene_trees_file(gene_trees_dir: str, tag: str, filtered: bool = False) -> str:
    base = os.path.join(gene_trees_dir, tag)
    return os.path.join(base, "all_trees_filtered.nwk" if filtered else "all_trees.nwk")

def parse_fasturec_output(fu_tsv: str) -> Tuple[float, str, int]:
    best_cost = float("inf")
    best_tree = ""
    n = 0
    with open(fu_tsv, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"^\s*([0-9]+(?:\.[0-9]+)?)\s+(.+)$", line)
            if not m:
                continue
            cost = float(m.group(1))
            tree = m.group(2).strip()
            if not tree.endswith(";"):
                tree += ";"
            n += 1
            if cost < best_cost:
                best_cost = cost
                best_tree = tree
    if n == 0 or not best_tree:
        raise ValueError(f"No parseable candidates in: {fu_tsv}")
    return best_cost, best_tree, n

def run_consensus_r(rscript: str, consensus_r: str, trees_file: str, out_tree: str,
                    strict_taxa: bool = False) -> bool:
    if not os.path.exists(trees_file):
        logger.error("Missing input for consensus: %s", trees_file)
        return False

    cmd = [
        rscript, consensus_r,
        "--trees_file", trees_file,
        "--out_tree", out_tree,
        "--min_support", "50.0",
        "--strict_taxa", "TRUE" if strict_taxa else "FALSE",
        "--use_intersection", "FALSE",
    ]

    logger.info("Running consensus (R): %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logger.info("Consensus OK: %s", out_tree)
        return True
    except subprocess.CalledProcessError as e:
        logger.warning("Consensus R script failed (often taxa mismatch).")
        logger.warning("R stderr:\n%s", e.stderr)
        return False

def detect_exe(user_path_or_name: str, candidates: List[str]) -> str:
    if user_path_or_name:
        if os.path.isfile(user_path_or_name) and os.access(user_path_or_name, os.X_OK):
            return user_path_or_name
        w = shutil.which(user_path_or_name)
        if w:
            return w
    for c in candidates:
        w = shutil.which(c)
        if w:
            return w
    return ""

def clean_tree_labels(in_path: str, out_path: str, strip_regex: str) -> None:
    rgx = re.compile(strip_regex) if strip_regex else None
    with open(in_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        for line in fin:
            s = line.rstrip("\n")
            if rgx:
                s = rgx.sub("", s)
            fout.write(s + "\n")

def extract_tip_labels_from_newick(newick: str) -> Set[str]:
    s = newick.strip()
    if not s:
        return set()
    s = re.sub(r":[^,)\s;]+", "", s)         # branch lengths
    s = re.sub(r"\)([^,);]+)", ")", s)       # internal node labels
    toks = re.split(r"[(),;]", s)
    tips = set()
    for t in toks:
        t = t.strip().strip("'\"")
        if t:
            tips.add(t)
    return tips

def filter_trees_to_max_taxa_set(in_path: str, out_path: str) -> Tuple[int, int, int]:
    trees: List[str] = []
    tipsets: List[Set[str]] = []

    with open(in_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if not line.endswith(";"):
                line += ";"
            tips = extract_tip_labels_from_newick(line)
            if len(tips) < 4:
                continue
            trees.append(line)
            tipsets.append(tips)

    n_total = len(trees)
    if n_total == 0:
        raise ValueError(f"No usable trees in: {in_path}")

    max_taxa = max(len(s) for s in tipsets)
    ref = next(s for s in tipsets if len(s) == max_taxa)
    kept = [t for t, tips in zip(trees, tipsets) if tips == ref]

    with open(out_path, "w", encoding="utf-8") as out:
        for t in kept:
            out.write(t + "\n")

    return n_total, len(kept), max_taxa

def run_iqtree_greedy_consensus(trees_file: str, out_tree: str, work_dir: str, cores: int,
                               minsup: float = 0.0, strip_regex: str = DEFAULT_STRIP_REGEX) -> bool:
    """
    IQ-TREE consensus:
      -con + -minsup 0.0 -> extended majority-rule ("greedy" w praktyce)
    """
    if not os.path.exists(trees_file):
        logger.error("Missing input for IQ-TREE consensus: %s", trees_file)
        return False

    exe = detect_exe("", ["iqtree2", "iqtree"])
    if not exe:
        logger.warning("IQ-TREE not found in PATH (iqtree2/iqtree). Skipping greedy consensus.")
        return False

    ensure_dir(work_dir)
    ensure_dir(os.path.dirname(out_tree))

    cleaned = os.path.join(work_dir, "gene_trees.cleaned.nwk")
    strict = os.path.join(work_dir, "gene_trees.cleaned.strict.nwk")
    log_path = os.path.join(work_dir, "iqtree_consensus.log")

    logger.info("Preparing trees for IQ-TREE consensus (strip_regex=%s)...", strip_regex if strip_regex else "NONE")
    clean_tree_labels(trees_file, cleaned, strip_regex)

    try:
        n_total, n_kept, max_taxa = filter_trees_to_max_taxa_set(cleaned, strict)
        logger.info("Strict taxa filter: kept %d/%d trees (max_taxa=%d)", n_kept, n_total, max_taxa)
        if n_kept < 2:
            logger.warning("Too few trees after strict taxa filtering, skipping IQ-TREE consensus.")
            return False
        input_for_iqtree = strict
    except Exception as e:
        logger.warning("Taxa filtering failed (%s). Trying IQ-TREE on cleaned trees without filtering.", e)
        input_for_iqtree = cleaned

    prefix = "iqtree_consensus"
    cmd = [
        exe,
        "-t", os.path.abspath(input_for_iqtree),
        "-con",
        "-minsup", str(minsup),
        "--prefix", prefix,
        "-T", str(cores),
        "-scale", "100",
        "-prec", "1",
    ]

    logger.info("Running IQ-TREE consensus: %s", " ".join(cmd))
    try:
        with open(log_path, "w") as lf:
            subprocess.run(cmd, check=True, cwd=work_dir, stdout=lf, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.error("IQ-TREE consensus failed. Check log: %s", log_path)
        return False

    contree = os.path.join(work_dir, prefix + ".contree")
    if not os.path.exists(contree):
        candidates = glob.glob(os.path.join(work_dir, "*.contree"))
        if candidates:
            contree = max(candidates, key=os.path.getmtime)

    if not os.path.exists(contree):
        logger.error("IQ-TREE finished but no .contree found in: %s (log: %s)", work_dir, log_path)
        return False

    shutil.copyfile(contree, out_tree)
    logger.info("Greedy consensus saved: %s", out_tree)
    return True

def run_supertree_fasturec(fasturec_bin: str, trees_file: str, out_tree: str, work_dir: str,
                          cores: int, strip_regex: str = DEFAULT_STRIP_REGEX) -> None:
    if not os.path.exists(trees_file):
        logger.error("Missing input for Fasturec: %s", trees_file)
        return

    exe = detect_exe(fasturec_bin, ["fasturec"])
    if not exe:
        logger.error("Fasturec executable not found: %s", fasturec_bin)
        return

    ensure_dir(work_dir)
    ensure_dir(os.path.dirname(out_tree))

    cleaned_trees = os.path.join(work_dir, "gene_trees.cleaned.nwk")
    logger.info("Preparing trees for Fasturec (strip_regex=%s)...", strip_regex if strip_regex else "NONE")
    clean_tree_labels(trees_file, cleaned_trees, strip_regex)

    run_log = os.path.join(work_dir, "fasturec_run.log")
    cmd = [exe, "-G", os.path.abspath(cleaned_trees), "-Z", "-c", "1"]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores)

    logger.info("Running Fasturec (Supertree): %s", " ".join(cmd))
    try:
        with open(run_log, "w") as lf:
            subprocess.run(cmd, check=True, cwd=work_dir, stdout=lf, stderr=subprocess.STDOUT, env=env)
    except subprocess.CalledProcessError:
        logger.error("Fasturec failed. Check log: %s", run_log)
        return

    candidates = glob.glob(os.path.join(work_dir, "*.tsv"))
    if not candidates:
        logger.error("No output .tsv found in %s", work_dir)
        return

    fu_tsv = max(candidates, key=os.path.getmtime)
    try:
        best_cost, best_tree, _ = parse_fasturec_output(fu_tsv)
        with open(out_tree, "w", encoding="utf-8") as f:
            f.write(best_tree + "\n")
        logger.info("Supertree saved: %s (Cost: %s)", out_tree, best_cost)
    except Exception as e:
        logger.error("Fasturec parse error: %s", e)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", choices=["orthologs", "paralogs"], required=True)
    parser.add_argument("--filtered", action="store_true", help="Use filtered trees (Bonus V.a)")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--fasturec_bin", default="fasturec")
    parser.add_argument("--skip_greedy_consensus", action="store_true", help="Skip IQ-TREE greedy consensus step")
    parser.add_argument("--greedy_minsup", type=float, default=0.0,
                        help="IQ-TREE -minsup (0=extended majority, 0.5=majority-rule)")
    parser.add_argument("--strip_regex", default=DEFAULT_STRIP_REGEX,
                        help="Regex to strip gene suffix in tip labels (default removes __<id>). Use empty string to disable.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    trees_file = pick_gene_trees_file(DEFAULT_GENE_TREES_BASE, args.tag, args.filtered)
    suffix = "_filtered" if args.filtered else ""
    strip_regex = args.strip_regex if args.strip_regex != "" else ""

    logger.info("Processing %s%s using input: %s", args.tag, suffix, trees_file)

    # 1) Consensus (R) + Greedy consensus (IQ-TREE) tylko dla orthologs
    if args.tag == "orthologs":
        out_cons = os.path.join(DEFAULT_OUT_DIR, f"{args.tag}{suffix}_consensus.nwk")
        run_consensus_r("Rscript", DEFAULT_CONS_R, trees_file, out_cons, strict_taxa=False)

        if not args.skip_greedy_consensus:
            out_greedy = os.path.join(DEFAULT_OUT_DIR, f"{args.tag}{suffix}_consensus_greedy_iqtree.nwk")
            work_greedy = os.path.join(DEFAULT_OUT_DIR, f"_{args.tag}{suffix}_iqtree_consensus_work")
            run_iqtree_greedy_consensus(
                trees_file=trees_file,
                out_tree=out_greedy,
                work_dir=work_greedy,
                cores=args.cores,
                minsup=args.greedy_minsup,
                strip_regex=strip_regex,
            )
    else:
        logger.info("Skipping Consensus for Paralogs (methodologically incorrect).")

    # 2) Supertree (always)
    out_super = os.path.join(DEFAULT_OUT_DIR, f"{args.tag}{suffix}_supertree_fasturec.nwk")
    work_dir = os.path.join(DEFAULT_OUT_DIR, f"_{args.tag}{suffix}_fasturec_work")

    run_supertree_fasturec(args.fasturec_bin, trees_file, out_super, work_dir, args.cores, strip_regex=strip_regex)

if __name__ == "__main__":
    main()
