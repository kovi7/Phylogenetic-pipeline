suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(phangorn)))

# ---- logger ----
log_line <- function(level, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("%s [%s] species_consensus: %s\n", ts, level, msg))
}
log_info <- function(msg) log_line("INFO", msg)
log_warn <- function(msg) log_line("WARNING", msg)
log_err  <- function(msg) { log_line("ERROR", msg); quit(status = 2) }

# ---- args parsing ----
args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(name, default = FALSE) {
  hit <- which(args == name)
  if (length(hit) == 0) return(default)
  if (hit < length(args) && !startsWith(args[hit + 1], "--")) {
    v <- tolower(args[hit + 1])
    return(v %in% c("1","true","t","yes","y"))
  }
  TRUE
}

get_value <- function(name, default = NA_character_) {
  hit <- which(args == name)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

trees_file <- get_value("--trees_file", NA_character_)
out_tree   <- get_value("--out_tree", NA_character_)
min_support_raw <- suppressWarnings(as.numeric(get_value("--min_support", "50")))
rooted <- get_flag("--rooted", default = FALSE)
strict_taxa <- get_flag("--strict_taxa", default = TRUE)
use_intersection <- get_flag("--use_intersection", default = FALSE)

if (is.na(trees_file) || is.na(out_tree)) {
  cat("Usage:\n")
  cat("  Rscript 06b_consensus_species_tree.R --trees_file <all_trees.nwk> --out_tree <out.nwk> [--min_support 50] [--rooted FALSE]\n")
  cat("Options:\n")
  cat("  --strict_taxa TRUE|FALSE      require same taxon set in each tree (default TRUE)\n")
  cat("  --use_intersection TRUE|FALSE if taxa differ, prune all trees to shared tips (default FALSE)\n")
  quit(status = 1)
}

if (!file.exists(trees_file)) log_err(paste0("Input tree file not found: ", trees_file))

# support threshold: accept 0.5 or 50
p <- min_support_raw
if (is.na(p)) log_err("min_support must be numeric")
if (p > 1) p <- p / 100
if (p <= 0 || p > 1) log_err("min_support must be in (0,1] or (0,100]")

log_info(paste0("Reading gene trees from: ", trees_file))
trs <- read.tree(trees_file)
if (is.null(trs)) log_err("read.tree returned NULL")
if (inherits(trs, "phylo")) trs <- list(trs)
if (length(trs) < 1) log_err("No trees loaded")

# normalize rooting
if (!rooted) {
  trs <- lapply(trs, function(x) unroot(x))
}

# taxon handling
tipsets <- lapply(trs, function(t) sort(t$tip.label))
ref_tips <- tipsets[[1]]
same_taxa <- vapply(tipsets, function(x) identical(x, ref_tips), logical(1))

if (!all(same_taxa)) {
  n_bad <- sum(!same_taxa)
  log_warn(paste0("Found ", n_bad, " tree(s) with different taxa set"))

  if (use_intersection) {
    shared <- Reduce(intersect, tipsets)
    if (length(shared) < 4) log_err("Intersection of taxa is too small (<4 tips)")
    log_info(paste0("Pruning all trees to shared taxa (n=", length(shared), ")"))
    trs <- lapply(trs, function(t) drop.tip(t, setdiff(t$tip.label, shared)))
  } else if (strict_taxa) {
    trs <- trs[same_taxa]
    log_warn(paste0("Keeping only trees with full taxa set: ", length(trs)))
    if (length(trs) < 1) log_err("No trees left after strict taxa filtering")
  } else {
    log_warn("Proceeding without taxa filtering (may fail or produce misleading consensus)")
  }
}

# compute consensus
log_info(paste0("Building consensus with p=", p))
cons <- consensus(trs, p = p)

# add clade supports (% of trees)
n <- length(trs)
cnt <- prop.clades(cons, trs)
if (!is.null(cnt)) cons$node.label <- as.character(round(100 * cnt / n, 1))

# write result
out_dir <- dirname(out_tree)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.tree(cons, file = out_tree)
log_info(paste0("Saved consensus tree: ", out_tree))
log_info(paste0("Trees used: ", n))