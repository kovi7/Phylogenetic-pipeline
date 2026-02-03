suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(phangorn)))
suppressMessages(suppressWarnings(library(TreeDist)))

# ---- logger ----
log_line <- function(level, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("%s [%s] tree_report: %s\n", ts, level, msg))
}
log_info <- function(msg) log_line("INFO", msg)
log_warn <- function(msg) log_line("WARNING", msg)
log_err  <- function(msg) { log_line("ERROR", msg); quit(status = 2) }

args <- commandArgs(trailingOnly = TRUE)

get_value <- function(name, default = NA_character_) {
  hit <- which(args == name)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

# Defaults
consensus_path        <- get_value("--consensus",        "trees/species_trees/orthologs_consensus.nwk")
consensus_greedy_path <- get_value("--consensus_greedy", "trees/species_trees/orthologs_consensus_greedy_iqtree.nwk")
super_ortho_path      <- get_value("--super_ortho",      "trees/species_trees/orthologs_supertree_fasturec.nwk")
super_para_path       <- get_value("--super_para",       "trees/species_trees/paralogs_supertree_fasturec.nwk")
timetree_path         <- get_value("--timetree",         "trees/timetree_tree.nwk")
article_path          <- get_value("--article",          "trees/article_tree.nwk")
out_dir               <- get_value("--out_dir",          "results/tree_report")

# ---- io checks ----
paths <- c(consensus_path, consensus_greedy_path, super_ortho_path, super_para_path, timetree_path, article_path)
names(paths) <- c("consensus","consensus_greedy","super_ortho","super_para","timetree","article")
for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) log_err(paste0("Missing file for --", nm, ": ", paths[[nm]]))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- read trees ----
read_one_tree <- function(p) {
  tr <- read.tree(p)
  if (is.null(tr)) log_err(paste0("Failed to read tree: ", p))
  if (inherits(tr, "multiPhylo")) {
    if (length(tr) < 1) log_err(paste0("No trees in file: ", p))
    tr <- tr[[1]]
  }
  tr
}

consensus        <- read_one_tree(consensus_path)
consensus_greedy <- read_one_tree(consensus_greedy_path)
super_ortho      <- read_one_tree(super_ortho_path)
super_para       <- read_one_tree(super_para_path)
timetree         <- read_one_tree(timetree_path)
article          <- read_one_tree(article_path)

# Unroot
consensus_u        <- unroot(consensus)
consensus_greedy_u <- unroot(consensus_greedy)
super_ortho_u      <- unroot(super_ortho)
super_para_u       <- unroot(super_para)
timetree_u         <- unroot(timetree)
article_u          <- unroot(article)

fmt <- function(x) sprintf("%.3f", x)

prune_to_common <- function(a, b) {
  common <- intersect(a$tip.label, b$tip.label)
  if (length(common) < 4) return(list(a = NULL, b = NULL, n = length(common)))
  a2 <- drop.tip(a, setdiff(a$tip.label, common))
  b2 <- drop.tip(b, setdiff(b$tip.label, common))
  list(a = a2, b = b2, n = length(common))
}

score_pair <- function(ref_tr, tr) {
  pr <- prune_to_common(ref_tr, tr)
  if (is.null(pr$a)) return(list(nrf = NA_real_, nmci = NA_real_, n = pr$n))
  nrf  <- RF.dist(pr$a, pr$b, normalize = TRUE)
  nmci <- MutualClusteringInfo(pr$a, pr$b, normalize = TRUE)
  list(nrf = nrf, nmci = nmci, n = pr$n)
}

# ---- Result trees (columns) ----
res_trees <- list(
  consensus_majority = consensus_u,
  consensus_greedy   = consensus_greedy_u,
  super_ortho        = super_ortho_u,
  super_para         = super_para_u
)

# ---- Reference trees (rows) ----
refs <- list(
  article  = article_u,
  timetree = timetree_u
)

# ---- summary table ----
tab <- data.frame(reference = c("article", "timetree"), stringsAsFactors = FALSE)

for (nm in names(res_trees)) {
  tab[[paste0(nm, "_nRF")]]   <- "-"
  tab[[paste0(nm, "_nMCI")]]  <- "-"
  tab[[paste0(nm, "_nTips")]] <- "-"
}

for (i in seq_len(nrow(tab))) {
  ref_name <- tab$reference[i]
  ref_tr <- refs[[ref_name]]

  for (nm in names(res_trees)) {
    tr <- res_trees[[nm]]
    sc <- score_pair(ref_tr, tr)
    tab[i, paste0(nm, "_nRF")]   <- ifelse(is.na(sc$nrf),  "-", fmt(sc$nrf))
    tab[i, paste0(nm, "_nMCI")]  <- ifelse(is.na(sc$nmci), "-", fmt(sc$nmci))
    tab[i, paste0(nm, "_nTips")] <- as.character(sc$n)
  }
}

# Extra: timetree vs article only in article row
tab[["timetree_vs_article_nRF"]]   <- "-"
tab[["timetree_vs_article_nMCI"]]  <- "-"
tab[["timetree_vs_article_nTips"]] <- "-"

sc_ta <- score_pair(article_u, timetree_u)
tab[tab$reference == "article", "timetree_vs_article_nRF"]   <- ifelse(is.na(sc_ta$nrf),  "-", fmt(sc_ta$nrf))
tab[tab$reference == "article", "timetree_vs_article_nMCI"]  <- ifelse(is.na(sc_ta$nmci), "-", fmt(sc_ta$nmci))
tab[tab$reference == "article", "timetree_vs_article_nTips"] <- as.character(sc_ta$n)

write.csv2(tab, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

# ---- plotting ----
plot_rect <- function(tr, main, show_support = FALSE) {
  plot(tr, type = "phylogram", use.edge.length = FALSE,
       main = "", cex = 0.75, label.offset = 0.01, no.margin = FALSE)

  mtext(main, side = 3, line = 1, font = 2, cex = 1)

  if (show_support && !is.null(tr$node.label)) {
    nodelabels(tr$node.label, frame = "rect", bg = "white",
               cex = 0.75, adj = c(0.5, -0.2))
  }
}

# Single-tree PDFs
pdf(file.path(out_dir, "tree_article.pdf"), width = 10, height = 7)
plot_rect(article_u, "Article tree")
dev.off()

pdf(file.path(out_dir, "tree_timetree.pdf"), width = 10, height = 7)
plot_rect(timetree_u, "TimeTree")
dev.off()

pdf(file.path(out_dir, "tree_consensus_majority.pdf"), width = 10, height = 7)
plot_rect(consensus_u, "Consensus (majority)", show_support = TRUE)
dev.off()

pdf(file.path(out_dir, "tree_consensus_greedy.pdf"), width = 10, height = 7)
plot_rect(consensus_greedy_u, "Consensus (greedy)", show_support = TRUE)
dev.off()

pdf(file.path(out_dir, "tree_supertree_ortho.pdf"), width = 10, height = 7)
plot_rect(super_ortho_u, "Supertree (orthologs)")
dev.off()

pdf(file.path(out_dir, "tree_supertree_paralogs.pdf"), width = 10, height = 7)
plot_rect(super_para_u, "Supertree (paralogs)")
dev.off()

# Panels: consensuses 
pdf(file.path(out_dir, "fig1_article_consensuses.pdf"), width = 16, height = 7)
par(mfrow = c(1, 3), mar = c(1, 1, 5, 1), xpd = NA)
plot_rect(article_u, "Article tree")
plot_rect(consensus_u, "Consensus (majority)", show_support = TRUE)
plot_rect(consensus_greedy_u, "Consensus (greedy)", show_support = TRUE)
dev.off()

pdf(file.path(out_dir, "fig2_timetree_consensuses.pdf"), width = 16, height = 7)
par(mfrow = c(1, 3), mar = c(1, 1, 5, 1), xpd = NA)
plot_rect(timetree_u, "TimeTree")
plot_rect(consensus_u, "Consensus (majority)", show_support = TRUE)
plot_rect(consensus_greedy_u, "Consensus (greedy)", show_support = TRUE)
dev.off()

# Panels: supertrees
pdf(file.path(out_dir, "fig3_article_supertrees.pdf"), width = 16, height = 7)
par(mfrow = c(1, 3), mar = c(1, 1, 5, 1), xpd = NA)
plot_rect(article_u, "Article tree")
plot_rect(super_ortho_u, "Supertree (orthologs)")
plot_rect(super_para_u, "Supertree (paralogs)")
dev.off()

pdf(file.path(out_dir, "fig4_timetree_supertrees.pdf"), width = 16, height = 7)
par(mfrow = c(1, 3), mar = c(1, 1, 5, 1), xpd = NA)
plot_rect(timetree_u, "TimeTree")
plot_rect(super_ortho_u, "Supertree (orthologs)")
plot_rect(super_para_u, "Supertree (paralogs)")
dev.off()

log_info(paste0("Done. Outputs saved to: ", normalizePath(out_dir)))