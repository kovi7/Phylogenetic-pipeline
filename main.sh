#!/usr/bin/env bash
set -euo pipefail

# --- KONFIGURACJA ---
FASTUREC_BIN="/mnt/c/Users/JK/Desktop/GP/projekt/sharma2022/fasturec/fasturec"
THREADS=8
# --------------------

# echo ">>> PIPELINE START"

# # 1. Pobieranie danych
# echo ">>> STEP 1: Downloading & Mapping"
# python scripts/01_proteoms_downloader.py -i data/accessions.txt
# # 1b. Mapowanie nazw
# python scripts/01b_mapper_unique.py -i data/accessions.txt --force

# # 2. Klastrowanie (MMseqs2)
# echo ">>> STEP 2: Clustering"
# rm -rf clustering
# bash scripts/02_mmseq_clustering.sh data/merged_fasta/merged_proteins.faa clustering

# # # 3. Tworzenie rodzin
# echo ">>> STEP 3: Building Families (Orthologs & Paralogs)"
# rm -rf families/orthologs families/paralogs
# python scripts/03_make_families.py \
#    --tag both \
#    --clusters_tsv clustering/clu_cluster.tsv \
#    --threads "$THREADS" \
#    --require_all_genomes \
#    --force

# 4. Multiuliniowienie (MAFFT)
# echo ">>> STEP 4: MSA (MAFFT)"
# rm -rf allignment
# python scripts/04_mafft_msa.py --tag orthologs --workers "$THREADS"
# python scripts/04_mafft_msa.py --tag paralogs --workers "$THREADS"

# # 5. Drzewa Genów (IQ-TREE)
# echo ">>> STEP 5: Gene Trees (IQ-TREE)"
# python scripts/05_iqtree_gene_trees.py --tag orthologs --workers "$THREADS" --threads_per_job 2 --force
# python scripts/05_iqtree_gene_trees.py --tag paralogs --workers "$THREADS" --threads_per_job 2 --force

# # 6. Drzewa Gatunków (Consensus & Supertree)
# echo ">>> STEP 6: Species Trees"

# --- ORTOLOGI ---
echo "--- Orthologs: Classic ---"
python scripts/06_species_trees.py --tag orthologs --cores "$THREADS" --fasturec_bin "$FASTUREC_BIN"

# # --- PARALOGI (V. b) ---
echo "--- Paralogs: Classic ---"
python scripts/06_species_trees.py --tag paralogs --cores "$THREADS" --fasturec_bin "$FASTUREC_BIN"

# # 7. Wizualizacja drzew + odległości 
echo "--- Visualization ---"
Rscript scripts/07_visualize.R

echo ">>> PIPELINE COMPLETED SUCCESSFULLY"