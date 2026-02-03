# Phylogenomic Pipeline

Ten pipeline automatyzuje budowę drzew filogenetycznych na podstawie proteomów genomowych. Proces obejmuje pobieranie danych z NCBI, klastrowanie sekwencji białkowych, identyfikację rodzin genów (ortologów i paralogów), wielosekwencyjne uliniowienie (MSA), budowę drzew genowych oraz inferencję drzew gatunków metodami konsensusowymi i superdrzewowymi.

Pipeline jest zaprojektowany dla zestawów genomów (np. Gammaproteobacteria), ale działa dla dowolnej listy proteomów NCBI. Obecna wersja zawiera ulepszenia w przetwarzaniu równoległym, obsłudze błędów i modularności.

## Wymagania systemowe

- **System operacyjny**: Linux/macOS (Windows via WSL)
- **Pamięć RAM**: Minimum 16 GB, zalecane 32 GB+ dla dużych zestawów danych
- **Dysk**: Kilka GB wolnego miejsca na dane wyjściowe
- **CPU**: Wielordzeniowy (zalecane 8+ wątków)

## Wymagania oprogramowania

### Narzędzia z Conda (środowisko `environment.yml`)
- **NCBI Datasets CLI**: Pobieranie proteomów z NCBI
- **MMseqs2**: Szybkie klastrowanie sekwencji białkowych
- **MAFFT**: Wielosekwencyjne uliniowienie (MSA)
- **IQ-TREE / IQ-TREE2**: Budowa drzew filogenetycznych (drzewa genowe + konsensus)
- **R** z pakietami: `ape`, `phangorn`, `TreeDist` (konsensus + wizualizacja)
- **Python** z pakietami: `biopython`, `tqdm` (skrypty pomocnicze)

### Narzędzia zewnętrzne (poza Conda)
- **Fasturec**: Algorytm do supertree (należy skompilować z C++). Ścieżka ustawiana w `main.sh` jako `FASTUREC_BIN`.

## Struktura katalogow
- `data/`
  - `accessions.txt` - lista identyfikatorów genomów NCBI
- `scripts/`:
  - `01_proteoms_downloader.py` - pobieranie proteomow (NCBI Datasets CLI)
  - `01b_mapper_unique.py` - unikalne ID bialek + mapowania
  - `02_mmseq_clustering.sh` - MMseqs2 clustering
  - `03_make_families.py` - budowa rodzin (orthologs / paralogs)
  - `04_mafft_msa.py` - MSA dla rodzin
  - `05_iqtree_gene_trees.py` - drzewa genowe + scalanie do `all_trees.nwk`
  - `06_species_trees.py` - drzewa gatunkow (consensus + supertree)
  - `06b_consensus_species_tree.R` - consensus majority (R)
  - `07_visualize.R` - porownania (nRF, nMCI) + PDFy
- `main.sh` - uruchamia caly pipeline krok po kroku

---

## Instalacja

1. **Utwórz środowisko Conda**:
   ```bash
   conda env create -f environment.yml
   conda activate genomics_pipeline
   ```

2. **Zainstaluj Fasturec**:
   - Pobierz i skompiluj Fasturec (np. https://bitbucket.org/pgor17/fasturec/downloads/).
   - Ustaw ścieżkę w `main.sh`: `FASTUREC_BIN="/path/to/fasturec"`

3. **Sprawdź instalację**:
   ```bash
   which datasets mmseqs mafft iqtree R python
   ```

## Przygotowanie danych wejściowych

1. **Plik accessions.txt**: Lista identyfikatorów genomów NCBI (GCF_... lub GCA_...), jeden na linię. Opcjonalnie etykieta po średniku.
   Przykład:
   ```
   GCF_000005845.2; Escherichia coli K-12
   GCF_000006945.2; Salmonella enterica
   ```

2. **Umieść plik** w `data/accessions.txt`.

## Uruchomienie pipeline'u

Pipeline jest uruchamiany za pomocą `main.sh`. Skrypt wykonuje kroki sekwencyjnie. Można uruchamiać pojedyncze skrypty dla debugowania.

### Pełne uruchomienie
```bash
./main.sh
```

### Konfiguracja
- Edytuj `main.sh` dla liczby wątków (`THREADS=8`) i ścieżki Fasturec.
- Większość kroków jest równoległa (używa `--threads`).

### Kroki pipeline'u

1. **Pobieranie proteomów** (`01_proteoms_downloader.py`):
   - Pobiera ZIP archiwa proteomów z NCBI.
   - Wyciąga pliki `protein.faa`.
   - Wyjście: `data/proteome_archives/`, `data/proteome_fasta/`

2. **Mapowanie i łączenie** (`01b_mapper_unique.py`):
   - Tworzy unikalne identyfikatory białek (genome__protein).
   - Łączy wszystkie sekwencje w jeden FASTA.
   - Wyjście: `data/maps/protein2genome.tsv`, `data/merged_fasta/merged_proteins.faa`

3. **Klastrowanie sekwencji** (`02_mmseq_clustering.sh`):
   - Grupuje podobne sekwencje za pomocą MMseqs2.
   - Wyjście: `clustering/clu_cluster.tsv`, `clustering/clu_rep_seq.fasta`

4. **Tworzenie rodzin genów** (`03_make_families.py`):
   - Identyfikuje rodziny ortologów (1 sekwencja/genom) i paralogów (wiele sekwencji/genom).
   - Filtruje rodziny na podstawie pokrycia genomów.
   - Wyjście: `families/orthologs/`, `families/paralogs/`, `families/families_stats.tsv`

5. **Wielosekwencyjne uliniowienie** (`04_mafft_msa.py`):
   - Uliniawia sekwencje w rodzinach za pomocą MAFFT.
   - Wyjście: `allignment/orthologs/`, `allignment/paralogs/`

6. **Budowa drzew genowych** (`05_iqtree_gene_trees.py`):
   - Inferuje drzewa filogenetyczne dla każdej rodziny za pomocą IQ-TREE.
   - Wyjście: `trees/gene_trees/`

7. **Budowa drzew gatunków** (`06_species_trees.py`, `06b_consensus_species_tree.R`):
   - Metoda konsensusowa (greedy consensus z IQ-TREE, majority-rule z ape).
   - Metoda superdrzewowa (Fasturec).
   - Wyjście: `trees/species_trees/`

8. **Wizualizacja** (`07_visualize.R`):
   - Porównuje drzewa gatunków, generuje wykresy.
   - Wyjście: Pliki PDF w `trees/`

## Wyjścia

- **Drzewa gatunków**: `trees/species_trees/` (formaty Newick .nwk)
- **Drzewa genowe**: `trees/gene_trees/` (.treefile)
- **Statystyki**: `families/families_stats.tsv`, logi w konsoli
- **Dane pośrednie**: FASTA, uliniowienia, klastry w odpowiednich katalogach

## Dostosowywanie

- **Liczba genomów**: Dostosuj `--min_genomes` w `03_make_families.py`.
- **Filtry**: Użyj `--require_all_genomes` dla ścisłych ortologów.
- **Tryb wyboru ortologów**: W `03_make_families.py` zmień stałą `ORTHO_PICK_MODE` na "closest" (wybór najbliższej sekwencji) lub "longest" (najdłuższa sekwencja) zamiast domyślnego "strict", żeby móc użyć metod generowania rodzin ortologicznych.
- **Bootstrap dla drzew genowych**: W `05_iqtree_gene_trees.py` użyj `--ufboot 1000` dla drzew z wsparciem bootstrap (np. 1000 replik). Następnie użyj tych drzew w `06_species_trees.py` dla budowy drzew gatunków z uwzględnieniem wsparcia.
- **Równoległość**: Zwiększ `--threads` dla szybszego przetwarzania.
- **Debugowanie**: Uruchamiaj skrypty pojedynczo z `--log_level DEBUG`.
