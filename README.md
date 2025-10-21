# InDelRecCoV - Indels and Recombination Analysis in Coronaviruses

This repository contains scripts to generate figures for the manuscript analyzing insertion-deletion events (indels) and recombination patterns in the Coronaviridae family.

## Repository Structure

```
InDelRecCoV/
├── scripts/                          # Python scripts to generate figures
│   ├── tree_apply_bootstrap.py       # Calculate bootstrap support values
│   ├── generate_figure_1.py          # Phylogenetic tree with bootstrap
│   ├── generate_figure_3.py          # Recombination rate plot
│   ├── generate_figure_6.py          # Ancestral deletion events
│   └── process_indelmip_results.py   # Process InDelMIP output
├── data/                             # Input data files
│   ├── treerooted.nwk                # Rooted phylogenetic tree
│   ├── RAxML_tree.raxml.bootstraps   # Bootstrap replicate trees
│   ├── names_correpondence.txt       # Bootstrap name mappings
│   ├── treerooted_with_bootstrap.nwk # Tree with bootstrap values
│   ├── outputmip.fasta               # InDelMIP binary sequences
│   ├── all_events_1_10_100.txt       # Deletion events data
│   ├── res_full.txt                  # LDhat recombination results
│   └── aln_coronaviridae.fasta       # Multiple sequence alignment
├── figures/                          # Output figures
├── requirements.txt                  # Python dependencies
└── README.md                         # This file
```

## Requirements

### Python Dependencies

```bash
pip install biopython matplotlib pandas numpy seaborn ete3 progress
```

Or install from requirements file:
```bash
pip install -r requirements.txt
```

Required packages:
- **Python 3.7+**
- `biopython` - For phylogenetic tree manipulation and sequence parsing
- `matplotlib` - For plotting
- `pandas` - For data manipulation
- `numpy` - For numerical operations
- `seaborn` - For enhanced plot styling
- `ete3` - For phylogenetic tree analysis (InDelMIP processing)
- `progress` - For progress bars (InDelMIP processing)

## Data Processing

### InDelMIP Results Processing

**Script:** `scripts/process_indelmip_results.py`

**Description:** Processes raw InDelMIP output to identify and classify deletion events across the phylogenetic tree. This interactive script:
- Preprocesses binary sequences (0s and 1s) representing deletions and nucleotide regions
- Detects continuous deletion regions (runs of 0s)
- Classifies deletions as ancestral (inherited) or new (lineage-specific)
- Performs sensitivity analysis for parameter selection
- Generates filtered output files for downstream analysis

**Input data:**
- `data/outputmip.fasta` - InDelMIP output in FASTA format
  - Sequences contain binary strings: 0 = deletion, 1 = conserved
  - One sequence per node in the phylogenetic tree
- `data/treerooted.nwk` - Phylogenetic tree

**Output files:**
- `all_events.txt` - All deletion events detected
- `new_events.txt` - Only new (non-ancestral) deletion events

**How to run:**
```bash
cd scripts
python3 process_indelmip_results.py
```

**Interactive parameters:**
The script will prompt you to select:
1. **Preprocessing error percentage** (0-100%): Determines how large regions of 1s within deletions can be to be ignored (merged into the deletion)
2. **Error tolerance percentage** (0-100%): Determines how similar two deletions must be to be considered the same event
3. **Minimum deletion size** (bp): Filters out small deletion events below this threshold

**Output format:**
Both output files are tab-separated with columns:
- `Node`: Node identifier in the phylogenetic tree
- `Deletion_Start`: Start position of deletion
- `Deletion_End`: End position of deletion
- `Size`: Size of deletion in base pairs
- `State`: Status (New, Ancestral)

---

## Figures

### Figure 1: Phylogenetic Tree of Coronavirinae

**Scripts:** 
1. `scripts/tree_apply_bootstrap.py` - Calculate bootstrap values
2. `scripts/generate_figure_1.py` - Generate tree figure

**Description:** Generates a phylogenetic tree with:
- Taxa colored by coronavirus genus (Alphacoronavirus, Betacoronavirus, Deltacoronavirus, Gammacoronavirus, Unclassified)
- **Real bootstrap support values** from RAxML analysis (0-100)
- Bootstrap values displayed and colored by confidence level:
  - ≥90: Green
  - 75-89: Blue
  - 50-74: Orange
  - 25-49: Red
  - <25: Purple

**Input data:**
- `data/treerooted.nwk` - Rooted phylogenetic tree in Newick format
- `data/RAxML_tree.raxml.bootstraps` - 100 bootstrap replicate trees from RAxML
- `data/names_correpondence.txt` - TSV file mapping bootstrap codes to full species names

**How to run (2 steps):**

**Step 1: Calculate bootstrap values**
```bash
cd data
python3 ../scripts/tree_apply_bootstrap.py \
  --ref treerooted.nwk \
  --boot RAxML_tree.raxml.bootstraps \
  --aliases names_correpondence.txt \
  --out treerooted_with_bootstrap.nwk \
  --log bootstrap_mapping.log
```

This creates `treerooted_with_bootstrap.nwk` with bootstrap percentages as internal node labels.

**Step 2: Generate the figure**
```bash
cd ../scripts
python3 generate_figure_1.py
```

This reads `treerooted_with_bootstrap.nwk` and generates the final figure.

**Output:**
- `data/treerooted_with_bootstrap.nwk` - Tree with bootstrap values
- `data/bootstrap_mapping.log` - Detailed log of bootstrap calculation
- `figures/Figure_1.pdf` - Final figure in PDF format
- `figures/Figure_1.png` - Final figure in PNG format

**Key features:**
- **Real bootstrap values** calculated from RAxML output (not simulated)
- All 64 internal nodes show bootstrap support
- Topology and branch lengths from reference tree preserved
- Automatic resolution of naming mismatches between bootstrap and reference trees

---

### Figure 3: Population Recombination Rate along the Genome

**Script:** `scripts/generate_figure_3.py`

**Description:** Visualizes population recombination rate (ρ) across the coronavirus genome with:
- Mean recombination rate line
- 95% confidence intervals (shaded area)
- Gene annotations bar showing all coding regions
- Consistent gene colors across figures

**Input data:**
- `data/res_full.txt` - LDhat output file with recombination rate estimates
  - Columns: Loci, Mean_rho, L95, U95
  - Tab-separated values
  - 71,920 genomic positions

**Output:**
- `figures/Figure_3.pdf`
- `figures/Figure_3.png`

**How to run:**
```bash
cd scripts
python3 generate_figure_3.py
```

---

### Figure 6: Ancestral Deletion Events in Coronavirinae

**Script:** `scripts/generate_figure_6.py`

**Description:** Visualizes ancestral deletion events mapped onto the phylogenetic tree with:
- Phylogenetic tree structure showing evolutionary relationships
- Deletion events plotted along genome positions
- Color-coded deletions by evolutionary state (ancestral, derived, terminal)
- Gene annotations showing where deletions occurred
- Two visualization modes: full tree and reduced genus-level view

**Input data:**
- `data/all_events_1_10_100.txt` - Deletion events data file
  - Columns: Node, Deletion_Start, Deletion_End, Size, State
  - States: ancestral, new, terminal
- `data/treerooted.nwk` - Phylogenetic tree

**Output:**
- `figures/Figure_6.pdf`
- `figures/Figure_6.png`

**How to run:**
```bash
cd scripts
python3 generate_figure_6.py
```

**Key features:**
- Combines phylogenetic tree with deletion event mapping
- Shows genome positions (0-72,000 bp) on horizontal axis
- Color-coded by deletion type:
  - Ancestral deletions: shared by multiple lineages
  - Derived deletions: lineage-specific
  - Terminal deletions: species-specific
- Gene structure overlay showing functional regions affected
- Optional reduced visualization focusing on genus-level patterns

---

## Data Files

### `data/treerooted.nwk`

Rooted phylogenetic tree of Coronavirinae subfamily in Newick format.
- 65 terminal taxa
- Rooted tree structure
- Taxa names format: `Genus-Species_Strain`
- Used by all scripts requiring phylogenetic relationships

### `data/outputmip.fasta`

InDelMIP output file containing binary sequences for deletion detection.
- Format: FASTA with binary sequences (0 = deletion, 1 = conserved)
- One sequence per node (both terminal and internal) in the phylogenetic tree
- Sequence length corresponds to the reference genome alignment
- Used for identifying and classifying ancestral vs. lineage-specific deletions
- Processed by `process_indelmip_results.py` to generate deletion event files

### `data/res_full.txt`

LDhat analysis results containing population recombination rate estimates.
- Format: Tab-separated values
- Columns:
  - `Loci`: Genomic position (1-71,920)
  - `Mean_rho`: Mean population recombination rate
  - `L95`: Lower 95% confidence bound
  - `U95`: Upper 95% confidence bound
- Total positions: 71,920
- Polymorphic sites: 47,226
- Mean ρ (polymorphic sites): 0.4417

### `data/all_events_1_10_100.txt`

Ancestral deletion events identified across the Coronavirinae phylogeny.
- Format: Tab-separated values
- Columns:
  - `Node`: Node identifier in the phylogenetic tree
  - `Deletion_Start`: Start position of deletion in genome
  - `Deletion_End`: End position of deletion in genome
  - `Size`: Length of deletion event (bp)
  - `State`: Evolutionary state (ancestral, new, terminal)
- Contains insertion-deletion (indel) events mapped to phylogenetic nodes
- Used to visualize evolutionary patterns of gene loss
- Generated by `process_indelmip_results.py` or provided pre-computed

### `data/RAxML_tree.raxml.bootstraps`

RAxML bootstrap replicate trees for phylogenetic support calculation.
- Format: Newick format (multiple trees separated by semicolons)
- Contains 100 bootstrap replicate trees
- Taxa names use short codes (e.g., AA, AB, AC...) that need mapping
- Used by `tree_apply_bootstrap.py` to calculate bootstrap support values
- Generated from RAxML phylogenetic analysis

### `data/names_correpondence.txt`

Mapping file between bootstrap tree codes and full species names.
- Format: Tab-separated values (TSV)
- Two columns: `bootstrap_code<TAB>full_species_name`
- Example: `AA	U-BatCoV1`, `CE	B-SARS-CoV-2`
- Required to correctly assign bootstrap values to the reference tree
- 65 mappings for all taxa in the phylogeny

### `data/aln_coronaviridae.fasta`

Multiple sequence alignment of Coronaviridae genomes.
- Format: FASTA
- Contains aligned genomic sequences for all 65 taxa
- Used as input for phylogenetic and recombination analyses
- Alignment length corresponds to reference genome (~72,000 bp)

### `data/treerooted_with_bootstrap.nwk` (generated)

Phylogenetic tree annotated with bootstrap support values.
- Format: Newick with bootstrap values as internal node labels
- Generated by `tree_apply_bootstrap.py`
- Contains same topology as `treerooted.nwk` but with bootstrap percentages
- Used by `generate_figure_1.py` to create the final tree figure

## Gene Annotations

The scripts use a consistent gene color scheme based on a custom palette. Gene positions are hardcoded in the scripts and include:

- **5'UTR** and **3'UTR**
- **ORF1ab** (major polyprotein)
- **Structural proteins**: S (Spike), E (Envelope), M (Membrane), N (Nucleocapsid)
- **Accessory proteins**: NS2, HE, NS4, NS6, etc.
- **Other ORFs**: Multiple ORF2-11 genes, some duplicated in different genome regions

## Citation

If you use these scripts or data, please cite:

[Manuscript in preparation - citation will be added upon publication]

## Authors

- Luis Daniel González-Vázquez
- David Ferreiro
- Filipe Pereira
- Miguel Arenas


## License

[To be determined]

## Contact

For questions or issues, please mail luisdaniel.gonzalez@uvigo.es

