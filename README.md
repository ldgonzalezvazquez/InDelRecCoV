# InDelRecCoV - Indels and Recombination Analysis in Coronaviruses

This repository contains scripts to generate figures for the manuscript analyzing insertion-deletion events (indels) and recombination patterns in the Coronaviridae family.

## Repository Structure

```
InDelRecCoV/
├── scripts/                                    # Python scripts to generate figures
│   ├── generate_figure_1.py                    # Phylogenetic tree with bootstrap
│   ├── generate_figure_3.py                    # Recombination rate plot
│   ├── generate_figure_6.py                    # Ancestral deletion events
│   └── process_indelmip_results.py             # Process InDelMIP output
├── data/                                       # Input data files
│   ├── treerooted.nwk                          # Rooted phylogenetic tree
│   ├── treerooted.raxml.support                # Tree with bootstrap support values
│   ├── outputmip.fasta                         # InDelMIP binary sequences
│   ├── all_events_1_10_100.txt                 # Deletion events data
│   ├── res_full.txt                            # LDhat recombination results
│   ├── aln_coronaviridae.fasta                 # Multiple sequence alignment
│   ├── Indel_ABC_estimation*/                  # R scripts and inputs needed to perform the indel rate estimation for each Coronaviridae and SARS-CoV-2 genome and spike data (link below)
│   │    ├── Coronaviridae_Genome/              
│   │    ├── Coronaviridae_Spike/              
│   │    ├── SARS-CoV-2_Genome/                 
│   │    ├── SARS-CoV-2_Spike/                  
│   │    └── Validation/                        # R scripts and inputs employed to test the method performance for each Coronaviridae and SARS-CoV-2 genome and spike data
│   │         ├── Model_selection/              # model selection of for each Coronaviridae and SARS-CoV-2 genome and spike data
│   │         |    ├── Coronaviridae_Genome/    
│   │         |    ├── Coronaviridae_Spike/    
│   │         |    ├── SARS-CoV-2_Genome/       
│   │         |    └── SARS-CoV-2_Spike/        
│   │         └── INDEL_estimation/             # and indel estimation of for each Coronaviridae and SARS-CoV-2 genome and spike data
│   │              ├── Coronaviridae_Genome/    
│   │              ├── Coronaviridae_Spike/     
│   │              ├── SARS-CoV-2_Genome/       
│   │              └── SARS-CoV-2_Spike/        
│   └──Indel_Simulations/                       # Simulated alignments obtained considering indel events under Binomial and Zipf truncated indel models for each Coronaviridae and SARS-CoV-2 genome and spike data
│       ├── Coronaviridae_Genome/
│       │    ├── Simulated_MSA_Coronaviridae_Genome_BN.tar.xz
│       │    └── Simulated_MSA_Coronaviridae_Genome_Zipf_T.tar.xz
│       ├── Coronaviridae_Spike/
│       │    ├── Simulated_MSA_Coronaviridae_Spike_BN.tar.xz
│       │    └── Simulated_MSA_Coronaviridae_Spike_Zipf_T.tar.xz
│       ├── SARS-CoV-2_Genome/
│       │    ├── Simulated_MSA_SARS-CoV-2_Genome_BN.tar.xz
│       │    └── Simulated_MSA_SARS-CoV-2_Genome_Zipf_T.tar.xz
│       └── SARS-CoV-2_Spike/
│            ├── Simulated_MSA_SARS-CoV-2_Spike_BN.tar.xz
│            └── Simulated_MSA_SARS-CoV-2_Spike_Zipf_T.tar.xz
├── figures/                                    # Output figures
├── requirements.txt                            # Python dependencies
└── README.md                                   # This file
```

* [Indel_ABC_estimation](https://universidadevigo-my.sharepoint.com/:f:/g/personal/luisdaniel_gonzalez_uvigo_gal/IgA_glwF291TS5bPb-CEwF6RAbEg6pFAEf0Qv4nliU9LnnA?e=AiQ8bH)

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

**Script:** `scripts/generate_figure_1.py`

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
- `data/treerooted.raxml.support` - Phylogenetic tree with bootstrap support values

**How to run:**
```bash
python3 scripts/generate_figure_1.py data/treerooted.raxml.support
```

**Output:**
- `figures/Figure_1.pdf` - Final figure in PDF format
- `figures/Figure_1.png` - Final figure in PNG format

**Key features:**
- **Real bootstrap values** from RAxML analysis
- All internal nodes show bootstrap support
- Clean, publication-ready visualization

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
- Used by scripts requiring phylogenetic relationships without bootstrap values

### `data/treerooted.raxml.support`

Phylogenetic tree of Coronavirinae subfamily with bootstrap support values.
- 65 terminal taxa
- Rooted tree structure
- Taxa names format: `Genus-Species_Strain`
- Bootstrap values as internal node labels (0-100)
- Used by `generate_figure_1.py` for tree visualization

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

### `data/aln_coronaviridae.fasta`

Multiple sequence alignment of Coronaviridae genomes.
- Format: FASTA
- Contains aligned genomic sequences for all 65 taxa
- Used as input for phylogenetic and recombination analyses
- Alignment length corresponds to reference genome (~72,000 bp)

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

