# mirMachine (2021): Genome-wide miRNA discovery and annotation pipeline for plant genomes  
[![DOI](https://img.shields.io/badge/PubMed-2021-blue)](https://pubmed.ncbi.nlm.nih.gov/33999024/)
![Perl](https://img.shields.io/badge/Perl-5.26+-purple)
![License](https://img.shields.io/badge/License-MIT-green)

> **📄 Publication:** mirMachine: A One-Stop Shop for Plant miRNA Annotation — [*Journal of Visualized Experiments (JoVE), 2021*](https://pubmed.ncbi.nlm.nih.gov/33999024/)<br/>
> **👤 Role:** First author and primary developer — designed algorithm and implemented full pipeline  
> **🎯 Impact:** Automated, genome-wide identification of known and novel miRNAs from plant genomes<br/>
> **Tech:** Perl • BLAST+ • RNAfold

## Overview
mirMachine is a fully automated, rule-based computational pipeline for **known and novel miRNA identification and annotation**. It combines sequence homology with hairpin structural validation to identify **high-confidence miRNA** candidates suitable for comparative genomics and functional analysis.

**Use cases:** Genome annotation projects • miRNA family expansion studies

## Methods
- **Algorithm:** Homology search (BLAST vs. miRBase) → Precursor extraction → Structural validation (RNAfold hairpin prediction) → Genomic mapping
- **Validation:** Benchmarked on *Arabidopsis thaliana* (TAIR10) and wheat (IWGSC RefSeq v2.0) genomes
- **Data:** Test data and example outputs available at [Releases/v1](https://github.com/hbusra/mirMachine/releases/tag/v1)

### Pipeline Overview
1. Homology-based candidate identification (BLAST vs. miRBase)
2. Precursor extraction and hairpin folding (RNAfold)
3. Rule-based filtering (structure + thermodynamics)
4. Optional: sRNA-seq integration for expression evidence and novel miRNA discovery
5. Genomic coordinate mapping and annotation output


## Installation
Step-by-step installation and usage are available in the JoVE video protocol.<br/>

**Dependencies:** BLAST+ (2.9+), ViennaRNA (RNAfold), Perl 5.26+

```bash
# Clone repository
git clone https://github.com/hbusra/mirMachine.git
cd mirMachine

# Add to PATH
export PATH=$PATH:$(pwd)

# Verify installation
bash mirMachine_submit.sh --help
```

### Conda (Recommended):
```bash
conda create -n mirmachine -c bioconda blast perl viennarna
conda activate mirmachine
```

## Quick Start

**Homology-based miRNA identification:**
```bash
bash mirMachine_submit.sh \
  -f genome.fasta \
  -i mature_high_conf_v22_1.fa.filtered.fasta \
  -n 10
```

**Novel miRNA Discovery with sRNA-seq:**
```bash
bash mirMachine_submit.sh \
  -f genome.fasta \
  -i smallRNA.fa \
  -sRNAseq \
  -lmax 24 -lmin 19 \
  -RPM 10
```

**Outputs:**
- `*.tbl` — Tabular annotations with genomic coordinates
- `*_mature.fa` — Mature miRNA sequences
- `*_hairpin.fa` — Precursor (pre-miRNA) sequences
- `*.log` — Pipeline diagnostics

## Notes
  - Designed specifically for plant genomes
  - Conservative filtering favors precision over exhaustive recall
  - Runtime scales with genome size and BLAST search space

## Citation
Cagirici et al., “mirMachine: a one-stop shop for plant miRNA annotation,” Journal of Visualized Experiments (JoVE), 2021. [PubMed](https://pubmed.ncbi.nlm.nih.gov/33999024/)

## Resources
- **Video protocol:** app.jove.com/v/62430/mirmachine-a-one-stop-shop-for-plant-mirna-annotation
- **Test data:** [Download from Releases/v1](https://github.com/hbusra/mirMachine/releases/tag/v1)
