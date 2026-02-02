# mirMachine
mirMachine (2021): Automated pipeline for the annotation of high-confidence miRNAs from genomic sequences

## Key Highlights
- Genome-wide microRNA (miRNA) discovery and annotation pipeline for plant genomes  
- Rule-based homology and structural validation for reproducibility  
- Publication: [JoVE (2021)](https://app.jove.com/t/62430/mirmachine-a-one-stop-shop-for-plant-mirna-annotation)  
- Produces interpretable miRNA annotations for downstream genomic and regulatory analysis

**Test data and outputs:** [Releases/v1 (TestData)](https://github.com/hbusra/mirMachine/releases/tag/v1)

## Overview
mirMachine is a fully automated, rule-based computational pipeline for **known and novel miRNA identification and annotation**. It combines sequence homology with hairpin structural validation to generate **high-confidence annotations** suitable for:
- Comparative genomics
- Functional analysis
- Genome annotation workflows

Benchmarked on Arabidopsis thaliana and wheat, mirMachine delivers high sensitivity and reduced false positives.

## Core Features
- Homology-based identification of known miRNAs  
- Rule-based discovery of novel miRNAs  
- RNA secondary structure validation (hairpin enforcement)  
- Genome-wide miRNA localization  
- Optional sRNA-seq integration for expression evidence and novel miRNA discovery  
- Fully automated execution

## Quick Start
1. Install dependencies: BLAST+, RNAfold, Perl  
2. Add mirMachine scripts to your PATH  
3. Run test data:

```bash
# Homology-based run
bash mirMachine_submit.sh -f genome.fasta -i mature_high_conf_v22_1.fa.filtered.fasta -n 10

# With sRNA-seq
bash mirMachine_submit.sh -f genome.fasta -i smallRNA.fa -sRNAseq -lmax 24 -lmin 19 -RPM 10
```

**Outputs**
- Tabular annotations with genomic coordinates (*.tbl)
- Mature miRNA and pre-miRNA (hairpin) FASTA files
- Log files reporting QC, warnings, and pipeline diagnostics

Step-by-step installation and usage are available in the [JoVE video](https://app.jove.com/v/62430/mirmachine-a-one-stop-shop-for-plant-mirna-annotation).


## Modes of Operation

**1. Homology-Based miRNA Identification**
- Requires only a genome FASTA
- Predicts genome-wide miRNA candidates
- Outputs mature and precursor sequences plus genomic coordinates

**2. Novel miRNA Discovery with sRNA-seq (Optional)**
- Pre-processing of sRNA-seq reads (adapter trimming and FASTQ→FASTA conversion) required
- Supports abundance filtering and mismatch control


## Benchmarking
- **Higher sensitivity** than miRDP2 in Arabidopsis
- **Improved performance** in wheat when combining homology and expression evidence
- **Reduced false positives** in complex plant genomes

## Citation
Cagirici et al., 
“mirMachine: a one-stop shop for plant miRNA annotation” 
Journal of Visualized Experiments (JoVE), 2021
