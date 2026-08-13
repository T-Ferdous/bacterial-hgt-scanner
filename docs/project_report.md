# Comprehensive Project Report: Bacterial Horizontal Gene Transfer (HGT) & Codon Usage Bias (CUB) Detection Pipeline

**Author:** Tasmia Ferdous  
**Date:** August 2026  
**Target Genome:** *Escherichia coli* CFT073 (NCBI RefSeq: NC_004431.1)  

---

## 1. Executive Summary

Horizontal Gene Transfer (HGT) is a major evolutionary mechanism in prokaryotes, driving rapid genomic adaptation, the dissemination of antimicrobial resistance (AMR), and the acquisition of pathogenic virulence factors. Foreign genetic elements introduced into a host bacterial genome often retain the distinct nucleotide composition and codon preferences of their donor species. 

This project developed an automated, end-to-end Python bioinformatics pipeline designed to scan bacterial genomes using sliding window analysis of local GC-content and Codon Usage Bias (CUB) dynamics. Applied to the uropathogenic strain *Escherichia coli* CFT073, the pipeline evaluated 5,227 overlapping genomic windows (5,000 bp window size, 1,000 bp step size) and successfully identified 276 candidate HGT regions meeting a statistical deviation threshold of |Z| >= 2.0.

---

## 2. Biological Background & Methodology

### 2.1 Genomic Composition and HGT Signatures
Bacterial species maintain characteristic average GC-content and codon usage profiles shaped by mutation bias and translational selection over long evolutionary timescales. When foreign DNA (e.g., genomic islands, prophages, or integrative conjugative elements) integrates into a recipient genome, its composition contrasts sharply with the host background. Over evolutionary time, this signal decays through "amelioration," making compositional filtering an effective technique for identifying relatively recent transfer events.

### 2.2 Algorithmic Workflow
1. **Data Acquisition:** The reference FASTA sequence for *E. coli* CFT073 (5.23 Mb) is automatically retrieved from NCBI Entrez APIs.
2. **Sliding-Window Filtering (`cub_filter.py`):** The sequence is traversed with a 5,000 bp window at 1,000 bp step intervals. Local GC-content percentage is calculated as:
   $$	ext{GC \%} = rac{n_G + n_C}{N_{total}} 	imes 100$$
3. **Statistical Anomaly Detection (`hgt_scanner.py`):** $Z$-scores are computed across all genomic windows relative to the host mean ($\mu$) and standard deviation ($\sigma$):
   $$Z = rac{x_i - \mu}{\sigma}$$
   Windows exhibiting $|Z| \ge 2.0$ are flagged as anomalous candidate islands.
4. **Data Export & Visualization (`plot_islands.py`):** Flagged genomic coordinates are formatted into standard BED files for genome browser inspection, and a genome-wide track plot is produced via Seaborn/Matplotlib.

---

## 3. Results & Findings

- **Total Genomic Windows Analyzed:** 5,227
- **Mean Host GC Content:** 50.5%
- **Identified HGT Candidate Windows:** 276 windows (|Z| >= 2.0)
- **Output Artifacts:**
  - Standard BED File: `data/processed/predicted_hgt_islands.bed`
  - Genome Track Plot: `docs/genome_gc_hgt_track.png`

The identified regions highlight prominent GC-content fluctuations across the *E. coli* CFT073 chromosome, consistent with known pathogenicity islands (PAIs) and prophage integrations typical of uropathogenic strains.

---

## 4. Software Architecture & Quality Assurance

The codebase was architected around modular Python practices and modern devops standards:
- **`src/hgt_cub/`**: Reusable Python package containing core filtering and scanning routines.
- **`tests/test_pipeline.py`**: Automated unit tests using `pytest` covering GC calculation edge cases, sliding window boundary limits, Z-score thresholds, and BED file output validation.
- **Version Control & CI Readiness**: Excluded compiled bytecode (`__pycache__`) and local raw sequence datasets via `.gitignore` to maintain a clean repository structure.

---

## 5. Limitations & Future Directions

While GC-content sliding window filtering provides a fast, unsupervised method for genome-wide screening, compositional methods have inherent limitations:
1. **Ameliorated Islands:** Ancient transfer events that have ameliorated to match the host GC composition will not be flagged.
2. **Same-Composition Donors:** Transfers between organisms with identical baseline GC ratios remain undetectable by nucleotide composition alone.

### Proposed Future Enhancements:
- **Codon Adaptation Index (CAI) & Relative Synonymous Codon Usage (RSCU):** Incorporate full codon frequency matrix analysis beyond GC content.
- **Comparative Phylogenomics:** Integrate BLAST-based ortholog profiling across related strains to validate candidate islands via sequence conservation.
- **Automated Annotation Integration:** Cross-reference BED coordinates with GFF3 annotation files to automatically output candidate HGT gene lists (e.g., transposases, integrases, resistance markers).

---

## 6. Conclusion

This project demonstrates an effective computational screening framework for identifying candidate genomic islands in bacterial pathogens. The software is modular, fully unit-tested, reproducible, and ready for integration into broader comparative genomics pipelines.
