# STD-NMR Saturation Mapping Framework

**Automated pipeline for atomic-resolution saturation transfer analysis in fragment-based drug discovery**

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?style=flat&logo=julia)](https://julialang.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat&logo=python)](https://python.org/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

This framework automates the quantitative analysis of Saturation Transfer Difference (STD) NMR spectroscopy, mapping saturation transfer at atomic resolution across protein 3D structures. The pipeline identifies systematic detection biases and optimizes experimental parameters for uniform coverage in fragment screening campaigns.

**Key Innovation:** Application of the Reciprocal Best Hit (RBH) algorithm—originally from bioinformatics—to NMR peak assignment, providing unbiased, high-confidence chemical shift matching with automated ambiguity detection.

### Core Features

- **Bioinformatics-Inspired Assignment:** RBH algorithm ensures mutual best-match verification between experimental and reference peak lists
- **Information-Theoretic Model Selection:** AICc-based kinetic fitting distinguishes mono-exponential from complex saturation buildup
- **PyMOL Integration:** Automated 3D visualization pipeline for publication-quality structural figures
- **Modular Architecture:** Clean separation between analysis logic, visualization, and configuration

### Demonstration: HEWL Case Study

The framework was validated using Hen Egg White Lysozyme (HEWL) natural abundance ¹³C data as a stress test:
- 8-fold kinetic heterogeneity detected (k_sat: 0.563–4.635 s⁻¹)
- Frequency-dependent blind spots identified in aromatic binding regions
- Quantitative pulse shape optimization for detection uniformity

**Note:** This demonstration used challenging low-SNR natural abundance data. Performance is significantly enhanced with isotopically labeled samples (recommended).

---

## Repository Structure

```
├── src/
│   ├── RBHAssign.jl          # Reciprocal Best Hit peak assignment algorithm
│   ├── ModelSelector.jl      # AICc-based kinetic model selection
│   ├── PyMOLGenerator.py     # 3D molecular visualization automation
│   └── NMRAnalysisUtils.jl   # Core NMR processing utilities
├── scripts/
│   └── run_pipeline.jl       # Master automation script
├── config/
│   └── hewl_settings.yaml    # Configuration (chemical shifts, tolerances, thresholds)
├── data/processed/
│   └── master-peak-list-rbh-optimized.tsv
├── results/
│   └── master-peak-list-rbh-optimized.tsv  # Sample output
└── Project.toml              # Julia environment specification
```

---

## Quick Start

### Prerequisites

- **Julia** ≥ 1.9
- **Python** ≥ 3.8 with PyMOL

### Installation

```bash
git clone https://github.com/yasirabdullah123/Automated-framework-for-STD-NMR-saturation-mapping.git
cd Automated-framework-for-STD-NMR-saturation-mapping

# Install Julia dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Configuration

All parameters are centralized in `config/hewl_settings.yaml`:

- **Chemical shift corrections:** System-specific offsets (¹H: -0.11 ppm, ¹³C: -3.05 ppm for HEWL)
- **Uncertainty thresholds:** 5% (isotopically labeled) or 25% (natural abundance)
- **Assignment tolerances:** ¹H (0.1 ppm), ¹³C (1.0 ppm)
- **Ambiguity detection:** Flags assignments where second-best match is within 20% of best hit

Adapt this file for your protein system and experimental conditions.

### Running the Pipeline

```bash
# Automated full pipeline
julia scripts/run_pipeline.jl

# Or run individual modules
julia src/RBHAssign.jl          # Peak assignment
julia src/ModelSelector.jl      # Kinetic analysis with AICc
python src/PyMOLGenerator.py    # 3D visualization
```

---

## Core Algorithms

### Reciprocal Best Hit (RBH) Assignment

**Implementation:** `src/RBHAssign.jl`

Adapted from ortholog identification in comparative genomics, this algorithm ensures mutual best-match verification:

1. **Forward assignment:** Each experimental peak finds its nearest reference match
2. **Reverse verification:** Each reference peak finds its nearest experimental match
3. **RBH criterion:** Assignments confirmed only when both directions agree

**Distance metric:** Normalized Euclidean distance accounting for dimension-specific tolerances:
```
d = sqrt((Δδ_H / 0.1 ppm)² + (Δδ_C / 1.0 ppm)²)
```

**Ambiguity detection:** Flags assignments where the second-best match is within 20% of the best match, indicating potential spectral overlap.

**Why this matters:** Unlike greedy assignment algorithms, RBH prevents false-positive assignments in crowded spectral regions, ensuring high-confidence peak identification.

### AICc-Based Model Selection

**Implementation:** `src/ModelSelector.jl`

Uses Akaike Information Criterion (corrected for small sample sizes) to select the optimal kinetic model:

- **Mono-exponential:** `I(t) = I₀(1 - exp(-k·t))`
- **Stretched exponential:** `I(t) = I₀(1 - exp(-(k·t)^β))`

The algorithm penalizes model complexity, preventing overfitting while identifying genuine multi-exponential behavior indicative of complex saturation pathways.

**Quality control:** Rejects fits with R² < 0.85 or parameter uncertainty > threshold (5% labeled, 25% natural abundance).

## Technical Highlights

This repository demonstrates:

- **Cross-disciplinary algorithm adaptation:** RBH from bioinformatics to NMR
- **Statistical rigor:** Information-theoretic model selection (AICc)
- **API integration:** Python-PyMOL bridge for automated 3D visualization
- **Modular design:** Clean separation of concerns with centralized configuration
- **Reproducibility:** Full dependency specification via `Project.toml`

---

## Use Cases

- **Method development:** Benchmark new assignment algorithms against RBH
- **Parameter optimization:** Adapt pipeline for custom saturation schemes
- **Quality control:** Detect systematic biases in STD-NMR experiments
- **Education:** Learn statistical model selection in biophysical data analysis

---

## Author

**Abdullah Yasir**
MSc Drug Discovery and Development
University College London

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- Dr. Christopher Waudby (UCL) - NMR methodology guidance
- BMRB Entry 4562 - Reference chemical shift database
