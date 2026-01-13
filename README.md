# NMR-STD-Analysis: Automated Framework for Atomic-Resolution Saturation Mapping in STD-NMR

**A modular Julia/Python pipeline for high-throughput saturation transfer analysis and parameter optimization**

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?style=flat&logo=julia)](https://julialang.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat&logo=python)](https://python.org/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

This repository provides a high-performance framework in Julia and Python to automate the quantitative analysis of Saturation Transfer Difference (STD) NMR spectroscopy. It enables researchers to map saturation transfer at atomic resolution across 3D protein structures, identifying systematic detection blind spots in fragment-based drug discovery.

The pipeline was developed using Hen Egg White Lysozyme (HEWL) as a model system, demonstrating robust performance even in challenging low-SNR environments (natural abundance ¹³C).

### Key Features

- **Reciprocal Best Hit (RBH) Assignment:** Bioinformatics-inspired algorithm for unbiased, high-confidence peak assignment with ambiguity detection
- **Statistical Model Selection:** Automated kinetic fitting using Akaike Information Criterion (AICc) to discriminate between mono-exponential and complex buildup kinetics
- **3D Structure Mapping:** Python-based PyMOL API integration for automated generation of publication-ready molecular visualizations
- **Parameter Optimization:** Systematic evaluation of saturation time, RF pulse shape, and frequency offset to maximize detection uniformity

### Case Study: HEWL Natural Abundance Analysis

Using a natural abundance ¹³C sample as a stress test, the pipeline demonstrated:
- 8-fold variation in saturation kinetics across protein sites (k_sat: 0.563–4.635 s⁻¹)
- Identification of frequency-dependent blind spots affecting aromatic binding regions
- Quantitative pulse shape performance comparison for uniformity optimization

**Note:** While this case study uses low-SNR natural abundance data, the framework is optimized for isotopically labeled samples and achieves maximum performance with ¹³C-labeled proteins.

---

## Repository Structure

```
├── src/                                    # Source code
│   ├── assignment/
│   │   └── peak_assignment_rbh.jl         # Reciprocal Best Hit peak assignment algorithm
│   ├── analysis/
│   │   ├── saturation_time_analysis.jl    # Temporal kinetics optimization (0.1-5.0 s)
│   │   ├── pulse_shape_analysis.jl        # Gaussian/E-BURP1/I-BURP2 comparison
│   │   └── frequency_offset_analysis.jl   # Blind spot mapping (-0.5 to 8.5 ppm)
│   ├── visualization/
│   │   ├── generate_publication_figures.jl # Master figure generation pipeline
│   │   ├── generate_frequency_heatmaps.jl  # Frequency-dependent saturation maps
│   │   ├── generate_pulse_heatmaps.jl      # Pulse shape comparison visualizations
│   │   └── generate_time_heatmaps.jl       # Temporal saturation kinetics maps
│   └── utils/
│       └── nmr_analysis_utils.jl           # Core NMR processing utilities
├── data/
│   └── *.tsv                               # Peak assignments and master lists
├── output/
│   ├── figures/                            # Publication-ready figures
│   └── pymol_scripts/
│       └── generate_std_nmr_optimizer.py   # 3D molecular visualization automation
└── docs/                                   # Documentation and manuscripts
```

---

## Quick Start

### Prerequisites

- **Julia** ≥ 1.9 with packages: `NMRTools`, `Plots`, `StatsPlots`, `LsqFit`, `Measurements`
- **Python** ≥ 3.8 with `PyMOL` for 3D visualization
- **ImageMagick** for figure processing

### Installation

```bash
git clone https://github.com/yasirabdullah123/An-HSQC-Based-Investigation-of-Saturation-Transfer-Heterogeneity-in-STD-NMR.git
cd An-HSQC-Based-Investigation-of-Saturation-Transfer-Heterogeneity-in-STD-NMR

# Install Julia dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Configuration

The pipeline uses `config/settings.yaml` for all configurable parameters:

- **Chemical shift corrections:** Adjust for systematic offsets between experimental and reference data
- **Uncertainty thresholds:** Default 5% for labeled samples, 25% for natural abundance
- **Assignment tolerances:** ¹H (0.1 ppm) and ¹³C (1.0 ppm) matching windows
- **Ambiguity detection:** Flags assignments where 2nd-best hit is within 20% of best

Edit this file to adapt the pipeline to your specific protein system and experimental conditions.

### Running Analyses

```bash
# Parameter optimization studies
cd src/analysis
julia saturation_time_analysis.jl      # saturation time analysis 
julia pulse_shape_analysis.jl          # Pulse shape analysis 
julia frequency_offset_analysis.jl     # saturation frequency analysis 

# Generate figures
cd ../visualization
julia generate_publication_figures.jl

# 3D molecular visualizations scripts
cd ../../output/pymol_scripts
python generate_std_nmr_optimizer.py
```

---

## Key Algorithms

### 1. Reciprocal Best Hit (RBH) Peak Assignment

**Location:** `src/assignment/peak_assignment_rbh.jl`

Novel bidirectional verification algorithm for robust HSQC peak assignment:
- Bidirectional verification ensuring >95% reliability
- Normalized Euclidean distance: `d = sqrt((Δδ_H/0.1)² + (Δδ_C/1.0)²)`
- Systematic chemical shift corrections for temperature/experimental conditions

### 2. Δ-Σ Diagnostic Framework


Quantitative blind spot identification:
- **Σ** = Total saturation accessibility (aliphatic + aromatic)
- **Δ** = Frequency-dependent bias (aliphatic - aromatic)
- **Critical threshold:** Σ ≥ 0.8 AND Δ < -0.15

### 3. Multi-Exponential Kinetic Fitting


Statistical model selection via AICc:
- Mono-exponential: `I(t) = I₀(1 - e^(-kt))`
- Stretched exponential: `I(t) = I₀(1 - e^(-(kt)^β))`
- Quality filters: R² > 0.85, uncertainty < 25%

---

## Data

### Raw NMR Experiments (`data/raw_experiments/`)

** systematically designed Bruker experiments:**
- **Temporal series:** 10 time points (0.1-5.0 s) → kinetic heterogeneity quantification
- **Frequency sweep:** 19 offsets (-0.5 to 8.5 ppm) → blind spot mapping
- **Pulse shapes:** Gaussian/E-BURP1/I-BURP2 → uniformity optimization
- **Temperature study:** 283K, 293K, 303K → thermal effects on spin diffusion


---

## Publication Figures

All figures generated by `src/visualization/generate_publication_figures.jl`:

- **Figure 1:** Global Kinetic Portrait (3-panel composite)
- **Figure 2:** Δ-Σ Diagnostic Plot (6 critical blind spots identified)
- **Figure 3:** Frequency Response Profile (average saturation behavior)
- **Figure 4:** Pulse Shape Performance (split raincloud visualization)

**Output:** High-resolution SVG (300 DPI equivalent), standardized I/I₀ color mapping

---


## Author

**Abdullah Yasir**
MSc Drug Discovery and Development, University College London

---

## Acknowledgments

- **Dr. Christopher Waudby** (UCL) - Supervision and NMR expertise
- **BMRB Entry 4562** - Reference chemical shift database
- **NMR community** - Open-source tools and methodologies
