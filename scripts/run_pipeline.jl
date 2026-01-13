#!/usr/bin/env julia
"""
STD-NMR Analysis Pipeline
=========================

This script orchestrates the complete analysis workflow for STD-NMR saturation mapping.

Usage:
    julia run_pipeline.jl

Steps:
    1. Load configuration from config/hewl_settings.yaml
    2. Assign peaks using Reciprocal Best Hit (RBH) algorithm
    3. Generate publication-quality figures with AICc model selection
    4. Export PyMOL visualization scripts
    5. Save results to results/ directory

Author: Abdullah Yasir
License: MIT
"""

using YAML

# Load configuration
config = YAML.load_file("../config/hewl_settings.yaml")

println("═"^70)
println("STD-NMR Saturation Mapping Pipeline")
println("═"^70)
println()

# Step 1: Peak Assignment using RBH
println("Step 1/4: Running Reciprocal Best Hit peak assignment...")
include("../src/RBHAssign.jl")
println("✓ Peak assignment complete")
println()

# Step 2: Model Selection with AICc
println("Step 2/4: Generating publication figures with AICc model selection...")
include("../src/ModelSelector.jl")
println("✓ Figure generation complete")
println()

# Step 3: PyMOL Script Generation
println("Step 3/4: Generating PyMOL visualization scripts...")
run(`python ../src/PyMOLGenerator.py`)
println("✓ PyMOL scripts generated")
println()

# Step 4: Summary
println("Step 4/4: Pipeline complete!")
println()
println("Output locations:")
println("  • Peak assignments: data/processed/master-peak-list-rbh-optimized.tsv")
println("  • Figures: results/figures/")
println("  • PyMOL scenes: results/pymol_scenes/")
println()
println("═"^70)
