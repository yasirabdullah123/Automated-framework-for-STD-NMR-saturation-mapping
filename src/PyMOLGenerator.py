#!/usr/bin/env python3
"""
STD-NMR Parameter Optimization Visualization Generator

This script creates PyMOL visualizations for optimizing STD-NMR parameters
based on systematic saturation transfer experiments with lysozyme.

The goal is to determine optimal conditions for STD-NMR experiments where
saturation at 1 ppm must efficiently propagate to aromatic regions (6-8 ppm)
for reliable ligand binding detection.

Experiment Design:
- Non-isotopically labeled lysozyme (10mM)
- 1H-13C HSQC for atomic-level precision
- Three parameter optimization studies:
  1. Saturation frequency (19 frequencies: -0.5 to 8.5 ppm)
  2. Saturation time (10 time points: 0.1 to 5.0 seconds)
  3. Pulse parameters (duration, shape, solvent: H2O vs D2O)

Output: Organized PyMOL scripts with scenes for parameter comparison
"""

import os
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np

def sanitize_atom_name(name: str) -> str:
    """
    Correctly clean atom names for PyMOL compatibility by removing only truly
    problematic characters like quotes or asterisks. Standard IUPAC atom names
    are left intact.
    """
    # Standard atom names like HD1, CG2, HDy, etc., are valid in PyMOL.
    # We only need to remove characters that can break the selection syntax.
    sanitized = name.replace("'", "").replace("*", "")
    return sanitized

def read_saturation_data(filepath: Path) -> list[tuple[str, str, float]]:
    """Read saturation data from pymol files"""
    data = []
    try:
        with filepath.open('r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    resi, atom, ratio = parts[0], parts[1], float(parts[2])
                    data.append((resi, atom, ratio))
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return data

def generate_std_nmr_master_script(directory_path: Path, output_dir: Path, experiment_type: str):
    """Generate comprehensive STD-NMR optimization script"""
    
    # Get all pymol files, but exclude old 283K format files
    all_pymol_files = sorted(directory_path.glob("pymol_*.txt"))
    # Filter out old format 283K files (pymol_283K_*.txt) to use the new format
    pymol_files = [f for f in all_pymol_files if not f.stem.startswith("pymol_283K_")]
    
    if not pymol_files:
        print(f"No pymol files found in {directory_path}")
        return

    # Create experiment-specific script
    script_name = f"std_nmr_{experiment_type}_optimization.pml"
    script_path = output_dir / script_name
    
    print(f"Generating STD-NMR optimization script: {script_name}")
    print(f"Experiment type: {experiment_type}")
    print(f"Processing {len(pymol_files)} conditions...")

    # Script header with experiment context
    commands = [
        "# STD-NMR Parameter Optimization Visualization",
        f"# Experiment: {experiment_type.replace('_', ' ').title()}",
        "# ",
        "# BACKGROUND:",
        "# Non-isotopically labeled lysozyme (10mM) with 1H-13C HSQC",
        "# Goal: Optimize STD-NMR parameters for ligand binding detection",
        "# Critical question: Does saturation at 1 ppm reach aromatic regions?",
        "# ",
        "# VISUALIZATION SCHEME:",
        "# Bright Yellow = Strong saturation (I/I0 < 0.7) - Optimal transfer",
        "# Orange = Moderate saturation (0.7 ≤ I/I0 < 0.85) - Good transfer", 
        "# Red = Weak saturation (0.85 ≤ I/I0 < 0.95) - Limited transfer",
        "# Dark Purple = No saturation (I/I0 ≥ 0.95) - Poor transfer",
        "# ",
        f"# Load with: @{script_name}",
        "",
        "# === INITIAL SETUP ===",
        "fetch 2vb1, async=0",
        "hide everything",
        "show cartoon", 
        "select backbone",
        "show lines, backbone",
        "set line_width, 1.5",
        "set cartoon_transparency, 0.6",
        "set sphere_scale, 0.4",
        "set sphere_transparency, 0.0",
        "bg_color white",
        "",
        "# Define STD-NMR color scheme (reverse plasma)",
        "set_color std_strong,   [0.992, 0.906, 0.145]  # Bright yellow - optimal",
        "set_color std_moderate, [0.988, 0.553, 0.349]  # Orange - good", 
        "set_color std_weak,     [0.843, 0.188, 0.153]  # Red - limited",
        "set_color std_none,     [0.267, 0.004, 0.329]  # Dark purple - poor",
        "",
        "# Color protein backbone gray",
        "color gray80, all",
        "",
    ]

    # Add experiment-specific context
    if experiment_type == "frequency_analysis":
        commands.extend([
            "# === FREQUENCY ANALYSIS CONTEXT ===",
            "# Testing 19 saturation frequencies (-0.5 to 8.5 ppm)",
            "# Key regions:",
            "#   Aliphatic: -0.5 to 2.0 ppm (typical STD saturation)",
            "#   Alpha-CH: 2.0 to 4.0 ppm (protein backbone)",
            "#   Aromatic: 6.0 to 9.0 ppm (ligand binding sites)",
            "# Question: Which frequencies provide best coverage?",
            "",
        ])
    elif experiment_type == "time_analysis":
        commands.extend([
            "# === SATURATION TIME ANALYSIS CONTEXT ===", 
            "# Testing 10 saturation times (0.1 to 5.0 seconds)",
            "# Shorter times: More selective saturation",
            "# Longer times: More complete saturation",
            "# Question: What's the optimal time for full protein coverage?",
            "",
        ])
    elif experiment_type == "pulse_analysis":
        commands.extend([
            "# === PULSE PARAMETER ANALYSIS CONTEXT ===",
            "# Testing different pulse conditions:",
            "#   Duration: 2.5ms, 5ms, 20ms, 50ms Gaussian",
            "#   Shape: Gaussian vs EBURP1",
            "#   Solvent: H2O vs D2O effects",
            "# Question: Which pulse parameters optimize saturation transfer?",
            "",
        ])
    elif experiment_type == "temperature_analysis_283K":
        commands.extend([
            "# === TEMPERATURE ANALYSIS CONTEXT (283K) ===",
            "# Temperature-dependent STD-NMR saturation transfer analysis",
            "# Testing 10 saturation times (0.1 to 5.0 seconds) at 283K (10°C)",
            "# On-resonance: 0.7824 ppm (temperature-adjusted from 1.0 ppm)",
            "# Off-resonance: 30 ppm control",
            "# Combined spectra: 57+59 (on-res) vs 56+58 (off-res)",
            "# Question: How does lower temperature affect saturation kinetics?",
            "#   Expected effects: Slower molecular tumbling, enhanced NOE",
            "#                     Slower exchange kinetics, different dynamics",
            "",
        ])
    elif experiment_type == "temperature_analysis_283K_fixed":
        commands.extend([
            "# === TEMPERATURE ANALYSIS CONTEXT (283K FIXED) ===",
            "# Temperature-dependent STD-NMR saturation transfer analysis (Fixed version)",
            "# Testing 10 saturation times (0.1 to 5.0 seconds) at 283K (10°C)",
            "# On-resonance: 0.7824 ppm (temperature-adjusted from 1.0 ppm)",
            "# Off-resonance: 30 ppm control",
            "# Combined spectra: 57+59 (on-res) vs 56+58 (off-res)",
            "# Improvements: Enhanced noise normalization, better error handling",
            "# Question: How does lower temperature affect saturation kinetics?",
            "#   Expected effects: Slower molecular tumbling, enhanced NOE",
            "#                     Slower exchange kinetics, different dynamics",
            "",
        ])

    # Process each condition
    all_scenes = []
    condition_data = {}
    
    for i, filepath in enumerate(pymol_files):
        condition_name = filepath.stem.replace("pymol_", "")
        scene_name = f"scene_{condition_name.replace('.', '_').replace('-', 'neg')}"
        all_scenes.append((scene_name, condition_name))
        
        print(f"  Processing condition {i+1}/{len(pymol_files)}: {condition_name}")
        
        data = read_saturation_data(filepath)
        if not data:
            print(f"    Warning: No data in {filepath.name}")
            continue
        
        condition_data[condition_name] = data
        
        # Analyze saturation efficiency
        ratios = [ratio for _, _, ratio in data]
        strong_count = sum(1 for r in ratios if r < 0.7)
        moderate_count = sum(1 for r in ratios if 0.7 <= r < 0.85)
        weak_count = sum(1 for r in ratios if 0.85 <= r < 0.95)
        none_count = sum(1 for r in ratios if r >= 0.95)
        
        commands.extend([
            f"# === CONDITION: {condition_name} ===",
            f"# Saturation efficiency: {strong_count} strong, {moderate_count} moderate,",
            f"#                        {weak_count} weak, {none_count} none",
            "",
        ])
        
        # Reset for new scene
        commands.extend([
            "color gray80, all",
            "hide spheres, all",
            "",
        ])
        
        # Categorize atoms by saturation strength
        atom_categories = defaultdict(list)
        residue_coverage = defaultdict(float)
        
        for resi, atom, ratio in data:
            sanitized_atom = sanitize_atom_name(atom)
            if sanitized_atom in ['CA', 'CB', 'CG', 'CD', 'CE', 'CZ', 'CH', 'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']:
                atom_sel = f"(resi {resi} and name {sanitized_atom})"
                
                # Update best saturation for this residue
                if resi not in residue_coverage or ratio < residue_coverage[resi]:
                    residue_coverage[resi] = ratio
                
                # Categorize by saturation strength
                if ratio < 0.7:
                    atom_categories['strong'].append(atom_sel)
                elif ratio < 0.85:
                    atom_categories['moderate'].append(atom_sel)
                elif ratio < 0.95:
                    atom_categories['weak'].append(atom_sel)
                else:
                    atom_categories['none'].append(atom_sel)
        
        # Create a single selection containing ALL atoms that will be shown as spheres
        all_visible_atoms = []
        for level in ['strong', 'moderate', 'weak', 'none']:
            all_visible_atoms.extend(atom_categories[level])
        
        if all_visible_atoms:
            # Important: Name the selection specifically for this scene
            scene_specific_selection_name = f"saturated_C_{scene_name}"
            full_selection_str = " or ".join(all_visible_atoms)
            commands.append("# Define a selection for all C atoms with saturation data in this scene")
            commands.append(f"select {scene_specific_selection_name}, {full_selection_str}")
            commands.append("")
        
        # Show and color atoms with chunking (atomic-level precision)
        for level, color in [('strong', 'std_strong'), ('moderate', 'std_moderate'), ('weak', 'std_weak'), ('none', 'std_none')]:
            atom_list = atom_categories[level]
            if atom_list:
                # Process in chunks to avoid command length limits
                chunk_size = 15
                for j in range(0, len(atom_list), chunk_size):
                    chunk = atom_list[j:j+chunk_size]
                    selection_str = " or ".join(chunk)
                    commands.append(f"show spheres, {selection_str}")
                    commands.append(f"color {color}, {selection_str}")
        
        # Highlight aromatic binding regions (structure only, no color override)
        commands.extend([
            "# Highlight aromatic binding regions (structure only)",
            "select aromatic_binding, resi 3+28+62+63+108+111+123",  # Phe/Tyr/Trp in active site
            "show sticks, aromatic_binding",
            "set stick_radius, 0.2, aromatic_binding",
            "set stick_transparency, 0.7, aromatic_binding",
            
            # COMPREHENSIVE ADDITION: Select ALL aromatic residues for complete analysis
            "# Select ALL aromatic residues on the protein for comprehensive analysis",
            "select all_aromatics, resn PHE+TYR+TRP+HIS",
            
            # Don't color override - let saturation colors show through
            "",
        ])
        
        # Add distance measurements from saturation center (unique per scene)
        commands.extend([
            "# Distance analysis (unique per scene)",
            "select sat_center, resi 1-10",  # N-terminal saturation region
            f"distance dist_{scene_name}, sat_center, aromatic_binding",
            f"hide labels, dist_{scene_name}",
            f"color yellow, dist_{scene_name}",
            "",
        ])
        
        # Add efficiency information (positioned off-protein)
        commands.extend([
            "# Add efficiency information (positioned off-protein)",
            f"pseudoatom efficiency_info_{scene_name}, pos=[30, 30, 30]",
            f"label efficiency_info_{scene_name}, 'Strong: {strong_count}, Moderate: {moderate_count}, Weak: {weak_count}, None: {none_count}'",
            f"hide spheres, efficiency_info_{scene_name}",
            "",
        ])
        
        # Store the scene
        commands.extend([
            f"scene {scene_name}, store",
            f"print 'Scene {scene_name}: {condition_name} - Coverage: {len(residue_coverage)} residues'",
            "",
        ])
        
        # Clean up scene-specific objects for next iteration
        commands.extend([
            "# Clean up scene-specific objects",
            f"delete dist_{scene_name}",
            f"delete efficiency_info_{scene_name}",
            "",
        ])

    # Add analysis commands
    commands.extend([
        "# === ANALYSIS COMMANDS ===",
        "# Use these commands to analyze saturation patterns:",
        "",
        "# Show only strongly saturated regions:",
        "# hide spheres, all",
        "# show spheres, all and std_strong",
        "",
        "# Highlight aromatic regions (key for ligand binding):",
        "# select aromatic_residues, resi 62+108+111+123  # Trp/Phe/Tyr",
        "# show sticks, aromatic_residues",
        "",
        "# Comprehensive aromatic analysis:",
        "# show sticks, all_aromatics",
        "# color red, all_aromatics and std_none  # Find aromatic 'blind spots'",
        "# color green, all_aromatics and std_strong  # Find efficiently saturated aromatics",
        "",
        "# Analyze saturation transfer efficiency:",
        "# select strong_sat, std_strong",
        "# select moderate_sat, std_moderate", 
        "# select weak_sat, std_weak",
        "# select no_sat, std_none",
        "",
        "# Distance analysis from saturation center:",
        "# show distance, dist_scene_name",
        "# hide distance, dist_scene_name",
        "",
        "# Compare saturation at different frequencies:",
        f"# scene {all_scenes[0][0] if all_scenes else 'scene_first'}",
        "",
        "# Quantitative analysis:",
        "# count_atoms strong_sat",
        "# count_atoms moderate_sat", 
        "# count_atoms weak_sat",
        "# count_atoms no_sat",
        "",
    ])

    # Add navigation help
    commands.extend([
        "# === SCENE NAVIGATION ===",
        f"# Total scenes: {len(all_scenes)}",
        "",
    ])
    
    for scene_name, condition_name in all_scenes:
        commands.append(f"# {scene_name} -> {condition_name}")
    
    # Add efficiency summary
    commands.extend([
        "",
        "# === EFFICIENCY SUMMARY ===",
        "# Overall saturation transfer efficiency across all conditions:",
        "",
    ])
    
    # Calculate overall statistics
    total_strong = sum(1 for data in condition_data.values() 
                      for _, _, ratio in data if ratio < 0.7)
    total_moderate = sum(1 for data in condition_data.values() 
                        for _, _, ratio in data if 0.7 <= ratio < 0.85)
    total_weak = sum(1 for data in condition_data.values() 
                    for _, _, ratio in data if 0.85 <= ratio < 0.95)
    total_none = sum(1 for data in condition_data.values() 
                    for _, _, ratio in data if ratio >= 0.95)
    total_atoms = total_strong + total_moderate + total_weak + total_none
    
    commands.extend([
        f"# Total atoms analyzed: {total_atoms}",
        f"# Strong saturation (< 0.7): {total_strong} ({total_strong/total_atoms*100:.1f}%)",
        f"# Moderate saturation (0.7-0.85): {total_moderate} ({total_moderate/total_atoms*100:.1f}%)",
        f"# Weak saturation (0.85-0.95): {total_weak} ({total_weak/total_atoms*100:.1f}%)",
        f"# No saturation (≥ 0.95): {total_none} ({total_none/total_atoms*100:.1f}%)",
        "",
        "# === FINAL SETUP ===",
        "orient",
        "zoom visible",
        "set scene_buttons, 1",
        "",
        f"print 'STD-NMR {experiment_type} optimization loaded!'",
        f"print 'Scenes: {len(all_scenes)} conditions'",
        "print 'Goal: Find optimal parameters for 1 ppm -> aromatic region transfer'",
        "print 'Use scene buttons or type: scene scene_name'",
        "",
    ])

    # Write the script
    with script_path.open('w') as f:
        f.write("\n".join(commands))
    
    print(f"✅ Created: {script_path}")
    return script_path

def main():
    """Main function with improved organization"""
    parser = argparse.ArgumentParser(
        description="STD-NMR Parameter Optimization Visualization Generator",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'directories',
        nargs='*', 
        default=["./sattimeanalysis", "./pulse_shape_analysis", "./31", "./sattimeanalysis_283K", "./sattimeanalysis_283K_fixed"],
        help="Directories to process"
    )
    parser.add_argument(
        '-o', '--output-dir',
        default="std_nmr_optimization_scenes",
        help="Output directory for organized scripts"
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("STD-NMR PARAMETER OPTIMIZATION VISUALIZATION")
    print("=" * 70)
    print("Goal: Optimize parameters for ligand binding detection")
    print("Critical: Ensure saturation from 1 ppm reaches aromatic regions")
    print("=" * 70)
    
    # Map directories to experiment types
    experiment_mapping = {
        "sattimeanalysis": "time_analysis",
        "pulse_shape_analysis": "pulse_analysis", 
        "31": "frequency_analysis",
        "sattimeanalysis_283K": "temperature_analysis_283K",
        "sattimeanalysis_283K_fixed": "temperature_analysis_283K_fixed"
    }
    
    generated_scripts = []
    
    for directory_str in args.directories:
        directory = Path(directory_str)
        if directory.exists() and directory.is_dir():
            exp_type = experiment_mapping.get(directory.name, directory.name)
            script_path = generate_std_nmr_master_script(directory, output_path, exp_type)
            if script_path:
                generated_scripts.append(script_path)
        else:
            print(f"⚠️  Directory not found: {directory}")
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Generated {len(generated_scripts)} optimization scripts:")
    for script in generated_scripts:
        print(f"  📄 {script.name}")
    print(f"\nAll scripts saved to: {output_path}")
    print("\nUsage in PyMOL:")
    print("  @std_nmr_frequency_analysis_optimization.pml")
    print("  @std_nmr_time_analysis_optimization.pml") 
    print("  @std_nmr_pulse_analysis_optimization.pml")
    print("=" * 70)

if __name__ == "__main__":
    main()