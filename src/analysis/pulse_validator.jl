using Plots
using NMRTools
using Measurements
using DelimitedFiles
using Printf
using Statistics
using Dates
using StatsPlots
using Clustering

# Include the centralized, verified alignment functions
include("../utils/nmr_analysis_utils.jl")

println("=== Pulse Shape Analysis: STD Saturation Transfer Comparison ===")
println("Comparing 4 different pulse shapes:")
println("Exp 34: 50ms Gaussian")
println("Exp 35: 5ms Gaussian") 
println("Exp 36: 20ms EBURP1")
println("Exp 41: 2.5ms Gaussian")

# Create output directory
output_dir = "pulse_shape_analysis"
if !isdir(output_dir)
    mkdir(output_dir)
    println("Created output directory: $(output_dir)")
end

# --- Use Centralized Peak Alignment ---
aligned_peaks = NMRAnalysisUtils.load_and_align_peaks("../data/master-peak-list-final.tsv", "../data/raw_experiments/33", x_shift=0.0, y_shift=0.0)
x = aligned_peaks.x
y = aligned_peaks.y
resi = aligned_peaks.resi
resn = aligned_peaks.resn
atomH = aligned_peaks.atomH
atomC = aligned_peaks.atomC

println("Loaded and aligned $(length(x)) peaks for analysis")

# --- ALIGNMENT VERIFICATION ---
println("
=== Verifying Peak Alignment ===")

# Load experiment 34 for plotting (first pulse shape experiment)
spec = loadnmr("../data/raw_experiments/34")
firstplane = spec[:,:,1]

# Simple alignment check plot (like prepare-peak-list.jl)
plot(firstplane)
scatter!(x, y, xflip=true, yflip=true, 
    xlabel="¹H Chemical Shift (ppm)", ylabel="¹³C Chemical Shift (ppm)", 
    title="Peak List Alignment Check (Exp 34)", legend=false, grid=false)
savefig(joinpath(output_dir, "alignment_check.svg"))
println("Alignment check plot saved to: $(joinpath(output_dir, "alignment_check.svg"))")

# Define experiments and their pulse shapes
experiments = [
    (exp="../data/raw_experiments/34", name="50ms Gaussian", duration="50ms", shape="Gaussian"),
    (exp="../data/raw_experiments/35", name="5ms Gaussian", duration="5ms", shape="Gaussian"),
    (exp="../data/raw_experiments/36", name="20ms EBURP1 ", duration="20ms", shape="EBURP1 "),
    (exp="../data/raw_experiments/41", name="2.5ms Gaussian", duration="2.5ms", shape="Gaussian")
]           

println("\nAnalyzing $(length(experiments)) pulse shape experiments:")
for exp_info in experiments
    println("Exp $(exp_info.exp): $(exp_info.name)")
end

# Load and process all experiments
println("\nLoading and processing experiments...")
all_specs = []
all_ratios = []

for exp_info in experiments
    println("Processing Experiment $(exp_info.exp): $(exp_info.name)")
    
    # Load experiment
    spec = loadnmr(exp_info.exp)
    spec = spec / spec[:noise]
    spec[:noise] = 1.0
    
    println("Experiment $(exp_info.exp) dimensions: $(size(spec))")
    
    # Calculate I/I0 ratios using internal reference with uncertainty filtering
    println("Using 25% uncertainty filtering (like analyze_exp31.jl)")
    
    II0 = map(1:length(x)) do i
        # Get reference intensity (slice 1 - 30 ppm off-resonance)
        ref_intensity = spec[Near(x[i]), Near(y[i]), 1] ± spec[:noise]
        # Get on-resonance intensity (slice 2)
        sat_intensity = spec[Near(x[i]), Near(y[i]), 2] ± spec[:noise]
        
        # Filter out peaks with poor reference intensity
        ref_value = Measurements.value(ref_intensity)
        if ref_value <= 0.1  # Reference intensity too small or negative
            return NaN ± 0.0  # Set to NaN if reference is poor
        end
        
        # Calculate ratio with automatic error propagation
        ratio = sat_intensity / ref_intensity
        
        # Apply uncertainty cutoff for I/I0 ratios
        # Note: 25% threshold used here due to low SNR in natural abundance 13C data
        # For isotopically labeled samples, 5% is recommended
        if Measurements.uncertainty(ratio) > 0.25
            return NaN ± 0.0  # Set to NaN if uncertainty exceeds threshold
        else
            return ratio
        end
    end
    
    push!(all_specs, spec)
    push!(all_ratios, II0)
    
    println("Calculated $(length(II0)) peak ratios for $(exp_info.name)")
end

# Convert ratios to matrix (peaks × experiments)
II0_matrix = hcat(all_ratios...)

println("Intensity ratio matrix shape: $(size(II0_matrix))")

# Generate PyMOL files for each pulse shape
println("\nGenerating PyMOL coloring files...")
for (idx, exp_info) in enumerate(experiments)
    filename = joinpath(output_dir, "pymol_$(exp_info.exp)_$(replace(exp_info.name, " " => "_")).txt")
    
    open(filename, "w") do f
        println(f, "## PyMOL coloring data for $(exp_info.name)")
        println(f, "# Lower values indicate stronger saturation transfer")
        println(f, "# Reference: Internal (slice 1 - 30 ppm off-resonance)")
        println(f, "# Uncertainty filtering: 25% threshold applied")
        println(f, "# Format: resi atomC ratio")
        
        for i in 1:length(x)
            n = resi[i]
            a = atomC[i]
            saturation = Measurements.value(II0_matrix[i, idx])
            
            # Only output non-NaN values (from uncertainty filtering only)
            if !isnan(saturation)
                println(f, "$n $a $saturation")
            end
        end
    end
    println("Saved $(filename)")
end



# Generate summary statistics
println("\n=== Pulse Shape Effectiveness Summary ===")
for (idx, exp_info) in enumerate(experiments)
    ratios = Measurements.value.(II0_matrix[:, idx])
    valid_ratios = ratios[.!isnan.(ratios)]
    
    if length(valid_ratios) > 0
        mean_ratio = mean(valid_ratios)
        median_ratio = median(valid_ratios)
        strong_sat_count = sum(valid_ratios .< 0.8)  # Count peaks with strong saturation
        
        println("$(exp_info.name):")
        println("  Valid peaks: $(length(valid_ratios))/$(length(ratios))")
        println("  Mean I/I₀: $(round(mean_ratio, digits=3))")
        println("  Median I/I₀: $(round(median_ratio, digits=3))")
        println("  Peaks with strong saturation (I/I₀ < 0.8): $(strong_sat_count)/$(length(valid_ratios))")
        println("  Saturation efficiency: $(round(100*(1-mean_ratio), digits=1))%")
    else
        println("$(exp_info.name): No valid peaks detected")
    end
end

# Save summary data to file
summary_file = joinpath(output_dir, "pulse_shape_summary.txt")
open(summary_file, "w") do f
    println(f, "# Pulse Shape Analysis Summary")
    println(f, "# Generated: $(now())")
    println(f, "# Experiments analyzed: $(length(experiments))")
    println(f, "# Peaks analyzed: $(length(x))")
    println(f, "#")
    println(f, "# Format: experiment_name mean_ratio median_ratio strong_sat_count efficiency_percent")
    
    for (idx, exp_info) in enumerate(experiments)
        ratios = Measurements.value.(II0_matrix[:, idx])
        valid_ratios = ratios[.!isnan.(ratios)]
        
        if length(valid_ratios) > 0
            mean_ratio = mean(valid_ratios)
            median_ratio = median(valid_ratios)
            strong_sat_count = sum(valid_ratios .< 0.8)
            efficiency = round(100*(1-mean_ratio), digits=1)
            
            println(f, "$(replace(exp_info.name, " " => "_")) $(round(mean_ratio, digits=3)) $(round(median_ratio, digits=3)) $(strong_sat_count) $(efficiency)")
        else
            println(f, "$(replace(exp_info.name, " " => "_")) NaN NaN 0 NaN")
        end
    end
end
println("Summary statistics saved to: $(summary_file)")

# Generate comprehensive optimization recommendations
println("\n=== STD NMR Optimization Recommendations ===")

# Find the best performing pulse shape
pulse_performance = []
for (idx, exp_info) in enumerate(experiments)
    ratios = Measurements.value.(II0_matrix[:, idx])
    valid_ratios = ratios[.!isnan.(ratios)]
    
    if length(valid_ratios) > 0
        mean_ratio = mean(valid_ratios)
        efficiency = round(100*(1-mean_ratio), digits=1)
        push!(pulse_performance, (exp_info.name, mean_ratio, efficiency))
    end
end

if length(pulse_performance) > 0
    sort!(pulse_performance, by=x->x[2])  # Sort by mean ratio (lower is better)
    best_pulse = pulse_performance[1]
    
    println("🏆 **Best Overall Pulse Shape**: $(best_pulse[1])")
    println("   Mean I/I₀: $(round(best_pulse[2], digits=3))")
    println("   Saturation Efficiency: $(best_pulse[3])%")
    
    # Provide specific recommendations
    println("\n📋 **Optimization Recommendations for STD NMR:**")
    println("1. **Primary Choice**: Use $(best_pulse[1]) for general STD experiments")
    
    println("\n2. **Experimental Considerations:**")
    println("   • Ensure saturation frequency is set to 1 ppm for optimal coverage")
    println("   • Use saturation times of 2-3 seconds for full buildup")
    println("   • Monitor data quality with uncertainty filtering")
    println("   • Consider ligand binding site when choosing pulse parameters")
    
    println("\n3. **Quality Control:**")
    println("   • Check PyMOL visualization files for spatial distribution")
    println("   • Verify saturation reaches aromatic regions if ligand binds there")
    println("   • Monitor for saturation transfer artifacts")
end

# Save comprehensive recommendations to file
recommendations_file = joinpath(output_dir, "std_optimization_recommendations.txt")
open(recommendations_file, "w") do f
    println(f, "# STD NMR Optimization Recommendations")
    println(f, "# Generated: $(now())")
    println(f, "# Based on Lysozyme Pulse Shape Analysis")
    println(f, "#")
    println(f, "# Key Findings:")
    
    if length(pulse_performance) > 0
        println(f, "# Best Overall: $(pulse_performance[1][1])")
        println(f, "# Mean I/I₀: $(pulse_performance[1][2])")
        println(f, "# Efficiency: $(pulse_performance[1][3])%")
        println(f, "#")
        println(f, "# Pulse Shape Rankings:")
        for (i, (name, ratio, efficiency)) in enumerate(pulse_performance)
            println(f, "# $(i). $(name) - I/I₀: $(ratio), Efficiency: $(efficiency)%")
        end
    end
end

println("\n✅ Pulse shape analysis complete!")
println("📊 PyMOL files generated for all $(length(experiments)) pulse shapes")
println("📋 Comprehensive optimization recommendations")
println("📁 All output files saved to directory: $(output_dir)/")
println("💡 Check 'std_optimization_recommendations.txt' for detailed guidance")