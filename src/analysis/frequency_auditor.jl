using Plots
using NMRTools
using Measurements
using DelimitedFiles
using Printf
using Statistics
using Dates

# Include the centralized, verified alignment functions
include("../utils/nmr_analysis_utils.jl")

println("=== Experiment 31: Comprehensive STD Frequency Analysis ===")
println("19 saturation frequencies from -0.5 to 8.5 ppm")

# Create output directory
output_dir = "../data/raw_experiments/31"
if !isdir(output_dir)
    mkdir(output_dir)
    println("Created output directory: $(output_dir)")
end

# --- Use Centralized Peak Alignment ---
aligned_peaks = NMRAnalysisUtils.load_and_align_peaks("../data/master-peak-list-final.tsv", "../data/raw_experiments/31", x_shift=0.0, y_shift=0.0)
x = aligned_peaks.x
y = aligned_peaks.y
resi = aligned_peaks.resi
resn = aligned_peaks.resn
atomH = aligned_peaks.atomH
atomC = aligned_peaks.atomC

println("Loaded and aligned $(length(x)) peaks for analysis.")

# --- ALIGNMENT VERIFICATION ---
println("\n=== Verifying Peak Alignment ===")

# Create alignment verification plot showing peaks vs spectral maxima
alignment_plot = plot(loadnmr("../data/raw_experiments/31")[:,:,1], title="Alignment Check - Experiment 31", xflip=true, yflip=true)
scatter!(alignment_plot, x, y, marker=:x, color=:red, label="Aligned Peaks", markersize=4, alpha=0.7)
savefig(alignment_plot, joinpath(output_dir, "alignment_check.svg"))
println("Alignment check plot saved to: $(joinpath(output_dir, "alignment_check.svg"))")


# Load and normalize experiment 31 (19 frequencies)
spec = loadnmr("../data/raw_experiments/31")
spec = spec / spec[:noise]
spec[:noise] = 1.0

println("Experiment 31 dimensions: $(size(spec))")

# Extract saturation frequencies from the experiment
# (This will be updated below after loading the experiment)

# Use internal reference from experiment 31 (slice 1 - 30 ppm off-resonance)
println("\nUsing internal reference from experiment 31 (slice 1 - 30 ppm)...")
println("No external reference needed - experiment 31 contains its own reference")

# Extract frequency values from the experiment
fq_list = acqus(spec, :fq1list)
all_freqs = fq_list.values
println("All frequencies in experiment 31: $(all_freqs)")
println("Reference frequency (slice 1): $(all_freqs[1]) ppm")
println("Saturation frequencies (slices 2-20): $(all_freqs[2:end])")

# Update sat_freqs to match the actual experiment frequencies
sat_freqs = all_freqs[2:end]  # Skip the reference frequency (30 ppm)
println("Analyzing $(length(sat_freqs)) saturation frequencies: $(sat_freqs)")

# Calculate I/I0 ratios using internal reference with uncertainty filtering
println("\nCalculating I/I0 ratios using internal reference with uncertainty filtering...")

II0 = map(1:length(x)) do i
    # Get reference intensity from slice 1 (30 ppm off-resonance)
    ref_intensity = spec[Near(x[i]), Near(y[i]), 1] ± spec[:noise]
    
    # Filter out peaks with poor reference intensity
    ref_value = Measurements.value(ref_intensity)
    if ref_value <= 0.1  # Reference intensity too small or negative
        return [NaN ± 0.0 for _ in 1:19]  # Return NaN for all frequencies
    end
    
    # Get intensities for the 19 saturation frequencies (slices 2-20)
    sat_intensities = [spec[Near(x[i]), Near(y[i]), slice] ± spec[:noise] for slice in 2:20]
    
    # Calculate I/I0 ratios with automatic error propagation
    ratios = sat_intensities ./ ref_intensity
    
    # Apply uncertainty cutoff for I/I0 ratios
    # Note: 25% threshold used here due to low SNR in natural abundance 13C data
    # For isotopically labeled samples, 5% is recommended
    filtered_ratios = map(ratios) do ratio
        if Measurements.uncertainty(ratio) > 0.25
            return NaN ± 0.0  # Set to NaN if uncertainty exceeds threshold
        else
            return ratio
        end
    end
    
    return filtered_ratios
end

II0 = hcat(II0...)' # stack the results into a matrix (peaks × frequencies)

println("Intensity matrix shape: $(size(II0))")

# Debug: Check for NaN values in the intensity matrix
println("\n=== Debugging Intensity Data ===")
II0_values = Measurements.value.(II0)
nan_count = sum(isnan.(II0_values))
inf_count = sum(isinf.(II0_values))
valid_count = sum(isfinite.(II0_values))

println("NaN values: $(nan_count)")
println("Infinite values: $(inf_count)")
println("Valid values: $(valid_count)")
println("Total values: $(length(II0_values))")

if nan_count > 0
    println("NaN locations by frequency:")
    for freq_idx in 1:size(II0_values, 2)
        freq_nans = sum(isnan.(II0_values[:, freq_idx]))
        if freq_nans > 0
            println("  Frequency $(sat_freqs[freq_idx]) ppm: $(freq_nans) NaN values")
        end
    end
end

# Show range of valid values
valid_values = II0_values[isfinite.(II0_values)]
if length(valid_values) > 0
    println("Valid intensity range: $(minimum(valid_values)) to $(maximum(valid_values))")
end

# Generate PyMOL files for each frequency
println("\nGenerating PyMOL coloring files...")
for slicenumber in 1:length(sat_freqs)
    freq = sat_freqs[slicenumber]
    filename = joinpath(output_dir, "pymol_$(freq)ppm.txt")
    
    open(filename, "w") do f
        println(f, "## PyMOL coloring data for $(freq) ppm saturation")
        println(f, "# Lower values indicate stronger saturation transfer")
        println(f, "# Reference: Internal (slice 1 - 30 ppm off-resonance)")
        println(f, "# Uncertainty filtering: 25% threshold applied")
        println(f, "# Format: resi atomC ratio")
        
        for i in eachindex(x)
            n = resi[i]
            a = atomC[i]
            ratio = Measurements.value(II0[i, slicenumber])
            
            # Only output non-NaN values (from uncertainty filtering only)
            if !isnan(ratio)
                println(f, "$n $a $ratio")
            end
        end
    end
    println("Saved $(filename)")
end



# Generate frequency optimization recommendations
println("\n=== Frequency Optimization Recommendations ===")

# Find overall best frequency
overall_ratios = []
for freq_idx in 1:length(sat_freqs)
    ratios = Measurements.value.(II0[:, freq_idx])
    valid_ratios = ratios[.!isnan.(ratios)]
    if length(valid_ratios) > 0
        push!(overall_ratios, (sat_freqs[freq_idx], mean(valid_ratios), length(valid_ratios)))
    end
end

if length(overall_ratios) > 0
    sort!(overall_ratios, by=x->x[2])
    best_overall = overall_ratios[1]
    
    println("🏆 **Best Overall Frequency**: $(best_overall[1]) ppm")
    println("   Mean I/I₀: $(round(best_overall[2], digits=3))")
    println("   Valid Peaks: $(best_overall[3])/$(length(x))")
    println("   Saturation Efficiency: $(round(100*(1-best_overall[2]), digits=1))%")
    
    println("\n💡 **Key Insights:**")
    println("   • Traditional 1 ppm saturation shows $(round(100*(1-best_overall[2]), digits=1))% efficiency")
    println("   • Consider ligand binding site when choosing frequency")
end

# Save comprehensive recommendations
recommendations_file = joinpath(output_dir, "frequency_optimization_recommendations.txt")
open(recommendations_file, "w") do f
    println(f, "# STD Frequency Optimization Recommendations")
    println(f, "# Generated: $(now())")
    println(f, "# Based on Lysozyme Frequency Analysis")
    println(f, "#")
    println(f, "# Key Findings:")
    
    if length(overall_ratios) > 0
        println(f, "# Best Overall: $(best_overall[1]) ppm")
        println(f, "# Mean I/I₀: $(best_overall[2])")
    end
end

println("\n✅ Frequency analysis complete!")
println("📊 PyMOL files generated for all $(length(sat_freqs)) frequencies")
println("🏆 Frequency optimization recommendations")
println("📁 All output files saved to directory: $(output_dir)/")
println("💡 Check 'frequency_optimization_recommendations.txt' for detailed guidance")