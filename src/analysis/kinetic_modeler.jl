using Plots
using NMRTools
using Measurements
using DelimitedFiles
using Printf
using Statistics
using Dates
using StatsPlots

include("../utils/nmr_analysis_utils.jl")

println("=== Saturation Time Analysis: STD Kinetics Study ===")
println("Analyzing saturation transfer buildup over time")
println("Experiment 33: Variable saturation times from 0.1 to 5.0 seconds")

# Create output directory
output_dir = "sattimeanalysis"
if !isdir(output_dir)
    mkdir(output_dir)
    println("Created output directory: $(output_dir)")
end

# --- Use Centralized Peak Alignment ---
# Load and align peaks using the verified function from nmr_analysis_utils.jl
# The reference for folding is now the target experiment itself (exp 33)
aligned_peaks = NMRAnalysisUtils.load_and_align_peaks("../data/master-peak-list-final.tsv", "../data/raw_experiments/33", x_shift=0.0, y_shift=0.0)
x = aligned_peaks.x
y = aligned_peaks.y
resi = aligned_peaks.resi
resn = aligned_peaks.resn
atomH = aligned_peaks.atomH
atomC = aligned_peaks.atomC

println("Loaded and aligned $(length(x)) peaks for analysis.")

# --- ALIGNMENT VERIFICATION ---
println("\n=== Verifying Peak Alignment ===")

# Load experiment 33 for plotting
spec = loadnmr("../data/raw_experiments/33")
firstplane = spec[:,:,1]

# Simple alignment check plot (like prepare-peak-list.jl)
plot(firstplane)
scatter!(x, y, xflip=true, yflip=true, 
    xlabel="¹H Chemical Shift (ppm)", ylabel="¹³C Chemical Shift (ppm)", 
    title="Peak List Alignment Check (Exp 33)", legend=false, grid=false)
savefig(joinpath(output_dir, "alignment_check.svg"))
println("Alignment check plot saved to: $(joinpath(output_dir, "alignment_check.svg"))")

# Load and normalize experiment 33 (saturation time series)
spec = loadnmr("../data/raw_experiments/33")
spec = spec / spec[:noise]
spec[:noise] = 1.0

println("Experiment 33 dimensions: $(size(spec))")

# Load reference spectrum (experiment 37/pdata/23 - 30 ppm off-resonance)
println("\nLoading reference spectrum (experiment 37/pdata/23)...")
ref_spec = loadnmr("37/pdata/23")
ref_spec = ref_spec / ref_spec[:noise]
ref_spec[:noise] = 1.0
println("Reference spectrum dimensions: $(size(ref_spec))")

# Define saturation times from vdlist
sat_times = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]  # vdlist times

println("Analyzing $(length(sat_times)) saturation time points:")
println("Reference: Experiment 37/pdata/23 (30 ppm off-resonance)")
println("Time points: $(sat_times) s (slices 1-$(length(sat_times)))")

# Measure peak heights for all time points
println("\nMeasuring peak intensities across all time points...")
println("Using ONLY 25% uncertainty filtering (no S/N filtering)")

II0 = map(1:length(x)) do i
    # Get reference intensity from experiment 37/pdata/23
    ref_intensity = ref_spec[Near(x[i]), Near(y[i])] ± ref_spec[:noise]
    
    # Filter out peaks with poor reference intensity
    ref_value = Measurements.value(ref_intensity)
    if ref_value <= 0.1  # Reference intensity too small or negative
        return [NaN ± 0.0 for _ in axes(spec, 3)]  # Return NaN for all time points
    end
    
    # Get data at all time points (slices 1 to end) with measurement errors
    sat_intensities = [spec[Near(x[i]), Near(y[i]), slice] ± spec[:noise] for slice in axes(spec, 3)]
    
    # Calculate ratios relative to experiment 37 reference with automatic error propagation
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

II0_matrix = hcat(II0...)' # stack the results into a matrix (peaks × time points)

println("Intensity ratio matrix shape: $(size(II0_matrix))")

# Generate PyMOL files for each time point
println("\nGenerating PyMOL coloring files...")
num_time_points = size(II0_matrix, 2)
for idx in 1:num_time_points
    sat_time = sat_times[idx]  # Direct indexing now that we analyze all time points
    filename = joinpath(output_dir, "pymol_$(sat_time)s.txt")
    
    open(filename, "w") do f
        println(f, "## PyMOL coloring data for $(sat_time) s saturation")
        println(f, "# Lower values indicate stronger saturation transfer")
        println(f, "# Reference: Experiment 37/pdata/23 (30 ppm off-resonance)")
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

