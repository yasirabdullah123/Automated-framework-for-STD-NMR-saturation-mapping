module NMRAnalysisUtils

export load_and_align_peaks, fold_coordinates_robust, test_and_visualize_alignment, estimate_noise

using Plots
using NMRTools
using DelimitedFiles
using Statistics

"""
    fold_coordinates_robust(coords, min_val, max_val, sw)

Folds spectral coordinates into the defined spectral window. This function is 
physically correct and handles any degree of folding by repeatedly adding or 
subtracting the spectral width (sw) until the coordinate falls within the 
window defined by min_val and max_val.

# Arguments
- `coords`: A vector of spectral coordinates to be folded.
- `min_val`: The minimum chemical shift of the spectral window.
- `max_val`: The maximum chemical shift of the spectral window.
- `sw`: The spectral width (in ppm).

# Returns
- A new vector of folded coordinates.
"""
function fold_coordinates_robust(coords, axis_start, axis_end, sw)
    # Correctly determine the min and max of the spectral window
    min_val = min(axis_start, axis_end)
    max_val = max(axis_start, axis_end)

    folded = copy(coords)
    for i in eachindex(folded)
        temp_coord = folded[i]
        # Use a while loop to handle any degree of folding
        if temp_coord < min_val
            while temp_coord < min_val
                temp_coord += sw
            end
        elseif temp_coord > max_val
            while temp_coord > max_val
                temp_coord -= sw
            end
        end
        # Final check in case folding overshoots
        if temp_coord > max_val
            temp_coord -= sw
        elseif temp_coord < min_val
            temp_coord += sw
        end
        folded[i] = temp_coord
    end
    return folded
end

"""
    load_and_align_peaks(peak_list_path, ref_exp_path; x_shift=0.0, y_shift=0.0)

Loads a master peak list, aligns it to a reference spectrum, and applies robust folding.

# Arguments
- `peak_list_path`: Path to the master peak list file (e.g., "../data/master-peak-list-final.tsv").
- `ref_exp_path`: Path to the reference NMR experiment directory (e.g., "../data/raw_experiments/4").
- `x_shift`: Optional systematic shift to apply to the ¹H dimension.
- `y_shift`: Optional systematic shift to apply to the ¹³C dimension.

# Returns
- A NamedTuple containing the aligned and folded peak data: `(x, y, resi, resn, atomH, atomC)`.
"""
function load_and_align_peaks(peak_list_path, ref_exp_path; x_shift=0.0, y_shift=0.0)
    println("--- Loading and Aligning Peaks (RBH-Optimized Format) ---")
    println("Peak List: $peak_list_path")
    println("Reference Experiment: $ref_exp_path")

    # BEST PRACTICE: Use readdlm with header=true to handle the new file format.
    # This automatically separates the data from the header row.
    data_matrix, header = readdlm(peak_list_path, '\t', header=true)

    # The new file has 7 columns. We only need columns 1-4 for this function.
    # We parse from `data_matrix`, which correctly contains only the data.
    x = float.(data_matrix[:, 3])  # ¹H Shift is still in Column 3
    y = float.(data_matrix[:, 4])  # ¹³C Shift is still in Column 4

    # Use `data_matrix` for assignments to ensure perfect alignment
    resi = map(data_matrix[:,1]) do assignment
        split(string(assignment), '.')[1]
    end
    resn = map(data_matrix[:,1]) do assignment
        split(string(assignment), '.')[2]
    end
    atomH = map(data_matrix[:,1]) do assignment
        split(string(assignment), '.')[3]
    end
    atomC = map(data_matrix[:,2]) do assignment
        split(string(assignment), '.')[3]
    end
    println("Loaded $(length(x)) high-confidence peaks from RBH-optimized master list.")

    # Apply systematic shifts if provided
    if x_shift != 0.0 || y_shift != 0.0
        x .+= x_shift
        y .+= y_shift
        println("Applied systematic shifts: H += $x_shift ppm, C += $y_shift ppm")
    end

    # Load reference spectrum to get folding parameters
    println("Loading reference spectrum for folding parameters...")
    ref_spec = loadnmr(ref_exp_path)

    # Get folding parameters for the ¹H dimension (F1)
    xsw = metadata(ref_spec, F1Dim, :swppm)
    x_axis = data(ref_spec, F1Dim)
    xmin, xmax = first(x_axis), last(x_axis)

    # Get folding parameters for the ¹³C dimension (F2)
    ysw = metadata(ref_spec, F2Dim, :swppm)
    y_axis = data(ref_spec, F2Dim)
    ymin, ymax = first(y_axis), last(y_axis)

    println("Folding Parameters (from Exp $ref_exp_path):")
    println("  H (F1): SW = $(round(xsw, digits=2)) ppm, Range = [$(round(xmin, digits=2)), $(round(xmax, digits=2))]")
    println("  C (F2): SW = $(round(ysw, digits=2)) ppm, Range = [$(round(ymin, digits=2)), $(round(ymax, digits=2))]")

    # Apply robust folding to both dimensions
    x_folded = fold_coordinates_robust(x, xmin, xmax, xsw)
    y_folded = fold_coordinates_robust(y, ymin, ymax, ysw)
    println("Applied robust folding to both dimensions.")

    return (x=x_folded, y=y_folded, resi=resi, resn=resn, atomH=atomH, atomC=atomC)
end

"""
    test_and_visualize_alignment(peak_list_path, target_exp_path, target_slice, output_path; x_shift=-0.11, y_shift=-3.05)

Tests and visualizes the alignment of a master peak list to a target spectrum.

# Arguments
- `peak_list_path`: Path to the master peak list file.
- `target_exp_path`: Path to the target NMR experiment directory.
- `target_slice`: The slice number of the target spectrum to use for alignment.
- `output_path`: Path to save the output visualization plot.
- `x_shift`: Systematic shift for the ¹H dimension.
- `y_shift`: Systematic shift for the ¹³C dimension.
"""
function test_and_visualize_alignment(peak_list_path, target_exp_path, target_slice, output_path; x_shift=-0.11, y_shift=-3.05)
    println("--- Running Alignment Test ---")
    println("Peak List: $peak_list_path")
    println("Target Spectrum: $target_exp_path, slice $target_slice")

    # Load master peak list
    results = readdlm(peak_list_path)
    x = float.(results[:, 3])
    y = float.(results[:, 4])
    println("Loaded $(length(x)) peaks from master list.")

    # Apply systematic shifts
    x .+= x_shift
    y .+= y_shift
    println("Applied systematic shifts: H += $x_shift ppm, C += $y_shift ppm")

    # Load target spectrum for folding parameters and visualization
    println("Loading target spectrum...")
    target_spec = loadnmr(target_exp_path)
    target_plane = target_spec[:,:,target_slice]

    # Get folding parameters from the target spectrum
    xsw = metadata(target_spec, F1Dim, :swppm)
    x_axis = data(target_spec, F1Dim)
    xmin, xmax = first(x_axis), last(x_axis)

    ysw = metadata(target_spec, F2Dim, :swppm)
    y_axis = data(target_spec, F2Dim)
    ymin, ymax = first(y_axis), last(y_axis)
    
    println("Folding Parameters (from Exp $target_exp_path):")
    println("  H (F1): SW = $(round(xsw, digits=2)) ppm, Range = [$(round(xmin, digits=2)), $(round(xmax, digits=2))]")
    println("  C (F2): SW = $(round(ysw, digits=2)) ppm, Range = [$(round(ymin, digits=2)), $(round(ymax, digits=2))]")

    # Apply robust folding
    x_folded = fold_coordinates_robust(x, xmin, xmax, xsw)
    y_folded = fold_coordinates_robust(y, ymin, ymax, ysw)
    println("Applied robust folding to both dimensions.")

    # Create visualization plot
    println("Generating alignment visualization plot...")
    p = plot(target_plane, title="Alignment Test: Exp $target_exp_path, Slice $target_slice", xflip=true, yflip=true)
    scatter!(p, x_folded, y_folded, marker=:x, color=:red, label="Aligned Peaks", markersize=4, alpha=0.7)
    
    # Save the plot
    savefig(p, output_path)
    println("Alignment plot saved to: $output_path")
    println("--- Test Complete ---")
end

end # end of module NMRAnalysisUtils