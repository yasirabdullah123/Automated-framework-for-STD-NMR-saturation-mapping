using Plots
using DelimitedFiles
using Printf
using Statistics
using Clustering
using Distances
using Colors

println("=== Creating Enhanced Pulse Shape Heatmaps ===")

# Function to read residue information from master peak list (no secondary structure)
function get_residue_info()
    println("Loading residue information from master peak list...")
    
    # Read the master peak list (no headers)
    data = readdlm("../data/master-peak-list-final.tsv", '\t')
    
    # Create residue info dictionary
    residue_info = Dict()
    
    for row in eachrow(data)
        # Parse the atom identifier (e.g., "112.ARG.HDy")
        atom_id = String(row[1])
        parts = split(atom_id, ".")
        
        if length(parts) >= 3
            resi = parse(Int, parts[1])
            resn = parts[2]
            carbon_id = String(row[2])
            carbon_parts = split(carbon_id, ".")
            atomC = carbon_parts[3]
            
            key = "$(resi)-$(atomC)"
            residue_info[key] = (resi=resi, resn=resn, atomC=atomC)
        end
    end
    
    println("Loaded information for $(length(residue_info)) atom-residue combinations")
    return residue_info
end

# Function to organize lysozyme carbons by type with proper residue-aware classification
function organize_lysozyme_carbons(data, residue_info)
    println("Organizing carbons by type with residue-aware classification...")
    
    # Define aromatic residue types for lysozyme
    aromatic_residues = Set(["PHE", "TYR", "TRP", "HIS"])
    
    carbon_groups = Dict(
        "backbone" => [],
        "sidechain_base" => [],
        "sidechain_mid" => [],
        "aromatic" => [],
        "sidechain_term" => []
    )
    
    # Update aromatic groupings to include non-standard names
    aromatic_names = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH2", "CDx", "CDy", "CGx", "CGy"]
    for entry in data
        resi, atomC = entry[1], entry[2]
        
        # Get residue name from residue_info
        key = "$(resi)-$(atomC)"
        resn = haskey(residue_info, key) ? residue_info[key].resn : "UNK"
        
        # Categorize by carbon type with residue awareness
        if atomC == "CA"
            push!(carbon_groups["backbone"], entry)
        elseif atomC == "CB"
            push!(carbon_groups["sidechain_base"], entry)
        elseif atomC in ["CG", "CG1", "CG2"]
            push!(carbon_groups["sidechain_mid"], entry)
        elseif resn in aromatic_residues && atomC in aromatic_names
            # Aromatic carbons (only for aromatic residues)
            push!(carbon_groups["aromatic"], entry)
        else
            # All other side chain carbons
            push!(carbon_groups["sidechain_term"], entry)
        end
    end
    
    # Print statistics
    for (carbon_type, entries) in carbon_groups
        println("$(carbon_type): $(length(entries)) atoms")
    end
    
    return carbon_groups
end

# Function to read PyMOL file
function read_pymol_file(filepath)
    data = []
    open(filepath, "r") do f
        for line in eachline(f)
            if !startswith(line, "#") && !isempty(strip(line))
                parts = split(strip(line))
                if length(parts) >= 3
                    resi = parse(Int, parts[1])
                    atomC = parts[2]
                    ratio = parse(Float64, parts[3])
                    push!(data, (resi, atomC, ratio))
                end
            end
        end
    end
    return data
end

# Function to create true atom-level matrix for pulse shape data
function create_atom_level_matrix_pulse_shape(carbon_data, pulse_shapes, pymol_files, pymol_dir)
    println("Creating atom-level matrix for pulse shape analysis...")
    
    # Get unique atoms
    unique_atoms = unique(["$(d[1])-$(d[2])" for d in carbon_data])
    sort!(unique_atoms)
    
    println("Found $(length(unique_atoms)) unique atoms")
    
    # Create atom-level matrix (atoms × pulse shapes)
    data_matrix = fill(NaN, length(unique_atoms), length(pulse_shapes))
    
    # Create atom labels
    atom_labels = []
    for atom_key in unique_atoms
        resi_str, atomC = split(atom_key, "-")
        resi = parse(Int, resi_str)
        push!(atom_labels, (resi, atomC))
    end
    
    # Fill matrix with data
    for (pulse_idx, pulse_info) in enumerate(pulse_shapes)
        # Robust file pattern matching
        pattern = Regex("$(pulse_info.exp)_.*\\.txt\$")
        matching_files = filter(x -> occursin(pattern, x), pymol_files)
        if !isempty(matching_files)
            filepath = joinpath(pymol_dir, matching_files[1])
            pulse_data = read_pymol_file(filepath)
            
            # Create lookup dictionary for this pulse shape (handle duplicates by averaging)
            pulse_dict = Dict()
            for d in pulse_data
                key = "$(d[1])-$(d[2])"
                if haskey(pulse_dict, key)
                    # Average with existing value
                    pulse_dict[key] = (pulse_dict[key] + d[3]) / 2.0
                else
                    pulse_dict[key] = d[3]
                end
            end
            
            # Fill matrix column
            for (atom_idx, atom_key) in enumerate(unique_atoms)
                if haskey(pulse_dict, atom_key)
                    data_matrix[atom_idx, pulse_idx] = pulse_dict[atom_key]
                end
            end
        end
    end
    
    println("Matrix created: $(size(data_matrix))")
    return data_matrix, atom_labels
end

# Function to create enhanced pulse shape heatmap
function create_true_atom_level_pulse_shape_heatmap(carbon_data, pulse_shapes, pymol_files, pymol_dir, carbon_type, title, output_filename, residue_info)
    println("\nCreating $(carbon_type) pulse shape heatmap...")
    
    # Create atom-level matrix
    data_matrix, atom_labels = create_atom_level_matrix_pulse_shape(carbon_data, pulse_shapes, pymol_files, pymol_dir)
    
    # Create enhanced labels with residue names
    enhanced_labels = []
    for (resi, atomC) in atom_labels
        key = "$(resi)-$(atomC)"
        if haskey(residue_info, key)
            info = residue_info[key]
            label = "$(info.resn)$(resi)-$(atomC)"
            push!(enhanced_labels, label)
        else
            label = "$(resi)-$(atomC)"
            push!(enhanced_labels, label)
        end
    end
    
    # Apply hierarchical clustering to atoms (rows)
    println("Applying hierarchical clustering to atoms...")
    valid_data = data_matrix[.!isnan.(data_matrix)]
    
    if length(valid_data) > 0
        # Use only complete rows for clustering
        complete_rows = [i for i in 1:size(data_matrix, 1) if !any(isnan.(data_matrix[i, :]))]
        
        if length(complete_rows) > 1
            clustering_data = data_matrix[complete_rows, :]
            
            # Calculate distance matrix (cluster atoms, not pulse shapes)
            # Calculate distance matrix for hierarchical clustering
            # Note: Euclidean distance clusters primarily by saturation magnitude.
            # For pathway-specific (shape) clustering, Pearson correlation distance
            # is recommended for future implementations.
            distances = pairwise(Euclidean(), clustering_data, dims=1)
            
            # Perform hierarchical clustering
            try
                hclust_result = hclust(distances, linkage=:ward)
                
                # Get cluster order
                cluster_order = hclust_result.order
                
                # Validate cluster order
                if length(cluster_order) == length(complete_rows) && all(1 .<= cluster_order .<= length(complete_rows))
                    # Reorder data and labels
                    data_matrix_clustered = data_matrix[complete_rows[cluster_order], :]
                    enhanced_labels_clustered = [enhanced_labels[i] for i in complete_rows[cluster_order]]
                    atoms_clustered = [atom_labels[i] for i in complete_rows[cluster_order]]
                    
                    println("Clustering successful: $(length(cluster_order)) atoms clustered")
                else
                    println("Invalid cluster order, using original order")
                    data_matrix_clustered = data_matrix
                    enhanced_labels_clustered = enhanced_labels
                    atoms_clustered = atom_labels
                end
            catch e
                println("Clustering failed: $(e), using original order")
                data_matrix_clustered = data_matrix
                enhanced_labels_clustered = enhanced_labels
                atoms_clustered = atom_labels
            end
        else
            println("Insufficient complete rows for clustering, using original order")
            data_matrix_clustered = data_matrix
            enhanced_labels_clustered = enhanced_labels
            atoms_clustered = atom_labels
        end
    else
        println("No valid data for clustering, using original order")
        data_matrix_clustered = data_matrix
        enhanced_labels_clustered = enhanced_labels
        atoms_clustered = atom_labels
    end
    
    # Create pulse shape labels
    pulse_labels = [pulse_info.name for pulse_info in pulse_shapes]
    
    # Prepare data for plotting
    valid_data = data_matrix_clustered[.!isnan.(data_matrix_clustered)]
    if length(valid_data) > 0
        data_min = minimum(valid_data)
        data_max = maximum(valid_data)
        println("Data range: $(round(data_min, digits=3)) to $(round(data_max, digits=3))")
    end
    
    # Create reverse plasma color scheme
    color_scheme = cgrad(:plasma, 100, rev=true)
    
    # Calculate optimal figure size based on number of atoms
    n_atoms = size(data_matrix_clustered, 1)
    n_pulse_shapes = size(data_matrix_clustered, 2)
    
    # Improved figure sizing - wider and taller
    fig_width = max(1600, n_pulse_shapes * 200)   # Wider for pulse shape data
    fig_height = max(1200, n_atoms * 10)          # Much taller
    
    # Create the heatmap with improved axis handling
    heatmap_plot = heatmap(
        data_matrix_clustered,
        xlabel="Pulse Shape",
        ylabel="Protein Atoms",
        title=title,
        color=color_scheme,
        clims=(0.0, 1.0),  # Fixed range for I/I₀ ratios
        size=(fig_width, fig_height),
        dpi=300,
        margin=20Plots.mm,
        left_margin=220Plots.mm,  # Much more space for labels
        bottom_margin=80Plots.mm,
        right_margin=20Plots.mm,
        top_margin=40Plots.mm,
        colorbar_title="I/I₀ Ratio",
        colorbar_titlefontsize=12,
        colorbar_tickfontsize=10,
        titlefontsize=16,
        xlabelfontsize=14,
        ylabelfontsize=14,
        xtickfontsize=12,
        ytickfontsize=7,  # Smaller font to fit better
        grid=false,
        framestyle=:box
    )
    
    # Add custom x-axis labels with proper rotation (all 4 pulse shapes)
    plot!(xticks=(1:length(pulse_labels), pulse_labels), xrotation=45)
    
    # IMPROVED y-axis labels - show as many as possible based on number of rows
    if n_atoms <= 20
        # Show all atom labels for small datasets
        y_tick_indices = 1:n_atoms
        y_tick_labels = enhanced_labels_clustered
    elseif n_atoms <= 50
        # Show every 2nd label for medium datasets
        y_tick_step = 2
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels_clustered[i] for i in y_tick_indices if i <= length(enhanced_labels_clustered)]
    elseif n_atoms <= 100
        # Show every 3rd label for larger datasets
        y_tick_step = 3
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels_clustered[i] for i in y_tick_indices if i <= length(enhanced_labels_clustered)]
    else
        # Show every 4th label for very large datasets
        y_tick_step = 4
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels_clustered[i] for i in y_tick_indices if i <= length(enhanced_labels_clustered)]
    end
    
    # Apply clean labels with no rotation
    plot!(yticks=(y_tick_indices, y_tick_labels))
    
    # Save the heatmap
    savefig(heatmap_plot, output_filename)
    println("Saved: $(output_filename)")
    
    # Also save as PNG
    png_filename = replace(output_filename, ".svg" => ".png")
    savefig(heatmap_plot, png_filename)
    println("Saved: $(png_filename)")
    
    return data_matrix_clustered, atoms_clustered
end

# Function to create comprehensive pulse shape heatmaps
function create_comprehensive_pulse_shape_heatmaps(pymol_dir::String="pulse_shape_analysis")
    println("\n=== Creating Comprehensive Pulse Shape Heatmaps ===")
    
    pymol_files = filter(x -> startswith(x, "pymol_") && endswith(x, ".txt"), readdir(pymol_dir))
    sort!(pymol_files)
    
    # Define pulse shapes
    pulse_shapes = [
        (exp="../data/raw_experiments/34", name="50ms Gaussian", duration="50ms", shape="Gaussian"),
        (exp="../data/raw_experiments/35", name="5ms Gaussian", duration="5ms", shape="Gaussian"),
        (exp="../data/raw_experiments/36", name="20ms EBURP1", duration="20ms", shape="EBURP1"),
        (exp="../data/raw_experiments/41", name="2.5ms Gaussian", duration="2.5ms", shape="Gaussian")
    ]
    
    println("Found $(length(pulse_shapes)) pulse shapes:")
    for pulse_info in pulse_shapes
        println("  $(pulse_info.name)")
    end
    
    # Load residue information
    residue_info = get_residue_info()
    
    # Read all data
    all_data = []
    for file in pymol_files
        filepath = joinpath(pymol_dir, file)
        data = read_pymol_file(filepath)
        append!(all_data, data)
    end
    
    # Remove duplicates and organize by carbon type
    unique_data = unique(all_data)
    carbon_groups = organize_lysozyme_carbons(unique_data, residue_info)
    
    # Print carbon type statistics
    println("\nCarbon type distribution:")
    for (carbon_type, data) in carbon_groups
        println("$(carbon_type): $(length(data)) atoms")
    end
    
    # Create output directory
    output_dir = "pulse_shape_analysis_enhanced"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    println("\n=== Creating 4 Targeted Pulse Shape Heatmaps ===")
    
    # 1. Backbone Network (individual CA atoms)
    if !isempty(carbon_groups["backbone"])
        create_true_atom_level_pulse_shape_heatmap(
            carbon_groups["backbone"], pulse_shapes, pymol_files, pymol_dir,
            "backbone", "Lysozyme Backbone Pulse Shape Response (CA atoms)",
            joinpath(output_dir, "1_backbone_CA_pulse_shape.svg"), residue_info
        )
    end
    
    # 2. Side Chain Base (individual CB atoms)
    if !isempty(carbon_groups["sidechain_base"])
        create_true_atom_level_pulse_shape_heatmap(
            carbon_groups["sidechain_base"], pulse_shapes, pymol_files, pymol_dir,
            "sidechain_base", "Side Chain Base Pulse Shape Response (CB atoms)",
            joinpath(output_dir, "2_sidechain_CB_pulse_shape.svg"), residue_info
        )
    end
    
    # 3. Aromatic Carbons (binding hotspots)
    if !isempty(carbon_groups["aromatic"])
        create_true_atom_level_pulse_shape_heatmap(
            carbon_groups["aromatic"], pulse_shapes, pymol_files, pymol_dir,
            "aromatic", "Aromatic Pulse Shape Response (Phe/Tyr/Trp)",
            joinpath(output_dir, "3_aromatic_pulse_shape.svg"), residue_info
        )
    end
    
    # 4. Terminal Carbons (deep side chain penetration)
    if !isempty(carbon_groups["sidechain_term"])
        create_true_atom_level_pulse_shape_heatmap(
            carbon_groups["sidechain_term"], pulse_shapes, pymol_files, pymol_dir,
            "sidechain_term", "Side Chain Terminal Pulse Shape Response (CD/CE/CZ)",
            joinpath(output_dir, "4_sidechain_terminals_pulse_shape.svg"), residue_info
        )
    end
    
    # Generate summary statistics
    println("\n=== Pulse Shape Analysis Summary ===")
    for (carbon_type, data) in carbon_groups
        if !isempty(data)
            unique_atoms = length(unique(["$(d[1])-$(d[2])" for d in data]))
            println("$(carbon_type): $(unique_atoms) unique atoms")
        end
    end
    
    println("\n✅ Pulse shape analysis complete!")
    println("📁 Results saved to: $(output_dir)/")
end

# Run the comprehensive pulse shape analysis
println("=== Starting Pulse Shape Heatmap Generation ===")
if abspath(PROGRAM_FILE) == @__FILE__
    dir = length(ARGS) > 0 ? ARGS[1] : "pulse_shape_analysis"
    create_comprehensive_pulse_shape_heatmaps(dir)
end

println("\n🎉 Pulse shape heatmap generation complete!")
println("📊 Generated 4 targeted pulse shape heatmaps:")
println("   • Backbone pulse shape response (CA atoms)")
println("   • Side chain base pulse shape response (CB atoms)")
println("   • Aromatic pulse shape response (Phe/Tyr/Trp)")
println("   • Side chain terminals pulse shape response (CD/CE/CZ)")
println("🔬 True atom-level pulse shape resolution achieved!")
println("📁 Files saved to: pulse_shape_analysis_enhanced/") 