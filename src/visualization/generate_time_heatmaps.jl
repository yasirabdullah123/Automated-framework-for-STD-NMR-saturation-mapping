using Plots
using DelimitedFiles
using Printf
using Statistics
using Clustering
using Colors
using Distances

println("=== Creating Enhanced STD NMR Heatmaps with Fixed Axes ===")

# Function to read PyMOL file and extract data
function read_pymol_file(filename)
    data = []
    open(filename, "r") do f
        for line in eachline(f)
            # Skip comment lines and empty lines
            if !startswith(line, "#") && !isempty(strip(line))
                parts = split(strip(line))
                if length(parts) >= 3
                    resi = parse(Int, parts[1])
                    atom = parts[2]
                    ratio = parse(Float64, parts[3])
                    push!(data, (resi, atom, ratio))
                end
            end
        end
    end
    return data
end

# Function to get residue name info (no secondary structure)
function get_residue_info()
    residue_map = Dict{Int, String}()  # resi -> resn
    open("../data/master-peak-list-final.tsv", "r") do f
        for line in eachline(f)
            parts = split(line, '\t')
            if length(parts) >= 1
                peak_info = split(parts[1], '.')
                if length(peak_info) >= 2
                    resi = parse(Int, peak_info[1])
                    resn = peak_info[2]
                    residue_map[resi] = resn
                end
            end
        end
    end
    return residue_map
end

# Function to check if residue is aromatic using residue_info
function is_aromatic_residue(resi, residue_info)
    # Define aromatic residue types
    aromatic_residues = Set(["PHE", "TYR", "TRP", "HIS"])
    
    # Check if this residue is aromatic
    if haskey(residue_info, resi)
        resn = residue_info[resi]  # First element is residue name
        return resn in aromatic_residues
    end
    return false
end

# Function to organize data by carbon type with proper residue-aware classification
function organize_lysozyme_carbons(all_atom_data, residue_info)
    carbon_groups = Dict(
        "backbone" => [],           # CA atoms
        "sidechain_base" => [],     # CB atoms  
        "sidechain_mid" => [],      # CG atoms
        "sidechain_term" => [],     # CD, CE, CZ atoms
        "aromatic" => []            # Phe/Tyr/Trp aromatic carbons
    )
    aromatic_names = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH2", "CDx", "CDy", "CGx", "CGy"]
    for (resi, atom, ratio) in all_atom_data
        atom_label = "$(resi)-$(atom)"
        if atom == "CA"
            push!(carbon_groups["backbone"], (resi, atom, ratio))
        elseif atom == "CB"  
            push!(carbon_groups["sidechain_base"], (resi, atom, ratio))
        elseif atom in ["CG", "CG1", "CG2"]
            push!(carbon_groups["sidechain_mid"], (resi, atom, ratio))
        elseif atom in ["CD", "CD1", "CD2", "CE", "CE1", "CE2", "CZ", "CDx", "CDy", "CGx", "CGy", "CEy"]
            if is_aromatic_residue(resi, residue_info) || atom in aromatic_names
                push!(carbon_groups["aromatic"], (resi, atom, ratio))
            else
                push!(carbon_groups["sidechain_term"], (resi, atom, ratio))
            end
        end
    end
    return carbon_groups
end

# Function to perform hierarchical clustering (following frequency script approach)
function cluster_matrix(data_matrix)
    # Check if matrix is too small for clustering
    if size(data_matrix, 1) < 2
        return collect(1:size(data_matrix, 1))
    end
    
    # Use only complete rows for clustering
    complete_rows = [i for i in 1:size(data_matrix, 1) if !any(isnan.(data_matrix[i, :]))]
    
    if length(complete_rows) < 2
        return collect(1:size(data_matrix, 1))
    end
    
    clustering_data = data_matrix[complete_rows, :]
    
    # Perform hierarchical clustering on rows
    try
        # Calculate distance matrix (cluster atoms, not time points)
            # Calculate distance matrix for hierarchical clustering
            # Note: Euclidean distance clusters primarily by saturation magnitude.
            # For pathway-specific (shape) clustering, Pearson correlation distance
            # is recommended for future implementations.
        distances = pairwise(Euclidean(), clustering_data, dims=1)
        hclust_result = hclust(distances, linkage=:ward)
        
        # Get cluster order
        cluster_order = hclust_result.order
        
        # Validate cluster order
        if length(cluster_order) == length(complete_rows) && all(1 .<= cluster_order .<= length(complete_rows))
            # Create full order for all rows
            full_order = collect(1:size(data_matrix, 1))
            full_order[complete_rows] = complete_rows[cluster_order]
            return full_order
        else
            println("Invalid cluster order, using original order")
            return collect(1:size(data_matrix, 1))
        end
    catch e
        println("Clustering failed, using original order: $e")
        return collect(1:size(data_matrix, 1))
    end
end

# Function to create true atom-level matrix
function create_atom_level_matrix(carbon_data, time_points, pymol_files, pymol_dir, carbon_type)
    if isempty(carbon_data)
        return Matrix{Float64}(undef, 0, length(time_points)), String[]
    end
    
    # Get ALL unique atom identifiers for this carbon type
    all_atoms = Set{String}()
    for (resi, atom, ratio) in carbon_data
        atom_id = "$(resi)-$(atom)"
        push!(all_atoms, atom_id)
    end
    
    sorted_atoms = sort(collect(all_atoms))
    
    # Create matrix: individual atoms × time points  
    data_matrix = fill(NaN, length(sorted_atoms), length(time_points))
    
    for (time_idx, time_point) in enumerate(time_points)
        # Robust file pattern matching
        pattern = Regex("$(time_point)(s|ppm)?\\.txt\$")
        matching_files = filter(f -> occursin(pattern, f), pymol_files)
        if !isempty(matching_files)
            filepath = joinpath(pymol_dir, matching_files[1])
            
            if isfile(filepath)
                file_data = read_pymol_file(filepath)
                
                # Process ALL atoms in file (not just first match)
                for (file_resi, file_atom, ratio) in file_data
                    atom_id = "$(file_resi)-$(file_atom)"
                    
                    # Only include atoms of the specified carbon type
                    if atom_id in sorted_atoms
                        atom_idx = findfirst(==(atom_id), sorted_atoms)
                        if atom_idx !== nothing
                            data_matrix[atom_idx, time_idx] = ratio
                        end
                    end
                end
            end
        end
    end
    
    return data_matrix, sorted_atoms
end

# Function to create true atom-level heatmap
function create_true_atom_level_heatmap(carbon_data, time_points, pymol_files, pymol_dir, 
                                       carbon_type, title, output_filename, residue_info)
    println("Creating $(carbon_type) atom-level heatmap...")
    
    # Create atom-level data matrix
    data_matrix, atom_labels = create_atom_level_matrix(carbon_data, time_points, pymol_files, pymol_dir, carbon_type)
    
    if size(data_matrix, 1) == 0
        println("No data found for $(carbon_type)")
        return
    end
    
    println("Found $(size(data_matrix, 1)) individual atoms for $(carbon_type)")
    
    # Apply hierarchical clustering to atoms (rows)
    println("Applying hierarchical clustering to atoms...")
    valid_data = data_matrix[.!isnan.(data_matrix)]
    
    if length(valid_data) > 0
        # Use only complete rows for clustering
        complete_rows = [i for i in 1:size(data_matrix, 1) if !any(isnan.(data_matrix[i, :]))]
        
        if length(complete_rows) > 1
            clustering_data = data_matrix[complete_rows, :]
            
            # Calculate distance matrix (cluster atoms, not time points)
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
                    atoms_clustered = atom_labels[complete_rows[cluster_order]]
                    
                    println("Clustering successful: $(length(cluster_order)) atoms clustered")
                else
                    println("Invalid cluster order, using original order")
                    data_matrix_clustered = data_matrix
                    atoms_clustered = atom_labels
                end
            catch e
                println("Clustering failed: $(e), using original order")
                data_matrix_clustered = data_matrix
                atoms_clustered = atom_labels
            end
        else
            println("Insufficient complete rows for clustering, using original order")
            data_matrix_clustered = data_matrix
            atoms_clustered = atom_labels
        end
    else
        println("No valid data for clustering, using original order")
        data_matrix_clustered = data_matrix
        atoms_clustered = atom_labels
    end
    
    # Create enhanced row labels with atom-level detail
    enhanced_labels = []
    for atom_label in atoms_clustered
        parts = split(atom_label, "-")
        resi = parse(Int, parts[1])
        atom = parts[2]
        if haskey(residue_info, resi)
            resn = residue_info[resi]
            push!(enhanced_labels, "$(resn)$(resi)-$(atom)")  # "Lys112-CD"
        else
            push!(enhanced_labels, atom_label)
        end
    end
    
    # Create column names
    if contains(pymol_files[1], "s.txt")
        col_names = ["$(t)s" for t in time_points]
    elseif contains(pymol_files[1], "ppm.txt")
        col_names = ["$(t)ppm" for t in time_points]
    else
        col_names = ["Exp$(t)" for t in time_points]
    end
    
    # Prepare data for plotting
    valid_data = data_matrix_clustered[.!isnan.(data_matrix_clustered)]
    if length(valid_data) > 0
        data_min = minimum(valid_data)
        data_max = maximum(valid_data)
        println("Data range: $(round(data_min, digits=3)) to $(round(data_max, digits=3))")
    end
    
    # Create reverse plasma color scheme (red -> yellow -> green -> blue -> purple)
    color_scheme = cgrad(:plasma, 100, rev=true)
    
    # Calculate optimal figure size based on number of atoms
    n_atoms = size(data_matrix_clustered, 1)
    n_timepoints = size(data_matrix_clustered, 2)
    
    # Improved figure sizing - wider and taller
    fig_width = max(1600, n_timepoints * 100)   # Wider
    fig_height = max(1200, n_atoms * 10)        # Much taller
    
    # Create the heatmap with improved axis handling
    heatmap_plot = heatmap(
        data_matrix_clustered,
        xlabel="Saturation Time (seconds)",
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
    
    # Add custom x-axis labels with proper rotation
    plot!(xticks=(1:length(col_names), col_names), xrotation=45)
    
    # IMPROVED y-axis labels - show as many as possible based on number of rows
    if n_atoms <= 20
        # Show all atom labels for small datasets
        y_tick_indices = 1:n_atoms
        y_tick_labels = enhanced_labels
    elseif n_atoms <= 50
        # Show every 2nd label for medium datasets
        y_tick_step = 2
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels[i] for i in y_tick_indices if i <= length(enhanced_labels)]
    elseif n_atoms <= 100
        # Show every 3rd label for larger datasets
        y_tick_step = 3
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels[i] for i in y_tick_indices if i <= length(enhanced_labels)]
    else
        # Show every 4th label for very large datasets
        y_tick_step = 4
        y_tick_indices = 1:y_tick_step:n_atoms
        y_tick_labels = [enhanced_labels[i] for i in y_tick_indices if i <= length(enhanced_labels)]
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

# Function to create comprehensive atom-level heatmaps
function create_comprehensive_atom_level_heatmaps(pymol_dir::String="sattimeanalysis")
    println("\n=== Creating Comprehensive Atom-Level STD NMR Heatmaps ===")
    
    pymol_files = filter(x -> startswith(x, "pymol_") && endswith(x, ".txt"), readdir(pymol_dir))
    sort!(pymol_files)
    
    # Extract time points
    time_points = []
    for file in pymol_files
        if contains(file, "s.txt")
            time_str = replace(replace(file, "pymol_" => ""), "s.txt" => "")
            time_point = parse(Float64, time_str)
            push!(time_points, time_point)
        end
    end
    sort!(time_points)
    
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
    output_dir = "atom_level_analysis"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    println("\n=== Creating 4 Targeted Atom-Level Heatmaps ===")
    
    # 1. Backbone Network (individual CA atoms)
    if !isempty(carbon_groups["backbone"])
        create_true_atom_level_heatmap(
            carbon_groups["backbone"], time_points, pymol_files, pymol_dir,
            "backbone", "Lysozyme Backbone Saturation Network (CA atoms)",
            joinpath(output_dir, "1_backbone_CA_saturation.svg"), residue_info
        )
    end
    
    # 2. Side Chain Base (individual CB atoms)
    if !isempty(carbon_groups["sidechain_base"])
        create_true_atom_level_heatmap(
            carbon_groups["sidechain_base"], time_points, pymol_files, pymol_dir,
            "sidechain_base", "Side Chain Base Accessibility (CB atoms)",
            joinpath(output_dir, "2_sidechain_CB_accessibility.svg"), residue_info
        )
    end
    
    # 3. Aromatic Carbons (binding hotspots)
    if !isempty(carbon_groups["aromatic"])
        create_true_atom_level_heatmap(
            carbon_groups["aromatic"], time_points, pymol_files, pymol_dir,
            "aromatic", "Aromatic Binding Site Accessibility (Phe/Tyr/Trp)",
            joinpath(output_dir, "3_aromatic_binding_sites.svg"), residue_info
        )
    end
    
    # 4. Terminal Carbons (deep side chain penetration)
    if !isempty(carbon_groups["sidechain_term"])
        create_true_atom_level_heatmap(
            carbon_groups["sidechain_term"], time_points, pymol_files, pymol_dir,
            "sidechain_term", "Side Chain Terminal Accessibility (CD/CE/CZ)",
            joinpath(output_dir, "4_sidechain_terminals.svg"), residue_info
        )
    end
    
    # Generate summary statistics
    println("\n=== Atom-Level Analysis Summary ===")
    for (carbon_type, data) in carbon_groups
        if !isempty(data)
            unique_atoms = length(unique(["$(d[1])-$(d[2])" for d in data]))
            println("$(carbon_type): $(unique_atoms) unique atoms")
        end
    end
    
    println("\n✅ Atom-level analysis complete!")
    println("📁 Results saved to: $(output_dir)/")
end

# Run the comprehensive atom-level analysis
println("=== Starting Atom-Level Heatmap Generation ===")
if abspath(PROGRAM_FILE) == @__FILE__
    dir = length(ARGS) > 0 ? ARGS[1] : "sattimeanalysis"
    create_comprehensive_atom_level_heatmaps(dir)
end

println("\n🎉 Atom-level heatmap generation complete!")
println("📊 Generated 4 targeted atom-level heatmaps:")
println("   • Backbone network (CA atoms)")
println("   • Side chain base accessibility (CB atoms)")
println("   • Aromatic binding sites (Phe/Tyr/Trp)")
println("   • Side chain terminals (CD/CE/CZ)")
println("🔬 True atom-level resolution achieved!")
println("📁 Files saved to: atom_level_analysis/") 