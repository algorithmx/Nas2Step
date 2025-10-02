using Nas2Step

# Alternative 6-tet decomposition of hex that ensures conformal interfaces
# Uses three diagonals through the hex in a consistent pattern
# Hex nodes: bottom face (1,2,3,4), top face (5,6,7,8)
# Grid positions: (i,j), (i+1,j), (i+1,j+1), (i,j+1) per layer
#   4---3    8---7
#   |   |    |   |
#   1---2    5---6
#
# Key insight: Use a single diagonal through the ENTIRE hex volume
# and decompose faces consistently relative to that diagonal
function hex_to_6tets(n1, n2, n3, n4, n5, n6, n7, n8)
    # Use the main body diagonal from n1 to n7
    # Decompose into 6 tets that all share this diagonal edge
    return [
        (n1, n2, n3, n7),  # tet containing face 1-2-3
        (n1, n3, n4, n7),  # tet containing face 1-3-4
        (n1, n5, n6, n7),  # tet containing face 1-5-6
        (n1, n6, n2, n7),  # tet containing face 1-6-2
        (n1, n4, n8, n7),  # tet containing face 1-4-8
        (n1, n8, n5, n7),  # tet containing face 1-8-5
    ]
end

function demo_box()
    nx, ny = 5, 5
    nz_per_region = 3
    xs = LinRange(0, 1, nx)
    ys = LinRange(0, 1, ny)

    # Create grid points with rough interface
    bottom_pts = [(x, y, 0.0) for y in ys for x in xs]
    # Add random perturbations to the interface z-coordinate
    common_pts = [(x, y, 1.0 + (rand() - 0.5) * 0.2) for y in ys for x in xs]
    top_pts = [(x, y, 2.0) for y in ys for x in xs]

    all_pts = []
    # Region 1 layers
    for k in 0:nz_per_region
        layer_pts = [
            (p_bot[1], p_bot[2], p_bot[3] * (1 - k/nz_per_region) + p_com[3] * (k/nz_per_region))
            for (p_bot, p_com) in zip(bottom_pts, common_pts)
        ]
        append!(all_pts, layer_pts)
    end
    # Region 2 layers (skip common layer)
    for k in 1:nz_per_region
        layer_pts = [
            (p_com[1], p_com[2], p_com[3] * (1 - k/nz_per_region) + p_top[3] * (k/nz_per_region))
            for (p_com, p_top) in zip(common_pts, top_pts)
        ]
        append!(all_pts, layer_pts)
    end

    function create_layer_tets(bottom_layer_idx)
        offset_bot = bottom_layer_idx * nx * ny
        offset_top = (bottom_layer_idx + 1) * nx * ny
        tets = []
        
        for j in 1:(ny-1), i in 1:(nx-1)
            # Hex corners in standard order
            n1 = offset_bot + (j - 1) * nx + i
            n2 = offset_bot + (j - 1) * nx + (i + 1)
            n3 = offset_bot + j * nx + (i + 1)
            n4 = offset_bot + j * nx + i
            n5 = offset_top + (j - 1) * nx + i
            n6 = offset_top + (j - 1) * nx + (i + 1)
            n7 = offset_top + j * nx + (i + 1)
            n8 = offset_top + j * nx + i
            
            # Create 6 tets
            hex_tets = hex_to_6tets(n1, n2, n3, n4, n5, n6, n7, n8)
            append!(tets, hex_tets)
        end
        return tets
    end

    # Region 1
    tets1 = []
    for k in 0:(nz_per_region-1)
        append!(tets1, create_layer_tets(k))
    end
    reg1 = Nas2Step.VolumeRegion("Region_1", all_pts, tets1)

    # Region 2
    tets2 = []
    for k in nz_per_region:(2*nz_per_region-1)
        append!(tets2, create_layer_tets(k))
    end
    reg2 = Nas2Step.VolumeRegion("Region_2", all_pts, tets2)

    out = Nas2Step.write_nas_volume("demo_box.nas", [reg1, reg2])
    println("Wrote $out with $(length(tets1) + length(tets2)) tetrahedra")
    println("  Region 1: $(length(tets1)) tetrahedra")
    println("  Region 2: $(length(tets2)) tetrahedra")
end

demo_box()
