"""
    test_interface_topology.jl

Unit tests for interface_topology.jl
Tests EdgeKey, Triangle, BoundingBox, build_interface_topology, and compute_interface_area.
"""

using Test

# Import CoordinateKeys functions for consistent EdgeKey handling
using Nas2Step: create_edge_key_int, coordinate_key_int

@testset "EdgeKey Construction and Hashing" begin
    @testset "EdgeKey Ordering" begin
        p1 = (0.0, 0.0, 0.0)
        p2 = (1.0, 1.0, 1.0)

        # Test that order doesn't matter
        ek1 = create_edge_key_int(p1, p2)
        ek2 = create_edge_key_int(p2, p1)
        @test ek1 == ek2
        @test hash(ek1) == hash(ek2)

        # Test internal ordering (smaller first)
        @test ek1.node1 <= ek1.node2
    end
    
    @testset "EdgeKey Hash Consistency" begin
        p1 = (0.0, 0.0, 0.0)
        p2 = (1.0, 1.0, 1.0)

        ek1 = create_edge_key_int(p1, p2)
        ek2 = create_edge_key_int(p1, p2)
        ek3 = create_edge_key_int(p2, p1)

        # All should have the same hash
        @test hash(ek1) == hash(ek2)
        @test hash(ek1) == hash(ek3)
        @test hash(ek2) == hash(ek3)

        # And all should be equal
        @test ek1 == ek2
        @test ek1 == ek3
        @test ek2 == ek3
    end
    
    @testset "EdgeKey Different Edges" begin
        p1 = (0.0, 0.0, 0.0)
        p2 = (1.0, 0.0, 0.0)
        p3 = (0.0, 1.0, 0.0)

        ek1 = create_edge_key_int(p1, p2)
        ek2 = create_edge_key_int(p1, p3)
        ek3 = create_edge_key_int(p2, p3)

        # Different edges should not be equal
        @test ek1 != ek2
        @test ek1 != ek3
        @test ek2 != ek3

        # And should have different hashes (with high probability)
        @test hash(ek1) != hash(ek2)
        @test hash(ek1) != hash(ek3)
    end
end

@testset "Triangle Geometric Calculations" begin
    @testset "Equilateral Triangle" begin
        c1 = (0.0, 0.0, 0.0)
        c2 = (1.0, 0.0, 0.0)
        c3 = (0.5, sqrt(3)/2, 0.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Test area (equilateral with side 1: area = sqrt(3)/4)
        expected_area = sqrt(3) / 4
        @test isapprox(tri.area, expected_area, rtol=1e-6)
        
        # Test centroid
        expected_centroid = (0.5, sqrt(3)/6, 0.0)
        @test isapprox(tri.centroid[1], expected_centroid[1], rtol=1e-6)
        @test isapprox(tri.centroid[2], expected_centroid[2], rtol=1e-6)
        @test isapprox(tri.centroid[3], expected_centroid[3], rtol=1e-6)
        
        # Test normal (should point in +Z direction for counter-clockwise vertices in XY plane)
        @test abs(tri.normal[3]) ≈ 1.0
        @test abs(tri.normal[1]) < 1e-10
        @test abs(tri.normal[2]) < 1e-10
    end
    
    @testset "Right Triangle" begin
        # Right triangle with legs of length 3 and 4, hypotenuse 5
        c1 = (0.0, 0.0, 0.0)
        c2 = (3.0, 0.0, 0.0)
        c3 = (0.0, 4.0, 0.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Area should be (3 * 4) / 2 = 6
        @test isapprox(tri.area, 6.0, rtol=1e-6)
        
        # Centroid at (1, 4/3, 0)
        @test isapprox(tri.centroid[1], 1.0, rtol=1e-6)
        @test isapprox(tri.centroid[2], 4.0/3.0, rtol=1e-6)
        @test isapprox(tri.centroid[3], 0.0, rtol=1e-6)
    end
    
    @testset "Degenerate Triangle (Collinear Points)" begin
        c1 = (0.0, 0.0, 0.0)
        c2 = (1.0, 0.0, 0.0)
        c3 = (2.0, 0.0, 0.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Area should be near zero
        @test tri.area < 1e-10
        
        # Normal magnitude should be near zero
        normal_mag = sqrt(tri.normal[1]^2 + tri.normal[2]^2 + tri.normal[3]^2)
        @test normal_mag < 1e-10 || normal_mag ≈ 0.0
    end
    
    @testset "Triangle in 3D Space" begin
        # Non-planar triangle (actually it's planar, but not in XY plane)
        c1 = (0.0, 0.0, 0.0)
        c2 = (1.0, 0.0, 1.0)
        c3 = (0.0, 1.0, 1.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Area should be positive
        @test tri.area > 0
        
        # Normal should be unit length (or zero for degenerate)
        normal_mag = sqrt(tri.normal[1]^2 + tri.normal[2]^2 + tri.normal[3]^2)
        @test isapprox(normal_mag, 1.0, rtol=1e-6) || normal_mag < 1e-10
    end
    
    @testset "Very Small Triangle" begin
        # Tiny triangle, but not degenerate
        scale = 1e-6
        c1 = (0.0, 0.0, 0.0)
        c2 = (scale, 0.0, 0.0)
        c3 = (scale/2, scale*sqrt(3)/2, 0.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Area should be positive but small
        expected_area = sqrt(3) / 4 * scale^2
        @test isapprox(tri.area, expected_area, rtol=1e-4)
    end
end

@testset "BoundingBox Construction" begin
    @testset "BoundingBox from Multiple Points" begin
        points = [
            (0.0, 0.0, 0.0),
            (1.0, 2.0, 3.0),
            (-1.0, 1.0, 2.0),
            (0.5, -0.5, 1.5)
        ]
        bbox = BoundingBox(points)
        
        # Check min corner
        @test bbox.min_corner[1] == -1.0
        @test bbox.min_corner[2] == -0.5
        @test bbox.min_corner[3] == 0.0
        
        # Check max corner
        @test bbox.max_corner[1] == 1.0
        @test bbox.max_corner[2] == 2.0
        @test bbox.max_corner[3] == 3.0
    end
    
    @testset "BoundingBox from Single Point" begin
        points = [(1.5, 2.5, 3.5)]
        bbox = BoundingBox(points)
        
        # Min and max should be the same
        @test bbox.min_corner == bbox.max_corner
        @test bbox.min_corner == (1.5, 2.5, 3.5)
    end
    
    @testset "BoundingBox from Empty Set" begin
        points = NTuple{3,Float64}[]
        bbox = BoundingBox(points)
        
        # Should return zero bounding box
        @test bbox.min_corner == (0.0, 0.0, 0.0)
        @test bbox.max_corner == (0.0, 0.0, 0.0)
    end
    
    @testset "BoundingBox Contains All Points" begin
    points = [Nas2StepTestUtils.random_point_3d() for _ in 1:20]
        bbox = BoundingBox(points)
        
        # All points should be inside the bounding box
        for p in points
            @test p[1] >= bbox.min_corner[1]
            @test p[1] <= bbox.max_corner[1]
            @test p[2] >= bbox.min_corner[2]
            @test p[2] <= bbox.max_corner[2]
            @test p[3] >= bbox.min_corner[3]
            @test p[3] <= bbox.max_corner[3]
        end
    end
end

@testset "InterfaceTopology Structure" begin
    @testset "Basic Topology Properties" begin
        # Create a minimal topology manually for testing
        shared_nodes = Set{NTuple{3,Float64}}([
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0)
        ])
        
        node_map = Dict{NTuple{3,Float64}, Tuple{Int,Int}}(
            (0.0, 0.0, 0.0) => (1, 101),
            (1.0, 0.0, 0.0) => (2, 102)
        )
        
        # Create some test triangles
        tri1 = Triangle(1, 2, 3, 1, (0.0,0.0,0.0), (1.0,0.0,0.0), (0.5,1.0,0.0))
        tri2 = Triangle(101, 102, 103, 2, (0.0,0.0,0.0), (1.0,0.0,0.0), (0.5,1.0,0.0))
        
        faces_A = [tri1]
        faces_B = [tri2]
        
        # Create edge maps
        edges_A = Dict{EdgeKey, Vector{Int}}()
        edges_B = Dict{EdgeKey, Vector{Int}}()
        
        ek1 = create_edge_key_int((0.0,0.0,0.0), (1.0,0.0,0.0))
        edges_A[ek1] = [1]
        edges_B[ek1] = [1]
        
        edges_shared = Set([ek1])
        edges_only_A = Set{EdgeKey}()
        edges_only_B = Set{EdgeKey}()
        
        bbox = BoundingBox([(0.0,0.0,0.0), (1.0,0.0,0.0), (0.5,1.0,0.0)])
        
        topology = InterfaceTopology(
            1, 2,
            shared_nodes, node_map,
            faces_A, faces_B,
            edges_A, edges_B,
            edges_only_A, edges_only_B, edges_shared,
            bbox,
            2, 1, 1, 1, 1, 1.0,
            0.0, 0.0, 0, 1.0  # consistency metrics: max_vertex_dist, mean_vertex_dist, edge_mismatch_count, triangulation_similarity
        )
        
        # Test basic properties
        @test topology.pidA == 1
        @test topology.pidB == 2
        @test topology.total_shared_nodes == 2
        @test topology.total_faces_A == 1
        @test topology.total_faces_B == 1
        @test topology.conformity_ratio == 1.0
    end
end

@testset "Coordinate Rounding and Tolerance" begin
    @testset "EdgeKey with Nearly Equal Coordinates" begin
        # Points that differ by less than tolerance
        p1 = (0.0, 0.0, 0.0)
        p2 = (0.0, 0.0, 1e-9)  # Very close to p1

        # After scaling to integer coordinates, these both become (0,0,0)
        # This is the intended behavior of the coordinate scaling system
        ek1 = create_edge_key_int(p1, p1)
        ek2 = create_edge_key_int(p2, p2)

        # They should be the same after scaling (both become (0,0,0))
        @test ek1 == ek2

        # Test with a larger difference that should remain different
        p3 = (0.0, 0.0, 1e-3)  # Still small, but should scale to (0,0,10)
        ek3 = create_edge_key_int(p3, p3)
        @test ek1 != ek3
    end
    
    @testset "Triangle Coordinate Consistency" begin
        # Test that triangle coordinates are stored correctly
        c1 = (0.123456789, 0.0, 0.0)
        c2 = (1.0, 0.0, 0.0)
        c3 = (0.5, 0.866025403, 0.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        @test tri.coord1 == c1
        @test tri.coord2 == c2
        @test tri.coord3 == c3
    end
end

@testset "Edge Cases and Error Handling" begin
    @testset "Triangle with Identical Nodes" begin
        # Triangle where all three vertices are the same
        c = (1.0, 1.0, 1.0)
        tri = Triangle(1, 1, 1, 1, c, c, c)
        
        # Area should be zero
        @test tri.area == 0.0
        
        # Centroid should be at the point
        @test tri.centroid == c
    end
    
    @testset "Triangle with Two Identical Nodes" begin
        # Triangle where two vertices are the same (degenerate)
        c1 = (0.0, 0.0, 0.0)
        c2 = (0.0, 0.0, 0.0)
        c3 = (1.0, 1.0, 1.0)
        tri = Triangle(1, 2, 3, 1, c1, c2, c3)
        
        # Area should be zero or very small
        @test tri.area < 1e-10
    end
end


