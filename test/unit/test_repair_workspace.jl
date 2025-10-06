"""
    test_repair_workspace.jl

Unit tests for repair_workspace.jl
Tests RepairWorkspace initialization, transaction lifecycle, face/node operations,
and rollback functionality.
"""

using Test

@testset "RepairWorkspace Initialization" begin
    @testset "Minimal Workspace Creation" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        @test ws.original_file == "synthetic_cube.nas"
        @test length(ws.working_nodes) > 0
        @test length(ws.working_faces) > 0
        @test ws.max_node_id > 0
        @test isempty(ws.modifications)
        @test ws.transaction_active == false
        @test ws.checkpoint_mod_count == 0
        @test ws.faces_added == 0
        @test ws.faces_deleted == 0
        @test ws.nodes_added == 0
    end
    
    @testset "Workspace Properties" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # Check that nodes are properly indexed
        for (node_id, coords) in ws.working_nodes
            @test node_id > 0
            @test length(coords) == 3
            @test all(isfinite.(coords))
        end
        
        # Check that faces reference valid nodes
        for (pid, faces) in ws.working_faces
            for face in faces
                @test length(face) == 3  # Triangular faces
                for node_id in face
                    @test haskey(ws.working_nodes, node_id)
                end
            end
        end
    end
end

@testset "Transaction Lifecycle - Commit" begin
    @testset "Begin Transaction" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        @test ws.transaction_active == false
        begin_transaction!(ws)
        @test ws.transaction_active == true
        @test ws.checkpoint_mod_count == 0
    end
    
    @testset "Commit Empty Transaction" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        initial_mod_count = length(ws.modifications)
        
        commit_transaction!(ws)
        @test ws.transaction_active == false
        @test ws.checkpoint_mod_count == initial_mod_count
    end
    
    @testset "Commit Transaction with Modifications" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        initial_mod_count = length(ws.modifications)
        
        # Make some modifications
        new_node_id = add_node!(ws, (10.0, 20.0, 30.0))
        @test ws.nodes_added == 1
        @test length(ws.modifications) == initial_mod_count + 1
        
        # Commit
        commit_transaction!(ws)
        @test ws.transaction_active == false
        @test ws.checkpoint_mod_count == length(ws.modifications)
        @test haskey(ws.working_nodes, new_node_id)
    end
    
    @testset "Commit Without Active Transaction" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        @test_throws ErrorException commit_transaction!(ws)
    end
end

@testset "Transaction Lifecycle - Rollback" begin
    @testset "Rollback Node Addition" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        initial_node_count = length(ws.working_nodes)
        initial_max_id = ws.max_node_id
        
        begin_transaction!(ws)
        checkpoint = ws.checkpoint_mod_count
        
        # Add a node
        new_node_id = add_node!(ws, (1.0, 2.0, 3.0))
        @test length(ws.working_nodes) == initial_node_count + 1
        @test ws.nodes_added == 1
        
        # Rollback
        rollback_transaction!(ws)
        
        # Verify state restored
        @test ws.transaction_active == false
        @test length(ws.working_nodes) == initial_node_count
        @test length(ws.modifications) == checkpoint
        @test ws.nodes_added == 0
        @test !haskey(ws.working_nodes, new_node_id)
    end
    
    @testset "Rollback Face Addition" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # Get a PID with faces
        pid = first(keys(ws.working_faces))
        initial_face_count = length(ws.working_faces[pid])
        
        # Get some valid node IDs
        node_ids = collect(keys(ws.working_nodes))[1:3]
        
        begin_transaction!(ws)
        checkpoint = ws.checkpoint_mod_count
        
        # Add a face
        add_face!(ws, pid, node_ids)
        @test length(ws.working_faces[pid]) == initial_face_count + 1
        @test ws.faces_added == 1
        
        # Rollback
        rollback_transaction!(ws)
        
        # Verify state restored
        @test length(ws.working_faces[pid]) == initial_face_count
        @test ws.faces_added == 0
        @test length(ws.modifications) == checkpoint
    end
    
    @testset "Rollback Face Deletion" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # Get a PID with faces
        pid = first(keys(ws.working_faces))
        initial_faces = deepcopy(ws.working_faces[pid])
        initial_face_count = length(initial_faces)
        
        begin_transaction!(ws)
        
        # Delete a face
        delete_face!(ws, pid, 1)
        @test length(ws.working_faces[pid]) == initial_face_count - 1
        @test ws.faces_deleted == 1
        
        # Rollback
        rollback_transaction!(ws)
        
        # Verify face restored
        @test length(ws.working_faces[pid]) == initial_face_count
        @test ws.faces_deleted == 0
        @test ws.working_faces[pid] == initial_faces
    end
    
    @testset "Rollback Multiple Modifications" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
    initial_state = Nas2StepTestUtils.capture_workspace_state(ws)
        
        begin_transaction!(ws)
        
        # Make multiple modifications
        node1 = add_node!(ws, (1.0, 0.0, 0.0))
        node2 = add_node!(ws, (2.0, 0.0, 0.0))
        node3 = add_node!(ws, (3.0, 0.0, 0.0))
        
        pid = first(keys(ws.working_faces))
        add_face!(ws, pid, [node1, node2, node3])
        
        @test ws.nodes_added == 3
        @test ws.faces_added == 1
        
        # Rollback all
        rollback_transaction!(ws)
        
        # Verify complete restoration
    @test Nas2StepTestUtils.workspace_state_equals(ws, initial_state)
    end
    
    @testset "Rollback Without Active Transaction" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        @test_throws ErrorException rollback_transaction!(ws)
    end
end

@testset "Nested Transaction Prevention" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
    
    begin_transaction!(ws)
    
    # Attempting nested transaction should error
    @test_throws ErrorException begin_transaction!(ws)
    
    # Clean up
    commit_transaction!(ws)
end

@testset "Node Operations" begin
    @testset "Add Node" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        initial_count = length(ws.working_nodes)
        initial_max_id = ws.max_node_id
        
        coords = (5.5, 6.6, 7.7)
        new_node_id = add_node!(ws, coords)
        
        @test new_node_id == initial_max_id + 1
        @test ws.max_node_id == new_node_id
        @test haskey(ws.working_nodes, new_node_id)
        @test ws.working_nodes[new_node_id] == coords
        @test length(ws.working_nodes) == initial_count + 1
        @test ws.nodes_added == 1
        
        commit_transaction!(ws)
    end
    
    @testset "Add Multiple Nodes" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        node_ids = Int[]
        for i in 1:5
            coords = (Float64(i), Float64(i*2), Float64(i*3))
            node_id = add_node!(ws, coords)
            push!(node_ids, node_id)
        end
        
        @test length(node_ids) == 5
        @test allunique(node_ids)
        @test ws.nodes_added == 5
        
        # Verify all nodes added
        for node_id in node_ids
            @test haskey(ws.working_nodes, node_id)
        end
        
        commit_transaction!(ws)
    end
end

@testset "Face Operations" begin
    @testset "Add Face with Valid Nodes" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        pid = first(keys(ws.working_faces))
        initial_count = length(ws.working_faces[pid])
        
        # Use existing nodes
        node_ids = collect(keys(ws.working_nodes))[1:3]
        add_face!(ws, pid, node_ids)
        
        @test length(ws.working_faces[pid]) == initial_count + 1
        @test ws.faces_added == 1
        @test ws.working_faces[pid][end] == node_ids
        
        commit_transaction!(ws)
    end
    
    @testset "Add Face with Missing Nodes" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        pid = first(keys(ws.working_faces))
        
        # Use invalid node IDs
        invalid_nodes = [9999, 10000, 10001]
        
        @test_throws ErrorException add_face!(ws, pid, invalid_nodes)
        
        rollback_transaction!(ws)
    end
    
    @testset "Delete Face - Valid Index" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        pid = first(keys(ws.working_faces))
        initial_count = length(ws.working_faces[pid])
        deleted_face = ws.working_faces[pid][1]
        
        delete_face!(ws, pid, 1)
        
        @test length(ws.working_faces[pid]) == initial_count - 1
        @test ws.faces_deleted == 1
        
        commit_transaction!(ws)
    end
    
    @testset "Delete Face - Out of Bounds" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        pid = first(keys(ws.working_faces))
        face_count = length(ws.working_faces[pid])
        
        # Try to delete non-existent face
        @test_throws ErrorException delete_face!(ws, pid, face_count + 10)
        @test_throws ErrorException delete_face!(ws, pid, 0)
        @test_throws ErrorException delete_face!(ws, pid, -1)
        
        rollback_transaction!(ws)
    end
end

@testset "Modification Tracking" begin
    @testset "Modification Count" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        @test length(ws.modifications) == 0
        
        add_node!(ws, (1.0, 2.0, 3.0))
        @test length(ws.modifications) == 1
        
        add_node!(ws, (4.0, 5.0, 6.0))
        @test length(ws.modifications) == 2
        
        pid = first(keys(ws.working_faces))
        node_ids = collect(keys(ws.working_nodes))[1:3]
        add_face!(ws, pid, node_ids)
        @test length(ws.modifications) == 3
        
        commit_transaction!(ws)
    end
    
    @testset "Modification Types" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        # Add node
        add_node!(ws, (1.0, 2.0, 3.0))
        @test ws.modifications[end].mod_type == NODE_ADDITION
        
        # Add face
        pid = first(keys(ws.working_faces))
        node_ids = collect(keys(ws.working_nodes))[1:3]
        add_face!(ws, pid, node_ids)
        @test ws.modifications[end].mod_type == FACE_ADDITION
        
        # Delete face
        delete_face!(ws, pid, 1)
        @test ws.modifications[end].mod_type == FACE_DELETION
        
        commit_transaction!(ws)
    end
end

@testset "Statistics Tracking" begin
    @testset "Face Statistics" begin
    ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        @test ws.faces_added == 0
        @test ws.faces_deleted == 0
        
        pid = first(keys(ws.working_faces))
        node_ids = collect(keys(ws.working_nodes))[1:3]
        
        # Add faces
        add_face!(ws, pid, node_ids)
        add_face!(ws, pid, node_ids)
        @test ws.faces_added == 2
        
        # Delete face
        delete_face!(ws, pid, 1)
        @test ws.faces_deleted == 1
        @test ws.faces_added == 2  # Should not change
        
        commit_transaction!(ws)
    end
    
    @testset "Node Statistics" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        begin_transaction!(ws)
        
        @test ws.nodes_added == 0
        
        for i in 1:5
            add_node!(ws, (Float64(i), 0.0, 0.0))
        end
        
        @test ws.nodes_added == 5
        
        commit_transaction!(ws)
    end
end

@testset "Checkpoint Management" begin
    @testset "Manual Checkpoint" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        @test ws.checkpoint_mod_count == 0
        
        begin_transaction!(ws)
        add_node!(ws, (1.0, 2.0, 3.0))
        commit_transaction!(ws)
        
        @test ws.checkpoint_mod_count == 1
    end
    
    @testset "Multiple Checkpoints" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # First transaction
        begin_transaction!(ws)
        add_node!(ws, (1.0, 0.0, 0.0))
        commit_transaction!(ws)
        checkpoint1 = ws.checkpoint_mod_count
        
        # Second transaction
        begin_transaction!(ws)
        add_node!(ws, (2.0, 0.0, 0.0))
        commit_transaction!(ws)
        checkpoint2 = ws.checkpoint_mod_count
        
        @test checkpoint2 > checkpoint1
        @test checkpoint2 == 2
    end
end

@testset "Complex Transaction Scenarios" begin
    @testset "Add and Delete Same Face" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        pid = first(keys(ws.working_faces))
        initial_count = length(ws.working_faces[pid])
        
        begin_transaction!(ws)
        
        node_ids = collect(keys(ws.working_nodes))[1:3]
        add_face!(ws, pid, node_ids)
        
        new_count = length(ws.working_faces[pid])
        @test new_count == initial_count + 1
        
        # Delete the face we just added
        delete_face!(ws, pid, new_count)
        
        @test length(ws.working_faces[pid]) == initial_count
        @test ws.faces_added == 1
        @test ws.faces_deleted == 1
        
        commit_transaction!(ws)
    end
    
    @testset "Sequential Transactions" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # Transaction 1
        begin_transaction!(ws)
        node1 = add_node!(ws, (1.0, 0.0, 0.0))
        commit_transaction!(ws)
        
        @test haskey(ws.working_nodes, node1)
        
        # Transaction 2
        begin_transaction!(ws)
        node2 = add_node!(ws, (2.0, 0.0, 0.0))
        commit_transaction!(ws)
        
        @test haskey(ws.working_nodes, node2)
        @test ws.nodes_added == 2
    end
    
    @testset "Rollback Then Commit" begin
        ws = Nas2StepTestUtils.create_minimal_workspace()
        
        # First transaction - rollback
        begin_transaction!(ws)
        add_node!(ws, (1.0, 0.0, 0.0))
        rollback_transaction!(ws)
        
        @test ws.nodes_added == 0
        
        # Second transaction - commit
        begin_transaction!(ws)
        node2 = add_node!(ws, (2.0, 0.0, 0.0))
        commit_transaction!(ws)
        
        @test ws.nodes_added == 1
        @test haskey(ws.working_nodes, node2)
    end
end


