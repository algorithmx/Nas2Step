Overview of Mismatch Classification
===

- Goal: Given an edge that exists on one side of an interface (A or B) but not on the other, determine why the meshes disagree and how they disagree topologically/geometrically.
- Inputs:
  - Edge: a geometric edge represented by EdgeKey((x1,y1,z1),(x2,y2,z2)) with canonical ordering (n1 <= n2).
  - InterfaceTopology: geometry/topology of both sides, including triangles, edge maps, and shared vertices.
- Output: An EdgeMismatch record with mismatch_type (one of T_JUNCTION, REFINEMENT, DIAGONAL, QUAD_MISMATCH, BOUNDARY_EDGE, NON_MANIFOLD, UNSHARED_ENDPOINT, UNKNOWN) and auxiliary data (e.g., hanging_nodes, quad_vertices).

Geometric primitives and tolerances
- EdgeKey canonicalizes endpoints so (a,b) == (b,a) for set/dict lookups.
- Coordinates are compared with a consistent small tolerance (default tol=1e-4).
- “Hanging node” detection: A target-side vertex that lies on the source edge segment (colinear and between endpoints within tol).

Data used in classification
- Source faces vs target faces: determined by present_in (:A or :B).
- Target nodes: all vertices from the opposite side; used to search for hanging nodes on the source edge.
- Edge incidence in source: number of source triangles that contain the edge.
- Shared vertex set: vertices common to both sides (after rounding/tolerance).

Decision flow (classify_edge_mismatch)
1) Hanging nodes on the opposite side
   - Find target-side vertices lying on the source edge.
   - If exactly 1 hanging node → T_JUNCTION.
   - If more than 1 hanging node → REFINEMENT.

2) No hanging nodes: topological checks on the source side
   - Count how many source triangles contain the edge:
     - Exactly 1 → BOUNDARY_EDGE (edge is on a boundary of the source mesh).
     - More than 2 → NON_MANIFOLD (invalid; edge shared by >2 triangles).
     - Zero → UNKNOWN (anomaly: edge not used by any source triangle).
     - Exactly 2 → proceed.

3) Endpoints shared?
   - Compute shared vertices across both sides.
   - If either endpoint is not in the shared set → UNSHARED_ENDPOINT.
     - Note on precedence: The classifier checks BOUNDARY_EDGE before UNSHARED_ENDPOINT. If an edge is both boundary and has unshared endpoints, it will be classified as BOUNDARY_EDGE. Tests should avoid this ambiguity (make the edge internal—used by 2 source triangles—when specifically testing UNSHARED_ENDPOINT).

4) Quad/diagonal analysis (when both endpoints are shared and the edge is internal on the source side)
   - Attempt to form a quad around the edge using the two adjacent source triangles.
   - Search for a corresponding pair of target-side triangles that use the opposite diagonal.
   - Compare the four boundary edges of the quads:
     - If boundaries match and a valid opposite diagonal is found → DIAGONAL.
     - If the four quad vertices are found but the boundary edges do not match → QUAD_MISMATCH.
     - If no consistent quad/opposite diagonal can be found → UNKNOWN.

Categories with first-principles examples (from tests)
- T_JUNCTION: One mesh has a vertex lying on the other mesh’s edge (exactly one hanging node). Example: A has diagonal (0,0)-(10,10); B has a center vertex (5,5) splitting that diagonal.
- REFINEMENT: Multiple hanging nodes on the edge (hierarchical subdivision). Example: A has long diagonal; B splits it with two internal points.
- DIAGONAL: Same quad boundary on both sides, but different diagonals (A uses AC, B uses BD).
- QUAD_MISMATCH: Same 4 vertices exist on both sides, but the boundary edges of the quad differ, so it’s not a pure diagonal flip.
- BOUNDARY_EDGE: Edge appears in only one triangle on the source side (not an internal interface edge).
- NON_MANIFOLD: Edge is shared by more than two triangles on the source side.
- UNSHARED_ENDPOINT: One or both edge endpoints are not in the shared vertex set (ensure the edge is internal to avoid boundary precedence).
- UNKNOWN: Catch-all for cases that do not fit the above patterns.

Practical notes from tests and debugging
- Building test topology: Ensure the InterfaceTopology fields are assigned correctly: edges_only_in_A, edges_only_in_B, then edges_shared (in that order). Misordering these leads to misleading membership tests.
- Ambiguities: Because classification prioritizes boundary detection before unshared endpoints, tests that assert UNSHARED_ENDPOINT should construct the source edge as internal (shared by exactly two source triangles) so it is not classified as BOUNDARY_EDGE first.
- Tolerance and rounding: Consistent rounding and isapprox checks are used throughout for robust geometric comparisons.

Minimal checklist for interpreting a mismatch
- Does the opposite side have 1 or more vertices lying on the edge? If yes → T_JUNCTION or REFINEMENT.
- If no, how many source triangles use the edge? 1 → BOUNDARY_EDGE; >2 → NON_MANIFOLD; 0 → UNKNOWN; 2 → continue.
- Are both endpoints shared between sides? If not → UNSHARED_ENDPOINT.
- Can you form a consistent quad and does the other side use the opposite diagonal with matching boundary? If yes → DIAGONAL; if boundaries differ → QUAD_MISMATCH; else → UNKNOWN.