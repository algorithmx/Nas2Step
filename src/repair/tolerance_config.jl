# Tolerance Configuration for Mesh Repair
# Centralized tolerance values to ensure consistency across all mesh repair operations

"""
Global tolerance configuration for mesh conformity analysis and repair.

All geometric comparisons (vertex matching, edge detection, etc.) should use
these standardized tolerance values to ensure consistency.
"""
module ToleranceConfig

# ============================================================================
# Primary Geometric Tolerance
# ============================================================================

"""
DEFAULT_GEOMETRIC_TOLERANCE: The standard tolerance for all geometric comparisons.

This value (1e-4) is used for:
- Vertex/node equality checks
- Edge endpoint matching  
- Hanging node detection
- Quad vertex identification
- All distance-based geometric operations

Corresponds to 0.0001 units in the mesh coordinate system.
For typical engineering meshes in mm or inches, this provides adequate precision
while being robust to floating-point round-off errors.
"""
const DEFAULT_GEOMETRIC_TOLERANCE = 1e-4

"""
COORDINATE_ROUNDING_DIGITS: Number of decimal places for coordinate rounding.

This value (4 digits) is used when coordinates need to be rounded for:
- Hash key generation
- Set membership tests using rounded coordinates
- JSON export formatting

Note: 4 digits corresponds to 1e-4 precision, matching DEFAULT_GEOMETRIC_TOLERANCE.
"""
const COORDINATE_ROUNDING_DIGITS = 4

# ============================================================================
# Quality and Angle Tolerances  
# ============================================================================

"""
DEFAULT_MIN_ANGLE_DEG: Minimum acceptable triangle angle in degrees.

Used for quality validation during repair planning.
"""
const DEFAULT_MIN_ANGLE_DEG = 10.0

"""
DEFAULT_MAX_ANGLE_DEG: Maximum acceptable triangle angle in degrees.

Used for quality validation during repair planning.
"""
const DEFAULT_MAX_ANGLE_DEG = 170.0

"""
DEFAULT_MAX_ASPECT_RATIO: Maximum acceptable triangle aspect ratio.

Used for quality validation during repair planning.
"""
const DEFAULT_MAX_ASPECT_RATIO = 20.0

"""
DEFAULT_FEATURE_ANGLE_DEG: Default feature angle for boundary detection.

Used to identify sharp features vs. smooth boundaries.
"""
const DEFAULT_FEATURE_ANGLE_DEG = 30.0

# ============================================================================
# Spatial Clustering
# ============================================================================

"""
DEFAULT_CLUSTER_DISTANCE: Maximum distance for grouping nearby mismatched vertices.

Used in vertex conformity analysis to identify spatial patterns in mismatches.
"""
const DEFAULT_CLUSTER_DISTANCE = 10.0

# ============================================================================
# Helper Functions
# ============================================================================

"""
    tolerance_squared()

Return the squared value of the default tolerance.
Useful for distance² comparisons to avoid sqrt() operations.
"""
tolerance_squared() = DEFAULT_GEOMETRIC_TOLERANCE * DEFAULT_GEOMETRIC_TOLERANCE

"""
    round_coordinate(value)

Round a single coordinate value using the standard precision.
"""
round_coordinate(value::Real) = round(value, digits=COORDINATE_ROUNDING_DIGITS)

"""
    round_coordinates(coords::NTuple{3,Float64})

Round a 3D coordinate tuple using the standard precision.
"""
function round_coordinates(coords::NTuple{3,Float64})::NTuple{3,Float64}
    return (round_coordinate(coords[1]), 
            round_coordinate(coords[2]), 
            round_coordinate(coords[3]))
end

end # module ToleranceConfig
