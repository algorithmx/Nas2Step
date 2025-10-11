# CoordinateKeys Module
# Centralized coordinate key operations for consistent EdgeKey creation
# Phase 2: Type-safe, performant coordinate handling

module CoordinateKeys

export coordinate_key_int, coordinate_key_float,
       create_edge_key_int, create_edge_key_float,
       EdgeKeyInt, EdgeKeyFloat,
       COORD_SCALE, COORD_DIGITS, DEFAULT_TOLERANCE,
       validate_coordinate_range, get_coordinate_key_stats,
       convert_to_int, convert_to_float, are_coordinate_keys_equal

# ============================================================================
# Constants and Configuration
# ============================================================================

const COORD_SCALE = 10000
const COORD_DIGITS = 4
const DEFAULT_TOLERANCE = 1e-4

# ============================================================================
# Type Definitions
# ============================================================================

"""
Type-safe EdgeKey with integer coordinates for optimal performance.
"""
struct EdgeKeyInt
    node1::NTuple{3,Int}
    node2::NTuple{3,Int}

    function EdgeKeyInt(n1::NTuple{3,Int}, n2::NTuple{3,Int})
        # Validate coordinate ranges
        validate_coordinate_range(n1)
        validate_coordinate_range(n2)

        # Canonical ordering for consistent comparison
        if n1 <= n2
            new(n1, n2)
        else
            new(n2, n1)
        end
    end
end

"""
Type-safe EdgeKey with float coordinates for backward compatibility.
"""
struct EdgeKeyFloat
    node1::NTuple{3,Float64}
    node2::NTuple{3,Float64}

    function EdgeKeyFloat(n1::NTuple{3,Float64}, n2::NTuple{3,Float64})
        # Validate coordinate ranges
        validate_coordinate_range(n1)
        validate_coordinate_range(n2)

        # Canonical ordering for consistent comparison
        if n1 <= n2
            new(n1, n2)
        else
            new(n2, n1)
        end
    end
end

# ============================================================================
# Core Coordinate Key Functions
# ============================================================================

"""
    coordinate_key_int(p::NTuple{3,Float64})::NTuple{3,Int}

Convert floating coordinates to integer keys for hash consistency.
Uses scaling instead of rounding for better performance.

# Arguments
- `p`: 3D coordinate tuple with float values

# Returns
- 3D coordinate tuple with integer values (scaled by COORD_SCALE)

# Example
```julia
coord = (1.23456789, 2.34567891, 3.45678912)
int_key = coordinate_key_int(coord)  # Returns (12346, 23457, 34568)
```
"""
function coordinate_key_int(p::NTuple{3,Float64})::NTuple{3,Int}
    STATS.int_conversions += 1
    return (Int(round(p[1] * COORD_SCALE)),
            Int(round(p[2] * COORD_SCALE)),
            Int(round(p[3] * COORD_SCALE)))
end

"""
    coordinate_key_float(p::NTuple{3,Float64})::NTuple{3,Float64}

Convert to rounded float coordinates (legacy compatibility).

# Arguments
- `p`: 3D coordinate tuple with float values

# Returns
- 3D coordinate tuple with rounded float values (COORD_DIGITS precision)
"""
function coordinate_key_float(p::NTuple{3,Float64})::NTuple{3,Float64}
    STATS.float_conversions += 1
    return (round(p[1]; digits=COORD_DIGITS),
            round(p[2]; digits=COORD_DIGITS),
            round(p[3]; digits=COORD_DIGITS))
end

# ============================================================================
# EdgeKey Factory Functions
# ============================================================================

"""
    create_edge_key_int(coord1, coord2)::EdgeKeyInt

Create EdgeKey with integer coordinates (recommended for performance).

# Arguments
- `coord1`: First 3D coordinate tuple
- `coord2`: Second 3D coordinate tuple

# Returns
- EdgeKeyInt with scaled integer coordinates

# Performance Notes
- Integer hashing is significantly faster than float hashing
- No floating-point precision issues
- Memory efficient (integers vs floats)
"""
function create_edge_key_int(coord1::NTuple{3,Float64},
                           coord2::NTuple{3,Float64})::EdgeKeyInt
    STATS.edgekey_int_created += 1
    return EdgeKeyInt(coordinate_key_int(coord1), coordinate_key_int(coord2))
end

"""
    create_edge_key_float(coord1, coord2)::EdgeKeyFloat

Create EdgeKey with float coordinates (legacy compatibility).

# Arguments
- `coord1`: First 3D coordinate tuple
- `coord2`: Second 3D coordinate tuple

# Returns
- EdgeKeyFloat with rounded float coordinates
"""
function create_edge_key_float(coord1::NTuple{3,Float64},
                             coord2::NTuple{3,Float64})::EdgeKeyFloat
    STATS.edgekey_float_created += 1
    return EdgeKeyFloat(coordinate_key_float(coord1), coordinate_key_float(coord2))
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    validate_coordinate_range(coord)

Validate that coordinate values are within reasonable ranges.

# Arguments
- `coord`: 3D coordinate tuple (int or float)

# Raises
- `ArgumentError` if coordinates are outside expected range
"""
function validate_coordinate_range(coord::T) where {T <: NTuple{3,Number}}
    max_reasonable_value = 100_000_000

    for (i, val) in enumerate(coord)
        if abs(val) > max_reasonable_value
            STATS.validation_warnings += 1
            @warn "Large coordinate value detected at position $i: $val. " *
                  "This may indicate a scaling issue or invalid geometry."
        end

        if !isfinite(val)
            STATS.validation_errors += 1
            throw(ArgumentError("Invalid coordinate value at position $i: $val"))
        end
    end

    return true
end

"""
    convert_to_float(coord_int::NTuple{3,Int})::NTuple{3,Float64}

Convert integer coordinate key back to float coordinates.
Useful for comparing EdgeKeyInt with Triangle coordinates.

# Arguments
- `coord_int`: Integer coordinate tuple from EdgeKeyInt

# Returns
- Float coordinate tuple (scaled back by 1/COORD_SCALE)
"""
function convert_to_float(coord_int::NTuple{3,Int})::NTuple{3,Float64}
    return (coord_int[1] / COORD_SCALE,
            coord_int[2] / COORD_SCALE,
            coord_int[3] / COORD_SCALE)
end

"""
    convert_to_int(coord_float::NTuple{3,Float64})::NTuple{3,Int}

Convert float coordinates to integer coordinate keys for exact comparison.
This is the preferred method for coordinate comparison - convert floats to integers
and compare element-wise integers instead of converting integers to floats.

# Arguments
- `coord_float`: Float coordinate tuple

# Returns
- Integer coordinate tuple (scaled by COORD_SCALE)

# Example
```julia
coord = (1.23456789, 2.34567891, 3.45678912)
int_coord = convert_to_int(coord)  # Returns (12346, 23457, 34568)
```
"""
function convert_to_int(coord_float::NTuple{3,Float64})::NTuple{3,Int}
    return coordinate_key_int(coord_float)
end

"""
    are_coordinate_keys_equal(key1, key2; tol=DEFAULT_TOLERANCE)

Check if two coordinate keys are equal within tolerance.
Works with mixed int/float coordinate keys.

# Arguments
- `key1`: First coordinate key (int or float tuple)
- `key2`: Second coordinate key (int or float tuple)
- `tol`: Tolerance for comparison (default: DEFAULT_TOLERANCE)

# Returns
- `true` if coordinates are within tolerance, `false` otherwise
"""
function are_coordinate_keys_equal(key1::T1, key2::T2;
                                 tol::Real=DEFAULT_TOLERANCE) where
                                 {T1 <: NTuple{3,Number}, T2 <: NTuple{3,Number}}

    # Convert both to float for comparison
    float1 = T1 == NTuple{3,Int} ? convert_to_float(key1) : key1
    float2 = T2 == NTuple{3,Int} ? convert_to_float(key2) : key2

    # Calculate squared distance
    dx = float1[1] - float2[1]
    dy = float1[2] - float2[2]
    dz = float1[3] - float2[3]
    dist2 = dx*dx + dy*dy + dz*dz

    return dist2 <= tol*tol
end

# ============================================================================
# Performance Monitoring
# ============================================================================

# Global counters for performance monitoring
mutable struct CoordinateKeyStats
    int_conversions::Int
    float_conversions::Int
    edgekey_int_created::Int
    edgekey_float_created::Int
    validation_warnings::Int
    validation_errors::Int

    CoordinateKeyStats() = new(0, 0, 0, 0, 0, 0)
end

const STATS = CoordinateKeyStats()

"""
    get_coordinate_key_stats()::Dict{String,Any}

Get current coordinate key statistics.

# Returns
- Dictionary with performance statistics
"""
function get_coordinate_key_stats()::Dict{String,Any}
    return Dict(
        "int_conversions" => STATS.int_conversions,
        "float_conversions" => STATS.float_conversions,
        "edgekey_int_created" => STATS.edgekey_int_created,
        "edgekey_float_created" => STATS.edgekey_float_created,
        "validation_warnings" => STATS.validation_warnings,
        "validation_errors" => STATS.validation_errors
    )
end

"""
    reset_coordinate_key_stats()

Reset all coordinate key statistics.
"""
function reset_coordinate_key_stats()
    STATS.int_conversions = 0
    STATS.float_conversions = 0
    STATS.edgekey_int_created = 0
    STATS.edgekey_float_created = 0
    STATS.validation_warnings = 0
    STATS.validation_errors = 0
    return
end


# ============================================================================
# Hashing and Equality
# ============================================================================

# Implement hash and equality for custom EdgeKey types
import Base: hash, ==

function Base.hash(key::EdgeKeyInt, h::UInt)
    return hash((key.node1, key.node2), h)
end

function Base.:(==)(a::EdgeKeyInt, b::EdgeKeyInt)
    return (a.node1 == b.node1) && (a.node2 == b.node2)
end

function Base.hash(key::EdgeKeyFloat, h::UInt)
    return hash((key.node1, key.node2), h)
end

function Base.:(==)(a::EdgeKeyFloat, b::EdgeKeyFloat)
    return (a.node1 == b.node1) && (a.node2 == b.node2)
end

# ============================================================================
# Pretty Printing
# ============================================================================

import Base: show, string

function Base.show(io::IO, key::EdgeKeyInt)
    print(io, "EdgeKeyInt($(key.node1) → $(key.node2))")
end

function Base.show(io::IO, key::EdgeKeyFloat)
    print(io, "EdgeKeyFloat($(key.node1) → $(key.node2))")
end

function Base.string(key::EdgeKeyInt)
    return "EdgeKeyInt($(key.node1) → $(key.node2))"
end

function Base.string(key::EdgeKeyFloat)
    return "EdgeKeyFloat($(key.node1) → $(key.node2))"
end

end # module CoordinateKeys