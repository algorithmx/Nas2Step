# Enhanced Feasibility Assessment Reporting

## What Was Improved

The feasibility assessment report has been significantly enhanced to provide more detailed and actionable information about why repair plans fail, with special focus on constraint violations.

## New Features

### 1. Detailed Constraint Blocking Analysis

**Before:**
```
Constraint violations breakdown:
  No affected triangles found for edge insertion: 16
```

**After:**
```
Constraint Blocking Analysis:
  Total plans blocked by constraints: 16/16 (100.0%)
  Breakdown:
    - Only constraints (quality OK): 0
    - Both constraints and quality: 16

  Top Constraint Violations (by frequency):

    1. "No affected triangles found for edge insertion"
       Affects: 16 plans (100.0% of all, 100.0% of infeasible)
       ℹ️  Note: These plans also have quality issues, so fixing constraint alone won't help
```

### 2. Impact Analysis for Each Constraint

For each constraint violation, the report now shows:

- **Frequency**: How many plans are affected
- **Percentage of all plans**: What fraction of total plans hit this constraint
- **Percentage of infeasible plans**: What fraction of failed plans have this issue
- **Would-be-feasible count**: How many plans would become feasible if ONLY this constraint were removed

### 3. Actionable Insights with Visual Indicators

The system now identifies **critical** vs **non-critical** constraints:

**Critical Constraint (Actionable):**
```
⚠️  CRITICAL: 15 plans would become feasible if this constraint were resolved
    (these plans have good quality but only this constraint issue)
```

**Non-Critical Constraint:**
```
ℹ️  Note: These plans also have quality issues, so fixing constraint alone won't help
```

### 4. Impact Summary

When actionable constraints are found:

```
⭐ Impact Summary:
   Removing 2 constraint type(s) could unlock 45 additional feasible plans
   🎯 This would be SUFFICIENT to meet the feasibility threshold!
```

Or if not sufficient:

```
⭐ Impact Summary:
   Removing 1 constraint type(s) could unlock 10 additional feasible plans
   ⚠️  Still need 3 more feasible plans to meet threshold
```

## Implementation Details

### New Data Structure: `constraint_impact`

```julia
constraint_impact[reason] = (
    count = count,                      # Number of plans with this violation
    pct_of_infeasible = pct,           # % of infeasible plans
    pct_of_all = pct,                  # % of all plans
    would_be_feasible = count,         # Plans that would become feasible
    is_blocking = bool                 # Is this truly blocking progress?
)
```

### Logic for Identifying Actionable Constraints

A constraint is considered **actionable** (critical) if:
1. It affects plans where quality is acceptable (`quality_acceptable == true`)
2. It's the ONLY constraint violation for those plans (`length(constraint_violations) == 1`)

This means removing that single constraint would immediately unlock those plans for repair.

### Enhanced Breakdown Categories

Plans are now categorized more precisely:
- **Only quality issues**: Quality fails but no constraints
- **Only constraints**: Constraints fail but quality is OK ← **These are actionable!**
- **Both**: Both quality and constraints fail ← **Fixing constraints alone won't help**

## Example Output Analysis

### Interface 7/14: PID 3 ↔ PID 4

```
Constraint Blocking Analysis:
  Total plans blocked by constraints: 170/170 (100.0%)
  Breakdown:
    - Only constraints (quality OK): 0
    - Both constraints and quality: 170

  Top Constraint Violations (by frequency):

    1. "No affected triangles found for edge insertion"
       Affects: 164 plans (96.5% of all, 96.5% of infeasible)
       ℹ️  Note: These plans also have quality issues, so fixing constraint alone won't help

    2. "Cannot plan quad retriangulation"
       Affects: 6 plans (3.5% of all, 3.5% of infeasible)
       ℹ️  Note: These plans also have quality issues, so fixing constraint alone won't help
```

**Key Insights from this report:**
1. 100% of plans are blocked by constraints
2. But ALL of them also have quality issues
3. This means the fundamental problem is not the constraints but the quality/geometry
4. No actionable constraints to remove - need to address underlying issues

### What This Tells Us

For the NC_Reduction_4.nas mesh:
- **Primary issue**: "No affected triangles found" (96-100% of plans across interfaces)
- **Secondary issue**: "Cannot plan quad retriangulation" (3-5% for diagonal mismatches)
- **Root cause**: Both issues are coupled with quality problems
- **Recommendation**: Need to investigate why affected triangles aren't found, not just remove constraints

## Benefits

### 1. Clearer Diagnosis
Users can now immediately see:
- Which constraints are blocking the most repairs
- Whether fixing a constraint would actually help
- How many repairs could be unlocked

### 2. Prioritization
The report clearly indicates:
- **Critical** constraints worth investigating (⚠️ marker)
- **Non-critical** constraints that are symptoms of deeper issues (ℹ️ marker)

### 3. Actionable Next Steps
Instead of just saying "170 plans violate constraints", users now know:
- Exactly which constraint types are the problem
- Whether removing them would help (and by how much)
- If there are quality issues that must be addressed first

### 4. Better Understanding
The percentage breakdowns help users understand:
- How widespread each issue is
- Whether it's a universal problem or specific to certain edges
- The relationship between constraints and quality issues

## Code Changes

### Modified File
- `src/repair/repair_planning.jl`

### Changes Made

1. **Enhanced `assess_interface_feasibility()` function** (lines 167-212):
   - Added tracking of which plans hit each constraint
   - Compute detailed impact statistics
   - Identify actionable vs non-actionable constraints
   - Added `constraint_impact` to return value

2. **Enhanced printing section** (lines 750-800):
   - New "Constraint Blocking Analysis" section
   - Per-constraint impact reporting with visual indicators
   - Impact summary showing potential unlocks
   - Guidance on whether threshold could be met

### Backward Compatibility
✅ All existing functionality preserved
✅ Existing tests still pass
✅ New fields are additions, not replacements

## Future Enhancements

Potential further improvements:

1. **Constraint dependency analysis**: Show if multiple constraints co-occur
2. **Historical tracking**: Compare constraint patterns across interfaces
3. **Recommendation engine**: Suggest specific actions based on constraint patterns
4. **Visualization**: Generate charts showing constraint distribution
5. **Export detailed reports**: JSON output with full constraint analysis

## Summary

The enhanced feasibility assessment provides:
- ✅ More detailed constraint analysis
- ✅ Clear identification of actionable vs non-actionable issues
- ✅ Quantified impact of each constraint type
- ✅ Visual indicators for critical issues
- ✅ Actionable guidance for next steps

**Result:** Users can now make informed decisions about which constraints to investigate and whether addressing them will actually improve repair success rates.
