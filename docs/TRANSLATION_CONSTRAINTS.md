# Time-Based Translation Distance Constraints Implementation

## Overview
This implementation adds support for maximum translation distance constraints between camera pairs during the translation averaging phase of the SfM pipeline. This allows enforcing motion priors such as speed limits based on capture time.

## Files Modified

### 1. Core Data Structures

#### `/home/yke/src/openMVG/src/openMVG/sfm/sfm_data.hpp`
- Added `Hash_Map<Pair, double> max_translation_distance_` to `SfM_Data` struct
- Stores pairwise maximum distance constraints (in meters)
- Maps (view_i, view_j) -> max_distance_meters

#### `/home/yke/src/openMVG/src/openMVG/sfm/sfm_data_io_cereal.cpp`
- Added serialization/deserialization for `max_translation_distance_` field
- Backward compatible - optional field with try-catch
- Logs number of constraints loaded

### 2. Translation Averaging Solver

#### `/home/yke/src/openMVG/src/openMVG/multiview/translation_averaging_solver.hpp`
- Added new function declaration: `solve_translations_problem_softl1_with_constraints()`
- Accepts `max_distance_constraints` map and `distance_constraint_weight` parameter

####`/home/yke/src/openMVG/src/openMVG/multiview/translation_averaging_solver_softl1.cpp`
- **New Cost Functor:** `MaxDistanceError`
  - Computes distance between two camera centers: `dist = ||t_i - t_j||`
  - One-sided penalty: `residual = weight * max(0, dist - max_distance)`
  - Only penalizes when distance exceeds the limit (allows stopping)

- **New Function:** `solve_translations_problem_softl1_with_constraints()`
  - Extends the standard SoftL1 solver
  - Adds distance constraint residual blocks for each constrained pair
  - Checks if both poses exist in the problem before adding constraint
  - Logs how many constraints were successfully added

### 3. Global SfM Pipeline Integration (NEEDS MANUAL COMPLETION)

#### `/home/yke/src/openMVG/src/openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.cpp`
**TODO:** The TRANSLATION_AVERAGING_SOFTL1 case needs to be updated to:
1. Check if `sfm_data.max_translation_distance_` is non-empty
2. If constraints exist, remap them from original pose IDs to reindexed IDs
3. Call `solve_translations_problem_softl1_with_constraints()` instead of the standard solver
4. Otherwise, use the standard solver

**Target location:** Line 247-266

**Required code pattern:**
```cpp
case TRANSLATION_AVERAGING_SOFTL1:
{
  std::vector<Vec3> vec_translations;
  
  // Check if we have max distance constraints
  if (!sfm_data.max_translation_distance_.empty())
  {
    // Remap constraints from original to reindexed pose IDs
    Hash_Map<Pair, double> reindexed_constraints;
    for (const auto & constraint : sfm_data.max_translation_distance_)
    {
      const Pair & orig_pair = constraint.first;
      if (reindex_forward.count(orig_pair.first) && 
          reindex_forward.count(orig_pair.second))
      {
        const IndexT new_i = reindex_forward.at(orig_pair.first);
        const IndexT new_j = reindex_forward.at(orig_pair.second);
        reindexed_constraints[Pair(new_i, new_j)] = constraint.second;
      }
    }
    
    OPENMVG_LOG_INFO << "Using SOFTL1 solver with " << reindexed_constraints.size()
                     << " max distance constraints.";
    
    if (!solve_translations_problem_softl1_with_constraints(
          vec_relative_motion_cpy, reindexed_constraints, vec_translations, 10.0))
    {
      OPENMVG_LOG_ERROR << "TRANSLATION_AVERAGING_SOFTL1 with constraints failed";
      return false;
    }
  }
  else
  {
    // Standard solver without constraints
    if (!solve_translations_problem_softl1(vec_relative_motion_cpy, vec_translations))
    {
      OPENMVG_LOG_ERROR << "TRANSLATION_AVERAGING_SOFTL1 failed";
      return false;
    }
  }

  // Update poses with the solution
  for (size_t i = 0; i < iNview; ++i)
  {
    const Vec3 & t = vec_translations[i];
    const IndexT pose_id = reindex_backward[i];
    const Mat3 & Ri = map_globalR.at(pose_id);
    sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
  }
}
break;
```

## Python Script

### `/home/yke/src/openMVG/scripts/populate_translation_constraints.py`
Offline preprocessing script to populate constraints from EXIF timestamps.

**Features:**
- Extracts millisecond-precision timestamps from EXIF (DateTimeOriginal + SubSecTimeOriginal)
- Computes `max_distance = max_velocity * delta_time` for temporally nearby frames
- Generates pairwise constraints within a configurable time window
- Outputs modified sfm_data.json with constraints

**Usage:**
```bash
python populate_translation_constraints.py \
  --sfm_data matches/sfm_data.json \
  --output matches/sfm_data_with_constraints.json \
  --max_velocity 1.0 \
  --time_window 10.0
```

**Parameters:**
- `--max_velocity`: Maximum speed in m/s (default: 1.0)
- `--time_window`: Maximum time between frames to constrain (default: 10.0s)

**Dependencies:**
```bash
pip install exifread
```

## Workflow

1. **Extract timestamps and compute constraints:**
   ```bash
   python scripts/populate_translation_constraints.py \
     --sfm_data matches/sfm_data.json \
     --output matches/sfm_data_constrained.json \
     --max_velocity 1.0  # 1 meter per second max speed
   ```

2. **Run global SfM with constraints:**
   ```bash
   openMVG_main_SfM \
     -i matches/sfm_data_constrained.json \
     -m matches \
     -o reconstruction \
     -t 2  # Use SOFTL1 translation averaging
   ```

3. The solver will automatically detect and use the constraints during translation averaging.

## Configuration

### Constraint Weight
The `distance_constraint_weight` parameter (default: 10.0) controls how strongly constraints are enforced relative to the relative translation residuals.

- **Higher weight** (e.g., 50.0): Stronger enforcement, may sacrifice fit to relative translations
- **Lower weight** (e.g., 1.0): Weaker enforcement, prioritizes relative translation fit
- **Recommended**: 10.0 provides good balance

### Time Window  
The Python script only generates constraints for pairs within the time window. This prevents:
- Unnecessary constraints between distant frames
- Overly constraining the problem

**Recommended:** 5-10 seconds for walking sequences

## Scale Handling

**Important:** The constraints assume the SfM reconstruction is in meters. This happens automatically if:

1. **GPS Priors exists:** The ViewPriors with pose_center already define metric scale
2. **No GPS:** The scale is arbitrary. Options:
   - Use the median speed approach (compute constraints relative to observed motion)
   - Add at least one known distance or GPS point
   - Scale thewhole reconstruction post-hoc

The constraints work in whatever unit system the reconstruction uses, but the `max_velocity` parameter must be specified in those units.

## Technical Notes

### Why Pairwise Constraints in SfM_Data?
- More general than timestamps (works for any motion prior)
- Decouples constraint computation from the solver
- Allows manual annotation or other constraint sources
- Can be pre-computed offline

### Why One-Sided Penalty?
- Allows the person to be stationary (no minimum distance)
- Only penalizes physically impossible motion (too fast)
- Soft constraint: can be violated if relative translations strongly disagree

### Reindexing
The translation averaging solver reindexes poses from [min_id, max_id] to [0, N-1] for efficiency. The constraints must be remapped accordingly before passing to the solver.

## Testing

TODO: Add unit tests for:
- MaxDistanceError cost functor
- solve_translations_problem_softl1_with_constraints()
- Python script timestamp extraction
- Constraint serialization/deserialization
