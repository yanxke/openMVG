# Elevation Check - DEBUGGING NOTES

## Issue Discovered

The elevation check is **rejecting 40% of pairs**, which is way too aggressive for flat-ground walking. Expected rejection rate should be **5-15%** max.

## Current Status: **DISABLED**

```cpp
double max_elevation_ratio = 0.0  // Default: OFF
```

To re-enable for testing:
```cpp
// In main_GeometricFilterWithPriors.cpp, line ~540
GeometricFilter_EMatrix_AC_Imu(4.0, 64, &map_imu_rotations, 
                               map_PutativeMatches.size(), 25.0, 
                               0.5)  // Try 0.5 (30°) or 0.7 (45°)
```

---

## Root Cause Analysis

### Potential Issues

1. **2-Point Estimate is Noisy**
   - Using only 2 points from 5 samples gives unreliable t direction
   - Small matching errors → large t direction errors
   - May not represent true translation

2. **Sample Selection is Poor**
   - Evenly-spaced sampling might pick outliers
   - Should use RANSAC consensus or better sampling

3. **Coordinate System Confusion**
   - Y-axis: Is it really "down" in the expected sense?
   - Need to verify camera orientation convention

4. **Walking Pattern Has Vertical Component**
   - Camera may naturally move up/down while walking
   - Phone in hand/pocket may not be level

---

## Debugging Steps

### 1. Add Logging to See Actual Values

```cpp
// In E_ACRobust_Imu.hpp, around line 305
const double elevation_ratio = std::abs(t(1)) / t_norm;

// Add this:
static int debug_count = 0;
if (debug_count++ < 10) {  // Log first 10 pairs
  OPENMVG_LOG_INFO << "Elevation debug: t=(" << t(0) << "," << t(1) << "," << t(2) 
                   << "), ratio=" << elevation_ratio;
}

if (elevation_ratio > m_max_elevation_ratio) {
  // ...
}
```

### 2. Use More Samples for Consensus

Instead of 2 points, try multiple 2-point solutions:

```cpp
// Try 10 different 2-point combinations from 5 samples
std::vector<double> ratios;
for (size_t i = 0; i < sample_indices.size()-1; ++i) {
  for (size_t j = i+1; j < sample_indices.size(); ++j) {
    // Solve with points i,j
    // Extract t, compute ratio
    ratios.push_back(ratio);
  }
}
// Use median ratio instead of single estimate
std::sort(ratios.begin(), ratios.end());
double median_ratio = ratios[ratios.size()/2];
```

### 3. Check Camera Orientation

```cpp
// Print R_relative to see if it's reasonable
static int rot_debug = 0;
if (rot_debug++ < 5) {
  OPENMVG_LOG_INFO << "R_relative:\n" << R_relative;
}
```

### 4. Validate Against Ground Truth

If you have a few pairs you know should pass:
- Manually check their elevation ratios
- See if 0.3 is actually too strict

---

## Better Implementation Ideas

### Option A: Use All Samples with Least Squares

```cpp
// Instead of 2-point, use all 5 samples with least squares
// Minimize: sum ||x2_i x (R * x1_i)||^2 subject to ||t|| = 1
// This gives a more robust estimate
```

### Option B: Use RANSAC on Samples

```cpp
// Run mini-RANSAC on the 5 samples
// Take best 2-point solution that has most inliers
// More robust to outliers
```

### Option C: Stricter Pre-filtering

```cpp
// Only check elevation if pre-filtering was very successful
if (filtered_indices.size() > 0.8 * putative_matches.size()) {
  // High inlier ratio → rotation is good → elevation check makes sense
  // Otherwise skip it
}
```

### Option D: Use Different Coordinate

```cpp
// Instead of |t_y| / ||t||, use |t_y| / sqrt(t_x^2 + t_z^2)
// This compares vertical to horizontal directly
const double horiz_norm = std::sqrt(t(0)*t(0) + t(2)*t(2));
if (horiz_norm > 1e-6) {
  const double elevation_ratio = std::abs(t(1)) / horiz_norm;
  // Now ratio=0.2 means 20cm vertical per 1m horizontal
}
```

---

## Recommended Fix

### Short-term: Increase Threshold

```cpp
max_elevation_ratio = 0.7  // ~35 degrees, more permissive
```

### Medium-term: Better Estimation

Use **Option B** (RANSAC on samples):

```cpp
// Try C(5,2)=10 combinations of 2 points
// Count how many give elevation_ratio < threshold
// If >70% agree it's flat, accept; else reject

int flat_count = 0;
for (auto [i,j] : all_pairs) {
  // Solve 2-point
  // Extract t, compute ratio
  if (ratio < threshold) flat_count++;
}
if (flat_count < 0.7 * num_pairs) reject();
```

### Long-term: IMU Accelerometer Integration

Use IMU accelerometer to **measure actual vertical displacement**:
- Integrate acceleration over time
- Compare against geometric estimate
- More physically accurate

---

## Why It's Rejecting 40%

Likely reasons:
1. **Threshold too strict**: 0.3 (17°) may be too low
   - Natural walking bounce
   - Phone/camera not perfectly level
   - Uneven ground

2. **Noisy 2-point estimates**: 
   - 2 points insufficient for reliable t direction
   - Small pixel errors → large direction errors

3. **Coordinate transform issue**:
   - Y-axis might not be "down" as expected
   - Landscape orientation handling?

---

## Testing Protocol

1. **Start with disabled** (current state)
2. **Enable with high threshold** (0.7 or 1.0)
3. **Add logging** to see actual t values
4. **Gradually decrease** threshold based on logs
5. **Validate** on known-good sequences

---

## Example: Expected vs Observed

For flat walking, camera ~1.5m high, 1m forward per step:

**Expected**:
- Horizontal: (1.0, 0.0, 0.0) m
- Vertical: (0.0, 0.05, 0.0) m (walking bounce)
- Ratio: 0.05 / 1.0 = **0.05** ✓ (should pass)

**If rejecting 40%**, actual ratios might be:
- Median ratio: ~0.4 (noisy estimates?)
- Or: Many pairs actually have elevation change
- Or: Coordinate system wrong (Y not vertical?)

---

## Action Items

- [ ] Add debug logging to see actual t values
- [ ] Check if Y-axis is really "down" as expected
- [ ] Try threshold 0.5-0.7 instead of 0.3
- [ ] Implement multi-sample consensus (Option B)
- [ ] Validate on 10-20 known-good pairs manually

---

## Conclusion

**Current recommendation**: **Keep disabled (0.0)** until properly debugged.

The pre-filtering + 2-point solver already provides significant speedup (~4x). The elevation check was meant as an additional optimization, but it's currently over-aggressive.

Focus on getting the main IMU-guided filtering working well first. The elevation check can be added back later after proper debugging.
