// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_FILTERS_HPP
#define OPENMVG_SFM_SFM_DATA_FILTERS_HPP

#include <set>

#include "openMVG/types.hpp"

namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

struct PixelResidualAdaptiveStats
{
  IndexT old_logic_outlier_observations = 0;
  IndexT rescued_observations = 0;
  IndexT removed_observations = 0;
};

struct AngleErrorAdaptiveStats
{
  IndexT old_logic_outlier_tracks = 0;
  IndexT rescued_tracks = 0;
  IndexT removed_tracks = 0;
};

/// List the view indexes that have valid camera intrinsic and pose.
std::set<IndexT> Get_Valid_Views
(
  const SfM_Data & sfm_data
);

/// Filter a list of pair: Keep only the pairs that are defined by the index list
template <typename IterablePairs, typename IterableIndex>
Pair_Set Pair_filter
(
  const IterablePairs & pairs,
  const IterableIndex & index
)
{
  Pair_Set kept_pairs;
  for (auto& it : pairs)
  {
    if (index.count(it.first) > 0 &&
      index.count(it.second) > 0)
    kept_pairs.insert(it);
  }
  return kept_pairs;
}

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
IndexT RemoveOutliers_PixelResidualError
(
  SfM_Data & sfm_data,
  const double dThresholdPixel,
  const unsigned int minTrackLength = 2
);

// Remove outlier observations by pixel residual with a relaxed threshold for tracks
// that contain at least one close-in-time observation pair.
IndexT RemoveOutliers_PixelResidualErrorAdaptive
(
  SfM_Data & sfm_data,
  const double dThresholdPixel,
  const double dThresholdPixelCloseTime,
  const double closeTimeSeconds,
  PixelResidualAdaptiveStats * stats = nullptr,
  const unsigned int minTrackLength = 2
);

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
IndexT RemoveOutliers_AngleError
(
  SfM_Data & sfm_data,
  const double dMinAcceptedAngle
);

// Remove tracks by triangulation angle with a relaxed threshold for close-in-time view pairs.
// A track is kept if:
// - max angle over all obs pairs >= dMinAcceptedAngle, OR
// - it contains at least one close-time obs pair and that pair angle >= dMinAcceptedAngleCloseTime.
IndexT RemoveOutliers_AngleErrorAdaptive
(
  SfM_Data & sfm_data,
  const double dMinAcceptedAngle,
  const double dMinAcceptedAngleCloseTime,
  const double closeTimeSeconds,
  AngleErrorAdaptiveStats * stats = nullptr
);

/// Erase pose with insufficient track observations
bool eraseMissingPoses
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose = 6
);

/// Erase observations with no defined pose
bool eraseObservationsWithMissingPoses
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_landmark = 2
);

/// Remove unstable content from analysis of the sfm_data structure
bool eraseUnstablePosesAndObservations
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose = 6,
  const IndexT min_points_per_landmark = 2
);

// Adaptive cleanup:
// - Uses min_points_per_landmark_default for normal tracks.
// - Uses min_points_per_landmark_close for tracks that contain at least one close-time view pair.
bool eraseUnstablePosesAndObservationsAdaptive
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose,
  const IndexT min_points_per_landmark_default,
  const IndexT min_points_per_landmark_close,
  const double closeTimeSeconds
);

/// Tell if the sfm_data structure is one CC or not
bool IsTracksOneCC
(
  const SfM_Data & sfm_data
);

/// Keep the largest connected component of tracks from the sfm_data structure
void KeepLargestViewCCTracks
(
  SfM_Data & sfm_data
);

/**
* @brief Implement a statistical Structure filter that remove 3D points that have:
* - a depth that is too large (threshold computed as factor * median ~= X84)
* @param sfm_data The sfm scene to filter (inplace filtering)
* @param k_factor The factor applied to the median depth per view
* @param k_min_point_per_pose Keep only poses that have at least this amount of points
* @param k_min_track_length Keep only tracks that have at least this length
* @return The min_median_value observed for all the view
*/
double DepthCleaning
(
  SfM_Data & sfm_data,
  const double k_factor =  5.2,    //5.2 * median ~= X84,
  const IndexT k_min_point_per_pose = 12,  // 6 min
  const IndexT k_min_track_length = 2      // 2 min
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_FILTERS_HPP
