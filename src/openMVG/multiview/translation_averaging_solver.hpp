// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_SOLVER_HPP
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_SOLVER_HPP

#include "openMVG/multiview/translation_averaging_common.hpp"

#include <vector>

//------------------
//-- Bibliography --
//------------------
//- [1] "Robust Global Translations with 1DSfM."
//- Authors: Kyle Wilson and Noah Snavely.
//- Date: September 2014.
//- Conference: ECCV.
//
//- [2] "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."
//- Authors: Pierre Moulon, Pascal Monasse and Renaud Marlet.
//- Date: December 2013.
//- Conference: ICCV.

namespace openMVG {

/// Implementation of [1] : "5. Solving the Translations Problem" equation (3)
/// Compute camera center positions from relative camera translations (translation directions).
bool
solve_translations_problem_l2_chordal
(
  const int* edges,
  const double* poses,
  const double* weights,
  int num_edges,
  double loss_width,
  double* X,
  double function_tolerance,
  double parameter_tolerance,
  int max_iterations
);

/**
* @brief Registration of relative translations to global translations. It implements LInf minimization of  [2]
*  as a SoftL1 minimization. It can use group of relative translation vectors (bearing, or n-uplets of translations).
*  All relative motions must be 1 connected component.
*
* @param[in] vec_initial_estimates group of relative motion information
*             Each group will have its own optimized scale
*             Bearing: 2 view estimates => essential matrices)
*             N-Uplets: N-view estimates => i.e. 3 view estimations means a triplet of relative motion

* @param[out] translations found global camera translations
* @param[in] d_l1_loss_threshold optional threshold for SoftL1 loss (-1: no loss function)
* @return True if the registration can be solved
*/
bool
solve_translations_problem_softl1
(
  const std::vector<openMVG::RelativeInfo_Vec > & vec_initial_estimates,
  std::vector<Eigen::Vector3d> & translations,
  const double d_l1_loss_threshold = 0.01
);

/**
* @brief Registration of relative translations to global translations with distance constraints.
*  Same as solve_translations_problem_softl1, but adds maximum distance constraints between
*  specified pairs of cameras (e.g., for enforcing speed limits based on time).
*
* @param[in] vec_initial_estimates group of relative motion information
* @param[in] global_rotations global rotations (aligned with pose indices) used to convert
*            translations to camera centers for distance constraints
* @param[in] max_distance_constraints map of (pose_i, pose_j) -> max_distance
*             Only penalizes if ||C_i - C_j|| > max_distance
* @param[in] distance_constraint_weight weight for distance constraint residuals (default: 10.0)
* @param[out] translations found global camera translations
* @param[in] d_l1_loss_threshold optional threshold for SoftL1 loss (-1: no loss function)
* @return True if the registration can be solved
*/
bool
solve_translations_problem_softl1_with_constraints
(
  const std::vector<openMVG::RelativeInfo_Vec > & vec_initial_estimates,
  const std::vector<Mat3> & global_rotations,
  const Hash_Map<Pair, double> & max_distance_constraints,
  std::vector<Eigen::Vector3d> & translations,
  const double distance_constraint_weight = 10.0,
  const double d_l1_loss_threshold = 0.01
);

} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_SOLVER_HPP
