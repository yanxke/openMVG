// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_TWO_POINT_IMU_HPP
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_TWO_POINT_IMU_HPP

#include "openMVG/numeric/numeric.h"
#include <vector>

namespace openMVG {
namespace essential {
namespace kernel {

/**
 * @brief Two-point Essential Matrix solver given known relative rotation from IMU
 * 
 * When the relative rotation R between two cameras is known (e.g., from IMU),
 * the essential matrix E = [t]_× · R requires only solving for translation t.
 * Since E has scale ambiguity, t has 2 DOF and can be solved with 2 point correspondences.
 * 
 * For each correspondence (x₁, x₂), the epipolar constraint gives:
 *   x₂ᵀ · E · x₁ = 0
 *   x₂ᵀ · [t]_× · R · x₁ = 0
 *   t · (x₂ × (R · x₁)) = 0
 * 
 * This is a linear constraint on t. With 2 correspondences, we get 2 constraints
 * and can solve for t (up to scale).
 * 
 * @note This solver assumes:
 *   - R is known and accurate (from IMU or other source)
 *   - Input points are bearing vectors (normalized)
 *   - The solver returns up to 1 solution
 */
struct TwoPointImuSolver
{
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 1 };

  /**
   * @brief Solve for Essential matrix given 2 bearing vector correspondences and known R
   * @param x1 3×2 bearing vectors from camera 1
   * @param x2 3×2 bearing vectors from camera 2
   * @param R_relative Known relative rotation: R = R₂ · R₁ᵀ
   * @param models Output Essential matrices (E = [t]_× · R)
   */
  static void Solve(
    const Mat3X & x1,
    const Mat3X & x2,
    const Mat3 & R_relative,
    std::vector<Mat3> * models)
  {
    assert(x1.cols() == 2);
    assert(x2.cols() == 2);
    assert(x1.rows() == 3);
    assert(x2.rows() == 3);

    models->clear();

    // Compute R · x₁ for both points
    const Vec3 Rx1_0 = R_relative * x1.col(0);
    const Vec3 Rx1_1 = R_relative * x1.col(1);

    // Compute v_i = x₂ᵢ × (R · x₁ᵢ)
    // t must be orthogonal to both v₀ and v₁
    const Vec3 v0 = x2.col(0).cross(Rx1_0);
    const Vec3 v1 = x2.col(1).cross(Rx1_1);

    // t = v₀ × v₁ (perpendicular to both constraints)
    const Vec3 t = v0.cross(v1);

    // Check if solution is degenerate (v₀ parallel to v₁)
    const double t_norm = t.norm();
    if (t_norm < 1e-10)
      return;  // Degenerate case: points don't constrain translation

    // Normalize t (essential matrix has scale ambiguity)
    const Vec3 t_normalized = t / t_norm;

    // Construct E = [t]_× · R
    Mat3 t_cross;
    t_cross <<  0,              -t_normalized(2),  t_normalized(1),
                t_normalized(2),  0,              -t_normalized(0),
               -t_normalized(1),  t_normalized(0),  0;

    const Mat3 E = t_cross * R_relative;

    models->push_back(E);
  }

  /**
   * @brief Solve for Essential matrix (convenience wrapper without R parameter)
   * 
   * @note This overload exists for compatibility but will abort.
   *       Use the version with R_relative parameter instead.
   */
  static void Solve(
    const Mat3X & x1,
    const Mat3X & x2,
    std::vector<Mat3> * models)
  {
    assert(false && "TwoPointImuSolver requires R_relative parameter");
    models->clear();
  }
};

} // namespace kernel
} // namespace essential
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_TWO_POINT_IMU_HPP
