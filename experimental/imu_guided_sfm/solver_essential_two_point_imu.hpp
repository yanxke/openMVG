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

struct TwoPointImuSolver
{
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 1 };

  static void Solve(
    const Mat3X & x1,
    const Mat3X & x2,
    const Mat3 & R_relative,
    std::vector<Mat3> * models)
  {
    assert(x1.cols() == 2);
    assert(x2.cols() == 2);
    models->clear();

    const Vec3 Rx1_0 = R_relative * x1.col(0);
    const Vec3 Rx1_1 = R_relative * x1.col(1);
    const Vec3 v0 = x2.col(0).cross(Rx1_0);
    const Vec3 v1 = x2.col(1).cross(Rx1_1);
    const Vec3 t = v0.cross(v1);

    const double t_norm = t.norm();
    if (t_norm < 1e-10)
      return;

    const Vec3 t_normalized = t / t_norm;
    Mat3 t_cross;
    t_cross <<  0,              -t_normalized(2),  t_normalized(1),
                t_normalized(2),  0,              -t_normalized(0),
               -t_normalized(1),  t_normalized(0),  0;

    models->push_back(t_cross * R_relative);
  }
};

} // namespace kernel
} // namespace essential
} // namespace openMVG

#endif
