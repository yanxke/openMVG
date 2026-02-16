// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PLY_HELPER_H
#define OPENMVG_SFM_PLY_HELPER_H

#include "openMVG/numeric/numeric.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace openMVG{
namespace plyHelper{

inline std::vector<unsigned char> ComputeKeepMaskByMedianDistance(const std::vector<Vec3> & vec_points)
{
  const size_t point_count = vec_points.size();
  std::vector<unsigned char> keep_mask(point_count, 1);
  const size_t remove_count = static_cast<size_t>(std::floor(point_count * 0.01));
  if (point_count == 0 || remove_count == 0)
  {
    return keep_mask;
  }

  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> zs;
  xs.reserve(point_count);
  ys.reserve(point_count);
  zs.reserve(point_count);
  for (const Vec3 & p : vec_points)
  {
    xs.push_back(p(0));
    ys.push_back(p(1));
    zs.push_back(p(2));
  }

  const auto median_of = [](std::vector<double> & v) -> double
  {
    const size_t n = v.size();
    const size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    if (n % 2 == 1)
    {
      return v[mid];
    }
    const double high = v[mid];
    std::nth_element(v.begin(), v.begin() + (mid - 1), v.begin() + mid);
    return 0.5 * (v[mid - 1] + high);
  };

  const Vec3 median_center(
    median_of(xs),
    median_of(ys),
    median_of(zs));

  std::vector<double> dist_sq;
  dist_sq.reserve(point_count);
  for (const Vec3 & p : vec_points)
  {
    dist_sq.push_back((p - median_center).squaredNorm());
  }

  const size_t keep_count = point_count - remove_count;
  std::vector<double> dist_sq_copy = dist_sq;
  std::nth_element(
    dist_sq_copy.begin(),
    dist_sq_copy.begin() + (keep_count - 1),
    dist_sq_copy.end());
  const double keep_threshold = dist_sq_copy[keep_count - 1];

  for (size_t i = 0; i < point_count; ++i)
  {
    keep_mask[i] = (dist_sq[i] <= keep_threshold) ? 1 : 0;
  }
  return keep_mask;
}

/// Export 3D point vector to PLY format
inline
bool
exportToPly
(
  const std::vector<Vec3> & vec_points,
  const std::string & sFileName
)
{
  const std::vector<unsigned char> keep_mask = ComputeKeepMaskByMedianDistance(vec_points);
  size_t kept_count = 0;
  for (unsigned char keep : keep_mask)
  {
    kept_count += keep;
  }

  std::ofstream outfile(sFileName.c_str());
  if (!outfile)
    return false;

  outfile << "ply"
    << "\n" << "format ascii 1.0"
    << "\n" << "element vertex " << kept_count
    << "\n" << "property double x"
    << "\n" << "property double y"
    << "\n" << "property double z"
    << "\n" << "property uchar red"
    << "\n" << "property uchar green"
    << "\n" << "property uchar blue"
    << "\n" << "end_header" << "\n";

  outfile << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  for (size_t i=0; i < vec_points.size(); ++i)
  {
    if (!keep_mask[i])
      continue;
    outfile
      << vec_points[i](0) << ' '
      << vec_points[i](1) << ' '
      << vec_points[i](2) << ' '
      << "255 255 255" << "\n";
  }
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

/// Export 3D point vector and camera position to PLY format
inline bool exportToPly
(
  const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName,
  const std::vector<Vec3> * vec_coloredPoints = nullptr
)
{
  const std::vector<unsigned char> keep_mask = ComputeKeepMaskByMedianDistance(vec_points);
  size_t kept_count = 0;
  for (unsigned char keep : keep_mask)
  {
    kept_count += keep;
  }

  std::ofstream outfile(sFileName.c_str());
  if (!outfile)
    return false;

  outfile << "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << kept_count+vec_camPos.size()
    << '\n' << "property double x"
    << '\n' << "property double y"
    << '\n' << "property double z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << "\n";

  outfile << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  for (size_t i=0; i < vec_points.size(); ++i)  {
    if (!keep_mask[i])
      continue;
    if (vec_coloredPoints == nullptr)
      outfile
        << vec_points[i](0) << ' '
        << vec_points[i](1) << ' '
        << vec_points[i](2) << ' '
        << "255 255 255\n";
    else
      outfile
        << vec_points[i](0) << ' '
        << vec_points[i](1) << ' '
        << vec_points[i](2) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](0)) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](1)) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](2))
        << "\n";
  }

  for (size_t i=0; i < vec_camPos.size(); ++i)  {
    outfile
      << vec_camPos[i](0) << ' '
      << vec_camPos[i](1) << ' '
      << vec_camPos[i](2) << ' '
      << "0 255 0\n";
  }
  outfile.flush();
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

} // namespace plyHelper
} // namespace openMVG

#endif // OPENMVG_SFM_PLY_HELPER_H
