// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// Copyright (c) 2026
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/istreamwrapper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <unordered_map>
#include <memory>
#include <numeric>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;

namespace
{

struct Args
{
  std::string input_dir;
  double min_time_delta = 30.0;
  double max_homography_error = 5.0;
  double min_imu_rotation_error = 90.0;
  double max_bbox_area_ratio = 0.10;
  double expansion_time_window = 2.0;
  bool dry_run = false;
};

struct ViewInfo
{
  IndexT id_intrinsic = UndefinedIndexT;
  std::string filename;
  int width = 0;
  int height = 0;
  bool has_timestamp = false;
  double timestamp_s = 0.0;
  bool has_imu_rotation = false;
  Mat3 imu_rotation = Mat3::Identity();
};

struct PairStats
{
  bool has_pts = false;
  bool temporal = false;
  bool bbox = false;
  bool hom = false;
  bool imu = false;
  bool exclude = false;

  double time_delta = std::numeric_limits<double>::quiet_NaN();
  double bbox_i = 0.0;
  double bbox_j = 0.0;
  double hom_err = std::numeric_limits<double>::infinity();
  double imu_err = std::numeric_limits<double>::quiet_NaN();
};

struct ExclusionsV2
{
  std::set<Pair> manual;
  std::set<Pair> auto_pairs;
};

struct PairHash
{
  std::size_t operator()(const Pair & p) const noexcept
  {
    const std::size_t h1 = std::hash<IndexT>{}(p.first);
    const std::size_t h2 = std::hash<IndexT>{}(p.second);
    return h1 ^ (h2 * 2654435761ULL + 0x9e3779b9ULL + (h1 << 6) + (h1 >> 2));
  }
};

Mat3 DeviceToCameraRotation()
{
  Mat3 R;
  R << 1.0, 0.0, 0.0,
       0.0, -1.0, 0.0,
       0.0, 0.0, -1.0;
  return R;
}

Mat3 DeviceLandscapeLeftRemap()
{
  Mat3 R;
  R << 0.0, -1.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 0.0, 1.0;
  return R;
}

Mat3 QuaternionToMat3(const double qw, const double qx, const double qy, const double qz)
{
  const double norm = std::sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
  if (norm < 1e-12)
  {
    return Mat3::Identity();
  }
  const double w = qw / norm;
  const double x = qx / norm;
  const double y = qy / norm;
  const double z = qz / norm;

  Mat3 R;
  R << 1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),       2.0 * (x * z + y * w),
       2.0 * (x * y + z * w),       1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w),
       2.0 * (x * z - y * w),       2.0 * (y * z + x * w),       1.0 - 2.0 * (x * x + y * y);
  return R;
}

double RotationErrorDeg(const Mat3 & R_a, const Mat3 & R_b)
{
  const Mat3 R_err = R_a * R_b.transpose();
  const double trace_value = (R_err.trace() - 1.0) * 0.5;
  const double clamped = std::max(-1.0, std::min(1.0, trace_value));
  return std::acos(clamped) * 180.0 / M_PI;
}

bool ParseTimestampFromFilename(const std::string & filename, double & epoch_seconds)
{
  static const std::regex pattern("_(\\d{4})-(\\d{2})-(\\d{2})-(\\d{2})-(\\d{2})-(\\d{2})-(\\d{3})");
  std::smatch m;
  if (!std::regex_search(filename, m, pattern))
  {
    return false;
  }

  std::tm tm = {};
  tm.tm_year = std::stoi(m[1]) - 1900;
  tm.tm_mon = std::stoi(m[2]) - 1;
  tm.tm_mday = std::stoi(m[3]);
  tm.tm_hour = std::stoi(m[4]);
  tm.tm_min = std::stoi(m[5]);
  tm.tm_sec = std::stoi(m[6]);

  std::time_t t = 0;
#if defined(_WIN32)
  t = _mkgmtime(&tm);
#else
  t = timegm(&tm);
#endif
  if (t == static_cast<std::time_t>(-1))
  {
    return false;
  }

  epoch_seconds = static_cast<double>(t) + static_cast<double>(std::stoi(m[7])) / 1000.0;
  return true;
}

double BBoxAreaRatio(const std::vector<Vec2> & pts, const int width, const int height)
{
  if (pts.size() < 2 || width <= 0 || height <= 0)
  {
    return 0.0;
  }
  double x_min = pts.front()(0), x_max = pts.front()(0);
  double y_min = pts.front()(1), y_max = pts.front()(1);
  for (const Vec2 & p : pts)
  {
    x_min = std::min(x_min, p(0));
    x_max = std::max(x_max, p(0));
    y_min = std::min(y_min, p(1));
    y_max = std::max(y_max, p(1));
  }
  const double bbox_area = (x_max - x_min) * (y_max - y_min);
  const double img_area = static_cast<double>(width) * static_cast<double>(height);
  return img_area > 0.0 ? bbox_area / img_area : 0.0;
}

bool ParsePairEntry(const rapidjson::Value & entry, Pair & pair)
{
  IndexT i = UndefinedIndexT;
  IndexT j = UndefinedIndexT;

  if (entry.IsArray() && entry.Size() == 2)
  {
    if (!entry[0].IsUint64() || !entry[1].IsUint64())
    {
      return false;
    }
    i = static_cast<IndexT>(entry[0].GetUint64());
    j = static_cast<IndexT>(entry[1].GetUint64());
  }
  else if (entry.IsObject() && entry.HasMember("i") && entry.HasMember("j"))
  {
    if (!entry["i"].IsUint64() || !entry["j"].IsUint64())
    {
      return false;
    }
    i = static_cast<IndexT>(entry["i"].GetUint64());
    j = static_cast<IndexT>(entry["j"].GetUint64());
  }
  else
  {
    return false;
  }

  if (i == j || i == UndefinedIndexT || j == UndefinedIndexT)
  {
    return false;
  }
  pair = Pair(std::min(i, j), std::max(i, j));
  return true;
}

std::set<Pair> ParsePairArray(const rapidjson::Value & array)
{
  std::set<Pair> pairs;
  if (!array.IsArray())
  {
    return pairs;
  }
  for (rapidjson::SizeType k = 0; k < array.Size(); ++k)
  {
    Pair p;
    if (ParsePairEntry(array[k], p))
    {
      pairs.insert(p);
    }
  }
  return pairs;
}

ExclusionsV2 LoadExclusions(const std::string & path)
{
  ExclusionsV2 out;
  if (!stlplus::is_file(path))
  {
    return out;
  }

  std::ifstream ifs(path.c_str());
  if (!ifs.is_open())
  {
    OPENMVG_LOG_WARNING << "Could not open exclusions file: " << path;
    return out;
  }

  rapidjson::IStreamWrapper isw(ifs);
  rapidjson::Document doc;
  doc.ParseStream(isw);
  if (doc.HasParseError())
  {
    OPENMVG_LOG_WARNING << "Could not parse exclusions file: " << path;
    return out;
  }

  if (doc.IsObject() && doc.HasMember("version") && doc["version"].IsInt() && doc["version"].GetInt() == 2)
  {
    if (doc.HasMember("manual"))
    {
      out.manual = ParsePairArray(doc["manual"]);
    }
    if (doc.HasMember("auto"))
    {
      out.auto_pairs = ParsePairArray(doc["auto"]);
    }
    return out;
  }

  if (doc.IsObject() && doc.HasMember("excluded_pairs"))
  {
    out.manual = ParsePairArray(doc["excluded_pairs"]);
  }
  else if (doc.IsArray())
  {
    out.manual = ParsePairArray(doc);
  }
  return out;
}

bool SaveExclusionsV2(
  const std::string & path,
  const std::set<Pair> & manual_pairs,
  const std::set<Pair> & auto_pairs)
{
  std::ofstream os(path.c_str());
  if (!os.is_open())
  {
    OPENMVG_LOG_ERROR << "Could not write exclusions file: " << path;
    return false;
  }

  std::set<Pair> merged = manual_pairs;
  merged.insert(auto_pairs.begin(), auto_pairs.end());

  os << "{\n";
  os << "  \"version\": 2,\n";

  auto writePairArray = [&os](const char * name, const std::set<Pair> & pairs)
  {
    os << "  \"" << name << "\": [";
    bool first = true;
    for (const Pair & p : pairs)
    {
      if (!first) os << ",";
      first = false;
      os << "\n    [" << p.first << ", " << p.second << "]";
    }
    if (!pairs.empty()) os << "\n  ";
    os << "]";
  };

  writePairArray("manual", manual_pairs);
  os << ",\n";
  writePairArray("auto", auto_pairs);
  os << ",\n";
  writePairArray("excluded_pairs", merged);
  os << "\n}\n";
  return true;
}

bool LoadFeatureCoords(const std::string & feat_path, std::vector<Vec2> & out)
{
  out.clear();
  std::ifstream ifs(feat_path.c_str());
  if (!ifs.is_open())
  {
    return false;
  }

  std::string line;
  while (std::getline(ifs, line))
  {
    std::istringstream iss(line);
    double x = 0.0, y = 0.0;
    if (!(iss >> x >> y))
    {
      continue;
    }
    out.emplace_back(x, y);
  }
  return !out.empty();
}

bool EstimateHomographyDLT(
  const std::vector<Vec2> & src,
  const std::vector<Vec2> & dst,
  const std::vector<size_t> & indices,
  Mat3 & H)
{
  if (indices.size() < 4)
  {
    return false;
  }

  auto normalizePoints = [](const std::vector<Vec2> & pts, const std::vector<size_t> & idx, std::vector<Vec2> & out_pts, Mat3 & T)
  {
    Vec2 centroid(0.0, 0.0);
    for (size_t id : idx)
    {
      centroid += pts[id];
    }
    centroid /= static_cast<double>(idx.size());

    double mean_dist = 0.0;
    for (size_t id : idx)
    {
      mean_dist += (pts[id] - centroid).norm();
    }
    mean_dist /= static_cast<double>(idx.size());
    const double scale = (mean_dist > 1e-12) ? (std::sqrt(2.0) / mean_dist) : 1.0;

    T = Mat3::Identity();
    T(0, 0) = scale;
    T(1, 1) = scale;
    T(0, 2) = -scale * centroid(0);
    T(1, 2) = -scale * centroid(1);

    out_pts.clear();
    out_pts.reserve(idx.size());
    for (size_t id : idx)
    {
      const Vec3 hp = T * Vec3(pts[id](0), pts[id](1), 1.0);
      out_pts.emplace_back(hp(0) / hp(2), hp(1) / hp(2));
    }
  };

  std::vector<Vec2> src_n, dst_n;
  Mat3 T_src = Mat3::Identity(), T_dst = Mat3::Identity();
  normalizePoints(src, indices, src_n, T_src);
  normalizePoints(dst, indices, dst_n, T_dst);

  Mat A(static_cast<int>(indices.size()) * 2, 9);
  A.setZero();
  for (size_t k = 0; k < indices.size(); ++k)
  {
    const double x = src_n[k](0);
    const double y = src_n[k](1);
    const double u = dst_n[k](0);
    const double v = dst_n[k](1);

    A(2 * static_cast<int>(k), 3) = -x;
    A(2 * static_cast<int>(k), 4) = -y;
    A(2 * static_cast<int>(k), 5) = -1.0;
    A(2 * static_cast<int>(k), 6) = v * x;
    A(2 * static_cast<int>(k), 7) = v * y;
    A(2 * static_cast<int>(k), 8) = v;

    A(2 * static_cast<int>(k) + 1, 0) = x;
    A(2 * static_cast<int>(k) + 1, 1) = y;
    A(2 * static_cast<int>(k) + 1, 2) = 1.0;
    A(2 * static_cast<int>(k) + 1, 6) = -u * x;
    A(2 * static_cast<int>(k) + 1, 7) = -u * y;
    A(2 * static_cast<int>(k) + 1, 8) = -u;
  }

  const Eigen::JacobiSVD<Mat> svd(A, Eigen::ComputeFullV);
  if (svd.matrixV().cols() < 9)
  {
    return false;
  }
  const Vec h = svd.matrixV().col(8);

  Mat3 Hn;
  Hn << h(0), h(1), h(2),
        h(3), h(4), h(5),
        h(6), h(7), h(8);

  H = T_dst.inverse() * Hn * T_src;
  if (std::abs(H(2, 2)) > 1e-12)
  {
    H /= H(2, 2);
  }
  return H.allFinite();
}

double ReprojectionErrorPx(const Mat3 & H, const Vec2 & p, const Vec2 & q)
{
  const Vec3 hp = H * Vec3(p(0), p(1), 1.0);
  if (std::abs(hp(2)) < 1e-12)
  {
    return std::numeric_limits<double>::infinity();
  }
  const Vec2 proj(hp(0) / hp(2), hp(1) / hp(2));
  return (proj - q).norm();
}

bool CheckHomographyPlanar(
  const std::vector<Vec2> & pts_i,
  const std::vector<Vec2> & pts_j,
  const double max_error_px,
  bool & fits,
  double & pct_error)
{
  fits = false;
  pct_error = std::numeric_limits<double>::infinity();

  if (pts_i.size() < 4 || pts_i.size() != pts_j.size())
  {
    return false;
  }

  const double ransac_threshold = 4.0;
  const int iterations = 200;
  std::vector<size_t> all_idx(pts_i.size());
  std::iota(all_idx.begin(), all_idx.end(), 0);

  std::mt19937 rng(42);
  Mat3 best_H = Mat3::Identity();
  size_t best_inliers = 0;

  auto sampleFour = [&](std::array<size_t, 4> & sample)
  {
    for (size_t k = 0; k < 4; ++k)
    {
      size_t idx = 0;
      bool unique = false;
      while (!unique)
      {
        idx = static_cast<size_t>(rng() % pts_i.size());
        unique = true;
        for (size_t m = 0; m < k; ++m)
        {
          if (sample[m] == idx)
          {
            unique = false;
            break;
          }
        }
      }
      sample[k] = idx;
    }
  };

  for (int it = 0; it < iterations; ++it)
  {
    std::array<size_t, 4> s = {};
    sampleFour(s);
    std::vector<size_t> subset(s.begin(), s.end());

    Mat3 H;
    if (!EstimateHomographyDLT(pts_i, pts_j, subset, H))
    {
      continue;
    }

    size_t inliers = 0;
    for (size_t k = 0; k < pts_i.size(); ++k)
    {
      const double e = ReprojectionErrorPx(H, pts_i[k], pts_j[k]);
      if (e < ransac_threshold)
      {
        ++inliers;
      }
    }
    if (inliers > best_inliers)
    {
      best_inliers = inliers;
      best_H = H;
      if (best_inliers == pts_i.size())
      {
        break;  // perfect model found, no need for more iterations
      }
    }
  }

  if (best_inliers < 4)
  {
    return false;
  }

  std::vector<double> errors;
  errors.reserve(pts_i.size());
  for (size_t k = 0; k < pts_i.size(); ++k)
  {
    errors.push_back(ReprojectionErrorPx(best_H, pts_i[k], pts_j[k]));
  }
  if (errors.empty())
  {
    return false;
  }

  std::sort(errors.begin(), errors.end());
  const size_t idx = std::min(errors.size() - 1, static_cast<size_t>(std::floor(0.95 * static_cast<double>(errors.size() - 1))));
  pct_error = errors[idx];
  fits = pct_error < max_error_px;
  return true;
}

bool RecoverRelativeRotationFromPoints(
  const std::vector<Vec2> & pts_i,
  const std::vector<Vec2> & pts_j,
  const cameras::IntrinsicBase * cam_i,
  const cameras::IntrinsicBase * cam_j,
  Mat3 & R_ij)
{
  if (!cam_i || !cam_j || pts_i.size() < 5 || pts_i.size() != pts_j.size())
  {
    return false;
  }

  Mat2X xI(2, pts_i.size()), xJ(2, pts_j.size());
  for (size_t k = 0; k < pts_i.size(); ++k)
  {
    xI.col(k) = pts_i[k];
    xJ.col(k) = pts_j[k];
  }

  const Mat3X bI = (*cam_i)(xI);
  const Mat3X bJ = (*cam_j)(xJ);
  std::vector<Mat3> E_vec;
  openMVG::EightPointRelativePoseSolver::Solve(bI, bJ, &E_vec);
  if (E_vec.empty())
  {
    return false;
  }

  std::vector<geometry::Pose3> poses;
  openMVG::MotionFromEssential(E_vec[0], &poses);
  if (poses.empty())
  {
    return false;
  }

  const geometry::Pose3 pose_identity(Mat3::Identity(), Vec3::Zero());
  int best_idx = 0;
  int best_count = -1;
  const size_t n_test = std::min<size_t>(20, static_cast<size_t>(bI.cols()));
  for (size_t k = 0; k < poses.size(); ++k)
  {
    int count = 0;
    for (size_t m = 0; m < n_test; ++m)
    {
      Vec3 X;
      if (openMVG::Triangulate2View(
            pose_identity.rotation(), pose_identity.translation(), bI.col(m),
            poses[k].rotation(), poses[k].translation(), bJ.col(m),
            X, openMVG::ETriangulationMethod::DEFAULT))
      {
        ++count;
      }
    }
    if (count > best_count)
    {
      best_count = count;
      best_idx = static_cast<int>(k);
    }
  }

  R_ij = poses[static_cast<size_t>(best_idx)].rotation();
  return R_ij.allFinite();
}

bool LoadGeometricRotations(
  const std::string & path,
  std::unordered_map<Pair, Mat3, PairHash> & pair_rotations)
{
  pair_rotations.clear();
  if (!stlplus::is_file(path))
  {
    return false;
  }

  std::ifstream ifs(path.c_str());
  if (!ifs.is_open())
  {
    return false;
  }

  rapidjson::IStreamWrapper isw(ifs);
  rapidjson::Document doc;
  doc.ParseStream(isw);
  if (doc.HasParseError() || !doc.IsObject() || !doc.HasMember("pairs") || !doc["pairs"].IsArray())
  {
    return false;
  }

  for (const auto & entry : doc["pairs"].GetArray())
  {
    if (!entry.IsObject() ||
        !entry.HasMember("i") || !entry["i"].IsUint64() ||
        !entry.HasMember("j") || !entry["j"].IsUint64() ||
        !entry.HasMember("qw") || !entry["qw"].IsNumber() ||
        !entry.HasMember("qx") || !entry["qx"].IsNumber() ||
        !entry.HasMember("qy") || !entry["qy"].IsNumber() ||
        !entry.HasMember("qz") || !entry["qz"].IsNumber())
    {
      continue;
    }
    const IndexT i = static_cast<IndexT>(entry["i"].GetUint64());
    const IndexT j = static_cast<IndexT>(entry["j"].GetUint64());
    const Pair key(std::min(i, j), std::max(i, j));
    pair_rotations[key] = QuaternionToMat3(
      entry["qw"].GetDouble(),
      entry["qx"].GetDouble(),
      entry["qy"].GetDouble(),
      entry["qz"].GetDouble());
  }
  return true;
}

bool ParseArgs(int argc, char ** argv, Args & args)
{
  CmdLine cmd;
  cmd.add(make_option('i', args.input_dir, "input_dir"));
  cmd.add(make_option('t', args.min_time_delta, "min_time_delta"));
  cmd.add(make_option('e', args.max_homography_error, "max_homography_error"));
  cmd.add(make_option('r', args.min_imu_rotation_error, "min_imu_rotation_error"));
  cmd.add(make_option('b', args.max_bbox_area_ratio, "max_bbox_area_ratio"));
  cmd.add(make_switch('d', "dry_run"));

  try
  {
    if (argc == 1)
    {
      throw std::string("Invalid command line parameter.");
    }
    cmd.process(argc, argv);
    args.dry_run = cmd.used('d');
  }
  catch (const std::string &)
  {
    OPENMVG_LOG_INFO
      << "Usage: " << argv[0] << '\n'
      << "[-i|--input_dir]             Input dir (task dir, matches dir, or openmvg dir)\n"
      << "[-t|--min_time_delta]        Min timestamp delta in seconds (default: 30)\n"
      << "[-e|--max_homography_error]  Max 95th-percentile homography reprojection error in px (default: 5)\n"
      << "[-r|--min_imu_rotation_error] Min IMU-vs-geom rotation error in deg (default: 90)\n"
      << "[-b|--max_bbox_area_ratio]   Max bbox area ratio in each image (default: 0.10)\n"
      << "[-d|--dry_run]               Do not write excluded_pairs.json\n";
    return false;
  }

  if (args.input_dir.empty())
  {
    OPENMVG_LOG_ERROR << "input_dir is required";
    return false;
  }
  return true;
}

std::string ResolveMatchesDir(const std::string & input_dir)
{
  const std::string d1 = stlplus::folder_append_separator(input_dir) + "openmvg/matches";
  const std::string d2 = stlplus::folder_append_separator(input_dir) + "matches";
  if (stlplus::folder_exists(d1))
  {
    return d1;
  }
  if (stlplus::folder_exists(d2))
  {
    return d2;
  }
  return input_dir;
}

} // namespace

int main(int argc, char ** argv)
{
  Args args;
  if (!ParseArgs(argc, argv, args))
  {
    return EXIT_FAILURE;
  }

  const std::string matches_dir = ResolveMatchesDir(args.input_dir);
  const std::string sfm_data_path = stlplus::create_filespec(matches_dir, "sfm_data", "json");
  const std::string matches_path = stlplus::create_filespec(matches_dir, "matches.e", "bin");
  const std::string rotations_path = stlplus::create_filespec(matches_dir, "matches.e.rotations", "json");
  const std::string exclusions_path = stlplus::create_filespec(matches_dir, "excluded_pairs", "json");

  OPENMVG_LOG_INFO << "=== openMVG_main_GenAutoExclusions ===";
  OPENMVG_LOG_INFO << "Matches directory: " << matches_dir;
  OPENMVG_LOG_INFO << "Exclusions file: " << exclusions_path;
  OPENMVG_LOG_INFO << "Criteria: min_time_delta=" << args.min_time_delta
                   << "s, max_bbox_area_ratio=" << args.max_bbox_area_ratio
                   << ", max_homography_error=" << args.max_homography_error
                   << "px, min_imu_rotation_error=" << args.min_imu_rotation_error
                   << "deg, expansion_time_window=" << args.expansion_time_window << "s";

  SfM_Data sfm_data;
  if (!Load(sfm_data, sfm_data_path, ESfM_Data(VIEWS | INTRINSICS)))
  {
    OPENMVG_LOG_ERROR << "Cannot read sfm_data file: " << sfm_data_path;
    return EXIT_FAILURE;
  }

  PairWiseMatches all_pairs;
  if (!Load(all_pairs, matches_path))
  {
    OPENMVG_LOG_ERROR << "Cannot read geometric matches file: " << matches_path;
    return EXIT_FAILURE;
  }

  std::map<IndexT, ViewInfo> views;
  const Mat3 R_cd = DeviceToCameraRotation();
  const Mat3 R_landscape = DeviceLandscapeLeftRemap();
  for (const auto & kv : sfm_data.GetViews())
  {
    const View * v = kv.second.get();
    if (!v)
    {
      continue;
    }
    ViewInfo info;
    info.id_intrinsic = v->id_intrinsic;
    info.filename = v->s_Img_path;
    info.width = static_cast<int>(v->ui_width);
    info.height = static_cast<int>(v->ui_height);
    if (v->b_has_capture_time_epoch_)
    {
      info.has_timestamp = true;
      info.timestamp_s = v->capture_time_epoch_;
    }
    else
    {
      info.has_timestamp = ParseTimestampFromFilename(v->s_Img_path, info.timestamp_s);
    }

    const ViewPriors * vp = dynamic_cast<const ViewPriors *>(v);
    if (vp && vp->b_has_imu_rotation_)
    {
      info.has_imu_rotation = true;
      info.imu_rotation = R_cd * R_landscape * vp->imu_rotation_.transpose();
    }

    views[kv.first] = info;
  }

  std::unordered_map<Pair, Mat3, PairHash> geometric_rotations;
  if (LoadGeometricRotations(rotations_path, geometric_rotations))
  {
    OPENMVG_LOG_INFO << "Loaded " << geometric_rotations.size()
                     << " geometric pair rotations from: " << rotations_path;
  }
  else
  {
    OPENMVG_LOG_INFO << "Rotations file not found or invalid, will recover from inlier points when needed.";
  }

  ExclusionsV2 existing = LoadExclusions(exclusions_path);
  OPENMVG_LOG_INFO << "Existing exclusions: manual=" << existing.manual.size()
                   << ", auto=" << existing.auto_pairs.size();

  std::unordered_map<IndexT, std::vector<Vec2>> feature_cache;
  std::set<Pair> new_auto_pairs;
  std::vector<std::pair<IndexT, double>> timestamped_views;
  timestamped_views.reserve(views.size());
  for (const auto & kv : views)
  {
    if (kv.second.has_timestamp)
    {
      timestamped_views.emplace_back(kv.first, kv.second.timestamp_s);
    }
  }

  std::sort(timestamped_views.begin(), timestamped_views.end(),
    [](const std::pair<IndexT, double> & a, const std::pair<IndexT, double> & b)
    { return a.second < b.second; });

  size_t n_checked = 0;
  size_t n_temporal = 0;
  size_t n_bbox = 0;
  size_t n_hom = 0;
  size_t n_imu = 0;
  const size_t total_pairs = all_pairs.size();
  size_t next_progress_pct = 10;

  for (const auto & kv : all_pairs)
  {
    const Pair pair = kv.first;
    const IndMatches & matches = kv.second;
    ++n_checked;

    if (total_pairs > 0)
    {
      const size_t pct_done = (n_checked * 100) / total_pairs;
      while (pct_done >= next_progress_pct && next_progress_pct <= 100)
      {
        OPENMVG_LOG_INFO << "Analysis progress: " << next_progress_pct << "% (" << n_checked
                         << "/" << total_pairs << " pairs)";
        next_progress_pct += 10;
      }
    }

    const auto it_i = views.find(pair.first);
    const auto it_j = views.find(pair.second);
    if (it_i == views.end() || it_j == views.end())
    {
      continue;
    }
    const ViewInfo & vi = it_i->second;
    const ViewInfo & vj = it_j->second;

    PairStats stats;

    if (vi.has_timestamp && vj.has_timestamp)
    {
      stats.time_delta = std::abs(vi.timestamp_s - vj.timestamp_s);
      if (stats.time_delta > args.min_time_delta)
      {
        stats.temporal = true;
        ++n_temporal;
      }
    }
    if (!stats.temporal)
    {
      continue;  // skip feature loading, bbox, homography RANSAC, and IMU check
    }

    auto loadFeaturesFor = [&](const IndexT view_id, const ViewInfo & info) -> const std::vector<Vec2> *
    {
      const auto result = feature_cache.emplace(view_id, std::vector<Vec2>());
      if (result.second)
      {
        const std::string feat_path = stlplus::create_filespec(
          matches_dir, stlplus::basename_part(info.filename), "feat");
        LoadFeatureCoords(feat_path, result.first->second);
      }
      return &result.first->second;
    };

    const std::vector<Vec2> * feats_i = loadFeaturesFor(pair.first, vi);
    const std::vector<Vec2> * feats_j = loadFeaturesFor(pair.second, vj);

    std::vector<Vec2> pts_i;
    std::vector<Vec2> pts_j;
    if (feats_i && feats_j && !feats_i->empty() && !feats_j->empty())
    {
      pts_i.reserve(matches.size());
      pts_j.reserve(matches.size());
      for (const IndMatch & m : matches)
      {
        const size_t idx_i = static_cast<size_t>(m.i_);
        const size_t idx_j = static_cast<size_t>(m.j_);
        if (idx_i < feats_i->size() && idx_j < feats_j->size())
        {
          pts_i.push_back((*feats_i)[idx_i]);
          pts_j.push_back((*feats_j)[idx_j]);
        }
      }
      stats.has_pts = !pts_i.empty() && pts_i.size() == pts_j.size();
    }

    const bool bbox_disabled = args.max_bbox_area_ratio >= 1.0;
    bool bbox_small = bbox_disabled;
    if (stats.has_pts)
    {
      stats.bbox_i = BBoxAreaRatio(pts_i, vi.width, vi.height);
      stats.bbox_j = BBoxAreaRatio(pts_j, vj.width, vj.height);
      if (!bbox_disabled &&
          stats.bbox_i <= args.max_bbox_area_ratio &&
          stats.bbox_j <= args.max_bbox_area_ratio)
      {
        stats.bbox = true;
        bbox_small = true;
        ++n_bbox;
      }
    }
    if (!bbox_small)
    {
      continue;  // skip expensive homography RANSAC and IMU check
    }

    const bool hom_disabled = args.max_homography_error >= 9999.0;
    bool hom_fits = hom_disabled;
    if (!hom_disabled && stats.has_pts && pts_i.size() >= 4)
    {
      bool fits = false;
      double p95_err = std::numeric_limits<double>::infinity();
      if (CheckHomographyPlanar(pts_i, pts_j, args.max_homography_error, fits, p95_err))
      {
        stats.hom_err = p95_err;
        if (fits)
        {
          stats.hom = true;
          hom_fits = true;
          ++n_hom;
        }
      }
    }
    if (!hom_fits)
    {
      continue;  // skip IMU check
    }

    Mat3 R_imu_rel = Mat3::Identity();
    bool has_imu_rel = false;
    if (vi.has_imu_rotation && vj.has_imu_rotation)
    {
      R_imu_rel = vj.imu_rotation * vi.imu_rotation.transpose();
      has_imu_rel = true;
    }

    Mat3 R_geom = Mat3::Identity();
    bool has_geom = false;
    auto it_rot = geometric_rotations.find(pair);
    if (it_rot != geometric_rotations.end())
    {
      R_geom = it_rot->second;
      has_geom = true;
    }
    else if (stats.has_pts && pts_i.size() >= 5)
    {
      const auto it_intr_i = sfm_data.GetIntrinsics().find(vi.id_intrinsic);
      const auto it_intr_j = sfm_data.GetIntrinsics().find(vj.id_intrinsic);
      if (it_intr_i != sfm_data.GetIntrinsics().end() &&
          it_intr_j != sfm_data.GetIntrinsics().end())
      {
        has_geom = RecoverRelativeRotationFromPoints(
          pts_i, pts_j, it_intr_i->second.get(), it_intr_j->second.get(), R_geom);
      }
    }

    bool imu_mismatch = false;
    if (has_imu_rel && has_geom)
    {
      stats.imu_err = RotationErrorDeg(R_imu_rel, R_geom);
      if (stats.imu_err > args.min_imu_rotation_error)
      {
        stats.imu = true;
        imu_mismatch = true;
        ++n_imu;
      }
    }

    stats.exclude = imu_mismatch;  // temporal, bbox, hom already confirmed above
    if (stats.exclude)
    {
      new_auto_pairs.insert(pair);
      std::ostringstream oss;
      oss << "Auto exclude pair (" << pair.first << ", " << pair.second << ")"
          << " dt=" << std::fixed << std::setprecision(1)
          << (std::isnan(stats.time_delta) ? 0.0 : stats.time_delta) << "s"
          << " bbox=(" << std::setprecision(3) << stats.bbox_i << ", " << stats.bbox_j << ")"
          << " hom95=" << std::setprecision(2)
          << (std::isfinite(stats.hom_err) ? stats.hom_err : -1.0)
          << " imu_err=" << std::setprecision(1)
          << (std::isnan(stats.imu_err) ? 0.0 : stats.imu_err) << "deg"
          << " [" << vi.filename << "] [" << vj.filename << "]";
      OPENMVG_LOG_INFO << oss.str();
    }
  }

  // Expansion pass:
  // For each seed excluded pair (i, j), find views within +/- expansion_time_window
  // around timestamps of i and j, then add cross-side combinations only if that pair
  // exists in geometrically-filtered matches.
  const std::set<Pair> seed_pairs = new_auto_pairs;
  size_t expanded_added = 0;
  for (const Pair & seed : seed_pairs)
  {
    const auto it_i = views.find(seed.first);
    const auto it_j = views.find(seed.second);
    if (it_i == views.end() || it_j == views.end())
    {
      continue;
    }

    std::vector<IndexT> side_i(1, seed.first);
    std::vector<IndexT> side_j(1, seed.second);

    auto collectInWindow = [&](const double t_center, std::vector<IndexT> & side)
    {
      const double lo = t_center - args.expansion_time_window;
      const double hi = t_center + args.expansion_time_window;
      const auto begin_it = std::lower_bound(
        timestamped_views.begin(), timestamped_views.end(), lo,
        [](const std::pair<IndexT, double> & a, double val) { return a.second < val; });
      const auto end_it = std::upper_bound(
        timestamped_views.begin(), timestamped_views.end(), hi,
        [](double val, const std::pair<IndexT, double> & a) { return val < a.second; });
      for (auto it = begin_it; it != end_it; ++it)
      {
        side.push_back(it->first);
      }
    };

    if (it_i->second.has_timestamp)
    {
      collectInWindow(it_i->second.timestamp_s, side_i);
    }

    if (it_j->second.has_timestamp)
    {
      collectInWindow(it_j->second.timestamp_s, side_j);
    }

    std::sort(side_i.begin(), side_i.end());
    side_i.erase(std::unique(side_i.begin(), side_i.end()), side_i.end());
    std::sort(side_j.begin(), side_j.end());
    side_j.erase(std::unique(side_j.begin(), side_j.end()), side_j.end());

    for (const IndexT a : side_i)
    {
      for (const IndexT b : side_j)
      {
        if (a == b)
        {
          continue;
        }
        const Pair p(std::min(a, b), std::max(a, b));
        if (all_pairs.find(p) == all_pairs.end())
        {
          continue;
        }
        const auto inserted = new_auto_pairs.insert(p);
        if (inserted.second)
        {
          ++expanded_added;
        }
      }
    }
  }

  OPENMVG_LOG_INFO << "Summary:";
  OPENMVG_LOG_INFO << "  Pairs analysed                   : " << n_checked;
  OPENMVG_LOG_INFO << "  Flagged by temporal check        : " << n_temporal;
  OPENMVG_LOG_INFO << "  Flagged by bounding-box check    : " << n_bbox;
  OPENMVG_LOG_INFO << "  Flagged by homography check      : " << n_hom;
  OPENMVG_LOG_INFO << "  Flagged by IMU mismatch check    : " << n_imu;
  OPENMVG_LOG_INFO << "  Pairs meeting all criteria (new) : " << new_auto_pairs.size();
  OPENMVG_LOG_INFO << "  Added by +/-" << args.expansion_time_window << "s expansion : " << expanded_added;
  OPENMVG_LOG_INFO << "  Existing manual exclusions       : " << existing.manual.size();

  if (args.dry_run)
  {
    OPENMVG_LOG_INFO << "dry_run enabled; not writing exclusions file.";
    return EXIT_SUCCESS;
  }

  if (!SaveExclusionsV2(exclusions_path, existing.manual, new_auto_pairs))
  {
    return EXIT_FAILURE;
  }

  std::set<Pair> added = new_auto_pairs;
  for (const Pair & p : existing.auto_pairs)
  {
    added.erase(p);
  }

  std::set<Pair> removed = existing.auto_pairs;
  for (const Pair & p : new_auto_pairs)
  {
    removed.erase(p);
  }

  OPENMVG_LOG_INFO << "Auto exclusions written: " << new_auto_pairs.size();
  OPENMVG_LOG_INFO << "Newly added: " << added.size();
  OPENMVG_LOG_INFO << "Removed: " << removed.size();
  OPENMVG_LOG_INFO << "Wrote: " << exclusions_path;

  return EXIT_SUCCESS;
}
