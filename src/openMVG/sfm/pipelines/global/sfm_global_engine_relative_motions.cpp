// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/types.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <ceres/types.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace {

bool ParseImuRotationFromUserComment(const std::string & comment, openMVG::Mat3 & rotation_dw)
{
  const std::string key = "Rotation:";
  const std::size_t pos = comment.find(key);
  if (pos == std::string::npos)
    return false;

  std::string rot_part = comment.substr(pos + key.size());
  for (char & c : rot_part)
  {
    if (c == '\n' || c == '\r')
      c = ' ';
  }

  std::vector<double> vals;
  vals.reserve(9);
  std::string token;
  std::stringstream ss(rot_part);
  while (std::getline(ss, token, ','))
  {
    std::stringstream t(token);
    double v = 0.0;
    if (t >> v)
      vals.push_back(v);
  }
  if (vals.size() < 9)
    return false;

  rotation_dw << vals[0], vals[1], vals[2],
                  vals[3], vals[4], vals[5],
                  vals[6], vals[7], vals[8];
  return true;
}

openMVG::Mat3 DeviceToCameraRotation()
{
  // Rear camera, device axes to camera axes (X right, Y down, Z forward).
  openMVG::Mat3 R;
  R << 1.0, 0.0, 0.0,
       0.0, -1.0, 0.0,
       0.0, 0.0, -1.0;
  return R;
}

openMVG::Mat3 DeviceLandscapeLeftRemap()
{
  // Landscape-left: phone rotated 90° CCW from portrait (top of phone on left).
  // This matrix transforms from portrait device axes to landscape device axes.
  // When multiplied as R_dc_portrait * R_landscape, gives R_dc_landscape.
  // Portrait -> Landscape: X -> -Y, Y -> X, Z -> Z (90° CW about Z).
  openMVG::Mat3 R;
  R << 0.0, -1.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 0.0, 1.0;
  return R;
}

double RotationAngularErrorDeg(const openMVG::Mat3 & R_est, const openMVG::Mat3 & R_prior)
{
  const openMVG::Mat3 R_err = R_est * R_prior.transpose();
  return openMVG::R2D(openMVG::getRotationMagnitude(R_err));
}

openMVG::Mat3 BestFitRotation(const std::vector<std::pair<openMVG::Mat3, openMVG::Mat3>> & pairs)
{
  if (pairs.empty())
    return openMVG::Mat3::Identity();

  openMVG::Mat3 M = openMVG::Mat3::Zero();
  for (const auto & pair : pairs)
  {
    const openMVG::Mat3 & R_global = pair.first;
    const openMVG::Mat3 & R_imu = pair.second;
    M += R_imu * R_global.transpose();
  }

  Eigen::JacobiSVD<openMVG::Mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  openMVG::Mat3 R = svd.matrixU() * svd.matrixV().transpose();
  if (R.determinant() < 0.0)
  {
    openMVG::Mat3 U = svd.matrixU();
    U.col(2) *= -1.0;
    R = U * svd.matrixV().transpose();
  }
  return R;
}

openMVG::Hash_Map<openMVG::IndexT, openMVG::IndexT>
BuildPoseToViewMap(const openMVG::sfm::SfM_Data & sfm_data)
{
  openMVG::Hash_Map<openMVG::IndexT, openMVG::IndexT> map_pose_to_view;
  for (const auto & view_it : sfm_data.GetViews())
  {
    const openMVG::sfm::View * view = view_it.second.get();
    if (!view || view->id_pose == openMVG::UndefinedIndexT)
      continue;
    if (map_pose_to_view.count(view->id_pose) == 0)
      map_pose_to_view[view->id_pose] = view->id_view;
  }
  return map_pose_to_view;
}

std::string ViewImagePath(const openMVG::sfm::SfM_Data & sfm_data, openMVG::IndexT view_id)
{
  const auto it = sfm_data.GetViews().find(view_id);
  if (it == sfm_data.GetViews().end() || !it->second)
    return std::string();
  return it->second->s_Img_path;
}

void WriteRotationsCsv(const std::string & path,
                       const openMVG::sfm::SfM_Data & sfm_data,
                       const openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3> & rotations)
{
  OPENMVG_LOG_INFO << "Writing CSV: " << path;
  std::ofstream csv(path.c_str());
  if (!csv)
    return;

  const auto pose_to_view = BuildPoseToViewMap(sfm_data);

  csv << "pose_id,view_id,image_path,qw,qx,qy,qz\n";
  csv << std::fixed << std::setprecision(10);
  for (const auto & it : rotations)
  {
    const openMVG::IndexT pose_id = it.first;
    const openMVG::Mat3 & R = it.second;
    openMVG::Quaternion q(R);
    q.normalize();

    openMVG::IndexT view_id = openMVG::UndefinedIndexT;
    std::string image_path;
    const auto view_it = pose_to_view.find(pose_id);
    if (view_it != pose_to_view.end())
    {
      view_id = view_it->second;
      image_path = ViewImagePath(sfm_data, view_id);
    }

    csv << pose_id << "," << view_id << "," << image_path << ","
        << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << "\n";
  }
}

void WritePosesCsv(const std::string & path, const openMVG::sfm::SfM_Data & sfm_data)
{
  OPENMVG_LOG_INFO << "Writing CSV: " << path;
  std::ofstream csv(path.c_str());
  if (!csv)
    return;

  const auto pose_to_view = BuildPoseToViewMap(sfm_data);

  csv << "pose_id,view_id,image_path,cx,cy,cz,qw,qx,qy,qz\n";
  csv << std::fixed << std::setprecision(10);
  for (const auto & it : sfm_data.GetPoses())
  {
    const openMVG::IndexT pose_id = it.first;
    const openMVG::geometry::Pose3 & pose = it.second;
    const openMVG::Vec3 & c = pose.center();
    openMVG::Quaternion q(pose.rotation());
    q.normalize();

    openMVG::IndexT view_id = openMVG::UndefinedIndexT;
    std::string image_path;
    const auto view_it = pose_to_view.find(pose_id);
    if (view_it != pose_to_view.end())
    {
      view_id = view_it->second;
      image_path = ViewImagePath(sfm_data, view_id);
    }

    csv << pose_id << "," << view_id << "," << image_path << ","
        << c(0) << "," << c(1) << "," << c(2) << ","
        << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << "\n";
  }
}

void WriteRelativeRotationsCsv(const std::string & path,
                               const openMVG::sfm::SfM_Data & sfm_data,
                               const openMVG::rotation_averaging::RelativeRotations & relatives_R)
{
  OPENMVG_LOG_INFO << "Writing CSV: " << path;
  std::ofstream csv(path.c_str());
  if (!csv)
    return;

  const auto pose_to_view = BuildPoseToViewMap(sfm_data);

  csv << "pose_i,pose_j,view_i,view_j,image_i,image_j,qw,qx,qy,qz,weight\n";
  csv << std::fixed << std::setprecision(10);
  for (const auto & rel : relatives_R)
  {
    openMVG::Quaternion q(rel.Rij);
    q.normalize();

    openMVG::IndexT view_i = openMVG::UndefinedIndexT;
    openMVG::IndexT view_j = openMVG::UndefinedIndexT;
    std::string image_i;
    std::string image_j;
    const auto it_i = pose_to_view.find(rel.i);
    const auto it_j = pose_to_view.find(rel.j);
    if (it_i != pose_to_view.end())
    {
      view_i = it_i->second;
      image_i = ViewImagePath(sfm_data, view_i);
    }
    if (it_j != pose_to_view.end())
    {
      view_j = it_j->second;
      image_j = ViewImagePath(sfm_data, view_j);
    }

    csv << rel.i << "," << rel.j << ","
        << view_i << "," << view_j << ","
        << image_i << "," << image_j << ","
        << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << ","
        << rel.weight << "\n";
  }
}

} // namespace

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

GlobalSfMReconstructionEngine_RelativeMotions::GlobalSfMReconstructionEngine_RelativeMotions(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), sLogging_file_(sloggingFile)
{

  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("GlobalSfMReconstructionEngine_RelativeMotions")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  eRotation_averaging_method_ = ROTATION_AVERAGING_L2;
  eTranslation_averaging_method_ = TRANSLATION_AVERAGING_L1;
}

GlobalSfMReconstructionEngine_RelativeMotions::~GlobalSfMReconstructionEngine_RelativeMotions()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_);
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  eRotation_averaging_method_ = eRotationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  eTranslation_averaging_method_ = eTranslationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetImuRotationPrior(
  double weight,
  bool filter_outliers,
  double max_error_deg,
  double histogram_bucket_deg)
{
  imu_rotation_weight_ = weight;
  imu_rotation_filter_outliers_ = filter_outliers;
  imu_rotation_max_error_deg_ = max_error_deg;
  imu_rotation_histogram_bucket_deg_ = histogram_bucket_deg;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (set_remainingIds.empty())
    {
      OPENMVG_LOG_WARNING << "Invalid input image graph for global SfM";
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
  }

  OPENMVG_LOG_INFO << "Starting relative motion estimation (solving for R and t per view pair)...";
  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Rotation Averaging failure!";
    return false;
  }

  matching::PairWiseMatches tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Translation Averaging failure!";
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Cannot initialize an initial structure!";
    return false;
  }
  if (!Adjust())
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Non-linear adjustment failure!";
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

/// Compute from relative rotations the global rotations of the camera poses
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations
)
{
  if (relatives_R.empty())
    return false;
  if (!sOut_directory_.empty())
  {
    const std::string csv_path =
      stlplus::create_filespec(sOut_directory_, "relative_rotations_before_rotation_averaging", "csv");
    WriteRelativeRotationsCsv(csv_path, sfm_data_, relatives_R);
  }
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    OPENMVG_LOG_INFO << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size();
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    //TRIPLET_ROTATION_INFERENCE_NONE;

  rotation_averaging::RelativeRotations relatives_with_priors = relatives_R;
  IndexT world_pose_id = UndefinedIndexT;
  size_t imu_prior_count = 0;
  size_t imu_filtered_count = 0;
  std::set<IndexT> imu_filtered_pose_ids;
  if (imu_rotation_weight_ > 0.0 && !imu_pose_rotations_.empty())
  {
    IndexT max_pose_id = 0;
    bool has_pose_id = false;
    for (const auto & view : sfm_data_.GetViews())
    {
      const IndexT pose_id = view.second->id_pose;
      if (pose_id != UndefinedIndexT)
      {
        has_pose_id = true;
        max_pose_id = std::max(max_pose_id, pose_id);
      }
    }
    if (has_pose_id)
    {
      world_pose_id = max_pose_id + 1;
      for (const auto & imu_it : imu_pose_rotations_)
      {
        relatives_with_priors.emplace_back(
          world_pose_id, imu_it.first, imu_it.second,
          static_cast<float>(imu_rotation_weight_));
        ++imu_prior_count;
      }
    }
    else
    {
      OPENMVG_LOG_WARNING << "IMU rotation priors skipped: no valid pose IDs.";
    }
  }

  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  auto run_rotation_averaging = [&](const rotation_averaging::RelativeRotations & rels) -> bool
  {
    system::Timer t;
    const bool ok = rotation_averaging_solver.Run(
      eRotation_averaging_method_, eRelativeRotationInferenceMethod,
      rels, global_rotations, world_pose_id);

    OPENMVG_LOG_INFO
      << "Found #global_rotations: " << global_rotations.size() << "\n"
      << "Timing: " << t.elapsed() << " seconds";
    return ok;
  };

  bool b_rotation_averaging = run_rotation_averaging(relatives_with_priors);
  if (b_rotation_averaging && imu_rotation_filter_outliers_ &&
      imu_rotation_weight_ > 0.0 && !imu_pose_rotations_.empty() &&
      imu_rotation_max_error_deg_ > 0.0)
  {
    std::vector<IndexT> outlier_pose_ids;
    for (const auto & imu_it : imu_pose_rotations_)
    {
      const auto it = global_rotations.find(imu_it.first);
      if (it == global_rotations.end())
        continue;
      const double err_deg = RotationAngularErrorDeg(it->second, imu_it.second);
      if (err_deg > imu_rotation_max_error_deg_)
        outlier_pose_ids.push_back(imu_it.first);
    }

    if (!outlier_pose_ids.empty())
    {
      imu_filtered_pose_ids.insert(outlier_pose_ids.begin(), outlier_pose_ids.end());
      imu_filtered_count = outlier_pose_ids.size();
      OPENMVG_LOG_WARNING << "IMU rotation prior filtering: removing "
                          << outlier_pose_ids.size() << " / "
                          << imu_pose_rotations_.size()
                          << " priors with error > "
                          << imu_rotation_max_error_deg_ << " deg.";

      rotation_averaging::RelativeRotations filtered_rels = relatives_R;
      for (const auto & imu_it : imu_pose_rotations_)
      {
        if (std::find(outlier_pose_ids.begin(), outlier_pose_ids.end(), imu_it.first) != outlier_pose_ids.end())
          continue;
        filtered_rels.emplace_back(
          world_pose_id, imu_it.first, imu_it.second,
          static_cast<float>(imu_rotation_weight_));
      }
      b_rotation_averaging = run_rotation_averaging(filtered_rels);
      imu_prior_count = imu_pose_rotations_.size() - outlier_pose_ids.size();
    }
    else
    {
      OPENMVG_LOG_INFO << "IMU rotation prior filtering: removed 0 / "
                       << imu_pose_rotations_.size()
                       << " priors with error > "
                       << imu_rotation_max_error_deg_ << " deg.";
    }
  }

  if (b_rotation_averaging && !imu_pose_rotations_.empty())
  {
    OPENMVG_LOG_INFO
      << "IMU rotation prior settings: weight=" << imu_rotation_weight_
      << ", max_error_deg=" << imu_rotation_max_error_deg_
      << ", filter_outliers=" << (imu_rotation_filter_outliers_ ? "true" : "false")
      << ", filtered=" << imu_filtered_count;

    std::vector<double> imu_errors_deg;
    imu_errors_deg.reserve(imu_pose_rotations_.size());
    size_t missing_global = 0;
    std::vector<std::pair<Mat3, Mat3>> imu_pairs;
    imu_pairs.reserve(imu_pose_rotations_.size());
    for (const auto & imu_it : imu_pose_rotations_)
    {
      if (imu_filtered_pose_ids.count(imu_it.first) > 0)
        continue;
      const auto it = global_rotations.find(imu_it.first);
      if (it == global_rotations.end())
      {
        ++missing_global;
        continue;
      }
      const double err_deg = RotationAngularErrorDeg(it->second, imu_it.second);
      imu_errors_deg.push_back(err_deg);
      imu_pairs.emplace_back(it->second, imu_it.second);
    }

    if (!imu_errors_deg.empty())
    {
      const double bucket = std::max(1.0, imu_rotation_histogram_bucket_deg_);
      const double max_err = *std::max_element(imu_errors_deg.cbegin(), imu_errors_deg.cend());
      const size_t bucket_count = static_cast<size_t>(std::ceil(max_err / bucket)) + 1;
      std::vector<size_t> buckets(bucket_count, 0);
      for (const double err : imu_errors_deg)
      {
        const size_t idx = std::min(static_cast<size_t>(err / bucket), bucket_count - 1);
        buckets[idx]++;
      }

      std::ostringstream os;
      os << "\nIMU rotation prior errors (deg):\n";
      minMaxMeanMedian<double>(imu_errors_deg.cbegin(), imu_errors_deg.cend(), os);
      os << "Histogram (" << bucket << " deg buckets):\n";
      for (size_t i = 0; i < buckets.size(); ++i)
      {
        const double start = i * bucket;
        const double end = (i + 1) * bucket;
        os << "  [" << start << "," << end << "): " << buckets[i] << "\n";
      }
      os << "IMU priors used: " << imu_prior_count
         << " (missing global rotations: " << missing_global << ")\n";
      OPENMVG_LOG_INFO << os.str();
    }

    if (imu_rotation_weight_ <= 0.0 && !imu_pairs.empty())
    {
      const Mat3 R_align = BestFitRotation(imu_pairs);
      std::vector<double> aligned_errors_deg;
      aligned_errors_deg.reserve(imu_pairs.size());
      for (const auto & pair : imu_pairs)
      {
        const Mat3 R_aligned = R_align * pair.first;
        aligned_errors_deg.push_back(RotationAngularErrorDeg(R_aligned, pair.second));
      }

      if (!aligned_errors_deg.empty())
      {
        const double bucket = std::max(1.0, imu_rotation_histogram_bucket_deg_);
        const double max_err = *std::max_element(aligned_errors_deg.cbegin(), aligned_errors_deg.cend());
        const size_t bucket_count = static_cast<size_t>(std::ceil(max_err / bucket)) + 1;
        std::vector<size_t> buckets(bucket_count, 0);
        for (const double err : aligned_errors_deg)
        {
          const size_t idx = std::min(static_cast<size_t>(err / bucket), bucket_count - 1);
          buckets[idx]++;
        }

        std::ostringstream os;
        os << "\nIMU rotation prior errors after best-fit alignment (deg):\n";
        minMaxMeanMedian<double>(aligned_errors_deg.cbegin(), aligned_errors_deg.cend(), os);
        os << "Histogram (" << bucket << " deg buckets):\n";
        for (size_t i = 0; i < buckets.size(); ++i)
        {
          const double start = i * bucket;
          const double end = (i + 1) * bucket;
          os << "  [" << start << "," << end << "): " << buckets[i] << "\n";
        }
        OPENMVG_LOG_INFO << os.str();
      }
    }
  }


  if (b_rotation_averaging)
  {
    // Compute & display rotation fitting residual errors
    std::vector<float> vec_rotation_fitting_error;
    vec_rotation_fitting_error.reserve(relatives_R.size());
    for (const auto & relative_R : relatives_R)
    {
      const Mat3 & Rij = relative_R.Rij;
      const IndexT i = relative_R.i;
      const IndexT j = relative_R.j;
      if (global_rotations.count(i)==0 || global_rotations.count(j)==0)
        continue;
      const Mat3 & Ri = global_rotations[i];
      const Mat3 & Rj = global_rotations[j];
      const Mat3 eRij(Rj.transpose()*Rij*Ri);
      const double angularErrorDegree = R2D(getRotationMagnitude(eRij));
      vec_rotation_fitting_error.push_back(angularErrorDegree);
    }

    if (!vec_rotation_fitting_error.empty())
    {
      const float error_max = *max_element(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend());
      Histogram<float> histo(0.0f,error_max, 20);
      histo.Add(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend());
      OPENMVG_LOG_INFO
        << "\nRelative/Global degree rotations residual errors {0," << error_max<< "}:\n"
        << histo.ToString();
      {
        Histogram<float> histo(0.0f, 5.0f, 20);
        histo.Add(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend());
        OPENMVG_LOG_INFO
          << "\nRelative/Global degree rotations residual errors {0,5}:\n"
          << histo.ToString();
      }
      std::ostringstream os;
      os << "\nStatistics about global rotation evaluation:\n";
      minMaxMeanMedian<float>(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend(), os);
      OPENMVG_LOG_INFO << os.str();
    }

    // Log input graph to the HTML report
    if (!sLogging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : sfm_data_.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
        graph::exportToGraphvizData(
          stlplus::create_filespec(sOut_directory_, sGraph_name),
          putativeGraph);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(sOut_directory_, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        html_doc_stream_->pushInfo(os.str());
      }
    }
  }
  if (b_rotation_averaging && !sOut_directory_.empty())
  {
    const std::string csv_path =
      stlplus::create_filespec(sOut_directory_, "camera_rotations_after_rotation_averaging", "csv");
    WriteRotationsCsv(csv_path, sfm_data_, global_rotations);
  }
  return b_rotation_averaging;
}

/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches,
    sOut_directory_);

  if (!sOut_directory_.empty())
  {
    const std::string csv_path =
      stlplus::create_filespec(sOut_directory_, "camera_poses_after_translation_averaging", "csv");
    WritePosesCsv(csv_path, sfm_data_);
  }

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

/// Compute the initial structure of the scene
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Initial_Structure
(
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
#if defined USE_ALL_VALID_MATCHES // not used by default
    matching::PairWiseMatches pose_supported_matches;
    for (const std::pair<Pair, IndMatches> & match_info :  matches_provider_->pairWise_matches_)
    {
      const View * vI = sfm_data_.GetViews().at(match_info.first.first).get();
      const View * vJ = sfm_data_.GetViews().at(match_info.first.second).get();
      if (sfm_data_.IsPoseAndIntrinsicDefined(vI) && sfm_data_.IsPoseAndIntrinsicDefined(vJ))
      {
        pose_supported_matches.insert(match_info);
      }
    }
    tracksBuilder.Build(pose_supported_matches);
#else
    // Use triplet validated matches
    tracksBuilder.Build(tripletWise_matches);
#endif
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data_.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = features_provider_->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    {
      std::ostringstream osTrack;
      //-- Display stats (suppress image id listing)
      osTrack
        << "\n------------------\n"
        << "-- Tracks Stats --\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n";

      std::map<uint32_t, uint32_t> map_Occurrence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurrence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (const auto & iter : map_Occurrence_TrackLength)  {
        osTrack << "\t" << iter.first << "\t" << iter.second << "\n";
      }
      osTrack << "\n";
      OPENMVG_LOG_INFO << osTrack.str();
    }
  }

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(sfm_data_);

    OPENMVG_LOG_INFO << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size());
    OPENMVG_LOG_INFO << "Triangulation took (s): " << timer.elapsed();

    // Export initial structure
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !sfm_data_.structure.empty();
}

// Adjust the scene (& remove outliers)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their uncertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  size_t view_priors_count = 0;
  size_t pose_center_prior_count = 0;
  for (const auto & view_it : sfm_data_.GetViews())
  {
    const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
    if (prior != nullptr)
    {
      ++view_priors_count;
      if (prior->b_use_pose_center_)
      {
        ++pose_center_prior_count;
      }
    }
  }
  // - refine only Structure and translations
  OPENMVG_LOG_INFO << "Bundle adjustment: refine translations + structure..."
                   << " (motion priors " << (this->b_use_motion_prior_ ? "enabled" : "disabled")
                   << ", view_priors=" << view_priors_count
                   << ", pose_center_prior=" << pose_center_prior_count << ")";
  bundle_adjustment_obj.ceres_options().progress_modulo_ = 5;
  bundle_adjustment_obj.ceres_options().progress_label_ = "T + X";
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  if (b_BA_Status)
  {
    if (!sOut_directory_.empty())
    {
      const std::string csv_path =
        stlplus::create_filespec(sOut_directory_, "camera_poses_after_ba_T_X", "csv");
      WritePosesCsv(csv_path, sfm_data_);
    }
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    OPENMVG_LOG_INFO << "Bundle adjustment: refine rotations + translations + structure..."
             << " (motion priors " << (this->b_use_motion_prior_ ? "enabled" : "disabled")
             << ", view_priors=" << view_priors_count
             << ", pose_center_prior=" << pose_center_prior_count << ")";
    bundle_adjustment_obj.ceres_options().progress_modulo_ = 5;
    bundle_adjustment_obj.ceres_options().progress_label_ = "R + T + X";
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    if (b_BA_Status && !sOut_directory_.empty())
    {
      const std::string csv_path =
        stlplus::create_filespec(sOut_directory_, "camera_poses_after_ba_RT_X", "csv");
      WritePosesCsv(csv_path, sfm_data_);
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    OPENMVG_LOG_INFO << "Bundle adjustment: refine intrinsics + motion + structure..."
             << " (motion priors " << (this->b_use_motion_prior_ ? "enabled" : "disabled")
             << ", view_priors=" << view_priors_count
             << ", pose_center_prior=" << pose_center_prior_count << ")";
    bundle_adjustment_obj.ceres_options().progress_modulo_ = 2;
    bundle_adjustment_obj.ceres_options().progress_label_ = "K + R + T + X";
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
    if (b_BA_Status && !sOut_directory_.empty())
    {
      const std::string csv_path =
        stlplus::create_filespec(sOut_directory_, "camera_poses_after_ba_KRT_X", "csv");
      WritePosesCsv(csv_path, sfm_data_);
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  OPENMVG_LOG_INFO << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    OPENMVG_LOG_INFO << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  OPENMVG_LOG_INFO << "Bundle adjustment: final refine after outlier removal..."
                   << " (motion priors " << (this->b_use_motion_prior_ ? "enabled" : "disabled")
                   << ", view_priors=" << view_priors_count
                   << ", pose_center_prior=" << pose_center_prior_count << ")";
  bundle_adjustment_obj.ceres_options().progress_modulo_ = 2;
  bundle_adjustment_obj.ceres_options().progress_label_ = "final";
  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }
  if (b_BA_Status && !sOut_directory_.empty())
  {
    const std::string csv_path =
      stlplus::create_filespec(sOut_directory_, "camera_poses_after_ba_final", "csv");
    WritePosesCsv(csv_path, sfm_data_);
  }
  return b_BA_Status;
}

void GlobalSfMReconstructionEngine_RelativeMotions::Compute_Relative_Rotations
(
  rotation_averaging::RelativeRotations & vec_relatives_R
)
{
  // Compute a relative pose for each edge of the pose pair graph
  const Relative_Pose_Engine::Relative_Pair_Poses relative_poses = [&]
  {
    Relative_Pose_Engine relative_pose_engine;
    if (!relative_pose_engine.Process(sfm_data_,
        matches_provider_,
        features_provider_))
      return Relative_Pose_Engine::Relative_Pair_Poses();
    else
      return relative_pose_engine.Get_Relative_Poses();
  }();

  // Export the rotation component from the computed relative poses
  OPENMVG_LOG_INFO << "Extracting relative rotations...";
  system::LoggerProgress extract_progress(
    static_cast<std::uint32_t>(relative_poses.size()),
    "- Relative rotation extraction -",
    25);
  for (const auto & relative_pose : relative_poses)
  {
    // Add the relative rotation to the relative 'rotation' pose graph
    vec_relatives_R.emplace_back(
      relative_pose.first.first, relative_pose.first.second,
      relative_pose.second.rotation(),
      1.f);
    ++extract_progress;
  }
  OPENMVG_LOG_INFO << "Relative rotations extracted.";

  if ((imu_rotation_weight_ > 0.0 || imu_rotation_histogram_bucket_deg_ > 0.0) && !imu_rotations_loaded_)
  {
    imu_pose_rotations_.clear();
    const Mat3 R_cd = DeviceToCameraRotation();
    const Mat3 R_landscape = DeviceLandscapeLeftRemap();

    // Debug: print transformation matrices once
    {
      std::ostringstream os;
      os << "\n=== IMU Debug: Transformation Matrices ===\n";
      os << "R_cd (DeviceToCameraRotation - portrait):\n" << R_cd << "\n\n";
      os << "R_landscape (DeviceLandscapeLeftRemap):\n" << R_landscape << "\n\n";
      os << "R_cd * R_landscape (combined device-to-camera for landscape-left):\n" << (R_cd * R_landscape) << "\n";
      OPENMVG_LOG_INFO << os.str();
    }

    size_t debug_count = 0;
    for (const auto & view_it : sfm_data_.GetViews())
    {
      const IndexT pose_id = view_it.second->id_pose;
      if (pose_id == UndefinedIndexT)
        continue;

      const std::string image_path =
        stlplus::folder_append_separator(sfm_data_.s_root_path) + view_it.second->s_Img_path;
      std::unique_ptr<openMVG::exif::Exif_IO> exifIO(
        new openMVG::exif::Exif_IO_EasyExif(image_path));

      std::string user_comment;
      if (!exifIO->UserComment(&user_comment))
        continue;

      Mat3 R_dw;
      if (!ParseImuRotationFromUserComment(user_comment, R_dw))
        continue;

      // R_dw: Device -> World (ENU). For landscape-left, we need:
      // R_dc_landscape = R_dc_portrait * R_landscape
      // R_wc = R_dc_landscape * R_dw^T = R_cd * R_landscape * R_dw^T
      // (NOT: R_cd * (R_dw * R_landscape)^T, which incorrectly transposes R_landscape)
      const Mat3 R_wc = R_cd * R_landscape * R_dw.transpose();
      imu_pose_rotations_[pose_id] = R_wc;

      // Debug: print first 3 IMU rotations
      if (debug_count < 0)
      {
        std::ostringstream os;
        os << "\n=== IMU Debug: View " << view_it.first << " (pose " << pose_id << ") ===\n";
        os << "Image: " << view_it.second->s_Img_path << "\n";
        os << "R_dw (Device->World from EXIF, row-major):\n" << R_dw << "\n\n";
        os << "R_dw columns (device axes in world coords):\n";
        os << "  Dev X in world: (" << R_dw(0,0) << ", " << R_dw(1,0) << ", " << R_dw(2,0) << ")\n";
        os << "  Dev Y in world: (" << R_dw(0,1) << ", " << R_dw(1,1) << ", " << R_dw(2,1) << ")\n";
        os << "  Dev Z in world: (" << R_dw(0,2) << ", " << R_dw(1,2) << ", " << R_dw(2,2) << ")\n";
        os << "  Camera looks at -Dev Z: (" << -R_dw(0,2) << ", " << -R_dw(1,2) << ", " << -R_dw(2,2) << ")\n\n";
        os << "R_wc (World->Camera, computed):\n" << R_wc << "\n\n";
        os << "R_wc rows (camera axes in world coords):\n";
        os << "  Cam X (right) in world: (" << R_wc(0,0) << ", " << R_wc(0,1) << ", " << R_wc(0,2) << ")\n";
        os << "  Cam Y (down) in world: (" << R_wc(1,0) << ", " << R_wc(1,1) << ", " << R_wc(1,2) << ")\n";
        os << "  Cam Z (forward) in world: (" << R_wc(2,0) << ", " << R_wc(2,1) << ", " << R_wc(2,2) << ")\n";
        OPENMVG_LOG_INFO << os.str();
        ++debug_count;
      }
    }
    imu_rotations_loaded_ = true;
    OPENMVG_LOG_INFO << "Loaded IMU rotation priors for "
                     << imu_pose_rotations_.size() << " / "
                     << sfm_data_.GetViews().size() << " views.";
  }

  // Log input graph to the HTML report
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    OPENMVG_LOG_INFO << "Exporting global relative rotation graphs...";
    system::Timer graph_timer;
    // Log a relative view graph
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data_.GetViews().cbegin(), sfm_data_.GetViews().cend(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(matches_provider_->pairWise_matches_));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, "global_relative_rotation_view_graph"),
        putativeGraph);
    }

    // Log a relative pose graph
    {
      std::set<IndexT> set_pose_ids;
      Pair_Set relative_pose_pairs;
      for (const auto & relative_R : vec_relatives_R)
      {
        const Pair relative_pose_indices(relative_R.i, relative_R.j);
        relative_pose_pairs.insert(relative_pose_indices);
        set_pose_ids.insert(relative_R.i);
        set_pose_ids.insert(relative_R.j);
      }
      const std::string sGraph_name = "global_relative_rotation_pose_graph";
      graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, sGraph_name),
        putativeGraph);
      using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(sOut_directory_, "global_relative_rotation_pose_graph", "svg")
         << "\" height=\"600\">\n";

      html_doc_stream_->pushInfo(os.str());
    }
    OPENMVG_LOG_INFO << "Graph export done in (s): " << graph_timer.elapsed();
  }
}

} // namespace sfm
} // namespace openMVG
