// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/graph/connectedComponent.hpp"
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

constexpr double kCloseTimeSeconds = 0.75;
// Adaptive angle filter thresholds used when close-time adaptive filtering is enabled.
// - kDefaultMinTriangulationAngleDeg: legacy minimum triangulation angle.
// - kCloseTimeMinTriangulationAngleDeg: relaxed minimum angle for close-time observation pairs.
constexpr double kDefaultMinTriangulationAngleDeg = 2.0;
constexpr double kCloseTimeMinTriangulationAngleDeg = 1.0;
constexpr double kDefaultPixelResidualThresholdPx = 4.0;
constexpr double kCloseTimePixelResidualThresholdPx = 8.0;

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

bool AreViewsCloseInTime(
  const openMVG::sfm::View * view_a,
  const openMVG::sfm::View * view_b)
{
  if (!view_a || !view_b ||
      !view_a->b_has_capture_time_epoch_ || !view_b->b_has_capture_time_epoch_)
  {
    return false;
  }
  return std::abs(view_a->capture_time_epoch_ - view_b->capture_time_epoch_) < kCloseTimeSeconds;
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

std::string EscapeCsvField(const std::string & value)
{
  if (value.find_first_of(",\"\n\r") == std::string::npos)
  {
    return value;
  }

  std::string escaped = "\"";
  for (const char c : value)
  {
    if (c == '"')
    {
      escaped += "\"\"";
    }
    else
    {
      escaped += c;
    }
  }
  escaped += "\"";
  return escaped;
}

void WriteNonMainComponentsCsv(
  const std::string & path,
  const openMVG::sfm::SfM_Data & sfm_data,
  const openMVG::Pair_Set & pairs)
{
  if (pairs.empty())
  {
    return;
  }

  openMVG::graph::indexedGraph putative_graph(pairs);
  using Graph = openMVG::graph::indexedGraph::GraphT;
  using EdgeMap = Graph::EdgeMap<bool>;
  EdgeMap cut_map(putative_graph.g);

  if (lemon::biEdgeConnectedCutEdges(putative_graph.g, cut_map) > 0)
  {
    using EdgeIterator = Graph::EdgeIt;
    EdgeIterator edge_it(putative_graph.g);
    for (EdgeMap::MapIt map_it(cut_map); map_it != lemon::INVALID; ++map_it, ++edge_it)
    {
      if (*map_it)
      {
        putative_graph.g.erase(edge_it);
      }
    }
  }

  const std::map<openMVG::IndexT, std::set<Graph::Node>> components =
    openMVG::graph::exportGraphToMapSubgraphs<Graph, openMVG::IndexT>(putative_graph.g);
  if (components.size() <= 1)
  {
    return;
  }

  openMVG::IndexT largest_component_id = openMVG::UndefinedIndexT;
  size_t largest_size = 0;
  for (const auto & component_it : components)
  {
    if (component_it.second.size() > largest_size)
    {
      largest_size = component_it.second.size();
      largest_component_id = component_it.first;
    }
  }

  OPENMVG_LOG_INFO << "Writing CSV: " << path;
  std::ofstream csv(path.c_str());
  if (!csv.is_open())
  {
    OPENMVG_LOG_WARNING << "Cannot write non-main-components CSV: " << path;
    return;
  }

  csv << "component_id,camera_id,image_id,image_name\n";
  for (const auto & component_it : components)
  {
    if (component_it.first == largest_component_id)
    {
      continue;
    }

    for (const auto & node : component_it.second)
    {
      const openMVG::IndexT image_id = (*putative_graph.node_map_id)[node];
      openMVG::IndexT camera_id = openMVG::UndefinedIndexT;
      std::string image_name;
      const auto view_it = sfm_data.GetViews().find(image_id);
      if (view_it != sfm_data.GetViews().end() && view_it->second)
      {
        camera_id = view_it->second->id_pose;
        image_name = view_it->second->s_Img_path;
      }

      csv << component_it.first << ","
          << camera_id << ","
          << image_id << ","
          << EscapeCsvField(image_name) << "\n";
    }
  }
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
                               const openMVG::rotation_averaging::RelativeRotations & relatives_R,
                               const openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3> * transformed_imu_pose_rotations = nullptr)
{
  OPENMVG_LOG_INFO << "Writing CSV: " << path;
  std::ofstream csv(path.c_str());
  if (!csv)
    return;

  openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3> imu_pose_rotations;
  for (const auto & view_it : sfm_data.GetViews())
  {
    const openMVG::sfm::View * view = view_it.second.get();
    if (!view || view->id_pose == openMVG::UndefinedIndexT)
      continue;

    const auto * view_priors =
      dynamic_cast<const openMVG::sfm::ViewPriors *>(view);
    if (!view_priors || !view_priors->b_has_imu_rotation_)
      continue;

    // Prefer transformed IMU rotations (same convention used by global SfM),
    // and only fallback to raw view priors when transformed rotations are unavailable.
    if (transformed_imu_pose_rotations)
    {
      const auto transformed_it = transformed_imu_pose_rotations->find(view->id_pose);
      if (transformed_it != transformed_imu_pose_rotations->end())
      {
        imu_pose_rotations[view->id_pose] = transformed_it->second;
        continue;
      }
    }
    imu_pose_rotations[view->id_pose] = view_priors->imu_rotation_;
  }

  const auto pose_to_view = BuildPoseToViewMap(sfm_data);

  csv << "pose_i,pose_j,view_i,view_j,image_i,image_j,qw,qx,qy,qz,weight,"
         "imu_relative_rotation_error_deg,imu_relative_rotation_error_mode\n";
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

    std::string err_deg_csv;
    std::string err_mode_csv;
    const auto imu_i_it = imu_pose_rotations.find(rel.i);
    const auto imu_j_it = imu_pose_rotations.find(rel.j);
    if (imu_i_it != imu_pose_rotations.end() && imu_j_it != imu_pose_rotations.end())
    {
      const openMVG::Mat3 rel_imu = imu_j_it->second * imu_i_it->second.transpose();
      const double err_direct_deg = RotationAngularErrorDeg(rel.Rij, rel_imu);
      const double err_transpose_deg = RotationAngularErrorDeg(rel.Rij.transpose(), rel_imu);
      const double err_deg = std::min(err_direct_deg, err_transpose_deg);
      const char * mode = (err_direct_deg <= err_transpose_deg) ? "direct" : "transpose";

      std::ostringstream err_ss;
      err_ss << std::fixed << std::setprecision(10) << err_deg;
      err_deg_csv = err_ss.str();
      err_mode_csv = mode;
    }

    csv << rel.i << "," << rel.j << ","
        << view_i << "," << view_j << ","
        << image_i << "," << image_j << ","
        << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << ","
        << rel.weight << ","
        << err_deg_csv << ","
        << err_mode_csv << "\n";
  }
}

void DumpPairWiseMatchesJson(
  const std::string & path,
  const openMVG::matching::PairWiseMatches & pairwise_matches)
{
  OPENMVG_LOG_INFO << "Writing debug JSON: " << path;
  std::ofstream os(path.c_str());
  if (!os.is_open())
  {
    OPENMVG_LOG_WARNING << "Cannot write debug matches JSON: " << path;
    return;
  }

  os << "{\n";
  os << "  \"num_pairs\": " << pairwise_matches.size() << ",\n";
  os << "  \"pairs\": [\n";
  bool first_pair = true;
  for (const auto & pair_entry : pairwise_matches)
  {
    if (!first_pair)
      os << ",\n";
    first_pair = false;
    const openMVG::Pair & pair = pair_entry.first;
    const openMVG::matching::IndMatches & matches = pair_entry.second;
    os << "    {\"i\": " << pair.first
       << ", \"j\": " << pair.second
       << ", \"num_matches\": " << matches.size()
       << ", \"matches\": [";
    bool first_match = true;
    for (const auto & m : matches)
    {
      if (!first_match)
        os << ", ";
      first_match = false;
      os << "[" << m.i_ << ", " << m.j_ << "]";
    }
    os << "]}";
  }
  os << "\n  ]\n";
  os << "}\n";
}

void DumpTracksJson(
  const std::string & path,
  const openMVG::tracks::STLMAPTracks & tracks)
{
  OPENMVG_LOG_INFO << "Writing debug JSON: " << path;
  std::ofstream os(path.c_str());
  if (!os.is_open())
  {
    OPENMVG_LOG_WARNING << "Cannot write debug tracks JSON: " << path;
    return;
  }

  os << "{\n";
  os << "  \"num_tracks\": " << tracks.size() << ",\n";
  os << "  \"tracks\": [\n";
  bool first_track = true;
  for (const auto & t : tracks)
  {
    if (!first_track)
      os << ",\n";
    first_track = false;
    os << "    {\"track_id\": " << t.first << ", \"observations\": [";
    bool first_obs = true;
    for (const auto & obs : t.second)
    {
      if (!first_obs)
        os << ", ";
      first_obs = false;
      os << "{\"view_id\": " << obs.first << ", \"feat_id\": " << obs.second << "}";
    }
    os << "]}";
  }
  os << "\n  ]\n";
  os << "}\n";
}

void DumpStructureObservationsJson(
  const std::string & path,
  const openMVG::sfm::SfM_Data & sfm_data)
{
  OPENMVG_LOG_INFO << "Writing debug JSON: " << path;
  std::ofstream os(path.c_str());
  if (!os.is_open())
  {
    OPENMVG_LOG_WARNING << "Cannot write debug structure JSON: " << path;
    return;
  }

  const openMVG::sfm::Landmarks & landmarks = sfm_data.GetLandmarks();
  os << "{\n";
  os << "  \"num_landmarks\": " << landmarks.size() << ",\n";
  os << "  \"landmarks\": [\n";
  bool first_landmark = true;
  for (const auto & lm : landmarks)
  {
    if (!first_landmark)
      os << ",\n";
    first_landmark = false;
    const openMVG::sfm::Observations & obs = lm.second.obs;
    os << "    {\"landmark_id\": " << lm.first
       << ", \"num_observations\": " << obs.size()
       << ", \"observations\": [";
    bool first_obs = true;
    for (const auto & o : obs)
    {
      if (!first_obs)
        os << ", ";
      first_obs = false;
      os << "{\"view_id\": " << o.first << ", \"feat_id\": " << o.second.id_feat << "}";
    }
    os << "]}";
  }
  os << "\n  ]\n";
  os << "}\n";
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

void GlobalSfMReconstructionEngine_RelativeMotions::SetCloseTimeAdaptiveFiltering(bool enabled)
{
  close_time_adaptive_filtering_ = enabled;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {

  //-------------------
  // Keep only the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = matches_provider_->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
    if (!sOut_directory_.empty())
    {
      WriteNonMainComponentsCsv(
        stlplus::create_filespec(sOut_directory_, "sfm_non_main_components", "csv"),
        sfm_data_,
        pairs);
    }
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
  if (!sOut_directory_.empty())
  {
    DumpPairWiseMatchesJson(
      stlplus::create_filespec(sOut_directory_, "sfm_debug_triplet_wise_matches", "json"),
      tripletWise_matches);
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
    WriteRelativeRotationsCsv(csv_path, sfm_data_, relatives_R, &imu_pose_rotations_);
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
  size_t imu_filtered_edge_count = 0;
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

  Hash_Map<IndexT, double> pose_timestamps;
  Hash_Map<IndexT, std::string> pose_img_names;
  for (const auto & view_it : sfm_data_.GetViews())
  {
    const auto & view = view_it.second;
    if (view->id_pose != UndefinedIndexT)
    {
      if (view->b_has_capture_time_epoch_)
        pose_timestamps[view->id_pose] = view->capture_time_epoch_;
      pose_img_names[view->id_pose] = view->s_Img_path;
    }
  }

  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  auto run_rotation_averaging = [&](const rotation_averaging::RelativeRotations & rels) -> bool
  {
    system::Timer t;
    const bool ok = rotation_averaging_solver.Run(
      eRotation_averaging_method_, eRelativeRotationInferenceMethod,
      rels, global_rotations, world_pose_id,
      &pose_timestamps, &pose_img_names);

    OPENMVG_LOG_INFO
      << "Found #global_rotations: " << global_rotations.size() << "\n"
      << "Timing: " << t.elapsed() << " seconds";
    return ok;
  };

  bool b_rotation_averaging = run_rotation_averaging(relatives_with_priors);
  if (b_rotation_averaging && imu_rotation_filter_outliers_ &&
      !imu_pose_rotations_.empty() &&
      imu_rotation_max_error_deg_ > 0.0)
  {
    if (imu_rotation_weight_ > 0.0)
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
    else
    {
      // Independent IMU filtering path (works even with zero IMU prior weight):
      // remove relative-rotation edges that strongly disagree with IMU pairwise rotation.
      rotation_averaging::RelativeRotations filtered_relatives;
      filtered_relatives.reserve(relatives_R.size());
      size_t comparable_edges = 0;

      for (const auto & rel : relatives_R)
      {
        bool keep = true;
        const auto imu_i_it = imu_pose_rotations_.find(rel.i);
        const auto imu_j_it = imu_pose_rotations_.find(rel.j);
        if (imu_i_it != imu_pose_rotations_.end() && imu_j_it != imu_pose_rotations_.end())
        {
          ++comparable_edges;
          const Mat3 rel_imu = imu_j_it->second * imu_i_it->second.transpose();
          const double err_direct_deg = RotationAngularErrorDeg(rel.Rij, rel_imu);
          const double err_transpose_deg = RotationAngularErrorDeg(rel.Rij.transpose(), rel_imu);
          const double err_deg = std::min(err_direct_deg, err_transpose_deg);
          if (err_deg > imu_rotation_max_error_deg_)
          {
            keep = false;
            ++imu_filtered_edge_count;
          }
        }
        if (keep)
          filtered_relatives.push_back(rel);
      }

      if (imu_filtered_edge_count > 0)
      {
        OPENMVG_LOG_WARNING
          << "IMU relative-rotation filtering (weight=0): removing "
          << imu_filtered_edge_count << " / " << comparable_edges
          << " comparable relative-rotation edges with error > "
          << imu_rotation_max_error_deg_ << " deg.";
        b_rotation_averaging = run_rotation_averaging(filtered_relatives);
      }
      else
      {
        OPENMVG_LOG_INFO
          << "IMU relative-rotation filtering (weight=0): removed 0 / "
          << comparable_edges << " comparable relative-rotation edges with error > "
          << imu_rotation_max_error_deg_ << " deg.";
      }
    }
  }

  if (b_rotation_averaging && !imu_pose_rotations_.empty())
  {
    OPENMVG_LOG_INFO
      << "IMU rotation prior settings: weight=" << imu_rotation_weight_
      << ", max_error_deg=" << imu_rotation_max_error_deg_
      << ", filter_outliers=" << (imu_rotation_filter_outliers_ ? "true" : "false")
      << ", filtered_priors=" << imu_filtered_count
      << ", filtered_edges=" << imu_filtered_edge_count;

    auto median_of = [](std::vector<double> values) -> double
    {
      if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();
      const size_t mid = values.size() / 2;
      std::nth_element(values.begin(), values.begin() + mid, values.end());
      double med = values[mid];
      if (values.size() % 2 == 0)
      {
        const auto max_low_it = std::max_element(values.begin(), values.begin() + mid);
        med = 0.5 * (med + *max_low_it);
      }
      return med;
    };

    std::vector<std::pair<Mat3, Mat3>> imu_pairs_base;
    imu_pairs_base.reserve(imu_pose_rotations_.size());
    std::vector<double> imu_errors_direct_deg;
    std::vector<double> imu_errors_transpose_deg;
    imu_errors_direct_deg.reserve(imu_pose_rotations_.size());
    imu_errors_transpose_deg.reserve(imu_pose_rotations_.size());
    size_t missing_global = 0;
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
      const Mat3 & R_global = it->second;
      const Mat3 & R_imu = imu_it.second;
      imu_pairs_base.emplace_back(R_global, R_imu);
      imu_errors_direct_deg.push_back(RotationAngularErrorDeg(R_global, R_imu));
      imu_errors_transpose_deg.push_back(RotationAngularErrorDeg(R_global, R_imu.transpose()));
    }

    const double direct_median = median_of(imu_errors_direct_deg);
    const double transpose_median = median_of(imu_errors_transpose_deg);
    const bool use_transpose_convention =
      !std::isnan(transpose_median) &&
      (std::isnan(direct_median) || transpose_median < direct_median);

    const std::vector<double> & imu_errors_deg =
      use_transpose_convention ? imu_errors_transpose_deg : imu_errors_direct_deg;

    std::vector<std::pair<Mat3, Mat3>> imu_pairs;
    imu_pairs.reserve(imu_pairs_base.size());
    for (const auto & pair : imu_pairs_base)
    {
      imu_pairs.emplace_back(
        pair.first,
        use_transpose_convention ? pair.second.transpose() : pair.second);
    }

    OPENMVG_LOG_INFO
      << "IMU absolute convention auto-selected: "
      << (use_transpose_convention ? "transpose" : "direct")
      << " (median_direct_deg=" << direct_median
      << ", median_transpose_deg=" << transpose_median << ")";

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
    sOut_directory_,
    close_time_adaptive_filtering_);

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
    if (!sOut_directory_.empty())
    {
      DumpTracksJson(
        stlplus::create_filespec(sOut_directory_, "sfm_debug_tracks_filtered_len3", "json"),
        map_selectedTracks);
    }

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
    if (!sOut_directory_.empty())
    {
      DumpStructureObservationsJson(
        stlplus::create_filespec(sOut_directory_, "sfm_debug_structure_after_triangulation", "json"),
        sfm_data_);
    }

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
  PixelResidualAdaptiveStats pixel_adaptive_stats;
  AngleErrorAdaptiveStats angle_adaptive_stats;
  if (close_time_adaptive_filtering_)
  {
    RemoveOutliers_PixelResidualErrorAdaptive(
      sfm_data_,
      kDefaultPixelResidualThresholdPx,
      kCloseTimePixelResidualThresholdPx,
      kCloseTimeSeconds,
      &pixel_adaptive_stats);
  }
  else
  {
    RemoveOutliers_PixelResidualError(sfm_data_, kDefaultPixelResidualThresholdPx);
  }
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  if (!sOut_directory_.empty())
  {
    DumpStructureObservationsJson(
      stlplus::create_filespec(sOut_directory_, "sfm_debug_structure_after_pixel_filter", "json"),
      sfm_data_);
  }
  if (close_time_adaptive_filtering_)
  {
    RemoveOutliers_AngleErrorAdaptive(
      sfm_data_,
      kDefaultMinTriangulationAngleDeg,
      kCloseTimeMinTriangulationAngleDeg,
      kCloseTimeSeconds,
      &angle_adaptive_stats);
  }
  else
  {
    RemoveOutliers_AngleError(sfm_data_, kDefaultMinTriangulationAngleDeg);
  }
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  if (!sOut_directory_.empty())
  {
    DumpStructureObservationsJson(
      stlplus::create_filespec(sOut_directory_, "sfm_debug_structure_after_angle_filter", "json"),
      sfm_data_);
  }
  OPENMVG_LOG_INFO << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter;
  if (close_time_adaptive_filtering_)
  {
    OPENMVG_LOG_INFO
      << "Close-time adaptive outlier summary:\n"
      << "-- pixel outliers by old logic: " << pixel_adaptive_stats.old_logic_outlier_observations << "\n"
      << "-- pixel outliers rescued by relaxed threshold: " << pixel_adaptive_stats.rescued_observations << "\n"
      << "-- pixel outliers removed by active logic: " << pixel_adaptive_stats.removed_observations << "\n"
      << "-- angular outlier tracks by old logic: " << angle_adaptive_stats.old_logic_outlier_tracks << "\n"
      << "-- angular outlier tracks rescued by relaxed threshold: " << angle_adaptive_stats.rescued_tracks << "\n"
      << "-- angular outlier tracks removed by active logic: " << angle_adaptive_stats.removed_tracks;
  }

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLengthDefault = 3; // normal
  bool cleaned = false;
  if (close_time_adaptive_filtering_)
  {
    const IndexT minTrackLengthCloseTime = minTrackLengthDefault;
    cleaned = eraseUnstablePosesAndObservationsAdaptive(
      sfm_data_,
      minPointPerPose,
      minTrackLengthDefault,
      minTrackLengthCloseTime,
      kCloseTimeSeconds);
  }
  else
  {
    cleaned = eraseUnstablePosesAndObservations(
      sfm_data_,
      minPointPerPose,
      minTrackLengthDefault);
  }
  if (cleaned)
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    OPENMVG_LOG_INFO << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }
  if (!sOut_directory_.empty())
  {
    DumpStructureObservationsJson(
      stlplus::create_filespec(sOut_directory_, "sfm_debug_structure_after_outlier_cleanup", "json"),
      sfm_data_);
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
    DumpStructureObservationsJson(
      stlplus::create_filespec(sOut_directory_, "sfm_debug_structure_after_final_ba", "json"),
      sfm_data_);
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

      const sfm::ViewPriors * view_priors =
        dynamic_cast<const sfm::ViewPriors *>(view_it.second.get());
      if (!view_priors || !view_priors->b_has_imu_rotation_)
        continue;

      const Mat3 & R_dw = view_priors->imu_rotation_;

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
        os << "R_dw (Device->World from sfm_data ViewPriors):\n" << R_dw << "\n\n";
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
