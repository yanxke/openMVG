// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_WITH_PRIORS_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_WITH_PRIORS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/solver_essential_five_point.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/robust_estimation/guided_matching.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include <atomic>
#include <iostream>

namespace openMVG {

namespace sfm {
  struct Regions_Provider;
}

namespace matching_image_collection {

/// Shared statistics for geometric filtering
struct GeometricFilterStats {
    std::atomic<size_t> total_pairs{0};
    std::atomic<size_t> rejected_by_heading{0};
    std::atomic<size_t> rejected_by_spot_check{0};
    std::atomic<size_t> rejected_by_min_inliers{0};
    std::atomic<size_t> ransac_attempted{0};
    std::atomic<size_t> total_features{0};
    std::atomic<size_t> ransac_inliers_total{0};
    std::atomic<size_t> matches_recovered{0};
    std::atomic<int> last_printed_percent{-1};
    // For auto-threshold tracking
    mutable std::mutex precision_mutex;
    double total_precision_accum{0.0};
    std::atomic<size_t> precision_count{0};

    // For Stage 1 inlier ratio
    std::atomic<size_t> putative_matches_total{0};
    std::atomic<size_t> stage1_matches_total{0};
    std::atomic<size_t> successful_pairs_count{0};
    // IMU rotation re-estimation (-U) diagnostics
    std::atomic<size_t> reestimate_pairs_attempted{0};
    std::atomic<size_t> reestimate_pairs_reduced{0};
    std::atomic<size_t> reestimate_pairs_rejected{0};
    std::atomic<size_t> reestimate_inliers_before_total{0};
    std::atomic<size_t> reestimate_inliers_after_total{0};
    std::atomic<size_t> reestimate_inliers_rejected_total{0};

    size_t total_expected{0};
    double heading_threshold = 0.0;
    size_t min_inliers_threshold = 0;

    void Print() const {
        if (total_pairs == 0) return;
        size_t h_rej = rejected_by_heading.load();
        size_t s_rej = rejected_by_spot_check.load();
        size_t m_rej = rejected_by_min_inliers.load();
        size_t r_att = ransac_attempted.load();
        size_t f_tot = total_features.load();
        size_t processed = total_pairs.load();
        int percent = (total_expected > 0) ? static_cast<int>(processed * 100 / total_expected) : 0;

        OPENMVG_LOG_INFO << "\n--- Geometric Filter Pipeline Statistics [" << percent << "%] ---";
        OPENMVG_LOG_INFO << "Total pairs processed:    " << processed << " / " << total_expected;
        if (heading_threshold > 0.0) {
            OPENMVG_LOG_INFO << "Rejected by Heading:      " << h_rej << " (" << (h_rej * 100.0 / processed) << "%) [Threshold: " << heading_threshold << " deg]";
        }
        if (s_rej > 0 || heading_threshold > 0.0) {
            OPENMVG_LOG_INFO << "Rejected by Spot Check:   " << s_rej << " (" << (s_rej * 100.0 / processed) << "%)";
        }
        if (min_inliers_threshold > 0) {
            OPENMVG_LOG_INFO << "Rejected by Min Inliers:  " << m_rej << " (" << (m_rej * 100.0 / processed) << "%) [Min count: " << min_inliers_threshold << "]";
        }
        OPENMVG_LOG_INFO << "RANSAC attempted:         " << r_att << " (" << (processed > 0 ? (r_att * 100.0 / processed) : 0.0) << "%)";
            size_t succ_pairs = successful_pairs_count.load();
            if (succ_pairs > 0) {
                size_t p_tot = putative_matches_total.load();
                size_t s1_tot = stage1_matches_total.load();
                size_t r_inl = ransac_inliers_total.load();
                size_t f_tot = matches_recovered.load(); // m_rec now stores FINAL (Stage 3)

                OPENMVG_LOG_INFO << "Geometric verification success (averages per SUCCESSFUL pair):";
                OPENMVG_LOG_INFO << "  Putative matches:          " << (static_cast<double>(p_tot) / succ_pairs);
                OPENMVG_LOG_INFO << "  Stage 1 (IMU filter) ratio: " << (p_tot > 0 ? (s1_tot * 100.0 / p_tot) : 0.0) << "%";
                OPENMVG_LOG_INFO << "  RANSAC Stage 1 inliers:    " << (static_cast<double>(r_inl) / succ_pairs);
                OPENMVG_LOG_INFO << "  Avg RANSAC Inlier ratio:   " << (s1_tot > 0 ? (r_inl * 100.0 / s1_tot) : 0.0) << "%";
                OPENMVG_LOG_INFO << "  Final matches:             " << (static_cast<double>(f_tot) / succ_pairs);
                
                {
                    std::lock_guard<std::mutex> lock(precision_mutex);
                    if (precision_count > 0) {
                        OPENMVG_LOG_INFO << "  Avg Auto Precision found:  " << std::sqrt(total_precision_accum / precision_count) << " px";
                    }
                }

                OPENMVG_LOG_INFO << "Geometric verification totals (across " << succ_pairs << " pairs):";
                OPENMVG_LOG_INFO << "  Total RANSAC Inliers: " << r_inl;
                OPENMVG_LOG_INFO << "  Total FINAL Matches:  " << f_tot;
            }
        size_t reest_attempted = reestimate_pairs_attempted.load();
        if (reest_attempted > 0) {
            const size_t reest_reduced = reestimate_pairs_reduced.load();
            const size_t reest_rejected_pairs = reestimate_pairs_rejected.load();
            const size_t reest_before = reestimate_inliers_before_total.load();
            const size_t reest_after = reestimate_inliers_after_total.load();
            const size_t reest_rejected = reestimate_inliers_rejected_total.load();

            OPENMVG_LOG_INFO << "IMU rotation re-estimation stats (-U):";
            OPENMVG_LOG_INFO << "  Pairs attempted:           " << reest_attempted;
            OPENMVG_LOG_INFO << "  Pairs with fewer inliers:  " << reest_reduced
                             << " (" << (reest_attempted > 0 ? (reest_reduced * 100.0 / reest_attempted) : 0.0) << "%)";
            OPENMVG_LOG_INFO << "  Pairs rejected at re-est.: " << reest_rejected_pairs
                             << " (" << (reest_attempted > 0 ? (reest_rejected_pairs * 100.0 / reest_attempted) : 0.0) << "%)";
            OPENMVG_LOG_INFO << "  Inliers before re-est.:    " << reest_before;
            OPENMVG_LOG_INFO << "  Inliers after re-est.:     " << reest_after;
            OPENMVG_LOG_INFO << "  Inliers further rejected:  " << reest_rejected
                             << " (" << (reest_before > 0 ? (reest_rejected * 100.0 / reest_before) : 0.0) << "%)";
        }
        OPENMVG_LOG_INFO << "------------------------------------------\n";
    }

    void PrintIfNewPercentage() {
        if (total_expected == 0) return;
        if (rejected_by_heading.load() == 0 && rejected_by_spot_check.load() == 0) return;
        int percent = static_cast<int>(total_pairs.load() * 100 / total_expected);
        int last_printed = last_printed_percent.load();
        
        // Calculate the current and last 5% milestone steps
        int current_milestone = percent / 5;
        int last_milestone = (last_printed < 0) ? -1 : last_printed / 5;

        if (current_milestone > last_milestone) {
            if (last_printed_percent.compare_exchange_weak(last_printed, percent)) {
                Print();
            }
        }
    }
};

/// Motion prior configuration for RANSAC-integrated prior scoring
struct MotionPriorConfig {
    double pitch_tol = 10.0;   // degrees - horizontal optical axis
    double roll_tol = 10.0;    // degrees - landscape mode
    double yaw_tol = 15.0;     // degrees - GPS compass consistency
    double alt_tol = 0.2;      // translation ty - flat surface
    double prior_weight = 1.0; // weight of prior penalty (0 = disabled)

    // Number of sample points for cheirality check during pose disambiguation.
    // Higher = more accurate pose recovery but slower.
    // Recommended ranges:
    //   8-16:   Fast, suitable for quick filtering passes
    //   32-64:  Good balance of speed and accuracy (default: 32)
    //   128+:   High accuracy, approaches full-point behavior
    //   0:      Use all points (same accuracy as original, slowest)
    size_t cheirality_samples = 32;

    // GPS headings for current pair (NaN if not available)
    double heading_i = std::numeric_limits<double>::quiet_NaN();
    double heading_j = std::numeric_limits<double>::quiet_NaN();

    // Fast rejection parameters
    double heading_max = 140.0;     // Max allowed yaw difference (degrees, 180 = disabled)
    int    spot_sample_size = 15;  // Number of points to spot-check (0 = disabled)
    int    spot_min_matches = 20;  // Min match count to trigger spot-check
    double spot_dot_thresh  = 0.2; // Parallax orthogonality tolerance
    double spot_min_ratio   = 0.6; // Required agreement ratio
};


/// Essential matrix Kernel adaptor with motion priors integrated into error computation
/// Priors influence model selection by penalizing Essential matrices that violate constraints
template <typename SolverArg,
  typename ErrorArg,
  typename ModelArg = Mat3>
class ACKernelAdaptorEssentialWithPriors
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptorEssentialWithPriors
  (
    const Mat2X &x1, const Mat3X & bearing1, int w1, int h1,
    const Mat2X &x2, const Mat3X & bearing2, int w2, int h2,
    const Mat3 & K1, const Mat3 & K2,
    const MotionPriorConfig & prior_config = MotionPriorConfig()
  ):x1_(x1),
    x2_(x2),
    bearing1_(bearing1),
    bearing2_(bearing2),
    N1_(Mat3::Identity()),
    N2_(Mat3::Identity()), logalpha0_(0.0),
    K1_(K1), K2_(K2),
    prior_config_(prior_config)
  {
    assert(2 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());

    assert(3 == bearing1_.rows());
    assert(bearing1_.rows() == bearing2_.rows());
    assert(bearing1_.cols() == bearing2_.cols());

    logalpha0_ = robust::ACParametrizationHelper<robust::AContrarioParametrizationType::POINT_TO_LINE>::LogAlpha0(w2, h2, 0.5);

    // Pre-select sample indices for cheirality check
    // See MotionPriorConfig::cheirality_samples for documentation on ranges
    const size_t num_points = bearing1_.cols();
    const size_t requested_samples = prior_config_.cheirality_samples;
    const size_t num_samples = (requested_samples == 0)
        ? num_points  // 0 means use all points
        : std::min(requested_samples, num_points);
    cheirality_samples_.reserve(num_samples);
    for (size_t i = 0; i < num_samples; ++i)
    {
      cheirality_samples_.push_back(static_cast<uint32_t>(i * num_points / num_samples));
    }
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<Model> *models
  ) const
  {
    const auto x1 = ExtractColumns(bearing1_, samples);
    const auto x2 = ExtractColumns(bearing2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error
  (
    uint32_t sample,
    const Model &model
  ) const
  {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    return ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }

  void Errors
  (
    const Model & model,
    std::vector<double> & vec_errors
  ) const
  {
    // Standard epipolar errors
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    vec_errors.resize(x1_.cols());
    for (uint32_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));

    // Skip prior penalty if disabled
    if (prior_config_.prior_weight <= 0.0)
      return;

    // Fast pose recovery using MotionFromEssential (O(1) SVD)
    std::vector<geometry::Pose3> poses;
    MotionFromEssential(model, &poses);

    if (poses.size() != 4)
      return;  // Unexpected, keep original errors

    // Determine correct pose via cheirality check on small sample
    const geometry::Pose3 pose1(Mat3::Identity(), Vec3::Zero());
    int best_pose_idx = -1;
    int best_count = 0;

    for (size_t pose_idx = 0; pose_idx < poses.size(); ++pose_idx)
    {
      const geometry::Pose3 &pose2 = poses[pose_idx];
      int count = 0;

      for (const uint32_t idx : cheirality_samples_)
      {
        Vec3 X;
        if (Triangulate2View(
          pose1.rotation(), pose1.translation(), bearing1_.col(idx),
          pose2.rotation(), pose2.translation(), bearing2_.col(idx),
          X, ETriangulationMethod::DEFAULT))
        {
          ++count;
        }
      }

      if (count > best_count)
      {
        best_count = count;
        best_pose_idx = static_cast<int>(pose_idx);
      }
    }

    // If no valid pose found, keep original errors
    if (best_pose_idx < 0 || best_count == 0)
      return;

    const geometry::Pose3 &best_pose = poses[best_pose_idx];

    // Extract Euler angles from rotation matrix
    // Convention: camera coords x-right, y-down, z-forward
    const Mat3 R = best_pose.rotation();
    const double pitch = std::asin(clamp(R(1, 2), -1.0, 1.0)) * 180.0 / M_PI;
    const double yaw = std::atan2(-R(0, 2), R(2, 2)) * 180.0 / M_PI;
    const double roll = std::atan2(-R(1, 0), R(1, 1)) * 180.0 / M_PI;
    const double ty = best_pose.translation()(1);

    // Compute prior penalty (normalized violations, soft penalty beyond tolerance)
    double penalty = 0.0;

    if (prior_config_.pitch_tol > 0.0)
    {
      double excess = std::max(0.0, std::abs(pitch) / prior_config_.pitch_tol - 1.0);
      penalty += excess * excess;  // Quadratic penalty for smooth gradient
    }

    if (prior_config_.roll_tol > 0.0)
    {
      double excess = std::max(0.0, std::abs(roll) / prior_config_.roll_tol - 1.0);
      penalty += excess * excess;
    }

    if (prior_config_.alt_tol > 0.0)
    {
      double excess = std::max(0.0, std::abs(ty) / prior_config_.alt_tol - 1.0);
      penalty += excess * excess;
    }

    // Yaw penalty (only if GPS headings available)
    if (!std::isnan(prior_config_.heading_i) &&
        !std::isnan(prior_config_.heading_j) &&
        prior_config_.yaw_tol > 0.0)
    {
      double expected_yaw = prior_config_.heading_j - prior_config_.heading_i;
      // Normalize to [-180, 180]
      while (expected_yaw > 180.0) expected_yaw -= 360.0;
      while (expected_yaw < -180.0) expected_yaw += 360.0;

      double yaw_diff = std::abs(yaw - expected_yaw);
      if (yaw_diff > 180.0) yaw_diff = 360.0 - yaw_diff;

      double excess = std::max(0.0, yaw_diff / prior_config_.yaw_tol - 1.0);
      penalty += excess * excess;
    }

    // Add prior penalty to all errors
    // This makes models violating priors less attractive to ACRANSAC
    if (penalty > 0.0)
    {
      const double penalty_term = prior_config_.prior_weight * penalty;
      for (auto& err : vec_errors)
        err += penalty_term;
    }
  }

  size_t NumSamples() const { return x1_.cols(); }
  void Unnormalize(Model * model) const {}
  double logalpha0() const {return logalpha0_;}
  double multError() const {return robust::ACParametrizationHelper<robust::AContrarioParametrizationType::POINT_TO_LINE>::MultError();}
  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return N2_;}
  double unormalizeError(double val) const { return val; }

private:
  Mat2X x1_, x2_;             // image points
  Mat3X bearing1_, bearing2_; // bearing vectors
  Mat3 N1_, N2_;              // Matrix used to normalize data
  double logalpha0_;          // Alpha0 is used to make the error adaptive to the image size
  Mat3 K1_, K2_;              // Intrinsic camera parameter
  MotionPriorConfig prior_config_;  // Motion prior tolerances
  std::vector<uint32_t> cheirality_samples_;  // Pre-selected sample indices for fast cheirality check
};


//-- A contrario essential matrix estimation with motion priors template functor
//   used for filtering pairs of putative correspondences
//
//   Strategy: Use custom kernel that integrates prior penalties into RANSAC
//   error computation, so ACRANSAC naturally prefers Essential matrices
//   that satisfy motion constraints.
struct GeometricFilter_EMatrix_AC_WithPriors
{
  // =========================================================================
  // FAST REJECTION PARAMETERS (TUNE HERE)
  // =========================================================================
  GeometricFilter_EMatrix_AC_WithPriors
  (
    double dPrecision = std::numeric_limits<double>::infinity(),
    uint32_t iteration = 1024,
    const MotionPriorConfig & default_prior_config = MotionPriorConfig(),
    const std::map<IndexT, double> * map_headings = nullptr,
    size_t total_expected = 0,
    size_t min_inliers = 0
  ):
    m_dPrecision(dPrecision),
    m_stIteration(iteration),
    m_E(Mat3::Identity()),
    m_dPrecision_robust(std::numeric_limits<double>::infinity()),
    m_default_prior_config(default_prior_config),
    m_map_headings(map_headings),
    m_min_inliers(min_inliers),
    m_stats(std::make_shared<GeometricFilterStats>())
  {
    m_stats->total_expected = total_expected;
    m_stats->heading_threshold = m_default_prior_config.heading_max;
    m_stats->min_inliers_threshold = min_inliers;
  }

  ~GeometricFilter_EMatrix_AC_WithPriors()
  {
    // Print stats only once (when the last shared copy is destroyed)
    if (m_stats.use_count() == 1)
    {
      m_stats->Print();
    }
  }

  /// Robust fitting of the ESSENTIAL matrix with motion priors integrated into RANSAC
  template<typename Regions_or_Features_ProviderT>
  bool Robust_estimation
  (
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<Regions_or_Features_ProviderT> & regions_provider,
    const Pair pairIndex,
    const matching::IndMatches & vec_PutativeMatches,
    matching::IndMatches & geometric_inliers)
  {
    geometric_inliers.clear();
    m_stats->total_pairs++;
    m_stats->PrintIfNewPercentage();

    // Get back corresponding view index
    const IndexT
      iIndex = pairIndex.first,
      jIndex = pairIndex.second;

    //--
    // Reject pair with missing Intrinsic information
    //--

    const sfm::View
      * view_I = sfm_data->views.at(iIndex).get(),
      * view_J = sfm_data->views.at(jIndex).get();

     // Check that valid cameras can be retrieved for the pair of views
    const cameras::IntrinsicBase
      * cam_I =
        sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr,
      * cam_J =
        sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

    if (!cam_I || !cam_J)
    {
      return false;
    }
    if (!isPinhole(cam_I->getType()) || !isPinhole(cam_J->getType()))
    {
      return false;
    }

    //--
    // Prepare prior config for this pair
    //--
    MotionPriorConfig pair_prior_config = m_default_prior_config;
    if (m_map_headings)
    {
      auto it_i = m_map_headings->find(iIndex);
      auto it_j = m_map_headings->find(jIndex);
      if (it_i != m_map_headings->end())
        pair_prior_config.heading_i = it_i->second;
      if (it_j != m_map_headings->end())
        pair_prior_config.heading_j = it_j->second;
    }

    //--
    // Fast Heading Check (Max Difference)
    //--
    if (pair_prior_config.heading_max < 180.0 &&
        !std::isnan(pair_prior_config.heading_i) && !std::isnan(pair_prior_config.heading_j))
    {
      double diff = std::abs(pair_prior_config.heading_i - pair_prior_config.heading_j);
      while (diff > 180.0) diff = 360.0 - diff;
      if (diff > pair_prior_config.heading_max)
      {
        m_stats->rejected_by_heading++;
        return false; // Skip pairs facing mostly opposite directions
      }
    }

    //--
    // Get corresponding point regions arrays
    //--

    Mat2X xI, xJ;
    MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

    const cameras::Pinhole_Intrinsic
      * ptrPinhole_I = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_I),
      * ptrPinhole_J = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_J);

    //--
    // Fast Rejection Spot-Check
    //--
    // If we have reasonably good priors, we can check if a small sample of matches
    // is consistent with the prior rotation. If not, the pair probably doesn't overlap.
    if (pair_prior_config.spot_sample_size > 0 &&
        xI.cols() > (size_t)pair_prior_config.spot_min_matches && 
        !std::isnan(pair_prior_config.heading_i))
    {
      // Assume R is roughly correct as per priors (Yaw only, level)
      double expected_yaw = (pair_prior_config.heading_j - pair_prior_config.heading_i) * M_PI / 180.0;
      Mat3 R_prior;
      R_prior << std::cos(expected_yaw), 0, -std::sin(expected_yaw),
                 0, 1, 0,
                 std::sin(expected_yaw), 0, std::cos(expected_yaw);
      
      const Mat3 K1_inv = ptrPinhole_I->K().inverse();
      const Mat3 K2_inv = ptrPinhole_J->K().inverse();

      // Check coplanarity ofparallax vectors: v_k = x' x (R x)
      // All v_k must be orthogonal to translation t.
      std::vector<Vec3> vs;
      const int sample_size = std::min<int>(pair_prior_config.spot_sample_size, (int)xI.cols());
      vs.reserve(sample_size);
      for (int k = 0; k < sample_size; ++k)
      {
        Vec3 p1 = K1_inv * Vec3(xI(0, k), xI(1, k), 1.0);
        Vec3 p2 = K2_inv * Vec3(xJ(0, k), xJ(1, k), 1.0);
        vs.push_back(p2.cross(R_prior * p1));
      }

      // 1-point test: assume t is cross product of pairs of vectors
      bool plausible = false;
      for (int i = 0; i < sample_size / 2 && !plausible; ++i)
      {
        for (int j = i + 1; j < sample_size / 2 + 1; ++j)
        {
          Vec3 t = vs[i].cross(vs[j]);
          if (t.norm() < 1e-6) continue;
          t.normalize();
          
          int inliers = 0;
          for (const auto& v : vs)
          {
            if (std::abs(v.dot(t)) < pair_prior_config.spot_dot_thresh * v.norm()) inliers++;
          }
          if (inliers >= static_cast<int>(pair_prior_config.spot_min_ratio * sample_size)) { plausible = true; break; }
        }
      }

      if (!plausible)
      {
        m_stats->rejected_by_spot_check++;
        return false;
      }
    }

    //--
    // Robust estimation using kernel with integrated priors
    //--

    using KernelType =
      ACKernelAdaptorEssentialWithPriors<
        openMVG::essential::kernel::FivePointSolver,
        openMVG::fundamental::kernel::EpipolarDistanceError,
        Mat3>;

    KernelType kernel(
      xI, (*cam_I)(xI),
      sfm_data->GetViews().at(iIndex)->ui_width, sfm_data->GetViews().at(iIndex)->ui_height,
      xJ, (*cam_J)(xJ),
      sfm_data->GetViews().at(jIndex)->ui_width, sfm_data->GetViews().at(jIndex)->ui_height,
      ptrPinhole_I->K(), ptrPinhole_J->K(),
      pair_prior_config);

    // Robustly estimate the Essential matrix with A Contrario ransac
    const double upper_bound_precision = Square(m_dPrecision);
    std::vector<uint32_t> vec_inliers;
    m_stats->ransac_attempted++;
    const auto ACRansacOut =
      openMVG::robust::ACRANSAC(kernel, vec_inliers, m_stIteration, &m_E, upper_bound_precision);

    if (vec_inliers.size() < std::max<size_t>(m_min_inliers, (size_t)(KernelType::MINIMUM_SAMPLES * 2.5)))
    {
      if (vec_inliers.size() > 0)
        m_stats->rejected_by_min_inliers++;
      vec_inliers.clear();
      return false;
    }

    m_dPrecision_robust = ACRansacOut.first;

    // Update geometric_inliers
    geometric_inliers.reserve(vec_inliers.size());
    for (const uint32_t & index : vec_inliers)
    {
      geometric_inliers.push_back( vec_PutativeMatches[index] );
    }
    return true;
  }

  bool Geometry_guided_matching
  (
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair pairIndex,
    const double dDistanceRatio,
    matching::IndMatches & matches
  )
  {
    if (m_dPrecision_robust != std::numeric_limits<double>::infinity())
    {
      // Get back corresponding view index
      const IndexT iIndex = pairIndex.first;
      const IndexT jIndex = pairIndex.second;

      const sfm::View * view_I = sfm_data->views.at(iIndex).get();
      const sfm::View * view_J = sfm_data->views.at(jIndex).get();

      // Check that valid cameras can be retrieved for the pair of views
      const cameras::IntrinsicBase * cam_I =
        sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr;
      const cameras::IntrinsicBase * cam_J =
        sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

      if (!cam_I || !cam_J)
      {
        return false;
      }

      if ( !isPinhole(cam_I->getType()) || !isPinhole(cam_J->getType()))
      {
        return false;
      }

      const auto * ptrPinhole_I = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_I);
      const auto * ptrPinhole_J = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_J);

      Mat3 F;
      FundamentalFromEssential(m_E, ptrPinhole_I->K(), ptrPinhole_J->K(), &F);

      const std::shared_ptr<features::Regions>
        regionsI = regions_provider->get(iIndex),
        regionsJ = regions_provider->get(jIndex);

      geometry_aware::GuidedMatching<
        Mat3,
        openMVG::fundamental::kernel::EpipolarDistanceError>(
          F,
          cam_I, *regionsI,
          cam_J, *regionsJ,
          Square(m_dPrecision_robust), Square(dDistanceRatio),
          matches);
    }
    return matches.size() != 0;
  }

  double m_dPrecision;    // upper_bound precision used for robust estimation
  uint32_t m_stIteration; // maximal number of iteration for robust estimation
  //
  //-- Stored data
  Mat3 m_E;
  double m_dPrecision_robust;
  MotionPriorConfig m_default_prior_config;
  const std::map<IndexT, double> * m_map_headings;
  size_t m_min_inliers;
  std::shared_ptr<GeometricFilterStats> m_stats;
};

} //namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_WITH_PRIORS_HPP
