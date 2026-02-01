// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_IMU_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_IMU_HPP

#include "openMVG/matching_image_collection/E_ACRobust_WithPriors.hpp"
#include "openMVG/multiview/solver_essential_two_point_imu.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include <map>

namespace openMVG {
namespace matching_image_collection {

/**
 * @brief Kernel adaptor for 2-point IMU-guided Essential matrix estimation
 * 
 * This kernel uses known IMU rotations to reduce Essential matrix estimation
 * from 5 points to 2 points, dramatically speeding up RANSAC.
 */
template <typename ErrorArg, typename ModelArg = Mat3>
class ACKernelAdaptorEssentialImu
{
public:
  using Solver = essential::kernel::TwoPointImuSolver;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptorEssentialImu(
    const Mat2X &x1, const Mat3X & bearing1, int w1, int h1,
    const Mat2X &x2, const Mat3X & bearing2, int w2, int h2,
    const Mat3 & K1, const Mat3 & K2,
    const Mat3 & R1_wc,  // World-to-Camera rotation for view 1 (from IMU)
    const Mat3 & R2_wc   // World-to-Camera rotation for view 2 (from IMU)
  ):x1_(x1),
    x2_(x2),
    bearing1_(bearing1),
    bearing2_(bearing2),
    N1_(Mat3::Identity()),
    N2_(Mat3::Identity()),
    logalpha0_(0.0),
    K1_(K1),
    K2_(K2),
    R_relative_(R2_wc * R1_wc.transpose())  // R = R2 · R1^T
  {
    assert(2 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());

    assert(3 == bearing1_.rows());
    assert(bearing1_.rows() == bearing2_.rows());
    assert(bearing1_.cols() == bearing2_.cols());

    logalpha0_ = robust::ACParametrizationHelper<robust::AContrarioParametrizationType::POINT_TO_LINE>::LogAlpha0(w2, h2, 0.5);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };  // 2
  enum { MAX_MODELS = Solver::MAX_MODELS };            // 1

  void Fit(
    const std::vector<uint32_t> &samples,
    std::vector<Model> *models) const
  {
    const auto x1 = ExtractColumns(bearing1_, samples);
    const auto x2 = ExtractColumns(bearing2_, samples);
    Solver::Solve(x1, x2, R_relative_, models);
  }

  double Error(
    uint32_t sample,
    const Model &model) const
  {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    return ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }

  void Errors(
    const Model & model,
    std::vector<double> & vec_errors) const
  {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    vec_errors.resize(x1_.cols());
    for (uint32_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));
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
  Mat3 N1_, N2_;              // Normalization matrices
  double logalpha0_;          // A Contrario parameter
  Mat3 K1_, K2_;              // Intrinsic camera parameters
  Mat3 R_relative_;           // Relative rotation R = R2 · R1^T from IMU
};


/**
 * @brief Geometric filter using IMU-guided 2-point essential matrix estimation
 * 
 * This filter uses IMU rotation data to dramatically accelerate essential matrix
 * estimation, reducing RANSAC from 5-point (~2000 iterations) to 2-point (~100 iterations).
 */
struct GeometricFilter_EMatrix_AC_Imu
{
  GeometricFilter_EMatrix_AC_Imu(
    double dPrecision = std::numeric_limits<double>::infinity(),
    uint32_t iteration = 256,  // Much fewer iterations needed with 2-point
    const std::map<IndexT, Mat3> * imu_rotations = nullptr,
    size_t total_expected = 0
  ):
    m_dPrecision(dPrecision),
    m_stIteration(iteration),
    m_E(Mat3::Identity()),
    m_dPrecision_robust(std::numeric_limits<double>::infinity()),
    m_imu_rotations(imu_rotations),
    m_stats(std::make_shared<GeometricFilterStats>())
  {
    m_stats->total_expected = total_expected;
  }

  ~GeometricFilter_EMatrix_AC_Imu()
  {
    // Print stats only once (when the last shared copy is destroyed)
    if (m_stats.use_count() == 1)
    {
      m_stats->Print();
    }
  }

  /// Robust estimation of Essential matrix using 2-point IMU-guided RANSAC
  template<typename Regions_or_Features_ProviderT>
  bool Robust_estimation(
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<Regions_or_Features_ProviderT> & regions_provider,
    const Pair pairIndex,
    const matching::IndMatches & vec_PutativeMatches,
    matching::IndMatches & geometric_inliers)
  {
    geometric_inliers.clear();
    m_stats->total_pairs++;
    m_stats->PrintIfNewPercentage();

    const IndexT iIndex = pairIndex.first;
    const IndexT jIndex = pairIndex.second;

    //-- Check for IMU rotation data
    if (!m_imu_rotations)
    {
      OPENMVG_LOG_WARNING << "IMU rotations not available, cannot use IMU-guided solver";
      return false;
    }

    const sfm::View
      * view_I = sfm_data->views.at(iIndex).get(),
      * view_J = sfm_data->views.at(jIndex).get();

    // Get IMU rotations for both views
    const IndexT pose_i = view_I->id_pose;
    const IndexT pose_j = view_J->id_pose;

    auto it_i = m_imu_rotations->find(pose_i);
    auto it_j = m_imu_rotations->find(pose_j);

    if (it_i == m_imu_rotations->end() || it_j == m_imu_rotations->end())
    {
      return false;  // Missing IMU data for one or both views
    }

    const Mat3 & R1_wc = it_i->second;  // World-to-Camera for view 1
    const Mat3 & R2_wc = it_j->second;  // World-to-Camera for view 2

    //-- Check intrinsics
    const cameras::IntrinsicBase
      * cam_I = sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr,
      * cam_J = sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

    if (!cam_I || !cam_J)
      return false;

    if (!isPinhole(cam_I->getType()) || !isPinhole(cam_J->getType()))
      return false;

    const cameras::Pinhole_Intrinsic
      * ptrPinhole_I = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_I),
      * ptrPinhole_J = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam_J);

    //-- Get corresponding point regions
    Mat2X xI, xJ;
    MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

    //-- Robust estimation using 2-point IMU kernel
    using KernelType =
      ACKernelAdaptorEssentialImu<
        openMVG::fundamental::kernel::EpipolarDistanceError,
        Mat3>;

    KernelType kernel(
      xI, (*cam_I)(xI),
      sfm_data->GetViews().at(iIndex)->ui_width, sfm_data->GetViews().at(iIndex)->ui_height,
      xJ, (*cam_J)(xJ),
      sfm_data->GetViews().at(jIndex)->ui_width, sfm_data->GetViews().at(jIndex)->ui_height,
      ptrPinhole_I->K(), ptrPinhole_J->K(),
      R1_wc, R2_wc);

    // Robustly estimate the Essential matrix with A Contrario RANSAC
    const double upper_bound_precision = Square(m_dPrecision);
    std::vector<uint32_t> vec_inliers;
    const auto ACRansacOut =
      openMVG::robust::ACRANSAC(kernel, vec_inliers, m_stIteration, &m_E, upper_bound_precision);

    if (vec_inliers.size() <= KernelType::MINIMUM_SAMPLES * 2.5)
    {
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

  bool Geometry_guided_matching(
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair pairIndex,
    const double dDistanceRatio,
    matching::IndMatches & matches)
  {
    if (m_dPrecision_robust != std::numeric_limits<double>::infinity())
    {
      const IndexT iIndex = pairIndex.first;
      const IndexT jIndex = pairIndex.second;

      const sfm::View * view_I = sfm_data->views.at(iIndex).get();
      const sfm::View * view_J = sfm_data->views.at(jIndex).get();

      const cameras::IntrinsicBase * cam_I =
        sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr;
      const cameras::IntrinsicBase * cam_J =
        sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

      if (!cam_I || !cam_J)
        return false;

      if ( !isPinhole(cam_I->getType()) || !isPinhole(cam_J->getType()))
        return false;

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

  double m_dPrecision;
  uint32_t m_stIteration;
  Mat3 m_E;
  double m_dPrecision_robust;
  const std::map<IndexT, Mat3> * m_imu_rotations;
  std::shared_ptr<GeometricFilterStats> m_stats;
};

} //namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_IMU_HPP
