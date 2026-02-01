// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_IMU_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_E_AC_ROBUST_IMU_HPP

#include "openMVG/matching_image_collection/E_ACRobust_WithPriors.hpp"
#include "solver_essential_two_point_imu.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <map>

namespace openMVG {
namespace matching_image_collection {

/**
 * @brief Kernel adaptor for 2-point IMU-guided Essential matrix estimation using Angular Error
 * 
 * This kernel uses known IMU rotations to solve for translation in 2-DOF.
 * Uses Angular Error (radians) for resolution-independent verification.
 */
class ACKernelAdaptorEssentialImuAngular
{
public:
  using Solver = essential::kernel::TwoPointImuSolver;
  using Model = Mat3;
  using ErrorT = openMVG::AngularError;

  ACKernelAdaptorEssentialImuAngular(
    const Mat3X & bearing1,
    const Mat3X & bearing2,
    const Mat3 & R1_wc,
    const Mat3 & R2_wc
  ):bearing1_(bearing1),
    bearing2_(bearing2),
    R_relative_(R2_wc * R1_wc.transpose())
  {
    // A Contrario parameter for angular error (standard 0.5)
    logalpha0_ = 0.5; 
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const auto x1 = ExtractColumns(bearing1_, samples);
    const auto x2 = ExtractColumns(bearing2_, samples);
    Solver::Solve(x1, x2, R_relative_, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    return ErrorT::Error(model, bearing1_.col(sample), bearing2_.col(sample));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const {
    vec_errors.resize(bearing1_.cols());
    for (uint32_t sample = 0; sample < bearing1_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(model, bearing1_.col(sample), bearing2_.col(sample));
  }

  size_t NumSamples() const { return bearing1_.cols(); }
  double logalpha0() const { return logalpha0_; }
  double multError() const { return 1.0; }
  double unormalizeError(double val) const { return val; }

private:
  Mat3X bearing1_, bearing2_;
  Mat3 R_relative_;
  double logalpha0_;
};

/**
 * @brief Geometric filter using IMU-guided 2-point essential matrix estimation
 * Uses Angular Error (degrees) for resolution-independence.
 */
struct GeometricFilter_EMatrix_AC_Imu
{
  GeometricFilter_EMatrix_AC_Imu(
    double dPrecisionDeg = 0.2, // Default: 0.2 degrees angular error
    uint32_t iteration = 256,
    const std::map<IndexT, Mat3> * imu_rotations = nullptr,
    size_t total_expected = 0
  ):
    m_dPrecisionDeg(dPrecisionDeg),
    m_stIteration(iteration),
    m_E(Mat3::Identity()),
    m_dPrecision_robust(std::numeric_limits<double>::infinity()),
    m_imu_rotations(imu_rotations),
    m_stats(std::make_shared<GeometricFilterStats>())
  {
    m_stats->total_expected = total_expected;
  }

  ~GeometricFilter_EMatrix_AC_Imu() {
    if (m_stats.use_count() == 1) m_stats->Print();
  }

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
    
    if (!m_imu_rotations) return false;

    const IndexT iIndex = pairIndex.first;
    const IndexT jIndex = pairIndex.second;

    const sfm::View * vI = sfm_data->views.at(iIndex).get();
    const sfm::View * vJ = sfm_data->views.at(jIndex).get();

    auto it_i = m_imu_rotations->find(vI->id_pose);
    auto it_j = m_imu_rotations->find(vJ->id_pose);
    if (it_i == m_imu_rotations->end() || it_j == m_imu_rotations->end()) return false;

    const cameras::IntrinsicBase
      * cam_I = sfm_data->GetIntrinsics().at(vI->id_intrinsic).get(),
      * cam_J = sfm_data->GetIntrinsics().at(vJ->id_intrinsic).get();

    Mat2X xI, xJ;
    MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

    ACKernelAdaptorEssentialImuAngular kernel((*cam_I)(xI), (*cam_J)(xJ), it_i->second, it_j->second);

    const double upper_bound_precision = D2R(m_dPrecisionDeg);
    std::vector<uint32_t> vec_inliers;
    const auto ACRansacOut =
      openMVG::robust::ACRANSAC(kernel, vec_inliers, m_stIteration, &m_E, upper_bound_precision);

    if (vec_inliers.size() <= 6) return false;

    m_dPrecision_robust = ACRansacOut.first;
    for (const uint32_t & index : vec_inliers) geometric_inliers.push_back( vec_PutativeMatches[index] );
    return true;
  }

  double m_dPrecisionDeg;
  uint32_t m_stIteration;
  Mat3 m_E;
  double m_dPrecision_robust;
  const std::map<IndexT, Mat3> * m_imu_rotations;
  std::shared_ptr<GeometricFilterStats> m_stats;
};

} //namespace matching_image_collection
} // namespace openMVG

#endif
